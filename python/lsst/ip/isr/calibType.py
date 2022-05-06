# This file is part of ip_isr.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
import abc
import datetime
import logging
import os.path
import warnings
import yaml
import numpy as np

from astropy.table import Table
from astropy.io import fits

from lsst.daf.base import PropertyList
from lsst.utils.introspection import get_full_type_name
from lsst.utils import doImport


__all__ = ["IsrCalib", "IsrProvenance"]


class IsrCalib(abc.ABC):
    """Generic calibration type.

    Subclasses must implement the toDict, fromDict, toTable, fromTable
    methods that allow the calibration information to be converted
    from dictionaries and afw tables.  This will allow the calibration
    to be persisted using the base class read/write methods.

    The validate method is intended to provide a common way to check
    that the calibration is valid (internally consistent) and
    appropriate (usable with the intended data).  The apply method is
    intended to allow the calibration to be applied in a consistent
    manner.

    Parameters
    ----------
    camera : `lsst.afw.cameraGeom.Camera`, optional
        Camera to extract metadata from.
    detector : `lsst.afw.cameraGeom.Detector`, optional
        Detector to extract metadata from.
    log : `logging.Logger`, optional
        Log for messages.
    """
    _OBSTYPE = "generic"
    _SCHEMA = "NO SCHEMA"
    _VERSION = 0

    def __init__(self, camera=None, detector=None, log=None, **kwargs):
        self._instrument = None
        self._raftName = None
        self._slotName = None
        self._detectorName = None
        self._detectorSerial = None
        self._detectorId = None
        self._filter = None
        self._calibId = None
        self._metadata = PropertyList()
        self.setMetadata(PropertyList())
        self.calibInfoFromDict(kwargs)

        # Define the required attributes for this calibration.
        self.requiredAttributes = set(["_OBSTYPE", "_SCHEMA", "_VERSION"])
        self.requiredAttributes.update(["_instrument", "_raftName", "_slotName",
                                        "_detectorName", "_detectorSerial", "_detectorId",
                                        "_filter", "_calibId", "_metadata"])

        self.log = log if log else logging.getLogger(__name__)

        if detector:
            self.fromDetector(detector)
        self.updateMetadata(camera=camera, detector=detector)

    def __str__(self):
        return f"{self.__class__.__name__}(obstype={self._OBSTYPE}, detector={self._detectorName}, )"

    def __eq__(self, other):
        """Calibration equivalence.

        Running ``calib.log.setLevel(0)`` enables debug statements to
        identify problematic fields.
        """
        if not isinstance(other, self.__class__):
            self.log.debug("Incorrect class type: %s %s", self.__class__, other.__class__)
            return False

        for attr in self._requiredAttributes:
            attrSelf = getattr(self, attr)
            attrOther = getattr(other, attr)

            if isinstance(attrSelf, dict):
                # Dictionary of arrays.
                if attrSelf.keys() != attrOther.keys():
                    self.log.debug("Dict Key Failure: %s %s %s", attr, attrSelf.keys(), attrOther.keys())
                    return False
                for key in attrSelf:
                    try:
                        if not np.allclose(attrSelf[key], attrOther[key], equal_nan=True):
                            self.log.debug("Array Failure: %s %s %s", key, attrSelf[key], attrOther[key])
                            return False
                    except TypeError:
                        # If it is not something numpy can handle
                        # (it's not a number or array of numbers),
                        # then it needs to have its own equivalence
                        # operator.
                        if attrSelf[key] != attrOther[key]:
                            return False
            elif isinstance(attrSelf, np.ndarray):
                # Bare array.
                if not np.allclose(attrSelf, attrOther, equal_nan=True):
                    self.log.debug("Array Failure: %s %s %s", attr, attrSelf, attrOther)
                    return False
            elif type(attrSelf) != type(attrOther):
                if set([attrSelf, attrOther]) == set([None, ""]):
                    # Fits converts None to "", but None is not "".
                    continue
                self.log.debug("Type Failure: %s %s %s %s %s", attr, type(attrSelf), type(attrOther),
                               attrSelf, attrOther)
                return False
            else:
                if attrSelf != attrOther:
                    self.log.debug("Value Failure: %s %s %s", attr, attrSelf, attrOther)
                    return False

        return True

    @property
    def requiredAttributes(self):
        return self._requiredAttributes

    @requiredAttributes.setter
    def requiredAttributes(self, value):
        self._requiredAttributes = value

    def getMetadata(self):
        """Retrieve metadata associated with this calibration.

        Returns
        -------
        meta : `lsst.daf.base.PropertyList`
            Metadata. The returned `~lsst.daf.base.PropertyList` can be
            modified by the caller and the changes will be written to
            external files.
        """
        return self._metadata

    def setMetadata(self, metadata):
        """Store a copy of the supplied metadata with this calibration.

        Parameters
        ----------
        metadata : `lsst.daf.base.PropertyList`
            Metadata to associate with the calibration.  Will be copied and
            overwrite existing metadata.
        """
        if metadata is not None:
            self._metadata.update(metadata)

        # Ensure that we have the obs type required by calibration ingest
        self._metadata["OBSTYPE"] = self._OBSTYPE
        self._metadata[self._OBSTYPE + "_SCHEMA"] = self._SCHEMA
        self._metadata[self._OBSTYPE + "_VERSION"] = self._VERSION

        if isinstance(metadata, dict):
            self.calibInfoFromDict(metadata)
        elif isinstance(metadata, PropertyList):
            self.calibInfoFromDict(metadata.toDict())

    def updateMetadata(self, camera=None, detector=None, filterName=None,
                       setCalibId=False, setCalibInfo=False, setDate=False,
                       **kwargs):
        """Update metadata keywords with new values.

        Parameters
        ----------
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Reference camera to use to set _instrument field.
        detector : `lsst.afw.cameraGeom.Detector`, optional
            Reference detector to use to set _detector* fields.
        filterName : `str`, optional
            Filter name to assign to this calibration.
        setCalibId : `bool`, optional
            Construct the _calibId field from other fields.
        setCalibInfo : `bool`, optional
            Set calibration parameters from metadata.
        setDate : `bool`, optional
            Ensure the metadata CALIBDATE fields are set to the current
            datetime.
        kwargs : `dict` or `collections.abc.Mapping`, optional
            Set of key=value pairs to assign to the metadata.
        """
        mdOriginal = self.getMetadata()
        mdSupplemental = dict()

        for k, v in kwargs.items():
            if isinstance(v, fits.card.Undefined):
                kwargs[k] = None

        if setCalibInfo:
            self.calibInfoFromDict(kwargs)

        if camera:
            self._instrument = camera.getName()

        if detector:
            self._detectorName = detector.getName()
            self._detectorSerial = detector.getSerial()
            self._detectorId = detector.getId()
            if "_" in self._detectorName:
                (self._raftName, self._slotName) = self._detectorName.split("_")

        if filterName:
            # TOD0 DM-28093: I think this whole comment can go away, if we
            # always use physicalLabel everywhere in ip_isr.
            # If set via:
            # exposure.getInfo().getFilter().getName()
            # then this will hold the abstract filter.
            self._filter = filterName

        if setDate:
            date = datetime.datetime.now()
            mdSupplemental["CALIBDATE"] = date.isoformat()
            mdSupplemental["CALIB_CREATION_DATE"] = date.date().isoformat()
            mdSupplemental["CALIB_CREATION_TIME"] = date.time().isoformat()

        if setCalibId:
            values = []
            values.append(f"instrument={self._instrument}") if self._instrument else None
            values.append(f"raftName={self._raftName}") if self._raftName else None
            values.append(f"detectorName={self._detectorName}") if self._detectorName else None
            values.append(f"detector={self._detectorId}") if self._detectorId else None
            values.append(f"filter={self._filter}") if self._filter else None

            calibDate = mdOriginal.get("CALIBDATE", mdSupplemental.get("CALIBDATE", None))
            values.append(f"calibDate={calibDate}") if calibDate else None

            self._calibId = " ".join(values)

        self._metadata["INSTRUME"] = self._instrument if self._instrument else None
        self._metadata["RAFTNAME"] = self._raftName if self._raftName else None
        self._metadata["SLOTNAME"] = self._slotName if self._slotName else None
        self._metadata["DETECTOR"] = self._detectorId
        self._metadata["DET_NAME"] = self._detectorName if self._detectorName else None
        self._metadata["DET_SER"] = self._detectorSerial if self._detectorSerial else None
        self._metadata["FILTER"] = self._filter if self._filter else None
        self._metadata["CALIB_ID"] = self._calibId if self._calibId else None
        self._metadata["CALIBCLS"] = get_full_type_name(self)

        mdSupplemental.update(kwargs)
        mdOriginal.update(mdSupplemental)

    def calibInfoFromDict(self, dictionary):
        """Handle common keywords.

        This isn't an ideal solution, but until all calibrations
        expect to find everything in the metadata, they still need to
        search through dictionaries.

        Parameters
        ----------
        dictionary : `dict` or `lsst.daf.base.PropertyList`
            Source for the common keywords.

        Raises
        ------
        RuntimeError :
            Raised if the dictionary does not match the expected OBSTYPE.

        """

        def search(haystack, needles):
            """Search dictionary 'haystack' for an entry in 'needles'
            """
            test = [haystack.get(x) for x in needles]
            test = set([x for x in test if x is not None])
            if len(test) == 0:
                if "metadata" in haystack:
                    return search(haystack["metadata"], needles)
                else:
                    return None
            elif len(test) == 1:
                value = list(test)[0]
                if value == "":
                    return None
                else:
                    return value
            else:
                raise ValueError(f"Too many values found: {len(test)} {test} {needles}")

        if "metadata" in dictionary:
            metadata = dictionary["metadata"]

            if self._OBSTYPE != metadata["OBSTYPE"]:
                raise RuntimeError(f"Incorrect calibration supplied.  Expected {self._OBSTYPE}, "
                                   f"found {metadata['OBSTYPE']}")

        self._instrument = search(dictionary, ["INSTRUME", "instrument"])
        self._raftName = search(dictionary, ["RAFTNAME"])
        self._slotName = search(dictionary, ["SLOTNAME"])
        self._detectorId = search(dictionary, ["DETECTOR", "detectorId"])
        self._detectorName = search(dictionary, ["DET_NAME", "DETECTOR_NAME", "detectorName"])
        self._detectorSerial = search(dictionary, ["DET_SER", "DETECTOR_SERIAL", "detectorSerial"])
        self._filter = search(dictionary, ["FILTER", "filterName"])
        self._calibId = search(dictionary, ["CALIB_ID"])

    @classmethod
    def determineCalibClass(cls, metadata, message):
        """Attempt to find calibration class in metadata.

        Parameters
        ----------
        metadata : `dict` or `lsst.daf.base.PropertyList`
            Metadata possibly containing a calibration class entry.
        message : `str`
            Message to include in any errors.

        Returns
        -------
        calibClass : `object`
            The class to use to read the file contents.  Should be an
            `lsst.ip.isr.IsrCalib` subclass.

        Raises
        ------
        ValueError :
            Raised if the resulting calibClass is the base
            `lsst.ip.isr.IsrClass` (which does not implement the
            content methods).
        """
        calibClassName = metadata.get("CALIBCLS")
        calibClass = doImport(calibClassName) if calibClassName is not None else cls
        if calibClass is IsrCalib:
            raise ValueError(f"Cannot use base class to read calibration data: {message}")
        return calibClass

    @classmethod
    def readText(cls, filename, **kwargs):
        """Read calibration representation from a yaml/ecsv file.

        Parameters
        ----------
        filename : `str`
            Name of the file containing the calibration definition.
        kwargs : `dict` or collections.abc.Mapping`, optional
            Set of key=value pairs to pass to the ``fromDict`` or
            ``fromTable`` methods.

        Returns
        -------
        calib : `~lsst.ip.isr.IsrCalibType`
            Calibration class.

        Raises
        ------
        RuntimeError :
            Raised if the filename does not end in ".ecsv" or ".yaml".
        """
        if filename.endswith((".ecsv", ".ECSV")):
            data = Table.read(filename, format="ascii.ecsv")
            calibClass = cls.determineCalibClass(data.meta, "readText/ECSV")
            return calibClass.fromTable([data], **kwargs)
        elif filename.endswith((".yaml", ".YAML")):
            with open(filename, "r") as f:
                data = yaml.load(f, Loader=yaml.CLoader)
            calibClass = cls.determineCalibClass(data["metadata"], "readText/YAML")
            return calibClass.fromDict(data, **kwargs)
        else:
            raise RuntimeError(f"Unknown filename extension: {filename}")

    def writeText(self, filename, format="auto"):
        """Write the calibration data to a text file.

        Parameters
        ----------
        filename : `str`
            Name of the file to write.
        format : `str`
            Format to write the file as.  Supported values are:
                ``"auto"`` : Determine filetype from filename.
                ``"yaml"`` : Write as yaml.
                ``"ecsv"`` : Write as ecsv.
        Returns
        -------
        used : `str`
            The name of the file used to write the data.  This may
            differ from the input if the format is explicitly chosen.

        Raises
        ------
        RuntimeError :
            Raised if filename does not end in a known extension, or
            if all information cannot be written.

        Notes
        -----
        The file is written to YAML/ECSV format and will include any
        associated metadata.
        """
        if format == "yaml" or (format == "auto" and filename.lower().endswith((".yaml", ".YAML"))):
            outDict = self.toDict()
            path, ext = os.path.splitext(filename)
            filename = path + ".yaml"
            with open(filename, "w") as f:
                yaml.dump(outDict, f)
        elif format == "ecsv" or (format == "auto" and filename.lower().endswith((".ecsv", ".ECSV"))):
            tableList = self.toTable()
            if len(tableList) > 1:
                # ECSV doesn't support multiple tables per file, so we
                # can only write the first table.
                raise RuntimeError(f"Unable to persist {len(tableList)}tables in ECSV format.")

            table = tableList[0]
            path, ext = os.path.splitext(filename)
            filename = path + ".ecsv"
            table.write(filename, format="ascii.ecsv")
        else:
            raise RuntimeError(f"Attempt to write to a file {filename} "
                               "that does not end in '.yaml' or '.ecsv'")

        return filename

    @classmethod
    def readFits(cls, filename, **kwargs):
        """Read calibration data from a FITS file.

        Parameters
        ----------
        filename : `str`
            Filename to read data from.
        kwargs : `dict` or collections.abc.Mapping`, optional
            Set of key=value pairs to pass to the ``fromTable``
            method.

        Returns
        -------
        calib : `lsst.ip.isr.IsrCalib`
            Calibration contained within the file.
        """
        tableList = []
        tableList.append(Table.read(filename, hdu=1))
        extNum = 2  # Fits indices start at 1, we've read one already.
        keepTrying = True

        while keepTrying:
            with warnings.catch_warnings():
                warnings.simplefilter("error")
                try:
                    newTable = Table.read(filename, hdu=extNum)
                    tableList.append(newTable)
                    extNum += 1
                except Exception:
                    keepTrying = False

        for table in tableList:
            for k, v in table.meta.items():
                if isinstance(v, fits.card.Undefined):
                    table.meta[k] = None

        calibClass = cls.determineCalibClass(tableList[0].meta, "readFits")
        return calibClass.fromTable(tableList, **kwargs)

    def writeFits(self, filename):
        """Write calibration data to a FITS file.

        Parameters
        ----------
        filename : `str`
            Filename to write data to.

        Returns
        -------
        used : `str`
            The name of the file used to write the data.

        """
        tableList = self.toTable()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=Warning, module="astropy.io")
            astropyList = [fits.table_to_hdu(table) for table in tableList]
            astropyList.insert(0, fits.PrimaryHDU())

            writer = fits.HDUList(astropyList)
            writer.writeto(filename, overwrite=True)
        return filename

    def fromDetector(self, detector):
        """Modify the calibration parameters to match the supplied detector.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            Detector to use to set parameters from.

        Raises
        ------
        NotImplementedError
            This needs to be implemented by subclasses for each
            calibration type.
        """
        raise NotImplementedError("Must be implemented by subclass.")

    @classmethod
    def fromDict(cls, dictionary, **kwargs):
        """Construct a calibration from a dictionary of properties.

        Must be implemented by the specific calibration subclasses.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties.
        kwargs : `dict` or collections.abc.Mapping`, optional
            Set of key=value options.

        Returns
        ------
        calib : `lsst.ip.isr.CalibType`
            Constructed calibration.

        Raises
        ------
        NotImplementedError :
            Raised if not implemented.
        """
        raise NotImplementedError("Must be implemented by subclass.")

    def toDict(self):
        """Return a dictionary containing the calibration properties.

        The dictionary should be able to be round-tripped through
        `fromDict`.

        Returns
        -------
        dictionary : `dict`
            Dictionary of properties.

        Raises
        ------
        NotImplementedError :
            Raised if not implemented.
        """
        raise NotImplementedError("Must be implemented by subclass.")

    @classmethod
    def fromTable(cls, tableList, **kwargs):
        """Construct a calibration from a dictionary of properties.

        Must be implemented by the specific calibration subclasses.

        Parameters
        ----------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables of properties.
        kwargs : `dict` or collections.abc.Mapping`, optional
            Set of key=value options.

        Returns
        ------
        calib : `lsst.ip.isr.CalibType`
            Constructed calibration.

        Raises
        ------
        NotImplementedError :
            Raised if not implemented.
        """
        raise NotImplementedError("Must be implemented by subclass.")

    def toTable(self):
        """Return a list of tables containing the calibration properties.

        The table list should be able to be round-tripped through
        `fromDict`.

        Returns
        -------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables of properties.

        Raises
        ------
        NotImplementedError :
            Raised if not implemented.
        """
        raise NotImplementedError("Must be implemented by subclass.")

    def validate(self, other=None):
        """Validate that this calibration is defined and can be used.

        Parameters
        ----------
        other : `object`, optional
            Thing to validate against.

        Returns
        -------
        valid : `bool`
            Returns true if the calibration is valid and appropriate.
        """
        return False

    def apply(self, target):
        """Method to apply the calibration to the target object.

        Parameters
        ----------
        target : `object`
            Thing to validate against.

        Returns
        -------
        valid : `bool`
            Returns true if the calibration was applied correctly.

        Raises
        ------
        NotImplementedError :
            Raised if not implemented.
        """
        raise NotImplementedError("Must be implemented by subclass.")


class IsrProvenance(IsrCalib):
    """Class for the provenance of data used to construct calibration.

    Provenance is not really a calibration, but we would like to
    record this when constructing the calibration, and it provides an
    example of the base calibration class.

    Parameters
    ----------
    instrument : `str`, optional
        Name of the instrument the data was taken with.
    calibType : `str`, optional
        Type of calibration this provenance was generated for.
    detectorName : `str`, optional
        Name of the detector this calibration is for.
    detectorSerial : `str`, optional
        Identifier for the detector.

    """
    _OBSTYPE = "IsrProvenance"

    def __init__(self, calibType="unknown",
                 **kwargs):
        self.calibType = calibType
        self.dimensions = set()
        self.dataIdList = list()

        super().__init__(**kwargs)

        self.requiredAttributes.update(["calibType", "dimensions", "dataIdList"])

    def __str__(self):
        return f"{self.__class__.__name__}(obstype={self._OBSTYPE}, calibType={self.calibType}, )"

    def __eq__(self, other):
        return super().__eq__(other)

    def updateMetadata(self, setDate=False, **kwargs):
        """Update calibration metadata.

        Parameters
        ----------
        setDate : `bool, optional
            Update the CALIBDATE fields in the metadata to the current
            time. Defaults to False.
        kwargs : `dict` or `collections.abc.Mapping`, optional
            Other keyword parameters to set in the metadata.
        """
        kwargs["calibType"] = self.calibType
        super().updateMetadata(setDate=setDate, **kwargs)

    def fromDataIds(self, dataIdList):
        """Update provenance from dataId List.

        Parameters
        ----------
        dataIdList : `list` [`lsst.daf.butler.DataId`]
            List of dataIds used in generating this calibration.
        """
        for dataId in dataIdList:
            for key in dataId:
                if key not in self.dimensions:
                    self.dimensions.add(key)
            self.dataIdList.append(dataId)

    @classmethod
    def fromTable(cls, tableList):
        """Construct provenance from table list.

        Parameters
        ----------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables to construct the provenance from.

        Returns
        -------
        provenance : `lsst.ip.isr.IsrProvenance`
            The provenance defined in the tables.
        """
        table = tableList[0]
        metadata = table.meta
        inDict = dict()
        inDict["metadata"] = metadata
        inDict["calibType"] = metadata["calibType"]
        inDict["dimensions"] = set()
        inDict["dataIdList"] = list()

        schema = dict()
        for colName in table.columns:
            schema[colName.lower()] = colName
            inDict["dimensions"].add(colName.lower())
        inDict["dimensions"] = sorted(inDict["dimensions"])

        for row in table:
            entry = dict()
            for dim in sorted(inDict["dimensions"]):
                entry[dim] = row[schema[dim]]
            inDict["dataIdList"].append(entry)

        return cls.fromDict(inDict)

    @classmethod
    def fromDict(cls, dictionary):
        """Construct provenance from a dictionary.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of provenance parameters.

        Returns
        -------
        provenance : `lsst.ip.isr.IsrProvenance`
            The provenance defined in the tables.
        """
        calib = cls()
        if calib._OBSTYPE != dictionary["metadata"]["OBSTYPE"]:
            raise RuntimeError(f"Incorrect calibration supplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['metadata']['OBSTYPE']}")
        calib.updateMetadata(setDate=False, setCalibInfo=True, **dictionary["metadata"])

        # These properties should be in the metadata, but occasionally
        # are found in the dictionary itself.  Check both places,
        # ending with `None` if neither contains the information.
        calib.calibType = dictionary["calibType"]
        calib.dimensions = set(dictionary["dimensions"])
        calib.dataIdList = dictionary["dataIdList"]

        calib.updateMetadata()
        return calib

    def toDict(self):
        """Return a dictionary containing the provenance information.

        Returns
        -------
        dictionary : `dict`
            Dictionary of provenance.
        """
        self.updateMetadata()

        outDict = {}

        metadata = self.getMetadata()
        outDict["metadata"] = metadata
        outDict["detectorName"] = self._detectorName
        outDict["detectorSerial"] = self._detectorSerial
        outDict["detectorId"] = self._detectorId
        outDict["instrument"] = self._instrument
        outDict["calibType"] = self.calibType
        outDict["dimensions"] = list(self.dimensions)
        outDict["dataIdList"] = self.dataIdList

        return outDict

    def toTable(self):
        """Return a list of tables containing the provenance.

        This seems inefficient and slow, so this may not be the best
        way to store the data.

        Returns
        -------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables containing the provenance information

        """
        tableList = []
        self.updateMetadata(setDate=True, setCalibInfo=True)

        catalog = Table(rows=self.dataIdList,
                        names=self.dimensions)
        filteredMetadata = {k: v for k, v in self.getMetadata().toDict().items() if v is not None}
        catalog.meta = filteredMetadata
        tableList.append(catalog)
        return tableList
