#
# LSST Data Management System
# Copyright 2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import abc
import copy
import datetime

import numpy as np
import yaml

import lsst.afw.table as afwTable
from lsst.daf.base import PropertyList
from lsst.pipe.base import Struct
from lsst.geom import Box2I, Point2I, Extent2I
from .applyLookupTable import applyLookupTable

__all__ = ["Linearizer",
           "LinearizeBase", "LinearizeLookupTable", "LinearizeSquared",
           "LinearizeProportional", "LinearizePolynomial", "LinearizeNone"]


class Linearizer(abc.ABC):
    """Parameter set for linearization.

    These parameters are included in cameraGeom.Amplifier, but
    should be accessible externally to allow for testing.

    Parameters
    ----------
    table : `numpy.array`, optional
        Lookup table; a 2-dimensional array of floats:
            - one row for each row index (value of coef[0] in the amplifier)
            - one column for each image value
        To avoid copying the table the last index should vary fastest
        (numpy default "C" order)
    detector : `lsst.afw.cameraGeom.Detector`
        Detector object
    override : `bool`, optional
        Override the parameters defined in the detector/amplifier.
    log : `lsst.log.Log`, optional
        Logger to handle messages.

    Raises
    ------
    RuntimeError :
        Raised if the supplied table is not 2D, or if the table has fewer
        columns than rows (indicating that the indices are swapped).
    """

    _OBSTYPE = "linearizer"
    """The dataset type name used for this class"""

    def __init__(self, table=None, detector=None, override=False, log=None):
        self._detectorName = None
        self._detectorSerial = None
        self._detectorId = None
        self._metadata = PropertyList()

        self.linearityCoeffs = dict()
        self.linearityType = dict()
        self.linearityThreshold = dict()
        self.linearityMaximum = dict()
        self.linearityUnits = dict()
        self.linearityBBox = dict()

        self.override = override
        self.populated = False
        self.log = log

        self.tableData = None
        if table is not None:
            if len(table.shape) != 2:
                raise RuntimeError("table shape = %s; must have two dimensions" % (table.shape,))
            if table.shape[1] < table.shape[0]:
                raise RuntimeError("table shape = %s; indices are switched" % (table.shape,))
            self.tableData = np.array(table, order="C")

        if detector:
            self.fromDetector(detector)

    def __call__(self, exposure):
        """Apply linearity, setting parameters if necessary.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to correct.

        Returns
        -------
        output : `lsst.pipe.base.Struct`
            Linearization results:
            ``"numAmps"``
                Number of amplifiers considered.
            ``"numLinearized"``
                Number of amplifiers linearized.
        """

    def fromDetector(self, detector):
        """Read linearity parameters from a detector.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.detector`
            Input detector with parameters to use.
        """
        self._detectorName = detector.getName()
        self._detectorSerial = detector.getSerial()
        self._detectorId = detector.getId()
        self.populated = True

        # Do not translate Threshold, Maximum, Units.
        for amp in detector.getAmplifiers():
            ampName = amp.getName()
            self.linearityCoeffs[ampName] = amp.getLinearityCoeffs()
            self.linearityType[ampName] = amp.getLinearityType()
            self.linearityBBox[ampName] = amp.getBBox()

    def fromYaml(self, yamlObject):
        """Read linearity parameters from a dict.

        Parameters
        ----------
        yamlObject : `dict`
            Dictionary containing detector and amplifier information.
        """
        self.setMetadata(metadata=yamlObject.get('metadata', None))
        self._detectorName = yamlObject['detectorName']
        self._detectorSerial = yamlObject['detectorSerial']
        self._detectorId = yamlObject['detectorId']
        self.populated = True
        self.override = True

        for ampName in yamlObject['amplifiers']:
            amp = yamlObject['amplifiers'][ampName]
            self.linearityCoeffs[ampName] = np.array(amp.get('linearityCoeffs', None), dtype=np.float64)
            self.linearityType[ampName] = amp.get('linearityType', 'None')
            self.linearityBBox[ampName] = amp.get('linearityBBox', None)

        if self.tableData is None:
            self.tableData = yamlObject.get('tableData', None)
            if self.tableData:
                self.tableData = np.array(self.tableData)

        return self

    def toDict(self):
        """Return linearity parameters as a dict.

        Returns
        -------
        outDict : `dict`:
        """
        # metadata copied from defects code
        now = datetime.datetime.utcnow()
        self.updateMetadata(date=now)

        outDict = {'metadata': self.getMetadata(),
                   'detectorName': self._detectorName,
                   'detectorSerial': self._detectorSerial,
                   'detectorId': self._detectorId,
                   'hasTable': self.tableData is not None,
                   'amplifiers': dict()}
        for ampName in self.linearityType:
            outDict['amplifiers'][ampName] = {'linearityType': self.linearityType[ampName],
                                              'linearityCoeffs': self.linearityCoeffs[ampName],
                                              'linearityBBox': self.linearityBBox[ampName]}
        if self.tableData is not None:
            outDict['tableData'] = self.tableData.tolist()

        return outDict

    @classmethod
    def readText(cls, filename):
        """Read linearity from text file.

        Parameters
        ----------
        filename : `str`
            Name of the file containing the linearity definition.
        Returns
        -------
        linearity : `~lsst.ip.isr.linearize.Linearizer``
            Linearity parameters.
        """
        data = ''
        with open(filename, 'r') as f:
            data = yaml.load(f, Loader=yaml.CLoader)
        return cls().fromYaml(data)

    def writeText(self, filename):
        """Write the linearity model to a text file.

        Parameters
        ----------
        filename : `str`
            Name of the file to write.

        Returns
        -------
        used : `str`
            The name of the file used to write the data.

        Raises
        ------
        RuntimeError :
            Raised if filename does not end in ".yaml".

        Notes
        -----
        The file is written to YAML format and will include any metadata
        associated with the `Linearity`.
        """
        outDict = self.toDict()
        if filename.lower().endswith((".yaml")):
            with open(filename, 'w') as f:
                yaml.dump(outDict, f)
        else:
            raise RuntimeError(f"Attempt to write to a file {filename} that does not end in '.yaml'")

        return filename

    @classmethod
    def fromTable(cls, table, tableExtTwo=None):
        """Read linearity from a FITS file.

        Parameters
        ----------
        table : `lsst.afw.table`
            afwTable read from input file name.
        tableExtTwo: `lsst.afw.table`, optional
            afwTable read from second extension of input file name

        Returns
        -------
        linearity : `~lsst.ip.isr.linearize.Linearizer``
            Linearity parameters.

        Notes
        -----
        The method reads a FITS file with 1 or 2 extensions. The metadata is read from the header of
        extension 1, which must exist.  Then the table is loaded, and  the ['AMPLIFIER_NAME', 'TYPE',
        'COEFFS', 'BBOX_X0', 'BBOX_Y0', 'BBOX_DX', 'BBOX_DY'] columns are read and used to
        set each dictionary by looping over rows.
        Eextension 2 is then attempted to read in the try block (which only exists for lookup tables).
        It has a column named 'LOOKUP_VALUES' that contains a vector of the lookup entries in each row.
        """
        metadata = table.getMetadata()
        schema = table.getSchema()

        linDict = dict()
        linDict['metadata'] = metadata
        linDict['detectorId'] = metadata['DETECTOR']
        linDict['detectorName'] = metadata['DETECTOR_NAME']
        try:
            linDict['detectorSerial'] = metadata['DETECTOR_SERIAL']
        except Exception:
            linDict['detectorSerial'] = 'NOT SET'
        linDict['amplifiers'] = dict()

        # Preselect the keys
        ampNameKey = schema['AMPLIFIER_NAME'].asKey()
        typeKey = schema['TYPE'].asKey()
        coeffsKey = schema['COEFFS'].asKey()
        x0Key = schema['BBOX_X0'].asKey()
        y0Key = schema['BBOX_Y0'].asKey()
        dxKey = schema['BBOX_DX'].asKey()
        dyKey = schema['BBOX_DY'].asKey()

        for record in table:
            ampName = record[ampNameKey]
            ampDict = dict()
            ampDict['linearityType'] = record[typeKey]
            ampDict['linearityCoeffs'] = record[coeffsKey]
            ampDict['linearityBBox'] = Box2I(Point2I(record[x0Key], record[y0Key]),
                                             Extent2I(record[dxKey], record[dyKey]))

            linDict['amplifiers'][ampName] = ampDict

        if tableExtTwo is not None:
            lookupValuesKey = 'LOOKUP_VALUES'
            linDict["tableData"] = [record[lookupValuesKey] for record in tableExtTwo]

        return cls().fromYaml(linDict)

    @classmethod
    def readFits(cls, filename):
        """Read linearity from a FITS file.

        Parameters
        ----------
        filename : `str`
            Name of the file containing the linearity definition.
        Returns
        -------
        linearity : `~lsst.ip.isr.linearize.Linearizer``
            Linearity parameters.

        Notes
        -----
        This method and `fromTable` read a FITS file with 1 or 2 extensions. The metadata is read from the
        header of extension 1, which must exist.  Then the table is loaded, and the ['AMPLIFIER_NAME',
        'TYPE', 'COEFFS', 'BBOX_X0', 'BBOX_Y0', 'BBOX_DX', 'BBOX_DY'] columns are read and used to
        set each dictionary by looping over rows.
        Extension 2 is then attempted to read in the try block (which only exists for lookup tables).
        It has a column named 'LOOKUP_VALUES' that contains a vector of the lookup entries in each row.
        """
        table = afwTable.BaseCatalog.readFits(filename)
        tableExtTwo = None
        try:
            tableExtTwo = afwTable.BaseCatalog.readFits(filename, 2)
        except Exception:
            pass
        return cls().fromTable(table, tableExtTwo=tableExtTwo)

    def toAmpTable(self, metadata):
        """Produce linearity catalog

        Parameters
        ----------
        metadata : `lsst.daf.base.PropertyList`
            Linearizer metadata

        Returns
        -------
        catalog : `lsst.afw.table.BaseCatalog`
            Catalog to write
        """
        metadata["LINEARITY_SCHEMA"] = "Linearity table"
        metadata["LINEARITY_VERSION"] = 1

        # Now pack it into a fits table.
        length = max([len(self.linearityCoeffs[x]) for x in self.linearityCoeffs.keys()])

        schema = afwTable.Schema()
        names = schema.addField("AMPLIFIER_NAME", type="String", size=16, doc="linearity amplifier name")
        types = schema.addField("TYPE", type="String", size=16, doc="linearity type names")
        coeffs = schema.addField("COEFFS", type="ArrayD", size=length, doc="linearity coefficients")
        boxX = schema.addField("BBOX_X0", type="I", doc="linearity bbox minimum x")
        boxY = schema.addField("BBOX_Y0", type="I", doc="linearity bbox minimum y")
        boxDx = schema.addField("BBOX_DX", type="I", doc="linearity bbox x dimension")
        boxDy = schema.addField("BBOX_DY", type="I", doc="linearity bbox y dimension")

        catalog = afwTable.BaseCatalog(schema)
        catalog.resize(len(self.linearityCoeffs.keys()))

        for ii, ampName in enumerate(self.linearityType):
            catalog[ii][names] = ampName
            catalog[ii][types] = self.linearityType[ampName]
            catalog[ii][coeffs] = np.array(self.linearityCoeffs[ampName], dtype=float)

            bbox = self.linearityBBox[ampName]
            catalog[ii][boxX], catalog[ii][boxY] = bbox.getMin()
            catalog[ii][boxDx], catalog[ii][boxDy] = bbox.getDimensions()
        catalog.setMetadata(metadata)

        return catalog

    def toTableDataTable(self, metadata):
        """Produce linearity catalog from table data

        Parameters
        ----------
        metadata : `lsst.daf.base.PropertyList`
            Linearizer metadata

        Returns
        -------
        catalog : `lsst.afw.table.BaseCatalog`
            Catalog to write
        """

        schema = afwTable.Schema()
        dimensions = self.tableData.shape
        lut = schema.addField("LOOKUP_VALUES", type='ArrayI', size=dimensions[1],
                              doc="linearity lookup data")
        catalog = afwTable.BaseCatalog(schema)
        catalog.resize(dimensions[0])

        for ii in range(dimensions[0]):
            catalog[ii][lut] = np.array(self.tableData[ii], dtype=np.intc)

        metadata["LINEARITY_LOOKUP"] = True
        catalog.setMetadata(metadata)

        return catalog

    def writeFits(self, filename):
        """Write the linearity model to a FITS file.

        Parameters
        ----------
        filename : `str`
            Name of the file to write.

        Notes
        -----
        The file is written to YAML format and will include any metadata
        associated with the `Linearity`.
        """
        now = datetime.datetime.utcnow()
        self.updateMetadata(date=now)
        metadata = copy.copy(self.getMetadata())
        catalog = self.toAmpTable(metadata)
        catalog.writeFits(filename)

        if self.tableData is not None:
            catalog = self.toTableDataTable(metadata)
            catalog.writeFits(filename, "a")

        return

    def getMetadata(self):
        """Retrieve metadata associated with this `Linearizer`.

        Returns
        -------
        meta : `lsst.daf.base.PropertyList`
            Metadata. The returned `~lsst.daf.base.PropertyList` can be
            modified by the caller and the changes will be written to
            external files.
        """
        return self._metadata

    def setMetadata(self, metadata=None):
        """Store a copy of the supplied metadata with the `Linearizer`.

        Parameters
        ----------
        metadata : `lsst.daf.base.PropertyList`, optional
            Metadata to associate with the linearizer.  Will be copied and
            overwrite existing metadata.  If not supplied the existing
            metadata will be reset.
        """
        if metadata is None:
            self._metadata = PropertyList()
        else:
            self._metadata = copy.copy(metadata)

        # Ensure that we have the obs type required by calibration ingest
        self._metadata["OBSTYPE"] = self._OBSTYPE

    def updateMetadata(self, date=None, detectorId=None, detectorName=None, instrumentName=None, calibId=None,
                       serial=None):
        """Update metadata keywords with new values.

        Parameters
        ----------
        date : `datetime.datetime`, optional
        detectorId : `int`, optional
        detectorName: `str`, optional
        instrumentName : `str`, optional
        calibId: `str`, optional
        serial: detector serial, `str`, optional

        """
        mdOriginal = self.getMetadata()
        mdSupplemental = dict()

        if date:
            mdSupplemental['CALIBDATE'] = date.isoformat()
            mdSupplemental['CALIB_CREATION_DATE'] = date.date().isoformat(),
            mdSupplemental['CALIB_CREATION_TIME'] = date.time().isoformat(),
        if detectorId:
            mdSupplemental['DETECTOR'] = f"{detectorId}"
        if detectorName:
            mdSupplemental['DETECTOR_NAME'] = detectorName
        if instrumentName:
            mdSupplemental['INSTRUME'] = instrumentName
        if calibId:
            mdSupplemental['CALIB_ID'] = calibId
        if serial:
            mdSupplemental['DETECTOR_SERIAL'] = serial

        mdOriginal.update(mdSupplemental)

    def getLinearityTypeByName(self, linearityTypeName):
        """Determine the linearity class to use from the type name.

        Parameters
        ----------
        linearityTypeName : str
            String name of the linearity type that is needed.

        Returns
        -------
        linearityType : `~lsst.ip.isr.linearize.LinearizeBase`
            The appropriate linearity class to use.  If no matching class
            is found, `None` is returned.
        """
        for t in [LinearizeLookupTable,
                  LinearizeSquared,
                  LinearizePolynomial,
                  LinearizeProportional,
                  LinearizeNone]:
            if t.LinearityType == linearityTypeName:
                return t
        return None

    def validate(self, detector=None, amplifier=None):
        """Validate linearity for a detector/amplifier.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`, optional
            Detector to validate, along with its amplifiers.
        amplifier : `lsst.afw.cameraGeom.Amplifier`, optional
            Single amplifier to validate.

        Raises
        ------
        RuntimeError :
            Raised if there is a mismatch in linearity parameters, and
            the cameraGeom parameters are not being overridden.
        """
        amplifiersToCheck = []
        if detector:
            if self._detectorName != detector.getName():
                raise RuntimeError("Detector names don't match: %s != %s" %
                                   (self._detectorName, detector.getName()))
            if int(self._detectorId) != int(detector.getId()):
                raise RuntimeError("Detector IDs don't match: %s != %s" %
                                   (int(self._detectorId), int(detector.getId())))
            if self._detectorSerial != detector.getSerial():
                raise RuntimeError("Detector serial numbers don't match: %s != %s" %
                                   (self._detectorSerial, detector.getSerial()))
            if len(detector.getAmplifiers()) != len(self.linearityCoeffs.keys()):
                raise RuntimeError("Detector number of amps = %s does not match saved value %s" %
                                   (len(detector.getAmplifiers()),
                                    len(self.linearityCoeffs.keys())))
            amplifiersToCheck.extend(detector.getAmplifiers())

        if amplifier:
            amplifiersToCheck.extend(amplifier)

        for amp in amplifiersToCheck:
            ampName = amp.getName()
            if ampName not in self.linearityCoeffs.keys():
                raise RuntimeError("Amplifier %s is not in linearity data" %
                                   (ampName, ))
            if amp.getLinearityType() != self.linearityType[ampName]:
                if self.override:
                    self.log.warn("Overriding amplifier defined linearityType (%s) for %s",
                                  self.linearityType[ampName], ampName)
                else:
                    raise RuntimeError("Amplifier %s type %s does not match saved value %s" %
                                       (ampName, amp.getLinearityType(), self.linearityType[ampName]))
            if (amp.getLinearityCoeffs().shape != self.linearityCoeffs[ampName].shape or not
                    np.allclose(amp.getLinearityCoeffs(), self.linearityCoeffs[ampName], equal_nan=True)):
                if self.override:
                    self.log.warn("Overriding amplifier defined linearityCoeffs (%s) for %s",
                                  self.linearityCoeffs[ampName], ampName)
                else:
                    raise RuntimeError("Amplifier %s coeffs %s does not match saved value %s" %
                                       (ampName, amp.getLinearityCoeffs(), self.linearityCoeffs[ampName]))

    def applyLinearity(self, image, detector=None, log=None):
        """Apply the linearity to an image.

        If the linearity parameters are populated, use those,
        otherwise use the values from the detector.

        Parameters
        ----------
        image : `~lsst.afw.image.image`
            Image to correct.
        detector : `~lsst.afw.cameraGeom.detector`
            Detector to use for linearity parameters if not already
            populated.
        log : `~lsst.log.Log`, optional
            Log object to use for logging.
        """
        if log is None:
            log = self.log
        if detector and not self.populated:
            self.fromDetector(detector)

        self.validate(detector)

        numAmps = 0
        numLinearized = 0
        numOutOfRange = 0
        for ampName in self.linearityType.keys():
            linearizer = self.getLinearityTypeByName(self.linearityType[ampName])
            numAmps += 1
            if linearizer is not None:
                ampView = image.Factory(image, self.linearityBBox[ampName])
                success, outOfRange = linearizer()(ampView, **{'coeffs': self.linearityCoeffs[ampName],
                                                               'table': self.tableData,
                                                               'log': self.log})
                numOutOfRange += outOfRange
                if success:
                    numLinearized += 1
                elif log is not None:
                    log.warn("Amplifier %s did not linearize.",
                             ampName)
        return Struct(
            numAmps=numAmps,
            numLinearized=numLinearized,
            numOutOfRange=numOutOfRange
        )


class LinearizeBase(metaclass=abc.ABCMeta):
    """Abstract base class functor for correcting non-linearity.

    Subclasses must define __call__ and set class variable
    LinearityType to a string that will be used for linearity type in
    the cameraGeom.Amplifier.linearityType field.

    All linearity corrections should be defined in terms of an
    additive correction, such that:

    corrected_value = uncorrected_value + f(uncorrected_value)
    """
    LinearityType = None  # linearity type, a string used for AmpInfoCatalogs

    @abc.abstractmethod
    def __call__(self, image, **kwargs):
        """Correct non-linearity.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Image to be corrected
        kwargs : `dict`
            Dictionary of parameter keywords:
            ``"coeffs"``
                Coefficient vector (`list` or `numpy.array`).
            ``"table"``
                Lookup table data (`numpy.array`).
            ``"log"``
                Logger to handle messages (`lsst.log.Log`).

        Returns
        -------
        output : `bool`
            If true, a correction was applied successfully.

        Raises
        ------
        RuntimeError:
            Raised if the linearity type listed in the
            detector does not match the class type.
        """
        pass


class LinearizeLookupTable(LinearizeBase):
    """Correct non-linearity with a persisted lookup table.

    The lookup table consists of entries such that given
    "coefficients" c0, c1:

    for each i,j of image:
        rowInd = int(c0)
        colInd = int(c1 + uncorrImage[i,j])
        corrImage[i,j] = uncorrImage[i,j] + table[rowInd, colInd]

    - c0: row index; used to identify which row of the table to use
            (typically one per amplifier, though one can have multiple
            amplifiers use the same table)
    - c1: column index offset; added to the uncorrected image value
            before truncation; this supports tables that can handle
            negative image values; also, if the c1 ends with .5 then
            the nearest index is used instead of truncating to the
            next smaller index
    """
    LinearityType = "LookupTable"

    def __call__(self, image, **kwargs):
        """Correct for non-linearity.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Image to be corrected
        kwargs : `dict`
            Dictionary of parameter keywords:
            ``"coeffs"``
                Columnation vector (`list` or `numpy.array`).
            ``"table"``
                Lookup table data (`numpy.array`).
            ``"log"``
                Logger to handle messages (`lsst.log.Log`).

        Returns
        -------
        output : `bool`
            If true, a correction was applied successfully.

        Raises
        ------
        RuntimeError:
            Raised if the requested row index is out of the table
            bounds.
        """
        numOutOfRange = 0

        rowInd, colIndOffset = kwargs['coeffs'][0:2]
        table = kwargs['table']
        log = kwargs['log']

        numTableRows = table.shape[0]
        rowInd = int(rowInd)
        if rowInd < 0 or rowInd > numTableRows:
            raise RuntimeError("LinearizeLookupTable rowInd=%s not in range[0, %s)" %
                               (rowInd, numTableRows))
        tableRow = table[rowInd, :]
        numOutOfRange += applyLookupTable(image, tableRow, colIndOffset)

        if numOutOfRange > 0 and log is not None:
            log.warn("%s pixels were out of range of the linearization table",
                     numOutOfRange)
        if numOutOfRange < image.getArray().size:
            return True, numOutOfRange
        else:
            return False, numOutOfRange


class LinearizePolynomial(LinearizeBase):
    """Correct non-linearity with a polynomial mode.

    corrImage = uncorrImage + sum_i c_i uncorrImage^(2 + i)

    where c_i are the linearity coefficients for each amplifier.
    Lower order coefficients are not included as they duplicate other
    calibration parameters:
        ``"k0"``
            A coefficient multiplied by uncorrImage**0 is equivalent to
            bias level.  Irrelevant for correcting non-linearity.
        ``"k1"``
            A coefficient multiplied by uncorrImage**1 is proportional
            to the gain.  Not necessary for correcting non-linearity.
    """
    LinearityType = "Polynomial"

    def __call__(self, image, **kwargs):
        """Correct non-linearity.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Image to be corrected
        kwargs : `dict`
            Dictionary of parameter keywords:
            ``"coeffs"``
                Coefficient vector (`list` or `numpy.array`).
                If the order of the polynomial is n, this list
                should have a length of n-1 ("k0" and "k1" are
                not needed for the correction).
            ``"log"``
                Logger to handle messages (`lsst.log.Log`).

        Returns
        -------
        output : `bool`
            If true, a correction was applied successfully.
        """
        if not np.any(np.isfinite(kwargs['coeffs'])):
            return False, 0
        if not np.any(kwargs['coeffs']):
            return False, 0

        ampArray = image.getArray()
        correction = np.zeros_like(ampArray)
        for order, coeff in enumerate(kwargs['coeffs'], start=2):
            correction += coeff * np.power(ampArray, order)
        ampArray += correction

        return True, 0


class LinearizeSquared(LinearizeBase):
    """Correct non-linearity with a squared model.

    corrImage = uncorrImage + c0*uncorrImage^2

    where c0 is linearity coefficient 0 for each amplifier.
    """
    LinearityType = "Squared"

    def __call__(self, image, **kwargs):
        """Correct for non-linearity.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Image to be corrected
        kwargs : `dict`
            Dictionary of parameter keywords:
            ``"coeffs"``
                Coefficient vector (`list` or `numpy.array`).
            ``"log"``
                Logger to handle messages (`lsst.log.Log`).

        Returns
        -------
        output : `bool`
            If true, a correction was applied successfully.
        """

        sqCoeff = kwargs['coeffs'][0]
        if sqCoeff != 0:
            ampArr = image.getArray()
            ampArr *= (1 + sqCoeff*ampArr)
            return True, 0
        else:
            return False, 0


class LinearizeProportional(LinearizeBase):
    """Do not correct non-linearity.
    """
    LinearityType = "Proportional"

    def __call__(self, image, **kwargs):
        """Do not correct for non-linearity.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Image to be corrected
        kwargs : `dict`
            Dictionary of parameter keywords:
            ``"coeffs"``
                Coefficient vector (`list` or `numpy.array`).
            ``"log"``
                Logger to handle messages (`lsst.log.Log`).

        Returns
        -------
        output : `bool`
            If true, a correction was applied successfully.
        """
        return True, 0


class LinearizeNone(LinearizeBase):
    """Do not correct non-linearity.
    """
    LinearityType = "None"

    def __call__(self, image, **kwargs):
        """Do not correct for non-linearity.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Image to be corrected
        kwargs : `dict`
            Dictionary of parameter keywords:
            ``"coeffs"``
                Coefficient vector (`list` or `numpy.array`).
            ``"log"``
                Logger to handle messages (`lsst.log.Log`).

        Returns
        -------
        output : `bool`
            If true, a correction was applied successfully.
        """
        return True, 0
