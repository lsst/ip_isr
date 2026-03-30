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
"""Support for image defects"""

__all__ = ("Defects", "DefectsReason",)

import logging
import itertools
import contextlib
import numpy as np
import math
import numbers
import astropy.table

import lsst.geom
import lsst.afw.table
import lsst.afw.detection
import lsst.afw.image
import lsst.afw.geom
from lsst.meas.algorithms import Defect
from .calibType import IsrCalib

log = logging.getLogger(__name__)

SCHEMA_NAME_KEY = "DEFECTS_SCHEMA"
SCHEMA_VERSION_KEY = "DEFECTS_SCHEMA_VERSION"


class Defects(IsrCalib):
    """Calibration handler for collections of `lsst.meas.algorithms.Defect`.

    Parameters
    ----------
    defectList : iterable, optional
        Collections of defects to apply to the image. Can be an iterable of
        `lsst.meas.algorithms.Defect` or `lsst.geom.BoxI`.
    metadata : `lsst.daf.base.PropertyList`, optional
        Metadata to associate with the defects.  Will be copied and
        overwrite existing metadata, if any. If not supplied the existing
        metadata will be reset.
    normalize_on_init : `bool`
        If True, normalization is applied to the defects in ``defectList`` to
        remove duplicates, eliminate overlaps, etc.

    Notes
    -----
    Defects are stored within this collection in a "reduced" or "normalized"
    form: rather than simply storing the bounding boxes which are added to the
    collection, we eliminate overlaps and duplicates. This normalization
    procedure may introduce overhead when adding many new defects; it may be
    temporarily disabled using the `Defects.bulk_update` context manager if
    necessary.

    The attributes stored in this calibration are:

    _defects : `list` [`lsst.meas.algorithms.Defect`]
        The collection of Defect objects.
    """

    """The calibration type used for ingest."""
    _OBSTYPE = "defects"
    _SCHEMA = ''
    _VERSION = 2.0

    def __init__(self, defectList=None, metadata=None, *, normalize_on_init=True, **kwargs):
        self._defects = []

        if defectList is not None:
            self._bulk_update = True
            for d in defectList:
                self.append(d)
        self._bulk_update = False

        if normalize_on_init:
            self._normalize()

        super().__init__(**kwargs)
        self.requiredAttributes.update(['_defects'])

    def _check_value(self, value):
        """Check that the supplied value is a `~lsst.meas.algorithms.Defect`
        or can be converted to one.

        Parameters
        ----------
        value : `object`
            Value to check.

        Returns
        -------
        new : `~lsst.meas.algorithms.Defect`
            Either the supplied value or a new object derived from it.

        Raises
        ------
        ValueError
            Raised if the supplied value can not be converted to
            `~lsst.meas.algorithms.Defect`
        """
        if isinstance(value, Defect):
            pass
        elif isinstance(value, lsst.geom.BoxI):
            value = Defect(value)
        elif isinstance(value, lsst.geom.PointI):
            value = Defect(lsst.geom.Box2I(value, lsst.geom.Extent2I(1, 1)))
        elif isinstance(value, lsst.afw.image.DefectBase):
            value = Defect(value.getBBox())
        else:
            raise ValueError(f"Defects must be of type Defect, BoxI, or PointI, not '{value!r}'")
        return value

    def __len__(self):
        return len(self._defects)

    def __getitem__(self, index):
        return self._defects[index]

    def __setitem__(self, index, value):
        """Can be given a `~lsst.meas.algorithms.Defect` or a `lsst.geom.BoxI`
        """
        self._defects[index] = self._check_value(value)
        self._normalize()

    def __iter__(self):
        return iter(self._defects)

    def __delitem__(self, index):
        del self._defects[index]

    def __eq__(self, other):
        """Compare if two `Defects` are equal.

        Two `Defects` are equal if their bounding boxes are equal and in
        the same order.  Metadata content is ignored.
        """
        super().__eq__(other)

        if not isinstance(other, self.__class__):
            return False

        # checking the bboxes with zip() only works if same length
        if len(self) != len(other):
            return False

        # Assume equal if bounding boxes are equal
        for d1, d2 in zip(self, other):
            if d1.getBBox() != d2.getBBox():
                return False

        return True

    def __str__(self):
        baseStr = super().__str__()
        return baseStr + ",".join(str(d.getBBox()) for d in self) + ")"

    def _normalize(self):
        """Recalculate defect bounding boxes for efficiency.

        Notes
        -----
        Ideally, this would generate the provably-minimal set of bounding
        boxes necessary to represent the defects. At present, however, that
        doesn't happen: see DM-24781. In the cases of substantial overlaps or
        duplication, though, this will produce a much reduced set.
        """
        # In bulk-update mode, normalization is a no-op.
        if self._bulk_update:
            return

        # If we have no defects, there is nothing to normalize.
        if len(self) == 0:
            return

        # work out the minimum and maximum bounds from all defect regions.
        minX, minY, maxX, maxY = float('inf'), float('inf'), float('-inf'), float('-inf')
        for defect in self:
            bbox = defect.getBBox()
            minX = min(minX, bbox.getMinX())
            minY = min(minY, bbox.getMinY())
            maxX = max(maxX, bbox.getMaxX())
            maxY = max(maxY, bbox.getMaxY())

        region = lsst.geom.Box2I(lsst.geom.Point2I(minX, minY),
                                 lsst.geom.Point2I(maxX, maxY))

        mask = lsst.afw.image.Mask(region)
        self.maskPixels(mask, maskName="BAD")
        self._defects = Defects.fromMask(mask, "BAD")._defects

    @contextlib.contextmanager
    def bulk_update(self):
        """Temporarily suspend normalization of the defect list.
        """
        self._bulk_update = True
        try:
            yield
        finally:
            self._bulk_update = False
            self._normalize()

    def append(self, value):
        self._defects.append(self._check_value(value))
        self._normalize()

    def insert(self, index, value):
        self._defects.insert(index, self._check_value(value))
        self._normalize()

    def copy(self):
        """Copy the defects to a new list, creating new defects from the
        bounding boxes.

        Returns
        -------
        new : `Defects`
            New list with new `Defect` entries.

        Notes
        -----
        This is not a shallow copy in that new `Defect` instances are
        created from the original bounding boxes.  It's also not a deep
        copy since the bounding boxes are not recreated.
        """
        return self.__class__(d.getBBox() for d in self)

    def transpose(self):
        """Make a transposed copy of this defect list.

        Returns
        -------
        retDefectList : `Defects`
            Transposed list of defects.
        """
        retDefectList = self.__class__()
        for defect in self:
            bbox = defect.getBBox()
            dimensions = bbox.getDimensions()
            nbbox = lsst.geom.Box2I(lsst.geom.Point2I(bbox.getMinY(), bbox.getMinX()),
                                    lsst.geom.Extent2I(dimensions[1], dimensions[0]))
            retDefectList.append(nbbox)
        return retDefectList

    def maskPixels(self, mask, maskName="BAD"):
        """Set mask plane based on these defects.

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage` or `lsst.afw.image.Mask`
            Image to process.  Only the mask plane is updated.
        maskName : str, optional
            Mask plane name to use.
        """
        # mask bad pixels
        if hasattr(mask, "getMask"):
            mask = mask.getMask()
        bitmask = mask.getPlaneBitMask(maskName)
        for defect in self:
            bbox = defect.getBBox()
            lsst.afw.geom.SpanSet(bbox).clippedTo(mask.getBBox()).setMask(mask, bitmask)

    def updateCounters(self, columns=None, hot=None, cold=None):
        """Update metadata with pixel and column counts.

        Parameters
        ----------
        columns : `int`, optional
            Number of full columns masked.
        hot : `dict` [`str`, `int`], optional
            Dictionary with the count of hot pixels, indexed by amplifier name.
        cold : `dict` [`str`, `int`], optional
            Dictionary with the count of hot pixels, indexed by amplifier name.
        """
        mdSupplemental = dict()
        if columns:
            mdSupplemental["LSST CALIB DEFECTS N_BAD_COLUMNS"] = columns
        if hot:
            for amp, count in hot.items():
                mdSupplemental[f"LSST CALIB DEFECTS {amp} N_HOT"] = count
        if cold:
            for amp, count in cold.items():
                mdSupplemental[f"LSST CALIB DEFECTS {amp} N_COLD"] = count
        self.getMetadata().update(mdSupplemental)

    def toFitsRegionTable(self):
        """Convert defect list to `~lsst.afw.table.BaseCatalog` using the
        FITS region standard.

        Returns
        -------
        table : `lsst.afw.table.BaseCatalog`
            Defects in tabular form.

        Notes
        -----
        The table created uses the
        `FITS regions <https://fits.gsfc.nasa.gov/registry/region.html>`_
        definition tabular format.  The ``X`` and ``Y`` coordinates are
        converted to FITS Physical coordinates that have origin pixel (1, 1)
        rather than the (0, 0) used in LSST software.
        """
        self.updateMetadata()
        nrows = len(self._defects)

        if nrows:
            # Adding entire columns is more efficient than adding
            # each element separately
            xCol = []
            yCol = []
            rCol = []
            shapes = []
            for i, defect in enumerate(self._defects):
                box = defect.getBBox()
                center = box.getCenter()
                # Correct for the FITS 1-based offset
                xCol.append(center.getX() + 1.0)
                yCol.append(center.getY() + 1.0)

                width = box.width
                height = box.height

                if width == 1 and height == 1:
                    # Call this a point
                    shapeType = "POINT"
                else:
                    shapeType = "BOX"

                # Strings have to be added per row
                shapes.append(shapeType)

                rCol.append(np.array([width, height], dtype=np.float64))

        table = astropy.table.Table({'X': xCol, 'Y': yCol, 'SHAPE': shapes,
                                     'R': rCol, 'ROTANG': np.zeros(nrows),
                                     'COMPONENT': np.arange(nrows)})
        table.meta = self.getMetadata().toDict()
        return table

    @classmethod
    def fromDict(cls, dictionary):
        """Construct a calibration from a dictionary of properties.

        Must be implemented by the specific calibration subclasses.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties.

        Returns
        -------
        calib : `lsst.ip.isr.CalibType`
            Constructed calibration.

        Raises
        ------
        RuntimeError
            Raised if the supplied dictionary is for a different
            calibration.
        """
        calib = cls()

        if calib._OBSTYPE != dictionary['metadata']['OBSTYPE']:
            raise RuntimeError(f"Incorrect crosstalk supplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['metadata']['OBSTYPE']}")

        calib.setMetadata(dictionary['metadata'])
        calib.calibInfoFromDict(dictionary)

        xCol = dictionary['x0']
        yCol = dictionary['y0']
        widthCol = dictionary['width']
        heightCol = dictionary['height']

        with calib.bulk_update:
            for x0, y0, width, height in zip(xCol, yCol, widthCol, heightCol):
                calib.append(lsst.geom.Box2I(lsst.geom.Point2I(x0, y0),
                                             lsst.geom.Extent2I(width, height)))
        return calib

    def toDict(self):
        """Return a dictionary containing the calibration properties.

        The dictionary should be able to be round-tripped through
        `fromDict`.

        Returns
        -------
        dictionary : `dict`
            Dictionary of properties.
        """
        self.updateMetadata()

        outDict = {}
        metadata = self.getMetadata()
        outDict['metadata'] = metadata

        xCol = []
        yCol = []
        widthCol = []
        heightCol = []

        nrows = len(self._defects)
        if nrows:
            for defect in self._defects:
                box = defect.getBBox()
                xCol.append(box.getBeginX())
                yCol.append(box.getBeginY())
                widthCol.append(box.getWidth())
                heightCol.append(box.getHeight())

        outDict['x0'] = xCol
        outDict['y0'] = yCol
        outDict['width'] = widthCol
        outDict['height'] = heightCol

        return outDict

    def toTable(self):
        """Convert defects to a simple table form that we use to write
        to text files.

        Returns
        -------
        table : `lsst.afw.table.BaseCatalog`
            Defects in simple tabular form.

        Notes
        -----
        These defect tables are used as the human readable definitions
        of defects in calibration data definition repositories.  The format
        is to use four columns defined as follows:

        x0 : `int`
            X coordinate of bottom left corner of box.
        y0 : `int`
            Y coordinate of bottom left corner of box.
        width : `int`
            X extent of the box.
        height : `int`
            Y extent of the box.
        """
        tableList = []
        self.updateMetadata()

        xCol = []
        yCol = []
        widthCol = []
        heightCol = []

        nrows = len(self._defects)
        if nrows:
            for defect in self._defects:
                box = defect.getBBox()
                xCol.append(box.getBeginX())
                yCol.append(box.getBeginY())
                widthCol.append(box.getWidth())
                heightCol.append(box.getHeight())

        catalog = astropy.table.Table({'x0': xCol, 'y0': yCol, 'width': widthCol, 'height': heightCol})
        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        catalog.meta = outMeta
        tableList.append(catalog)

        return tableList

    @staticmethod
    def _get_values(values, n=1):
        """Retrieve N values from the supplied values.

        Parameters
        ----------
        values : `numbers.Number` or `list` or `np.array`
            Input values.
        n : `int`
            Number of values to retrieve.

        Returns
        -------
        vals : `list` or `np.array` or `numbers.Number`
            Single value from supplied list if ``n`` is 1, or `list`
            containing first ``n`` values from supplied values.

        Notes
        -----
        Some supplied tables have vectors in some columns that can also
        be scalars.  This method can be used to get the first number as
        a scalar or the first N items from a vector as a vector.
        """
        if n == 1:
            if isinstance(values, numbers.Number):
                return values
            else:
                return values[0]

        return values[:n]

    @classmethod
    def fromTable(cls, tableList, normalize_on_init=True):
        """Construct a `Defects` from the contents of a
        `~lsst.afw.table.BaseCatalog`.

        Parameters
        ----------
        table : `lsst.afw.table.BaseCatalog`
            Table with one row per defect.
        normalize_on_init : `bool`, optional
            If `True`, normalization is applied to the defects listed in the
            table to remove duplicates, eliminate overlaps, etc. Otherwise
            the defects in the returned object exactly match those in the
            table.

        Returns
        -------
        defects : `Defects`
            A `Defects` list.

        Notes
        -----
        Two table formats are recognized.  The first is the
        `FITS regions <https://fits.gsfc.nasa.gov/registry/region.html>`_
        definition tabular format written by `toFitsRegionTable` where the
        pixel origin is corrected from FITS 1-based to a 0-based origin.
        The second is the legacy defects format using columns ``x0``, ``y0``
        (bottom left hand pixel of box in 0-based coordinates), ``width``
        and ``height``.

        The FITS standard regions can only read BOX, POINT, or ROTBOX with
        a zero degree rotation.
        """
        table = tableList[0]
        defectList = []

        schema = table.columns
        # Check schema to see which definitions we have
        if "X" in schema and "Y" in schema and "R" in schema and "SHAPE" in schema:
            # This is a FITS region style table
            isFitsRegion = True
        elif "x0" in schema and "y0" in schema and "width" in schema and "height" in schema:
            # This is a classic LSST-style defect table
            isFitsRegion = False
        else:
            raise ValueError("Unsupported schema for defects extraction")

        for record in table:
            if isFitsRegion:
                # Coordinates can be arrays (some shapes in the standard
                # require this)
                # Correct for FITS 1-based origin
                xcen = cls._get_values(record['X']) - 1.0
                ycen = cls._get_values(record['Y']) - 1.0
                shape = record['SHAPE'].upper().rstrip()
                if shape == "BOX":
                    box = lsst.geom.Box2I.makeCenteredBox(lsst.geom.Point2D(xcen, ycen),
                                                          lsst.geom.Extent2I(cls._get_values(record['R'],
                                                                                             n=2)))
                elif shape == "POINT":
                    # Handle the case where we have an externally created
                    # FITS file.
                    box = lsst.geom.Point2I(xcen, ycen)
                elif shape == "ROTBOX":
                    # Astropy regions always writes ROTBOX
                    rotang = cls._get_values(record['ROTANG'])
                    # We can support 0 or 90 deg
                    if math.isclose(rotang % 90.0, 0.0):
                        # Two values required
                        r = cls._get_values(record['R'], n=2)
                        if math.isclose(rotang % 180.0, 0.0):
                            width = r[0]
                            height = r[1]
                        else:
                            width = r[1]
                            height = r[0]
                        box = lsst.geom.Box2I.makeCenteredBox(lsst.geom.Point2D(xcen, ycen),
                                                              lsst.geom.Extent2I(width, height))
                    else:
                        log.warning("Defect can not be defined using ROTBOX with non-aligned rotation angle")
                        continue
                else:
                    log.warning("Defect lists can only be defined using BOX or POINT not %s", shape)
                    continue

            else:
                # This is a classic LSST-style defect table
                box = lsst.geom.Box2I(lsst.geom.Point2I(record['x0'], record['y0']),
                                      lsst.geom.Extent2I(record['width'], record['height']))

            defectList.append(box)

        defects = cls(defectList, normalize_on_init=normalize_on_init)
        newMeta = dict(table.meta)
        defects.updateMetadata(setCalibInfo=True, **newMeta)

        return defects

    @classmethod
    def readLsstDefectsFile(cls, filename, normalize_on_init=False):
        """Read defects information from a legacy LSST format text file.

        Parameters
        ----------
        filename : `str`
            Name of text file containing the defect information.

        normalize_on_init : `bool`, optional
            If `True`, normalization is applied to the defects listed in the
            table to remove duplicates, eliminate overlaps, etc. Otherwise
            the defects in the returned object exactly match those in the
            table.

        Returns
        -------
        defects : `Defects`
            The defects.

        Notes
        -----
        These defect text files are used as the human readable definitions
        of defects in calibration data definition repositories.  The format
        is to use four columns defined as follows:

        x0 : `int`
            X coordinate of bottom left corner of box.
        y0 : `int`
            Y coordinate of bottom left corner of box.
        width : `int`
            X extent of the box.
        height : `int`
            Y extent of the box.

        Files of this format were used historically to represent defects
        in simple text form.  Use `Defects.readText` and `Defects.writeText`
        to use the more modern format.
        """
        # Use loadtxt so that ValueError is thrown if the file contains a
        # non-integer value. genfromtxt converts bad values to -1.
        defect_array = np.loadtxt(filename,
                                  dtype=[("x0", "int"), ("y0", "int"),
                                         ("x_extent", "int"), ("y_extent", "int")])

        defects = (lsst.geom.Box2I(lsst.geom.Point2I(row["x0"], row["y0"]),
                                   lsst.geom.Extent2I(row["x_extent"], row["y_extent"]))
                   for row in defect_array)

        return cls(defects, normalize_on_init=normalize_on_init)

    @classmethod
    def fromFootprintList(cls, fpList):
        """Compute a defect list from a footprint list, optionally growing
        the footprints.

        Parameters
        ----------
        fpList : `list` of `lsst.afw.detection.Footprint`
            Footprint list to process.

        Returns
        -------
        defects : `Defects`
            List of defects.
        """
        # normalize_on_init is set to False to avoid recursively calling
        # fromMask/fromFootprintList in Defects.__init__.
        return cls(itertools.chain.from_iterable(lsst.afw.detection.footprintToBBoxList(fp)
                                                 for fp in fpList), normalize_on_init=False)

    @classmethod
    def fromMask(cls, mask, maskName):
        """Compute a defect list from a specified mask plane.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask` or `lsst.afw.image.MaskedImage`
            Image to process.
        maskName : `str` or `list`
            Mask plane name, or list of names to convert.

        Returns
        -------
        defects : `Defects`
            Defect list constructed from masked pixels.
        """
        if hasattr(mask, "getMask"):
            mask = mask.getMask()
        thresh = lsst.afw.detection.Threshold(mask.getPlaneBitMask(maskName),
                                              lsst.afw.detection.Threshold.BITMASK)
        fpList = lsst.afw.detection.FootprintSet(mask, thresh).getFootprints()
        return cls.fromFootprintList(fpList)


class DefectsReason(Defects):
    """Handles the definition of the reason(s) for pixels in defects.

    Parameters
    ----------
        exp : `lsst.afw.image.exposure.Exposure`
            Exposure where defects are being found.
    """

    def __init__(self, exp):
        # Image will be populated with defects reason bits.
        self.defectsReasonImage = lsst.afw.image.ImageI(exp.getBBox())

        # List of possible reasons for defects.
        self.reasonList = ['UNLABELED', 'HOT_PIXEL', 'COLD_PIXEL',
                           'HOT_COLUMN', 'COLD_COLUMN', 'MANUAL_DEFECT',
                           'EDGE', 'VAMPIRE_PIXEL']

    def setReasonValue(self, reasonDefect, fpSet):
        """Set value of the defect reason image to the defect reason bit
        at the location of the pixels in a defect.

        Parameters
        ----------
        reasonDefect : `str`
            Reason for the defect.
        fpSet : `lsst.afw.detection._detection.FootprintSet`
            Footprint set of the defects.
        """

        # Check the reason exists
        if reasonDefect in self.reasonList:
            # Get the bit corresponding to the reason
            bitmask = 2**self.bitFromReason(reasonDefect)
            # Set the defect reason image to the corresponding bit value
            lsst.afw.geom.SpanSet(fpSet).setImage(self.defectsReasonImage, bitmask)
        else:
            raise RuntimeError(f"Reason {reasonDefect} not existent.")

    def getReasonMappingDict(self, reasonList):
        """Get the mappings between a reason and a bit and vice versa.

        Parameters
        ----------
        reasonList : `List`
            A list of reason.

        Returns
        -------
        mapping : `str` `int`
            Mapping between a reason and a bit.
        mapping_reverse : `int` `str`
            Mapping between a bit and a reason.
        """

        mapping = {reason: bit for bit, reason in enumerate(reasonList)}
        mapping_reverse = {bit: reason for bit, reason in enumerate(reasonList)}

        return mapping, mapping_reverse

    def reasonFromBit(self, bit):
        """Get the reason from a bit.

        Parameters
        ----------
        bit : `int`
            Reason bit.
        """
        return self.getReasonMappingDict(self.reasonList)[1][bit]

    def bitFromReason(self, reason):
        """Get the bit from a reason.

        Parameters
        ----------
        reason : `str`
            Reason for defect.
        """
        return self.getReasonMappingDict(self.reasonList)[0][reason]

    def addReasonDictMetadata(self):
        """Add the mappings between reason and bit and vice versa in the
        metadata of the defect reason image.
        """

        # The reason image metadata have all the mapping for all the possible
        # reasons wether they are in the image or not.
        exposureMetadata = self.defectsReasonImage.metadata

        # Add the reason to bit correspondance to
        # the defect reason image metadata
        exposureMetadata.update({f"DEFECT REASON {reason}": bit for reason, bit
                                 in self.getReasonMappingDict(self.reasonList)[0].items()})
        exposureMetadata.update({f"DEFECT REASON BIT {bit}": reason for reason, bit
                                 in self.getReasonMappingDict(self.reasonList)[1].items()})

    def interpret(self, value):
        """Outputs the reasons corresponding to a value in a pixel of
        the defect reason image.

        Parameters
        ----------
        value : `int`
            Value in a pixel of the defect reason image.

        Returns
        -------
        reasonsString : `str`
            The reasons corresponding to the value.
        """
        reasonsString = str()
        # loop over possible reasons
        for reason in self.reasonList:
            reasonBit = np.bitwise_and(value, 2**self.bitFromReason(reason))
            # if the resulting bit is different from 0 then append the reason
            if reasonBit != 0:
                if not reasonsString:
                    # if the reason string is empty, we start it
                    reasonsString = reason
                else:
                    # if it is not empty, we append the defects reason
                    # after a comma
                    reasonsString += "," + reason

        return reasonsString

    def getAsString(self, x, y):
        """Outputs the reasons corresponding to a pixel location in
        the defect reason image.

        Parameters
        ----------
        x : `int`
            x coordinate of a pixel in the defect reason image.
        y : `int`
            y coordinate of a pixel in the defect reason image.

        Returns
        -------
        reasonsString : `str`
            The reasons corresponding to the pixel location.
        """
        # Get the value at the x,y location of the exposure
        xyValue = self.defectsReasonImage.image.array[y, x]

        # Return the interpretation of that value
        return self.interpret(xyValue)

    def writeFits(self, path):
        """Write  to a FITS file.

        Parameters
        ----------
        filename : `str`
            Filename to write data to.
        """
        hdr = astropy.io.fits.Header()
        hdr.update(self.defectsReasonImage.metadata)
        hdu = astropy.io.fits.PrimaryHDU(
            data=self.defectsReasonImage, header=hdr
        )
        hdl = astropy.io.fits.HDUList([hdu])
        hdl.writeto(path)

    @classmethod
    def readFits(cls, path):
        """Read a reason image from disk."

        Parameters
        ----------
        path : `str`
            The file to read

        Returns
        -------
        reasonImage : `~lsst.ip.isr.DefectsReason`
        """

        fitsfile = astropy.io.fits.open(path)
        md = fitsfile[0].header

        reasonImage = fitsfile[0].data
        # Get the list of reasons from the fits header
        reasonList = []
        for i in range(len(md)):
            if isinstance(md[i], str):
                reasonList.append(md[i])

        return cls(reasonImage, reasonList)
