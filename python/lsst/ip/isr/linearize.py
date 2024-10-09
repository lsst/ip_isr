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

__all__ = ["Linearizer",
           "LinearizeBase", "LinearizeLookupTable", "LinearizeSquared",
           "LinearizeProportional", "LinearizePolynomial", "LinearizeSpline", "LinearizeNone"]

import abc
import numpy as np

from astropy.table import Table

import lsst.afw.math as afwMath
from lsst.pipe.base import Struct
from lsst.geom import Box2I, Point2I, Extent2I
from .applyLookupTable import applyLookupTable
from .calibType import IsrCalib


class Linearizer(IsrCalib):
    """Parameter set for linearization.

    These parameters are included in `lsst.afw.cameraGeom.Amplifier`, but
    should be accessible externally to allow for testing.

    Parameters
    ----------
    table : `numpy.array`, optional
        Lookup table; a 2-dimensional array of floats:

        - one row for each row index (value of coef[0] in the amplifier)
        - one column for each image value

        To avoid copying the table the last index should vary fastest
        (numpy default "C" order)
    detector : `lsst.afw.cameraGeom.Detector`, optional
        Detector object.  Passed to self.fromDetector() on init.
    log : `logging.Logger`, optional
        Logger to handle messages.
    kwargs : `dict`, optional
        Other keyword arguments to pass to the parent init.

    Raises
    ------
    RuntimeError
        Raised if the supplied table is not 2D, or if the table has fewer
        columns than rows (indicating that the indices are swapped).

    Notes
    -----
    The linearizer attributes stored are:

    hasLinearity : `bool`
        Whether a linearity correction is defined for this detector.
    override : `bool`
        Whether the detector parameters should be overridden.
    ampNames : `list` [`str`]
        List of amplifier names to correct.
    linearityCoeffs : `dict` [`str`, `numpy.array`]
        Coefficients to use in correction.  Indexed by amplifier
        names.  The format of the array depends on the type of
        correction to apply.
    linearityType : `dict` [`str`, `str`]
        Type of correction to use, indexed by amplifier names.
    linearityBBox : `dict` [`str`, `lsst.geom.Box2I`]
        Bounding box the correction is valid over, indexed by
        amplifier names.
    fitParams : `dict` [`str`, `numpy.array`], optional
        Linearity fit parameters used to construct the correction
        coefficients, indexed as above.
    fitParamsErr : `dict` [`str`, `numpy.array`], optional
        Uncertainty values of the linearity fit parameters used to
        construct the correction coefficients, indexed as above.
    fitChiSq : `dict` [`str`, `float`], optional
        Chi-squared value of the linearity fit, indexed as above.
    fitResiduals : `dict` [`str`, `numpy.array`], optional
        Residuals of the fit, indexed as above. Used for
        calculating photdiode corrections
    fitResidualsSigmaMad : `dict` [`str`, `float`], optional
        Robust median-absolute-deviation of fit residuals, scaled
        by the signal level.
    linearFit : The linear fit to the low flux region of the curve.
        [intercept, slope].
    tableData : `numpy.array`, optional
        Lookup table data for the linearity correction.

    Notes
    -----
    Version 1.4 adds ``linearityTurnoff`` and ``linearityMaxSignal``.
    """
    _OBSTYPE = "LINEARIZER"
    _SCHEMA = 'Gen3 Linearizer'
    _VERSION = 1.4

    def __init__(self, table=None, **kwargs):
        self.hasLinearity = False
        self.override = False

        self.ampNames = list()
        self.linearityCoeffs = dict()
        self.linearityType = dict()
        self.linearityBBox = dict()
        self.fitParams = dict()
        self.fitParamsErr = dict()
        self.fitChiSq = dict()
        self.fitResiduals = dict()
        self.fitResidualsSigmaMad = dict()
        self.linearFit = dict()
        self.linearityTurnoff = dict()
        self.linearityMaxSignal = dict()
        self.tableData = None
        if table is not None:
            if len(table.shape) != 2:
                raise RuntimeError("table shape = %s; must have two dimensions" % (table.shape,))
            if table.shape[1] < table.shape[0]:
                raise RuntimeError("table shape = %s; indices are switched" % (table.shape,))
            self.tableData = np.array(table, order="C")

        self.linearityUnits = 'adu'

        super().__init__(**kwargs)
        self.requiredAttributes.update(['hasLinearity', 'override',
                                        'ampNames',
                                        'linearityCoeffs', 'linearityType', 'linearityBBox',
                                        'fitParams', 'fitParamsErr', 'fitChiSq',
                                        'fitResiduals', 'fitResidualsSigmaMad', 'linearFit', 'tableData',
                                        'units', 'linearityTurnoff', 'linearityMaxSignal'])

    def updateMetadata(self, setDate=False, **kwargs):
        """Update metadata keywords with new values.

        This calls the base class's method after ensuring the required
        calibration keywords will be saved.

        Parameters
        ----------
        setDate : `bool`, optional
            Update the CALIBDATE fields in the metadata to the current
            time. Defaults to False.
        kwargs :
            Other keyword parameters to set in the metadata.
        """
        kwargs['HAS_LINEARITY'] = self.hasLinearity
        kwargs['OVERRIDE'] = self.override
        kwargs['HAS_TABLE'] = self.tableData is not None
        kwargs['LINEARITY_UNITS'] = self.linearityUnits

        super().updateMetadata(setDate=setDate, **kwargs)

    def fromDetector(self, detector):
        """Read linearity parameters from a detector.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.detector`
            Input detector with parameters to use.

        Returns
        -------
        calib : `lsst.ip.isr.Linearizer`
            The calibration constructed from the detector.
        """
        self._detectorName = detector.getName()
        self._detectorSerial = detector.getSerial()
        self._detectorId = detector.getId()
        self.hasLinearity = True

        # Do not translate Threshold, Maximum, Units.
        for amp in detector.getAmplifiers():
            ampName = amp.getName()
            self.ampNames.append(ampName)
            self.linearityType[ampName] = amp.getLinearityType()
            self.linearityCoeffs[ampName] = amp.getLinearityCoeffs()
            self.linearityBBox[ampName] = amp.getBBox()

        self.linearityUnits = 'adu'

        return self

    @classmethod
    def fromDict(cls, dictionary):
        """Construct a calibration from a dictionary of properties

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties

        Returns
        -------
        calib : `lsst.ip.isr.Linearity`
            Constructed calibration.

        Raises
        ------
        RuntimeError
            Raised if the supplied dictionary is for a different
            calibration.
        """

        calib = cls()

        if calib._OBSTYPE != dictionary['metadata']['OBSTYPE']:
            raise RuntimeError(f"Incorrect linearity supplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['metadata']['OBSTYPE']}")

        calib.setMetadata(dictionary['metadata'])

        calib.hasLinearity = dictionary.get('hasLinearity',
                                            dictionary['metadata'].get('HAS_LINEARITY', False))
        calib.override = dictionary.get('override', True)

        calib.linearityUnits = dictionary.get('linearityUnits', 'adu')

        if calib.hasLinearity:
            for ampName in dictionary['amplifiers']:
                amp = dictionary['amplifiers'][ampName]
                calib.ampNames.append(ampName)
                calib.linearityCoeffs[ampName] = np.array(amp.get('linearityCoeffs', [0.0]))
                calib.linearityType[ampName] = amp.get('linearityType', 'None')
                calib.linearityBBox[ampName] = amp.get('linearityBBox', None)

                calib.fitParams[ampName] = np.array(amp.get('fitParams', [0.0]))
                calib.fitParamsErr[ampName] = np.array(amp.get('fitParamsErr', [0.0]))
                calib.fitChiSq[ampName] = amp.get('fitChiSq', np.nan)
                calib.fitResiduals[ampName] = np.array(amp.get('fitResiduals', [0.0]))
                calib.fitResidualsSigmaMad[ampName] = np.array(amp.get('fitResidualsSigmaMad', np.nan))
                calib.linearFit[ampName] = np.array(amp.get('linearFit', [0.0]))

                calib.linearityTurnoff[ampName] = np.array(amp.get('linearityTurnoff', np.nan))
                calib.linearityMaxSignal[ampName] = np.array(amp.get('linearityMaxSignal', np.nan))

            calib.tableData = dictionary.get('tableData', None)
            if calib.tableData:
                calib.tableData = np.array(calib.tableData)

        return calib

    def toDict(self):
        """Return linearity parameters as a dict.

        Returns
        -------
        outDict : `dict`:
        """
        self.updateMetadata()

        outDict = {'metadata': self.getMetadata(),
                   'detectorName': self._detectorName,
                   'detectorSerial': self._detectorSerial,
                   'detectorId': self._detectorId,
                   'hasTable': self.tableData is not None,
                   'amplifiers': dict(),
                   'linearityUnits': self.linearityUnits,
                   }
        for ampName in self.linearityType:
            outDict['amplifiers'][ampName] = {
                'linearityType': self.linearityType[ampName],
                'linearityCoeffs': self.linearityCoeffs[ampName].tolist(),
                'linearityBBox': self.linearityBBox[ampName],
                'fitParams': self.fitParams[ampName].tolist(),
                'fitParamsErr': self.fitParamsErr[ampName].tolist(),
                'fitChiSq': self.fitChiSq[ampName],
                'fitResiduals': self.fitResiduals[ampName].tolist(),
                'fitResidualsSigmaMad': self.fitResiduals[ampName],
                'linearFit': self.linearFit[ampName].tolist(),
                'linearityTurnoff': self.linearityTurnoff[ampName],
                'linearityMaxSignal': self.linearityMaxSignal[ampName],
            }
        if self.tableData is not None:
            outDict['tableData'] = self.tableData.tolist()

        return outDict

    @classmethod
    def fromTable(cls, tableList):
        """Read linearity from a FITS file.

        This method uses the `fromDict` method to create the
        calibration, after constructing an appropriate dictionary from
        the input tables.

        Parameters
        ----------
        tableList : `list` [`astropy.table.Table`]
            afwTable read from input file name.

        Returns
        -------
        linearity : `~lsst.ip.isr.linearize.Linearizer``
            Linearity parameters.

        Notes
        -----
        The method reads a FITS file with 1 or 2 extensions. The metadata is
        read from the header of extension 1, which must exist.  Then the table
        is loaded, and  the ['AMPLIFIER_NAME', 'TYPE', 'COEFFS', 'BBOX_X0',
        'BBOX_Y0', 'BBOX_DX', 'BBOX_DY'] columns are read and used to set each
        dictionary by looping over rows.
        Extension 2 is then attempted to read in the try block (which only
        exists for lookup tables). It has a column named 'LOOKUP_VALUES' that
        contains a vector of the lookup entries in each row.
        """
        coeffTable = tableList[0]

        metadata = coeffTable.meta
        inDict = dict()
        inDict['metadata'] = metadata
        inDict['hasLinearity'] = metadata.get('HAS_LINEARITY', False)
        inDict['amplifiers'] = dict()
        inDict['linearityUnits'] = metadata.get('LINEARITY_UNITS', 'adu')

        for record in coeffTable:
            ampName = record['AMPLIFIER_NAME']

            fitParams = record['FIT_PARAMS'] if 'FIT_PARAMS' in record.columns else np.array([0.0])
            fitParamsErr = record['FIT_PARAMS_ERR'] if 'FIT_PARAMS_ERR' in record.columns else np.array([0.0])
            fitChiSq = record['RED_CHI_SQ'] if 'RED_CHI_SQ' in record.columns else np.nan
            fitResiduals = record['FIT_RES'] if 'FIT_RES' in record.columns else np.array([0.0])
            fitResidualsSigmaMad = record['FIT_RES_SIGMAD'] if 'FIT_RES_SIGMAD' in record.columns else np.nan
            linearFit = record['LIN_FIT'] if 'LIN_FIT' in record.columns else np.array([0.0])

            linearityTurnoff = record['LINEARITY_TURNOFF'] if 'LINEARITY_TURNOFF' in record.columns \
                else np.nan
            linearityMaxSignal = record['LINEARITY_MAX_SIGNAL'] if 'LINEARITY_MAX_SIGNAL' in record.columns \
                else np.nan

            inDict['amplifiers'][ampName] = {
                'linearityType': record['TYPE'],
                'linearityCoeffs': record['COEFFS'],
                'linearityBBox': Box2I(Point2I(record['BBOX_X0'], record['BBOX_Y0']),
                                       Extent2I(record['BBOX_DX'], record['BBOX_DY'])),
                'fitParams': fitParams,
                'fitParamsErr': fitParamsErr,
                'fitChiSq': fitChiSq,
                'fitResiduals': fitResiduals,
                'fitResidualsSigmaMad': fitResidualsSigmaMad,
                'linearFit': linearFit,
                'linearityTurnoff': linearityTurnoff,
                'linearityMaxSignal': linearityMaxSignal,
            }

        if len(tableList) > 1:
            tableData = tableList[1]
            inDict['tableData'] = [record['LOOKUP_VALUES'] for record in tableData]

        return cls().fromDict(inDict)

    def toTable(self):
        """Construct a list of tables containing the information in this
        calibration.

        The list of tables should create an identical calibration
        after being passed to this class's fromTable method.

        Returns
        -------
        tableList : `list` [`astropy.table.Table`]
            List of tables containing the linearity calibration
            information.
        """

        tableList = []
        self.updateMetadata()
        catalog = Table([{'AMPLIFIER_NAME': ampName,
                          'TYPE': self.linearityType[ampName],
                          'COEFFS': self.linearityCoeffs[ampName],
                          'BBOX_X0': self.linearityBBox[ampName].getMinX(),
                          'BBOX_Y0': self.linearityBBox[ampName].getMinY(),
                          'BBOX_DX': self.linearityBBox[ampName].getWidth(),
                          'BBOX_DY': self.linearityBBox[ampName].getHeight(),
                          'FIT_PARAMS': self.fitParams[ampName],
                          'FIT_PARAMS_ERR': self.fitParamsErr[ampName],
                          'RED_CHI_SQ': self.fitChiSq[ampName],
                          'FIT_RES': self.fitResiduals[ampName],
                          'FIT_RES_SIGMAD': self.fitResidualsSigmaMad[ampName],
                          'LIN_FIT': self.linearFit[ampName],
                          'LINEARITY_TURNOFF': self.linearityTurnoff[ampName],
                          'LINEARITY_MAX_SIGNAL': self.linearityMaxSignal[ampName],
                          } for ampName in self.ampNames])
        catalog.meta = self.getMetadata().toDict()
        tableList.append(catalog)

        if self.tableData is not None:
            catalog = Table([{'LOOKUP_VALUES': value} for value in self.tableData])
            tableList.append(catalog)
        return tableList

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
                  LinearizeSpline,
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
        RuntimeError
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
            # TODO: DM-38778: This check fails on LATISS due to an
            #                 error in the camera configuration.
            # if self._detectorSerial != detector.getSerial():
            #      raise RuntimeError(
            #          "Detector serial numbers don't match: %s != %s" %
            #          (self._detectorSerial, detector.getSerial()))
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
                    self.log.debug("Overriding amplifier defined linearityType (%s) for %s",
                                   self.linearityType[ampName], ampName)
                else:
                    raise RuntimeError("Amplifier %s type %s does not match saved value %s" %
                                       (ampName, amp.getLinearityType(), self.linearityType[ampName]))
            if (amp.getLinearityCoeffs().shape != self.linearityCoeffs[ampName].shape or not
                    np.allclose(amp.getLinearityCoeffs(), self.linearityCoeffs[ampName], equal_nan=True)):
                if self.override:
                    self.log.debug("Overriding amplifier defined linearityCoeffs (%s) for %s",
                                   self.linearityCoeffs[ampName], ampName)
                else:
                    raise RuntimeError("Amplifier %s coeffs %s does not match saved value %s" %
                                       (ampName, amp.getLinearityCoeffs(), self.linearityCoeffs[ampName]))

    def applyLinearity(self, image, detector=None, log=None, gains=None):
        """Apply the linearity to an image.

        If the linearity parameters are populated, use those,
        otherwise use the values from the detector.

        Parameters
        ----------
        image : `~lsst.afw.image.image`
            Image to correct.
        detector : `~lsst.afw.cameraGeom.detector`, optional
            Detector to use to determine exposure trimmed state.  If
            supplied, but no other linearity information is provided
            by the calibration, then the static solution stored in the
            detector will be used.
        log : `~logging.Logger`, optional
            Log object to use for logging.
        gains : `dict` [`str`, `float`], optional
            Dictionary of amp name to gain. If this is provided then
            linearity terms will be converted from adu to electrons.
            Only used for Spline linearity corrections.
        """
        if log is None:
            log = self.log
        if detector and not self.hasLinearity:
            self.fromDetector(detector)

        self.validate(detector)

        isTrimmed = None
        if detector:
            if detector.getBBox() == image.getBBox():
                isTrimmed = True
            else:
                isTrimmed = False

        numAmps = 0
        numLinearized = 0
        numOutOfRange = 0
        for ampName in self.linearityType.keys():
            linearizer = self.getLinearityTypeByName(self.linearityType[ampName])
            numAmps += 1

            if gains and self.linearityUnits == 'adu':
                gainValue = gains[ampName]
            else:
                gainValue = 1.0

            if linearizer is not None:
                match isTrimmed:
                    case True:
                        bbox = detector[ampName].getBBox()
                    case False:
                        bbox = detector[ampName].getRawBBox()
                    case None:
                        bbox = self.linearityBBox[ampName]

                ampView = image.Factory(image, bbox)
                success, outOfRange = linearizer()(ampView, **{'coeffs': self.linearityCoeffs[ampName],
                                                               'table': self.tableData,
                                                               'log': self.log,
                                                               'gain': gainValue})
                numOutOfRange += outOfRange
                if success:
                    numLinearized += 1
                elif log is not None:
                    log.warning("Amplifier %s did not linearize.",
                                ampName)
        return Struct(
            numAmps=numAmps,
            numLinearized=numLinearized,
            numOutOfRange=numOutOfRange
        )


class LinearizeBase(metaclass=abc.ABCMeta):
    """Abstract base class functor for correcting non-linearity.

    Subclasses must define ``__call__`` and set class variable
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

            ``coeffs``
                Coefficient vector (`list` or `numpy.array`).
            ``table``
                Lookup table data (`numpy.array`).
            ``log``
                Logger to handle messages (`logging.Logger`).

        Returns
        -------
        output : `bool`
            If `True`, a correction was applied successfully.

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

            ``coeffs``
                Columnation vector (`list` or `numpy.array`).
            ``table``
                Lookup table data (`numpy.array`).
            ``log``
                Logger to handle messages (`logging.Logger`).

        Returns
        -------
        output : `tuple` [`bool`, `int`]
            If true, a correction was applied successfully.  The
            integer indicates the number of pixels that were
            uncorrectable by being out of range.

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
        tableRow = np.array(table[rowInd, :], dtype=image.getArray().dtype)

        numOutOfRange += applyLookupTable(image, tableRow, colIndOffset)

        if numOutOfRange > 0 and log is not None:
            log.warning("%s pixels were out of range of the linearization table",
                        numOutOfRange)
        if numOutOfRange < image.getArray().size:
            return True, numOutOfRange
        else:
            return False, numOutOfRange


class LinearizePolynomial(LinearizeBase):
    """Correct non-linearity with a polynomial mode.

    .. code-block::

        corrImage = uncorrImage + sum_i c_i uncorrImage^(2 + i)

    where ``c_i`` are the linearity coefficients for each amplifier.
    Lower order coefficients are not included as they duplicate other
    calibration parameters:

    ``k0``
        A coefficient multiplied by ``uncorrImage**0`` is equivalent to
        bias level.  Irrelevant for correcting non-linearity.
    ``k1``
        A coefficient multiplied by ``uncorrImage**1`` is proportional
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

            ``coeffs``
                Coefficient vector (`list` or `numpy.array`).
                If the order of the polynomial is n, this list
                should have a length of n-1 ("k0" and "k1" are
                not needed for the correction).
            ``log``
                Logger to handle messages (`logging.Logger`).

        Returns
        -------
        output : `tuple` [`bool`, `int`]
            If true, a correction was applied successfully.  The
            integer indicates the number of pixels that were
            uncorrectable by being out of range.
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

            ``coeffs``
                Coefficient vector (`list` or `numpy.array`).
            ``log``
                Logger to handle messages (`logging.Logger`).

        Returns
        -------
        output : `tuple` [`bool`, `int`]
            If true, a correction was applied successfully.  The
            integer indicates the number of pixels that were
            uncorrectable by being out of range.
        """

        sqCoeff = kwargs['coeffs'][0]
        if sqCoeff != 0:
            ampArr = image.getArray()
            ampArr *= (1 + sqCoeff*ampArr)
            return True, 0
        else:
            return False, 0


class LinearizeSpline(LinearizeBase):
    """Correct non-linearity with a spline model.

    corrImage = uncorrImage - Spline(coeffs, uncorrImage)

    Notes
    -----

    The spline fit calculates a correction as a function of the
    expected linear flux term.  Because of this, the correction needs
    to be subtracted from the observed flux.

    """
    LinearityType = "Spline"

    def __call__(self, image, **kwargs):
        """Correct for non-linearity.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Image to be corrected
        kwargs : `dict`
            Dictionary of parameter keywords:

            ``coeffs``
                Coefficient vector (`list` or `numpy.array`).
            ``log``
                Logger to handle messages (`logging.Logger`).
            ``gain``
                Gain value to apply.

        Returns
        -------
        output : `tuple` [`bool`, `int`]
            If true, a correction was applied successfully.  The
            integer indicates the number of pixels that were
            uncorrectable by being out of range.
        """
        splineCoeff = kwargs['coeffs']
        gain = kwargs.get('gain', 1.0)
        centers, values = np.split(splineCoeff, 2)
        values = values*gain
        # If the spline is not anchored at zero, remove the offset
        # found at the lowest flux available, and add an anchor at
        # flux=0.0 if there's no entry at that point.
        if values[0] != 0:
            offset = values[0]
            values -= offset
        if centers[0] != 0.0:
            centers = np.concatenate(([0.0], centers))
            values = np.concatenate(([0.0], values))

        interp = afwMath.makeInterpolate(centers.tolist(), values.tolist(),
                                         afwMath.stringToInterpStyle("AKIMA_SPLINE"))

        ampArr = image.getArray()
        delta = interp.interpolate(ampArr.ravel())
        ampArr -= np.array(delta).reshape(ampArr.shape)

        return True, 0


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

            ``coeffs``
                Coefficient vector (`list` or `numpy.array`).
            ``log``
                Logger to handle messages (`logging.Logger`).

        Returns
        -------
        output : `tuple` [`bool`, `int`]
            If true, a correction was applied successfully.  The
            integer indicates the number of pixels that were
            uncorrectable by being out of range.
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

            ``coeffs``
                Coefficient vector (`list` or `numpy.array`).
            ``log``
                Logger to handle messages (`logging.Logger`).

        Returns
        -------
        output : `tuple` [`bool`, `int`]
            If true, a correction was applied successfully.  The
            integer indicates the number of pixels that were
            uncorrectable by being out of range.
        """
        return True, 0
