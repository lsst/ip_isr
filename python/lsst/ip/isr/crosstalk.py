#
# LSST Data Management System
# Copyright 2008-2017 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
"""
Apply intra-detector crosstalk corrections
"""

__all__ = ["CrosstalkCalib", "CrosstalkConfig", "CrosstalkTask",
           "NullCrosstalkTask"]

import numpy as np
from astropy.table import Table

import lsst.afw.math
import lsst.afw.detection
import lsst.daf.butler
from lsst.pex.config import Config, Field, ChoiceField, ListField
from lsst.pipe.base import Task

from lsst.ip.isr import IsrCalib


class CrosstalkCalib(IsrCalib):
    """Calibration of amp-to-amp crosstalk coefficients.

    Parameters
    ----------
    detector : `lsst.afw.cameraGeom.Detector`, optional
        Detector to use to pull coefficients from.
    nAmp : `int`, optional
        Number of amplifiers to initialize.
    log : `logging.Logger`, optional
        Log to write messages to.
    **kwargs :
        Parameters to pass to parent constructor.

    Notes
    -----
    The crosstalk attributes stored are:

    hasCrosstalk : `bool`
        Whether there is crosstalk defined for this detector.
    nAmp : `int`
        Number of amplifiers in this detector.
    crosstalkShape : `tuple` [`int`, `int`]
        A tuple containing the shape of the ``coeffs`` matrix.  This
        should be equivalent to (``nAmp``, ``nAmp``).
    coeffs : `numpy.ndarray`
        A matrix containing the crosstalk coefficients.  coeff[i][j]
        contains the coefficients to calculate the contribution
        amplifier_j has on amplifier_i (each row[i] contains the
        corrections for detector_i).
    coeffErr : `numpy.ndarray`, optional
        A matrix (as defined by ``coeffs``) containing the standard
        distribution of the crosstalk measurements.
    coeffNum : `numpy.ndarray`, optional
        A matrix containing the number of pixel pairs used to measure
        the ``coeffs`` and ``coeffErr``.
    coeffValid : `numpy.ndarray`, optional
        A matrix of Boolean values indicating if the coefficient is
        valid, defined as abs(coeff) > coeffErr / sqrt(coeffNum).
    coeffsSqr : `numpy.ndarray`, optional
        A matrix containing potential quadratic crosstalk coefficients
        (see e.g., Snyder+21, 2001.03223). coeffsSqr[i][j]
        contains the coefficients to calculate the contribution
        amplifier_j has on amplifier_i (each row[i] contains the
        corrections for detector_i).
    coeffErrSqr : `numpy.ndarray`, optional
        A matrix (as defined by ``coeffsSqr``) containing the standard
        distribution of the quadratic term of the crosstalk measurements.
    interChip : `dict` [`numpy.ndarray`]
        A dictionary keyed by detectorName containing ``coeffs``
        matrices used to correct for inter-chip crosstalk with a
        source on the detector indicated.

    Version 1.1 adds quadratic coefficients, a matrix with the ratios
    of amplifiers gains per detector, and a field to indicate the units
    of the numerator and denominator of the source and target signals, with
    "adu" meaning "ADU / ADU" and "electron" meaning "e- / e-".
    """
    _OBSTYPE = 'CROSSTALK'
    _SCHEMA = 'Gen3 Crosstalk'
    _VERSION = 1.1

    def __init__(self, detector=None, nAmp=0, **kwargs):
        self.hasCrosstalk = False
        self.nAmp = nAmp if nAmp else 0
        self.crosstalkShape = (self.nAmp, self.nAmp)

        self.coeffs = np.zeros(self.crosstalkShape) if self.nAmp else None
        self.coeffErr = np.zeros(self.crosstalkShape) if self.nAmp else None
        self.coeffNum = np.zeros(self.crosstalkShape,
                                 dtype=int) if self.nAmp else None
        self.coeffValid = np.zeros(self.crosstalkShape,
                                   dtype=bool) if self.nAmp else None
        # Quadratic terms, if any.
        self.coeffsSqr = np.zeros(self.crosstalkShape) if self.nAmp else None
        self.coeffErrSqr = np.zeros(self.crosstalkShape) if self.nAmp else None

        # Gain ratios
        self.ampGainRatios = np.zeros(self.crosstalkShape) if self.nAmp else None

        # Units
        self.crosstalkRatiosUnits = 'adu' if self.nAmp else None

        self.interChip = {}

        super().__init__(**kwargs)
        self.requiredAttributes.update(['hasCrosstalk', 'nAmp', 'coeffs',
                                        'coeffErr', 'coeffNum', 'coeffValid',
                                        'coeffsSqr', 'coeffErrSqr',
                                        'ampGainRatios', 'crosstalkRatiosUnits',
                                        'interChip'])
        if detector:
            self.fromDetector(detector)

    def updateMetadata(self, setDate=False, **kwargs):
        """Update calibration metadata.

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
        kwargs['DETECTOR'] = self._detectorId
        kwargs['DETECTOR_NAME'] = self._detectorName
        kwargs['DETECTOR_SERIAL'] = self._detectorSerial
        kwargs['HAS_CROSSTALK'] = self.hasCrosstalk
        kwargs['NAMP'] = self.nAmp
        self.crosstalkShape = (self.nAmp, self.nAmp)
        kwargs['CROSSTALK_SHAPE'] = self.crosstalkShape
        kwargs['CROSSTALK_RATIOS_UNITS'] = self.crosstalkRatiosUnits

        super().updateMetadata(setDate=setDate, **kwargs)

    def fromDetector(self, detector, coeffVector=None, coeffSqrVector=None):
        """Set calibration parameters from the detector.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            Detector to use to set parameters from.
        coeffVector : `numpy.array`, optional
            Use the detector geometry (bounding boxes and flip
            information), but use ``coeffVector`` instead of the
            output of ``detector.getCrosstalk()``.
        coeffSqrVector : `numpy.array`, optional
            Quadratic crosstalk coefficients.

        Returns
        -------
        calib : `lsst.ip.isr.CrosstalkCalib`
            The calibration constructed from the detector.

        """
        self._detectorId = detector.getId()
        self._detectorName = detector.getName()
        self._detectorSerial = detector.getSerial()

        self.nAmp = len(detector)
        self.crosstalkShape = (self.nAmp, self.nAmp)

        if coeffVector is not None:
            crosstalkCoeffs = coeffVector
        else:
            crosstalkCoeffs = detector.getCrosstalk()
        if coeffSqrVector is not None:
            self.coeffsSqr = coeffSqrVector
        else:
            self.coeffsSqr = np.zeros(self.crosstalkShape)
        if len(crosstalkCoeffs) == 1 and crosstalkCoeffs[0] == 0.0:
            return self
        self.coeffs = np.array(crosstalkCoeffs).reshape(self.crosstalkShape)

        if self.coeffs.shape != self.crosstalkShape:
            raise RuntimeError("Crosstalk coefficients do not match detector shape. "
                               f"{self.crosstalkShape} {self.nAmp}")
        # Set default as in __init__
        self.coeffErr = np.zeros(self.crosstalkShape)
        self.coeffNum = np.zeros(self.crosstalkShape, dtype=int)
        self.coeffValid = np.ones(self.crosstalkShape, dtype=bool)
        self.coeffErrSqr = np.zeros(self.crosstalkShape)
        self.ampGainRatios = np.zeros(self.crosstalkShape)
        self.crosstalkRatiosUnits = 'adu'

        self.interChip = {}

        self.hasCrosstalk = True
        self.updateMetadata()

        return self

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

        if 'detectorName' in dictionary:
            calib._detectorName = dictionary.get('detectorName')
        elif 'DETECTOR_NAME' in dictionary:
            calib._detectorName = dictionary.get('DETECTOR_NAME')
        elif 'DET_NAME' in dictionary['metadata']:
            calib._detectorName = dictionary['metadata']['DET_NAME']
        else:
            calib._detectorName = None

        if 'detectorSerial' in dictionary:
            calib._detectorSerial = dictionary.get('detectorSerial')
        elif 'DETECTOR_SERIAL' in dictionary:
            calib._detectorSerial = dictionary.get('DETECTOR_SERIAL')
        elif 'DET_SER' in dictionary['metadata']:
            calib._detectorSerial = dictionary['metadata']['DET_SER']
        else:
            calib._detectorSerial = None

        if 'detectorId' in dictionary:
            calib._detectorId = dictionary.get('detectorId')
        elif 'DETECTOR' in dictionary:
            calib._detectorId = dictionary.get('DETECTOR')
        elif 'DETECTOR' in dictionary['metadata']:
            calib._detectorId = dictionary['metadata']['DETECTOR']
        elif calib._detectorSerial:
            calib._detectorId = calib._detectorSerial
        else:
            calib._detectorId = None

        if 'instrument' in dictionary:
            calib._instrument = dictionary.get('instrument')
        elif 'INSTRUME' in dictionary['metadata']:
            calib._instrument = dictionary['metadata']['INSTRUME']
        else:
            calib._instrument = None

        calib.hasCrosstalk = dictionary.get('hasCrosstalk',
                                            dictionary['metadata'].get('HAS_CROSSTALK', False))
        if calib.hasCrosstalk:
            calib.nAmp = dictionary.get('nAmp', dictionary['metadata'].get('NAMP', 0))
            calib.crosstalkShape = (calib.nAmp, calib.nAmp)
            calib.coeffs = np.array(dictionary['coeffs']).reshape(calib.crosstalkShape)
            calib.crosstalkRatiosUnits = dictionary.get(
                'crosstalkRatiosUnits',
                dictionary['metadata'].get('CROSSTALK_RATIOS_UNITS', None))
            if 'coeffErr' in dictionary:
                calib.coeffErr = np.array(dictionary['coeffErr']).reshape(calib.crosstalkShape)
            else:
                calib.coeffErr = np.zeros_like(calib.coeffs)
            if 'coeffNum' in dictionary:
                calib.coeffNum = np.array(dictionary['coeffNum']).reshape(calib.crosstalkShape)
            else:
                calib.coeffNum = np.zeros_like(calib.coeffs, dtype=int)
            if 'coeffValid' in dictionary:
                calib.coeffValid = np.array(dictionary['coeffValid']).reshape(calib.crosstalkShape)
            else:
                calib.coeffValid = np.ones_like(calib.coeffs, dtype=bool)
            if 'coeffsSqr' in dictionary:
                calib.coeffsSqr = np.array(dictionary['coeffsSqr']).reshape(calib.crosstalkShape)
            else:
                calib.coeffsSqr = np.zeros_like(calib.coeffs)
            if 'coeffErrSqr' in dictionary:
                calib.coeffErrSqr = np.array(dictionary['coeffErrSqr']).reshape(calib.crosstalkShape)
            else:
                calib.coeffErrSqr = np.zeros_like(calib.coeffs)
            if 'ampGainRatios' in dictionary:
                calib.ampGainRatios = np.array(dictionary['ampGainRatios']).reshape(calib.crosstalkShape)
            else:
                calib.ampGainRatios = np.zeros_like(calib.coeffs)
            if 'crosstalkRatiosUnits' in dictionary:
                calib.crosstalkRatiosUnits = dictionary['crosstalkRatiosUnits']
            else:
                calib.crosstalkRatiosUnits = None

            calib.interChip = dictionary.get('interChip', None)
            if calib.interChip:
                for detector in calib.interChip:
                    coeffVector = calib.interChip[detector]
                    calib.interChip[detector] = np.array(coeffVector).reshape(calib.crosstalkShape)

        calib.updateMetadata()
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

        outDict['hasCrosstalk'] = self.hasCrosstalk
        outDict['nAmp'] = self.nAmp
        outDict['crosstalkShape'] = self.crosstalkShape
        outDict['crosstalkRatiosUnits'] = self.crosstalkRatiosUnits

        ctLength = self.nAmp*self.nAmp
        outDict['coeffs'] = self.coeffs.reshape(ctLength).tolist()

        if self.coeffErr is not None:
            outDict['coeffErr'] = self.coeffErr.reshape(ctLength).tolist()
        if self.coeffNum is not None:
            outDict['coeffNum'] = self.coeffNum.reshape(ctLength).tolist()
        if self.coeffValid is not None:
            outDict['coeffValid'] = self.coeffValid.reshape(ctLength).tolist()
        if self.coeffsSqr is not None:
            outDict['coeffsSqr'] = self.coeffsSqr.reshape(ctLength).tolist()
        if self.coeffErrSqr is not None:
            outDict['coeffErrSqr'] = self.coeffErrSqr.reshape(ctLength).tolist()
        if self.ampGainRatios is not None:
            outDict['ampGainRatios'] = self.ampGainRatios.reshape(ctLength).tolist()

        if self.interChip:
            outDict['interChip'] = dict()
            for detector in self.interChip:
                outDict['interChip'][detector] = self.interChip[detector].reshape(ctLength).tolist()

        return outDict

    @classmethod
    def fromTable(cls, tableList):
        """Construct calibration from a list of tables.

        This method uses the `fromDict` method to create the
        calibration, after constructing an appropriate dictionary from
        the input tables.

        Parameters
        ----------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables to use to construct the crosstalk
            calibration.

        Returns
        -------
        calib : `lsst.ip.isr.CrosstalkCalib`
            The calibration defined in the tables.

        """
        coeffTable = tableList[0]

        metadata = coeffTable.meta
        inDict = dict()
        inDict['metadata'] = metadata
        inDict['hasCrosstalk'] = metadata['HAS_CROSSTALK']
        inDict['nAmp'] = metadata['NAMP']
        calibVersion = metadata['CROSSTALK_VERSION']
        if calibVersion < 1.1:
            inDict['crosstalkRatiosUnits'] = ''
        else:
            inDict['crosstalkRatiosUnits'] = metadata['CROSSTALK_RATIOS_UNITS']
        inDict['coeffs'] = coeffTable['CT_COEFFS']
        if 'CT_ERRORS' in coeffTable.columns:
            inDict['coeffErr'] = coeffTable['CT_ERRORS']
        if 'CT_COUNTS' in coeffTable.columns:
            inDict['coeffNum'] = coeffTable['CT_COUNTS']
        if 'CT_VALID' in coeffTable.columns:
            inDict['coeffValid'] = coeffTable['CT_VALID']
        if 'CT_COEFFS_SQR' in coeffTable.columns:
            inDict['coeffsSqr'] = coeffTable['CT_COEFFS_SQR']
        if 'CT_ERRORS_SQR' in coeffTable.columns:
            inDict['coeffErrSqr'] = coeffTable['CT_ERRORS_SQR']
        if 'CT_AMP_GAIN_RATIOS' in coeffTable.columns:
            inDict['ampGainRatios'] = coeffTable['CT_AMP_GAIN_RATIOS']

        if len(tableList) > 1:
            inDict['interChip'] = dict()
            interChipTable = tableList[1]
            for record in interChipTable:
                inDict['interChip'][record['IC_SOURCE_DET']] = record['IC_COEFFS']

        return cls().fromDict(inDict)

    def toTable(self):
        """Construct a list of tables containing the information in this
        calibration.

        The list of tables should create an identical calibration
        after being passed to this class's fromTable method.

        Returns
        -------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables containing the crosstalk calibration
            information.

        """
        tableList = []
        self.updateMetadata()
        catalog = Table([{'CT_COEFFS': self.coeffs.reshape(self.nAmp*self.nAmp),
                          'CT_ERRORS': self.coeffErr.reshape(self.nAmp*self.nAmp),
                          'CT_COUNTS': self.coeffNum.reshape(self.nAmp*self.nAmp),
                          'CT_VALID': self.coeffValid.reshape(self.nAmp*self.nAmp),
                          'CT_COEFFS_SQR': self.coeffsSqr.reshape(self.nAmp*self.nAmp),
                          'CT_ERRORS_SQR': self.coeffErrSqr.reshape(self.nAmp*self.nAmp),
                          'CT_AMP_GAIN_RATIOS': self.ampGainRatios.reshape(self.nAmp*self.nAmp),
                          }])
        # filter None, because astropy can't deal.
        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        catalog.meta = outMeta
        tableList.append(catalog)

        if self.interChip:
            interChipTable = Table([{'IC_SOURCE_DET': sourceDet,
                                     'IC_COEFFS': self.interChip[sourceDet].reshape(self.nAmp*self.nAmp)}
                                   for sourceDet in self.interChip.keys()])
            tableList.append(interChipTable)
        return tableList

    # Implementation methods.
    @staticmethod
    def extractAmp(image, amp, ampTarget, isTrimmed=False):
        """Extract the image data from an amp, flipped to match ampTarget.

        Parameters
        ----------
        image : `lsst.afw.image.Image` or `lsst.afw.image.MaskedImage`
            Image containing the amplifier of interest.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Amplifier on image to extract.
        ampTarget : `lsst.afw.cameraGeom.Amplifier`
            Target amplifier that the extracted image will be flipped
            to match.
        isTrimmed : `bool`
            The image is already trimmed.
            TODO : DM-15409 will resolve this.

        Returns
        -------
        output : `lsst.afw.image.Image`
            Image of the amplifier in the desired configuration.
        """
        X_FLIP = {lsst.afw.cameraGeom.ReadoutCorner.LL: False,
                  lsst.afw.cameraGeom.ReadoutCorner.LR: True,
                  lsst.afw.cameraGeom.ReadoutCorner.UL: False,
                  lsst.afw.cameraGeom.ReadoutCorner.UR: True}
        Y_FLIP = {lsst.afw.cameraGeom.ReadoutCorner.LL: False,
                  lsst.afw.cameraGeom.ReadoutCorner.LR: False,
                  lsst.afw.cameraGeom.ReadoutCorner.UL: True,
                  lsst.afw.cameraGeom.ReadoutCorner.UR: True}

        output = image[amp.getBBox() if isTrimmed else amp.getRawDataBBox()]
        thisAmpCorner = amp.getReadoutCorner()
        targetAmpCorner = ampTarget.getReadoutCorner()

        # Flipping is necessary only if the desired configuration doesn't match
        # what we currently have.
        xFlip = X_FLIP[targetAmpCorner] ^ X_FLIP[thisAmpCorner]
        yFlip = Y_FLIP[targetAmpCorner] ^ Y_FLIP[thisAmpCorner]
        return lsst.afw.math.flipImage(output, xFlip, yFlip)

    @staticmethod
    def calculateBackground(mi, badPixels=["BAD"]):
        """Estimate median background in image.

        Getting a great background model isn't important for crosstalk
        correction, since the crosstalk is at a low level. The median should
        be sufficient.

        Parameters
        ----------
        mi : `lsst.afw.image.MaskedImage`
            MaskedImage for which to measure background.
        badPixels : `list` of `str`
            Mask planes to ignore.
        Returns
        -------
        bg : `float`
            Median background level.
        """
        mask = mi.getMask()
        stats = lsst.afw.math.StatisticsControl()
        stats.setAndMask(mask.getPlaneBitMask(badPixels))
        return lsst.afw.math.makeStatistics(mi, lsst.afw.math.MEDIAN, stats).getValue()

    def subtractCrosstalk(self, thisExposure, sourceExposure=None, crosstalkCoeffs=None,
                          crosstalkCoeffsSqr=None,
                          badPixels=["BAD"], minPixelToMask=45000,
                          crosstalkStr="CROSSTALK", isTrimmed=False,
                          backgroundMethod="None", doSqrCrosstalk=False):
        """Subtract the crosstalk from thisExposure, optionally using a
        different source.

        We set the mask plane indicated by ``crosstalkStr`` in a target
        amplifier for pixels in a source amplifier that exceed
        ``minPixelToMask``. Note that the correction is applied to all pixels
        in the amplifier, but only those that have a substantial crosstalk
        are masked with ``crosstalkStr``.

        The uncorrected image is used as a template for correction. This is
        good enough if the crosstalk is small (e.g., coefficients < ~ 1e-3),
        but if it's larger you may want to iterate.

        Parameters
        ----------
        thisExposure : `lsst.afw.image.Exposure`
            Exposure for which to subtract crosstalk.
        sourceExposure : `lsst.afw.image.Exposure`, optional
            Exposure to use as the source of the crosstalk.  If not set,
            thisExposure is used as the source (intra-detector crosstalk).
        crosstalkCoeffs : `numpy.ndarray`, optional.
            Coefficients to use to correct crosstalk.
        crosstalkCoeffsSqr : `numpy.ndarray`, optional.
            Quadratic coefficients to use to correct crosstalk.
        badPixels : `list` of `str`, optional
            Mask planes to ignore.
        minPixelToMask : `float`, optional
            Minimum pixel value (relative to the background level) in
            source amplifier for which to set ``crosstalkStr`` mask plane
            in target amplifier.
        crosstalkStr : `str`, optional
            Mask plane name for pixels greatly modified by crosstalk
            (above minPixelToMask).
        isTrimmed : `bool`, optional
            The image is already trimmed.
            This should no longer be needed once DM-15409 is resolved.
        backgroundMethod : `str`, optional
            Method used to subtract the background.  "AMP" uses
            amplifier-by-amplifier background levels, "DETECTOR" uses full
            exposure/maskedImage levels.  Any other value results in no
            background subtraction.
        doSqrCrosstalk: `bool`, optional
            Should the quadratic crosstalk coefficients be used for the
            crosstalk correction?

        Notes
        -----

        For a given image I, we want to find the crosstalk subtrahend
        image CT, such that
                         I_corrected = I - CT
        The subtrahend image is the sum of all crosstalk contributions
        that appear in I, so we can build it up by amplifier. Each
        amplifier A in image I sees the contributions from all other
        amplifiers B_v != A.  For the current linear model, we set `sImage`
        equal to the segment of the subtrahend image CT corresponding to
        amplifier A, and then build it up as:
        simage_linear = sum_v coeffsA_v * (B_v - bkg_v) where coeffsA_v
        is the vector of crosstalk coefficients for sources that cause
        images in amplifier A.  The bkg_v term in this equation is
        identically 0.0 for all cameras except obs_subaru (and is only
        non-zero there for historical reasons).
        To include the non-linear term, we can again add to the subtrahend
        image using the same loop, as:

        simage_nonlinear = sum_v (coeffsA_v * B_v) + (NLcoeffsA_v * B_v * B_v)
                         = sum_v   linear_term_v   +    nonlinear_term_v

        where coeffsA_v is the linear term, and NLcoeffsA_v are the quadratic
        component. For LSSTCam, it has been observed that the linear_term_v >>
        nonlinear_term_v.
        """
        mi = thisExposure.getMaskedImage()
        mask = mi.getMask()
        detector = thisExposure.getDetector()
        if self.hasCrosstalk is False:
            self.fromDetector(detector, coeffVector=crosstalkCoeffs)

        numAmps = len(detector)
        if numAmps != self.nAmp:
            raise RuntimeError(f"Crosstalk built for {self.nAmp} in {self._detectorName}, received "
                               f"{numAmps} in {detector.getName()}")

        if doSqrCrosstalk and crosstalkCoeffsSqr is None:
            raise RuntimeError("Attempted to perform NL crosstalk correction without NL "
                               "crosstalk coefficients.")

        if sourceExposure:
            source = sourceExposure.getMaskedImage()
            sourceDetector = sourceExposure.getDetector()
        else:
            source = mi
            sourceDetector = detector

        if crosstalkCoeffs is not None:
            coeffs = crosstalkCoeffs
        else:
            coeffs = self.coeffs
        self.log.debug("CT COEFF: %s", coeffs)

        if doSqrCrosstalk:
            if crosstalkCoeffsSqr is not None:
                coeffsSqr = crosstalkCoeffsSqr
            else:
                coeffsSqr = self.coeffsSqr
            self.log.debug("CT COEFF SQR: %s", coeffsSqr)
        # Set background level based on the requested method.  The
        # thresholdBackground holds the offset needed so that we only mask
        # pixels high relative to the background, not in an absolute
        # sense.
        thresholdBackground = self.calculateBackground(source, badPixels)

        backgrounds = [0.0 for amp in sourceDetector]
        if backgroundMethod is None:
            pass
        elif backgroundMethod == "AMP":
            backgrounds = [self.calculateBackground(source[amp.getBBox()], badPixels)
                           for amp in sourceDetector]
        elif backgroundMethod == "DETECTOR":
            backgrounds = [self.calculateBackground(source, badPixels) for amp in sourceDetector]

        # Set the crosstalkStr bit for the bright pixels (those which will have
        # significant crosstalk correction)
        crosstalkPlane = mask.addMaskPlane(crosstalkStr)
        footprints = lsst.afw.detection.FootprintSet(source,
                                                     lsst.afw.detection.Threshold(minPixelToMask
                                                                                  + thresholdBackground))
        footprints.setMask(mask, crosstalkStr)
        crosstalk = mask.getPlaneBitMask(crosstalkStr)

        # Define a subtrahend image to contain all the scaled crosstalk signals
        subtrahend = source.Factory(source.getBBox())
        subtrahend.set((0, 0, 0))

        coeffs = coeffs.transpose()
        # Apply NL coefficients
        if doSqrCrosstalk:
            coeffsSqr = coeffsSqr.transpose()
            mi2 = mi.clone()
            mi2.scaledMultiplies(1.0, mi)
        for ss, sAmp in enumerate(sourceDetector):
            sImage = subtrahend[sAmp.getBBox() if isTrimmed else sAmp.getRawDataBBox()]
            for tt, tAmp in enumerate(detector):
                if coeffs[ss, tt] == 0.0:
                    continue
                tImage = self.extractAmp(mi, tAmp, sAmp, isTrimmed)
                tImage.getMask().getArray()[:] &= crosstalk  # Remove all other masks
                tImage -= backgrounds[tt]
                sImage.scaledPlus(coeffs[ss, tt], tImage)
                # Add the nonlinear term
                if doSqrCrosstalk:
                    tImageSqr = self.extractAmp(mi2, tAmp, sAmp, isTrimmed)
                    sImage.scaledPlus(coeffsSqr[ss, tt], tImageSqr)

        # Set crosstalkStr bit only for those pixels that have been
        # significantly modified (i.e., those masked as such in 'subtrahend'),
        # not necessarily those that are bright originally.
        mask.clearMaskPlane(crosstalkPlane)
        mi -= subtrahend  # also sets crosstalkStr bit for bright pixels

    def subtractCrosstalkParallelOverscanRegion(self, thisExposure, crosstalkCoeffs=None,
                                                crosstalkCoeffsSqr=None,
                                                badPixels=["BAD"], crosstalkStr="CROSSTALK",
                                                detectorConfig=None, doSqrCrosstalk=False):
        """Subtract crosstalk just from the parallel overscan region.

        This assumes that serial overscan has been previously subtracted.

        Parameters
        ----------
        thisExposure : `lsst.afw.image.Exposure`
            Exposure for which to subtract crosstalk.
        crosstalkCoeffs : `numpy.ndarray`, optional.
            Coefficients to use to correct crosstalk.
        crosstalkCoeffsSqr : `numpy.ndarray`, optional.
            Quadratic coefficients to use to correct crosstalk.
        badPixels : `list` of `str`, optional
            Mask planes to ignore.
        crosstalkStr : `str`, optional
            Mask plane name for pixels greatly modified by crosstalk
            (above minPixelToMask).
        detectorConfig : `lsst.ip.isr.overscanDetectorConfig`, optional
            Per-amplifier configs to use.
        doSqrCrosstalk: `bool`, optional
            Should the quadratic crosstalk coefficients be used for the
            crosstalk correction?
        """
        mi = thisExposure.getMaskedImage()
        mask = mi.getMask()
        detector = thisExposure.getDetector()
        if self.hasCrosstalk is False:
            self.fromDetector(detector, coeffVector=crosstalkCoeffs)

        numAmps = len(detector)
        if numAmps != self.nAmp:
            raise RuntimeError(f"Crosstalk built for {self.nAmp} in {self._detectorName}, received "
                               f"{numAmps} in {detector.getName()}")

        if doSqrCrosstalk and crosstalkCoeffsSqr is None:
            raise RuntimeError("Attempted to perform NL crosstalk correction without NL "
                               "crosstalk coefficients.")

        source = mi
        sourceDetector = detector

        if crosstalkCoeffs is not None:
            coeffs = crosstalkCoeffs
        else:
            coeffs = self.coeffs
        if doSqrCrosstalk:
            if crosstalkCoeffsSqr is not None:
                coeffsSqr = crosstalkCoeffsSqr
            else:
                coeffsSqr = self.coeffsSqr
            self.log.debug("CT COEFF SQR: %s", coeffsSqr)

        crosstalkPlane = mask.addMaskPlane(crosstalkStr)
        crosstalk = mask.getPlaneBitMask(crosstalkStr)

        subtrahend = source.Factory(source.getBBox())
        subtrahend.set((0, 0, 0))

        coeffs = coeffs.transpose()
        # Apply NL coefficients
        if doSqrCrosstalk:
            coeffsSqr = coeffsSqr.transpose()
            mi2 = mi.clone()
            mi2.scaledMultiplies(1.0, mi)
        for ss, sAmp in enumerate(sourceDetector):
            if detectorConfig is not None:
                ampConfig = detectorConfig.getOverscanAmpconfig(sAmp.getName())
                if not ampConfig.doParallelOverscanCrosstalk:
                    # Skip crosstalk correction for this amplifier.
                    continue

            sImage = subtrahend[sAmp.getRawParallelOverscanBBox()]
            for tt, tAmp in enumerate(detector):
                if coeffs[ss, tt] == 0.0:
                    continue
                tImage = self.extractAmp(mi, tAmp, sAmp, False, parallelOverscan=True)
                tImage.getMask().getArray()[:] &= crosstalk  # Remove all other masks
                sImage.scaledPlus(coeffs[ss, tt], tImage)
                # Add the nonlinear term, if any.
                if doSqrCrosstalk:
                    tImageSqr = self.extractAmp(mi2, tAmp, sAmp, False, parallelOverscan=True)
                    sImage.scaledPlus(coeffsSqr[ss, tt], tImageSqr)
        # Set crosstalkStr bit only for those pixels that have been
        # significantly modified (i.e., those masked as such in 'subtrahend'),
        # not necessarily those that are bright originally.
        mask.clearMaskPlane(crosstalkPlane)
        mi -= subtrahend  # also sets crosstalkStr bit for bright pixels


class CrosstalkConfig(Config):
    """Configuration for intra-detector crosstalk removal."""
    minPixelToMask = Field(
        dtype=float,
        doc="Set crosstalk mask plane for pixels over this value.",
        default=45000
    )
    crosstalkMaskPlane = Field(
        dtype=str,
        doc="Name for crosstalk mask plane.",
        default="CROSSTALK"
    )
    crosstalkBackgroundMethod = ChoiceField(
        dtype=str,
        doc="Type of background subtraction to use when applying correction.",
        default="None",
        allowed={
            "None": "Do no background subtraction.",
            "AMP": "Subtract amplifier-by-amplifier background levels.",
            "DETECTOR": "Subtract detector level background."
        },
    )
    useConfigCoefficients = Field(
        dtype=bool,
        doc="Ignore the detector crosstalk information in favor of CrosstalkConfig values?",
        default=False,
    )
    crosstalkValues = ListField(
        dtype=float,
        doc=("Amplifier-indexed crosstalk coefficients to use.  This should be arranged as a 1 x nAmp**2 "
             "list of coefficients, such that when reshaped by crosstalkShape, the result is nAmp x nAmp. "
             "This matrix should be structured so CT * [amp0 amp1 amp2 ...]^T returns the column "
             "vector [corr0 corr1 corr2 ...]^T."),
        default=[0.0],
    )
    crosstalkShape = ListField(
        dtype=int,
        doc="Shape of the coefficient array.  This should be equal to [nAmp, nAmp].",
        default=[1],
    )
    doQuadraticCrosstalkCorrection = Field(
        dtype=bool,
        doc="Use quadratic crosstalk coefficients in the crosstalk correction",
        default=False,
    )

    def getCrosstalk(self, detector=None):
        """Return a 2-D numpy array of crosstalk coefficients in the proper
        shape.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.detector`
            Detector that is to be crosstalk corrected.

        Returns
        -------
        coeffs : `numpy.ndarray`
            Crosstalk coefficients that can be used to correct the detector.

        Raises
        ------
        RuntimeError
            Raised if no coefficients could be generated from this
            detector/configuration.
        """
        if self.useConfigCoefficients is True:
            coeffs = np.array(self.crosstalkValues).reshape(self.crosstalkShape)
            if detector is not None:
                nAmp = len(detector)
                if coeffs.shape != (nAmp, nAmp):
                    raise RuntimeError("Constructed crosstalk coeffients do not match detector shape. "
                                       f"{coeffs.shape} {nAmp}")
            return coeffs
        elif detector is not None and detector.hasCrosstalk() is True:
            # Assume the detector defines itself consistently.
            return detector.getCrosstalk()
        else:
            raise RuntimeError("Attempted to correct crosstalk without crosstalk coefficients")

    def hasCrosstalk(self, detector=None):
        """Return a boolean indicating if crosstalk coefficients exist.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.detector`
            Detector that is to be crosstalk corrected.

        Returns
        -------
        hasCrosstalk : `bool`
            True if this detector/configuration has crosstalk coefficients
            defined.
        """
        if self.useConfigCoefficients is True and self.crosstalkValues is not None:
            return True
        elif detector is not None and detector.hasCrosstalk() is True:
            return True
        else:
            return False


class CrosstalkTask(Task):
    """Apply intra-detector crosstalk correction."""
    ConfigClass = CrosstalkConfig
    _DefaultName = 'isrCrosstalk'

    def run(self,
            exposure, crosstalk=None,
            crosstalkSources=None, isTrimmed=False, camera=None, parallelOverscanRegion=False,
            detectorConfig=None,
            ):
        """Apply intra-detector crosstalk correction

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure for which to remove crosstalk.
        crosstalkCalib : `lsst.ip.isr.CrosstalkCalib`, optional
            External crosstalk calibration to apply.  Constructed from
            detector if not found.
        crosstalkSources : `defaultdict`, optional
            Image data for other detectors that are sources of
            crosstalk in exposure.  The keys are expected to be names
            of the other detectors, with the values containing
            `lsst.afw.image.Exposure` at the same level of processing
            as ``exposure``.
            The default for intra-detector crosstalk here is None.
        isTrimmed : `bool`, optional
            The image is already trimmed.
            This should no longer be needed once DM-15409 is resolved.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Camera associated with this exposure.  Only used for
            inter-chip matching.
        parallelOverscanRegion : `bool`, optional
            Do subtraction in parallel overscan region (only)?
        detectorConfig : `lsst.ip.isr.OverscanDetectorConfig`, optional
            Per-amplifier configs used when parallelOverscanRegion=True.

        Raises
        ------
        RuntimeError
            Raised if called for a detector that does not have a
            crosstalk correction.  Also raised if the crosstalkSource
            is not an expected type.
        """
        if not crosstalk:
            crosstalk = CrosstalkCalib(log=self.log)
            crosstalk = crosstalk.fromDetector(exposure.getDetector(),
                                               coeffVector=self.config.crosstalkValues)
        if not crosstalk.log:
            crosstalk.log = self.log

        doSqrCrosstalk = self.config.doQuadraticCrosstalkCorrection
        if doSqrCrosstalk and crosstalk.coeffsSqr is None:
            raise RuntimeError("Attempted to perform NL crosstalk correction without NL "
                               "crosstalk coefficients.")
        if doSqrCrosstalk:
            crosstalkCoeffsSqr = crosstalk.coeffsSqr
        else:
            crosstalkCoeffsSqr = None

        if not crosstalk.hasCrosstalk:
            raise RuntimeError("Attempted to correct crosstalk without crosstalk coefficients.")
        elif parallelOverscanRegion:
            self.log.info("Applying crosstalk correction to parallel overscan region.")
            crosstalk.subtractCrosstalkParallelOverscanRegion(
                exposure,
                crosstalkCoeffs=crosstalk.coeffs,
                crosstalkCoeffsSqr=crosstalkCoeffsSqr,
                detectorConfig=detectorConfig,
                doSqrCrosstalk=doSqrCrosstalk,
            )
        else:
            self.log.info("Applying crosstalk correction.")
            crosstalk.subtractCrosstalk(exposure, crosstalkCoeffs=crosstalk.coeffs,
                                        crosstalkCoeffsSqr=crosstalkCoeffsSqr,
                                        minPixelToMask=self.config.minPixelToMask,
                                        crosstalkStr=self.config.crosstalkMaskPlane, isTrimmed=isTrimmed,
                                        backgroundMethod=self.config.crosstalkBackgroundMethod,
                                        doSqrCrosstalk=doSqrCrosstalk)

            if crosstalk.interChip:
                if crosstalkSources:
                    # Parse crosstalkSources: Identify which detectors we have
                    # available
                    if isinstance(crosstalkSources[0], lsst.afw.image.Exposure):
                        # Received afwImage.Exposure
                        sourceNames = [exp.getDetector().getName() for exp in crosstalkSources]
                    elif isinstance(crosstalkSources[0], lsst.daf.butler.DeferredDatasetHandle):
                        # Received dafButler.DeferredDatasetHandle
                        detectorList = [source.dataId['detector'] for source in crosstalkSources]
                        sourceNames = [camera[detector].getName() for detector in detectorList]
                    else:
                        raise RuntimeError("Unknown object passed as crosstalk sources.",
                                           type(crosstalkSources[0]))

                    for detName in crosstalk.interChip:
                        if detName not in sourceNames:
                            self.log.warning("Crosstalk lists %s, not found in sources: %s",
                                             detName, sourceNames)
                            continue
                        # Get the coefficients.
                        interChipCoeffs = crosstalk.interChip[detName]

                        sourceExposure = crosstalkSources[sourceNames.index(detName)]
                        if isinstance(sourceExposure, lsst.daf.butler.DeferredDatasetHandle):
                            # Dereference the dafButler.DeferredDatasetHandle.
                            sourceExposure = sourceExposure.get()
                        if not isinstance(sourceExposure, lsst.afw.image.Exposure):
                            raise RuntimeError("Unknown object passed as crosstalk sources.",
                                               type(sourceExposure))

                        self.log.info("Correcting detector %s with ctSource %s",
                                      exposure.getDetector().getName(),
                                      sourceExposure.getDetector().getName())
                        crosstalk.subtractCrosstalk(exposure, sourceExposure=sourceExposure,
                                                    crosstalkCoeffs=interChipCoeffs,
                                                    minPixelToMask=self.config.minPixelToMask,
                                                    crosstalkStr=self.config.crosstalkMaskPlane,
                                                    isTrimmed=isTrimmed,
                                                    backgroundMethod=self.config.crosstalkBackgroundMethod)
                else:
                    self.log.warning("Crosstalk contains interChip coefficients, but no sources found!")


class NullCrosstalkTask(CrosstalkTask):
    def run(self, exposure, crosstalkSources=None):
        self.log.info("Not performing any crosstalk correction")
