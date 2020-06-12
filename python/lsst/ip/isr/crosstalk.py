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
import numpy as np
from astropy.table import Table

import lsst.afw.math
import lsst.afw.detection
from lsst.pex.config import Config, Field, ChoiceField, ListField
from lsst.pipe.base import Task

from lsst.ip.isr import IsrCalib


__all__ = ["CrosstalkCalib", "CrosstalkConfig", "CrosstalkTask",
           "NullCrosstalkTask"]


class CrosstalkCalib(IsrCalib):
    """Calibration of amp-to-amp crosstalk coefficients.

    Parameters
    ----------
    detector : `lsst.afw.cameraGeom.Detector`, optional
        Detector to use to pull coefficients from.
    log : `lsst.log.Log`, optional
        Log to write messages to.
    **kwargs :
        Parameters to pass to parent constructor.
    """
    _OBSTYPE = 'CROSSTALK'
    _SCHEMA = 'Gen3 Crosstalk'
    _VERSION = 1.0

    def __init__(self, detector=None, **kwargs):
        self.hasCrosstalk = False
        self.nAmp = 0
        self.crosstalkShape = (self.nAmp, self.nAmp)

        self.coeffs = None
        self.coeffErr = None
        self.coeffNum = None
        self.interChip = None

        super().__init__(**kwargs)
        self.requiredAttributes.update(['hasCrosstalk', 'nAmp', 'coeffs',
                                        'coeffErr', 'coeffNum', 'interChip'])
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
        kwargs['DETECTOR'] = self._detectorName
        kwargs['DETECTOR_SERIAL'] = self._detectorSerial
        kwargs['HAS_CROSSTALK'] = self.hasCrosstalk
        kwargs['NAMP'] = self.nAmp
        self.crosstalkShape = (self.nAmp, self.nAmp)
        kwargs['CROSSTALK_SHAPE'] = self.crosstalkShape

        super().updateMetadata(setDate=setDate, **kwargs)

    def fromDetector(self, detector, coeffVector=None):
        """Set calibration parameters from the detector.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            Detector to use to set parameters from.
        coeffVector : `numpy.array`, optional
            Use the detector geometry (bounding boxes and flip
            information), but use ``coeffVector`` instead of the
            output of ``detector.getCrosstalk()``.

        Returns
        -------
        calib : `lsst.ip.isr.CrosstalkCalib`
            The calibration constructed from the detector.

        """
        if detector.hasCrosstalk() or coeffVector:
            self._detectorName = detector.getName()
            self._detectorSerial = detector.getSerial()

            self.nAmp = len(detector)
            self.crosstalkShape = (self.nAmp, self.nAmp)

            if coeffVector is not None:
                crosstalkCoeffs = coeffVector
            else:
                crosstalkCoeffs = detector.getCrosstalk()
            if len(crosstalkCoeffs) == 1 and crosstalkCoeffs[0] == 0.0:
                return self
            self.coeffs = np.array(crosstalkCoeffs).reshape(self.crosstalkShape)

            if self.coeffs.shape != self.crosstalkShape:
                raise RuntimeError("Crosstalk coefficients do not match detector shape. "
                                   f"{self.crosstalkShape} {self.nAmp}")

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
        RuntimeError :
            Raised if the supplied dictionary is for a different
            calibration.
        """
        calib = cls()

        if calib._OBSTYPE != dictionary['metadata']['OBSTYPE']:
            raise RuntimeError(f"Incorrect crosstalk supplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['metadata']['OBSTYPE']}")

        calib.setMetadata(dictionary['metadata'])
        calib._detectorName = dictionary['metadata']['DETECTOR']
        calib._detectorSerial = dictionary['metadata']['DETECTOR_SERIAL']

        calib.hasCrosstalk = dictionary.get('hasCrosstalk',
                                            dictionary['metadata'].get('HAS_CROSSTALK', False))
        if calib.hasCrosstalk:
            calib.nAmp = dictionary.get('nAmp', dictionary['metadata'].get('NAMP', 0))
            calib.crosstalkShape = (calib.nAmp, calib.nAmp)
            calib.coeffs = np.array(dictionary['coeffs']).reshape(calib.crosstalkShape)

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

        ctLength = self.nAmp*self.nAmp
        outDict['coeffs'] = self.coeffs.reshape(ctLength).tolist()

        if self.coeffErr is not None:
            outDict['coeffErr'] = self.coeffErr.reshape(ctLength).tolist()
        if self.coeffNum is not None:
            outDict['coeffNum'] = self.coeffNum.reshape(ctLength).tolist()

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

        inDict['coeffs'] = coeffTable['CT_COEFFS']
        if len(tableList) > 1:
            inDict['interChip'] = dict()
            interChipTable = tableList[1]
            for record in interChipTable:
                inDict['interChip'][record['IC_SOURCE_DET']] = record['IC_COEFFS']

        return cls().fromDict(inDict)

    def toTable(self):
        """Construct a list of tables containing the information in this calibration.

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
        catalog = Table([{'CT_COEFFS': self.coeffs.reshape(self.nAmp*self.nAmp)}])
        catalog.meta = self.getMetadata().toDict()
        tableList.append(catalog)

        if self.interChip:
            interChipTable = Table([{'IC_SOURCE_DET': sourceDet,
                                     'IC_COEFFS': self.interChip[sourceDet].reshape(self.nAmp*self.nAmp)}
                                    for sourceDet in self.interChip.keys()])
            tableList.append(interChipTable)
        return tableList

    # Implementation methods.
    def extractAmp(self, image, amp, ampTarget, isTrimmed=False):
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

        # Flipping is necessary only if the desired configuration doesn't match what we currently have
        xFlip = X_FLIP[targetAmpCorner] ^ X_FLIP[thisAmpCorner]
        yFlip = Y_FLIP[targetAmpCorner] ^ Y_FLIP[thisAmpCorner]
        self.log.debug("Extract amp: %s %s %s %s",
                       amp.getName(), ampTarget.getName(), thisAmpCorner, targetAmpCorner)
        return lsst.afw.math.flipImage(output, xFlip, yFlip)

    def calculateBackground(self, mi, badPixels=["BAD"]):
        """Estimate median background in image.

        Getting a great background model isn't important for crosstalk correction,
        since the crosstalk is at a low level. The median should be sufficient.

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
                          badPixels=["BAD"], minPixelToMask=45000,
                          crosstalkStr="CROSSTALK", isTrimmed=False,
                          backgroundMethod="None"):
        """Subtract the crosstalk from thisExposure, optionally using a different source.

        We set the mask plane indicated by ``crosstalkStr`` in a target amplifier
        for pixels in a source amplifier that exceed ``minPixelToMask``. Note that
        the correction is applied to all pixels in the amplifier, but only those
        that have a substantial crosstalk are masked with ``crosstalkStr``.

        The uncorrected image is used as a template for correction. This is good
        enough if the crosstalk is small (e.g., coefficients < ~ 1e-3), but if it's
        larger you may want to iterate.

        Parameters
        ----------
        thisExposure : `lsst.afw.image.Exposure`
            Exposure for which to subtract crosstalk.
        sourceExposure : `lsst.afw.image.Exposure`, optional
            Exposure to use as the source of the crosstalk.  If not set,
            thisExposure is used as the source (intra-detector crosstalk).
        crosstalkCoeffs : `numpy.ndarray`, optional.
            Coefficients to use to correct crosstalk.
        badPixels : `list` of `str`
            Mask planes to ignore.
        minPixelToMask : `float`
            Minimum pixel value (relative to the background level) in
            source amplifier for which to set ``crosstalkStr`` mask plane
            in target amplifier.
        crosstalkStr : `str`
            Mask plane name for pixels greatly modified by crosstalk
            (above minPixelToMask).
        isTrimmed : `bool`
            The image is already trimmed.
            This should no longer be needed once DM-15409 is resolved.
        backgroundMethod : `str`
            Method used to subtract the background.  "AMP" uses
            amplifier-by-amplifier background levels, "DETECTOR" uses full
            exposure/maskedImage levels.  Any other value results in no
            background subtraction.
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
        for ii, iAmp in enumerate(sourceDetector):
            iImage = subtrahend[iAmp.getBBox() if isTrimmed else iAmp.getRawDataBBox()]
            for jj, jAmp in enumerate(detector):
                if coeffs[ii, jj] == 0.0:
                    continue
                jImage = self.extractAmp(mi, jAmp, iAmp, isTrimmed)
                jImage.getMask().getArray()[:] &= crosstalk  # Remove all other masks
                jImage -= backgrounds[jj]
                iImage.scaledPlus(coeffs[ii, jj], jImage)

        # Set crosstalkStr bit only for those pixels that have been significantly modified (i.e., those
        # masked as such in 'subtrahend'), not necessarily those that are bright originally.
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

    def getCrosstalk(self, detector=None):
        """Return a 2-D numpy array of crosstalk coefficients in the proper shape.

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
            Raised if no coefficients could be generated from this detector/configuration.
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
            True if this detector/configuration has crosstalk coefficients defined.
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

    def prepCrosstalk(self, dataRef, crosstalk=None):
        """Placeholder for crosstalk preparation method, e.g., for inter-detector crosstalk.

        Parameters
        ----------
        dataRef : `daf.persistence.butlerSubset.ButlerDataRef`
            Butler reference of the detector data to be processed.
        crosstalk : `~lsst.ip.isr.CrosstalkConfig`
            Crosstalk calibration that will be used.

        See also
        --------
        lsst.obs.decam.crosstalk.DecamCrosstalkTask.prepCrosstalk
        """
        return

    def run(self, exposure, crosstalk=None,
            crosstalkSources=None, isTrimmed=False):
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
        isTrimmed : `bool`
            The image is already trimmed.
            This should no longer be needed once DM-15409 is resolved.

        Raises
        ------
        RuntimeError
            Raised if called for a detector that does not have a
            crosstalk correction.
        """
        if not crosstalk:
            crosstalk = CrosstalkCalib(log=self.log)
            crosstalk = crosstalk.fromDetector(exposure.getDetector(),
                                               coeffVector=self.config.crosstalkValues)
        if not crosstalk.log:
            crosstalk.log = self.log
        if not crosstalk.hasCrosstalk:
            raise RuntimeError("Attempted to correct crosstalk without crosstalk coefficients.")

        else:
            self.log.info("Applying crosstalk correction.")
            crosstalk.subtractCrosstalk(exposure, crosstalkCoeffs=crosstalk.coeffs,
                                        minPixelToMask=self.config.minPixelToMask,
                                        crosstalkStr=self.config.crosstalkMaskPlane, isTrimmed=isTrimmed,
                                        backgroundMethod=self.config.crosstalkBackgroundMethod)

            if crosstalk.interChip:
                if crosstalkSources:
                    for detName in crosstalk.interChip:
                        if isinstance(crosstalkSources[0], 'lsst.afw.image.Exposure'):
                            # Received afwImage.Exposure
                            sourceNames = [exp.getDetector().getName() for exp in crosstalkSources]
                        else:
                            # Received dafButler.DeferredDatasetHandle
                            sourceNames = [expRef.get(datasetType='isrOscanCorr').getDetector().getName()
                                           for expRef in crosstalkSources]
                        if detName not in sourceNames:
                            self.log.warn("Crosstalk lists %s, not found in sources: %s",
                                          detName, sourceNames)
                            continue
                        interChipCoeffs = crosstalk.interChip[detName]
                        sourceExposure = crosstalkSources[sourceNames.index(detName)]
                        crosstalk.subtractCrosstalk(exposure, sourceExposure=sourceExposure,
                                                    crosstalkCoeffs=interChipCoeffs,
                                                    minPixelToMask=self.config.minPixelToMask,
                                                    crosstalkStr=self.config.crosstalkMaskPlane,
                                                    isTrimmed=isTrimmed,
                                                    backgroundMethod=self.config.crosstalkBackgroundMethod)
                else:
                    self.log.warn("Crosstalk contains interChip coefficients, but no sources found!")


class NullCrosstalkTask(CrosstalkTask):
    def run(self, exposure, crosstalkSources=None):
        self.log.info("Not performing any crosstalk correction")
