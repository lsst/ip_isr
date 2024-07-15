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

__all__ = ["IsrMockLSSTConfig", "IsrMockLSST", "RawMockLSST",
           "CalibratedRawMockLSST", "ReferenceMockLSST",
           "BiasMockLSST", "DarkMockLSST", "FlatMockLSST", "FringeMockLSST",
           "BfKernelMockLSST", "DefectMockLSST", "CrosstalkCoeffMockLSST",
           "TransmissionMockLSST"]
import numpy as np

import lsst.geom as geom
import lsst.pex.config as pexConfig
from .crosstalk import CrosstalkCalib
from .isrMock import IsrMockConfig, IsrMock


# FIXME: go through the units on this , it should not be so hard.
# And getting it right is *important* now.


class IsrMockLSSTConfig(IsrMockConfig):
    """Configuration parameters for isrMockLSST.
    """
    # Detector parameters and "Exposure" parameters,
    # mostly inherited from IsrMockConfig.
    isLsstLike = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="If True, products have one raw image per amplifier, otherwise, one raw image per detector.",
    )
    calibMode = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Set to true to produce mock calibration products, e.g. combined bias, dark, flat, etc.",
    )
    doAdd2DBias = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add 2D bias residual frame to data.",
    )
    doAddClockInjectedOffset = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add clock-injected offset to data (on-chip bias level).",
    )
    clockInjectedOffsetLevel = pexConfig.Field(
        dtype=float,
        default=8500.0,
        doc="Clock-injected offset (on-chip bias level), in electrons.",
    )
    noise2DBias = pexConfig.Field(
        dtype=float,
        default=2.0,
        doc="Noise (in electrons) to generate a 2D bias residual frame.",
    )
    doAddDarkNoiseOnly = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Add only dark current noise, for testing consistency.",
    )
    doAddParallelOverscanRamp = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add overscan ramp to parallel overscan and data regions.",
    )
    doAddSerialOverscanRamp = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add overscan ramp to serial overscan and data regions.",
    )
    doApplyGain = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add gain to data.",
    )
    doRoundADU = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Round ADU values to nearest integer.",
    )

    def setDefaults(self):
        super().setDefaults()

        self.gain = 1.7
        self.skyLevel = 1700.0  # e-
        self.sourceFlux = [50_000.0]  # e-
        self.overscanScale = 170.0  # e-
        self.biasLevel = 20_000.0  # ADU
        self.doAddCrosstalk = True


class IsrMockLSST(IsrMock):
    """Class to generate consistent mock images for ISR testing.
    """
    ConfigClass = IsrMockLSSTConfig
    _DefaultName = "isrMockLSST"

    def __init__(self, **kwargs):
        # cross-talk coeffs, bf kernel are defined in the parent class.
        super().__init__(**kwargs)

    def run(self):
        """Generate a mock ISR product following LSSTCam ISR, and return it.

        Returns
        -------
        image : `lsst.afw.image.Exposure`
            Simulated ISR image with signals added.
        dataProduct :
            Simulated ISR data products.
        None :
            Returned if no valid configuration was found.

        Raises
        ------
        RuntimeError
            Raised if both doGenerateImage and doGenerateData are specified.
        """
        if self.config.doGenerateImage and self.config.doGenerateData:
            raise RuntimeError("Only one of doGenerateImage and doGenerateData may be specified.")
        elif self.config.doGenerateImage:
            return self.makeImage()
        elif self.config.doGenerateData:
            return self.makeData()
        else:
            return None

    def makeImage(self):
        """Generate a simulated ISR LSST image.

        Returns
        -------
        exposure : `lsst.afw.image.Exposure` or `dict`
            Simulated ISR image data.

        Notes
        -----
        This method constructs a "raw" data image.
        """
        exposure = self.getExposure()

        # Set up random number generators for consistency of components,
        # no matter the group that are configured.
        rngSky = np.random.RandomState(seed=self.config.rngSeed + 1)
        rngDark = np.random.RandomState(seed=self.config.rngSeed + 2)
        rng2DBias = np.random.RandomState(seed=self.config.rngSeed + 3)
        rngOverscan = np.random.RandomState(seed=self.config.rngSeed + 4)
        rngReadNoise = np.random.RandomState(seed=self.config.rngSeed + 5)

        # We introduce effects as they happen from a source to the signal,
        # so the effects go from electrons to ADU.
        # The ISR steps will then correct these effects in the reverse order.
        for idx, amp in enumerate(exposure.getDetector()):

            # Get image bbox and data
            bbox = None
            if self.config.isTrimmed:
                bbox = amp.getBBox()
                bboxFull = bbox
            else:
                bbox = amp.getRawDataBBox()
                bboxFull = amp.getRawBBox()

            # This is the image data (excluding pre/overscans).
            ampImageData = exposure.image[bbox]
            # This is the full data (including pre/overscans if untrimmed).
            ampFullData = exposure.image[bboxFull]

            # Astrophysical signals are all in electrons (e-).
            # These are only applied to the imaging portion of the
            # amplifier (ampImageData)

            if self.config.doAddSky:
                # The sky effects are in electrons.
                self.amplifierAddNoise(
                    ampImageData,
                    self.config.skyLevel,
                    np.sqrt(self.config.skyLevel),
                    rng=rngSky,
                )

            if self.config.doAddSource:
                for sourceAmp, sourceFlux, sourceX, sourceY in zip(self.config.sourceAmp,
                                                                   self.config.sourceFlux,
                                                                   self.config.sourceX,
                                                                   self.config.sourceY):
                    if idx == sourceAmp:
                        # The source flux is in electrons,
                        # but the sourceFlux is configured in ADU (*** NOT NOW)
                        # TODO: DM-42880 to set configs to correct units
                        self.amplifierAddSource(ampImageData, sourceFlux, sourceX, sourceY)

            if self.config.doAddFringe:
                # Fringes are added in electrons.
                self.amplifierAddFringe(amp,
                                        ampImageData,
                                        np.array(self.config.fringeScale),
                                        x0=np.array(self.config.fringeX0),
                                        y0=np.array(self.config.fringeY0))

            if self.config.doAddFlat:
                if self.config.calibMode:
                    # In case we are making a combined flat,
                    # add a non-zero signal so the mock flat can be multiplied
                    self.amplifierAddNoise(ampImageData, 1.0, 0.0)
                # Multiply each amplifier by a Gaussian centered on u0 and v0
                u0 = exposure.getDetector().getBBox().getDimensions().getX()/2.
                v0 = exposure.getDetector().getBBox().getDimensions().getY()/2.
                self.amplifierMultiplyFlat(amp, ampImageData, self.config.flatDrop, u0=u0, v0=v0)

            # On-chip electronic effects.

            # TODO
            # 1. Add BF effect (e-)

            # 2. Add dark current (e-) to imaging portion of the amplifier.
            if self.config.doAddDark or self.config.doAddDarkNoiseOnly:
                if self.config.doAddDarkNoiseOnly:
                    darkLevel = 0.0
                else:
                    darkLevel = self.config.darkRate * self.config.darkTime
                if self.config.calibMode:
                    darkNoise = 0.0
                else:
                    darkNoise = np.sqrt(self.config.darkRate * self.config.darkTime)

                self.amplifierAddNoise(ampImageData, darkLevel, darkNoise, rng=rngDark)

            # 3. Add 2D bias residual (e-) to imaging portion of the amplifier.
            if self.config.doAdd2DBias:
                self.amplifierAddNoise(
                    ampImageData,
                    0.0,
                    self.config.noise2DBias,
                    rng=rng2DBias,
                )

            # 3. Add serial CTI (e-) to amplifier (imaging + overscan).
            # TODO

            # 4. Add clock-injected offset (e-) to amplifier
            #    (imaging + overscan).
            # This is just an offset that will be crosstalked and modified by
            # the gain, and does not have a noise associated with it.
            if self.config.doAddClockInjectedOffset:
                self.amplifierAddNoise(
                    ampFullData,
                    self.config.clockInjectedOffsetLevel,
                    0.0,
                )

            # 5. Add non-linearity (e-) to amplifier (imaging + overscan).
            # TODO

            # 6./7. Add serial and parallel overscan slopes (e-)
            #       (overscan regions).
            if (self.config.doAddParallelOverscanRamp or self.config.doAddSerialOverscanRamp) and \
               not self.config.isTrimmed:
                parallelOverscanBBox = amp.getRawParallelOverscanBBox()
                parallelOverscanData = exposure.image[parallelOverscanBBox]

                serialOverscanBBox = self.getFullSerialOverscanBBox(amp)
                serialOverscanData = exposure.image[serialOverscanBBox]

                # Add read noise of mean 0
                # to the parallel and serial overscan regions.
                self.amplifierAddNoise(
                    parallelOverscanData,
                    0.0,
                    self.config.readNoise,
                    rng=rngOverscan,
                )
                self.amplifierAddNoise(
                    serialOverscanData,
                    0.0,
                    self.config.readNoise,
                    rng=rngOverscan,
                )

                if self.config.doAddParallelOverscanRamp:
                    # Apply gradient along the X axis.
                    self.amplifierAddXGradient(ampFullData, -1.0 * self.config.overscanScale,
                                               1.0 * self.config.overscanScale)

                if self.config.doAddSerialOverscanRamp:
                    # Apply the gradient along the Y axis.
                    self.amplifierAddYGradient(ampFullData, -1.0 * self.config.overscanScale,
                                               1.0 * self.config.overscanScale)

            # 8. Add read noise (e-) to the amplifier (imaging + overscan).
            # (Unsure if this should be before or after crosstalk)
            # (Probably some of both)
            # (Probably doesn't matter)
            if not self.config.calibMode:
                self.amplifierAddNoise(
                    ampImageData,
                    0.0,
                    self.config.readNoise,
                    rng=rngReadNoise,
                )

        # 9. Add crosstalk (e-) to all the amplifiers (imaging + overscan).
        if self.config.doAddCrosstalk:
            ctCalib = CrosstalkCalib()
            exposureClean = exposure.clone()
            for idxS, ampS in enumerate(exposure.getDetector()):
                for idxT, ampT in enumerate(exposure.getDetector()):
                    ampDataTarget = exposure.image[ampT.getBBox() if self.config.isTrimmed
                                                   else ampT.getRawBBox()]
                    ampDataSource = ctCalib.extractAmp(exposureClean.image, ampS, ampT,
                                                       isTrimmed=self.config.isTrimmed,
                                                       fullAmplifier=True)
                    self.amplifierAddCT(ampDataSource, ampDataTarget, self.crosstalkCoeffs[idxS][idxT])

        for amp in exposure.getDetector():
            # Get image bbox and data (again).
            bbox = None
            if self.config.isTrimmed:
                bbox = amp.getBBox()
                bboxFull = bbox
            else:
                bbox = amp.getRawDataBBox()
                bboxFull = amp.getRawBBox()

            # This is the image data (excluding pre/overscans).
            ampImageData = exposure.image[bbox]
            # This is the full data (including pre/overscans if untrimmed).
            ampFullData = exposure.image[bboxFull]

            # 10. Gain un-normalize (from e- to floating point ADU)
            if self.config.doApplyGain:
                self.applyGain(ampFullData, self.config.gain)

            # 11. Add overall bias level (ADU) to the amplifier
            #    (imaging + overscan)
            if self.config.doAddBias:
                self.addBiasLevel(ampFullData, self.config.biasLevel)

            # 12. Round/Truncate to integers (ADU)
            if self.config.doRoundADU:
                self.roundADU(ampFullData)

        if self.config.doGenerateAmpDict:
            expDict = dict()
            for amp in exposure.getDetector():
                expDict[amp.getName()] = exposure
            return expDict
        else:
            return exposure

    def addBiasLevel(self, ampData, biasLevel):
        """Add bias level to an amplifier's image data.

        Parameters
        ----------
        ampData : `lsst.afw.image.ImageF`
            Amplifier image to operate on.
        biasLevel : `float`
            Bias level to be added to the image.
        """
        ampArr = ampData.array
        ampArr[:] = ampArr[:] + biasLevel

    def amplifierMultiplyFlat(self, amp, ampData, fracDrop, u0=100.0, v0=100.0):
        """Multiply an amplifier's image data by a flat-like pattern.

        Parameters
        ----------
        amp : `lsst.afw.ampInfo.AmpInfoRecord`
            Amplifier to operate on. Needed for amp<->exp coordinate
            transforms.
        ampData : `lsst.afw.image.ImageF`
            Amplifier image to operate on.
        fracDrop : `float`
            Fractional drop from center to edge of detector along x-axis.
        u0 : `float`
            Peak location in detector coordinates.
        v0 : `float`
            Peak location in detector coordinates.
        """
        if fracDrop >= 1.0:
            raise RuntimeError("Flat fractional drop cannot be greater than 1.0")

        sigma = u0 / np.sqrt(2.0 * fracDrop)

        for x in range(0, ampData.getDimensions().getX()):
            for y in range(0, ampData.getDimensions().getY()):
                (u, v) = self.localCoordToExpCoord(amp, x, y)
                f = np.exp(-0.5 * ((u - u0)**2 + (v - v0)**2) / sigma**2)
                ampData.array[y][x] = (ampData.array[y][x] * f)

    def applyGain(self, ampData, gain):
        """Apply gain to the amplifier's data.
        This method divides the data by the gain
        because the mocks need to convert the data in electron to ADU,
        so it does the inverse operation to applyGains in isrFunctions.

        Parameters
        ----------
        ampData : `lsst.afw.image.ImageF`
            Amplifier image to operate on.
        gain : `float`
            Gain value in e^-/DN.
        """
        ampArr = ampData.array
        ampArr[:] = ampArr[:] / gain

    def roundADU(self, ampData):
        """Round ADU to nearest integer.

        Parameters
        ----------
        ampData : `lsst.afw.image.ImageF`
            Amplifier image to operate on.
        """
        ampArr = ampData.array
        ampArr[:] = np.around(ampArr)

    def amplifierAddXGradient(self, ampData, start, end):
        """Add a x-axis linear gradient to an amplifier's image data.

         This method operates in the amplifier coordinate frame.

        Parameters
        ----------
        ampData : `lsst.afw.image.ImageF`
            Amplifier image to operate on.
        start : `float`
            Start value of the gradient (at x=0).
        end : `float`
            End value of the gradient (at x=xmax).
        """
        nPixX = ampData.getDimensions().getX()
        ampArr = ampData.array
        ampArr[:] = ampArr[:] + (np.interp(range(nPixX), (0, nPixX - 1), (start, end)).reshape(1, nPixX)
                                 + np.zeros(ampData.getDimensions()).transpose())

    def getFullSerialOverscanBBox(self, amp):
        """Get the full serial overscan bounding box from an amplifier.

        This includes the serial/parallel overscan region.

        Parameters
        ----------
        amp : `blah`

        Returns
        -------
        bbox : `lsst.geom.Box2I`
        """
        # This only works for untrimmed data.
        bbox = amp.getRawDataBBox()

        parallelOverscanBBox = amp.getRawParallelOverscanBBox()
        grownImageBBox = bbox.expandedTo(parallelOverscanBBox)

        serialOverscanBBox = amp.getRawSerialOverscanBBox()
        # Extend the serial overscan bbox to include corners
        serialOverscanBBox = geom.Box2I(
            geom.Point2I(serialOverscanBBox.getMinX(),
                         grownImageBBox.getMinY()),
            geom.Extent2I(serialOverscanBBox.getWidth(),
                          grownImageBBox.getHeight()))

        return serialOverscanBBox


class RawMockLSST(IsrMockLSST):
    """Generate a raw exposure suitable for ISR.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.isTrimmed = False
        self.config.doGenerateImage = True
        self.config.doGenerateAmpDict = False

        # Add astro effects
        self.config.doAddSky = True
        self.config.doAddSource = True

        # Add optical effects
        self.config.doAddFringe = True

        # Add instrument effects
        self.config.doAddParallelOverscanRamp = True
        self.config.doAddSerialOverscanRamp = True
        self.config.doAddCrosstalk = True
        self.config.doAddBias = True
        self.config.doAddDark = True

        self.config.doAddFlat = True


class TrimmedRawMockLSST(RawMockLSST):
    """Generate a trimmed raw exposure.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.isTrimmed = True
        self.config.doAddParallelOverscanRamp = False
        self.config.doAddSerialOverscanRamp = False


# Question: what is this used for?
class CalibratedRawMockLSST(RawMockLSST):
    """Generate a trimmed raw exposure.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.isTrimmed = True
        self.config.doGenerateImage = True

        self.config.doAddSky = True
        self.config.doAddSource = True

        self.config.doAddFringe = True

        self.config.doAddParallelOverscanRamp = False
        self.config.doAddSerialOverscanRamp = False
        self.config.doAddCrosstalk = False
        self.config.doAddBias = False
        self.config.doAdd2DBias = False
        self.config.doAddDark = False
        self.config.doApplyGain = False
        self.config.doAddFlat = False
        self.config.doAddClockInjectedOffset = False

        self.config.biasLevel = 0.0
        # Assume combined calibrations are made with 16 inputs.
        self.config.readNoise *= 0.25

        self.config.doRoundADU = False


class ReferenceMockLSST(IsrMockLSST):
    """Parent class for those that make reference calibrations.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.isTrimmed = True
        self.config.doGenerateImage = True

        self.config.calibMode = True

        self.config.doAddSky = False
        self.config.doAddSource = False

        self.config.doAddFringe = False

        self.config.doAddParallelOverscanRamp = False
        self.config.doAddSerialOverscanRamp = False
        self.config.doAddCrosstalk = False
        self.config.doAddBias = False
        self.config.doAdd2DBias = False
        self.config.doAddDark = False
        self.config.doApplyGain = False
        self.config.doAddFlat = False
        self.config.doAddClockInjectedOffset = False

        # Reference calibrations are not integerized.
        self.config.doRoundADU = False


# Classes to generate calibration products mocks.
class DarkMockLSST(ReferenceMockLSST):
    """Simulated reference dark calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doAddDark = True
        self.config.darkTime = 1.0


class BiasMockLSST(ReferenceMockLSST):
    """Simulated combined bias calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # We assume a perfect noiseless bias frame.
        # A combined bias has mean 0
        # so we set its bias level to 0.
        # This is equivalent to doAddBias = False
        # but we do the following instead to be consistent
        # with any other bias products we might want to produce.
        self.config.doAddBias = True
        self.config.biasLevel = 0.0
        self.config.doApplyGain = True
        self.config.doAdd2DBias = True


class FlatMockLSST(ReferenceMockLSST):
    """Simulated reference flat calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doAddFlat = True


class FringeMockLSST(ReferenceMockLSST):
    """Simulated reference fringe calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doAddFringe = True


class BfKernelMockLSST(IsrMockLSST):
    """Simulated brighter-fatter kernel.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doGenerateImage = False
        self.config.doGenerateData = True

        self.config.doBrighterFatter = True
        self.config.doDefects = False
        self.config.doCrosstalkCoeffs = False
        self.config.doTransmissionCurve = False


class DefectMockLSST(IsrMockLSST):
    """Simulated defect list.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doGenerateImage = False
        self.config.doGenerateData = True

        self.config.doBrighterFatter = False
        self.config.doDefects = True
        self.config.doCrosstalkCoeffs = False
        self.config.doTransmissionCurve = False


class CrosstalkCoeffMockLSST(IsrMockLSST):
    """Simulated crosstalk coefficient matrix.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doGenerateImage = False
        self.config.doGenerateData = True

        self.config.doBrighterFatter = False
        self.config.doDefects = False
        self.config.doCrosstalkCoeffs = True
        self.config.doTransmissionCurve = False


class TransmissionMockLSST(IsrMockLSST):
    """Simulated transmission curve.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doGenerateImage = False
        self.config.doGenerateData = True

        self.config.doBrighterFatter = False
        self.config.doDefects = False
        self.config.doCrosstalkCoeffs = False
        self.config.doTransmissionCurve = True
