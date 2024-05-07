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
    # Signal parameters.
    # Most of them are inherited from isrMockConfig but we update
    # some to LSSTcam expected values.
    # TODO: DM-42880 Update values to what is expected in LSSTCam.
    biasLevel = pexConfig.Field(
        dtype=float,
        default=25000.0,
        doc="Background contribution to be generated from the bias offset in ADU.",
    )
    flatMode = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Set to true for producing mock flats.",
    )
    # Inclusion parameters are inherited from isrMock.
    doAddParallelOverscan = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add overscan ramp to parallel overscan and data regions.",
    )
    doAddSerialOverscan = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add overscan ramp to serial overscan and data regions.",
    )
    doApplyGain = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add gain to data.",
    )


class IsrMockLSST(IsrMock):
    """Class to generate consistent mock images for ISR testing.

    ISR testing currently relies on one-off fake images that do not
    accurately mimic the full set of detector effects. This class
    uses the test camera/detector/amplifier structure defined in
    `lsst.afw.cameraGeom.testUtils` to avoid making the test data
    dependent on any of the actual obs package formats.
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

        # We introduce effects as they happen from a source to the signal,
        # so the effects go from electrons to ADU.
        # The ISR steps will then correct these effects in the reverse order.
        for idx, amp in enumerate(exposure.getDetector()):

            # Get image bbox and data
            bbox = None
            if self.config.isTrimmed:
                bbox = amp.getBBox()
            else:
                bbox = amp.getRawDataBBox()

            ampData = exposure.image[bbox]

            # Sky effects in e-
            if self.config.doAddSky:
                self.amplifierAddNoise(ampData, self.config.skyLevel, np.sqrt(self.config.skyLevel))

            if self.config.doAddSource:
                for sourceAmp, sourceFlux, sourceX, sourceY in zip(self.config.sourceAmp,
                                                                   self.config.sourceFlux,
                                                                   self.config.sourceX,
                                                                   self.config.sourceY):
                    if idx == sourceAmp:
                        self.amplifierAddSource(ampData, sourceFlux, sourceX, sourceY)

            # Other effects in e-
            if self.config.doAddFringe:
                self.amplifierAddFringe(amp, ampData, np.array(self.config.fringeScale),
                                        x0=np.array(self.config.fringeX0),
                                        y0=np.array(self.config.fringeY0))

            if self.config.doAddFlat:
                if ampData.getArray().sum() == 0.0:
                    # Add noise
                    self.amplifierAddNoise(ampData, 1.0, 0.0)
                # Multiply each amplifiers by a Gaussian centered on u0 and v0
                u0 = exposure.getDetector().getBBox().getDimensions().getX()/2.
                v0 = exposure.getDetector().getBBox().getDimensions().getY()/2.
                self.amplifierMultiplyFlat(amp, ampData, self.config.flatDrop, u0=u0, v0=v0)

            # ISR effects
            # 1. Add dark in e- (different from isrMock which does it in ADU)
            if self.config.doAddDark:
                self.amplifierAddNoise(ampData,
                                       self.config.darkRate * self.config.darkTime,
                                       np.sqrt(self.config.darkRate * self.config.darkTime))

            # 2. Gain normalize  (from e- to ADU)
            # TODO: DM-43601 gain from PTC per amplifier
            # TODO: DM-36639 gain with temperature dependence
            if self.config.doApplyGain:
                self.applyGain(ampData, self.config.gain)

            # 3. Add read noise to the image region in ADU.
            if not self.config.flatMode:
                self.amplifierAddNoise(ampData, 0.0,
                                    self.config.readNoise / self.config.gain)

        # 4. Apply cross-talk in ADU
        if self.config.doAddCrosstalk:
            ctCalib = CrosstalkCalib()
            exposureClean = exposure.clone()
            for idxS, ampS in enumerate(exposure.getDetector()):
                for idxT, ampT in enumerate(exposure.getDetector()):
                    ampDataTarget = exposure.image[ampT.getBBox() if self.config.isTrimmed
                                                   else ampT.getRawDataBBox()]
                    ampDataSource = ctCalib.extractAmp(exposureClean.image, ampS, ampT,
                                                       isTrimmed=self.config.isTrimmed)
                    self.amplifierAddCT(ampDataSource, ampDataTarget, self.crosstalkCoeffs[idxS][idxT])


        # We now apply parallel and serial overscans
        for amp in exposure.getDetector():
            # Get image bbox and data
            bbox = None
            if self.config.isTrimmed:
                bbox = amp.getBBox()
            else:
                bbox = amp.getRawDataBBox()
            ampData = exposure.image[bbox]

            if self.config.doAddParallelOverscan or self.config.doAddSerialOverscan or self.config.doAddBias:

                allData = ampData

                if self.config.doAddParallelOverscan or self.config.doAddSerialOverscan:
                    # 5. Apply parallel overscan in ADU
                    # First get the parallel and serial overscan bbox
                    # and corresponding data
                    parallelOscanBBox = amp.getRawParallelOverscanBBox()
                    parallelOscanData = exposure.image[parallelOscanBBox]

                    grownImageBBox = bbox.expandedTo(parallelOscanBBox)

                    serialOscanBBox = amp.getRawSerialOverscanBBox()
                    serialOscanBBox = geom.Box2I(
                            geom.Point2I(serialOscanBBox.getMinX(),
                                        grownImageBBox.getMinY()),
                            geom.Extent2I(serialOscanBBox.getWidth(),
                                        grownImageBBox.getHeight()),
                        )
                    serialOscanData = exposure.image[serialOscanBBox]

                    # Add read noise with or without a bias level
                    # to the parallel and serial overscan regions.
                    self.amplifierAddNoise(parallelOscanData, 0.0,
                                            self.config.readNoise / self.config.gain)

                    self.amplifierAddNoise(serialOscanData, 0.0,
                                            self.config.readNoise / self.config.gain)

                    grownImageBBoxAll = grownImageBBox.expandedTo(serialOscanBBox)
                    allData = exposure.image[grownImageBBoxAll]

                    if self.config.doAddParallelOverscan:
                        # Apply gradient along the Y axis
                        self.amplifierAddXGradient(allData, -1.0 * self.config.overscanScale,
                                            1.0 * self.config.overscanScale)

        # 6. Add Parallel overscan xtalk.
        # TODO: DM-43286

                # Add bias level to the whole image
                # (science and overscan regions if any)
                self.addBiasLevel(allData, self.config.biasLevel if self.config.doAddBias else 0.0)

                if self.config.doAddSerialOverscan:
                    # Apply gradient along the Y axis
                    self.amplifierAddYGradient(allData, -1.0 * self.config.overscanScale,
                                               1.0 * self.config.overscanScale)

        if self.config.doGenerateAmpDict:
            expDict = dict()
            for amp in exposure.getDetector():
                expDict[amp.getName()] = exposure
            return expDict
        else:
            return exposure

    # Simple data values.
    def addBiasLevel(self, ampData, biasLevel):
        """Add Gaussian noise to an amplifier's image data.

         This method operates in the amplifier coordinate frame.

        Parameters
        ----------
        ampData : `lsst.afw.image.ImageF`
            Amplifier image to operate on.
        biasLevel : `float`
            Bias level to be added to the image
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

        Notes
        -----
        This uses a 2-d Gaussian to simulate an illumination pattern
        that falls off towards the edge of the detector. The (x, y)
        coordinates are in the frame of the amplifier, and (u, v) in
        the frame of the full trimmed image.
        """
        if fracDrop >= 1.0:
            raise RuntimeError("Flat fractional drop cannot be greater than 1.0")

        sigma = u0 / np.sqrt(2.0 *fracDrop)

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

    def amplifierAddXGradient(self, ampData, start, end):
        """Add a x-axis linear gradient to an amplifier's image data.

         This method operates in the amplifier coordinate frame.

        Parameters
        ----------
        ampData : `lsst.afw.image.ImageF`
            Amplifier image to operate on.
        start : `float`
            Start value of the gradient (at y=0).
        end : `float`
            End value of the gradient (at y=ymax).
        """
        nPixX = ampData.getDimensions().getX()
        ampArr = ampData.array
        ampArr[:] = ampArr[:] + (np.interp(range(nPixX), (0, nPixX - 1), (start, end)).reshape(1, nPixX)
                                 + np.zeros(ampData.getDimensions()).transpose())


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

        # Add instru effects
        self.config.doAddParallelOverscan = True
        self.config.doAddSerialOverscan = True
        self.config.doAddCrosstalk = False
        self.config.doAddBias = True
        self.config.doAddDark = True

        self.config.doAddFlat = True


class TrimmedRawMockLSST(RawMockLSST):
    """Generate a trimmed raw exposure.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.isTrimmed = True
        self.config.doAddParallelOverscan = False
        self.config.doAddSerialOverscan = False


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

        self.config.doAddParallelOverscan = False
        self.config.doAddSerialOverscan = False
        self.config.doAddCrosstalk = False
        self.config.doAddBias = False
        self.config.doAddDark = False
        self.config.doApplyGain = False
        self.config.doAddFlat = False

        self.config.biasLevel = 0.0
        # Assume combined calibrations are made with 16 inputs.
        self.config.readNoise *= 0.25


class ReferenceMockLSST(IsrMockLSST):
    """Parent class for those that make reference calibrations.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.isTrimmed = True
        self.config.doGenerateImage = True

        self.config.doAddSky = False
        self.config.doAddSource = False

        self.config.doAddFringe = False

        self.config.doAddParallelOverscan = False
        self.config.doAddSerialOverscan = False
        self.config.doAddCrosstalk = False
        self.config.doAddBias = False
        self.config.doAddDark = False
        self.config.doApplyGain = False
        self.config.doAddFlat = False


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
        # A combined bias has mean 0
        # so we set its bias level to 0.
        # This is equivalent to doAddBias = False
        # but we do the following instead to be consistent
        # with any other bias products we might want to produce.
        self.config.doAddBias = True
        self.config.biasLevel = 0.0
        self.config.doApplyGain = True
        # Assume combined calibrations are made with 16 inputs.
        self.config.readNoise = 10.0*0.25


class FlatMockLSST(ReferenceMockLSST):
    """Simulated reference flat calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doAddFlat = True
        self.config.flatMode = True


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

        # calibration products configs
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
