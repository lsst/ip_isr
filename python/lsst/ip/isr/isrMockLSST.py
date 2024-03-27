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
        default=30000.0,
        doc="Background contribution to be generated from the bias offset in ADU.",
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
                bbox = amp.getRawDataBBox().shiftedBy(amp.getRawXYOffset())

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
                    self.amplifierAddNoise(ampData, 1.0, 0.0)
                u0 = exposure.getDimensions().getX()
                v0 = exposure.getDimensions().getY()
                self.amplifierMultiplyFlat(amp, ampData, self.config.flatDrop, u0=u0, v0=v0)

            # ISR effects
            # 1. Add dark in e- (different from isrMock which does it in ADU)
            if self.config.doAddDark:
                self.amplifierAddNoise(ampData,
                                       self.config.darkRate * self.config.darkTime,
                                       np.sqrt(self.config.darkRate * self.config.darkTime))

            # 2. Gain normalize  (from e- to ADU)
            # TODO: DM-??? gain from PTC per amplifier
            # TODO: DM-??? gain with temperature dependence
            if self.config.doApplyGain:
                self.applyGain(ampData, self.config.gain)

            # 3. We either add the bias frame or only the read noise
            # to the image region in ADU.
            self.amplifierAddNoise(ampData, self.config.biasLevel if self.config.doAddBias else 0.0,
                                   self.config.readNoise / self.config.gain)

        # 4. Apply cross-talk in ADU
        if self.config.doAddCrosstalk:
            ctCalib = CrosstalkCalib()
            for idxS, ampS in enumerate(exposure.getDetector()):
                for idxT, ampT in enumerate(exposure.getDetector()):
                    ampDataT = exposure.image[ampT.getBBox() if self.config.isTrimmed
                                              else ampT.getRawDataBBox().shiftedBy(ampT.getRawXYOffset())]
                    outAmp = ctCalib.extractAmp(exposure.getImage(), ampS, ampT,
                                                isTrimmed=self.config.isTrimmed)
                    self.amplifierAddCT(outAmp, ampDataT, self.crosstalkCoeffs[idxS][idxT])

        # We now apply parallel and serial overscans
        for amp in exposure.getDetector():
            # Get image bbox and data
            bbox = None
            if self.config.isTrimmed:
                bbox = amp.getBBox()
            else:
                bbox = amp.getRawDataBBox().shiftedBy(amp.getRawXYOffset())
            ampData = exposure.image[bbox]

            # Get overscan bbox and data
            if not self.config.isTrimmed:
                parallelOscanBBox = amp.getRawParallelOverscanBBox().shiftedBy(amp.getRawXYOffset())
                parallelOscanData = exposure.image[parallelOscanBBox]

                serialOscanBBox = amp.getRawSerialOverscanBBox().shiftedBy(amp.getRawXYOffset())

            # 5. Apply parallel overscan in ADU
            if self.config.doAddParallelOverscan:
                if not self.config.isTrimmed:
                    # Add bias frame or read noise
                    # to the parallel overscan region.
                    self.amplifierAddNoise(parallelOscanData, self.config.biasLevel
                                           if self.config.doAddBias else 0.0,
                                           self.config.readNoise / self.config.gain)
                    # Apply gradient along the Y axis
                    # to the parallel overscan region.
                    self.amplifierAddYGradient(parallelOscanData, -1.0 * self.config.overscanScale,
                                               1.0 * self.config.overscanScale)

                # Apply gradient along the Y axis to the image region
                self.amplifierAddYGradient(ampData, -1.0 * self.config.overscanScale,
                                           1.0 * self.config.overscanScale)

        # 6. Add Parallel overscan xtalk.
        # TODO: DM-43286

            if self.config.doAddSerialOverscan:
                if not self.config.isTrimmed:
                    # We grow the image to the parallel overscan region
                    # (we do this instead of using the whole raw region
                    # in case there are prescan regions)
                    grownImageBBox = bbox.expandedTo(parallelOscanBBox)
                    # Now we grow the serial overscan region
                    # to include the corners
                    serialOscanBBox = geom.Box2I(
                        geom.Point2I(serialOscanBBox.getMinX(),
                                     grownImageBBox.getMinY()),
                        geom.Extent2I(serialOscanBBox.getWidth(),
                                      grownImageBBox.getHeight()),
                    )
                    serialOscanData = exposure.image[serialOscanBBox]

                    # Add bias frame or read noise
                    # to the serial overscan region.
                    self.amplifierAddNoise(serialOscanData, self.config.biasLevel
                                           if self.config.doAddBias else 0.0,
                                           self.config.readNoise / self.config.gain)

                    # 7. Apply serial overscan in ADU
                    # Apply gradient along the X axis to both overscan regions.
                    self.amplifierAddXGradient(serialOscanData, -1.0 * self.config.overscanScale,
                                               1.0 * self.config.overscanScale)
                    self.amplifierAddXGradient(parallelOscanData, -1.0 * self.config.overscanScale,
                                               1.0 * self.config.overscanScale)

                # Apply gradient along the X axis to the image region.
                self.amplifierAddXGradient(ampData, -1.0 * self.config.overscanScale,
                                           1.0 * self.config.overscanScale)

        if self.config.doGenerateAmpDict:
            expDict = dict()
            for amp in exposure.getDetector():
                expDict[amp.getName()] = exposure
            return expDict
        else:
            return exposure

    def applyGain(self, ampData, gain):
        """Apply gain to the amplifier's data.

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
        self.config.doAddBias = True
        self.config.doApplyGain = True
        # Assume combined calibrations are made with 16 inputs.
        self.config.readNoise = 10.0*0.25


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
