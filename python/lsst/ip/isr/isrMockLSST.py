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

__all__ = ["IsrMockLSSTConfig", "IsrMockLSST"]

import copy
import numpy as np
import tempfile

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
    doAddGain = pexConfig.Field(
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
        This method currently constructs a "raw" data image by:

        * Generating a simulated sky with noise
        * Adding a single Gaussian "star"
        * Adding the fringe signal
        * Multiplying the frame by the simulated flat
        * Adding dark current (and noise)
        * Adding a bias offset (and noise)
        * Adding an overscan gradient parallel to the pixel y-axis
        * Simulating crosstalk by adding a scaled version of each
          amplifier to each other amplifier.

        The exposure with image data constructed this way is in one of
        three formats.

        * A single image, with overscan and prescan regions retained
        * A single image, with overscan and prescan regions trimmed
        * A `dict`, containing the amplifer data indexed by the
          amplifier name.

        The nonlinearity, CTE, and brighter fatter are currently not
        implemented.

        Note that this method generates an image in the reverse
        direction as the ISR processing, as the output image here has
        had a series of instrument effects added to an idealized
        exposure.
        """
        exposure = self.getExposure()

        for idx, amp in enumerate(exposure.getDetector()):

            # Get image bbox and data
            imageBBox = amp.getRawDataBBox()
            ampData = exposure.image[imageBBox]

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

            # Post ISR effects in e-
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
                                       np.sqrt(self.config.darkRate
                                               * self.config.darkTime))

            # 2. Gain normalize  (from e- to ADU)
            # TODO: DM-??? gain from PTC per amplifier
            # TODO: DM-??? gain with temperature dependence
            if self.config.doAddGain:
                self.addGain(ampData, self.config.gain)

        # 3. TODO: Add bias frame (make fake bias frame - could be 0)
        # 4. Apply cross-talk in ADU
        if self.config.doAddCrosstalk:
            ctCalib = CrosstalkCalib()
            for idxS, ampS in enumerate(exposure.getDetector()):
                for idxT, ampT in enumerate(exposure.getDetector()):
                    ampDataT = exposure.image[ampT.getRawDataBBox()]
                    outAmp = ctCalib.extractAmp(exposure.getImage(), ampS, ampT,
                                                isTrimmed=False)
                    self.amplifierAddCT(outAmp, ampDataT, self.crosstalkCoeffs[idxS][idxT])

        for amp in exposure.getDetector():
            parallelOscanBBox = amp.getRawParallelOverscanBBox()
            parallelOscanData = exposure.image[parallelOscanBBox]

            serialOscanBBox = amp.getRawSerialOverscanBBox()

            # 5. Apply parallel overscan in ADU
            if self.config.doAddParallelOverscan:
                # Add noise to the parallel overscan region.
                self.amplifierAddNoise(parallelOscanData, 0.0,
                                       self.config.readNoise / self.config.gain)

                # Apply gradient along Y axis to the image
                # and parallel overscan regions.
                self.amplifierAddYGradient(ampData, -1.0 * self.config.overscanScale,
                                           1.0 * self.config.overscanScale)
                self.amplifierAddYGradient(parallelOscanData, -1.0 * self.config.overscanScale,
                                           1.0 * self.config.overscanScale)

        # 6. Add Parallel overscan xtalk.
        # TODO: DM-???

            # 7. Add bias level to each amplifier in ADU
            # The bias level is in ADU and readNoise in electrons.
            # We always want to do this when adding overscans.
                self.amplifierAddNoise(ampData, self.config.biasLevel,
                                       self.config.readNoise / self.config.gain)

            # 8. Apply serial overscan in ADU
            if self.config.doAddSerialOverscan:
                # 1. We grow the image to the parallel oversan region
                # (we do this in case there are prescan regions)
                grownImageBBox = imageBBox.expandedTo(parallelOscanBBox)
                # 2. Now we grow the serial overscan region
                # to include the corners
                serialOscanBBox = geom.Box2I(
                    geom.Point2I(serialOscanBBox.getMinX(),
                                 grownImageBBox.getMinY()),
                    geom.Extent2I(serialOscanBBox.getWidth(),
                                  grownImageBBox.getHeight()),
                )
                serialOscanData = exposure.image[serialOscanBBox]

                # Add noise to the serial overscan region.
                self.amplifierAddNoise(serialOscanData, 0.0,
                                       self.config.readNoise / self.config.gain)

                # Apply gradient along X axis to the whole amp
                # First the serial overscan and corners region
                self.amplifierAddXGradient(serialOscanData, -1.0 * self.config.overscanScale,
                                           1.0 * self.config.overscanScale)
                # And second the image and parallel overscan region
                self.amplifierAddXGradient(ampData, -1.0 * self.config.overscanScale,
                                           1.0 * self.config.overscanScale)
                self.amplifierAddXGradient(parallelOscanData, -1.0 * self.config.overscanScale,
                                           1.0 * self.config.overscanScale)

        if self.config.doGenerateAmpDict:
            expDict = dict()
            for amp in exposure.getDetector():
                expDict[amp.getName()] = exposure
            return expDict
        else:
            return exposure

    def addGain(self, ampData, gain):
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
        self.config.doAddCrosstalk = True
        self.config.doAddBias = True
        self.config.doAddDark = True

        self.config.doAddFlat = True


class CalibratedRawMockLSST(RawMockLSST):
    """Generate a trimmed raw exposure.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
#        self.config.isTrimmed = True
        self.config.doGenerateImage = True

        self.config.doAddSky = True
        self.config.doAddSource = True

        self.config.doAddFringe = True

        self.config.doAddParallelOverscan = False
        self.config.doAddSerialOverscan = False
        self.config.doAddCrosstalk = False
        self.config.doAddBias = False
        self.config.doAddDark = False
        self.config.doAddGain = False
        self.config.doAddFlat = False

        self.config.biasLevel = 0.0
        self.config.readNoise = 10.0


class ReferenceMockLSST(IsrMockLSST):
    """Parent class for those that make reference calibrations.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doGenerateImage = True

        self.config.doAddSky = False
        self.config.doAddSource = False

        self.config.doAddFringe = False

        self.config.doAddParallelOverscan = False
        self.config.doAddSerialOverscan = False
        self.config.doAddCrosstalk = False
        self.config.doAddBias = False
        self.config.doAddDark = False
        self.config.doAddGain = False
        self.config.doAddFlat = False


# Classes to generate calibration products mocks.
class DarkMockLSST(ReferenceMockLSST):
    """Simulated master dark calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doAddDark = True
        self.config.darkTime = 1.0


class BiasMockLSST(ReferenceMockLSST):
    """Simulated master bias calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doAddBias = True
        self.config.doAddGain = True
        self.config.readNoise = 10.0


class FlatMockLSST(ReferenceMockLSST):
    """Simulated master flat calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doAddFlat = True


class FringeMockLSST(ReferenceMockLSST):
    """Simulated master fringe calibration.
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


class MockLSSTDataContainer(object):
    """Container for holding ISR LSST mock objects.
    """
    dataId = "isrMockLSST Fake Data"
    darkval = 2.  # e-/sec
    oscan = 250.  # DN
    gradient = .10
    exptime = 15.0  # seconds
    darkexptime = 15.0  # seconds

    def __init__(self, **kwargs):
        if 'config' in kwargs.keys():
            self.config = kwargs['config']
        else:
            self.config = None

    def expectImage(self):
        if self.config is None:
            self.config = IsrMockLSSTConfig()
        self.config.doGenerateImage = True
        self.config.doGenerateData = False

    def expectData(self):
        if self.config is None:
            self.config = IsrMockLSSTConfig()
        self.config.doGenerateImage = False
        self.config.doGenerateData = True

    def get(self, dataType, **kwargs):
        """Return an appropriate data product.

        Parameters
        ----------
        dataType : `str`
            Type of data product to return.

        Returns
        -------
        mock : IsrMockLSST.run() result
            The output product.
        """
        if "_filename" in dataType:
            self.expectData()
            return tempfile.mktemp(), "mock"
        elif 'transmission_' in dataType:
            self.expectData()
            return TransmissionMockLSST(config=self.config).run()
        elif dataType == 'ccdExposureId':
            self.expectData()
            return 20090913
        elif dataType == 'camera':
            self.expectData()
            return IsrMockLSST(config=self.config).getCamera()
        elif dataType == 'raw':
            self.expectImage()
            return RawMockLSST(config=self.config).run()
        elif dataType == 'bias':
            self.expectImage()
            return BiasMockLSST(config=self.config).run()
        elif dataType == 'dark':
            self.expectImage()
            return DarkMockLSST(config=self.config).run()
        elif dataType == 'flat':
            self.expectImage()
            return FlatMockLSST(config=self.config).run()
        elif dataType == 'fringe':
            self.expectImage()
            return FringeMockLSST(config=self.config).run()
        elif dataType == 'defects':
            self.expectData()
            return DefectMockLSST(config=self.config).run()
        elif dataType == 'bfKernel':
            self.expectData()
            return BfKernelMockLSST(config=self.config).run()
        elif dataType == 'linearizer':
            return None
        elif dataType == 'crosstalkSources':
            return None
        else:
            raise RuntimeError("ISR DataRefMock cannot return %s.", dataType)


class MockFringeLSSTContainer(object):
    """Container for mock fringe data.
    """
    dataId = "isrMockLSST Fake Data"
    darkval = 2.  # e-/sec
    oscan = 250.  # DN
    gradient = .10
    exptime = 15  # seconds
    darkexptime = 40.  # seconds

    def __init__(self, **kwargs):
        if 'config' in kwargs.keys():
            self.config = kwargs['config']
        else:
            self.config = IsrMockLSSTConfig()
            self.config.doAddFringe = True
            self.config.readNoise = 10.0

    def get(self, dataType, **kwargs):
        """Return an appropriate data product.

        Parameters
        ----------
        dataType : `str`
            Type of data product to return.

        Returns
        -------
        mock : IsrMockLSST.run() result
            The output product.
        """
        if "_filename" in dataType:
            return tempfile.mktemp(), "mock"
        elif 'transmission_' in dataType:
            return TransmissionMockLSST(config=self.config).run()
        elif dataType == 'ccdExposureId':
            return 20090913
        elif dataType == 'camera':
            return IsrMockLSST(config=self.config).getCamera()
        elif dataType == 'raw':
            return CalibratedRawMockLSST(config=self.config).run()
        elif dataType == 'bias':
            return BiasMockLSST(config=self.config).run()
        elif dataType == 'dark':
            return DarkMockLSST(config=self.config).run()
        elif dataType == 'flat':
            return FlatMockLSST(config=self.config).run()
        elif dataType == 'fringe':
            fringes = []
            configCopy = copy.deepcopy(self.config)
            for scale, x, y in zip(self.config.fringeScale, self.config.fringeX0, self.config.fringeY0):
                configCopy.fringeScale = [1.0]
                configCopy.fringeX0 = [x]
                configCopy.fringeY0 = [y]
                fringes.append(FringeMockLSST(config=configCopy).run())
            return fringes
        elif dataType == 'defects':
            return DefectMockLSST(config=self.config).run()
        elif dataType == 'bfKernel':
            return BfKernelMockLSST(config=self.config).run()
        elif dataType == 'linearizer':
            return None
        elif dataType == 'crosstalkSources':
            return None
        else:
            return None
