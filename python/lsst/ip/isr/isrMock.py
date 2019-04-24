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

import copy
import numpy as np
import tempfile

import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.cameraGeom.utils as afwUtils
import lsst.afw.cameraGeom.testUtils as afwTestUtils
from lsst.meas.algorithms import Defects
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .crosstalk import X_FLIP, Y_FLIP

__all__ = ["IsrMockConfig", "IsrMock", "RawMock", "TrimmedRawMock", "RawDictMock",
           "CalibratedRawMock", "MasterMock",
           "BiasMock", "DarkMock", "FlatMock", "FringeMock", "UntrimmedFringeMock",
           "BfKernelMock", "DefectMock", "CrosstalkCoeffMock", "TransmissionMock",
           "DataRefMock"]


class IsrMockConfig(pexConfig.Config):
    """Configuration parameters for isrMock.

    These parameters produce generic fixed position signals from
    various sources, and combine them in a way that matches how those
    signals are combined to create real data. The camera used is the
    test camera defined by the afwUtils code.
    """
    # Detector parameters. "Exposure" parameters.
    isLsstLike = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="If True, products have one raw image per amplifier, otherwise, one raw image per detector.",
    )
    isTrimmed = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="If True, amplifiers have been trimmed and mosaicked to remove regions outside the data BBox.",
    )
    detectorIndex = pexConfig.Field(
        dtype=int,
        default=20,
        doc="Index for the detector to use. The default value uses a standard 2x4 array of amps.",
    )
    rngSeed = pexConfig.Field(
        dtype=int,
        default=20000913,
        doc="Seed for random number generator used to add noise.",
    )
    # TODO: DM-18345 Check that mocks scale correctly when gain != 1.0
    gain = pexConfig.Field(
        dtype=float,
        default=1.0,
        doc="Gain for simulated data in e^-/DN.",
    )
    readNoise = pexConfig.Field(
        dtype=float,
        default=5.0,
        doc="Read noise of the detector in e-.",
    )
    expTime = pexConfig.Field(
        dtype=float,
        default=5.0,
        doc="Exposure time for simulated data.",
    )

    # Signal parameters
    skyLevel = pexConfig.Field(
        dtype=float,
        default=1000.0,
        doc="Background contribution to be generated from 'the sky' in DN.",
    )
    sourceFlux = pexConfig.ListField(
        dtype=float,
        default=[45000.0],
        doc="Peak flux level (in DN) of simulated 'astronomical sources'.",
    )
    sourceAmp = pexConfig.ListField(
        dtype=int,
        default=[0],
        doc="Amplifier to place simulated 'astronomical sources'.",
    )
    sourceX = pexConfig.ListField(
        dtype=float,
        default=[50.0],
        doc="Peak position (in amplifier coordinates) of simulated 'astronomical sources'.",
    )
    sourceY = pexConfig.ListField(
        dtype=float,
        default=[25.0],
        doc="Peak position (in amplifier coordinates) of simulated 'astronomical sources'.",
    )
    overscanScale = pexConfig.Field(
        dtype=float,
        default=100.0,
        doc="Amplitude (in DN) of the ramp function to add to overscan data.",
    )
    biasLevel = pexConfig.Field(
        dtype=float,
        default=8000.0,
        doc="Background contribution to be generated from the bias offset in DN.",
    )
    darkRate = pexConfig.Field(
        dtype=float,
        default=5.0,
        doc="Background level contribution (in e-/s) to be generated from dark current.",
    )
    darkTime = pexConfig.Field(
        dtype=float,
        default=5.0,
        doc="Exposure time for the dark current contribution.",
    )
    flatDrop = pexConfig.Field(
        dtype=float,
        default=0.1,
        doc="Fractional flux drop due to flat from center to edge of detector along x-axis.",
    )
    fringeScale = pexConfig.ListField(
        dtype=float,
        default=[200.0],
        doc="Peak fluxes for the components of the fringe ripple in DN.",
    )
    fringeX0 = pexConfig.ListField(
        dtype=float,
        default=[-100],
        doc="Center position for the fringe ripples.",
    )
    fringeY0 = pexConfig.ListField(
        dtype=float,
        default=[-0],
        doc="Center position for the fringe ripples.",
    )

    # Inclusion parameters
    doAddSky = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Apply 'sky' signal to output image.",
    )
    doAddSource = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add simulated source to output image.",
    )
    doAddCrosstalk = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Apply simulated crosstalk to output image. This cannot be corrected by ISR, "
        "as detector.hasCrosstalk()==False.",
    )
    doAddOverscan = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="If untrimmed, add overscan ramp to overscan and data regions.",
    )
    doAddBias = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add bias signal to data.",
    )
    doAddDark = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add dark signal to data.",
    )
    doAddFlat = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add flat signal to data.",
    )
    doAddFringe = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add fringe signal to data.",
    )

    # Datasets to create and return instead of generating an image.
    doTransmissionCurve = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Return a simulated transmission curve.",
    )
    doDefects = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Return a simulated defect list.",
    )
    doBrighterFatter = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Return a simulated brighter-fatter kernel.",
    )
    doCrosstalkCoeffs = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Return the matrix of crosstalk coefficients.",
    )
    doDataRef = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Return a simulated gen2 butler dataRef.",
    )
    doGenerateImage = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Return the generated output image if True.",
    )
    doGenerateData = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Return a non-image data structure if True.",
    )
    doGenerateAmpDict = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Return a dict of exposure amplifiers instead of an afwImage.Exposure.",
    )


class IsrMock(pipeBase.Task):
    """Class to generate consistent mock images for ISR testing.

    ISR testing currently relies on one-off fake images that do not
    accurately mimic the full set of detector effects. This class
    uses the test camera/detector/amplifier structure defined in
    `lsst.afw.cameraGeom.testUtils` to avoid making the test data
    dependent on any of the actual obs package formats.
    """
    ConfigClass = IsrMockConfig
    _DefaultName = "isrMock"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.rng = np.random.RandomState(self.config.rngSeed)
        self.crosstalkCoeffs = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, -1e-3, 0.0, 0.0],
                                         [1e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                         [1e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                         [1e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                         [1e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                         [1e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                         [1e-2, 0.0, 0.0, 2.2e-2, 0.0, 0.0, 0.0, 0.0],
                                         [1e-2, 5e-3, 5e-4, 3e-3, 4e-2, 5e-3, 5e-3, 0.0]])

        self.bfKernel = np.array([[1., 4., 7., 4., 1.],
                                  [4., 16., 26., 16., 4.],
                                  [7., 26., 41., 26., 7.],
                                  [4., 16., 26., 16., 4.],
                                  [1., 4., 7., 4., 1.]]) / 273.0

    def run(self):
        """Generate a mock ISR product, and return it.

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

    def makeData(self):
        """Generate simulated ISR data.

        Currently, only the class defined crosstalk coefficient
        matrix, brighter-fatter kernel, a constant unity transmission
        curve, or a simple single-entry defect list can be generated.

        Returns
        -------
        dataProduct :
            Simulated ISR data product.
        """
        if sum(map(bool, [self.config.doBrighterFatter,
                          self.config.doDefects,
                          self.config.doTransmissionCurve,
                          self.config.doCrosstalkCoeffs])) != 1:
            raise RuntimeError("Only one data product can be generated at a time.")
        elif self.config.doBrighterFatter is True:
            return self.makeBfKernel()
        elif self.config.doDefects is True:
            return self.makeDefectList()
        elif self.config.doTransmissionCurve is True:
            return self.makeTransmissionCurve()
        elif self.config.doCrosstalkCoeffs is True:
            return self.crosstalkCoeffs
        else:
            return None

    def makeBfKernel(self):
        """Generate a simple Gaussian brighter-fatter kernel.

        Returns
        -------
        kernel : `numpy.ndarray`
            Simulated brighter-fatter kernel.
        """
        return self.bfKernel

    def makeDefectList(self):
        """Generate a simple single-entry defect list.

        Returns
        -------
        defectList : `lsst.meas.algorithms.Defects`
            Simulated defect list
        """
        return Defects([lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                        lsst.geom.Extent2I(40, 50))])

    def makeCrosstalkCoeff(self):
        """Generate the simulated crosstalk coefficients.

        Returns
        -------
        coeffs : `numpy.ndarray`
            Simulated crosstalk coefficients.
        """

        return self.crosstalkCoeffs

    def makeTransmissionCurve(self):
        """Generate a simulated flat transmission curve.

        Returns
        -------
        transmission : `lsst.afw.image.TransmissionCurve`
            Simulated transmission curve.
        """

        return afwImage.TransmissionCurve.makeIdentity()

    def makeImage(self):
        """Generate a simulated ISR image.

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
            bbox = None
            if self.config.isTrimmed is True:
                bbox = amp.getBBox()
            else:
                bbox = amp.getRawDataBBox()

            ampData = exposure.image[bbox]

            if self.config.doAddSky is True:
                self.amplifierAddNoise(ampData, self.config.skyLevel, np.sqrt(self.config.skyLevel))

            if self.config.doAddSource is True:
                for sourceAmp, sourceFlux, sourceX, sourceY in zip(self.config.sourceAmp,
                                                                   self.config.sourceFlux,
                                                                   self.config.sourceX,
                                                                   self.config.sourceY):
                    if idx == sourceAmp:
                        self.amplifierAddSource(ampData, sourceFlux, sourceX, sourceY)

            if self.config.doAddFringe is True:
                self.amplifierAddFringe(amp, ampData, np.array(self.config.fringeScale),
                                        x0=np.array(self.config.fringeX0),
                                        y0=np.array(self.config.fringeY0))

            if self.config.doAddFlat is True:
                if ampData.getArray().sum() == 0.0:
                    self.amplifierAddNoise(ampData, 1.0, 0.0)
                u0 = exposure.getDimensions().getX()
                v0 = exposure.getDimensions().getY()
                self.amplifierMultiplyFlat(amp, ampData, self.config.flatDrop, u0=u0, v0=v0)

            if self.config.doAddDark is True:
                self.amplifierAddNoise(ampData,
                                       self.config.darkRate * self.config.darkTime / self.config.gain,
                                       np.sqrt(self.config.darkRate *
                                               self.config.darkTime / self.config.gain))

        if self.config.doAddCrosstalk is True:

            for idxS, ampS in enumerate(exposure.getDetector()):
                for idxT, ampT in enumerate(exposure.getDetector()):
                    ampDataS = exposure.image[ampS.getBBox() if self.config.isTrimmed
                                              else ampS.getRawDataBBox()]
                    ampDataT = exposure.image[ampT.getBBox() if self.config.isTrimmed
                                              else ampT.getRawDataBBox()]
                    ampDataS = afwMath.flipImage(ampDataS,
                                                 (X_FLIP[ampS.getReadoutCorner()] ^
                                                  X_FLIP[ampT.getReadoutCorner()]),
                                                 (Y_FLIP[ampS.getReadoutCorner()] ^
                                                  Y_FLIP[ampT.getReadoutCorner()]))
                    self.amplifierAddCT(ampDataS, ampDataT, self.crosstalkCoeffs[idxT][idxS])

        for amp in exposure.getDetector():
            bbox = None
            if self.config.isTrimmed is True:
                bbox = amp.getBBox()
            else:
                bbox = amp.getRawDataBBox()

            ampData = exposure.image[bbox]

            if self.config.doAddBias is True:
                self.amplifierAddNoise(ampData, self.config.biasLevel,
                                       self.config.readNoise / self.config.gain)

            if self.config.doAddOverscan is True:
                oscanBBox = amp.getRawHorizontalOverscanBBox()
                oscanData = exposure.image[oscanBBox]
                self.amplifierAddNoise(oscanData, self.config.biasLevel,
                                       self.config.readNoise / self.config.gain)

                self.amplifierAddYGradient(ampData, -1.0 * self.config.overscanScale,
                                           1.0 * self.config.overscanScale)
                self.amplifierAddYGradient(oscanData, -1.0 * self.config.overscanScale,
                                           1.0 * self.config.overscanScale)

        if self.config.doGenerateAmpDict is True:
            expDict = dict()
            for amp in exposure.getDetector():
                expDict[amp.getName()] = exposure
            return expDict
        else:
            return exposure

    # afw primatives to construct the image structure
    def getCamera(self):
        """Construct a test camera object.

        Returns
        -------
        camera : `lsst.afw.cameraGeom.camera`
            Test camera.
        """
        cameraWrapper = afwTestUtils.CameraWrapper(self.config.isLsstLike)
        camera = cameraWrapper.camera
        return camera

    def getExposure(self):
        """Construct a test exposure.

        The test exposure has a simple WCS set, as well as a list of
        unlikely header keywords that can be removed during ISR
        processing to exercise that code.

        Returns
        -------
        exposure : `lsst.afw.exposure.Exposure`
            Construct exposure containing masked image of the
            appropriate size.
        """
        camera = self.getCamera()
        detector = camera[self.config.detectorIndex]
        image = afwUtils.makeImageFromCcd(detector,
                                          isTrimmed=self.config.isTrimmed,
                                          showAmpGain=False,
                                          rcMarkSize=0,
                                          binSize=1,
                                          imageFactory=afwImage.ImageF)

        var = afwImage.ImageF(image.getDimensions())
        mask = afwImage.Mask(image.getDimensions())
        image.assign(0.0)

        maskedImage = afwImage.makeMaskedImage(image, mask, var)
        exposure = afwImage.makeExposure(maskedImage)
        exposure.setDetector(detector)
        exposure.setWcs(self.getWcs())

        visitInfo = afwImage.VisitInfo(exposureTime=self.config.expTime, darkTime=self.config.darkTime)
        exposure.getInfo().setVisitInfo(visitInfo)

        metadata = exposure.getMetadata()
        metadata.add("SHEEP", 7.3, "number of sheep on farm")
        metadata.add("MONKEYS", 155, "monkeys per tree")
        metadata.add("VAMPIRES", 4, "How scary are vampires.")

        for amp in exposure.getDetector():
            amp.setLinearityCoeffs((0., 1., 0., 0.))
            amp.setLinearityType("Polynomial")
            amp.setGain(self.config.gain)

        exposure.image.array[:] = np.zeros(exposure.getImage().getDimensions()).transpose()
        exposure.mask.array[:] = np.zeros(exposure.getMask().getDimensions()).transpose()
        exposure.variance.array[:] = np.zeros(exposure.getVariance().getDimensions()).transpose()

        return exposure

    def getWcs(self):
        """Construct a dummy WCS object.

        Taken from the deprecated ip_isr/examples/exampleUtils.py.

        This is not guaranteed, given the distortion and pixel scale
        listed in the afwTestUtils camera definition.

        Returns
        -------
        wcs : `lsst.afw.geom.SkyWcs`
            Test WCS transform.
        """
        return afwGeom.makeSkyWcs(crpix=lsst.geom.Point2D(0.0, 100.0),
                                  crval=lsst.geom.SpherePoint(45.0, 25.0, lsst.geom.degrees),
                                  cdMatrix=afwGeom.makeCdMatrix(scale=1.0*lsst.geom.degrees))

    def localCoordToExpCoord(self, ampData, x, y):
        """Convert between a local amplifier coordinate and the full
        exposure coordinate.

        Parameters
        ----------
        ampData : `lsst.afw.image.ImageF`
            Amplifier image to use for conversions.
        x : `int`
            X-coordinate of the point to transform.
        y : `int`
            Y-coordinate of the point to transform.

        Returns
        -------
        u : `int`
            Transformed x-coordinate.
        v : `int`
            Transformed y-coordinate.

        Notes
        -----
        The output is transposed intentionally here, to match the
        internal transpose between numpy and afw.image coordinates.
        """
        u = x + ampData.getBBox().getBeginX()
        v = y + ampData.getBBox().getBeginY()

        return (v, u)

    # Simple data values.
    def amplifierAddNoise(self, ampData, mean, sigma):
        """Add Gaussian noise to an amplifier's image data.

         This method operates in the amplifier coordinate frame.

        Parameters
        ----------
        ampData : `lsst.afw.image.ImageF`
            Amplifier image to operate on.
        mean : `float`
            Mean value of the Gaussian noise.
        sigma : `float`
            Sigma of the Gaussian noise.
        """
        ampArr = ampData.array
        ampArr[:] = ampArr[:] + self.rng.normal(mean, sigma,
                                                size=ampData.getDimensions()).transpose()

    def amplifierAddYGradient(self, ampData, start, end):
        """Add a y-axis linear gradient to an amplifier's image data.

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
        nPixY = ampData.getDimensions().getY()
        ampArr = ampData.array
        ampArr[:] = ampArr[:] + (np.interp(range(nPixY), (0, nPixY - 1), (start, end)).reshape(nPixY, 1) +
                                 np.zeros(ampData.getDimensions()).transpose())

    def amplifierAddSource(self, ampData, scale, x0, y0):
        """Add a single Gaussian source to an amplifier.

         This method operates in the amplifier coordinate frame.

        Parameters
        ----------
        ampData : `lsst.afw.image.ImageF`
            Amplifier image to operate on.
        scale : `float`
            Peak flux of the source to add.
        x0 : `float`
            X-coordinate of the source peak.
        y0 : `float`
            Y-coordinate of the source peak.
        """
        for x in range(0, ampData.getDimensions().getX()):
            for y in range(0, ampData.getDimensions().getY()):
                ampData.array[y][x] = (ampData.array[y][x] +
                                       scale * np.exp(-0.5 * ((x - x0)**2 + (y - y0)**2) / 3.0**2))

    def amplifierAddCT(self, ampDataSource, ampDataTarget, scale):
        """Add a scaled copy of an amplifier to another, simulating crosstalk.

         This method operates in the amplifier coordinate frame.

        Parameters
        ----------
        ampDataSource : `lsst.afw.image.ImageF`
            Amplifier image to add scaled copy from.
        ampDataTarget : `lsst.afw.image.ImageF`
            Amplifier image to add scaled copy to.
        scale : `float`
            Flux scale of the copy to add to the target.

        Notes
        -----
        This simulates simple crosstalk between amplifiers.
        """
        ampDataTarget.array[:] = (ampDataTarget.array[:] +
                                  scale * ampDataSource.array[:])

    # Functional form data values.
    def amplifierAddFringe(self, amp, ampData, scale, x0=100, y0=0):
        """Add a fringe-like ripple pattern to an amplifier's image data.

        Parameters
        ----------
        amp : `~lsst.afw.ampInfo.AmpInfoRecord`
            Amplifier to operate on. Needed for amp<->exp coordinate transforms.
        ampData : `lsst.afw.image.ImageF`
            Amplifier image to operate on.
        scale : `numpy.array` or `float`
            Peak intensity scaling for the ripple.
        x0 : `numpy.array` or `float`, optional
            Fringe center
        y0 : `numpy.array` or `float`, optional
            Fringe center

        Notes
        -----
        This uses an offset sinc function to generate a ripple
        pattern. True fringes have much finer structure, but this
        pattern should be visually identifiable. The (x, y)
        coordinates are in the frame of the amplifier, and (u, v) in
        the frame of the full trimmed image.
        """
        for x in range(0, ampData.getDimensions().getX()):
            for y in range(0, ampData.getDimensions().getY()):
                (u, v) = self.localCoordToExpCoord(amp, x, y)
                ampData.getArray()[y][x] = np.sum((ampData.getArray()[y][x] +
                                                   scale *
                                                   np.sinc(((u - x0) / 50)**2 +
                                                           ((v - y0) / 50)**2)))

    def amplifierMultiplyFlat(self, amp, ampData, fracDrop, u0=100.0, v0=100.0):
        """Multiply an amplifier's image data by a flat-like pattern.

        Parameters
        ----------
        amp : `lsst.afw.ampInfo.AmpInfoRecord`
            Amplifier to operate on. Needed for amp<->exp coordinate transforms.
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

        sigma = u0 / np.sqrt(-2.0 * np.log(fracDrop))

        for x in range(0, ampData.getDimensions().getX()):
            for y in range(0, ampData.getDimensions().getY()):
                (u, v) = self.localCoordToExpCoord(amp, x, y)
                f = np.exp(-0.5 * ((u - u0)**2 + (v - v0)**2) / sigma**2)
                ampData.array[y][x] = (ampData.array[y][x] * f)


class RawMock(IsrMock):
    """Generate a raw exposure suitable for ISR.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.isTrimmed = False
        self.config.doGenerateImage = True
        self.config.doGenerateAmpDict = False
        self.config.doAddOverscan = True
        self.config.doAddSky = True
        self.config.doAddSource = True
        self.config.doAddCrosstalk = False
        self.config.doAddBias = True
        self.config.doAddDark = True


class TrimmedRawMock(RawMock):
    """Generate a trimmed raw exposure.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.isTrimmed = True
        self.config.doAddOverscan = False


class CalibratedRawMock(RawMock):
    """Generate a trimmed raw exposure.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.isTrimmed = True
        self.config.doGenerateImage = True
        self.config.doAddOverscan = False
        self.config.doAddSky = True
        self.config.doAddSource = True
        self.config.doAddCrosstalk = False

        self.config.doAddBias = False
        self.config.doAddDark = False
        self.config.doAddFlat = False
        self.config.doAddFringe = True

        self.config.biasLevel = 0.0
        self.config.readNoise = 10.0


class RawDictMock(RawMock):
    """Generate a raw exposure dict suitable for ISR.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doGenerateAmpDict = True


class MasterMock(IsrMock):
    """Parent class for those that make master calibrations.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.isTrimmed = True
        self.config.doGenerateImage = True
        self.config.doAddOverscan = False
        self.config.doAddSky = False
        self.config.doAddSource = False
        self.config.doAddCrosstalk = False

        self.config.doAddBias = False
        self.config.doAddDark = False
        self.config.doAddFlat = False
        self.config.doAddFringe = False


class BiasMock(MasterMock):
    """Simulated master bias calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doAddBias = True
        self.config.readNoise = 10.0


class DarkMock(MasterMock):
    """Simulated master dark calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doAddDark = True
        self.config.darkTime = 1.0


class FlatMock(MasterMock):
    """Simulated master flat calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doAddFlat = True


class FringeMock(MasterMock):
    """Simulated master fringe calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.doAddFringe = True


class UntrimmedFringeMock(FringeMock):
    """Simulated untrimmed master fringe calibration.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config.isTrimmed = False


class BfKernelMock(IsrMock):
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


class DefectMock(IsrMock):
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


class CrosstalkCoeffMock(IsrMock):
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


class TransmissionMock(IsrMock):
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


class DataRefMock(object):
    """Simulated gen2 butler data ref.

    Currently only supports get and put operations, which are most
    likely to be called for data in ISR processing.

    """
    dataId = "isrMock Fake Data"
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
            self.config = IsrMockConfig()
        self.config.doGenerateImage = True
        self.config.doGenerateData = False

    def expectData(self):
        if self.config is None:
            self.config = IsrMockConfig()
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
        mock : IsrMock.run() result
            The output product.
        """
        if "_filename" in dataType:
            self.expectData()
            return tempfile.mktemp(), "mock"
        elif 'transmission_' in dataType:
            self.expectData()
            return TransmissionMock(config=self.config).run()
        elif dataType == 'ccdExposureId':
            self.expectData()
            return 20090913
        elif dataType == 'camera':
            self.expectData()
            return IsrMock(config=self.config).getCamera()
        elif dataType == 'raw':
            self.expectImage()
            return RawMock(config=self.config).run()
        elif dataType == 'bias':
            self.expectImage()
            return BiasMock(config=self.config).run()
        elif dataType == 'dark':
            self.expectImage()
            return DarkMock(config=self.config).run()
        elif dataType == 'flat':
            self.expectImage()
            return FlatMock(config=self.config).run()
        elif dataType == 'fringe':
            self.expectImage()
            return FringeMock(config=self.config).run()
        elif dataType == 'defects':
            self.expectData()
            return DefectMock(config=self.config).run()
        elif dataType == 'bfKernel':
            self.expectData()
            return BfKernelMock(config=self.config).run()
        elif dataType == 'linearizer':
            return None
        elif dataType == 'crosstalkSources':
            return None
        else:
            raise RuntimeError("ISR DataRefMock cannot return %s.", dataType)

    def put(self, exposure, filename):
        """Write an exposure to a FITS file.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Image data to write out.
        filename : `str`
            Base name of the output file.
        """
        exposure.writeFits(filename+".fits")


class FringeDataRefMock(object):
    """Simulated gen2 butler data ref.

    Currently only supports get and put operations, which are most
    likely to be called for data in ISR processing.

    """
    dataId = "isrMock Fake Data"
    darkval = 2.  # e-/sec
    oscan = 250.  # DN
    gradient = .10
    exptime = 15  # seconds
    darkexptime = 40.  # seconds

    def __init__(self, **kwargs):
        if 'config' in kwargs.keys():
            self.config = kwargs['config']
        else:
            self.config = IsrMockConfig()
            self.config.isTrimmed = True
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
        mock : IsrMock.run() result
            The output product.
        """
        if "_filename" in dataType:
            return tempfile.mktemp(), "mock"
        elif 'transmission_' in dataType:
            return TransmissionMock(config=self.config).run()
        elif dataType == 'ccdExposureId':
            return 20090913
        elif dataType == 'camera':
            return IsrMock(config=self.config).getCamera()
        elif dataType == 'raw':
            return CalibratedRawMock(config=self.config).run()
        elif dataType == 'bias':
            return BiasMock(config=self.config).run()
        elif dataType == 'dark':
            return DarkMock(config=self.config).run()
        elif dataType == 'flat':
            return FlatMock(config=self.config).run()
        elif dataType == 'fringe':
            fringes = []
            configCopy = copy.deepcopy(self.config)
            for scale, x, y in zip(self.config.fringeScale, self.config.fringeX0, self.config.fringeY0):
                configCopy.fringeScale = [1.0]
                configCopy.fringeX0 = [x]
                configCopy.fringeY0 = [y]
                fringes.append(FringeMock(config=configCopy).run())
            return fringes
        elif dataType == 'defects':
            return DefectMock(config=self.config).run()
        elif dataType == 'bfKernel':
            return BfKernelMock(config=self.config).run()
        elif dataType == 'linearizer':
            return None
        elif dataType == 'crosstalkSources':
            return None
        else:
            return None

    def put(self, exposure, filename):
        """Write an exposure to a FITS file.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Image data to write out.
        filename : `str`
            Base name of the output file.
        """
        exposure.writeFits(filename+".fits")
