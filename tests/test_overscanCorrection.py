#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

import unittest
import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as cameraGeom
import lsst.ip.isr as ipIsr
import lsst.pipe.base as pipeBase


def computeImageMedianAndStd(image):
    """Function to calculate median and std of image data.

    Parameters
    ----------
    image : `lsst.afw.image.Image`
        Image to measure statistics on.

    Returns
    -------
    median : `float`
        Image median.
    std : `float`
        Image stddev.
    """
    median = np.nanmedian(image.getArray())
    std = np.nanstd(image.getArray())

    return (median, std)


class IsrTestCases(lsst.utils.tests.TestCase):

    def updateConfigFromKwargs(self, config, **kwargs):
        """Common config from keywords.
        """
        fitType = kwargs.get('fitType', None)
        if fitType:
            config.overscan.fitType = fitType

        order = kwargs.get('order', None)
        if order:
            config.overscan.order = order

    def makeExposure(self, isTransposed=False):
        # Define the camera geometry we'll use.
        cameraBuilder = cameraGeom.Camera.Builder("Fake Camera")
        detectorBuilder = cameraBuilder.add("Fake amp", 0)

        ampBuilder = cameraGeom.Amplifier.Builder()

        dataBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                   lsst.geom.Extent2I(10, 10))

        if isTransposed is True:
            fullBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                       lsst.geom.Point2I(9, 12))
            overscanBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 10),
                                           lsst.geom.Point2I(9, 12))
        else:
            fullBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                       lsst.geom.Point2I(12, 9))

            overscanBBox = lsst.geom.Box2I(lsst.geom.Point2I(10, 0),
                                           lsst.geom.Point2I(12, 9))

        ampBuilder.setRawBBox(fullBBox)
        ampBuilder.setRawSerialOverscanBBox(overscanBBox)
        ampBuilder.setRawDataBBox(dataBBox)

        detectorBuilder.append(ampBuilder)
        camera = cameraBuilder.finish()
        detector = camera[0]

        # Define image data.
        maskedImage = afwImage.MaskedImageF(fullBBox)
        maskedImage.set(10, 0x0, 1)

        overscan = afwImage.MaskedImageF(maskedImage, overscanBBox)
        overscan.set(2, 0x0, 1)

        exposure = afwImage.ExposureF(maskedImage, None)
        exposure.setDetector(detector)
        return exposure

    def checkOverscanCorrectionY(self, **kwargs):
        exposure = self.makeExposure(isTransposed=True)
        detector = exposure.getDetector()

        # These subimages are needed below.
        overscan = exposure[detector.getAmplifiers()[0].getRawSerialOverscanBBox()]
        maskedImage = exposure[detector.getAmplifiers()[0].getRawBBox()]

        config = ipIsr.IsrTask.ConfigClass()
        self.updateConfigFromKwargs(config, **kwargs)

        if kwargs['fitType'] == "MEDIAN_PER_ROW":
            # Add a bad point to test outlier rejection.
            overscan.getImage().getArray()[0, 0] = 12345

            # Shrink the sigma clipping limit to handle the fact that the
            # bad point is not be rejected at higher thresholds (2/0.74).
            config.overscan.numSigmaClip = 2.7

        isrTask = ipIsr.IsrTask(config=config)
        isrTask.overscan.run(exposure, detector.getAmplifiers()[0], isTransposed=True)

        height = maskedImage.getHeight()
        width = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                if j == 10 and i == 0 and kwargs['fitType'] == "MEDIAN_PER_ROW":
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 12343)
                elif j >= 10:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 0)
                else:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 8)

    def checkOverscanCorrectionX(self, **kwargs):
        exposure = self.makeExposure(isTransposed=False)
        detector = exposure.getDetector()

        # These subimages are needed below.
        maskedImage = exposure[detector.getAmplifiers()[0].getRawBBox()]

        config = ipIsr.IsrTask.ConfigClass()
        self.updateConfigFromKwargs(config, **kwargs)

        isrTask = ipIsr.IsrTask(config=config)
        isrTask.overscan.run(exposure, detector.getAmplifiers()[0], isTransposed=False)

        height = maskedImage.getHeight()
        width = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                if i >= 10:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 0)
                else:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 8)

    def checkOverscanCorrectionSineWave(self, **kwargs):
        """vertical sine wave along long direction"""
        # Define the camera geometry we'll use.
        cameraBuilder = cameraGeom.Camera.Builder("Fake Camera")
        detectorBuilder = cameraBuilder.add("Fake amp", 0)

        ampBuilder = cameraGeom.Amplifier.Builder()

        dataBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                   lsst.geom.Extent2I(70, 500))

        fullBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                   lsst.geom.Extent2I(100, 500))

        overscanBBox = lsst.geom.Box2I(lsst.geom.Point2I(70, 0),
                                       lsst.geom.Extent2I(30, 500))

        ampBuilder.setRawBBox(fullBBox)
        ampBuilder.setRawSerialOverscanBBox(overscanBBox)
        ampBuilder.setRawDataBBox(dataBBox)

        detectorBuilder.append(ampBuilder)
        camera = cameraBuilder.finish()
        detector = camera[0]

        # Define image data.
        maskedImage = afwImage.MaskedImageF(fullBBox)
        maskedImage.set(50, 0x0, 1)

        overscan = afwImage.MaskedImageF(maskedImage, overscanBBox)
        overscan.set(0, 0x0, 1)

        exposure = afwImage.ExposureF(maskedImage, None)
        exposure.setDetector(detector)

        # vertical sine wave along long direction
        x = np.linspace(0, 2*3.14159, 500)
        a, w = 15, 5*3.14159
        sineWave = 20 + a*np.sin(w*x)
        sineWave = sineWave.astype(int)

        fullImage = np.repeat(sineWave, 100).reshape((500, 100))
        maskedImage.image.array += fullImage

        config = ipIsr.IsrTask.ConfigClass()
        self.updateConfigFromKwargs(config, **kwargs)

        isrTask = ipIsr.IsrTask(config=config)
        isrTask.overscan.run(exposure, detector.getAmplifiers()[0])

        height = maskedImage.getHeight()
        width = maskedImage.getWidth()

        for j in range(height):
            for i in range(width):
                if i >= 70:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 0.0)
                else:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 50.0)

    def test_MedianPerRowOverscanCorrection(self):
        self.checkOverscanCorrectionY(fitType="MEDIAN_PER_ROW")
        self.checkOverscanCorrectionX(fitType="MEDIAN_PER_ROW")
        self.checkOverscanCorrectionSineWave(fitType="MEDIAN_PER_ROW")

    def test_MedianOverscanCorrection(self):
        self.checkOverscanCorrectionY(fitType="MEDIAN")
        self.checkOverscanCorrectionX(fitType="MEDIAN")

    def checkPolyOverscanCorrectionX(self, **kwargs):
        exposure = self.makeExposure(isTransposed=False)
        detector = exposure.getDetector()

        # These subimages are needed below.
        overscan = exposure[detector.getAmplifiers()[0].getRawSerialOverscanBBox()]
        maskedImage = exposure[detector.getAmplifiers()[0].getRawBBox()]

        bbox = detector.getAmplifiers()[0].getRawSerialOverscanBBox()
        overscan.getMaskedImage().set(2, 0x0, 1)
        for i in range(bbox.getDimensions()[1]):
            for j, off in enumerate([-0.5, 0.0, 0.5]):
                overscan.image[j, i, afwImage.LOCAL] = 2+i+off

        config = ipIsr.IsrTask.ConfigClass()
        self.updateConfigFromKwargs(config, **kwargs)

        isrTask = ipIsr.IsrTask(config=config)
        isrTask.overscan.run(exposure, detector.getAmplifiers()[0], isTransposed=False)

        height = maskedImage.getHeight()
        width = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                if i == 10:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], -0.5)
                elif i == 11:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 0)
                elif i == 12:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 0.5)
                else:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 10 - 2 - j)

    def checkPolyOverscanCorrectionY(self, **kwargs):
        exposure = self.makeExposure(isTransposed=True)
        detector = exposure.getDetector()

        # These subimages are needed below.
        overscan = exposure[detector.getAmplifiers()[0].getRawSerialOverscanBBox()]
        maskedImage = exposure[detector.getAmplifiers()[0].getRawBBox()]

        bbox = detector.getAmplifiers()[0].getRawSerialOverscanBBox()
        overscan.getMaskedImage().set(2, 0x0, 1)
        for i in range(bbox.getDimensions()[0]):
            for j, off in enumerate([-0.5, 0.0, 0.5]):
                overscan.image[i, j, afwImage.LOCAL] = 2+i+off
        # maskedImage.getMaskedImage().set(10, 0x0, 1)

        config = ipIsr.IsrTask.ConfigClass()
        self.updateConfigFromKwargs(config, **kwargs)

        isrTask = ipIsr.IsrTask(config=config)
        isrTask.overscan.run(exposure, detector.getAmplifiers()[0], isTransposed=True)

        height = maskedImage.getHeight()
        width = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                if j == 10:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], -0.5)
                elif j == 11:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 0)
                elif j == 12:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 0.5)
                else:
                    self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 10 - 2 - i)

    def test_PolyOverscanCorrection(self):
        for fitType in ("POLY", "CHEB", "LEG"):
            self.checkPolyOverscanCorrectionX(fitType=fitType, order=5)
            self.checkPolyOverscanCorrectionY(fitType=fitType, order=5)

    def test_SplineOverscanCorrection(self):
        for fitType in ("NATURAL_SPLINE", "CUBIC_SPLINE", "AKIMA_SPLINE"):
            self.checkPolyOverscanCorrectionX(fitType=fitType, order=5)
            self.checkPolyOverscanCorrectionY(fitType=fitType, order=5)

    def test_overscanCorrection(self):
        """Expect that this should reduce the image variance with a full fit.
        The default fitType of MEDIAN will reduce the median value.

        This needs to operate on a RawMock() to have overscan data to use.

        The output types may be different when fitType != MEDIAN.
        """
        exposure = self.makeExposure(isTransposed=False)
        detector = exposure.getDetector()
        amp = detector.getAmplifiers()[0]

        statBefore = computeImageMedianAndStd(exposure.image[amp.getRawDataBBox()])

        config = ipIsr.IsrTask.ConfigClass()
        isrTask = ipIsr.IsrTask(config=config)
        oscanResults = isrTask.overscan.run(exposure, amp)

        self.assertIsInstance(oscanResults, pipeBase.Struct)
        self.assertIsInstance(oscanResults.imageFit, float)
        self.assertIsInstance(oscanResults.overscanFit, float)
        self.assertIsInstance(oscanResults.overscanImage, afwImage.ExposureF)

        statAfter = computeImageMedianAndStd(exposure.image[amp.getRawDataBBox()])
        self.assertLess(statAfter[0], statBefore[0])

    def test_overscanCorrection_isNotInt(self):
        """Expect smaller median/smaller std after.
        Expect exception if overscan fit type isn't known.
        """
        exposure = self.makeExposure(isTransposed=False)
        detector = exposure.getDetector()
        amp = detector.getAmplifiers()[0]

        for fitType in ('MEAN', 'MEDIAN', 'MEDIAN_PER_ROW', 'MEANCLIP', 'POLY', 'CHEB',
                        'NATURAL_SPLINE', 'CUBIC_SPLINE'):
            if fitType in ('NATURAL_SPLINE', 'CUBIC_SPLINE'):
                order = 3
            else:
                order = 1
                config = ipIsr.IsrTask.ConfigClass()
                config.overscan.order = order
                config.overscan.fitType = fitType
                isrTask = ipIsr.IsrTask(config=config)

            response = isrTask.overscan.run(exposure, amp)

            self.assertIsInstance(response, pipeBase.Struct,
                                  msg=f"overscanCorrection overscanIsNotInt Bad response: {fitType}")
            self.assertIsNotNone(response.imageFit,
                                 msg=f"overscanCorrection overscanIsNotInt Bad imageFit: {fitType}")
            self.assertIsNotNone(response.overscanFit,
                                 msg=f"overscanCorrection overscanIsNotInt Bad overscanFit: {fitType}")
            self.assertIsInstance(response.overscanImage, afwImage.ExposureF,
                                  msg=f"overscanCorrection overscanIsNotInt Bad overscanImage: {fitType}")


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
