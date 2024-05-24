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

    def updateOverscanConfigFromKwargs(self, config, **kwargs):
        """Common config from keywords.
        """
        fitType = kwargs.get('fitType', None)
        if fitType:
            config.fitType = fitType

        order = kwargs.get('order', None)
        if order:
            config.order = order

    def makeExposure(self, addRamp=False, isTransposed=False):
        # Define the camera geometry we'll use.
        cameraBuilder = cameraGeom.Camera.Builder("Fake Camera")
        detectorBuilder = cameraBuilder.add("Fake amp", 0)

        ampBuilder = cameraGeom.Amplifier.Builder()

        dataBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                   lsst.geom.Extent2I(10, 10))

        if isTransposed is True:
            fullBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                       lsst.geom.Point2I(12, 12))
            serialOverscanBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 10),
                                                 lsst.geom.Point2I(9, 12))
            parallelOverscanBBox = lsst.geom.Box2I(lsst.geom.Point2I(10, 0),
                                                   lsst.geom.Point2I(12, 9))
        else:
            fullBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                       lsst.geom.Point2I(12, 12))
            serialOverscanBBox = lsst.geom.Box2I(lsst.geom.Point2I(10, 0),
                                                 lsst.geom.Point2I(12, 9))
            parallelOverscanBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 10),
                                                   lsst.geom.Point2I(9, 12))

        ampBuilder.setRawBBox(fullBBox)
        ampBuilder.setRawSerialOverscanBBox(serialOverscanBBox)
        ampBuilder.setRawParallelOverscanBBox(parallelOverscanBBox)
        ampBuilder.setRawDataBBox(dataBBox)

        detectorBuilder.append(ampBuilder)
        camera = cameraBuilder.finish()
        detector = camera[0]

        # Define image data.
        maskedImage = afwImage.MaskedImageF(fullBBox)
        maskedImage.set(2, 0x0, 1)

        dataImage = afwImage.MaskedImageF(maskedImage, dataBBox)
        dataImage.set(10, 0x0, 1)

        if addRamp:
            for column in range(dataBBox.getWidth()):
                maskedImage.image.array[:, column] += column

        exposure = afwImage.ExposureF(maskedImage, None)
        exposure.setDetector(detector)
        return exposure

    def checkOverscanCorrectionY(self, **kwargs):
        # We check serial overscan with the "old" and "new" tasks.
        for taskClass in (ipIsr.overscan.OverscanCorrectionTask, ipIsr.overscan.SerialOverscanCorrectionTask):
            exposure = self.makeExposure(isTransposed=True)
            detector = exposure.getDetector()

            # These subimages are needed below.
            overscan = exposure[detector.getAmplifiers()[0].getRawSerialOverscanBBox()]
            maskedImage = exposure[detector.getAmplifiers()[0].getRawBBox()]

            config = taskClass.ConfigClass()
            self.updateOverscanConfigFromKwargs(config, **kwargs)

            if kwargs['fitType'] == "MEDIAN_PER_ROW":
                # Add a bad point to test outlier rejection.
                overscan.getImage().getArray()[0, 0] = 12345

                # Shrink the sigma clipping limit to handle the fact that the
                # bad point is not be rejected at higher thresholds (2/0.74).
                config.numSigmaClip = 2.7

            overscanTask = taskClass(config=config)
            _ = overscanTask.run(exposure, detector.getAmplifiers()[0], isTransposed=True)

            height = maskedImage.getHeight()
            width = maskedImage.getWidth()
            for j in range(height):
                for i in range(width):
                    if j == 10 and i == 0 and kwargs['fitType'] == "MEDIAN_PER_ROW":
                        self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 12343)
                    elif j >= 10 and i < 10:
                        self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 0)
                    elif i < 10:
                        self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 8)

    def checkOverscanCorrectionX(self, **kwargs):
        # We check serial ovsercan with "old" and "new" tasks.
        for taskClass in (ipIsr.overscan.OverscanCorrectionTask, ipIsr.overscan.SerialOverscanCorrectionTask):
            exposure = self.makeExposure(isTransposed=False)
            detector = exposure.getDetector()

            # These subimages are needed below.
            maskedImage = exposure[detector.getAmplifiers()[0].getRawBBox()]

            config = taskClass.ConfigClass()
            self.updateOverscanConfigFromKwargs(config, **kwargs)

            overscanTask = taskClass(config=config)
            _ = overscanTask.run(exposure, detector.getAmplifiers()[0], isTransposed=False)

            height = maskedImage.getHeight()
            width = maskedImage.getWidth()
            for j in range(height):
                for i in range(width):
                    if i >= 10 and j < 10:
                        self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 0)
                    elif j < 10:
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

        for taskClass in (ipIsr.overscan.OverscanCorrectionTask, ipIsr.overscan.SerialOverscanCorrectionTask):
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

            config = taskClass.ConfigClass()
            self.updateOverscanConfigFromKwargs(config, **kwargs)

            overscanTask = taskClass(config=config)
            _ = overscanTask.run(exposure, detector.getAmplifiers()[0])

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

    def test_MeanPerRowOverscanCorrection(self):
        self.checkOverscanCorrectionY(fitType="MEAN_PER_ROW")
        self.checkOverscanCorrectionX(fitType="MEAN_PER_ROW")
        self.checkOverscanCorrectionSineWave(fitType="MEAN_PER_ROW")

    def test_MedianOverscanCorrection(self):
        self.checkOverscanCorrectionY(fitType="MEDIAN")
        self.checkOverscanCorrectionX(fitType="MEDIAN")

    def checkPolyOverscanCorrectionX(self, **kwargs):
        for taskClass in (ipIsr.overscan.OverscanCorrectionTask, ipIsr.overscan.SerialOverscanCorrectionTask):
            exposure = self.makeExposure(isTransposed=False)
            detector = exposure.getDetector()

            # Fill the full serial overscan region with a polynomial,
            # all the way into the parallel overscan region.
            amp = detector.getAmplifiers()[0]
            serialOverscanBBox = amp.getRawSerialOverscanBBox()
            imageBBox = amp.getRawDataBBox()
            parallelOverscanBBox = amp.getRawParallelOverscanBBox()
            imageBBox = imageBBox.expandedTo(parallelOverscanBBox)

            serialOverscanBBox = lsst.geom.Box2I(
                lsst.geom.Point2I(serialOverscanBBox.getMinX(),
                                  imageBBox.getMinY()),
                lsst.geom.Extent2I(serialOverscanBBox.getWidth(),
                                   imageBBox.getHeight()),
            )

            overscan = exposure[serialOverscanBBox]
            maskedImage = exposure[detector.getAmplifiers()[0].getRawBBox()]

            bbox = serialOverscanBBox
            overscan.getMaskedImage().set(2, 0x0, 1)
            for i in range(bbox.getDimensions()[1]):
                for j, off in enumerate([-0.5, 0.0, 0.5]):
                    overscan.image[j, i, afwImage.LOCAL] = 2+i+off

            config = taskClass.ConfigClass()
            self.updateOverscanConfigFromKwargs(config, **kwargs)

            overscanTask = taskClass(config=config)
            _ = overscanTask.run(exposure, detector.getAmplifiers()[0], isTransposed=False)

            height = maskedImage.getHeight()
            width = maskedImage.getWidth()
            for j in range(height):
                for i in range(width):
                    if j < 10:
                        if i == 10:
                            self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], -0.5)
                        elif i == 11:
                            self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 0)
                        elif i == 12:
                            self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 0.5)
                        else:
                            self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 10 - 2 - j)

    def checkPolyOverscanCorrectionY(self, **kwargs):
        for taskClass in (ipIsr.overscan.OverscanCorrectionTask, ipIsr.overscan.SerialOverscanCorrectionTask):
            exposure = self.makeExposure(isTransposed=True)
            detector = exposure.getDetector()

            # Fill the full serial overscan region with a polynomial,
            # all the way into the parallel overscan region.
            amp = detector.getAmplifiers()[0]
            serialOverscanBBox = amp.getRawSerialOverscanBBox()
            imageBBox = amp.getRawDataBBox()
            parallelOverscanBBox = amp.getRawParallelOverscanBBox()
            imageBBox = imageBBox.expandedTo(parallelOverscanBBox)

            serialOverscanBBox = lsst.geom.Box2I(
                lsst.geom.Point2I(serialOverscanBBox.getMinX(), imageBBox.getEndY()),
                lsst.geom.Extent2I(imageBBox.getWidth(), serialOverscanBBox.getHeight()),
            )

            overscan = exposure[serialOverscanBBox]
            maskedImage = exposure[detector.getAmplifiers()[0].getRawBBox()]

            bbox = serialOverscanBBox
            overscan.getMaskedImage().set(2, 0x0, 1)
            for i in range(bbox.getDimensions()[0]):
                for j, off in enumerate([-0.5, 0.0, 0.5]):
                    overscan.image[i, j, afwImage.LOCAL] = 2+i+off

            config = taskClass.ConfigClass()
            self.updateOverscanConfigFromKwargs(config, **kwargs)

            overscanTask = taskClass(config=config)
            _ = overscanTask.run(exposure, detector.getAmplifiers()[0], isTransposed=True)

            height = maskedImage.getHeight()
            width = maskedImage.getWidth()
            for j in range(height):
                for i in range(width):
                    if i < 10:
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
        for taskClass in (ipIsr.overscan.OverscanCorrectionTask, ipIsr.overscan.SerialOverscanCorrectionTask):
            exposure = self.makeExposure(isTransposed=False)
            detector = exposure.getDetector()
            amp = detector.getAmplifiers()[0]

            statBefore = computeImageMedianAndStd(exposure.image[amp.getRawDataBBox()])

            config = taskClass.ConfigClass()
            overscanTask = taskClass(config=config)
            oscanResults = overscanTask.run(exposure, amp)

            self.assertIsInstance(oscanResults, pipeBase.Struct)
            self.assertIsInstance(oscanResults.imageFit, float)
            self.assertIsInstance(oscanResults.overscanFit, float)
            self.assertIsInstance(oscanResults.overscanImage, afwImage.ExposureF)

            statAfter = computeImageMedianAndStd(exposure.image[amp.getRawDataBBox()])
            self.assertLess(statAfter[0], statBefore[0])

    def test_parallelOverscanCorrection(self):
        """Expect that this should reduce the image variance with a full fit.
        The default fitType of MEDIAN will reduce the median value.

        This needs to operate on a RawMock() to have overscan data to use.

        This test checks that the outputs match, and that the serial
        overscan is the trivial value (2.0), and that the parallel
        overscan is the median of the ramp inserted (4.5)
        """
        for taskType in ("combined", "separate"):
            exposure = self.makeExposure(addRamp=True, isTransposed=False)
            detector = exposure.getDetector()
            amp = detector.getAmplifiers()[0]

            statBefore = computeImageMedianAndStd(exposure.image[amp.getRawDataBBox()])

            for fitType in ('MEDIAN', 'MEDIAN_PER_ROW'):
                # This tests these two types to cover scalar and vector
                # calculations.
                exposureCopy = exposure.clone()

                if taskType == "combined":
                    config = ipIsr.overscan.OverscanCorrectionTask.ConfigClass()
                    config.doParallelOverscan = True
                    config.fitType = fitType

                    overscanTask = ipIsr.overscan.OverscanCorrectionTask(config=config)
                    oscanResults = overscanTask.run(exposureCopy, amp)
                else:
                    configSerial = ipIsr.overscan.SerialOverscanCorrectionTask.ConfigClass()
                    configSerial.fitType = fitType

                    serialOverscanTask = ipIsr.overscan.SerialOverscanCorrectionTask(config=configSerial)
                    serialResults = serialOverscanTask.run(exposureCopy, amp)

                    configParallel = ipIsr.overscan.ParallelOverscanCorrectionTask.ConfigClass()
                    configParallel.fitType = fitType

                    parallelOverscanTask = ipIsr.overscan.ParallelOverscanCorrectionTask(
                        config=configParallel,
                    )
                    oscanResults = parallelOverscanTask.run(exposureCopy, amp)

                self.assertIsInstance(oscanResults, pipeBase.Struct)
                if fitType == 'MEDIAN':
                    self.assertIsInstance(oscanResults.imageFit, float)
                    self.assertIsInstance(oscanResults.overscanFit, float)
                else:
                    self.assertIsInstance(oscanResults.imageFit, np.ndarray)
                    self.assertIsInstance(oscanResults.overscanFit, np.ndarray)
                self.assertIsInstance(oscanResults.overscanImage, afwImage.ExposureF)

                statAfter = computeImageMedianAndStd(exposureCopy.image[amp.getRawDataBBox()])
                self.assertLess(statAfter[0], statBefore[0])

                # Test the output value for the serial and parallel overscans
                if taskType == "combined":
                    self.assertAlmostEqual(oscanResults.overscanMean[0], 2.0, delta=0.001)
                    self.assertAlmostEqual(oscanResults.overscanMean[1], 4.5, delta=0.001)
                else:
                    self.assertAlmostEqual(serialResults.overscanMean, 2.0, delta=0.001)
                    self.assertAlmostEqual(oscanResults.overscanMean, 4.5, delta=0.001)

                if fitType != 'MEDIAN':
                    # The ramp that has been inserted should be fully
                    # removed by the overscan fit, removing all of the
                    # signal.  This isn't true of the constant fit, so do
                    # not test that here.
                    self.assertLess(statAfter[1], statBefore[1])
                    self.assertAlmostEqual(statAfter[1], 0.0, delta=0.001)

    def test_bleedParallelOverscanCorrection(self):
        """Expect that this should reduce the image variance with a full fit.
        The default fitType of MEDIAN will reduce the median value.

        This needs to operate on a RawMock() to have overscan data to use.

        This test adds a large artificial bleed to the overscan region,
        which should be masked and patched with the median of the
        other pixels.
        """
        for taskType in ("combined", "separate"):
            exposure = self.makeExposure(addRamp=True, isTransposed=False)
            detector = exposure.getDetector()
            amp = detector.getAmplifiers()[0]

            maskedImage = exposure.getMaskedImage()
            overscanBleedBox = lsst.geom.Box2I(lsst.geom.Point2I(4, 10),
                                               lsst.geom.Extent2I(2, 3))
            overscanBleed = afwImage.MaskedImageF(maskedImage, overscanBleedBox)
            overscanBleed.set(110000, 0x0, 1)

            statBefore = computeImageMedianAndStd(exposure.image[amp.getRawDataBBox()])

            for fitType in ('MEDIAN', 'MEDIAN_PER_ROW', 'POLY'):
                # We only test these three types as this should cover the
                # scalar calculations, the generic vector calculations,
                # and the specific C++ MEDIAN_PER_ROW case.
                exposureCopy = exposure.clone()

                if taskType == "combined":
                    config = ipIsr.overscan.OverscanCorrectionTask.ConfigClass()
                    config.doParallelOverscan = True
                    config.parallelOverscanMaskGrowSize = 1
                    config.fitType = fitType

                    overscanTask = ipIsr.overscan.OverscanCorrectionTask(config=config)
                    # This next line is usually run as part of IsrTask
                    overscanTask.maskParallelOverscan(exposureCopy, detector)
                    oscanResults = overscanTask.run(exposureCopy, amp)
                else:
                    configSerial = ipIsr.overscan.SerialOverscanCorrectionTask.ConfigClass()
                    configSerial.fitType = fitType

                    serialOverscanTask = ipIsr.overscan.SerialOverscanCorrectionTask(config=configSerial)
                    serialResults = serialOverscanTask.run(exposureCopy, amp)

                    configParallel = ipIsr.overscan.ParallelOverscanCorrectionTask.ConfigClass()
                    configParallel.parallelOverscanMaskGrowSize = 1
                    configParallel.fitType = fitType

                    parallelOverscanTask = ipIsr.overscan.ParallelOverscanCorrectionTask(
                        config=configParallel,
                    )
                    # This next line is usually run as part of IsrTask
                    parallelOverscanTask.maskParallelOverscan(exposureCopy, detector, saturationLevel=100000.)
                    oscanResults = parallelOverscanTask.run(exposureCopy, amp)

                self.assertIsInstance(oscanResults, pipeBase.Struct)
                if fitType == 'MEDIAN':
                    self.assertIsInstance(oscanResults.imageFit, float)
                    self.assertIsInstance(oscanResults.overscanFit, float)
                else:
                    self.assertIsInstance(oscanResults.imageFit, np.ndarray)
                    self.assertIsInstance(oscanResults.overscanFit, np.ndarray)
                self.assertIsInstance(oscanResults.overscanImage, afwImage.ExposureF)

                statAfter = computeImageMedianAndStd(exposureCopy.image[amp.getRawDataBBox()])
                self.assertLess(statAfter[0], statBefore[0])

                # Test the output value for the serial and parallel
                # overscans.
                if taskType == "combined":
                    self.assertAlmostEqual(oscanResults.overscanMean[0], 2.0, delta=0.001)
                    self.assertAlmostEqual(oscanResults.overscanMean[1], 4.5, delta=0.001)
                    self.assertAlmostEqual(oscanResults.residualMean[1], 0.0, delta=0.001)
                else:
                    self.assertAlmostEqual(serialResults.overscanMean, 2.0, delta=0.001)
                    self.assertAlmostEqual(oscanResults.overscanMean, 4.5, delta=0.001)
                    self.assertAlmostEqual(oscanResults.residualMean, 0.0, delta=0.001)

                if fitType != 'MEDIAN':
                    # Check the bleed isn't oversubtracted.  This is the
                    # average of the two mid-bleed pixels as the patching
                    # uses the median correction value there, and there is
                    # still a residual ramp in this region.  The large
                    # delta allows the POLY fit to pass, which has sub-ADU
                    # differences.
                    self.assertAlmostEqual(exposureCopy.image.array[5][0],
                                           0.5 * (exposureCopy.image.array[5][4]
                                                  + exposureCopy.image.array[5][5]), delta=0.3)
                    # These fits should also reduce the image stdev, as
                    # they are modeling the ramp.
                    self.assertLess(statAfter[1], statBefore[1])

    def test_bleedParallelOverscanCorrectionFailure(self):
        """Expect that this should reduce the image variance with a full fit.
        The default fitType of MEDIAN will reduce the median value.

        This needs to operate on a RawMock() to have overscan data to use.

        This adds a large artificial bleed to the overscan region,
        which should be masked and patched with the median of the
        other pixels.
        """
        for taskType in ("combined", "separate"):
            exposure = self.makeExposure(addRamp=True, isTransposed=False)
            detector = exposure.getDetector()
            amp = detector.getAmplifiers()[0]

            maskedImage = exposure.getMaskedImage()
            overscanBleedBox = lsst.geom.Box2I(lsst.geom.Point2I(4, 10),
                                               lsst.geom.Extent2I(2, 3))
            overscanBleed = afwImage.MaskedImageF(maskedImage, overscanBleedBox)
            overscanBleed.set(10000, 0x0, 1)  # This level is below the mask threshold.

            statBefore = computeImageMedianAndStd(exposure.image[amp.getRawDataBBox()])

            for fitType in ('MEDIAN', 'MEDIAN_PER_ROW'):
                # We only test these three types as this should cover the
                # scalar calculations, the generic vector calculations,
                # and the specific C++ MEDIAN_PER_ROW case.
                exposureCopy = exposure.clone()

                if taskType == "combined":
                    config = ipIsr.overscan.OverscanCorrectionTask.ConfigClass()
                    config.doParallelOverscan = True
                    config.parallelOverscanMaskGrowSize = 1
                    # Ensure we don't mask anything
                    config.maxDeviation = 100000
                    config.fitType = fitType

                    overscanTask = ipIsr.overscan.OverscanCorrectionTask(config=config)
                    oscanResults = overscanTask.run(exposureCopy, amp)

                    oscanMeanSerial, oscanMeanParallel = oscanResults.overscanMean
                    oscanMedianParallel = oscanResults.overscanMedian[1]
                else:
                    configSerial = ipIsr.overscan.SerialOverscanCorrectionTask.ConfigClass()
                    # Ensure we don't mask anything
                    configSerial.maxDeviation = 100000
                    configSerial.fitType = fitType

                    serialOverscanTask = ipIsr.overscan.SerialOverscanCorrectionTask(config=configSerial)
                    serialResults = serialOverscanTask.run(exposureCopy, amp)

                    configParallel = ipIsr.overscan.ParallelOverscanCorrectionTask.ConfigClass()
                    configParallel.maxDeviation = 100000
                    configParallel.parallelOverscanMaskGrowSize = 1
                    configParallel.fitType = fitType

                    parallelOverscanTask = ipIsr.overscan.ParallelOverscanCorrectionTask(
                        config=configParallel,
                    )
                    oscanResults = parallelOverscanTask.run(exposureCopy, amp)

                    oscanMeanSerial = serialResults.overscanMean
                    oscanMeanParallel = oscanResults.overscanMean
                    oscanMedianParallel = oscanResults.overscanMedian

                self.assertIsInstance(oscanResults, pipeBase.Struct)
                if fitType == 'MEDIAN':
                    self.assertIsInstance(oscanResults.imageFit, float)
                    self.assertIsInstance(oscanResults.overscanFit, float)
                else:
                    self.assertIsInstance(oscanResults.imageFit, np.ndarray)
                    self.assertIsInstance(oscanResults.overscanFit, np.ndarray)
                self.assertIsInstance(oscanResults.overscanImage, afwImage.ExposureF)

                statAfter = computeImageMedianAndStd(exposureCopy.image[amp.getRawDataBBox()])
                self.assertLess(statAfter[0], statBefore[0])

                # Test the output value for the serial and parallel
                # overscans.
                self.assertAlmostEqual(oscanMeanSerial, 2.0, delta=0.001)
                # These are the wrong values:
                if fitType == 'MEDIAN':
                    # Check that the constant case is now biased, at 6.5
                    # instead of 4.5:
                    self.assertAlmostEqual(oscanMeanParallel, 6.5, delta=0.001)
                else:
                    # This is not correcting the bleed, so it will be printed
                    # onto the image, making the stdev after correction worse
                    # than before.
                    self.assertGreater(statAfter[1], statBefore[1])

                    # Check that the median overscan value matches the
                    # constant fit:
                    self.assertAlmostEqual(oscanMedianParallel, 6.5, delta=0.001)
                    # Check that the mean isn't what we found before, and
                    # is larger:
                    self.assertNotEqual(oscanMeanParallel, 4.5)
                    self.assertGreater(oscanMeanParallel, 4.5)
                    self.assertGreater(exposureCopy.image.array[5][0],
                                       0.5 * (exposureCopy.image.array[5][4]
                                              + exposureCopy.image.array[5][5]))

    def test_overscanCorrection_isNotInt(self):
        """Expect smaller median/smaller std after.
        Expect exception if overscan fit type isn't known.
        """
        for taskClass in (ipIsr.overscan.OverscanCorrectionTask, ipIsr.overscan.SerialOverscanCorrectionTask):
            exposure = self.makeExposure(isTransposed=False)
            detector = exposure.getDetector()
            amp = detector.getAmplifiers()[0]

            for fitType in ('MEAN', 'MEDIAN', 'MEDIAN_PER_ROW', 'MEANCLIP', 'POLY', 'CHEB',
                            'NATURAL_SPLINE', 'CUBIC_SPLINE'):
                if fitType in ('NATURAL_SPLINE', 'CUBIC_SPLINE'):
                    order = 3
                else:
                    order = 1

                config = taskClass.ConfigClass()
                config.order = order
                config.fitType = fitType

                overscanTask = taskClass(config=config)

                response = overscanTask.run(exposure, amp)

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
