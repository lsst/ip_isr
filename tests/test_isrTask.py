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

import unittest
import numpy as np

import lsst.afw.image as afwImage
import lsst.ip.isr.isrMock as isrMock
import lsst.utils.tests
from lsst.ip.isr.isrTask import (IsrTask, IsrTaskConfig)
from lsst.ip.isr.isrQa import IsrQaConfig
from lsst.pipe.base import Struct


def countMaskedPixels(maskedImage, maskPlane):
    """Function to count the number of masked pixels of a given type.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to measure the mask on.
    maskPlane : `str`
        Name of the mask plane to count

    Returns
    -------
    nMask : `int`
        Number of masked pixels.
    """
    bitMask = maskedImage.getMask().getPlaneBitMask(maskPlane)
    isBit = maskedImage.getMask().getArray() & bitMask > 0
    numBit = np.sum(isBit)
    return numBit


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


class IsrTaskTestCases(lsst.utils.tests.TestCase):
    """Test IsrTask methods with trimmed raw data.
    """
    def setUp(self):
        self.config = IsrTaskConfig()
        self.config.qa = IsrQaConfig()
        self.task = IsrTask(config=self.config)
        self.dataRef = isrMock.DataRefMock()
        self.camera = isrMock.IsrMock().getCamera()

        self.inputExp = isrMock.TrimmedRawMock().run()
        self.amp = self.inputExp.getDetector()[0]
        self.mi = self.inputExp.getMaskedImage()

    def validateIsrData(self, results):
        """results should be a struct with components that are
        not None if included in the configuration file.
        """
        self.assertIsInstance(results, Struct)
        if self.config.doBias is True:
            self.assertIsNotNone(results.bias)
        if self.config.doDark is True:
            self.assertIsNotNone(results.dark)
        if self.config.doFlat is True:
            self.assertIsNotNone(results.flat)
        if self.config.doFringe is True:
            self.assertIsNotNone(results.fringes)
        if self.config.doDefect is True:
            self.assertIsNotNone(results.defects)
        if self.config.doBrighterFatter is True:
            self.assertIsNotNone(results.bfKernel)
        if self.config.doAttachTransmissionCurve is True:
            self.assertIsNotNone(results.opticsTransmission)
            self.assertIsNotNone(results.filterTransmission)
            self.assertIsNotNone(results.sensorTransmission)
            self.assertIsNotNone(results.atmosphereTransmission)

    def test_readIsrData_noTrans(self):
        """Test that all necessary calibration frames are retrieved.
        """
        self.config.doAttachTransmissionCurve = False
        self.task = IsrTask(config=self.config)
        results = self.task.readIsrData(self.dataRef, self.inputExp)
        self.validateIsrData(results)

    def test_readIsrData_withTrans(self):
        """Test that all necessary calibration frames are retrieved.
        """
        self.config.doAttachTransmissionCurve = True
        self.task = IsrTask(config=self.config)
        results = self.task.readIsrData(self.dataRef, self.inputExp)
        self.validateIsrData(results)

    def test_ensureExposure(self):
        """Test that an exposure has a usable instance class.
        """
        self.assertIsInstance(self.task.ensureExposure(self.inputExp, self.camera, 0),
                              afwImage.Exposure)

    def test_convertItoF(self):
        """Test conversion from integer to floating point pixels.
        """
        self.assertIsInstance(self.task.convertIntToFloat(self.inputExp).getImage()[1, 1],
                              float)

    def test_updateVariance(self):
        """Expect The variance image should have a larger median value after
        this operation.
        """
        statBefore = computeImageMedianAndStd(self.inputExp.variance[self.amp.getBBox()])
        self.task.updateVariance(self.inputExp, self.amp)
        statAfter = computeImageMedianAndStd(self.inputExp.variance[self.amp.getBBox()])
        self.assertGreater(statAfter[0], statBefore[0])

    def test_darkCorrection(self):
        """Expect the median image value should decrease after this operation.
        """
        darkIm = isrMock.DarkMock().run()

        statBefore = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.task.darkCorrection(self.inputExp, darkIm)
        statAfter = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.assertLess(statAfter[0], statBefore[0])

    def test_darkCorrection_noVisitInfo(self):
        """Expect the median image value should decrease after this operation.
        """
        darkIm = isrMock.DarkMock().run()
        darkIm.getInfo().setVisitInfo(None)

        statBefore = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.task.darkCorrection(self.inputExp, darkIm)
        statAfter = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.assertLess(statAfter[0], statBefore[0])

    def test_flatCorrection(self):
        """Expect the image median should increase (divide by < 1).
        """
        flatIm = isrMock.FlatMock().run()

        statBefore = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.task.flatCorrection(self.inputExp, flatIm)
        statAfter = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.assertGreater(statAfter[1], statBefore[1])

    def test_saturationDetection(self):
        """Expect the saturation level detection/masking to scale with
        threshold.
        """
        self.amp.setSaturation(1000.0)
        self.task.saturationDetection(self.inputExp, self.amp)
        countBefore = countMaskedPixels(self.mi, "SAT")

        self.amp.setSaturation(25.0)
        self.task.saturationDetection(self.inputExp, self.amp)
        countAfter = countMaskedPixels(self.mi, "SAT")

        self.assertLessEqual(countBefore, countAfter)

    def test_measureBackground(self):
        """Expect the background measurement runs successfully and to save
        metadata values.
        """
        self.config.qa.flatness.meshX = 20
        self.config.qa.flatness.meshY = 20
        self.task.measureBackground(self.inputExp, self.config.qa)
        self.assertIsNotNone(self.inputExp.getMetadata().getScalar('SKYLEVEL'))

    def test_flatContext(self):
        """Expect the flat context manager runs successfully and leaves the
        image data the same.
        """
        darkExp = isrMock.DarkMock().run()
        flatExp = isrMock.FlatMock().run()

        mi = self.inputExp.getMaskedImage().clone()
        with self.task.flatContext(self.inputExp, flatExp, darkExp):
            self.assertTrue(True)

        self.assertMaskedImagesAlmostEqual(mi, self.inputExp.getMaskedImage(), -3)


class IsrTaskUnTrimmedTestCases(lsst.utils.tests.TestCase):
    """Test IsrTask methods using untrimmed raw data.
    """
    def setUp(self):
        self.config = IsrTaskConfig()
        self.config.qa = IsrQaConfig()
        self.task = IsrTask(config=self.config)

        self.mockConfig = isrMock.IsrMockConfig()
        self.mockConfig.isTrimmed = False
        self.doGenerateImage = True
        self.dataRef = isrMock.DataRefMock(config=self.mockConfig)
        self.camera = isrMock.IsrMock(config=self.mockConfig).getCamera()

        self.inputExp = isrMock.RawMock(config=self.mockConfig).run()
        self.amp = self.inputExp.getDetector()[0]
        self.mi = self.inputExp.getMaskedImage()

    def batchSetConfiguration(self, value):
        """Set the configuration state to a consistent value.

        Disable options we do not need as well.

        Parameters
        ----------
        value : `bool`
            Value to switch common ISR configuration options to.
        """
        self.config.qa.flatness.meshX = 20
        self.config.qa.flatness.meshY = 20
        self.config.doWrite = False
        self.config.doLinearize = False
        self.config.doCrosstalk = False

        self.config.doConvertIntToFloat = value
        self.config.doSaturation = value
        self.config.doSuspect = value
        self.config.doSetBadRegions = value
        self.config.doOverscan = value
        self.config.doBias = value
        self.config.doVariance = value
        self.config.doWidenSaturationTrails = value
        self.config.doBrighterFatter = value
        self.config.doDefect = value
        self.config.doSaturationInterpolation = value
        self.config.doDark = value
        self.config.doStrayLight = value
        self.config.doFlat = value
        self.config.doFringe = value
        self.config.doAddDistortionModel = value
        self.config.doMeasureBackground = value
        self.config.doVignette = value
        self.config.doAttachTransmissionCurve = value
        self.config.doUseOpticsTransmission = value
        self.config.doUseFilterTransmission = value
        self.config.doUseSensorTransmission = value
        self.config.doUseAtmosphereTransmission = value
        self.config.qa.saveStats = value
        self.config.qa.doThumbnailOss = value
        self.config.qa.doThumbnailFlattened = value

        self.config.doApplyGains = not value
        self.config.doCameraSpecificMasking = value
        self.config.vignette.doWriteVignettePolygon = value

    def validateIsrResults(self):
        """results should be a struct with components that are
        not None if included in the configuration file.

        Returns
        -------
        results : `pipeBase.Struct`
            Results struct generated from the current ISR configuration.
        """
        self.task = IsrTask(config=self.config)
        results = self.task.run(self.inputExp,
                                camera=self.camera,
                                bias=self.dataRef.get("bias"),
                                dark=self.dataRef.get("dark"),
                                flat=self.dataRef.get("flat"),
                                bfKernel=self.dataRef.get("bfKernel"),
                                defects=self.dataRef.get("defects"),
                                fringes=Struct(fringes=self.dataRef.get("fringe"), seed=1234),
                                opticsTransmission=self.dataRef.get("transmission_"),
                                filterTransmission=self.dataRef.get("transmission_"),
                                sensorTransmission=self.dataRef.get("transmission_"),
                                atmosphereTransmission=self.dataRef.get("transmission_")
                                )

        self.assertIsInstance(results, Struct)
        self.assertIsInstance(results.exposure, afwImage.Exposure)
        return results

    def test_overscanCorrection(self):
        """Expect that this should reduce the image variance with a full fit.
        The default fitType of MEDIAN will reduce the median value.

        This needs to operate on a RawMock() to have overscan data to use.

        The output types may be different when fitType != MEDIAN.
        """
        statBefore = computeImageMedianAndStd(self.inputExp.image[self.amp.getRawDataBBox()])

        oscanResults = self.task.overscanCorrection(self.inputExp, self.amp)
        self.assertIsInstance(oscanResults, Struct)
        self.assertIsInstance(oscanResults.imageFit, float)
        self.assertIsInstance(oscanResults.overscanFit, float)
        self.assertIsInstance(oscanResults.overscanImage, afwImage.MaskedImageF)

        statAfter = computeImageMedianAndStd(self.inputExp.image[self.amp.getRawDataBBox()])
        self.assertLess(statAfter[0], statBefore[0])

    def test_runDataRef(self):
        """Expect a dataRef to be handled correctly.
        """
        self.config.doLinearize = False
        self.config.doWrite = False
        self.task = IsrTask(config=self.config)
        print(self.task.config)
        results = self.task.runDataRef(self.dataRef)

        self.assertIsInstance(results, Struct)
        self.assertIsInstance(results.exposure, afwImage.Exposure)

    def test_run_allTrue(self):
        """Expect successful run with expected outputs when all non-exclusive
        configuration options are on.

        Output results should be tested more precisely by the
        individual function tests.

        """
        self.batchSetConfiguration(True)
        self.validateIsrResults()

    def test_run_allFalse(self):
        """Expect successful run with expected outputs when all non-exclusive
        configuration options are off.

        Output results should be tested more precisely by the
        individual function tests.

        """
        self.batchSetConfiguration(False)
        self.validateIsrResults()

    def test_failCases(self):
        """Expect failure with crosstalk enabled.

        Output results should be tested more precisely by the
        individual function tests.
        """
        self.batchSetConfiguration(True)

        # This breaks it
        self.config.doCrosstalk = True

        with self.assertRaises(RuntimeError):
            self.validateIsrResults()

    def test_maskingCase_noMasking(self):
        """Test masking cases of configuration parameters.
        """
        self.batchSetConfiguration(True)
        self.config.overscanFitType = "POLY"
        self.config.overscanOrder = 1

        self.config.doSaturation = False
        self.config.doWidenSaturationTrails = False
        self.config.doSaturationInterpolation = False
        self.config.doSuspect = False
        self.config.doSetBadRegions = False
        self.config.doDefect = False
        self.config.doBrighterFatter = False

        results = self.validateIsrResults()

        self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 0)

    def test_maskingCase_satMasking(self):
        """Test masking cases of configuration parameters.
        """
        self.batchSetConfiguration(True)
        self.config.overscanFitType = "POLY"
        self.config.overscanOrder = 1

        self.config.saturation = 20000.0
        self.config.doSaturation = True
        self.config.doWidenSaturationTrails = True

        self.config.doSaturationInterpolation = False
        self.config.doSuspect = False
        self.config.doSetBadRegions = False
        self.config.doDefect = False
        self.config.doBrighterFatter = False

        results = self.validateIsrResults()

        self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 0)

    def test_maskingCase_satMaskingAndInterp(self):
        """Test masking cases of configuration parameters.
        """
        self.batchSetConfiguration(True)
        self.config.overscanFitType = "POLY"
        self.config.overscanOrder = 1

        self.config.saturation = 20000.0
        self.config.doSaturation = True
        self.config.doWidenSaturationTrails = True
        self.config.doSaturationInterpolation = True

        self.config.doSuspect = False
        self.config.doSetBadRegions = False
        self.config.doDefect = False
        self.config.doBrighterFatter = False

        results = self.validateIsrResults()

        self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 0)

    def test_maskingCase_throughEdge(self):
        """Test masking cases of configuration parameters.
        """
        self.batchSetConfiguration(True)
        self.config.overscanFitType = "POLY"
        self.config.overscanOrder = 1

        self.config.saturation = 20000.0
        self.config.doSaturation = True
        self.config.doWidenSaturationTrails = True
        self.config.doSaturationInterpolation = True
        self.config.numEdgeSuspect = 5
        self.config.doSuspect = True

        self.config.doSetBadRegions = False
        self.config.doDefect = False
        self.config.doBrighterFatter = False

        results = self.validateIsrResults()

        self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 0)

    def test_maskingCase_throughDefects(self):
        """Test masking cases of configuration parameters.
        """
        self.batchSetConfiguration(True)
        self.config.overscanFitType = "POLY"
        self.config.overscanOrder = 1

        self.config.saturation = 20000.0
        self.config.doSaturation = True
        self.config.doWidenSaturationTrails = True
        self.config.doSaturationInterpolation = True
        self.config.numEdgeSuspect = 5
        self.config.doSuspect = True
        self.config.doDefect = True

        self.config.doSetBadRegions = False
        self.config.doBrighterFatter = False

        results = self.validateIsrResults()

        self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 2000)
        self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 3940)
        self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 2000)

    def test_maskingCase_throughBad(self):
        """Test masking cases of configuration parameters.
        """
        self.batchSetConfiguration(True)
        self.config.overscanFitType = "POLY"
        self.config.overscanOrder = 1

        self.config.saturation = 20000.0
        self.config.doSaturation = True
        self.config.doWidenSaturationTrails = True
        self.config.doSaturationInterpolation = True

        self.config.doSuspect = True
        self.config.doDefect = True
        self.config.doSetBadRegions = True
        self.config.doBrighterFatter = False

        results = self.validateIsrResults()

        self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 2000)
        self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 0)
        self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 2000)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
