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
from lsst.ip.isr.isrTaskLSST import (IsrTaskLSST, IsrTaskLSSTConfig)
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


class IsrTaskLSSTTestCases(lsst.utils.tests.TestCase):
    """Test IsrTaskLSST methods with trimmed raw data.
    """
    def setUp(self):
        print('New instance')
        self.config = IsrTaskLSSTConfig()
        self.config.overscan.doParallelOverscan = True
        self.config.qa = IsrQaConfig()
        self.task = IsrTaskLSST(config=self.config)
        self.camera = isrMock.IsrMock().getCamera()

        self.inputExp = isrMock.TrimmedRawMock().run()
        self.amp = self.inputExp.getDetector()[0]
        self.mi = self.inputExp.getMaskedImage()
        self.detector = self.inputExp.getDetector()

    def validateIsrData(self, results):
        """results should be a struct with components that are
        not None if included in the configuration file.
        """
        self.assertIsInstance(results, Struct)
        if self.config.doBias is True:
            self.assertIsNotNone(results.bias)
        if self.config.doDark is True:
            self.assertIsNotNone(results.dark)
        if self.config.doDefect is True:
            self.assertIsNotNone(results.defects)
        if self.config.doBrighterFatter is True:
            self.assertIsNotNone(results.bfKernel)

    def test_updateVariance(self):
        """Expect The variance image should have a larger median value after
        this operation.
        """
        statBefore = computeImageMedianAndStd(self.inputExp.variance[self.amp.getBBox()])
        self.task.updateVariance(self.inputExp, self.amp)
        statAfter = computeImageMedianAndStd(self.inputExp.variance[self.amp.getBBox()])
        self.assertGreater(statAfter[0], statBefore[0])
        self.assertFloatsAlmostEqual(statBefore[0], 0.0, atol=1e-2)
        self.assertFloatsAlmostEqual(statAfter[0], 8170.0195, atol=1e-2)

    def test_darkCorrection(self):
        """Expect the median image value should decrease after this operation.
        """
        darkIm = isrMock.DarkMock().run()

        statBefore = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.task.darkCorrection(self.inputExp, darkIm)
        statAfter = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.assertLess(statAfter[0], statBefore[0])
        self.assertFloatsAlmostEqual(statBefore[0], 8070.0195, atol=1e-2)
        self.assertFloatsAlmostEqual(statAfter[0], 8045.7773, atol=1e-2)

    def test_darkCorrection_noVisitInfo(self):
        """Expect the median image value should decrease after this operation.
        """
        darkIm = isrMock.DarkMock().run()
        darkIm.getInfo().setVisitInfo(None)

        statBefore = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.task.darkCorrection(self.inputExp, darkIm)
        statAfter = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.assertLess(statAfter[0], statBefore[0])
        self.assertFloatsAlmostEqual(statBefore[0], 8070.0195, atol=1e-2)
        self.assertFloatsAlmostEqual(statAfter[0], 8045.7773, atol=1e-2)

    def test_flatCorrection(self):
        """Expect the image median should increase (divide by < 1).
        """
        flatIm = isrMock.FlatMock().run()

        statBefore = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.task.flatCorrection(self.inputExp, flatIm)
        statAfter = computeImageMedianAndStd(self.inputExp.image[self.amp.getBBox()])
        self.assertGreater(statAfter[1], statBefore[1])
        self.assertFloatsAlmostEqual(statAfter[1], 147407.02, atol=1e-2)
        self.assertFloatsAlmostEqual(statBefore[1], 147.55304, atol=1e-2)

    # def test_saturationDetection(self):
    #     """Expect the saturation level detection/masking to scale with
    #     threshold.
    #     """
    #     ampB = self.amp.rebuild()
    #     ampB.setSaturation(9000.0)
    #     self.task.saturationDetection(self.inputExp, ampB.finish())
    #     countBefore = countMaskedPixels(self.mi, "SAT")

    #     ampB.setSaturation(8250.0)
    #     self.task.saturationDetection(self.inputExp, ampB.finish())
    #     countAfter = countMaskedPixels(self.mi, "SAT")

    #     self.assertLessEqual(countBefore, countAfter)
    #     self.assertEqual(countBefore, 43)
    #     self.assertEqual(countAfter, 136)

    def test_flatContext(self):
        """Expect the flat context manager runs successfully (applying both
        flat and dark within the context), and results in the same
        image data after completion.
        """
        darkExp = isrMock.DarkMock().run()
        flatExp = isrMock.FlatMock().run()

        mi = self.inputExp.getMaskedImage().clone()
        with self.task.flatContext(self.inputExp, flatExp, darkExp):
            contextStat = computeImageMedianAndStd(self.inputExp.getMaskedImage().getImage())
            self.assertFloatsAlmostEqual(contextStat[0], 37165.594, atol=1e-2)

        self.assertMaskedImagesAlmostEqual(mi, self.inputExp.getMaskedImage())

    # def test_failCases(self):
    #     """Expect failure with crosstalk enabled.

    #     Output results should be tested more precisely by the
    #     individual function tests.
    #     """
    #     self.batchSetConfiguration(True)

    #     # This breaks it
    #     self.config.doCrosstalk = True

    #     with self.assertRaises(RuntimeError):
    #         self.validateIsrResults()

    # def test_maskingCase_negativeVariance(self):
    #     """Test masking cases of configuration parameters.
    #     """
    #     self.batchSetConfiguration(True)
    #     self.config.overscan.doParallelOverscan = False
    #     self.config.overscan.fitType = "POLY"
    #     self.config.overscan.order = 1

    #     self.config.doSaturation = False
    #     self.config.doWidenSaturationTrails = False
    #     self.config.doSaturationInterpolation = False
    #     self.config.doSuspect = False
    #     self.config.doSetBadRegions = False
    #     self.config.doDefect = False
    #     self.config.doBrighterFatter = False

    #     self.config.maskNegativeVariance = True
    #     self.config.doInterpolate = False

    #     results = self.validateIsrResults()

    #     self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 40800)

    # def test_maskingCase_noMasking(self):
    #     """Test masking cases of configuration parameters.
    #     """
    #     self.batchSetConfiguration(True)
    #     self.config.overscan.fitType = "POLY"
    #     self.config.overscan.order = 1

    #     self.config.doSaturation = False
    #     self.config.doWidenSaturationTrails = False
    #     self.config.doSaturationInterpolation = False
    #     self.config.doSuspect = False
    #     self.config.doSetBadRegions = False
    #     self.config.doDefect = False
    #     self.config.doBrighterFatter = False

    #     self.config.maskNegativeVariance = False
    #     self.config.doInterpolate = False

    #     results = self.validateIsrResults()

    #     self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 0)

    # def test_maskingCase_satMasking(self):
    #     """Test masking cases of configuration parameters.
    #     """
    #     self.batchSetConfiguration(True)
    #     self.config.overscan.fitType = "POLY"
    #     self.config.overscan.order = 1

    #     self.config.saturation = 20000.0
    #     self.config.doSaturation = True
    #     self.config.doWidenSaturationTrails = True

    #     self.config.doSaturationInterpolation = False
    #     self.config.doSuspect = False
    #     self.config.doSetBadRegions = False
    #     self.config.doDefect = False
    #     self.config.doBrighterFatter = False

    #     self.config.maskNegativeVariance = False  # These are mock images.

    #     results = self.validateIsrResults()

    #     self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 0)

    # def test_maskingCase_satMaskingAndInterp(self):
    #     """Test masking cases of configuration parameters.
    #     """
    #     self.batchSetConfiguration(True)
    #     self.config.overscan.fitType = "POLY"
    #     self.config.overscan.order = 1

    #     self.config.saturation = 20000.0
    #     self.config.doSaturation = True
    #     self.config.doWidenSaturationTrails = True
    #     self.config.doSaturationInterpolation = True

    #     self.config.doSuspect = False
    #     self.config.doSetBadRegions = False
    #     self.config.doDefect = False
    #     self.config.doBrighterFatter = False

    #     self.config.maskNegativeVariance = False  # These are mock images.

    #     results = self.validateIsrResults()

    #     self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 0)

    # def test_maskingCase_throughEdge(self):
    #     """Test masking cases of configuration parameters.
    #     """
    #     self.batchSetConfiguration(True)
    #     self.config.overscan.fitType = "POLY"
    #     self.config.overscan.order = 1

    #     self.config.saturation = 20000.0
    #     self.config.doSaturation = True
    #     self.config.doWidenSaturationTrails = True
    #     self.config.doSaturationInterpolation = True
    #     self.config.numEdgeSuspect = 5
    #     self.config.doSuspect = True

    #     self.config.doSetBadRegions = False
    #     self.config.doDefect = False
    #     self.config.doBrighterFatter = False

    #     self.config.maskNegativeVariance = False  # These are mock images.

    #     results = self.validateIsrResults()

    #     self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 0)

    # def test_maskingCase_throughDefects(self):
    #     """Test masking cases of configuration parameters.
    #     """
    #     self.batchSetConfiguration(True)
    #     self.config.overscan.fitType = "POLY"
    #     self.config.overscan.order = 1

    #     self.config.saturation = 20000.0
    #     self.config.doSaturation = True
    #     self.config.doWidenSaturationTrails = True
    #     self.config.doSaturationInterpolation = True
    #     self.config.numEdgeSuspect = 5
    #     self.config.doSuspect = True
    #     self.config.doDefect = True

    #     self.config.doSetBadRegions = False
    #     self.config.doBrighterFatter = False

    #     self.config.maskNegativeVariance = False  # These are mock images.

    #     results = self.validateIsrResults()

    #     self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 2000)
    #     self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 3940)
    #     self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 2000)

    # def test_maskingCase_throughDefectsAmpEdges(self):
    #     """Test masking cases of configuration parameters.
    #     """
    #     self.batchSetConfiguration(True)
    #     self.config.overscan.fitType = "POLY"
    #     self.config.overscan.order = 1

    #     self.config.saturation = 20000.0
    #     self.config.doSaturation = True
    #     self.config.doWidenSaturationTrails = True
    #     self.config.doSaturationInterpolation = True
    #     self.config.numEdgeSuspect = 5
    #     self.config.doSuspect = True
    #     self.config.doDefect = True
    #     self.config.edgeMaskLevel = 'AMP'

    #     self.config.doSetBadRegions = False
    #     self.config.doBrighterFatter = False

    #     self.config.maskNegativeVariance = False  # These are mock images.

    #     results = self.validateIsrResults()

    #     self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 2000)
    #     self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 11280)
    #     self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 2000)

    # def test_maskingCase_throughBad(self):
    #     """Test masking cases of configuration parameters.
    #     """
    #     self.batchSetConfiguration(True)
    #     self.config.overscan.fitType = "POLY"
    #     self.config.overscan.order = 1

    #     self.config.saturation = 20000.0
    #     self.config.doSaturation = True
    #     self.config.doWidenSaturationTrails = True
    #     self.config.doSaturationInterpolation = True

    #     self.config.doSuspect = True
    #     self.config.doDefect = True
    #     self.config.doSetBadRegions = True
    #     self.config.doBrighterFatter = False

    #     self.config.maskNegativeVariance = False  # These are mock images.

    #     results = self.validateIsrResults()

    #     self.assertEqual(countMaskedPixels(results.exposure, "SAT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "INTRP"), 2000)
    #     self.assertEqual(countMaskedPixels(results.exposure, "SUSPECT"), 0)
    #     self.assertEqual(countMaskedPixels(results.exposure, "BAD"), 2000)

    # def test_binnedExposures(self):
    #     """Ensure that binned exposures have correct sizes."""
    #     self.batchSetConfiguration(True)
    #     self.config.doBinnedExposures = True
    #     self.config.binFactor1 = 8
    #     self.config.binFactor2 = 64

    #     results = self.validateIsrResults()

    #     original = results.exposure.image.array.shape
    #     bin1 = results.outputBin1Exposure.image.array.shape
    #     bin2 = results.outputBin2Exposure.image.array.shape

    #     # Binning truncates, so check that the original and obtained
    #     # binned images have the correct offset.
    #     self.assertEqual(original[0] - bin1[0] * 8, 4)
    #     self.assertEqual(original[1] - bin1[1] * 8, 0)
    #     self.assertEqual(original[0] - bin2[0] * 64, 12)
    #     self.assertEqual(original[1] - bin2[1] * 64, 8)


class IsrTaskLSSTUnTrimmedTestCases(lsst.utils.tests.TestCase):
    """Test IsrTask methods using untrimmed raw data,
    testing pre-assembly and assembly
    """
    def setUp(self):
        self.config = IsrTaskLSSTConfig()
        self.config.overscan.doParallelOverscan = True
        self.config.qa = IsrQaConfig()
        self.task = IsrTaskLSST(config=self.config)

        self.mockConfig = isrMock.IsrMockConfig()
        self.mockConfig.isTrimmed = False
        self.doGenerateImage = True
        self.dataContainer = isrMock.MockDataContainer(config=self.mockConfig)
        self.camera = isrMock.IsrMock(config=self.mockConfig).getCamera()

        self.inputExp = isrMock.RawMock(config=self.mockConfig).run()
        self.amp = self.inputExp.getDetector()[0]
        self.mi = self.inputExp.getMaskedImage()
        #TO DO initiate ptc data sets and fill gains with 1.5

    def batchSetConfiguration(self, value):
        """Set the configuration state to a consistent value.

        Disable options we do not need as well.

        Parameters
        ----------
        value : `bool`
            Value to switch common ISR configuration options to.
        """
        self.config.doDiffNonLinearCorrection = value
        self.config.doOverscan = value
        self.config.doAssembleCcd = value

        self.config.doSetBadRegions = False
        self.config.doBias = False
        self.config.doVariance = False
        self.config.doWidenSaturationTrails = False
        self.config.doBrighterFatter = False
        self.config.doDefect = False
        self.config.doSaturationInterpolation = False
        self.config.doDark = False
        self.config.qa.saveStats = False

    def validateIsrResults(self):
        """results should be a struct with components that are
        not None if included in the configuration file.

        Returns
        -------
        results : `pipeBase.Struct`
            Results struct generated from the current ISR configuration.
        """
        self.task = IsrTaskLSST(config=self.config)
        results = self.task.run(self.inputExp,
                                camera=self.camera,
                                bias=self.dataContainer.get("bias"),
                                dark=self.dataContainer.get("dark"),
                                flat=self.dataContainer.get("flat"),
                                bfKernel=self.dataContainer.get("bfKernel"),
                                defects=self.dataContainer.get("defects")
                                )

        self.assertIsInstance(results, Struct)
        self.assertIsInstance(results.exposure, afwImage.Exposure)
        return results

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


    def test_overscanCorrection(self):
        self.task.overscanCorrection(self.detector,self.inputExp)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main(failfast=True)
