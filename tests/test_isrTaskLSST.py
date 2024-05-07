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
import lsst.ip.isr.isrMockLSST as isrMockLSST
import lsst.utils.tests
from lsst.ip.isr.overscanAmpConfig import OverscanAmpConfig
from lsst.ip.isr.overscanAmpConfig import OverscanDetectorConfig
from lsst.ip.isr.isrTaskLSST import (IsrTaskLSST, IsrTaskLSSTConfig)
from lsst.ip.isr.crosstalk import CrosstalkCalib
from lsst.ip.isr.deferredCharge import DeferredChargeCalib
from lsst.pipe.base import Struct
from lsst.ip.isr import PhotonTransferCurveDataset


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
    """Test IsrTaskLSST pipeline step-by-step
    to produce calibrations and to correct a mock image.
    """

    def doSetUp_mock(self,config):

        # Create a mock image
        self.mock = isrMockLSST.IsrMockLSST(config=config)
        self.inputExp = self.mock.run()
        self.camera = self.mock.getCamera()
        self.detector = self.inputExp.getDetector()

        # Get number of amps
        self.namp = len(self.detector.getAmplifiers())

        # Create a mock bias
        self.bias = isrMockLSST.BiasMockLSST().run()

        # Create a mock CTI, by instatianting the class with default settings
        # so we can test the pipeline would run
        # TODO: Update with some mock CTI data DM-????
        # self.cti = DeferredChargeCalib()

        # Create a mock flat
        self.flat = isrMockLSST.FlatMockLSST().run()

        # Create a mock dark
        self.dark = isrMockLSST.DarkMockLSST().run()

        # Create a mock brighter-fatter kernel
        self.bfkernel = isrMockLSST.BfKernelMockLSST().run()

        # Create a crosstalk calib with default coefficients from the mock data
        self.crosstalk = CrosstalkCalib(nAmp=self.namp)
        self.crosstalk.hasCrosstalk = True
        self.crosstalk.coeffs = isrMockLSST.CrosstalkCoeffMockLSST().run()

        # Create mock defects
        self.defect = isrMockLSST.DefectMockLSST().run()

        # We need a mock PTC
        ampNames = [x.getName() for x in self.detector.getAmplifiers()]
        print(ampNames)
        # Make PTC.
        # The arguments other than the ampNames aren't important for this.
        self.ptc = PhotonTransferCurveDataset(ampNames,
                                              ptcFitType='DUMMY_PTC',
                                              covMatrixSide=1)
        for ampName in ampNames:
            self.ptc.gain[ampName] = 3.5  # gain in e-/ADU
            self.ptc.noise[ampName] = 8.5  # read noise in ADU

    def setMockConfigFalse(self):
        """Set all configs corresponding to ISR steps to False.
        """
        self.mockConfig.isTrimmed = False
        self.mockConfig.doGenerateImage = True
        self.mockConfig.doGenerateData = False
        self.mockConfig.doGenerateAmpDict = False

        self.mockConfig.doAddSky = True
        self.mockConfig.doAddSource = True

        self.mockConfig.doAddBias = False
        self.mockConfig.doAddFringe = False
        self.mockConfig.doAddFlat = False
        self.mockConfig.doAddDark = False
        self.mockConfig.doApplyGain = False
        self.mockConfig.doAddCrosstalk = False
        self.mockConfig.doAddParallelOverscan = False
        self.mockConfig.doAddSerialOverscan = False

    def setIsrConfig(self):
        """Set all configs corresponding to ISR steps to False.
        """

        # Set up ISR task configs
        self.defaultAmpConfig = self.config.overscanCamera.getOverscanDetectorConfig(self.detector).defaultAmpConfig

        self.config.doDiffNonLinearCorrection = False
        overscanDetectorConfig = self.config.overscanCamera.getOverscanDetectorConfig(self.detector)
        overscanDetectorConfig.defaultAmpConfig.doSerialOverscan = False
        overscanDetectorConfig.defaultAmpConfig.doParallelOverscanCrosstalk = False
        overscanDetectorConfig.defaultAmpConfig.doParallelOverscan = False
        self.config.doAssembleCcd = False
        self.config.doLinearize = False
        self.config.doCrosstalk = False
        self.config.doBias = False
        self.config.doGainsCorrection = False
        self.config.doApplyGains = False
        self.config.doDeferredCharge = False
        self.config.doVariance = False
        self.config.doDefect = False
        self.config.doNanMasking = False
        self.config.doWidenSaturationTrails = False
        self.config.doSaveInterpPixels = False
        self.config.doSetBadRegions = False
        self.config.doInterpolate = False
        self.config.doDark = False
        self.config.doBrighterFatter = False
        self.config.doFlat = False

    def validateIsrResults(self):
        """results should be a struct with components that are
        not None if included in the configuration file.

        Returns
        -------
        results : `pipeBase.Struct`
            Results struct generated from the current ISR configuration.
        """


        self.task = IsrTaskLSST(config=self.config)

        mockMean = 0.
        for amp in self.inputExp.getDetector():
            bbox = amp.getRawDataBBox()
            ampData = self.inputExp.image[bbox]
            mockMean += np.mean(ampData.array)

        # Not testing dnlLUT (not existant yet), deferred Charge,
        # linearizer, crosstalk, bfgains
        results = self.task.run(self.inputExp,
                                camera=self.camera,
                                bias=self.bias,
                                ptc=self.ptc,
                                crosstalk=self.crosstalk,
                                # deferredChargeCalib=self.cti,
                                defects=self.defect,
                                bfKernel=self.bfkernel,
                                dark=self.dark,
                                flat=self.flat
                                )

        outputMean = np.mean(results.outputExposure.image.array)
        # Test that the output has a smaller mean than the input mock
        self.assertLess(outputMean, mockMean)
        # Test that the output is a struct
        self.assertIsInstance(results, Struct)
        # Test that the output has an exposure with expected format
        self.assertIsInstance(results.exposure, afwImage.Exposure)
        return results

    def test_run_serialOverscanCorrection(self):
        """Test up to serial overscan correction.
        """
        self.mockConfig = isrMockLSST.IsrMockLSSTConfig()
        self.setMockConfigFalse()
        self.mockConfig.doAddSerialOverscan = True
        self.doSetUp_mock(config=self.mockConfig)

        self.config = IsrTaskLSSTConfig()
        self.setIsrConfig()
        self.defaultAmpConfig.doSerialOverscan = True

        results = self.validateIsrResults()

    def test_run_parallelOverscanCrosstalkCorrection(self):
        """Test up to parallel overscan crosstalk correction.
        """
        # TODO: DM-43286
        pass

    def test_run_parallelOverscanCorrection(self):
        """Test up to parallel overscan correction.
        """

        self.mockConfig = isrMockLSST.IsrMockLSSTConfig()
        self.setMockConfigFalse()
        self.mockConfig.doAddSerialOverscan = True
        self.mockConfig.doAddParallelOverscan = True
        self.doSetUp_mock(config=self.mockConfig)

        self.config = IsrTaskLSSTConfig()
        self.setIsrConfig()
        self.defaultAmpConfig.doSerialOverscan = True
        self.defaultAmpConfig.doParallelOverscan = True

        results = self.validateIsrResults()

    def test_run_linearize(self):
        """Test up to linearizer.
        """
        # TODO DM-???
        pass

    def test_run_crosstalkCorrection(self):
        """Test up to crosstalk correction.
        """
        self.mockConfig = isrMockLSST.IsrMockLSSTConfig()
        self.setMockConfigFalse()
        self.mockConfig.doAddSerialOverscan = True
        self.mockConfig.doAddParallelOverscan = True
        self.mockConfig.doAddCrosstalk = True
        self.doSetUp_mock(config=self.mockConfig)

        self.config = IsrTaskLSSTConfig()
        self.setIsrConfig()
        self.defaultAmpConfig.doSerialOverscan = True
        self.defaultAmpConfig.doParallelOverscan = True
        self.config.doAssembleCcd = True
        self.config.doCrosstalk = True

        results = self.validateIsrResults()

    def test_run_biasCorrection(self):
        """Test up to crosstalk correction.
        """
        self.mockConfig = isrMockLSST.IsrMockLSSTConfig()
        self.setMockConfigFalse()
        self.mockConfig.doAddSerialOverscan = True
        self.mockConfig.doAddParallelOverscan = True
        self.mockConfig.doAddCrosstalk = True
        self.mockConfig.doAddBias = True
        self.doSetUp_mock(config=self.mockConfig)

        self.config = IsrTaskLSSTConfig()
        self.setIsrConfig()
        self.defaultAmpConfig.doSerialOverscan = True
        self.defaultAmpConfig.doParallelOverscan = True
        self.config.doAssembleCcd = True
        self.config.doCrosstalk = True
        self.config.doBias = True

        results = self.validateIsrResults()

    def test_run_correctGains(self):
        """Test up to gains correction
        # TODO DM-??
        """
        pass

    def test_run_applyGains(self):
        """Test up to crosstalk correction.
        """
        self.mockConfig = isrMockLSST.IsrMockLSSTConfig()
        self.setMockConfigFalse()
        self.mockConfig.doAddSerialOverscan = True
        self.mockConfig.doAddParallelOverscan = True
        self.mockConfig.doAddCrosstalk = True
        self.mockConfig.doAddBias = True
        self.mockConfig.doApplyGain = True
        self.mockConfig.gain = 3.5
        self.doSetUp_mock(config=self.mockConfig)

        self.config = IsrTaskLSSTConfig()
        self.setIsrConfig()
        self.defaultAmpConfig.doSerialOverscan = True
        self.defaultAmpConfig.doParallelOverscan = True
        self.config.doAssembleCcd = True
        self.config.doCrosstalk = True
        self.config.doBias = True
        self.config.doApplyGains = True

        results = self.validateIsrResults()

    def test_run_doVarianceDefectsMasking(self):
        """Test up to masking of bad pixels.
        """
        self.mockConfig = isrMockLSST.IsrMockLSSTConfig()
        self.setMockConfigFalse()
        self.mockConfig.doAddSerialOverscan = True
        self.mockConfig.doAddParallelOverscan = True
        self.mockConfig.doAddCrosstalk = True
        self.mockConfig.doAddBias = True
        self.mockConfig.doApplyGain = True
        self.mockConfig.gain = 3.5
        self.doSetUp_mock(config=self.mockConfig)

        self.config = IsrTaskLSSTConfig()
        self.setIsrConfig()
        self.defaultAmpConfig.doSerialOverscan = True
        self.defaultAmpConfig.doParallelOverscan = True
        self.config.doAssembleCcd = True
        self.config.doCrosstalk = True
        self.config.doBias = True
        self.config.doApplyGains = True
        self.config.doVariance = True
        self.config.doDefect = True
        self.config.doNanMasking = True
        self.config.doWidenSaturationTrails = True

        results = self.validateIsrResults()

    def test_run_doDark(self):
        """Test up to dark correction.
        """
        self.mockConfig = isrMockLSST.IsrMockLSSTConfig()
        self.setMockConfigFalse()
        self.mockConfig.doAddSerialOverscan = True
        self.mockConfig.doAddParallelOverscan = True
        self.mockConfig.doAddCrosstalk = True
        self.mockConfig.doAddBias = True
        self.mockConfig.doApplyGain = True
        self.mockConfig.gain = 3.5
        self.mockConfig.doAddDark = True
        self.doSetUp_mock(config=self.mockConfig)

        self.config = IsrTaskLSSTConfig()
        self.setIsrConfig()
        self.defaultAmpConfig.doSerialOverscan = True
        self.defaultAmpConfig.doParallelOverscan = True
        self.config.doAssembleCcd = True
        self.config.doCrosstalk = True
        self.config.doBias = True
        self.config.doApplyGains = True
        self.config.doVariance = True
        self.config.doDefect = True
        self.config.doNanMasking = True
        self.config.doWidenSaturationTrails = True
        self.config.doDark = True

        results = self.validateIsrResults()

    def test_run_doBFcorrection(self):
        """Test up to do BF correction
        # TODO DM-??
        """
        pass

    def test_run_doFlat(self):
        """Test up to flat correction
        """
        self.mockConfig = isrMockLSST.IsrMockLSSTConfig()
        self.setMockConfigFalse()
        self.mockConfig.doAddSerialOverscan = True
        self.mockConfig.doAddParallelOverscan = True
        self.mockConfig.doAddCrosstalk = True
        self.mockConfig.doAddBias = True
        self.mockConfig.doApplyGain = True
        self.mockConfig.gain = 3.5
        self.mockConfig.doAddDark = True
        self.mockConfig.doAddFlat = True
        self.doSetUp_mock(config=self.mockConfig)

        self.config = IsrTaskLSSTConfig()
        self.setIsrConfig()
        self.defaultAmpConfig.doSerialOverscan = True
        self.defaultAmpConfig.doParallelOverscan = True
        self.config.doAssembleCcd = True
        self.config.doCrosstalk = True
        self.config.doBias = True
        self.config.doApplyGains = True
        self.config.doVariance = True
        self.config.doDefect = True
        self.config.doNanMasking = True
        self.config.doWidenSaturationTrails = True
        self.config.doDark = True
        self.config.doFlat = True

        results = self.validateIsrResults()

    def test_run_setPixelValues(self):
        """Test up to flat correction
        """
        self.mockConfig = isrMockLSST.IsrMockLSSTConfig()
        self.setMockConfigFalse()
        self.mockConfig.doAddSerialOverscan = True
        self.mockConfig.doAddParallelOverscan = True
        self.mockConfig.doAddCrosstalk = True
        self.mockConfig.doAddBias = True
        self.mockConfig.doApplyGain = True
        self.mockConfig.gain = 3.5
        self.mockConfig.doAddDark = True
        self.mockConfig.doAddFlat = True
        self.doSetUp_mock(config=self.mockConfig)

        self.config = IsrTaskLSSTConfig()
        self.setIsrConfig()
        self.defaultAmpConfig.doSerialOverscan = True
        self.defaultAmpConfig.doParallelOverscan = True
        self.config.doAssembleCcd = True
        self.config.doCrosstalk = True
        self.config.doBias = True
        self.config.doApplyGains = True
        self.config.doVariance = True
        self.config.doDefect = True
        self.config.doNanMasking = True
        self.config.doWidenSaturationTrails = True
        self.config.doDark = True
        self.config.doFlat = True
        self.config.doSetBadRegions = True
        self.config.doInterpolate = True

        results = self.validateIsrResults()


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main(failfast=True)
