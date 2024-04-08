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
    def setUp(self):

        # Set up ISR task configs
        self.config = IsrTaskLSSTConfig()

        # Create a mock image with ISR effects in
        self.mockConfig = isrMockLSST.IsrMockLSSTConfig()

        self.camera = isrMockLSST.IsrMockLSST(config=self.mockConfig).getCamera()
        self.inputExp = isrMockLSST.RawMockLSST(config=self.mockConfig).run()
        self.detector = self.inputExp.getDetector()

        # Create a mock bias
        self.bias = isrMockLSST.BiasMockLSST(config=self.mockConfig).run()

        # Create a mock flat
        self.flat = isrMockLSST.FlatMockLSST(config=self.mockConfig).run()

        # Create a mock dark
        self.dark = isrMockLSST.DarkMockLSST(config=self.mockConfig).run()

        # Create a mock brighter-fatter kernel
        self.bfkernel = isrMockLSST.BfKernelMockLSST(config=self.mockConfig).run()

        # Create a mock crosstalk coeff matrix
        self.crosstalkCoeff = isrMockLSST.CrosstalkCoeffMockLSST(config=self.mockConfig).run()

        # Create mock defects
        self.defect = isrMockLSST.DefectMockLSST(config=self.mockConfig).run()

        # We need a mock PTC
        ampNames = [x.getName() for x in self.detector.getAmplifiers()]
        print(ampNames)
        # Make PTC.
        # The arguments other than the ampNames aren't important for this.
        self.ptc = PhotonTransferCurveDataset(ampNames,
                                              ptcFitType='DUMMY_PTC',
                                              covMatrixSide=1)
        for ampName in ampNames:
            self.ptc.gain[ampName] = 1.5  # gain in e-/ADU
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

        self.mockConfig.doAddFringe = False
        self.mockConfig.doAddFlat = False
        self.mockConfig.doAddDark = False
        self.mockConfig.doApplyGain = False
        self.mockConfig.doAddCrosstalk = False
        self.mockConfig.doAddParallelOverscan = False
        self.mockConfig.doAddSerialOverscan = False

    def setIsrConfigStepsFalse(self):
        """Set all configs corresponding to ISR steps to False.
        """

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
        # Not testing dnlLUT (not existant yet), deferred Charge,
        # linearizer, crosstalk, bfgains
        results = self.task.run(self.inputExp,
                                camera=self.camera,
                                bias=self.bias,
                                ptc=self.ptc,
                                crosstalk=self.crosstalkCoeff,
                                defects=self.defect,
                                bfKernel=self.bfkernel,
                                dark=self.dark,
                                flat=self.flat
                                )
        import IPython
        IPython.embed()

        self.assertIsInstance(results, Struct)
        self.assertIsInstance(results.exposure, afwImage.Exposure)
        return results

    def test_run_serialOverscanCorrection(self):
        """Test up to serial overscan correction.
        """
        self.setIsrConfigStepsFalse()
        self.setMockConfigFalse()
        self.inputExp = isrMockLSST.RawMockLSST(config=self.mockConfig).run()
        self.mockConfig.doAddSerialOverscan = True
        self.config.overscanCamera.getOverscanDetectorConfig(self.detector).defaultAmpConfig.doSerialOverscan = True
        results = self.validateIsrResults()


    def notest_run_parallelOverscanCrosstalkCorrection(self):
        """Test up to parallel overscan crosstalk correction.
        """
        # TODO: DM-43286
        pass

    def notest_run_parallelOverscanCorrection(self):
        """Test up to parallel overscan correction.
        """
        self.setIsrConfigStepsFalse()
        self.config.overscanCamera.getOverscanDetectorConfig(self.detector).defaultAmpConfig.doSerialOverscan = True
        self.config.overscanCamera.getOverscanDetectorConfig(self.detector).defaultAmpConfig.doParallelOverscan = True
        results = self.validateIsrResults()

    def notest_run_assembleCcd(self):
        """Test up to assembly.
        """
        self.setIsrConfigStepsFalse()
        self.config.overscanCamera.getOverscanDetectorConfig(self.detector).defaultAmpConfig.doSerialOverscan = True
        self.config.overscanCamera.getOverscanDetectorConfig(self.detector).defaultAmpConfig.doParallelOverscan = True
        self.config.doAssembleCcd = True
        results = self.validateIsrResults()

    def notest_run_linearize(self):
        """Test up to linearizer.
        """
        # TODO DM-???

    def notest_run_crosstalkCorrection(self):
        """Test up to crosstalk correction.
        """
        self.setIsrConfigStepsFalse()
        self.config.overscanCamera.getOverscanDetectorConfig(self.detector).defaultAmpConfig.doSerialOverscan = True
        self.config.overscanCamera.getOverscanDetectorConfig(self.detector).defaultAmpConfig.doParallelOverscan = True
        self.config.doAssembleCcd = True
        self.config.doCrosstalk = True
        results = self.validateIsrResults()


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main(failfast=True)
