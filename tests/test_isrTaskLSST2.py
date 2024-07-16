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

import lsst.ip.isr.isrMockLSST as isrMockLSST
import lsst.utils.tests
from lsst.ip.isr.isrTaskLSST import (IsrTaskLSST, IsrTaskLSSTConfig)
from lsst.ip.isr.crosstalk import CrosstalkCalib
from lsst.ip.isr import PhotonTransferCurveDataset


class IsrTaskLSSTTestCase(lsst.utils.tests.TestCase):
    """Test IsrTaskLSST"""
    def setUp(self):
        mock = isrMockLSST.IsrMockLSST()
        self.camera = mock.getCamera()
        self.detector = self.camera[mock.config.detectorIndex]
        self.namp = len(self.detector)

        # Create calibration frames
        self.bias = isrMockLSST.BiasMockLSST().run()
        self.dark = isrMockLSST.DarkMockLSST().run()
        self.flat = isrMockLSST.FlatMockLSST().run()
        self.bf_kernel = isrMockLSST.BfKernelMockLSST().run()

        self.crosstalk = CrosstalkCalib(nAmp=self.namp)
        self.crosstalk.hasCrosstalk = True
        self.crosstalk.coeffs = isrMockLSST.CrosstalkCoeffMockLSST().run()

        self.defects = isrMockLSST.DefectMockLSST().run()

        amp_names = [x.getName() for x in self.detector.getAmplifiers()]
        self.ptc = PhotonTransferCurveDataset(amp_names,
                                              ptcFitType='DUMMY_PTC',
                                              covMatrixSide=1)

        # TODO: check units of ptc noise
        for amp_name in amp_names:
            self.ptc.gain[amp_name] = mock.config.gainDict.get(amp_name, mock.config.gain)
            self.ptc.noise[amp_name] = mock.config.readNoise * mock.config.gain

        # TODO:
        # self.cti = isrMockLSST.DeferredChargeMockLSST().run()

        # TODO:
        # self.linearizer = ???

    def test_isrBootstrapBias(self):
        """Test processing of a ``bootstrap`` bias frame."""
        mock_config = self.get_mock_config_no_signal()

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = self.get_isr_config_minimal_corrections()
        isr_config.doBias = True

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(input_exp.clone(), bias=self.bias)

        # Rerun without doing the bias correction.
        isr_config.doBias = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(input_exp.clone())

        good_pixels = self.get_non_defect_pixels(result.exposure.mask)

        self.assertLess(
            np.mean(result.exposure.image.array[good_pixels]),
            np.mean(result2.exposure.image.array[good_pixels]),
        )
        self.assertLess(
            np.std(result.exposure.image.array[good_pixels]),
            np.std(result2.exposure.image.array[good_pixels]),
        )

        delta = result2.exposure.image.array - result.exposure.image.array
        self.assertFloatsAlmostEqual(delta[good_pixels], self.bias.image.array[good_pixels], atol=1e-5)

    def test_isrBootstrapDark(self):
        """Test processing of a ``bootstrap`` dark frame."""
        mock_config = self.get_mock_config_no_signal()
        mock_config.doAddDark = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = self.get_isr_config_minimal_corrections()
        isr_config.doBias = True
        isr_config.doDark = True

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(input_exp.clone(), bias=self.bias, dark=self.dark)

        # Rerun without doing the dark correction.
        isr_config.doDark = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(input_exp.clone(), bias=self.bias)

        good_pixels = self.get_non_defect_pixels(result.exposure.mask)

        self.assertLess(
            np.mean(result.exposure.image.array[good_pixels]),
            np.mean(result2.exposure.image.array[good_pixels]),
        )
        # The mock dark has no noise, so these should be equal.
        self.assertFloatsAlmostEqual(
            np.std(result.exposure.image.array[good_pixels]),
            np.std(result2.exposure.image.array[good_pixels]),
            atol=1e-6,
        )

        delta = result2.exposure.image.array - result.exposure.image.array
        exp_time = input_exp.getInfo().getVisitInfo().getExposureTime()
        self.assertFloatsAlmostEqual(
            delta[good_pixels],
            self.dark.image.array[good_pixels] * exp_time,
            atol=1e-5,
        )

    def test_isrBootstrapFlat(self):
        """Test processing of a ``bootstrap`` flat frame."""
        mock_config = self.get_mock_config_no_signal()
        mock_config.doAddDark = True
        mock_config.doAddFlat = True
        # The doAddSky option adds the equivalent of flat-field flux.
        mock_config.doAddSky = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = self.get_isr_config_minimal_corrections()
        isr_config.doBias = True
        isr_config.doDark = True
        isr_config.doFlat = True

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(input_exp.clone(), bias=self.bias, dark=self.dark, flat=self.flat)

        # Rerun without doing the flat correction.
        isr_config.doFlat = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(input_exp.clone(), bias=self.bias, dark=self.dark)

        good_pixels = self.get_non_defect_pixels(result.exposure.mask)

        # Applying the flat will increase the counts.
        self.assertGreater(
            np.mean(result.exposure.image.array[good_pixels]),
            np.mean(result2.exposure.image.array[good_pixels]),
        )
        # And will decrease the sigma.
        self.assertLess(
            np.std(result.exposure.image.array[good_pixels]),
            np.std(result2.exposure.image.array[good_pixels]),
        )

        ratio = result2.exposure.image.array / result.exposure.image.array
        self.assertFloatsAlmostEqual(ratio[good_pixels], self.flat.image.array[good_pixels], atol=1e-5)

    def test_isrBias(self):
        """Test processing of a bias frame."""
        mock_config = self.get_mock_config_no_signal()

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = self.get_isr_config_electronic_corrections()
        isr_config.doBias = True
        # We do not do defect correction when processing biases.
        isr_config.doDefect = False

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(input_exp.clone(), bias=self.bias, crosstalk=self.crosstalk)

        # Rerun without doing the bias correction.
        isr_config.doBias = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(input_exp.clone(), crosstalk=self.crosstalk)

        good_pixels = self.get_non_defect_pixels(result.exposure.mask)

        self.assertLess(
            np.mean(result.exposure.image.array[good_pixels]),
            np.mean(result2.exposure.image.array[good_pixels]),
        )
        self.assertLess(
            np.std(result.exposure.image.array[good_pixels]),
            np.std(result2.exposure.image.array[good_pixels]),
        )

        delta = result2.exposure.image.array - result.exposure.image.array
        self.assertFloatsAlmostEqual(delta[good_pixels], self.bias.image.array[good_pixels], atol=1e-5)

    def test_isrDark(self):
        """Test processing of a dark frame."""
        mock_config = self.get_mock_config_no_signal()
        mock_config.doAddDark = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = self.get_isr_config_electronic_corrections()
        isr_config.doBias = True
        isr_config.doDark = True
        # We do not do defect correction when processing darks.
        isr_config.doDefect = False

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(input_exp.clone(), bias=self.bias, dark=self.dark, crosstalk=self.crosstalk)

        # Rerun without doing the dark correction.
        isr_config.doDark = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(input_exp.clone(), bias=self.bias, crosstalk=self.crosstalk)

        good_pixels = self.get_non_defect_pixels(result.exposure.mask)

        self.assertLess(
            np.mean(result.exposure.image.array[good_pixels]),
            np.mean(result2.exposure.image.array[good_pixels]),
        )
        # The mock dark has no noise, so these should be equal.
        self.assertFloatsAlmostEqual(
            np.std(result.exposure.image.array[good_pixels]),
            np.std(result2.exposure.image.array[good_pixels]),
            atol=1e-6,
        )

        delta = result2.exposure.image.array - result.exposure.image.array
        exp_time = input_exp.getInfo().getVisitInfo().getExposureTime()
        self.assertFloatsAlmostEqual(
            delta[good_pixels],
            self.dark.image.array[good_pixels] * exp_time,
            atol=1e-5,
        )

    def test_isrFlat(self):
        """Test processing of a flat frame."""
        mock_config = self.get_mock_config_no_signal()
        mock_config.doAddDark = True
        mock_config.doAddFlat = True
        # The doAddSky option adds the equivalent of flat-field flux.
        mock_config.doAddSky = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = self.get_isr_config_electronic_corrections()
        isr_config.doBias = True
        isr_config.doDark = True
        isr_config.doFlat = True
        # Although we usually do not do defect interpolation when
        # processing flats, this is a good test of the interpolation.
        isr_config.doDefect = True

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(
            input_exp.clone(),
            bias=self.bias,
            dark=self.dark,
            flat=self.flat,
            crosstalk=self.crosstalk,
            defects=self.defects,
        )

        # Rerun without doing the bias correction.
        isr_config.doFlat = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(
            input_exp.clone(),
            bias=self.bias,
            dark=self.dark,
            crosstalk=self.crosstalk,
            defects=self.defects,
        )

        # With defect correction, we should not need to filter out bad
        # pixels.

        # Applying the flat will increase the counts.
        self.assertGreater(
            np.mean(result.exposure.image.array),
            np.mean(result2.exposure.image.array),
        )
        # And will decrease the sigma.
        self.assertLess(
            np.std(result.exposure.image.array),
            np.std(result2.exposure.image.array),
        )

        # Generate a flat without any defects for comparison
        # (including interpolation)
        flat_nodefect_config = isrMockLSST.FlatMockLSST.ConfigClass()
        flat_nodefect_config.doAddBrightDefects = False
        flat_nodefects = isrMockLSST.FlatMockLSST(config=flat_nodefect_config).run()

        ratio = result2.exposure.image.array / result.exposure.image.array
        self.assertFloatsAlmostEqual(ratio, flat_nodefects.image.array, atol=1e-4)

    def test_isrSkyImage(self):
        """Test processing of a sky image."""
        mock_config = self.get_mock_config_no_signal()
        mock_config.doAddDark = True
        mock_config.doAddFlat = True
        # Set this to False until we have fringe correction.
        mock_config.doAddFringe = False
        mock_config.doAddSky = True
        mock_config.doAddSource = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = self.get_isr_config_electronic_corrections()
        isr_config.doBias = True
        isr_config.doDark = True
        isr_config.doFlat = True

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(
            input_exp.clone(),
            bias=self.bias,
            dark=self.dark,
            flat=self.flat,
            crosstalk=self.crosstalk,
            defects=self.defects,
            ptc=self.ptc,
        )

        clean_mock_config = self.get_mock_config_clean()
        # We want the dark noise for more direct comparison.
        clean_mock_config.doAddDarkNoiseOnly = True
        clean_mock_config.doAddSky = True
        clean_mock_config.doAddSource = True

        clean_mock = isrMockLSST.IsrMockLSST(config=clean_mock_config)
        clean_exp = clean_mock.run()

        delta = result.exposure.image.array - clean_exp.image.array

        good_pixels = self.get_non_defect_pixels(result.exposure.mask)

        # We compare the good pixels in the entirety.
        self.assertLess(np.std(delta[good_pixels]), 5.0)
        self.assertLess(np.max(np.abs(delta[good_pixels])), 5.0*5)

        # And overall where the interpolation is a bit worse but
        # the statistics are still fine.
        self.assertLess(np.std(delta), 5.1)

    def test_isrSkyImageSaturated(self):
        """Test processing of a sky image.

        This variation uses saturated pixels instead of defects.
        """
        mock_config = self.get_mock_config_no_signal()
        mock_config.doAddDark = True
        mock_config.doAddFlat = True
        # Set this to False until we have fringe correction.
        mock_config.doAddFringe = False
        mock_config.doAddSky = True
        mock_config.doAddSource = True
        mock_config.brightDefectLevel = 50000.0  # Above saturation.

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = self.get_isr_config_electronic_corrections()
        isr_config.doBias = True
        isr_config.doDark = True
        isr_config.doFlat = True
        # We turn off defect masking to test the saturation code.
        # However, the same pixels below should be masked/interpolated.
        isr_config.doDefect = False

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(
            input_exp.clone(),
            bias=self.bias,
            dark=self.dark,
            flat=self.flat,
            crosstalk=self.crosstalk,
            defects=self.defects,
            ptc=self.ptc,
        )

        clean_mock_config = self.get_mock_config_clean()
        # We want the dark noise for more direct comparison.
        clean_mock_config.doAddDarkNoiseOnly = True
        clean_mock_config.doAddSky = True
        clean_mock_config.doAddSource = True

        clean_mock = isrMockLSST.IsrMockLSST(config=clean_mock_config)
        clean_exp = clean_mock.run()

        delta = result.exposure.image.array - clean_exp.image.array

        good_pixels = self.get_non_defect_pixels(result.exposure.mask)

        # We compare the good pixels in the entirety.
        self.assertLess(np.std(delta[good_pixels]), 5.0)
        self.assertLess(np.max(np.abs(delta[good_pixels])), 5.0*5)

        # And overall where the interpolation is a bit worse but
        # the statistics are still fine.
        self.assertLess(np.std(delta), 5.1)

    def get_mock_config_no_signal(self):
        """Get an IsrMockLSSTConfig with all signal set to False.

        This will have all the electronic effects turned on (including
        2D bias).
        """
        mock_config = isrMockLSST.IsrMockLSSTConfig()
        mock_config.isTrimmed = False
        mock_config.doAddDark = False
        mock_config.doAddFlat = False
        mock_config.doAddFringe = False
        mock_config.doAddSky = False
        mock_config.doAddSource = False

        mock_config.doAdd2DBias = True
        mock_config.doAddBias = True
        mock_config.doAddCrosstalk = True
        mock_config.doAddBrightDefects = True
        mock_config.doAddClockInjectedOffset = True
        mock_config.doAddParallelOverscanRamp = True
        mock_config.doAddSerialOverscanRamp = True
        mock_config.doApplyGain = True
        mock_config.doRoundADU = True
        # NOTE: additional electronic effects (BF, CTI, Linearity) should
        # be added here when they are supported.

        # We always want to generate the image with these configs.
        mock_config.doGenerateImage = True

        return mock_config

    def get_mock_config_clean(self):
        """Get an IsrMockLSSTConfig trimmed with all electronic signatures
        turned off.
        """
        mock_config = isrMockLSST.IsrMockLSSTConfig()
        mock_config.doAddBias = False
        mock_config.doAdd2DBias = False
        mock_config.doAddClockInjectedOffset = False
        mock_config.doAddDark = False
        mock_config.doAddDarkNoiseOnly = False
        mock_config.doAddFlat = False
        mock_config.doAddFringe = False
        mock_config.doAddSky = False
        mock_config.doAddSource = False
        mock_config.doRoundADU = False
        mock_config.doApplyGain = False
        mock_config.doAddCrosstalk = False
        mock_config.doAddBrightDefects = False
        mock_config.doAddParallelOverscanRamp = False
        mock_config.doAddSerialOverscanRamp = False

        mock_config.isTrimmed = True
        mock_config.doGenerateImage = True

        return mock_config

    def get_isr_config_minimal_corrections(self):
        """Get an IsrTaskLSSTConfig with minimal corrections.
        """
        isr_config = IsrTaskLSSTConfig()
        isr_config.doBias = False
        isr_config.doDark = False
        isr_config.doDeferredCharge = False
        isr_config.doLinearize = False
        isr_config.doCorrectGains = False
        isr_config.doCrosstalk = False
        isr_config.doDefect = False
        isr_config.doBrighterFatter = False
        isr_config.doFlat = False
        # We override the leading/trailing to skip here because of the limited
        # size of the test camera overscan regions.
        defaultAmpConfig = isr_config.overscanCamera.getOverscanDetectorConfig(self.detector).defaultAmpConfig
        defaultAmpConfig.doSerialOverscan = True
        defaultAmpConfig.serialOverscanConfig.leadingToSkip = 0
        defaultAmpConfig.serialOverscanConfig.trailingToSkip = 0
        defaultAmpConfig.doParallelOverscan = True
        defaultAmpConfig.parallelOverscanConfig.leadingToSkip = 0
        defaultAmpConfig.parallelOverscanConfig.trailingToSkip = 0

        isr_config.doAssembleCcd = True

        return isr_config

    def get_isr_config_electronic_corrections(self):
        """Get an IsrTaskLSSTConfig with electronic corrections.

        This tests all the corrections that we support in the mocks/ISR.
        """
        isr_config = IsrTaskLSSTConfig()
        # We add these as appropriate in the tests.
        isr_config.doBias = False
        isr_config.doDark = False
        isr_config.doFlat = False

        # These are the electronic effects the tests support (in addition
        # to overscan).
        isr_config.doCrosstalk = True
        isr_config.doDefect = True

        # These are the electronic effects we do not support in tetss yet.
        isr_config.doDeferredCharge = False
        isr_config.doLinearize = False
        isr_config.doCorrectGains = False
        isr_config.doBrighterFatter = False

        # We override the leading/trailing to skip here because of the limited
        # size of the test camera overscan regions.
        defaultAmpConfig = isr_config.overscanCamera.getOverscanDetectorConfig(self.detector).defaultAmpConfig
        defaultAmpConfig.doSerialOverscan = True
        defaultAmpConfig.serialOverscanConfig.leadingToSkip = 0
        defaultAmpConfig.serialOverscanConfig.trailingToSkip = 0
        defaultAmpConfig.doParallelOverscan = True
        defaultAmpConfig.parallelOverscanConfig.leadingToSkip = 0
        defaultAmpConfig.parallelOverscanConfig.trailingToSkip = 0

        isr_config.doAssembleCcd = True

        return isr_config

    def get_non_defect_pixels(self, mask_origin):
        """Get the non-defect pixels to compare.

        Parameters
        ----------
        mask_origin : `lsst.afw.image.MaskX`
            The origin mask (for shape and type).

        Returns
        -------
        pix_x, pix_y : `tuple` [`np.ndarray`]
            x and y values of good pixels.
        """
        mask_temp = mask_origin.clone()
        mask_temp[:, :] = 0

        for defect in self.defects:
            mask_temp[defect.getBBox()] = 1

        return np.where(mask_temp.array == 0)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
