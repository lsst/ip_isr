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

        self.defect = isrMockLSST.DefectMockLSST().run()

        amp_names = [x.getName() for x in self.detector.getAmplifiers()]
        self.ptc = PhotonTransferCurveDataset(amp_names,
                                              ptcFitType='DUMMY_PTC',
                                              covMatrixSide=1)
        # TODO: different gains ...
        # Can we simulate different gains?  Probably not yet.  FIX THIS.
        for amp_name in amp_names:
            self.ptc.gain[amp_name] = mock.config.gain
            self.ptc.noise[amp_name] = mock.config.readNoise * mock.config.gain

        # TODO:
        # self.cti = isrMockLSST.DeferredChargeMockLSST().run()

        # TODO:
        # self.linearizer = ???

    def test_isrBootstrapBias(self):
        """Test processing of a ``bootstrap`` bias frame."""
        mock_config = isrMockLSST.IsrMockLSSTConfig()
        mock_config.isTrimmed = False
        mock_config.doAddBias = True
        mock_config.doAdd2DBias = True
        mock_config.doAddDark = False
        mock_config.doAddFlat = False
        mock_config.doAddFringe = False
        mock_config.doAddSky = False
        mock_config.doAddSource = False
        mock_config.doGenerateImage = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = IsrTaskLSSTConfig()
        isr_config.doBias = True
        isr_config.doDark = False
        isr_config.doDeferredCharge = False
        isr_config.doLinearize = False
        isr_config.doCorrectGains = False
        isr_config.doCrosstalk = False
        isr_config.doDefect = False
        isr_config.doBrighterFatter = False
        isr_config.doFlat = False

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(input_exp.clone(), bias=self.bias)

        # Rerun without doing the bias correction.
        isr_config.doBias = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(input_exp.clone())

        self.assertLess(
            np.mean(result.exposure.image.array),
            np.mean(result2.exposure.image.array),
        )
        self.assertLess(
            np.std(result.exposure.image.array),
            np.std(result2.exposure.image.array),
        )

        delta = result2.exposure.image.array - result.exposure.image.array
        self.assertFloatsAlmostEqual(delta, self.bias.image.array, atol=1e-5)

    def test_isrBootstrapDark(self):
        """Test processing of a ``bootstrap`` dark frame."""
        mock_config = isrMockLSST.IsrMockLSSTConfig()
        mock_config.isTrimmed = False
        mock_config.doAddBias = True
        mock_config.doAdd2DBias = True
        mock_config.doAddDark = True
        mock_config.doAddFlat = False
        mock_config.doAddFringe = False
        mock_config.doAddSky = False
        mock_config.doAddSource = False
        mock_config.doGenerateImage = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = IsrTaskLSSTConfig()
        isr_config.doBias = True
        isr_config.doDark = True
        isr_config.doDeferredCharge = False
        isr_config.doLinearize = False
        isr_config.doCorrectGains = False
        isr_config.doCrosstalk = False
        isr_config.doDefect = False
        isr_config.doBrighterFatter = False
        isr_config.doFlat = False

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(input_exp.clone(), bias=self.bias, dark=self.dark)

        # Rerun without doing the bias correction.
        isr_config.doDark = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(input_exp.clone(), bias=self.bias)

        self.assertLess(
            np.mean(result.exposure.image.array),
            np.mean(result2.exposure.image.array),
        )
        # The mock dark has no noise, so these should be equal.
        self.assertFloatsAlmostEqual(
            np.std(result.exposure.image.array),
            np.std(result2.exposure.image.array),
        )

        delta = result2.exposure.image.array - result.exposure.image.array
        exp_time = input_exp.getInfo().getVisitInfo().getExposureTime()
        self.assertFloatsAlmostEqual(delta, self.dark.image.array * exp_time, atol=1e-5)

    def test_isrBootstrapFlat(self):
        """Test processing of a ``bootstrap`` flat frame."""
        mock_config = isrMockLSST.IsrMockLSSTConfig()
        mock_config.isTrimmed = False
        mock_config.doAddBias = True
        mock_config.doAdd2DBias = True
        mock_config.doAddDark = True
        mock_config.doAddFlat = True
        mock_config.doAddFringe = False
        # The doAddSky option adds the equivalent of flat-field flux.
        mock_config.doAddSky = True
        mock_config.doAddSource = False
        mock_config.doGenerateImage = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = IsrTaskLSSTConfig()
        isr_config.doBias = True
        isr_config.doDark = True
        isr_config.doFlat = True
        isr_config.doDeferredCharge = False
        isr_config.doLinearize = False
        isr_config.doCorrectGains = False
        isr_config.doCrosstalk = False
        isr_config.doDefect = False
        isr_config.doBrighterFatter = False

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(input_exp.clone(), bias=self.bias, dark=self.dark, flat=self.flat)

        # Rerun without doing the bias correction.
        isr_config.doFlat = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(input_exp.clone(), bias=self.bias, dark=self.dark)

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

        ratio = result2.exposure.image.array / result.exposure.image.array
        self.assertFloatsAlmostEqual(ratio, self.flat.image.array, atol=1e-5)

    def test_isrBias(self):
        """Test processing of a bias frame."""
        mock_config = isrMockLSST.IsrMockLSSTConfig()
        mock_config.isTrimmed = False
        mock_config.doAddBias = True
        mock_config.doAdd2DBias = True
        mock_config.doAddDark = False
        mock_config.doAddFlat = False
        mock_config.doAddFringe = False
        mock_config.doAddSky = False
        mock_config.doAddSource = False
        mock_config.doGenerateImage = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = IsrTaskLSSTConfig()
        isr_config.doBias = True
        isr_config.doDark = False
        isr_config.doDeferredCharge = False
        isr_config.doLinearize = False
        isr_config.doCorrectGains = False
        isr_config.doCrosstalk = True
        isr_config.doDefect = False
        isr_config.doBrighterFatter = False
        isr_config.doFlat = False

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(input_exp.clone(), bias=self.bias, crosstalk=self.crosstalk)

        # Rerun without doing the bias correction.
        isr_config.doBias = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(input_exp.clone(), crosstalk=self.crosstalk)

        self.assertLess(
            np.mean(result.exposure.image.array),
            np.mean(result2.exposure.image.array),
        )
        self.assertLess(
            np.std(result.exposure.image.array),
            np.std(result2.exposure.image.array),
        )

        delta = result2.exposure.image.array - result.exposure.image.array
        self.assertFloatsAlmostEqual(delta, self.bias.image.array, atol=1e-5)

    def test_isrDark(self):
        """Test processing of a dark frame."""
        mock_config = isrMockLSST.IsrMockLSSTConfig()
        mock_config.isTrimmed = False
        mock_config.doAddBias = True
        mock_config.doAdd2DBias = True
        mock_config.doAddDark = True
        mock_config.doAddFlat = False
        mock_config.doAddFringe = False
        mock_config.doAddSky = False
        mock_config.doAddSource = False
        mock_config.doGenerateImage = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = IsrTaskLSSTConfig()
        isr_config.doBias = True
        isr_config.doDark = True
        isr_config.doDeferredCharge = False
        isr_config.doLinearize = False
        isr_config.doCorrectGains = False
        isr_config.doCrosstalk = True
        isr_config.doDefect = False
        isr_config.doBrighterFatter = False
        isr_config.doFlat = False

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(input_exp.clone(), bias=self.bias, dark=self.dark, crosstalk=self.crosstalk)

        # Rerun without doing the bias correction.
        isr_config.doDark = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(input_exp.clone(), bias=self.bias, crosstalk=self.crosstalk)

        self.assertLess(
            np.mean(result.exposure.image.array),
            np.mean(result2.exposure.image.array),
        )
        # The mock dark has no noise, so these should be equal.
        self.assertFloatsAlmostEqual(
            np.std(result.exposure.image.array),
            np.std(result2.exposure.image.array),
            atol=1e-6,
        )

        delta = result2.exposure.image.array - result.exposure.image.array
        exp_time = input_exp.getInfo().getVisitInfo().getExposureTime()
        self.assertFloatsAlmostEqual(delta, self.dark.image.array * exp_time, atol=1e-5)

    def test_isrFlat(self):
        """Test processing of a flat frame."""
        mock_config = isrMockLSST.IsrMockLSSTConfig()
        mock_config.isTrimmed = False
        mock_config.doAddBias = True
        mock_config.doAdd2DBias = True
        mock_config.doAddDark = True
        mock_config.doAddFlat = True
        mock_config.doAddFringe = False
        # The doAddSky option adds the equivalent of flat-field flux.
        mock_config.doAddSky = True
        mock_config.doAddSource = False
        mock_config.doGenerateImage = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = IsrTaskLSSTConfig()
        isr_config.doBias = True
        isr_config.doDark = True
        isr_config.doFlat = True
        # These should be set to true when we support them in tests.
        isr_config.doDeferredCharge = False
        isr_config.doLinearize = False
        isr_config.doCorrectGains = False
        isr_config.doCrosstalk = True
        isr_config.doDefect = False
        isr_config.doBrighterFatter = False

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(
            input_exp.clone(),
            bias=self.bias,
            dark=self.dark,
            flat=self.flat,
            crosstalk=self.crosstalk,
        )

        # Rerun without doing the bias correction.
        isr_config.doFlat = False
        isr_task2 = IsrTaskLSST(config=isr_config)
        result2 = isr_task2.run(input_exp.clone(), bias=self.bias, dark=self.dark, crosstalk=self.crosstalk)

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

        ratio = result2.exposure.image.array / result.exposure.image.array
        self.assertFloatsAlmostEqual(ratio, self.flat.image.array, atol=1e-5)

    def test_isrSkyImage(self):
        """Test processing of a sky image."""
        mock_config = isrMockLSST.IsrMockLSSTConfig()
        mock_config.isTrimmed = False
        mock_config.doAddBias = True
        mock_config.doAdd2DBias = True
        mock_config.doAddDark = True
        mock_config.doAddFlat = True
        mock_config.doAddCrosstalk = True
        # Set this to False until we have fringe correction.
        mock_config.doAddFringe = False
        mock_config.doAddSky = True
        mock_config.doAddSource = True
        mock_config.doGenerateImage = True

        mock = isrMockLSST.IsrMockLSST(config=mock_config)
        input_exp = mock.run()

        isr_config = IsrTaskLSSTConfig()
        isr_config.doBias = True
        isr_config.doDark = True
        isr_config.doFlat = True
        isr_config.doCrosstalk = True
        # This makes the one region look bad ...
        isr_config.doDefect = False

        # These should be set to true when we support them in tests.
        isr_config.doDeferredCharge = False
        isr_config.doLinearize = False
        isr_config.doCorrectGains = False
        isr_config.doCrosstalk = True
        isr_config.doBrighterFatter = False

        isr_task = IsrTaskLSST(config=isr_config)
        result = isr_task.run(
            input_exp.clone(),
            bias=self.bias,
            dark=self.dark,
            flat=self.flat,
            crosstalk=self.crosstalk,
            defects=self.defect,
            ptc=self.ptc,
        )

        clean_mock_config = isrMockLSST.IsrMockLSSTConfig()
        clean_mock_config.isTrimmed = True
        clean_mock_config.doAddBias = False
        clean_mock_config.doAdd2DBias = False
        clean_mock_config.doAddDark = False
        clean_mock_config.doAddDarkNoiseOnly = True
        clean_mock_config.doAddFlat = False
        clean_mock_config.doAddFringe = False
        clean_mock_config.doAddSky = True
        clean_mock_config.doAddSource = True
        clean_mock_config.doGenerateImage = True
        clean_mock_config.doRoundADU = False
        clean_mock_config.doApplyGain = False
        clean_mock_config.doAddCrosstalk = False

        clean_mock = isrMockLSST.IsrMockLSST(config=clean_mock_config)
        clean_exp = clean_mock.run()

        delta = result.exposure.image.array - clean_exp.image.array

        # TODO:
        # * Fix documented units and such.
        # * Add a defect bad column, consistently.
        # * Add a saturated set of pixels; check that they are masked
        #   and exclude from comparison.

        # There is noise from the overscan correction, given the
        # small overscan regions. This can/should be improved.
        self.assertLess(np.std(delta), 7.0)
        self.assertLess(np.max(np.abs(delta)), 7.0*4)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
