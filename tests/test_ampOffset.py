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
import unittest
import numpy as np

import lsst.utils.tests
from lsst.utils.tests import methodParameters
from lsst.ip.isr.ampOffset import AmpOffsetConfig, AmpOffsetTask
from lsst.ip.isr.isrMock import IsrMock


class AmpOffsetTest(lsst.utils.tests.TestCase):
    def setUp(self):
        # Test with a single detector featuring 4x2 amplifiers to ensure
        # functionality in a general 2-dimensional scenario.
        config = IsrMock.ConfigClass()
        config.isLsstLike = True
        config.doAddBias = False
        config.doAddDark = False
        config.doAddFlat = False
        config.doAddFringe = False
        config.doGenerateImage = True
        config.doGenerateData = True
        config.doGenerateAmpDict = True
        self.mock = IsrMock(config=config)
        self.measuredValuesBackground = [
            0.13336708,
            0.06299729,
            0.34589145,
            0.27551733,
            -0.27551731,
            -0.34589142,
            -0.06299731,
            -0.1333671,
        ]
        self.measuredValuesRampBackground = [
            0.19094621,
            0.21447723,
            0.19441773,
            0.21793666,
            -0.21794104,
            -0.1944167,
            -0.21447818,
            -0.19094191,
        ]
        self.measuredSigmaBackground = 1.7861559064582404
        self.measuredSigmaRampBackground = 1.8148713749702268

    def tearDown(self):
        del self.mock

    def buildExposure(self, valueType="constant_interval", addBackground=False):
        # constant_interval - like what we had before
        # random
        # artificially make it worse so that the weightig can show its advantage -- weight-sensitive
        # E F G H
        # A B C D
        exp = self.mock.getExposure()
        detector = exp.getDetector()
        amps = detector.getAmplifiers()
        values = {
            "constant_interval": np.linspace(-2.5, 2.5, len(amps)),
            "random": np.random.RandomState(seed=1746).uniform(-2.5, 2.5, len(amps)),
            "artificial": np.array([5, 1, 3, -4, 30, 25, 22, 27])
        }
        self.values = values[valueType]
        for amp, value in zip(amps, self.values):
            exp.image[amp.getBBox()] = value
        if addBackground:
            exp.image.array += 100
        return exp

    @methodParameters(valueType=["constant_interval", "random", "artificial"])
    def testAmpOffset(self, valueType):
        exp = self.buildExposure(addBackground=False, valueType=valueType)
        config = AmpOffsetConfig()
        config.doBackground = False
        config.doDetection = False
        config.ampEdgeWidth = 12  # Given 100x51 amps in our mock detector.
        if valueType=="artificial":
            config.ampEdgeMaxOffset = 50
        config.applyWeights = True
        task = AmpOffsetTask(config=config)
        pedestals = task.run(exp).pedestals
        # self.assertEqual(np.sum(exp.image.array), 0)
        # for pedestal, value in zip(pedestals, self.values):
        #     self.assertAlmostEqual(pedestal, value, 6)
        print(f"applyWeights = {config.applyWeights}")
        print(f"valueType: {valueType} -> self.values = {self.values}")
        print(f"pedestals = {pedestals}")
        print(f"np.std(pedestals) = {np.std(pedestals)}")
        print(f"np.std(pedestals-self.values) = {np.std(pedestals-self.values)}")
        print("-"*80)

    # def testAmpOffsetWithBackground(self):
    #     exp = self.buildExposure(addBackground=True)
    #     amps = exp.getDetector().getAmplifiers()
    #     config = AmpOffsetConfig()
    #     config.doBackground = True
    #     config.doDetection = True
    #     config.ampEdgeWidth = 12
    #     task = AmpOffsetTask(config=config)
    #     pedestals = task.run(exp).pedestals
    #     nAmps = len(amps)
    #     for i in range(nAmps // 2):
    #         self.assertAlmostEqual(pedestals[i], -pedestals[nAmps - i - 1], 5)
    #     breakpoint()
    #     for pedestal, value in zip(pedestals, self.measuredValuesBackground):
    #         self.assertAlmostEqual(pedestal, value, 5)
    #     # If we are getting it wrong, let's not get it wrong by more than some
    #     # specified DN.
    #     self.assertAlmostEqual(np.std(pedestals - self.values), self.measuredSigmaBackground, 5)

    # def testAmpOffsetWithRampBackground(self):
    #     exp = self.buildExposure(addBackground=True)
    #     amps = exp.getDetector().getAmplifiers()
    #     yscale = 100.0
    #     xscale = 0.0
    #     # Add a gradient.
    #     self.amplifierAddYGradient(exp.image, 0.0, yscale)
    #     # Add another gradient to the other direction.
    #     self.amplifierAddXGradient(exp.image, 0.0, xscale)
    #     config = AmpOffsetConfig()
    #     config.doBackground = True
    #     config.doDetection = True
    #     config.ampEdgeWidth = 12
    #     task = AmpOffsetTask(config=config)
    #     pedestals = task.run(exp).pedestals
    #     nAmps = len(amps)
    #     for i in range(nAmps // 2):
    #         self.assertAlmostEqual(pedestals[i], -pedestals[nAmps - i - 1], 5)
    #     breakpoint()
    #     for pedestal, value in zip(pedestals, self.measuredValuesRampBackground):
    #         self.assertAlmostEqual(pedestal, value, 5)
    #     self.assertAlmostEqual(np.std(pedestals - self.values), self.measuredSigmaRampBackground)

    # The two static methods below are taken from ip_isr/isrMock.
    @staticmethod
    def amplifierAddYGradient(ampData, start, end):
        nPixY = ampData.getDimensions().getY()
        ampArr = ampData.array
        ampArr[:] = ampArr[:] + (
            np.interp(range(nPixY), (0, nPixY - 1), (start, end)).reshape(nPixY, 1)
            + np.zeros(ampData.getDimensions()).transpose()
        )

    @staticmethod
    def amplifierAddXGradient(ampData, start, end):
        nPixX = ampData.getDimensions().getX()
        ampArr = ampData.array
        ampArr[:] = ampArr[:] + (
            np.interp(range(nPixX), (0, nPixX - 1), (start, end)).reshape(1, nPixX)
            + np.zeros(ampData.getDimensions()).transpose()
        )


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
