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

# The following values are used to test the AmpOffsetTask.
BACKGROUND_VALUE = 100
RAMP_XSCALE = 0.0
RAMP_YSCALE = 100.0


class AmpOffsetTest(lsst.utils.tests.TestCase):
    def setUp(self):
        # Testing with a single detector that has 8 amplifiers in a 4x2
        # configuration to ensure functionality in a general 2-dimensional
        # scenario. Each amplifier measures 100x51 in dimensions.
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
        self.measuredPedestalsConstantBackground = {
            "unweighted": {
                "symmetric": [
                    -1.36344239,
                    -0.997554,
                    -0.4337341,
                    -0.06784741,
                    0.06784741,
                    0.43373411,
                    0.99755399,
                    1.36344238,
                ],
                "random": [
                    1.2545466,
                    -1.42373527,
                    -0.32187415,
                    -0.61296588,
                    0.44789211,
                    0.21339304,
                    -0.6362251,
                    1.07896865,
                ],
                "artificial": [
                    -3.96260237,
                    -6.31372346,
                    -5.54533972,
                    -9.8481738,
                    8.68735117,
                    5.97670723,
                    3.85524471,
                    7.15053624,
                ],
            },
            "weighted": {
                "symmetric": [
                    -1.36344235,
                    -0.99755403,
                    -0.43373414,
                    -0.06784737,
                    0.06784738,
                    0.43373415,
                    0.99755403,
                    1.36344234,
                ],
                "random": [
                    1.28420419,
                    -1.42457154,
                    -0.34101654,
                    -0.62264482,
                    0.41823451,
                    0.21422931,
                    -0.61708271,
                    1.0886476,
                ],
                "artificial": [
                    -3.96861485,
                    -6.30713176,
                    -5.5651146,
                    -9.82897816,
                    8.69336365,
                    5.97011553,
                    3.87501958,
                    7.13134059,
                ],
            },
        }
        self.measuredPedestalsRampBackground = {
            "unweighted": {
                "symmetric": [
                    -1.35330781,
                    -0.92859869,
                    -0.50268509,
                    -0.07798556,
                    0.0779834,
                    0.50268579,
                    0.92859927,
                    1.35330869,
                ],
                "random": [
                    1.31054124,
                    -1.64586418,
                    -0.15561936,
                    -0.54440476,
                    0.30516155,
                    0.25678385,
                    -0.78329776,
                    1.25669942,
                ],
                "artificial": [
                    -3.92486412,
                    -6.64218097,
                    -4.85473319,
                    -9.58046338,
                    8.67801534,
                    5.66513899,
                    3.59097321,
                    7.06811412,
                ],
            },
            "weighted": {
                "symmetric": [
                    -1.35330818,
                    -0.92859835,
                    -0.5026847,
                    -0.07798592,
                    0.07798378,
                    0.50268544,
                    0.92859888,
                    1.35330905,
                ],
                "random": [
                    1.31103334,
                    -1.63924511,
                    -0.16414238,
                    -0.54299292,
                    0.30466945,
                    0.25016478,
                    -0.77477473,
                    1.25528757,
                ],
                "artificial": [
                    -3.93099297,
                    -6.63242338,
                    -4.86026077,
                    -9.57856454,
                    8.68414419,
                    5.6553814,
                    3.5965008,
                    7.06621528,
                ],
            },
        }
        self.measuredSigma = {
            "unweighted": {
                "symmetric": 3.550909965491515e-08,
                "random": 9.271480035637249e-08,
                "artificial": 3.9430895759341275e-14,
            },
            "weighted": {
                "symmetric": 3.6093886632617167e-08,
                "random": 9.280982569622311e-08,
                "artificial": 8.690997162288794e-15,
            },
        }
        self.measuredSigmaConstantBackground = {
            "unweighted": {
                "symmetric": 0.06679099233901276,
                "random": 0.12836424357658588,
                "artificial": 0.5766923710558869,
            },
            "weighted": {
                "symmetric": 0.06679099233983563,
                "random": 0.12854047710163782,
                "artificial": 0.579062202931159,
            },
        }
        self.measuredSigmaRampBackground = {
            "unweighted": {
                "symmetric": 0.8049174282700369,
                "random": 0.8128787501698868,
                "artificial": 6.090700238137678,
            },
            "weighted": {
                "symmetric": 0.8049174282702547,
                "random": 0.8089364786360683,
                "artificial": 6.089900455983161,
            },
        }

    def tearDown(self):
        del self.mock

    def buildExposure(self, valueType, addBackground=False, rampBackground=False):
        """
        Build and return an exposure with different types of value
        distributions across its amplifiers.

        Parameters
        ----------
        valueType : `str`
            Determines the distribution type of values across the amplifiers.
            - "symmetric": Creates a symmetric constant interval distribution
            of values.
            - "random": Generates a random distribution of values.
            - "artificial": Uses a predefined array of values to simulate a
            weight-sensitive condition for the output pedestals. This set of
            values designed to show more change across the short interface and
            will not solve exactly.

        addBackground : `bool`, optional
            If True, adds a background value to the entire exposure.

        rampBackground : `bool`, optional
            Whether the added background should be a ramp.

        Returns
        -------
        exp : `~lsst.afw.image.Exposure`
            An exposure object modified according to the specified valueType
            and background addition.

        Notes
        -----
        This method is used to generate different scenarios of exposure data
        for testing and analysis. The 'artificial' valueType is particularly
        useful for testing the robustness of algorithms under non-ideal or
        challenging data conditions.
        """
        exp = self.mock.getExposure()
        detector = exp.getDetector()
        amps = detector.getAmplifiers()
        values = {
            "symmetric": np.linspace(-2.5, 2.5, len(amps)),
            "random": np.random.RandomState(seed=1746).uniform(-2.5, 2.5, len(amps)),
            "artificial": np.array([5, 1, 3, -4, 30, 25, 22, 27]),
        }
        self.values = values[valueType]
        for amp, value in zip(amps, self.values):
            exp.image[amp.getBBox()] = value
        if addBackground:
            exp.image.array += BACKGROUND_VALUE
            if rampBackground:
                # Add a gradient.
                self.amplifierAddYGradient(exp.image, 0.0, RAMP_YSCALE)
                # Add another gradient to the other direction.
                self.amplifierAddXGradient(exp.image, 0.0, RAMP_XSCALE)
        else:
            assert not rampBackground, "rampBackground requires addBackground=True"
        return exp

    def runAmpOffsetWithBackground(self, valueType, rampBackground=False):
        """
        Tests the AmpOffsetTask on an exposure with a background added.

        Parameters
        ----------
        valueType : `str`
            Determines the distribution type of values across the amplifiers.
            See `buildExposure` for details.

        rampBackground : `bool`
            Whether the added background should be a ramp.
        """
        if rampBackground:
            measuredPedestals = self.measuredPedestalsRampBackground
            measuredSigma = self.measuredSigmaRampBackground
        else:
            measuredPedestals = self.measuredPedestalsConstantBackground
            measuredSigma = self.measuredSigmaConstantBackground

        for applyWeights in [False, True]:
            exp = self.buildExposure(valueType, addBackground=True, rampBackground=rampBackground)
            amps = exp.getDetector().getAmplifiers()
            config = AmpOffsetConfig()
            config.doBackground = True
            config.doDetection = True
            config.ampEdgeWidth = 12
            config.applyWeights = applyWeights
            config.doApplyAmpOffset = True  # Updates the exposure in place.
            if valueType == "random":
                # For this specific case, the fraction of unmasked pixels for
                # amp interface 01 is unusually small.
                config.ampEdgeMinFrac = 0.1
            if valueType == "artificial":
                # For this extreme case, we expect the interface offsets to be
                # unusually large.
                config.ampEdgeMaxOffset = 50
            task = AmpOffsetTask(config=config)
            pedestals = task.run(exp).pedestals
            nAmps = len(amps)
            if valueType == "symmetric":
                for i in range(nAmps // 2):
                    self.assertAlmostEqual(pedestals[i], -pedestals[nAmps - i - 1], 5)

            ampBBoxes = [amp.getBBox() for amp in amps]
            maskedImage = exp.getMaskedImage()
            nX = exp.getWidth() // (task.shortAmpSide * config.backgroundFractionSample) + 1
            nY = exp.getHeight() // (task.shortAmpSide * config.backgroundFractionSample) + 1
            bg = task.background.fitBackground(maskedImage, nx=int(nX), ny=int(nY))
            bgSubtractedValues = []
            for i, bbox in enumerate(ampBBoxes):
                ampBgImage = bg.getImageF(
                    interpStyle=task.background.config.algorithm,
                    undersampleStyle=task.background.config.undersampleStyle,
                    bbox=bbox,
                )
                if not rampBackground:
                    bgSubtractedValues.append(self.values[i] + BACKGROUND_VALUE - np.mean(ampBgImage.array))
                else:
                    # With the added gradient, averaging is required for a
                    # proper approximation.
                    meanValueWithBg = exp.image[bbox].array.mean()
                    bgSubtractedValues.append(meanValueWithBg - np.mean(ampBgImage.array))

            approximatePedestals = np.array(bgSubtractedValues) - np.mean(bgSubtractedValues)

            weightType = "weighted" if applyWeights else "unweighted"
            for pedestal, value in zip(pedestals, measuredPedestals[weightType][valueType]):
                self.assertAlmostEqual(pedestal, value, 5)
            # If we are getting it wrong, let's not get it wrong by more than
            # some specified DN.
            self.assertAlmostEqual(
                np.std(pedestals - approximatePedestals),
                measuredSigma[weightType][valueType],
                4 if rampBackground else 12,
            )
            if valueType == "artificial":
                if not applyWeights:
                    sigmaUnweighted = np.std(pedestals - approximatePedestals)
                else:
                    # Verify that the weighted sigma differs from the
                    # unweighted sigma. It's not anticipated for the weighted
                    # sigma to be consistently smaller, given that the
                    # estimated background isn't uniform and our expected
                    # pedestals are approximations.
                    sigmaWeighted = np.std(pedestals - approximatePedestals)
                    self.assertNotEqual(sigmaWeighted, sigmaUnweighted)

    def testAmpOffsetEffectOnExposure(self):
        exp0 = self.buildExposure("random", addBackground=True, rampBackground=True)
        exp = exp0.clone()
        config = AmpOffsetConfig()
        config.doBackground = True
        config.doDetection = True
        config.ampEdgeWidth = 12
        config.applyWeights = True

        # Configure to not apply amp offset to the exposure and run the task.
        # Verify that the exposure remains unchanged.
        config.doApplyAmpOffset = False
        AmpOffsetTask(config=config).run(exp)
        self.assertFloatsEqual(exp0.image.array, exp.image.array)

        # Configure to apply amp offset to the exposure and run the task.
        # Verify that the exposure is updated.
        config.doApplyAmpOffset = True
        AmpOffsetTask(config=config).run(exp)
        self.assertFloatsNotEqual(exp0.image.array, exp.image.array)

    @methodParameters(valueType=["symmetric", "random", "artificial"])
    def testAmpOffset(self, valueType):
        for applyWeights in [False, True]:
            exp = self.buildExposure(valueType, addBackground=False)
            config = AmpOffsetConfig()
            config.doBackground = False
            config.doDetection = False
            config.ampEdgeWidth = 12  # Given 100x51 amps in our mock detector.
            config.doApplyAmpOffset = True  # Updates the exposure in place.
            if valueType == "artificial":
                # For this extreme case, we expect the interface offsets to be
                # unusually large.
                config.ampEdgeMaxOffset = 50
            config.applyWeights = applyWeights
            task = AmpOffsetTask(config=config)
            pedestals = task.run(exp).pedestals
            if valueType == "symmetric":
                self.assertAlmostEqual(np.sum(exp.image.array), 0, 6)
            truePedestals = self.values - np.mean(self.values)
            for pedestal, value in zip(pedestals, truePedestals):
                self.assertAlmostEqual(pedestal, value, 6)
            weightType = "weighted" if applyWeights else "unweighted"
            self.assertAlmostEqual(
                np.std(pedestals - truePedestals), self.measuredSigma[weightType][valueType], 12
            )
            if valueType == "artificial":
                if not applyWeights:
                    sigmaUnweighted = np.std(pedestals - truePedestals)
                else:
                    # Verify that the weighted sigma differs from the
                    # unweighted sigma. It's not anticipated for the weighted
                    # sigma to be consistently smaller, given the numerical
                    # noise from exceedingly small value calculations and the
                    # variations in library versions and operating systems.
                    sigmaWeighted = np.std(pedestals - truePedestals)
                    self.assertNotEqual(sigmaWeighted, sigmaUnweighted)

    @methodParameters(valueType=["symmetric", "random", "artificial"])
    def testAmpOffsetWithConstantBackground(self, valueType):
        self.runAmpOffsetWithBackground(valueType, rampBackground=False)

    @methodParameters(valueType=["symmetric", "random", "artificial"])
    def testAmpOffsetWithRampBackground(self, valueType):
        self.runAmpOffsetWithBackground(valueType, rampBackground=True)

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
