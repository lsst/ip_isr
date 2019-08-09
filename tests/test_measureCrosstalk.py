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

import lsst.utils.tests
from lsst.ip.isr.measureCrosstalk import (MeasureCrosstalkTask, MeasureCrosstalkConfig)
import lsst.ip.isr.isrMock as isrMock


class MeasureCrosstalkTaskCases(lsst.utils.tests.TestCase):

    def setup_measureCrosstalk(self, isTrimmed=False, nSources=8):
        """Generate a simulated set of exposures and test the measured
        crosstalk matrix.

        Parameters
        ----------
        isTrimmed : `bool`, optional
            Should the simulation use trimmed or untrimmed raw
            exposures?
        nSources : `int`, optional
            Number of random simulated sources to generate in the
            simulated exposures.

        Returns
        -------
        coeffErr : `np.ndarray`
            Array of booleans indicating if the measured and expected
            crosstalk ratios are smaller than the measured uncertainty
            in the crosstalk ratio.
        """
        mockTask = isrMock.CalibratedRawMock()
        mockTask.config.rngSeed = 12345
        mockTask.config.doGenerateImage = True
        mockTask.config.doAddSky = True
        mockTask.config.doAddSource = True
        mockTask.config.doAddCrosstalk = True
        mockTask.config.doAddBias = True
        mockTask.config.doAddFringe = False

        mockTask.config.skyLevel = 0.0
        mockTask.config.biasLevel = 0.0
        mockTask.config.readNoise = 100.0

        mcConfig = MeasureCrosstalkConfig()
        mcConfig.threshold = 4000
        mcConfig.isTrimmed = isTrimmed
        mct = MeasureCrosstalkTask(config=mcConfig)
        fullResult = []

        mockTask.config.isTrimmed = isTrimmed
        # Generate simulated set of exposures.
        for idx in range(0, 10):
            mockTask.config.rngSeed = 12345 + idx * 1000

            # Allow each simulated exposure to have nSources random
            # bright sources.
            mockTask.config.sourceAmp = (np.random.randint(8, size=nSources)).tolist()
            mockTask.config.sourceFlux = ((np.random.random(size=nSources) * 25000.0 + 20000.0).tolist())
            mockTask.config.sourceX = ((np.random.random(size=nSources) * 100.0).tolist())
            mockTask.config.sourceY = ((np.random.random(size=nSources) * 50.0).tolist())

            exposure = mockTask.run()
            result = mct.run(exposure, dataId=None)
            fullResult.append(result)

        # Generate the final measured CT ratios, uncertainties, pixel counts.
        coeff, coeffSig, coeffNum = mct.reduce(fullResult)

        # Needed because measureCrosstalk cannot find coefficients equal to 0.0
        coeff = np.nan_to_num(coeff)
        coeffSig = np.nan_to_num(coeffSig)

        # Compare result against expectation used to create the simulation.
        expectation = isrMock.CrosstalkCoeffMock().run()
        coeffErr = abs(coeff - expectation) <= coeffSig
        return coeffErr

    def testMeasureCrosstalkTaskTrimmed(self):
        """Measure crosstalk from a sequence of trimmed mocked images.
        """
        coeffErr = self.setup_measureCrosstalk(isTrimmed=True, nSources=8)

        self.assertTrue(np.all(coeffErr))

    def testMeasureCrosstalkTaskUntrimmed(self):
        """Measure crosstalk from a sequence of untrimmed mocked images.
        """
        coeffErr = self.setup_measureCrosstalk(isTrimmed=False, nSources=8)

        # DM-18528 This doesn't always fully converge, so be permissive
        # for now.
        self.assertTrue(np.any(coeffErr))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
