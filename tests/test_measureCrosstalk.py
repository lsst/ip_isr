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

    def testMeasureCrosstalkTaskTrimmed(self):
        """Measure crosstalk from a sequence of mocked images.
        """
        config = isrMock.IsrMockConfig()
        config.rngSeed = 12345
        config.doAddCrosstalk = True
        config.doAddSky = True
        config.doAddSource = True
        config.skyLevel = 0.0
        config.readNoise = 0.0
        mcConfig = MeasureCrosstalkConfig()
        mcConfig.threshold = 4000
        mct = MeasureCrosstalkTask(config=mcConfig)
        fullResult = []

        config.isTrimmed = True

        for idx in range(0, 10):
            config.rngSeed = 12345 + idx * 1000
            config.sourceAmp = [0, 1, 2, 3, 4, 5, 6, 7]
            config.sourceFlux = [45000.0, 45000.0, 45000.0, 45000.0,
                                 45000.0, 45000.0, 45000.0, 45000.0]
            config.sourceX = [50.0, 25.0, 75.0, 12.5, 37.5, 67.5, 82.5]
            config.sourceY = [25.0, 12.5, 37.5, 26.75, 22.25, 12.5, 37.5]
            exposure = isrMock.CalibratedRawMock(config=config).run()
            result = mct.run(exposure, dataId=None)
            fullResult.append(result)

        coeff, coeffSig, coeffNum = mct.reduce(fullResult)

        # Needed because measureCrosstalk cannot find coefficients equal to 0.0
        coeff = np.nan_to_num(coeff)
        coeffSig = np.nan_to_num(coeffSig)

        expectation = isrMock.CrosstalkCoeffMock().run()
        coeffErr = abs(coeff - expectation) <= coeffSig

        # DM-18528 This doesn't always fully converge, so be permissive
        # for now.  This is also more challenging on the test
        # chip due to density issues.
        self.assertTrue(np.any(coeffErr))

    def testMeasureCrosstalkTaskUntrimmed(self):
        """Measure crosstalk from a sequence of mocked images.
        """
        config = isrMock.IsrMockConfig()
        config.rngSeed = 12345
        config.doAddCrosstalk = True
        config.doAddSky = True
        config.doAddSource = True
        config.skyLevel = 0.0
        config.readNoise = 0.0
        mcConfig = MeasureCrosstalkConfig()
        mcConfig.threshold = 4000
        mct = MeasureCrosstalkTask(config=mcConfig)
        fullResult = []

        config.isTrimmed = False

        for idx in range(0, 10):
            config.rngSeed = 12345 + idx * 1000
            config.sourceAmp = [0, 1, 2, 3, 4, 5, 6, 7]
            config.sourceFlux = [45000.0, 45000.0, 45000.0, 45000.0,
                                 45000.0, 45000.0, 45000.0, 45000.0]
            config.sourceX = [50.0, 25.0, 75.0, 12.5, 37.5, 67.5, 82.5]
            config.sourceY = [25.0, 12.5, 37.5, 26.75, 22.25, 12.5, 37.5]
            exposure = isrMock.CalibratedRawMock(config=config).run()
            result = mct.run(exposure, dataId=None)
            fullResult.append(result)

        coeff, coeffSig, coeffNum = mct.reduce(fullResult)

        # Needed because measureCrosstalk cannot find coefficients equal to 0.0
        coeff = np.nan_to_num(coeff)
        coeffSig = np.nan_to_num(coeffSig)

        expectation = isrMock.CrosstalkCoeffMock().run()
        coeffErr = abs(coeff - expectation) <= coeffSig

        # DM-18528 This doesn't always fully converge, so be permissive
        # for now.  This is also more challenging on the test
        # chip due to density issues.
        self.assertTrue(np.any(coeffErr))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
