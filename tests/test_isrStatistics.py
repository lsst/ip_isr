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

import lsst.utils.tests

from lsst.ip.isr.isrMock import RawMock
from lsst.ip.isr import IsrTask


class IsrStatisticsTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        config = IsrTask.ConfigClass()
        config.doAssembleCcd = True
        config.doOverscan = True
        config.overscan.fitType = 'MEDIAN_PER_ROW'
        config.doBias = False
        config.doDark = False
        config.doFlat = False
        config.doLinearize = False
        config.doCalculateStatistics = True
        config.doFringe = False
        config.isrStats.doBandingStatistics = True

        self.task = IsrTask(config=config)
        self.exposure = RawMock().run()

        self.overscans = []
        for amp in self.exposure.getDetector().getAmplifiers():
            overscanResults = self.task.overscanCorrection(self.exposure, amp)
            self.overscans.append(overscanResults)

    def test_bandingStatistics(self):
        """Look at the banding on a raw mock image.

        The value obtained is far larger than we expect in real data.
        """
        bandingResults = self.task.isrStats.measureBanding(self.exposure, self.overscans)

        self.assertFloatsAlmostEqual(bandingResults['DET_BANDING'], 153.6335, atol=1e2)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
