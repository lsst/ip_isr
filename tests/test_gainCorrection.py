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
import tempfile

import numpy as np

import lsst.utils.tests

from lsst.ip.isr import GainCorrection


class GainCorrectionTest(lsst.utils.tests.TestCase):
    """Test the GainCorrection dataset."""
    def setUp(self):
        self.gainCorrection = GainCorrection(
            ampNames=["C01", "C02"],
            gainAdjustments=[1.0, 1.1],
        )

    def _checkEqual(self, a, b):
        self.assertEqual(b.metadata, a.metadata)
        np.testing.assert_array_equal(b.ampNames, a.ampNames)
        self.assertEqual(b.gainAdjustments.dtype, a.gainAdjustments.dtype)
        np.testing.assert_array_almost_equal(b.gainAdjustments, a.gainAdjustments)

    def testRoundTrip(self):
        """Test persistence round-tripping."""

        with tempfile.NamedTemporaryFile(suffix=".yaml") as f:
            usedFilename = self.gainCorrection.writeText(f.name)
            fromText = GainCorrection.readText(usedFilename)
        self._checkEqual(fromText, self.gainCorrection)

        with tempfile.NamedTemporaryFile(suffix=".fits") as f:
            usedFilename = self.gainCorrection.writeFits(f.name)
            fromFits = GainCorrection.readFits(usedFilename)
        self._checkEqual(fromFits, self.gainCorrection)

    def testCorrectGains(self):
        gains = {
            ampName: 1.5
            for ampName in self.gainCorrection.ampNames
        }

        self.gainCorrection.correctGains(gains)

        for i, ampName in enumerate(self.gainCorrection.ampNames):
            np.testing.assert_almost_equal(1.5*self.gainCorrection.gainAdjustments[i], gains[ampName])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
