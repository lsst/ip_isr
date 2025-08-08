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

import os
import tempfile
import unittest

import lsst.utils.tests

from lsst.ip.isr import ShutterMotionProfile

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class ShutterMotionV1TestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        filename = os.path.join(TESTDIR, "data", "MC_O_20250526_000095_shutterMotionProfileOpen.json")
        self.calib = ShutterMotionProfile.readText(filename)

    def testRoundTrip(self):
        outPath = tempfile.mktemp() + '.yaml'
        self.calib.writeText(outPath)
        newCalib = ShutterMotionProfile.readText(outPath)
        self.assertEqual(self.calib, newCalib)

        outPath = tempfile.mktemp() + '.fits'
        self.calib.writeFits(outPath)
        newCalib = ShutterMotionProfile.readFits(outPath)
        self.assertEqual(self.calib, newCalib)

    def testMidpointCalculation(self):
        mid_accel, mid_position = self.calib.calculateMidpoint(modelName="hallSensorFit")
        self.assertFloatsAlmostEqual(mid_accel, 0.45052044442394185)
        self.assertFloatsAlmostEqual(mid_position, 0.44985220648315977)

        mid_accel, mid_position = self.calib.calculateMidpoint(modelName="motorEncoderFit")
        self.assertFloatsAlmostEqual(mid_accel, 0.45040944835186697)
        self.assertFloatsAlmostEqual(mid_position, 0.4497698788310154)


class ShutterMotionV2TestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        # We do not yet have true v2 shutter motion profiles, so this
        # is a v1 file modified to match the expected changes.
        filename = os.path.join(TESTDIR, "data", "MC_O_20250720_000554_shutterMotionProfileClose.json")
        self.calib = ShutterMotionProfile.readText(filename)

    def testRoundTrip(self):
        outPath = tempfile.mktemp() + '.yaml'
        self.calib.writeText(outPath)
        newCalib = ShutterMotionProfile.readText(outPath)
        self.assertEqual(self.calib, newCalib)

        outPath = tempfile.mktemp() + '.fits'
        self.calib.writeFits(outPath)
        newCalib = ShutterMotionProfile.readFits(outPath)
        self.assertEqual(self.calib, newCalib)

    def testMidpointCalculation(self):
        mid_accel, mid_position = self.calib.calculateMidpoint(modelName="hallSensorFit")
        self.assertFloatsAlmostEqual(mid_accel, 0.4503675064864023)
        self.assertFloatsAlmostEqual(mid_position, 0.45007949121581464)

        mid_accel, mid_position = self.calib.calculateMidpoint(modelName="motorEncoderFit")
        self.assertFloatsAlmostEqual(mid_accel, 0.4503689163377142)
        self.assertFloatsAlmostEqual(mid_position, 0.45008433680071347)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
