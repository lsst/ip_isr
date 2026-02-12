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

from lsst.afw.image import ExposureF
from lsst.ip.isr import ShutterMotionProfile, ShutterMotionProfileFull

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


class ShutterMotionFullTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        # Get the profiles from the text files.
        file_open = os.path.join(TESTDIR, "data",
                                 "MC_O_20251203_000317_shutterMotionProfileOpen.json")
        file_close = os.path.join(TESTDIR, "data",
                                  "MC_O_20251203_000317_shutterMotionProfileClose.json")
        self.profile_open = ShutterMotionProfile.readText(file_open)
        self.profile_close = ShutterMotionProfile.readText(file_close)

        # Simulate the image that would match these profiles, using
        # the as-written header of this exposure.
        self.exposure = ExposureF(1, 1)  # Image size doesn't matter.

        # Fill the exposure header with the important information
        self.exposure.metadata['SHUTTER OPEN STARTTIME TAI ISOT'] = '2025-12-04T08:46:56.592'
        self.exposure.metadata['SHUTTER OPEN STARTTIME TAI MJD'] = 61013.365932771936
        self.exposure.metadata['SHUTTER OPEN SIDE'] = 'MINUSX'
        self.exposure.metadata['SHUTTER OPEN MODEL'] = 'ThreeJerksModelv1'
        self.exposure.metadata['SHUTTER OPEN HALLSENSORFIT MODELSTARTTIME'] = 0.0005325114037250868
        self.exposure.metadata['SHUTTER OPEN HALLSENSORFIT PIVOTPOINT1'] = 0.2249632907921461
        self.exposure.metadata['SHUTTER OPEN HALLSENSORFIT PIVOTPOINT2'] = 0.6764472494529934
        self.exposure.metadata['SHUTTER OPEN HALLSENSORFIT JERK0'] = 32998.841790483784
        self.exposure.metadata['SHUTTER OPEN HALLSENSORFIT JERK1'] = -32961.436183162106
        self.exposure.metadata['SHUTTER OPEN HALLSENSORFIT JERK2'] = 33455.024219797444
        self.exposure.metadata['SHUTTER CLOSE STARTTIME TAI ISOT'] = '2025-12-04T08:47:06.593'
        self.exposure.metadata['SHUTTER CLOSE STARTTIME TAI MJD'] = 61013.366048524156
        self.exposure.metadata['SHUTTER CLOSE SIDE'] = 'PLUSX'
        self.exposure.metadata['SHUTTER CLOSE MODEL'] = 'ThreeJerksModelv1'
        self.exposure.metadata['SHUTTER CLOSE HALLSENSORFIT MODELSTARTTIME'] = 0.000900220338163209
        self.exposure.metadata['SHUTTER CLOSE HALLSENSORFIT PIVOTPOINT1'] = 0.2239073050036253
        self.exposure.metadata['SHUTTER CLOSE HALLSENSORFIT PIVOTPOINT2'] = 0.6804015164010819
        self.exposure.metadata['SHUTTER CLOSE HALLSENSORFIT JERK0'] = 33134.34774228922
        self.exposure.metadata['SHUTTER CLOSE HALLSENSORFIT JERK1'] = -32799.5059339621
        self.exposure.metadata['SHUTTER CLOSE HALLSENSORFIT JERK2'] = 36183.41782883381

    def testMidpointCalculations(self):
        # A "full" profile should yield the same results as the json
        # profiles.
        profile_full = ShutterMotionProfileFull.fromExposure(self.exposure)

        open_mid, close_mid = profile_full.calculateMidpoints()

        open_acc, _ = self.profile_open.calculateMidpoint()
        close_acc, _ = self.profile_close.calculateMidpoint()

        self.assertFloatsAlmostEqual(open_mid, open_acc)
        self.assertFloatsAlmostEqual(close_mid, close_acc)

    def testEdgeCases(self):
        # These are options that will not work with the exposure based
        # profiles, as they only have the Hall sensor fit, and no
        # position data.
        profile_full = ShutterMotionProfileFull.fromExposure(self.exposure)
        with self.assertRaises(KeyError):
            profile_full.profile_open.calculateMidpoint(skipPosition=False)

        with self.assertRaises(KeyError):
            profile_full.profile_close.calculateMidpoint(skipPosition=False)

        with self.assertRaises(RuntimeError):
            profile_full.profile_open.calculateMidpoint(modelName="motorEncoderFit", skipPosition=True)

        with self.assertRaises(RuntimeError):
            profile_full.profile_close.calculateMidpoint(modelName="motorEncoderFit", skipPosition=True)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
