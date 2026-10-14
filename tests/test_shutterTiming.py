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
"""Tests of `lsst.ip.isr.shutterTiming` against the reference results in
tests/data/shutterTiming (see its README.md).
"""


import glob
import json
import math
import os
import unittest

import numpy as np

import lsst.geom
import lsst.utils.tests
from lsst.afw.cameraGeom import Orientation
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
from lsst.daf.base import PropertyList
from lsst.ip.isr.shutterTiming import (
    ShutterTimingConfig,
    ShutterTimingFlag,
    ShutterTimingStatus,
    computeShutterTiming,
)

TESTDIR = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(TESTDIR, "data", "shutterTiming")
BASE = "MC_O_20260712_000100"

#: Agreement target (s): centre, focal-plane and per-source times, T_eff.
TOL_S = 10e-6
SECONDS_PER_DAY = 86400.0
#: The reference's flag bits that this implementation keeps.
KEPT_FLAGS = int(ShutterTimingFlag.NO_PROFILE | ShutterTimingFlag.PARAM_RANGE
                 | ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN)


def readFixture(name):
    with open(os.path.join(DATADIR, name + ".json")) as f:
        return json.load(f)


def makeDetector(detectorId, geometry):
    """Build an afw detector with 10 um pixels and no rotation, as every
    LSSTCam detector.

    Parameters
    ----------
    detectorId : `int`
        The detector ID.
    geometry : `dict`
        ``center_mm`` (focal-plane position of the centre), ``center_pix``
        (its pixel coordinates) and ``size`` (pixels), as in the fixtures.

    Returns
    -------
    detector : `lsst.afw.cameraGeom.Detector`
        The detector.
    """
    orientation = Orientation(lsst.geom.Point2D(*geometry["center_mm"]),
                              lsst.geom.Point2D(*geometry["center_pix"]),
                              lsst.geom.Angle(0.0, lsst.geom.degrees))
    bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(*geometry["size"]))
    return DetectorWrapper(id=detectorId, bbox=bbox, pixelSize=(0.01, 0.01),
                           orientation=orientation).detector


class ShutterTimingTestCase(lsst.utils.tests.TestCase):
    """Agreement with the reference results, and the cases that occur in
    processing.
    """

    @classmethod
    def setUpClass(cls):
        cls.config = ShutterTimingConfig()
        with open(os.path.join(DATADIR, "mutations.json")) as f:
            cls.mutations = json.load(f)["mutations"]
        doc = readFixture(BASE)
        cls.base = doc["metadata"]
        cls.detectors = {int(det): makeDetector(int(det), expected["geometry"])
                         for det, expected in doc["detectors"].items()}

    def compute(self, md, det=94):
        return computeShutterTiming(md, self.detectors[det], md["EXPTIME"], self.config)

    def mutated(self, name):
        """Return the base metadata with a mutation's cards set."""
        md = dict(self.base)
        md.update(self.mutations[name]["set"])
        return md

    def checkAgainstReference(self, doc, metadata, label):
        """Check all detectors and samples of a fixture against the
        reference; its flags are compared on the kept bits.
        """
        # Every detector: status, flags and its times.
        timings = {}
        for det, expected in doc["detectors"].items():
            t = timings[int(det)] = computeShutterTiming(
                metadata, makeDetector(int(det), expected["geometry"]), metadata["EXPTIME"], self.config
            )
            label2 = f"{label} {det}"

            self.assertEqual(t.status, ShutterTimingStatus.OK, label2)
            self.assertEqual(t.message, "", label2)
            self.assertEqual(int(t.flags), expected["qc_flags"] & KEPT_FLAGS, label2)

            diffs = dict(
                center=abs(t.centerMidpointMjdTai - expected["center_mjd_tai"]) * SECONDS_PER_DAY,
                teff=abs(t.effectiveExposureTime - expected["effective_exposure_time_s"]),
                focalPlane=abs(t.focalPlaneMidpointMjdTai - doc["visit_mjd_tai"]) * SECONDS_PER_DAY,
            )
            for k, v in diffs.items():
                self.assertLessEqual(v, TOL_S, f"{label2}: {k}")

        # Per-source times at the sample positions, on the detector and up
        # to 150 pixels off it.
        for det, x, y, expected in doc["samples"]:
            got = float(timings[det].midpointMjdTai(np.array(x), np.array(y)))
            self.assertLessEqual(abs(got - expected)*SECONDS_PER_DAY, TOL_S, f"{label} {det} ({x}, {y})")

    def testFixtures(self):
        """Check real exposures, one per blade direction, against the
        reference on detectors from the centre to the corners of the focal
        plane.
        """
        names = sorted(os.path.basename(p)[:-5] for p in glob.glob(os.path.join(DATADIR, "MC_O_*.json")))
        self.assertEqual(len(names), 2)

        for name in names:
            doc = readFixture(name)
            self.assertEqual(len(doc["detectors"]), 7)
            self.checkAgainstReference(doc, doc["metadata"], name)

    def testNoCards(self):
        """Check that exposures before the shutter cards are UNAVAILABLE,
        with no times.
        """
        md = {k: v for k, v in self.base.items() if not k.startswith("SHUTTER")}
        xy = np.array([0.0, 2000.0])

        for det in (0, 94, 188):
            t = self.compute(md, det)
            self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE)
            self.assertEqual(t.flags, ShutterTimingFlag.NO_PROFILE)
            self.assertTrue(math.isnan(t.centerMidpointMjdTai))
            self.assertTrue(math.isnan(t.focalPlaneMidpointMjdTai))
            self.assertTrue(np.all(np.isnan(t.midpointMjdTai(xy, xy))))

    def testReanchoring(self):
        """Check that shutter clocks that disagree are re-anchored to the
        header, so the close-card shift does not move the result.
        """
        timings = []
        for name in ("close_start_plus_10ms", "close_start_minus_5ms"):
            md = self.mutated(name)
            self.checkAgainstReference(self.mutations[name], md, name)
            t = self.compute(md)
            self.assertTrue(t.flags & ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN)
            timings.append(t)
        self.assertEqual(timings[0].centerMidpointMjdTai, timings[1].centerMidpointMjdTai)

        # Clocks that agree are not re-anchored.
        self.assertFalse(self.compute(self.base).flags & ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN)

    def testParamRange(self):
        """Check that Hall fits outside the nominal ranges, including a
        negative JERK0 with |JERK0| in range, are UNAVAILABLE.
        """
        for name in ("open_pivot1_out_of_range", "close_jerk2_out_of_range", "open_jerk0_negative"):
            for det in (0, 94, 188):
                t = self.compute(self.mutated(name), det)
                self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE, name)
                self.assertEqual(t.flags, ShutterTimingFlag.PARAM_RANGE, name)
                self.assertTrue(np.isnan(t.midpointMjdTai(np.array(2000.0), np.array(2000.0))))

    def testPropertyList(self):
        """Check that a `PropertyList` gives the same result as a `dict`."""
        doc = readFixture("MC_O_20260105_000250")
        pl = PropertyList()
        for k, v in doc["metadata"].items():
            pl.set(k, v)

        a, b = self.compute(pl), self.compute(doc["metadata"])
        xy = np.array([0.0, 2000.0, 4000.0])
        self.assertEqual(a.centerMidpointMjdTai, b.centerMidpointMjdTai)
        np.testing.assert_array_equal(a.midpointMjdTai(xy, xy), b.midpointMjdTai(xy, xy))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
