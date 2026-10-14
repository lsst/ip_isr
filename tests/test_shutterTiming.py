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
from lsst.afw.cameraGeom import FOCAL_PLANE, Orientation
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
from lsst.daf.base import PropertyList
from lsst.ip.isr.shutterTiming import (
    ShutterTimingConfig,
    ShutterTimingFlag,
    ShutterTimingStatus,
    _DetectorGeometry,
    computeShutterTiming,
)

TESTDIR = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(TESTDIR, "data", "shutterTiming")
BEAM_FILE = os.path.join(DATADIR, "beam_at_L3S1_z9.618_rot0_evaluated.tnt")
BASE = "MC_O_20260712_000100"

#: Agreement target (s): centre times, lever-arm terms, residuals, T_eff.
TOL_S = 10e-6
#: Lever arm (pixels) of the quadratic terms.
LEVER = 2000.0
LEVER_POWERS = np.array([LEVER, LEVER**2, LEVER, LEVER**2, LEVER**2])
SECONDS_PER_DAY = 86400.0
#: The reference's flag bits that this implementation keeps.
KEPT_FLAGS = int(ShutterTimingFlag.NO_PROFILE | ShutterTimingFlag.PARAM_RANGE
                 | ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN | ShutterTimingFlag.BEAM_EXTRAPOLATED)


def readFixture(name):
    with open(os.path.join(DATADIR, name + ".json")) as f:
        return json.load(f)


def makeDetector(detectorId, geometry, yaw=0.0):
    """Build an afw detector with 10 um pixels.

    Parameters
    ----------
    detectorId : `int`
        The detector ID.
    geometry : `dict`
        ``center_mm`` (focal-plane position of the centre), ``center_pix``
        (its pixel coordinates) and ``size`` (pixels), as in the fixtures.
    yaw : `float`, optional
        Rotation (degrees); every LSSTCam detector has 0.

    Returns
    -------
    detector : `lsst.afw.cameraGeom.Detector`
        The detector.
    """
    orientation = Orientation(lsst.geom.Point2D(*geometry["center_mm"]),
                              lsst.geom.Point2D(*geometry["center_pix"]),
                              lsst.geom.Angle(yaw, lsst.geom.degrees))
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
        cls.config.beamFile = BEAM_FILE
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
        reference.

        The reference's DEGRADED is OK here, and its flags are compared on
        the kept bits.  The reference gives no time 150 px off the detector,
        where this implementation evaluates the quadratic like anywhere else;
        callers only pass on-detector positions, so those samples are only
        checked to be finite.
        """
        # Every detector: status, flags and the fitted quantities.
        timings = {}
        for det, expected in doc["detectors"].items():
            t = timings[int(det)] = computeShutterTiming(
                metadata, makeDetector(int(det), expected["geometry"]), metadata["EXPTIME"], self.config
            )
            label2 = f"{label} {det}"

            self.assertEqual(t.status, ShutterTimingStatus.OK, label2)
            self.assertEqual(t.message, "", label2)
            self.assertEqual(t.axis, expected["axis"], label2)
            self.assertEqual(int(t.flags), expected["qc_flags"] & KEPT_FLAGS, label2)

            diffs = dict(
                center=abs(t.centerMidpointMjdTai - expected["center_mjd_tai"]) * SECONDS_PER_DAY,
                lever=np.max(np.abs(np.array(t.coefficients) - expected["coefficients_s"]) * LEVER_POWERS),
                residual=abs(t.maxAbsResidual - expected["max_abs_residual_s"]),
                teff=abs(t.effectiveExposureTime - expected["effective_exposure_time_s"]),
                focalPlane=abs(t.focalPlaneMidpointMjdTai - doc["visit_mjd_tai"]) * SECONDS_PER_DAY,
            )
            for k, v in diffs.items():
                self.assertLessEqual(v, TOL_S, f"{label2}: {k}")

        # Per-source times at the sample positions.
        for det, x, y, expected, _ in doc["samples"]:
            got = float(timings[det].midpointMjdTai(np.array(x), np.array(y)))
            self.assertTrue(math.isfinite(got), f"{label} {det} ({x}, {y})")
            if expected is not None:
                self.assertLessEqual(abs(got - expected)*SECONDS_PER_DAY, TOL_S,
                                     f"{label} {det} ({x}, {y})")

        return timings

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
            timings = self.checkAgainstReference(doc, doc["metadata"], name)

            # The outer detectors extend past the beam table.
            self.assertTrue(any(t.flags & ShutterTimingFlag.BEAM_EXTRAPOLATED for t in timings.values()))

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
        self.assertEqual(self.compute(pl), self.compute(doc["metadata"]))

    def testFromDetector(self):
        """Check the linearized detector geometry: corners and blade axis at
        every yaw.
        """
        geometry = dict(center_mm=(-127.0, 211.5), center_pix=(2035.5, 1999.5), size=(4072, 4000))
        for yaw in (0.0, 90.0, 180.0, 270.0):
            detector = makeDetector(17, geometry, yaw)
            geom = _DetectorGeometry.fromDetector(detector)
            self.assertEqual((geom.detectorId, geom.nx, geom.ny), (17, 4072, 4000))

            # CCS = (DVCS y, DVCS x); afw's corner order is its own.
            pix = geom.pixelCorners()
            yd, xd = geom.pixelToCcs(pix[:, 0], pix[:, 1])
            for p in detector.getCorners(FOCAL_PLANE):
                self.assertLess(np.min(np.hypot(xd - p.getX(), yd - p.getY())), 1e-6)

            self.assertEqual(geom.bladeAxis(), "y" if yaw in (0.0, 180.0) else "x")


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
