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

import csv
import dataclasses
import glob
import json
import math
import os
import time
import unittest

import numpy as np

import lsst.geom
import lsst.utils.tests
from lsst.afw.cameraGeom import FOCAL_PLANE, DetectorType, Orientation
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
from lsst.daf.base import PropertyList
from lsst.ip.isr.shutterTiming import (
    DEGRADED_FLAGS,
    DetectorGeometry,
    ShutterBeamModel,
    ShutterTimingConfig,
    ShutterTimingFlag,
    ShutterTimingStatus,
    computeShutterTiming,
    loadShutterBeam,
)

TESTDIR = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(TESTDIR, "data", "shutterTiming")
BEAM_FILE = os.path.join(DATADIR, "beam_at_L3S1_z9.618_rot0_evaluated.tnt")
GEOMETRY_FILE = os.path.join(DATADIR, "lsstcam_detector_geometry.csv")
NO_CARDS = "MC_O_20250810_000030"

#: Agreement target (s): centre times, lever-arm terms, residuals, T_eff.
TOL_S = 10e-6
#: Lever arm (pixels) of the quadratic terms.
LEVER = 2000.0
LEVER_POWERS = np.array([LEVER, LEVER**2, LEVER, LEVER**2, LEVER**2])
SECONDS_PER_DAY = 86400.0


def readGeometries():
    """The detectors of the fixture CSV, as DetectorGeometry."""
    out = {}
    with open(GEOMETRY_FILE, newline="") as f:
        for r in csv.DictReader(f):
            v = {k: float(x) for k, x in r.items() if k not in ("name", "type", "physical_type")}
            det = int(v["detector"])
            out[det] = DetectorGeometry(
                detectorId=det, nx=int(v["nx"]), ny=int(v["ny"]),
                centerPixel=(v["center_x_pix"], v["center_y_pix"]),
                centerMm=(v["fp_center_x_mm"], v["fp_center_y_mm"]),
                jacobian=((v["dfp_dxpix_x"], v["dfp_dypix_x"]), (v["dfp_dxpix_y"], v["dfp_dypix_y"])),
                isScience=(r["type"] == "SCIENCE"),
            )
    return out


def readFixture(name):
    with open(os.path.join(DATADIR, name + ".json")) as f:
        return json.load(f)


class ReferenceChecks:
    """Comparisons with the reference results (mixin)."""

    def assertAgrees(self, timing, expected, label):
        """Compare one detector's result with a fixture row; return the
        absolute differences (s).
        """
        self.assertEqual(timing.axis, expected["axis"], label)
        self.assertEqual(int(timing.flags), expected["qc_flags"], label)
        diffs = dict(
            center=abs(timing.centerMjdTai - expected["center_mjd_tai"]) * SECONDS_PER_DAY,
            lever=np.max(np.abs(np.array(timing.coefficients) - expected["coefficients_s"]) * LEVER_POWERS),
            residual=abs(timing.maxAbsResidual - expected["max_abs_residual_s"]),
            teff=abs(timing.effectiveExposureTime - expected["effective_exposure_time_s"]),
        )
        for k, v in diffs.items():
            self.assertLessEqual(v, TOL_S, f"{label}: {k}")
        expStatus = (ShutterTimingStatus.DEGRADED
                     if (expected["qc_flags"] & int(DEGRADED_FLAGS)
                         or not expected["max_abs_residual_s"] <= self.config.degradedResidual)
                     else ShutterTimingStatus.OK)
        self.assertEqual(timing.status, expStatus, label)
        self.assertEqual(timing.message == "", expStatus == ShutterTimingStatus.OK, label)
        return diffs

    def assertSamples(self, timings, samples, label):
        worst = 0.0
        for det, x, y, t, status in samples:
            timing = timings[det]
            got = timing.tMidMjdTai(x, y)
            self.assertEqual(int(timing.sourceStatus(x, y)), status, f"{label} {det} ({x}, {y})")
            if t is None:
                self.assertTrue(np.isnan(got))
            else:
                worst = max(worst, abs(float(got) - t) * SECONDS_PER_DAY)
        self.assertLessEqual(worst, TOL_S, label)
        return worst


class ShutterTimingFixtureTestCase(ReferenceChecks, lsst.utils.tests.TestCase):
    """Agreement with the reference results on the fixture exposures."""

    @classmethod
    def setUpClass(cls):
        cls.geometries = readGeometries()
        cls.beam = loadShutterBeam(BEAM_FILE)
        cls.config = ShutterTimingConfig()
        cls.config.beamFile = BEAM_FILE

    def testFixtures(self):
        names = sorted(os.path.basename(p)[:-5] for p in glob.glob(os.path.join(DATADIR, "MC_O_*.json")))
        self.assertEqual(len(names), 5)
        worst = {}
        for name in names:
            if name == NO_CARDS:
                continue
            doc = readFixture(name)
            self.assertEqual(len(doc["detectors"]), 189)
            timings = {}
            for det, expected in doc["detectors"].items():
                timing = computeShutterTiming(doc["metadata"], self.geometries[int(det)], self.config)
                timings[int(det)] = timing
                for k, v in self.assertAgrees(timing, expected, f"{name} {det}").items():
                    worst[k] = max(worst.get(k, 0.0), v)
                fp = abs(timing.focalPlaneMjdTai - doc["visit_mjd_tai"]) * SECONDS_PER_DAY
                self.assertLessEqual(fp, TOL_S)
                worst["focalPlane"] = max(worst.get("focalPlane", 0.0), fp)
                self.assertEqual(timing.policy, doc["policy"])
                self.assertAlmostEqual(timing.headerMidMjdTai, doc["header_mid_mjd_tai"], delta=1e-9)
            worst["samples"] = max(worst.get("samples", 0.0),
                                   self.assertSamples(timings, doc["samples"], name))
        print("\nworst |computed - reference| (us): "
              + ", ".join(f"{k} {v * 1e6:.2g}" for k, v in worst.items()))

    def testNoCards(self):
        doc = readFixture(NO_CARDS)
        for det in (0, 94, 188):
            timing = computeShutterTiming(doc["metadata"], self.geometries[det], self.config)
            self.assertEqual(timing.status, ShutterTimingStatus.UNAVAILABLE)
            self.assertEqual(timing.flags, ShutterTimingFlag.NO_PROFILE)
            self.assertIn("STARTTIME", timing.message)
            self.assertTrue(math.isnan(timing.centerMjdTai))
            self.assertTrue(all(math.isnan(c) for c in timing.coefficients))
            self.assertTrue(math.isnan(timing.focalPlaneMjdTai))
            x = np.array([0.0, 2000.0, 4000.0])
            y = np.array([0.0, 2000.0, 3999.0])
            self.assertTrue(np.all(np.isnan(timing.tMidMjdTai(x, y))))
            np.testing.assert_array_equal(timing.sourceStatus(x, y), ShutterTimingStatus.UNAVAILABLE)
            self.assertEqual(timing.summary()["status"], "UNAVAILABLE")

    def testEmptyBeamFile(self):
        doc = readFixture("MC_O_20260712_000100")
        config = ShutterTimingConfig()
        self.assertEqual(config.beamFile, "")
        timing = computeShutterTiming(doc["metadata"], self.geometries[94], config)
        self.assertEqual(timing.status, ShutterTimingStatus.UNAVAILABLE)
        self.assertIn("beam", timing.message)
        self.assertTrue(np.isnan(timing.tMidMjdTai(2000.0, 2000.0)))
        # An explicit beam overrides the empty beamFile.
        timing = computeShutterTiming(doc["metadata"], self.geometries[94], config, beam=self.beam)
        self.assertNotEqual(timing.status, ShutterTimingStatus.UNAVAILABLE)

    def testPropertyList(self):
        doc = readFixture("MC_O_20260105_000250")
        pl = PropertyList()
        for k, v in doc["metadata"].items():
            pl.set(k, v)
        det = self.geometries[94]
        a = computeShutterTiming(doc["metadata"], det, self.config)
        b = computeShutterTiming(pl, det, self.config)
        self.assertEqual(a.status, b.status)
        self.assertEqual(a.flags, b.flags)
        self.assertEqual(a.centerMjdTai, b.centerMjdTai)
        self.assertEqual(a.coefficients, b.coefficients)
        self.assertEqual(a.focalPlaneMjdTai, b.focalPlaneMjdTai)
        self.assertAgrees(b, doc["detectors"]["94"], "PropertyList")
        # HIERARCH-prefixed keys are accepted too.
        c = computeShutterTiming({"HIERARCH " + k: v for k, v in doc["metadata"].items()}, det, self.config)
        self.assertEqual(a.centerMjdTai, c.centerMjdTai)

    def testMultiValuedCards(self):
        """A card with several values is malformed: the same flags and result
        as for that card missing (PropertyList) or non-numeric (mapping).
        """
        doc = readFixture("MC_O_20260712_000100")
        det = self.geometries[94]
        for card in ("SHUTTER CLOSE HALLSENSORFIT JERK1", "SHUTTER OPEN STARTTIME TAI MJD",
                     "SHUTTER OPEN SIDE", "EXPTIME"):
            missing = dict(doc["metadata"])
            del missing[card]
            expected = computeShutterTiming(missing, det, self.config)
            pl = PropertyList()
            for k, v in doc["metadata"].items():
                pl.set(k, v)
            pl.add(card, doc["metadata"][card])
            self.assertEqual(pl.valueCount(card), 2)
            md = dict(doc["metadata"])
            md[card] = [doc["metadata"][card]] * 2
            for metadata in (pl, md):
                t = computeShutterTiming(metadata, det, self.config)
                self.assertEqual(t.status, expected.status, card)
                self.assertEqual(t.flags, expected.flags, card)
                np.testing.assert_array_equal(t.centerMjdTai, expected.centerMjdTai)
            if card != "EXPTIME":  # a missing EXPTIME only skips a clock check
                self.assertNotEqual(expected.flags, ShutterTimingFlag.NONE, card)

    def testSummary(self):
        doc = readFixture("MC_O_20260105_000249")
        timing = computeShutterTiming(doc["metadata"], self.geometries[94], self.config)
        s = timing.summary()
        self.assertEqual(set(s), {"status", "flags", "message", "centerMjdTai", "focalPlaneMjdTai",
                                  "centerMinusHeaderMid", "maxAbsResidual", "policy"})
        self.assertEqual(s["status"], "OK")
        self.assertEqual(s["message"], "")
        self.assertAlmostEqual(s["centerMinusHeaderMid"],
                               (timing.centerMjdTai - doc["header_mid_mjd_tai"]) * SECONDS_PER_DAY)
        for v in s.values():
            self.assertIsInstance(v, (str, int, float))

    def testVectorized(self):
        doc = readFixture("MC_O_20260712_000100")
        timing = computeShutterTiming(doc["metadata"], self.geometries[94], self.config)
        x = np.linspace(-200, 4300, 30).reshape(5, 6)
        y = np.linspace(-150, 4200, 6)
        t = timing.tMidMjdTai(x, y)
        st = timing.sourceStatus(x, y)
        self.assertEqual(t.shape, (5, 6))
        self.assertEqual(st.dtype, np.uint8)
        for i in range(5):
            for j in range(6):
                self.assertEqual(st[i, j], timing.sourceStatus(x[i, j], y[j]))
                np.testing.assert_array_equal(t[i, j], timing.tMidMjdTai(x[i, j], y[j]))
        st = timing.sourceStatus([np.nan, 10.0, np.inf], [10.0, np.nan, 10.0])
        np.testing.assert_array_equal(st, ShutterTimingStatus.UNAVAILABLE)

    def testMaskedArrays(self):
        """Masked entries of masked arrays are UNAVAILABLE with a NaN time;
        unmasked entries are as for plain arrays.
        """
        doc = readFixture("MC_O_20260712_000100")
        timing = computeShutterTiming(doc["metadata"], self.geometries[94], self.config)
        x = np.array([10.0, 2000.0, 3000.0, 4000.0])
        y = np.array([20.0, 1000.0, 2000.0, 3900.0])
        mask = np.array([False, True, False, True])
        plainT = timing.tMidMjdTai(x, y)
        plainS = timing.sourceStatus(x, y)
        self.assertTrue(np.all(np.isfinite(plainT)))
        for mx, my in ((np.ma.array(x, mask=mask), y), (x, np.ma.array(y, mask=mask)),
                       (np.ma.array(x.astype(int), mask=mask), np.ma.array(y, mask=False))):
            t = timing.tMidMjdTai(mx, my)
            st = timing.sourceStatus(mx, my)
            self.assertNotIsInstance(t, np.ma.MaskedArray)
            self.assertTrue(np.all(np.isnan(t[mask])))
            np.testing.assert_array_equal(st[mask], ShutterTimingStatus.UNAVAILABLE)
            np.testing.assert_array_equal(st[~mask], plainS[~mask])
            np.testing.assert_allclose(t[~mask], timing.tMidMjdTai(np.asarray(mx, dtype=float)[~mask],
                                                                   np.asarray(my, dtype=float)[~mask]))


class ShutterTimingMutationTestCase(ReferenceChecks, lsst.utils.tests.TestCase):
    """Mutated header cards, against the reference results in
    tests/data/shutterTiming/mutations.json.
    """

    @classmethod
    def setUpClass(cls):
        cls.geometries = readGeometries()
        cls.config = ShutterTimingConfig()
        cls.config.beamFile = BEAM_FILE
        with open(os.path.join(DATADIR, "mutations.json")) as f:
            cls.mutations = json.load(f)["mutations"]
        cls.base = readFixture("MC_O_20260712_000100")["metadata"]

    def compute(self, md, det=94, config=None):
        return computeShutterTiming(md, self.geometries[det], config or self.config)

    def checkAgainstReference(self, name, config=None):
        """All detectors and samples of a mutation agree with the reference.
        """
        doc = self.mutations[name]
        timings = {}
        saved = self.config
        if config is not None:
            self.config = config
        try:
            self._checkAgainstReference(doc, name, timings)
        finally:
            self.config = saved
        return timings

    def _checkAgainstReference(self, doc, name, timings):
        for det, expected in doc["detectors"].items():
            t = self.compute(doc["metadata"], int(det))
            timings[int(det)] = t
            self.assertAgrees(t, expected, f"{name} {det}")
            self.assertEqual(t.policy, doc["policy"])
            self.assertLessEqual(abs(t.focalPlaneMjdTai - doc["visit_mjd_tai"]) * SECONDS_PER_DAY, TOL_S)
        self.assertSamples(timings, doc["samples"], name)
        return timings

    def testCloseStartShifted(self):
        for name in ("close_start_plus_10ms", "close_start_minus_5ms"):
            timings = self.checkAgainstReference(name)
            t = timings[94]
            self.assertTrue(t.flags & ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN)
            self.assertTrue(t.flags & ShutterTimingFlag.PRE_CLOCK_EPOCH)
            self.assertEqual(t.policy, "header_anchor")
            self.assertEqual(t.status, ShutterTimingStatus.DEGRADED)
            self.assertIn("CLOCK_CLOSE_VS_OPEN", t.message)
        # Re-anchored: the shift of the close card does not move the result.
        a = self.compute(self.mutations["close_start_plus_10ms"]["metadata"])
        b = self.compute(self.mutations["close_start_minus_5ms"]["metadata"])
        self.assertEqual(a.centerMjdTai, b.centerMjdTai)
        # Re-anchoring needs MJD-BEG.
        md = dict(self.mutations["close_start_plus_10ms"]["metadata"])
        del md["MJD-BEG"]
        t = self.compute(md)
        self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE)
        self.assertTrue(t.flags & ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN)

    def testA1(self):
        """A nonzero A1 per direction, against the reference with the same
        ``a1_mm`` (a dropped or sign-flipped A1 moves times by ~0.5 ms).
        """
        for name in ("a1_decreasing", "a1_increasing"):
            doc = self.mutations[name]
            config = ShutterTimingConfig()
            config.beamFile = BEAM_FILE
            config.a1Decreasing = doc["a1_mm"]["-1"]
            config.a1Increasing = doc["a1_mm"]["1"]
            self.assertNotEqual(config.a1Decreasing + config.a1Increasing, 0.0)
            timings = self.checkAgainstReference(name, config)
            # The A1 of the exposure's direction is applied (~0.5 ms).
            t0 = self.compute(readFixture(doc["base"])["metadata"]).centerMjdTai
            self.assertGreater(abs(timings[94].centerMjdTai - t0) * SECONDS_PER_DAY, 100e-6, name)

    def testNegativeJerk0(self):
        """|JERK0| in range but negative: only the displacement-at-0.9 s range
        check catches it; the mean profile is used, as the reference.
        """
        timings = self.checkAgainstReference("open_jerk0_negative")
        t = timings[94]
        self.assertTrue(t.flags & ShutterTimingFlag.PARAM_RANGE)
        self.assertTrue(t.flags & ShutterTimingFlag.MEAN_PROFILE)
        self.assertEqual(t.status, ShutterTimingStatus.DEGRADED)

    def testBladeOverlap(self):
        """EXPTIME 0 or tiny: the blades overlap, UNAVAILABLE (the reference
        rejects the exposure).
        """
        doc = self.mutations["exptime_zero"]
        self.assertEqual(doc["status"], ShutterTimingStatus.UNAVAILABLE)
        self.assertIn("Overlap", doc["reason"])
        for exptime in (0.0, 1e-6, 0.01):
            md = dict(self.base)
            md["EXPTIME"] = exptime
            for det in (0, 94, 188):
                t = self.compute(md, det)
                self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE, (exptime, det))
                self.assertIn("overlapping", t.message)
        # Close start == open start, with EXPTIME consistent (no
        # re-anchoring).
        md = dict(self.base)
        md["SHUTTER CLOSE STARTTIME TAI MJD"] = md["SHUTTER OPEN STARTTIME TAI MJD"]
        md["EXPTIME"] = 0.0
        self.assertEqual(self.compute(md).status, ShutterTimingStatus.UNAVAILABLE)

    def testOffDetectorLimit(self):
        """Exactly ``offDetectorLimit`` pixels off: DEGRADED; beyond:
        UNAVAILABLE.
        """
        t = self.compute(self.base)
        self.assertEqual(t.status, ShutterTimingStatus.OK)
        limit = self.config.offDetectorLimit
        nx, ny = t.geometry.nx, t.geometry.ny
        cx, cy = t.geometry.centerPixel
        x = np.array([-0.5 - limit, -0.5 - limit - 1e-6, nx - 0.5 + limit, nx - 0.5 + limit + 1e-6,
                      cx, cx, cx, cx])
        y = np.array([cy, cy, cy, cy, -0.5 - limit, -0.5 - limit - 1e-6, ny - 0.5 + limit,
                      ny - 0.5 + limit + 1e-6])
        expected = [ShutterTimingStatus.DEGRADED, ShutterTimingStatus.UNAVAILABLE] * 4
        np.testing.assert_array_equal(t.sourceStatus(x, y), expected)
        tt = t.tMidMjdTai(x, y)
        self.assertTrue(np.all(np.isfinite(tt[0::2])))
        self.assertTrue(np.all(np.isnan(tt[1::2])))

    def testParamRange(self):
        for name in ("open_pivot1_out_of_range", "close_jerk2_out_of_range"):
            timings = self.checkAgainstReference(name)
            t = timings[94]
            self.assertTrue(t.flags & ShutterTimingFlag.PARAM_RANGE)
            self.assertTrue(t.flags & ShutterTimingFlag.MEAN_PROFILE)
            self.assertEqual(t.status, ShutterTimingStatus.DEGRADED)
            self.assertEqual(t.policy, "profile")

    def testOpenBegFlagOnly(self):
        timings = self.checkAgainstReference("open_beg_minus_20ms")
        t = timings[94]
        self.assertEqual(t.flags, ShutterTimingFlag.CLOCK_OPEN_VS_BEG)
        self.assertEqual(t.status, ShutterTimingStatus.OK)

    def testUnreadableStartTime(self):
        doc = self.mutations["open_start_string"]
        self.assertEqual(doc["status"], ShutterTimingStatus.UNAVAILABLE)
        t = self.compute(doc["metadata"])
        self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE)
        self.assertEqual(t.flags, ShutterTimingFlag.NO_PROFILE)
        for value in (math.nan, None, "", True, math.inf):
            md = dict(self.base)
            md["SHUTTER CLOSE STARTTIME TAI MJD"] = value
            self.assertEqual(self.compute(md).status, ShutterTimingStatus.UNAVAILABLE, repr(value))
        md = dict(self.base)
        md["SHUTTER OPEN SIDE"] = "SIDEWAYS"
        self.assertEqual(self.compute(md).flags, ShutterTimingFlag.NO_PROFILE)

    def testUnusableHallFit(self):
        """A missing or malformed Hall-fit card: the mean profile
        (NO_PROFILE | MEAN_PROFILE, DEGRADED), as the reference.
        """
        both = ShutterTimingFlag.NO_PROFILE | ShutterTimingFlag.MEAN_PROFILE
        for name in ("close_jerk1_string", "open_jerk0_missing"):
            timings = self.checkAgainstReference(name)
            t = timings[94]
            self.assertEqual(t.flags & both, both, name)
            self.assertEqual(t.status, ShutterTimingStatus.DEGRADED, name)
            self.assertIn("NO_PROFILE", t.message)
        ref = self.compute(self.mutations["close_jerk1_string"]["metadata"])
        for card, value in (("SHUTTER CLOSE HALLSENSORFIT JERK1", math.nan),
                            ("SHUTTER CLOSE HALLSENSORFIT JERK1", math.inf),
                            ("SHUTTER CLOSE HALLSENSORFIT JERK1", -math.inf),
                            ("SHUTTER CLOSE HALLSENSORFIT JERK1", "35000 mm/s3"),
                            ("SHUTTER CLOSE MODEL", "FourJerksModel")):
            md = dict(self.base)
            md[card] = value
            t = self.compute(md)
            self.assertEqual(t.status, ShutterTimingStatus.DEGRADED, (card, value))
            self.assertEqual(t.flags, ref.flags, (card, value))
            self.assertEqual(t.centerMjdTai, ref.centerMjdTai, (card, value))
        # Numeric strings are numbers.
        md = dict(self.base)
        md["SHUTTER CLOSE HALLSENSORFIT JERK2"] = str(md["SHUTTER CLOSE HALLSENSORFIT JERK2"])
        self.assertEqual(self.compute(md).centerMjdTai, self.compute(self.base).centerMjdTai)

    def testGarbageNeverRaises(self):
        """Absurd but numeric values give a result, never an exception."""
        rng = np.random.default_rng(42)
        keys = [k for k, v in self.base.items() if k.startswith("SHUTTER") and isinstance(v, float)]
        for _ in range(40):
            md = dict(self.base)
            for k in rng.choice(keys, size=3, replace=False):
                md[k] = float(rng.choice([0.0, -1.0, 1e30, -1e30, 1e-30, md[k] * 1.5, md[k] + 1.0]))
            t = self.compute(md, det=int(rng.choice([0, 94, 188])))
            self.assertIn(t.status, set(ShutterTimingStatus))
            if t.status == ShutterTimingStatus.UNAVAILABLE:
                self.assertNotEqual(t.message, "")
            else:
                self.assertTrue(np.isfinite(t.centerMjdTai))


class ShutterTimingGeometryTestCase(lsst.utils.tests.TestCase):
    """DetectorGeometry and the beam model."""

    def testFromDetector(self):
        for yaw in (0.0, 90.0, 180.0, 270.0):
            orientation = Orientation(lsst.geom.Point2D(-127.0, 211.5), lsst.geom.Point2D(2035.5, 1999.5),
                                      lsst.geom.Angle(yaw, lsst.geom.degrees))
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(4072, 4000))
            detector = DetectorWrapper(id=17, bbox=bbox, pixelSize=(0.01, 0.01),
                                       orientation=orientation).detector
            geom = DetectorGeometry.fromDetector(detector)
            self.assertEqual((geom.detectorId, geom.nx, geom.ny), (17, 4072, 4000))
            self.assertFloatsAlmostEqual(np.array(geom.centerPixel), np.array([2035.5, 1999.5]))
            self.assertFloatsAlmostEqual(np.array(geom.centerMm), np.array([-127.0, 211.5]), atol=1e-9)
            corners = np.array([[p.getX(), p.getY()] for p in detector.getCorners(FOCAL_PLANE)])
            pix = geom._pixelCorners()
            for c in corners:  # afw's corner order is its own: match each corner
                xd, yd = geom.pixelToDvcs(pix[:, 0], pix[:, 1])
                self.assertLess(np.min(np.hypot(xd - c[0], yd - c[1])), 1e-6)
            xc, yc = geom.pixelToCcs(pix[:, 0], pix[:, 1], -1)
            xd, yd = geom.pixelToDvcs(pix[:, 0], pix[:, 1])
            self.assertFloatsAlmostEqual(xc, -yd)
            self.assertFloatsAlmostEqual(yc, -xd)
            self.assertEqual(geom._bladeAxis(), "y" if yaw in (0.0, 180.0) else "x")

    def testFromDetectorTiming(self):
        """An afw detector equal to a fixture detector gives its timing."""
        g = readGeometries()[94]
        orientation = Orientation(lsst.geom.Point2D(*g.centerMm), lsst.geom.Point2D(*g.centerPixel),
                                  lsst.geom.Angle(0.0, lsst.geom.degrees))
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(g.nx, g.ny))
        detector = DetectorWrapper(id=94, bbox=bbox, pixelSize=(0.01, 0.01), orientation=orientation).detector
        config = ShutterTimingConfig()
        config.beamFile = BEAM_FILE
        doc = readFixture("MC_O_20260712_000100")
        a = computeShutterTiming(doc["metadata"], detector, config)
        b = computeShutterTiming(doc["metadata"], g, config)
        self.assertLess(abs(a.centerMjdTai - b.centerMjdTai) * SECONDS_PER_DAY, 1e-6)
        self.assertFloatsAlmostEqual(np.array(a.coefficients), np.array(b.coefficients), rtol=1e-6)

    def testNonScienceDetector(self):
        """An afw detector that is not SCIENCE is UNAVAILABLE (the beam table
        covers only the science beams), even where the geometry is covered.
        """
        g = readGeometries()[94]
        orientation = Orientation(lsst.geom.Point2D(*g.centerMm), lsst.geom.Point2D(*g.centerPixel),
                                  lsst.geom.Angle(0.0, lsst.geom.degrees))
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(g.nx, g.ny))
        config = ShutterTimingConfig()
        config.beamFile = BEAM_FILE
        md = readFixture("MC_O_20260712_000100")["metadata"]
        for detType in (DetectorType.WAVEFRONT, DetectorType.GUIDER, DetectorType.FOCUS):
            detector = DetectorWrapper(id=94, bbox=bbox, pixelSize=(0.01, 0.01), orientation=orientation,
                                       detType=detType).detector
            t = computeShutterTiming(md, detector, config)
            self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE, detType)
            self.assertIn("not SCIENCE", t.message)
            self.assertTrue(np.isnan(t.tMidMjdTai(2000.0, 2000.0)))
            self.assertFalse(DetectorGeometry.fromDetector(detector).isScience)
        detector = DetectorWrapper(id=94, bbox=bbox, pixelSize=(0.01, 0.01), orientation=orientation,
                                   detType=DetectorType.SCIENCE).detector
        self.assertTrue(DetectorGeometry.fromDetector(detector).isScience)
        self.assertEqual(computeShutterTiming(md, detector, config).status, ShutterTimingStatus.OK)

    def testNonScienceGeometry(self):
        """A DetectorGeometry with ``isScience=False`` is UNAVAILABLE, as the
        afw detector it describes: the LSSTCam guider and wavefront detectors
        (ids 189-204), and a science geometry flagged as non-science.
        """
        config = ShutterTimingConfig()
        config.beamFile = BEAM_FILE
        md = readFixture("MC_O_20260712_000100")["metadata"]
        geometries = readGeometries()
        nonScience = sorted(d for d, g in geometries.items() if not g.isScience)
        self.assertEqual(nonScience, list(range(189, 205)))
        for det in nonScience:
            t = computeShutterTiming(md, geometries[det], config)
            self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE, det)
            self.assertIn("not SCIENCE", t.message)
            self.assertTrue(np.isnan(t.centerMjdTai))
            np.testing.assert_array_equal(t.sourceStatus([0.0, 2000.0], [0.0, 1000.0]),
                                          ShutterTimingStatus.UNAVAILABLE)
        g = geometries[94]
        self.assertTrue(g.isScience)
        self.assertEqual(computeShutterTiming(md, g, config).status, ShutterTimingStatus.OK)
        t = computeShutterTiming(md, dataclasses.replace(g, isScience=False), config)
        self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE)

    def testOffDetector(self):
        g = readGeometries()[0]
        off = g.offDetector([-0.5, 4071.5, -10.0, 2000.0, 4100.0], [-0.5, 3999.5, 5.0, 4030.0, -20.0])
        np.testing.assert_array_equal(off, [0.0, 0.0, 9.5, 30.5, 28.5])

    def testBeam(self):
        beam = ShutterBeamModel.fromFile(BEAM_FILE)
        np.testing.assert_array_equal(beam.levels, [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99])
        self.assertIs(loadShutterBeam(BEAM_FILE), loadShutterBeam(BEAM_FILE))
        d = beam.hullDistance(np.array([0.0, 296.25, 400.0]), np.array([0.0, 0.0, 0.0]))
        self.assertEqual(d[0], 0.0)
        self.assertEqual(d[1], 0.0)
        self.assertAlmostEqual(d[2], 400.0 - 296.25)
        # A URI works as well as a path.
        uriBeam = ShutterBeamModel.fromFile("file://" + BEAM_FILE)
        np.testing.assert_array_equal(uriBeam.hullDistance(350.0, 10.0), beam.hullDistance(350.0, 10.0))

    def testLoadShutterBeamCache(self):
        """The cache keys on the absolute URI and holds a read-only model."""
        cwd = os.getcwd()
        try:
            os.chdir(DATADIR)
            relative = loadShutterBeam(os.path.basename(BEAM_FILE))
            os.chdir(TESTDIR)
            self.assertIs(loadShutterBeam(BEAM_FILE), relative)
            self.assertIs(loadShutterBeam("file://" + BEAM_FILE), relative)
            self.assertIs(loadShutterBeam(os.path.join("data", "shutterTiming", os.path.basename(BEAM_FILE))),
                          relative)
            # The bare name no longer resolves after the chdir.
            with self.assertRaises(Exception):
                loadShutterBeam(os.path.basename(BEAM_FILE))
        finally:
            os.chdir(cwd)
        for name in ("_grid", "_r", "_points", "_hullEq", "_nodes", "_weights"):
            arr = getattr(relative, name)
            self.assertFalse(arr.flags.writeable, name)
            with self.assertRaises(ValueError):
                arr.flat[0] = 0.0
        levels = relative.levels
        levels[0] = 0.5  # a copy
        self.assertEqual(relative.levels[0], 0.01)


class ShutterTimingSpeedTestCase(lsst.utils.tests.TestCase):
    def testSpeed(self):
        config = ShutterTimingConfig()
        config.beamFile = BEAM_FILE
        loadShutterBeam(BEAM_FILE)
        geometries = readGeometries()
        md = readFixture("MC_O_20260712_000100")["metadata"]
        computeShutterTiming(md, geometries[94], config)
        best = math.inf
        for det in (0, 30, 94, 120, 188):
            t0 = time.perf_counter()
            computeShutterTiming(md, geometries[det], config)
            best = min(best, time.perf_counter() - t0)
        print(f"\ncomputeShutterTiming: {best * 1e3:.1f} ms per detector (best of 5)")
        self.assertLess(best, 0.05)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
