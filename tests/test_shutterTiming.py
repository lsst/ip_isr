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
import unittest

import numpy as np

import lsst.geom
import lsst.utils.tests
from lsst.afw.cameraGeom import FOCAL_PLANE, DetectorType, Orientation
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
from lsst.daf.base import PropertyList
from lsst.ip.isr.shutterTiming import (
    _DEGRADED_FLAGS,
    ShutterTimingConfig,
    ShutterTimingFlag,
    ShutterTimingStatus,
    _DetectorGeometry,
    _ShutterBeamModel,
    computeShutterTiming,
    loadShutterBeam,
)

TESTDIR = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(TESTDIR, "data", "shutterTiming")
BEAM_FILE = os.path.join(DATADIR, "beam_at_L3S1_z9.618_rot0_evaluated.tnt")
GEOMETRY_FILE = os.path.join(DATADIR, "lsstcam_detector_geometry.csv")
NO_CARDS = "MC_O_20250810_000030"
BASE = "MC_O_20260712_000100"

#: Agreement target (s): centre times, lever-arm terms, residuals, T_eff.
TOL_S = 10e-6
#: Lever arm (pixels) of the quadratic terms.
LEVER = 2000.0
LEVER_POWERS = np.array([LEVER, LEVER**2, LEVER, LEVER**2, LEVER**2])
SECONDS_PER_DAY = 86400.0


def readGeometries():
    """The detectors of the fixture CSV."""
    out = {}
    with open(GEOMETRY_FILE, newline="") as f:
        for r in csv.DictReader(f):
            v = {k: float(x) for k, x in r.items() if k not in ("name", "type", "physical_type")}
            det = int(v["detector"])
            out[det] = _DetectorGeometry(
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


def makeConfig():
    config = ShutterTimingConfig()
    config.beamFile = BEAM_FILE
    return config


def makeDetector(g, detType=DetectorType.SCIENCE, yaw=0.0):
    """An afw detector with the geometry of ``g`` (yaw 0)."""
    orientation = Orientation(lsst.geom.Point2D(*g.centerMm), lsst.geom.Point2D(*g.centerPixel),
                              lsst.geom.Angle(yaw, lsst.geom.degrees))
    bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(g.nx, g.ny))
    return DetectorWrapper(id=g.detectorId, bbox=bbox, pixelSize=(0.01, 0.01), orientation=orientation,
                           detType=detType).detector


class ShutterTimingTestBase(lsst.utils.tests.TestCase):
    """Shared fixtures and comparisons with the reference results."""

    @classmethod
    def setUpClass(cls):
        cls.geometries = readGeometries()
        cls.config = makeConfig()
        with open(os.path.join(DATADIR, "mutations.json")) as f:
            cls.mutations = json.load(f)["mutations"]
        cls.base = readFixture(BASE)["metadata"]

    def compute(self, md, det=94, config=None):
        return computeShutterTiming(md, self.geometries[det], config or self.config)

    def mutated(self, name):
        """The base metadata with a mutation's cards set and deleted."""
        m = self.mutations[name]
        md = dict(self.base)
        md.update(m.get("set", {}))
        for k in m.get("delete", []):
            del md[k]
        return md

    def assertAgrees(self, timing, expected, label, config):
        """Compare one detector's result with a fixture row."""
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
                     if (expected["qc_flags"] & int(_DEGRADED_FLAGS)
                         or not expected["max_abs_residual_s"] <= config.degradedResidual)
                     else ShutterTimingStatus.OK)
        self.assertEqual(timing.status, expStatus, label)
        self.assertEqual(timing.message == "", expStatus == ShutterTimingStatus.OK, label)

    def assertSamples(self, timings, samples, label):
        """A per-source time exactly where the sample status is not 2, and
        agreeing with the reference.
        """
        for det, x, y, t, status in samples:
            got = float(timings[det].tMidMjdTai(x, y))
            self.assertEqual(math.isfinite(got), status != ShutterTimingStatus.UNAVAILABLE,
                             f"{label} {det} ({x}, {y})")
            if t is not None:
                self.assertLessEqual(abs(got - t) * SECONDS_PER_DAY, TOL_S, f"{label} {det} ({x}, {y})")

    def checkAgainstReference(self, doc, metadata, label, config):
        """All detectors and samples of a fixture agree with the reference."""
        timings = {}
        for det, expected in doc["detectors"].items():
            t = computeShutterTiming(metadata, self.geometries[int(det)], config)
            timings[int(det)] = t
            self.assertAgrees(t, expected, f"{label} {det}", config)
            self.assertLessEqual(abs(t.focalPlaneMjdTai - doc["visit_mjd_tai"]) * SECONDS_PER_DAY, TOL_S)
        self.assertSamples(timings, doc["samples"], label)
        return timings


class ShutterTimingFixtureTestCase(ShutterTimingTestBase):
    """Agreement with the reference results on the fixture exposures."""

    def testFixtures(self):
        names = sorted(os.path.basename(p)[:-5] for p in glob.glob(os.path.join(DATADIR, "MC_O_*.json")))
        self.assertEqual(len(names), 5)
        for name in names:
            if name == NO_CARDS:
                continue
            doc = readFixture(name)
            self.assertEqual(len(doc["detectors"]), 189)
            timings = self.checkAgainstReference(doc, doc["metadata"], name, self.config)
            for t in timings.values():
                self.assertAlmostEqual(t.headerMidMjdTai, doc["header_mid_mjd_tai"], delta=1e-9)

    def testNoCards(self):
        doc = readFixture(NO_CARDS)
        for det in (0, 94, 188):
            timing = self.compute(doc["metadata"], det)
            self.assertEqual(timing.status, ShutterTimingStatus.UNAVAILABLE)
            self.assertEqual(timing.flags, ShutterTimingFlag.NO_PROFILE)
            self.assertIn("STARTTIME", timing.message)
            self.assertTrue(math.isnan(timing.centerMjdTai))
            self.assertTrue(all(math.isnan(c) for c in timing.coefficients))
            self.assertTrue(math.isnan(timing.focalPlaneMjdTai))
            self.assertTrue(math.isfinite(timing.headerMidMjdTai))
            x = np.array([0.0, 2000.0, 4000.0])
            self.assertTrue(np.all(np.isnan(timing.tMidMjdTai(x, x))))
            self.assertEqual(timing.summary()["status"], "UNAVAILABLE")

    def testBadBeamFile(self):
        """An empty or unreadable beamFile raises."""
        config = ShutterTimingConfig()
        self.assertEqual(config.beamFile, "")
        with self.assertRaises(ValueError):
            self.compute(self.base, config=config)
        config.beamFile = os.path.join(DATADIR, "no_such_beam.tnt")
        with self.assertRaises(Exception):
            self.compute(self.base, config=config)

    def testPropertyList(self):
        doc = readFixture("MC_O_20260105_000250")
        pl = PropertyList()
        for k, v in doc["metadata"].items():
            pl.set(k, v)
        a = self.compute(doc["metadata"])
        b = self.compute(pl)
        self.assertEqual(a, b)
        self.assertAgrees(b, doc["detectors"]["94"], "PropertyList", self.config)

    def testSummary(self):
        doc = readFixture("MC_O_20260105_000249")
        timing = self.compute(doc["metadata"])
        s = timing.summary()
        self.assertEqual(set(s), {"status", "flags", "message", "centerMjdTai", "focalPlaneMjdTai",
                                  "centerMinusHeaderMid", "maxAbsResidual"})
        self.assertEqual(s["status"], "OK")
        self.assertEqual(s["message"], "")
        self.assertAlmostEqual(s["centerMinusHeaderMid"],
                               (timing.centerMjdTai - doc["header_mid_mjd_tai"]) * SECONDS_PER_DAY)
        for v in s.values():
            self.assertIsInstance(v, (str, int, float))

    def testVectorized(self):
        timing = self.compute(self.base)
        x = np.linspace(-200, 4300, 30).reshape(5, 6)
        y = np.linspace(-150, 4200, 6)
        t = timing.tMidMjdTai(x, y)
        self.assertEqual(t.shape, (5, 6))
        for i in range(5):
            for j in range(6):
                np.testing.assert_array_equal(t[i, j], timing.tMidMjdTai(x[i, j], y[j]))
        t = timing.tMidMjdTai([np.nan, 10.0, np.inf, 10], [10.0, np.nan, 10.0, 10])
        np.testing.assert_array_equal(np.isfinite(t), [False, False, False, True])
        self.assertEqual(timing.tMidMjdTai(10, 10).shape, ())


class ShutterTimingMutationTestCase(ShutterTimingTestBase):
    """Mutated header cards, against the reference results in
    tests/data/shutterTiming/mutations.json.
    """

    def testMutations(self):
        for name, doc in self.mutations.items():
            md = self.mutated(name)
            if doc["status"] == ShutterTimingStatus.UNAVAILABLE:
                for det in (0, 94, 188):
                    self.assertEqual(self.compute(md, det).status, ShutterTimingStatus.UNAVAILABLE, name)
            else:
                self.checkAgainstReference(doc, md, name, self.config)

    def testCloseStartShifted(self):
        for name in ("close_start_plus_10ms", "close_start_minus_5ms"):
            t = self.compute(self.mutated(name))
            self.assertTrue(t.flags & ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN)
            self.assertTrue(t.flags & ShutterTimingFlag.PRE_CLOCK_EPOCH)
            self.assertEqual(t.status, ShutterTimingStatus.DEGRADED)
            self.assertIn("CLOCK_CLOSE_VS_OPEN", t.message)
        # Re-anchored: the shift of the close card does not move the result.
        a = self.compute(self.mutated("close_start_plus_10ms"))
        b = self.compute(self.mutated("close_start_minus_5ms"))
        self.assertEqual(a.centerMjdTai, b.centerMjdTai)
        # Re-anchoring needs MJD-BEG.
        md = self.mutated("close_start_plus_10ms")
        del md["MJD-BEG"]
        t = self.compute(md)
        self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE)
        self.assertTrue(t.flags & ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN)

    def testParamRange(self):
        """Out-of-range fits, including a negative JERK0 with |JERK0| in range
        (only the displacement-at-0.9 s check catches it): the mean profile.
        """
        for name in ("open_pivot1_out_of_range", "close_jerk2_out_of_range", "open_jerk0_negative"):
            t = self.compute(self.mutated(name))
            self.assertTrue(t.flags & ShutterTimingFlag.PARAM_RANGE, name)
            self.assertTrue(t.flags & ShutterTimingFlag.MEAN_PROFILE, name)
            self.assertEqual(t.status, ShutterTimingStatus.DEGRADED, name)

    def testOpenBegFlagOnly(self):
        t = self.compute(self.mutated("open_beg_minus_20ms"))
        self.assertEqual(t.flags, ShutterTimingFlag.CLOCK_OPEN_VS_BEG)
        self.assertEqual(t.status, ShutterTimingStatus.OK)

    def testBladeOverlap(self):
        """EXPTIME 0 or tiny: the blades overlap, UNAVAILABLE."""
        self.assertIn("Overlap", self.mutations["exptime_zero"]["reason"])
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

    def testPartialBladeOverlap(self):
        """EXPTIME 0.05 s: the blades overlap over part of the focal plane.

        No time at the detector centre: UNAVAILABLE.  A centre time but not
        at every grid point: DEGRADED, no quadratic and no per-source times.
        """
        md = dict(self.base)
        md["EXPTIME"] = 0.05
        md["SHUTTER CLOSE STARTTIME TAI MJD"] = (md["SHUTTER OPEN STARTTIME TAI MJD"]
                                                 + (0.05 + 0.00052) / SECONDS_PER_DAY)
        counts = dict(full=0, partial=0, none=0)
        for det, g in self.geometries.items():
            if not g.isScience:
                continue
            t = self.compute(md, det)
            if t.status == ShutterTimingStatus.UNAVAILABLE:
                counts["none"] += 1
                self.assertIn("detector centre", t.message)
                continue
            self.assertTrue(math.isfinite(t.centerMjdTai))
            self.assertTrue(math.isfinite(t.effectiveExposureTime))
            x = np.array([g.centerPixel[0], -0.5, g.nx - 0.5, 100.0])
            y = np.array([g.centerPixel[1], -0.5, g.ny - 0.5, 3000.0])
            if np.all(np.isfinite(t.coefficients)):
                counts["full"] += 1
                self.assertTrue(np.all(np.isfinite(t.tMidMjdTai(x, y))))
                continue
            counts["partial"] += 1
            self.assertEqual(t.status, ShutterTimingStatus.DEGRADED, det)
            self.assertTrue(math.isnan(t.maxAbsResidual))
            self.assertIn("no per-source times", t.message)
            self.assertTrue(np.all(np.isnan(t.tMidMjdTai(x, y))))
        self.assertGreater(min(counts.values()), 0, counts)

    def testOffDetectorLimit(self):
        """Times up to exactly ``offDetectorLimit`` pixels off; NaN beyond."""
        t = self.compute(self.base)
        self.assertEqual(t.status, ShutterTimingStatus.OK)
        limit = self.config.offDetectorLimit
        nx, ny = t.geometry.nx, t.geometry.ny
        cx, cy = t.geometry.centerPixel
        x = np.array([-0.5 - limit, -0.5 - limit - 1e-6, nx - 0.5 + limit, nx - 0.5 + limit + 1e-6,
                      cx, cx, cx, cx])
        y = np.array([cy, cy, cy, cy, -0.5 - limit, -0.5 - limit - 1e-6, ny - 0.5 + limit,
                      ny - 0.5 + limit + 1e-6])
        np.testing.assert_array_equal(np.isfinite(t.tMidMjdTai(x, y)), [True, False] * 4)

    def testUnreadableStartTime(self):
        t = self.compute(self.mutated("open_start_string"))
        self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE)
        self.assertEqual(t.flags, ShutterTimingFlag.NO_PROFILE)
        for card, values in (("SHUTTER CLOSE STARTTIME TAI MJD", (math.nan, None, "", True, math.inf, "1.0")),
                             ("SHUTTER OPEN SIDE", ("SIDEWAYS", None, 1.0))):
            for value in values:
                md = dict(self.base)
                md[card] = value
                t = self.compute(md)
                self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE, (card, value))
                self.assertEqual(t.flags, ShutterTimingFlag.NO_PROFILE, (card, value))

    def testUnusableHallFit(self):
        """A missing or non-numeric Hall-fit card, or another model: the mean
        profile (NO_PROFILE | MEAN_PROFILE, DEGRADED).
        """
        both = ShutterTimingFlag.NO_PROFILE | ShutterTimingFlag.MEAN_PROFILE
        for name in ("close_jerk1_string", "open_jerk0_missing"):
            t = self.compute(self.mutated(name))
            self.assertEqual(t.flags & both, both, name)
            self.assertEqual(t.status, ShutterTimingStatus.DEGRADED, name)
            self.assertIn("NO_PROFILE", t.message)
        ref = self.compute(self.mutated("close_jerk1_string"))
        for card, value in (("SHUTTER CLOSE HALLSENSORFIT JERK1", math.nan),
                            ("SHUTTER CLOSE HALLSENSORFIT JERK1", -math.inf),
                            ("SHUTTER CLOSE HALLSENSORFIT JERK1", "35000"),
                            ("SHUTTER CLOSE HALLSENSORFIT JERK1", False),
                            ("SHUTTER CLOSE MODEL", "FourJerksModel")):
            md = dict(self.base)
            md[card] = value
            t = self.compute(md)
            self.assertEqual(t.status, ShutterTimingStatus.DEGRADED, (card, value))
            self.assertEqual(t.flags, ref.flags, (card, value))
            self.assertEqual(t.centerMjdTai, ref.centerMjdTai, (card, value))
        # Ints are numbers.
        md = dict(self.base)
        md["SHUTTER CLOSE HALLSENSORFIT JERK2"] = int(md["SHUTTER CLOSE HALLSENSORFIT JERK2"])
        self.assertEqual(self.compute(md).status, ShutterTimingStatus.OK)

    def testGarbageNeverRaises(self):
        """Absurd card values give a result, never an exception."""
        rng = np.random.default_rng(42)
        keys = [k for k in self.base if k.startswith("SHUTTER") or k in ("MJD-BEG", "MJD-END", "EXPTIME")]
        garbage = [0.0, -1.0, 1e30, -1e30, 1e-30, math.nan, math.inf, None, True, "x", "", [1.0]]
        for _ in range(60):
            md = dict(self.base)
            for k in rng.choice(keys, size=3, replace=False):
                scaled = [md[k] * 1.5, md[k] + 1.0] if isinstance(md[k], float) else []
                choices = garbage + scaled
                md[k] = choices[rng.integers(len(choices))]
            if rng.random() < 0.2:
                del md[keys[rng.integers(len(keys))]]
            t = self.compute(md, det=int(rng.choice([0, 94, 188])))
            if t.status == ShutterTimingStatus.UNAVAILABLE:
                self.assertNotEqual(t.message, "")
            else:
                self.assertTrue(np.isfinite(t.centerMjdTai))


class ShutterTimingGeometryTestCase(ShutterTimingTestBase):
    """Detector geometry and the beam model."""

    def testFromDetector(self):
        for yaw in (0.0, 90.0, 180.0, 270.0):
            orientation = Orientation(lsst.geom.Point2D(-127.0, 211.5), lsst.geom.Point2D(2035.5, 1999.5),
                                      lsst.geom.Angle(yaw, lsst.geom.degrees))
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(4072, 4000))
            detector = DetectorWrapper(id=17, bbox=bbox, pixelSize=(0.01, 0.01),
                                       orientation=orientation).detector
            geom = _DetectorGeometry.fromDetector(detector)
            self.assertEqual((geom.detectorId, geom.nx, geom.ny), (17, 4072, 4000))
            self.assertFloatsAlmostEqual(np.array(geom.centerPixel), np.array([2035.5, 1999.5]))
            self.assertFloatsAlmostEqual(np.array(geom.centerMm), np.array([-127.0, 211.5]), atol=1e-9)
            # CCS = (DVCS y, DVCS x); afw's corner order is its own.
            pix = geom.pixelCorners()
            yd, xd = geom.pixelToCcs(pix[:, 0], pix[:, 1])
            for p in detector.getCorners(FOCAL_PLANE):
                self.assertLess(np.min(np.hypot(xd - p.getX(), yd - p.getY())), 1e-6)
            self.assertEqual(geom.bladeAxis(), "y" if yaw in (0.0, 180.0) else "x")

    def testFromDetectorTiming(self):
        """An afw detector equal to a fixture detector gives its timing."""
        g = self.geometries[94]
        a = self.compute(self.base, config=self.config, det=94)
        b = computeShutterTiming(self.base, makeDetector(g), self.config)
        self.assertLess(abs(a.centerMjdTai - b.centerMjdTai) * SECONDS_PER_DAY, 1e-6)
        self.assertFloatsAlmostEqual(np.array(a.coefficients), np.array(b.coefficients), rtol=1e-6)

    def testNonScienceDetector(self):
        """Non-science detectors are UNAVAILABLE (the beam table covers only
        the science beams), even where the geometry is covered: afw
        detectors by type, and the LSSTCam guider and wavefront geometries.
        """
        g = self.geometries[94]
        for detType in (DetectorType.WAVEFRONT, DetectorType.GUIDER, DetectorType.FOCUS):
            detector = makeDetector(g, detType)
            self.assertFalse(_DetectorGeometry.fromDetector(detector).isScience)
            t = computeShutterTiming(self.base, detector, self.config)
            self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE, detType)
            self.assertIn("not SCIENCE", t.message)
            self.assertTrue(np.isnan(t.tMidMjdTai(2000.0, 2000.0)))
        self.assertEqual(computeShutterTiming(self.base, makeDetector(g), self.config).status,
                         ShutterTimingStatus.OK)
        nonScience = sorted(d for d, geom in self.geometries.items() if not geom.isScience)
        self.assertEqual(nonScience, list(range(189, 205)))
        for det in nonScience + [None]:
            geom = self.geometries[det] if det is not None else dataclasses.replace(g, isScience=False)
            t = computeShutterTiming(self.base, geom, self.config)
            self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE, det)
            self.assertIn("not SCIENCE", t.message)
            self.assertTrue(np.isnan(t.centerMjdTai))

    def testOffDetector(self):
        g = self.geometries[0]
        off = g.offDetector([-0.5, 4071.5, -10.0, 2000.0, 4100.0], [-0.5, 3999.5, 5.0, 4030.0, -20.0])
        np.testing.assert_array_equal(off, [0.0, 0.0, 9.5, 30.5, 28.5])

    def testBeam(self):
        """Reading, the hull distance, and the per-path cache."""
        beam = loadShutterBeam(BEAM_FILE)
        self.assertIs(loadShutterBeam(BEAM_FILE), beam)
        p = np.array([[0.0, 0.0], [296.25, 0.0], [400.0, 0.0]])
        np.testing.assert_array_equal(beam.isOutside(p), [False, False, True])
        q, dist = beam._hullProject(p[2:])
        np.testing.assert_allclose(q, [[296.25, 0.0]], atol=1e-9)
        self.assertAlmostEqual(dist[0], 400.0 - 296.25)
        # A URI works as well as a path (and is cached separately).
        uriBeam = loadShutterBeam("file://" + BEAM_FILE)
        self.assertIsNot(uriBeam, beam)
        self.assertIsInstance(uriBeam, _ShutterBeamModel)
        x, y = np.array([0.0, 100.0, 300.0]), np.array([0.0, -50.0, 280.0])
        np.testing.assert_array_equal(uriBeam.quadratureNodes(x, y)[0], beam.quadratureNodes(x, y)[0])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
