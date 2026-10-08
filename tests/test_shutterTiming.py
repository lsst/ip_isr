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
from lsst.ip.isr._shutterTrajectory import ThreeJerksParams, ThreeJerksTrajectory
from lsst.ip.isr.shutterTiming import (
    _BEAM_MARGIN_MM,
    _FIT_CARDS,
    ShutterTimingConfig,
    ShutterTimingFlag,
    ShutterTimingStatus,
    _DetectorGeometry,
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
#: The reference's flag bits that this implementation keeps.
KEPT_FLAGS = int(ShutterTimingFlag.NO_PROFILE | ShutterTimingFlag.PARAM_RANGE
                 | ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN | ShutterTimingFlag.BEAM_EXTRAPOLATED)


def readGeometries():
    """The science detectors of the fixture CSV."""
    out = {}
    with open(GEOMETRY_FILE, newline="") as f:
        for r in csv.DictReader(f):
            if r["type"] != "SCIENCE":
                continue
            v = {k: float(x) for k, x in r.items() if k not in ("name", "type", "physical_type")}
            det = int(v["detector"])
            out[det] = _DetectorGeometry(
                detectorId=det, nx=int(v["nx"]), ny=int(v["ny"]),
                centerPixel=(v["center_x_pix"], v["center_y_pix"]),
                centerMm=(v["fp_center_x_mm"], v["fp_center_y_mm"]),
                jacobian=((v["dfp_dxpix_x"], v["dfp_dypix_x"]), (v["dfp_dxpix_y"], v["dfp_dypix_y"])),
            )
    return out


def readFixture(name):
    with open(os.path.join(DATADIR, name + ".json")) as f:
        return json.load(f)


def makeDetector(detectorId, centerMm, centerPixel, size, yaw=0.0):
    orientation = Orientation(lsst.geom.Point2D(*centerMm), lsst.geom.Point2D(*centerPixel),
                              lsst.geom.Angle(yaw, lsst.geom.degrees))
    bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(*size))
    return DetectorWrapper(id=detectorId, bbox=bbox, pixelSize=(0.01, 0.01),
                           orientation=orientation).detector


class ShutterTimingTestBase(lsst.utils.tests.TestCase):
    """Shared fixtures and comparisons with the reference results."""

    @classmethod
    def setUpClass(cls):
        cls.geometries = readGeometries()
        cls.config = ShutterTimingConfig()
        cls.config.beamFile = BEAM_FILE
        with open(os.path.join(DATADIR, "mutations.json")) as f:
            cls.mutations = json.load(f)["mutations"]
        cls.base = readFixture(BASE)["metadata"]

    def compute(self, md, det=94):
        return computeShutterTiming(md, self.geometries[det], md["EXPTIME"], self.config)

    def mutated(self, name):
        """The base metadata with a mutation's cards set."""
        md = dict(self.base)
        md.update(self.mutations[name]["set"])
        return md

    def checkAgainstReference(self, doc, metadata, label):
        """All detectors and samples of a fixture agree with the reference.

        The reference's DEGRADED is OK here, and its flags are compared on
        the kept bits.  The reference gives no time 150 px off the detector,
        where this implementation evaluates the quadratic like anywhere else;
        callers only pass on-detector positions, so those samples are only
        checked to be finite.
        """
        timings = {}
        for det, expected in doc["detectors"].items():
            t = timings[int(det)] = computeShutterTiming(
                metadata, self.geometries[int(det)], metadata["EXPTIME"], self.config
            )
            label2 = f"{label} {det}"
            self.assertEqual(t.status, ShutterTimingStatus.OK, label2)
            self.assertEqual(t.message, "", label2)
            self.assertEqual(t.axis, expected["axis"], label2)
            self.assertEqual(int(t.flags), expected["qc_flags"] & KEPT_FLAGS, label2)
            diffs = dict(
                center=abs(t.centerMjdTai - expected["center_mjd_tai"]) * SECONDS_PER_DAY,
                lever=np.max(np.abs(np.array(t.coefficients) - expected["coefficients_s"]) * LEVER_POWERS),
                residual=abs(t.maxAbsResidual - expected["max_abs_residual_s"]),
                teff=abs(t.effectiveExposureTime - expected["effective_exposure_time_s"]),
                focalPlane=abs(t.focalPlaneMjdTai - doc["visit_mjd_tai"]) * SECONDS_PER_DAY,
            )
            for k, v in diffs.items():
                self.assertLessEqual(v, TOL_S, f"{label2}: {k}")
        for det, x, y, expected, _ in doc["samples"]:
            got = float(timings[det].tMidMjdTai(x, y))
            self.assertTrue(math.isfinite(got), f"{label} {det} ({x}, {y})")
            if expected is not None:
                self.assertLessEqual(abs(got - expected) * SECONDS_PER_DAY, TOL_S,
                                     f"{label} {det} ({x}, {y})")
        return timings


class ShutterTimingTestCase(ShutterTimingTestBase):
    """Agreement with the reference results, and the cases that occur in
    processing.
    """

    def testFixtures(self):
        names = sorted(os.path.basename(p)[:-5] for p in glob.glob(os.path.join(DATADIR, "MC_O_*.json")))
        self.assertEqual(len(names), 5)
        for name in names:
            if name == NO_CARDS:
                continue
            doc = readFixture(name)
            self.assertEqual(len(doc["detectors"]), 189)
            timings = self.checkAgainstReference(doc, doc["metadata"], name)
            self.assertTrue(any(t.flags & ShutterTimingFlag.BEAM_EXTRAPOLATED for t in timings.values()))
            for t in timings.values():
                self.assertAlmostEqual(t.headerMidMjdTai, doc["header_mid_mjd_tai"], delta=1e-9)

    def testNoCards(self):
        """Exposures before the shutter cards: UNAVAILABLE, no times."""
        doc = readFixture(NO_CARDS)
        for det in (0, 94, 188):
            t = self.compute(doc["metadata"], det)
            self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE)
            self.assertEqual(t.flags, ShutterTimingFlag.NO_PROFILE)
            self.assertTrue(math.isnan(t.centerMjdTai))
            self.assertTrue(math.isnan(t.focalPlaneMjdTai))
            self.assertTrue(np.all(np.isnan(t.tMidMjdTai([0.0, 2000.0], [0.0, 2000.0]))))

    def testReanchoring(self):
        """Shutter clocks that disagree: re-anchored to the header, so the
        close-card shift does not move the result.
        """
        timings = []
        for name in ("close_start_plus_10ms", "close_start_minus_5ms"):
            md = self.mutated(name)
            self.checkAgainstReference(self.mutations[name], md, name)
            t = self.compute(md)
            self.assertTrue(t.flags & ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN)
            timings.append(t)
        self.assertEqual(timings[0].centerMjdTai, timings[1].centerMjdTai)
        self.assertFalse(self.compute(self.base).flags & ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN)

    def testParamRange(self):
        """Hall fits outside the nominal ranges, including a negative JERK0
        with |JERK0| in range: UNAVAILABLE.
        """
        for name in ("open_pivot1_out_of_range", "close_jerk2_out_of_range", "open_jerk0_negative"):
            for det in (0, 94, 188):
                t = self.compute(self.mutated(name), det)
                self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE, name)
                self.assertEqual(t.flags, ShutterTimingFlag.PARAM_RANGE, name)
                self.assertTrue(np.isnan(t.tMidMjdTai(2000.0, 2000.0)))

    def testBladeOverlap(self):
        """EXPTIME 0: the blades overlap, UNAVAILABLE."""
        md = self.mutated("exptime_zero")
        for det in (0, 94, 188):
            t = self.compute(md, det)
            self.assertEqual(t.status, ShutterTimingStatus.UNAVAILABLE, det)
            self.assertTrue(math.isnan(t.centerMjdTai))

    def testPropertyList(self):
        doc = readFixture("MC_O_20260105_000250")
        pl = PropertyList()
        for k, v in doc["metadata"].items():
            pl.set(k, v)
        self.assertEqual(self.compute(pl), self.compute(doc["metadata"]))
        one = dict(doc, detectors={"94": doc["detectors"]["94"]},
                   samples=[s for s in doc["samples"] if s[0] == 94])
        self.checkAgainstReference(one, pl, "PropertyList")

    def testVectorized(self):
        t = self.compute(self.base)
        x = np.linspace(-0.5, 4071.5, 30).reshape(5, 6)
        y = np.linspace(-0.5, 3999.5, 6)
        tt = t.tMidMjdTai(x, y)
        self.assertEqual(tt.shape, (5, 6))
        for i in range(5):
            for j in range(6):
                self.assertEqual(tt[i, j], t.tMidMjdTai(x[i, j], y[j]))
        self.assertEqual(t.tMidMjdTai(10, 10).shape, ())
        got = t.tMidMjdTai([np.nan, 10.0, np.inf, 10.0], [10.0, np.nan, 10.0, 10.0])
        np.testing.assert_array_equal(np.isfinite(got), [False, False, False, True])

    def testFromDetector(self):
        """afw detector geometry: corners, blade axis, and the same timing as
        the equivalent fixture geometry.
        """
        for yaw in (0.0, 90.0, 180.0, 270.0):
            detector = makeDetector(17, (-127.0, 211.5), (2035.5, 1999.5), (4072, 4000), yaw)
            geom = _DetectorGeometry.fromDetector(detector)
            self.assertEqual((geom.detectorId, geom.nx, geom.ny), (17, 4072, 4000))
            # CCS = (DVCS y, DVCS x); afw's corner order is its own.
            pix = geom.pixelCorners()
            yd, xd = geom.pixelToCcs(pix[:, 0], pix[:, 1])
            for p in detector.getCorners(FOCAL_PLANE):
                self.assertLess(np.min(np.hypot(xd - p.getX(), yd - p.getY())), 1e-6)
            self.assertEqual(geom.bladeAxis(), "y" if yaw in (0.0, 180.0) else "x")
        g = self.geometries[94]
        a = self.compute(self.base, 94)
        b = computeShutterTiming(self.base, makeDetector(94, g.centerMm, g.centerPixel, (g.nx, g.ny)),
                                 self.base["EXPTIME"], self.config)
        self.assertLess(abs(a.centerMjdTai - b.centerMjdTai) * SECONDS_PER_DAY, 1e-6)
        self.assertFloatsAlmostEqual(np.array(a.coefficients), np.array(b.coefficients), rtol=1e-6)


class ShutterBeamTestCase(lsst.utils.tests.TestCase):
    """Beam-table interpolation (the corner detectors use the triangulated
    and extrapolated regions; the fixtures check them end to end).
    """

    def setUp(self):
        self.beam = loadShutterBeam(BEAM_FILE)
        rows = np.loadtxt(BEAM_FILE, skiprows=2)
        self.nodes = {}
        for x, y, q, _, s in rows:
            self.nodes.setdefault((x, y), {})[q] = s - x
        self.points = np.array(list(self.nodes))
        self.r = np.array([[v[q] for q in sorted(v)] for v in self.nodes.values()])

    def testNodes(self):
        """The table is reproduced at every node."""
        r = self.beam._residual(self.points[:, 0], self.points[:, 1])
        np.testing.assert_allclose(r, self.r, rtol=0, atol=1e-9)

    def testBilinear(self):
        """The centre of a fully tabulated cell is the mean of its corners."""
        x0, y0, h = 0.0, 0.0, 42.25
        corners = [self.nodes[(x0 + i * h, y0 + j * h)] for i in (0, 1) for j in (0, 1)]
        expected = np.mean([[c[q] for q in sorted(c)] for c in corners], axis=0)
        r = self.beam._residual(np.array([x0 + h / 2]), np.array([y0 + h / 2]))[0]
        np.testing.assert_allclose(r, expected, rtol=0, atol=1e-9)

    def testOutsideHull(self):
        """Outside the table's hull: held at the nearest hull point, up to
        the margin; NaN beyond.
        """
        k = np.argmax(self.points.sum(axis=1))  # the hull vertex farthest along (1, 1)
        p = self.points[k]
        out = np.array([p + d / np.sqrt(2.0) for d in (5.0, _BEAM_MARGIN_MM - 1.0, _BEAM_MARGIN_MM + 1.0)])
        self.assertTrue(np.all(self.beam.isOutside(out)))
        r = self.beam._residual(out[:, 0], out[:, 1])
        np.testing.assert_allclose(r[:2], np.array([self.r[k]] * 2), rtol=0, atol=1e-5)
        self.assertTrue(np.all(np.isnan(r[2])))


class ShutterTrajectoryTestCase(lsst.utils.tests.TestCase):
    """Inversion of the ThreeJerksModelv1 trajectory."""

    def testInversion(self):
        md = readFixture(BASE)["metadata"]
        for which in ("OPEN", "CLOSE"):
            fit = ThreeJerksParams(*(md[f"SHUTTER {which} HALLSENSORFIT {c}"] for c in _FIT_CARDS))
            for sign, start in ((-1, 750.76), (1, -0.05)):
                traj = ThreeJerksTrajectory(fit, start, sign)
                self.assertGreater(traj.tStop, 0.85)
                t = np.concatenate([np.linspace(0.0, traj.tStop, 1001), [fit.pivot1, fit.pivot2]])
                d = traj.sva(t)[0]
                self.assertTrue(np.all(np.diff(d[:1001]) > 0))
                inv = traj.timeOfDisplacement(d)
                np.testing.assert_allclose(traj.sva(inv)[0], d, rtol=0, atol=1e-9)
                # The time is ill-conditioned only where the blade stops.
                moving = t < 0.99 * traj.tStop
                np.testing.assert_allclose(inv[moving], t[moving], rtol=0, atol=1e-12)
                self.assertTrue(np.all(np.isnan(traj.timeOfDisplacement([-1e-6, traj.dMax + 1e-6]))))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
