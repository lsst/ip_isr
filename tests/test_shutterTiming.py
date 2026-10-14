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
import os
import unittest

import numpy as np

import lsst.geom
import lsst.utils.tests
from lsst.afw.cameraGeom import Orientation
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
from lsst.afw.image import ExposureF, ExposureSummaryStats, VisitInfo
from lsst.daf.base import PropertyList
from lsst.ip.isr import IsrTaskLSST, IsrTaskLSSTConfig, ShutterMotionProfile
from lsst.ip.isr.shutterTiming import (
    ShutterTiming,
    ShutterTimingConfig,
    computeShutterTiming,
    getShutterTiming,
)

TESTDIR = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(TESTDIR, "data", "shutterTiming")
BASE = "MC_O_20260712_000100"

#: Agreement target (s): centre and per-source times.
TOL_S = 10e-6
SECONDS_PER_DAY = 86400.0
#: The reference's flag for re-anchored start times (CLOCK_CLOSE_VS_OPEN).
REANCHORED_FLAG = 16


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


def makeExposure(metadata, detector, requestedExposureTime=True, exposureId=0):
    """Build a 1x1 exposure (pixels are not used) with the given metadata
    (`dict`) and detector.

    Unless ``requestedExposureTime`` is False, the metadata also gets the
    ``EXPTIME REQUESTED`` card from ``EXPTIME``, as the obs_lsst raw
    formatter adds it.  ``exposureId`` sets the VisitInfo ID.
    """
    exposure = ExposureF(1, 1)
    exposure.setDetector(detector)
    exposure.info.setVisitInfo(VisitInfo(id=exposureId))
    for key, value in metadata.items():
        exposure.metadata.set(key, value)
    if requestedExposureTime:
        exposure.metadata.set("EXPTIME REQUESTED", metadata["EXPTIME"])
    return exposure


class ShutterTimingTestCase(lsst.utils.tests.TestCase):
    """Agreement with the reference results, and the cases that occur in
    processing.
    """

    @classmethod
    def setUpClass(cls):
        cls.config = ShutterTimingConfig()
        doc = cls.baseDoc = readFixture(BASE)
        cls.base = doc["metadata"]
        cls.detectors = {int(det): makeDetector(int(det), expected["geometry"])
                         for det, expected in doc["detectors"].items()}

    def compute(self, md, det=94):
        return getShutterTiming(makeExposure(md, self.detectors[det]), self.config)

    def changed(self, card, value):
        """Return the base metadata with one card set to ``value``."""
        md = dict(self.base)
        md[card] = value
        return md

    def checkAgainstReference(self, doc, metadata, label):
        """Check all detectors and samples of a fixture against the
        reference.
        """
        # Every detector: re-anchoring and the centre time.
        timings = {}
        for det, expected in doc["detectors"].items():
            exposure = makeExposure(metadata, makeDetector(int(det), expected["geometry"]))
            t = timings[int(det)] = getShutterTiming(exposure, self.config)
            label2 = f"{label} {det}"

            self.assertEqual(t.reanchored, bool(expected["qc_flags"] & REANCHORED_FLAG), label2)
            self.assertLessEqual(abs(t.centerMidpointMjdTai - expected["center_mjd_tai"])*SECONDS_PER_DAY,
                                 TOL_S, label2)

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

    def testRecord(self):
        """Check that the times survive a round trip through metadata cards
        and summary statistics, and that those without them give none.
        """
        timing = self.compute(self.base)
        y = np.array([-100.0, 0.0, 2000.0, 3999.0])

        metadata = PropertyList()
        timing.toMetadata(metadata)
        readBack = ShutterTiming.fromMetadata(metadata)
        self.assertFloatsEqual(readBack.midpointMjdTai(y, y), timing.midpointMjdTai(y, y))
        for field in ("centerExposureTime", "reanchored", "source"):
            self.assertEqual(getattr(readBack, field), getattr(timing, field), field)
        self.assertIsNone(ShutterTiming.fromMetadata(PropertyList()))

        summary = ExposureSummaryStats()
        timing.updateSummaryStats(summary)
        readBack = ShutterTiming.fromSummaryStats(summary)
        self.assertFloatsEqual(readBack.midpointMjdTai(y, y), timing.midpointMjdTai(y, y))
        self.assertEqual(readBack.centerExposureTime, timing.centerExposureTime)
        self.assertIsNone(ShutterTiming.fromSummaryStats(ExposureSummaryStats()))

    def testProfiles(self):
        """Check that the camera's motion profiles give the times of the
        header cards of the same exposure, to within the correction for the
        measured start positions (<= 31 us; here they are nominal to 0.2 um).
        """
        def profile(direction):
            return ShutterMotionProfile.readText(os.path.join(
                TESTDIR, "data", f"MC_O_20251203_000317_shutterMotionProfile{direction}.json"))

        # The header cards of that exposure, as written by the camera.
        md = {k: v for k, v in self.base.items() if not k.startswith("SHUTTER")}
        for direction in ("Open", "Close"):
            p = profile(direction)
            i = p.fit_name.index("hallSensorFit")
            prefix = f"SHUTTER {direction.upper()} "
            md.update({prefix + "STARTTIME TAI ISOT": p.metadata["startTime_tai"],
                       prefix + "STARTTIME TAI MJD": p.metadata["startTime_mjd"],
                       prefix + "SIDE": p.metadata["side"],
                       prefix + "MODEL": p.fit_model,
                       prefix + "HALLSENSORFIT MODELSTARTTIME": p.fit_start_time[i],
                       prefix + "HALLSENSORFIT PIVOTPOINT1": p.fit_pivot1[i],
                       prefix + "HALLSENSORFIT PIVOTPOINT2": p.fit_pivot2[i],
                       prefix + "HALLSENSORFIT JERK0": p.fit_jerk0[i],
                       prefix + "HALLSENSORFIT JERK1": p.fit_jerk1[i],
                       prefix + "HALLSENSORFIT JERK2": p.fit_jerk2[i]})
        md["EXPTIME"] = 10.0
        fromHeader = self.compute(md, 0)
        fromProfiles = computeShutterTiming(profile("Open"), profile("Close"), self.detectors[0],
                                            md["MJD-BEG"], 10.0, self.config)

        self.assertEqual((fromHeader.source, fromProfiles.source), ("header", "profile"))
        y = np.array([0.0, 2000.0, 3999.0])
        dt = (fromProfiles.midpointMjdTai(y, y) - fromHeader.midpointMjdTai(y, y))*SECONDS_PER_DAY
        self.assertTrue(np.all(np.abs(dt) <= 31e-6), dt)

    def testIsrTask(self):
        """Check that ISR records the times of the header cards in the
        metadata.
        """
        config = IsrTaskLSSTConfig()
        config.doShutterTiming = True
        task = IsrTaskLSST(config=config)

        exposure = makeExposure(self.base, self.detectors[94], requestedExposureTime=False)
        task.recordShutterTiming(exposure, None, None, self.base["EXPTIME"])
        recorded = ShutterTiming.fromMetadata(exposure.metadata)
        expected = self.compute(self.base)

        self.assertEqual(recorded.centerMidpointMjdTai, expected.centerMidpointMjdTai)
        self.assertFloatsEqual(recorded.chebyshevCoefficients, expected.chebyshevCoefficients)

        # Later steps get the recorded times, not ones from the header.
        exposure.metadata.set("SHUTTER OPEN HALLSENSORFIT JERK0", 0.0)
        self.assertEqual(getShutterTiming(exposure).centerMidpointMjdTai, expected.centerMidpointMjdTai)

    def testNoCards(self):
        """Check that exposures before the shutter cards have no corrected
        times.
        """
        md = {k: v for k, v in self.base.items() if not k.startswith("SHUTTER")}
        for det in (0, 94, 188):
            self.assertIsNone(self.compute(md, det))

    def testReanchoring(self):
        """Check that shutter timestamps that disagree with the requested
        exposure time are re-anchored to the header, so the shift of the
        close start time does not move the result.
        """
        card = "SHUTTER CLOSE STARTTIME TAI MJD"
        timings = []
        for shiftMs in (10.0, -5.0):
            md = self.changed(card, self.base[card] + shiftMs*1e-3/SECONDS_PER_DAY)
            self.checkAgainstReference(self.baseDoc["reanchored"], md, f"close start {shiftMs:+} ms")
            t = self.compute(md)
            self.assertTrue(t.reanchored)
            timings.append(t)
        self.assertEqual(timings[0].centerMidpointMjdTai, timings[1].centerMidpointMjdTai)

        # Clocks that agree are not re-anchored.
        self.assertFalse(self.compute(self.base).reanchored)

        # Images made before the requested exposure time was recorded: the
        # timestamps are used as they are, except for the exposures known to
        # need re-anchoring.
        exposure = makeExposure(md, self.detectors[94], requestedExposureTime=False)
        unchecked = getShutterTiming(exposure, self.config)
        self.assertFalse(unchecked.reanchored)
        self.assertNotEqual(unchecked.centerMidpointMjdTai, timings[1].centerMidpointMjdTai)

        exposure = makeExposure(md, self.detectors[94], requestedExposureTime=False, exposureId=2026031900172)
        listed = getShutterTiming(exposure, self.config)
        self.assertTrue(listed.reanchored)
        self.assertEqual(listed.centerMidpointMjdTai, timings[1].centerMidpointMjdTai)

    def testParamRange(self):
        """Check that Hall fits outside the nominal ranges, including a
        negative JERK0 with |JERK0| in range, have no corrected times.
        """
        jerk0 = "SHUTTER OPEN HALLSENSORFIT JERK0"
        changes = (
            ("SHUTTER OPEN HALLSENSORFIT PIVOTPOINT1", 0.30),
            ("SHUTTER CLOSE HALLSENSORFIT JERK2", 50000.0),
            (jerk0, -self.base[jerk0]),
        )
        for card, value in changes:
            for det in (0, 94, 188):
                self.assertIsNone(self.compute(self.changed(card, value), det), card)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
