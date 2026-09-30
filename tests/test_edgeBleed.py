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

import numpy as np

import lsst.geom as geom
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.image as afwImage
import lsst.utils.tests
import lsst.ip.isr.isrFunctions as ipIsrFunctions
from lsst.ip.isr.isrTask import IsrTaskConfig
from lsst.ip.isr.masking import DECamEdgeBleedMaskTask

AMP_WIDTH = 200
AMP_HEIGHT = 800
SKY = 1000.0
NOISE = 10.0
DIP = 200.0


def makeDetector():
    """Two side-by-side amplifiers: A reads out at the top, B at the bottom."""
    camBuilder = cameraGeom.Camera.Builder("testCam")
    detBuilder = camBuilder.add("testDet", 0)
    detBuilder.setBBox(geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(2*AMP_WIDTH, AMP_HEIGHT)))
    detBuilder.setPixelSize(geom.Extent2D(0.015, 0.015))
    detBuilder.setOrientation(cameraGeom.Orientation())
    for name, x0, corner in (("A", 0, cameraGeom.ReadoutCorner.UL),
                             ("B", AMP_WIDTH, cameraGeom.ReadoutCorner.LL)):
        bbox = geom.Box2I(geom.Point2I(x0, 0), geom.Extent2I(AMP_WIDTH, AMP_HEIGHT))
        amp = cameraGeom.Amplifier.Builder()
        amp.setName(name)
        amp.setBBox(bbox)
        amp.setRawBBox(bbox)
        amp.setRawDataBBox(bbox)
        amp.setReadoutCorner(corner)
        amp.setGain(1.0)
        amp.setReadNoise(5.0)
        amp.setSaturation(60000.0)
        detBuilder.append(amp)
    return camBuilder.finish()["testDet"]


class MaskDECamEdgeBleedTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.detector = makeDetector()
        self.exposure = afwImage.ExposureF(self.detector.getBBox())
        self.exposure.setDetector(self.detector)
        rng = np.random.default_rng(12345)
        self.exposure.image.array[:] = SKY + rng.normal(0.0, NOISE, self.exposure.image.array.shape)
        self.exposure.variance.array[:] = NOISE**2
        self.satBit = self.exposure.mask.getPlaneBitMask("SAT")

    def addSatBlock(self, x0, y0, width, height):
        """Mark a width x height block as SAT (image value irrelevant)."""
        self.exposure.mask.array[y0:y0 + height, x0:x0 + width] |= self.satBit
        self.exposure.image.array[y0:y0 + height, x0:x0 + width] = 60000.0

    def addDip(self, x0, y0, width, height):
        """Depress a width x height block below sky."""
        self.exposure.image.array[y0:y0 + height, x0:x0 + width] -= DIP

    def runMasking(self):
        before = self.exposure.mask.array.copy()
        ipIsrFunctions.maskDECamEdgeBleed(self.exposure)
        return before, self.exposure.mask.array

    def test_topReadEdgeBleedIsMasked(self):
        # 100 x 200 = 20000 SAT pixels touching the top edge of amp A.
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)
        # Dip in the 100 rows nearest the top edge, full width of amp A.
        self.addDip(0, AMP_HEIGHT - 100, AMP_WIDTH, 100)

        before, after = self.runMasking()

        # height 100 -> margin int(100*0.125) + 1 = 13 -> 113 rows masked.
        expected = before.copy()
        expected[AMP_HEIGHT - 113:, 0:AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(after, expected)

    def test_noDipIsNotMasked(self):
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)

        before, after = self.runMasking()

        np.testing.assert_array_equal(after, before)

    def test_smallFootprintIsNotMasked(self):
        # 50 x 100 = 5000 SAT pixels: below satMinArea.
        self.addSatBlock(50, AMP_HEIGHT - 100, 50, 100)
        self.addDip(0, AMP_HEIGHT - 100, AMP_WIDTH, 100)

        before, after = self.runMasking()

        np.testing.assert_array_equal(after, before)

    def test_largeFootprintIsMasked(self):
        # 300 x 400 = 120000 SAT pixels: tall DECam bleeds have footprints of
        # 1e5-2e5 pixels, so there is no upper area cut.  Only x 0-49 of the
        # dip rows are unsaturated, still enough low pixels per row to be
        # detectable.
        self.addSatBlock(50, 400, 300, 400)
        self.addDip(0, AMP_HEIGHT - 100, AMP_WIDTH, 100)

        before, after = self.runMasking()

        expected = before.copy()
        expected[AMP_HEIGHT - 113:, 0:AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(after, expected)

    def test_dipInFewUsablePixelsIsMasked(self):
        # Saturation covers 160 of the 200 columns in the 100 rows nearest
        # the read edge, and the dip is in 25 of the remaining 40 columns:
        # 25 low pixels per row is far below the old count threshold of 30,
        # but 62 per cent of the usable pixels are low.
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)
        self.addSatBlock(40, AMP_HEIGHT - 100, 160, 100)
        self.addDip(0, AMP_HEIGHT - 100, 25, 100)

        before, after = self.runMasking()

        expected = before.copy()
        expected[AMP_HEIGHT - 113:, 0:AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(after, expected)

    def test_sparseLowPixelsAreNotADip(self):
        # 2 low pixels per row, 2 per cent of the 100 usable columns in the
        # check rows, is below minLowFraction: not an edge bleed.
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)
        self.addDip(0, AMP_HEIGHT - 100, 2, 100)

        before, after = self.runMasking()

        np.testing.assert_array_equal(after, before)

    def test_shallowDipInFirstRowsIsMasked(self):
        # A dip in 10 of the 100 usable columns over the 4 rows nearest the
        # read edge: 2 per cent averaged over the 20 check rows, but 8 per
        # cent over the 5-row short window.  4 rows -> int(4*0.125) + 1 = 1
        # -> 5 rows masked.
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)
        self.addDip(0, AMP_HEIGHT - 4, 10, 4)

        before, after = self.runMasking()

        expected = before.copy()
        expected[AMP_HEIGHT - 5:, 0:AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(after, expected)

    def test_shortWindowCanBeDisabled(self):
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)
        self.addDip(0, AMP_HEIGHT - 4, 10, 4)

        before = self.exposure.mask.array.copy()
        ipIsrFunctions.maskDECamEdgeBleed(self.exposure, nRowsCheckShort=0)

        np.testing.assert_array_equal(self.exposure.mask.array, before)

    def test_blockedRowsBeyondBorderAreMasked(self):
        # Both edges of amp A carry a 10-row BAD border (as DECam defects
        # do).  A drained read register leaves rows 10-59 from the read edge
        # fully BAD as well, with no measurable dip beyond them: those 50
        # blocked rows are themselves the edge bleed, so the mask covers
        # them plus the margin (60 -> int(60*0.125) + 1 = 8 -> 68 rows).
        badBit = self.exposure.mask.getPlaneBitMask("BAD")
        self.exposure.mask.array[:10, 0:AMP_WIDTH] |= badBit
        self.exposure.mask.array[AMP_HEIGHT - 60:, 0:AMP_WIDTH] |= badBit
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)

        before, after = self.runMasking()

        expected = before.copy()
        expected[AMP_HEIGHT - 68:, 0:AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(after, expected)

    def test_normalBorderIsNotABleed(self):
        # The same 10-row BAD border at both edges, but nothing beyond it:
        # the read-edge border is not counted as blocked rows.
        badBit = self.exposure.mask.getPlaneBitMask("BAD")
        self.exposure.mask.array[:10, 0:AMP_WIDTH] |= badBit
        self.exposure.mask.array[AMP_HEIGHT - 10:, 0:AMP_WIDTH] |= badBit
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)

        before, after = self.runMasking()

        np.testing.assert_array_equal(after, before)

    def test_fewBlockedRowsNeedADip(self):
        # 3 blocked rows beyond the border are below minBlockedRows; with no
        # dip beyond them nothing is masked.
        badBit = self.exposure.mask.getPlaneBitMask("BAD")
        self.exposure.mask.array[:10, 0:AMP_WIDTH] |= badBit
        self.exposure.mask.array[AMP_HEIGHT - 13:, 0:AMP_WIDTH] |= badBit
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)

        before, after = self.runMasking()

        np.testing.assert_array_equal(after, before)

    def test_blockedRowsExtendAMeasuredDip(self):
        # 30 blocked rows beyond a 10-row border, then a 40-row dip: the
        # height is measured from the physical edge through the dip
        # (80 -> int(80*0.125) + 1 = 11 -> 91 rows).
        badBit = self.exposure.mask.getPlaneBitMask("BAD")
        self.exposure.mask.array[:10, 0:AMP_WIDTH] |= badBit
        self.exposure.mask.array[AMP_HEIGHT - 40:, 0:AMP_WIDTH] |= badBit
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)
        self.addDip(0, AMP_HEIGHT - 80, AMP_WIDTH, 40)

        before, after = self.runMasking()

        expected = before.copy()
        expected[AMP_HEIGHT - 91:, 0:AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(after, expected)

    def test_footprintFarFromEdgeIsNotMasked(self):
        # Footprint top at row 599, 200 rows short of the read edge.
        self.addSatBlock(50, 400, 100, 200)
        self.addDip(0, AMP_HEIGHT - 100, AMP_WIDTH, 100)

        before, after = self.runMasking()

        np.testing.assert_array_equal(after, before)

    def test_bottomReadEdgeBleedIsMasked(self):
        # Amp B reads out at the bottom.
        self.addSatBlock(AMP_WIDTH + 50, 0, 100, 200)
        self.addDip(AMP_WIDTH, 0, AMP_WIDTH, 100)

        before, after = self.runMasking()

        expected = before.copy()
        # height 100 -> margin int(100*0.125) + 1 = 13 -> 113 rows masked.
        expected[:113, AMP_WIDTH:2*AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(after, expected)

    def test_dipInOtherAmpIsNotMasked(self):
        # Footprint in amp A (top read edge); dip only in amp B's top rows,
        # which is not B's read edge and has no footprint.
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)
        self.addDip(AMP_WIDTH, AMP_HEIGHT - 100, AMP_WIDTH, 100)

        before, after = self.runMasking()

        np.testing.assert_array_equal(after, before)

    def test_unusableEdgeRowsAreSkippedButMasked(self):
        # The 20 rows nearest the read edge are NO_DATA (as after trimming
        # on DECam), so the dip can only be seen from row 20 inward.  The
        # measured height still counts from the physical edge.
        noDataBit = self.exposure.mask.getPlaneBitMask("NO_DATA")
        self.exposure.mask.array[AMP_HEIGHT - 20:, 0:AMP_WIDTH] |= noDataBit
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)
        self.addDip(0, AMP_HEIGHT - 100, AMP_WIDTH, 80)

        before, after = self.runMasking()

        expected = before.copy()
        expected[AMP_HEIGHT - 113:, 0:AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(after, expected)

    def test_suspectPixelsCountAsLow(self):
        # DECam flags the rows nearest the read edge SUSPECT; those pixels
        # must still be counted when confirming the dip.
        suspectBit = self.exposure.mask.getPlaneBitMask("SUSPECT")
        self.exposure.mask.array[AMP_HEIGHT - 35:, 0:AMP_WIDTH] |= suspectBit
        self.addSatBlock(50, AMP_HEIGHT - 200, 100, 200)
        self.addDip(0, AMP_HEIGHT - 100, AMP_WIDTH, 100)

        before, after = self.runMasking()

        expected = before.copy()
        expected[AMP_HEIGHT - 113:, 0:AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(after, expected)

    def test_nonDefaultParameters(self):
        # A 100-row footprint ending 150 rows short of the edge is only a
        # candidate with a larger approachRows; a 50 ADU dip is only seen
        # with a smaller nSigma; a zero marginFraction masks height + 1 rows.
        self.addSatBlock(50, AMP_HEIGHT - 250, 200, 100)
        self.exposure.image.array[AMP_HEIGHT - 100:, 0:AMP_WIDTH] -= 50.0

        before = self.exposure.mask.array.copy()
        ipIsrFunctions.maskDECamEdgeBleed(self.exposure, approachRows=150, nSigma=3.0,
                                          marginFraction=0.0)

        expected = before.copy()
        expected[AMP_HEIGHT - 101:, 0:AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(self.exposure.mask.array, expected)

    def test_alternateSaturatedMaskName(self):
        self.exposure.mask.addMaskPlane("MYSAT")
        mySatBit = self.exposure.mask.getPlaneBitMask("MYSAT")
        self.exposure.mask.array[AMP_HEIGHT - 200:, 50:150] |= mySatBit
        self.addDip(0, AMP_HEIGHT - 100, AMP_WIDTH, 100)

        before = self.exposure.mask.array.copy()
        ipIsrFunctions.maskDECamEdgeBleed(self.exposure, saturatedMaskName="MYSAT")

        expected = before.copy()
        expected[AMP_HEIGHT - 113:, 0:AMP_WIDTH] |= mySatBit
        np.testing.assert_array_equal(self.exposure.mask.array, expected)

    def test_partialAmpIsSkipped(self):
        # An exposure covering only the bottom 600 rows of amp B has a
        # qualifying footprint and dip, but the amp is not fully contained
        # so nothing is masked (rather than failing on the sub-image).
        self.addSatBlock(AMP_WIDTH + 50, 0, 100, 200)
        self.addDip(AMP_WIDTH, 0, AMP_WIDTH, 100)
        partial = geom.Box2I(geom.Point2I(AMP_WIDTH, 0), geom.Extent2I(AMP_WIDTH, 600))
        subExposure = self.exposure[partial]

        before = subExposure.mask.array.copy()
        ipIsrFunctions.maskDECamEdgeBleed(subExposure)

        np.testing.assert_array_equal(subExposure.mask.array, before)

    def test_noDetectorRaises(self):
        exposure = afwImage.ExposureF(geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(50, 50)))
        with self.assertRaises(RuntimeError):
            ipIsrFunctions.maskDECamEdgeBleed(exposure)


class DECamEdgeBleedMaskTaskTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.detector = makeDetector()
        self.exposure = afwImage.ExposureF(self.detector.getBBox())
        self.exposure.setDetector(self.detector)
        rng = np.random.default_rng(12345)
        self.exposure.image.array[:] = SKY + rng.normal(0.0, NOISE, self.exposure.image.array.shape)
        self.exposure.variance.array[:] = NOISE**2
        self.satBit = self.exposure.mask.getPlaneBitMask("SAT")

    def test_defaults(self):
        config = DECamEdgeBleedMaskTask.ConfigClass()
        self.assertEqual(config.satMinArea, 10000)
        self.assertEqual(config.approachRows, 20)
        self.assertEqual(config.nSigma, 5.0)
        self.assertEqual(config.nRowsCheck, 20)
        self.assertEqual(config.nRowsCheckShort, 5)
        self.assertEqual(config.minUsablePixelsPerRow, 30)
        self.assertEqual(config.minLowFraction, 0.03)
        self.assertEqual(config.minLowFractionExtent, 0.01)
        self.assertEqual(config.minBlockedRows, 5)
        self.assertEqual(config.marginFraction, 0.125)
        self.assertEqual(config.saturatedMaskName, "SAT")
        config.validate()

    def test_validateExtentFractionNotAboveDipFraction(self):
        config = DECamEdgeBleedMaskTask.ConfigClass()
        config.minLowFractionExtent = 2*config.minLowFraction
        with self.assertRaises(ValueError):
            config.validate()

    def test_validateFractionsInRange(self):
        config = DECamEdgeBleedMaskTask.ConfigClass()
        config.minLowFraction = 1.5
        with self.assertRaises(ValueError):
            config.validate()

    def test_runPassesNewParameters(self):
        # minLowFraction above the 62 per cent dip fraction of a 25-column dip
        # in 40 usable columns must switch masking off through the task.
        self.exposure.mask.array[AMP_HEIGHT - 200:, 50:150] |= self.satBit
        self.exposure.mask.array[AMP_HEIGHT - 100:, 40:AMP_WIDTH] |= self.satBit
        self.exposure.image.array[AMP_HEIGHT - 100:, 0:25] -= DIP
        before = self.exposure.mask.array.copy()

        config = DECamEdgeBleedMaskTask.ConfigClass()
        config.minLowFraction = 0.9
        config.minLowFractionExtent = 0.5
        DECamEdgeBleedMaskTask(config=config).run(self.exposure)
        np.testing.assert_array_equal(self.exposure.mask.array, before)

        DECamEdgeBleedMaskTask().run(self.exposure)
        expected = before.copy()
        expected[AMP_HEIGHT - 113:, 0:AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(self.exposure.mask.array, expected)

    def test_runMasksEdgeBleed(self):
        self.exposure.mask.array[AMP_HEIGHT - 200:, 50:150] |= self.satBit
        self.exposure.image.array[AMP_HEIGHT - 100:, 0:AMP_WIDTH] -= DIP
        before = self.exposure.mask.array.copy()

        task = DECamEdgeBleedMaskTask()
        task.run(self.exposure)

        expected = before.copy()
        expected[AMP_HEIGHT - 113:, 0:AMP_WIDTH] |= self.satBit
        np.testing.assert_array_equal(self.exposure.mask.array, expected)

    def test_retargetIntoIsrTaskConfig(self):
        config = IsrTaskConfig()
        config.doCameraSpecificMasking = True
        config.masking.retarget(DECamEdgeBleedMaskTask)
        config.masking.nSigma = 3.0
        config.validate()
        self.assertEqual(config.masking.nSigma, 3.0)
        self.assertIs(config.masking.target, DECamEdgeBleedMaskTask)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
