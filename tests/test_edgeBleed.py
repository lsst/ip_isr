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

    def test_largeFootprintIsNotMasked(self):
        # 300 x 400 = 120000 SAT pixels: above satMaxArea. Only x 0-49 of the
        # dip rows are unsaturated, still enough low pixels per row to be
        # detectable, so only the area cut prevents masking.
        self.addSatBlock(50, 400, 300, 400)
        self.addDip(0, AMP_HEIGHT - 100, AMP_WIDTH, 100)

        before, after = self.runMasking()

        np.testing.assert_array_equal(after, before)

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

    def test_noDetectorLogsWarningAndReturns(self):
        exposure = afwImage.ExposureF(geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(50, 50)))
        with self.assertLogs("lsst.ip.isr.isrFunctions", level="WARNING") as cm:
            ipIsrFunctions.maskDECamEdgeBleed(exposure)
        self.assertEqual(len(cm.output), 1)
        self.assertIn("no detector", cm.output[0])


class IsrTaskConfigEdgeBleedTestCase(lsst.utils.tests.TestCase):
    def test_defaults(self):
        config = IsrTaskConfig()
        self.assertFalse(config.doDECamEdgeBleedMask)
        self.assertEqual(config.decamEdgeBleedSatMinArea, 10000)
        self.assertEqual(config.decamEdgeBleedSatMaxArea, 100000)
        self.assertEqual(config.decamEdgeBleedApproachRows, 20)
        self.assertEqual(config.decamEdgeBleedNSigma, 5.0)
        self.assertEqual(config.decamEdgeBleedNRowsCheck, 20)
        self.assertEqual(config.decamEdgeBleedMinLowPixelsPerRow, 30)
        self.assertEqual(config.decamEdgeBleedMinLowPixelsExtent, 10)
        self.assertEqual(config.decamEdgeBleedMarginFraction, 0.125)
        config.doDECamEdgeBleedMask = True
        config.validate()

    def test_requiresSaturation(self):
        config = IsrTaskConfig()
        config.doDECamEdgeBleedMask = True
        config.doSaturation = False
        with self.assertRaises(ValueError):
            config.validate()


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
