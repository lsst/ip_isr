#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import unittest

import lsst.utils.tests
import lsst.geom
import lsst.afw.image as afwImage
import lsst.ip.isr as ipIsr


class IsrTestCases(lsst.utils.tests.TestCase):

    def testSaturation(self):
        """Test saturation threshold masking and interpolation.

        The test image used here is a simulated 20x20 square with a
        10-pixel long defect in the y-direction.
        """
        saturation = 1000

        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(19, 19))
        maskedImage = afwImage.MaskedImageF(bbox)
        maskedImage.set(100, 0x0, 1)

        bbox = lsst.geom.Box2I(lsst.geom.Point2I(9, 5),
                               lsst.geom.Point2I(9, 15))
        submi = afwImage.MaskedImageF(maskedImage, bbox, afwImage.PARENT, False)
        submi.set(saturation, 0x0, 1)

        ipIsr.makeThresholdMask(
            maskedImage=maskedImage,
            threshold=saturation,
            growFootprints=0,
            maskName='SAT',
        )
        ipIsr.interpolateFromMask(
            maskedImage=maskedImage,
            fwhm=5.0,
            growSaturatedFootprints=1,
            maskNameList=['SAT'],
        )

        mask = maskedImage.getMask()
        bitmaskSat = mask.getPlaneBitMask('SAT')
        bitmaskInterp = mask.getPlaneBitMask('INTRP')
        height = maskedImage.getHeight()
        width = maskedImage.getWidth()

        for j in range(height):
            for i in range(width):
                # Grown saturation mask; one around the mask at 9
                if i >= 8 and i <= 10:
                    if (i, j) in [(8, 4), (8, 16), (10, 4), (10, 16)]:
                        # Should not be saturated or interpolated at all
                        self.assertEqual(mask[i, j, afwImage.LOCAL] & bitmaskInterp, 0)
                        self.assertEqual(mask[i, j, afwImage.LOCAL] & bitmaskSat, 0)
                    elif (j > 4 and j < 16) and (i == 8 or i == 10):
                        # Not saturated but interpolated over
                        self.assertEqual(mask[i, j, afwImage.LOCAL] & bitmaskInterp, bitmaskInterp)
                    elif (j == 4 or j == 16):
                        # Interpolated over; bottom/top
                        self.assertEqual(mask[i, j, afwImage.LOCAL] & bitmaskInterp, bitmaskInterp)
                    elif (j > 4 and j < 16 and i == 9):
                        # Both saturated and interpolated over; guts of it
                        self.assertEqual(mask[i, j, afwImage.LOCAL] & bitmaskInterp, bitmaskInterp)
                        self.assertEqual(mask[i, j, afwImage.LOCAL] & bitmaskSat, bitmaskSat)
                    else:
                        # Neither; above or below the mask
                        self.assertEqual(mask[i, j, afwImage.LOCAL], 0)
                else:
                    self.assertEqual(mask[i, j, afwImage.LOCAL], 0)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
