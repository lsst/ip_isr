#!/usr/bin/env python

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
import os
import unittest

import lsst.utils.tests as tests
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.ip.isr as ipIsr

class IsrTestCases(unittest.TestCase):
    def testSaturation(self):
        saturation = 1000

        bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Point2I(19,19))
        maskedImage       = afwImage.MaskedImageF(bbox)
        maskedImage.set(100, 0x0, 1)
        exposure = afwImage.ExposureF(maskedImage, None)

        bbox     = afwGeom.Box2I(afwGeom.Point2I(9,5),
                                 afwGeom.Point2I(9,15))
        submi    = afwImage.MaskedImageF(maskedImage, bbox, afwImage.PARENT, False)
        submi.set(saturation, 0x0, 1)

        ipIsr.makeThresholdMask(
            maskedImage = maskedImage,
            threshold = saturation,
            growFootprints = 0,
            maskName = 'SAT',
        )
        ipIsr.interpolateFromMask(
            maskedImage = maskedImage,
            fwhm = 5.0,
            growFootprints = 1,
            maskName = 'SAT',
        )

        mask = maskedImage.getMask()
        bitmaskBad    = mask.getPlaneBitMask('BAD')
        bitmaskSat    = mask.getPlaneBitMask('SAT')
        bitmaskInterp = mask.getPlaneBitMask('INTRP')
        height        = maskedImage.getHeight()
        width         = maskedImage.getWidth()

        for j in range(height):
            for i in range(width):
                # Grown saturation mask; one around the mask at 9
                if i >= 8 and i <= 10:
                    if (i,j) in [(8,4),(8,16),(10,4),(10,16)]:
                        #Should not be saturated or interpolated at all
                        self.assertEqual(mask.get(i,j) & bitmaskInterp, 0)
                        self.assertEqual(mask.get(i,j) & bitmaskSat, 0)
                    elif (j >4 and j < 16) and (i == 8 or i == 10):
                        # Not saturated but interpolated over
                        self.assertEqual(mask.get(i,j) & bitmaskInterp, bitmaskInterp)
                    elif (j == 4 or j == 16):
                        # Interpolated over; bottom/top
                        self.assertEqual(mask.get(i,j) & bitmaskInterp, bitmaskInterp)
                    elif (j > 4 and j < 16 and i == 9):
                        # Both saturated and interpolated over; guts of it
                        self.assertEqual(mask.get(i,j) & bitmaskInterp, bitmaskInterp)
                        self.assertEqual(mask.get(i,j) & bitmaskSat,    bitmaskSat)
                    else:
                        # Neither; above or below the mask
                        self.assertEqual(mask.get(i,j), 0)
                else:
                    self.assertEqual(mask.get(i,j), 0)


def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(IsrTestCases)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
