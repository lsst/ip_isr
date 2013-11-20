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
    def setUp(self):
        self.overscanKeyword = "BIASSEC"

    def tearDown(self):
        del self.overscanKeyword

    def testOverscanCorrectionY(self):
        bbox = afwGeom.Box2I(afwGeom.Point2I(0,0),
                            afwGeom.Point2I(9,12))
        maskedImage = afwImage.MaskedImageF(bbox)
        maskedImage.set(10, 0x0, 1)

        # these should be functionally equivalent
        bbox     = afwGeom.Box2I(afwGeom.Point2I(0,10),
                                 afwGeom.Point2I(9,12))
        biassec  = '[1:10,11:13]'
        overscan = afwImage.MaskedImageF(maskedImage, bbox, afwImage.PARENT)
        overscan.set(2, 0x0, 1)
        
        exposure = afwImage.ExposureF(maskedImage, None)
        metadata = exposure.getMetadata()
        metadata.setString(self.overscanKeyword, biassec)

        ipIsr.overscanCorrection(maskedImage, overscan.getImage(), fitType = "MEDIAN")

        height        = maskedImage.getHeight()
        width         = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                if j >= 10:
                    self.assertEqual(maskedImage.getImage().get(i,j), 0)
                else:
                    self.assertEqual(maskedImage.getImage().get(i,j), 8)

    def testOverscanCorrectionX(self):
        bbox = afwGeom.Box2I(afwGeom.Point2I(0,0),
                            afwGeom.Point2I(12,9))
        maskedImage = afwImage.MaskedImageF(bbox)
        maskedImage.set(10, 0x0, 1)

        # these should be functionally equivalent
        bbox     = afwGeom.Box2I(afwGeom.Point2I(10,0),
                                 afwGeom.Point2I(12,9))
        biassec  = '[11:13,1:10]'
        overscan = afwImage.MaskedImageF(maskedImage, bbox, afwImage.PARENT)
        overscan.set(2, 0x0, 1)
        
        exposure = afwImage.ExposureF(maskedImage, None)
        metadata = exposure.getMetadata()
        metadata.setString(self.overscanKeyword, biassec)

        ipIsr.overscanCorrection(maskedImage, overscan.getImage(), fitType = "MEDIAN")

        height        = maskedImage.getHeight()
        width         = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                if i >= 10:
                    self.assertEqual(maskedImage.getImage().get(i,j), 0)
                else:
                    self.assertEqual(maskedImage.getImage().get(i,j), 8)

    def checkPolyOverscanCorrectionX(self, fitType):
        bbox = afwGeom.Box2I(afwGeom.Point2I(0,0),
                            afwGeom.Point2I(12,9))
        maskedImage = afwImage.MaskedImageF(bbox)
        maskedImage.set(10, 0x0, 1)

        # these should be functionally equivalent
        bbox     = afwGeom.Box2I(afwGeom.Point2I(10,0),
                                 afwGeom.Point2I(12,9))
        biassec  = '[11:13,1:10]'
        overscan = afwImage.MaskedImageF(maskedImage, bbox, afwImage.PARENT)
        overscan.set(2, 0x0, 1)
        for i in range(bbox.getDimensions()[1]):
            for j,off in enumerate([-0.5, 0.0, 0.5]):
                overscan.getImage().set(j,i,2+i+off)
        
        exposure = afwImage.ExposureF(maskedImage, None)

        ipIsr.overscanCorrection(maskedImage, overscan.getImage(), fitType=fitType)

        height        = maskedImage.getHeight()
        width         = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                if i == 10:
                    self.assertEqual(maskedImage.getImage().get(i,j), -0.5)
                elif i == 11:
                    self.assertEqual(maskedImage.getImage().get(i,j), 0)
                elif i == 12:
                    self.assertEqual(maskedImage.getImage().get(i,j), 0.5)
                else:
                    self.assertEqual(maskedImage.getImage().get(i,j), 10 - 2 - j)

    def checkPolyOverscanCorrectionY(self, fitType):
        bbox = afwGeom.Box2I(afwGeom.Point2I(0,0),
                            afwGeom.Point2I(9,12))
        maskedImage = afwImage.MaskedImageF(bbox)
        maskedImage.set(10, 0x0, 1)

        # these should be functionally equivalent
        bbox     = afwGeom.Box2I(afwGeom.Point2I(0,10),
                                 afwGeom.Point2I(9,12))
        biassec  = '[11:13,1:10]'
        overscan = afwImage.MaskedImageF(maskedImage, bbox, afwImage.PARENT)
        overscan.set(2, 0x0, 1)
        for i in range(bbox.getDimensions()[0]):
            for j,off in enumerate([-0.5, 0.0, 0.5]):
                overscan.getImage().set(i,j,2+i+off)
        
        exposure = afwImage.ExposureF(maskedImage, None)

        ipIsr.overscanCorrection(maskedImage, overscan.getImage(), fitType=fitType)

        height        = maskedImage.getHeight()
        width         = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                if j == 10:
                    self.assertEqual(maskedImage.getImage().get(i,j), -0.5)
                elif j == 11:
                    self.assertEqual(maskedImage.getImage().get(i,j), 0)
                elif j == 12:
                    self.assertEqual(maskedImage.getImage().get(i,j), 0.5)
                else:
                    self.assertEqual(maskedImage.getImage().get(i,j), 10 - 2 - i)

    def testPolyOverscanCorrection(self):
        for fitType in ("POLY", "CHEB", "LEG"):
            self.checkPolyOverscanCorrectionX(fitType)
            self.checkPolyOverscanCorrectionY(fitType)

        
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
