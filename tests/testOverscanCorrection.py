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

import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.pex.policy as pexPolicy
from lsst.ip.isr import Isr
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.isr', Verbosity)

isrDir     = eups.productDir('ip_isr')


# Policy file
InputIsrPolicy = os.path.join(isrDir, 'pipeline', 'isrPolicy.paf')

class IsrTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = pexPolicy.Policy.createPolicy(InputIsrPolicy)
        self.isr = Isr()

    def tearDown(self):
        del self.policy

    def testOverscanCorrectionY(self):
	bbox = afwGeom.Box2I(afwGeom.Point2I(0,0),
			    afwGeom.Point2I(9,12))
        mi = afwImage.MaskedImageF(bbox)
        mi.set(10, 0x0, 1)

        # these should be functionally equivalent
        bbox     = afwGeom.Box2I(afwGeom.Point2I(0,10),
                                 afwGeom.Point2I(9,12))
        biassec  = '[1:10,11:13]'
        overscan = afwImage.MaskedImageF(mi, bbox, afwImage.PARENT)
        overscan.set(2, 0x0, 1)
        
        overscanKeyword = self.policy.getString('overscanPolicy.overscanKeyword')
        fitType = self.policy.getPolicy('overscanPolicy').getString('overscanFitType')

        exposure = afwImage.ExposureF(mi, afwImage.Wcs())
        metadata = exposure.getMetadata()
        metadata.setString(overscanKeyword, biassec)

        self.isr.overscanCorrection(exposure, bbox, fitType)

        height        = mi.getHeight()
        width         = mi.getWidth()
        for j in range(height):
            for i in range(width):
                if j >= 10:
                    self.assertEqual(mi.getImage().get(i,j), 0)
                else:
                    self.assertEqual(mi.getImage().get(i,j), 8)

    def testOverscanCorrectionX(self):
	bbox = afwGeom.Box2I(afwGeom.Point2I(0,0),
			    afwGeom.Point2I(12,9))
        mi = afwImage.MaskedImageF(bbox)
        mi.set(10, 0x0, 1)

        # these should be functionally equivalent
        bbox     = afwGeom.Box2I(afwGeom.Point2I(10,0),
                                 afwGeom.Point2I(12,9))
        biassec  = '[11:13,1:10]'
        overscan = afwImage.MaskedImageF(mi, bbox, afwImage.PARENT)
        overscan.set(2, 0x0, 1)
        
        overscanKeyword = self.policy.getString('overscanPolicy.overscanKeyword')
        fitType = self.policy.getPolicy('overscanPolicy').getString('overscanFitType')
        exposure = afwImage.ExposureF(mi, afwImage.Wcs())
        metadata = exposure.getMetadata()
        metadata.setString(overscanKeyword, biassec)

        self.isr.overscanCorrection(exposure, bbox, fitType)

        height        = mi.getHeight()
        width         = mi.getWidth()
        for j in range(height):
            for i in range(width):
                if i >= 10:
                    self.assertEqual(mi.getImage().get(i,j), 0)
                else:
                    self.assertEqual(mi.getImage().get(i,j), 8)

    def testPolyOverscanCorrectionX(self):
	bbox = afwGeom.Box2I(afwGeom.Point2I(0,0),
			    afwGeom.Point2I(12,9))
        mi = afwImage.MaskedImageF(bbox)
        mi.set(10, 0x0, 1)

        # these should be functionally equivalent
        bbox     = afwGeom.Box2I(afwGeom.Point2I(10,0),
                                 afwGeom.Point2I(12,9))
        biassec  = '[11:13,1:10]'
        overscan = afwImage.MaskedImageF(mi, bbox, afwImage.PARENT)
        overscan.set(2, 0x0, 1)
        for i in range(bbox.getDimensions()[1]):
            for j,off in enumerate([-0.5, 0.0, 0.5]):
                overscan.getImage().set(j,i,2+i+off)
        
        exposure = afwImage.ExposureF(mi, afwImage.Wcs())

        self.isr.overscanCorrection(exposure, bbox, fittype="POLY")

        height        = mi.getHeight()
        width         = mi.getWidth()
        for j in range(height):
            for i in range(width):
                if i == 10:
                    self.assertEqual(mi.getImage().get(i,j), -0.5)
                elif i == 11:
                    self.assertEqual(mi.getImage().get(i,j), 0)
                elif i == 12:
                    self.assertEqual(mi.getImage().get(i,j), 0.5)
                else:
                    self.assertEqual(mi.getImage().get(i,j), 10 - 2 - j)

    def testPolyOverscanCorrectionY(self):
	bbox = afwGeom.Box2I(afwGeom.Point2I(0,0),
			    afwGeom.Point2I(9,12))
        mi = afwImage.MaskedImageF(bbox)
        mi.set(10, 0x0, 1)

        # these should be functionally equivalent
        bbox     = afwGeom.Box2I(afwGeom.Point2I(0,10),
                                 afwGeom.Point2I(9,12))
        biassec  = '[11:13,1:10]'
        overscan = afwImage.MaskedImageF(mi, bbox, afwImage.PARENT)
        overscan.set(2, 0x0, 1)
        for i in range(bbox.getDimensions()[0]):
            for j,off in enumerate([-0.5, 0.0, 0.5]):
                overscan.getImage().set(i,j,2+i+off)
        
        exposure = afwImage.ExposureF(mi, afwImage.Wcs())

        self.isr.overscanCorrection(exposure, bbox, fittype="POLY")

        height        = mi.getHeight()
        width         = mi.getWidth()
        for j in range(height):
            for i in range(width):
                if j == 10:
                    self.assertEqual(mi.getImage().get(i,j), -0.5)
                elif j == 11:
                    self.assertEqual(mi.getImage().get(i,j), 0)
                elif j == 12:
                    self.assertEqual(mi.getImage().get(i,j), 0.5)
                else:
                    self.assertEqual(mi.getImage().get(i,j), 10 - 2 - i)

#####
        
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
