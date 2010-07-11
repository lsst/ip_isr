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
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr
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
        
    def tearDown(self):
        del self.policy

    def doFlat(self, scaling):
        flatScaleKeyword = self.policy.getString('flatPolicy.flatScaleKeyword')
        filenameKeyword  = self.policy.getString('filenameKeyword')
        
        mi = afwImage.MaskedImageF(10,10)
        mi.getImage().set(10)
        exposure = afwImage.ExposureF(mi, afwImage.Wcs())
        
        flat = afwImage.MaskedImageF(10,10)
        flat.getImage().set(1)
        flatexposure = afwImage.ExposureF(flat, afwImage.Wcs())
        dmetadata = flatexposure.getMetadata()
        dmetadata.setString(filenameKeyword, 'Unittest Flat')

        ipIsr.flatCorrection(exposure, flatexposure, 'USER', scaling)

        height        = mi.getHeight()
        width         = mi.getWidth()
        for j in range(height):
            for i in range(width):
                self.assertAlmostEqual(mi.getImage().get(i,j), 10 / (1./scaling), 5)

    def testFlat1(self):
        self.doFlat(scaling=10)

    def testFlat2(self):
        self.doFlat(scaling=0.1)

    def testFlat3(self):
        self.doFlat(scaling=3.7)

    def doIllum(self, scaling):
        flatScaleKeyword = self.policy.getString('flatPolicy.flatScaleKeyword')
        filenameKeyword  = self.policy.getString('filenameKeyword')
        
        mi = afwImage.MaskedImageF(10,10)
        mi.getImage().set(10)
        exposure = afwImage.ExposureF(mi, afwImage.Wcs())
        
        illum = afwImage.MaskedImageF(10,10)
        illum.getImage().set(1)
        illumexposure = afwImage.ExposureF(illum, afwImage.Wcs())
        dmetadata = illumexposure.getMetadata()
        dmetadata.setString(filenameKeyword, 'Unittest Illum')

        ipIsr.illuminationCorrection(exposure, illumexposure, scaling)

        height        = mi.getHeight()
        width         = mi.getWidth()
        for j in range(height):
            for i in range(width):
                self.assertAlmostEqual(mi.getImage().get(i,j), 10 / (1./scaling), 5)

    def testIllum1(self):
        self.doIllum(scaling=10)

    def testIllum2(self):
        self.doIllum(scaling=0.1)

    def testIllum3(self):
        self.doIllum(scaling=3.7)
    
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
