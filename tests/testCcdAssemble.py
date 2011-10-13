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
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr
import lsst.pex.logging as logging
import lsst.daf.base as dafBase

import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils

import lsst.afw.display.ds9 as ds9
Verbosity = 4
display = False
logging.Trace_setVerbosity('lsst.ip.isr', Verbosity)

afwDir = eups.productDir('afw')

# Policy file
cameraPolicy = os.path.join(afwDir, 'tests', 'TestCameraGeom.paf')

class IsrTestCases(unittest.TestCase):
    
    def setUp(self):
        self.cameraPolicy = cameraGeomUtils.getGeomPolicy(cameraPolicy)
        afwImage.Filter.reset()
        afwImage.FilterProperty.reset()

        filterPolicy = pexPolicy.Policy()
        filterPolicy.add("lambdaEff", 470.0)
        afwImage.Filter.define(afwImage.FilterProperty("g", filterPolicy))
        
    def tearDown(self):
        del self.cameraPolicy

    def testCcdAssemble(self):
        ccd = cameraGeomUtils.makeCcd(self.cameraPolicy,
                cameraGeom.Id(1234))
        exposureList = []
        for n,a in enumerate(ccd):
            mi = afwImage.MaskedImageF(a.getElectronicAllPixels())
            mi.getImage().set(n)
            mi.getVariance().set(n)
            exp  = afwImage.ExposureF(mi)
            exp.setDetector(a)
            exp.setFilter(afwImage.Filter("g"))
            exp.getCalib().setExptime(15)
            exposureList.append(exp)
        aexp = ipIsr.assembleCcd(exposureList, ccd, reNorm = False,
                isTrimmed=True, inCcs=False)
        xpos = [50,150]
        ypos = [25,76,137,188]
        ind = 0
        for ix in xpos:
            for iy in ypos:
                self.assertEqual(aexp.getMaskedImage().getImage().get(ix,iy), ind)
                ind += 1
        if display:
            ds9.mtv(aexp)

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
