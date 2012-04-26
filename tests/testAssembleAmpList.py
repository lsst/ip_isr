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

import eups
import lsst.utils.tests as tests
import lsst.afw.image as afwImage
import lsst.pex.policy as pexPolicy
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
from lsst.ip.isr import AssembleCcdTask

afwDir = eups.productDir('afw')

# Policy file
CameraPolicyPath = os.path.join(afwDir, 'tests', 'TestCameraGeom.paf')

class AssembleCcdTestCase(unittest.TestCase):
    def setUp(self):
        self.cameraPolicy = cameraGeomUtils.getGeomPolicy(CameraPolicyPath)
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
            mi = afwImage.MaskedImageF(a.getDiskAllPixels())
            mi.getImage().set(n)
            mi.getVariance().set(n)
            exp  = afwImage.ExposureF(mi)
            exp.setDetector(a)
            exp.setFilter(afwImage.Filter("g"))
            exp.getCalib().setExptime(15)
            exposureList.append(exp)
        
        assemblerConfig = AssembleCcdTask.ConfigClass()
        assemblerConfig.doRenorm = False
        assembler = AssembleCcdTask(config=assemblerConfig)
        
        assembledExposure = assembler.assembleAmpList(exposureList)
        xpos = [50,150]
        ypos = [25,76,137,188]
        ind = 0
        for ix in xpos:
            for iy in ypos:
                self.assertEqual(assembledExposure.getMaskedImage().getImage().get(ix,iy), ind)
                ind += 1
        
def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(AssembleCcdTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
