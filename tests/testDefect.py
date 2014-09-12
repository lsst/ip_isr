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
import unittest

import lsst.utils.tests as tests
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils
import lsst.ip.isr as ipIsr

try:
    type(display)
except NameError:
    display = False

class DefectTestCases(unittest.TestCase):
    def setUp(self):
        self.setVal = 10.

    def tearDown(self):
        del self.setVal

    def testDefectBase(self):
        """Test DefectBases"""

        defectList = measAlg.DefectListT()
        ccdImage = afwImage.MaskedImageF(250,225)  
        ccdImage.set(self.setVal, 0, self.setVal)
        #
        # Insert some defects into the Ccd
        #
        for x0, y0, x1, y1 in [
            (34,  0,   35,  80 ),
            (34,  81,  34,  100),
            (180, 100, 182, 130),
            ]:
            bbox = afwGeom.Box2I(afwGeom.Point2I(x0, y0), afwGeom.Point2I(x1, y1))
            defectList.append(measAlg.Defect(bbox))
            bad = ccdImage.Factory(ccdImage, bbox, afwImage.PARENT)
            bad.set(100)

        ipIsr.maskPixelsFromDefectList(ccdImage, defectList, maskName='BAD')
        mask = ccdImage.getMask()
        bitMask = mask.getPlaneBitMask('BAD')
        for d in defectList:
            bad = mask.Factory(mask, d.getBBox(), afwImage.PARENT)
            self.assertTrue((bad.getArray()&bitMask == bitMask).all())

        if display:
            ds9.mtv(ccdImage.getImage(), title="Defects")
            for d in defectList:
                displayUtils.drawBBox(d.getBBox(), ctype=ds9.CYAN, borderWidth=.5)
            ds9.incrDefaultFrame()

        ipIsr.interpolateDefectList(ccdImage, defectList, 2.)
        im = ccdImage.getImage()
        for d in defectList:
            intrp = im.Factory(im, d.getBBox())
            self.assertTrue((intrp.getArray() == self.setVal).all())

        if display:
            ds9.mtv(ccdImage.getImage(), title="Defects Interpolated")
            for d in defectList:
                displayUtils.drawBBox(d.getBBox(), ctype=ds9.CYAN, borderWidth=.5)
            ds9.incrDefaultFrame()
        
def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(DefectTestCases)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
