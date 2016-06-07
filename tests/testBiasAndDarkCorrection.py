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
import lsst.afw.geom as afwGeom
import lsst.ip.isr as ipIsr

class IsrTestCases(unittest.TestCase):
    def setUp(self):
        self.pmin = afwGeom.Point2I(1,1)
        self.pmax = afwGeom.Point2I(10,10)
        self.meanCountsKeyword = "IMMODE"
        self.filenameKeyword = "filename"

    def tearDown(self):
        del self.pmin
        del self.pmax
        del self.meanCountsKeyword
        del self.filenameKeyword

    def testBias(self):
        maskedImage = afwImage.MaskedImageF(afwGeom.Box2I(self.pmin, self.pmax))
        maskedImage.getImage().set(10)

        bias = afwImage.MaskedImageF(afwGeom.Box2I(self.pmin, self.pmax))
        bias.getImage().set(1)
        biasexposure = afwImage.ExposureF(bias, None)
        bmetadata = biasexposure.getMetadata()
        bmetadata.setDouble(self.meanCountsKeyword, 1.)
        bmetadata.setString(self.filenameKeyword, 'Unittest Bias')

        ipIsr.biasCorrection(maskedImage, biasexposure.getMaskedImage())

        height = maskedImage.getHeight()
        width  = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                self.assertEqual(maskedImage.getImage().get(i,j), 9)

    def doDark(self, scaling):
        maskedImage = afwImage.MaskedImageF(afwGeom.Box2I(self.pmin, self.pmax))
        maskedImage.getImage().set(10)

        dark = afwImage.MaskedImageF(afwGeom.Box2I(self.pmin, self.pmax))
        dark.getImage().set(1)
        darkexposure = afwImage.ExposureF(dark, None)
        dmetadata = darkexposure.getMetadata()
        dmetadata.setString(self.filenameKeyword, 'Unittest Dark')

        ipIsr.darkCorrection(maskedImage, darkexposure.getMaskedImage(), 1., scaling)

        height        = maskedImage.getHeight()
        width         = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                self.assertAlmostEqual(maskedImage.getImage().get(i,j), 10 - 1./scaling, 5)

    def testDark1(self):
        self.doDark(scaling=10)

    def testDark2(self):
        self.doDark(scaling=0.1)

    def testDark3(self):
        self.doDark(scaling=3.7)


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
