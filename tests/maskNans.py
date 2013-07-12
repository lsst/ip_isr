#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2012-2013 LSST Corporation.
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
import sys
import unittest

import numpy
import lsst.utils.tests as tests
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9
from lsst.ip.isr import maskNans

try:
    debug
except NameError:
    debug = False

class MaskNansTestCase(unittest.TestCase):
    def setUp(self):
        self.size = 100
        self.freqImage = 34
        self.freqMask = 12
        self.freqVariance = 23
        self.afterMask = 2
        self.allowMask = 1

    def check(self, ImageClass):
        image = ImageClass(self.size, self.size)
        x,y = numpy.indices((self.size, self.size))
        image.getImage().getArray()[y,x] = numpy.where(x*y % self.freqImage, 0, numpy.nan)
        image.getMask().getArray()[y,x] = numpy.where(x*y % self.freqMask, 0, self.allowMask)
        image.getVariance().getArray()[y,x] = numpy.where(x*y % self.freqVariance, 0, numpy.nan)

        if debug:
            ds9.mtv(image.getImage(), frame=1, title="Image")
            ds9.mtv(image.getVariance(), frame=2, title="Variance")
            ds9.mtv(image.getMask(), frame=3, title="Original mask")

        isUnmasked = numpy.logical_not(image.getMask().getArray() & self.allowMask)
        isNan = numpy.logical_or(numpy.isnan(image.getImage().getArray()),
                                 numpy.isnan(image.getVariance().getArray()))

        maskExpected = numpy.where(numpy.logical_and(isUnmasked, isNan), self.afterMask,
                                   image.getMask().getArray())
        numExpected = numpy.count_nonzero(maskExpected & self.afterMask)

        numNans = maskNans(image, self.afterMask, self.allowMask)

        if debug:
            ds9.mtv(image.getMask(), frame=4, title="UNC-ed mask")

        self.assertEqual(numNans, numExpected)
        self.assertTrue(numpy.all(image.getMask().getArray() == maskExpected))

    def test(self):
        for ImageClass in (afwImage.MaskedImageF, afwImage.MaskedImageD):
            self.check(ImageClass)

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(MaskNansTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--debug":
        debug = True
    run(True)
