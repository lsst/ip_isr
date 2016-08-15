#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import unittest
import pickle
import os
import numpy

import lsst.utils.tests as tests
import lsst.afw.image as afwImage
import lsst.ip.isr as ipIsr


class BrighterFatterTestCases(unittest.TestCase):

    def setUp(self):
        self.filename = "bf_kernel.pkl"
        kernel = afwImage.ImageF(17, 17)
        kernel.set(9, 9, 1)
        kernel.getArray().dump(self.filename)

    def tearDown(self):
        os.unlink(self.filename)

    def testBrighterFatterInterface(self):
        """Test brighter fatter correction interface using a delta function kernel on a flat image"""

        image = afwImage.ImageF(100, 100)
        image.set(100)
        ref_image = afwImage.ImageF(image, True)

        mi = afwImage.makeMaskedImage(image)
        exp = afwImage.makeExposure(mi)

        isrTask = ipIsr.IsrTask()
        with open(self.filename) as f:
            bfKernel = pickle.load(f)

        isrTask.brighterFatterCorrection(exp, bfKernel, 5, 100, False)
        self.assertTrue(numpy.all(ref_image.getArray() == image.getArray()))


def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(BrighterFatterTestCases)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
