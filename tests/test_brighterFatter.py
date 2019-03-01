#
# LSST Data Management System
# Copyright 2008-2017 AURA/LSST.
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

import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.ip.isr.isrFunctions as isrFunctions


class BrighterFatterTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        self.filename = "bf_kernel.pkl"
        kernel = afwImage.ImageF(17, 17)
        kernel[9, 9, afwImage.LOCAL] = 1
        kernelPickleString = kernel.getArray().dumps()
        # kernel.getArray().dump(self.filename) triggers an "unclosed file" warning with numpy 1.13.1
        with open(self.filename, 'wb') as f:
            f.write(kernelPickleString)

    def tearDown(self):
        os.unlink(self.filename)

    def testBrighterFatterInterface(self):
        """Test brighter fatter correction interface using a delta function kernel on a flat image"""

        image = afwImage.ImageF(100, 100)
        image.set(100)
        ref_image = afwImage.ImageF(image, True)

        mi = afwImage.makeMaskedImage(image)
        exp = afwImage.makeExposure(mi)

        with open(self.filename, 'rb') as f:
            bfKernel = pickle.load(f)

        isrFunctions.brighterFatterCorrection(exp, bfKernel, 5, 100, False)
        self.assertImagesEqual(ref_image, image)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
