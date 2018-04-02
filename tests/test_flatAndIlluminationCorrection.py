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
from __future__ import absolute_import, division, print_function

from builtins import range
import unittest

import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.ip.isr as ipIsr


class IsrTestCases(unittest.TestCase):

    def setUp(self):
        self.pmin = afwGeom.Point2I(1, 1)
        self.pmax = afwGeom.Point2I(10, 10)
        self.flatScaleKeyword = "IMMODE"
        self.filenameKeyword = "filename"

    def tearDown(self):
        del self.pmin
        del self.pmax
        del self.flatScaleKeyword
        del self.filenameKeyword

    def doFlat(self, scaling):
        maskedImage = afwImage.MaskedImageF(afwGeom.Box2I(self.pmin, self.pmax))
        maskedImage.getImage().set(10)

        flat = afwImage.MaskedImageF(afwGeom.Box2I(self.pmin, self.pmax))
        flat.getImage().set(1)
        flatexposure = afwImage.ExposureF(flat, None)
        dmetadata = flatexposure.getMetadata()
        dmetadata.setString(self.filenameKeyword, 'Unittest Flat')

        ipIsr.flatCorrection(maskedImage, flatexposure.getMaskedImage(), 'USER', scaling)

        height = maskedImage.getHeight()
        width = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                self.assertAlmostEqual(maskedImage.getImage().get(i, j), 10 / (1./scaling), 5)

    def testFlat1(self):
        self.doFlat(scaling=10)

    def testFlat2(self):
        self.doFlat(scaling=0.1)

    def testFlat3(self):
        self.doFlat(scaling=3.7)

    def doIllum(self, scaling):
        maskedImage = afwImage.MaskedImageF(afwGeom.Box2I(self.pmin, self.pmax))
        maskedImage.getImage().set(10)

        illum = afwImage.MaskedImageF(afwGeom.Box2I(self.pmin, self.pmax))
        illum.getImage().set(1)
        illumexposure = afwImage.ExposureF(illum, None)
        dmetadata = illumexposure.getMetadata()
        dmetadata.setString(self.filenameKeyword, 'Unittest Illum')

        ipIsr.illuminationCorrection(maskedImage, illumexposure.getMaskedImage(), scaling)

        height = maskedImage.getHeight()
        width = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                self.assertAlmostEqual(maskedImage.getImage().get(i, j), 10 / (1./scaling), 5)

    def testIllum1(self):
        self.doIllum(scaling=10)

    def testIllum2(self):
        self.doIllum(scaling=0.1)

    def testIllum3(self):
        self.doIllum(scaling=3.7)

    def testGainAndReadnoise(self):
        import lsst.afw.image as afwImage
        from lsst.afw.cameraGeom.testUtils import DetectorWrapper
        from lsst.ip.isr import IsrTask

        isrTask = IsrTask()

        detector = DetectorWrapper().detector
        raw = afwImage.ExposureF(detector.getBBox())

        level = 10
        readNoise = 1
        raw.image.set(level)

        amp = detector[0]
        amp.setReadNoise(readNoise)

        for gain in [-1, 0, 0.1, 1]:
            amp.setGain(gain)
            isrTask.updateVariance(raw, amp)
            if gain <= 0:               # behave the same way as amp.setGain
                gain = 1

            self.assertEqual(raw.variance.get(0, 0), level/gain + readNoise**2)

class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
