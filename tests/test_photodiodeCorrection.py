# This file is part of ip_isr.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest
import tempfile

import lsst.afw.cameraGeom as cameraGeom
import lsst.geom as geom
import lsst.utils.tests

from lsst.ip.isr import PhotodiodeCorrection


class PhotodiodeCorrectionTestCase(lsst.utils.tests.TestCase):
    def setUp(self):

        self.pairs = ['[3021120600588, 3021120600589]',
                      '[3021120600639, 3021120600640]',
                      '[3021120600711, 3021120600712]',
                      '[3021120600774, 3021120600775]',
                      '[3021120600612, 3021120600613]',
                      '[3021120600675, 3021120600676]']

        self.corrections = [0.0008682616635145324, 0.0008634151174990838,
                            0.0008343013803181451, 0.0007904133024002547,
                            0.0007386119393599833, 0.000780670380485576]

        camBuilder = cameraGeom.Camera.Builder("fakeCam")
        detBuilder = camBuilder.add('det_a', 1)
        detBuilder.setSerial("123")

        bbox = geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(100, 100))
        orientation = cameraGeom.Orientation()
        pixelSize = lsst.geom.Extent2D(1, 1)
        detBuilder.setBBox(bbox)
        detBuilder.setOrientation(orientation)
        detBuilder.setPixelSize(pixelSize)

        self.camera = camBuilder
        self.detector = detBuilder

    def testDataOnly(self):
        calib = PhotodiodeCorrection()
        for i, pair in enumerate(self.pairs):
            calib.abscissaCorrections[pair] = self.corrections[i]

        outPath = tempfile.mktemp() + '.yaml'
        calib.writeText(outPath)
        newCorrection = PhotodiodeCorrection().readText(outPath)
        self.assertEqual(calib, newCorrection)

        outPath = tempfile.mktemp() + '.fits'
        calib.writeFits(outPath)
        newCorrection = PhotodiodeCorrection().readFits(outPath)
        self.assertEqual(calib, newCorrection)

    def testFullyPopulated(self):
        calib = PhotodiodeCorrection(detector=self.detector,
                                     camera=self.camera)

        for i, pair in enumerate(self.pairs):
            calib.abscissaCorrections[pair] = self.corrections[i]

        outPath = tempfile.mktemp() + '.yaml'
        calib.writeText(outPath)
        newCorrection = PhotodiodeCorrection().readText(outPath)
        self.assertEqual(calib, newCorrection)

        outPath = tempfile.mktemp() + '.fits'
        calib.writeFits(outPath)
        newCorrection = PhotodiodeCorrection().readFits(outPath)
        self.assertEqual(calib, newCorrection)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
