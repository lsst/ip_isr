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

import numpy as np

import lsst.afw.cameraGeom as cameraGeom
import lsst.geom as geom
import lsst.utils.tests

from lsst.ip.isr import PhotodiodeCalib


class PhotodiodeTestCase(lsst.utils.tests.TestCase):
    def setUp(self):

        self.timeSeries = [0.0, 0.08496094, 0.1689453, 0.2529297,
                           0.3378906, 0.421875, 0.5068359, 0.5908203, 0.6757813,
                           0.7597656, 0.84375, 0.9287109, 1.012695, 1.097656, 1.181641,
                           1.266602]
        self.currentSeries = [-1.350031E-12, 1.598721E-13, -1.456613E-13,
                              -1.385558E-13, -3.517187E-13, -3.375078E-13, 3.597549E-10,
                              3.591616E-10, 3.599823E-10, 3.597158E-10, 3.606893E-10,
                              3.602736E-10, 3.582272E-10, 3.59293E-10, 3.602878E-10,
                              3.588703E-10]

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
        calib = PhotodiodeCalib(timeSamples=self.timeSeries,
                                currentSamples=self.currentSeries)

        self.assertFloatsAlmostEqual(calib.integrate(), 2.88414e-10, rtol=1e-14)
        self.assertFloatsAlmostEqual(calib.integrateDirectSum(), 2.88414e-10, rtol=1e-14)
        self.assertFloatsAlmostEqual(calib.integrateTrimmedSum(), 2.88871e-10, rtol=1e-14)
        self.assertFloatsAlmostEqual(calib.integrateChargeSum(), 3.5960649e-09, rtol=1e-14)

        self.assertEqual(calib.timeSamples.shape, (16, ))
        self.assertEqual(calib.currentSamples.shape, (16, ))

        outPath = tempfile.mktemp() + '.yaml'
        calib.writeText(outPath)
        newPhotodiode = PhotodiodeCalib().readText(outPath)
        self.assertEqual(calib, newPhotodiode)

        outPath = tempfile.mktemp() + '.fits'
        calib.writeFits(outPath)
        newPhotodiode = PhotodiodeCalib().readFits(outPath)
        self.assertEqual(calib, newPhotodiode)

    def testFullyPopulated(self):
        calib = PhotodiodeCalib(timeSamples=self.timeSeries,
                                currentSamples=self.currentSeries,
                                detector=self.detector,
                                camera=self.camera)

        self.assertFloatsAlmostEqual(calib.integrate(), 2.88414e-10, rtol=1e-14)
        self.assertFloatsAlmostEqual(calib.integrateDirectSum(), 2.88414e-10, rtol=1e-14)
        self.assertFloatsAlmostEqual(calib.integrateTrimmedSum(), 2.88871e-10, rtol=1e-14)
        self.assertFloatsAlmostEqual(calib.integrateChargeSum(), 3.5960649e-09, rtol=1e-14)

        self.assertEqual(calib.timeSamples.shape, (16, ))
        self.assertEqual(calib.currentSamples.shape, (16, ))

        outPath = tempfile.mktemp() + '.yaml'
        calib.writeText(outPath)
        newPhotodiode = PhotodiodeCalib().readText(outPath)
        self.assertEqual(calib, newPhotodiode)

        outPath = tempfile.mktemp() + '.fits'
        calib.writeFits(outPath)
        newPhotodiode = PhotodiodeCalib().readFits(outPath)
        self.assertEqual(calib, newPhotodiode)

    def testIntegrateChargeSum(self):
        """
        This tests the .integrateChargeSum method in the presence of
        an above-baseline point which is still below the initial
        threshold for identifying points to exclude in the baseline
        current estimate.
        """
        # Compute a value to insert at time 0.47.
        test_time = 0.47
        test_current = ((max(self.currentSeries) - min(self.currentSeries))/40.
                        + min(self.currentSeries))
        new_times = np.array(self.timeSeries + [test_time])
        new_currents = np.array(self.currentSeries + [test_current])
        # Put these arrays in time order, though this isn't strictly needed.
        index = np.argsort(new_times)
        calib = PhotodiodeCalib(timeSamples=new_times[index],
                                currentSamples=new_currents[index])
        self.assertFloatsAlmostEqual(calib.integrateChargeSum(), 3.6049277e-09, rtol=1e-14)

    def testIntegrateMean(self):
        rng = np.random.RandomState(12345)

        timeSamples = np.linspace(0, 5, 100)
        currentSamples = rng.normal(loc=-1e-9, scale=1e-11, size=len(timeSamples))
        calib = PhotodiodeCalib(timeSamples=timeSamples, currentSamples=currentSamples)
        calib.currentScale = -1.0
        expTime = 5.0
        expected = 1e-9 * expTime

        self.assertFloatsAlmostEqual(calib.integrateMean(expTime), expected, rtol=1e-3)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
