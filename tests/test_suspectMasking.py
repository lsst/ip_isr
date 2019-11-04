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

import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as cameraGeom
import lsst.ip.isr as ipIsr


class IsrTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        ampInfo = cameraGeom.Amplifier.Builder()

        ampInfo.setRawBBox(lsst.geom.Box2I(lsst.geom.Point2I(-5, 7), lsst.geom.Extent2I(53, 104)))
        ampInfo.setSuspectLevel(25000)
        self.ampInfo = ampInfo
        self.isrTask = ipIsr.IsrTask()

    def tearDown(self):
        self.ampInfo = None
        self.isrTask = None

    def testBasicMasking(self):
        """Test that masking works
        """
        maxVal = 32000
        fracSuspect = 0.3
        suspectLevel = maxVal*(1 - fracSuspect)

        bbox = self.ampInfo.getRawBBox()
        self.ampInfo.setSuspectLevel(suspectLevel)
        maskedImage = makeRampMaskedImage(bbox, 0, maxVal)
        imArr = maskedImage.getImage().getArray()
        desSetArr = imArr >= suspectLevel
        exposure = afwImage.ExposureF(maskedImage)
        inMaskedImage = maskedImage.Factory(maskedImage, True)  # deep copy
        self.isrTask.suspectDetection(exposure, self.ampInfo)
        maskArr = maskedImage.getMask().getArray()
        suspectMask = maskedImage.getMask().getPlaneBitMask("SUSPECT")
        measSetArr = maskArr == suspectMask
        self.assertImagesEqual(desSetArr, measSetArr)
        self.assertMaskedImagesAlmostEqual(inMaskedImage, maskedImage, doMask=False)

    def testNanLevel(self):
        """Test that setting the suspect level to nan disables masking
        """
        bbox = self.ampInfo.getRawBBox()
        self.ampInfo.setSuspectLevel(float("nan"))
        maskedImage = makeRampMaskedImage(bbox, 0, 32000)
        exposure = afwImage.ExposureF(maskedImage)
        inMaskedImage = maskedImage.Factory(maskedImage, True)  # deep copy
        self.isrTask.suspectDetection(exposure, self.ampInfo)
        self.assertMaskedImagesAlmostEqual(inMaskedImage, maskedImage)

    def testRenamedMasking(self):
        """Test that masking works using some other mask name instead of the default
        """
        AltMaskName = "BAD"  # pick something that exists for simplicity

        isrConfig = ipIsr.IsrTask.ConfigClass()
        isrConfig.suspectMaskName = AltMaskName
        isrTask = ipIsr.IsrTask(config=isrConfig)

        maxVal = 32000
        fracSuspect = 0.3
        suspectLevel = maxVal*(1 - fracSuspect)

        bbox = self.ampInfo.getRawBBox()
        self.ampInfo.setSuspectLevel(suspectLevel)
        maskedImage = makeRampMaskedImage(bbox, 0, maxVal)
        imArr = maskedImage.getImage().getArray()
        desSetArr = imArr >= suspectLevel
        exposure = afwImage.ExposureF(maskedImage)
        inMaskedImage = maskedImage.Factory(maskedImage, True)  # deep copy
        isrTask.suspectDetection(exposure, self.ampInfo)
        maskArr = maskedImage.getMask().getArray()
        suspectMask = maskedImage.getMask().getPlaneBitMask(AltMaskName)
        measSetArr = maskArr == suspectMask
        self.assertImagesEqual(desSetArr, measSetArr)
        self.assertMaskedImagesAlmostEqual(inMaskedImage, maskedImage, doMask=False)


def makeRampMaskedImage(bbox, minVal, maxVal, imgClass=afwImage.MaskedImageF):
    """Make a ramp image of the specified size and image class

    Image values start from 0 at the lower left corner and increase by 1 along rows
    Variance values equal image values + 100
    Mask values equal 0
    """
    mi = imgClass(bbox)
    imageArr = mi.getImage().getArray()
    varianceArr = mi.getVariance().getArray()
    maskArr = mi.getMask().getArray()
    imData = np.linspace(minVal, maxVal, imageArr.size)
    imData.shape = (bbox.getHeight(), bbox.getWidth())
    imageArr[:] = imData
    varianceArr[:] = 100 + imData
    maskArr[:] = 0
    return mi


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
