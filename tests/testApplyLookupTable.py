#
# LSST Data Management System
# Copyright 2017 LSST Corporation.
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
import unittest

import numpy as np

import lsst.utils.tests
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.afw.image.testUtils import makeRampImage
from lsst.ip.isr import applyLookupTable


def referenceApply(image, table, indOffset):
    """!Reference implementation of applyLookupTable

    The algorithm is as follows:
        numOutOfRange = 0
        For each i,j of the image:
            lookupInd = int(indOffset + image[i,j])
            if lookupInd not in range [0, table.size() - 1]:
                set lookupInd to nearest edge and increment numOutOfRange
            image[i,j] += table[lookupInd]
        return numOutOfRange

    @param[in,out] image  image to which to add the values; modified in place
    @param[in] table  lookup table
    @param[in] indOffset  scalar added to image value before truncating to lookup column

    @return the number of pixels whose values were out of range
    """
    imArr = image.getArray()
    indArr = np.array(imArr + indOffset, dtype=int)
    maxInd = len(table) - 1
    numBadPoints = np.sum(indArr < 0)
    numBadPoints += np.sum(indArr > maxInd)
    indArr = np.where(indArr < 0, 0, indArr)
    indArr = np.where(indArr >= maxInd, maxInd, indArr)
    imArr += table[indArr]
    return numBadPoints


class ApplyLookupTableTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        np.random.seed(42)

    def testBasics(self):
        """!Test basic functionality of applyLookupTable
        """
        bbox = afwGeom.Box2I(afwGeom.Point2I(-31, 22), afwGeom.Extent2I(100, 85))
        imMin = -5
        imMax = 2500
        tableLen = 2000
        tableSigma = 55
        for indOffset in (0, -50, 234):
            for imageClass in (afwImage.ImageF, afwImage.ImageD):
                inImage = makeRampImage(bbox=bbox, start=imMin, stop=imMax, imageClass=imageClass)
                table = np.random.normal(scale=tableSigma, size=tableLen)
                table = np.array(table, dtype=inImage.getArray().dtype)

                refImage = imageClass(inImage, True)
                refNumBad = referenceApply(image=refImage, table=table, indOffset=indOffset)

                measImage = imageClass(inImage, True)
                measNumBad = applyLookupTable(measImage, table, indOffset)

                self.assertEqual(refNumBad, measNumBad)
                self.assertImagesAlmostEqual(refImage, measImage)

    def testKnown(self):
        """Test that a given image and lookup table produce the known answer

        Apply a negative ramp table to a positive ramp image to get a flat image,
        but have one value out of range at each end, to offset each end point by one
        """
        # generate a small ramp image with ascending integer values
        # starting at some small negative value going positive
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(3, 4))
        numPix = bbox.getWidth()*bbox.getHeight()
        start = -3
        stop = start + numPix - 1
        im = makeRampImage(bbox=bbox, start=start, stop=stop, imageClass=afwImage.ImageF)
        # generate a ramp lookup table with descending integer values,
        # with a range offset by a small arbitrary value from the image ramp
        # make it two elements too short so we can have one value out of range at each end
        numOutOfRangePerEnd = 1
        numOutOfRange = 2*numOutOfRangePerEnd
        tableOffset = -2
        table = np.linspace(
            start=stop + tableOffset,
            stop=numOutOfRange + start + tableOffset,
            num=numPix - numOutOfRange)
        table = np.array(table, dtype=im.getArray().dtype)
        # apply the table with the first and last image value out of range by one
        indOffset = -(start + numOutOfRangePerEnd)
        measNumOutOfRange = applyLookupTable(im, table, indOffset)
        self.assertEqual(numOutOfRange, measNumOutOfRange)
        # at this point the image should all have the same value
        # except the first point will be one less and the last one more
        imArr = im.getArray()
        desVal = start + numOutOfRangePerEnd + table[0]
        desImArr = np.zeros(numPix, dtype=im.getArray().dtype)
        desImArr[:] = desVal
        desImArr[0] -= 1
        desImArr[-1] += 1
        desImArr.shape = imArr.shape
        self.assertClose(desImArr, imArr)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
