import unittest

import numpy as np

import lsst.utils.tests as tests
import lsst.afw.image as afwImage
from lsst.ip.isr import applyLookupTable

def referenceApply(image, table, indOffset):
    """Reference implementation of applyLookupTable

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

class ApplyLookupTableTestCase(unittest.TestCase):
    def setUp(self):
        pass
        
    def tearDown(self):
        pass

    def testBasics(self):
        """Test basic functionality of applyLookupTable
        """
        imDim = (100, 85)
        imMin = -5
        imMax = 2500
        tableLen = 2000
        tableSigma = 55
        for indOffset in (0, -50, 234):
            for imageClass in (afwImage.ImageF, afwImage.ImageD):
                inImage = makeRampImage(imageClass, imDim=imDim, minVal=imMin, maxVal=imMax)
                table = np.random.normal(scale=tableSigma, size=tableLen)

                refImage = imageClass(inImage, True)
                refNumBad = referenceApply(image=refImage, table=table, indOffset=indOffset)

                measImage = imageClass(inImage, True)
                measNumBad = applyLookupTable(refImage, table, indOffset)

                self.assertEqual(refNumBad, measNumBad)
                self.assertImagesNearlyEqual(refImage, measImage)

def makeRampImage(imClass, imDim, minVal, maxVal):
    im = imClass(*imDim)
    imArr = im.getArray()
    rampArr = np.linspace(start=minVal, stop=maxVal, num=imDim[0]*imDim[1])
    imArr[:] = np.reshape(rampArr, (imDim[1], imDim[0])) # numpy arrays are transposed w.r.t. afwImage
    return im

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(ApplyLookupTableTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
