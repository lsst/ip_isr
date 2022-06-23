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
import numpy as np

import lsst.utils.tests
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.image as afwImage
import lsst.ip.isr.isrFunctions as isrFunctions
from lsst.ip.isr import BrighterFatterKernel


class BrighterFatterTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        """Set up a no-op BFK dataset
        """
        cameraBuilder = cameraGeom.Camera.Builder('fake camera')
        detectorWrapper = cameraGeom.testUtils.DetectorWrapper(numAmps=4, cameraBuilder=cameraBuilder)
        self.detector = detectorWrapper.detector
        camera = cameraBuilder.finish()

        self.bfk = BrighterFatterKernel(level='AMP', camera=camera, detectorId=1)
        self.bfk.shape = (17, 17)
        self.bfk.badAmps = ['amp 3']

        covar = np.zeros((8, 8))
        covar[0, 0] = 1.0

        kernel = np.zeros(self.bfk.shape)
        kernel[8, 8] = 1.0

        for amp in self.detector:
            ampName = amp.getName()
            if amp in self.bfk.badAmps:
                self.bfk.expIdMask[ampName] = [False, False, False, False, False, False, False, False, False,
                                               False]
            else:
                self.bfk.expIdMask[ampName] = [True, True, True, True, True, True, True, True, False, False]
            self.bfk.rawMeans[ampName] = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
            self.bfk.rawVariances[ampName] = np.array(self.bfk.rawMeans[ampName], dtype=float)
            self.bfk.rawXcorrs[ampName] = [covar for _ in self.bfk.rawMeans[ampName]]
            self.bfk.gain[ampName] = 1.0
            self.bfk.noise[ampName] = 5.0

            self.bfk.meanXcorrs[ampName] = kernel
            self.bfk.valid[ampName] = (ampName != 'amp 3')

            self.bfk.ampKernels[ampName] = kernel

    def test_BrighterFatterInterface(self):
        """Test brighter fatter correction interface using a delta function
        kernel on a flat image"""

        image = afwImage.ImageF(100, 100)
        image.set(100)
        ref_image = afwImage.ImageF(image, True)

        mi = afwImage.makeMaskedImage(image)
        exp = afwImage.makeExposure(mi)

        self.bfk.makeDetectorKernelFromAmpwiseKernels(self.detector.getName())
        kernelToUse = self.bfk.detKernels[self.detector.getName()]

        isrFunctions.brighterFatterCorrection(exp, kernelToUse, 5, 100, False)
        self.assertImagesEqual(ref_image, image)

    def test_BrighterFatterIO(self):
        dictionary = self.bfk.toDict()
        newBfk = BrighterFatterKernel().fromDict(dictionary)
        self.assertEqual(self.bfk, newBfk)

        tables = self.bfk.toTable()
        newBfk = BrighterFatterKernel().fromTable(tables)
        self.assertEqual(self.bfk, newBfk)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
