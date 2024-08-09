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

import lsst.afw.image as afwImage
import lsst.ip.isr.isrMockLSST as isrMockLSST
from lsst.ip.isr.brighterFatterKernel import BrighterFatterKernel


class IsrMockLSSTCases(lsst.utils.tests.TestCase):
    """Test the generation of IsrMockLSST data.
    """
    def setUp(self):
        self.inputExp = isrMockLSST.TrimmedRawMockLSST().run()
        self.mi = self.inputExp.maskedImage

        self.defects = isrMockLSST.DefectMockLSST().run()
        self.nonDefectPixels = self.getNonDefectPixels(self.mi.mask)

    def test_simple(self):
        """Check trimmed raw data are generated as expected,
        taking the same approach as in test_isrMock.
        """
        initialMean = np.median(self.mi.image.array[self.nonDefectPixels])
        initialStd = np.std(self.mi.image.array[self.nonDefectPixels])

        # Build and subtract a bias calibration
        biasAdu = isrMockLSST.BiasMockLSST(adu=True).run()
        self.mi.image.array[:] = (self.mi.image.array[:] - biasAdu.image.array[:])
        newMean = np.median(self.mi.image.array[self.nonDefectPixels])
        newStd = np.std(self.mi.image.array[self.nonDefectPixels])

        self.assertFloatsAlmostEqual(newStd, initialStd, rtol=1e-5)

        initialMean = newMean
        initialStd = newStd

        darkAdu = isrMockLSST.DarkMockLSST(adu=True).run()
        self.mi.image.array[:] = (self.mi.image.array[:]
                                  - darkAdu.image.array[:])
        newMean = np.median(self.mi.image.array[self.nonDefectPixels])
        newStd = np.std(self.mi.image.array[self.nonDefectPixels])

        self.assertLess(newMean, initialMean)

        initialMean = newMean
        initialStd = newStd

        flatAdu = isrMockLSST.FlatMockLSST(adu=True).run()
        self.mi.image.array[:] = (self.mi.image.array[:]
                                  - flatAdu.image.array[:])
        newMean = np.median(self.mi.image.array[self.nonDefectPixels])
        newStd = np.std(self.mi.image.array[self.nonDefectPixels])

        self.assertAlmostEqual(newMean, initialMean, -2)
        self.assertLess(newStd, initialStd)

        initialMean = newMean
        initialStd = newStd

        fringeAdu = isrMockLSST.FringeMockLSST(adu=True).run()
        self.mi.image.array[:] = (self.mi.image.array[:]
                                  - fringeAdu.image.array[:])
        newMean = np.median(self.mi.image.array[self.nonDefectPixels])
        newStd = np.std(self.mi.image.array[self.nonDefectPixels])

        self.assertAlmostEqual(newMean, initialMean, -1)
        self.assertAlmostEqual(newStd, initialStd, -1)

    # TODO: add tests that bias is consistent, etc.

    def test_untrimmedSimple(self):
        """Test untrimmed mocks are genetared.
        """
        exposureLowNoise = isrMockLSST.RawMockLSST().run()

        rawMock = isrMockLSST.RawMockLSST()
        rawMock.config.readNoise = 100.
        exposureHighNoise = rawMock.run()

        lowNoiseStd = np.std(exposureLowNoise.image.array[:])
        highNoiseStd = np.std(exposureHighNoise.image.array[:])

        self.assertLess(lowNoiseStd, highNoiseStd)

    def test_productTypes(self):
        """Tests non-image data are returned as the expected type,
        taking the same approach as in test_isrMock.
        """
        self.assertIsInstance(isrMockLSST.BfKernelMockLSST().run(), BrighterFatterKernel)
        self.assertIsInstance(isrMockLSST.CrosstalkCoeffMockLSST().run(), np.ndarray)

        self.assertIsInstance(isrMockLSST.DefectMockLSST().run()[0], lsst.meas.algorithms.Defect)
        self.assertIsInstance(isrMockLSST.TransmissionMockLSST().run(), afwImage.TransmissionCurve)

    def test_edgeCases(self):
        """Tests that improperly specified configurations do not return data,
        taking the same approach as in test_isrMock
        """
        config = isrMockLSST.IsrMockLSSTConfig()
        self.assertIsNone(isrMockLSST.IsrMockLSST(config=config).run())

        with self.assertRaises(RuntimeError):
            config.doGenerateData = True
            isrMockLSST.IsrMockLSST(config=config).run()

    def getNonDefectPixels(self, maskOrigin):
        """Get the non-defect pixels to compare.

        Parameters
        ----------
        maskOrigin : `lsst.afw.image.MaskX`
            The origin mask (for shape and type).

        Returns
        -------
        pix_x, pix_y : `tuple` [`np.ndarray`]
            x and y values of good pixels.
        """
        maskTemp = maskOrigin.clone()
        maskTemp[:, :] = 0

        for defect in self.defects:
            maskTemp[defect.getBBox()] = 1

        return np.where(maskTemp.array == 0)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
