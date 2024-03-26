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


class IsrMockLSSTCases(lsst.utils.tests.TestCase):
    """Test the generation of IsrMockLSST data.
    """
    def setUp(self):
        self.inputExp = isrMockLSST.CalibratedRawMockLSST().run()
        self.mi = self.inputExp.getMaskedImage()

    def test_simple(self):
        """Check the raw data is generated as expected,
        taking the same approach as in test_isrMock.
        """

        initialMean = np.median(self.mi.getImage().getArray()[:])
        initialStd = np.std(self.mi.getImage().getArray()[:])

        bias = isrMockLSST.BiasMockLSST().run()
        self.mi.getImage().getArray()[:] = (self.mi.getImage().getArray()[:]
                                            - bias.getMaskedImage().getImage().getArray()[:])
        newMean = np.median(self.mi.getImage().getArray()[:])
        newStd = np.std(self.mi.getImage().getArray()[:])

        # here we do a check on the standard deviation, instead of the mean
        # as done in other tests because the bias mock has mean 0
        self.assertLess(newStd, initialStd)

        initialMean = newMean
        initialStd = newStd

        dark = isrMockLSST.DarkMockLSST().run()
        self.mi.getImage().getArray()[:] = (self.mi.getImage().getArray()[:]
                                            - dark.getMaskedImage().getImage().getArray()[:])
        newMean = np.median(self.mi.getImage().getArray()[:])
        newStd = np.std(self.mi.getImage().getArray()[:])

        self.assertLess(newMean, initialMean)

        initialMean = newMean
        initialStd = newStd

        flat = isrMockLSST.FlatMockLSST().run()
        self.mi.getImage().getArray()[:] = (self.mi.getImage().getArray()[:]
                                            - flat.getMaskedImage().getImage().getArray()[:])
        newMean = np.median(self.mi.getImage().getArray()[:])
        newStd = np.std(self.mi.getImage().getArray()[:])

        self.assertAlmostEqual(newMean, initialMean, -2)
        self.assertLess(newStd, initialStd)

        initialMean = newMean
        initialStd = newStd

        fringe = isrMockLSST.FringeMockLSST().run()
        self.mi.getImage().getArray()[:] = (self.mi.getImage().getArray()[:]
                                            - fringe.getMaskedImage().getImage().getArray()[:])
        newMean = np.median(self.mi.getImage().getArray()[:])
        newStd = np.std(self.mi.getImage().getArray()[:])

        self.assertLess(newMean, initialMean)


    def test_untrimmedSimple(self):
        """Confirm untrimmed data classes are generated consistently.
        """
        exposureLowNoise = isrMockLSST.RawMockLSST().run()

        rawMock = isrMockLSST.RawMockLSST()
        rawMock.config.readNoise = 10.
        exposureHighNoise = rawMock.run()


        initialStd = np.std(exposureLowNoise.getMaskedImage().getImage().getArray()[:])

        diff = (exposureHighNoise.getMaskedImage().getImage().getArray()[:]
                - exposureLowNoise.getMaskedImage().getImage().getArray()[:])

        newStd = np.std(diff[:])

        self.assertLess(newStd, initialStd)

    def test_productTypes(self):
        """Tests non-image data are returned as the expected type,
        taking the same approach as in test_isrMock.
        """
        self.assertIsInstance(isrMockLSST.BfKernelMockLSST().run(), np.ndarray)
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


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
