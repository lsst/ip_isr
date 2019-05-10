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
import lsst.ip.isr.isrMock as isrMock


class IsrMockCases(lsst.utils.tests.TestCase):
    """Test the generation of IsrMock data.
    """
    def setUp(self):
        self.inputExp = isrMock.TrimmedRawMock().run()
        self.mi = self.inputExp.getMaskedImage()

    def test_simple(self):
        """Chain raw and calibration mock data.

        This test should confirm the raw data is generated as expected.
        """

        initialMean = np.median(self.mi.getImage().getArray()[:])
        initialStd = np.std(self.mi.getImage().getArray()[:])

        bias = isrMock.BiasMock().run()
        self.mi.getImage().getArray()[:] = (self.mi.getImage().getArray()[:] -
                                            bias.getMaskedImage().getImage().getArray()[:])
        newMean = np.median(self.mi.getImage().getArray()[:])
        newStd = np.std(self.mi.getImage().getArray()[:])

        self.assertLess(newMean, initialMean)

        initialMean = newMean
        initialStd = newStd

        dark = isrMock.DarkMock().run()
        self.mi.getImage().getArray()[:] = (self.mi.getImage().getArray()[:] -
                                            dark.getMaskedImage().getImage().getArray()[:])
        newMean = np.median(self.mi.getImage().getArray()[:])
        newStd = np.std(self.mi.getImage().getArray()[:])

        self.assertLess(newMean, initialMean)

        initialMean = newMean
        initialStd = newStd

        flat = isrMock.FlatMock().run()
        self.mi.getImage().getArray()[:] = (self.mi.getImage().getArray()[:] -
                                            flat.getMaskedImage().getImage().getArray()[:])
        newMean = np.median(self.mi.getImage().getArray()[:])
        newStd = np.std(self.mi.getImage().getArray()[:])

        self.assertAlmostEqual(newMean, initialMean, -2)
        self.assertLess(newStd, initialStd)

        initialMean = newMean
        initialStd = newStd

        fringe = isrMock.FringeMock().run()
        self.mi.getImage().getArray()[:] = (self.mi.getImage().getArray()[:] -
                                            fringe.getMaskedImage().getImage().getArray()[:])
        newMean = np.median(self.mi.getImage().getArray()[:])
        newStd = np.std(self.mi.getImage().getArray()[:])

        self.assertLess(newMean, initialMean)

    def test_untrimmedSimple(self):
        """Confirm untrimmed data classes are generated consistently.
        """
        exposure = isrMock.RawMock().run()
        fringe = isrMock.UntrimmedFringeMock().run()

        initialStd = np.std(exposure.getMaskedImage().getImage().getArray()[:])

        diff = (exposure.getMaskedImage().getImage().getArray()[:] -
                fringe.getMaskedImage().getImage().getArray()[:])

        newStd = np.std(diff[:])

        self.assertLess(newStd, initialStd)

    def test_productTypes(self):
        """Test non-image data is returned as the expected type.
        """
        self.assertIsInstance(isrMock.BfKernelMock().run(), np.ndarray)
        self.assertIsInstance(isrMock.CrosstalkCoeffMock().run(), np.ndarray)

        self.assertIsInstance(isrMock.DefectMock().run()[0], lsst.meas.algorithms.Defect)
        self.assertIsInstance(isrMock.TransmissionMock().run(), afwImage.TransmissionCurve)

    def test_edgeCases(self):
        """Test that improperly specified configurations do not return data.
        """
        config = isrMock.IsrMockConfig()
        self.assertIsNone(isrMock.IsrMock(config=config).run())

        with self.assertRaises(RuntimeError):
            config.doGenerateData = True
            isrMock.IsrMock(config=config).run()


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
