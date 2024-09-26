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

import lsst.utils.tests
import lsst.afw.image as afwImage
from lsst.ip.isr.binExposureTask import binExposure


class BinExposureTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        '''Set up both ExposureF and ExposureI inputs,
        plus reference images that represent the expected
        output from an 8x8 binning.
        '''
        self.imageF = afwImage.ImageF(2048, 2048)
        self.imageF.set(100)
        self.miF = afwImage.makeMaskedImage(self.imageF)
        self.expF = afwImage.makeExposure(self.miF)

        self.imageI = afwImage.ImageI(2048, 2048)
        self.imageI.set(100)
        self.miI = afwImage.makeMaskedImage(self.imageI)
        self.expI = afwImage.makeExposure(self.miI)

        refImageF = afwImage.ImageF(256, 256)
        refImageF.set(100)
        refMiF = afwImage.makeMaskedImage(refImageF)
        self.refExpF = afwImage.makeExposure(refMiF)

        refImageI = afwImage.ImageI(256, 256)
        refImageI.set(100)
        refMiI = afwImage.makeMaskedImage(refImageI)
        self.refExpI = afwImage.makeExposure(refMiI)

    def test_ExposureBinning(self):
        '''Test an 8x8 binning.'''
        binnedExposure = binExposure(self.expF, binFactor=8)
        self.assertMaskedImagesEqual(
            self.refExpF.getMaskedImage(),
            binnedExposure.getMaskedImage(),
        )

        binnedExposure = binExposure(self.expI, binFactor=8)
        self.assertMaskedImagesEqual(
            self.refExpI.getMaskedImage(),
            binnedExposure.getMaskedImage(),
        )

    def test_BinExsposureDataTypes(self):
        '''Test that a TypeError is raised should the input
        not be of type `lsst.afw.image.Exposure or one of its sub-types.
        '''
        with self.assertRaises(TypeError):
            _ = binExposure(None)
        with self.assertRaises(TypeError):
            _ = binExposure(self.imageF)
        with self.assertRaises(TypeError):
            _ = binExposure(self.miF)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
