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
from lsst.ip.isr.binImageDataTask import binImageData


class BinImageDataTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        """Set up the various image and image-like data type inputs,
        plus reference images that represent the expected output from
        an 8x8 binning.
        """
        self.imageF = afwImage.ImageF(2048, 2048)
        self.imageF.set(100)
        self.miF = afwImage.makeMaskedImage(self.imageF)
        self.expF = afwImage.makeExposure(self.miF)

        self.imageI = afwImage.ImageI(2048, 2048)
        self.imageI.set(100)
        self.miI = afwImage.makeMaskedImage(self.imageI)
        self.expI = afwImage.makeExposure(self.miI)

        self.refImageF = afwImage.ImageF(256, 256)
        self.refImageF.set(100)
        self.refMiF = afwImage.makeMaskedImage(self.refImageF)
        self.refExpF = afwImage.makeExposure(self.refMiF)

        self.refImageI = afwImage.ImageI(256, 256)
        self.refImageI.set(100)
        self.refMiI = afwImage.makeMaskedImage(self.refImageI)
        self.refExpI = afwImage.makeExposure(self.refMiI)

    def test_ImageDataBinning(self):
        """Test an 8x8 binning."""
        images = [(self.imageF, self.refImageF), (self.imageI, self.refImageI)]
        maskedImages = [(self.miF, self.refMiF), (self.miI, self.refMiI)]
        exposures = [(self.expF, self.refExpF), (self.expI, self.refExpI)]

        for input, ref in images:
            binnedData = binImageData(input, binFactor=8)
            self.assertImagesEqual(ref, binnedData)
        for input, ref in maskedImages:
            binnedData = binImageData(input, binFactor=8)
            self.assertMaskedImagesEqual(ref, binnedData)
        for input, ref in exposures:
            binnedData = binImageData(input, binFactor=8)
            self.assertMaskedImagesEqual(
                ref.getMaskedImage(), binnedData.getMaskedImage()
            )

    def test_BinImageDataTypes(self):
        """Test that a TypeError is raised should the input
        not be of type `lsst.afw.image.Image`, `lsst.afw.image.MaskedImage`, or
        `lsst.afw.image.Exposure, or one of their sub-types.
        """
        with self.assertRaises(TypeError):
            _ = binImageData(None)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys

    setup_module(sys.modules[__name__])
    unittest.main()
