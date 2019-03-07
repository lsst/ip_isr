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

import numpy as np

import lsst.utils.tests
import lsst.afw.image as afwImage
from lsst.ip.isr import maskNans

display = False
if display:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)


class MaskNansTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.size = 100
        self.freqImage = 34
        self.freqMask = 12
        self.freqVariance = 23
        self.afterMask = 2
        self.allowMask = 1

    def check(self, ImageClass):
        image = ImageClass(self.size, self.size)
        x, y = np.indices((self.size, self.size))
        image.getImage().getArray()[y, x] = np.where(x*y % self.freqImage, 0, np.nan)
        image.getMask().getArray()[y, x] = np.where(x*y % self.freqMask, 0, self.allowMask)
        image.getVariance().getArray()[y, x] = np.where(x*y % self.freqVariance, 0, np.nan)

        if display:
            afwDisplay.Display(frame=1).mtv(image.getImage(), title=self._testMethodName + ": Image")
            afwDisplay.Display(frame=2).mtv(image.getVariance(), title=self._testMethodName + ": Variance")
            afwDisplay.Display(frame=3).mtv(image.getMask(), title=self._testMethodName + ": Original mask")

        isUnmasked = np.logical_not(image.getMask().getArray() & self.allowMask)
        isNan = np.logical_or(np.isnan(image.getImage().getArray()),
                              np.isnan(image.getVariance().getArray()))

        maskExpected = np.where(np.logical_and(isUnmasked, isNan), self.afterMask,
                                image.getMask().getArray())
        numExpected = np.count_nonzero(maskExpected & self.afterMask)

        numNans = maskNans(image, self.afterMask, self.allowMask)

        if display:
            afwDisplay.Display(frame=4).mtv(image.getMask(), title=self._testMethodName + ": UNC-ed mask")

        self.assertEqual(numNans, numExpected)
        self.assertMasksEqual(image.getMask(), maskExpected)

    def testMaskImageF(self):
        self.check(afwImage.MaskedImageF)

    def testMaskImageD(self):
        self.check(afwImage.MaskedImageD)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
