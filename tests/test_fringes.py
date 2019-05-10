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
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.pipe.base as pipeBase
from lsst.ip.isr.fringe import FringeTask

import lsst.ip.isr.isrMock as isrMock

try:
    display
except NameError:
    display = False
else:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)


class FringeDataRef(object):
    """Quacks like a ButlerDataRef, so we can provide an in-memory fringe frame.
    """
    def __init__(self, fringe):
        self.fringe = fringe
        self.dataId = {'test': True}

    def get(self, name="fringe", immediate=False):
        if name == "fringe":
            return self.fringe
        if name == "ccdExposureId":
            return 1000


def createFringe(width, height, xFreq, xOffset, yFreq, yOffset):
    """Create a fringe frame.

    Parameters
    ----------
    width, height : `int`
       Size of image.
    xFreq, yFreq : `float`
       Frequency of sinusoids in x and y.
    xOffset, yOffset : `float`
       Phase of sinusoids in x and y.

    Returns
    -------
    exp : `lsst.afw.image.ExposureF`
       Fringe frame.
    """
    image = afwImage.ImageF(width, height)
    array = image.getArray()
    x, y = np.indices(array.shape)
    array[x, y] = np.sin(xFreq*x + xOffset) + np.sin(yFreq*y + yOffset)
    mi = afwImage.makeMaskedImage(image)
    exp = afwImage.makeExposure(mi)
    exp.setFilter(afwImage.Filter('FILTER'))
    return exp


class FringeTestCase(lsst.utils.tests.TestCase):
    """Tests of the FringeTask.
    """
    def setUp(self):
        self.size = 512
        self.config = FringeTask.ConfigClass()
        self.config.filters = ['FILTER']
        self.config.num = 5000
        self.config.small = 1
        self.config.large = 128
        self.config.pedestal = False
        self.config.iterations = 10
        afwImageUtils.defineFilter('FILTER', lambdaEff=0)

    def tearDown(self):
        afwImageUtils.resetFilters()

    def checkFringe(self, task, exp, fringes, stddevMax):
        """Run fringe subtraction and verify.

        Parameters
        ----------
        task : `lsst.ip.isr.fringe.FringeTask`
           Task to run.
        exp : `lsst.afw.image.ExposureF`
           Science exposure.
        fringes : `list` of `lsst.afw.image.ExposureF`
           Data reference that will provide the fringes.
        stddevMax : `float`
           Maximum allowable standard deviation.
        """
        if display:
            frame = 0
            afwDisplay.Display(frame=frame).mtv(exp, title=self._testMethodName + ": Science exposure")
            frame += 1
            if not isinstance(fringes, list):
                fringe = [fringes]
            else:
                fringe = fringes
            for i, f in enumerate(fringe):
                afwDisplay.Display(frame=frame).mtv(f, title=self._testMethodName +
                                                    ": Fringe frame %d" % (i + 1))
                frame += 1

        task.run(exp, fringes)

        mi = exp.getMaskedImage()

        if display:
            afwDisplay.Display(frame=frame).mtv(exp, title=self._testMethodName + ": Subtracted")
            frame += 1

        mi -= afwMath.makeStatistics(mi, afwMath.MEAN).getValue()
        self.assertLess(afwMath.makeStatistics(mi, afwMath.STDEV).getValue(), stddevMax)

    def testSingle(self, pedestal=0.0, stddevMax=1.0e-4):
        """Test subtraction of a single fringe frame.

        Parameters
        ----------
        pedestal : `float`, optional
           Pedestal to add into fringe frame
        stddevMax : `float`, optional
           Maximum allowable standard deviation.
        """
        xFreq = np.pi/10.0
        xOffset = 1.0
        yFreq = np.pi/15.0
        yOffset = 0.5
        scale = 1.0
        fringe = createFringe(self.size, self.size, xFreq, xOffset, yFreq, yOffset)
        fMi = fringe.getMaskedImage()
        fMi += pedestal
        exp = createFringe(self.size, self.size, xFreq, xOffset, yFreq, yOffset)
        eMi = exp.getMaskedImage()
        eMi *= scale

        task = FringeTask(name="fringe", config=self.config)
        self.checkFringe(task, exp, fringe, stddevMax)

    def testBad(self, bad="BAD"):
        """Test fringe subtraction with bad inputs.

        Parameters
        ----------
        bad : `str`, optional
           Mask plane to use.
        """
        xFreq = np.pi/10.0
        xOffset = 1.0
        yFreq = np.pi/15.0
        yOffset = 0.5
        fringe = createFringe(self.size, self.size, xFreq, xOffset, yFreq, yOffset)
        exp = createFringe(self.size, self.size, xFreq, xOffset, yFreq, yOffset)

        # This is a bad CCD: entirely masked
        exp.maskedImage.image.set(0.0)
        mask = exp.maskedImage.mask
        mask.set(mask.getPlaneBitMask(bad))

        self.config.stats.badMaskPlanes = [bad]
        task = FringeTask(name="fringe", config=self.config)
        task.run(exp, fringe)
        self.assertFloatsEqual(exp.maskedImage.image.array, 0.0)

    def testPedestal(self):
        """Test subtraction of a fringe frame with a pedestal.
        """
        self.config.pedestal = True
        self.testSingle(pedestal=10000.0, stddevMax=1.0e-3)  # Not sure why this produces worse sttdev
        self.testMultiple(pedestal=10000.0)

    def testMultiple(self, pedestal=0.0):
        """Test subtraction of multiple fringe frames

        Paramters
        ---------
        pedestal : `float`, optional
           Pedestal to add into fringe frame.
        """
        xFreqList = [0.1, 0.13, 0.06]
        xOffsetList = [0.0, 0.1, 0.2]
        yFreqList = [0.09, 0.12, 0.07]
        yOffsetList = [0.3, 0.2, 0.1]
        fringeList = [createFringe(self.size, self.size, xFreq, xOffset, yFreq, yOffset)
                      for xFreq, xOffset, yFreq, yOffset in
                      zip(xFreqList, xOffsetList, yFreqList, yOffsetList)]

        for fringe in fringeList:
            fMi = fringe.getMaskedImage()
            fMi += pedestal
        # Generate science frame
        scales = [0.33, 0.33, 0.33]
        image = afwImage.ImageF(self.size, self.size)
        image.set(0)
        for s, f in zip(scales, fringeList):
            image.scaledPlus(s, f.getMaskedImage().getImage())
        mi = afwImage.makeMaskedImage(image)
        exp = afwImage.makeExposure(mi)
        exp.setFilter(afwImage.Filter('FILTER'))

        task = FringeTask(name="multiFringe", config=self.config)
        self.checkFringe(task, exp, fringeList, stddevMax=1.0e-2)

    def testRunDataRef(self, pedestal=0.0, stddevMax=1.0e-4):
        """Test the .runDataRef method for complete test converage.

        Paramters
        ---------
        pedestal : `float`, optional
           Pedestal to add into fringe frame.
        stddevMax : `float`, optional
           Maximum allowable standard deviation.
        """
        xFreq = np.pi/10.0
        xOffset = 1.0
        yFreq = np.pi/15.0
        yOffset = 0.5
        scale = 1.0
        fringe = createFringe(self.size, self.size, xFreq, xOffset, yFreq, yOffset)
        fMi = fringe.getMaskedImage()
        fMi += pedestal
        exp = createFringe(self.size, self.size, xFreq, xOffset, yFreq, yOffset)
        eMi = exp.getMaskedImage()
        eMi *= scale

        task = FringeTask(name="fringe", config=self.config)
        dataRef = FringeDataRef(fringe)
        task.runDataRef(exp, dataRef)

        mi = exp.getMaskedImage()
        mi -= afwMath.makeStatistics(mi, afwMath.MEAN).getValue()
        self.assertLess(afwMath.makeStatistics(mi, afwMath.STDEV).getValue(), stddevMax)

    def test_readFringes(self):
        """Test that fringes can be successfully accessed from the butler.
        """
        task = FringeTask()
        dataRef = isrMock.DataRefMock()

        result = task.readFringes(dataRef, assembler=None)
        self.assertIsInstance(result, pipeBase.Struct)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
