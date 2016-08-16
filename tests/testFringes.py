#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2012 LSST Corporation.
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
import unittest

import numpy
import lsst.utils.tests
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9
import lsst.afw.image.utils as afwImageUtils
from lsst.ip.isr.fringe import FringeTask

try:
    debug
except NameError:
    debug = False


def checkDebug():
    """Turn on Task debugging if desired"""
    if debug:
        import lsstDebug
        print "Importing debug settings..."

        def DebugInfo(name):
            di = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
            if name == "lsst.ip.isr.fringe":
                di.plot = True
            return di
        lsstDebug.Info = DebugInfo


class FringeDataRef(object):
    """Quacks like a ButlerDataRef, so we can provide an in-memory fringe frame"""

    def __init__(self, fringe):
        self.fringe = fringe
        self.dataId = {'test': True}

    def get(self, name="fringe", immediate=False):
        if name == "fringe":
            return self.fringe
        if name == "ccdExposureId":
            return 1000


def createFringe(width, height, xFreq, xOffset, yFreq, yOffset):
    """Create a fringe frame

    @param width, height    Size of image
    @param xFreq, yFreq     Frequency of sinusoids in x and y
    @param xOffset, yOffset Phase of sinusoids in x and y
    @return Fringe frame
    """
    image = afwImage.ImageF(width, height)
    array = image.getArray()
    x, y = numpy.indices(array.shape)
    array[x, y] = numpy.sin(xFreq*x + xOffset) + numpy.sin(yFreq*y + yOffset)
    mi = afwImage.makeMaskedImage(image)
    exp = afwImage.makeExposure(mi)
    exp.setFilter(afwImage.Filter('FILTER'))
    return exp

frame = 1  # ds9 frame


class FringeTestCase(lsst.utils.tests.TestCase):
    """Tests of the FringeTask"""

    def setUp(self):
        self.size = 512
        self.config = FringeTask.ConfigClass()
        self.config.filters = ['FILTER']
        self.config.num = 5000
        self.config.small = 1
        self.config.large = 128
        self.config.pedestal = False
        afwImageUtils.defineFilter('FILTER', lambdaEff=0)

    def tearDown(self):
        afwImageUtils.resetFilters()

    def checkFringe(self, task, exp, fringes, stddevMax):
        """Run fringe subtraction and verify

        @param task         Task to run
        @param exp          Science exposure
        @param dataRef      Data reference that will provide the fringes
        @param stddevMax    Maximum allowable standard deviation
        """
        if debug:
            global frame
            ds9.mtv(exp, frame=frame, title="Science exposure")
            frame += 1
            if not isinstance(fringes, list):
                fringe = [fringes]
            for i, f in enumerate(fringe):
                ds9.mtv(f, frame=frame, title="Fringe frame %d" % (i+1))
                frame += 1

        task.run(exp, fringes)

        mi = exp.getMaskedImage()

        if debug:
            ds9.mtv(exp, frame=frame, title="Subtracted")
            frame += 1

        mi -= afwMath.makeStatistics(mi, afwMath.MEAN).getValue()
        self.assertLess(afwMath.makeStatistics(mi, afwMath.STDEV).getValue(), stddevMax)

    def testSingle(self, pedestal=0.0, stddevMax=1.0e-4):
        """Test subtraction of a single fringe frame

        @param pedestal    Pedestal to add into fringe frame
        @param stddevMax    Maximum allowable standard deviation
        """
        xFreq = numpy.pi / 10.0
        xOffset = 1.0
        yFreq = numpy.pi / 15.0
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

    def testPedestal(self):
        """Test subtraction of a fringe frame with a pedestal"""
        self.config.pedestal = True
        self.testSingle(pedestal=10000.0, stddevMax=1.0e-3)  # Not sure why this produces worse sttdev
        self.testMultiple(pedestal=10000.0)

    def testMultiple(self, pedestal=0.0):
        """Test subtraction of multiple fringe frames

        @param pedestal    Pedestal to add into fringe frame
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
        """Test the .runDataRef method for complete test converage

        @param pedestal    Pedestal to add into fringe frame
        @param stddevMax   Maximum allowable standard deviation
        """
        xFreq = numpy.pi / 10.0
        xOffset = 1.0
        yFreq = numpy.pi / 15.0
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


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
