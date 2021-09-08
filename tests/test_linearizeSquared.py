#
# LSST Data Management System
# Copyright 2017 LSST Corporation.
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
import pickle
import logging

import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as cameraGeom
from lsst.afw.geom.testUtils import BoxGrid
from lsst.afw.image.testUtils import makeRampImage
from lsst.ip.isr import Linearizer


def refLinearizeSquared(image, detector):
    """!Basic implementation of squared non-linearization correction

    corr = uncorr + coeff[0]*uncorr^2

    @param[in,out] image  image to correct in place (an lsst.afw.image.Image of some type)
    @param[in] detector  detector info (an lsst.afw.cameraGeom.Detector)
    """
    ampInfoCat = detector.getAmplifiers()
    for ampInfo in ampInfoCat:
        bbox = ampInfo.getBBox()
        sqCoeff = ampInfo.getLinearityCoeffs()[0]
        viewArr = image.Factory(image, bbox).getArray()
        viewArr[:] = viewArr + sqCoeff*viewArr**2


class LinearizeSquaredTestCase(lsst.utils.tests.TestCase):
    """!Unit tests for LinearizeSquared"""

    def setUp(self):
        # the following values are all arbitrary, but sane and varied
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-31, 22), lsst.geom.Extent2I(100, 85))
        self.numAmps = (2, 3)
        self.sqCoeffs = np.array([[0, 5e-6, 2.5e-5], [1e-5, 1.1e-6, 2.1e-5]], dtype=float)
        self.detector = self.makeDetector()

    def tearDown(self):
        # destroy LSST objects so memory test passes
        self.bbox = None
        self.detector = None

    def testBasics(self):
        """!Test basic functionality of LinearizeSquared
        """
        for imageClass in (afwImage.ImageF, afwImage.ImageD):
            inImage = makeRampImage(bbox=self.bbox, start=-5, stop=2500, imageClass=imageClass)

            measImage = inImage.Factory(inImage, True)
            linCorr = Linearizer(detector=self.detector)
            linRes = linCorr.applyLinearity(image=measImage, detector=self.detector)
            desNumLinearized = np.sum(self.sqCoeffs.flatten() > 0)
            self.assertEqual(linRes.numLinearized, desNumLinearized)
            self.assertEqual(linRes.numAmps, len(self.detector.getAmplifiers()))

            refImage = inImage.Factory(inImage, True)
            refLinearizeSquared(image=refImage, detector=self.detector)

            self.assertImagesAlmostEqual(refImage, measImage)

            # make sure logging is accepted
            log = logging.getLogger("ip.isr.LinearizeSquared")
            linRes = linCorr.applyLinearity(image=measImage, detector=self.detector, log=log)

    def testKnown(self):
        """!Test a few known values
        """
        numAmps = (2, 2)
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(4, 4))
        # make a 4x4 image with 4 identical 2x2 subregions that flatten to -1, 0, 1, 2
        im = afwImage.ImageF(bbox)
        imArr = im.getArray()
        imArr[:, :] = np.array(((-1, 0, -1, 0),
                                (1, 2, 1, 2),
                                (-1, 0, -1, 0),
                                (1, 2, 1, 2)), dtype=imArr.dtype)

        sqCoeffs = np.array(((0, 0.11), (-0.15, -12)))
        detector = self.makeDetector(bbox=bbox, numAmps=numAmps, sqCoeffs=sqCoeffs)
        ampInfoCat = detector.getAmplifiers()

        linSq = Linearizer(detector=detector)
        linSq.applyLinearity(im, detector=detector)

        # amp 0 has 0 squared coefficient and so makes no correction
        imArr0 = im.Factory(im, ampInfoCat[0].getBBox()).getArray()
        linCoeff0 = ampInfoCat[0].getLinearityCoeffs()[0]
        self.assertEqual(0, linCoeff0)
        self.assertFloatsAlmostEqual(imArr0.flatten(), (-1, 0, 1, 2))

        # test all amps
        for ampInfo in ampInfoCat:
            imArr = im.Factory(im, ampInfo.getBBox()).getArray()
            linCoeff = ampInfo.getLinearityCoeffs()[0]
            expect = np.array((-1 + linCoeff, 0, 1 + linCoeff, 2 + 4*linCoeff), dtype=imArr.dtype)
            self.assertFloatsAlmostEqual(imArr.flatten(), expect)

    def testPickle(self):
        """!Test that a LinearizeSquared can be pickled and unpickled
        """
        inImage = makeRampImage(bbox=self.bbox, start=-5, stop=2500)
        linSq = Linearizer(detector=self.detector)

        refImage = inImage.Factory(inImage, True)
        refNumOutOfRange = linSq.applyLinearity(refImage, detector=self.detector)

        pickledStr = pickle.dumps(linSq)
        restoredLlt = pickle.loads(pickledStr)

        measImage = inImage.Factory(inImage, True)
        measNumOutOfRange = restoredLlt.applyLinearity(measImage, detector=self.detector)

        self.assertEqual(refNumOutOfRange, measNumOutOfRange)
        self.assertImagesAlmostEqual(refImage, measImage)

    def makeDetector(self, bbox=None, numAmps=None, sqCoeffs=None, linearityType="Squared"):
        """!Make a detector

        @param[in] bbox  bounding box for image
        @param[n] numAmps  x,y number of amplifiers (pair of int)
        @param[in] sqCoeffs  square coefficient for each amplifier (2D array of float)
        @param[in] detName  detector name (a str)
        @param[in] detID  detector ID (an int)
        @param[in] detSerial  detector serial numbe (a str)
        @param[in] linearityType  name of linearity type (a str)

        @return a detector (an lsst.afw.cameraGeom.Detector)
        """
        bbox = bbox if bbox is not None else self.bbox
        numAmps = numAmps if numAmps is not None else self.numAmps
        sqCoeffs = sqCoeffs if sqCoeffs is not None else self.sqCoeffs

        detName = "det_a"
        detId = 1
        detSerial = "123"
        orientation = cameraGeom.Orientation()
        pixelSize = lsst.geom.Extent2D(1, 1)

        camBuilder = cameraGeom.Camera.Builder("fakeCam")
        detBuilder = camBuilder.add(detName, detId)
        detBuilder.setSerial(detSerial)
        detBuilder.setBBox(bbox)
        detBuilder.setOrientation(orientation)
        detBuilder.setPixelSize(pixelSize)

        boxArr = BoxGrid(box=bbox, numColRow=numAmps)
        for i in range(numAmps[0]):
            for j in range(numAmps[1]):
                ampInfo = cameraGeom.Amplifier.Builder()
                ampInfo.setName("amp %d_%d" % (i + 1, j + 1))
                ampInfo.setBBox(boxArr[i, j])
                ampInfo.setLinearityType(linearityType)
                ampInfo.setLinearityCoeffs([sqCoeffs[i, j]])
                detBuilder.append(ampInfo)

        return detBuilder


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
