from future import standard_library
standard_library.install_aliases()
from builtins import range
import unittest
import pickle

import numpy as np

import lsst.utils.tests
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.cameraGeom as cameraGeom
from lsst.afw.geom.testUtils import BoxGrid
from lsst.afw.image.testUtils import makeRampImage
from lsst.ip.isr import LinearizeSquared
from lsst.log import Log


def refLinearizeSquared(image, detector):
    """!Basic implementation of squared non-linearization correction

    corr = uncorr + coeff[0]*uncorr^2

    @param[in,out] image  image to correct in place (an lsst.afw.image.Image of some type)
    @param[in] detector  detector info (an lsst.afw.cameraGeom.Detector)
    """
    ampInfoCat = detector.getAmpInfoCatalog()
    for ampInfo in ampInfoCat:
        bbox = ampInfo.getBBox()
        sqCoeff = ampInfo.getLinearityCoeffs()[0]
        viewArr = image.Factory(image, bbox).getArray()
        viewArr[:] = viewArr + sqCoeff*viewArr**2


class LinearizeSquaredTestCase(lsst.utils.tests.TestCase):
    """!Unit tests for LinearizeSquared"""

    def setUp(self):
        # the following values are all arbitrary, but sane and varied
        self.bbox = afwGeom.Box2I(afwGeom.Point2I(-31, 22), afwGeom.Extent2I(100, 85))
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
            linSq = LinearizeSquared()
            linRes = linSq(image=measImage, detector=self.detector)
            desNumLinearized = np.sum(self.sqCoeffs.flatten() > 0)
            self.assertEqual(linRes.numLinearized, desNumLinearized)
            self.assertEqual(linRes.numAmps, len(self.detector.getAmpInfoCatalog()))

            refImage = inImage.Factory(inImage, True)
            refLinearizeSquared(image=refImage, detector=self.detector)

            self.assertImagesNearlyEqual(refImage, measImage)

            # make sure logging is accepted
            log = Log.getLogger("ip.isr.LinearizeSquared")
            linRes = linSq(image=measImage, detector=self.detector, log=log)

    def testKnown(self):
        """!Test a few known values
        """
        numAmps = (2, 2)
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(4, 4))
        # make a 4x4 image with 4 identical 2x2 subregions that flatten to -1, 0, 1, 2
        im = afwImage.ImageF(bbox)
        imArr = im.getArray()
        imArr[:, :] = np.array(((-1, 0, -1, 0),
                                (1, 2, 1, 2),
                                (-1, 0, -1, 0),
                                (1, 2, 1, 2)), dtype=imArr.dtype)

        sqCoeffs = np.array(((0, 0.11), (-0.15, -12)))
        detector = self.makeDetector(bbox=bbox, numAmps=numAmps, sqCoeffs=sqCoeffs)
        ampInfoCat = detector.getAmpInfoCatalog()

        linSq = LinearizeSquared()
        linSq(im, detector=detector)

        # amp 0 has 0 squared coefficient and so makes no correction
        imArr0 = im.Factory(im, ampInfoCat[0].getBBox()).getArray()
        linCoeff0 = ampInfoCat[0].getLinearityCoeffs()[0]
        self.assertEqual(0, linCoeff0)
        self.assertClose(imArr0.flatten(), (-1, 0, 1, 2))

        # test all amps
        for ampInfo in ampInfoCat:
            imArr = im.Factory(im, ampInfo.getBBox()).getArray()
            linCoeff = ampInfo.getLinearityCoeffs()[0]
            expect = np.array((-1 + linCoeff, 0, 1 + linCoeff, 2 + 4*linCoeff), dtype=imArr.dtype)
            self.assertClose(imArr.flatten(), expect)

    def testPickle(self):
        """!Test that a LinearizeSquared can be pickled and unpickled
        """
        inImage = makeRampImage(bbox=self.bbox, start=-5, stop=2500)
        linSq = LinearizeSquared()

        refImage = inImage.Factory(inImage, True)
        refNumOutOfRange = linSq(refImage, self.detector)

        pickledStr = pickle.dumps(linSq)
        restoredLlt = pickle.loads(pickledStr)

        measImage = inImage.Factory(inImage, True)
        measNumOutOfRange = restoredLlt(measImage, self.detector)

        self.assertEqual(refNumOutOfRange, measNumOutOfRange)
        self.assertImagesNearlyEqual(refImage, measImage)

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
        schema = afwTable.AmpInfoTable.makeMinimalSchema()
        ampInfoCat = afwTable.AmpInfoCatalog(schema)
        boxArr = BoxGrid(box=bbox, numColRow=numAmps)
        for i in range(numAmps[0]):
            for j in range(numAmps[1]):
                ampInfo = ampInfoCat.addNew()
                ampInfo.setName("amp %d_%d" % (i + 1, j + 1))
                ampInfo.setBBox(boxArr[i, j])
                ampInfo.setLinearityType(linearityType)
                ampInfo.setLinearityCoeffs([sqCoeffs[i, j]])
        detName = "det_a"
        detId = 1
        detSerial = "123"
        orientation = cameraGeom.Orientation()
        pixelSize = afwGeom.Extent2D(1, 1)
        transMap = {}
        return cameraGeom.Detector(
            detName,
            detId,
            cameraGeom.SCIENCE,
            detSerial,
            bbox,
            ampInfoCat,
            orientation,
            pixelSize,
            transMap,
        )


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
