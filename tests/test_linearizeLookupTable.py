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

import numpy as np

import lsst.utils.tests
import lsst.utils
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.cameraGeom as cameraGeom
from lsst.afw.geom.testUtils import BoxGrid
from lsst.afw.image.testUtils import makeRampImage
from lsst.ip.isr import applyLookupTable, LinearizeLookupTable
from lsst.log import Log


def refLinearize(image, detector, table):
    """!Basic implementation of lookup table based non-linearity correction

    @param[in,out] image  image to correct in place (an lsst.afw.image.Image of some type)
    @param[in] detector  detector info (an lsst.afw.cameraGeom.Detector)
    @param[in] table  lookup table: a 2D array of values of the same type as image;
            - one row for each row index (value of coef[0] in the amp info catalog)
            - one column for each image value

    @return the number of pixels whose values were out of range of the lookup table
    """
    ampInfoCat = detector.getAmpInfoCatalog()
    numOutOfRange = 0
    for ampInfo in ampInfoCat:
        bbox = ampInfo.getBBox()
        rowInd, colIndOffset = ampInfo.getLinearityCoeffs()[0:2]
        rowInd = int(rowInd)
        tableRow = table[rowInd, :]
        imView = image.Factory(image, bbox)
        numOutOfRange += applyLookupTable(imView, tableRow, colIndOffset)
    return numOutOfRange


class LinearizeLookupTableTestCase(lsst.utils.tests.TestCase):
    """!Unit tests for LinearizeLookupTable"""

    def setUp(self):
        # the following values are all arbitrary, but sane and varied
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-31, 22), lsst.geom.Extent2I(100, 85))
        self.numAmps = (2, 3)
        self.colIndOffsets = np.array([[0, -50, 2.5], [37, 1, -3]], dtype=float)
        self.rowInds = np.array([[0, 1, 4], [3, 5, 2]])
        numCols = self.numAmps[0]*self.numAmps[1]
        self.assertLess(np.max(self.rowInds), numCols, "error in test conditions; invalid row index")
        self.detector = self.makeDetector()

    def tearDown(self):
        # destroy LSST objects so memory test passes
        self.bbox = None
        self.detector = None

    def testBasics(self):
        """!Test basic functionality of LinearizeLookupTable
        """
        for imageClass in (afwImage.ImageF, afwImage.ImageD):
            inImage = makeRampImage(bbox=self.bbox, start=-5, stop=250, imageClass=imageClass)
            table = self.makeTable(inImage)

            measImage = inImage.Factory(inImage, True)
            llt = LinearizeLookupTable(table=table, detector=self.detector)
            linRes = llt(measImage, self.detector)

            refImage = inImage.Factory(inImage, True)
            refNumOutOfRange = refLinearize(image=refImage, detector=self.detector, table=table)

            self.assertEqual(linRes.numAmps, len(self.detector.getAmpInfoCatalog()))
            self.assertEqual(linRes.numAmps, linRes.numLinearized)
            self.assertEqual(linRes.numOutOfRange, refNumOutOfRange)
            self.assertImagesAlmostEqual(refImage, measImage)

            # make sure logging is accepted
            log = Log.getLogger("ip.isr.LinearizeLookupTable")
            linRes = llt(image=measImage, detector=self.detector, log=log)

    def testErrorHandling(self):
        """!Test error handling in LinearizeLookupTable
        """
        image = makeRampImage(bbox=self.bbox, start=-5, stop=250)
        table = self.makeTable(image)
        llt = LinearizeLookupTable(table=table, detector=self.detector)

        # bad name
        detBadName = self.makeDetector(detName="bad_detector_name")
        with self.assertRaises(RuntimeError):
            llt(image, detBadName)

        # bad serial
        detBadSerial = self.makeDetector(detSerial="bad_detector_serial")
        with self.assertRaises(RuntimeError):
            llt(image, detBadSerial)

        # bad number of amplifiers
        badNumAmps = (self.numAmps[0]-1, self.numAmps[1])
        detBadNumMaps = self.makeDetector(numAmps=badNumAmps)
        with self.assertRaises(RuntimeError):
            llt(image, detBadNumMaps)

        # bad linearity type
        detBadLinType = self.makeDetector(linearityType="bad_linearity_type")
        with self.assertRaises(RuntimeError):
            llt(image, detBadLinType)

        # wrong dimension
        badTable = table[..., np.newaxis]
        with self.assertRaises(RuntimeError):
            LinearizeLookupTable(table=badTable, detector=self.detector)

        # wrong size
        badTable = np.resize(table, (2, 8))
        with self.assertRaises(RuntimeError):
            LinearizeLookupTable(table=badTable, detector=self.detector)

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

        def castAndReshape(arr):
            arr = np.array(arr, dtype=float)
            arr.shape = numAmps
            return arr

        rowInds = castAndReshape((3, 2, 1, 0))  # avoid the trivial mapping to exercise more of the code
        colIndOffsets = castAndReshape((0, 0, 1, 1))
        detector = self.makeDetector(bbox=bbox, numAmps=numAmps, rowInds=rowInds, colIndOffsets=colIndOffsets)
        ampInfoCat = detector.getAmpInfoCatalog()

        # note: table rows are reversed relative to amplifier order because rowInds is a descending ramp
        table = np.array(((7, 6, 5, 4), (1, 1, 1, 1), (5, 4, 3, 2), (0, 0, 0, 0)), dtype=imArr.dtype)

        llt = LinearizeLookupTable(table=table, detector=detector)

        lltRes = llt(image=im, detector=detector)
        self.assertEqual(lltRes.numOutOfRange, 2)

        # amp 0 is a constant correction of 0; one image value is out of range, but it doesn't matter
        imArr0 = im.Factory(im, ampInfoCat[0].getBBox()).getArray()
        self.assertFloatsAlmostEqual(imArr0.flatten(), (-1, 0, 1, 2))

        # amp 1 is a correction of (5, 4, 3, 2), but the first image value is under range
        imArr1 = im.Factory(im, ampInfoCat[1].getBBox()).getArray()
        self.assertFloatsAlmostEqual(imArr1.flatten(), (4, 5, 5, 5))

        # amp 2 is a constant correction of +1; all image values are in range, but it doesn't matter
        imArr2 = im.Factory(im, ampInfoCat[2].getBBox()).getArray()
        self.assertFloatsAlmostEqual(imArr2.flatten(), (0, 1, 2, 3))

        # amp 3 is a correction of (7, 6, 5, 4); all image values in range
        imArr1 = im.Factory(im, ampInfoCat[3].getBBox()).getArray()
        self.assertFloatsAlmostEqual(imArr1.flatten(), (6, 6, 6, 6))

    def testPickle(self):
        """!Test that a LinearizeLookupTable can be pickled and unpickled
        """
        inImage = makeRampImage(bbox=self.bbox, start=-5, stop=2500)
        table = self.makeTable(inImage)
        llt = LinearizeLookupTable(table=table, detector=self.detector)

        refImage = inImage.Factory(inImage, True)
        refNumOutOfRange = llt(refImage, self.detector)

        pickledStr = pickle.dumps(llt)
        restoredLlt = pickle.loads(pickledStr)

        measImage = inImage.Factory(inImage, True)
        measNumOutOfRange = restoredLlt(measImage, self.detector)

        self.assertEqual(refNumOutOfRange, measNumOutOfRange)
        self.assertImagesAlmostEqual(refImage, measImage)

    def makeDetector(self, bbox=None, numAmps=None, rowInds=None, colIndOffsets=None,
                     detName="det_a", detSerial="123", linearityType="LookupTable"):
        """!Make a detector

        @param[in] bbox  bounding box for image
        @param[n] numAmps  x,y number of amplifiers (pair of int)
        @param[in] rowInds  index of lookup table for each amplifier (array of shape numAmps)
        @param[in] colIndOffsets  column index offset for each amplifier (array of shape numAmps)
        @param[in] detName  detector name (a str)
        @param[in] detSerial  detector serial numbe (a str)
        @param[in] linearityType  name of linearity type (a str)

        @return a detector (an lsst.afw.cameraGeom.Detector)
        """
        bbox = bbox if bbox is not None else self.bbox
        numAmps = numAmps if numAmps is not None else self.numAmps
        rowInds = rowInds if rowInds is not None else self.rowInds
        colIndOffsets = colIndOffsets if colIndOffsets is not None else self.colIndOffsets

        schema = afwTable.AmpInfoTable.makeMinimalSchema()
        ampInfoCat = afwTable.AmpInfoCatalog(schema)
        boxArr = BoxGrid(box=bbox, numColRow=numAmps)
        for i in range(numAmps[0]):
            for j in range(numAmps[1]):
                ampInfo = ampInfoCat.addNew()
                ampInfo.setName("amp %d_%d" % (i + 1, j + 1))
                ampInfo.setBBox(boxArr[i, j])
                ampInfo.setLinearityType(linearityType)
                # setLinearityCoeffs is picky about getting a mixed int/float list.
                ampInfo.setLinearityCoeffs(np.array([rowInds[i, j], colIndOffsets[i, j], 0, 0], dtype=float))
        detId = 1
        orientation = cameraGeom.Orientation()
        pixelSize = lsst.geom.Extent2D(1, 1)
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

    def makeTable(self, image, numCols=None, numRows=2500, sigma=55):
        """!Make a 2D lookup table

        @param[in] image  image whose type is used for the table
        @param[in] numCols  number of columns for table; defaults to self.numCols
        @param[in] numRows  number of rows for the table
        @param[in] sigma  standard deviation of normal distribution
        """
        numCols = numCols or self.numAmps[0]*self.numAmps[1]
        dtype = image.getArray().dtype
        table = np.random.normal(scale=sigma, size=(numCols, numRows))
        return np.array(table, dtype=dtype)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
