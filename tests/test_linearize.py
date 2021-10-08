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

import logging
import numpy as np

import lsst.utils.tests
import lsst.utils
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.cameraGeom as cameraGeom
from lsst.afw.geom.testUtils import BoxGrid
from lsst.afw.image.testUtils import makeRampImage
from lsst.ip.isr import applyLookupTable, Linearizer


def referenceImage(image, detector, linearityType, inputData, table=None):
    """Generate a reference linearization.

    Parameters
    ----------
    image: `lsst.afw.image.Image`
        Image to linearize.
    detector: `lsst.afw.cameraGeom.Detector`
        Detector this image is from.
    linearityType: `str`
        Type of linearity to apply.
    inputData: `numpy.array`
        An array of values for the linearity correction.
    table: `numpy.array`, optional
        An optional lookup table to use.

    Returns
    -------
    outImage: `lsst.afw.image.Image`
        The output linearized image.
    numOutOfRange: `int`
        The number of values that could not be linearized.

    Raises
    ------
    RuntimeError :
        Raised if an invalid linearityType is supplied.
    """
    numOutOfRange = 0
    for ampIdx, amp in enumerate(detector.getAmplifiers()):
        ampIdx = (ampIdx // 3, ampIdx % 3)
        bbox = amp.getBBox()
        imageView = image.Factory(image, bbox)

        if linearityType == 'Squared':
            sqCoeff = inputData[ampIdx]
            array = imageView.getArray()

            array[:] = array + sqCoeff*array**2
        elif linearityType == 'LookupTable':
            rowInd, colIndOffset = inputData[ampIdx]
            rowInd = int(rowInd)
            tableRow = table[rowInd, :]
            numOutOfRange += applyLookupTable(imageView, tableRow, colIndOffset)
        elif linearityType == 'Polynomial':
            coeffs = inputData[ampIdx]
            array = imageView.getArray()
            summation = np.zeros_like(array)
            for index, coeff in enumerate(coeffs):
                summation += coeff*np.power(array, (index + 2))
            array += summation
        elif linearityType == 'Spline':
            centers, values = np.split(inputData, 2)  # This uses the full data
            interp = afwMath.makeInterpolate(centers.tolist(), values.tolist(),
                                             afwMath.stringToInterpStyle('AKIMA_SPLINE'))
            array = imageView.getArray()
            delta = interp.interpolate(array.flatten())
            array -= np.array(delta).reshape(array.shape)
        else:
            raise RuntimeError(f"Unknown linearity: {linearityType}")
    return image, numOutOfRange


class LinearizeTestCase(lsst.utils.tests.TestCase):
    """Unit tests for linearizers.
    """

    def setUp(self):
        # This uses the same arbitrary values used in previous tests.
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-31, 22), lsst.geom.Extent2I(100, 85))
        self.ampArrangement = (2, 3)
        self.numAmps = self.ampArrangement[0]*self.ampArrangement[1]
        # Squared Parameters
        self.sqCoeffs = np.array([[0, 5e-6, 2.5e-5], [1e-5, 1.1e-6, 2.1e-6]], dtype=float)

        # Lookup Table Parameters
        self.colIndOffsets = np.array([[0, -50, 2.5], [37, 1, -3]], dtype=float)
        self.rowInds = np.array([[0, 1, 4], [3, 5, 2]])
        # This creates a 2x3 array (matching the amplifiers) that contains a
        # 2x1 array containing [colIndOffset_i, rowInd_i].
        self.lookupIndices = np.transpose(np.stack((self.rowInds, self.colIndOffsets), axis=0),
                                          axes=[1, 2, 0])

        self.table = np.random.normal(scale=55, size=(self.numAmps, 2500))
        self.assertLess(np.max(self.rowInds), self.numAmps, "error in test conditions; invalid row index")

        # Polynomial Parameters: small perturbation on Squared
        self.polyCoeffs = np.array([[[0, 1e-7], [5e-6, 1e-7], [2.5e-5, 1e-7]],
                                    [[1e-5, 1e-7], [1.1e-6, 1e-7], [2.1e-6, 1e-7]]], dtype=float)

        # Spline coefficients: should match a 1e-6 Squared solution
        self.splineCoeffs = np.array([-100, 0.0, 1000, 2000, 3000, 4000, 5000,
                                      0.0, 0.0, 1.0, 4.0, 9.0, 16.0, 25.0])
        self.log = logging.getLogger("lsst.ip.isr.testLinearizer")

    def tearDown(self):
        # destroy LSST objects so memory test passes.
        self.bbox = None
        self.detector = None

    def compareResults(self, linearizedImage, linearizedOutOfRange, linearizedCount, linearizedAmps,
                       referenceImage, referenceOutOfRange, referenceCount, referenceAmps):
        """Run assert tests on results.

        Parameters
        ----------
        linearizedImage : `lsst.afw.image.Image`
            Corrected image.
        linearizedOutOfRange : `int`
            Number of measured out-of-range pixels.
        linearizedCount : `int`
            Number of amplifiers that should be linearized.
        linearizedAmps : `int`
            Total number of amplifiers checked.
        referenceImage : `lsst.afw.image.Image`
            Truth image to compare against.
        referenceOutOfRange : `int`
            Number of expected out-of-range-pixels.
        referenceCount : `int`
            Number of amplifiers that are expected to be linearized.
        referenceAmps : `int`
            Expected number of amplifiers checked.
        """
        self.assertImagesAlmostEqual(linearizedImage, referenceImage)
        self.assertEqual(linearizedOutOfRange, referenceOutOfRange)
        self.assertEqual(linearizedCount, referenceCount)
        self.assertEqual(linearizedAmps, referenceAmps)

    def testBasics(self):
        """Test basic linearization functionality.
        """
        for imageClass in (afwImage.ImageF, afwImage.ImageD):
            inImage = makeRampImage(bbox=self.bbox, start=-5, stop=2500, imageClass=imageClass)

            for linearityType in ('Squared', 'LookupTable', 'Polynomial', 'Spline'):
                detector = self.makeDetector(linearityType)
                table = None
                inputData = {'Squared': self.sqCoeffs,
                             'LookupTable': self.lookupIndices,
                             'Polynomial': self.polyCoeffs,
                             'Spline': self.splineCoeffs}[linearityType]
                if linearityType == 'LookupTable':
                    table = np.array(self.table, dtype=inImage.getArray().dtype)
                linearizer = Linearizer(detector=detector, table=table)

                measImage = inImage.Factory(inImage, True)
                result = linearizer.applyLinearity(measImage, detector=detector, log=self.log)
                refImage, refNumOutOfRange = referenceImage(inImage.Factory(inImage, True),
                                                            detector, linearityType, inputData, table)

                # This is necessary for the same tests to be used on
                # all types.  The first amplifier has 0.0 for the
                # coefficient, which should be tested (it has a log
                # message), but we are not linearizing an amplifier
                # with no correction, so it fails the test that
                # numLinearized == numAmps.
                zeroLinearity = 1 if linearityType == 'Squared' else 0

                self.compareResults(measImage, result.numOutOfRange, result.numLinearized, result.numAmps,
                                    refImage, refNumOutOfRange, self.numAmps - zeroLinearity, self.numAmps)

                # Test a stand alone linearizer.  This ignores validate checks.
                measImage = inImage.Factory(inImage, True)
                storedLinearizer = self.makeLinearizer(linearityType)
                storedResult = storedLinearizer.applyLinearity(measImage, log=self.log)

                self.compareResults(measImage, storedResult.numOutOfRange, storedResult.numLinearized,
                                    storedResult.numAmps,
                                    refImage, refNumOutOfRange, self.numAmps - zeroLinearity, self.numAmps)

                # "Save to yaml" and test again
                storedDict = storedLinearizer.toDict()
                storedLinearizer = Linearizer().fromDict(storedDict)

                measImage = inImage.Factory(inImage, True)
                storedLinearizer = self.makeLinearizer(linearityType)
                storedResult = storedLinearizer.applyLinearity(measImage, log=self.log)

                self.compareResults(measImage, storedResult.numOutOfRange, storedResult.numLinearized,
                                    storedResult.numAmps,
                                    refImage, refNumOutOfRange, self.numAmps - zeroLinearity, self.numAmps)

                # "Save to fits" and test again
                storedTable = storedLinearizer.toTable()
                storedLinearizer = Linearizer().fromTable(storedTable)

                measImage = inImage.Factory(inImage, True)
                storedLinearizer = self.makeLinearizer(linearityType)
                storedResult = storedLinearizer.applyLinearity(measImage, log=self.log)

                self.compareResults(measImage, storedResult.numOutOfRange, storedResult.numLinearized,
                                    storedResult.numAmps,
                                    refImage, refNumOutOfRange, self.numAmps - zeroLinearity, self.numAmps)

    def makeDetector(self, linearityType, bbox=None):
        """Generate a fake detector for the test.

        Parameters
        ----------
        linearityType : `str`
            Which linearity to assign to the detector's cameraGeom.
        bbox : `lsst.geom.Box2I`, optional
            Bounding box to use for the detector.

        Returns
        -------
        detBuilder : `lsst.afw.cameraGeom.Detector`
            The fake detector.
        """
        bbox = bbox if bbox is not None else self.bbox
        numAmps = self.ampArrangement

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
                if linearityType == 'Squared':
                    ampInfo.setLinearityCoeffs([self.sqCoeffs[i, j]])
                elif linearityType == 'LookupTable':
                    # setLinearityCoeffs is picky about getting a mixed int/float list.
                    ampInfo.setLinearityCoeffs(np.array([self.rowInds[i, j], self.colIndOffsets[i, j],
                                                         0, 0], dtype=float))
                elif linearityType == 'Polynomial':
                    ampInfo.setLinearityCoeffs(self.polyCoeffs[i, j])
                elif linearityType == 'Spline':
                    ampInfo.setLinearityCoeffs(self.splineCoeffs)
                detBuilder.append(ampInfo)

        return detBuilder

    def makeLinearizer(self, linearityType, bbox=None):
        """Construct a linearizer with the test coefficients.

        Parameters
        ----------
        linearityType : `str`
            Type of linearity to use.  The coefficients are set by the
            setUp method.
        bbox : `lsst.geom.Box2I`
            Bounding box for the full detector.  Used to assign
            amp-based bounding boxes.

        Returns
        -------
        linearizer : `lsst.ip.isr.Linearizer`
            A fully constructed, persistable linearizer.
        """
        bbox = bbox if bbox is not None else self.bbox
        numAmps = self.ampArrangement
        boxArr = BoxGrid(box=bbox, numColRow=numAmps)
        linearizer = Linearizer()
        linearizer.hasLinearity = True

        for i in range(numAmps[0]):
            for j in range(numAmps[1]):
                ampName = f"amp {i+1}_{j+1}"
                ampBox = boxArr[i, j]
                linearizer.ampNames.append(ampName)

                if linearityType == 'Squared':
                    linearizer.linearityCoeffs[ampName] = np.array([self.sqCoeffs[i, j]])
                elif linearityType == 'LookupTable':
                    linearizer.linearityCoeffs[ampName] = np.array(self.lookupIndices[i, j])
                    linearizer.tableData = self.table
                elif linearityType == 'Polynomial':
                    linearizer.linearityCoeffs[ampName] = np.array(self.polyCoeffs[i, j])
                elif linearityType == 'Spline':
                    linearizer.linearityCoeffs[ampName] = np.array(self.splineCoeffs)

                linearizer.linearityType[ampName] = linearityType
                linearizer.linearityBBox[ampName] = ampBox
                linearizer.fitParams[ampName] = np.array([])
                linearizer.fitParamsErr[ampName] = np.array([])
                linearizer.fitChiSq[ampName] = np.nan

        return linearizer


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
