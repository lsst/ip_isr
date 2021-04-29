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
import itertools
import tempfile

import numpy as np

import lsst.geom
import lsst.utils.tests
import lsst.afw.image
import lsst.afw.table
import lsst.afw.cameraGeom as cameraGeom

from lsst.pipe.base import Struct
from lsst.ip.isr import IsrTask, CrosstalkCalib, CrosstalkTask, NullCrosstalkTask

try:
    display
except NameError:
    display = False
else:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)


outputName = None    # specify a name (as a string) to save the output crosstalk coeffs.


class CrosstalkTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        width, height = 250, 500
        self.numAmps = 4
        numPixelsPerAmp = 1000
        # crosstalk[i][j] is the fraction of the j-th amp present on the i-th amp.
        self.crosstalk = [[0.0, 1e-4, 2e-4, 3e-4],
                          [3e-4, 0.0, 2e-4, 1e-4],
                          [4e-4, 5e-4, 0.0, 6e-4],
                          [7e-4, 8e-4, 9e-4, 0.0]]
        self.value = 12345
        self.crosstalkStr = "XTLK"

        # A bit of noise is important, because otherwise the pixel distributions are razor-thin
        # and then rejection doesn't work
        rng = np.random.RandomState(12345)
        self.noise = rng.normal(0.0, 0.1, (2*height, 2*width))

        # Create amp images
        withoutCrosstalk = [lsst.afw.image.ImageF(width, height) for _ in range(self.numAmps)]
        for image in withoutCrosstalk:
            image.set(0)
            xx = rng.randint(0, width, numPixelsPerAmp)
            yy = rng.randint(0, height, numPixelsPerAmp)
            image.getArray()[yy, xx] = self.value

        # Add in crosstalk
        withCrosstalk = [image.Factory(image, True) for image in withoutCrosstalk]
        for ii, iImage in enumerate(withCrosstalk):
            for jj, jImage in enumerate(withoutCrosstalk):
                value = self.crosstalk[ii][jj]
                iImage.scaledPlus(value, jImage)

        # Put amp images together
        def construct(imageList):
            image = lsst.afw.image.ImageF(2*width, 2*height)
            image.getArray()[:height, :width] = imageList[0].getArray()
            image.getArray()[:height, width:] = imageList[1].getArray()[:, ::-1]  # flip in x
            image.getArray()[height:, :width] = imageList[2].getArray()[::-1, :]  # flip in y
            image.getArray()[height:, width:] = imageList[3].getArray()[::-1, ::-1]  # flip in x and y
            image.getArray()[:] += self.noise
            return image

        # Construct detector
        detName = 'detector 1'
        detId = 1
        detSerial = 'serial 1'
        orientation = cameraGeom.Orientation()
        pixelSize = lsst.geom.Extent2D(1, 1)
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                               lsst.geom.Extent2I(2*width, 2*height))
        crosstalk = np.array(self.crosstalk, dtype=np.float32)

        camBuilder = cameraGeom.Camera.Builder("fakeCam")
        detBuilder = camBuilder.add(detName, detId)
        detBuilder.setSerial(detSerial)
        detBuilder.setBBox(bbox)
        detBuilder.setOrientation(orientation)
        detBuilder.setPixelSize(pixelSize)
        detBuilder.setCrosstalk(crosstalk)

        # Construct second detector in this fake camera
        detName = 'detector 2'
        detId = 2
        detSerial = 'serial 2'
        orientation = cameraGeom.Orientation()
        pixelSize = lsst.geom.Extent2D(1, 1)
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                               lsst.geom.Extent2I(2*width, 2*height))
        crosstalk = np.array(self.crosstalk, dtype=np.float32)

        detBuilder2 = camBuilder.add(detName, detId)
        detBuilder2.setSerial(detSerial)
        detBuilder2.setBBox(bbox)
        detBuilder2.setOrientation(orientation)
        detBuilder2.setPixelSize(pixelSize)
        detBuilder2.setCrosstalk(crosstalk)

        # Create amp info
        for ii, (xx, yy, corner) in enumerate([(0, 0, lsst.afw.cameraGeom.ReadoutCorner.LL),
                                               (width, 0, lsst.afw.cameraGeom.ReadoutCorner.LR),
                                               (0, height, lsst.afw.cameraGeom.ReadoutCorner.UL),
                                               (width, height, lsst.afw.cameraGeom.ReadoutCorner.UR)]):

            amp = cameraGeom.Amplifier.Builder()
            amp.setName("amp %d" % ii)
            amp.setBBox(lsst.geom.Box2I(lsst.geom.Point2I(xx, yy),
                                        lsst.geom.Extent2I(width, height)))
            amp.setRawDataBBox(lsst.geom.Box2I(lsst.geom.Point2I(xx, yy),
                                               lsst.geom.Extent2I(width, height)))
            amp.setReadoutCorner(corner)
            detBuilder.append(amp)
            detBuilder2.append(amp)

        cam = camBuilder.finish()
        ccd1 = cam.get('detector 1')
        ccd2 = cam.get('detector 2')

        self.exposure = lsst.afw.image.makeExposure(lsst.afw.image.makeMaskedImage(construct(withCrosstalk)))
        self.exposure.setDetector(ccd1)

        # Create a single ctSource that will be used for interChip CT correction.
        self.ctSource = lsst.afw.image.makeExposure(lsst.afw.image.makeMaskedImage(construct(withCrosstalk)))
        self.ctSource.setDetector(ccd2)

        self.corrected = construct(withoutCrosstalk)

        if display:
            disp = lsst.afw.display.Display(frame=1)
            disp.mtv(self.exposure, title="exposure")
            disp = lsst.afw.display.Display(frame=0)
            disp.mtv(self.corrected, title="corrected exposure")

    def tearDown(self):
        del self.exposure
        del self.corrected

    def checkCoefficients(self, coeff, coeffErr, coeffNum):
        """Check that coefficients are as expected

        Parameters
        ----------
        coeff : `numpy.ndarray`
            Crosstalk coefficients.
        coeffErr : `numpy.ndarray`
            Crosstalk coefficient errors.
        coeffNum : `numpy.ndarray`
            Number of pixels to produce each coefficient.
        """
        for matrix in (coeff, coeffErr, coeffNum):
            self.assertEqual(matrix.shape, (self.numAmps, self.numAmps))
        self.assertFloatsAlmostEqual(coeff, np.array(self.crosstalk), atol=1.0e-6)

        for ii in range(self.numAmps):
            self.assertEqual(coeff[ii, ii], 0.0)
            self.assertTrue(np.isnan(coeffErr[ii, ii]))
            self.assertEqual(coeffNum[ii, ii], 1)

        self.assertTrue(np.all(coeffErr[ii, jj] > 0 for ii, jj in
                               itertools.product(range(self.numAmps), range(self.numAmps)) if ii != jj))
        self.assertTrue(np.all(coeffNum[ii, jj] > 0 for ii, jj in
                               itertools.product(range(self.numAmps), range(self.numAmps)) if ii != jj))

    def checkSubtracted(self, exposure):
        """Check that the subtracted image is as expected

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Crosstalk-subtracted exposure.
        """
        image = exposure.getMaskedImage().getImage()
        mask = exposure.getMaskedImage().getMask()
        self.assertFloatsAlmostEqual(image.getArray(), self.corrected.getArray(), atol=2.0e-2)
        self.assertIn(self.crosstalkStr, mask.getMaskPlaneDict())
        self.assertGreater((mask.getArray() & mask.getPlaneBitMask(self.crosstalkStr) > 0).sum(), 0)

    def testDirectAPI(self):
        """Test that individual function calls work"""
        calib = CrosstalkCalib()
        calib.coeffs = np.array(self.crosstalk).transpose()
        calib.subtractCrosstalk(self.exposure, crosstalkCoeffs=calib.coeffs,
                                minPixelToMask=self.value - 1,
                                crosstalkStr=self.crosstalkStr)
        self.checkSubtracted(self.exposure)

        outPath = tempfile.mktemp() if outputName is None else "{}-isrCrosstalk".format(outputName)
        outPath += '.yaml'
        calib.writeText(outPath)

    def testTaskAPI(self):
        """Test that the Tasks work

        Checks both MeasureCrosstalkTask and the CrosstalkTask.
        """
        # make exposure available to NullIsrTask
        # without NullIsrTask's `self` hiding this test class's `self`
        exposure = self.exposure

        class NullIsrTask(IsrTask):
            def runDataRef(self, dataRef):
                return Struct(exposure=exposure)

        coeff = np.array(self.crosstalk).transpose()
        config = IsrTask.ConfigClass()
        config.crosstalk.minPixelToMask = self.value - 1
        config.crosstalk.crosstalkMaskPlane = self.crosstalkStr
        isr = IsrTask(config=config)
        calib = CrosstalkCalib().fromDetector(self.exposure.getDetector(), coeffVector=coeff)
        isr.crosstalk.run(self.exposure, crosstalk=calib)
        self.checkSubtracted(self.exposure)

    def test_prepCrosstalk(self):
        """Test that prep crosstalk does not error when given a dataRef with no
        crosstalkSources to find.
        """
        dataRef = Struct(dataId={'fake': 1})
        task = CrosstalkTask()
        result = task.prepCrosstalk(dataRef)
        self.assertIsNone(result)

    def test_nullCrosstalkTask(self):
        """Test that the null crosstalk task does not create an error.
        """
        exposure = self.exposure
        task = NullCrosstalkTask()
        result = task.run(exposure, crosstalkSources=None)
        self.assertIsNone(result)

    def test_interChip(self):
        """Test that passing an external exposure as the crosstalk source works.
        """
        exposure = self.exposure
        ctSources = [self.ctSource]

        coeff = np.array(self.crosstalk).transpose()
        calib = CrosstalkCalib().fromDetector(exposure.getDetector(), coeffVector=coeff)
        # Now convert this into zero intra-chip, full inter-chip:
        calib.interChip['detector 2'] = coeff
        calib.coeffs = np.zeros_like(coeff)

        # Process and check as above
        config = IsrTask.ConfigClass()
        config.crosstalk.minPixelToMask = self.value - 1
        config.crosstalk.crosstalkMaskPlane = self.crosstalkStr
        isr = IsrTask(config=config)
        isr.crosstalk.run(exposure, crosstalk=calib, crosstalkSources=ctSources)
        self.checkSubtracted(exposure)

    def test_crosstalkIO(self):
        """Test that crosstalk doesn't change on being converted to persistable formats."""

        # Add the interchip crosstalk as in the previous test.
        exposure = self.exposure

        coeff = np.array(self.crosstalk).transpose()
        calib = CrosstalkCalib().fromDetector(exposure.getDetector(), coeffVector=coeff)
        # Now convert this into zero intra-chip, full inter-chip:
        calib.interChip['detector 2'] = coeff

        outPath = tempfile.mktemp() + '.yaml'
        calib.writeText(outPath)
        newCrosstalk = CrosstalkCalib().readText(outPath)
        self.assertEqual(calib, newCrosstalk)

        outPath = tempfile.mktemp() + '.fits'
        calib.writeFits(outPath)
        newCrosstalk = CrosstalkCalib().readFits(outPath)
        self.assertEqual(calib, newCrosstalk)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
