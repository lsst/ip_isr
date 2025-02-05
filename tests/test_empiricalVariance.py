#
# LSST Data Management System
# Copyright 2018 LSST Corporation.
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
import numpy as np
import unittest

import lsst.utils.tests

from lsst.daf.base import PropertyList
import lsst.afw.cameraGeom as cameraGeom
from lsst.geom import Point2I, Extent2I, Box2I, Extent2D
from lsst.afw.image import ExposureF, VisitInfo

from lsst.ip.isr.isrTask import IsrTask


def makeAmplifier(name, bbox, rawImageBox, overscanBox, gain, readNoise, saturation):
    amp = cameraGeom.Amplifier.Builder()
    amp.setName(name)
    amp.setBBox(bbox)
    amp.setRawDataBBox(rawImageBox)
    amp.setRawHorizontalOverscanBBox(overscanBox)
    amp.setGain(gain)
    amp.setReadNoise(readNoise)
    amp.setSaturation(saturation)
    amp.setSuspectLevel(np.nan)
    return amp


class EmpiricalVarianceTestCast(lsst.utils.tests.TestCase):
    def setUp(self):
        """Constructs a CCD with two amplifiers and prepares for ISR"""
        rng = np.random.Generator(np.random.MT19937(42))

        baseValue = 100.0
        gain = 1.0
        readNoise = 123456789.0
        saturation = 987654321.0
        height = 234
        imageSize = Extent2I(123, height)
        overscanSize = Extent2I(16, height)
        self.sigma = 1.234

        # Set up the various regions
        overscan1 = Box2I(Point2I(0, 0), overscanSize)
        image1 = Box2I(Point2I(overscanSize[0], 0), imageSize)
        image2 = Box2I(Point2I(overscanSize[0] + imageSize[0], 0), imageSize)
        overscan2 = Box2I(Point2I(overscanSize[0] + 2*imageSize[0], 0), overscanSize)

        leftBox = Box2I(overscan1.getMin(), Extent2I(overscan1.getWidth() + image1.getWidth(), height))
        rightBox = Box2I(image2.getMin(), Extent2I(image2.getWidth() + overscan2.getWidth(), height))

        target1 = Box2I(Point2I(0, 0), imageSize)
        target2 = Box2I(Point2I(image1.getWidth(), 0), imageSize)

        # Set the pixels
        exposure = ExposureF(Box2I(Point2I(0, 0), Extent2I(imageSize[0]*2 + overscanSize[0]*2, height)))
        yy = np.arange(0, height, 1, dtype=np.float32)
        leftImage = ExposureF(exposure, leftBox)
        leftImage.image.array[:] = baseValue + yy[:, np.newaxis]
        rightImage = ExposureF(exposure, rightBox)
        rightImage.image.array[:] = baseValue - yy[:, np.newaxis]

        leftOverscan = ExposureF(exposure, overscan1)
        leftOverscan.image.array += rng.normal(0.0, self.sigma, leftOverscan.image.array.shape)
        rightOverscan = ExposureF(exposure, overscan2)
        rightOverscan.image.array += rng.normal(0.0, self.sigma, leftOverscan.image.array.shape)
        exposure.mask.array[:] = 0.0
        exposure.variance.array[:] = np.nan

        # Construct the detectors
        amp1 = makeAmplifier("left", target1, image1, overscan1, gain, readNoise, saturation)
        amp2 = makeAmplifier("right", target2, image2, overscan2, gain, readNoise, saturation)
        ccdBox = Box2I(Point2I(0, 0), Extent2I(image1.getWidth() + image2.getWidth(), height))
        camBuilder = cameraGeom.Camera.Builder("fakeCam")
        detBuilder = camBuilder.add("detector", 1)
        detBuilder.setSerial("det1")
        detBuilder.setBBox(ccdBox)
        detBuilder.setPixelSize(Extent2D(1.0, 1.0))
        detBuilder.setOrientation(cameraGeom.Orientation())
        detBuilder.append(amp1)
        detBuilder.append(amp2)
        cam = camBuilder.finish()
        exposure.setDetector(cam.get('detector'))

        header = PropertyList()
        header.add("EXPTIME", 0.0)
        exposure.getInfo().setVisitInfo(VisitInfo(header))

        self.exposure = exposure
        self.config = IsrTask.ConfigClass()

        # Disable everything we don't care about
        self.config.doBias = False
        self.config.doDark = False
        self.config.doFlat = False
        self.config.doFringe = False
        self.config.doDefect = False
        self.config.doWrite = False
        self.config.expectWcs = False
        self.config.doLinearize = False
        self.config.doCrosstalk = False
        self.config.doBrighterFatter = False
        self.config.doAttachTransmissionCurve = False
        self.config.doAssembleCcd = False
        self.config.doNanMasking = False
        self.config.doInterpolate = False

        self.config.maskNegativeVariance = False  # This runs on mocks.
        # Set the things that match our test setup
        self.config.overscan.fitType = "CHEB"
        self.config.overscan.order = 1
        self.config.doEmpiricalReadNoise = True

        self.task = IsrTask(config=self.config)

    def tearDown(self):
        del self.exposure

    def testEmpiricalVariance(self):
        results = self.task.run(self.exposure)
        postIsr = results.exposure

        self.assertFloatsEqual(postIsr.mask.array, 0)
        # Image is not exactly zero because the noise in the overscan (required
        # to be able to set the empirical variance) leads to a slight
        # misestimate in the polynomial fit.
        self.assertFloatsAlmostEqual(np.median(postIsr.image.array), 0.0, atol=5.0e-2)
        self.assertFloatsAlmostEqual(np.nanmedian(postIsr.variance.array), self.sigma**2, rtol=5.0e-2)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
