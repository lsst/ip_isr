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
import unittest

import lsst.utils.tests
import lsst.geom
import lsst.afw.geom as afwGeom
from lsst.afw.image import ExposureF, ExposureInfo, PhotoCalib, VisitInfo
from lsst.afw.geom.wcsUtils import makeDistortedTanWcs
from lsst.afw.cameraGeom import FIELD_ANGLE, FOCAL_PLANE, PIXELS
from lsst.afw.cameraGeom.testUtils import DetectorWrapper, CameraWrapper
from lsst.afw.geom.utils import wcsAlmostEqualOverBBox
from lsst.ip.isr import IsrTask
import lsst.ip.isr.isrFunctions as isrFunctions


class AddDistortionModelTestCase(lsst.utils.tests.TestCase):
    """Test IsrTask.addDistortionModel.

    DEPRECATED: to be removed with addDistortionModel
    """

    def setUp(self):
        self.camera = CameraWrapper().camera
        self.detector = DetectorWrapper().detector
        self.crpix = lsst.geom.Point2D(50, 100)
        self.crval = lsst.geom.SpherePoint(36, 71, lsst.geom.degrees)
        scale = 1.0*lsst.geom.arcseconds
        self.cdMatrix = afwGeom.makeCdMatrix(scale=scale)
        self.wcs = afwGeom.makeSkyWcs(crpix=self.crpix, crval=self.crval, cdMatrix=self.cdMatrix)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-10, 10), lsst.geom.Extent2I(1000, 1022))
        self.exposure = ExposureF(self.bbox)

        # set the few items of ExposureInfo needed by IsrTask.run
        # when only adding a distortion model
        exposureInfo = ExposureInfo(photoCalib=PhotoCalib(1.0),
                                    detector=self.detector,
                                    visitInfo=VisitInfo(exposureTime=1.0),
                                    wcs=self.wcs)

        self.exposure.setInfo(exposureInfo)

    def tearDown(self):
        self.detector = None
        self.exposure = None

    def testAddDistortionMethod(self):
        """Call IsrTask.addDistortionModel directly"""
        isrFunctions.addDistortionModel(self.exposure, self.camera)
        self.assertFalse(wcsAlmostEqualOverBBox(self.wcs, self.exposure.getWcs(), self.bbox))

        desiredWcs = self.makeDesiredDistortedWcs()
        self.assertWcsAlmostEqualOverBBox(desiredWcs, self.exposure.getWcs(), self.bbox)

    def makeMinimalIsrConfig(self):
        """Return an IsrConfig with all boolean flags disabled"""
        isrConfig = IsrTask.ConfigClass()
        for name in isrConfig:
            if name.startswith("do"):
                setattr(isrConfig, name, False)
        return isrConfig

    def makeDesiredDistortedWcs(self):
        """Make the expected distorted WCS"""
        pixelToFocalPlane = self.detector.getTransform(PIXELS, FOCAL_PLANE)
        focalPlaneToFieldAngle = self.camera.getTransformMap().getTransform(FOCAL_PLANE, FIELD_ANGLE)
        return makeDistortedTanWcs(self.wcs, pixelToFocalPlane, focalPlaneToFieldAngle)

    def testRunWithoutAddDistortionModel(self):
        """Test IsrTask.run with config.doAddDistortionModel false"""
        isrConfig = self.makeMinimalIsrConfig()
        isrTask = IsrTask(config=isrConfig)

        # the camera argument is not needed
        exposure = isrTask.run(ccdExposure=self.exposure).exposure
        self.assertEqual(self.wcs, exposure.getWcs())

        # and the camera argument is ignored if provided
        exposure2 = isrTask.run(ccdExposure=self.exposure, camera=self.camera).exposure
        self.assertEqual(self.wcs, exposure2.getWcs())


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
