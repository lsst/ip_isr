#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#

"""
Tests setting of valid focal plane polygon intersection with ccd corners

Run with:
   python testSetValidPolygonIntersect.py
or
   python
   >>> import testSetValidPolygonIntersect; testSetValidPolygonIntersect.run()
"""

import unittest

import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.ip.isr.vignette import setValidPolygonCcdIntersect
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
from lsst.afw.cameraGeom import PIXELS, FOCAL_PLANE


def makeCircularPolygon(fpCenterX, fpCenterY, fpRadius, numPolygonPoints):
    theta = np.linspace(0, 2*np.pi, num=numPolygonPoints, endpoint=False)
    xx = fpRadius*np.cos(theta) + fpCenterX
    yy = fpRadius*np.sin(theta) + fpCenterY
    points = np.array([xx, yy]).transpose()
    polygon = afwGeom.Polygon([lsst.geom.Point2D(x, y) for x, y in reversed(points)])
    return polygon


def makeSquarePolygon(fpX0, fpY0, fpSize):
    xx = [fpX0, fpX0, fpX0 + fpSize - 1, fpX0 + fpSize - 1, fpX0]
    yy = [fpY0, fpY0 + fpSize - 1, fpY0 + fpSize - 1, fpY0, fpY0]
    points = np.array([xx, yy]).transpose()
    polygon = afwGeom.Polygon([lsst.geom.Point2D(x, y) for x, y in points])
    return polygon


class SetValidPolygonIntersectTestCase(lsst.utils.tests.TestCase):
    """A test case for setting of valid focal plane polygon intersection with
    ccd corners.
    """

    def testSetPolygonIntersect(self):
        # Create a detector
        detector = DetectorWrapper().detector
        numPolygonPoints = 50
        # Create an exposure with bounding box defined by detector
        exposure = afwImage.ExposureF(detector.getBBox())
        exposure.setDetector(detector)

        pixelSizeMm = exposure.getDetector().getPixelSize()[0]

        pixX0 = exposure.getX0()
        pixY0 = exposure.getY0()
        pixX1 = pixX0 + exposure.getWidth() - 1
        pixY1 = pixY0 + exposure.getHeight() - 1

        fpCenter = exposure.getDetector().getCenter(FOCAL_PLANE)
        fpCenterX = fpCenter[0]
        fpCenterY = fpCenter[1]
        pixCenter = exposure.getDetector().getCenter(PIXELS)

        # Make a polygon that encompases entire ccd (radius of 2*max of
        # width/height)
        fpRadius = 2.0*max(exposure.getWidth()*pixelSizeMm, exposure.getHeight()*pixelSizeMm)
        fpPolygon = makeCircularPolygon(fpCenterX, fpCenterY, fpRadius, numPolygonPoints)
        # Set the polygon that is the intersection of fpPolygon and ccd
        setValidPolygonCcdIntersect(exposure, fpPolygon)
        # Since the ccd is fully contained in the fpPolygon, the intersection
        # should be the ccdPolygon itself
        ccdPolygonPix = afwGeom.Polygon(exposure.getDetector().getCorners(PIXELS))
        self.assertEqual(exposure.getInfo().getValidPolygon(), ccdPolygonPix)

        # Make a polygon that is entirely within, but smaller than, the ccd
        # (radius of 0.2*min of width/height)
        fpRadius = 0.2*min(exposure.getWidth()*pixelSizeMm, exposure.getHeight()*pixelSizeMm)
        fpPolygon = makeCircularPolygon(fpCenterX, fpCenterY, fpRadius, numPolygonPoints)
        # Set the polygon that is the intersection of fpPolygon and ccd
        setValidPolygonCcdIntersect(exposure, fpPolygon)
        # all vertices of polygon should be contained within the ccd
        for x in exposure.getInfo().getValidPolygon():
            self.assertTrue(ccdPolygonPix.contains(lsst.geom.Point2D(x)))
        # intersection is smaller than the ccd
        self.assertNotEqual(exposure.getInfo().getValidPolygon(), ccdPolygonPix)

        # make a simple square polygon that partly intersects the ccd, centered
        # at ccd center
        fpPolygonSize = max(exposure.getWidth()*pixelSizeMm, exposure.getHeight()*pixelSizeMm)
        fpPolygon = makeSquarePolygon(fpCenterX, fpCenterY, fpPolygonSize)
        setValidPolygonCcdIntersect(exposure, fpPolygon)
        # Check that the polygon contains the central pixel (offset by one to
        # actually be "contained")
        pixCenterPlusOne = lsst.geom.Point2D(pixCenter[0] + 1, pixCenter[1] + 1)
        self.assertTrue(exposure.getInfo().getValidPolygon().contains(lsst.geom.Point2D(pixCenterPlusOne)))
        # Check that the polygon contains the upper right ccd edge
        self.assertTrue(exposure.getInfo().getValidPolygon().contains(lsst.geom.Point2D(pixX1, pixY1)))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
