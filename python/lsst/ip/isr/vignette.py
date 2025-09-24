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

__all__ = ('VignetteConfig', 'VignetteTask', 'setValidPolygonCcdIntersect',
           'maskVignettedRegion')

import logging
import numpy as np

import lsst.geom as geom
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.geom as afwGeom

from lsst.pex.config import Config, Field
from lsst.pipe.base import Task


class VignetteConfig(Config):
    """Settings to define vignetting pattern.
    """
    xCenter = Field(
        dtype=float,
        doc="Center of vignetting pattern, in focal plane x coordinates.",
        default=0.0,
    )
    yCenter = Field(
        dtype=float,
        doc="Center of vignetting pattern, in focal plane y coordinates.",
        default=0.0,
    )
    radius = Field(
        dtype=float,
        doc="Radius of vignetting pattern, in focal plane coordinates.",
        default=100.0,
        check=lambda x: x >= 0
    )
    numPolygonPoints = Field(
        dtype=int,
        doc="Number of points used to define the vignette polygon.",
        default=100,
    )


class VignetteTask(Task):
    """Define a simple circular vignette pattern and optionally update mask
    plane.
    """
    ConfigClass = VignetteConfig
    _DefaultName = "isrVignette"

    def run(
        self,
        exposure=None,
        doUpdateMask=True,
        maskPlane="NO_DATA",
        vignetteValue=None,
        log=None,
        doUpdatePolygon=True,
    ):
        """Generate circular vignette pattern.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`, optional
            Exposure to construct, apply, and optionally mask vignette for.
        doUpdateMask : `bool`, optional
            If true, the mask will be updated to mask the vignetted region.
        maskPlane : `str`, optional
            Mask plane to assign vignetted pixels to.
        vignetteValue : `float` or `None`, optional
            Value to assign to the image array pixels within the ``polygon``
            region.  If `None`, image pixel values are not replaced.
        log : `logging.Logger`, optional
            Log object to write to.
        doUpdatePolygon : `bool`, optional
            If true, the valid polygon will be updated.

        Returns
        -------
        polygon : `lsst.afw.geom.Polygon`
            Polygon defining the boundary of the vignetted region.
        """
        theta = np.linspace(0, 2*np.pi, num=self.config.numPolygonPoints, endpoint=False)
        x = self.config.radius*np.cos(theta) + self.config.xCenter
        y = self.config.radius*np.sin(theta) + self.config.yCenter
        points = np.array([x, y]).transpose()
        fpPolygon = afwGeom.Polygon([geom.Point2D(x1, y1) for x1, y1 in reversed(points)])
        if exposure is None:
            return fpPolygon

        # Exposure was provided, so attach the validPolygon associated with the
        # vignetted region.
        validPolygon = setValidPolygonCcdIntersect(exposure, fpPolygon, log=log, setPolygon=doUpdatePolygon)

        if doUpdateMask:
            maskVignettedRegion(
                exposure,
                validPolygon,
                maskPlane="NO_DATA",
                vignetteValue=vignetteValue,
                log=log,
            )

        return fpPolygon


def setValidPolygonCcdIntersect(ccdExposure, fpPolygon, log=None, setPolygon=True):
    """Set valid polygon on ccdExposure associated with focal plane polygon.

    The ccd exposure's valid polygon is the intersection of fpPolygon,
    a valid polygon in focal plane coordinates, and the ccd corners,
    in ccd pixel coordinates.

    Parameters
    ----------
    ccdExposure : `lsst.afw.image.Exposure`
        Exposure to process.
    fpPolygon : `lsst.afw.geom.Polygon`
        Polygon in focal plane coordinates.
    log : `logging.Logger`, optional
        Log object to write to.
    setPolygon : `bool`, optional
        Actually set the polygon, or just return it?

    Returns
    -------
    validPolygon : `lsst.afw.geom.Polygon`
        Valid polygon with the intersection.
    """
    # Get ccd corners in focal plane coordinates
    ccd = ccdExposure.getDetector()
    fpCorners = ccd.getCorners(cameraGeom.FOCAL_PLANE)
    ccdPolygon = afwGeom.Polygon(fpCorners)
    # Get intersection of ccd corners with fpPolygon
    try:
        intersect = ccdPolygon.intersectionSingle(fpPolygon)
    except afwGeom.SinglePolygonException:
        intersect = None
    if intersect is not None:
        # Transform back to pixel positions and build new polygon
        ccdPoints = ccd.transform(intersect, cameraGeom.FOCAL_PLANE, cameraGeom.PIXELS)
        validPolygon = afwGeom.Polygon(ccdPoints)
        if setPolygon:
            ccdExposure.getInfo().setValidPolygon(validPolygon)
    else:
        if log is not None:
            log.info("Ccd exposure does not overlap with focal plane polygon.  Not setting validPolygon.")
        validPolygon = None

    return validPolygon


def maskVignettedRegion(exposure, polygon, maskPlane="NO_DATA", vignetteValue=None, log=None):
    """Add mask bit to image pixels according to vignetted polygon region.

    NOTE: this function could be used to mask and replace pixels associated
    with any polygon in the exposure pixel coordinates.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Image whose mask plane is to be updated.
    polygon : `lsst.afw.geom.Polygon`
        Polygon region defining the vignetted region in the pixel coordinates
        of ``exposure``.
    maskPlane : `str`, optional
        Mask plane to assign vignetted pixels to.
    vignetteValue : `float` or `None`, optional
        Value to assign to the image array pixels within the ``polygon``
        region.  If `None`, image pixel values are not replaced.
    log : `logging.Logger`, optional
        Log object to write to.

    Raises
    ------
    RuntimeError
        Raised if no valid polygon exists.
    """
    log = log if log else logging.getLogger(__name__)
    if not polygon:
        log.info("No polygon provided.  Masking nothing.")
        return

    fullyIlluminated = True
    if not all(polygon.contains(exposure.getBBox().getCorners())):
        fullyIlluminated = False
    log.info("Exposure is fully illuminated? %s", fullyIlluminated)

    if not fullyIlluminated:
        # Scan pixels.
        mask = exposure.getMask()
        xx, yy = np.meshgrid(np.arange(0, mask.getWidth(), dtype=float),
                             np.arange(0, mask.getHeight(), dtype=float))

        vignMask = ~(polygon.contains(xx, yy))

        bitMask = mask.getPlaneBitMask(maskPlane)
        maskArray = mask.getArray()
        maskArray[vignMask] |= bitMask
        log.info("Exposure contains {} vignetted pixels which are now masked with mask plane {}.".
                 format(np.count_nonzero(vignMask), maskPlane))

        if vignetteValue is not None:
            imageArray = exposure.getImage().getArray()
            imageArray[vignMask] = vignetteValue
            log.info("Vignetted pixels in image array have been replaced with {}.".format(vignetteValue))
