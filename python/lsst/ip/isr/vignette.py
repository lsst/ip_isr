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
import numpy

import lsst.afw.geom as afwGeom

from lsst.pex.config import Config, Field
from lsst.pipe.base import Task


class VignetteConfig(Config):
    """
    Settings to define vignetteing pattern
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
    doWriteVignettePolygon = Field(
        dtype=bool,
        doc="Persist polygon used to define vignetted region?",
        default=False,
    )


class VignetteTask(Task):
    """Define a simple circular vignette pattern.
    """
    ConfigClass = VignetteConfig
    _DefaultName = "isrVignette"

    def run(self, exposure):
        """Generate circular vignette pattern.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to construct vignette for.

        Returns
        -------
        polygon : `lsst.afw.Geom.Polygon`
            Polygon defining the boundary of the vignetted region.
        """

        if self.config.doWriteVignettePolygon:
            theta = numpy.linspace(0, 2*numpy.pi, num=self.config.numPolygonPoints, endpoint=False)
            x = self.config.radius*numpy.cos(theta) + self.config.xCenter
            y = self.config.radius*numpy.sin(theta) + self.config.yCenter
            points = numpy.array([x, y]).transpose()
            polygon = afwGeom.Polygon([afwGeom.Point2D(x1, y1) for x1, y1 in reversed(points)])
            return polygon
        else:
            return None
