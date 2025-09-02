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
"""
Flat gradient fit storage class.
"""

__all__ = ["FlatGradient"]

from astropy.table import Table
import numpy as np
from scipy.interpolate import Akima1DInterpolator

from lsst.ip.isr import IsrCalib


class FlatGradient(IsrCalib):
    """Flat gradient measurements.

    Parameters
    ----------
    log : `logging.Logger`, optional
        Log to write messages to. If `None` a default logger will be used.
    **kwargs :
        Additional parameters.
    """

    _OBSTYPE = "flatGradient"
    _SCHEMA = "FlatGradient"
    _VERSION = 1.0

    def __init__(self, **kwargs):

        self.radialSplineNodes = np.zeros(1)
        self.radialSplineValues = np.zeros(1)
        self.itlRatio = 1.0
        self.centroidX = 0.0
        self.centroidY = 0.0
        self.centroidDeltaX = 0.0
        self.centroidDeltaY = 0.0
        self.gradientX = 0.0
        self.gradientY = 0.0
        self.normalizationFactor = 1.0

        super().__init__(**kwargs)

        self.requiredAttributes.update(
            [
                "radialSplineNodes",
                "radialSplineValues",
                "itlRatio",
                "centroidX",
                "centroidY",
                "centroidDeltaX",
                "centroidDeltaY",
                "gradientX",
                "gradientY",
                "normalizationFactor",
            ],
        )

        self.updateMetadata(setCalibInfo=True, setCalibId=True, **kwargs)

    def setParameters(
        self,
        *,
        radialSplineNodes,
        radialSplineValues,
        itlRatio=1.0,
        centroidX=0.0,
        centroidY=0.0,
        centroidDeltaX=0.0,
        centroidDeltaY=0.0,
        gradientX=0.0,
        gradientY=0.0,
        normalizationFactor=1.0,
    ):
        """Set the parameters for the gradient model.

        Parameters
        ----------
        radialSplineNodes : `np.ndarray`
            Array of spline nodes.
        radialSplineValues : `np.ndarray`
            Array of spline values (same length as ``radialSplineNodes``).
        itlRatio : `float`, optional
            Ratio of flat for ITL detectors to E2V detectors.
        centroidX : `float`, optional
            X centroid of the focal plane (mm). This will be used as the
            pivot for the gradient plane.
        centroidY : `float`, optional
            Y centroid of the focal plane (mm). This will be used as the
            pivot for the gradient plane.
        centroidDeltaX : `float`, optional
            Centroid offset (mm). This is used in the radial function to
            allow for mis-centering in the illumination gradient.
        centroidDeltaY : `float`, optional
            Centroid offset (mm). This is used in the radial function to
            allow for mis-centering in the illumination gradient.
        gradientX : `float`, optional
            Slope of gradient in x direction (throughput/mm).
        gradientY : `float`, optional
            Slope of gradient in y direction (throughput/mm).
        normalizationFactor : `float`, optional
            Overall normalization factor (used to, e.g. make the
            center of the focal plane equal to 1.0 vs. a focal-plane
            average.
        """
        if len(radialSplineNodes) != len(radialSplineValues):
            raise ValueError("The number of spline nodes and values must be equal.")

        self.radialSplineNodes = radialSplineNodes
        self.radialSplineValues = radialSplineValues
        self.itlRatio = itlRatio
        self.centroidX = centroidX
        self.centroidY = centroidY
        self.centroidDeltaX = centroidDeltaX
        self.centroidDeltaY = centroidDeltaY
        self.gradientX = gradientX
        self.gradientY = gradientY
        self.normalizationFactor = normalizationFactor

    def computeRadialSplineModelXY(self, x, y):
        """Compute the radial spline model values from x/y.

        The spline model is a 1D Akima spline. When computed, the values
        from the model describe the radial function of the full focal
        plane flat-field. Dividing by this model will yield a radially
        flattened flat-field.

        Parameters
        ----------
        x : `np.ndarray`
            Array of focal plane x values (mm).
        y : `np.ndarray`
            Array of focal plane y values (mm).

        Returns
        -------
        splineModel : `np.ndarray`
            Spline model values at the x/y positions.
        """
        centroidX = self.centroidX + self.centroidDeltaX
        centroidY = self.centroidY + self.centroidDeltaY

        radius = np.sqrt((x - centroidX)**2. + (y - centroidY)**2.)

        return self.computeRadialSplineModel(radius)

    def computeRadialSplineModel(self, radius):
        """Compute the radial spline model values from radii.

        The spline model is a 1D Akima spline. When computed, the values
        from the model describe the radial function of the full focal
        plane flat-field. Dividing by this model will yield a radially
        flattened flat-field.

        Parameters
        ----------
        radius : `np.ndarray`
            Array of focal plane radii (mm).

        Returns
        -------
        splineModel : `np.ndarray`
            Spline model values at the radius positions.
        """
        spl = Akima1DInterpolator(self.radialSplineNodes, self.radialSplineValues)

        return spl(np.clip(radius, self.radialSplineNodes[0], self.radialSplineNodes[-1]))

    def computeGradientModel(self, x, y):
        """Compute the gradient model values.

        The gradient model is a plane constrained to be 1.0 at the
        ``centroidX``, ``centroidY`` values. Dividing by this model will
        remove the planar gradient in a flat field. Note that the planar
        gradient pivot is always at the same position, and does not
        move with the radial gradient centroid so as to keep the
        model fit more stable.

        Parameters
        ----------
        x : `np.ndarray`
            Array of focal plane x values (mm).
        y : `np.ndarray`
            Array of focal plane y values (mm).

        Returns
        -------
        gradientModel : `np.ndarray`
            Gradient model values at the x/y positions.
        """
        gradient = 1 + self.gradientX*(x - self.centroidX) + self.gradientY*(y - self.centroidY)

        return gradient

    def computeFullModel(self, x, y, is_itl):
        """Compute the full gradient model given x/y and itl booleans.

        This returns the full model that can be applied directly
        to data that was used in a fit.

        Parameters
        ----------
        x : `np.ndarray`
            Array of focal plane x values (mm).
        y : `np.ndarray`
            Array of focal plane y values (mm).
        is_itl : `np.ndarray`
            Boolean array of whether each point is from an ITL detector.

        Returns
        -------
        model : `np.ndarray`
            Model values at each position.
        """
        model = self.computeRadialSplineModelXY(x, y) / self.computeGradientModel(x, y)
        model[is_itl] *= self.itlRatio

        return model

    @classmethod
    def fromDict(cls, dictionary):
        """Construct a FlatGradient from a dictionary of properties.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties.

        Returns
        -------
        calib : `lsst.ip.isr.FlatGradient`
            Constructed calibration.
        """
        calib = cls()

        calib.setMetadata(dictionary["metadata"])

        calib.radialSplineNodes = np.asarray(dictionary["radialSplineNodes"])
        calib.radialSplineValues = np.asarray(dictionary["radialSplineValues"])
        calib.itlRatio = dictionary["itlRatio"]
        calib.centroidX = dictionary["centroidX"]
        calib.centroidY = dictionary["centroidY"]
        calib.centroidDeltaX = dictionary["centroidDeltaX"]
        calib.centroidDeltaY = dictionary["centroidDeltaY"]
        calib.gradientX = dictionary["gradientX"]
        calib.gradientY = dictionary["gradientY"]
        calib.normalizationFactor = dictionary["normalizationFactor"]

        calib.updateMetadata()
        return calib

    def toDict(self):
        """Return a dictionary containing the calibration properties.

        Returns
        -------
        dictionary : `dict`
            Dictionary of properties.
        """
        self.updateMetadata()

        outDict = dict()
        metadata = self.getMetadata()
        outDict["metadata"] = metadata

        outDict["radialSplineNodes"] = self.radialSplineNodes.tolist()
        outDict["radialSplineValues"] = self.radialSplineValues.tolist()
        outDict["itlRatio"] = float(self.itlRatio)
        outDict["centroidX"] = float(self.centroidX)
        outDict["centroidY"] = float(self.centroidY)
        outDict["centroidDeltaX"] = float(self.centroidDeltaX)
        outDict["centroidDeltaY"] = float(self.centroidDeltaY)
        outDict["gradientX"] = float(self.gradientX)
        outDict["gradientY"] = float(self.gradientY)
        outDict["normalizationFactor"] = float(self.normalizationFactor)

        return outDict

    @classmethod
    def fromTable(cls, tableList):
        """Construct a calibration from a list of tables.

        Parameters
        ----------
        tableList : `list` [`astropy.table.Table`]
            List of table(s) to use to construct the FlatGradient.

        Returns
        -------
        calib : `lsst.ip.isr.FlatGradient`
            The calibration defined in the table(s).
        """
        gradientTable = tableList[0]

        metadata = gradientTable.meta
        inDict = dict()
        inDict["metadata"] = metadata
        inDict["radialSplineNodes"] = np.array(gradientTable[0]["RADIAL_SPLINE_NODES"], dtype=np.float64)
        inDict["radialSplineValues"] = np.array(gradientTable[0]["RADIAL_SPLINE_VALUES"], dtype=np.float64)
        inDict["itlRatio"] = float(gradientTable[0]["ITL_RATIO"][0])
        inDict["centroidX"] = float(gradientTable[0]["CENTROID_X"][0])
        inDict["centroidY"] = float(gradientTable[0]["CENTROID_Y"][0])
        inDict["centroidDeltaX"] = float(gradientTable[0]["CENTROID_DELTA_X"][0])
        inDict["centroidDeltaY"] = float(gradientTable[0]["CENTROID_DELTA_Y"][0])
        inDict["gradientX"] = float(gradientTable[0]["GRADIENT_X"][0])
        inDict["gradientY"] = float(gradientTable[0]["GRADIENT_Y"][0])
        inDict["normalizationFactor"] = float(gradientTable[0]["NORMALIZATION_FACTOR"][0])

        return cls().fromDict(inDict)

    def toTable(self):
        """Construct a list of table(s) containing the FlatGradient data.

        Returns
        -------
        tableList : `list` [`astropy.table.Table`]
            List of tables containing the FlatGradient information.
        """
        tableList = []
        self.updateMetadata()

        catalog = Table(
            data=({
                "RADIAL_SPLINE_NODES": self.radialSplineNodes,
                "RADIAL_SPLINE_VALUES": self.radialSplineValues,
                "ITL_RATIO": np.array([self.itlRatio]),
                "CENTROID_X": np.array([self.centroidX]),
                "CENTROID_Y": np.array([self.centroidY]),
                "CENTROID_DELTA_X": np.array([self.centroidDeltaX]),
                "CENTROID_DELTA_Y": np.array([self.centroidDeltaY]),
                "GRADIENT_X": np.array([self.gradientX]),
                "GRADIENT_Y": np.array([self.gradientY]),
                "NORMALIZATION_FACTOR": np.array([self.normalizationFactor]),
            },)
        )

        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        catalog.meta = outMeta
        tableList.append(catalog)

        return tableList
