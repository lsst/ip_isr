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
Intrinsic Zernikes storage class.
"""

__all__ = ["IntrinsicZernikes"]

import numpy as np
from astropy import units as u
from astropy.table import Table
from scipy.interpolate import LinearNDInterpolator

from lsst.ip.isr import IsrCalib


class IntrinsicZernikes(IsrCalib):
    """Intrinsic Zernike coefficients.

    Stores Zernike wavefront-error coefficients sampled at a set of
    focal-plane field angles.  At query time the coefficients are
    interpolated to an arbitrary field position.

    Field angles are expressed in the Camera Coordinate System (CCS), also
    known as the Engineering Diagram Coordinate System.  See
    `LSE-349 <https://ls.st/LSE-349>`_ for the definition.

    Parameters
    ----------
    table : `astropy.table.Table`, optional
        Source table.  Must contain columns:

        ``"x"``
            Field x positions (in CCS) with angular units (e.g. ``u.deg``).
        ``"y"``
            Field y positions (in CCS) with angular units (e.g. ``u.deg``).
        ``"Z{j}"``
            One column per Noll index *j*, with length units
            (e.g. ``u.um``).

    Attributes
    ----------
    field_x : `numpy.ndarray`
        CCS x field positions in degrees for all sample points,
        shape ``(n_points,)``.
    field_y : `numpy.ndarray`
        CCS y field positions in degrees for all sample points,
        shape ``(n_points,)``.
    noll_indices : `numpy.ndarray`
        Noll indices of the stored Zernike terms, shape ``(n_zernikes,)``.
    values : `numpy.ndarray`
        Zernike coefficients in microns, shape
        ``(n_points, n_zernikes)``.
    interpolator : `scipy.interpolate.LinearNDInterpolator` or `None`
        Interpolator built from ``field_x``, ``field_y``, and
        ``values``.  ``None`` until the calibration is populated.
    """

    _OBSTYPE = "INTRINSIC_ZERNIKES"
    _SCHEMA = "Intrinsic Zernikes"
    _VERSION = 1.0

    def __init__(self, table=None, **kwargs):
        self.field_x = np.array([])
        self.field_y = np.array([])
        self.values = np.array([])
        self.noll_indices = np.array([])
        self.interpolator = None

        super().__init__(**kwargs)

        if table is not None:
            self.field_x = table["x"].to("deg").value
            self.field_y = table["y"].to("deg").value
            zcols = [col for col in table.colnames if col.startswith("Z")]
            self.noll_indices = np.array(sorted([int(col[1:]) for col in zcols]))
            zks = np.column_stack(
                [
                    table[col].to("um").value for col in zcols
                ]
            )
            self.values = zks
            self._createInterpolator()

        self.requiredAttributes.update(["field_x", "field_y", "values", "noll_indices"])

    def _createInterpolator(self):
        self.interpolator = LinearNDInterpolator(
            np.column_stack((self.field_x, self.field_y)),
            self.values
        )

    @classmethod
    def fromDict(cls, dictionary):
        """Construct an IntrinsicZernikes from dictionary of properties.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties.

        Returns
        -------
        calib : `lsst.ip.isr.IntrinsicZernikes`
            Constructed calibration.

        Raises
        ------
        RuntimeError
            Raised if the supplied dictionary is for a different
            calibration type.
        """
        calib = cls()

        if calib._OBSTYPE != dictionary["metadata"]["OBSTYPE"]:
            raise RuntimeError(
                f"Incorrect intrinsic zernikes supplied. "
                f"Expected {calib._OBSTYPE}, found {dictionary['metadata']['OBSTYPE']}"
            )

        calib.setMetadata(dictionary["metadata"])
        calib.field_x = np.array(dictionary["field_x"])
        calib.field_y = np.array(dictionary["field_y"])
        calib.values = np.array(dictionary["values"])
        calib.noll_indices = np.array(dictionary["noll_indices"])
        calib._createInterpolator()

        calib.updateMetadata()
        return calib

    def toDict(self):
        """Return a dictionary containing the calibration properties.

        The dictionary should be able to be round-tripped through
        `fromDict`.

        Returns
        -------
        dictionary : `dict`
            Dictionary of properties.
        """
        self.updateMetadata()

        outDict = {}
        outDict["metadata"] = self.getMetadata()
        outDict["field_x"] = self.field_x.tolist()
        outDict["field_y"] = self.field_y.tolist()
        outDict["values"] = self.values.tolist()
        outDict["noll_indices"] = self.noll_indices.tolist()

        return outDict

    @classmethod
    def fromTable(cls, tableList):
        """Construct calibration from a list of tables.

        Parameters
        ----------
        tableList : `list` [`astropy.table.Table`]
            List of tables to use to construct the intrinsic zernikes
            calibration.

        Returns
        -------
        calib : `lsst.ip.isr.IntrinsicZernikes`
            The calibration defined in the tables.
        """
        table = tableList[0]
        calib = cls(table=table)
        calib.setMetadata(table.meta)
        calib.updateMetadata()
        return calib

    def toTable(self):
        """Construct a list of tables containing the information in this
        calibration.

        The list of tables should be able to be round-tripped through
        `fromTable`.

        Returns
        -------
        tableList : `list` [`astropy.table.Table`]
            List of tables containing the intrinsic zernikes calibration
            information.
        """
        self.updateMetadata()

        data = {
            "x": self.field_x * u.deg,
            "y": self.field_y * u.deg,
        }
        for i, j in enumerate(self.noll_indices):
            data[f"Z{j}"] = self.values[:, i] * u.um

        table = Table(data)

        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        table.meta = outMeta

        return [table]

    def getIntrinsicZernikes(self, field_x, field_y, noll_indices=None):
        """
        Get the intrinsic Zernike coefficients at a given field position.

        Parameters
        ----------
        field_x : `array-like`
            CCS x-field positions in degrees.
        field_y : `array-like`
            CCS y-field positions in degrees.
        noll_indices : `list` [`int`], optional
            List of Noll indices to return. If None, return all.

        Returns
        -------
        zernikes : `array-like`
            Array of Zernike coefficient values in microns corresponding to the
            requested Noll indices and field positions.
        """
        if noll_indices is None:
            noll_indices = self.noll_indices

        point = np.array([field_x, field_y]).T
        interpolated_values = self.interpolator(point)

        noll_indices = np.array(noll_indices)
        noll_mask = np.isin(self.noll_indices, noll_indices)
        return interpolated_values[..., noll_mask]
