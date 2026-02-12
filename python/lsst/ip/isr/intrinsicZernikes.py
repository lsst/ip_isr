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
from astropy.table import Table
from scipy.interpolate import RegularGridInterpolator

from lsst.ip.isr import IsrCalib


class IntrinsicZernikes(IsrCalib):
    """
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
            self.field_x = np.unique(table["x"].to("deg").value)
            self.field_y = np.unique(table["y"].to("deg").value)
            zcols = [col for col in table.colnames if col.startswith("Z")]
            self.noll_indices = np.array(sorted([int(col[1:]) for col in zcols]))
            zks = np.column_stack(
                [
                    table[col].to("um").value for col in zcols
                ]
            )
            self.values = zks.reshape(self.field_y.size, self.field_x.size, -1)
            self._createInterpolator()

        self.requiredAttributes.update(["field_x", "field_y", "values", "noll_indices"])

    def _createInterpolator(self):
        self.interpolator = RegularGridInterpolator(
            (self.field_y, self.field_x),
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
        inDict = {}
        inDict["metadata"] = table.meta
        inDict["field_x"] = table["FIELD_X"][0]
        inDict["field_y"] = table["FIELD_Y"][0]
        inDict["noll_indices"] = table["NOLL_INDICES"][0]
        inDict["values"] = table["VALUES"][0].reshape(
            inDict["field_y"].size,
            inDict["field_x"].size,
            -1
        )
        return cls.fromDict(inDict)

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
        tableList = []
        self.updateMetadata()

        table = Table({
            "FIELD_X": [self.field_x],
            "FIELD_Y": [self.field_y],
            "NOLL_INDICES": [self.noll_indices],
            "VALUES": [self.values.ravel()],
        })

        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        table.meta = outMeta

        tableList.append(table)
        return tableList

    def getIntrinsicZernikes(self, field_x, field_y, noll_indices=None):
        """
        Get the intrinsic Zernike coefficients at a given field position.

        Parameters
        ----------
        field_x : `array-like`
            The x-field positions in degrees.
        field_y : `array-like`
            The y-field positions in degrees.
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

        point = np.array([field_y, field_x]).T
        interpolated_values = self.interpolator(point)

        noll_indices = np.array(noll_indices)
        noll_mask = np.isin(self.noll_indices, noll_indices)
        return interpolated_values[..., noll_mask]
