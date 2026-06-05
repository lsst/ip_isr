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
    interpolated to an arbitrary field position provided in CCS.

    The coefficients are stored in two coordinate systems, each as an
    independent set of sample points and values:

    - the Camera Coordinate System (CCS), which corresponds to the
     focal plane heights and any other CCS contribution and
    - the Optical Coordinate System (OCS), corresponding to
    the intrinsics, which are by nature defined in the optical
    coordinates.

    See `LSE-349 <https://ls.st/LSE-349>`_ for the definitions.

    The CCS sample points and values are stored on the un-suffixed
    ``field_x``, ``field_y`` and ``values`` attributes (and serialized
    under those same keys).  These names are kept unchanged from version
    1, which only stored the CCS system, so that version 1 calibrations
    round-trip unchanged.  The OCS system is stored alongside on the
    ``*_ocs`` attributes/keys.

    Parameters
    ----------
    table : `astropy.table.Table`, optional
        Source table in the CCS.  Must contain columns:

        ``"x"``
            Field x positions with angular units (e.g. ``u.deg``).
        ``"y"``
            Field y positions with angular units (e.g. ``u.deg``).
        ``"Z{j}"``
            One column per Noll index *j*, with length units
            (e.g. ``u.um``).
    table_ocs : `astropy.table.Table`, optional
        Source table in the OCS, with the same column layout as
        ``table``.

    Attributes
    ----------
    field_x, field_y : `numpy.ndarray`
        CCS x/y field positions in degrees for all sample points,
        shape ``(n_points_ccs,)``.
    field_x_ocs, field_y_ocs : `numpy.ndarray`
        OCS x/y field positions in degrees for all sample points,
        shape ``(n_points_ocs,)``.
    noll_indices : `numpy.ndarray`
        Noll indices of the stored Zernike terms, shape ``(n_zernikes,)``.
    values, values_ocs : `numpy.ndarray`
        Zernike coefficients in microns for the CCS and OCS sample
        points, shape ``(n_points, n_zernikes)``.
    interpolator, interpolator_ocs : `scipy.interpolate.LinearNDInterpolator`
    or `None`
        Interpolators built from the CCS and OCS sample points and
        values.  ``None`` until the corresponding system is populated.

    Version 1.1 adds the OCS coordinate system alongside the CCS system
    stored by version 1.
    """

    _OBSTYPE = "INTRINSIC_ZERNIKES"
    _SCHEMA = "Intrinsic Zernikes"
    _VERSION = 1.1

    def __init__(self, table=None, table_ocs=None, **kwargs):
        # CCS uses the un-suffixed names (field_x/field_y/values/
        # interpolator) for backwards compatibility with version 1, which
        # only stored the CCS system under those names.
        self.field_x = np.array([])
        self.field_y = np.array([])
        self.values = np.array([])
        self.field_x_ocs = np.array([])
        self.field_y_ocs = np.array([])
        self.values_ocs = np.array([])
        self.noll_indices = np.array([])
        self.noll_indices_ocs = np.array([])
        self.interpolator = None
        self.interpolator_ocs = None

        super().__init__(**kwargs)

        if table is not None:
            (self.field_x, self.field_y, self.values, self.noll_indices) = (
                self._unpackTable(table)
            )
            self.interpolator = self._makeInterpolator(
                self.field_x, self.field_y, self.values
            )
        if table_ocs is not None:
            (self.field_x_ocs, self.field_y_ocs, self.values_ocs, self.noll_indices_ocs) = (
                self._unpackTable(table_ocs)
            )
            self.interpolator_ocs = self._makeInterpolator(
                self.field_x_ocs, self.field_y_ocs, self.values_ocs
            )

        self.requiredAttributes.update(
            [
                "field_x",
                "field_y",
                "values",
                "field_x_ocs",
                "field_y_ocs",
                "values_ocs",
                "noll_indices",
                "noll_indices_ocs",
            ]
        )

    @staticmethod
    def _unpackTable(table):
        """Unpack a source table into field positions, values, and Noll
        indices.

        Parameters
        ----------
        table : `astropy.table.Table`
            Source table with ``"x"``, ``"y"``, and ``"Z{j}"`` columns.

        Returns
        -------
        field_x, field_y : `numpy.ndarray`
            Field positions in degrees.
        values : `numpy.ndarray`
            Zernike coefficients in microns, shape
            ``(n_points, n_zernikes)`` ordered by ascending Noll index.
        noll_indices : `numpy.ndarray`
            Sorted Noll indices.
        """
        field_x = table["x"].to("deg").value
        field_y = table["y"].to("deg").value
        zcols = [col for col in table.colnames if col.startswith("Z")]
        noll_indices = np.array(sorted(int(col[1:]) for col in zcols))
        values = np.column_stack([table[f"Z{j}"].to("um").value for j in noll_indices])
        return field_x, field_y, values, noll_indices

    @staticmethod
    def _makeInterpolator(field_x, field_y, values):
        """Build a field-position interpolator, or `None` if there are no
        sample points."""
        if np.asarray(field_x).size == 0:
            return None
        return LinearNDInterpolator(np.column_stack((field_x, field_y)), values)

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
        # CCS keys (field_x/field_y/values) are always present, including
        # in version 1 dictionaries.  The OCS keys are optional, so version
        # 1 dictionaries (CCS only) load with an empty OCS system.
        calib.field_x = np.array(dictionary["field_x"])
        calib.field_y = np.array(dictionary["field_y"])
        calib.values = np.array(dictionary["values"])
        calib.field_x_ocs = np.array(dictionary.get("field_x_ocs", []))
        calib.field_y_ocs = np.array(dictionary.get("field_y_ocs", []))
        calib.values_ocs = np.array(dictionary.get("values_ocs", []))
        calib.noll_indices = np.array(dictionary["noll_indices"])
        calib.noll_indices_ocs = np.array(dictionary["noll_indices_ocs"])
        calib.interpolator = cls._makeInterpolator(
            calib.field_x, calib.field_y, calib.values
        )
        calib.interpolator_ocs = cls._makeInterpolator(
            calib.field_x_ocs, calib.field_y_ocs, calib.values_ocs
        )

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
        outDict["field_x_ocs"] = self.field_x_ocs.tolist()
        outDict["field_y_ocs"] = self.field_y_ocs.tolist()
        outDict["values_ocs"] = self.values_ocs.tolist()
        outDict["noll_indices"] = self.noll_indices.tolist()
        outDict["noll_indices_ocs"] = self.noll_indices_ocs.tolist()

        return outDict

    @classmethod
    def fromTable(cls, tableList):
        """Construct calibration from a list of tables.

        Parameters
        ----------
        tableList : `list` [`astropy.table.Table`]
            List of tables to use to construct the intrinsic zernikes
            calibration.  Each table is dispatched to the CCS or OCS
            coordinate system according to its ``coord_sys`` metadata
            entry (defaulting to ``"CCS"``).  A version 1 single-table
            calibration, which has no ``coord_sys`` entry, is therefore
            read as the CCS system.

        Returns
        -------
        calib : `lsst.ip.isr.IntrinsicZernikes`
            The calibration defined in the tables.
        """
        tables = {}
        for table in tableList:
            coord_sys = table.meta.get("coord_sys", "CCS")
            if coord_sys not in ("CCS", "OCS"):
                raise RuntimeError(
                    f"Invalid coordinate system {coord_sys} in table metadata; "
                    f"expected 'CCS' or 'OCS'"
                )
            tables[coord_sys] = table

        calib = cls(table=tables.get("CCS"), table_ocs=tables.get("OCS", None))
        # ``coord_sys`` is a per-table annotation used only to dispatch each
        # table above; drop it so it does not leak into the calibration
        # metadata (which must match across a toTable/fromTable round-trip).
        meta = dict(tableList[0].meta)
        meta.pop("coord_sys", None)
        calib.setMetadata(meta)
        calib.updateMetadata()
        return calib

    def toTable(self):
        """Construct a list of tables containing the information in this
        calibration.

        One table is produced per populated coordinate system.  The CCS
        table is always emitted; the OCS table is only emitted when the
        OCS system holds sample points, so a CCS-only calibration (e.g.
        one read from a version 1 file) round-trips to a single table,
        exactly as in version 1.  The list of tables should be able to be
        round-tripped through `fromTable`.

        Returns
        -------
        tableList : `list` [`astropy.table.Table`]
            List of tables containing the intrinsic zernikes calibration
            information, one per populated coordinate system.
        """
        self.updateMetadata()

        inMeta = self.getMetadata().toDict()
        baseMeta = {k: v for k, v in inMeta.items() if v is not None}
        baseMeta.update({k: "" for k, v in inMeta.items() if v is None})

        systems = [("CCS", self.field_x, self.field_y, self.values)]
        if np.asarray(self.field_x_ocs).size > 0:
            systems.append(("OCS", self.field_x_ocs, self.field_y_ocs, self.values_ocs))

        tableList = []
        for coord_sys, field_x, field_y, values in systems:
            data = {
                "x": field_x * u.deg,
                "y": field_y * u.deg,
            }
            for i, j in enumerate(self.noll_indices):
                column = values[:, i] if values.ndim == 2 else np.array([])
                data[f"Z{j}"] = column * u.um

            table = Table(data)
            meta = dict(baseMeta)
            meta["coord_sys"] = coord_sys
            table.meta = meta
            tableList.append(table)

        return tableList

    def getIntrinsicZernikes(
        self, field_x, field_y, rotTelPos=0.0, noll_indices=None
    ):
        """
        Get the intrinsic Zernike coefficients at a given field position.

        The returned coefficients are the sum of the CCS contribution
        (heights_ccs), interpolated at the requested field position,
        and the OCS contribution (measured_intrinsics),
        interpolated at the field position rotated by
        ``rotTelPos``.  For calibrations that only store one
        coordinate system (e.g. version 1 files, which only carry CCS),
        the missing OCS contribution is simply omitted from the sum.

        Parameters
        ----------
        field_x : `array-like`
            x-field positions in degrees (CCS).
        field_y : `array-like`
            y-field positions in degrees (CCS).
        rotTelPos : `float`, optional
            Rotation angle in degrees applied to the query point before
            interpolating the OCS contribution.  Defaults to 0.
        noll_indices : `list` [`int`], optional
            List of Noll indices to return. If None, return all.

        Returns
        -------
        zernikes : `array-like`
            Array of Zernike coefficient values in microns corresponding to the
            requested Noll indices and field positions.
        """
        if noll_indices is None:
            noll_indices = sorted(
                set(self.noll_indices).intersection(self.noll_indices_ocs)
            )
        noll_indices = np.array(noll_indices)
        noll_mask = np.isin(self.noll_indices, noll_indices)

        field_x = np.asarray(field_x)
        field_y = np.asarray(field_y)

        total = None
        point = np.array([field_x, field_y]).T
        total = self.interpolator(point)

        if self.interpolator_ocs is not None:
            # Rotate the query point into the OCS frame before interpolating.
            theta = np.deg2rad(rotTelPos)
            cos_a, sin_a = np.cos(theta), np.sin(theta)
            x_ocs = cos_a * field_x - sin_a * field_y
            y_ocs = sin_a * field_x + cos_a * field_y
            point_ocs = np.array([x_ocs, y_ocs]).T
            ocs_values = self.interpolator_ocs(point_ocs)
            total = ocs_values if total is None else total + ocs_values

        return total[..., noll_mask]
