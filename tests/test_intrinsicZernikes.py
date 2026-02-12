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
import tempfile

import numpy as np
from astropy.table import Table
import astropy.units as u

import lsst.utils.tests

from lsst.ip.isr import IntrinsicZernikes


class IntrinsicZernikesTestCase(lsst.utils.tests.TestCase):
    """Test the IntrinsicZernikes calibration class."""

    def setUp(self):
        """Create test data for intrinsic Zernikes."""
        self.field_x = np.array([-1.0, 0.0, 1.0])
        self.field_y = np.array([-0.5, 0.0, 0.5])

        # Create test Zernike coefficients for Noll indices 4, 5, 6
        # (defocus, astigmatism)
        self.noll_indices = np.array([4, 5, 6])

        # Create values: shape (n_y, n_x, n_zernikes)
        rng = np.random.default_rng(seed=57721)
        self.values = rng.normal(
            scale=0.1,
            size=(len(self.field_y), len(self.field_x), len(self.noll_indices))
        )  # microns

        # Create an astropy table in the format expected by __init__
        x_grid, y_grid = np.meshgrid(self.field_x, self.field_y)
        self.inputTable = Table()
        self.inputTable["x"] = (x_grid.ravel() * u.deg)
        self.inputTable["y"] = (y_grid.ravel() * u.deg)

        # Add Zernike columns
        for i, noll in enumerate(self.noll_indices):
            self.inputTable[f"Z{noll}"] = (self.values[:, :, i].ravel() * u.um)

        # Create the calibration object
        self.calib = IntrinsicZernikes(table=self.inputTable)

    def test_initialization_with_table(self):
        """Test that IntrinsicZernikes initializes correctly from a table."""
        np.testing.assert_array_equal(self.calib.field_x, self.field_x)
        np.testing.assert_array_equal(self.calib.field_y, self.field_y)
        np.testing.assert_array_equal(self.calib.noll_indices, self.noll_indices)
        np.testing.assert_array_equal(self.calib.values, self.values)
        self.assertIsNotNone(self.calib.interpolator)

    def test_metadata(self):
        """Test that metadata is properly set."""
        metadata = self.calib.getMetadata()
        self.assertEqual(metadata["OBSTYPE"], "INTRINSIC_ZERNIKES")
        self.assertEqual(metadata["INTRINSIC_ZERNIKES_SCHEMA"], "Intrinsic Zernikes")
        self.assertEqual(metadata["INTRINSIC_ZERNIKES_VERSION"], 1.0)

    def test_dict_roundtrip(self):
        """Test round-tripping through dictionary."""
        newCalib = IntrinsicZernikes.fromDict(self.calib.toDict())
        self.assertEqual(newCalib, self.calib)

    def test_table_roundtrip(self):
        """Test round-tripping through table."""
        newCalib = IntrinsicZernikes.fromTable(self.calib.toTable())
        self.assertEqual(newCalib, self.calib)

    def test_yaml_roundtrip(self):
        """Test round-tripping through YAML file."""
        with tempfile.NamedTemporaryFile(suffix='.yaml', delete=False) as f:
            filename = f.name

        try:
            self.calib.writeText(filename)
            newCalib = IntrinsicZernikes.readText(filename)
            self.assertEqual(newCalib, self.calib)
        finally:
            import os
            if os.path.exists(filename):
                os.remove(filename)

    def test_fits_roundtrip(self):
        """Test round-tripping through FITS file."""
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
            filename = f.name

        try:
            self.calib.writeFits(filename)
            newCalib = IntrinsicZernikes.readFits(filename)
            self.assertEqual(newCalib, self.calib)
        finally:
            import os
            if os.path.exists(filename):
                os.remove(filename)

    def test_fromDict_wrong_obstype(self):
        """Test that fromDict raises error for wrong OBSTYPE."""
        outDict = self.calib.toDict()
        outDict["metadata"]["OBSTYPE"] = "WRONG_TYPE"

        with self.assertRaises(RuntimeError) as context:
            IntrinsicZernikes.fromDict(outDict)

        self.assertIn("Incorrect intrinsic zernikes supplied", str(context.exception))
        self.assertIn("INTRINSIC_ZERNIKES", str(context.exception))
        self.assertIn("WRONG_TYPE", str(context.exception))

    def test_getIntrinsicZernikes(self):
        """Test interpolation of Zernike coefficients."""
        # Test at a grid point
        field_x_test = 0.0
        field_y_test = 0.0

        zernikes = self.calib.getIntrinsicZernikes(field_x_test, field_y_test)

        # In this case, we're on a grid point
        self.assertFloatsEqual(zernikes, self.values[1, 1, :])

        # Test with specific Noll indices
        zernikes_subset = self.calib.getIntrinsicZernikes(
            field_x_test, field_y_test,
            noll_indices=[4]
        )
        self.assertFloatsEqual(zernikes_subset, zernikes[:, 0])

    def test_getIntrinsicZernikes_array(self):
        """Test interpolation with array inputs."""
        field_x_test = np.array([0.0, 0.5])
        field_y_test = np.array([0.0, 0.25])

        zernikes = self.calib.getIntrinsicZernikes(field_x_test, field_y_test)
        z0 = self.calib.getIntrinsicZernikes(field_x_test[0], field_y_test[0])
        z1 = self.calib.getIntrinsicZernikes(field_x_test[1], field_y_test[1])
        np.testing.assert_array_almost_equal(zernikes[[0]], z0)
        np.testing.assert_array_almost_equal(zernikes[[1]], z1)

        # Should return shape (n_points, n_zernikes)
        self.assertEqual(zernikes.shape, (2, len(self.noll_indices)))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
