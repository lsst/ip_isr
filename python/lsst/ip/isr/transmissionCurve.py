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

__all__ = ["IntermediateTransmissionCurve"]

import numpy as np
from astropy import units as u
from astropy.units import cds

from lsst.afw.image import TransmissionCurve
from .calibType import IsrCalib


class IntermediateTransmissionCurve(IsrCalib):
    """Definition for the TransmissionCurve format used as inputs.

    Parameters
    ----------
    filename : `str`
        Filename of a transmission curve dataset.
    """

    _OBSTYPE = "TRANSMISSION CURVE"
    _SCHEMA = ""
    _VERSION = 1.0

    def __init__(self, filename=None):
        super().__init__()
        if filename:
            self.readText(filename)

        self.requiredAttributes.update([''])
        # self.setTransmissionCurveRepresentation()
        self.isSpatiallyConstant = True

    @classmethod
    def fromTable(cls, tableList):
        """
        """
        calib = cls()

        metadata = tableList[0].meta
        calib.setMetadata(metadata)
        calib.updateMetadata()
        calib.data = tableList[0]

        calib.setTransmissionCurveRepresentation()
        return calib

    def setTransmissionCurveRepresentation(self):
        """Construct transmission curve representation from the data that was
        read.
        """
        doAverageCurves = False
        if 'wavelength' not in self.data.columns or 'throughput' not in self.data.columns:
            raise RuntimeError("Expected columns not found.")
        if 'amp_name' in self.data.columns:
            doAverageCurves = True

        # These will be used to construct the
        # ``lsst.afw.image.TransmissionCurve`` object.
        wavelengths = None
        throughput = None

        if doAverageCurves:
            curveStack = []
            amplifierNames = set(self.data['amp_name'])
            comparisonIndices = np.where(self.data['amp_name'] == next(iter(amplifierNames)))
            wavelengths = self.data[comparisonIndices]['wavelength']

            for amplifier in amplifierNames:
                indices = np.where(self.data['amp_name'] == amplifier)
                if len(self.data[indices]) != len(self.data[comparisonIndices]):
                    raise RuntimeError("Incompatible lengths in average.")
                if not np.array_equal(self.data[indices]['wavelength'], wavelengths):
                    raise RuntimeError("Mismatch in wavelength samples.")

                if curveStack is not None:
                    curveStack = np.column_curveStack((curveStack, self.data[indices]['efficiency']))
                else:
                    curveStack = self.data[indices]['efficiency']
            throughput = np.mean(curveStack, 1)
        else:
            wavelengths = self.data['wavelength']
            throughput = self.data['throughput']

        # Convert units:
        with cds.enable():
            wavelengths = wavelengths.to(u.Angstrom)  # These need to be in Angstroms, for consistency.
            throughput = throughput.to(u.dimensionless_unscaled)  # These need to be fractions, not percent.
        tt = throughput.to_value().astype(np.float64)
        ww = wavelengths.to_value().astype(np.float64)
        self.transmissionCurve = TransmissionCurve.makeSpatiallyConstant(tt, ww,
                                                                         throughputAtMin=0.0,
                                                                         throughputAtMax=0.0)

    def writeFits(self, outputFilename):
        return self.transmissionCurve.writeFits(outputFilename)
