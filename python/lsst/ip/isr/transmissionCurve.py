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

__all__ = ["IntermediateTransmissionCurve",
           "IntermediateOpticsTransmissionCurve",
           "IntermediateFilterTransmissionCurve",
           "IntermediateSensorTransmissionCurve",
           "IntermediateAtmosphereTransmissionCurve",
           "IntermediateSystemTransmissionCurve",
]

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

    _OBSTYPE = "transmission_curve"
    _SCHEMA = ""
    _VERSION = 1.0

    def __init__(self, filename=None):
        self.data = None
        self.transmissionCurve = None
        self.isSpatiallyConstant = True
        super().__init__()

        # Because we are not persisting this calib as itself, we
        # should skip adding any other attributes.
        self.requiredAttributes.update(['isSpatiallyConstant'])

    def setMetadata(self, metadata):
        # Inherits from lsst.ip.isr.IsrCalib.setMetadata.
        super().setMetadata(metadata)
        if 'OBSTYPE' in metadata:
            _OBSTYPE = metadata['OBSTYPE']

    @classmethod
    def fromTable(cls, tableList):
        """Construct intermediate transmission curve from a list of input
        tables.  Only the first table is used.

        Parameters
        ----------
        tableList : `list` [`astropy.table.Table`]
            List containing input tables.

        Returns
        -------
        calib : `lsst.ip.isr.IntermediateTransmissionCurve`
            The final calibration.
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

        Raises
        ------
        RuntimeError
            This is raised if no table data exists in the calibration,
            if there are array length mismatches, or if the wavelength
            sampling for multi-amp tables differ.
        """
        if self.data is None:
            raise RuntimeError("No table data was found to convert to a transmission curve!")

        throughputKey = None
        if 'wavelength' not in self.data.columns:
            raise RuntimeError("Expected column [wavelength] not found.")
        if 'efficiency' in self.data.columns:
            throughputKey = 'efficiency'
        elif 'throughput' in self.data.columns:
            throughputKey = 'throughput'
        else:
            raise RuntimeError("Expected columns [throughput|efficiency] not found.")

        doAverageCurves = False
        if 'amp_name' in self.data.columns:
            doAverageCurves = True

        # These will be used to construct the
        # ``lsst.afw.image.TransmissionCurve`` object.
        wavelengths = None
        throughput = None

        if doAverageCurves:
            curveStack = None
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
                    curveStack = np.column_stack((curveStack, self.data[indices][throughputKey]))
                else:
                    curveStack = self.data[indices][throughputKey]
            throughput = np.mean(curveStack, 1)

            # This averaging operation has stripped units.
            throughput = throughput * self.data[throughputKey].unit
        else:
            wavelengths = self.data['wavelength']
            throughput = self.data[throughputKey]

        # Convert units:
        # import pdb; pdb.set_trace()
        with cds.enable():
            # These need to be in Angstroms, for consistency.
            wavelengths = wavelengths.to(u.Angstrom).to_value()
            # This is ugly.  Fix.
            if throughput.unit != u.dimensionless_unscaled and throughput.unit != u.UnrecognizedUnit('-'):
                # These need to be fractions, not percent.
                throughput = throughput.to(u.dimensionless_unscaled).to_value()

        self.transmissionCurve = TransmissionCurve.makeSpatiallyConstant(
            throughput.astype(np.float64),
            wavelengths.astype(np.float64),
            throughputAtMin=0.0,
            throughputAtMax=0.0
        )

    def getTransmissionCurve(self):
        if self.transmissionCurve is None and self.data is None:
            raise RuntimeError("No transmission curve data found.")
        if self.transmissionCurve is None:
            self.setTransmissionCurveRepresentation()
        return self.transmissionCurve

    def writeFits(self, outputFilename):
        """Write the transmission curve data to a file.

        Parameters
        ----------
        outputFilename : `str`
            Destination filename.

        Returns
        -------
        outputFilename : `str`
            The output filename actually used.

        Raises
        ------
        RuntimeError
            Raised if no transmission curve can be created.
        """
        if self.transmissionCurve is None and self.data is None:
            raise RuntimeError("No transmission curve data found.")
        if self.transmissionCurve is None:
            self.setTransmissionCurveRepresentation()

        return self.transmissionCurve.writeFits(outputFilename)


class IntermediateSensorTransmissionCurve(IntermediateTransmissionCurve):
    _OBSTYPE = 'transmission_sensor'


class IntermediateFilterTransmissionCurve(IntermediateTransmissionCurve):
    _OBSTYPE = 'transmission_filter'


class IntermediateOpticsTransmissionCurve(IntermediateTransmissionCurve):
    _OBSTYPE = 'transmission_optics'


class IntermediateAtmosphereTransmissionCurve(IntermediateTransmissionCurve):
    _OBSTYPE = 'transmission_atmosphere'


class IntermediateSystemTransmissionCurve(IntermediateTransmissionCurve):
    _OBSTYPE = 'transmission_system'
