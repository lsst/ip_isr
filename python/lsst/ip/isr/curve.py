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
QE Curve datatype
"""
import numpy as np
from astropy.table import Table
from lsst.ip.isr import IsrCalib


__all__ = ["QeCurve"]


class QeCurve(IsrCalib):
    """Quantum Efficiency Curve Calibration type

    Parameters
    ----------
    detector : `lsst.afw.cameraGeom.Detector`, optional
        Detector to use to pull coefficients from.
    log : `logging.Logger`, optional
        Log to write messages to.
    **kwargs :
        Parameters to pass to parent constructor.

    Notes
    -----
    """
    _OBSTYPE = 'QE'
    _SCHEMA = 'qe_curve'
    _VERSION = 1.0

    def __init__(self, detector=None, **kwargs):
        super().__init__(**kwargs)
        self.requiredAttributes.update(['wavelength', 'efficiency'])
        self.wavelength = {}
        self.efficiency = {}

        if detector:
            self.fromDetector(detector)

    def fromDetector(self, detector):
        """Set calibration parameters from the detector.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            Detector to use to set parameters from.

        Returns
        -------
        calib : `lsst.ip.isr.QeCurve`
            The calibration constructed from the detector.
        """
        self._detectorId = detector.getId()
        self._detectorName = detector.getName()
        self._detectorSerial = detector.getSerial()

        for amp in detector.getAmplifiers():
            ampName = amp.getName()
            self.wavelength[ampName] = np.array([])
            self.efficiency[ampName] = np.array([])

        self.updateMetadata()
        return self

    @classmethod
    def fromDict(cls, dictionary):
        """Construct a calibration from a dictionary of properties.

        Must be implemented by the specific calibration subclasses.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties.

        Returns
        -------
        calib : `lsst.ip.isr.CalibType`
            Constructed calibration.

        Raises
        ------
        RuntimeError :
            Raised if the supplied dictionary is for a different
            calibration.
        """
        calib = cls()

        if calib._OBSTYPE != dictionary['metadata']['OBSTYPE']:
            raise RuntimeError(f"Incorrect datatype supplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['metadata']['OBSTYPE']}")

        calib.setMetadata(dictionary['metadata'])

        for ampName in dictionary['wavelength'].keys():
            calib.wavelength[ampName] = np.array(dictionary['wavelength'][ampName])
            calib.efficiency[ampName] = np.array(dictionary['efficiency'][ampName])

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
        metadata = self.getMetadata()
        outDict['metadata'] = metadata

        for ampName in self.wavelength.keys():
            outDict['wavelength'][ampName] = self.wavelength[ampName].tolist()
            outDict['efficiency'][ampName] = self.efficiency[ampName].tolist()

        return outDict

    @classmethod
    def fromTable(cls, tableList):
        """Construct calibration from a list of tables.

        This method uses the `fromDict` method to create the
        calibration, after constructing an appropriate dictionary from
        the input tables.

        Parameters
        ----------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables to use to construct the crosstalk
            calibration.

        Returns
        -------
        calib : `lsst.ip.isr.CrosstalkCalib`
            The calibration defined in the tables.

        """
        curveTable = tableList[0]

        metadata = curveTable.meta
        inDict = dict()
        inDict['metadata'] = metadata

        ampNames = set(curveTable['amp_name'])

        for ampName in ampNames:
            indicies = np.argwhere(curveTable['amp_name'] == ampName)
            inDict['wavelength']['ampName'] = curveTable['wavelength'][indicies]
            inDict['efficiency']['ampName'] = curveTable['efficiency'][indicies]

        return cls().fromDict(inDict)

    def toTable(self):
        """Construct a list of tables containing the information in this
        calibration.

        The list of tables should create an identical calibration
        after being passed to this class's fromTable method.

        Returns
        -------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables containing the crosstalk calibration
            information.

        """
        self.updateMetadata()

        ampNames = []
        wavelengths = []
        efficiencies = []

        for ampName in self.wavelength.keys():
            ampNames.extend([ampName for _ in self.wavelength[ampName].tolist()])
            wavelengths.extend(self.wavelength[ampName].tolist())
            efficiencies.extend(self.efficincy[ampName].tolist())

        catalog = Table([{'amp_name': ampNames,
                          'wavelength': wavelengths,
                          'efficiency': efficiencies}])
        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        catalog.meta = outMeta

        return([catalog])
