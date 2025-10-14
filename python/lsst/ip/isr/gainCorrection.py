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
Gain correction storage class.
"""

__all__ = ["GainCorrection"]

from astropy.table import Table
import numpy as np

from lsst.ip.isr import IsrCalib


class GainCorrection(IsrCalib):
    """Gain correction parameters.

    Parameters
    ----------
    ampNames : `list` [`str`]
        List of amplifier names.
    gainAdjustments : `list` [`float`]
        List of gain adjustment parameters.
    **kwargs :
        Additional parameters.
    """

    _OBSTYPE = "gainCorrection"
    _SCHEMA = "GainCorrection"
    _VERSION = 1.0

    def __init__(self, ampNames=[], gainAdjustments=[], **kwargs):
        if len(ampNames) != len(gainAdjustments):
            raise ValueError("Number of ampNames must be the same as number of gainAdjustments.")

        self.ampNames = ampNames
        self.gainAdjustments = np.asarray(gainAdjustments)

        super().__init__(**kwargs)
        self.requiredAttributes.update(["ampNames", "gainAdjustments"])

        self.updateMetadata(setCalibInfo=True, setCalibId=True, **kwargs)

    def setParameters(
        self,
        *,
        ampNames=[],
        gainAdjustments=[],
    ):
        """Set the parameters for the gain correction model.

        Parameters
        ----------
        ampNames : `list` [`str`]
            List of amplifier names.
        gainAdjustments : `list` [`float`]
            List of gain adjustment parameters.
        """
        if len(ampNames) != len(gainAdjustments):
            raise ValueError("Number of ampNames must be the same as number of gainAdjustments.")

        self.ampNames = ampNames
        self.gainAdjustments = gainAdjustments

    @classmethod
    def fromDict(cls, dictionary):
        """Construct a GainCorrection from a dictionary of properties.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties.

        Returns
        -------
        calib : `lsst.ip.isr.GainCorrection`
            Constructed calibration.
        """
        calib = cls()

        calib.setMetadata(dictionary["metadata"])

        calib.ampNames = dictionary["ampNames"]
        calib.gainAdjustments = np.asarray(dictionary["gainAdjustments"])

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

        outDict["ampNames"] = self.ampNames
        outDict["gainAdjustments"] = self.gainAdjustments.tolist()

        return outDict

    @classmethod
    def fromTable(cls, tableList):
        """Construct a calibration from a list of tables.

        Parameters
        ----------
        tableList : `list` [`astropy.table.Table`]
            List of table(s) to use to construct the GainCorrection.

        Returns
        -------
        calib : `lsst.ip.isr.GainCorrection`
            The calibration defined in the table(s).
        """
        gainCorrectionTable = tableList[0]

        inDict = dict()

        inDict["metadata"] = gainCorrectionTable.meta
        inDict["ampNames"] = list(gainCorrectionTable["AMP_NAME"])
        inDict["gainAdjustments"] = np.asarray(gainCorrectionTable["GAIN_ADJUSTMENT"], dtype=np.float64)

        return cls().fromDict(inDict)

    def toTable(self):
        """Construct a list of table(s) containing the GainCorrection data.

        Returns
        -------
        tableList : `list` [`astropy.table.Table`]
            List of tables containing the GainCorrection information.
        """
        tableList = []
        self.updateMetadata()

        catalog = Table()
        catalog["AMP_NAME"] = self.ampNames
        catalog["GAIN_ADJUSTMENT"] = self.gainAdjustments

        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        catalog.meta = outMeta
        tableList.append(catalog)

        return tableList

    def correctGains(self, gains, exposure=None):
        """Correct a dictionary of gains (in place).

        Parameters
        ----------
        gains : `dict` [`str`, `float`]
            Array of gains to correct.
        exposure : `lsst.afw.image.Exposure`, optional
            Exposure with additional metadata for correction.
        """
        for i, ampName in enumerate(self.ampNames):
            gains[ampName] *= self.gainAdjustments[i]
