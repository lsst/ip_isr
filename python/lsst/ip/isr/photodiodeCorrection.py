#
# LSST Data Management System
# Copyright 2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
from .calibType import IsrCalib

__all__ = ["PhotodiodeCorrection"]


class PhotodiodeCorrection(IsrCalib):
    """Parameter set for linearization.

    These parameters are included in cameraGeom.Amplifier, but
    should be accessible externally to allow for testing.

    Parameters
    ----------
    log : `logging.Logger`, optional
        Logger to handle messages.
    kwargs : `dict`, optional
        Other keyword arguments to pass to the parent init.

    Raises
    ------
    RuntimeError :
        Raised if the supplied table is not 2D, or if the table has fewer
        columns than rows (indicating that the indices are swapped).

    Notes
    -----
    The photodiode correction attributes stored are:

    abscissaCorrections : 'dict' 'float'
        Correction value indexed by exposure pair
    """
    _OBSTYPE = "PHOTODIODE_CORRECTION"
    _SCHEMA = 'PhotodiodeCorrection'
    _VERSION = 1.1

    def __init__(self, table=None, **kwargs):
        self.abscissaCorrections = dict()

        super().__init__(**kwargs)
        self.requiredAttributes.update(['abscissaCorrections'])

    def updateMetadata(self, setDate=False, **kwargs):
        """Update metadata keywords with new values.

        This calls the base class's method after ensuring the required
        calibration keywords will be saved.

        Parameters
        ----------
        setDate : `bool`, optional
            Update the CALIBDATE fields in the metadata to the current
            time. Defaults to False.
        kwargs :
            Other keyword parameters to set in the metadata.
        """

        super().updateMetadata(setDate=setDate, **kwargs)

    @classmethod
    def fromDict(cls, dictionary):
        """Construct a PhotodiodeCorrection from a dictionary of properties.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties.

        Returns
        -------
        calib : `lsst.ip.isr.PhotodiodeCorrection`
            Constructed photodiode data.

        Raises
        ------
        RuntimeError :
            Raised if the supplied dictionary is for a different
            calibration type.
        """
        calib = cls()

        if calib._OBSTYPE != dictionary['metadata']['OBSTYPE']:
            raise RuntimeError(f"Incorrect photodiode correctionsupplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['metadata']['OBSTYPE']}")

        calib.setMetadata(dictionary['metadata'])

        calib.abscissaCorrections = np.array(dictionary['abscissaCorrections'])

        calib.updateMetadata()
        return calib

    def toDict(self):
        """Return a dictionary containing the photodiode correction properties.

        The dictionary should be able to be round-tripped through.
        `fromDict`.

        Returns
        -------
        dictionary : `dict`
            Dictionary of properties.
        """
        self.updateMetadata()

        outDict = {}
        outDict['metadata'] = self.getMetadata()

        outDict['abscissaCorrections'] = self.abscissaCorrections.tolist()

        return outDict

    @classmethod
    def fromTable(cls, tableList):
        """Construct calibration from a list of tables.

        This method uses the `fromDict` method to create the
        calibration after constructing an appropriate dictionary from
        the input tables.

        Parameters
        ----------
        tableList : `list` [`astropy.table.Table`]
            List of tables to use to construct the crosstalk
            calibration.

        Returns
        -------
        calib : `lsst.ip.isr.PhotodiodeCorrection`
            The calibration defined in the tables.
        """
        dataTable = tableList[0]

        metadata = dataTable.meta
        inDict = {}
        inDict['metadata'] = metadata
        inDict['abscissaCorrections'] = dataTable['PD_CORR']

        return cls().fromDict(inDict)

    def toTable(self):
        """Construct a list of tables containing the information in this
        calibration.

        The list of tables should create an identical calibration
        after being passed to this class's fromTable method.

        Returns
        -------
        tableList : `list` [`astropy.table.Table`]
            List of tables containing the photodiode correction
            information.
        """
        self.updateMetadata()
        catalog = Table([{'PD_CORR': self.absciccaCorrections}])
        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        catalog.meta = outMeta

        return([catalog])

    def validate(self):
        """Validate photodiode correction"""
        return
