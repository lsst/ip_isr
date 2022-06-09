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
Photodiode storage class.
"""

__all__ = ["PhotodiodeCalib"]

import numpy as np
from astropy.table import Table

from lsst.ip.isr import IsrCalib


class PhotodiodeCalib(IsrCalib):
    """Independent current measurements from photodiode for linearity
    calculations.

    Parameters
    ----------
    timeSamples : `list` or `numpy.ndarray`, optional
        List of samples the photodiode was measured at.
    currentSamples : `list` or `numpy.ndarray`, optional
        List of current measurements at each time sample.
    **kwargs :
        Additional parameters. These will be passed to the parent
        constructor with the exception of:

        ``"integrationMethod"``
            Name of the algorithm to use to integrate the current
            samples. Allowed values are ``DIRECT_SUM``,
            ``TRIMMED_SUM``, and ``CHARGE_SUM`` (`str`).
        ``"currentScale"``
            Scale factor to apply to the current samples for the
            ``CHARGE_SUM`` integration method. A typical value
            would be `-1`, to flip the sign of the integrated charge.
    """

    _OBSTYPE = 'PHOTODIODE'
    _SCHEMA = 'Photodiode'
    _VERSION = 1.0

    def __init__(self, timeSamples=None, currentSamples=None, **kwargs):
        if timeSamples is not None and currentSamples is not None:
            if len(timeSamples) != len(currentSamples):
                raise RuntimeError(f"Inconsitent vector lengths: {len(timeSamples)} vs {len(currentSamples)}")
            else:
                self.timeSamples = np.array(timeSamples).ravel()
                self.currentSamples = np.array(currentSamples).ravel()
        else:
            self.timeSamples = np.array([]).ravel()
            self.currentSamples = np.array([]).ravel()

        super().__init__(**kwargs)

        if 'integrationMethod' in kwargs:
            self.integrationMethod = kwargs.pop('integrationMethod')
        else:
            self.integrationMethod = 'DIRECT_SUM'

        if 'currentScale' in kwargs:
            self.currentScale = kwargs.pop('currentScale')
        else:
            self.currentScale = 1.0

        if 'day_obs' in kwargs:
            self.updateMetadata(day_obs=kwargs['day_obs'])
        if 'seq_num' in kwargs:
            self.updateMetadata(seq_num=kwargs['seq_num'])

        self.requiredAttributes.update(['timeSamples', 'currentSamples', 'integrationMethod'])

    @classmethod
    def fromDict(cls, dictionary):
        """Construct a PhotodiodeCalib from a dictionary of properties.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties.

        Returns
        -------
        calib : `lsst.ip.isr.PhotodiodeCalib`
            Constructed photodiode data.

        Raises
        ------
        RuntimeError
            Raised if the supplied dictionary is for a different
            calibration type.
        """
        calib = cls()

        if calib._OBSTYPE != dictionary['metadata']['OBSTYPE']:
            raise RuntimeError(f"Incorrect photodiode supplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['metadata']['OBSTYPE']}")

        calib.setMetadata(dictionary['metadata'])

        calib.timeSamples = np.array(dictionary['timeSamples']).ravel()
        calib.currentSamples = np.array(dictionary['currentSamples']).ravel()
        calib.integrationMethod = dictionary.get('integrationMethod', "DIRECT_SUM")

        calib.updateMetadata()
        return calib

    def toDict(self):
        """Return a dictionary containing the photodiode properties.

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

        outDict['timeSamples'] = self.timeSamples.tolist()
        outDict['currentSamples'] = self.currentSamples.tolist()

        outDict['integrationMethod'] = self.integrationMethod

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
        calib : `lsst.ip.isr.PhotodiodeCalib`
            The calibration defined in the tables.
        """
        dataTable = tableList[0]

        metadata = dataTable.meta
        inDict = {}
        inDict['metadata'] = metadata
        if 'OBSTYPE' not in metadata:
            inDict['metadata']['OBSTYPE'] = cls._OBSTYPE
        inDict['integrationMethod'] = metadata.pop('INTEGRATION_METHOD', 'DIRECT_SUM')

        for key in ('TIME', 'Elapsed Time', ):
            if key in dataTable.columns:
                inDict['timeSamples'] = dataTable[key]
        for key in ('CURRENT', 'Signal', ):
            if key in dataTable.columns:
                inDict['currentSamples'] = dataTable[key]

        return cls().fromDict(inDict)

    def toTable(self):
        """Construct a list of tables containing the information in this
        calibration.

        The list of tables should create an identical calibration
        after being passed to this class's fromTable method.

        Returns
        -------
        tableList : `list` [`astropy.table.Table`]
            List of tables containing the photodiode calibration
            information.
        """
        self.updateMetadata()
        catalog = Table([{'TIME': self.timeSamples,
                          'CURRENT': self.currentSamples}])
        inMeta = self.getMetadata().toDict()
        outMeta = {k: v for k, v in inMeta.items() if v is not None}
        outMeta.update({k: "" for k, v in inMeta.items() if v is None})
        outMeta['INTEGRATION_METHOD'] = self.integrationMethod
        catalog.meta = outMeta

        return [catalog]

    @classmethod
    def readTwoColumnPhotodiodeData(cls, filename):
        """Construct a PhotodiodeCalib by reading the simple column format.

        Parameters
        ----------
        filename : `str`
            File to read samples from.

        Returns
        -------
        calib : `lsst.ip.isr.PhotodiodeCalib`
            The calibration defined in the file.
        """
        import os.path

        rawData = np.loadtxt(filename, dtype=[('time', 'float'), ('current', 'float')])

        basename = os.path.basename(filename)
        cleaned = os.path.splitext(basename)[0]
        _, _, day_obs, seq_num = cleaned.split("_")

        return cls(timeSamples=rawData['time'], currentSamples=rawData['current'],
                   day_obs=int(day_obs), seq_num=int(seq_num))

    def integrate(self):
        """Integrate the current.

        Raises
        ------
        RuntimeError
            Raised if the integration method is not known.
        """
        if self.integrationMethod == 'DIRECT_SUM':
            return self.integrateDirectSum()
        elif self.integrationMethod == 'TRIMMED_SUM':
            return self.integrateTrimmedSum()
        elif self.integrationMethod == 'CHARGE_SUM':
            return self.integrateChargeSum()
        else:
            raise RuntimeError(f"Unknown integration method {self.integrationMethod}")

    def integrateDirectSum(self):
        """Integrate points.

        This uses numpy's trapezoidal integrator.

        Returns
        -------
        sum : `float`
            Total charge measured.
        """
        return np.trapz(self.currentSamples, x=self.timeSamples)

    def integrateTrimmedSum(self):
        """Integrate points with a baseline level subtracted.

        This uses numpy's trapezoidal integrator.

        Returns
        -------
        sum : `float`
            Total charge measured.

        See Also
        --------
        lsst.eotask.gen3.eoPtc
        """
        currentThreshold = ((max(self.currentSamples) - min(self.currentSamples))/5.0
                            + min(self.currentSamples))
        lowValueIndices = np.where(self.currentSamples < currentThreshold)
        baseline = np.median(self.currentSamples[lowValueIndices])
        return np.trapz(self.currentSamples - baseline, self.timeSamples)

    def integrateChargeSum(self):
        """For this method, the values in .currentSamples are actually the
        integrated charge values as measured by the ammeter for each
        sampling interval.  We need to do a baseline subtraction,
        based on the charge values when the LED is off, then sum up
        the corrected signals.

        Returns
        -------
        sum : `float`
            Total charge measured.
        """
        dt = self.timeSamples[1:] - self.timeSamples[:-1]
        # The .currentSamples values are the current integrals over
        # the interval preceding the current time stamp, so omit the
        # first value.
        charge = self.currentScale*self.currentSamples[1:]
        # The current per interval to use for baseline subtraction
        # without assuming all of the dt values are the same:
        current = charge/dt
        # To determine the baseline current level, exclude points with
        # signal levels > 5% of the maximum (measured relative to the
        # overall minimum), and extend that selection 2 entries on
        # either side to avoid otherwise low-valued points that sample
        # the signal ramp and which should not be included in the
        # baseline estimate.
        dy = np.max(current) - np.min(current)
        signal, = np.where(current > dy/20. + np.min(current))
        imin = signal[0] - 2
        imax = signal[-1] + 2
        bg = np.concatenate([np.arange(0, imin), np.arange(imax, len(current))])
        bg_current = np.sum(charge[bg])/np.sum(dt[bg])
        # Return the background-subtracted total charge.
        return np.sum(charge - bg_current*dt)
