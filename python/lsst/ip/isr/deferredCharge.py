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
import numpy as np
from astropy.table import Table

from lsst.pex.config import Config, Field
from lsst.pipe.base import Task
from .isrFunctions import gainContext
from .calibType import IsrCalib


__all__ = ('DeferredChargeConfig', 'DeferredChargeTask')


class SerialTrap():
    """Represents a serial register trap.

    Parameters
    ----------
    size : `float`
        Size of the charge trap, in electrons.
    emission_time : `float`
        Trap emission time constant, in inverse transfers.
    pixel : `int`
        Serial pixel location of the trap.
    trapType : `str`
        Type of trap capture to use.
    coeffs : `list` [`float`]
        Coefficients for the capture process.
    """

    def __init__(self, size, emission_time, pixel, trapType, coeffs):
        if size < 0.0:
            raise ValueError('Trap size must be greater than or equal to 0.')
        self.size = size

        if emission_time <= 0.0:
            raise ValueError('Emission time must be greater than 0.')
        if np.isnan(emission_time):
            raise ValueError('Emission time must be real-valued, not NaN')
        self.emission_time = emission_time

        self.pixel = int(pixel)

        self.trapType = trapType
        self.coeffs = coeffs

        self._trap_array = None
        self._trapped_charge = None

    @property
    def trap_array(self):
        return self._trap_array

    @property
    def trapped_charge(self):
        return self._trapped_charge

    def initialize(self, ny, nx, prescan_width):
        """Initialize trapping arrays for simulated readout.

        Parameters
        ----------
        ny : `int`
            Number of rows to simulate.
        nx : `int`
            Number of columns to simulate.
        prescan_width : `int`
            Additional transfers due to prescan.
        """
        if self.pixel > nx+prescan_width:
            raise ValueError('Trap location {0} must be less than {1}'.format(self.pixel,
                                                                              nx+prescan_width))

        self._trap_array = np.zeros((ny, nx+prescan_width))
        self._trap_array[:, self.pixel] = self.size
        self._trapped_charge = np.zeros((ny, nx+prescan_width))

    def release_charge(self):
        """Release charge through exponential decay.

        Returns
        -------
        released_charge : `float`
            Charge released.
        """
        released_charge = self._trapped_charge*(1-np.exp(-1./self.emission_time))
        self._trapped_charge -= released_charge

        return released_charge

    def trap_charge(self, free_charge):
        """Perform charge capture using a logistic function.

        Parameters
        ----------
        free_charge : `float`
            Charge available to be trapped.

        Returns
        -------
        captured_charge : `float`
            Amount of charge actually trapped.
        """
        captured_charge = (np.clip(self.capture(free_charge), self.trapped_charge, self._trap_array)
                           - self.trapped_charge)
        self._trapped_charge += captured_charge

        return captured_charge

    def capture(self, pixel_signals):
        """Trap capture function.

        Parameters
        ----------
        pixel_signals : `list` [`float`]
            Input pixel values.

        Returns
        -------
        captured_charge : `list` [`float`]
            Amount of charge captured from each pixel.
        """
        if self.trapType == 'linear':
            scaling = self.coeffs[0]
            return np.minimum(self.size, pixel_signals*scaling)
        elif self.trapType == 'logistic':
            f0, k = (self.coeffs[0], self.coeffs[1])
            return self.size/(1.+np.exp(-k*(pixel_signals-f0)))
        elif self.trapType == 'spline':
            raise NotImplementedError("Spline currently not implemented.")
            # super().__init__(200000., emission_time, pixel)
            # self.f = interpolant


class DeferredChargeCalib(IsrCalib):
    """Calibration of deferred charge/CTI parameters.
    """
    _OBSTYPE = 'CTI'
    _SCHEMA = 'Deferred Charge'
    _VERSION = 1.0

    def __init__(self, detector=None, **kwargs):
        self.driftScale = {}
        self.decayTime = {}
        self.globalCti = {}
        self.serialTraps = {}

        super().__init__(**kwargs)
        self.requiredAttributes.update(['driftScale', 'decayTime', 'globalCti', 'serialTraps'])

    @classmethod
    def fromDict(cls, dictionary):
        calib = cls()

        if calib._OBSTYPE != dictionary['metadata']['OBSTYPE']:
            raise RuntimeError(f"Incorrect CTI supplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['metadata']['OBSTYPE']}")

        calib.setMetadata(dictionary['metadata'])

        calib.driftScale = dictionary['driftScale']
        calib.decayTime = dictionary['decayTime']
        calib.globalCti = dictionary['globalCti']

        for ampName in dictionary['serialTraps']:
            ampTraps = dictionary['serialTraps'][ampName]
            calib.serialTraps[ampName] = SerialTrap(ampTraps['size'], ampTraps['emissionTime'],
                                                    ampTraps['pixel'], ampTraps['trapType'],
                                                    ampTraps['coeffs'])
        calib.updateMetadata()
        return calib

    def toDict(self):
        self.updateMetadata()
        outDict = {}
        outDict['metadata'] = self.getMetadata()

        outDict['driftScale'] = self.driftScale
        outDict['decayTime'] = self.decayTime
        outDict['globalCti'] = self.globalCti

        outDict['serialTraps'] = {}
        for ampName in self.serialTraps:
            ampTrap = {'size': self.serialTraps[ampName].size,
                       'emissionTime': self.serialTraps[ampName].emission_time,
                       'pixel': self.serialTraps[ampName].pixel,
                       'trapType': self.serialTraps[ampName].trapType,
                       'coeffs': self.serialTraps[ampName].coeffs}
            outDict['serialTraps'][ampName] = ampTrap

        return outDict

    @classmethod
    def fromTable(cls, tableList):
        ampTable = tableList[0]

        inDict = {}
        inDict['metadata'] = ampTable.meta

        amps = ampTable['AMPLIFIER']
        driftScales = ampTable['DRIFT_SCALE']
        decayTimes = ampTable['DECAY_TIME']
        globalCti = ampTable['GLOBAL_CTI']

        inDict['driftScale'] = {amp: value for amp, value in zip(amps, driftScales)}
        inDict['decayTimes'] = {amp: value for amp, value in zip(amps, decayTimes)}
        inDict['globalCti'] = {amp: value for amp, value in zip(amps, globalCti)}

        trapTable = tableList[1]

        amps = trapTable['AMPLIFIER']
        sizes = trapTable['SIZE']
        emissionTimes = trapTable['EMISSION_TIME']
        pixels = trapTable['PIXEL']
        trapType = trapTable['TYPE']
        coeffs = trapTable['coeffs']

        for index, amp in enumerate(amps):
            ampTrap = {}
            ampTrap['size'] = sizes[index]
            ampTrap['emissionTime'] = emissionTimes[index]
            ampTrap['pixel'] = pixels[index]
            ampTrap['trapType'] = trapType[index]
            ampTrap['coeffs'] = coeffs[index]

            inDict['serialTraps'][amp] = ampTrap

        return cls.fromDict(inDict)

    def toTable(self):
        tableList = []
        self.updateMetadata()

        ampList = []
        driftScales = []
        decayTimes = []
        globalCti = []

        for amp in self.driftScale.keys():
            ampList.append(amp)
            driftScales.append(self.driftScale[amp])
            decayTimes.append(self.decayTimes[amp])
            globalCti.append(self.globalCti[amp])

        ampTable = Table({'AMPLIFIER': ampList,
                          'DRIFT_SCALE': driftScales,
                          'DECAY_TIME': decayTimes,
                          'GLOBAL_CTI': globalCti,
                          })

        ampTable.meta = self.getMetadata().toDict()
        tableList.append(ampTable)

        # serial Traps??
        return tableList


class DeferredChargeConfig(Config):
    """Settings for deferred charge correction.
    """
    nPixelOffsetCorrection = Field(
        dtype=int,
        doc="Number of prior pixels to CZW DOC.",
        default=15,
    )
    nPixelTrapCorrection = Field(
        dtype=int,
        doc="Number of prior pixels to CZW DOC.",
        default=6,
    )


class DeferredChargeTask(Task):
    """CZW DOC.
    """
    ConfigClass = DeferredChargeConfig
    _DefaultName = 'isrDeferredCharge'

    def run(self, exposure, overscans, ctiCalib, gains=None):
        """Correct deferred charge/CTI issues.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to correct the deferred charge on.
        gains : `dict` [`str`, `float`]
            A dictionary, keyed by amplifier name, of the gains to
            use.  If gains is None, the nominal gains in the amplifier
            object are used
        ctiCalib : `lsst.ip.isr.DeferredChargeCalib`
            Calibration object containing the charge transfer
            inefficiency model.

        Returns
        -------
        exposure : `lsst.afw.image.Exposure`
            The corrected exposure.
        """
        image = exposure.getMaskedImage().getImage()
        detector = exposure.getDetector()
        with gainContext(exposure, image, True, gains):
            for amp in detector.getAmplifiers():
                ampImage = image[amp.getBBox()]

                if ctiCalib.driftScale[amp] > 0.0:
                    correctedAmpImage = self.local_offset_inverse(ampImage.getArray(),
                                                                  ctiCalib.driftScale[amp],
                                                                  ctiCalib.decayTime[amp],
                                                                  self.config.nPixelOffsetCorrection)
                else:
                    correctedAmpImage = ampImage.clone()

                correctedAmpImage = self.local_trap_inverse(correctedAmpImage.getArray(),
                                                            ctiCalib.serialTraps[amp],
                                                            ctiCalib.globalCti[amp],
                                                            self.config.nPixelTrapCorrection)
            image.getArray()[:, :] = correctedAmpImage[:, :]

        return exposure

    @staticmethod
    def local_offset_inverse(inputArr, scale, decay_time, num_previous_pixels=4):
        """
        """
        r = np.exp(-1/decay_time)
        Ny, Nx = inputArr.shape

        offset = np.zeros((num_previous_pixels, Ny, Nx))
        offset[0, :, :] = scale*np.maximum(0, inputArr)

        for n in range(1, num_previous_pixels):
            offset[n, :, n:] = scale*np.maximum(0, inputArr[:, :-n])*(r**n)

        L = np.amax(offset, axis=0)

        # This probably should just return L, the correction
        outputArr = inputArr - L

        return outputArr

    @staticmethod
    def local_trap_inverse(inputArr, trap, global_cti=0.0, num_previous_pixels=4):
        """Apply localized trapping inverse operator to pixel signals."""

        Ny, Nx = inputArr.shape
        a = 1 - global_cti
        r = np.exp(-1/trap.emission_time)

        # Estimate trap occupancies during readout
        trap_occupancy = np.zeros((num_previous_pixels, Ny, Nx))
        for n in range(num_previous_pixels):
            trap_occupancy[n, :, n+1:] = trap.capture(np.maximum(0, inputArr))[:, :-(n+1)]*(r**n)
        trap_occupancy = np.amax(trap_occupancy, axis=0)

        # Estimate captured charge
        C = trap.capture(np.maximum(0, inputArr)) - trap_occupancy*r
        C[C < 0] = 0.

        # Estimate released charge
        R = np.zeros(inputArr.shape)
        R[:, 1:] = trap_occupancy[:, 1:]*(1-r)
        T = R - C

        # This probably should just return a*T, the correction amount
        outputArr = inputArr - a*T

        return outputArr
