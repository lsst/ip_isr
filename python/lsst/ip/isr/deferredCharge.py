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

# import lsst.afw.math as afwMath
from lsst.afw.cameraGeom import ReadoutCorner
from lsst.pex.config import Config, Field
from lsst.pipe.base import Task
from .isrFunctions import gainContext
from .calibType import IsrCalib

import scipy.interpolate as interp


__all__ = ('DeferredChargeConfig', 'DeferredChargeTask', 'SerialTrap', 'DeferredChargeCalib')


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
    trap_type : `str`
        Type of trap capture to use.
    coeffs : `list` [`float`]
        Coefficients for the capture process.
    """

    def __init__(self, size, emission_time, pixel, trap_type, coeffs):
        if size < 0.0:
            raise ValueError('Trap size must be greater than or equal to 0.')
        self.size = size

        if emission_time <= 0.0:
            raise ValueError('Emission time must be greater than 0.')
        if np.isnan(emission_time):
            raise ValueError('Emission time must be real-valued, not NaN')
        self.emission_time = emission_time

        self.pixel = int(pixel)

        self.trap_type = trap_type
        self.coeffs = coeffs

        if self.trap_type == 'spline':
            centers, values = np.split(np.array(self.coeffs), 2)
            self.interp = interp.interp1d(centers, values)

        self._trap_array = None
        self._trapped_charge = None

    def __eq__(self, other):
        # A trap is equal to another trap if all of the initialization
        # parameters are equal.  All other properties are only filled
        # during use, and are not persisted into the calibration.
        if self.size != other.size:
            return False
        if self.emission_time != other.emission_time:
            return False
        if self.pixel != other.pixel:
            return False
        if self.trap_type != other.trap_type:
            return False
        if self.coeffs != other.coeffs:
            return False
        return True

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
        if self.trap_type == 'linear':
            scaling = self.coeffs[0]
            return np.minimum(self.size, pixel_signals*scaling)
        elif self.trap_type == 'logistic':
            f0, k = (self.coeffs[0], self.coeffs[1])
            return self.size/(1.+np.exp(-k*(pixel_signals-f0)))
        elif self.trap_type == 'spline':
            return self.interp(pixel_signals)


class DeferredChargeCalib(IsrCalib):
    """Calibration containing deferred charge/CTI parameters.

    Parameters
    ----------
    detector : `lsst.afw.cameraGeom.Detector`, optional
        Detector to use for metadata properties.
    **kwargs :
        Additional parameters to pass to parent constructor.

    Notes
    -----
    The charge transfer inefficiency attributes stored are:

    driftScale : `dict` [`str`, `float`]
        A dictionary, keyed by amplifier name, of the local electronic
        offset drift scale parameter, A_L in Snyder+2021.
    decayTime : `dict` [`str`, `float`]
        A dictionary, keyed by amplifier name, of the local electronic
        offset decay time, \tau_L in Snyder+2021.
    globalCti : `dict` [`str`, `float`]
        A dictionary, keyed by amplifier name, of the mean global CTI
        paramter, b in Snyder+2021.
    serialTraps : `dict` [`str`, `lsst.ip.isr.SerialTrap`]
        A dictionary, keyed by amplifier name, containing a single
        serial trap for each amplifier.
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
        """Construct a calibration from a dictionary of properties.

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
            raise RuntimeError(f"Incorrect CTI supplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['metadata']['OBSTYPE']}")

        calib.setMetadata(dictionary['metadata'])

        calib.driftScale = dictionary['driftScale']
        calib.decayTime = dictionary['decayTime']
        calib.globalCti = dictionary['globalCti']

        for ampName in dictionary['serialTraps']:
            ampTraps = dictionary['serialTraps'][ampName]
            calib.serialTraps[ampName] = SerialTrap(ampTraps['size'], ampTraps['emissionTime'],
                                                    ampTraps['pixel'], ampTraps['trap_type'],
                                                    ampTraps['coeffs'])
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
        outDict['metadata'] = self.getMetadata()

        outDict['driftScale'] = self.driftScale
        outDict['decayTime'] = self.decayTime
        outDict['globalCti'] = self.globalCti

        outDict['serialTraps'] = {}
        for ampName in self.serialTraps:
            ampTrap = {'size': self.serialTraps[ampName].size,
                       'emissionTime': self.serialTraps[ampName].emission_time,
                       'pixel': self.serialTraps[ampName].pixel,
                       'trap_type': self.serialTraps[ampName].trap_type,
                       'coeffs': self.serialTraps[ampName].coeffs}
            outDict['serialTraps'][ampName] = ampTrap

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
            calibration.  Two tables are expected in this list, the
            first containing the per-amplifier CTI parameters, and the
            second containing the parameters for serial traps.

        Returns
        -------
        calib : `lsst.ip.isr.CrosstalkCalib`
            The calibration defined in the tables.
        """
        ampTable = tableList[0]

        inDict = {}
        inDict['metadata'] = ampTable.meta

        amps = ampTable['AMPLIFIER']
        driftScale = ampTable['DRIFT_SCALE']
        decayTime = ampTable['DECAY_TIME']
        globalCti = ampTable['GLOBAL_CTI']

        inDict['driftScale'] = {amp: value for amp, value in zip(amps, driftScale)}
        inDict['decayTime'] = {amp: value for amp, value in zip(amps, decayTime)}
        inDict['globalCti'] = {amp: value for amp, value in zip(amps, globalCti)}

        inDict['serialTraps'] = {}
        trapTable = tableList[1]

        amps = trapTable['AMPLIFIER']
        sizes = trapTable['SIZE']
        emissionTimes = trapTable['EMISSION_TIME']
        pixels = trapTable['PIXEL']
        trap_type = trapTable['TYPE']
        coeffs = trapTable['COEFFS']

        for index, amp in enumerate(amps):
            ampTrap = {}
            ampTrap['size'] = sizes[index]
            ampTrap['emissionTime'] = emissionTimes[index]
            ampTrap['pixel'] = pixels[index]
            ampTrap['trap_type'] = trap_type[index]
            ampTrap['coeffs'] = np.array(coeffs[index])[~np.isnan(coeffs[index])].tolist()

            inDict['serialTraps'][amp] = ampTrap

        return cls.fromDict(inDict)

    def toTable(self):
        """Construct a list of tables containing the information in this
        calibration.

        The list of tables should create an identical calibration
        after being passed to this class's fromTable method.

        Returns
        -------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables containing the crosstalk calibration
            information.  Two tables are generated for this list, the
            first containing the per-amplifier CTI parameters, and the
            second containing the parameters for serial traps.
        """
        tableList = []
        self.updateMetadata()

        ampList = []
        driftScale = []
        decayTime = []
        globalCti = []

        for amp in self.driftScale.keys():
            ampList.append(amp)
            driftScale.append(self.driftScale[amp])
            decayTime.append(self.decayTime[amp])
            globalCti.append(self.globalCti[amp])

        ampTable = Table({'AMPLIFIER': ampList,
                          'DRIFT_SCALE': driftScale,
                          'DECAY_TIME': decayTime,
                          'GLOBAL_CTI': globalCti,
                          })

        ampTable.meta = self.getMetadata().toDict()
        tableList.append(ampTable)

        ampList = []
        sizeList = []
        timeList = []
        pixelList = []
        typeList = []
        coeffList = []

        # Get maximum coeff length
        maxCoeffLength = 0
        for trap in self.serialTraps.values():
            maxCoeffLength = np.maximum(maxCoeffLength, len(trap.coeffs))

        # Pack and pad the end of the coefficients with NaN values.
        for amp, trap in self.serialTraps.items():
            ampList.append(amp)
            sizeList.append(trap.size)
            timeList.append(trap.emission_time)
            pixelList.append(trap.pixel)
            typeList.append(trap.trap_type)

            coeffs = trap.coeffs
            if len(coeffs) != maxCoeffLength:
                coeffs = np.pad(coeffs, (0, maxCoeffLength - len(coeffs)),
                                constant_values=np.nan).tolist()
            coeffList.append(coeffs)

        trapTable = Table({'AMPLIFIER': ampList,
                           'SIZE': sizeList,
                           'EMISSION_TIME': timeList,
                           'PIXEL': pixelList,
                           'TYPE': typeList,
                           'COEFFS': coeffList})

        tableList.append(trapTable)

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
    useGains = Field(
        dtype=bool,
        doc="If true, scale by the gain.",
        default=False,
    )
    zeroUnusedPixels = Field(
        dtype=bool,
        doc="If true, set serial prescan and parallel overscan to zero before correction.",
        default=False,
    )


class DeferredChargeTask(Task):
    """Task to correct an exposure for charge transfer inefficiency.

    This uses the methods described by Snyder et al. 2021, Journal of
    Astronimcal Telescopes, Instruments, and Systems, 7,
    048002. doi:10.1117/1.JATIS.7.4.048002 (Snyder+21).
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
        if self.config.useGains:
            if gains is None:
                gains = {amp.getName(): amp.getGain() for amp in detector.getAmplifiers()}
        else:
            gains = {amp.getName(): 1.0 for amp in detector.getAmplifiers()}

        with gainContext(exposure, image, True, gains):
            for amp in detector.getAmplifiers():
                ampName = amp.getName()

                ampImage = image[amp.getRawBBox()]
                if self.config.zeroUnusedPixels:
                    # We don't apply overscan subtraction, so zero these
                    # out for now.
                    ampImage[amp.getRawParallelOverscanBBox()].getArray()[:, :] = 0.0
                    ampImage[amp.getRawSerialPrescanBBox()].getArray()[:, :] = 0.0

                # The algorithm expects that the readout corner is in
                # the lower left corner.  Flip it to be so:

                ampData = self.flipData(ampImage.getArray(), amp)

                if ctiCalib.driftScale[ampName] > 0.0:
                    correctedAmpData = self.local_offset_inverse(ampData,
                                                                 ctiCalib.driftScale[ampName],
                                                                 ctiCalib.decayTime[ampName],
                                                                 self.config.nPixelOffsetCorrection)
                else:
                    correctedAmpData = ampData.copy()

                correctedAmpData = self.local_trap_inverse(correctedAmpData,
                                                           ctiCalib.serialTraps[ampName],
                                                           ctiCalib.globalCti[ampName],
                                                           self.config.nPixelTrapCorrection)

                # Undo flips here.  The method is symmetric.
                correctedAmpData = self.flipData(correctedAmpData, amp)
                image[amp.getBBox()].getArray()[:, :] = correctedAmpData[:, :]

        return exposure

    @staticmethod
    def flipData(ampData, amp):
        """Flip data array such that readout corner is at lower-left.

        Parameters
        ----------
        ampData : `np.ndarray`, (nx, ny)
            Image data to flip.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Amplifier to get readout corner information.

        Returns
        -------
        ampData : `np.ndarray`, (nx, ny)
            Flipped image data.
        """
        X_FLIP = {ReadoutCorner.LL: False,
                  ReadoutCorner.LR: True,
                  ReadoutCorner.UL: False,
                  ReadoutCorner.UR: True}
        Y_FLIP = {ReadoutCorner.LL: False,
                  ReadoutCorner.LR: False,
                  ReadoutCorner.UL: True,
                  ReadoutCorner.UR: True}

        if X_FLIP(amp.getReadoutCorner()):
            ampData = np.fliplr(ampData)
        if Y_FLIP(amp.getReadoutCorner()):
            ampData = np.flipud(ampData)

        return ampData

    @staticmethod
    def local_offset_inverse(inputArr, drift_scale, decay_time, num_previous_pixels=15):
        """Remove CTI effects from local offsets.

        This implements equation 10 of Snyder+21.  For an image with
        CTI, s'(m, n), the correction factor is equal to the maximum
        value of the set of:
            {A_L s'(m, n - j) exp(-j t / \tau_L)}_j=0^jmax

        Parameters
        ----------
        inputArr : `np.ndarray`, (nx, ny)
            Input image data to correct.
        drift_scale : `float`
            Drift scale (Snyder+21 A_L value) to use in correction.
        decay_time : `float`
            Decay time (Snyder+21 \tau_L) of the correction.
        num_previous_pixels : `int`, optional
            Number of previous pixels to use for correction.  As the
            CTI has an exponential decay, this essentially truncates
            the correction where that decay scales the input charge to
            near zero.

        Returns
        -------
        outputArr : `np.ndarray`, (nx, ny)
            Corrected image data.
        """
        r = np.exp(-1/decay_time)
        Ny, Nx = inputArr.shape

        # j = 0 term:
        offset = np.zeros((num_previous_pixels, Ny, Nx))
        offset[0, :, :] = drift_scale*np.maximum(0, inputArr)

        # j = 1..jmax terms:
        for n in range(1, num_previous_pixels):
            offset[n, :, n:] = drift_scale*np.maximum(0, inputArr[:, :-n])*(r**n)

        Linv = np.amax(offset, axis=0)
        outputArr = inputArr - Linv

        return outputArr

    @staticmethod
    def local_trap_inverse(inputArr, trap, global_cti=0.0, num_previous_pixels=6):
        """Apply localized trapping inverse operator to pixel signals.

        This implements equation 13 of Snyder+21.  For an image with
        CTI, s'(m, n), the correction factor is equal to the maximum
        value of the set of:
            {A_L s'(m, n - j) exp(-j t / \tau_L)}_j=0^jmax

        Parameters
        ----------
        inputArr : `np.ndarray`, (nx, ny)
            Input image data to correct.
        trap : `lsst.ip.isr.SerialTrap`
            Serial trap describing the capture and release of charge.
        global_cti: `float`
            Mean charge transfer inefficiency, b from Snyder+21.
        num_previous_pixels : `int`, optional
            Number of previous pixels to use for correction.

        Returns
        -------
        outputArr : `np.ndarray`, (nx, ny)
            Corrected image data.

        """
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

        outputArr = inputArr - a*T

        return outputArr
