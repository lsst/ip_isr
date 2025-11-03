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

__all__ = ('DeferredChargeConfig',
           'DeferredChargeTask',
           'SerialTrap',
           'OverscanModel',
           'SimpleModel',
           'SimulatedModel',
           'SegmentSimulator',
           'FloatingOutputAmplifier',
           'DeferredChargeCalib',
           )

import copy
import numpy as np
import warnings
from astropy.table import Table

from lsst.afw.cameraGeom import ReadoutCorner
from lsst.pex.config import Config, Field
from lsst.pipe.base import Task
from .isrFunctions import gainContext
from .calibType import IsrCalib

import scipy.interpolate as interp


class SerialTrap():
    """Represents a serial register trap.

    Parameters
    ----------
    size : `float`
        Size of the charge trap, in electrons.
    emission_time : `float`
        Trap emission time constant, in inverse transfers.
    pixel : `int`
        Serial pixel location of the trap, including the prescan.
    trap_type : `str`
        Type of trap capture to use.  Should be one of ``linear``,
        ``logistic``, or ``spline``.
    coeffs : `list` [`float`]
        Coefficients for the capture process.  Linear traps need one
        coefficient, logistic traps need two, and spline based traps
        need to have an even number of coefficients that can be split
        into their spline locations and values.

    Raises
    ------
    ValueError
        Raised if the specified parameters are out of expected range.
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

        if int(pixel) != pixel:
            raise ValueError('Fraction value for pixel not allowed.')
        self.pixel = int(pixel)

        self.trap_type = trap_type
        self.coeffs = coeffs

        if self.trap_type not in ('linear', 'logistic', 'spline'):
            raise ValueError('Unknown trap type: %s', self.trap_type)

        if self.trap_type == 'spline':
            # Note that ``spline`` is actually a piecewise linear interpolation
            # in the model and the application, and not a true spline.
            centers, values = np.split(np.array(self.coeffs, dtype=np.float64), 2)
            # Ensure all NaN values are stripped out
            values = values[~np.isnan(centers)]
            centers = centers[~np.isnan(centers)]
            centers = centers[~np.isnan(values)]
            values = values[~np.isnan(values)]
            self.interp = interp.interp1d(
                centers,
                values,
                bounds_error=False,
                fill_value=(values[0], values[-1]),
            )

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

        Raises
        ------
        ValueError
            Raised if the trap falls outside of the image.
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

        Raises
        ------
        RuntimeError
            Raised if the trap type is invalid.
        """
        if self.trap_type == 'linear':
            scaling = self.coeffs[0]
            return np.minimum(self.size, pixel_signals*scaling)
        elif self.trap_type == 'logistic':
            f0, k = (self.coeffs[0], self.coeffs[1])
            return self.size/(1.+np.exp(-k*(pixel_signals-f0)))
        elif self.trap_type == 'spline':
            return self.interp(pixel_signals)
        else:
            raise RuntimeError(f"Invalid trap capture type: {self.trap_type}.")


class OverscanModel:
    """Base class for handling model/data fit comparisons.
    This handles all of the methods needed for the lmfit Minimizer to
    run.
    """

    @staticmethod
    def model_results(params, signal, num_transfers, start=1, stop=10):
        """Generate a realization of the overscan model, using the specified
        fit parameters and input signal.

        Parameters
        ----------
        params : `lmfit.Parameters`
            Object containing the model parameters.
        signal : `np.ndarray`, (nMeasurements)
            Array of image means.
        num_transfers : `int`
            Number of serial transfers that the charge undergoes.
        start : `int`, optional
            First overscan column to fit. This number includes the
            last imaging column, and needs to be adjusted by one when
            using the overscan bounding box.
        stop : `int`, optional
            Last overscan column to fit. This number includes the
            last imaging column, and needs to be adjusted by one when
            using the overscan bounding box.

            Returns
        -------
        results : `np.ndarray`, (nMeasurements, nCols)
            Model results.
        """
        raise NotImplementedError("Subclasses must implement the model calculation.")

    def loglikelihood(self, params, signal, data, error, *args, **kwargs):
        """Calculate log likelihood of the model.

        Parameters
        ----------
        params : `lmfit.Parameters`
            Object containing the model parameters.
        signal : `np.ndarray`, (nMeasurements)
            Array of image means.
        data : `np.ndarray`, (nMeasurements, nCols)
            Array of overscan column means from each measurement.
        error : `float`
            Fixed error value.
        *args :
            Additional position arguments.
        **kwargs :
            Additional keyword arguments.

        Returns
        -------
        logL : `float`
            The log-likelihood of the observed data given the model
            parameters.
        """
        model_results = self.model_results(params, signal, *args, **kwargs)

        inv_sigma2 = 1.0/(error**2.0)
        diff = model_results - data

        return -0.5*(np.sum(inv_sigma2*(diff)**2.))

    def negative_loglikelihood(self, params, signal, data, error, *args, **kwargs):
        """Calculate negative log likelihood of the model.

        Parameters
        ----------
        params : `lmfit.Parameters`
            Object containing the model parameters.
        signal : `np.ndarray`, (nMeasurements)
            Array of image means.
        data : `np.ndarray`, (nMeasurements, nCols)
            Array of overscan column means from each measurement.
        error : `float`
            Fixed error value.
        *args :
            Additional position arguments.
        **kwargs :
            Additional keyword arguments.

        Returns
        -------
        negativelogL : `float`
            The negative log-likelihood of the observed data given the
            model parameters.
        """
        ll = self.loglikelihood(params, signal, data, error, *args, **kwargs)

        return -ll

    def rms_error(self, params, signal, data, error, *args, **kwargs):
        """Calculate RMS error between model and data.

        Parameters
        ----------
        params : `lmfit.Parameters`
            Object containing the model parameters.
        signal : `np.ndarray`, (nMeasurements)
            Array of image means.
        data : `np.ndarray`, (nMeasurements, nCols)
            Array of overscan column means from each measurement.
        error : `float`
            Fixed error value.
        *args :
            Additional position arguments.
        **kwargs :
            Additional keyword arguments.

        Returns
        -------
        rms : `float`
            The rms error between the model and input data.
        """
        model_results = self.model_results(params, signal, *args, **kwargs)

        diff = model_results - data
        rms = np.sqrt(np.mean(np.square(diff)))

        return rms

    def difference(self, params, signal, data, error, *args, **kwargs):
        """Calculate the flattened difference array between model and data.

        Parameters
        ----------
        params : `lmfit.Parameters`
            Object containing the model parameters.
        signal : `np.ndarray`, (nMeasurements)
            Array of image means.
        data : `np.ndarray`, (nMeasurements, nCols)
            Array of overscan column means from each measurement.
        error : `float`
            Fixed error value.
        *args :
            Additional position arguments.
        **kwargs :
            Additional keyword arguments.

        Returns
        -------
        difference : `np.ndarray`, (nMeasurements*nCols)
            The rms error between the model and input data.
        """
        model_results = self.model_results(params, signal, *args, **kwargs)
        diff = (model_results-data).flatten()

        return diff


class SimpleModel(OverscanModel):
    """Simple analytic overscan model."""

    @staticmethod
    def model_results(params, signal, num_transfers, start=1, stop=10):
        """Generate a realization of the overscan model, using the specified
        fit parameters and input signal.

        Parameters
        ----------
        params : `lmfit.Parameters`
            Object containing the model parameters.
        signal : `np.ndarray`, (nMeasurements)
            Array of image means.
        num_transfers : `int`
            Number of serial transfers that the charge undergoes.
        start : `int`, optional
            First overscan column to fit. This number includes the
            last imaging column, and needs to be adjusted by one when
            using the overscan bounding box.
        stop : `int`, optional
            Last overscan column to fit. This number includes the
            last imaging column, and needs to be adjusted by one when
            using the overscan bounding box.

        Returns
        -------
        res : `np.ndarray`, (nMeasurements, nCols)
            Model results.
        """
        v = params.valuesdict()
        v['cti'] = 10**v['ctiexp']

        # Adjust column numbering to match DM overscan bbox.
        start += 1
        stop += 1

        x = np.arange(start, stop+1)
        res = np.zeros((signal.shape[0], x.shape[0]))

        for i, s in enumerate(signal):
            # This is largely equivalent to equation 2.  The minimum
            # indicates that a trap cannot emit more charge than is
            # available, nor can it emit more charge than it can hold.
            # This scales the exponential release of charge from the
            # trap.  The next term defines the contribution from the
            # global CTI at each pixel transfer, and the final term
            # includes the contribution from local CTI effects.
            res[i, :] = (np.minimum(v['trapsize'], s*v['scaling'])
                         * (np.exp(1/v['emissiontime']) - 1.0)
                         * np.exp(-x/v['emissiontime'])
                         + s*num_transfers*v['cti']**x
                         + v['driftscale']*s*np.exp(-x/float(v['decaytime'])))

        return res


class SimulatedModel(OverscanModel):
    """Simulated overscan model."""

    @staticmethod
    def model_results(params, signal, num_transfers, amp, start=1, stop=10, trap_type=None):
        """Generate a realization of the overscan model, using the specified
        fit parameters and input signal.

        Parameters
        ----------
        params : `lmfit.Parameters`
            Object containing the model parameters.
        signal : `np.ndarray`, (nMeasurements)
            Array of image means.
        num_transfers : `int`
            Number of serial transfers that the charge undergoes.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Amplifier to use for geometry information.
        start : `int`, optional
            First overscan column to fit. This number includes the
            last imaging column, and needs to be adjusted by one when
            using the overscan bounding box.
        stop : `int`, optional
            Last overscan column to fit. This number includes the
            last imaging column, and needs to be adjusted by one when
            using the overscan bounding box.
        trap_type : `str`, optional
            Type of trap model to use.

        Returns
        -------
        results : `np.ndarray`, (nMeasurements, nCols)
            Model results.
        """
        v = params.valuesdict()

        # Adjust column numbering to match DM overscan bbox.
        start += 1
        stop += 1

        # Electronics effect optimization
        output_amplifier = FloatingOutputAmplifier(1.0, v['driftscale'], v['decaytime'])

        # CTI optimization
        v['cti'] = 10**v['ctiexp']

        # Trap type for optimization
        if trap_type is None:
            trap = None
        elif trap_type == 'linear':
            trap = SerialTrap(v['trapsize'], v['emissiontime'], 1, 'linear',
                              [v['scaling']])
        elif trap_type == 'logistic':
            trap = SerialTrap(v['trapsize'], v['emissiontime'], 1, 'logistic',
                              [v['f0'], v['k']])
        else:
            raise ValueError('Trap type must be linear or logistic or None')

        # Simulate ramp readout
        imarr = np.zeros((signal.shape[0], amp.getRawDataBBox().getWidth()))
        ramp = SegmentSimulator(imarr, amp.getRawSerialPrescanBBox().getWidth(), output_amplifier,
                                cti=v['cti'], traps=trap)
        ramp.ramp_exp(signal)
        model_results = ramp.readout(serial_overscan_width=amp.getRawSerialOverscanBBox().getWidth(),
                                     parallel_overscan_width=0)

        ncols = amp.getRawSerialPrescanBBox().getWidth() + amp.getRawDataBBox().getWidth()

        return model_results[:, ncols+start-1:ncols+stop]


class SegmentSimulator:
    """Controls the creation of simulated segment images.

    Parameters
    ----------
    imarr : `np.ndarray` (nx, ny)
        Image data array.
    prescan_width : `int`
        Number of serial prescan columns.
    output_amplifier : `lsst.cp.pipe.FloatingOutputAmplifier`
        An object holding some deferred charge parameters.
    cti : `float`
        Global CTI value.
    traps : `list` [`lsst.ip.isr.SerialTrap`]
        Serial traps to simulate.
    """

    def __init__(self, imarr, prescan_width, output_amplifier, cti=0.0, traps=None):
        # Image array geometry
        self.prescan_width = prescan_width
        self.ny, self.nx = imarr.shape

        self.segarr = np.zeros((self.ny, self.nx+prescan_width))
        self.segarr[:, prescan_width:] = imarr

        # Serial readout information
        self.output_amplifier = output_amplifier
        if isinstance(cti, np.ndarray):
            raise ValueError("cti must be single value, not an array.")
        self.cti = cti

        self.serial_traps = None
        self.do_trapping = False
        if traps is not None:
            if not isinstance(traps, list):
                traps = [traps]
            for trap in traps:
                self.add_trap(trap)

    def add_trap(self, serial_trap):
        """Add a trap to the serial register.

        Parameters
        ----------
        serial_trap : `lsst.ip.isr.SerialTrap`
            The trap to add.
        """
        try:
            self.serial_traps.append(serial_trap)
        except AttributeError:
            self.serial_traps = [serial_trap]
            self.do_trapping = True

    def ramp_exp(self, signal_list):
        """Simulate an image with varying flux illumination per row.

        This method simulates a segment image where the signal level
        increases along the horizontal direction, according to the
        provided list of signal levels.

        Parameters
        ----------
        signal_list : `list` [`float`]
            List of signal levels.

        Raises
        ------
        ValueError
            Raised if the length of the signal list does not equal the
            number of rows.
        """
        if len(signal_list) != self.ny:
            raise ValueError("Signal list does not match row count.")

        ramp = np.tile(signal_list, (self.nx, 1)).T
        self.segarr[:, self.prescan_width:] += ramp

    def readout(self, serial_overscan_width=10, parallel_overscan_width=0):
        """Simulate serial readout of the segment image.

        This method performs the serial readout of a segment image
        given the appropriate SerialRegister object and the properties
        of the ReadoutAmplifier.  Additional arguments can be provided
        to account for the number of desired overscan transfers. The
        result is a simulated final segment image, in ADU.

        Parameters
        ----------
        serial_overscan_width : `int`, optional
            Number of serial overscan columns.
        parallel_overscan_width : `int`, optional
            Number of parallel overscan rows.

        Returns
        -------
        result : `np.ndarray` (nx, ny)
            Simulated image, including serial prescan, serial
            overscan, and parallel overscan regions. Result in electrons.
        """
        # Create output array
        iy = int(self.ny + parallel_overscan_width)
        ix = int(self.nx + self.prescan_width + serial_overscan_width)

        image = np.random.default_rng().normal(
            loc=self.output_amplifier.global_offset,
            scale=self.output_amplifier.noise,
            size=(iy, ix),
        )

        free_charge = copy.deepcopy(self.segarr)

        # Set flow control parameters
        do_trapping = self.do_trapping
        cti = self.cti

        offset = np.zeros(self.ny)
        cte = 1 - cti
        if do_trapping:
            for trap in self.serial_traps:
                trap.initialize(self.ny, self.nx, self.prescan_width)

        for i in range(ix):
            # Trap capture
            if do_trapping:
                for trap in self.serial_traps:
                    captured_charge = trap.trap_charge(free_charge)
                    free_charge -= captured_charge

            # Pixel-to-pixel proportional loss
            transferred_charge = free_charge*cte
            deferred_charge = free_charge*cti

            # Pixel transfer and readout
            offset = self.output_amplifier.local_offset(offset,
                                                        transferred_charge[:, 0])
            image[:iy-parallel_overscan_width, i] += transferred_charge[:, 0] + offset

            free_charge = np.pad(transferred_charge, ((0, 0), (0, 1)),
                                 mode='constant')[:, 1:] + deferred_charge

            # Trap emission
            if do_trapping:
                for trap in self.serial_traps:
                    released_charge = trap.release_charge()
                    free_charge += released_charge

        return image


class FloatingOutputAmplifier:
    """Object representing the readout amplifier of a single channel.

    Parameters
    ----------
    gain : `float`
        Gain of the amplifier. Currently not used.
    scale : `float`
        Drift scale for the amplifier.
    decay_time : `float`
        Decay time for the bias drift.
    noise : `float`, optional
        Amplifier read noise.
    offset : `float`, optional
        Global CTI offset.
    """

    def __init__(self, gain, scale, decay_time, noise=0.0, offset=0.0):

        self.gain = gain
        self.noise = noise
        self.global_offset = offset

        self.update_parameters(scale, decay_time)

    def local_offset(self, old, signal):
        """Calculate local offset hysteresis.

        Parameters
        ----------
        old : `np.ndarray`, (,)
            Previous iteration.
        signal : `np.ndarray`, (,)
            Current column measurements.
        Returns
        -------
        offset : `np.ndarray`
            Local offset.
        """
        new = self.scale*signal

        return np.maximum(new, old*np.exp(-1/self.decay_time))

    def update_parameters(self, scale, decay_time):
        """Update parameter values, if within acceptable values.

        Parameters
        ----------
        scale : `float`
            Drift scale for the amplifier.
        decay_time : `float`
            Decay time for the bias drift.

        Raises
        ------
        ValueError
            Raised if the input parameters are out of range.
        """
        if scale < 0.0:
            raise ValueError("Scale must be greater than or equal to 0.")
        if np.isnan(scale):
            raise ValueError("Scale must be real-valued number, not NaN.")
        self.scale = scale
        if decay_time <= 0.0:
            raise ValueError("Decay time must be greater than 0.")
        if np.isnan(decay_time):
            raise ValueError("Decay time must be real-valued number, not NaN.")
        self.decay_time = decay_time


class DeferredChargeCalib(IsrCalib):
    r"""Calibration containing deferred charge/CTI parameters.

    This includes, parameters from Snyder+2021 and exstimates of
    the serial and parallel CTI using the extended pixel edge
    response (EPER) method (also defined in Snyder+2021).

    Parameters
    ----------
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
    signals : `dict` [`str`, `np.ndarray`]
        A dictionary, keyed by amplifier name, of the mean signal
        level for each input measurement.
    inputGain : `dict` [`str`, `float`]
        A dictionary, keyed by amplifier name of the input gain used
        to calculate the overscan statistics and produce this calib.
    serialEper : `dict` [`str`, `np.ndarray`, `float`]
        A dictionary, keyed by amplifier name, of the serial EPER
        estimator of serial CTI, given in a list for each input
        measurement.
    parallelEper : `dict` [`str`, `np.ndarray`, `float`]
        A dictionary, keyed by amplifier name, of the parallel
        EPER estimator of parallel CTI, given in a list for each
        input measurement.
    serialCtiTurnoff : `dict` [`str`, `float`]
        A dictionary, keyed by amplifier name, of the serial CTI
        turnoff (unit: electrons).
    parallelCtiTurnoff : `dict` [`str`, `float`]
        A dictionary, keyed by amplifier name, of the parallel CTI
        turnoff (unit: electrons).
    serialCtiTurnoffSamplingErr : `dict` [`str`, `float`]
        A dictionary, keyed by amplifier name, of the serial CTI
        turnoff sampling error (unit: electrons).
    parallelCtiTurnoffSamplingErr : `dict` [`str`, `float`]
        A dictionary, keyed by amplifier name, of the parallel CTI
        turnoff sampling error (unit: electrons).

    Also, the values contained in this calibration are all derived
    from and image and overscan in units of electron as these are
    the most natural units in which to compute deferred charge.
    However, this means the the user should supply a reliable set
    of gains when computing the CTI statistics during ISR.

    Version 1.1 deprecates the USEGAINS attribute and standardizes
        everything to electron units.
    Version 1.2 adds the ``signal``, ``serialEper``, ``parallelEper``,
        ``serialCtiTurnoff``, ``parallelCtiTurnoff``,
        ``serialCtiTurnoffSamplingErr``, ``parallelCtiTurnoffSamplingErr``
        attributes.
    Version 1.3 adds the `inputGain` attribute.
    """
    _OBSTYPE = 'CTI'
    _SCHEMA = 'Deferred Charge'
    _VERSION = 1.3

    def __init__(self, **kwargs):
        self.driftScale = {}
        self.decayTime = {}
        self.globalCti = {}
        self.serialTraps = {}
        self.signals = {}
        self.inputGain = {}
        self.serialEper = {}
        self.parallelEper = {}
        self.serialCtiTurnoff = {}
        self.parallelCtiTurnoff = {}
        self.serialCtiTurnoffSamplingErr = {}
        self.parallelCtiTurnoffSamplingErr = {}

        # Check for deprecated kwargs
        if kwargs.pop("useGains", None) is not None:
            warnings.warn("useGains is deprecated, and will be removed "
                          "after v28.", FutureWarning)

        super().__init__(**kwargs)

        # Units are always in electron.
        self.updateMetadata(UNITS='electron')

        self.requiredAttributes.update(['driftScale', 'decayTime', 'globalCti', 'serialTraps',
                                        'inputGain', 'signals', 'serialEper', 'parallelEper',
                                        'serialCtiTurnoff', 'parallelCtiTurnoff',
                                        'serialCtiTurnoffSamplingErr',
                                        'parallelCtiTurnoffSamplingErr'])

    def fromDetector(self, detector):
        """Read metadata parameters from a detector.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.detector`
            Input detector with parameters to use.

        Returns
        -------
        calib : `lsst.ip.isr.Linearizer`
            The calibration constructed from the detector.
        """

        pass

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
        RuntimeError
            Raised if the supplied dictionary is for a different
            calibration.
        """
        calib = cls()

        if calib._OBSTYPE != dictionary['metadata']['OBSTYPE']:
            raise RuntimeError(f"Incorrect CTI supplied. Expected {calib._OBSTYPE}, "
                               f"found {dictionary['metadata']['OBSTYPE']}")

        calib.setMetadata(dictionary['metadata'])

        calib.inputGain = dictionary['inputGain']
        calib.driftScale = dictionary['driftScale']
        calib.decayTime = dictionary['decayTime']
        calib.globalCti = dictionary['globalCti']
        calib.serialCtiTurnoff = dictionary['serialCtiTurnoff']
        calib.parallelCtiTurnoff = dictionary['parallelCtiTurnoff']
        calib.serialCtiTurnoffSamplingErr = dictionary['serialCtiTurnoffSamplingErr']
        calib.parallelCtiTurnoffSamplingErr = dictionary['parallelCtiTurnoffSamplingErr']

        allAmpNames = dictionary['driftScale'].keys()

        # Some amps might not have a serial trap solution, so
        # dictionary['serialTraps'].keys() might not be equal
        # to dictionary['driftScale'].keys()
        for ampName in dictionary['serialTraps']:
            ampTraps = dictionary['serialTraps'][ampName]
            calib.serialTraps[ampName] = SerialTrap(ampTraps['size'], ampTraps['emissionTime'],
                                                    ampTraps['pixel'], ampTraps['trap_type'],
                                                    ampTraps['coeffs'])

        for ampName in allAmpNames:
            calib.signals[ampName] = np.array(dictionary['signals'][ampName], dtype=np.float64)
            calib.serialEper[ampName] = np.array(dictionary['serialEper'][ampName], dtype=np.float64)
            calib.parallelEper[ampName] = np.array(dictionary['parallelEper'][ampName], dtype=np.float64)

        calib.updateMetadata()
        return calib

    def toDict(self):
        """Return a dictionary containing the calibration properties.
        The dictionary should be able to be round-tripped through
        ``fromDict``.

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
        outDict['signals'] = self.signals
        outDict['inputGain'] = self.inputGain
        outDict['serialEper'] = self.serialEper
        outDict['parallelEper'] = self.parallelEper
        outDict['serialCtiTurnoff'] = self.serialCtiTurnoff
        outDict['parallelCtiTurnoff'] = self.parallelCtiTurnoff
        outDict['serialCtiTurnoffSamplingErr'] = self.serialCtiTurnoffSamplingErr
        outDict['parallelCtiTurnoffSamplingErr'] = self.parallelCtiTurnoffSamplingErr

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

        This method uses the ``fromDict`` method to create the
        calibration, after constructing an appropriate dictionary from
        the input tables.

        Parameters
        ----------
        tableList : `list` [`lsst.afw.table.Table`]
            List of tables to use to construct the CTI
            calibration.  Two tables are expected in this list, the
            first containing the per-amplifier CTI parameters, and the
            second containing the parameters for serial traps.

        Returns
        -------
        calib : `lsst.ip.isr.DeferredChargeCalib`
            The calibration defined in the tables.

        Raises
        ------
        ValueError
            Raised if the trap type or trap coefficients are not
            defined properly.
        """
        ampTable = tableList[0]

        inDict = {}
        inDict['metadata'] = ampTable.meta
        calibVersion = inDict['metadata']['CTI_VERSION']

        amps = ampTable['AMPLIFIER']
        driftScale = ampTable['DRIFT_SCALE']
        decayTime = ampTable['DECAY_TIME']
        globalCti = ampTable['GLOBAL_CTI']

        inDict['driftScale'] = {amp: value for amp, value in zip(amps, driftScale)}
        inDict['decayTime'] = {amp: value for amp, value in zip(amps, decayTime)}
        inDict['globalCti'] = {amp: value for amp, value in zip(amps, globalCti)}

        # Version check
        if calibVersion < 1.1:
            # This version might be in the wrong units (not electron),
            # and does not contain the gain information to convert
            # into a new calibration version.
            raise RuntimeError(f"Using old version of CTI calibration (ver. {calibVersion} < 1.1), "
                               "which is no longer supported.")
        elif calibVersion < 1.2:
            inDict['signals'] = {amp: np.array([np.nan]) for amp in amps}
            inDict['serialEper'] = {amp: np.array([np.nan]) for amp in amps}
            inDict['parallelEper'] = {amp: np.array([np.nan]) for amp in amps}
            inDict['serialCtiTurnoff'] = {amp: np.nan for amp in amps}
            inDict['parallelCtiTurnoff'] = {amp: np.nan for amp in amps}
            inDict['serialCtiTurnoffSamplingErr'] = {amp: np.nan for amp in amps}
            inDict['parallelCtiTurnoffSamplingErr'] = {amp: np.nan for amp in amps}
        else:
            signals = ampTable['SIGNALS']
            serialEper = ampTable['SERIAL_EPER']
            parallelEper = ampTable['PARALLEL_EPER']
            serialCtiTurnoff = ampTable['SERIAL_CTI_TURNOFF']
            parallelCtiTurnoff = ampTable['PARALLEL_CTI_TURNOFF']
            serialCtiTurnoffSamplingErr = ampTable['SERIAL_CTI_TURNOFF_SAMPLING_ERR']
            parallelCtiTurnoffSamplingErr = ampTable['PARALLEL_CTI_TURNOFF_SAMPLING_ERR']
            inDict['signals'] = {amp: value for amp, value in zip(amps, signals)}
            inDict['serialEper'] = {amp: value for amp, value in zip(amps, serialEper)}
            inDict['parallelEper'] = {amp: value for amp, value in zip(amps, parallelEper)}
            inDict['serialCtiTurnoff'] = {amp: value for amp, value in zip(amps, serialCtiTurnoff)}
            inDict['parallelCtiTurnoff'] = {amp: value for amp, value in zip(amps, parallelCtiTurnoff)}
            inDict['serialCtiTurnoffSamplingErr'] = {
                amp: value for amp, value in zip(amps, serialCtiTurnoffSamplingErr)
            }
            inDict['parallelCtiTurnoffSamplingErr'] = {
                amp: value for amp, value in zip(amps, parallelCtiTurnoffSamplingErr)
            }
        if calibVersion < 1.3:
            inDict['inputGain'] = {amp: np.nan for amp in amps}
        else:
            inputGain = ampTable['INPUT_GAIN']
            inDict['inputGain'] = {amp: value for amp, value in zip(amps, inputGain)}

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

            # Unpad any trailing NaN values: find the continuous array
            # of NaNs at the end of the coefficients, and remove them.
            inCoeffs = coeffs[index]
            breakIndex = 1
            nanValues = np.where(np.isnan(inCoeffs))[0]
            if nanValues is not None:
                coeffLength = len(inCoeffs)
                while breakIndex < coeffLength:
                    if coeffLength - breakIndex in nanValues:
                        breakIndex += 1
                    else:
                        break
            breakIndex -= 1  # Remove the fixed offset.
            if breakIndex != 0:
                outCoeffs = inCoeffs[0: coeffLength - breakIndex]
            else:
                outCoeffs = inCoeffs
            ampTrap['coeffs'] = outCoeffs.tolist()

            if ampTrap['trap_type'] == 'linear':
                if len(ampTrap['coeffs']) < 1:
                    raise ValueError("CTI Amplifier %s coefficients for trap has illegal length %d.",
                                     amp, len(ampTrap['coeffs']))
            elif ampTrap['trap_type'] == 'logistic':
                if len(ampTrap['coeffs']) < 2:
                    raise ValueError("CTI Amplifier %s coefficients for trap has illegal length %d.",
                                     amp, len(ampTrap['coeffs']))
            elif ampTrap['trap_type'] == 'spline':
                if len(ampTrap['coeffs']) % 2 != 0:
                    raise ValueError("CTI Amplifier %s coefficients for trap has illegal length %d.",
                                     amp, len(ampTrap['coeffs']))
            else:
                raise ValueError('Unknown trap type: %s', ampTrap['trap_type'])

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
        signals = []
        inputGain = []
        serialEper = []
        parallelEper = []
        serialCtiTurnoff = []
        parallelCtiTurnoff = []
        serialCtiTurnoffSamplingErr = []
        parallelCtiTurnoffSamplingErr = []

        for amp in self.driftScale.keys():
            ampList.append(amp)
            driftScale.append(self.driftScale[amp])
            decayTime.append(self.decayTime[amp])
            globalCti.append(self.globalCti[amp])
            signals.append(self.signals[amp])
            inputGain.append(self.inputGain[amp])
            serialEper.append(self.serialEper[amp])
            parallelEper.append(self.parallelEper[amp])
            serialCtiTurnoff.append(self.serialCtiTurnoff[amp])
            parallelCtiTurnoff.append(self.parallelCtiTurnoff[amp])
            serialCtiTurnoffSamplingErr.append(
                self.serialCtiTurnoffSamplingErr[amp]
            )
            parallelCtiTurnoffSamplingErr.append(
                self.parallelCtiTurnoffSamplingErr[amp]
            )

        ampTable = Table({
            'AMPLIFIER': ampList,
            'DRIFT_SCALE': driftScale,
            'DECAY_TIME': decayTime,
            'GLOBAL_CTI': globalCti,
            'SIGNALS': signals,
            'INPUT_GAIN': inputGain,
            'SERIAL_EPER': serialEper,
            'PARALLEL_EPER': parallelEper,
            'SERIAL_CTI_TURNOFF': serialCtiTurnoff,
            'PARALLEL_CTI_TURNOFF': parallelCtiTurnoff,
            'SERIAL_CTI_TURNOFF_SAMPLING_ERR': serialCtiTurnoffSamplingErr,
            'PARALLEL_CTI_TURNOFF_SAMPLING_ERR': parallelCtiTurnoffSamplingErr,
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
        doc="Number of prior pixels to use for local offset correction.",
        default=15,
    )
    nPixelTrapCorrection = Field(
        dtype=int,
        doc="Number of prior pixels to use for trap correction.",
        default=6,
    )
    zeroUnusedPixels = Field(
        dtype=bool,
        doc="If true, set serial prescan and parallel overscan to zero before correction.",
        default=True,
    )


class DeferredChargeTask(Task):
    """Task to correct an exposure for charge transfer inefficiency.

    This uses the methods described by Snyder et al. 2021, Journal of
    Astronimcal Telescopes, Instruments, and Systems, 7,
    048002. doi:10.1117/1.JATIS.7.4.048002 (Snyder+21).
    """
    ConfigClass = DeferredChargeConfig
    _DefaultName = 'isrDeferredCharge'

    def run(self, exposure, ctiCalib, gains=None):
        """Correct deferred charge/CTI issues.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to correct the deferred charge on.
        ctiCalib : `lsst.ip.isr.DeferredChargeCalib`
            Calibration object containing the charge transfer
            inefficiency model.
        gains : `dict` [`str`, `float`]
            A dictionary, keyed by amplifier name, of the gains to
            use.  If gains is None, the nominal gains in the amplifier
            object are used.

        Returns
        -------
        exposure : `lsst.afw.image.Exposure`
            The corrected exposure.

        Notes
        -------
        This task will read the exposure metadata and determine if
        applying gains if necessary. The correction takes place in
        units of electrons. If bootstrapping, the gains used
        will just be 1.0. and the input/output units will stay in
        adu. If the input image is in adu, the output image will be
        in units of electrons. If the input image is in electron,
        the output image will be in electron.
        """
        image = exposure.getMaskedImage().image
        detector = exposure.getDetector()

        # Get the image and overscan units.
        imageUnits = exposure.getMetadata().get("LSST ISR UNITS")

        # The deferred charge correction assumes that everything is in
        # electron units. Make it so:
        applyGains = False
        if imageUnits == "adu":
            applyGains = True

        # If we need to convert the image to electrons, check that gains
        # were supplied. CTI should not be solved or corrected without
        # supplied gains.
        if applyGains and gains is None:
            raise RuntimeError("No gains supplied for deferred charge correction.")

        with gainContext(exposure, image, apply=applyGains, gains=gains, isTrimmed=False):
            # Both the image and the overscan are in electron units.
            for amp in detector.getAmplifiers():
                ampName = amp.getName()

                ampImage = image[amp.getRawBBox()]
                if self.config.zeroUnusedPixels:
                    # We don't apply overscan subtraction, so zero these
                    # out for now.
                    ampImage[amp.getRawParallelOverscanBBox()].array[:, :] = 0.0
                    ampImage[amp.getRawSerialPrescanBBox()].array[:, :] = 0.0

                # The algorithm expects that the readout corner is in
                # the lower left corner.  Flip it to be so:
                ampData = self.flipData(ampImage.array, amp)

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
                image[amp.getRawBBox()].array[:, :] = correctedAmpData[:, :]

        return exposure

    @staticmethod
    def flipData(ampData, amp):
        """Flip data array such that readout corner is at lower-left.

        Parameters
        ----------
        ampData : `numpy.ndarray`, (nx, ny)
            Image data to flip.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Amplifier to get readout corner information.

        Returns
        -------
        ampData : `numpy.ndarray`, (nx, ny)
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

        if X_FLIP[amp.getReadoutCorner()]:
            ampData = np.fliplr(ampData)
        if Y_FLIP[amp.getReadoutCorner()]:
            ampData = np.flipud(ampData)

        return ampData

    @staticmethod
    def local_offset_inverse(inputArr, drift_scale, decay_time, num_previous_pixels=15):
        r"""Remove CTI effects from local offsets.

        This implements equation 10 of Snyder+21.  For an image with
        CTI, s'(m, n), the correction factor is equal to the maximum
        value of the set of:

        .. code-block::

            {A_L s'(m, n - j) exp(-j t / \tau_L)}_j=0^jmax

        Parameters
        ----------
        inputArr : `numpy.ndarray`, (nx, ny)
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
        outputArr : `numpy.ndarray`, (nx, ny)
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
        r"""Apply localized trapping inverse operator to pixel signals.

        This implements equation 13 of Snyder+21.  For an image with
        CTI, s'(m, n), the correction factor is equal to the maximum
        value of the set of:

        .. code-block::

            {A_L s'(m, n - j) exp(-j t / \tau_L)}_j=0^jmax

        Parameters
        ----------
        inputArr : `numpy.ndarray`, (nx, ny)
            Input image data to correct.
        trap : `lsst.ip.isr.SerialTrap`
            Serial trap describing the capture and release of charge.
        global_cti: `float`
            Mean charge transfer inefficiency, b from Snyder+21.
        num_previous_pixels : `int`, optional
            Number of previous pixels to use for correction.

        Returns
        -------
        outputArr : `numpy.ndarray`, (nx, ny)
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
