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

__all__ = ["IsrStatisticsTaskConfig", "IsrStatisticsTask"]

import numpy as np
import astropy.stats
from scipy.signal.windows import hamming, hann, gaussian
from scipy.signal import butter, filtfilt
from scipy.stats import linregress

import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig

from lsst.afw.cameraGeom import ReadoutCorner
from .isrFunctions import gainContext, isTrimmedExposure


class IsrStatisticsTaskConfig(pexConfig.Config):
    """Image statistics options.
    """
    doCtiStatistics = pexConfig.Field(
        dtype=bool,
        doc="Measure CTI statistics from image and overscans?",
        default=False,
    )
    doApplyGainsForCtiStatistics = pexConfig.Field(
        dtype=bool,
        doc="Apply gain to the overscan region when measuring CTI statistics?",
        default=True,
        deprecated="This field is no longer used. Will be removed after v28.",
    )

    doBandingStatistics = pexConfig.Field(
        dtype=bool,
        doc="Measure image banding metric?",
        default=False,
    )
    bandingKernelSize = pexConfig.Field(
        dtype=int,
        doc="Width of box for boxcar smoothing for banding metric.",
        default=3,
        check=lambda x: x == 0 or x % 2 != 0,
    )
    bandingFractionLow = pexConfig.Field(
        dtype=float,
        doc="Fraction of values to exclude from low samples.",
        default=0.1,
        check=lambda x: x >= 0.0 and x <= 1.0
    )
    bandingFractionHigh = pexConfig.Field(
        dtype=float,
        doc="Fraction of values to exclude from high samples.",
        default=0.9,
        check=lambda x: x >= 0.0 and x <= 1.0,
    )
    bandingUseHalfDetector = pexConfig.Field(
        dtype=float,
        doc="Use only the first half set of amplifiers.",
        default=True,
    )

    doProjectionStatistics = pexConfig.Field(
        dtype=bool,
        doc="Measure projection metric?",
        default=False,
    )
    projectionKernelSize = pexConfig.Field(
        dtype=int,
        doc="Width of box for boxcar smoothing of projections.",
        default=0,
        check=lambda x: x == 0 or x % 2 != 0,
    )
    doProjectionFft = pexConfig.Field(
        dtype=bool,
        doc="Generate FFTs from the image projections?",
        default=False,
    )
    projectionFftWindow = pexConfig.ChoiceField(
        dtype=str,
        doc="Type of windowing to use prior to calculating FFT.",
        default="HAMMING",
        allowed={
            "HAMMING": "Hamming window.",
            "HANN": "Hann window.",
            "GAUSSIAN": "Gaussian window.",
            "NONE": "No window."
        }
    )

    doDivisaderoStatistics = pexConfig.Field(
        dtype=bool,
        doc="Measure divisadero tearing statistics?",
        default=False,
    )
    divisaderoEdgePixels = pexConfig.Field(
        dtype=int,
        doc="Number of edge pixels excluded from divisadero linear fit.",
        default=25,
    )
    divisaderoNumImpactPixels = pexConfig.Field(
        dtype=int,
        doc="Number of edge pixels to examine for divisadero tearing.",
        default=2,
    )
    divisaderoProjectionMinimum = pexConfig.Field(
        dtype=int,
        doc="Minimum row to consider when taking robust mean of columns.",
        default=10,
    )
    divisaderoProjectionMaximum = pexConfig.Field(
        dtype=int,
        doc="Maximum row to consider when taking robust mean of columns",
        default=210,
    )

    doCopyCalibDistributionStatistics = pexConfig.Field(
        dtype=bool,
        doc="Copy calibration distribution statistics to output?",
        default=False,
    )
    expectedDistributionLevels = pexConfig.ListField(
        dtype=float,
        doc="Percentile levels expected in the calibration header.",
        default=[0, 5, 16, 50, 84, 95, 100],
    )

    doBiasShiftStatistics = pexConfig.Field(
        dtype=bool,
        doc="Measure number of image shifts in overscan?",
        default=False,
    )
    biasShiftFilterOrder = pexConfig.Field(
        dtype=int,
        doc="Filter order for Butterworth highpass filter.",
        default=5,
    )
    biasShiftCutoff = pexConfig.Field(
        dtype=float,
        doc="Cutoff frequency for highpass filter.",
        default=1.0/15.0,
    )
    biasShiftWindow = pexConfig.Field(
        dtype=int,
        doc="Filter window size in pixels for highpass filter.",
        default=30,
    )
    biasShiftThreshold = pexConfig.Field(
        dtype=float,
        doc="S/N threshold for bias shift detection.",
        default=3.0,
    )
    biasShiftRowSkip = pexConfig.Field(
        dtype=int,
        doc="Number of rows to skip for the bias shift detection.",
        default=30,
    )
    biasShiftColumnSkip = pexConfig.Field(
        dtype=int,
        doc="Number of columns to skip when averaging the overscan region.",
        default=3,
    )

    doAmplifierCorrelationStatistics = pexConfig.Field(
        dtype=bool,
        doc="Measure amplifier correlations?",
        default=False,
    )

    stat = pexConfig.Field(
        dtype=str,
        default="MEANCLIP",
        doc="Statistic name to use to measure regions.",
    )
    nSigmaClip = pexConfig.Field(
        dtype=float,
        default=3.0,
        doc="Clipping threshold for background",
    )
    nIter = pexConfig.Field(
        dtype=int,
        default=3,
        doc="Clipping iterations for background",
    )
    badMask = pexConfig.ListField(
        dtype=str,
        default=["BAD", "INTRP", "SAT"],
        doc="Mask planes to ignore when identifying source pixels."
    )


class IsrStatisticsTask(pipeBase.Task):
    """Task to measure arbitrary statistics on ISR processed exposures.

    The goal is to wrap a number of optional measurements that are
    useful for calibration production and detector stability.
    """
    ConfigClass = IsrStatisticsTaskConfig
    _DefaultName = "isrStatistics"

    def __init__(self, statControl=None, **kwargs):
        super().__init__(**kwargs)
        self.statControl = afwMath.StatisticsControl(self.config.nSigmaClip, self.config.nIter,
                                                     afwImage.Mask.getPlaneBitMask(self.config.badMask))
        self.statType = afwMath.stringToStatisticsProperty(self.config.stat)

    def run(self, inputExp, ptc=None, overscanResults=None, **kwargs):
        """Task to run arbitrary statistics.

        The statistics should be measured by individual methods, and
        add to the dictionary in the return struct.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            The exposure to measure.
        ptc : `lsst.ip.isr.PtcDataset`, optional
            A PTC object containing gains to use.
        overscanResults : `list` [`lsst.pipe.base.Struct`], optional
            List of overscan results.  Expected fields are:

            ``imageFit``
                Value or fit subtracted from the amplifier image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanFit``
                Value or fit subtracted from the overscan image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanImage``
                Image of the overscan region with the overscan
                correction applied (`lsst.afw.image.Image`). This
                quantity is used to estimate the amplifier read noise
                empirically.
        **kwargs :
             Keyword arguments.  Calibrations being passed in should
             have an entry here.

        Returns
        -------
        resultStruct : `lsst.pipe.base.Struct`
            Contains the measured statistics as a dict stored in a
            field named ``results``.

        Raises
        ------
        RuntimeError
            Raised if the amplifier gains could not be found.
        """
        # Find gains.
        detector = inputExp.getDetector()
        if ptc is not None:
            gains = ptc.gain
        elif detector is not None:
            gains = {amp.getName(): amp.getGain() for amp in detector.getAmplifiers()}
        else:
            raise RuntimeError("No source of gains provided.")

        ctiResults = None
        if self.config.doCtiStatistics:
            ctiResults = self.measureCti(inputExp, overscanResults, gains)

        bandingResults = None
        if self.config.doBandingStatistics:
            bandingResults = self.measureBanding(inputExp, overscanResults)

        projectionResults = None
        if self.config.doProjectionStatistics:
            projectionResults = self.measureProjectionStatistics(inputExp, overscanResults)

        divisaderoResults = None
        if self.config.doDivisaderoStatistics:
            divisaderoResults = self.measureDivisaderoStatistics(inputExp, **kwargs)

        calibDistributionResults = None
        if self.config.doCopyCalibDistributionStatistics:
            calibDistributionResults = self.copyCalibDistributionStatistics(inputExp, **kwargs)

        biasShiftResults = None
        if self.config.doBiasShiftStatistics:
            biasShiftResults = self.measureBiasShifts(inputExp, overscanResults)

        ampCorrelationResults = None
        if self.config.doAmplifierCorrelationStatistics:
            ampCorrelationResults = self.measureAmpCorrelations(inputExp, overscanResults)

        mjd = inputExp.getMetadata().get("MJD", None)

        return pipeBase.Struct(
            results={"CTI": ctiResults,
                     "BANDING": bandingResults,
                     "PROJECTION": projectionResults,
                     "CALIBDIST": calibDistributionResults,
                     "BIASSHIFT": biasShiftResults,
                     "AMPCORR": ampCorrelationResults,
                     "MJD": mjd,
                     'DIVISADERO': divisaderoResults,
                     },
        )

    def measureCti(self, inputExp, overscans, gains):
        """Task to measure CTI statistics.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            Exposure to measure.
        overscans : `list` [`lsst.pipe.base.Struct`]
            List of overscan results (expects base units of adu).
            Expected fields are:

            ``imageFit``
                Value or fit subtracted from the amplifier image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanFit``
                Value or fit subtracted from the overscan image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanImage``
                Image of the overscan region with the overscan
                correction applied (`lsst.afw.image.Image`). This
                quantity is used to estimate the amplifier read noise
                empirically.
        gains : `dict` [`str` `float`]
            Dictionary of per-amplifier gains, indexed by amplifier name.

        Returns
        -------
        outputStats : `dict` [`str`, [`dict` [`str`, `float`]]]
            Dictionary of measurements, keyed by amplifier name and
            statistics segment. Everything in units based on electron.
        """
        outputStats = {}

        detector = inputExp.getDetector()
        image = inputExp.image

        # It only makes sense to measure CTI in electron units.
        # Make it so.
        imageUnits = inputExp.getMetadata().get("LSST ISR UNITS")
        applyGain = False
        if imageUnits == "adu":
            applyGain = True

        # Check if the image is trimmed.
        isTrimmed = isTrimmedExposure(inputExp)

        # Ensure we have the same number of overscans as amplifiers.
        assert len(overscans) == len(detector.getAmplifiers())

        with gainContext(inputExp, image, applyGain, gains, isTrimmed=isTrimmed:
            for ampIter, amp in enumerate(detector.getAmplifiers()):
                ampStats = {}
                readoutCorner = amp.getReadoutCorner()

                # Full data region.
                dataRegion = image[amp.getBBox() if isTrimmed else amp.getRawDataBBox()]

                # Get the mean of the image
                ampStats["IMAGE_MEAN"] = afwMath.makeStatistics(dataRegion, self.statType,
                                                                self.statControl).getValue()

                # First and last image columns.
                pixelA = afwMath.makeStatistics(dataRegion.array[:, 0],
                                                self.statType,
                                                self.statControl).getValue()
                pixelZ = afwMath.makeStatistics(dataRegion.array[:, -1],
                                                self.statType,
                                                self.statControl).getValue()

                # We want these relative to the readout corner.  If that's
                # on the right side, we need to swap them.
                if readoutCorner in (ReadoutCorner.LR, ReadoutCorner.UR):
                    ampStats["FIRST_MEAN"] = pixelZ
                    ampStats["LAST_MEAN"] = pixelA
                else:
                    ampStats["FIRST_MEAN"] = pixelA
                    ampStats["LAST_MEAN"] = pixelZ

                # Measure the columns of the overscan.
                if overscans[ampIter] is None:
                    # The amplifier is likely entirely bad, and needs to
                    # be skipped.
                    self.log.warning("No overscan information available for ISR statistics for amp %s.",
                                     amp.getName())
                    nCols = amp.getRawSerialOverscanBBox().getWidth()
                    ampStats["OVERSCAN_COLUMNS"] = np.full((nCols, ), np.nan)
                    ampStats["OVERSCAN_VALUES"] = np.full((nCols, ), np.nan)
                else:
                    overscanImage = overscans[ampIter].overscanImage

                    columns = []
                    values = []
                    for column in range(0, overscanImage.getWidth()):
                        # If overscan.doParallelOverscan=True, the
                        # overscanImage will contain both the serial
                        # and parallel overscan regions.
                        # Only the serial overscan correction is
                        # implemented, so we must select only the
                        # serial overscan rows for a given column.
                        nRows = amp.getRawSerialOverscanBBox().getHeight()
                        osMean = afwMath.makeStatistics(overscanImage.image.array[:nRows, column],
                                                        self.statType, self.statControl).getValue()
                        columns.append(column)
                        # The overscan input is always in adu, but it only
                        # makes sense to measure CTI in electron units.
                        values.append(osMean * gains[amp.getName()])

                    # We want these relative to the readout corner.  If that's
                    # on the right side, we need to swap them.
                    if readoutCorner in (ReadoutCorner.LR, ReadoutCorner.UR):
                        ampStats["OVERSCAN_COLUMNS"] = list(reversed(columns))
                        ampStats["OVERSCAN_VALUES"] = list(reversed(values))
                    else:
                        ampStats["OVERSCAN_COLUMNS"] = columns
                        ampStats["OVERSCAN_VALUES"] = values

                outputStats[amp.getName()] = ampStats

        return outputStats

    @staticmethod
    def makeKernel(kernelSize):
        """Make a boxcar smoothing kernel.

        Parameters
        ----------
        kernelSize : `int`
            Size of the kernel in pixels.

        Returns
        -------
        kernel : `np.array`
            Kernel for boxcar smoothing.
        """
        if kernelSize > 0:
            kernel = np.full(kernelSize, 1.0 / kernelSize)
        else:
            kernel = np.array([1.0])
        return kernel

    def measureBanding(self, inputExp, overscans):
        """Task to measure banding statistics.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            Exposure to measure.
        overscans : `list` [`lsst.pipe.base.Struct`]
            List of overscan results.  Expected fields are:

            ``imageFit``
                Value or fit subtracted from the amplifier image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanFit``
                Value or fit subtracted from the overscan image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanImage``
                Image of the overscan region with the overscan
                correction applied (`lsst.afw.image.Image`). This
                quantity is used to estimate the amplifier read noise
                empirically.

        Returns
        -------
        outputStats : `dict` [`str`, [`dict` [`str`, `float`]]]
            Dictionary of measurements, keyed by amplifier name and
            statistics segment.
        """
        outputStats = {}

        detector = inputExp.getDetector()
        kernel = self.makeKernel(self.config.bandingKernelSize)

        outputStats["AMP_BANDING"] = []
        for amp, overscanData in zip(detector.getAmplifiers(), overscans):
            overscanFit = np.array(overscanData.overscanFit)
            overscanArray = overscanData.overscanImage.image.array
            rawOverscan = np.mean(overscanArray + overscanFit, axis=1)

            smoothedOverscan = np.convolve(rawOverscan, kernel, mode="valid")

            low, high = np.quantile(smoothedOverscan, [self.config.bandingFractionLow,
                                                       self.config.bandingFractionHigh])
            outputStats["AMP_BANDING"].append(float(high - low))

        if self.config.bandingUseHalfDetector:
            fullLength = len(outputStats["AMP_BANDING"])
            outputStats["DET_BANDING"] = float(np.nanmedian(outputStats["AMP_BANDING"][0:fullLength//2]))
        else:
            outputStats["DET_BANDING"] = float(np.nanmedian(outputStats["AMP_BANDING"]))

        return outputStats

    def measureProjectionStatistics(self, inputExp, overscans):
        """Task to measure metrics from image slicing.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            Exposure to measure.
        overscans : `list` [`lsst.pipe.base.Struct`]
            List of overscan results.  Expected fields are:

            ``imageFit``
                Value or fit subtracted from the amplifier image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanFit``
                Value or fit subtracted from the overscan image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanImage``
                Image of the overscan region with the overscan
                correction applied (`lsst.afw.image.Image`). This
                quantity is used to estimate the amplifier read noise
                empirically.

        Returns
        -------
        outputStats : `dict` [`str`, [`dict` [`str`, `float`]]]
            Dictionary of measurements, keyed by amplifier name and
            statistics segment.
        """
        outputStats = {}

        detector = inputExp.getDetector()
        kernel = self.makeKernel(self.config.projectionKernelSize)

        outputStats["AMP_VPROJECTION"] = {}
        outputStats["AMP_HPROJECTION"] = {}
        convolveMode = "valid"
        if self.config.doProjectionFft:
            outputStats["AMP_VFFT_REAL"] = {}
            outputStats["AMP_VFFT_IMAG"] = {}
            outputStats["AMP_HFFT_REAL"] = {}
            outputStats["AMP_HFFT_IMAG"] = {}
            convolveMode = "same"

        for amp in detector.getAmplifiers():
            ampArray = inputExp.image[amp.getBBox()].array

            horizontalProjection = np.mean(ampArray, axis=0)
            verticalProjection = np.mean(ampArray, axis=1)

            horizontalProjection = np.convolve(horizontalProjection, kernel, mode=convolveMode)
            verticalProjection = np.convolve(verticalProjection, kernel, mode=convolveMode)

            outputStats["AMP_HPROJECTION"][amp.getName()] = horizontalProjection.tolist()
            outputStats["AMP_VPROJECTION"][amp.getName()] = verticalProjection.tolist()

            if self.config.doProjectionFft:
                horizontalWindow = np.ones_like(horizontalProjection)
                verticalWindow = np.ones_like(verticalProjection)
                if self.config.projectionFftWindow == "NONE":
                    pass
                elif self.config.projectionFftWindow == "HAMMING":
                    horizontalWindow = hamming(len(horizontalProjection))
                    verticalWindow = hamming(len(verticalProjection))
                elif self.config.projectionFftWindow == "HANN":
                    horizontalWindow = hann(len(horizontalProjection))
                    verticalWindow = hann(len(verticalProjection))
                elif self.config.projectionFftWindow == "GAUSSIAN":
                    horizontalWindow = gaussian(len(horizontalProjection))
                    verticalWindow = gaussian(len(verticalProjection))
                else:
                    raise RuntimeError(f"Invalid window function: {self.config.projectionFftWindow}")

                horizontalFFT = np.fft.rfft(np.multiply(horizontalProjection, horizontalWindow))
                verticalFFT = np.fft.rfft(np.multiply(verticalProjection, verticalWindow))

                outputStats["AMP_HFFT_REAL"][amp.getName()] = np.real(horizontalFFT).tolist()
                outputStats["AMP_HFFT_IMAG"][amp.getName()] = np.imag(horizontalFFT).tolist()
                outputStats["AMP_VFFT_REAL"][amp.getName()] = np.real(verticalFFT).tolist()
                outputStats["AMP_VFFT_IMAG"][amp.getName()] = np.imag(verticalFFT).tolist()

        return outputStats

    def copyCalibDistributionStatistics(self, inputExp, **kwargs):
        """Copy calibration statistics for this exposure.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            The exposure being processed.
        **kwargs :
            Keyword arguments with calibrations.

        Returns
        -------
        outputStats : `dict` [`str`, [`dict` [`str`, `float`]]]
            Dictionary of measurements, keyed by amplifier name and
            statistics segment.
        """
        outputStats = {}

        # Amp level elements
        for amp in inputExp.getDetector():
            ampStats = {}

            for calibType in ("bias", "dark", "flat"):
                if kwargs.get(calibType, None) is not None:
                    metadata = kwargs[calibType].getMetadata()
                    for pct in self.config.expectedDistributionLevels:
                        key = f"LSST CALIB {calibType.upper()} {amp.getName()} DISTRIBUTION {pct}-PCT"
                        ampStats[key] = metadata.get(key, np.nan)

            for calibType in ("defects"):
                if kwargs.get(calibType, None) is not None:
                    metadata = kwargs[calibType].getMetadata()
                    for key in (f"LSST CALIB {calibType.upper()} {amp.getName()} N_HOT",
                                f"LSST CALIB {calibType.upper()} {amp.getName()} N_COLD"):
                        ampStats[key] = metadata.get(key, np.nan)
            outputStats[amp.getName()] = ampStats

        # Detector level elements
        for calibType in ("defects"):
            if kwargs.get(calibType, None) is not None:
                metadata = kwargs[calibType].getMetadata()
                for key in (f"LSST CALIB {calibType.upper()} N_BAD_COLUMNS"):
                    outputStats["detector"][key] = metadata.get(key, np.nan)

        return outputStats

    def measureBiasShifts(self, inputExp, overscanResults):
        """Measure number of bias shifts from overscan data.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            Exposure to measure.
        overscans : `list` [`lsst.pipe.base.Struct`]
            List of overscan results.  Expected fields are:

            ``imageFit``
                Value or fit subtracted from the amplifier image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanFit``
                Value or fit subtracted from the overscan image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanImage``
                Image of the overscan region with the overscan
                correction applied (`lsst.afw.image.Image`). This
                quantity is used to estimate the amplifier read noise
                empirically.

        Returns
        -------
        outputStats : `dict` [`str`, [`dict` [`str`, `float`]]]
            Dictionary of measurements, keyed by amplifier name and
            statistics segment.

        Notes
        -----
        Based on eop_pipe implementation:
        https://github.com/lsst-camera-dh/eo_pipe/blob/main/python/lsst/eo/pipe/biasShiftsTask.py  # noqa: E501 W505
        """
        outputStats = {}

        detector = inputExp.getDetector()
        for amp, overscans in zip(detector, overscanResults):
            ampStats = {}
            # Add fit back to data
            rawOverscan = overscans.overscanImage.image.array + overscans.overscanFit

            # Collapse array, skipping first three columns
            rawOverscan = np.mean(rawOverscan[:, self.config.biasShiftColumnSkip:], axis=1)

            # Scan for shifts
            noise, shift_peaks = self._scan_for_shifts(rawOverscan)
            ampStats["LOCAL_NOISE"] = float(noise)
            ampStats["BIAS_SHIFTS"] = shift_peaks

            outputStats[amp.getName()] = ampStats
        return outputStats

    def _scan_for_shifts(self, overscanData):
        """Scan overscan data for shifts.

        Parameters
        ----------
        overscanData : `list` [`float`]
             Overscan data to search for shifts.

        Returns
        -------
        noise : `float`
            Noise estimated from Butterworth filtered overscan data.
        peaks : `list` [`float`, `float`, `int`, `int`]
            Shift peak information, containing the convolved peak
            value, the raw peak value, and the lower and upper bounds
            of the region checked.
        """
        numerator, denominator = butter(self.config.biasShiftFilterOrder,
                                        self.config.biasShiftCutoff,
                                        btype="high", analog=False)
        noise = np.std(filtfilt(numerator, denominator, overscanData))
        kernel = np.concatenate([np.arange(self.config.biasShiftWindow),
                                 np.arange(-self.config.biasShiftWindow + 1, 0)])
        kernel = kernel/np.sum(kernel[:self.config.biasShiftWindow])

        convolved = np.convolve(overscanData, kernel, mode="valid")
        convolved = np.pad(convolved, (self.config.biasShiftWindow - 1, self.config.biasShiftWindow))

        shift_check = np.abs(convolved)/noise
        shift_mask = shift_check > self.config.biasShiftThreshold
        shift_mask[:self.config.biasShiftRowSkip] = False

        shift_regions = np.flatnonzero(np.diff(np.r_[np.int8(0),
                                                     shift_mask.view(np.int8),
                                                     np.int8(0)])).reshape(-1, 2)
        shift_peaks = []
        for region in shift_regions:
            region_peak = np.argmax(shift_check[region[0]:region[1]]) + region[0]
            if self._satisfies_flatness(region_peak, convolved[region_peak], overscanData):
                shift_peaks.append(
                    [float(convolved[region_peak]), float(region_peak),
                     int(region[0]), int(region[1])])
        return noise, shift_peaks

    def _satisfies_flatness(self, shiftRow, shiftPeak, overscanData):
        """Determine if a region is flat.

        Parameters
        ----------
        shiftRow : `int`
            Row with possible peak.
        shiftPeak : `float`
            Value at the possible peak.
        overscanData : `list` [`float`]
            Overscan data used to fit around the possible peak.

        Returns
        -------
        isFlat : `bool`
            Indicates if the region is flat, and so the peak is valid.
        """
        prerange = np.arange(shiftRow - self.config.biasShiftWindow, shiftRow)
        postrange = np.arange(shiftRow, shiftRow + self.config.biasShiftWindow)

        preFit = linregress(prerange, overscanData[prerange])
        postFit = linregress(postrange, overscanData[postrange])

        if shiftPeak > 0:
            preTrend = (2*preFit[0]*len(prerange) < shiftPeak)
            postTrend = (2*postFit[0]*len(postrange) < shiftPeak)
        else:
            preTrend = (2*preFit[0]*len(prerange) > shiftPeak)
            postTrend = (2*postFit[0]*len(postrange) > shiftPeak)

        return (preTrend and postTrend)

    def measureAmpCorrelations(self, inputExp, overscanResults):
        """Measure correlations between amplifier segments.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            Exposure to measure.
        overscans : `list` [`lsst.pipe.base.Struct`]
            List of overscan results.  Expected fields are:

            ``imageFit``
                Value or fit subtracted from the amplifier image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanFit``
                Value or fit subtracted from the overscan image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanImage``
                Image of the overscan region with the overscan
                correction applied (`lsst.afw.image.Image`). This
                quantity is used to estimate the amplifier read noise
                empirically.

        Returns
        -------
        outputStats : `dict` [`str`, [`dict` [`str`, `float`]]]
            Dictionary of measurements, keyed by amplifier name and
            statistics segment.

        Notes
        -----
        Based on eo_pipe implementation:
        https://github.com/lsst-camera-dh/eo_pipe/blob/main/python/lsst/eo/pipe/raft_level_correlations.py  # noqa: E501 W505
        """
        detector = inputExp.getDetector()
        rows = []

        for ampId, overscan in enumerate(overscanResults):
            rawOverscan = overscan.overscanImage.image.array + overscan.overscanFit
            rawOverscan = rawOverscan.ravel()

            ampImage = inputExp[detector[ampId].getBBox()]
            ampImage = ampImage.image.array.ravel()

            for ampId2, overscan2 in enumerate(overscanResults):
                osC = 1.0
                imC = 1.0
                if ampId2 != ampId:
                    rawOverscan2 = overscan2.overscanImage.image.array + overscan2.overscanFit
                    rawOverscan2 = rawOverscan2.ravel()

                    osC = np.corrcoef(rawOverscan, rawOverscan2)[0, 1]

                    ampImage2 = inputExp[detector[ampId2].getBBox()]
                    ampImage2 = ampImage2.image.array.ravel()

                    imC = np.corrcoef(ampImage, ampImage2)[0, 1]
                rows.append(
                    {'detector': detector.getId(),
                     'detectorComp': detector.getId(),
                     'ampName': detector[ampId].getName(),
                     'ampComp': detector[ampId2].getName(),
                     'imageCorr': float(imC),
                     'overscanCorr': float(osC),
                     }
                )
        return rows

    def measureDivisaderoStatistics(self, inputExp, **kwargs):
        """Measure Max Divisadero Tearing effect per amp.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            Exposure to measure. Usually a flat.
        **kwargs :
            The flat will be selected from here.

        Returns
        -------
        outputStats : `dict` [`str`, [`dict` [`str`, `float`]]]
            Dictionary of measurements, keyed by amplifier name and
            statistics segment.
            Measurements include

            - DIVISADERO_PROFILE: Robust mean of rows between
              divisaderoProjection<Maximum|Minumum> on readout edge of ccd
              normalized by a linear fit to the same rows.
            - DIVISADERO_MAX_PAIR: Tuple of maximum of the absolute values of
              the DIVISADERO_PROFILE, for number of pixels (specified by
              divisaderoNumImpactPixels on left and right side of amp.
            - DIVISADERO_MAX: Maximum of the absolute values of the
              the DIVISADERO_PROFILE, for the divisaderoNumImpactPixels on
              boundaries of neighboring amps (including the pixels in those
              neighborboring amps).
        """
        outputStats = {}

        for amp in inputExp.getDetector():
            # Copy unneeded if we do not ever modify the array by flipping
            ampArray = inputExp.image[amp.getBBox()].array
            # slice the outer top or bottom of the amp: the readout side
            if amp.getReadoutCorner().name in ('UL', 'UR'):
                minRow = amp.getBBox().getHeight() - self.config.divisaderoProjectionMaximum
                maxRow = amp.getBBox().getHeight() - self.config.divisaderoProjectionMinimum
            else:
                minRow = self.config.divisaderoProjectionMinimum
                maxRow = self.config.divisaderoProjectionMaximum

            segment = slice(minRow, maxRow)
            projection, _, _ = astropy.stats.sigma_clipped_stats(ampArray[segment, :], axis=0)

            ampStats = {}
            projection /= np.median(projection)
            columns = np.arange(len(projection))

            segment = slice(self.config.divisaderoEdgePixels, -self.config.divisaderoEdgePixels)
            model = np.polyfit(columns[segment], projection[segment], 1)
            modelProjection = model[0] * columns + model[1]
            divisaderoProfile = projection / modelProjection

            # look for max at the edges:
            leftMax = np.nanmax(np.abs(divisaderoProfile[0:self.config.divisaderoNumImpactPixels] - 1.0))
            rightMax = np.nanmax(np.abs(divisaderoProfile[-self.config.divisaderoNumImpactPixels:] - 1.0))

            ampStats['DIVISADERO_PROFILE'] = np.array(divisaderoProfile).tolist()
            ampStats['DIVISADERO_MAX_PAIR'] = [leftMax, rightMax]
            outputStats[amp.getName()] = ampStats

        detector = inputExp.getDetector()
        xCenters = [amp.getBBox().getCenterX() for amp in detector]
        yCenters = [amp.getBBox().getCenterY() for amp in detector]
        xIndices = np.ceil(xCenters / np.min(xCenters) / 2).astype(int) - 1
        yIndices = np.ceil(yCenters / np.min(yCenters) / 2).astype(int) - 1
        ampIds = np.zeros((len(set(yIndices)), len(set(xIndices))), dtype=int)
        for ampId, xIndex, yIndex in zip(np.arange(len(detector)), xIndices, yIndices):
            ampIds[yIndex, xIndex] = ampId

        # Loop over amps again because the DIVISIDERO_MAX will be the max
        # of the profile on its boundary with its neighboring amps
        for i, amp in enumerate(detector):
            y, x = np.where(ampIds == i)
            end = ampIds.shape[1] - 1
            xInd = x[0]
            yInd = y[0]
            thisAmpsPair = outputStats[amp.getName()]['DIVISADERO_MAX_PAIR']

            if x == 0:
                # leftmost amp: take the max of your right side and
                myMax = thisAmpsPair[1]
                # your neighbor's left side
                neighborMax = outputStats[detector[ampIds[yInd, 1]].getName()]['DIVISADERO_MAX_PAIR'][0]
            elif x == end:
                # rightmost amp: take the max of your left side and
                myMax = thisAmpsPair[0]
                # your neighbor's right side
                neighborMax = outputStats[detector[ampIds[yInd, end - 1]].getName()]['DIVISADERO_MAX_PAIR'][1]
            else:
                # Middle amp: take the max of both your own sides and the
                myMax = max(thisAmpsPair)
                leftName = detector[ampIds[yInd, max(xInd - 1, 0)]].getName()
                rightName = detector[ampIds[yInd, min(xInd + 1, ampIds.shape[1] - 1)]].getName()
                # right side of the neighbor to your left
                # and left side of your neighbor to your right
                neighborMax = max(outputStats[leftName]['DIVISADERO_MAX_PAIR'][1],
                                  outputStats[rightName]['DIVISADERO_MAX_PAIR'][0])

            divisaderoMax = max([myMax, neighborMax])
            outputStats[amp.getName()]['DIVISADERO_MAX'] = divisaderoMax

        return outputStats
