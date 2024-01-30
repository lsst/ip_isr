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

import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig

from lsst.afw.cameraGeom import ReadoutCorner


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
    doProjectionInCameraCoords = pexConfig.Field(
        dtype=bool,
        doc="Flip pixels to match readout before measurement?",
        default=False,
    )
    projectionMinimum = pexConfig.Field(
        dtype=int,
        doc="Minimum coordinate to consider.",
        default=0,
    )
    projectionMaximum = pexConfig.Field(
        dtype=int,
        doc="Maximum coordinate to consider.",
        default=-1,
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
    divisaderoImpactPixels = pexConfig.Field(
        dtype=int,
        doc="Number of edge pixels to examine for divisadero tearing.",
        default=2,
    )
    divisaderoProjectionMinimum = pexConfig.Field(
        dtype=int,
        doc="Minimum row to consider when taking mean of columns.",
        default=10,
    )
    divisaderoProjectionMaximum = pexConfig.Field(
        dtype=int,
        doc="Maximum row to consider when takening mean of columns",
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

    def run(self, inputExp, overscanResults=None, **kwargs):
        """Task to run arbitrary statistics.

        The statistics should be measured by individual methods, and
        add to the dictionary in the return struct.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            The exposure to measure.
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
             Keyword arguments.  Calibrations being passed in shoule
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
        if kwargs['ptc'] is not None:
            gains = kwargs['ptc'].gain
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
            projectionResults = self.measureProjectionStatistics(inputExp)

        divisaderoResults = None
        if self.config.doDivisaderoStatistics:
            divisaderoResults = self.measureDivisaderoStatistics(inputExp, **kwargs)

        calibDistributionResults = None
        if self.config.doCopyCalibDistributionStatistics:
            calibDistributionResults = self.copyCalibDistributionStatistics(inputExp, **kwargs)

        return pipeBase.Struct(
            results={"CTI": ctiResults,
                     "BANDING": bandingResults,
                     "PROJECTION": projectionResults,
                     "CALIBDIST": calibDistributionResults,
                     "DIVISADERO": divisaderoResults,
                     },
        )

    def measureCti(self, inputExp, overscans, gains):
        """Task to measure CTI statistics.

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
        gains : `dict` [`str` `float`]
            Dictionary of per-amplifier gains, indexed by amplifier name.

        Returns
        -------
        outputStats : `dict` [`str`, [`dict` [`str`,`float]]
            Dictionary of measurements, keyed by amplifier name and
            statistics segment.
        """
        outputStats = {}

        detector = inputExp.getDetector()
        image = inputExp.image

        # Ensure we have the same number of overscans as amplifiers.
        assert len(overscans) == len(detector.getAmplifiers())

        for ampIter, amp in enumerate(detector.getAmplifiers()):
            ampStats = {}
            gain = gains[amp.getName()]
            readoutCorner = amp.getReadoutCorner()
            # Full data region.
            dataRegion = image[amp.getBBox()]
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
                nCols = amp.getSerialOverscanBBox().getWidth()
                ampStats["OVERSCAN_COLUMNS"] = np.full((nCols, ), np.nan)
                ampStats["OVERSCAN_VALUES"] = np.full((nCols, ), np.nan)
            else:
                overscanImage = overscans[ampIter].overscanImage
                columns = []
                values = []
                for column in range(0, overscanImage.getWidth()):
                    osMean = afwMath.makeStatistics(overscanImage.image.array[:, column],
                                                    self.statType, self.statControl).getValue()
                    columns.append(column)
                    if self.config.doApplyGainsForCtiStatistics:
                        values.append(gain * osMean)
                    else:
                        values.append(osMean)

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
        outputStats : `dict` [`str`, [`dict` [`str`,`float]]
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
        outputStats : `dict` [`str`, [`dict` [`str`,`float]]
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

            if self.config.doProjectionInCameraCoords:
                if amp.getRawFlipX():
                    ampArray[:] = ampArray[:, ::-1]
                if amp.getRawFlipY():
                    ampArray[:] = ampArray[::-1, :]

            segment = slice(self.config.projectionMinimum, self.config.projectionMaximum)
            horizontalProjection = np.mean(ampArray[:, segment], axis=1)
            verticalProjection = np.mean(ampArray[segment, :], axis=0)

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
        outputStats : `dict` [`str`, [`dict` [`str`,`float]]
            Dictionary of measurements, keyed by amplifier name and
            statistics segment.
        """
        outputStats = {}

        for amp in inputExp.getDetector():
            ampStats = {}

            for calibType in ("bias", "dark", "flat"):
                if kwargs.get(calibType, None) is not None:
                    metadata = kwargs[calibType].getMetadata()
                    for pct in self.config.expectedDistributionLevels:
                        key = f"LSST CALIB {calibType.upper()} {amp.getName()} DISTRIBUTION {pct}-PCT"
                        ampStats[key] = metadata.get(key, np.nan)
                outputStats[amp.getName()] = ampStats

        return outputStats

    def measureDivisaderoStatistics(self, inputExp, **kwargs):
        """Task to measure metrics from image slicing.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            Exposure to measure.
        **kwargs :
            The flat will be selected from here.

        Returns
        -------
        outputStats : `dict` [`str`, [`dict` [`str`,`float]]
            Dictionary of measurements, keyed by amplifier name and
            statistics segment.
        """
        outputStats = {}

        for amp in inputExp.getDetector():
            ampArray = inputExp.image[amp.getBBox()].array.copy()
            # slice the top or the bottom of the amp, which is the readout side
            if amp.getReadoutCorner().name in ('UL', 'UR'):
                minRow = amp.getBBox().getHeight() - self.config.divisaderoProjectionMaximum
                maxRow = amp.getBBox().getHeight() - self.config.divisaderoProjectionMinimum
            else:
                minRow = self.config.divisaderoProjectionMinimum
                maxRow = self.config.divisaderoProjectionMaximum

            segment = slice(minRow, maxRow)
            projection, _, _ = astropy.stats.sigma_clipped_stats(ampArray[segment, :], axis=0)

            ampStats = {}
            projection = projection
            projection /= np.median(projection)
            columns = np.arange(len(projection))

            segment = slice(self.config.divisaderoEdgePixels, - self.config.divisaderoEdgePixels)
            model = np.polyfit(columns[segment], projection[segment], 1)
            modelProjection = model[0] * columns + model[1]
            divisaderoProfile = projection / modelProjection

            # look for max at the edges:
            leftMax = np.nanmax(np.abs(divisaderoProfile[0:self.config.divisaderoImpactPixels] - 1.0))
            rightMax = np.nanmax(np.abs(divisaderoProfile[-self.config.divisaderoImpactPixels:] - 1.0))

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
