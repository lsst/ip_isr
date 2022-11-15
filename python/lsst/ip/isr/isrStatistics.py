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

    doBandingStatistics = pexConfig.Field(
        dtype=bool,
        doc="Measure image banding metric?",
        default=False,
    )
    bandingKernelSize = pexConfig.Field(
        dtype=int,
        doc="Width of box for boxcar smoothing.",
        default=3,
    )
    bandingFraction = pexConfig.Field(
        dtype=float,
        doc="Fraction of values to exclude from both high and low samples.",
        default=0.1,
    )
    bandingUseHalfDetector = pexConfig.Field(
        dtype=float,
        doc="Use only the first half set of amplifiers.",
        default=True,
    )

    doSliceStatistics = pexConfig.Field(
        dtype=bool,
        doc="Measure slice metric.",
        default=False,
    )

    stat = pexConfig.Field(
        dtype=str,
        default='MEANCLIP',
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

        sliceResults = None
        if self.config.doSliceStatistics:
            sliceResults = self.measureSliceStatistics(inputExp, overscanResults)

        return pipeBase.Struct(
            results={'CTI': ctiResults,
                     'BANDING': bandingResults,
                     'SLICE': sliceResults,
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
            ampStats['IMAGE_MEAN'] = afwMath.makeStatistics(dataRegion, self.statType,
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
                ampStats['FIRST_MEAN'] = pixelZ
                ampStats['LAST_MEAN'] = pixelA
            else:
                ampStats['FIRST_MEAN'] = pixelA
                ampStats['LAST_MEAN'] = pixelZ

            # Measure the columns of the overscan.
            if overscans[ampIter] is None:
                # The amplifier is likely entirely bad, and needs to
                # be skipped.
                self.log.warn("No overscan information available for ISR statistics for amp %s.",
                              amp.getName())
                nCols = amp.getSerialOverscanBBox().getWidth()
                ampStats['OVERSCAN_COLUMNS'] = np.full((nCols, ), np.nan)
                ampStats['OVERSCAN_VALUES'] = np.full((nCols, ), np.nan)
            else:
                overscanImage = overscans[ampIter].overscanImage
                columns = []
                values = []
                for column in range(0, overscanImage.getWidth()):
                    osMean = afwMath.makeStatistics(overscanImage.image.array[:, column],
                                                    self.statType, self.statControl).getValue()
                    columns.append(column)
                    values.append(gain * osMean)

                # We want these relative to the readout corner.  If that's
                # on the right side, we need to swap them.
                if readoutCorner in (ReadoutCorner.LR, ReadoutCorner.UR):
                    ampStats['OVERSCAN_COLUMNS'] = list(reversed(columns))
                    ampStats['OVERSCAN_VALUES'] = list(reversed(values))
                else:
                    ampStats['OVERSCAN_COLUMNS'] = columns
                    ampStats['OVERSCAN_VALUES'] = values

            outputStats[amp.getName()] = ampStats

        return outputStats

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
        kernel = np.full(self.config.bandingKernelSize, 1.0 / self.config.bandingKernelSize)

        outputStats['AMP_BANDING'] = []
        for amp, overscanData in zip(detector.getAmplifiers(), overscans):
            overscanFit = np.array(overscanData.overscanFit)
            overscanArray = overscanData.overscanImage.image.array
            rawOverscan = np.mean(overscanArray + overscanFit, axis=1)

            smoothedOverscan = np.convolve(rawOverscan, kernel, mode='valid')

            low, high = np.quantile(smoothedOverscan, [self.config.bandingFraction,
                                                       1.0 - self.config.bandingFraction])
            outputStats['AMP_BANDING'].append(float(high - low))

        if self.config.bandingUseHalfDetector:
            fullLength = len(outputStats['AMP_BANDING'])
            outputStats['DET_BANDING'] = float(np.median(outputStats['AMP_BANDING'][0:fullLength//2]))
        else:
            outputStats['DET_BANDING'] = float(np.median(outputStats['AMP_BANDING']))

        return outputStats

    def measureSliceStatistics(self, inputExp, overscans):
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

        outputStats['AMP_VSLICE'] = {}
        outputStats['AMP_HSLICE'] = {}
        for amp in detector.getAmplifiers():
            ampArray = inputExp.image[amp.getBBox()].array
            horizontalSlice = np.mean(ampArray, axis=0)
            verticalSlice = np.mean(ampArray, axis=1)
            outputStats['AMP_HSLICE'][amp.getName()] = horizontalSlice.tolist()
            outputStats['AMP_VSLICE'][amp.getName()] = verticalSlice.tolist()
            # import pdb; pdb.set_trace()

        return outputStats
