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

__all__ = ["AmpOffsetConfig", "AmpOffsetTask"]

import warnings

import numpy as np
from lsst.afw.math import MEANCLIP, StatisticsControl, makeStatistics
from lsst.afw.table import SourceTable
from lsst.meas.algorithms import SourceDetectionTask, SubtractBackgroundTask
from lsst.pex.config import Config, ConfigurableField, Field
from lsst.pipe.base import Task


class AmpOffsetConfig(Config):
    """Configuration parameters for AmpOffsetTask."""

    ampEdgeInset = Field(
        doc="Number of pixels the amp edge strip is inset from the amp edge. A thin strip of pixels running "
        "parallel to the edge of the amp is used to characterize the average flux level at the amp edge.",
        dtype=int,
        default=5,
    )
    ampEdgeWidth = Field(
        doc="Pixel width of the amp edge strip, starting at ampEdgeInset and extending inwards.",
        dtype=int,
        default=64,
    )
    ampEdgeMinFrac = Field(
        doc="Minimum allowed fraction of viable pixel rows along an amp edge. No amp offset estimate will be "
        "generated for amp edges that do not have at least this fraction of unmasked pixel rows.",
        dtype=float,
        default=0.5,
    )
    ampEdgeMaxOffset = Field(
        doc="Maximum allowed amp offset ADU value. If a measured amp offset value is larger than this, the "
        "result will be discarded and therefore not used to determine amp pedestal corrections.",
        dtype=float,
        default=5.0,
    )
    ampEdgeWindow = Field(
        doc="Pixel size of the sliding window used to generate rolling average amp offset values.",
        dtype=int,
        default=512,
    )
    doBackground = Field(
        doc="Estimate and subtract background prior to amp offset estimation?",
        dtype=bool,
        default=True,
    )
    background = ConfigurableField(
        doc="An initial background estimation step run prior to amp offset calculation.",
        target=SubtractBackgroundTask,
    )
    doDetection = Field(
        doc="Detect sources and update cloned exposure prior to amp offset estimation?",
        dtype=bool,
        default=True,
    )
    detection = ConfigurableField(
        doc="Source detection to add temporary detection footprints prior to amp offset calculation.",
        target=SourceDetectionTask,
    )


class AmpOffsetTask(Task):
    """Calculate and apply amp offset corrections to an exposure."""

    ConfigClass = AmpOffsetConfig
    _DefaultName = "isrAmpOffset"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Always load background subtask, even if doBackground=False;
        # this allows for default plane bit masks to be defined.
        self.makeSubtask("background")
        if self.config.doDetection:
            self.makeSubtask("detection")

    def run(self, exposure):
        """Calculate amp offset values, determine corrective pedestals for each
        amp, and update the input exposure in-place.

        Parameters
        ----------
        exposure: `lsst.afw.image.Exposure`
            Exposure to be corrected for amp offsets.
        """

        # Generate an exposure clone to work on and establish the bit mask.
        exp = exposure.clone()
        bitMask = exp.mask.getPlaneBitMask(self.background.config.ignoredPixelMask)
        self.log.info(
            "Ignored mask planes for amp offset estimation: [%s].",
            ", ".join(self.background.config.ignoredPixelMask),
        )

        # Fit and subtract background.
        if self.config.doBackground:
            maskedImage = exp.getMaskedImage()
            bg = self.background.fitBackground(maskedImage)
            bgImage = bg.getImageF(self.background.config.algorithm, self.background.config.undersampleStyle)
            maskedImage -= bgImage

        # Detect sources and update cloned exposure mask planes in-place.
        if self.config.doDetection:
            schema = SourceTable.makeMinimalSchema()
            table = SourceTable.make(schema)
            # Detection sigma, used for smoothing and to grow detections, is
            # normally measured from the PSF of the exposure. As the PSF hasn't
            # been measured at this stage of processing, sigma is instead
            # set to an approximate value here (which should be sufficient).
            _ = self.detection.run(table=table, exposure=exp, sigma=2)

        # Safety check: do any pixels remain for amp offset estimation?
        if (exp.mask.array & bitMask).all():
            self.log.warning("All pixels masked: cannot calculate any amp offset corrections.")
        else:
            # Set up amp offset inputs.
            im = exp.image
            im.array[(exp.mask.array & bitMask) > 0] = np.nan
            amps = exp.getDetector().getAmplifiers()
            ampEdgeOuter = self.config.ampEdgeInset + self.config.ampEdgeWidth
            sctrl = StatisticsControl()

            # Determine amplifier interface geometry.
            # RAISE WARNING IF ALL AMPS NOT SAME SIZE
            ampInterfaces = self.getAmpInterfaces(amps)
            print(ampInterfaces)

            # loop over each amp edge boundary to extract amp offset values
            ampOffsets = []
            for ii in range(1, len(amps)):
                ampA = im[amps[ii - 1].getBBox()].array
                ampB = im[amps[ii].getBBox()].array
                stripA = ampA[:, -ampEdgeOuter: -self.config.ampEdgeInset]
                stripB = ampB[:, self.config.ampEdgeInset: ampEdgeOuter]

                # catch warnings to prevent all-NaN slice RuntimeWarning
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", r"All-NaN (slice|axis) encountered")
                    edgeA = np.nanmedian(stripA, axis=1)
                    edgeB = np.nanmedian(stripB, axis=1)
                edgeDiff = edgeB - edgeA
                # compute rolling averages
                edgeDiffSum = np.convolve(np.nan_to_num(edgeDiff), np.ones(self.config.ampEdgeWindow), "same")
                edgeDiffNum = np.convolve(~np.isnan(edgeDiff), np.ones(self.config.ampEdgeWindow), "same")
                edgeDiffAvg = edgeDiffSum / np.clip(edgeDiffNum, 1, None)
                edgeDiffAvg[np.isnan(edgeDiff)] = np.nan
                # take clipped mean of rolling average data as amp offset value
                ampOffset = makeStatistics(edgeDiffAvg, MEANCLIP, sctrl).getValue()
                # perform a couple of do-no-harm safety checks:
                # a) the fraction of unmasked pixel rows is > ampEdgeMinFrac,
                # b) the absolute offset ADU value is < ampEdgeMaxOffset
                ampEdgeGoodFrac = 1 - (np.sum(np.isnan(edgeDiffAvg)) / len(edgeDiffAvg))
                minFracFail = ampEdgeGoodFrac < self.config.ampEdgeMinFrac
                maxOffsetFail = np.abs(ampOffset) > self.config.ampEdgeMaxOffset
                if minFracFail or maxOffsetFail:
                    ampOffset = 0
                ampOffsets.append(ampOffset)
                self.log.debug(
                    f"amp edge {ii}{ii+1} : "
                    f"viable edge frac = {ampEdgeGoodFrac}, "
                    f"edge offset = {ampOffset:.3f}"
                )

            # solve for pedestal values and update original exposure in-place
            A = np.array(
                [[-1.0, 1.0, 0.0, 0.0], [1.0, -2.0, 1.0, 0.0], [0.0, 1.0, -2.0, 1.0], [0.0, 0.0, 1.0, -1.0]]
            )
            B = np.array(
                [ampOffsets[0], ampOffsets[1] - ampOffsets[0], ampOffsets[2] - ampOffsets[1], -ampOffsets[2]]
            )
            # if least-squares minimization fails, convert NaNs to zeroes,
            # ensuring that no values are erroneously added/subtracted
            pedestals = np.nan_to_num(np.linalg.lstsq(A, B, rcond=None)[0])
            metadata = exposure.getMetadata()
            for ii, (amp, pedestal) in enumerate(zip(amps, pedestals)):
                ampIm = exposure.image[amp.getBBox()].array
                ampIm -= pedestal
                metadata.set(
                    f"PEDESTAL{ii + 1}", float(pedestal), f"Pedestal level subtracted from amp {ii + 1}"
                )
            self.log.info(f"amp pedestal values: {', '.join([f'{x:.2f}' for x in pedestals])}")

        # raise NotImplementedError("Amp offset task should be retargeted
        # by a camera specific version.")

    def getAmpInterfaces(self, amps):
        """Determine amp geometry and amp interfaces from a list of amplifiers.

        Parse an input list of amplifiers to determine the layout of amps
        within a detector, and identify all amp interfaces (i.e., the
        horizontal and vertical junctions between amps).

        Returns a matrix with a shape corresponding to the geometry of the amps
        in the detector, with the value of each element corresponding to a
        binary bit value indicating the type of amp interfaces relevant to that
        amp. Possible binary bit values are: 1 (top-edge amp interface), 2
        (right-edge amp interface), 4 (bottom-edge amp interface), and 8
        (left-edge amp interface).

        Parameters
        ----------
        amps: `list` [`lsst.afw.cameraGeom.Amplifier`]
            List of amplifier objects.

        Returns
        -------
        ampInterfaces: `numpy.ndarray`
            Matrix with amp geometry and binary bit interface information.
        """
        xCenters = [amp.getBBox().getCenterX() for amp in amps]
        yCenters = [amp.getBBox().getCenterY() for amp in amps]

        xIndices = np.ceil(xCenters / (2 * np.min(xCenters))).astype(int) - 1
        yIndices = np.ceil(yCenters / (2 * np.min(yCenters))).astype(int) - 1
        ampInterfaces = np.zeros((len(np.unique(yIndices)), len(np.unique(xIndices))), dtype=int)

        for xIndex, yIndex, amp in zip(xIndices, yIndices, amps):
            if yIndex < np.max(yIndices):
                ampInterfaces[yIndex, xIndex] += 1
            if xIndex < np.max(xIndices):
                ampInterfaces[yIndex, xIndex] += 2
            if yIndex > 0:
                ampInterfaces[yIndex, xIndex] += 4
            if xIndex > 0:
                ampInterfaces[yIndex, xIndex] += 8

        return ampInterfaces
