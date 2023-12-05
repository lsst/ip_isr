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
from lsst.pipe.base import Struct, Task


class AmpOffsetConfig(Config):
    """Configuration parameters for AmpOffsetTask."""

    def setDefaults(self):
        self.background.algorithm = "AKIMA_SPLINE"
        self.background.useApprox = False
        self.background.ignoredPixelMask = [
            "BAD",
            "SAT",
            "INTRP",
            "CR",
            "EDGE",
            "DETECTED",
            "DETECTED_NEGATIVE",
            "SUSPECT",
            "NO_DATA",
        ]
        self.detection.reEstimateBackground = False

        # This maintains existing behavior and test values after DM-39796.
        self.detection.thresholdType = "stdev"

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
    ampEdgeWindowFrac = Field(
        doc="Fraction of the amp edge lengths utilized as the sliding window for generating rolling average "
        "amp offset values. It should be reconfigured for every instrument (HSC, LSSTCam, etc.) and should "
        "not exceed 1. If not provided, it defaults to the fraction that recovers the pixel size of the "
        "sliding window used in obs_subaru for compatibility with existing HSC data.",
        dtype=float,
        default=512 / 4176,
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
    backgroundFractionSample = Field(
        doc="The fraction of the shorter side of the amplifier used for background binning.",
        dtype=float,
        default=1.0,
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
    applyWeights = Field(
        doc="Weights the amp offset calculation by the length of the interface between amplifiers. Applying "
        "weights does not affect outcomes for amplifiers in a 2D grid with square-shaped amplifiers or in "
        "any 1D layout on a detector, regardless of whether the amplifiers are square.",
        dtype=bool,
        default=True,
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
        # Initialize all of the instance variables here.
        self.shortAmpSide = 0

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
        amps = exp.getDetector().getAmplifiers()

        # Check that all amps have the same gemotry.
        ampDims = [amp.getBBox().getDimensions() for amp in amps]
        if not all(dim == ampDims[0] for dim in ampDims):
            raise RuntimeError("All amps should have the same geometry.")
        else:
            # The zeroth amp is representative of all amps in the detector.
            self.ampDims = ampDims[0]
            # Dictionary mapping side numbers to interface lengths.
            # See `getAmpAssociations()` for details about sides.
            self.interfaceLengthLookupBySide = {i: self.ampDims[i % 2] for i in range(4)}

        # Determine amplifier geometry.
        ampWidths = {amp.getBBox().getWidth() for amp in amps}
        ampHeights = {amp.getBBox().getHeight() for amp in amps}
        if len(ampWidths) > 1 or len(ampHeights) > 1:
            raise NotImplementedError(
                "Amp offset correction is not yet implemented for detectors with differing amp sizes."
            )

        # Assuming all the amps have the same geometry.
        self.shortAmpSide = np.min(ampDims[0])

        # Check that the edge width and inset are not too large.
        if self.config.ampEdgeWidth >= self.shortAmpSide - 2 * self.config.ampEdgeInset:
            raise RuntimeError(
                f"The edge width ({self.config.ampEdgeWidth}) plus insets ({self.config.ampEdgeInset}) "
                f"exceed the amp's short side ({self.shortAmpSide}). This setup leads to incorrect results."
            )

        # Fit and subtract background.
        if self.config.doBackground:
            maskedImage = exp.getMaskedImage()
            # Assuming all the detectors are the same.
            nX = exp.getWidth() // (self.shortAmpSide * self.config.backgroundFractionSample) + 1
            nY = exp.getHeight() // (self.shortAmpSide * self.config.backgroundFractionSample) + 1
            # This ensures that the `binSize` is as large as possible,
            # preventing background subtraction from inadvertently removing the
            # amp offset signature. Here it's set to the shorter dimension of
            # the amplifier by default (`backgroundFractionSample` = 1), which
            # seems reasonable.
            bg = self.background.fitBackground(maskedImage, nx=int(nX), ny=int(nY))
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
            self.log.warning(
                "All pixels masked: cannot calculate any amp offset corrections. All pedestals are being set "
                "to zero."
            )
            pedestals = np.zeros(len(amps))
        else:
            # Set up amp offset inputs.
            im = exp.image
            im.array[(exp.mask.array & bitMask) > 0] = np.nan

            if self.config.ampEdgeWindowFrac > 1:
                raise RuntimeError(
                    f"The specified fraction (`ampEdgeWindowFrac`={self.config.ampEdgeWindowFrac}) of the "
                    "edge length exceeds 1. This leads to complications downstream, after convolution in "
                    "the `getInterfaceOffset()` method. Please modify the `ampEdgeWindowFrac` value in the "
                    "config to be 1 or less and rerun."
                )

            # Obtain association and offset matrices.
            A, sides = self.getAmpAssociations(amps)
            B = self.getAmpOffsets(im, amps, A, sides)

            # If least-squares minimization fails, convert NaNs to zeroes,
            # ensuring that no values are erroneously added/subtracted.
            pedestals = np.nan_to_num(np.linalg.lstsq(A, B, rcond=None)[0])

        metadata = exposure.getMetadata()
        for amp, pedestal in zip(amps, pedestals):
            ampIm = exposure.image[amp.getBBox()].array
            ampIm -= pedestal
            ampName = amp.getName()
            metadata.set(
                f"LSST ISR AMPOFFSET PEDESTAL {ampName}",
                float(pedestal),
                f"Pedestal level subtracted from amp {ampName}",
            )
        self.log.info(f"amp pedestal values: {', '.join([f'{x:.4f}' for x in pedestals])}")

        return Struct(pedestals=pedestals)

    def getAmpAssociations(self, amps):
        """Determine amp geometry and amp associations from a list of
        amplifiers.

        Parse an input list of amplifiers to determine the layout of amps
        within a detector, and identify all amp sides (i.e., the
        horizontal and vertical junctions between amps).

        Returns a matrix with a shape corresponding to the geometry of the amps
        in the detector.

        Parameters
        ----------
        amps : `list` [`lsst.afw.cameraGeom.Amplifier`]
            List of amplifier objects used to deduce associations.

        Returns
        -------
        ampAssociations : `numpy.ndarray`
            An N x N matrix (N = number of amplifiers) that illustrates the
            connections between amplifiers within the detector layout. Each row
            and column index corresponds to the ampIds of a specific pair of
            amplifiers, and the matrix elements indicate their associations as
            follows:
            0: No association
            -1: Association exists (direction specified in the ampSides matrix)
            n >= 1: Diagonal elements indicate the number of neighboring
                    amplifiers for the corresponding ampId==row==column number.

        ampSides : `numpy.ndarray`
            An N x N matrix (N = the number of amplifiers) representing the amp
            side information corresponding to the `ampAssociations`
            matrix. The elements are integers defined as below:
            -1: No side due to no association or the same amp (diagonals)
            0: Side on the bottom
            1: Side on the right
            2: Side on the top
            3: Side on the left
        """
        xCenters = [amp.getBBox().getCenterX() for amp in amps]
        yCenters = [amp.getBBox().getCenterY() for amp in amps]
        xIndices = np.ceil(xCenters / np.min(xCenters) / 2).astype(int) - 1
        yIndices = np.ceil(yCenters / np.min(yCenters) / 2).astype(int) - 1

        nAmps = len(amps)
        ampIds = np.zeros((len(set(yIndices)), len(set(xIndices))), dtype=int)

        for ampId, xIndex, yIndex in zip(np.arange(nAmps), xIndices, yIndices):
            ampIds[yIndex, xIndex] = ampId

        ampAssociations = np.zeros((nAmps, nAmps), dtype=int)
        ampSides = np.full_like(ampAssociations, -1)

        for ampId in ampIds.ravel():
            neighbors, sides = self.getNeighbors(ampIds, ampId)
            interfaceWeights = (
                1
                if not self.config.applyWeights
                else np.array([self.interfaceLengthLookupBySide[side] for side in sides])
            )
            ampAssociations[ampId, neighbors] = -1 * interfaceWeights
            ampSides[ampId, neighbors] = sides
            ampAssociations[ampId, ampId] = -ampAssociations[ampId].sum()

        if ampAssociations.sum() != 0:
            raise RuntimeError("The `ampAssociations` array does not sum to zero.")

        if not np.all(ampAssociations == ampAssociations.T):
            raise RuntimeError("The `ampAssociations` is not symmetric about the diagonal.")

        self.log.debug("amp associations:\n%s", ampAssociations)
        self.log.debug("amp sides:\n%s", ampSides)

        return ampAssociations, ampSides

    def getNeighbors(self, ampIds, ampId):
        """Get the neighbor amplifiers and their sides for a given
        amplifier.

        Parameters
        ----------
        ampIds : `numpy.ndarray`
            Matrix with amp side association information.
        ampId : `int`
            The amplifier ID for which neighbor amplifiers and side IDs
            are to be found.

        Returns
        -------
        neighbors : `list` [`int`]
            List of neighbor amplifier IDs.
        sides : `list` [`int`]
            List of side IDs, with each ID corresponding to its respective
            neighbor amplifier.
        """
        m, n = ampIds.shape
        r, c = np.ravel(np.where(ampIds == ampId))
        neighbors, sides = [], []
        sideLookup = {
            0: (r + 1, c),
            1: (r, c + 1),
            2: (r - 1, c),
            3: (r, c - 1),
        }
        for side, (row, column) in sideLookup.items():
            if 0 <= row < m and 0 <= column < n:
                neighbors.append(ampIds[row][column])
                sides.append(side)
        return neighbors, sides

    def getAmpOffsets(self, im, amps, associations, sides):
        """Calculate the amp offsets for all amplifiers.

        Parameters
        ----------
        im : `lsst.afw.image._image.ImageF`
            Amplifier image to extract data from.
        amps : `list` [`lsst.afw.cameraGeom.Amplifier`]
            List of amplifier objects.
        associations : numpy.ndarray
            An N x N matrix containing amp association information, where N is
            the number of amplifiers.
        sides : numpy.ndarray
            An N x N matrix containing amp side information, where N is the
            number of amplifiers.

        Returns
        -------
        ampsOffsets : `numpy.ndarray`
            1D float array containing the calculated amp offsets for all
            amplifiers.
        """
        ampsOffsets = np.zeros(len(amps))
        ampsEdges = self.getAmpEdges(im, amps, sides)
        interfaceOffsetLookup = {}

        for ampId, ampAssociations in enumerate(associations):
            ampNeighbors = np.ravel(np.where(ampAssociations < 0))
            for ampNeighbor in ampNeighbors:
                ampSide = sides[ampId][ampNeighbor]
                interfaceWeight = (
                    1 if not self.config.applyWeights else self.interfaceLengthLookupBySide[ampSide]
                )
                edgeA = ampsEdges[ampId][ampSide]
                edgeB = ampsEdges[ampNeighbor][(ampSide + 2) % 4]
                if ampId < ampNeighbor:
                    interfaceOffset = self.getInterfaceOffset(ampId, ampNeighbor, edgeA, edgeB)
                    interfaceOffsetLookup[f"{ampId}{ampNeighbor}"] = interfaceOffset
                else:
                    interfaceOffset = -interfaceOffsetLookup[f"{ampNeighbor}{ampId}"]
                ampsOffsets[ampId] += interfaceWeight * interfaceOffset
        return ampsOffsets

    def getAmpEdges(self, im, amps, ampSides):
        """Calculate the amp edges for all amplifiers.

        Parameters
        ----------
        im : `lsst.afw.image._image.ImageF`
            Amplifier image to extract data from.
        amps : `list` [`lsst.afw.cameraGeom.Amplifier`]
            List of amplifier objects.
        ampSides : `numpy.ndarray`
            An N x N matrix containing amp side information, where N is the
            number of amplifiers.

        Returns
        -------
        ampEdges : `dict` [`int`, `dict` [`int`, `numpy.ndarray`]]
            A dictionary containing amp edge(s) for each amplifier,
            corresponding to one or more potential sides, where each edge is
            associated with a side. The outer dictionary has integer keys
            representing amplifier IDs, and the inner dictionary has integer
            keys representing side IDs for each amplifier and values that are
            1D arrays of floats representing the 1D medianified strips from the
            amp image, referred to as "amp edge":
            {ampID: {sideID: numpy.ndarray}, ...}
        """
        ampEdgeOuter = self.config.ampEdgeInset + self.config.ampEdgeWidth
        ampEdges = {}
        slice_map = {
            0: (slice(-ampEdgeOuter, -self.config.ampEdgeInset), slice(None)),
            1: (slice(None), slice(-ampEdgeOuter, -self.config.ampEdgeInset)),
            2: (slice(self.config.ampEdgeInset, ampEdgeOuter), slice(None)),
            3: (slice(None), slice(self.config.ampEdgeInset, ampEdgeOuter)),
        }
        for ampId, (amp, ampSides) in enumerate(zip(amps, ampSides)):
            ampEdges[ampId] = {}
            ampIm = im[amp.getBBox()].array
            # Loop over identified sides.
            for ampSide in ampSides:
                if ampSide < 0:
                    continue
                strip = ampIm[slice_map[ampSide]]
                # Catch warnings to prevent all-NaN slice RuntimeWarning.
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", r"All-NaN (slice|axis) encountered")
                    ampEdges[ampId][ampSide] = np.nanmedian(strip, axis=ampSide % 2)  # 1D medianified strip
        return ampEdges

    def getInterfaceOffset(self, ampIdA, ampIdB, edgeA, edgeB):
        """Calculate the amp offset for a given interface between two
        amplifiers.

        Parameters
        ----------
        ampIdA : int
            ID of the first amplifier.
        ampIdB : int
            ID of the second amplifier.
        edgeA : numpy.ndarray
            Amp edge for the first amplifier.
        edgeB : numpy.ndarray
            Amp edge for the second amplifier.

        Returns
        -------
        interfaceOffset : float
            The calculated amp offset value for the given interface between
            amps A and B.
        """
        interfaceId = f"{ampIdA}{ampIdB}"
        sctrl = StatisticsControl()
        # NOTE: Taking the difference with the order below fixes the sign flip
        # in the B matrix.
        edgeDiff = edgeA - edgeB
        window = int(self.config.ampEdgeWindowFrac * len(edgeDiff))
        # Compute rolling averages.
        edgeDiffSum = np.convolve(np.nan_to_num(edgeDiff), np.ones(window), "same")
        edgeDiffNum = np.convolve(~np.isnan(edgeDiff), np.ones(window), "same")
        edgeDiffAvg = edgeDiffSum / np.clip(edgeDiffNum, 1, None)
        edgeDiffAvg[np.isnan(edgeDiff)] = np.nan
        # Take clipped mean of rolling average data as amp offset value.
        interfaceOffset = makeStatistics(edgeDiffAvg, MEANCLIP, sctrl).getValue()
        # Perform a couple of do-no-harm safety checks:
        # a) The fraction of unmasked pixel rows is > ampEdgeMinFrac,
        # b) The absolute offset ADU value is < ampEdgeMaxOffset.
        ampEdgeGoodFrac = 1 - (np.sum(np.isnan(edgeDiffAvg)) / len(edgeDiffAvg))
        minFracFail = ampEdgeGoodFrac < self.config.ampEdgeMinFrac
        maxOffsetFail = np.abs(interfaceOffset) > self.config.ampEdgeMaxOffset
        if minFracFail or maxOffsetFail:
            interfaceOffset = 0
            if minFracFail:
                self.log.warning(
                    f"The fraction of unmasked pixels for amp interface {interfaceId} is below the threshold "
                    f"({ampEdgeGoodFrac:.2f} < {self.config.ampEdgeMinFrac}). Setting the interface offset "
                    f"to {interfaceOffset}."
                )
            if maxOffsetFail:
                self.log.warning(
                    "The absolute offset value exceeds the limit "
                    f"({np.abs(interfaceOffset):.2f} > {self.config.ampEdgeMaxOffset} ADU). Setting the "
                    f"interface offset to {interfaceOffset}."
                )
        self.log.debug(
            f"amp interface {interfaceId} : "
            f"viable edge difference frac = {ampEdgeGoodFrac}, "
            f"interface offset = {interfaceOffset:.3f}"
        )
        return interfaceOffset
