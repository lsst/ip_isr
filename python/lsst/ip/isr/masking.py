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
# import os

__all__ = ["MaskingConfig", "MaskingTask", "DECamEdgeBleedMaskConfig", "DECamEdgeBleedMaskTask"]

from lsst.pex.config import Config, Field
from lsst.pipe.base import Task

from . import isrFunctions


class MaskingConfig(Config):
    doSpecificMasking = Field(
        dtype=bool,
        doc="Masking configuration.",
        default=False,
    )


class MaskingTask(Task):
    """Perform extra masking for detector issues such as ghosts and glints.
    """
    ConfigClass = MaskingConfig
    _DefaultName = "isrMasking"

    def run(self, exposure):
        """Mask a known bad region of an exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to construct detector-specific masks for.

        Returns
        -------
        status : scalar
            This task is currently not implemented, and should be
            retargeted by a camera specific version.
        """
        return


class DECamEdgeBleedMaskConfig(Config):
    satMinArea = Field(
        dtype=int,
        doc="Minimum area (pixels) of a saturated footprint to be an edge bleed candidate.",
        default=10000,
    )
    satMaxArea = Field(
        dtype=int,
        doc="Maximum area (pixels, exclusive) of a saturated footprint to be an edge bleed candidate.",
        default=100000,
    )
    approachRows = Field(
        dtype=int,
        doc="Saturated footprint must come within this many rows of the read edge.",
        default=20,
    )
    nSigma = Field(
        dtype=float,
        doc="A pixel is counted as depressed if it is below sky by more than this many sigma.",
        default=5.0,
    )
    nRowsCheck = Field(
        dtype=int,
        doc="Number of rows from the read edge used to confirm an edge bleed.",
        default=20,
    )
    minLowPixelsPerRow = Field(
        dtype=int,
        doc=("Mean depressed pixels per row over the check rows required (exceeded) to confirm an "
             "edge bleed. Read-edge rows with fewer usable pixels than this are skipped before the "
             "check rows start."),
        default=30,
    )
    minLowPixelsExtent = Field(
        dtype=int,
        doc=("Depressed pixels a row must have (more than this) to count toward the edge bleed "
             "height. Should be smaller than minLowPixelsPerRow."),
        default=10,
    )
    marginFraction = Field(
        dtype=float,
        doc=("Extra rows masked beyond the measured edge bleed height, as a fraction of that height "
             "(plus one row)."),
        default=0.125,
    )
    saturatedMaskName = Field(
        dtype=str,
        doc="Name of mask plane holding saturated pixels; must match the parent ISR task.",
        default="SAT",
    )

    def validate(self):
        super().validate()
        if self.minLowPixelsExtent >= self.minLowPixelsPerRow:
            raise ValueError("minLowPixelsExtent must be smaller than minLowPixelsPerRow.")


class DECamEdgeBleedMaskTask(MaskingTask):
    """Mask DECam-style edge bleeds: depressed rows next to the read register
    below a large saturated star.

    Intended as a retarget for ``IsrTaskConfig.masking`` (enabled with
    ``doCameraSpecificMasking``).  Must run after saturation masking.
    """
    ConfigClass = DECamEdgeBleedMaskConfig
    _DefaultName = "decamEdgeBleedMask"

    def run(self, exposure):
        """Mask edge bleeds in every amplifier of an exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Assembled exposure with saturation already masked.  The mask
            plane is modified in place.
        """
        isrFunctions.maskDECamEdgeBleed(
            exposure,
            satMinArea=self.config.satMinArea,
            satMaxArea=self.config.satMaxArea,
            approachRows=self.config.approachRows,
            nSigma=self.config.nSigma,
            nRowsCheck=self.config.nRowsCheck,
            minLowPixelsPerRow=self.config.minLowPixelsPerRow,
            minLowPixelsExtent=self.config.minLowPixelsExtent,
            marginFraction=self.config.marginFraction,
            saturatedMaskName=self.config.saturatedMaskName,
            log=self.log,
        )
