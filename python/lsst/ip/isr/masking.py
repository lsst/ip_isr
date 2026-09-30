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
        doc="Minimum area (pixels) of a saturated footprint to be considered.",
        default=10000,
    )
    approachRows = Field(
        dtype=int,
        doc=("A saturated footprint must come within this many rows of the read edge to trigger the "
             "edge bleed check."),
        default=20,
    )
    nSigma = Field(
        dtype=float,
        doc="A pixel is \"low\" if it is more than this many sigma below sky.",
        default=5.0,
    )
    nRowsCheck = Field(
        dtype=int,
        doc="Number of rows from the read edge to check for an edge bleed dip.",
        default=20,
    )
    nRowsCheckShort = Field(
        dtype=int,
        doc="A second, shorter window from the read edge to check for a short edge bleed dip.",
        default=5,
    )
    minUsablePixelsPerRow = Field(
        dtype=int,
        doc=("Minimum number of pixels not saturated, bad or missing to search for dips. "
             "Rows with fewer pixels may be \"blocked\" due to a drained register."),
        default=30,
    )
    minLowFraction = Field(
        dtype=float,
        doc=("Mean fraction of usable pixels that are low, over either nRowsCheck or nRowsCheckShort, "
             "to confirm a dip."),
        default=0.03,
    )
    minLowFractionExtent = Field(
        dtype=float,
        doc=("Fraction of usable pixels that are low that each row must exceed to count toward the "
             "bleed height; the scan stops after five consecutive rows that do not. "
             "Should not exceed minLowFraction."),
        default=0.01,
    )
    minBlockedRows = Field(
        dtype=int,
        doc="Minimum number of blocked rows beyond the normal border at the read edge to be masked if found.",
        default=5,
    )
    marginFraction = Field(
        dtype=float,
        doc=("Extra rows masked beyond the measured height, as a fraction of that height "
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
        if not 0.0 < self.minLowFraction < 1.0 or not 0.0 < self.minLowFractionExtent < 1.0:
            raise ValueError("minLowFraction and minLowFractionExtent must be between 0 and 1.")
        if self.minLowFractionExtent > self.minLowFraction:
            raise ValueError("minLowFractionExtent must not exceed minLowFraction.")


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
            approachRows=self.config.approachRows,
            nSigma=self.config.nSigma,
            nRowsCheck=self.config.nRowsCheck,
            nRowsCheckShort=self.config.nRowsCheckShort,
            minUsablePixelsPerRow=self.config.minUsablePixelsPerRow,
            minLowFraction=self.config.minLowFraction,
            minLowFractionExtent=self.config.minLowFractionExtent,
            minBlockedRows=self.config.minBlockedRows,
            marginFraction=self.config.marginFraction,
            saturatedMaskName=self.config.saturatedMaskName,
            log=self.log,
        )
