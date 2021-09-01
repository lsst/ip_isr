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

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.meas.algorithms import (SubtractBackgroundTask, SourceDetectionTask)


class AmpOffsetConfig(pexConfig.Config):
    """Configuration parameters for AmpOffsetTask.
    """
    ampEdgeInset = pexConfig.Field(
        doc="Number of pixels the amp edge strip is inset from the amp edge. A thin strip of pixels running "
        "parallel to the edge of the amp is used to characterize the average flux level at the amp edge.",
        dtype=int,
        default=5,
    )
    ampEdgeWidth = pexConfig.Field(
        doc="Pixel width of the amp edge strip, starting at ampEdgeInset and extending inwards.",
        dtype=int,
        default=64,
    )
    ampEdgeMinFrac = pexConfig.Field(
        doc="Minimum allowed fraction of viable pixel rows along an amp edge. No amp offset estimate will be "
        "generated for amp edges that do not have at least this fraction of unmasked pixel rows.",
        dtype=float,
        default=0.5,
    )
    ampEdgeMaxOffset = pexConfig.Field(
        doc="Maximum allowed amp offset ADU value. If a measured amp offset value is larger than this, the "
        "result will be discarded and therefore not used to determine amp pedestal corrections.",
        dtype=float,
        default=5.0,
    )
    ampEdgeWindow = pexConfig.Field(
        doc="Pixel size of the sliding window used to generate rolling average amp offset values.",
        dtype=int,
        default=512,
    )
    doBackground = pexConfig.Field(
        doc="Estimate and subtract background prior to amp offset estimation?",
        dtype=bool,
        default=True,
    )
    background = pexConfig.ConfigurableField(
        doc="An initial background estimation step run prior to amp offset calculation.",
        target=SubtractBackgroundTask,
    )
    doDetection = pexConfig.Field(
        doc="Detect sources and update cloned exposure prior to amp offset estimation?",
        dtype=bool,
        default=True,
    )
    detection = pexConfig.ConfigurableField(
        doc="Source detection to add temporary detection footprints prior to amp offset calculation.",
        target=SourceDetectionTask,
    )


class AmpOffsetTask(pipeBase.Task):
    """Calculate and apply amp offset corrections to an exposure.
    """
    ConfigClass = AmpOffsetConfig
    _DefaultName = "isrAmpOffset"

    def __init__(self, *args, **kwargs):
        super().__init__(**kwargs)
        # always load background subtask, even if doBackground=False;
        # this allows for default plane bit masks to be defined
        self.makeSubtask("background")
        if self.config.doDetection:
            self.makeSubtask("detection")

    def run(self, exposure):
        """Calculate amp offset values, determine corrective pedestals for each
        amp, and update the input exposure in-place. This task is currently not
        implemented, and should be retargeted by a camera specific version.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to be corrected for any amp offsets.
        """
        raise NotImplementedError("Amp offset task should be retargeted by a camera specific version.")
