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
from abc import ABC, abstractmethod
from typing import Optional

from lsst.pex.config import Config, Field, ListField
from lsst.pipe.base import Task
from lsst.geom import Angle


class StrayLightConfig(Config):
    doRotatorAngleCorrection = Field(
        dtype=bool,
        doc="",
        default=False,
    )
    # TODO DM-28093: change the doc to specify that these are physical labels
    filters = ListField(
        dtype=str,
        doc="Filters that need straylight correction.",
        default=[],
    )


class StrayLightTask(Task):
    """Remove stray light from instruments.

    This is a dummy task to be retargeted with an camera-specific version.
    """
    ConfigClass = StrayLightConfig
    _DefaultName = "isrStrayLight"

    def readIsrData(self, dataRef, rawExposure):
        """Read and return calibration products relevant for correcting
        stray light in the given exposure.

        Parameters
        ----------
        dataRef : `daf.persistence.butlerSubset.ButlerDataRef`
            Butler reference of the detector data to be processed
        rawExposure : `afw.image.Exposure`
            The raw exposure that will later be corrected with the
            retrieved calibration data; should not be modified in this
            method.

        Returns
        -------
        straylightData : `object`, optional
            An opaque object that should be passed as the second argument to
            the `run` method.  If `None`, no stray light correction will be
            performed for the given image.  Any other object (e.g. `True`)
            may be used to signal that stray light correction should be
            performed even if there is nothing to read.

        Notes
        -----
        This method will be called only when `IsrTask` is run by the Gen2
        Middleware (i.e. CmdLineTask).
        """
        return None

    def check(self, exposure):
        """Check if stray light correction should be run.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to correct.
        """
        return False

    def run(self, exposure, strayLightData):
        """Correct stray light.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
           Exposure to correct.
        strayLightData : `object`, optional
            An opaque object that contains any calibration data used to
            correct for stray light.
        """
        raise NotImplementedError("Must be implemented by subclasses.")

    def checkFilter(self, exposure):
        """Check whether we should fringe-subtract the science exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to check the filter of.

        Returns
        -------
        needsFringe : `bool`
            If True, then the exposure has a filter listed in the
            configuration, and should have the fringe applied.
        """
        return exposure.getFilter().getName() in self.config.filters


class StrayLightData(ABC):
    """An abstract base class for rotator-dependent stray light information.
    """

    @abstractmethod
    def evaluate(self, angle_start: Angle, angle_end: Optional[Angle] = None):
        """Get a stray light array for a range of rotator angles.

        Parameters
        ----------
        angle_begin : `float`
            Instrument rotation angle at the start of the exposure.
        angle_end : `float`, optional
            Instrument rotation angle at the end of the exposure.
            If not provided, the returned array will reflect a snapshot at
            `angle_start`.

        Returns
        -------
        array : `numpy.ndarray`
            A stray-light background image for this exposure.
        """
        raise NotImplementedError("Must be implemented by subclasses.")
