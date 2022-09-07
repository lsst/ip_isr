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

__all__ = ["MaskingConfig", "MaskingTask"]

from lsst.pex.config import Config, Field
from lsst.pipe.base import Task


class MaskingConfig(Config):
    doSpecificMasking = Field(
        dtype=bool,
        doc="",
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
