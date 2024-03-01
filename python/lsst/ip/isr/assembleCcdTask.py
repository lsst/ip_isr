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

__all__ = ["AssembleCcdTask", "AssembleCcdConfig"]

import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display as afwDisplay
import lsst.afw.image as afwImage
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsstDebug import getDebugFrame


class AssembleCcdConfig(pexConfig.Config):
    doTrim = pexConfig.Field(
        doc="trim out non-data regions?",
        dtype=bool,
        default=True,
    )
    keysToRemove = pexConfig.ListField(
        doc="FITS headers to remove (in addition to DATASEC, BIASSEC, TRIMSEC and perhaps GAIN)",
        dtype=str,
        default=(),
    )


class AssembleCcdTask(pipeBase.Task):
    """Assemble a set of amplifier images into a full detector size set of
    pixels.

    The keys for removal specified in
    `lsst.ip.isr.AssembleCcdConfig.keysToRemove` are added to a default set:
    ('DATASEC', 'BIASSEC', 'TRIMSEC', 'GAIN').
    """
    ConfigClass = AssembleCcdConfig
    _DefaultName = "assembleCcd"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)

        self.allKeysToRemove = ('DATASEC', 'BIASSEC', 'TRIMSEC', 'GAIN') + tuple(self.config.keysToRemove)

    def assembleCcd(self, assembleInput):
        """Assemble a set of amps into a single CCD size image.

        Parameters
        ----------
        assembleInput : `dict` [`str`, `lsst.afw.image.Exposure`] or \
                        `lsst.afw.image.Exposure`
            Either a dictionary of amp exposures, or a single exposure
            containing all raw amps. If a dictionary of amp exposures, the key
            should be the amp name.

        Returns
        -------
        assembledCcd : `lsst.afw.image.Exposure`
            An exposure of the assembled amp sections.

        Raises
        ------
        TypeError
            Raised if the input exposures to be assembled do not adhere to the
            required format.
        RuntimeError
            Raised if the detector set on the input exposure is not set.
        """
        ccd = None
        if isinstance(assembleInput, dict):
            # assembleInput is a dictionary of amp name: amp exposure

            # Assume all amps have the same detector, so get the detector from
            # an arbitrary amp.
            ccd = next(iter(assembleInput.values())).getDetector()

            def getNextExposure(amp):
                return assembleInput[amp.getName()]
        elif hasattr(assembleInput, "getMaskedImage"):
            # assembleInput is a single exposure
            ccd = assembleInput.getDetector()

            def getNextExposure(amp):
                return assembleInput
        else:
            raise TypeError("Expected either a dictionary of amp exposures or a single raw exposure")

        if ccd is None:
            raise RuntimeError("No ccd detector found")

        if not self.config.doTrim:
            outBox = cameraGeomUtils.calcRawCcdBBox(ccd)
        else:
            outBox = ccd.getBBox()
        outExposure = afwImage.ExposureF(outBox)
        outMI = outExposure.getMaskedImage()

        if self.config.doTrim:
            assemble = cameraGeom.assembleAmplifierImage
        else:
            assemble = cameraGeom.assembleAmplifierRawImage

        for amp in ccd:
            inMI = getNextExposure(amp).getMaskedImage()
            assemble(outMI, inMI, amp)
        #
        # If we are returning an "untrimmed" image (with overscans and
        # extended register) we need to update the ampInfo table in the
        # Detector as we've moved the amp images into
        # place in a single Detector image
        #
        if not self.config.doTrim:
            ccd = cameraGeom.makeUpdatedDetector(ccd)

        outExposure.setDetector(ccd)
        self.postprocessExposure(outExposure=outExposure, inExposure=getNextExposure(ccd[0]))

        return outExposure

    def postprocessExposure(self, outExposure, inExposure):
        """Set exposure non-image attributes, including wcs and metadata and
        display exposure (if requested).

        Call after assembling the pixels.

        Parameters
        ----------
        outExposure : `lsst.afw.image.Exposure`
            The exposure to modify by copying metadata (after removing unwanted
            keywords), wcs, filter, and detector from ``inExposure``.
        inExposure : `lsst.afw.image.Exposure`
            The input exposure providing metadata, wcs, filter, and detector.
        """
        if inExposure.hasWcs():
            outExposure.setWcs(inExposure.getWcs())

        exposureMetadata = inExposure.getMetadata()
        for key in self.allKeysToRemove:
            if exposureMetadata.exists(key):
                exposureMetadata.remove(key)
        outExposure.setMetadata(exposureMetadata)

        # note: don't copy PhotoCalib, because it is assumed to be unknown in
        # raw data.
        outExposure.info.id = inExposure.info.id
        outExposure.setFilter(inExposure.getFilter())
        outExposure.getInfo().setVisitInfo(inExposure.getInfo().getVisitInfo())

        frame = getDebugFrame(self._display, "assembledExposure")
        if frame:
            afwDisplay.Display(frame=frame).mtv(outExposure, title="postprocessExposure")
