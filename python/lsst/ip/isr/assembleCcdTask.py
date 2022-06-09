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

    This task assembles sections of an image into a larger mosaic.
    The sub-sections are typically amplifier sections and are to be
    assembled into a detector size pixel grid. The assembly is driven
    by the entries in the raw amp information.  The task can be
    configured to return a detector image with non-data
    (e.g. overscan) pixels included.  The task can also renormalize
    the pixel values to a nominal gain of 1.  The task also removes
    exposure metadata that has context in raw amps, but not in trimmed
    detectors.  This includes the set ('DATASEC', 'BIASSEC',
    'TRIMSEC', 'GAIN'), as well as any other keys specificied in the
    task config..

    Parameters
    ----------
    **kwargs : Any
        Keyword parameters.

    Notes
    -----

    Debug plotting of the assembled image can be enabled by setting
    the ``assembledExposure`` keyword in the lsstDebug `debug.display`
    dictionary to True.
    """

    ConfigClass = AssembleCcdConfig
    _DefaultName = "assembleCcd"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)

        self.allKeysToRemove = ('DATASEC', 'BIASSEC', 'TRIMSEC', 'GAIN') + tuple(self.config.keysToRemove)

    def assembleCcd(self, assembleInput):
        """Assemble a set of amps into a single CCD size image

        Parameters
        ----------
        assembleInput : `lsst.afw.image.Exposure` or `dict`
            If an `~lsst.afw.image.Exposure` is passed, it is assumed
            to contain all the raw amps as regions within that
            exposure.  If the input is a `dict`, then the keys are the
            amplifier names, and the values are
            `~lsst.afw.image.Exposure` entries containing the data for
            that amplifier

        Returns
        -------
        outExposure : `lsst.afw.image.Exposure`
            The assembled exposure.

        Raises
        ------
        TypeError
            Raised if the input is not an `lsst.afw.image.Exposure` or a `dict`.
        RuntimeError
            Raised if no input detector can be identified.
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

        # If we are returning an "untrimmed" image (with overscans and
        # extended register) we need to update the ampInfo table in the
        # Detector as we've moved the amp images into
        # place in a single Detector image
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
