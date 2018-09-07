#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.image as afwImage
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsstDebug import getDebugFrame
from lsst.afw.display import getDisplay

__all__ = ["AssembleCcdTask"]


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

## @addtogroup LSST_task_documentation
## @{
## @page AssembleCcdTask
## @ref AssembleCcdTask_ "AssembleCcdTask"
## @copybrief AssembleCcdTask
## @}


class AssembleCcdTask(pipeBase.Task):
    """Assemble a set of amplifier images into a full detector size set of pixels.

    This task assembles sections of an image into a larger mosaic.  The sub-sections
    are typically amplifier sections and are to be assembled into a detector size pixel grid.
    The assembly is driven by the entries in the raw amp information.  The task can be configured
    to return a detector image with non-data (e.g. overscan) pixels included.  The task can also
    renormalize the pixel values to a nominal gain of 1.  The task also removes exposure metadata that
    has context in raw amps, but not in trimmed detectors (e.g. 'BIASSEC').

    The `lsst.pipe.base.cmdLineTask.CmdLineTask` command line task interface supports a
    flag -d to import debug.py from your PYTHONPATH; see
    "http://lsst-web.ncsa.illinois.edu/~buildbot/doxygen/x_masterDoxyDoc/base_debug.html"
    Using lsstDebug to control debugging output for more about debug.py files.

    The available debug variables in AssembleCcdTask are:

    display :
        A dictionary containing debug point names as keys with frame number as value. Valid keys are:

    assembledExposure :
        display assembled exposure

    Examples
    --------

    This code is in runAssembleTask.py in the examples directory, and can be run as e.g.

    >>> python examples/runAssembleTask.py

    Import the task.  There are other imports.  Read the source file for more info.

    >>> from lsst.ip.isr import AssembleCcdTask

    The above numbers can be changed.  The assumption that the readout corner is flipped on every other amp is
    hardcoded in createDetector.

    Run the assembler task

    .. code-block :: none

        def runAssembler():
        '''Run the task to assemble amps into a ccd'''
        # Create the assemble task with default config
        assembleConfig = AssembleCcdTask.ConfigClass()
        assembleTask = AssembleCcdTask(config=assembleConfig)
        frame = 0
        # The assemble task can operate in two ways:
        # 1. On a dictionary of amp size images
        # 2. On a single amp mosaic image
        # Both methods should produce the same output.
        for isPerAmp in (True, False):
            assemblyInput = exampleUtils.makeAssemblyInput(isPerAmp)
            assembledExposure = assembleTask.assembleCcd(assemblyInput)
            ds9.mtv(assembledExposure.getMaskedImage(), frame=frame, title="Per amp input is %s"%(isPerAmp))
            frame += 1

    To investigate the ip_isr_assemble_Debug, put something like

    .. code-block:: none

        import lsstDebug
        def DebugInfo(name):
            di = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
            if name == "lsst.ip.isr.assembleCcdTask":
                di.display = {'assembledExposure':2}
            return di

        lsstDebug.Info = DebugInfo

    into your debug.py file and run runAssembleTask.py with the @c --debug flag.

    Notes
    -----
    Conversion notes:
        Display code should be updated once we settle on a standard way of controlling what is displayed.
    """
    ConfigClass = AssembleCcdConfig
    _DefaultName = "assembleCcd"

    def __init__(self, **kwargs):
        """Initialize the AssembleCcdTask

        The keys for removal specified in the config are added to a default set:
        ('DATASEC', 'BIASSEC', 'TRIMSEC', 'GAIN')
        """
        pipeBase.Task.__init__(self, **kwargs)

        self.allKeysToRemove = ('DATASEC', 'BIASSEC', 'TRIMSEC', 'GAIN') + tuple(self.config.keysToRemove)

    def assembleCcd(self, assembleInput):
        """Assemble a set of amps into a single CCD size image

        Parameters
        ----------
        assembleInput : `lsst.afw.image.Exposures`
            Either a dictionary of amp lsst.afw.image.Exposures or a single
            lsst.afw.image.Exposure containing all raw
            amps.  If a dictionary of amp exposures,
            the key should be the amp name.

        Returns
        -------
        assembledCcd : `lsst.afw.image.Exposure`
            An lsst.afw.image.Exposure of the assembled amp sections.

        Raises
        ------
        TypeError
            with the following string:

            .. code-block:: none

                Expected either a dictionary of amp exposures or a single raw exposure
                The input exposures to be assembled do not adhere to the required format.


        RuntimeError
            with the following string:

            .. code-block:: none

                No ccd detector found
                The detector set on the input exposure is not set.

        """
        ccd = None
        if isinstance(assembleInput, dict):
            # assembleInput is a dictionary of amp name: amp exposure

            # Assume all amps have the same detector, so get the detector from an arbitrary amp
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
        # If we are returning an "untrimmed" image (with overscans and extended register) we
        # need to update the ampInfo table in the Detector as we've moved the amp images into
        # place in a single Detector image
        #
        if not self.config.doTrim:
            ccd = cameraGeom.makeUpdatedDetector(ccd)

        outExposure.setDetector(ccd)
        self.postprocessExposure(outExposure=outExposure, inExposure=getNextExposure(ccd[0]))

        return outExposure

    def postprocessExposure(self, outExposure, inExposure):
        """Set exposure non-image attributes, including wcs and metadata and display exposure (if requested)

        Call after assembling the pixels

        Parameters
        ----------
        outExposure :
            assembled exposure:

            - removes unwanted keywords
            - sets calib, filter, and detector

        inExposure :
            input exposure
        """
        self.setWcs(outExposure=outExposure, inExposure=inExposure)

        exposureMetadata = inExposure.getMetadata()
        for key in self.allKeysToRemove:
            if exposureMetadata.exists(key):
                exposureMetadata.remove(key)
        outExposure.setMetadata(exposureMetadata)

        # note: Calib is not copied, presumably because it is assumed unknown in raw data
        outExposure.setFilter(inExposure.getFilter())
        outExposure.getInfo().setVisitInfo(inExposure.getInfo().getVisitInfo())

        frame = getDebugFrame(self._display, "assembledExposure")
        if frame:
            getDisplay(frame).mtv(outExposure)

    def setWcs(self, outExposure, inExposure):
        """Set output WCS = input WCS, adjusted as required for datasecs not starting at lower left corner

        Parameters
        ----------
        outExposure :
            assembled exposure; wcs is set
        inExposure :
            input exposure
        """
        if inExposure.hasWcs():
            wcs = inExposure.getWcs()
            ccd = outExposure.getDetector()
            amp0 = ccd[0]
            if amp0 is None:
                raise RuntimeError("No amplifier detector information found")
            adjustedWcs = cameraGeomUtils.prepareWcsData(wcs, amp0, isTrimmed=self.config.doTrim)
            outExposure.setWcs(adjustedWcs)
