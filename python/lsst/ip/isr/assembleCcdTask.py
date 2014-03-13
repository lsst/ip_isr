# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
from .isr import calcEffectiveGain

__all__ = ["AssembleCcdTask"]

class AssembleCcdConfig(pexConfig.Config):
    setGain = pexConfig.Field(
        doc = "set gain?",
        dtype = bool,
        default = True,
    )
    doRenorm = pexConfig.Field(
        doc = "renormalize to a gain of 1? (ignored if setGain false)",
        dtype = bool,
        default = True,
    )
    doTrim = pexConfig.Field(
        doc = "trim out non-data regions?",
        dtype = bool,
        default = True,
    )
    keysToRemove = pexConfig.ListField(
        doc = "FITS headers to remove (in addition to DATASEC, BIASSEC, TRIMSEC and perhaps GAIN)",
        dtype = str,
        default = (),
    )


class AssembleCcdTask(pipeBase.Task):
    """Assemble a CCD
    """
    ConfigClass = AssembleCcdConfig
    _DefaultName = "assembleCcd"
    
    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        
        self.allKeysToRemove = ('DATASEC', 'BIASSEC', 'TRIMSEC', 'GAIN') + tuple(self.config.keysToRemove)
    
    def assembleCcd(self, assembleInput):
        """Assemble a set of amps into a single CCD size image
        @param assembleInput -- Either a dictionary of amp exposures or a single exposure containing all raw
                                amps.  If a dictionary of amp exposures, the key should be the amp name.
        """
        ccd = None
        if hasattr(assembleInput, "has_key"):
            # Get a detector object for this set of amps
            ccd = assembleInput.itervalues().next().getDetector()
            # Sent a dictionary of input exposures, assume one amp per key keyed on amp name
            def getNextExposure(amp):
                return assembleInput[amp.getName()]
        elif hasattr(assembleInput, "getMaskedImage"):
            ccd = assembleInput.getDetector()
            # A single exposure was sent.  Use this to assemble.
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
        outExposure.setDetector(ccd)
        self.postprocessExposure(outExposure=outExposure, inExposure=getNextExposure(ccd[0]))
    
        return outExposure
    
    def postprocessExposure(self, outExposure, inExposure):
        """Set exposure non-image attributes, including wcs and metadata and display exposure (if requested)
        
        Call after assembling the pixels
        
        @param[in,out]  outExposure assembled exposure:
                                    - removes unwanted keywords
                                    - sets calib, filter, and detector
        @param[in]      inExposure  input exposure
        """
        self.setWcs(outExposure = outExposure, inExposure = inExposure)

        exposureMetadata = inExposure.getMetadata()
        for key in self.allKeysToRemove:
            if exposureMetadata.exists(key):
                exposureMetadata.remove(key)
        outExposure.setMetadata(exposureMetadata)

        if self.config.setGain:
            self.setGain(outExposure = outExposure)

        inCalib = inExposure.getCalib()
        outCalib = outExposure.getCalib()
        outCalib.setExptime(inCalib.getExptime())
        outCalib.setMidTime(inCalib.getMidTime())

        outExposure.setFilter(inExposure.getFilter())

        self.display("assembledExposure", exposure=outExposure)

    def setWcs(self, outExposure, inExposure):
        """Set output WCS = input WCS, adjusted as required for datasecs not starting at lower left corner

        @param[in,out]  outExposure     assembled exposure; wcs is set
        @param[in]      inExposure      input exposure
        """
        if inExposure.hasWcs():
            wcs = inExposure.getWcs()
            ccd = outExposure.getDetector()
            amp0 = ccd[0]
            if amp0 is None:
                raise RuntimeError("No amplifier detector information found")
            cameraGeomUtils.prepareWcsData(wcs, amp0)
            outExposure.setWcs(wcs)
        else:
            self.log.log(self.log.WARN, "No WCS found in input exposure")

    def setGain(self, outExposure):
        """Renormalize, if requested, and set gain metadata

        @param[in,out]  outExposure     assembled exposure:
                                        - scales the pixels if config.doRenorm is true
                                        - adds some gain keywords to the metadata
        """
        if outExposure.getMaskedImage().getVariance().getArray().max() == 0:
            raise RuntimeError("Can't calculate the effective gain since the variance plane is set to zero")
        ccd = outExposure.getDetector()
        exposureMetadata = outExposure.getMetadata()
        gain = 0.
        namps = 0
        for amp in ccd:
            gain += amp.getGain()
            namps += 1
        gain /= namps
        if self.config.doRenorm:
            mi = outExposure.getMaskedImage()
            mi *= gain
            exposureMetadata.set("GAIN", 1.0)
        medgain, meangain = calcEffectiveGain(outExposure.getMaskedImage())
        exposureMetadata.add("MEDGAIN", medgain)
        exposureMetadata.add("MEANGAIN", meangain)
        exposureMetadata.add("GAINEFF", medgain)
