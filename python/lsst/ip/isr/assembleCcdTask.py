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
    
    def run(self, inExposure):
        """Assemble a CCD by trimming non-data areas

        @param[in,out] inExposure   input exposure; the setTrimmed flag of the ccd device info may be modified
        @return a pipe_base Struct with one field:
        - exposure: assembled exposure
        """
        outExposure = self.assemblePixels(inExposure=inExposure)

        self.setExposureComponents(outExposure=outExposure, inExposure=inExposure)

        self.display("assembledExposure", exposure=outExposure)
    
        return pipeBase.Struct(exposure=outExposure)
    
    def assemblePixels(self, inExposure):
        """Assemble CCD pixels

        @param[in,out]  inExposure  input exposure; the setTrimmed flag of the ccd device info may be modified
        @return         outExposure assembled exposure: just the pixel data and detector are set
        """
        ccd = cameraGeom.cast_Ccd(inExposure.getDetector())
        if ccd is None:
            raise RuntimeError("Detector not a ccd")
        ccd.setTrimmed(self.config.doTrim)

        outExposure = afwImage.ExposureF(ccd.getAllPixels(self.config.doTrim))
        outMI = outExposure.getMaskedImage()
        inMI = inExposure.getMaskedImage()
        for amp in ccd:
            outView = outMI.Factory(outMI, amp.getAllPixels(self.config.doTrim), afwImage.LOCAL)
            if self.config.doTrim:
                inBBox = amp.getDiskDataSec()
            else:
                inBBox = amp.getDiskAllPixels()
            inView = inMI.Factory(inMI, inBBox, afwImage.PARENT)
            outView <<= amp.prepareAmpData(inView)

        outExposure.setDetector(ccd)
        return outExposure

    def makeWcs(self, inExposure):
        """Create output WCS = input WCS offset for the datasec of the lower left amplifier.

        @param[in]      inExposure      input exposure
        """
        if inExposure.hasWcs():
            wcs = inExposure.getWcs()
            ccd = cameraGeom.cast_Ccd(inExposure.getDetector())
            amp0 = cameraGeom.cast_Amp(ccd[0])
            if amp0 is None:
                raise RuntimeError("No amplifier detector information found")
            amp0.prepareWcsData(wcs)
        else:
            wcs = None
        return wcs

    def setGain(self, outExposure):
        """Renormalize, if requested, and set gain metadata

        @param[in,out]  outExposure     assembled exposure:
                                        - scales the pixels if config.doRenorm is true
                                        - adds some gain keywords to the metadata
        """
        if outExposure.getMaskedImage().getVariance().getArray().max() == 0:
            raise RuntimeError("Can't calculate the effective gain since the variance plane is set to zero")
        ccd = cameraGeom.cast_Ccd(outExposure.getDetector())
        exposureMetadata = outExposure.getMetadata()
        gain = 0
        namps = 0
        for amp in ccd:
            gain += cameraGeom.cast_Amp(amp).getElectronicParams().getGain()
            namps += 1.
        gain /= float(namps)
        if self.config.doRenorm:
            mi = outExposure.getMaskedImage()
            mi *= gain
            exposureMetadata.set("GAIN", 1.0)
        medgain, meangain = calcEffectiveGain(outExposure.getMaskedImage())
        exposureMetadata.add("MEDGAIN", medgain)
        exposureMetadata.add("MEANGAIN", meangain)
        exposureMetadata.add("GAINEFF", medgain)
    
    def setExposureComponents(self, outExposure, inExposure):
        """Set exposure non-image attributes, including wcs and metadata
        
        @param[in,out]  outExposure assembled exposure:
                                    - removes unwanted keywords
                                    - sets calib, filter, and detector
        @param[in]      inExposure  input exposure
        """
        wcs = self.makeWcs(inExposure = inExposure)
        if wcs is not None:
            outExposure.setWcs(wcs)
        else:
            self.log.log(self.log.WARN, "No WCS found in input exposure")

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
