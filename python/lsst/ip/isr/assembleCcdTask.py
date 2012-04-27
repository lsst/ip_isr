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
    
    def assembleCcd(self, inExposure):
        """Assemble a CCD by copying data sections and (usually) trimming out the rest

        @param[in,out] inExposure   input exposure; the setTrimmed flag of the ccd device info may be modified
        @return assembled exposure
        """
        ccd = cameraGeom.cast_Ccd(inExposure.getDetector())
        if ccd is None:
            raise RuntimeError("Detector not a ccd")
        ccd.setTrimmed(self.config.doTrim)

        outExposure = afwImage.ExposureF(ccd.getAllPixels(self.config.doTrim))
        outMI = outExposure.getMaskedImage()
        inMI = inExposure.getMaskedImage()

        # Precompute amplifier information to avoid recomputing it for each image plane.
        ampInfoList = []
        for amp in ccd:
            outBBox = amp.getAllPixels(self.config.doTrim)
            if self.config.doTrim:
                inBBox = amp.getDiskDataSec()
            else:
                inBBox = amp.getDiskAllPixels()
            ampInfoList.append(pipeBase.Struct(
                amp = amp,
                outBBox = outBBox,
                inBBox = inBBox,
            ))
        
        # Process one image plane at a time, since that is probably better for cache performance
        # and is a requirement of amp.prepareAmpData.
        # Unfortunately image planes cannot be accessed by index or name,
        # so use getattr to access each image plane getter function
        for imagePlaneGetterName in ("getImage", "getMask", "getVariance"):
            outImage = getattr(outMI, imagePlaneGetterName)()
            inImage = getattr(inMI, imagePlaneGetterName)()
            for ampInfo in ampInfoList:
                outView = outImage.Factory(outImage, ampInfo.outBBox, afwImage.LOCAL)
                inView = inImage.Factory(inImage, ampInfo.inBBox, afwImage.PARENT)
                outView <<= ampInfo.amp.prepareAmpData(inView)

        outExposure.setDetector(ccd)

        self.postprocessExposure(outExposure=outExposure, inExposure=inExposure)
        
        return outExposure
    
    def assembleAmpList(self, ampExposureList):
        """Assemble a collection of amplifier exposures into a CCD

        @param[in,out]  ampExposureList collection of amp exposures to assemble;
                                        the setTrimmed flag of the ccd device info may be modified
        @return assembled exposure
        """
        ampExp0 = ampExposureList[0]
        amp0 = cameraGeom.cast_Amp(ampExp0.getDetector())
        if amp0 is None:
            raise RuntimeError("No amp detector found in first amp exposure")
        ccd = cameraGeom.cast_Ccd(amp0.getParent())
        if ccd is None:
            raise RuntimeError("No ccd detector found in amp detector")
        ccd.setTrimmed(self.config.doTrim)

        outExposure = afwImage.ExposureF(ccd.getAllPixels(self.config.doTrim))
        outMI = outExposure.getMaskedImage()
        
        # Precompute amplifier information to avoid recomputing it for each image plane.
        ampInfoList = []
        for ampExp in ampExposureList:
            amp = cameraGeom.cast_Amp(ampExp.getDetector())
            outBBox = amp.getAllPixels(self.config.doTrim)
            if self.config.doTrim:
                inBBox = amp.getDiskDataSec()
            else:
                inBBox = amp.getDiskAllPixels()
            ampInfoList.append(pipeBase.Struct(
                amp = amp,
                outBBox = outBBox,
                inBBox = inBBox,
                ampMaskedImage = ampExp.getMaskedImage(),
            ))
        
        # Process one image plane at a time, since that is probably better for cache performance
        # and since it is the only way amp.prepareAmpData works.
        # Unfortunately image planes cannot be accessed by index or name,
        # so use getattr to access each image plane getter function
        for imagePlaneGetterName in ("getImage", "getMask", "getVariance"):
            outImage = getattr(outMI, imagePlaneGetterName)()
            for ampInfo in ampInfoList:
                outView = outImage.Factory(outImage, ampInfo.outBBox, afwImage.LOCAL)
                inImage = getattr(ampInfo.ampMaskedImage, imagePlaneGetterName)()
                inView = inImage.Factory(inImage, ampInfo.inBBox, afwImage.PARENT)
                outView <<= ampInfo.amp.prepareAmpData(inView)

        outExposure.setDetector(ccd)

        self.postprocessExposure(outExposure=outExposure, inExposure=ampExposureList[0])
    
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
            ccd = cameraGeom.cast_Ccd(outExposure.getDetector())
            amp0 = cameraGeom.cast_Amp(ccd[0])
            if amp0 is None:
                raise RuntimeError("No amplifier detector information found")
            amp0.prepareWcsData(wcs)
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
