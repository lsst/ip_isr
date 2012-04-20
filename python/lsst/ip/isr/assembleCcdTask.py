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
    setGain = pexConfig.ConfigField(dtype = bool, default = True,
        doc = "set gain?")
    doRenorm = pexConfig.ConfigField(dtype = bool, default = True,
        doc = "renormalize to a gain of 1? (ignored if setGain false)")
    keysToRemove = pexConfig.ListField(dtype = str, default = (),
        doc = "FITS headers to remove (in addition to DATASEC, BIASSEC, TRIMSEC and perhaps GAIN)")


class AssembleCcdTask(pipeBase.Task):
    """Assemble a CCD
    """
    ConfigClass = AssembleCcdConfig
    _DefaultName = "assembleCcd"
    
    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        
        if self.config.setGain and self.config.doRenorm:
            #Don't want to remove GAIN keyword as it is set by the renormalization.
            self.ktr = ['DATASEC', 'BIASSEC', 'TRIMSEC']
        else:
            self.ktr = ['DATASEC', 'BIASSEC', 'TRIMSEC', 'GAIN']
        for k in self.config.keysToRemove:
            self.ktr.append(k)
    
    def run(self, inExposure):
        """Assemble a CCD by trimming non-data areas

        @param[in]      inExposure   input exposure
        """
        ccd = cameraGeom.cast_Ccd(inExposure.getDetector())
        llAmp = cameraGeom.cast_Amp(ccd[0]) # lower left amplifier

        if ccd is None or not isinstance(ccd, cameraGeom.Ccd) or \
               llAmp is None or not isinstance(llAmp, cameraGeom.Amp):
            raise RuntimeError("Detector in exposure does not match calling pattern")
        
        # convert CCD for assembled exposure
        ccd.setTrimmed(True)
        
        outExposure = self.assemblePixels(
            inExposure = inExposure,
            ccd = ccd,
        )

        self.setExposureComponents(
            outExposure = outExposure,
            inExposure = inExposure,
            ccd = ccd,
            llAmp = llAmp,
        )

        if self.config.setGain:
            if ccdVariance.getArray().max() == 0:
                raise("Can't calculate the effective gain since the variance plane is set to zero")
            self.setGain(
                outExposure = outExposure,
                ccd = ccd,
            )

        self.display("assembledExposure", exposure = outExposure)
    
        return pipeBase.Struct(
            exposure = outExposure
        )
    
    def assemblePixels(self, inExposure, ccd):
        """Assemble CCD pixels

        @param[in]      inExposure  input exposure
        @param[in]      ccd         device info for assembled exposure
        @return         outExposure assembled exposure (just the pixels data is set)
        """
        ccdImage = cameraGeomUtils.makeImageFromCcd(
            ccd = ccd,
            imageSource = GetCcdImageData(inMaskedImage.getImage()),
            imageFactory = inMaskedImage.getImage().Factory,
            bin = False,
        )
        ccdVariance = cameraGeomUtils.makeImageFromCcd(
            ccd = ccd,
            imageSource = GetCcdImageData(inMaskedImage.getVariance()),
            imageFactory = afwImage.ImageF,
            bin = False,
        )
        ccdMask = cameraGeomUtils.makeImageFromCcd(
            ccd = ccd,
            imageSource = GetCcdImageData(inMaskedImage.getMask()),
            imageFactory = afwImage.MaskU,
            bin = False,
        )
        mi = afwImage.makeMaskedImage(ccdImage, ccdMask, ccdVariance)
        return afwImage.makeExposure(mi)

    def makeWcs(self, inExposure, llAmp):
        """Create output WCS = input WCS offset for the datasec of the lower left amplifier.

        @param[in]      inExposure   input exposure
        @param[in]      llAmp        device info for lower left amplifier of input exposure
        """
        if inExposure.hasWcs():
            wcs = inExposure.getWcs()
            llAmp.prepareWcsData(wcs)
            #shift the reference pixel to account for the location of amp 0 in the ccd.
            dataSec = llAmp.getDataSec()
            wcs.shiftReferencePixel(dataSec.getMinX(), dataSec.getMinY())
        else:
            wcs = None
        return wcs

    def setGain(self, outExposure, ccd):
        """Renormalize, if requested, and set gain metadata

        @param[in,out]  outExposure     assembled exposure:
                                        - modifies the pixels if config.doRenorm is true
                                        - adds some gain keywords to the metadata
        @param[in]      llAmp           device info for lower left amplifier of input exposure
        """
        exposureMetadata = outExposure.getMetadata()
        gain = 0
        namps = 0
        for a in ccd:
            gain += cameraGeom.cast_Amp(a).getElectronicParams().getGain()
            namps += 1.
        gain /= float(namps)
        if self.config.doRenorm:
            mi = outExposure.getMaskedImage()
            mi *= gain
            exposureMetadata.set("GAIN", 1.0)
        (medgain, meangain) = calcEffectiveGain(outExposure.getMaskedImage())
        exposureMetadata.add("MEDGAIN", medgain)
        exposureMetadata.add("MEANGAIN", meangain)
        exposureMetadata.add("GAINEFF", medgain)
    
    def setExposureComponents(self, outExposure, inExposure, ccd, llAmp):
        """Set exposure non-image attributes, such as wcs and metadata
        
        @param[in,out]  outExposure assembled exposure:
                                    - removes unwanted keywords
                                    - sets calib, filter, and detector
        @param[in]      inExposure  input exposure
        @param[in]      ccd         device info for assembled exposure
        @param[in]      llAmp       device info for lower left amplifier of input exposure
        """
        wcs = self.makeWcs(inExposure, llAmp)
        if wcs is not None:
            outExposure.setWcs(wcs)

        exposureMetadata = inExposure.getMetadata().copy()
        for k in self.ktr:
            if exposureMetadata.exists(k):
                exposureMetadata.remove(k)
        outExposure.setMetadata(exposureMetadata)

        inCalib = inExposure.getCalib()
        outCalib = outExposure.getCalib()
        outCalib.setExptime(calib.getExptime())
        outCalib.setMidTime(calib.getMidTime())

        outExposure.setFilter(inExposure.getFilter())

        outExposure.setDetector(ccd)

class GetCcdImageData(cameraGeomUtils.GetCcdImage):
    def __init__(self, image, isTrimmed=True):
        self.image = image
        self.isRaw = True
        self.isTrimmed = isTrimmed

    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.ImageF):
        if self.isTrimmed:
            bbox = amp.getDiskDataSec()
        else:
            bbox = amp.getDiskAllPixels()
        subImage = imageFactory(self.image, bbox, afwImage.PARENT)
        return amp.prepareAmpData(subImage)
