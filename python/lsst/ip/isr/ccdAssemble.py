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
import lsstDebug

import lsst.pex.policy as pexPolicy
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display.ds9 as ds9
import isr
import os,sys,eups,math

class singleImageFactory(cameraGeomUtils.GetCcdImage):
    def __init__(self, exposure, isTrimmed=True):
        self.exposure = exposure
        self.isRaw = True
        self.isTrimmed = isTrimmed
    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.ImageF):
        if self.isTrimmed:
            img = imageFactory(self.exposure.getMaskedImage().getImage(),
                    amp.getDiskDataSec(), afwImage.PARENT)
        else:
            img = imageFactory(self.exposure.getMaskedImage().getImage(),
                    amp.getDiskAllPixels(), afwImage.PARENT)
        return amp.prepareAmpData(img)

class singleMaskFactory(cameraGeomUtils.GetCcdImage):
    def __init__(self, exposure, isTrimmed=True):
        self.exposure = exposure
        self.isRaw = True
        self.isTrimmed = isTrimmed
    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.MaskU):
	if self.isTrimmed is True:
	    img = imageFactory(self.exposure.getMaskedImage().getMask(),
		    amp.getDiskDataSec(), afwImage.PARENT)
	else:
	    img = imageFactory(self.exposuree.getMaskedImage().getMask(),
		    amp.getDiskAllPixels(), afwImage.PARENT)
	return amp.prepareAmpData(img)

class singleVarianceFactory(cameraGeomUtils.GetCcdImage):
    def __init__(self, exposure, isTrimmed=True):
        self.exposure = exposure
        self.isRaw = True
        self.isTrimmed = isTrimmed
    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.ImageF,
            isTrimmed=True):
	if self.isTrimmed:
	    img = imageFactory(self.exposure.getMaskedImage().getVariance(),
		    amp.getDiskDataSec(), afwImage.PARENT)
	else:
	    img = imageFactory(self.exposure.getMaskedImage().getVariance(),
		    amp.getDiskAllPixels(), afwImage.PARENT)
	return amp.prepareAmpData(img)

class listImageFactory(cameraGeomUtils.GetCcdImage):
    def __init__(self, exposures, isTrimmed=True):
        self.exposures = exposures
        self.isRaw = True
        self.isTrimmed = isTrimmed
    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.ImageF):
        for e in self.exposures:             
            if e.getDetector().getId() == amp.getId():
                if self.isTrimmed:
                    img = imageFactory(e.getMaskedImage().getImage(),
                            amp.getDiskDataSec(), afwImage.PARENT)
                else:
                    img = imageFactory(e.getMaskedImage().getImage(),
                            amp.getDiskAllPixels(), afwImage.PARENT)
                return amp.prepareAmpData(img)
        return None

class listMaskFactory(cameraGeomUtils.GetCcdImage):
    def __init__(self, exposures, isTrimmed=True):
        self.exposures = exposures
        self.isRaw = True
        self.isTrimmed = isTrimmed
    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.MaskU):
        for e in self.exposures:             
            if e.getDetector().getId() == amp.getId():
                if self.isTrimmed is True:
                    img = imageFactory(e.getMaskedImage().getMask(),
                            amp.getDiskDataSec(), afwImage.PARENT)
                else:
                    img = imageFactory(e.getMaskedImage().getMask(),
                            amp.getDiskAllPixels(), afwImage.PARENT)
                return amp.prepareAmpData(img)
        return None

class listVarianceFactory(cameraGeomUtils.GetCcdImage):
    def __init__(self, exposures, isTrimmed=True):
        self.exposures = exposures
        self.isRaw = True
        self.isTrimmed = isTrimmed
    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.ImageF,
            isTrimmed=True):
        for e in self.exposures:             
            if e.getDetector().getId() == amp.getId():
                if self.isTrimmed:
                    img = imageFactory(e.getMaskedImage().getVariance(),
                            amp.getDiskDataSec(), afwImage.PARENT)
                else:
                    img = imageFactory(e.getMaskedImage().getVariance(),
                            amp.getDiskAllPixels(), afwImage.PARENT)
                return amp.prepareAmpData(img)
        return None

def getFixedWcs(exposures, amp):
    if len(exposures) > 1:
        for exp in exposures:
            if exp.getDetector().getId() == amp.getId():
                if exp.hasWcs():
                    wcs = exp.getWcs()
                    amp.prepareWcsData(wcs)
                else:
                    wcs = None
    else:
        if exposures[0].hasWcs():
            wcs = exposures[0].getWcs()
            amp.prepareWcsData(wcs)
        else:
            wcs = None
    origin = amp.getDataSec()
    #shift the reference pixel to account for the location of the amp in the 
    #ccd.
    wcs.shiftReferencePixel(origin.getMinX(), origin.getMinY())
    return wcs

def assembleCcd(exposures, ccd, reNorm=True, isTrimmed=True, keysToRemove=[], imageFactory=afwImage.ImageF):
    display = lsstDebug.Info(__name__).display 
    ccd.setTrimmed(isTrimmed)
    wcs = getFixedWcs(exposures, cameraGeom.cast_Amp(ccd[15]))
    filter = exposures[0].getFilter()
    metadata = exposures[0].getMetadata()
    calib = exposures[0].getCalib()
    for k in keysToRemove:
        if metadata.exists(k):
            metadata.remove(k)
    #detector = cameraGeom.cast_Ccd(exposures[0].getDetector().getParent())
    #dl = detector.getDefects()
    gain = 0
    namps = 0
    #for a in detector:
    for a in ccd:
        gain += cameraGeom.cast_Amp(a).getElectronicParams().getGain()
        namps += 1.
    gain /= namps
    if len(exposures) > 1:
        lif = listImageFactory(exposures)
        lmf = listMaskFactory(exposures)
        lvf = listVarianceFactory(exposures)
    else:
        lif = singleImageFactory(exposures[0])
        lmf = singleMaskFactory(exposures[0])
        lvf = singleVarianceFactory(exposures[0])
    ccdImage = cameraGeomUtils.makeImageFromCcd(ccd, imageSource = lif,
            imageFactory = imageFactory, bin=False)
    ccdVariance = cameraGeomUtils.makeImageFromCcd(ccd, imageSource = lvf,
            imageFactory = afwImage.ImageF, bin=False)
    ccdMask = cameraGeomUtils.makeImageFromCcd(ccd, imageSource = lmf,
            imageFactory = afwImage.MaskU, bin=False)
    mi = afwImage.makeMaskedImage(ccdImage,
        ccdMask, ccdVariance)
    if reNorm:
        mi *= gain
        metadata.set("GAIN", 1.0)
    ccdExposure = afwImage.makeExposure(mi)
    if wcs is not None:
        ccdExposure.setWcs(wcs)
    ccdExposure.setMetadata(metadata)
    ccdExposure.setFilter(filter)
    #ccdExposure.setDetector(detector)
    ccdExposure.setDetector(ccd)
    ccdExposure.getCalib().setExptime(calib.getExptime())
    ccdExposure.getCalib().setMidTime(calib.getMidTime())
    if not ccdVariance.getArray().max() == 0:
        (medgain, meangain) = isr.calcEffectiveGain(ccdExposure)
        metadata.add("MEDGAIN", medgain)
        metadata.add("MEANGAIN", meangain)
        metadata.add("GAINEFF", medgain)
    else:
        metadata.add("MEDGAIN", 0.)
        metadata.add("MEANGAIN", 0.)
        metadata.add("GAINEFF", 0.)

    if display:
        ds9.mtv(ccdExposure)
    
    return ccdExposure
