import lsst.pex.policy as pexPolicy
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display.ds9 as ds9
import Isr
import os,sys,eups

class listImageFactory(cameraGeomUtils.GetCcdImage):
    def __init__(self, exposures):
        self.exposures = exposures
        self.isRaw = True
    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.ImageF):
        for e in self.exposures:             
            if e.getDetector().getId() == amp.getId():
              img = imageFactory(e.getMaskedImage().getImage(),
                      amp.getDiskDataSec())
              return img
        return None

class listMaskFactory(cameraGeomUtils.GetCcdImage):
    def __init__(self, exposures):
        self.exposures = exposures
        self.isRaw = True
    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.ImageF):
        for e in self.exposures:             
            if e.getDetector().getId() == amp.getId():
              img = imageFactory(e.getMaskedImage().getMask(),
                      amp.getDiskDataSec())
              return img
        return None

class listVarianceFactory(cameraGeomUtils.GetCcdImage):
    def __init__(self, exposures):
        self.exposures = exposures
        self.isRaw = True
    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.ImageF):
        for e in self.exposures:             
            if e.getDetector().getId() == amp.getId():
              img = imageFactory(e.getMaskedImage().getVariance(),
                      amp.getDiskDataSec())
              return img
        return None

def assembleCcd(exposures, ccd, isTrimmed = True, isOnDisk = True, keysToRemove = []):
    wcs = exposures[0].getWcs()
    filter = exposures[0].getFilter()
    metadata = exposures[0].getMetadata()
    for k in keysToRemove:
        if metadata.exists(k):
            metadata.remove(k)
    detector = cameraGeom.cast_Ccd(exposures[0].getDetector().getParent())
    dl = detector.getDefects()
    gain = 0
    namps = 0
    for a in detector:
        gain += cameraGeom.cast_Amp(a).getElectronicParams().getGain()
        namps += 1.
    gain /= namps
    lif = listImageFactory(exposures)
    lmf = listMaskFactory(exposures)
    lvf = listVarianceFactory(exposures)
    ccdImage = cameraGeomUtils.makeImageFromCcd(ccd, imageSource = lif,
            isTrimmed = isTrimmed, imageFactory = afwImage.ImageF, bin=False)
    ccdVariance = cameraGeomUtils.makeImageFromCcd(ccd, imageSource = lvf,
            isTrimmed = isTrimmed, imageFactory = afwImage.ImageF, bin=False)
    ccdMask = cameraGeomUtils.makeImageFromCcd(ccd, imageSource = lmf,
            isTrimmed = isTrimmed, imageFactory = afwImage.MaskU, bin=False)
    mi = afwImage.makeMaskedImage(ccdImage,
        ccdMask, ccdVariance)
    mi *= gain
    metadata.set("GAIN", 1.0)
    ccdExposure = afwImage.makeExposure(mi, wcs)
    ccdExposure.setWcs(wcs)
    ccdExposure.setMetadata(metadata)
    ccdExposure.setFilter(filter)
    ccdExposure.setDetector(detector)
    (medgain, meangain) = Isr.calcEffectiveGain(ccdExposure)
    metadata.add("MEDGAIN", medgain)
    metadata.add("MEANGAIN", meangain)
    metadata.add("GAINEFF", medgain)
    
    return ccdExposure
