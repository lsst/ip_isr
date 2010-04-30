import lsst.pex.policy as pexPolicy
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display.ds9 as ds9
import os,sys,eups

class listImageFactory(cameraGeomUtils.GetCcdImage):
    def __init__(self, exposures):
        self.exposures = exposures
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
    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.ImageF):
        for e in self.exposures:             
            if e.getDetector().getId() == amp.getId():
              img = imageFactory(e.getMaskedImage().getVariance(),
                      amp.getDiskDataSec())
              return img
        return None

def assembleCcd(exposures, ccd, isTrimmed = True, isOnDisk = True):
    wcs = exposures[0].getWcs()
    filter = exposures[0].getFilter()
    metadata = exposures[0].getMetadata()
    metadata.remove("BIASSEC")
    metadata.remove("DATASEC")
    detector = cameraGeom.cast_Ccd(exposures[0].getDetector().getParent())

    lif = listImageFactory(exposures)
    lmf = listMaskFactory(exposures)
    lvf = listVarianceFactory(exposures)
    ccdImage = cameraGeomUtils.makeImageFromCcd(ccd, imageSource = lif,
            isTrimmed = isTrimmed, imageFactory = afwImage.ImageF)
    ccdVariance = cameraGeomUtils.makeImageFromCcd(ccd, imageSource = lvf,
            isTrimmed = isTrimmed, imageFactory = afwImage.ImageF)
    #ccdMask = afwImage.MaskU(ccdImage.getWidth(), ccdImage.getHeight())
    #ccdMask.set(0)
    ccdMask = cameraGeomUtils.makeImageFromCcd(ccd, imageSource = lmf,
            isTrimmed = isTrimmed, imageFactory = afwImage.MaskU)
    ccdExposure = afwImage.makeExposure(afwImage.makeMaskedImage(ccdImage,
        ccdMask, ccdVariance), wcs)
    ccdExposure.setWcs(wcs)
    ccdExposure.setMetadata(metadata)
    ccdExposure.setFilter(filter)
    ccdExposure.setDetector(detector)
    
    return ccdExposure
