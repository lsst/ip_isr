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

import lsst.pex.policy as pexPolicy
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display.ds9 as ds9
import isr
import os,sys,eups,math

class listImageFactory(cameraGeomUtils.GetCcdImage):
    def __init__(self, exposures):
        self.exposures = exposures
        self.isRaw = True
    def getImage(self, ccd, amp, expType=None, imageFactory=afwImage.ImageF):
        for e in self.exposures:             
            if e.getDetector().getId() == amp.getId():
              img = imageFactory(e.getMaskedImage().getImage(),
                      amp.getDiskDataSec(), afwImage.PARENT)
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
                      amp.getDiskDataSec(), afwImage.PARENT)
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
                      amp.getDiskDataSec(), afwImage.PARENT)
              return img
        return None

def assembleCcd(exposures, ccd, reNorm = True, isOnDisk = True, keysToRemove = []):
    if exposures[0].hasWcs():
        wcs = exposures[0].getWcs()
    else:
        wcs = False
    filter = exposures[0].getFilter()
    metadata = exposures[0].getMetadata()
    calib = exposures[0].getCalib()
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
            isTrimmed = True, imageFactory = afwImage.ImageF, bin=False)
    ccdVariance = cameraGeomUtils.makeImageFromCcd(ccd, imageSource = lvf,
            isTrimmed = True, imageFactory = afwImage.ImageF, bin=False)
    ccdMask = cameraGeomUtils.makeImageFromCcd(ccd, imageSource = lmf,
            isTrimmed = True, imageFactory = afwImage.MaskU, bin=False)
    mi = afwImage.makeMaskedImage(ccdImage,
        ccdMask, ccdVariance)
    if reNorm:
        mi *= gain
        metadata.set("GAIN", 1.0)
    ccdExposure = afwImage.makeExposure(mi)
    if wcs is not False:
        print "WCS: ",wcs
        ccdExposure.setWcs(wcs)
    ccdExposure.setMetadata(metadata)
    ccdExposure.setFilter(filter)
    ccdExposure.setDetector(detector)
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
    
    return ccdExposure
