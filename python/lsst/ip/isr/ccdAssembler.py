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

import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display.ds9 as ds9
from .isr import Isr


class CcdAssembler(object):
    def __init__(self, exposure, ccd, reNorm=True, setGain=True, keysToRemove=[], isTrimmed=True, display=False):
        from .imageFactories import SingleImageFactory
        from .imageFactories import SingleMaskFactory
        from .imageFactories import SingleVarianceFactory
        #set the reference exposure from the list.  Assume the first is fine...
        self.exposure = exposure
        #set the reference amp for doing the WCS manipulation
        amp = cameraGeom.cast_Amp(ccd[0])
        if amp is None:
            raise("ccd does not contain amp level detectors")
        else:
            self.amp = amp
        self.detector = ccd
        self.detector.setTrimmed(isTrimmed)
        self.reNorm = reNorm
        self.ktr = keysToRemove
        self.outputImageFactory = self.exposure.getMaskedImage().getImage().Factory
        self.ifactory = SingleImageFactory(self.exposure)
        self.mfactory = SingleMaskFactory(self.exposure)
        self.vfactory = SingleVarianceFactory(self.exposure)
        self.filter = self.exposure.getFilter()
        self.metadata = self.exposure.getMetadata()
        self.calib = self.exposure.getCalib()
        self.display

    def getFixedWcs(self):
        if self.exposure.hasWcs():
            wcs = self.exposure.getWcs()
            self.amp.prepareWcsData(wcs)
            #shift the reference pixel to account for the location of the amp in the ccd.
            origin = self.amp.getDataSec()
            wcs.shiftReferencePixel(origin.getMinX(), origin.getMinY())
        else:
            wcs = None
        return wcs

    def setGain(self, ccdExposure):
        gain = 0
        namps = 0
        for a in self.detector:
            gain += cameraGeom.cast_Amp(a).getElectronicParams().getGain()
            namps += 1.
        gain /= namps
        if self.reNorm:
            mi *= gain
            self.metadata.set("GAIN", 1.0)
        isr = Isr()
        (medgain, meangain) = isr.calcEffectiveGain(ccdExposure)
        self.metadata.add("MEDGAIN", medgain)
        self.metadata.add("MEANGAIN", meangain)
        self.metadata.add("GAINEFF", medgain)
        
    def setExposureComponents(self, ccdExposure):
        wcs = self.getFixedWcs()
        if wcs is not None:
            ccdExposure.setWcs(wcs)
        for k in self.ktr:
            if metadata.exists(k):
                metadata.remove(k)
        ccdExposure.setMetadata(self.metadata)
        ccdExposure.setFilter(self.filter)
        ccdExposure.setDetector(self.detector)
        ccdExposure.getCalib().setExptime(self.calib.getExptime())
        ccdExposure.getCalib().setMidTime(self.calib.getMidTime())


    def assembleCcd(self):
        ccdImage = cameraGeomUtils.makeImageFromCcd(self.detector, imageSource = self.ifactory,
            imageFactory = self.outputImageFactory, bin=False)
        ccdVariance = cameraGeomUtils.makeImageFromCcd(self.detector, imageSource = self.vfactory,
            imageFactory = afwImage.ImageF, bin=False)
        ccdMask = cameraGeomUtils.makeImageFromCcd(self.detector, imageSource = self.mfactory,
            imageFactory = afwImage.MaskU, bin=False)
        mi = afwImage.makeMaskedImage(ccdImage, ccdMask, ccdVariance)
        ccdExposure = afwImage.makeExposure(mi)
        if self.setGain:
            if not ccdVariance.getArray().max() == 0:
                self.setGain(ccdExposure)
            else:
                raise("Can't calculate the effective gain since the variance plane is set to zero")
        self.setExposureComponents(ccdExposure)

        if self.display:
            ds9.mtv(ccdExposure)
    
        return ccdExposure
