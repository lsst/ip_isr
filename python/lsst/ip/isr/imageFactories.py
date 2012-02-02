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
import lsst.afw.cameraGeom.utils as cameraGeomUtils
class SingleImageFactory(cameraGeomUtils.GetCcdImage):
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

class SingleMaskFactory(cameraGeomUtils.GetCcdImage):
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

class SingleVarianceFactory(cameraGeomUtils.GetCcdImage):
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

class ListImageFactory(cameraGeomUtils.GetCcdImage):
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

class ListMaskFactory(cameraGeomUtils.GetCcdImage):
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

class ListVarianceFactory(cameraGeomUtils.GetCcdImage):
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
