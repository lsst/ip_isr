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
from .ccdAssembler import CcdAssembler


class CcdListAssembler(CcdAssembler):
    def __init__(self, exposureList, ccd, reNorm=True, setGain=True, keysToRemove=[], isTrimmed=True, display=False):
        from .imageFactories import ListImageFactory
        from .imageFactories import ListMaskFactory
        from .imageFactories import ListVarianceFactory
        #set the reference exposure from the list.  Assume the first is fine...
        self.exposure = exposureList[0]
        #set the reference amp for doing the WCS manipulation
        amp = cameraGeom.cast_Amp(self.exposure.getDetector())
        if amp is None:
            raise("First exposure does not have an amp level detector")
        else:
            self.amp = amp
        self.detector = ccd
        self.detector.setTrimmed(isTrimmed)
        self.reNorm = reNorm
        self.ktr = keysToRemove
        self.outputImageFactory = self.exposure.getMaskedImage().getImage().Factory
        self.ifactory = ListImageFactory(exposureList)
        self.mfactory = ListMaskFactory(exposureList)
        self.vfactory = ListVarianceFactory(exposureList)
        self.filter = self.exposure.getFilter()
        self.metadata = self.exposure.getMetadata()
        self.calib = self.exposure.getCalib()
        self.display = display

