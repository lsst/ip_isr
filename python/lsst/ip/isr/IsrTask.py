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

import math, re
import lsst.afw.detection   as afwDetection
import lsst.afw.image       as afwImage
import lsst.afw.geom        as afwGeom
import lsst.afw.cameraGeom  as cameraGeom
import lsst.afw.math        as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.meas.algorithms as algorithms
import lsst.pex.logging     as pexLog
import lsst.pex.exceptions  as pexExcept
from copy import deep_copy

# relative imports
import isrLib

import lsst.pex.policy as pexPolicy
import lsst.pipe.base as pipeBase
class IsrTask(pipeBase.Task):
    def __init__(self, **keyArgs):
        pipeBase.Task.__init__(self, **keyArgs)
        self.isr = Isr()
        self.methodList = []
        self.lookupTable = lookupTable
        for methodname in policy.getArray("methodList"):
            self.methodList.append(getattr(self, methodname))

    def run(self, exposure, calibSet):
        #The ISR routines operate in place.  A copy of the original exposure
        #will be made and the reduced exposure will be returned.
        workingExposure = deep_copy(exposure)
        for m in self.methodList:
            workingExposure = m(workingExposure, calibSet)
        return workingExposure #I guess this should be a struct instead...

    @staticmethod    
    def getPolicy():
        policy = pexPolicy.Policy()
        policy.set("Fwhm", 4.) #In pixels
        policy.set("GrowSatruationFootprints", True)
        policy.set("SatruatedMaskName", "SAT")
        policy.set("FlatScalingType", "USER")
        policy.set("scalingValue", 1.)
        policy.set("methodList", ["doSaturationDetection", "doOverscanCorrection", "doBiasSubtraction", \
                   "doDarkCorrection", "doFlatCorrection"])
        return policy
    
    def doCrosstalkCorrection(self, exposure, calibSet):
        pass

    def doSaturationDetection(self, exposure, calibSet):
        amp = cameraGeom.cast_Amp(exposure.getDetector())
        ep = amp.getElectronicParams()
        satvalue = ep.getSaturationLevel()
        fwhm = self.policy.get("Fwhm")
        grow = self.policy.get("GrowSatruationFootprints")
        maskname = self.policy.get("SaturatedMaskName")
        self.isr.saturationDetection(exposure, satvalue, doMask=True, maskName=maskname)
        return exposure

    def doLinearization(self, exposure, calibSet):
        self.isr.linearization(exposure, self.lookupTable)
        return exposure

    def doOverscanCorrection(self, exposure, calibSet):
        amp = cameraGeom.cast_Amp(exposure.getDetector())
        overscanBbox = amp.getDiskBiasSec()
        self.isr.overscanCorrection(exposure, overscanBBox, fittype=fittype, polyorder, imageFactory=afwImage.ImageF)
        return exposure

    def doBiasSubtraction(self, exposure, calibSet):
        bias = calibSet['bias']
        self.isr.biasCorrection(exposure, bias)
        return exposure 

    def doDarkCorrection(self, exposure, calibSet):
        dark = calibSet['dark']
        darkscaling = darkexposure.getCalib().getExptime()
        expscaling = exposure.getCalib().getExptime()
        self.isr.darkCorrection(exposure, dark, expscaling, darkscaling)
        return exposure

    def doFringeCorrection(self, exposure, calibSet):
        pass

    def doFlatCorrection(self, exposure, calibSet):
        flat = calibSet['flat']
        scalingtype = self.policy.get("flatScalingType")
        scalingvalue = self.policy.get("scalingValue")
        self.isr.flatCorrection(exposure, flat, scalingtype, scaling = scalingvalue)   
        return exposure

    def doIlluminationCorrection(self, exposure, calibSet):
        pass
