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

import lsst.afw.image       as afwImage
import lsst.afw.cameraGeom  as cameraGeom
from copy import deep_copy
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .isr import Linearization
from .isr import Isr

import isrLib
class IsrConfig(pexConfig.Config):
    fwhm = pexConfig.Field(
        dtype = float,
        doc = "FWHM of PSF (arcsec)",
        default = 1.0,
    )
    #This is needed for both the detection and correction aspects
    saturatedMaskName = pexConfig.Field(
        dtype = str,
        doc = "Name of mask plane to use in saturation detection",
        default = "SAT",
    )
    flatScalingType = pexConfig.ChoiceField(
        dtype = str,
        doc = "The method for scaling the flat on the fly.",
        default = 'USER',
        allowed = {"USER": "User defined scaling",
            "MEAN": "Scale by the inverse of the mean",
            "MEDIAN": "Scale by the inverse of the median",
        },
    )
    flatScalingValue = pexConfig.Field(
        dtype = float,
        doc = "If scaling type is USER, a value for the scaling must be provided",
        default = 1.0,
    )
    nGrowSaturated = pexConfig.Field(
        dtype = int,
        doc = "Number of pixels by which to grow the saturated footprints",
        default = 1,
    )
    growSaturationFootprints = pexConfig.Field(
        dtype = bool,
        doc = "Should the saturated footprints be grown?",
        default = True,
    )
    methodList = pexConfig.ListField(
        dtype = float,
        doc = "The list of ISR corrections to apply in the order they should be applied",
        default = ["doSaturationDetection", "doLinearization", "doOverscanCorrection", "doBiasSubtraction", \
                   "doDarkCorrection", "doFlatCorrection"],
    )
    
class IsrTask(pipeBase.Task):
    ConfigClass = IsrConfig
    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.isr = Isr()
        self.methodList = []
        for methodname in config.methodList:
            self.methodList.append(getattr(self, methodname))

    def run(self, exposure, calibSet):
        """Do instrument signature removal on an exposure: saturation, bias, overscan, dark, flat, fringe correction

        @param exposure Apply ISR to this Exposure
        @param calibSet Dictionary of calibration products (bias/zero, dark, flat, fringe, linearization information)
        @return a pipeBase.Struct with fields:
        - postIsrExposure: the exposure after application of ISR
        - metadata: metadata about the ISR process
        """

        #The ISR routines operate in place.  A copy of the original exposure
        #will be made and the reduced exposure will be returned.
        workingExposure = deep_copy(exposure)
        for m in self.methodList:
            workingExposure = m(workingExposure, calibSet)
        return workingExposure #I guess this should be a struct instead...

    def doCrosstalkCorrection(self, exposure, calibSet):
        pass

    def doSaturationCorrection(self, exposure, calibSet):
        amp = cameraGeom.cast_Amp(exposure.getDetector())
        ep = amp.getElectronicParams()
        satvalue = ep.getSaturationLevel()
        fwhm = self.config.fwhm
        grow = self.config.growSatruationFootprints
        maskname = self.config.saturatedMaskName
        if grow:
            grow = self.config.nGrowSaturated
        else:
            grow = 0
        self.isr.saturationCorrection(exposure, satvalue, fwhm, growFootprints=grow, maskName=maskname)
        return exposure

    def doSaturationDetection(self, exposure, calibSet):
        amp = cameraGeom.cast_Amp(exposure.getDetector())
        ep = amp.getElectronicParams()
        satvalue = ep.getSaturationLevel()
        maskname = self.config.saturatedMaskName
        self.isr.saturationDetection(exposure, satvalue, doMask=True, maskName=maskname)
        return exposure

    def doSaturationInterpolation(self, exposure, calibSet):
        maskname = self.config.saturatedMaskName
        fwhm = self.config.fwhm
        grow = self.config.growSatruationFootprints
        if grow:
            grow = self.config.nGrowSaturated
        else:
            grow = 0
        self.isr.saturationInterpolation(exposure, fwhm, growFootprints=grow, maskName=maskname)
        return exposure

    def doLinearization(self, exposure, calibSet):
        linearizer = Linearization(calibSet['linearityFile'])
        linearizer.apply(exposure)
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
        scalingtype = self.config.flatScalingType
        scalingvalue = self.config.flatScalingValue
        self.isr.flatCorrection(exposure, flat, scalingtype, scaling = scalingvalue)   
        return exposure

    def doIlluminationCorrection(self, exposure, calibSet):
        pass
