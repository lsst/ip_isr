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
import lsst.meas.algorithms as measAlg
import lsst.afw.cameraGeom  as cameraGeom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
#from .isr import Linearization
from .isr import Isr
from .ccdAssembler import CcdAssembler
from . import isrLib

class IsrTaskConfig(pexConfig.Config):
    doWrite = pexConfig.Field(dtype=bool, doc="Write output?", default=True)
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
    overscanFitType = pexConfig.ChoiceField(
        dtype = str,
        doc = "The method for fitting the overscan bias level.",
        default = 'MEDIAN',
        allowed = {"POLY": "Fit polynomial to the longest axis of the overscan region",
            "MEAN": "Correct using the mean of the overscan region",
            "MEDIAN": "Correct using the median of the overscan region",
        },
    )
    overscanPolyOrder = pexConfig.Field(
        dtype = int,
        doc = "Order of polynomial to fit if overscan fit type is POLY",
        default = 1,
    )
    growSaturationFootprintSize = pexConfig.Field(
        dtype = int,
        doc = "Number of pixels by which to grow the saturation footprints",
        default = 1,
    )
    growDefectFootprintSize = pexConfig.Field(
        dtype = int,
        doc = "Number of pixels by which to grow the defect (bad and nan) footprints",
        default = 1,
    )
    fluxMag0T1 = pexConfig.Field(
        dtype = float,
        doc = "The approximate flux of a zero-magnitude object in a one-second exposure",
        default = 1e10,
    )
    setGainAssembledCcd = pexConfig.Field(
        dtype = bool,
        doc = "update exposure metadata in the assembled ccd to reflect the effective gain of the assembled chip",
        default = True,
    )
    keysToRemoveFromAssembledCcd = pexConfig.ListField(
        dtype = str,
        doc = "fields to remove from the metadata of the assembled ccd.",
        default = [],
    )
    reNormAssembledCcd = pexConfig.Field(
        dtype = bool,
        doc = "renormalize the assembled chips to have unity gain.  False if setGain is False",
        default = True,
    )
    methodList = pexConfig.ListField(
        dtype = str,   
        doc = "The list of ISR corrections to apply in the order they should be applied",
        default = ["doConversionForIsr", "doSaturationDetection", "doOverscanCorrection", "doBiasSubtraction", "doVariance", "doDarkCorrection", "doFlatCorrection"],
    )
    
class IsrTask(pipeBase.Task):
    ConfigClass = IsrTaskConfig
    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.isr = Isr()
        self.methodList = []
        for methodname in self.config.methodList:
            self.methodList.append(getattr(self, methodname))

    def run(self, exposure, calibSet):
        """Do instrument signature removal on an exposure: saturation, bias, overscan, dark, flat, fringe correction

        @param exposure Apply ISR to this Exposure
        @param calibSet Dictionary of calibration products (bias/zero, dark, flat, fringe, linearization information)
        @return a pipeBase.Struct with fields:
        - postIsrExposure: the exposure after application of ISR
        """

        #The ISR routines operate in place.  A copy of the original exposure
        #will be made and the reduced exposure will be returned.
        workingExposure = exposure.Factory(exposure, True)
        for m in self.methodList:
            workingExposure = m(workingExposure, calibSet)
        return pipeBase.Struct(postIsrExposure=workingExposure)

    def runButler(self, butler, dataid):
        """Run the ISR given a butler
        @param butler Butler describing the data repository
        @param dataid A data identifier of the amp to process
        @return a pieBase.Struct see self.run for returned fields
        """
        calibSet = self.makeCalibDict(butler, dataid)
        output = self.run(butler.get("raw", dataid), calibSet)
        if self.config.doWrite:
            butler.put(output.postIsrExposure, "postISR", dataId=dataid)
        return output
        
    def makeCalibDict(self, butler, dataId):
        ret = {}
        required = {"doBiasSubtraction": "bias",
                    "doDarkCorrection": "dark",
                    "doFlatCorrection": "flat",
                    }
        for method in required.keys():
            if method in self.config.methodList:
                calib = required[method]
                ret[calib] = butler.get(calib, dataId)
        return ret

    def doConversionForIsr(self, exposure, calibSet):
        if not isinstance(exposure, afwImage.ExposureU):
            raise Exception("ipIsr.convertImageForIsr: Expecting Uint16 image. Got\
                %s."%(exposure.__repr__()))

        newexposure = exposure.convertF()
        mi = newexposure.getMaskedImage()
        var = afwImage.ImageF(mi.getBBox(afwImage.PARENT))
        mask = afwImage.MaskU(mi.getBBox(afwImage.PARENT))
        mask.set(0)
        newexposure.setMaskedImage(afwImage.MaskedImageF(mi.getImage(), mask, var))
        return newexposure

    def doVariance(self, exposure, calibSet):
        for amp in self._getAmplifiers(exposure):
            exp = exposure.Factory(exposure, amp.getDiskDataSec())
            self.isr.updateVariance(exp.getMaskedImage(), amp.getElectronicParams().getGain())
        return exposure

    def doCrosstalkCorrection(self, exposure, calibSet):
        pass

    def doSaturationCorrection(self, exposure, calibSet):
        fwhm = self.config.fwhm
        grow = self.config.growSaturationFootprintSize
        maskname = self.config.saturatedMaskName
        for amp in self._getAmplifiers(exposure):
            ep = amp.getElectronicParams()
            satvalue = ep.getSaturationLevel()
            exp = exposure.Factory(exposure, amp.getDiskDataSec())
            self.isr.saturationCorrection(exp.getMaskedImage(), satvalue, fwhm, growFootprints=grow, maskName=maskname)
        return exposure

    def doSaturationDetection(self, exposure, calibSet):
        for amp in self._getAmplifiers(exposure):
            datasec = amp.getDiskDataSec()
            exp = exposure.Factory(exposure, datasec)
            ep = amp.getElectronicParams()
            satvalue = ep.getSaturationLevel()
            maskname = self.config.saturatedMaskName
            self.isr.makeThresholdMask(exp.getMaskedImage(), satvalue, growFootprints=0, maskName=maskname)
        return exposure

    def doSaturationInterpolation(self, exposure, calibSet):
        #Don't loop over amps since saturation can cross amp boundaries
        maskname = self.config.saturatedMaskName
        fwhm = self.config.fwhm
        grow = self.config.growSaturationFootprintSize
        self.isr.interpolateFromMask(exposure.getMaskedImage(), fwhm, growFootprints=grow, maskName=maskname)
        return exposure
    
    def doMaskAndInterpDefect(self, exposure, calibSet):
        #Don't loop over amps since defects could cross amp boundaries
        fwhm = self.config.fwhm
        grow = self.config.growDefectFootprintSize
        defectBaseList = cameraGeom.cast_Ccd(exposure.getDetector()).getDefects()
        defectList = measAlg.DefectListT()
        #mask bad pixels in the camera class
        #create master list of defects and add those from the camera class
        for d in defectBaseList:
            bbox = d.getBBox()
            nd = measAlg.Defect(bbox)
            defectList.append(nd)
        self.isr.maskPixelsFromDefectList(exposure.getMaskedImage(), defectList, maskName='BAD')
        defectList = self.isr.getDefectListFromMask(exposure.getMaskedImage(), maskName='BAD', growFootprints=grow)
        self.isr.interpolateDefectList(exposure.getMaskedImage(), defectList, fwhm)
        return exposure

    def doMaskAndInterpNan(self, exposure, calibSet):
        #Don't loop over amps since nans could cross amp boundaries
        fwhm = self.config.fwhm
        grow = self.config.growDefectFootprintSize
        #find unmasked bad pixels and mask them
        exposure.getMaskedImage().getMask().addMaskPlane("UNMASKEDNAN") 
        unc = isrLib.UnmaskedNanCounterF()
        unc.apply(exposure.getMaskedImage())
        nnans = unc.getNpix()
        expmeta = exposure.getMetadata()
        expmeta.set("NUMNANS", nnans)
        if not nnans == 0:
		raise RuntimeError("There were %i unmasked NaNs"%(nnans))
        #get footprints of bad pixels not in the camera class
        undefects = self.isr.getDefectListFromMask(exposure.getMaskedImage(), maskName='UNMASKEDNAN', growFootprints=grow)
        #interpolate all bad pixels
        self.isr.interpolateDefectList(exposure.getMaskedImage(), undefects, fwhm)
        return exposure

    '''
    def doLinearization(self, exposure, calibSet):
        linearizer = Linearization(calibSet['linearityFile'])
        linearizer.apply(exposure)
        return exposure
    '''
    
    def doOverscanCorrection(self, exposure, calibSet):
        fittype = self.config.overscanFitType
        polyorder = self.config.overscanPolyOrder
        for amp in self._getAmplifiers(exposure):
            expImage = exposure.getMaskedImage().getImage()
            overscan = expImage.Factory(expImage, amp.getDiskBiasSec())
            exp = exposure.Factory(exposure, amp.getDiskDataSec())
            self.isr.overscanCorrection(exp.getMaskedImage(), overscan, fittype=fittype, polyorder=polyorder,
                                        imageFactory=afwImage.ImageF)
        return exposure

    def doBiasSubtraction(self, exposure, calibSet):
        biasExposure = calibSet['bias']
        for amp in self._getAmplifiers(exposure):
            exp, bias = self._getCalibration(exposure, biasExposure, amp)
            self.isr.biasCorrection(exp.getMaskedImage(), bias.getMaskedImage())
        
        return exposure 

    def doDarkCorrection(self, exposure, calibSet):
        darkexposure = calibSet['dark']
        darkscaling = darkexposure.getCalib().getExptime()
        expscaling = exposure.getCalib().getExptime()
        
        for amp in self._getAmplifiers(exposure):
            exp, dark = self._getCalibration(exposure, darkexposure, amp)
            self.isr.darkCorrection(exp.getMaskedImage(), dark.getMaskedImage(), expscaling, darkscaling)
        return exposure

    def doFringeCorrection(self, exposure, calibSet):
        pass

    def _getAmplifiers(self, exposure):
        """Return list of all amplifiers in an Exposure"""
        amp = cameraGeom.cast_Amp(exposure.getDetector())
        if amp is not None:
            return [amp]
        ccd = cameraGeom.cast_Ccd(exposure.getDetector())
        assert ccd is not None
        return [cameraGeom.cast_Amp(a) for a in ccd]

    def _getCalibration(self, exposure, calibration, amp):
        """Get a suitably-sized calibration exposure"""
        exp = exposure
        calib = calibration
        if exp.getDimensions() != amp.getDataSec().getDimensions():
            # Just the amp of interest
            if exp.getDimensions() == amp.getParent().getAllPixels().getDimensions():
                # We have a full CCD; just cut out the amp
                exp = exp.Factory(exp, amp.getDataSec())
            else:
                # Probably dealing with a disk image
                exp = exp.Factory(exp, amp.getDiskDataSec())
        if exp.getDimensions() == calib.getDimensions():
            return exp, calib
        # Try just the calibration's pixels of interest
        try:
            calib = calib.Factory(calib, amp.getDataSec(True)) # Calib is likely trimmed and assembled
        except:
            pass
        if exp.getDimensions() == calib.getDimensions():
            return exp, calib
        raise RuntimeError("Dimensions for exposure (%s) and calibration (%s) don't match" % \
                           (exposure.getDimensions(), calibration.getDimensions()))

    def doFlatCorrection(self, exposure, calibSet):
        flatfield = calibSet['flat']
        scalingtype = self.config.flatScalingType
        scalingvalue = self.config.flatScalingValue

        for amp in self._getAmplifiers(exposure):
            exp, flat = self._getCalibration(exposure, flatfield, amp)
            self.isr.flatCorrection(exp.getMaskedImage(), flat.getMaskedImage(), scalingtype, scaling = scalingvalue)   
        return exposure

    def doIlluminationCorrection(self, exposure, calibSet):
        pass

    def doCcdAssembly(self, exposureList):
        renorm = self.config.reNormAssembledCcd
        setgain = self.config.setGainAssembledCcd
        k2rm = self.config.keysToRemoveFromAssembledCcd
        fluxMag0T1 = self.config.fluxMag0T1
        assembler = CcdAssembler(exposureList, reNorm=renorm, setGain=setgain, keysToRemove=k2rm,
                                 fluxMag0T1=fluxMag0T1)
        return assembler.assembleCcd()
