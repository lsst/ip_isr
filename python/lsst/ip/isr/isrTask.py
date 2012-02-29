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
    calibsAreTrimmed = pexConfig.Field(
        dtype = bool,
        doc = "Have calibs been trimmed of non-science pixels?",
        default = True,
    )
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
    doConversionForIsr = pexConfig.ChoiceField(
        dtype = str,
        doc = "Do the conversion from int to float for ISR",
        allowed = {"pre":"Do conversion before assembly", "None":"Don't do conversion"},
        default = "pre",
    )
    doSaturationDetection = pexConfig.ChoiceField(
        dtype = str,
        doc = "Detect saturated pixels and mask them",
        allowed = {"pre":"Do saturation correction before assembly", "None":"Don't do detection"},
        default = "pre",
    )
    doLinearization = pexConfig.ChoiceField(
        dtype = str,
        doc = "Do linearization of sensor response",
        allowed = {"pre":"Do linearization before assembly", "None":"Don't do conversion"},
        default = "None",
    )
    doOverscanCorrection = pexConfig.ChoiceField(
        dtype = str,
        doc = "Correct image for overscan",
        allowed = {"pre":"Do overscan correction before assembly.", "None":"Don't do correction"},
        default = "pre",
    )
    doTrimExposure = pexConfig.ChoiceField(
        dtype = str,
        doc = "Trim exposure to just imaging pixels.",
        allowed = {"pre":"Do trimming before assembly.", "None":"Don't do trimming."},
        default = "None",
    )
    doBiasSubtraction = pexConfig.ChoiceField(
        dtype = str,
        doc = "Do bias correction",
        allowed = {"pre":"Do bias correction before assembly", "post":"Do bias correction after assembly", "None":"Don't do correction"},
        default = "pre",
    )
    doVariance = pexConfig.ChoiceField(
        dtype = str,
        doc = "Update variance plane",
        allowed = {"pre":"Update variance plane before assembly", "post":"Update variance plane after assembly", "None":"Don't do update"},
        default = "pre",
    )
    doDarkCorrection = pexConfig.ChoiceField(
        dtype = str,
        doc = "Correct for sensor dark current.",
        allowed = {"pre":"Do dark current before assembly", "post":"Do dark current correction after assembly", "None":"Don't do correction"},
        default = "pre",
    )
    doFringeCorrection = pexConfig.ChoiceField(
        dtype = str,
        doc = "Do correction for sky line fringing",
        allowed = {"pre":"Do fringe correction before assembly", "post":"Do fringe correction after assembly", "None":"Don't do correction"},
        default = "None",
    )
    doFlatCorrection = pexConfig.ChoiceField(
        dtype = str,
        doc = "Do the flat field correction",
        allowed = {"pre":"Do flat correction before assembly", "post":"Do flat correction after assembly", "None":"Don't do correction"},
        default = "pre",
    )
    doAssembly = pexConfig.Field(
        dtype = bool,
        doc = "Do assembly of channels into a single sensor",
        default = True,
    )
    doSaturationInterpolation = pexConfig.ChoiceField(
        dtype = str,
        doc = "Interpolate pixels masked as saturated",
        allowed = {"pre":"Do flat correction before assembly", "post":"Do the interpolation after assembly", "None":"Don't interpolate saturated pixels."},
        default = "post",
    )
    doMaskAndInterpDefect = pexConfig.ChoiceField(
        dtype = str,
        doc = "Mask and interpolate over known defects.",
        allowed = {"pre":"Do flat correction before assembly", "post":"Do the masking and interpolation after assembly", "None":"Don't find and interpolate defects"},
        default = "post",
    )
    doMaskAndInterpNan = pexConfig.ChoiceField(
        dtype = str,
        doc = "Do look for unmasked NaN/Infs in the image.  If found mask and interpolate them.",
        allowed = {"pre":"Do flat correction before assembly", "post":"Do the masking and interpolation after assembly.", "None":"Don't look for nans/infs"},
        default = "post",
    )
    transposeForInterpolation = pexConfig.Field(
        dtype = bool,
        doc = "Since defect interpolation happens along rows, some data may need to be transposed before interpolation",
        default = False,
    )
    
    
class IsrTask(pipeBase.Task):
    ConfigClass = IsrTaskConfig
    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.isr = Isr()
        self.methodList = ["doConversionForIsr", "doSaturationDetection", "doLinearization", "doOverscanCorrection", "doTrimExposure", "doBiasSubtraction", "doVariance", "doDarkCorrection", "doFringeCorrection", "doFlatCorrection", "doMaskAndInterpDefect", "doSaturationInterpolation", "doMaskAndInterpNan"]
        self.preList = []
        self.postList = []
        for methodname in self.methodList:
            mvalue = getattr(self.config, methodname)
            assembleReached = False
            if mvalue == "pre" and not assembleReached:
                self.preList.append(methodname)
            elif mvalue == "post" and self.config.doAssembly:
                assembleReached = True
                self.postList.append(methodname)
            elif mvalue == "None":
                continue
            else:
                raise ValueError("Got invalid value for %s: %s\n Either you have tagged a pre method after assembly or you have post methods with assembly turned off."%(methodname, mvalue))

    def run(self, sDataRef):
        """Do instrument signature removal on an exposure: saturation, bias, overscan, dark, flat, fringe correction

        @param sDataRef daf.persistence.butlerSubset.ButlerDataRef of the data to be processed
        @return a pipeBase.Struct with fields:
        - exposure: the exposure after application of ISR
        """

        isChannel = sDataRef.butlerSubset.level == 'channel'
        if isChannel and (not self.postList or self.config.doAssembly):
            raise RuntimeError("Must provide more than one channel to do assembly and post assembly processing")
        exposureList = []
        #Get the list of amps and the associated calibration products for doing pre assembly processing
        ampList, calibList = self.makeAmpList(sDataRef)
        for ampExp, calibDict in zip(ampList, calibList):
            amp = cameraGeom.cast_Amp(ampExp.getDetector())
            workingExposure = ampExp.Factory(ampExp, True)
            #Run each correction for each amp
            for m in self.preList:
                workingExposure = getattr(self, m)(workingExposure, calibDict, amp)
            exposureList.append(workingExposure)

        if not isChannel and self.config.doAssembly:
            workingExposure = self.doCcdAssembly(exposureList)
            ccd = cameraGeom.cast_Ccd(workingExposure.getDetector())

            #Get calibs for the post assembly corrections
            calibDict = self.makeCalibDict(sDataRef, ccd, self.postList)

            #Run post assembly corrections
            for m in self.postList:
                workingExposure = getattr(self,m)(workingExposure, calibDict, ccd)
            if self.config.doWrite:
                #Persist data
                sDataRef.put(workingExposure, "postISRCCD")

        else:
            if self.config.doWrite:
                #Persist each of the corrected amps if assembly is not done
                for exp in exposureList:
                    sDataRef.put(exp, "postISR")

        return pipeBase.Struct(exposure=workingExposure)

    def checkIsAmp(self, detector):
        """Check if a detector is of type cameraGeom.Amp

        @param Detector cameraGeom.Detector to be checked
        @return boolean True if Amp, else False
        """
        if not isinstance(detector, cameraGeom.Amp):
            return False
        else:
            return True

    def makeAmpList(self, dataRef):
        """Make a list of amps and the associated calibration products

        @param dataRef ButlerDataRef from which to retrieve relevant data 
        @return ampList, calibList list of amp level raw data and a list of dictionaries of the calibration products for each amp
        """
        ampList = []
        calibList = []
        if dataRef.butlerSubset.level == 'channel':
            ampExp = dataRef.get('raw')
            amp = cameraGeom.cast_Amp(ampExp.getDetector())
            ampList.append(ampExp)
            calibList.append(self.makeCalibDict(dataRef, amp, self.preList))
        elif dataRef.butlerSubset.level == 'sensor':
            if 'channel' in dataRef.subLevels():
                for c in dataRef.subItems('channel'):
                    ampExp = c.get('raw')
                    amp = cameraGeom.cast_Amp(ampExp.getDetector())
                    ampList.append(ampExp)
                    calibList.append(self.makeCalibDict(c, amp, self.preList))
            else:
                ccdExp = dataRef.get('raw')
                ccd = cameraGeom.cast_Ccd(ccdExp.getDetector())
                for amp in ccd:
                    ampList.append(ccdExp.Factory(ccdExp, amp.getDiskAllPixels()))
                    calibList.append(self.makeCalibDict(dataRef, amp, self.preList))
        else:
            raise RuntimeException("The butler data reference must be at the channel or sensor level.  The supplied data reference is at the %s level."%(dataRef.butlerSubset.level))
        return ampList, calibList 

    def makeCalibDict(self, dataRef, detector, methodList):
        """Make a dictionary of calibration products for a detector given the methods that will be run on it

        @param dataRef ButlerDataRef from which to retrieve relevant data 
        @param detector cameraGeom.Detector to use for getting bounding boxes  
        @param methodList list of method names for retrieving calibration products
        @return calibList dictionary of calibration products
        """
        ret = {}
        required = {"doBiasSubtraction": "bias",
                    "doDarkCorrection": "dark",
                    "doFlatCorrection": "flat",
                    "doFringeCorrection": "fringe",
                    }
        for method in required.keys():
            if method in methodList:
                calib = required[method]
                calibExp = dataRef.get(calib)
                if not self.config.calibsAreTrimmed:
                    ret[calib] = calibExp.Factory(calibExp, detector.getDiskDataSec())
                else:
                    ret[calib] = calibExp
        return ret

    def doConversionForIsr(self, exposure, calibSet, detector):
        """Convert from int to float image for ISR processing

        @param exposure afwImage.Exposure to operate on
        @param calibSet dictionary of calibration products
        @param detector cameraGeom.Detector for the exposure
        @return newexposure afwImage.Exposure corrected exposure
        """
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

    def doVariance(self, exposure, calibSet, detector):
        if not self.checkIsAmp(detector):
            raise RuntimeError("This method must be executed on an amp.")
        self.isr.updateVariance(exposure.getMaskedImage(), detector.getElectronicParams().getGain())
        return exposure

    def doCrosstalkCorrection(self, exposure, calibSet, detector):
        if not self.checkIsAmp(detector):
            raise RuntimeError("This method must be executed on an amp.")
        pass

    def doSaturationCorrection(self, exposure, calibSet, detector):
        if not self.checkIsAmp(detector):
            raise RuntimeError("This method must be executed on a single channel.")
        fwhm = self.config.fwhm
        grow = self.config.growSaturationFootprintSize
        maskname = self.config.saturatedMaskName
        ep = detector.getElectronicParams()
        satvalue = ep.getSaturationLevel()
        self.isr.saturationCorrection(exposure.getMaskedImage(), satvalue, fwhm, growFootprints=grow, maskName=maskname)
        return exposure

    def doSaturationDetection(self, exposure, calibSet, detector):
        if not self.checkIsAmp(detector):
            raise RuntimeError("This method must be executed on an amp.")
        ep = detector.getElectronicParams()
        satvalue = ep.getSaturationLevel()
        maskname = self.config.saturatedMaskName
        self.isr.makeThresholdMask(exposure.getMaskedImage(), satvalue, growFootprints=0, maskName=maskname)
        return exposure

    def doSaturationInterpolation(self, exposure, calibSet, detector):
        #Don't loop over amps since saturation can cross amp boundaries
        maskname = self.config.saturatedMaskName
        fwhm = self.config.fwhm
        grow = self.config.growSaturationFootprintSize
        if self.config.transposeForInterpolation:
            mi = self.isr.transposeMaskedImage(exposure.getMaskedImage())
            self.isr.interpolateFromMask(mi, fwhm, growFootprints=grow, maskName=maskname)
            mi = self.isr.transposeMaskedImage(mi)
            exposure.setMaskedImage(mi)
        else:
            self.isr.interpolateFromMask(exposure.getMaskedImage(), fwhm, growFootprints=grow, maskName=maskname)
        return exposure
    
    def doMaskAndInterpDefect(self, exposure, calibSet, detector):
        #Don't loop over amps since defects could cross amp boundaries
        fwhm = self.config.fwhm
        grow = self.config.growDefectFootprintSize
        defectBaseList = detector.getDefects()
        defectList = measAlg.DefectListT()
        #mask bad pixels in the camera class
        #create master list of defects and add those from the camera class
        for d in defectBaseList:
            bbox = d.getBBox()
            nd = measAlg.Defect(bbox)
            defectList.append(nd)
        self.isr.maskPixelsFromDefectList(exposure.getMaskedImage(), defectList, maskName='BAD')
        defectList = self.isr.getDefectListFromMask(exposure.getMaskedImage(), maskName='BAD', growFootprints=grow)
        if self.config.transposeForInterpolation:
            mi = self.isr.transposeMaskedImage(exposure.getMaskedImage())
            defectList = self.isr.transposeDefectList(defectList)
            self.isr.interpolateDefectList(mi, defectList, fwhm)
            mi = self.isr.transposeMaskedImage(mi)
            exposure.setMaskedImage(mi)
        else:
            self.isr.interpolateDefectList(exposure.getMaskedImage(), defectList, fwhm)
        return exposure

    def doMaskAndInterpNan(self, exposure, calibSet, detector):
        #Don't loop over amps since nans could cross amp boundaries
        fwhm = self.config.fwhm
        grow = self.config.growDefectFootprintSize
        #find unmasked bad pixels and mask them
        exposure.getMaskedImage().getMask().addMaskPlane("UNMASKEDNAN") 
        unc = isrLib.UnmaskedNanCounterF()
        unc.apply(exposure.getMaskedImage())
        nnans = unc.getNpix()
        self.metadata.set("NUMNANS", nnans)
        if not nnans == 0:
		raise RuntimeError("There were %i unmasked NaNs"%(nnans))
        #get footprints of bad pixels not in the camera class
        if nnans > 0:
            undefects = self.isr.getDefectListFromMask(exposure.getMaskedImage(), maskName='UNMASKEDNAN', growFootprints=grow)
            #interpolate all bad pixels
            if self.config.transposeForInterpolation:
                mi = self.isr.transposeMaskedImage(exposure.getMaskedImage())
                defectList = self.isr.transposeDefectList(uudefects)
                self.isr.interpolateDefectList(mi, defectList, fwhm)
                mi = self.isr.transposeMaskedImage(mi)
                exposure.setMaskedImage(mi)
            else:
                self.isr.interpolateDefectList(exposure.getMaskedImage(), uudefects, fwhm)
        return exposure

    '''
    def doLinearization(self, exposure, calibSet, detector):
        linearizer = Linearization(calibSet['linearityFile'])
        linearizer.apply(exposure)
        return exposure
    '''
    
    def doOverscanCorrection(self, exposure, calibSet, detector):
        if not self.checkIsAmp(detector):
            raise RuntimeError("This method must be executed on an amp.")
        fittype = self.config.overscanFitType
        polyorder = self.config.overscanPolyOrder
        expImage = exposure.getMaskedImage().getImage()
        overscan = expImage.Factory(expImage, detector.getDiskBiasSec())
        self.isr.overscanCorrection(exposure.getMaskedImage(), overscan, fittype=fittype, polyorder=polyorder,
                                        imageFactory=overscan.Factory)
        return exposure
 
    def doTrimExposure(self, exposure, calibSet, detector):
        if not self.checkIsAmp(detector):
            raise RuntimeError("This method must be executed on an amp.")
        return exposure.Factory(exposure, detector.getDiskDataSec())

    def doBiasSubtraction(self, exposure, calibSet, detector):
        biasExposure = calibSet['bias']
        self.isr.biasCorrection(exposure.getMaskedImage(), biasExposure.getMaskedImage())
        return exposure 

    def doDarkCorrection(self, exposure, calibSet, detector):
        darkexposure = calibSet['dark']
        darkscaling = darkexposure.getCalib().getExptime()
        expscaling = exposure.getCalib().getExptime()
        
        self.isr.darkCorrection(exposure.getMaskedImage(), darkexposure.getMaskedImage(), expscaling, darkscaling)
        return exposure

    def doFringeCorrection(self, exposure, calibSet, detector):
        pass


    def doFlatCorrection(self, exposure, calibSet, detector):
        flatfield = calibSet['flat']
        scalingtype = self.config.flatScalingType
        scalingvalue = self.config.flatScalingValue

        self.isr.flatCorrection(exposure.getMaskedImage(), flatfield.getMaskedImage(), scalingtype, scaling = scalingvalue)   
        return exposure

    def doIlluminationCorrection(self, exposure, calibSet, detector):
        pass

    def doCcdAssembly(self, exposureList):
        renorm = self.config.reNormAssembledCcd
        setgain = self.config.setGainAssembledCcd
        k2rm = self.config.keysToRemoveFromAssembledCcd
        assembler = CcdAssembler(exposureList, reNorm=renorm, setGain=setgain, keysToRemove=k2rm)
        return assembler.assembleCcd()
