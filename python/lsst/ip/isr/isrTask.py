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
import lsst.meas.algorithms as measAlg
import lsst.afw.cameraGeom as cameraGeom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from . import isr
from .isrLib import UnmaskedNanCounterF
from .assembleCcdTask import AssembleCcdTask

class IsrTaskConfig(pexConfig.Config):
    doBias = pexConfig.Field(
        dtype = bool,
        doc = "Apply bias frame correction?",
        default = True,
    )
    doDark = pexConfig.Field(
        dtype = bool,
        doc = "Apply dark frame correction?",
        default = True,
    )
    doFlat = pexConfig.Field(
        dtype = bool,
        doc = "Apply flat field correction?",
        default = True,
    )
    doWrite = pexConfig.Field(
        dtype = bool,
        doc = "Persist visitCCD?",
        default = True,
    )
    assembleCcd = pexConfig.ConfigurableField(
        target = AssembleCcdTask,
        doc = "CCD assembly task",
    )
    fwhm = pexConfig.Field(
        dtype = float,
        doc = "FWHM of PSF (arcsec)",
        default = 1.0,
    )
    saturatedMaskName = pexConfig.Field(
        dtype = str,
        doc = "Name of mask plane to use in saturation detection and interpolation",
        default = "SAT",
    )
    flatScalingType = pexConfig.ChoiceField(
        dtype = str,
        doc = "The method for scaling the flat on the fly.",
        default = 'USER',
        allowed = {
            "USER":   "Scale by flatUserScale",
            "MEAN":   "Scale by the inverse of the mean",
            "MEDIAN": "Scale by the inverse of the median",
        },
    )
    flatUserScale = pexConfig.Field(
        dtype = float,
        doc = "If flatScalingType is 'USER' then scale flat by this amount; ignored otherwise",
        default = 1.0,
    )
    overscanFitType = pexConfig.ChoiceField(
        dtype = str,
        doc = "The method for fitting the overscan bias level.",
        default = 'MEDIAN',
        allowed = {
            "POLY": "Fit polynomial to the longest axis of the overscan region",
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
        default = 0,
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
    
    
class IsrTask(pipeBase.Task):
    ConfigClass = IsrTaskConfig
    _DefaultName = "isr"

    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.makeSubtask("assembleCcd")
        self.transposeForInterpolation = False

    @pipeBase.timeMethod
    def run(self, sensorRef):
        """Perform instrument signature removal on an exposure
        
        Steps include:
        - Detect saturation, apply overscan correction, bias, dark and flat
        - Perform CCD assembly
        - Interpolate over defects, saturated pixels and any remaining NaNs
        - Persist the ISR-corrected exposure as "visitCcd" if config.doWrite is True

        @param sensorRef daf.persistence.butlerSubset.ButlerDataRef of the data to be processed
        @return a pipeBase.Struct with fields:
        - exposure: the exposure after application of ISR
        """
        self.log.log(self.log.INFO, "Performing ISR on sensor %s" % (sensorRef.dataId))
        ccdExp = sensorRef.get('raw')
        ccd = cameraGeom.cast_Ccd(ccdExp.getDetector())

        ccdExp = isr.floatImageFromInt(ccdExp)
        
        for amp in ccd:
            ampExp = ccdExp.Factory(ccdExp, amp.getDiskAllPixels())
            amp = cameraGeom.cast_Amp(ampExp.getDetector())
    
            self.saturationDetection(ampExp, amp)

            self.overscanCorrection(ampExp, amp)

        if self.config.doBias:
            self.biasCorrection(ccdExp, sensorRef)
        
        if self.config.doDark:
            self.darkCorrection(ccdExp, sensorRef)
        
        self.updateVariance(ccdExp, ccd)
        
        if self.config.doFlat:
            self.flatCorrection(ccdExp, sensorRef)
        
        ccdExp = self.assembleCcd.assembleCcd(ccdExp)
        
        self.maskAndInterpDefect(ccdExp)
        
        self.saturationInterpolation(ccdExp)
        
        self.maskAndInterpNan(ccdExp)

        if self.config.doWrite:
            sensorRef.put(ccdExp, "visitCCD")
        
        self.display("visitCCD", ccdExp)

        return pipeBase.Struct(
            exposure = ccdExp,
        )

    def checkIsAmp(self, amp):
        """Check if a amp is of type cameraGeom.Amp

        @param Detector cameraGeom.Detector to be checked
        @return True if Amp, else False
        """
        return isinstance(amp, cameraGeom.Amp)
    
    def biasCorrection(self, exposure, dataRef):
        """Apply bias correction in place
    
        @param[in,out]  exposure        exposure to process
        @param[in]      dataRef         data reference at same level as exposure
        """
        biasMaskedImage = dataRef.get("bias").getMaskedImage()
        isr.biasCorrection(ccdExp.getMaskedImage(), biasMaskedImage)

    def darkCorrection(self, exposure, dataRef):
        """Apply dark correction in place
    
        @param[in,out]  exposure        exposure to process
        @param[in]      dataRef         data reference at same level as exposure
        """
        darkExposure = dataRef.get("dark")
        darkCalib = darkExposure.getCalib()
        isr.darkCorrection(
            maskedImage = exposure.getMaskedImage(),
            darkMaskedImage = darkExposure.getMaskedImage(),
            expScale = darkCalib.getExptime(),
            darkScale = darkCalib.getExptime(),
        )
    
    def updateVariance(self, ccdExposure, ccd):
        """Set the variance plane based on the image plane
        
        @param[in,out]  ccdExposure     exposure to process
        @param[in]      ccd             CCD detector information
        """
        isr.updateVariance(
            maskedImage = ccdExp.getMaskedImage(),
            gain = ccd.getElectronicParams().getGain(),
            readNoise = ccd.getElectronicParams().getReadNoise(),
        )

    def flatCorrection(self, exposure, dataRef):
        """Apply flat correction in place
    
        @param[in,out]  exposure        exposure to process
        @param[in]      dataRef         data reference at same level as exposure
        """
        flatfield = dataRef['flat']
        isr.flatCorrection(
            maskedImage = exposure.getMaskedImage(),
            flatMaskedImage = flatfield.getMaskedImage(),
            scalingType = self.config.flatScalingType,
            userScale = self.config.flatUserScale,
        )

    def saturationCorrection(self, ampExposure, amp):
        """Perform saturation detection and interpolation as a single step.
        
        parm[in,out]    ampExposure exposure to process
        @param[in]      amp         amplifier device data

        @warning Only suitable for CCDs with a single amplifier; otherwise you don't correctly handle
        bleed trails that cross amplifier boundaries. For more complex CCDs: call saturationDetection
        for each amplifier, then saturationInterpolation on the assembled CCD.
        """
        if not self.checkIsAmp(amp):
            raise RuntimeError("This method must be executed on an amp.")
        isr.saturationCorrection(
            maskedImage = ampExposure.getMaskedImage(),
            saturation = amp.getElectronicParams().getSaturationLevel(),
            fwhm = self.config.fwhm,
            growFootprints = self.config.growSaturationFootprintSize,
            maskName = self.config.saturatedMaskName,
        )
        return ampExposure

    def saturationDetection(self, ampExposure, amp):
        """Detect saturated pixels, in place
        
        Saturated pixels are identified by the SAT mask plane (e.g. by saturationDetection).
        
        @param[in,out]  ampExposure amplifier ampExposure to process
        @param[in]      amp         amplifier device data
        """
        if not self.checkIsAmp(amp):
            raise RuntimeError("This method must be executed on an amp.")
        isr.makeThresholdMask(
            maskedImage = ampExposure.getMaskedImage(),
            saturation = amp.getElectronicParams().getSaturationLevel(),
            growFootprints = 0,
            maskName = self.config.saturatedMaskName,
        )
        return ampExposure

    def saturationInterpolation(self, ccdExposure):
        """Interpolate over saturated pixels, in place
        
        @param[in,out]  ccdExposure     exposure to process

        @warning:
        - Call saturationDetection first, so that saturated pixels have been identified in the mask.
        - Call this after CCD assembly, since saturated regions may cross amplifier boundaries
        """
        if self.transposeForInterpolation:
            mi = isr.transposeMaskedImage(ccdExposure.getMaskedImage())
            isr.interpolateFromMask(
                maskedImage = mi,
                fwhm = self.config.fwhm,
                growFootprints = self.config.growSaturationFootprintSize,
                maskName = self.config.saturatedMaskName,
            )
            mi = isr.transposeMaskedImage(mi)
            ccdExposure.setMaskedImage(mi)
        else:
            isr.interpolateFromMask(
                maskedImage = ccdExposure.getMaskedImage(),
                fwhm = self.config.fwhm,
                growFootprints = self.config.growSaturationFootprintSize,
                maskName = self.config.saturatedMaskName,
            )
    
    def maskAndInterpDefect(self, ccdExposure):
        """Mask defects and interpolate over them, in place

        @param[in,out]  ccdExposure     exposure to process
        
        @warning: call this after CCD assembly, since defects may cross amplifier boundaries
        """
        ccd = cameraGeom.cast_Ccd(ccdExp.getDetector())
        fwhm = self.config.fwhm
        grow = self.config.growDefectFootprintSize
        defectBaseList = ccd.getDefects()
        defectList = measAlg.DefectListT()
        #mask bad pixels in the camera class
        #create master list of defects and add those from the camera class
        for d in defectBaseList:
            bbox = d.getBBox()
            nd = measAlg.Defect(bbox)
            defectList.append(nd)
        isr.maskPixelsFromDefectList(ccdExposure.getMaskedImage(), defectList, maskName='BAD')
        defectList = isr.getDefectListFromMask(ccdExposure.getMaskedImage(), maskName='BAD', growFootprints=grow)
        if self.transposeForInterpolation:
            mi = isr.transposeMaskedImage(ccdExposure.getMaskedImage())
            defectList = isr.transposeDefectList(defectList)
            isr.interpolateDefectList(mi, defectList, fwhm)
            mi = isr.transposeMaskedImage(mi)
            ccdExposure.setMaskedImage(mi)
        else:
            isr.interpolateDefectList(ccdExposure.getMaskedImage(), defectList, fwhm)

    def maskAndInterpNan(self, exposure):
        """Mask NaNs and interpolate over them, in place

        @param[in,out]  exposure        exposure to process
        """
        fwhm = self.config.fwhm
        grow = self.config.growDefectFootprintSize
        #find unmasked bad pixels and mask them
        exposure.getMaskedImage().getMask().addMaskPlane("UNMASKEDNAN") 
        unc = UnmaskedNanCounterF()
        unc.apply(exposure.getMaskedImage())
        nnans = unc.getNpix()
        self.metadata.set("NUMNANS", nnans)

        #get footprints of bad pixels not in the camera class
        if nnans > 0:
            self.log.log(self.log.WARN, "There were %i unmasked NaNs" % (nnans,))
            undefects = isr.getDefectListFromMask(exposure.getMaskedImage(), maskName='UNMASKEDNAN', growFootprints=grow)
            #interpolate all bad pixels
            if self.transposeForInterpolation:
                mi = isr.transposeMaskedImage(exposure.getMaskedImage())
                defectList = isr.transposeDefectList(uudefects)
                isr.interpolateDefectList(mi, defectList, fwhm)
                mi = isr.transposeMaskedImage(mi)
                exposure.setMaskedImage(mi)
            else:
                isr.interpolateDefectList(exposure.getMaskedImage(), uudefects, fwhm)

    def overscanCorrection(self, ampExposure, amp):
        """Apply overscan correction, in place

        @param[in,out]  ampExposure amplifier ampExposure to process
        @param[in]      amp         amplifier device data
        """
        if not self.checkIsAmp(amp):
            raise RuntimeError("This method must be executed on an amp.")
        expImage = ampExposure.getMaskedImage().getImage()
        overscan = expImage.Factory(expImage, amp.getDiskBiasSec())
        isr.overscanCorrection(
            maskedImage = ampExposure.getMaskedImage(),
            overscanData = ovescan,
            fitType = self.config.overscanFitType,
            polyOrder = self.config.overscanPolyOrder,
            imageFactory = overscan.Factory,
        )
