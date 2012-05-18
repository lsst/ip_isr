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
        doc = "Persist postISRCCD?",
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
    
    
class IsrTask(pipeBase.CmdLineTask):
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
        - Persist the ISR-corrected exposure as "postISRCCD" if config.doWrite is True

        @param sensorRef daf.persistence.butlerSubset.ButlerDataRef of the data to be processed
        @return a pipeBase.Struct with fields:
        - exposure: the exposure after application of ISR
        """
        self.log.log(self.log.INFO, "Performing ISR on sensor %s" % (sensorRef.dataId))
        ccdExposure = sensorRef.get('raw')
        ccd = cameraGeom.cast_Ccd(ccdExposure.getDetector())
    
        ccdExposure = self.convertIntToFloat(ccdExposure)
        
        for amp in ccd:
            self.saturationDetection(ccdExposure, amp)

            self.overscanCorrection(ccdExposure, amp)
        
        ccdExposure = self.assembleCcd.assembleCcd(ccdExposure)
        ccd = cameraGeom.cast_Ccd(ccdExposure.getDetector())

        if self.config.doBias:
            self.biasCorrection(ccdExposure, sensorRef)
        
        if self.config.doDark:
            self.darkCorrection(ccdExposure, sensorRef)
        
        for amp in ccd:
            ampExposure = ccdExposure.Factory(ccdExposure, amp.getAllPixels(True), afwImage.PARENT)

            self.updateVariance(ampExposure, amp)
        
        if self.config.doFlat:
            self.flatCorrection(ccdExposure, sensorRef)
        
        self.maskAndInterpDefect(ccdExposure)
        
        self.saturationInterpolation(ccdExposure)
        
        self.maskAndInterpNan(ccdExposure)

        if self.config.doWrite:
            sensorRef.put(ccdExposure, "postISRCCD")
        
        self.display("postISRCCD", ccdExposure)

        return pipeBase.Struct(
            exposure = ccdExposure,
        )

    def checkIsAmp(self, amp):
        """Check if a amp is of type cameraGeom.Amp

        @param Detector cameraGeom.Detector to be checked
        @return True if Amp, else False
        """
        return isinstance(amp, cameraGeom.Amp)

    def convertIntToFloat(self, exposure):
        """Convert an exposure from uint16 to float, set variance plane to 1 and mask plane to 0
        """
        if not isinstance(exposure, afwImage.ExposureU):
            raise Exception("ipIsr.convertImageForIsr: Expecting Uint16 image. Got %r" % (exposure,))

        newexposure = exposure.convertF()
        maskedImage = newexposure.getMaskedImage()
        varArray = maskedImage.getVariance().getArray()
        varArray[:,:] = 1
        maskArray = maskedImage.getMask().getArray()
        maskArray[:,:] = 0
        return newexposure
    
    def biasCorrection(self, exposure, dataRef):
        """Apply bias correction in place
    
        @param[in,out]  exposure        exposure to process
        @param[in]      dataRef         data reference at same level as exposure
        """
        biasMaskedImage = dataRef.get("bias").getMaskedImage()
        isr.biasCorrection(exposure.getMaskedImage(), biasMaskedImage)

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
    
    def updateVariance(self, ampExposure, amp):
        """Set the variance plane based on the image plane, plus amplifier gain and read noise
        
        @param[in,out]  ampExposure     exposure to process
        @param[in]      amp             amplifier detector information
        """
        if not self.checkIsAmp(amp):
            raise RuntimeError("This method must be executed on an amp.")
        isr.updateVariance(
            maskedImage = ampExposure.getMaskedImage(),
            gain = amp.getElectronicParams().getGain(),
            readNoise = amp.getElectronicParams().getReadNoise(),
        )

    def flatCorrection(self, exposure, dataRef):
        """Apply flat correction in place
    
        @param[in,out]  exposure        exposure to process
        @param[in]      dataRef         data reference at same level as exposure
        """
        flatfield = dataRef.get("flat")
        isr.flatCorrection(
            maskedImage = exposure.getMaskedImage(),
            flatMaskedImage = flatfield.getMaskedImage(),
            scalingType = self.config.flatScalingType,
            userScale = self.config.flatUserScale,
        )

    def saturationDetection(self, exposure, amp):
        """Detect saturated pixels and mask them using mask plane "SAT", in place
        
        @param[in,out]  exposure    exposure to process; only the amp DataSec is processed
        @param[in]      amp         amplifier device data
        """
        if not self.checkIsAmp(amp):
            raise RuntimeError("This method must be executed on an amp.")
        maskedImage = exposure.getMaskedImage()
        dataView = maskedImage.Factory(maskedImage, amp.getDiskDataSec(), afwImage.PARENT)
        isr.makeThresholdMask(
            maskedImage = dataView,
            threshold = amp.getElectronicParams().getSaturationLevel(),
            growFootprints = 0,
            maskName = self.config.saturatedMaskName,
        )

    def saturationInterpolation(self, ccdExposure):
        """Interpolate over saturated pixels, in place
        
        @param[in,out]  ccdExposure     exposure to process

        @warning:
        - Call saturationDetection first, so that saturated pixels have been identified in the "SAT" mask.
        - Call this after CCD assembly, since saturated regions may cross amplifier boundaries
        """
        if self.transposeForInterpolation:
            maskedImage = isr.transposeMaskedImage(ccdExposure.getMaskedImage())
            isr.interpolateFromMask(
                maskedImage = maskedImage,
                fwhm = self.config.fwhm,
                growFootprints = self.config.growSaturationFootprintSize,
                maskName = self.config.saturatedMaskName,
            )
            maskedImage = isr.transposeMaskedImage(maskedImage)
            ccdExposure.setMaskedImage(maskedImage)
        else:
            isr.interpolateFromMask(
                maskedImage = ccdExposure.getMaskedImage(),
                fwhm = self.config.fwhm,
                growFootprints = self.config.growSaturationFootprintSize,
                maskName = self.config.saturatedMaskName,
            )
    
    def maskAndInterpDefect(self, ccdExposure):
        """Mask defects using mask plane "BAD" and interpolate over them, in place

        @param[in,out]  ccdExposure     exposure to process
        
        @warning: call this after CCD assembly, since defects may cross amplifier boundaries
        """
        maskedImage = ccdExposure.getMaskedImage()
        ccd = cameraGeom.cast_Ccd(ccdExposure.getDetector())
        defectBaseList = ccd.getDefects()
        defectList = measAlg.DefectListT()
        # mask bad pixels in the camera class
        # create master list of defects and add those from the camera class
        for d in defectBaseList:
            bbox = d.getBBox()
            nd = measAlg.Defect(bbox)
            defectList.append(nd)
        isr.maskPixelsFromDefectList(maskedImage, defectList, maskName='BAD')
        if self.transposeForInterpolation:
            maskedImage = isr.transposeMaskedImage(maskedImage)
            defectList = isr.transposeDefectList(defectList)
            isr.interpolateDefectList(
                maskedImage = maskedImage,
                defectList = defectList,
                fwhm = self.config.fwhm,
            )
            maskedImage = isr.transposeMaskedImage(maskedImage)
            ccdExposure.setMaskedImage(maskedImage)
        else:
            isr.interpolateDefectList(
                maskedImage = maskedImage,
                defectList = defectList,
                fwhm = self.config.fwhm,
            )

    def maskAndInterpNan(self, exposure):
        """Mask unmasked NaNs using mask plane "UNMAKSEDNAN" and interpolate over them, in place

        @param[in,out]  exposure        exposure to process
        """
        maskedImage = exposure.getMaskedImage()
        
        # Find and mask unmasked NaNs
        maskedImage.getMask().addMaskPlane("UNMASKEDNAN") 
        unc = UnmaskedNanCounterF()
        unc.apply(exposure.getMaskedImage())
        numNans = unc.getNpix()
        self.metadata.set("NUMNANS", numNans)

        # Interpolate over these previously-unmasked NaNs
        if numNans > 0:
            self.log.log(self.log.WARN, "There were %i unmasked NaNs" % (numNans,))
            nanDefectList = isr.getDefectListFromMask(
                maskedImage = maskedImage,
                maskName = 'UNMASKEDNAN',
                growFootprints = 0,
            )
            if self.transposeForInterpolation:
                maskedImage = isr.transposeMaskedImage(exposure.getMaskedImage())
                defectList = isr.transposeDefectList(nanDefectList)
                isr.interpolateDefectList(
                    maskedImage = maskedImage,
                    defectList = nanDefectList,
                    fwhm = self.config.fwhm,
                )
                maskedImage = isr.transposeMaskedImage(maskedImage)
                exposure.setMaskedImage(maskedImage)
            else:
                isr.interpolateDefectList(
                    maskedImage = exposure.getMaskedImage(),
                    defectList = nanDefectList,
                    fwhm = self.config.fwhm,
                )

    def overscanCorrection(self, exposure, amp):
        """Apply overscan correction, in place

        @param[in,out]  exposure    exposure to process; must include both DataSec and BiasSec pixels
        @param[in]      amp         amplifier device data
        """
        if not self.checkIsAmp(amp):
            raise RuntimeError("This method must be executed on an amp.")
        maskedImage = exposure.getMaskedImage()
        dataView = maskedImage.Factory(maskedImage, amp.getDiskDataSec(), afwImage.PARENT)

        expImage = exposure.getMaskedImage().getImage()
        overscanImage = expImage.Factory(expImage, amp.getDiskBiasSec(), afwImage.PARENT)

        isr.overscanCorrection(
            ampMaskedImage = dataView,
            overscanImage = overscanImage,
            fitType = self.config.overscanFitType,
            polyOrder = self.config.overscanPolyOrder,
        )
