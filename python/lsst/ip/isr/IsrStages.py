import time, re, os, math, pdb
import lsst.utils           as utils
import lsst.afw.detection   as afwDetection
import lsst.afw.image       as afwImage
import lsst.afw.math        as afwMath
import lsst.meas.algorithms as algorithms
import lsst.pex.logging     as pexLog
from lsst.pex.logging import Log
import lsst.pex.policy      as pexPolicy
import lsst.pex.exceptions  as pexExcept
import lsst.meas.algorithms.defects as measDefects
import lsst.sdqa            as sdqa

import lsst.afw.display.ds9 as ds9

# relative imports
import isrLib

#
### Actual Pipeline stage to run the ISR
#
from lsst.pex.harness.Stage import Stage

class IsrStage(Stage):

    def __init__(self, stageId = -1, policy = None):
        # call base constructor
        Stage.__init__(self,stageId, policy)
        # initialize a log
        self.log = Log(Log.getDefaultLog(), 
                "lsst.ip.isr.IsrStages")

    def process(self):
        self.activeClipboard = self.inputQueue.getNextDataset()

        # 
        inputImageKey    = self._policy.get('inputImageKey')
        inputMetadataKey = self._policy.get('inputMetadataKey')
        calibDataKey     = self._policy.get('calibDataKey')
        inputImage       = self.activeClipboard.get(inputImageKey)
        inputMetadata    = self.activeClipboard.get(inputMetadataKey)
        calibData        = self.activeClipboard.get(calibDataKey)

        # Grab the necessary calibration data products
        bias             = self.activeClipboard.get('biasExposure')
        dark             = self.activeClipboard.get('darkExposure')
        flat             = self.activeClipboard.get('flatExposure')

        # Get the image's (e.g. Amp's) origin on the master (e.g. CCD)
        # image.  
        #ampBBox          = self.activeClipboard.get('ampBBox')

        # Step 1 : create an exposure
        inputExposure    = ExposureFromInputData(inputImage, inputMetadata, False)
        mi      = inputExposure.getMaskedImage()
        self.log.log(Log.DEBUG, "Input MaskedImage X0, Y0: %d %d" % (mi.getX0(), mi.getY0()))

        ###
        # Isr Substages
        #

        isrPolicy = self._policy.get('isrPolicy')
        
        # Linearize
        # Note - a replacement lookup table requires an integer image

        if keyIsSet('linearize', isrPolicy):
            linearizePath    = calibData.get('linearizePath')
            linearityPolicy = pexPolicy.Policy.createPolicy(linearizePath)
            linearityTable  = LookupTableFromPolicy(linearityPolicy)
            Linearization(inputExposure, isrPolicy, linearityTable)

        # Saturation correction
        if keyIsSet('saturationCorrection', isrPolicy):
            SaturationCorrection(inputExposure, isrPolicy)
        
        # Overscan correction
        if keyIsSet('overscanCorrection', isrPolicy):
            OverscanCorrection(inputExposure, isrPolicy)

        # Trim; yields new exposure
        if keyIsSet('trim', isrPolicy):
            calibratedExposure = TrimNew(inputExposure, isrPolicy)

        mi      =  calibratedExposure.getMaskedImage()
        self.log.log(Log.DEBUG, "Trimmed MaskedImage X0, Y0: %d %d" % (mi.getX0(), mi.getY0()))
        
        # Merge the metadata
        calibratedExposure.getMetadata().combine(inputExposure.getMetadata())

        mi      =  calibratedExposure.getMaskedImage()
        self.log.log(Log.DEBUG, "md combined MaskedImage X0, Y0: %d %d" % (mi.getX0(), mi.getY0()))
        

        # Bias correct
        if keyIsSet('biasCorrection', isrPolicy):
            BiasCorrection(calibratedExposure, bias, isrPolicy)

        # Dark correct
        if keyIsSet('darkCorrection', isrPolicy):
            DarkCorrection(calibratedExposure, dark, isrPolicy)

        # Flat field
        if keyIsSet('flatCorrection', isrPolicy):
            FlatCorrection(calibratedExposure, flat, isrPolicy)

        # Fringe; not for DC3a
        # fringeImage    = afwImage.ImageF(fringePath)
        # fringeMetadata = afwImage.readMetadata(fringePath)
        # fringe = ExposureFromInputData(fringeImage, fringeMetadata, makeWcs=False, policy=isrPolicy)
        
        # Finally, mask bad pixels
        # Finally, grab the additional information in calibDataKey
        if keyIsSet('maskBadPixels', isrPolicy):
            defectPath       = calibData.get('defectPath')
            defectList = measDefects.policyToBadRegionList(defectPath)
            #
            # The defects file is in amp coordinates, and we need to shift to the CCD frame
            #
            dx, dy = calibratedExposure.getMaskedImage().getXY0()
            for defect in defectList:
                defect.shift(dx, dy)

            MaskBadPixelsDef(calibratedExposure, isrPolicy, defectList)

        # And cosmic rays
        if keyIsSet('cosmicRayCorrection', isrPolicy):
            CrRejection(calibratedExposure, isrPolicy)
        mi      =  calibratedExposure.getMaskedImage()
        self.log.log(Log.DEBUG, "Output MaskedImage X0, Y0: %d %d" % (mi.getX0(), mi.getY0()))
        
        calibratedExposureKey = self._policy.get('calibratedExposureKey')
        self.activeClipboard.put(calibratedExposureKey, calibratedExposure)
        self.outputQueue.addDataset(self.activeClipboard)

        

        #
        # Isr Substages
        ###


    

#
### STAGE : Assemble Exposure from input Image
#

def CalculateSdqaRatings(exposure,
                         stageName   = 'lsst.ip.isr.calculatesdqaratings'):
    sdqaRatingSet = sdqa.SdqaRatingSet()
    counter       = isrLib.CountMaskedPixelsF()
    mi   = exposure.getMaskedImage()

    bitmaskSat = mi.getMask().getPlaneBitMask('SAT')
    counter.apply( mi, bitmaskSat )
    satRating  = sdqa.SdqaRating('ip.isr.numSaturatedPixels', counter.getCount(), 0, sdqa.SdqaRating.AMP)
    sdqaRatingSet.push_back(satRating)
    pexLog.Trace(stageName, 4, 'Found %d saturated pixels' % (counter.getCount()))
    
    bitmaskCr  = mi.getMask().getPlaneBitMask('CR')
    counter.apply( mi, bitmaskCr )
    crRating   = sdqa.SdqaRating('ip.isr.numCosmicRayPixels', counter.getCount(), 0, sdqa.SdqaRating.AMP)
    sdqaRatingSet.push_back(crRating)
    pexLog.Trace(stageName, 4, 'Found %d pixels with cosmic rays' % (counter.getCount()))

    return sdqaRatingSet

#
### STAGE : Assemble Exposure from input Image
#

def ExposureFromInputData(image, metadata,
                          makeWcs     = True,
                          policy      = None,
                          defaultGain = 1.0,
                          stageName   = 'lsst.ip.isr.exposurefrominputdata'):

    # Generate an empty mask
    mask = afwImage.MaskU(image.getDimensions())
    mask.set(0)

    # Generate a variance from the image pixels and gain
    var  = afwImage.ImageF(image, True)
    
    if metadata.exists('gain'):
        gain = metadata.get('gain')
    elif policy:
        filenameKeyword = policy.get('filenameKeyword')
        filename        = metadata.get(filenameKeyword)
        if policy.exists('defaultGainKeyword'):
            gainKeyword = policy.get('defaultGainKeyword')
            if metadata.exists(gainKeyword):
                gain = metadata.get(gainKeyword)
            else:
                pexLog.Trace(stageName, 4, 'Using default gain=%f for %s' % (defaultGain, filename))
                gain = defaultGain
        else:
            pexLog.Trace(stageName, 4, 'Using default gain=%f for %s' % (defaultGain, filename))
            gain = defaultGain
    else:
        pexLog.Trace(stageName, 4, 'Using default gain=%f' % (defaultGain))
        gain = defaultGain
    # Normalize by the gain
    var /= gain

    # makeMaskedImage() will make a MaskedImage with the same type as Image
    mi   = afwImage.makeMaskedImage(image, mask, var)
#    mi.setXY0(ampBBox.getX0(), ampBBox.getY0())

    if makeWcs:
        # Extract the Wcs info from the input metadata
        wcs      = afwImage.Wcs(metadata)
    else:
        wcs      = afwImage.Wcs()
        
    # makeExposure will make an Exposure with the same type as MaskedImage
    exposure = afwImage.makeExposure(mi, wcs)
    exposure.setMetadata(metadata)

    return exposure

#
### STAGE : Validation of the image sizes, contents, etc.
#

def ValidateCalibration(exposure, calibration, policy):
    """
    Make sure that the images are the same size, were derived from
    the same chunk of the focal plane, etc

    Things to check are :
     * Image Size (all)
     * From the same piece of the focal plane (illum, flat)
     * From the same piece of silicon (bad pixel mask, bias)
     * Through the same filter (dome flat)
     * Appropriate for the date range (anything time variable; dflats, etc)
    """
    
    pass

#
### STAGE : Bad pixel correction
#

def MaskFromDefects(dimensions, fpList):
    # the output LSST Mask image
    mask    = afwImage.MaskU(dimensions)
    mask.set(0)
    bitmask = mask.getPlaneBitMask('BAD')

    # set the bits
    afwDetection.setMaskFromFootprintList(mask, fpList, bitmask)

    return mask
    

def DefectsFromCfhtImage(fitsfile):
    # input bad pixel image
    image   = afwImage.ImageF(fitsfile)
    image  *= -1

    # turn into masked image for detection
    mi = afwImage.MaskedImageF(image)

    # find bad regions
    thresh    = afwDetection.Threshold(-0.5)
    ds        = afwDetection.DetectionSetF(mi, thresh)
    fpList    = ds.getFootprints()

    return fpList


def MaskBadPixelsFp(exposure, policy, fpList,
                    interpolate = True,
                    maskName    = 'BAD',
                    stageSig    = isrLib.ISR_BADP,
                    stageName   = 'lsst.ip.isr.maskbadpixels'):
                  
    raise RuntimeError, "Do not call me; use MaskBadPixelsDef.  Talk to RHL if you disagree"

    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    # mask bad pixels
    mi      = exposure.getMaskedImage()
    mask    = mi.getMask()
    bitmask = mask.getPlaneBitMask(maskName)
    afwDetection.setMaskFromFootprintList(mask, fpList, bitmask)

    if interpolate:
        # and interpolate over them
        defaultFwhm = policy.get('defaultFwhm')
        psf = algorithms.createPSF('DoubleGaussian', 0, 0, defaultFwhm/(2*math.sqrt(2*math.log(2))))
        for fp in fpList:
            defect = afwDetection.Defect(fp.getBBox())
            algorithms.interpolateOverDefects(mi, psf, defect)

        stageSummary = 'with interpolation'
    else:
        stageSummary = 'without interpolation'

        
    # common outputs
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))


def MaskBadPixelsDef(exposure, policy, defectList,
                     interpolate = True,
                     maskName    = 'BAD',
                     stageSig    = isrLib.ISR_BADP,
                     stageName   = 'lsst.ip.isr.maskbadpixels'):
                  
    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    # mask bad pixels
    mi      = exposure.getMaskedImage()
    mask    = mi.getMask()
    bitmask = mask.getPlaneBitMask(maskName)
    
    for defect in defectList:
        bbox = defect.getBBox()
        afwDetection.setMaskFromFootprint(mask, afwDetection.Footprint(bbox), bitmask)    

    if interpolate:
        # and interpolate over them
        defaultFwhm = policy.get('defaultFwhm')
        psf = algorithms.createPSF('DoubleGaussian', 0, 0, defaultFwhm/(2*math.sqrt(2*math.log(2))))

        fallbackValue = afwMath.makeStatistics(mi.getImage(), afwMath.MEANCLIP).getValue()

        algorithms.interpolateOverDefects(mi, psf, defectList, fallbackValue)

        stageSummary = 'with interpolation'
    else:
        stageSummary = 'without interpolation'

        
    # common outputs
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))

#
### STAGE : Linearization
#

def LookupTableFromPolicy(tablePolicy,
                          stageName = 'lsst.ip.isr.lookuptablefrompolicy'):
    tableType   = tablePolicy.getString('type')
    tableLength = tablePolicy.getInt('length')
    tableValues = tablePolicy.getArray('value')
    assert len(tableValues) == tableLength
    tableValues = afwMath.vectorD(tableValues)

    if tableType == 'Replace':
        lookupTable = isrLib.LookupTableReplaceF(tableValues)
    elif tableType == 'Multiplicative':
        lookupTable = isrLib.LookupTableMultiplicativeF(tableValues)
    else:
        pexLog.Trace(stageName, 4, 'Unknown table type : %s' % (tableType))
        return None
    
    return lookupTable


def Linearization(exposure, policy,
                  lookupTable = None,
                  stageSig    = isrLib.ISR_LIN,
                  stageName   = 'lsst.ip.isr.linearization',
                  policyPath  = os.path.join(utils.productDir('ip_isr'), 'pipeline')):

    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    if lookupTable == None:
        lookupTableName = policy.getPolicy('linearizePolicy').getString('lookupTableName')
        lookupTable     = LookupTableFromPolicy(os.path.join(policyPath, lookupTableName))
    else:
        lookupTableName = 'provided to ipIsr.Linearization'

    gain = metadata.get('gain')
    mi   = exposure.getMaskedImage()
    lookupTable.apply(mi, gain)
    
    # common outputs
    stageSummary = 'using table %s' % (lookupTableName)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))

#
### STAGE : Cosmic Ray Rejection
#

def CrRejection(exposure, policy,
                stageSig      = isrLib.ISR_CRREJ,
                stageName     = 'lsst.ip.isr.crreject',
                subBackground = True):
    
    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    crPolicy    = policy.getPolicy('crRejectionPolicy')
    # gain is LSST norm
    gainKeyword = 'gain' 
    gain        = metadata.get(gainKeyword)
    # needed for CR
    crPolicy.set('e_per_dn', gain)

    mi = exposure.getMaskedImage()
    if subBackground:
        # how much of this do we put in policy?
        bctrl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
        bctrl.setNxSample(max(2, int(mi.getWidth()/256) + 1))
        bctrl.setNySample(max(2, int(mi.getHeight()/256) + 1))
        bctrl.sctrl.setNumSigmaClip(3)
        bctrl.sctrl.setNumIter(3)
        
        im      = mi.getImage()
        backobj = afwMath.makeBackground(im, bctrl)
        im     -= backobj.getImageF()

    # NOTE - this background issue needs to be resolved
    bg = 0.
    
    defaultFwhm = policy.get('defaultFwhm')
    psf         = algorithms.createPSF('DoubleGaussian', 0, 0, defaultFwhm/(2*math.sqrt(2*math.log(2))))
    crs         = algorithms.findCosmicRays(mi, psf, bg, crPolicy, False)    
    
    if subBackground:
        im     += backobj.getImageF() 
    
    # common outputs
    stageSummary = 'with background subtraction = %s; found %d CRs' % (str(subBackground),
                                                                       len(crs))
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))

#
### STAGE : Saturation correction
#

def SaturationCorrection(exposure, policy,
                         interpolate = True,
                         maskName    = 'SAT',
                         stageSig    = isrLib.ISR_SAT,
                         stageName   = 'lsst.ip.isr.saturationcorrection'):

    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    try:
        satKeyword = policy.getPolicy('saturationPolicy').getString('saturationKeyword')
        saturation = metadata.get(satKeyword)
    except:
        saturation = policy.getPolicy('saturationPolicy').get('defaultSaturation')
        pexLog.Trace(stageName, 4, 'Unable to read %s, using default saturation of %s' % (satKeyword, saturation))    
        
    mi         = exposure.getMaskedImage()
    mask       = mi.getMask()
    bitmask    = mask.getPlaneBitMask(maskName)

    # find saturated regions
    thresh     = afwDetection.Threshold(saturation)
    ds         = afwDetection.DetectionSetF(mi, thresh)
    fpList     = ds.getFootprints()
    # we will turn them into defects for interpolating
    defectList = algorithms.DefectListT()
    
    # grow them
    growSaturated = policy.getPolicy('saturationPolicy').getInt('growSaturated')
    for fp in fpList:
        # if "True", growing requires a convolution
        # if "False", its faster
        fpGrow = afwDetection.growFootprint(fp, growSaturated, False)
        afwDetection.setMaskFromFootprint(mask, fpGrow, bitmask)

        if interpolate:
            defect = algorithms.Defect(fpGrow.getBBox())
            defectList.push_back(defect)

    # interpolate over them
    if interpolate:
        mask.addMaskPlane('INTERP')
        defaultFwhm   = policy.get('defaultFwhm')
        psf = algorithms.createPSF('DoubleGaussian', 0, 0, defaultFwhm/(2*math.sqrt(2*math.log(2))))
        algorithms.interpolateOverDefects(mi, psf, defectList)
    
    # common outputs
    stageSummary = 'using %s=%.2f' % (satKeyword, saturation)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))

#
### STAGE : Bias / dark correction
#

def BiasCorrection(exposure, bias, policy,
                   stageSig  = isrLib.ISR_BIAS,
                   stageName = 'lsst.ip.isr.biascorrection'):

    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    bmetadata         = bias.getMetadata()
#    filename          = bmetadata.getString('filename')
# TEMPORARY 
    filename = 'Zero'
    mi  = exposure.getMaskedImage()

    bmi = bias.getMaskedImage()
    mi -= bmi

    # common outputs
    stageSummary = 'using %s' % (filename)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))
#
# should be removing BIASSEC here.   Where did we use it?  Did we ignore it?
#  

def DarkCorrection(exposure, dark, policy,
                   stageSig  = isrLib.ISR_DARK,
                   stageName = 'lsst.ip.isr.darkcorrection'):
    
    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    dmetadata         = dark.getMetadata()
    filename          = dmetadata.getString('filename')

    # expTime is LSST norm
    scalingKeyword    = 'expTime'
    expscaling        = metadata.get(scalingKeyword)
    darkscaling       = dmetadata.get(scalingKeyword)
    scale             = expscaling / darkscaling

    mi  = exposure.getMaskedImage()
    mi.scaledMinus(scale, dark.getMaskedImage())

    # common outputs
    stageSummary = 'using %s with scale=%.2f' % (filename, scale)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))
    
                   
#
### STAGE : Flat / illum correction
#

def FlatCorrection(exposure, flat, policy,
                   stageSig  = isrLib.ISR_DFLAT,
                   stageName = 'lsst.ip.isr.flatcorrection'):

    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    fmetadata         = flat.getMetadata()
#    filename          = fmetadata.getString('filename')
# TEMPORARY 
    filename = 'Flat'

    scalingKeyword    = policy.getPolicy('flatPolicy').getString('flatScaleKeyword') # e.g. MEAN
    flatscaling = 1.0
    if fmetadata.exists(scalingKeyword):
        flatscaling = fmetadata.get(scalingKeyword)
    else:
        # Figure it out from the data
        scalingKeyword = policy.getPolicy('flatPolicy').getString('flatScaleFitType')
        if scalingKeyword == 'MEAN':
            flatscaling = afwMath.makeStatistics(flat.getMaskedImage().getImage(), afwMath.MEAN).getValue(afwMath.MEAN)
        elif scalingKeyword == 'MEDIAN':
            flatscaling = afwMath.makeStatistics(flat.getMaskedImage().getImage(), afwMath.MEDIAN).getValue(afwMath.MEDIAN)
        else:
            raise pexExcept.LsstException, '%s : %s not implemented' % (stageName, flatscalingType)            

    mi   = exposure.getMaskedImage()
    fmi  = flat.getMaskedImage()
    mi.scaledDivides(1./flatscaling, fmi)
    
    # common outputs
    stageSummary = 'using %s with %s=%.2f' % (filename, scalingKeyword, flatscaling)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))


def IlluminationCorrection(exposure, illum, policy,
                           stageSig  = isrLib.ISR_ILLUM,
                           stageName = 'lsst.ip.isr.illuminationcorrection'):

    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    imetadata         = illum.getMetadata()
    filename          = imetadata.getString('filename')

    scalingKeyword    = policy.getPolicy('illuminationPolicy').getString('illumScaleKeyword')
    illumscaling      = imetadata.get(scalingKeyword)

    mi   = exposure.getMaskedImage()
    mi.scaledDivides(1./illumscaling, illum.getMaskedImage())
    
    # common outputs
    stageSummary = 'using %s with scale=%.2f' % (filename, illumscaling)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))


#
### STAGE : Trim / overscan correction
#

# Now implemented in C++
#
#def BBoxFromDatasec(string,
#                    stageName = 'lsst.ip.isr.bboxfromdatasec'):
#    
#    c = re.compile('^\[(\d+):(\d+),(\d+):(\d+)\]$')
#    m = c.match(string)
#    if m:
#        startCol, endCol, startRow, endRow = m.groups()
#        # Beware the FITS convention
#        startCol -= floor((1 + 0.5 - afwImage.PixelZeroPos))
#        startRow -= floor((1 + 0.5 - afwImage.PixelZeroPos))
#    else:
#        raise pexExcept.LsstException, '%s : Cannot parse %s' % (stageName, string)
#
#    bbox = afwImage.BBox(afwImage.PointI(startCol, startRow),
#                         endCol-startCol,
#                         endRow-startRow)
#    return bbox

def TrimNew(exposure, policy,
            stageSig  = isrLib.ISR_TRIM,
            stageName = 'lsst.ip.isr.trim'):
    """
    This returns a new Exposure that is a subsection of the input exposure.
    
    NOTE : do we need to deal with the WCS in any way, shape, or form?
    """
    
    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    trimsec = None
    if metadata.exists('trimsec'):
        trimsecKeyword = 'trimsec'
        trimsec = metadata.getString(trimsecKeyword)
    else:
        trimsecKeyword = policy.getPolicy('trimPolicy').getString('trimsecKeyword')
        if metadata.exists(trimsecKeyword):
            trimsec = metadata.getString(trimsecKeyword)
    if trimsec == None:
        raise pexExcept.LsstException, '%s : cannot find trimsec' % (stageName)        
    trimsecBBox = isrLib.BBoxFromDatasec(trimsec)

    #
    # In this case, we do NOT want the X0, Y0 of exposure to be adjusted, so save it and restore it
    #
    saveXY0 = exposure.getMaskedImage().getXY0()

    trimmedExposure = afwImage.ExposureF(exposure, trimsecBBox)

    trimmedExposure.getMaskedImage().setXY0(saveXY0)

    # remove trimsec from metadata
    metadata.remove(trimsecKeyword)

    # shift datasec to be relative to trimsec, which defines new origin
    # remove biassec

    datasec = None
    if metadata.exists('datasec'):
        datasecKeyword = 'datasec'
        datasec = metadata.getString(datasecKeyword)
    else:
        datasecKeyword = policy.getPolicy('trimPolicy').getString('datasecKeyword')
        if metadata.exists(datasecKeyword):
            datasec = metadata.getString(datasecKeyword)
    if datasec == None:
        raise pexExcept.LsstException, '%s : cannot find datasec' % (stageName)        
    datasecBBox = isrLib.BBoxFromDatasec(datasec)

    datasecBBox.shift(-trimsecBBox.getX0(), -trimsecBBox.getY0())

    datasec = isrLib.DatasecFromBBox(datasecBBox)

    metadata.set(datasecKeyword, datasec)

    # Any changes needed for wcs??
    
    # common outputs
    stageSummary = 'using trimsec %s' % (trimsec)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))

    return trimmedExposure


def OverscanCorrection(exposure, policy,
                       stageSig  = isrLib.ISR_OSCAN,
                       stageName = 'lsst.ip.isr.overscancorrection'):
    """
    This returns a new Exposure that is a subsection of the input exposure.
    
    NOTE : do we need to deal with the WCS in any way, shape, or form?
    """

    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    mi = exposure.getMaskedImage()

    overscan = None
    if metadata.exists('overscan'):
        overscanKeyword = 'overscan'
        overscan = metadata.getString(overscanKeyword)
    else:
        overscanKeyword = policy.getPolicy('overscanPolicy').getString('overscanKeyword')
        if metadata.exists(overscanKeyword):
            overscan = metadata.getString(overscanKeyword) 
    if overscan == None:
        raise pexExcept.LsstException, '%s : cannot find overscan region' % (stageName)     
    overscanBBox    = isrLib.BBoxFromDatasec(overscan)

    # if "True", do a deep copy
    overscanData    = afwImage.ImageF(exposure.getMaskedImage().getImage(), overscanBBox, False)

    # what type of overscan modeling?
    overscanFitType = policy.getPolicy('overscanPolicy').getString('overscanFitType')
    if overscanFitType == 'MEAN':
        offset = afwMath.makeStatistics(overscanData, afwMath.MEAN).getValue(afwMath.MEAN)
        mi    -= offset
    elif overscanFitType == 'MEDIAN':
        offset = afwMath.makeStatistics(overscanData, afwMath.MEDIAN).getValue(afwMath.MEDIAN)
        mi    -= offset
    elif overscanFitType == 'POLY':
        polyOrder = policy.getPolicy('overscanPolicy').getInt('polyOrder')
        raise pexExcept.LsstException, '%s : %s not implemented' % (stageName, overscanFitType)
    else:
        raise pexExcept.LsstException, '%s : %s an invalid overscan type' % (stageName, overscanFitType)

    # remove overscan from metadata
    exposure.getMetadata().remove(overscanKeyword)

    # common outputs
    stageSummary = 'using overscan section %s with %s=%f' % (overscan, overscanFitType, offset)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))

#
### STAGE : Fringe correction
#

def FringeCorrection(exposure, fringe, policy,
                     stageSig  = isrLib.ISR_FRING,
                     stageName = 'lsst.ip.isr.fringecorrection'):

    fringeSkyKeyword   = policy.getPolicy('fringePolicy').getString('fringeSkyKeyword')
    fringeScaleKeyword = policy.getPolicy('fringePolicy').getString('fringeScaleKeyword')
    
    raise pexExcept.LsstException, '%s not implemented' % (stageName)

#
### STAGE : Pupil correction
#

def PupilCorrection(exposure, fringe, policy,
                    stageSig  = isrLib.ISR_PUPIL,
                    stageName = 'lsst.ip.isr.pupilcorrection'):

    raise pexExcept.LsstException, '%s not implemented' % (stageName)

def keyIsSet(key, policy):
    if not policy.exists(key):
        return False
    
    return(policy.getBool(key))
