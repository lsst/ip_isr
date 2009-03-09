import time, re, os, math
import eups
import lsst.afw.detection   as afwDetection
import lsst.afw.image       as afwImage
import lsst.afw.math        as afwMath
import lsst.meas.algorithms as algorithms
import lsst.pex.logging     as pexLog
import lsst.pex.exceptions  as pexExcept

# relative imports
import isrLib

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
        defaultFwhm = policy.getDouble('defaultFwhm')
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
        afwDetection.setMaskFromFootprint(mask, defect.getFootprint(), bitmask)    

    if interpolate:
        # and interpolate over them
        defaultFwhm = policy.getDouble('defaultFwhm')
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
        lookupTable = isrLib.LookupTableReplaceI(tableValues)
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
                  policyPath  = os.path.join(eups.productDir('ip_isr'), 'pipeline')):

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

    mi = exposure.getMaskedImage()
    lookupTable.apply(mi)
    
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
                subBackground = False):
    
    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    crPolicy    = policy.getPolicy('crRejectionPolicy')
    gainKeyword = crPolicy.getString('gainKeyword')
    gain        = metadata.getDouble(gainKeyword)
    # needed for CR
    crPolicy.set('e_per_dn', gain)

    mi = exposure.getMaskedImage()
    import pdb
    pdb.set_trace()
    if subBackground:
        # how much of this do we put in policy?
        bctrl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
        bctrl.setNxSample(int(mi.getWidth()/256) + 1)
        bctrl.setNySample(int(mi.getHeight()/256) + 1)
        bctrl.sctrl.setNumSigmaClip(3.0)
        bctrl.sctrl.setNumIter(2)
        
        im      = mi.getImage()
        backobj = afwMath.makeBackground(im, bctrl)
        im     -= backobj.getImageF()

    # NOTE - this background issue needs to be resolved
    bg = 0.
    
    defaultFwhm = policy.getDouble('defaultFwhm')
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
        saturation = metadata.getDouble(satKeyword)
    except:
        saturation = policy.getPolicy('saturationPolicy').getDouble('defaultSaturation')
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
        fpGrow = afwDetection.growFootprint(fp, growSaturated)
        afwDetection.setMaskFromFootprint(mask, fpGrow, bitmask)

        if interpolate:
            defect = algorithms.Defect(fpGrow.getBBox())
            defectList.push_back(defect)

    # interpolate over them
    if interpolate:
        mask.addMaskPlane('INTERP')
        defaultFwhm   = policy.getDouble('defaultFwhm')
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
    
    filenameKeyword   = policy.getString('filenameKeyword')
    filename          = bmetadata.getString(filenameKeyword)

    meanCountsKeyword = policy.getPolicy('biasPolicy').getString('meanCountsKeyword')
    meanCounts        = bmetadata.getDouble(meanCountsKeyword)

    mi  = exposure.getMaskedImage()
    bmi = bias.getMaskedImage()
    mi -= bmi

    # common outputs
    stageSummary = 'using %s with mean=%.2f' % (filename, meanCounts)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))
    

def DarkCorrection(exposure, dark, policy,
                   stageSig  = isrLib.ISR_DARK,
                   stageName = 'lsst.ip.isr.darkcorrection'):
    
    # common input test
    metadata   = exposure.getMetadata()
    if metadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    dmetadata         = dark.getMetadata()
    filenameKeyword   = policy.getString('filenameKeyword')
    filename          = dmetadata.getString(filenameKeyword)

    scalingKeyword    = policy.getPolicy('darkPolicy').getString('darkScaleKeyword') # e.g. EXPTIME
    expscaling        = metadata.getDouble(scalingKeyword)
    darkscaling       = dmetadata.getDouble(scalingKeyword)
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
    filenameKeyword   = policy.getString('filenameKeyword')
    filename          = fmetadata.getString(filenameKeyword)

    scalingKeyword    = policy.getPolicy('flatPolicy').getString('flatScaleKeyword') # e.g. MEAN
    flatscaling       = fmetadata.getDouble(scalingKeyword)

    mi   = exposure.getMaskedImage()
    mi.scaledDivides(1./flatscaling, flat.getMaskedImage())
    
    # common outputs
    stageSummary = 'using %s with scale=%.2f' % (filename, flatscaling)
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
    filenameKeyword   = policy.getString('filenameKeyword')
    filename          = imetadata.getString(filenameKeyword)

    scalingKeyword    = policy.getPolicy('illuminationPolicy').getString('illumScaleKeyword')
    illumscaling      = imetadata.getDouble(scalingKeyword)

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

    trimsecKeyword  = policy.getPolicy('trimPolicy').getString('trimsecKeyword')
    trimsec         = metadata.getString(trimsecKeyword)
    trimsecBBox     = isrLib.BBoxFromDatasec(trimsec)

    # if "True", do a deep copy
    trimmedExposure = afwImage.ExposureF(exposure, trimsecBBox, False)
    llc = trimsecBBox.getLLC()
    trimmedExposure.setXY0(llc)

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
    
    overscanKeyword = policy.getPolicy('overscanPolicy').getString('overscanKeyword')
    overscan        = metadata.getString(overscanKeyword)
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
