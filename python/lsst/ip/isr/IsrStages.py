import time, re

import lsst.afw.detection   as afwDetection
import lsst.afw.image       as afwImage
import lsst.afw.math        as afwMath
import lsst.meas.algorithms as algorithms
import lsst.pex.logging     as pexLog
import lsst.pex.exceptions  as pexExcept

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


def MaskBadPixels(exposure, policy, fpList,
                  interpolate = True,
                  maskName    = 'BAD',
                  stageSig    = ipIsr.ISR_BADP,
                  stageName   = 'lsst.ip.isr.maskbadpixels'):
                  
    # common input test
    metadata   = exposure.getMetadata()
    if chunkMetadata.exists(stageSig):
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
        psf = algorithms.createPSF('DGPSF', 0, defaultFwhm/(2*sqrt(2*log(2))))
        mask.addMaskPlane('INTERP')
        for fp in fpList:
            defect = afwDetection.Defect(fp.getBbox())
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
        lookupTable = ipIsr.LookupTableReplaceF(tableValues)
    elif tableType == 'Multiplicative':
        lookupTable = ipIsr.LookupTableMultiplicativeF(tableValues)
    else:
        pexLog.Trace(stageName, 4, 'Unknown table type : %s' % (tableType))
        return None
    
    return lookupTable


def Linearization(exposure, policy,
                  lookupTable = None,
                  stageSig    = ipIsr.ISR_LIN,
                  stageName   = 'lsst.ip.isr.linearization'):

    # common input test
    metadata   = exposure.getMetadata()
    if chunkMetadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    if lookupTable == None:
        lookupTableName = policy.getString('lookupTableName')
        lookupTable = LookupTableFromPolicy(lookupTableName)

    mi = exposure.getMaskedImage()
    lookupTable.apply(mi)
    
    # common outputs
    stageSummary = 'using table %s' % (satKeyword, saturation))    
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))

#
### STAGE : Saturation correction
#

def SaturationCorrection(exposure, policy,
                         maskName  = 'SAT',
                         stageSig  = ipIsr.ISR_SAT,
                         stageName = 'lsst.ip.isr.saturationcorrection'):

    # common input test
    metadata   = exposure.getMetadata()
    if chunkMetadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return
    
    satKeyword = policy.getString('saturationKeyword')
    saturation = metadata.getDouble(satKeyword)

    mi         = exposure.getMaskedImage()
    mask       = mi.getMask()
    bitmask    = mask.getPlaneBitMask(maskName)
    
    # find saturated regions
    thresh    = afwDetection.Threshold(saturation)
    ds        = afwDetection.DetectionSetF(mi, thresh)
    fpList    = ds.getFootprints()
    
    # grow them
    growSaturated = satPolicy.getInt('growSaturated')
    for fp in fpList:
        fpGrow = afwDetection.growFootprint(fp, growSaturated)
        afwDetection.setMaskFromFootprint(mask, fpGrow, bitmask)

    # interpolate over them
    defaultFwhm   = satPolicy.getDouble('defaultFwhm')
    psf = algorithms.createPSF('DGPSF', 0, defaultFwhm/(2*sqrt(2*log(2))))
    mask.addMaskPlane('INTERP')
    for fp in fpList:
        defect = afwDetection.Defect(fp.getBbox())
        algorithms.interpolateOverDefects(mi, psf, defect)
    
    # common outputs
    stageSummary = 'using %s=%.2f' % (satKeyword, saturation))    
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))

#
### STAGE : Bias / dark correction
#

def BiasCorrection(exposure, bias, policy,
                   stageSig  = ipIsr.ISR_BIAS,
                   stageName = 'lsst.ip.isr.biascorrection'):

    # common input test
    metadata   = exposure.getMetadata()
    if chunkMetadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    bmetadata         = bias.getMetadata()
    filenameKeyword   = policy.getString('filenameKeyword')
    meanCountsKeyword = policy.getString('meanCountsKeyword')
    filename          = bmetadata(filenameKeyword)
    meanCounts        = bmetadata(meanCountsKeyword)

    mi  = exposure.getMaskedImage()
    bmi = bias.getMaskedImage()
    mi -= bmi

    # common outputs
    stageSummary = 'using %s with mean=%.2f' % (filename, meanCounts)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))
    

def DarkCorrection(exposure, dark, policy,
                   stageSig  = ipIsr.ISR_DARK,
                   stageName = 'lsst.ip.isr.darkcorrection'):
    
    # common input test
    metadata   = exposure.getMetadata()
    if chunkMetadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    dmetadata         = dark.getMetadata()
    filenameKeyword   = policy.getString('filenameKeyword')
    filename          = dmetadata(filenameKeyword)

    scalingKeyword    = policy.getString('darkScaleKeyword') # e.g. EXPTIME
    expscaling        = metadata.getDouble(scalingKeyword)
    darkscaling       = dmetadata.getDouble(scalingKeyword)
    scale             = expscaling / darkscaling

    #mi   = exposure.getMaskedImage()
    #bmi  = afwImage.MaskedImageF(bias.getMaskedImage(), True)
    #bmi *= scale
    mi.scaledMinus(scale, bias.getMaskedImage())

    # common outputs
    stageSummary = 'using %s with scale=%.2f' % (filename, scale)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))
    
                   
#
### STAGE : Flat / illum correction
#

def FlatCorrection(exposure, flat, policy,
                   stageSig  = ipIsr.ISR_FLAT,
                   stageName = 'lsst.ip.isr.flatcorrection'):

    # common input test
    metadata   = exposure.getMetadata()
    if chunkMetadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    fmetadata         = flat.getMetadata()
    filenameKeyword   = policy.getString('filenameKeyword')
    filename          = fmetadata(filenameKeyword)

    scalingKeyword    = policy.getString('flatScaleKeyword') # e.g. MEAN
    flatscaling       = dmetadata.getDouble(scalingKeyword)

    mi   = exposure.getMaskedImage()
    #fmi  = afwImage.MaskedImageF(flat.getMaskedImage(), True)
    #if flatscaling != 1.:
    #    fmi /= flatscaling
    #mi  /= fmi
    mi.scaledDivides(1./flatscaling, flat.getMaskedImage())
    
    # common outputs
    stageSummary = 'using %s with scale=%.2f' % (filename, flatscaling)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))


def IlluminationCorrection(exposure, illum, policy,
                           stageSig  = ipIsr.ISR_ILLUM,
                           stageName = 'lsst.ip.isr.illuminationcorrection'):

    # common input test
    metadata   = exposure.getMetadata()
    if chunkMetadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    imetadata         = illum.getMetadata()
    filenameKeyword   = policy.getString('filenameKeyword')
    filename          = imetadata(filenameKeyword)

    scalingKeyword    = policy.getString('illumScaleKeyword') # e.g. MEAN
    illumscaling      = dmetadata.getDouble(scalingKeyword)

    mi   = exposure.getMaskedImage()
    #imi  = afwImage.MaskedImageF(illum.getMaskedImage(), True)
    #if illumscaling != 1.:
    #    imi /= illumscaling
    #mi  /= imi
    mi.scaledDivides(1./illumscaling, illum.getMaskedImage())
    
    # common outputs
    stageSummary = 'using %s with scale=%.2f' % (filename, flatscaling)
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))


#
### STAGE : Trim / overscan correction
#

def BboxFromDatasec(string,
                    stageName = 'lsst.ip.isr.bboxfromdatasec'):
    
    c = re.compile('^(\d)\:(\d)\,(\d)\:(\d)$')
    m = c.match(string)
    if m:
        startCol, endCol, startRow, endRow = m.groups()
        # Beware the FITS convention
        startCol -= floor((1 + 0.5 - afwImage.PixelZeroPos))
        startRow -= floor((1 + 0.5 - afwImage.PixelZeroPos))
    else:
        raise pexExcept.LsstException, '%s : Cannot parse %s' % (stageName, string)

    bbox = afwImage.BBox(afwImage.PointI(startCol, startRow),
                         endCol-startCol,
                         endRow-startRow)
    return bbox


def TrimNew(exposure, policy,
            stageSig  = ipIsr.ISR_TRIM,
            stageName = 'lsst.ip.isr.trim'):
    """
    This returns a new Exposure that is a subsection of the input exposure.
    
    NOTE : do we need to deal with the WCS in any way, shape, or form?
    """
    
    # common input test
    metadata   = exposure.getMetadata()
    if chunkMetadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    trimsecKeyword  = policy.getString('trimsecKeyword')
    trimsec         = metadata.getString(trimsecKeyword)
    trimsecBbox     = BboxFromDatasec(trimsec)

    # if "True", do a deep copy
    trimmedExposure = afwImage.ExposureF(exposure, trimsecBbox, False)

    # common outputs
    stageSummary = 'using trimsec %s' % (trimsec))    
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))

    return trimmedExposure


def OverscanCorrection(exposure, policy,
                       stageSig  = ipIsr.ISR_OSCAN,
                       stageName = 'lsst.ip.isr.overscancorrection'):
    """
    This returns a new Exposure that is a subsection of the input exposure.
    
    NOTE : do we need to deal with the WCS in any way, shape, or form?
    """

    # common input test
    metadata   = exposure.getMetadata()
    if chunkMetadata.exists(stageSig):
        pexLog.Trace(stageName, 4, '%s has already been run' % (stageSig))
        return

    mi = exposure.getMaskedImage()
    
    overscanKeyword = policy.getString('overscanKeyword')
    overscan        = metadata.getString(overscanKeyword)
    overscanBbox    = BboxFromDatasec(overscan)

    # if "True", do a deep copy
    overscanData    = afwImage.MaskedImage(exposure.getMaskedImage(), overscanBbox, False)

    # what type of overscan modeling?
    overscanFitType = policy.getString('overscanFitType')
    if overscanFitType == 'MEAN':
        mean   = afwMath.makeStatistics(overscanData, afwMath.MEAN).getValue(afwMath.MEAN)
        mi    -= mean
    elif overscanFitType == 'MEDIAN':
        median = afwMath.makeStatistics(overscanData, afwMath.MEDIAN).getValue(afwMath.MEDIAN)
        mi    -= median
    else:
        raise pexExcept.LsstException, '%s : %s not implemented' % (stageName, overscanFitType)

    # common outputs
    stageSummary = 'using overscan section %s' % (overscan))    
    pexLog.Trace(stageName, 4, '%s %s' % (stageSig, stageSummary))    
    metadata.setString(stageSig, '%s; %s' % (stageSummary, time.asctime()))

#
### STAGE : Fringe correction
#

def FringeCorrection(exposure, fringe, policy,
                     stageSig  = ipIsr.ISR_FRING,
                     stageName = 'lsst.ip.isr.fringecorrection'):

    raise pexExcept.LsstException, '%s not implemented' % (stageName)
