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

import time, os, math
import lsst.utils           as utils
import lsst.afw.detection   as afwDetection
import lsst.afw.image       as afwImage
import lsst.afw.geom        as afwGeom
import lsst.afw.cameraGeom  as cameraGeom
import lsst.afw.math        as afwMath
import lsst.afw.coord       as afwCoord
import lsst.meas.algorithms as algorithms
import lsst.pex.logging     as pexLog
import lsst.pex.policy      as pexPolicy
import lsst.pex.exceptions  as pexExcept
import lsst.daf.base        as dafBase

# relative imports
import isrLib
class Bbox(object):
    def __init__(self, x0, y0, dx, dy):
        self.x0 = x0
	self.y0 = y0
	self.width = dx
	self.height = dy
	self.x1 = self.x0 + self.width - 1
	self.y1 = self.y0 + self.height - 1
    def shift(self, dx, dy):
	self.x0 += dx
	self.x1 += dx
	self.y0 += dy
	self.y1 += dy
    def __str__(self):
        return "(%i,%i) -- (%i,%i)"%(self.x0,self.y0,self.x1,self.y1)

def createPsf(fwhm):
    """Make a PSF"""
    ksize = 4*int(fwhm) + 1
    return afwDetection.createPsf('DoubleGaussian', ksize, ksize, fwhm/(2*math.sqrt(2*math.log(2))))

def calcEffectiveGain(exposure):
    mi = exposure.getMaskedImage()
    im = afwImage.ImageF(mi.getImage(), True)
    var = mi.getVariance()
    im /= var
    medgain = afwMath.makeStatistics(im, afwMath.MEDIAN).getValue()
    meangain = afwMath.makeStatistics(im, afwMath.MEANCLIP).getValue()
    return medgain, meangain

def convertImageForIsr(exposure):
    if not isinstance(exposure, afwImage.ExposureU):
        raise Exception("ipIsr.convertImageForIsr: Expecting Uint16 image. Got\
                %s."%(exposure.__repr__()))

    newexposure = exposure.convertF()
    amp = cameraGeom.cast_Amp(exposure.getDetector())
    mi = newexposure.getMaskedImage()
    var = afwImage.ImageF(mi.getDimensions())
    mask = afwImage.MaskU(newexposure.getWidth(), newexposure.getHeight())
    mask.set(0)
    newexposure.setMaskedImage(afwImage.MaskedImageF(mi.getImage(), mask, var))
    return newexposure

def calculateSdqaCcdRatings(exposure):
    metrics = {}
    metrics['nSaturatePix'] = 0
    metrics['nBadCalibPix'] = 0
    metrics['imageClipMean4Sig3Pass'] = None
    metrics['imageSigma'] = None
    metrics['imageMedian'] = None
    metrics['imageMin'] = None
    metrics['imageMax'] = None
    metadata = exposure.getMetadata()
    mi = exposure.getMaskedImage()
    mask = mi.getMask()
    badbitmask = mask.getPlaneBitMask('BAD')
    satbitmask = mask.getPlaneBitMask('SAT')
    intrpbitmask = mask.getPlaneBitMask('INTRP')
    sctrl = afwMath.StatisticsControl()
    sctrl.setNumIter(4)
    sctrl.setNumSigmaClip(3)
    satmask = afwImage.MaskU(mask, True)
    badmask = afwImage.MaskU(mask, True)
    satmask &= satbitmask
    badmask &= badbitmask
    satmaskim = afwImage.ImageU(satmask.getDimensions())
    satmaskim <<= satmask
    badmaskim = afwImage.ImageU(badmask.getDimensions())
    badmaskim <<= badmask
    thresh = afwDetection.Threshold(0.5)
    fs = afwDetection.makeFootprintSet(satmaskim, thresh)
    for f in fs.getFootprints():
        metrics['nSaturatePix'] += f.getNpix()
    fs = afwDetection.makeFootprintSet(badmaskim, thresh)
    for f in fs.getFootprints():
        metrics['nBadCalibPix'] += f.getNpix()
    metrics['imageClipMean4Sig3Pass'] = afwMath.makeStatistics(mi, afwMath.MEANCLIP, sctrl).getValue()
    metrics['imageSigma'] = afwMath.makeStatistics(mi, afwMath.STDEVCLIP, sctrl).getValue()
    metrics['imageMedian'] = afwMath.makeStatistics(mi, afwMath.MEDIAN, sctrl).getValue()
    metrics['imageMin'] = afwMath.makeStatistics(mi, afwMath.MIN, sctrl).getValue()
    metrics['imageMax'] = afwMath.makeStatistics(mi, afwMath.MAX, sctrl).getValue()
    for k in metrics.keys():
        metadata.set(k, metrics[k])

def calculateSdqaAmpRatings(exposure, biasBBox, dataBBox):
    metrics = {}
    metrics['nSaturatePix'] = 0
    metrics['overscanMean'] = None
    metrics['overscanStdDev'] = None
    metrics['overscanMedian'] = None
    metrics['overscanMin'] = None
    metrics['overscanMax'] = None
    mi = exposure.getMaskedImage()
    metadata = exposure.getMetadata()
    trimmi = afwImage.MaskedImageF(mi, dataBBox)
    biasmi = afwImage.MaskedImageF(mi, biasBBox)
    mask = mi.getMask()
    satbitmask = mask.getPlaneBitMask('SAT')
    satmask = trimmi.getMask()
    satmask &= satbitmask
    satmaskim = afwImage.ImageU(satmask.getDimensions())
    satmaskim <<= satmask
    thresh = afwDetection.Threshold(0.5)
    fs = afwDetection.makeFootprintSet(satmaskim, thresh)
    for f in fs.getFootprints():
        metrics['nSaturatePix'] += f.getNpix()
    sctl = afwMath.StatisticsControl()
    metrics['overscanMean'] = afwMath.makeStatistics(biasmi, afwMath.MEAN, sctl).getValue() 
    metrics['overscanStdDev'] = afwMath.makeStatistics(biasmi, afwMath.STDEV, sctl).getValue() 
    metrics['overscanMedian'] = afwMath.makeStatistics(biasmi, afwMath.MEDIAN, sctl).getValue() 
    metrics['overscanMin'] = afwMath.makeStatistics(biasmi, afwMath.MIN, sctl).getValue() 
    metrics['overscanMax'] = afwMath.makeStatistics(biasmi, afwMath.MAX, sctl).getValue() 
    for k in metrics.keys():
        metadata.set(k, metrics[k])

def exposureFromInputData(image, metadata, ampBBox,
                          makeWcs     = True,
                          policy      = None,
                          defaultGain = 1.0,
                          imageSource = afwImage.ImageF):
    methodName = "isr.exposureFromInputData"

    # Generate an empty mask
    mask = afwImage.MaskU(image.getDimensions())
    mask.set(0)

    # Generate a variance from the image pixels and gain
    var  = afwImage.ImageF(image, True)
    #var = afwImage.ImageF(image)
    
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
                pexLog.Trace(methodName, 4, 'Using default gain=%f for %s' % (defaultGain, filename))
                gain = defaultGain
        else:
            pexLog.Trace(methodName, 4, 'Using default gain=%f for %s' % (defaultGain, filename))
            gain = defaultGain
    else:
        pexLog.Trace(methodName, 4, 'Using default gain=%f' % (defaultGain))
        gain = defaultGain
    # Normalize by the gain
    var /= gain

    # makeMaskedImage() will make a MaskedImage with the same type as Image
    mi   = afwImage.makeMaskedImage(image, mask, var)
    mi.setXY0(ampBBox.getX0(), ampBBox.getY0())

    if makeWcs:
        # Extract the Wcs info from the input metadata
        wcs      = afwImage.makeWcs(metadata)
    else:
        wcs      = afwImage.Wcs()
        
    # makeExposure will make an Exposure with the same type as MaskedImage
    exposure = afwImage.makeExposure(mi, wcs)
    exposure.setMetadata(metadata)

    return exposure


def validateCalibration(exposure, calibration, policy):
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


def maskFromDefects(dimensions, fpList):
    # the output LSST Mask image
    mask = afwImage.MaskU(dimensions)
    mask.set(0)
    bitmask = mask.getPlaneBitMask('BAD')

    # set the bits
    afwDetection.setMaskFromFootprintList(mask, fpList, bitmask)

    return mask
    

def defectsFromBoolImage(fitsfile, invert=False):
    # input bad pixel image
    # This assumes an image with ones and zeros
    image = afwImage.ImageF(fitsfile)
    if invert:
        image *= -1
        thresh = afwDetection.Threshold(-0.5)
    else:
        thresh = afwDetection.Threshold(0.5)
   
    # turn into masked image for detection
    mi = afwImage.MaskedImageF(image)

    # find bad regions
    ds = afwDetection.FootprintSetF(mi, thresh)
    fpList = ds.getFootprints()

    return fpList


def maskBadPixelsFp(exposure, policy, fpList,
                    interpolate = True,
                    maskName    = 'BAD'):
                  
    raise RuntimeError, "Do not call me; use MaskBadPixelsDef.  Talk to RHL if you disagree"

    # common input test
    metadata   = exposure.getMetadata()

    # mask bad pixels
    mi      = exposure.getMaskedImage()
    mask    = mi.getMask()
    bitmask = mask.getPlaneBitMask(maskName)
    afwDetection.setMaskFromFootprintList(mask, fpList, bitmask)

    if interpolate:
        # and interpolate over them
        psf = createPsf(policy.get('defaultFwhm'))

        for fp in fpList:
            defect = afwDetection.Defect(fp.getBBox())
            algorithms.interpolateOverDefects(mi, psf, defect)

        
def interpolateDefectList(exposure, defectList, fwhm, fallbackValue=None):
    mi = exposure.getMaskedImage()
    psf = createPsf(fwhm)
    if fallbackValue is None:
        fallbackValue = afwMath.makeStatistics(mi.getImage(), afwMath.MEANCLIP).getValue()
    algorithms.interpolateOverDefects(mi, psf, defectList, fallbackValue)
    

def maskBadPixelsDef(exposure, defectList, fwhm,
                     interpolate = True,
                     maskName    = 'BAD'):
                  
    # mask bad pixels
    mi      = exposure.getMaskedImage()
    mask    = mi.getMask()
    bitmask = mask.getPlaneBitMask(maskName)
    for defect in defectList:
        bbox = defect.getBBox()
        afwDetection.setMaskFromFootprint(mask, afwDetection.Footprint(bbox), bitmask)    

    if interpolate:
        # and interpolate over them
        psf = createPsf(fwhm)
        fallbackValue = afwMath.makeStatistics(mi.getImage(), afwMath.MEANCLIP).getValue()
        algorithms.interpolateOverDefects(mi, psf, defectList, fallbackValue)

        

def lookupTableFromPolicy(tablePolicy):
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
        return None
    
    return lookupTable


def linearization(exposure, lookupTable):

    # common input test
    metadata   = exposure.getMetadata()
    gain = metadata.get('gain')
    mi   = exposure.getMaskedImage()
    lookupTable.apply(mi, gain)
    
#Get rid of this: SIMON
def backgroundSubtraction(exposure, gridsize="32",
        interptype="AKIMA_SPLINE", nsigma=3.0, niter=3.0):

    metadata = exposure.getMetadata()
    mi = exposure.getMaskedImage()
    bctrl = None
    try:
        bctrl = {
            'LINEAR'               :
            afwMath.BackgroundControl(afwMath.Interpolate.LINEAR), 
            'NATURAL_SPLINE'       :
            afwMath.BackgroundControl(afwMath.Interpolate.NATURAL_SPLINE), 
            'CUBIC_SPLINE'         :
            afwMath.BackgroundControl(afwMath.Interpolate.CUBIC_SPLINE), 
            'CUBIC_SPLINE_PERIODIC':
            afwMath.BackgroundControl(afwMath.Interpolate.CUBIC_SPLINE_PERIODIC), 
            'AKIMA_SPLINE'         :
            afwMath.BackgroundControl(afwMath.Interpolate.AKIMA_SPLINE), 
            'AKIMA_SPLINE_PERIODIC':
            afwMath.BackgroundControl(afwMath.Interpolate.CUBIC_SPLINE_PERIODIC) 
        }[interptype]
    except:
        bctrl = afwMath.BackgroundControl(afwMath.Interpolate.AKIMA_SPLINE)
    bctrl.setNxSample(max(2, int(mi.getWidth()/gridsize) + 1))
    bctrl.setNySample(max(2, int(mi.getHeight()/gridsize) + 1))
    bctrl.sctrl.setNumSigmaClip(nsigma)
    bctrl.sctrl.setNumIter(niter)
    im      = mi.getImage()
    backobj = afwMath.makeBackground(im, bctrl)
    im     -= backobj.getImageF()

    
def saturationDetection(exposure, saturation, doMask = True,
                         maskName = 'SAT'):

    mi         = exposure.getMaskedImage()
    if False:
        ds9.mtv(mi, frame=0)
    

    # find saturated regions
    thresh     = afwDetection.Threshold(saturation)
    ds         = afwDetection.makeFootprintSet(mi, thresh)
    fpList     = ds.getFootprints()
    # we will turn them into defects for interpolating
    defectList = algorithms.DefectListT()
    
    # grow them
    bboxes = []
    for fp in fpList:
        if doMask:
            mask       = mi.getMask()
            bitmask    = mask.getPlaneBitMask(maskName)
            afwDetection.setMaskFromFootprint(mask, fp, bitmask)
        for bbox in afwDetection.footprintToBBoxList(fp):
            bboxes.append(Bbox(bbox.getX0(), bbox.getY0(), bbox.getWidth(),
                bbox.getHeight()))
    return bboxes

def defectListFromMask(exposure, growFootprints = 1, maskName = 'SAT'):
    mi = exposure.getMaskedImage()
    mask = mi.getMask()
    satmask = afwImage.MaskU(mask, True)
    satmask &= mask.getPlaneBitMask(maskName)
    thresh = afwDetection.Threshold(0.5)
    maskimg = afwImage.ImageU(satmask.getDimensions())
    maskimg <<= satmask
    ds = afwDetection.makeFootprintSet(maskimg, thresh)
    fpList = ds.getFootprints()
    satDefectList = algorithms.DefectListT()
    for fp in fpList:
        if growFootprints > 0:
            # if "True", growing requires a convolution
            # if "False", its faster
            fpGrow = afwDetection.growFootprint(fp, growFootprints, False)
        else:
            fpGrow = fp
        for bbox in afwDetection.footprintToBBoxList(fpGrow):
            defect = algorithms.Defect(bbox)
            satDefectList.push_back(defect)
    if 'INTRP' not in mask.getMaskPlaneDict().keys():
        mask.addMaskPlane('INTRP')
    return satDefectList

def saturationInterpolation(exposure, fwhm, growFootprints = 1, maskName = 'SAT'):
    mi = exposure.getMaskedImage()
    mask = mi.getMask()
    satmask = afwImage.MaskU(mask, True)
    satmask &= mask.getPlaneBitMask(maskName)
    thresh = afwDetection.Threshold(0.5)
    maskimg = afwImage.ImageU(satmask.getDimensions())
    maskimg <<= satmask
    ds = afwDetection.makeFootprintSet(maskimg, thresh)
    fpList = ds.getFootprints()
    satDefectList = algorithms.DefectListT()
    for fp in fpList:
        if growFootprints > 0:
            # if "True", growing requires a convolution
            # if "False", its faster
            fpGrow = afwDetection.growFootprint(fp, growFootprints, False)
        else:
            fpGrow = fp
        for bbox in afwDetection.footprintToBBoxList(fpGrow):
            defect = algorithms.Defect(bbox)
            satDefectList.push_back(defect)
    if 'INTRP' not in mask.getMaskPlaneDict().keys():
        mask.addMaskPlane('INTRP')
    psf = createPsf(fwhm)

    algorithms.interpolateOverDefects(mi, psf, satDefectList)


def saturationCorrection(exposure, saturation, fwhm, growSaturated = False,
                         interpolate = True,
                         maskName    = 'SAT'):

    mi         = exposure.getMaskedImage()
    mask       = mi.getMask()
    bitmask    = mask.getPlaneBitMask(maskName)
    if False:
        ds9.mtv(mi, frame=0)
    

    # find saturated regions
    thresh     = afwDetection.Threshold(saturation)
    ds         = afwDetection.makeFootprintSet(mi, thresh)
    fpList     = ds.getFootprints()
    # we will turn them into defects for interpolating
    defectList = algorithms.DefectListT()
    
    # grow them
    for fp in fpList:
        # if "True", growing requires a convolution
        # if "False", its faster
        fpGrow = afwDetection.growFootprint(fp, growSaturated, False)
        afwDetection.setMaskFromFootprint(mask, fpGrow, bitmask)

        if interpolate:
            for bbox in afwDetection.footprintToBBoxList(fpGrow):
                defect = algorithms.Defect(bbox)
                defectList.push_back(defect)

    # interpolate over them
    if interpolate:
        mask.addMaskPlane('INTRP')
        psf = createPsf(fwhm)
        algorithms.interpolateOverDefects(mi, psf, defectList)
    

def biasCorrection(exposure, bias):

    mi  = exposure.getMaskedImage()
    bmi = bias.getMaskedImage()
    mi -= bmi

def darkCorrection(exposure, dark, expscaling, darkscaling):
    
    scale             = expscaling / darkscaling
    mi  = exposure.getMaskedImage()
    mi.scaledMinus(scale, dark.getMaskedImage())

def updateVariance(exposure):
    mi = exposure.getMaskedImage()
    var = afwImage.ImageF(mi.getImage(), True)
    amp = cameraGeom.cast_Amp(exposure.getDetector())
    ep = amp.getElectronicParams()
    gain = ep.getGain()
    var /= gain
    mi = afwImage.makeMaskedImage(mi.getImage(), mi.getMask(), var)
    exposure.setMaskedImage(mi)
    #stdev = afwMath.makeStatistics(mi, afwMath.STDEVCLIP).getValue()
    #varmean = afwMath.makeStatistics(var, afwMath.MEANCLIP).getValue()
    
    
                   
def flatCorrection(exposure, flat, scalingtype, scaling = 1.0):

    flatscaling = 1.0
    # Figure out scaling from the data
    if scalingtype == 'MEAN':
        flatscaling = afwMath.makeStatistics(flat.getMaskedImage().getImage(), afwMath.MEAN).getValue(afwMath.MEAN)
    elif scalingtype == 'MEDIAN':
        flatscaling = afwMath.makeStatistics(flat.getMaskedImage().getImage(), afwMath.MEDIAN).getValue(afwMath.MEDIAN)
    elif scalingtype == 'USER':
        flatscaling = scaling
    else:
        raise pexExcept.LsstException, '%s : %s not implemented' % ("flatCorrection", scalingtype)            
    mi   = exposure.getMaskedImage()
    fmi  = flat.getMaskedImage()
    mi.scaledDivides(1./flatscaling, fmi)
    


def illuminationCorrection(exposure, illum, illumscaling):

    # common input test

    mi   = exposure.getMaskedImage()
    mi.scaledDivides(1./illumscaling, illum.getMaskedImage())
    


def trimNew(exposure, ampBBox, trimsec=None, trimsecKeyword='trimsec'):
    """
    This returns a new Exposure that is a subsection of the input exposure.
    
    NOTE : do we need to deal with the WCS in any way, shape, or form?
    """
    methodName = 'trimNew'
    
    # common input test

    metadata   = exposure.getMetadata()

    if metadata.exists(trimsecKeyword) and trimsec==None:
        trimsec = metadata.getString(trimsecKeyword)


    if trimsec == None:
        raise pexExcept.LsstException, '%s : cannot find trimsec' % (methodName)        

    trimsecBBox = isrLib.BBoxFromDatasec(trimsec)

    if not (trimsecBBox.getDimensions() == ampBBox.getDimensions()):
        raise pexException.LsstException, '%s : amp bounding box not same as\
        trim section'%(methodName)

    trimmedExposure = afwImage.ExposureF(exposure, trimsecBBox)
    trimmedExposure.getMaskedImage().setXY0(ampBBox.getLLC())

    # remove trimsec from metadata
    trimmedExposure.getMetadata().remove(trimsecKeyword)
    # n.b. what other changes are needed here?
    # e.g. wcs info, overscan, etc
    

    return trimmedExposure



def overscanCorrection(exposure, overscanBBox, fittype, overscanKeyword =
        'overscan', polyorder = 1):
    """
    """

    # common input test
    mi = exposure.getMaskedImage()

    # if "True", do a deep copy
    overscanData    = afwImage.ImageF(exposure.getMaskedImage().getImage(), overscanBBox, False)

    # what type of overscan modeling?
    offset = 0
    if fittype == 'MEAN':
        offset = afwMath.makeStatistics(overscanData, afwMath.MEAN).getValue(afwMath.MEAN)
        mi    -= offset
    elif fittype == 'MEDIAN':
        offset = afwMath.makeStatistics(overscanData, afwMath.MEDIAN).getValue(afwMath.MEDIAN)
        mi    -= offset
    elif fittype == 'POLY':
        raise pexExcept.LsstException, '%s : %s not implemented' % ("overscanCorrection", fittype)
    else:
        raise pexExcept.LsstException, '%s : %s an invalid overscan type' % ("overscanCorrection", fittype)

    return 


def fringeCorrection(exposure, fringe):

    raise pexExcept.LsstException, '%s not implemented' % ("ipIsr.fringCorrection")


def pupilCorrection(exposure, pupil):

    raise pexExcept.LsstException, '%s not implemented' % (stageName)
