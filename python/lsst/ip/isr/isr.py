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
import math

import numpy

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg

def createPsf(fwhm):
    """Make a PSF"""
    ksize = 4*int(fwhm) + 1
    return afwDetection.createPsf('DoubleGaussian', ksize, ksize, fwhm/(2*math.sqrt(2*math.log(2))))

def calcEffectiveGain(maskedImage):
    im = afwImage.ImageF(maskedImage.getImage(), True)
    var = maskedImage.getVariance()
    im /= var
    medgain = afwMath.makeStatistics(im, afwMath.MEDIAN).getValue()
    meangain = afwMath.makeStatistics(im, afwMath.MEANCLIP).getValue()
    return medgain, meangain

def transposeMaskedImage(maskedImage):
    imarr = maskedImage.getImage().getArray().T.__copy__()
    vararr = maskedImage.getVariance().getArray().T.__copy__()
    maskarr = maskedImage.getMask().getArray().T.__copy__()
    return afwImage.makeMaskedImageFromArrays(imarr, maskarr, vararr)


def calculateSdqaCcdRatings(maskedImage, metadata):
    metrics = {}
    metrics['nSaturatePix'] = 0
    metrics['nBadCalibPix'] = 0
    metrics['imageClipMean4Sig3Pass'] = None
    metrics['imageSigma'] = None
    metrics['imageMedian'] = None
    metrics['imageMin'] = None
    metrics['imageMax'] = None
    mask = maskedImage.getMask()
    badbitmask = mask.getPlaneBitMask('BAD')
    satbitmask = mask.getPlaneBitMask('SAT')
    intrpbitmask = mask.getPlaneBitMask('INTRP')
    sctrl = afwMath.StatisticsControl()
    sctrl.setNumIter(3)
    sctrl.setNumSigmaClip(4)
    sctrl.setAndMask(satbitmask | badbitmask | intrpbitmask)
    satmask = afwImage.MaskU(mask, True)
    badmask = afwImage.MaskU(mask, True)
    satmask &= satbitmask
    badmask &= badbitmask
    satmaskim = afwImage.ImageU(satmask.getBBox(afwImage.PARENT))
    satmaskim <<= satmask
    badmaskim = afwImage.ImageU(badmask.getBBox(afwImage.PARENT))
    badmaskim <<= badmask
    thresh = afwDetection.Threshold(0.5)
    fs = afwDetection.makeFootprintSet(satmaskim, thresh)
    for f in fs.getFootprints():
        metrics['nSaturatePix'] += f.getNpix()
    fs = afwDetection.makeFootprintSet(badmaskim, thresh)
    for f in fs.getFootprints():
        metrics['nBadCalibPix'] += f.getNpix()
    stats = afwMath.makeStatistics(maskedImage, afwMath.MEANCLIP | \
        afwMath.STDEVCLIP | afwMath.MEDIAN | afwMath.MIN |\
        afwMath.MAX, sctrl)
    metrics['imageClipMean4Sig3Pass'] = stats.getValue(afwMath.MEANCLIP)
    metrics['imageSigma'] = stats.getValue(afwMath.STDEVCLIP)
    metrics['imageMedian'] = stats.getValue(afwMath.MEDIAN)
    metrics['imageMin'] = stats.getValue(afwMath.MIN)
    metrics['imageMax'] = stats.getValue(afwMath.MAX)
    for k in metrics.keys():
        metadata.set(k, metrics[k])

def calculateSdqaAmpRatings(maskedImage, metadata, biasBBox, dataBBox):
    metrics = {}
    metrics['nSaturatePix'] = 0
    metrics['overscanMean'] = None
    metrics['overscanStdDev'] = None
    metrics['overscanMedian'] = None
    metrics['overscanMin'] = None
    metrics['overscanMax'] = None
    trimmi = afwImage.MaskedImageF(maskedImage, dataBBox, False)
    biasmi = afwImage.MaskedImageF(maskedImage, biasBBox, False)
    mask = maskedImage.getMask()
    satbitmask = mask.getPlaneBitMask('SAT')
    sctrl = afwMath.StatisticsControl()
    sctrl.setAndMask(satbitmask)
    satmask = trimmi.getMask()
    satmask &= satbitmask
    satmaskim = afwImage.ImageU(satmask.getBBox(afwImage.PARENT))
    satmaskim <<= satmask
    thresh = afwDetection.Threshold(0.5)
    fs = afwDetection.makeFootprintSet(satmaskim, thresh)
    for f in fs.getFootprints():
        metrics['nSaturatePix'] += f.getNpix()
    stats = afwMath.makeStatistics(biasmi, afwMath.MEAN | \
        afwMath.STDEV | afwMath.MEDIAN | afwMath.MIN |\
        afwMath.MAX,sctrl)
    metrics['overscanMean'] = stats.getValue(afwMath.MEAN)
    metrics['overscanStdDev'] = stats.getValue(afwMath.STDEV)
    metrics['overscanMedian'] = stats.getValue(afwMath.MEDIAN)
    metrics['overscanMin'] = stats.getValue(afwMath.MIN)
    metrics['overscanMax'] = stats.getValue(afwMath.MAX)
    for k in metrics.keys():
        metadata.set(k, metrics[k])

def interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=None):
    psf = createPsf(fwhm)
    if fallbackValue is None:
        fallbackValue = afwMath.makeStatistics(maskedImage.getImage(), afwMath.MEANCLIP).getValue()
    measAlg.interpolateOverDefects(maskedImage, psf, defectList, fallbackValue)

def defectListFromFootprintList(fpList, growFootprints=1):
    defectList = measAlg.DefectListT()
    for fp in fpList:
        if growFootprints > 0:
            # if "True", growing requires a convolution
            # if "False", its faster
            fpGrow = afwDetection.growFootprint(fp, growFootprints, False)
        else:
            fpGrow = fp
        for bbox in afwDetection.footprintToBBoxList(fpGrow):
            defect = measAlg.Defect(bbox)
            defectList.push_back(defect)
    return defectList

def transposeDefectList(defectList):
    retDefectList = measAlg.DefectListT()
    for defect in defectList:
        bbox = defect.getBBox()
        nbbox = afwGeom.Box2I(afwGeom.Point2I(bbox.getMinY(), bbox.getMinX()), 
             afwGeom.Extent2I(bbox.getDimensions()[1], bbox.getDimensions()[0]))
        retDefectList.push_back(measAlg.Defect(nbbox))
    return retDefectList

def maskPixelsFromDefectList(maskedImage, defectList, maskName='BAD'):
    # mask bad pixels
    mask = maskedImage.getMask()
    bitmask = mask.getPlaneBitMask(maskName)
    for defect in defectList:
        bbox = defect.getBBox()
        afwDetection.setMaskFromFootprint(mask, afwDetection.Footprint(bbox), bitmask)

def getDefectListFromMask(maskedImage, maskName, growFootprints=1):
    mask = maskedImage.getMask()
    workmask = afwImage.MaskU(mask, True)
    workmask &= mask.getPlaneBitMask(maskName)
    thresh = afwDetection.Threshold(0.5)
    maskimg = afwImage.ImageU(workmask.getBBox(afwImage.PARENT))
    maskimg <<= workmask
    ds = afwDetection.makeFootprintSet(maskimg, thresh)
    fpList = ds.getFootprints()
    return defectListFromFootprintList(fpList, growFootprints)

def makeThresholdMask(maskedImage, threshold, growFootprints=1, maskName = 'SAT'):
    # find saturated regions
    thresh = afwDetection.Threshold(threshold)
    ds = afwDetection.makeFootprintSet(maskedImage, thresh)
    fpList = ds.getFootprints()
    # set mask
    mask = maskedImage.getMask()
    bitmask = mask.getPlaneBitMask(maskName)
    if growFootprints > 0:
        for fp in fpList:
            fp = afwDetection.growFootprint(fp, growFootprints)
    afwDetection.setMaskFromFootprintList(mask, fpList, bitmask)

    return defectListFromFootprintList(fpList, growFootprints=0)

def interpolateFromMask(maskedImage, fwhm, growFootprints = 1, maskName = 'SAT'):
    defectList = getDefectListFromMask(maskedImage, maskName, growFootprints)
    if 'INTRP' not in maskedImage.getMask().getMaskPlaneDict().keys():
        maskedImage.getMask.addMaskPlane('INTRP')
    psf = createPsf(fwhm)
    measAlg.interpolateOverDefects(maskedImage, psf, defectList)

def saturationCorrection(maskedImage, saturation, fwhm, growFootprints=1, interpolate = True, maskName = 'SAT'):
    defectList = makeThresholdMask(maskedImage, saturation, grwoFootprints=growFootprints, maskName=maskName)
    if interpolate:
        measAlg.interpolateOverDefects(maskedImage, createPsf(fwhm), defectList)

def biasCorrection(maskedImage, biasMaskedImage):
    maskedImage -= biasMaskedImage

def darkCorrection(maskedImage, darkMaskedImage, expscaling, darkscaling):
    scale = expscaling / darkscaling
    maskedImage.scaledMinus(scale, darkMaskedImage)

def updateVariance(maskedImage, gain):
    var = maskedImage.getVariance()
    var <<= maskedImage.getImage()
    var /= gain

def flatCorrection(maskedImage, flatMaskedImage, scalingtype, scaling = 1.0):
    flatscaling = 1.0
    # Figure out scaling from the data
    # I'm not sure we should be doing this here, but maybe
    if scalingtype == 'MEAN':
        flatscaling = afwMath.makeStatistics(flatMaskedImage.getImage(), afwMath.MEAN).getValue(afwMath.MEAN)
    elif scalingtype == 'MEDIAN':
        flatscaling = afwMath.makeStatistics(flatMaskedImage.getImage(), afwMath.MEDIAN).getValue(afwMath.MEDIAN)
    elif scalingtype == 'USER':
        flatscaling = scaling
    else:
        raise pexExcept.LsstException, '%s : %s not implemented' % ("flatCorrection", scalingtype)
    
    maskedImage.scaledDivides(1./flatscaling, flatMaskedImage)

def illuminationCorrection(maskedImage, illumMaskedImage, illumscaling):
    # common input test
    maskedImage.scaledDivides(1./illumscaling, illumMaskedImage)


def trimAmp(exposure, trimBbox=None):
    """Return a new Exposure that is a subsection of the input exposure.

    NOTE: do we need to deal with the WCS in any way, shape, or form?
    """
    if trimBbox is not None:
        return exposureFactory(exposure, trimBbox, LOCAL)
    else:
        amp = cameraGeom.cast_Amp(exposure.getDetector())
        return exposureFactory(exposure, amp.getDiskDataSec(false), LOCAL)
    # n.b. what other changes are needed here?
    # e.g. wcs info, overscan, etc

def overscanCorrection(maskedImage, overscanData, fittype='MEDIAN', polyorder=1, imageFactory=afwImage.ImageF):
    """
    """
    typemap = {afwImage.ImageU:numpy.uint16, afwImage.ImageI:numpy.int32, afwImage.ImageF:numpy.float32, afwImage.ImageD:numpy.float64}

    # what type of overscan modeling?
    offset = 0
    if fittype == 'MEAN':
        offset = afwMath.makeStatistics(overscanData, afwMath.MEAN).getValue(afwMath.MEAN)
        maskedImage -= offset
    elif fittype == 'MEDIAN':
        offset = afwMath.makeStatistics(overscanData, afwMath.MEDIAN).getValue(afwMath.MEDIAN)
        maskedImage -= offset
    elif fittype == 'POLY':
        biasArray = overscanData.getArray()
        #Assume we want to fit along the long axis
        aind = numpy.argmin(biasArray.shape)
        find = numpy.argmin(biasArray.shape)
        fitarr = numpy.median(biasArray, axis=aind)
        coeffs = numpy.polyfit(range(len(fitarr)), fitarr, deg=polyorder)
        offsets = numpy.polyval(coeffs, range(len(fitarr)))
        width, height = maskedImage.getDimensions()
        offarr = numpy.zeros((height, width), dtype = typemap[imageFactory])
        if aind == 1:
            for i in range(len(offsets)):
                offarr[i] = offsets[i]
        elif aind == 0:
            offarr = offarr.T
            for i in range(len(offsets)):
                offarr[i] = offsets[i]
            offarr = offarr.T
        else:
            raise pexExcept.LsstException, "Non-2D array returned from MaskedImage.getArray()"
        im = afwImage.makeImageFromArray(offarr)
        maskedImage -= im 
    else:
        raise pexExcept.LsstException, '%s : %s an invalid overscan type' % ("overscanCorrection", fittype)

def fringeCorrection(maskedImage, fringe):
    raise pexExcept.LsstException, '%s not implemented' % ("ipIsr.fringCorrection")


def pupilCorrection(maskedImage, pupil):

    raise pexExcept.LsstException, '%s not implemented' % (stageName)

# class Linearization(object):
#     def __init__(linearityFile=None):
#         if linearityFile is not None:
#             self.readFile(linearityFile)
#         else:
#             self.makeLinearReplace()
#     def getImageFactoryFromExposure(exposure):
#         return exposure.getMaskedImage().getImage().Factory
#     def apply(exposure):
#         self.getImageFactoryFromExposure(exposure)
#         if type is "LUT":
#             imageData = exposure.getMaskedImage(). 
#         mi = exposure.getMaskedImage()
#         	
#     def readFile(filename):
#     def writeFile(filename):
#     def makeLinearReplace(self):
#         self.type = 
#     def makeLinearMult(self):

