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
    """Make a double Gaussian PSF
    
    @param[in]      fwhm            FWHM of double Gaussian smoothing kernel
    @return psf
    """
    ksize = 4*int(fwhm) + 1
    return measAlg.DoubleGaussianPsf(ksize, ksize, fwhm/(2*math.sqrt(2*math.log(2))))

def calcEffectiveGain(maskedImage):
    """Calculate effective gain

    @param[in]      maskedImage     masked image to process
    @return (median gain, mean gain) in e-/ADU
    """
    im = afwImage.ImageF(maskedImage.getImage(), True)
    var = maskedImage.getVariance()
    im /= var
    medgain = afwMath.makeStatistics(im, afwMath.MEDIAN).getValue()
    meangain = afwMath.makeStatistics(im, afwMath.MEANCLIP).getValue()
    return medgain, meangain

def transposeMaskedImage(maskedImage):
    """Make a transposed copy of a masked image

    @param[in]      maskedImage     masked image to process
    @return transposed masked image
    """
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
    fs = afwDetection.FootprintSet(satmaskim, thresh)
    for f in fs.getFootprints():
        metrics['nSaturatePix'] += f.getNpix()
    fs = afwDetection.FootprintSet(badmaskim, thresh)
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
    fs = afwDetection.FootprintSet(satmaskim, thresh)
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
    """Interpolate over defects specified in a defect list

    @param[in,out]  maskedImage     masked image to process
    @param[in]      defectList      defect list
    @param[in]      fwhm            FWHM of double Gaussian smoothing kernel
    @param[in]      fallbackValue   fallback value if an interpolated value cannot be determined;
                                    if None then use clipped mean image value
    """
    psf = createPsf(fwhm)
    if fallbackValue is None:
        fallbackValue = afwMath.makeStatistics(maskedImage.getImage(), afwMath.MEANCLIP).getValue()
    if 'INTRP' not in maskedImage.getMask().getMaskPlaneDict().keys():
        maskedImage.getMask.addMaskPlane('INTRP')
    measAlg.interpolateOverDefects(maskedImage, psf, defectList, fallbackValue)

def defectListFromFootprintList(fpList, growFootprints=1):
    """Compute a defect list from a footprint list, optionally growing the footprints
    
    @param[in]      fpList          footprint list
    @param[in]      growFootprints  amount by which to grow footprints of detected regions
    @return defect list
    """
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
    """Make a transposed copy of a defect list
    
    @param[in]      defectList      defect list
    @return defect list with transposed defects
    """
    retDefectList = measAlg.DefectListT()
    for defect in defectList:
        bbox = defect.getBBox()
        nbbox = afwGeom.Box2I(afwGeom.Point2I(bbox.getMinY(), bbox.getMinX()), 
             afwGeom.Extent2I(bbox.getDimensions()[1], bbox.getDimensions()[0]))
        retDefectList.push_back(measAlg.Defect(nbbox))
    return retDefectList

def maskPixelsFromDefectList(maskedImage, defectList, maskName='BAD'):
    """Set mask plane based on a defect list

    @param[in,out]  maskedImage     masked image to process; mask plane is updated
    @param[in]      defectList      defect list
    @param[in]      maskName        mask plane name
    """
    # mask bad pixels
    mask = maskedImage.getMask()
    bitmask = mask.getPlaneBitMask(maskName)
    for defect in defectList:
        bbox = defect.getBBox()
        afwDetection.setMaskFromFootprint(mask, afwDetection.Footprint(bbox), bitmask)

def getDefectListFromMask(maskedImage, maskName, growFootprints=1):
    """Compute a defect list from a specified mask plane

    @param[in]      maskedImage     masked image to process
    @param[in]      maskName        mask plane name
    @param[in]      growFootprints  amount by which to grow footprints of detected regions
    """
    mask = maskedImage.getMask()
    workmask = afwImage.MaskU(mask, True)
    workmask &= mask.getPlaneBitMask(maskName)
    thresh = afwDetection.Threshold(0.5)
    maskimg = afwImage.ImageU(workmask.getBBox(afwImage.PARENT))
    maskimg <<= workmask
    ds = afwDetection.FootprintSet(maskimg, thresh)
    fpList = ds.getFootprints()
    return defectListFromFootprintList(fpList, growFootprints)

def makeThresholdMask(maskedImage, threshold, growFootprints=1, maskName = 'SAT'):
    """Mask pixels based on threshold detection
    
    @param[in,out]  maskedImage     masked image to process; the mask is altered
    @param[in]      threshold       detection threshold
    @param[in]      growFootprints  amount by which to grow footprints of detected regions
    @param[in]      maskName        mask plane name
    """
    # find saturated regions
    thresh = afwDetection.Threshold(threshold)
    fs = afwDetection.FootprintSet(maskedImage, thresh)

    if growFootprints > 0:
        fs = afwDetection.FootprintSet(fs, growFootprints)

    fpList = fs.getFootprints()
    # set mask
    mask = maskedImage.getMask()
    bitmask = mask.getPlaneBitMask(maskName)
    afwDetection.setMaskFromFootprintList(mask, fpList, bitmask)

    return defectListFromFootprintList(fpList, growFootprints=0)

def interpolateFromMask(maskedImage, fwhm, growFootprints=1, maskName='SAT', fallbackValue=None):
    """Interpolate over defects identified by a particular mask plane
    
    @param[in,out]  maskedImage     masked image to process
    @param[in]      fwhm            FWHM of double Gaussian smoothing kernel
    @param[in]      growFootprints  amount by which to grow footprints of detected regions
    @param[in]      maskName        mask plane name
    @param[in]      fallbackValue   value of last resort for interpolation
    """
    defectList = getDefectListFromMask(maskedImage, maskName, growFootprints)
    interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=fallbackValue)

def saturationCorrection(maskedImage, saturation, fwhm, growFootprints=1, interpolate=True, maskName='SAT',
                         fallbackValue=None):
    """Mark saturated pixels and optionally interpolate over them

    @param[in,out]  maskedImage     masked image to process
    @param[in]      saturation      saturation level (used as a detection threshold)
    @param[in]      fwhm            FWHM of double Gaussian smoothing kernel
    @param[in]      growFootprints  amount by which to grow footprints of detected regions
    @param[in]      interpolate     interpolate over saturated pixels?
    @param[in]      maskName        mask plane name
    @param[in]      fallbackValue   value of last resort for interpolation
    """
    defectList = makeThresholdMask(
        maskedImage = maskedImage,
        threshold = saturation,
        growFootprints = growFootprints,
        maskName = maskName,
    )
    if interpolate:
        interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=fallbackValue)

def biasCorrection(maskedImage, biasMaskedImage):
    """Apply bias correction in place

    @param[in,out]  maskedImage     masked image to correct
    @param[in]      biasMaskedImage bias, as a masked image
    """
    maskedImage -= biasMaskedImage

def darkCorrection(maskedImage, darkMaskedImage, expScale, darkScale):
    """Apply dark correction in place
    
    maskedImage -= dark * expScaling / darkScaling

    @param[in,out]  maskedImage     masked image to correct
    @param[in]      darkMaskedImage dark masked image
    @param[in]      expScale        exposure scale
    @param[in]      darkScale       dark scale
    """
    if maskedImage.getBBox() != darkMaskedImage.getBBox():
        raise RuntimeError("maskedImage bbox %s != darkMaskedImage bbox %s" % \
            (maskedImage.getBBox(), darkMaskedImage.getBBox()))

    scale = expScale / darkScale
    maskedImage.scaledMinus(scale, darkMaskedImage)

def updateVariance(maskedImage, gain, readNoise):
    """Set the variance plane based on the image plane

    @param[in,out]  maskedImage     masked image; image plane is read and variance plane is written
    @param[in]      gain            amplifier gain (e-/ADU)
    @param[in]      readNoise       amplifier read noise (ADU/pixel)
    """
    var = maskedImage.getVariance()
    var <<= maskedImage.getImage()
    var /= gain
    var += readNoise**2

def flatCorrection(maskedImage, flatMaskedImage, scalingType, userScale=1.0):
    """Apply flat correction in place

    @param[in,out]  maskedImage     masked image to correct
    @param[in]      flatMaskedImage flat field masked image
    @param[in]      scalingType     how to compute flat scale; one of 'MEAN', 'MEDIAN' or 'USER'
    @param[in]      userScale       scale to use if scalingType is 'USER', else ignored
    """
    if maskedImage.getBBox() != flatMaskedImage.getBBox():
        raise RuntimeError("maskedImage bbox %s != flatMaskedImage bbox %s" % \
            (maskedImage.getBBox(), flatMaskedImage.getBBox()))

    # Figure out scale from the data
    # I'm not sure we should be doing this here, but maybe
    if scalingType == 'MEAN':
        flatScale = afwMath.makeStatistics(flatMaskedImage.getImage(), afwMath.MEAN).getValue(afwMath.MEAN)
    elif scalingType == 'MEDIAN':
        flatScale = afwMath.makeStatistics(flatMaskedImage.getImage(), afwMath.MEDIAN).getValue(afwMath.MEDIAN)
    elif scalingType == 'USER':
        flatScale = userScale
    else:
        raise pexExcept.LsstException, '%s : %s not implemented' % ("flatCorrection", scalingType)
    
    maskedImage.scaledDivides(1.0/flatScale, flatMaskedImage)

def illuminationCorrection(maskedImage, illumMaskedImage, illumScale):
    """Apply illumination correction in place

    @param[in,out]  maskedImage     masked image to correct
    @param[in]      illumMaskedImage illumination correction masked image
    @param[in]      illumScale      scale value for illumination correction
    """
    if maskedImage.getBBox() != illumMaskedImage.getBBox():
        raise RuntimeError("maskedImage bbox %s != illumMaskedImage bbox %s" % \
            (maskedImage.getBBox(), illumMaskedImage.getBBox()))

    maskedImage.scaledDivides(1./illumScale, illumMaskedImage)

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

def overscanCorrection(ampMaskedImage, overscanImage, fitType='MEDIAN', polyOrder=1):
    """Apply overscan correction in place

    @param[in,out]  ampMaskedImage  masked image to correct
    @param[in]      overscanImage   overscan data as an image
    @param[in]      fitType         type of fit for overscan correction; one of:
                                    - 'MEAN'
                                    - 'MEDIAN'
                                    - 'POLY'
    @param[in]      polyOrder       polynomial order (ignored unless fitType='POLY')
    """
    ampImage = ampMaskedImage.getImage()
    if fitType == 'MEAN':
        offImage = afwMath.makeStatistics(overscanImage, afwMath.MEAN).getValue(afwMath.MEAN)
    elif fitType == 'MEDIAN':
        offImage = afwMath.makeStatistics(overscanImage, afwMath.MEDIAN).getValue(afwMath.MEDIAN)
    elif fitType == 'POLY':
        biasArray = overscanImage.getArray()
        # Fit along the long axis, so take median of each short row and fit the resulting array
        shortInd = numpy.argmin(biasArray.shape)
        medianBiasArr = numpy.median(biasArray, axis=shortInd)
        coeffs = numpy.polyfit(range(len(medianBiasArr)), medianBiasArr, deg=polyOrder)
        fitBiasArr = numpy.polyval(coeffs, range(len(medianBiasArr)))
        offImage = ampImage.Factory(ampImage.getDimensions())
        offArray = offImage.getArray()
        if shortInd == 1:
            print "shortInd=1"
            offArray[:,:] = fitBiasArr[:,numpy.newaxis]
        else:
            offArray[:,:] = fitBiasArr[numpy.newaxis,:]
    else:
        raise pexExcept.LsstException, '%s : %s an invalid overscan type' % \
            ("overscanCorrection", fitType)
    ampImage -= offImage

# def fringeCorrection(maskedImage, fringe):
#     raise NotImplementedError()

# def pupilCorrection(maskedImage, pupil):
#     raise NotImplementedError()

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

