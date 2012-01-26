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

import math, re
import numpy
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
class Isr(object):
    def __init__(self, display=False):
        self.display = display
    def createPsf(self, fwhm):
        """Make a PSF"""
        ksize = 4*int(fwhm) + 1
        return afwDetection.createPsf('DoubleGaussian', ksize, ksize, fwhm/(2*math.sqrt(2*math.log(2))))

    def calcEffectiveGain(self, exposure):
        mi = exposure.getMaskedImage()
        im = afwImage.ImageF(mi.getImage(), True)
        var = mi.getVariance()
        im /= var
        medgain = afwMath.makeStatistics(im, afwMath.MEDIAN).getValue()
        meangain = afwMath.makeStatistics(im, afwMath.MEANCLIP).getValue()
        return medgain, meangain

    def convertImageForIsr(self, exposure):
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

    def calculateSdqaCcdRatings(self, exposure):
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
        stats = afwMath.makeStatistics(mi, afwMath.MEANCLIP | \
            afwMath.STDEVCLIP | afwMath.MEDIAN | afwMath.MIN |\
            afwMath.MAX, sctrl)
        metrics['imageClipMean4Sig3Pass'] = stats.getValue(afwMath.MEANCLIP)
        metrics['imageSigma'] = stats.getValue(afwMath.STDEVCLIP)
        metrics['imageMedian'] = stats.getValue(afwMath.MEDIAN)
        metrics['imageMin'] = stats.getValue(afwMath.MIN)
        metrics['imageMax'] = stats.getValue(afwMath.MAX)
        for k in metrics.keys():
            metadata.set(k, metrics[k])

    def calculateSdqaAmpRatings(self, exposure, biasBBox, dataBBox):
        metrics = {}
        metrics['nSaturatePix'] = 0
        metrics['overscanMean'] = None
        metrics['overscanStdDev'] = None
        metrics['overscanMedian'] = None
        metrics['overscanMin'] = None
        metrics['overscanMax'] = None
        mi = exposure.getMaskedImage()
        metadata = exposure.getMetadata()
        trimmi = afwImage.MaskedImageF(mi, dataBBox, False)
        biasmi = afwImage.MaskedImageF(mi, biasBBox, False)
        mask = mi.getMask()
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

    def maskFromDefects(self, dimensions, fpList):
        # the output LSST Mask image
        minpt = afwGeom.Point2I(0,0)
        maxpt = afwGeom.Point2I(dimensions[0], dimensions[1])
        mask = afwImage.MaskU(afwGeom.Box2I(minpt, maxpt), afwImage.LOCAL)
        mask.set(0)
        bitmask = mask.getPlaneBitMask('BAD')

        # set the bits
        afwDetection.setMaskFromFootprintList(mask, fpList, bitmask)

        return mask


    def defectsFromBoolImage(self, fitsfile, invert=False):
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

    def interpolateDefectList(self, exposure, defectList, fwhm, fallbackValue=None):
        mi = exposure.getMaskedImage()
        psf = self.createPsf(fwhm)
        if fallbackValue is None:
            fallbackValue = afwMath.makeStatistics(mi.getImage(), afwMath.MEANCLIP).getValue()
        measAlg.interpolateOverDefects(mi, psf, defectList, fallbackValue)

    def maskBadPixelsDef(self, exposure, defectList, fwhm=None,
                 interpolate=True,
                 maskName='BAD'):

        # mask bad pixels
        mi = exposure.getMaskedImage()
        mask = mi.getMask()
        bitmask = mask.getPlaneBitMask(maskName)
        for defect in defectList:
            bbox = defect.getBBox()
            afwDetection.setMaskFromFootprint(mask, afwDetection.Footprint(bbox), bitmask)

        if interpolate:
            # and interpolate over them
            assert fwhm and fwhm > 0, "FWHM not provided for interpolation"
            interpolateDefectList(exposure, defectList, fwhm)

    def linearization(self, exposure, lookupTable):

        # common input test
        metadata   = exposure.getMetadata()
        gain = metadata.get('gain')
        mi   = exposure.getMaskedImage()
        lookupTable.apply(mi, gain)

    def saturationDetection(self, exposure, saturation, doMask = True, maskName = 'SAT'):

        mi = exposure.getMaskedImage()
        if self.display:
            ds9.mtv(mi, frame=0)


        # find saturated regions
        thresh = afwDetection.Threshold(saturation)
        ds = afwDetection.makeFootprintSet(mi, thresh)
        fpList = ds.getFootprints()
        # we will turn them into defects for interpolating
        defectList = measAlg.DefectListT()

        # grow them
        bboxes = []
        for fp in fpList:
            if doMask:
                mask = mi.getMask()
                bitmask = mask.getPlaneBitMask(maskName)
                afwDetection.setMaskFromFootprint(mask, fp, bitmask)
            bboxes.append(afwDetection.footprintToBBoxList(fp))
        return bboxes

    def saturationInterpolation(self, exposure, fwhm, growFootprints = 1, maskName = 'SAT'):
        mi = exposure.getMaskedImage()
        mask = mi.getMask()
        satmask = afwImage.MaskU(mask, True)
        satmask &= mask.getPlaneBitMask(maskName)
        thresh = afwDetection.Threshold(0.5)
        maskimg = afwImage.ImageU(satmask.getBBox(afwImage.PARENT))
        maskimg <<= satmask
        ds = afwDetection.makeFootprintSet(maskimg, thresh)
        fpList = ds.getFootprints()
        satDefectList = measAlg.DefectListT()
        for fp in fpList:
            if growFootprints > 0:
                # if "True", growing requires a convolution
                # if "False", its faster
                fpGrow = afwDetection.growFootprint(fp, growFootprints, False)
            else:
                fpGrow = fp
            for bbox in afwDetection.footprintToBBoxList(fpGrow):
                defect = measAlg.Defect(bbox)
                satDefectList.push_back(defect)
        if 'INTRP' not in mask.getMaskPlaneDict().keys():
            mask.addMaskPlane('INTRP')
        psf = self.createPsf(fwhm)

        measAlg.interpolateOverDefects(mi, psf, satDefectList)


    def saturationCorrection(self, exposure, saturation, fwhm, growFootprints=1,
                 interpolate = True,
                 maskName    = 'SAT'):
        #This will be slower than necessary.  Should write another method that does detection and interpolation at the same time.

        saturationDetection(exposure, saturation, doMask=True, maskName=maskName)
        if interpolate:
            saturationInterpolation(exposure, fwhm, growFootprints=growFootprints, maskName=maskName)
        if self.display:
            ds9.mtv(exposure.getmaskedImage(), frame=0)


    def biasCorrection(self, exposure, bias):

        mi = exposure.getMaskedImage()
        bmi = bias.getMaskedImage()
        mi -= bmi

    def darkCorrection(self, exposure, dark, expscaling, darkscaling):

        scale = expscaling / darkscaling
        mi  = exposure.getMaskedImage()
        mi.scaledMinus(scale, dark.getMaskedImage())

    def updateVariance(self, exposure):
        mi = exposure.getMaskedImage()
        var = afwImage.ImageF(mi.getImage(), True)
        amp = cameraGeom.cast_Amp(exposure.getDetector())
        ep = amp.getElectronicParams()
        gain = ep.getGain()
        var /= gain
        mi = afwImage.makeMaskedImage(mi.getImage(), mi.getMask(), var)
        exposure.setMaskedImage(mi)



    def flatCorrection(self, exposure, flat, scalingtype, scaling = 1.0):

        flatscaling = 1.0
        # Figure out scaling from the data
        # I'm not sure we should be doing this here, but maybe
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

        if self.display:
            ds9.mtv(mi, title="Flattened")

    def illuminationCorrection(self, exposure, illum, illumscaling):

        # common input test

        mi   = exposure.getMaskedImage()
        mi.scaledDivides(1./illumscaling, illum.getMaskedImage())



    def trimAmp(self, exposure, trimBbox=None):
        """
        This returns a new Exposure that is a subsection of the input exposure.

        NOTE : do we need to deal with the WCS in any way, shape, or form?
        """
        if trimBbox is not None:
            return exposureFactory(exposure, trimBbox, LOCAL)
        else:
            amp = cameraGeom.cast_Amp(exposure.getDetector())
            return exposureFactory(exposure, amp.getDiskDataSec(false), LOCAL)
        # n.b. what other changes are needed here?
        # e.g. wcs info, overscan, etc

    def overscanCorrection(self, exposure, overscanData, fittype='MEDIAN', polyorder=1, imageFactory=afwImage.ImageF):
        """
        """
        typemap = {afwImage.ImageU:numpy.uint16, afwImage.ImageI:numpy.int32, afwImage.ImageF:numpy.float32, afwImage.ImageD:numpy.float64}

        # common input test
        mi = exposure.getMaskedImage()

        # what type of overscan modeling?
        offset = 0
        if fittype == 'MEAN':
            offset = afwMath.makeStatistics(overscanData, afwMath.MEAN).getValue(afwMath.MEAN)
            mi    -= offset
        elif fittype == 'MEDIAN':
            offset = afwMath.makeStatistics(overscanData, afwMath.MEDIAN).getValue(afwMath.MEDIAN)
            mi    -= offset
        elif fittype == 'POLY':
            biasArray = overscanData.getArray()
            #Assume we want to fit along the long axis
            aind = numpy.argmin(biasArray.shape)
            find = numpy.argmin(biasArray.shape)
            fitarr = numpy.median(biasArray, axis=aind)
            coeffs = numpy.polyfit(range(len(fitarr)), fitarr, deg=polyorder)
            offsets = numpy.polyval(coeffs, range(len(fitarr)))
            width, height = mi.getDimensions()
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
            mi -= im 
        else:
            raise pexExcept.LsstException, '%s : %s an invalid overscan type' % ("overscanCorrection", fittype)

    def fringeCorrection(self, exposure, fringe):

        raise pexExcept.LsstException, '%s not implemented' % ("ipIsr.fringCorrection")


    def pupilCorrection(self, exposure, pupil):

        raise pexExcept.LsstException, '%s not implemented' % (stageName)

'''
class Linearization(object):
    def __init__(self, linearityFile):
    
    def apply(self, exposure):
        mi = exposure.getMaskedImage()
        	
    def readFile(self):
    def writeFile(self):
'''    
