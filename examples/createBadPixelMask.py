import lsst.afw.detection as afwDetection
import lsst.afw.image     as afwImage

def maskFromImage(fitsfile):
    # input bad pixel image
    maskImage   = afwImage.ImageF(fitsfile)
    maskImage  *= -1

    # turn into masked image for detection
    maskedImage = afwImage.MaskedImageF(maskImage)

    # find bad regions
    thresh    = afwDetection.Threshold(-0.5)
    ds        = afwDetection.DetectionSetF(maskedImage, thresh)
    fpList    = ds.getFootprints()

    # the output LSST Mask image
    maskMask  = afwImage.MaskU(maskImage.getDimensions())
    bitmask   = maskMask.getPlaneBitMask('BAD')

    # set the bits
    for fp in fpList:
        #for s in fp.getSpans():
        #    print s.toString()
        afwDetection.setMaskFromFootprint(maskMask, fp, bitmask)

    return maskMask

if __name__ == '__main__':
    import sys
    mask = maskFromImage(sys.argv[1])
    mask.writeFits('Mask.fits')
