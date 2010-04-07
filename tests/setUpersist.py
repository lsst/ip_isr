#!/usr/bin/env python
import eups
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import os

def PersistImageU(imagename):
    mymeta = dafBase.PropertySet()
    imf = afwImage.ImageF(imagename, 0, mymeta)
    mask = afwImage.MaskU(imf.getDimensions())
    mask.set(0)
    var = afwImage.ImageF(imf)
    mi = afwImage.makeMaskedImage(imf, mask, var)
    exp = afwImage.ExposureF(mi, afwImage.Wcs())
    exp.setMetadata(mymeta)
    exp.writeFits("test_out.fits")

    md = afwImage.readMetadata("test_out.fits")
    if md.exists("BZERO"):
        print "Found BZERO in file"

if __name__ == "__main__":
    imsimdir = eups.productDir('ip_isr')
    PersistImageU(os.path.join(imsimdir,'tests','imsim_85751839_R23_S11_C00_E000.fits.gz'))
