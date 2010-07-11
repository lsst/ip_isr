#!/usr/bin/env python

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
