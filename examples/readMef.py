#!/usr/bin/env python
# Read and get header information from a MEF file.

import os
import re
import sys
import pyfits

import lsst.afw.image as afwImage

if __name__ == "__main__":
    
    if len(sys.argv) != 2:
        print 'Usage : readMef.py infile.fits'
        sys.exit(1)
       
    inImage = sys.argv[1]
    if not os.path.isfile(inImage):
        print 'ERROR: Cannot locate', inImage
        sys.exit(1)
        
##     try:
##         print 'Reading fits file with afwImage.ImageD.readfits'        
##         inFile = afwImage.ImageD().readfits(inImage)
##     except:
##         print 'ERROR: Cannot Read', inImage
##         sys.exit(1)
        
    try:
        print 'Reading fits file with pyfits.open'
        inFile  = pyfits.open(inImage)
    except:
        print 'ERROR: Cannot open', inImage
        sys.exit(1)
        
    nExt = len(inFile)
    print 'Number of Extensions in %s :' % inImage, nExt
    print 'First extension [00] is a header'

    for extension in range(1,nExt):	
        header = inFile[extension].header
        gain = header['GAIN']
        print 'Gain for extension %d :' % extension, gain
        rdNoise = header['RDNOISE']
        print 'RDNoise for extension %d :' % extension, rdNoise

#    outFile = re.sub('.fits', '_%d_image.fits' % extension, inFile)
#    print 'Writing ', outFile
#    afwImage.ImageD().writefits(outile)


