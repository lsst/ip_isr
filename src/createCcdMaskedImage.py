# file createCcdMaskedImage.py
#
# author Nicole Silvestri, University of Washington
#
# contact nms@astro.washington.edu
#
#

# This code is designed to be a pre-stage to the nightly IPP.  It will
# read (currently with PyFits) a list of multi-extension FITS [MEF]
# files from any telescope/camera configuration and break it into its
# individual CCD MaskedImages (in the form of
# imType_mjd_filter_image#_ccd#_img|var|msk.fits). The CCD
# MaskedImages are then written to FitsStorage (the destinationDir).

# The code requires two input files:
#
# 1. inputFileList.txt: the list of MEF image and bad pixel mask files, one
#    pair per line, in the form of 'imageFileName.fits badPixelMaskName.fits' 
# 2. policyFile.paf: the file containing information specific to the data
#    being split (see cfhtDataPolicy.paf as an example)

import gzip
import numarray
import os
import pyfits
import re
import string
import sys

import eups

import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.detection.detectionLib as detectionLib
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy

dataDir = eups.productDir("afwdata")
if not dataDir:
    raise RuntimeError("Must set up afwdata to run this program.")
cfhtDataDir = os.path.join(dataDir, "/CFHT/D4/")
destinationDir = os.path.join(cfhtDataDir, "/testfiles")
ccdNumList = (0,1,2,3,4,5,6,7,8) # Specific CCDs for tests    
currDir = os.path.abspath(os.path.dirname(__file__))
     
def main():
   
    if len(sys.argv) != 3:
        print 'Usage : python createCcdMaskedImage.py inputFileList.txt policyFile.paf'
        sys.exit(1)

# Need to implement this if all input files are Gzipped        
#    inImageGzip = sys.argv[1]
#    inImage = re.sub('.gz', '', inImageGzip)

#    if not os.path.isfile(inImageGzip):
#    print 'Uncompressing MEF Image', inImageGzip
#    inF = gzip.GzipFile(inImageGzip, 'rb');
#    inTemp = inF.read()
#    inF.close()

#    print 'Writing Uncompressed MEF File', inImage
#    outF = file(inImage, 'wb');
#    outF.write(inTemp)
#    outF.close()

    inputFileList = sys.argv[1]
    if not os.path.isfile(inputFileList):
        print 'ERROR: Cannot locate', inputFileList
        sys.exit(1)

    policyFile = sys.argv[2]
    if not os.path.isfile(policyFile):
        print 'ERROR: Cannot locate', policyFile
        sys.exit(1)
    policy = pexPolicy.Policy.createPolicy(policyFile)   

    fileList = open(inputFileList, "rU")  
    files = fileList.readlines()
    numFiles = len(files)
    print 'Number of files to process: ', numFiles
        
    for image in files:
        # strip trailing whitespace, returns, etc.
        image = image.strip()
        # ignore blank lines
        if not image:
            continue
        # ignore comment lines
        if image.startswith("#"):
            continue
        imageList = image.split()
        if len(imageList) < 2 or len(imageList) > 3:
            print "Cannot parse: " image
        inImage, badPixFile = imageList[0:2]
        
        # afwImage readFits can not open MEF files. Using PyFits...
        try:
            print 'Reading Input Image Fits File with PyFits.'
            inFile  = open(inImage)
        except:
            print 'ERROR: Cannot open', inImage
            sys.exit(1)

        try:
            print 'Reading Bad Pixel Mask Fits File with PyFits.'
            badPixelMask = pyfits.open(badPixFile)
        except:
            print 'ERROR: Cannot open', badPixFile
            sys.exit(1)
     
        nExt = len(inFile)
        print 'Number of CCDs in %s :' % inImage, nExt-1

        nExtMask = len(badPixelMask)
        print 'Number of CCDs in %s :' % badPixFile, nExtMask-1

        # First extension [00] of an MEF is general header information
        generalHeader = inFile[0].header
    
        #for extension in range(1,nExt):
        for ccdNum in ccdNumList:     # split only CCDs given above for tests 
            extension = ccdNum + 1
            ccdHeader = inFile[extension].header

            # set up the output file names based on header info
            # NOTE: this is CFHTLS header-specific.
            # will need to generalize this...
            filterName = ccdHeader['FILTER']
            filt = filterName.split('.')
            filter = filt[0]
            mjdate = ccdHeader['MJDATE']
            mdate = int(mjdate)
            mjdDate = repr(mdate)
            imageRoot = re.sub('.fits', '', inImage)
            imageNumber = imageRoot[0:6]
            imageExt = re.sub(imageRoot, '_%d_img.fits' % (extension), imageRoot)

            # Write the individual Science Images to disk
            #outImageFileName = 'raw-' + mjdDate + '-' + filter + '-' + imageNumber + imageExt 
            #outImageFileName = re.sub('.fits', '_%d_img.fits' % (ccdNum), inImage)
            outImageFileName = re.sub('.fits', '_%d_img.fits' % (extension), inImage)
            outImgDir = os.path.join(destinationDir,outImageFileName);
            newInFile = pyfits.PrimaryHDU(inFile[extension].data, ccdHeader)
            newInFile.writeto(outImgDir)
            print 'Wrote:', outImageFileName
            scienceImage = afwImage.ImageF()
            scienceImage.readFits(outImgDir)
            print 'Reading: ', outImageFileName  

            # Synthesize the Varaince Image using the Science Image, gain,
            # and rdNoise and write it to disk
        
            gain = ccdHeader['GAIN']
            print 'Gain: ', gain
            rdNoise = ccdHeader['RDNOISE']
            print 'RdNoise: ', rdNoise
            varianceImage = afwImage.ImageF()
            varianceImage.readFits(outImgDir)

            # The varaince: sigma^2= DN/gain + (rdNoise^2/gain^2)
            scale = rdNoise**2 / gain**2
            varianceImage /= gain
            varianceImage += scale
            imageExt = re.sub(imageRoot, '_%d_var.fits' % (extension), imageRoot)
            #outVarianceFileName = 'raw-' + mjdDate + '-' + filter + '-' + imageNumber + imageExt
            #outVarianceFileName = re.sub('.fits', '_%d_var.fits' % (ccdNum), inImage)
            outVarianceFileName = re.sub('.fits', '_%d_var.fits' % (extension), inImage)
            outVarDir = os.path.join(destinationDir,outVarianceFileName)
            varianceImage.writeFits(outVarDir)
            print 'Wrote:', outVarianceFileName
        
           # Create and write the Mask. The following uses code from AB's
           # createCfhtMask.py which is specific to CFHT badPixelMask data.  This will
           # have to be modified to use badPixelMasks from other data sets later.
      
           maskFormat = re.compile('^\[(\d+):(\d+),(\d+):(\d+)\]$')
           maskPlane  = 'BAD' 
      
	   # Check that the bad pixel mask and the science image are the same size
           numRows = scienceImage.getRows()
           numCols = scienceImage.getCols()  
           numMaskRows, numMaskCols = badPixelMask[1].data.shape

           if numRows != numMaskRows or numCols != numMaskCols:
               print 'ERROR: Image and Bad Pixel Mask are different sizes.'
               sys.exit(1)

           header = badPixelMask[extension].header
	
           # make a new Mask image
           mask = afwImage.MaskUPtr(afwImage.MaskU(numMaskCols, numMaskRows))
   
           # note that this will have its mask planes initialized by default
           # we want to modify BAD
           badBitMask = mask.getPlaneBitMask(maskPlane)
       
           # put them all in a list and do at once
           footprintList = detectionLib.FootprintContainerT()
       
           for card in header.ascardlist().keys():
               if card.startswith('MASK_'):
                   maskArea = header[card]
                   # the convention of the cards is
                   # [col_min:row_min,col_max:row_max]
                   # inclusive
                   match = maskFormat.match(maskArea)
                   if match == None:
                       # unable to match mask area!
                       print '# WARNING: Extension', extension, 'unable to parse', maskArea
                       continue

                   group = map(int, match.groups())

                   # turn into a Footprint!
                   # we have to account for 2 things
                   #
                   # the regions in the fits file are 1-indexed
                   # while the footprints are 0-indexed
                   # therefore subtract 1 from the origin

                   maskBBox2i = afwImage.BBox2i(
                       group[0]-1,             # col min
                       group[1]-1,             # row min
                       group[2]-group[0]+1,    # col span
                       group[3]-group[1]+1     # row span
                   )
                   maskFootprint = detectionLib.FootprintPtrT( detectionLib.Footprint(maskBBox2i) )
                   footprintList.push_back(maskFootprint)

           # set all the bad masks at once
           detectionLib.setMaskFromFootprintList(mask, footprintList, badBitMask)
                
           # Write the Mask to disk
           imageExt = re.sub(imageRoot, '_%d_msk.fits' % (extension), imageRoot)
           #outMaskFileName = 'raw-' + mjdDate + '-' + filter + '-' + imageNumber + imageExt
       
           #outMaskFileName = re.sub('.fits', '_%d_msk.fits' % (ccdNum), inImage)
           outMaskFileName = re.sub('.fits', '_%d_msk.fits' % (extension), inImage)
           outMskDir = os.path.join(destinationDir,outMaskFileName);
           mask.writeFits(outMskDir)
           print 'Wrote: ', outMaskFileName
       
    inFile.close()
#  os.remove(inImage)
 
if __name__ == "__main__":
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)


        


