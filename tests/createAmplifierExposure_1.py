# file createAmplifierImage.py
#
# author Nicole Silvestri, University of Washington
#
# contact nms@astro.washington.edu
#
# This code is designed to be a pre-stage to the nightly
# IPP.  It will process a file (currently multi or single extension FITS file)
# from any telescope/camera configuration and break it into its
# individual amplifier images.  The amplifier images are then written
# to FitsStorage.
#
# A Policy file with information specific to the data set is required
# for this stage.
#
import os
import re
import sys
import pyfits
import string
import numarray

import eups

import lsst.detection.detectionLib as detectionLib
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy

## dataDir = eups.productDir("afwdata")
## if not dataDir:
##     raise RuntimeError("You must set up DMS/testdata/afwdata to run this code.")
     
def main():

    destinationDir = '/astro/users/nms/code/isrtrunk/examples/'

    if len(sys.argv) != 3:
        print 'Usage : createAmplifierImage.py inputFile.fits badPixelMaskFile.fits'
        sys.exit(1)
       
    ccdNumList = (12, 13, 21, 22)

    inImage = sys.argv[1]
    if not os.path.isfile(inImage):
        print 'ERROR: Cannot locate', inImage
        sys.exit(1)

    # afwImage.readfits can not open MEF files. Using PyFits...
    try:
        print 'Reading Input Image Fits File.'
        inFile  = pyfits.open(inImage)
    except:
        print 'ERROR: Cannot open', inImage
        sys.exit(1)

    badPixFile = sys.argv[2]
    if not os.path.isfile(badPixFile):
        print 'ERROR: Cannot locate', badPixFile
        sys.exit(1)
 
    try:
        print 'Reading Bad Pixel Mask Fits File.'
        badPixelMask = pyfits.open(badPixFile)
    except:
        print 'ERROR: Cannot open', badPixFile
        sys.exit(1)
     
    nExt = len(inFile)
    print 'Number of CCDs in %s :' % inImage, nExt-1

    nExtMask = len(badPixelMask)
    print 'Number of CCDs in %s :' %badPixelMask, nExtMask-1

    # First extension [00] of an MEF is general header information
    # which will be concatenated with the individual image or
    # headers.
    
    generalHeader = inFile[0].header
    
    #for extension in range(1,nExt):
    for extension in range(ccdNumList+1):  # Images Russell wants for Kaiser co-add tests
        ccdHeader = inFile[extension].header

 #       print 'Concatinating first extension [00] header with individual CCD headers.'
 #       ccdHeader = inFile[0].header.append(ccdHeader)


        # Write the individual Science Images to disk
    
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

        #newVariance = pyfits.PrimaryHDU(newVarFile, ccdHeader)
        outVarianceFileName = re.sub('.fits', '_%d_var.fits' % (extension), inImage)
        outVarDir = os.path.join(destinationDir,outVarianceFileName)
        #newVaraince.writeto(outVarDir)    
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
                    group[3]-group[1]+1
                )    # row span
                maskFootprint = detectionLib.FootprintPtrT( detectionLib.Footprint(maskBBox2i) )
                footprintList.push_back(maskFootprint)

        # set all the bad masks at once!
        detectionLib.setMaskFromFootprintList(mask, footprintList, badBitMask)
                
        # Write the Mask to disk
        outMaskFileName = re.sub('.fits', '_%d_msk.fits' % (extension), inImage)
        outMskDir = os.path.join(destinationDir,outMaskFileName);
        mask.writeFits(outMskDir)
        print 'Wrote: ',outMaskFileName
        
    inName = re.sub('.fits', '_%d' % (extension), inImage)
    inCcdMaskedImage = afwImage.MaskedImageF()
    inCcdMaskedImage.readFits(inName)
    inCcdWCS = afwImage.Wcs(inCcdMaskedImage.getImage().getMetaData())
    inCcdExposure = afwImage.ExposureF(inCcdMaskedImage, inCcdWCS)

    inMiRows = inCcdMaskedImage.getRows()
    inMiCols = inCcdMaskedImage.getCols()
    
    policyFile = sys.argv[2]
    policy = pexPolicy.Policy.createPolicy(policyFile)

    # get the number of amplifiers in each row/col of the CCD Image           
    numRowAmps = policy.getInt('numRowAmps')
    numColAmps = policy.getInt('numColAmps')

    # get the number of rows/cols/ in each amplifier image
    nRowSubImages = inMiRows/numRowAmps
    nColSubImages = inMiCols/numColAmps

    amp = 1
    for row in range(numRowAmps):
        for col in range(numColAmps):

            # Create the bounding box for each amplifier
            exposureBbox = afwImage.BBox2i(col * nColSubImages,
                                           row * nRowSubImages,
                                           nColSubImages,
                                           nRowSubImages)

            # Write the new amplifier exposure to disk
            outAmpExposure = inCcdExposure.getSubExposure(exposureBbox)
            outFileName = re.sub('.fits', '_%d_%d_img.fits' % (extension, amp), inImage)
            outdir = os.path.join(destinationDir,outFileName);
            outAmpExposure.writeFits(outdir)
            print 'Wrote ', outFileName
            amp += 1

if __name__ == "__main__":
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)


        


