# file dataFormatConversion.py
#
# author Nicole Silvestri, University of Washington
#
# contact nms@astro.washington.edu
#
# This code is designed to process any format file (currently MEF or
# single FITS file) from any telescope/camera configuration and
# convert it into a standard LSST Exposure.  A Policy file with
# information specific to the data set is required for this stage.
#
# This is written following the DataFormatConversionStage proposal on
# Trac (/trac/wiki/DataFormatConversionStage).  To fully implement
# this, RO's Metsdata Proposal must be accepted
# (trac/wiki/ImageMetadataProposal)

import os
import re
import sys
import pyfits
import eups
import lsst.detection.detectionLib as detectionLib
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase

import lsst.ip.isr as ipIsr
import lsst.pex.logging as pexLogging
import lsst.pex.policy pexPolicy

import createMaskForCFHT as createMask

## dataDir = eups.productDir("afwdata")
## if not dataDir:
##     raise RuntimeError("You must set up afwdata to run these tests")

def main():

    # Created Policy file named "dataFormatConversionStageDictionary.paf"
    if len(sys.argv) != 3:
        print 'Usage : dataFormatConversion.py inFile.fits policyFile.paf'
        sys.exit(1)
       
    inImage = sys.argv[1]
    if not os.path.isfile(inImage):
        print 'ERROR: Cannot locate', inImage
        sys.exit(1)

    # afwImage.readfits can not open MEF files. Using PyFits...
    try:
        print 'Reading fits file with pyfits'
        inFile  = pyfits.open(inImage)
    except:
        print 'ERROR: Cannot open', inImage
        sys.exit(1)

    policyFile = sys.argv[2]
    policy = pexPolicy.Policy.createPolicy(policyFile)

    def verbosity = 0
	   
    if options.verbosity > 0:
        print "Verbosity =", options.verbosity
        pexLogging.Trace_setVerbosity("lsst.ip.isr", options.verbosity)
               
    nExt = len(inFile)
    print 'Number of Extensions in %s :' % inImage, nExt
    if (len > 1):
        print 'First extension [00] header will be concatenated with individual image headers.'
    
    # For now (07/21/08) extract keyword values needed for the ISR.  If
    # Metadata proposal is accepted, these values will become member
    # variables. Query policy file for gain and rdNoise values or keywords.

    numRowAmps = pexPolicy.Policy.getDouble(policyFile.numRowAmps)
    numColAmps = pexPolicy.Policy.getDouble(policyFile.numColAmps)

    for extension in range(1,nExt):
        inputImage = afwImage.ImageF(inFile[extension])

        nImageRows = inputImage.getRows()
        nImageCols = inputImage.getCols()
        
        # get the Mask Image - will need to check for masks on input eventually
##                 try:
##                     inputMaskImage = inputMaskedImage.getMask()
##                 else:          
                       # if no Mask Image, create one from the 'BAD'
                       # pixel info.
                    
        # Modified Andy Becker's code slightly to extract Masks from
        # CFHT images. This code will need to be generalized in the
        # future to accommodate other datasets.
        
        inputMask = createMask.getMask(inputImage, extension)

        nRowSubImages = nImageRows/numRowAmps
        nColSubImages = nImageCols/numColAmps
	
        for row in range(nRowAmps):
            for col in range(nColAmps):

                bbox = afwImage.BBox2i(col * nColSubImages,
                                       row * nRowSubImages,
                                       nColSubImages,
                                       nRowSubImages)
                
                header = inFile[extension].header
                
                # Get the member variable 'gain' from the header.
                # Determine if a keyword or numerical value exists in the Policy for gain.
                
                def is_a_number(pexPolicy.Policy.getDouble(policyFile.gain)):
                    try:
                        gain = double(pexPolicy.Policy.getDouble(policyFile.gain))
                    else:
                        gainKeyword = pexPolicy.Policy.getString(policyFile.gain)
                        gain = header[gainKeyword]
               
                print 'The gain for extension %d and amplifier %d :' % extension % col, gain

                # Get the member variable 'rdNoise' from the header.
                # Determine if a keyword or numerical value exists in the Policy for rdNoise. 

                def is_a_number(pexPolicy.Policy.getDouble(policyFile.rdNoise)):
                    try:
                        rdNoise = double(pexPolicy.Policy.getDouble(policyFile.rdNoise))
                    else:
                        rdNoiseKeyword = pexPolicy.Policy.getString(policyFile.rdNoise)
                        rdNoise = header[rdNoiseKeyword]
                
                print 'The rdNoise for extension %d and amplifier %d :' % extension % col, rdNoise

                # Get the rest of the required member variables from
                # the header while we are at it
                
                datasec = header['DATASEC']
                trimsec = header['FILTER']
                biassec = header['BIASSEC']
                exptime = header['EXPTIME']
                taiMjd = header['TAI']
                binning = header['BIN']
                exptype = header['EXPTYPE']

                # Remove all of the above member variables (gain,
                # rdNoise, datasec, trimsec, biassec, exptime) from
                # the header to avoid future duplication issues

                del header['DATASEC']
                del header['FILTER']
                del header['BIASSEC']
                del header['EXPTIME']
                del header['TAI']
                del header['BIN']
                del header['EXPTYPE']
                del header[gainKeyword]
                del header[rdNoiseKeyword]
                
                # Get the Variance image. - will need to check for variance on input eventually
##                 try:
##                     inputVarianceImage = inputMaskedImage.getVariance()
##                 else:
                    # if no variance image, create the variance image
                    # from the gain, rdnoise and pixel values of the
                    # input image
                afwImage.Image<ImageT> inputVarianceImage(nImageCols, nImageRows)
                inputVarianceImage = inputImage/gain + 1/gain**2 * rdNoise**2

                # Get the WCS information from the header.                
                # This extracts all of the metadata?! Shouldn't we
                # extract only the WCS keywords into the WCS object
                # here (eg. CRPix, etc)?
                
                inputWCS = afwImage.Wcs(inputImage.getMetaData())

                # Remove the WCS information from the header to avoid
                # future duplication issues.  Persistence needs to fix
                # this issue with WCS later.

                del header['PIXSCAL1']
                del header['PIXSCAL2']
                del header['CRPIX1']
                del header['CRPIX2']
                del header['CD1_1']
                del header['CD1_2']
                del header['CD2_1']
                del header['CD2_2']
                del header['RADECSYS']
                del header['CRVAL1']
                del header['CRVAL2']
                del header['CTYPE1']
                del header['CTYPE2']
                del header['EQUINOX']
                
                # For now, create the MaskedImage
                newInputMaskedImage = afwImage.MaskedImageF(inputImage, inputVaraince, inputMask)

                # Create the Exposure, using the MaskedImage
                inputExposure = afwImage.ExposureF(inputMaskedImage, inputWCS)
                
                # In the future...create the Exposure from the
                # individual Image, Varaince, Mask, member variables,
                # and extra metadata
                # inputExposure = afwImage.ExposureF(inputImage, inputVaraince, InputMask, 

                outputExposure = inputExposure.getSubExposure(bbox)
                               
                outFileName = re.sub('.fits', '_%d_%d_image.fits' % (extension, col), inFile)
                outputExposure.writeFits(outFileName)
                print 'Wrote ', outFileName

if __name__ == "__main__":
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)





