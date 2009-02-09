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

import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy

import createCfhtMask as createMask

dataDir = eups.productDir("afwdata")
if not dataDir:
    raise RuntimeError("You must set up afwdata to run these tests")

def main():

    # Created a policy file named 'pipeline/datasetPolicy.paf'.  Each
    # dataset will need a Policy File like this with information
    # specific to the input images.
    
    if len(sys.argv) != 3:
        print 'Usage : createChunkExposure.py inFile.fits datasetPolicy.paf'
        sys.exit(1)
       
    inImage = sys.argv[1]
    if not os.path.isfile(inImage):
        print 'ERROR: Cannot locate', inImage
        sys.exit(1)

    # afwImage.readfits can not open MEF files. Using PyFits...
    try:
        print 'Reading fits file.'
        inFile  = pyfits.open(inImage)
    except:
        print 'ERROR: Cannot open', inImage
        sys.exit(1)

    policyFile = sys.argv[2]
    policy = pexPolicy.Policy.createPolicy(policyFile)
               
    nExt = len(inFile)
    print 'Number of Extensions in %s :' % inImage, nExt
            
    # For now (07/21/08) extract keyword values needed for the ISR.
    # If Metadata proposal is accepted, these values will become
    # member variables. Query policy file for gain and rdNoise values
    # or keywords and the number of apms/row and amps/column.

    numRowAmps = policy.getInt('numRowAmps')
    numColAmps = policy.getInt('numColAmps')

    # First extension of an MEF is general header information which
    # will be concatenated with the individual image or amplifier
    # headers.
    
    generalHeader = inFile[0].header
    
    for extension in range(1,nExt):

        ampHeader = inFile[extension].header

        # Create temporary images so we can read them into a MaskedImage
        tempImgFileName = re.sub('.fits', '_%d_img.fits' % (extension), inImage)
        inFile[extension].writeto(tempImgFileName)
        print 'Wrote:', tempImgFileName
        nCol = ampHeader['NAXIS1']
        nRow = ampHeader['NAXIS2']

        # Create an empty temporary Variance Image that is the same size as the Image
        tempVarFileName = re.sub('.fits', '_%d_var.fits' % (extension), inImage)
        varImage = afwImage.ImageF
        tempVariance = varImage(nCol, nRow)
        tempVariance.writeFits(tempVarFileName)
        print 'Wrote:', tempVarFileName

        # Create an empty temporary Mask Image that is the same size as the Image 
        tempMskFileName = re.sub('.fits', '_%d_msk.fits' % (extension), inImage)
        mskImage = afwImage.MaskU
        tempMask = mskImage(nCol, nRow)
        tempMask.writeFits(tempMskFileName)
        print 'Wrote:', tempMskFileName
                                                         
        tempName = re.sub('.fits', '_%d' % (extension), inImage)
        inputMaskedImage = afwImage.MaskedImageF(tempName)
        
        nImageRows = inputMaskedImage.getRows()
        nImageCols = inputMaskedImage.getCols()

        nRowSubImages = nImageRows/numRowAmps
        nColSubImages = nImageCols/numColAmps
	
        for row in range(nRowAmps):
            for col in range(nColAmps):

                exposureBbox = afwImage.BBox2i(col * nColSubImages,
                                               row * nRowSubImages,
                                               nColSubImages,
                                               nRowSubImages)

                # Concatinate the first extension [00] of an MEF file
                # with the individual headers of each amplifier.
                
                if (nExt > 1):
                    print 'Concatinating first extension [00] header with individual amplifier headers.'
                    header = generalHeader.append(ampHeader)
                
                # Get the member variable 'gain' from the header.
                # Determine if a keyword or numerical value exists in the Policy for gain.

                gain = policy.getDouble(gain)
                
                if type(gain) == type(1.0):
                    print 'Numerical value for gain found in Policy:' , gain
                else:
                    gainKeyword = policy.getString('gain')
                    gain = header[gainKeyword]
                    print 'Gain recovered from FITS keyword.'
               
                print 'The gain for extension %d and amplifier %d :' % extension % col, gain

                # Get the member variable 'rdNoise' from the header.
                # Determine if a keyword or numerical value exists in the Policy for rdNoise. 

                rdNoise = policy.getDouble('rdNoise')
                if type(gain) == type(1.0):
                   print 'Numerical value for rdNoise found in Policy:', rdNoise
                else:
                    rdNoiseKeyword = policy.getString('rdNoise')
                    rdNoise = header[rdNoiseKeyword]
                    print 'RdNoise recovered from FITS keyword.'
                
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
                
                # Get or create the Variance Image.
                try:
                    inputVariance = inputMaskedImage.getVariance()
                except noVariance:
                    
                    # No Variance Image, create the variance image
                    # from the gain, rdnoise and pixel values of the
                    # input image
                    
                    # afwImage.Image<ImageT> inputVariance(nImageCols, nImageRows)
                    inputVariance = inputImage/gain + 1/gain**2 * rdNoise**2

                # Get or create the Chunk Mask Image.
                # Check the Policy for a bad pixel mask if a Mask does not already exist.
                try:
                    inputMask = inputMaskedImage.getMask()
                except noMask:    
                    bpmAvailable = policy.getDouble('bpm')
                    if (bpmAvailable == 0):
            
                        # No Bad Pixel Mask, create one from the 'BAD' pixel
                        # info in the image header cards.
                    
                        # Modified Andy Becker's code slightly to extract Masks from
                        # CFHT images. This code will need to be generalized in the
                        # future to accommodate other datasets.
        
                        inputMask = createMask.getMask(inputImage, extension)
                    else:
                        bpmName = policy.getString('bpmName')
                        badPixelMask = afwImage.MaskF(bpmName)
                        ##  maskPlaneDict = badPixelMask.getMaskPlaneDict()
                        ##  inputMask = afwImage.MaskF(nImageCols, nImageRows, maskPlaneDict)
                        inputMask = badPixelMask.copy()
                    
                # Get the WCS information from the header.                
                # This extracts all of the metadata?! We should
                # extract only the WCS keywords into the WCS object
                # here (eg. CRPix, etc)...
                
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
                
                # For now, create the new MaskedImage
                newInputMaskedImage = afwImage.MaskedImageF(inputImage, inputVaraince, inputMask)

                # Create the Exposure, using the MaskedImage
                inputExposure = afwImage.ExposureF(newInputMaskedImage, inputWCS)
                
                # In the future...maybe create the Exposure from the
                # individual Image, Varaince, Mask, member variables,
                # and extra metadata:
                #inputExposure = afwImage.ExposureF(inputImage, inputVaraince, inputMask, datasec,filter, biassec, tiaMjd, exptime, binning, exptype, gain, rdnoise, inputWcs, header)

                ampExposure = inputExposure.getSubExposure(exposureBbox)

                # Need to remove overscan strip region(s) from the
                # amplifier image and save them as individual images

                biassecFormat = re.compile('^\[(\d+):(\d+),(\d+):(\d+)\]$')
                match = biassecFormat.match(biassec)
                if match == None:
                    # unable to match biassec area!
                    print '# WARNING: Extn', nExt, 'unable to parse', biassec
                    continue
                
                group = map(int, match.groups())
                overscanBbox = afwImage.BBox2i(group[0],group[1],             
                                               group[2]-group[0],    
                                               group[3]-group[1])
                overscanExposure = ampExposure.getSubExposure(overscanBbox)

                # Create a subimage of the Chunk Exposure that
                # excludes the overscan region and write out the final
                # Chunk Exposure that should now represent on
                # amplifier of the CCD and have no overscan region.
                
                outputExposureBbox = afwImage.BBox2i((col * nColSubImages)-group[0],
                                                     (row * nRowSubImages)-group[1],
                                                     nColSubImages-(group[2]-group[0]),
                                                     nRowSubImages-(group[3]-group[1]))
                outputExposure = ampExposure.getSubExposure(outputExposureBbox)

                outFileName = re.sub('.fits', '_%d_%d_img.fits' % (extension, col), inImage)
                outputExposure.writeFits(outFileName)
                print 'Wrote ', outFileName

                # Clean up intermediate files
                os.remove(tempImgFileName)
                os.remove(tempVarFileName)
                os.remove(tempMskFileName)


if __name__ == "__main__":
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)





