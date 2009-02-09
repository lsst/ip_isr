# file createChunkExposure.py
#
# author Nicole Silvestri, University of Washington
#
# contact nms@astro.washington.edu
#
# This code is designed to be a pre-stage to the nightly IPP.  It
# takes as input a list of CCD MaskedImages (excluding the _img.fits
# extension) and splits them into one LSST Exposure per chunk
# requested of the CCD Image.  This 'Chunk' Exposure is then written
# to FitsStorage.  Currently, the code assumes that the chunk is
# defined as one amplifier but may be moodified to accept different
# chunk sizes.
#
# The code also takes an input a Policy file (policyFile.paf) with
# information specific to the data set being processed.  This file is
# used to map the dataset-specific header information to standard LSST
# metadata names.

import os
import re
import sys
import pyfits
import string
import numarray

import eups

import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy

dataDir = eups.productDir("afwdata")
if not dataDir:
    raise RuntimeError("You must set up afwdata to run this code.")
cfhtDataDir = os.path.join(dataDir, "/CFHT/D4/")
destinationDir = os.path.join(cfhtDataDir, "/tempfiles/")

def main():

    if len(sys.argv) != 3:
        print 'Usage : createChunkExposure.py inputList.txt policyFile.paf'
        sys.exit(1)

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
        inCcdMaskedImage = afwImage.MaskedImageF(image)
        print 'Read: ', image + '_img.fits'
        inCcdWCS = afwImage.Wcs(inCcdMaskedImage.getImage().getMetaData())
        inCcdExposure = afwImage.ExposureF(inCcdMaskedImage, inCcdWCS)
        inCcdMetadata = inCcdMaskedImage.getImage().getMetaData()

        inMiRows = inCcdMaskedImage.getRows()
        inMiCols = inCcdMaskedImage.getCols()
            
        # get the number of amplifiers in each row/col of the CCD Image
        numRowAmps = policy.getInt('numRowAmps')
        numColAmps = policy.getInt('numColAmps')

        # get the number of rows/cols in each amplifier image
        nRowSubImages = inMiRows/numRowAmps
        nColSubImages = inMiCols/numColAmps

        ampNum = 1
        for row in range(numRowAmps):
            for col in range(numColAmps):

                # Create the bounding box for each amplifier
                exposureBbox = afwImage.BBox2i(col * nColSubImages,
                                               row * nRowSubImages,
                                               nColSubImages,
                                               nRowSubImages)
                ccdNumField = policy.getString('ccdField')
                print 'CCD Number Field: ', ccdNumField
                ccdNumber = inCcdMetadata.findUnique(ccdNumField)
                ccdNum = ccdNumber.getValue()
                print 'CCD Number: ', ccdNum
                # Write the new amplifier exposure to disk
                outAmpExposure = inCcdExposure.getSubExposure(exposureBbox)
                outFileName = re.sub('_%d' % ccdNum, '_%d_%d' % (ccdNum, ampNum), image)
                outdir = os.path.join(destinationDir,outFileName);
                outAmpExposure.writeFits(outdir)
                print 'Wrote ', outFileName
                ampNum += 1

if __name__ == "__main__":
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)


        


