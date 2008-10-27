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

dataDir = eups.productDir("afwdata")
if not dataDir:
    raise RuntimeError("You must set up afwdata to run this code.")
destinationDir = '/astro/users/nms/code/isrtrunk/examples/tempfiles/'
    
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
        inCcdMaskedImage = afwImage.MaskedImageF()
        inCcdMaskedImage.readFits(image)
        print 'Read: ', image + '_img.fits'
        inCcdWCS = afwImage.Wcs(inCcdMaskedImage.getImage().getMetaData())
        inCcdExposure = afwImage.ExposureF(inCcdMaskedImage, inCcdWCS)
        inCcdMetadata = inCcdMaskedImage.getImage().getMetaData()
	ccdNumField = policy.getString('ccdField')
        print 'CCD Number Field: ', ccdNumField
        ccdNumber = inCcdMetadata.findUnique(ccdNumField)
        ccdNum = ccdNumber.getValue()
        print 'CCD Number: ', ccdNum

if __name__ == "__main__":
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)


