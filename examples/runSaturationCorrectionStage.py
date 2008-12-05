"""
@brief Example code to run the nightly ISR Saturation Correction Stage.

@author Nicole M. Silvestri
        University of Washington
        nms@astro.washington.edu

        created: Mon Nov 24, 2008 

@file

"""
import eups
import os
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr
import lsst.ip.isr.SaturationCorrection as ipIsrSat 

def main():

    pexLog.Trace.setVerbosity("lsst.ip.isr", 4)
    
    dataDir = eups.productDir("afwdata")
    if not dataDir:
        raise RuntimeError("Must set up afwdata to run this program.")
    isrDir = eups.productDir("ip_isr")
    if not isrDir:
        raise RuntimeError("Must set up ip_isr to run this program.")
    
    chunkExposureInPath = os.path.join(dataDir, "CFHT", "D4", "raw-53535-i-797722_1")
    isrPolicyPath = os.path.join(isrDir, "pipeline", "isrPolicy.paf")
    lookupTablePath = os.path.join(isrDir, "pipeline", "saturationLookupTable.tab")
    chunkExposureOutPath = os.path.join(dataDir, "satStageTestExposure_1")
    
    chunkExposure = afwImage.ExposureD()
    chunkExposure.readFits(chunkExposureInPath)

    isrPolicy = pexPolicy.Policy.createPolicy(isrPolicyPath)
    
    chunkMaskedImage = chunkExposure.getMaskedImage()
    numCols = chunkMaskedImage.getCols()
    numRows = chunkMaskedImage.getRows()
    numPixels = numCols * numRows
        
    lookupTable = open(lookupTablePath, "rU")  
    pixelValues = lookupTable.readlines()
    numPix = len(pixelValues)
    print 'Number of pixels: ', numPix
    for pixels in pixelValues:
        # strip trailing whitespace, returns, etc.
        pixels = pixels.strip()
        # ignore blank lines
        if not pixels:
            continue
        # ignore comment lines
        if pixels.startswith("#"):
            continue
        lookupList = pixels.split()
        if len(lookupList) < numPix or len(lookupList) > numPix:
            print "Cannot parse: ", pixels

    ipIsrSat.saturationCorrection(chunkExposure, isrPolicy, lookupList)

    pexLog.Trace("lsst.ip.isr.saturationCorrection", 4, "Writing chunkExposure to %s [_img|var|msk.fits]" % (chunkExposureOutPath,))
    chunkExposure.writeFits(chunkExposureOutPath)    

if __name__ == "__main__":
   
    memId0 = dafBase.Citizen_getNextMemId()
    main()
   
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)
