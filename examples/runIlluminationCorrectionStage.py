"""
@brief Example code to run the Data Release and nightly ISR
Illumination Correction Stage.

@author Nicole M. Silvestri
        University of Washington
        nms@astro.washington.edu

 file created Tue Nov 25, 2008 

@file

"""
import eups
import os
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr.IlluminationCorrection as ipIsrIllum

def main():

    #Set verbosity - increase from zero to see messages
    pexLog.Trace.setVerbosity("lsst.ip.isr", 4)

    #Check for data and ISR directories
    dataDir = eups.productDir("isrdata")
    if not dataDir:
        raise RuntimeError("Must set up isrdata to run this program.")
    isrDir = eups.productDir("ip_isr")
    if not isrDir:
        raise RuntimeError("Must set up ip_isr to run this program.")

    #INPUTS:
    masterChunkExposureInPath = os.path.join(dataDir, "CFHT", "D4", "05Am05.flat.0.36.00_1_img.fits")  
    # going to have to make the night sky flat
    masterSFChunkExposureInPath = os.path.join(dataDir, "CFHT", "D4", "05Am05.NSflat.0.36.00_1_img.fits")
    isrPolicyPath = os.path.join(isrDir, "pipeline", "isrPolicy.paf")
    masterIcpChunkExposureInPath = os.path.join(dataDir, "CFHT", "D4", "05Am05.IC.0.36.00_1_img.fits") ##WRONG NAME
    masterDfpChunkExposureInPath = os.path.join(dataDir, "CFHT", "D4", "05Am05.flat.0.36.00_1_img.fits")

    #OUTPUTS:
    masterChunkExposureOutPath = os.path.join(dataDir, "IllumStageTestExposure_1")
    masterChunkExposureOutPathDR = os.path.join(dataDir, "IllumStageDRTestExposure_1")
   
    masterChunkExposure = afwImage.ExposureD(masterChunkExposureInPath)
    masterSFChunkExposure = afwImage.ExposureD(masterSFChunkExposureInPath)
    masterIcpChunkExposure = afwImage.ExposureD(masterIcpChunkExposureInPath)
    masterDfpChunkExposure = afwImage.ExposureD(masterDfpChunkExposureInPath)

    isrPolicy = pexPolicy.Policy.createPolicy(isrPolicyPath)

    # Call the Data Release version
    ipIsrIllum.illuminationCorrectionDR(masterChunkExposure, masterSFChunkExposure, isrPolicy)
    pexLog.Trace("lsst.ip.isr.illuminationCorrection", 4, "Writing chunkExposure to %s [_img|var|msk.fits]" % (masterChunkExposurOutPathDR,))
    masterChunkExposure.writeFits(masterChunkExposureOutPathDR)

    # Call the nightly processing version
    ipIsrIllum.illuminationCorrection(masterChunkExposure, masterDfpChunkExposure, masterIcpChunkExposure, isrPolicy)
    pexLog.Trace("lsst.ip.isr.illuminationCorrection", 4, "Writing chunkExposure to %s [_img|var|msk.fits]" % (masterChunkExposureOutPath,))
    masterChunkExposure.writeFits(masterChunkExposureOutPath)
    
if __name__ == "__main__":
    
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)
