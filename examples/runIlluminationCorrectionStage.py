"""
@brief Example code to run the nightly ISR's 'Illumination Correction Stage'.

@author Nicole M. Silvestri
        University of Washington
        nms@astro.washington.edu

 file created Tue Nov 25, 2008 

@file

"""

if __name__ == "__main__":
    import eups
    import os
    import lsst.afw.image as afwImage
    import lsst.daf.base as dafBase
    import lsst.pex.logging as pexLog
    import lsst.pex.policy as pexPolicy
    import lsst.ip.isr.IlluminationCorrection as ipIsrIllum 

    pexLog.Trace.setVerbosity("lsst.ip.isr", 4)
    
    dataDir = eups.productDir("afwdata")
    if not dataDir:
        raise RuntimeError("Must set up afwdata to run this program.")
    isrDir = eups.productDir("ip_isr")
    if not isrDir:
        raise RuntimeError("Must set up ip_isr to run this program.")
    
    chunkExposureInPath = os.path.join(dataDir, "linearStageTestExposure_1")
    masterChunkExposureInPath = os.path.join(dataDir, "CFHT", "D4", "05Am05.flat.0.36.00_1_img.fits")
    # going to have to make the night sky flat
    masterSFChunkExposureInPath = os.path.join(dataDir, "CFHT", "D4", "05Am05.NSflat.0.36.00_1_img.fits")
    isrPolicyPath = os.path.join(isrDir, "pipeline", "isrPolicy.paf")
    chunkExposureOutPath = os.path.join(dataDir, "IllumStageTestExposure_1")
    
    memId0 = dafBase.Citizen_getNextMemId()

    masterChunkExposure = afwImage.ExposureD()
    masterChunkExposure.readFits(masterChunkExposureInPath)
    masterSFChunkExposure = afwImage.ExposureD()
    masterSFChunkExposure.readFits(masterSFChunkExposureInPath) 

    isrPolicy = pexPolicy.Policy.createPolicy(isrPolicyPath)
    
    ipIsrIllum.illuminationCorrection(masterChunkExposure, masterSFChunkExposure, isrPolicy)

    pexLog.Trace("lsst.ip.isr.illuminationCorrection", 4, "Writing chunkExposure to %s [_img|var|msk.fits]" % (chunkExposureOutPath,))
    chunkExposure.writeFits(masterChunkExposureOutPath)
    
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)
