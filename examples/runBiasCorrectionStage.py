"""
@brief Example code to run the nightly ISR's 'Bias Correction Stage'.

@author Nicole M. Silvestri, University of Washington
 contact nms@astro.washington.edu

 file created Mon Nov 24, 2008 

@file

"""

if __name__ == "__main__":
    import eups
    import os
    import lsst.afw.image as afwImage
    import lsst.daf.base as dafBase
    import lsst.pex.logging as pexLog
    import lsst.pex.policy as pexPolicy
    import lsst.ip.isr.BiasCorrection as ipIsrBias 

    pexLog.Trace.setVerbosity("lsst.ip.isr", 4)
    
    dataDir = eups.productDir("afwdata")
    if not dataDir:
        raise RuntimeError("Must set up afwdata to run this program.")
    isrDir = eups.productDir("ip_isr")
    if not isrDir:
        raise RuntimeError("Must set up ip_isr to run this program.")
    
    chunkExposureInPath = os.path.join(dataDir, "overStageTestExposure_1")
    masterChunkExposureInPath = os.path.join(dataDir, "CFHT", "D4", "05Am05.bias.0.36.00_1_img.fits")
    isrPolicyPath = os.path.join(isrDir, "pipeline", "isrPolicy.paf")
    chunkExposureOutPath = os.path.join(dataDir, "biasStageTestExposure_1")
    
    memId0 = dafBase.Citizen_getNextMemId()

    chunkExposure = afwImage.ExposureD()
    chunkExposure.readFits(chunkExposureInPath)
    masterChunkExposure = afwImage.ExposureD()
    masterChunkExposure.readFits(masterChunkExposureInPath)

    isrPolicy = pexPolicy.Policy.createPolicy(isrPolicyPath)
    
    ipIsrBias.biasCorrection(chunkExposure, masterChunkExposure, isrPolicy)

    pexLog.Trace("lsst.ip.isr.biasCorrection", 4, "Writing chunkExposure to %s [_img|var|msk.fits]" % (chunkExposureOutPath,))
    chunkExposure.writeFits(chunkExposureOutPath)
    
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)
