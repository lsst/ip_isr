"""
@brief Example code to run the nightly ISR Bias Correction Stage.

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
import lsst.ip.isr.BiasCorrection as ipIsrBias
from lsst.ip.isr.VerifyMasterFile import VerifyMasterFile

def main():

    pexLog.Trace.setVerbosity("lsst.ip.isr", 4)
    
    dataDir = eups.productDir("isrdata")
    if not dataDir:
        raise RuntimeError("Must set up isrdata to run this program.")
    isrDir = eups.productDir("ip_isr")
    if not isrDir:
        raise RuntimeError("Must set up ip_isr to run this program.")
    
    chunkExposureInPath       = os.path.join(dataDir, "overStageTestExposure_1")
    masterChunkExposureInPath = os.path.join(dataDir, "CFHT", "D4", "05Am05.bias.0.36.00_1")
    isrPolicyPath             = os.path.join(isrDir, "pipeline", "isrPolicy.paf")
    chunkExposureOutPath      = os.path.join(dataDir, "biasStageTestExposure_1")

    chunkExposure       = afwImage.ExposureF(chunkExposureInPath)
    masterChunkExposure = afwImage.ExposureF(masterChunkExposureInPath)
    isrPolicy           = pexPolicy.Policy.createPolicy(isrPolicyPath)

    if not VerifyMasterFile(chunkExposure, masterChunkExposure, 'lsst.ip.isr.biasCorrection', isrPolicy.getPolicy("biasPolicy")):
        raise pexExcept.LsstException, "Can not verify Master file for Bias correction"
    
    ipIsrBias.biasCorrection(chunkExposure, masterChunkExposure, isrPolicy)

    pexLog.Trace("lsst.ip.isr.biasCorrection", 4, "Writing chunkExposure to %s [_img|var|msk.fits]" % (chunkExposureOutPath,))
    chunkExposure.writeFits(chunkExposureOutPath)

if __name__ == "__main__":
   
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)
