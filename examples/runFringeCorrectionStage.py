"""
@brief Example code to run the nightly ISR's 'Fringe Correction Stage'.

@author Nicole M. Silvestri
        University of Washington
        nms@astro.washington.edu

        created: Tue Nov 25, 2008 

@file

"""
import eups
import os
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.ip.isr.FringeCorrection as ipIsrFringe

def main():

    # Set verbosity - increase from zero to see messages
    pexLog.Trace.setVerbosity("lsst.ip.isr", 4)

    # Verify that data and ISR directories are setup
    dataDir = eups.productDir("afwdata")
    if not dataDir:
        raise RuntimeError("Must set up afwdata to run this program.")
    isrDir = eups.productDir("ip_isr")
    if not isrDir:
        raise RuntimeError("Must set up ip_isr to run this program.")

    #INPUTS
    chunkExposureInPath = os.path.join(dataDir, "FlatStageTestExposure_1")
    masterChunkExposureInPath = os.path.join(dataDir, "CFHT", "D4", "05Am05.fringe.i.36.01_1")
    isrPolicyPath = os.path.join(isrDir, "pipeline", "isrPolicy.paf")

    #OUTPUTS
    chunkExposureOutPath = os.path.join(dataDir, "FringeStageTestExposure_1")

    chunkExposure = afwImage.ExposureD()
    chunkExposure.readFits(chunkExposureInPath)
    masterChunkExposure = afwImage.ExposureD()
    masterChunkExposure.readFits(masterChunkExposureInPath)

    isrPolicy = pexPolicy.Policy.createPolicy(isrPolicyPath)

    # run it
    ipIsrFringe.fringeCorrection(chunkExposure, masterChunkExposure, isrPolicy)

    # write it
    pexLog.Trace("lsst.ip.isr.fringeCorrection", 4, "Writing chunkExposure to %s [_img|var|msk.fits]" % (chunkExposureOutPath,))
    chunkExposure.writeFits(chunkExposureOutPath)

if __name__ == "__main__":
    
    memId0 = dafBase.Citizen_getNextMemId()
    main()
   
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)
