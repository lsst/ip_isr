"""
Test ISR sub-stages.
"""

import os
import math
import pdb                                   # we may want to say pdb.set_trace()
import unittest

import numpy

import eups

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.utils.tests as utilsTests
import lsst.pex.logging as pexLog
import lsst.pex.exceptions as pexEx
import lsst.pex.policy as pexPolicy
import lsst.ip.isr.isrLib as isrTrunk

Verbosity = 0 # increase to see trace
pexLog.Trace_setVerbosity("lsst.ip.isr", Verbosity)

dataDir = eups.productDir("afwdata/trunk/CFHT/D4/")
if not dataDir:
    raise RuntimeError("Must set up afwdata to run these tests!")

## INPUT IMAGE AND PATH NAMES

InputChunkExposure = "raw-53535-i-797722_1" # The Chunk Exposure to be calibrated

InputBiasName = "05Am05.bias.0.36.00_1"     # Master Calibration Image Names
#InputDarkName = ""
InputFlatName = "05Am05.flat.i.36.01_1"
InputFringeName = "05Am05.fringe.i.36.00_1"
#InputIllumName = ""
#InputPupilName = ""

InputDatasetPolicy = "cfhtDataPolicy.paf"   # Contains information specific to the data set being processed
InputIsrPolicy =  "isrPolicy.paf"           # Contains information directing ISR pipeline behavior
SatLookupTable = "satLookUpTable"           # Use Lookup Table in lieu of a function
LinLookupTable = "linLookupTable"           # Use Lookup Table in lieu of a function
#FringeFile = "fringeFile"                  # For scaling the fringe image

currDir = os.path.abspath(os.path.dirname(__file__))
isrDir = "~/code/isrtrunk/src/"
inFilePath = os.path.join(dataDir, InputChunkExposure)
biasPath = os.path.join(dataDir, InputBiasName)
#darkPath = os.path.join(dataDir, InputDarkName)
flatPath = os.path.join(dataDir, InputFlatName)
fringePath = os.path.join(dataDir, InputFringeName)
#illumPath = os.path.join(dataDir, InputIllumName)
#pupilPath = os.path.join(dataDir, InputPupilName)

datasetPolicyPath = os.path.join(isrDir, InputDatasetPolicy)
isrPolicyPath = os.path.join(isrDir, InputIsrPolicy)
satLookUpTablePath = (isrDir, SatLookupTable)
linLookupTablePath = (isrDir, LinLookupTable)
#fringeFilePath = (isrDir, FringeFile)

## OUTPUT IMAGE AND PATH NAMES

OutputSatName = "satCorExposure"
OutputOverName = "overCorName"
OutputBiasName = "biasCorExposure"

outSatPath = os.path.join(dataDir, OutputSatName)
outOverPath = os.path.join(dataDir, OutputOverName)
outBiasPath = os.path.join(dataDir, OutputBiasName)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class isrTestCases(unittest.TestCase):
    """
    Tests for the ISR sub-stages.
    """
    def setUp(self):
       
        chunkExposure = afwImage.ExposureF(inFilePath)
        biasExposure = afwImage.ExposureF(biasPath)
        #darkExposure = afwImage.ExposureF(darkPath)
        flatExposure = afwImage.ExposureF(flatPath)
        fringeExposure = afwImage.ExposureF(fringePath)
        illumExposure = afwImage.ExposureF(illumPath)
        pupilExposure = afwImage.ExposureF(pupilPath)

        isrPolicy = pexPolicy.Policy.createPolicy(isrPolicyPath)
        datasetPolicy = pexPolicy.Policy.createPolicy(datasetPolicyPath)

    def tearDown(self):

        del self.chunkExposure
        del self.biasExposure
        del self.darkExposure
        del self.flatExpoure
        del self.fringeExposure
        del self.illumExposure
        del self.pupilExposure
        del self.isrPolicy
        del self.datasetPolicy

    def testSaturationCorrection(self):
        
        satCor = isrTrunk.saturationCorrectionForChunkExposure(chunkExposure, isrPolicy, datasetPolicy)

        chunkExposure.writeFits(outSatPath)
        
        #satCor = isrTrunk.saturationCorrectionForChunkExposure(chunkExposure, isrPolicy, datasetPolicy, satLookupTablePath)

    def testOverscanTrimCorrection(self):
        
        overscanCor = isrTrunk.overscanCorrectAndTrimChunkExposure(chunkExposure, isrPolicy, datasetPolicy)

##     def testBiasCorrection(self):
        
##         biasCor = isrTrunk.biasCorrectChunkExposure(chunkExposure, biasExposure, isrPolicy, datasetPolicy)

##     def xtestDarkCurrentCorrection(self):
        
##         #darkCor = isrTrunk.darkCurrentCorrectChunkExposure(chunkExposure, darkExposure, isrPolicy, datasetPolicy)

##     def xtestLinearizationCorrection(self):
        
##         linearize = isrTrunk.linearizeChunkExposure(chunkExposure, isrPolicy, datasetPolicy, linLookupTablePath)

##     def xtestFlatFieldCorrection(self):

##         flatCor = isrTrunk.flatFieldCorrectChunkExposure(chunkExposure, flatExposure, isrPolicy, datasetPolicy)

##     def xtestIlluminationCorrection(self):
        
##         #illumCor = isrTrunk.illuminationCorrection(chunkExposure, illumExposure, isrPolicy, datasetPolicy)

##     def xtestFringeCorrection(self):

##         #fringeCor = isrTrunk.defringeChunkExposure(chunkExposure, fringeExposure, isrPolicy, datasetPolicy)
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """
    Returns a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(isrTestCases)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
	   

if __name__ == "__main__":
    #   utilsTests.run(suite())
    run(True)
