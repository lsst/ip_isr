"""
Test ISR sub-stage, 'Saturation Correction for Chunk Exposure'.
"""

import os
import math
import pdb  # we may want to say pdb.set_trace()
import unittest

import numpy

import eups

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.utils.tests as utilsTests
import lsst.pex.logging as pexLog
import lsst.pex.exceptions as pexEx
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as isrTrunk

Verbosity = 0 # increase to see trace
pexLog.Trace_setVerbosity("lsst.ip.isr", Verbosity)

dataDir = eups.productDir("afwdata/trunk/CFHT/D4/")
if not dataDir:
    raise RuntimeError("Must set up afwdata to run these tests!")

## INPUT IMAGE AND PATH NAMES

InputChunkExposure = "raw-53535-i-797722_1" # The Chunk Exposure to be calibrated

InputDatasetPolicy = "cfhtDataPolicy.paf"   # Contains information specific to the data set being processed
InputIsrPolicy =  "isrPolicy.paf"           # Contains information directing ISR pipeline behavior
SatLookupTable = "satLookUpTable"           # Use Lookup Table in lieu of a function

currDir = os.path.abspath(os.path.dirname(__file__))
isrDir = "~/code/isrtrunk/src/"
inFilePath = os.path.join(dataDir, InputChunkExposure)

datasetPolicyPath = os.path.join(isrDir, InputDatasetPolicy)
isrPolicyPath = os.path.join(isrDir, InputIsrPolicy)
satLookUpTablePath = (isrDir, SatLookupTable)

## OUTPUT IMAGE AND PATH NAMES

OutputSatName = "satCorExposure"

outSatPath = os.path.join(dataDir, OutputSatName)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class isrTestCases(unittest.TestCase):
    """
    Tests for the ISR sub-stage, 'Saturation Correction for Chunk Exposure'.
    """
    def setUp(self):
       
        chunkExposure = afwImage.ExposureF()
        chunkExposure.readFits(inFilePath)

        isrPolicy = pexPolicy.Policy.createPolicy(isrPolicyPath)
        datasetPolicy = pexPolicy.Policy.createPolicy(datasetPolicyPath)

    def tearDown(self):

        del self.chunkExposure
        del self.isrPolicy
        del self.datasetPolicy

    def testSaturationCorrection(self):
        
        satCor = isrTrunk.saturationCorrectionForChunkExposure(chunkExposure, isrPolicy, datasetPolicy)

        chunkExposure.writeFits(outSatPath)
        
        #satCor = isrTrunk.saturationCorrectionForChunkExposure(chunkExposure, isrPolicy, datasetPolicy, satLookupTablePath)
        
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
