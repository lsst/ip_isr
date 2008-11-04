"""
Purpose: A simple test for the ISR sub-stage,
         'Saturation Correction for Chunk Exposure'.

Author: Nicole M. Silvestri,
        University of Washington
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
import lsst.afw.image.testUtils as imUtilsTests
import lsst.pex.logging as pexLog
import lsst.pex.exceptions as pexEx
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr

Verbosity = 0 # increase to see trace
pexLog.Trace_setVerbosity("lsst.ip.isr", Verbosity)

dataDir = eups.productDir("afwdata")
if not dataDir:
    raise RuntimeError("Must set up afwdata to run these tests!")

## INPUT IMAGE AND PATH NAMES

InputChunkExposure = "raw-53535-i-797722_1" # The Chunk Exposure to be calibrated

InputDatasetPolicy = "cfhtDataPolicy.paf"   # Contains information specific to the data set being processed
InputIsrPolicy =  "isrPolicy.paf"           # Contains information directing ISR pipeline behavior
SatLookupTable = "satLookUpTable"           # Use Lookup Table in lieu of a function

currDir = os.path.abspath(os.path.dirname(__file__))
isrDir = "../src/"
pafDir = "../pipeline/"
inFilePath = os.path.join(dataDir, "CFHT", "D4", InputChunkExposure)

datasetPolicyPath = os.path.join(pafDir, InputDatasetPolicy)
isrPolicyPath = os.path.join(pafDir, InputIsrPolicy)
satLookUpTablePath = (pafDir, SatLookupTable)

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
##         isrPolicy = pexPolicy.Policy(isrPolicyPath)
##         datasetPolicy = pexPolicy.Policy(datasetPolicyPath)

        isrPolicy = pexPolicy.Policy.createPolicy(isrPolicyPath)
        datasetPolicy = pexPolicy.Policy.createPolicy(datasetPolicyPath)

    def tearDown(self):

        del self.chunkExposure
        del self.isrPolicy
        del self.datasetPolicy

    def testSaturationCorrection(self):
        
        satCor = ipIsr.saturationCorrectionForChunkExposure(chunkExposure, isrPolicy, datasetPolicy)

        chunkExposure.writeFits(outSatPath)
        
        #satCor = ipIsr.saturationCorrectionForChunkExposure(chunkExposure, isrPolicy, datasetPolicy, satLookupTablePath)
        
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
