"""
Purpose: A simple test for the ISR sub-stage,
         'Overscan Correct and Trim Chunk Exposure'.

Author: Nicole M. Silvestri,
        University of Washington
        nms@astro.washington.edu
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

currDir = os.path.abspath(os.path.dirname(__file__))
inFilePath = os.path.join(dataDir, "CFHT", "D4", InputChunkExposure)

datasetPolicyPath = os.path.join(currDir, "../pipeline", InputDatasetPolicy)
isrPolicyPath = os.path.join(currDir, "../pipeline", InputIsrPolicy)

## OUTPUT IMAGE AND PATH NAMES

OutputName = "testOverCorExposure"

outputPath = os.path.join(dataDir, OutputName)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class isrTestCases(unittest.TestCase):
    """
    Tests for the ISR stage, 'Overscan Correct and Trim Chunk Exposure'.
    """
    def setUp(self):
        self.chunkExposure = afwImage.ExposureF()
        self.chunkExposure.readFits(inFilePath)
        self.isrPolicy = pexPolicy.Policy.createPolicy(isrPolicyPath)
        self.datasetPolicy = pexPolicy.Policy.createPolicy(datasetPolicyPath)

    def tearDown(self):
        del self.chunkExposure
        del self.isrPolicy
        del self.datasetPolicy

    def testOverscanTrimCorrection(self):
        
        ipIsr.overscanCorrectAndTrimChunkExposure(self.chunkExposure, self.isrPolicy, self.datasetPolicy)

        self.chunkExposure.writeFits(outputPath)
        
       
        
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
