#!/usr/bin/env python
import os

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.isr', Verbosity)

isrDataDir = eups.productDir('isrdata')
isrDir     = eups.productDir('ip_isr')

# The Chunk Exposure to be calibrated
InputExposure  = os.path.join(isrDataDir, 'CFHT/D4', 'raw-53535-i-797722_1')

# Master Calibration Image Names
InputBias      = os.path.join(isrDataDir, 'CFHT/D4', '05Am05.bias.0.36.00_1')    
InputFlat      = os.path.join(isrDataDir, 'CFHT/D4', '05Am05.flat.i.36.01_1')
InputFringe    = os.path.join(isrDataDir, 'CFHT/D4', '05Am05.fringe.i.36.00_1')

# Policy file
InputIsrPolicy = os.path.join(isrDir, 'pipeline', 'isrPolicy.paf')

class IsrTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = pexPolicy.Policy.createPolicy(InputIsrPolicy)
        
    def tearDown(self):
        del self.policy

    def doFlat(self, scaling):
        flatScaleKeyword = self.policy.getString('flatPolicy.flatScaleKeyword')
        filenameKeyword  = self.policy.getString('filenameKeyword')
        
        mi = afwImage.MaskedImageF(10,10)
        mi.getImage().set(10)
        exposure = afwImage.ExposureF(mi)
        metadata = exposure.getMetadata()
        metadata.setDouble(flatScaleKeyword, 1.)
        
        flat = afwImage.MaskedImageF(10,10)
        flat.getImage().set(1)
        flatexposure = afwImage.ExposureF(flat)
        dmetadata = flatexposure.getMetadata()
        dmetadata.setDouble(flatScaleKeyword, scaling)
        dmetadata.setString(filenameKeyword, 'Unittest Flat')

        ipIsr.FlatCorrection(exposure, flatexposure, self.policy)

        height        = mi.getHeight()
        width         = mi.getWidth()
        for j in range(height):
            for i in range(width):
                self.assertAlmostEqual(mi.getImage().get(i,j), 10 / (1./scaling), 5)

    def testFlat1(self):
        self.doFlat(scaling=10)

    def testFlat2(self):
        self.doFlat(scaling=0.1)

    def testFlat3(self):
        self.doFlat(scaling=3.7)

    def doIllum(self, scaling):
        illumScaleKeyword = self.policy.getString('illuminationPolicy.illumScaleKeyword')
        filenameKeyword  = self.policy.getString('filenameKeyword')
        
        mi = afwImage.MaskedImageF(10,10)
        mi.getImage().set(10)
        exposure = afwImage.ExposureF(mi)
        metadata = exposure.getMetadata()
        metadata.setDouble(illumScaleKeyword, 1.)
        
        illum = afwImage.MaskedImageF(10,10)
        illum.getImage().set(1)
        illumexposure = afwImage.ExposureF(illum)
        dmetadata = illumexposure.getMetadata()
        dmetadata.setDouble(illumScaleKeyword, scaling)
        dmetadata.setString(filenameKeyword, 'Unittest Illum')

        ipIsr.IlluminationCorrection(exposure, illumexposure, self.policy)

        height        = mi.getHeight()
        width         = mi.getWidth()
        for j in range(height):
            for i in range(width):
                self.assertAlmostEqual(mi.getImage().get(i,j), 10 / (1./scaling), 5)

    def testIllum1(self):
        self.doIllum(scaling=10)

    def testIllum2(self):
        self.doIllum(scaling=0.1)

    def testIllum3(self):
        self.doIllum(scaling=3.7)
    
#####
        
def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(IsrTestCases)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
