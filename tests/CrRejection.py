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
#InputExposure  = os.path.join(isrDataDir, 'CFHT/D4', 'raw-53535-i-797722_1')

# Master Calibration Image Names
#InputBias      = os.path.join(isrDataDir, 'CFHT/D4', '05Am05.bias.0.36.00_1')    
#InputFlat      = os.path.join(isrDataDir, 'CFHT/D4', '05Am05.flat.i.36.01_1')
#InputFringe    = os.path.join(isrDataDir, 'CFHT/D4', '05Am05.fringe.i.36.00_1')

# Policy file
InputIsrPolicy = os.path.join(isrDir, 'pipeline', 'isrPolicy.paf')

class IsrTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = pexPolicy.Policy.createPolicy(InputIsrPolicy)
        
    def tearDown(self):
        del self.policy

    def testCrRejectionVariance(self):

        # RHL debiases the interpolation; is off by a small factor
        # when all the pixel values are equal and the bias is non-zero
        mi       = afwImage.MaskedImageF(20,20)
        mi.set(100, 0x0, 1)
        exposure = afwImage.ExposureF(mi, afwImage.Wcs())
        metadata = exposure.getMetadata()
        metadata.set('gain', 1.0)
        
        mi.set(7, 7, (1000, 0x0, 1))

        ipIsr.CrRejection(exposure, self.policy)

        bitmaskCr     = mi.getMask().getPlaneBitMask('CR')
        bitmaskInterp = mi.getMask().getPlaneBitMask('INTRP')
        height        = mi.getHeight()
        width         = mi.getWidth()

        for j in range(height):
            for i in range(width):
                if i == 7 and j == 7:
                    self.assertEqual(mi.getMask().get(i,j) & bitmaskInterp, bitmaskInterp)
                    self.assertEqual(mi.getMask().get(i,j) & bitmaskCr, bitmaskCr)
                    self.assertAlmostEqual(mi.getImage().get(i,j)/100., 1, 1)
                else:
                    self.assertEqual(mi.getMask().get(i,j), 0)
                    self.assertEqual(mi.getImage().get(i,j), 100)

    def testCrRejectionNoVariance(self):
        mi       = afwImage.MaskedImageF(20,20)
        mi.set(100, 0x0, 0)
        exposure = afwImage.ExposureF(mi, afwImage.Wcs())
        metadata = exposure.getMetadata()
        metadata.set('gain', 1.0)

        mi.set(7, 7, (1000, 0x0, 0))

        ipIsr.CrRejection(exposure, self.policy)

        bitmaskCr     = mi.getMask().getPlaneBitMask('CR')
        bitmaskInterp = mi.getMask().getPlaneBitMask('INTRP')
        height        = mi.getHeight()
        width         = mi.getWidth()

        for j in range(height):
            for i in range(width):
                if i == 7 and j == 7:
                    self.assertEqual(mi.getMask().get(i,j) & bitmaskInterp, bitmaskInterp)
                    self.assertEqual(mi.getMask().get(i,j) & bitmaskCr, bitmaskCr)
                    self.assertEqual(mi.getImage().get(i,j), 100)
                else:
                    self.assertEqual(mi.getMask().get(i,j), 0)
                    self.assertEqual(mi.getImage().get(i,j), 100)
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
