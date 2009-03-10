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

    def testSaturation(self):
        saturation = 1000
        
        saturationKeyword = self.policy.getString('saturationPolicy.saturationKeyword')
        growSaturated     = self.policy.getInt('saturationPolicy.growSaturated')

        mi       = afwImage.MaskedImageF(20,20)
        mi.set(100, 0x0, 1)
        exposure = afwImage.ExposureF(mi)
        metadata = exposure.getMetadata()
        metadata.setDouble(saturationKeyword, saturation)
        
        bbox     = afwImage.BBox(afwImage.PointI(9,5),
                                 afwImage.PointI(9,15))
        submi    = afwImage.MaskedImageF(mi, bbox)
        submi.set(saturation)

        ipIsr.SaturationCorrection(exposure, self.policy)

        bitmaskBad    = mi.getMask().getPlaneBitMask('BAD')
        bitmaskSat    = mi.getMask().getPlaneBitMask('SAT')
        bitmaskInterp = mi.getMask().getPlaneBitMask('INTRP')
        height        = mi.getHeight()
        width         = mi.getWidth()

        for j in range(height):
            for i in range(width):
                # Grown saturation mask; one around the mask at 9
                if i >= 8 and i <= 10:
                    if (j == 4 or j == 16) and (i == 8 or i == 10):
                        # Not saturated but interpolated over
                        self.assertEqual(mi.getMask().get(i,j) & bitmaskInterp, bitmaskInterp)
                    elif (j == 4 or j == 16):
                        # Both saturated and interpolated over; bottom/top
                        self.assertEqual(mi.getMask().get(i,j) & bitmaskInterp, bitmaskInterp)
                        self.assertEqual(mi.getMask().get(i,j) & bitmaskSat,    bitmaskSat)
                    elif (j > 4 and j < 16):
                        # Both saturated and interpolated over; guts of it
                        self.assertEqual(mi.getMask().get(i,j) & bitmaskInterp, bitmaskInterp)
                        self.assertEqual(mi.getMask().get(i,j) & bitmaskSat,    bitmaskSat)
                    else:
                        # Neither; above or below the mask
                        self.assertEqual(mi.getMask().get(i,j), 0)
                else:
                    self.assertEqual(mi.getMask().get(i,j), 0)
                           

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
