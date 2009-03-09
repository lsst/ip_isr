#!/usr/bin/env python
import os

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr
import lsst.pex.logging as logging

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

    def testLinearizationReplace(self):
        # create a basic lookup table
        lookupPolicy = pexPolicy.Policy()
        lookupPolicy.set('type', 'Replace')
        lookupPolicy.set('length', 10)
        lookupPolicy.add('value', 0)
        lookupPolicy.add('value', 1)
        lookupPolicy.add('value', 2)
        lookupPolicy.add('value', 3)
        lookupPolicy.add('value', 4)
        lookupPolicy.add('value', 5)
        lookupPolicy.add('value', 7)
        lookupPolicy.add('value', 7)
        lookupPolicy.add('value', 7)
        lookupPolicy.add('value', 9)
        mi       = afwImage.MaskedImageI(10,10)
        exposure = afwImage.ExposureI(mi)
        for i in range(10):
            for j in range(10):
                exposure.getMaskedImage().getImage().set(i, j, i)
                exposure.getMaskedImage().getVariance().set(i, j, i)

        lookupTable = ipIsr.LookupTableFromPolicy(lookupPolicy)
        ipIsr.Linearization(exposure, self.policy, lookupTable=lookupTable)

        mi = exposure.getMaskedImage()
        for i in range(10):
            for j in range(10):
                if (i != 6) and (i != 8):
                    self.assertEqual(mi.getImage().get(i,j), i)
                    self.assertEqual(mi.getVariance().get(i,j), i)
                else:
                    self.assertEqual(mi.getImage().get(i,j), 7)
                    self.assertEqual(mi.getVariance().get(i,j), 7)

    def testLinearizationMultiplicative(self):
        # create a basic lookup table
        lookupPolicy = pexPolicy.Policy()
        lookupPolicy.set('type', 'Multiplicative')
        lookupPolicy.set('length', 10)
        lookupPolicy.add('value', 1.)
        lookupPolicy.add('value', 1.)
        lookupPolicy.add('value', 1.)
        lookupPolicy.add('value', 1.)
        lookupPolicy.add('value', 1.)
        lookupPolicy.add('value', 1.)
        lookupPolicy.add('value', 1.)
        lookupPolicy.add('value', 1.1)
        lookupPolicy.add('value', 1.2)
        lookupPolicy.add('value', 1.)
        mi       = afwImage.MaskedImageF(10,10)
        exposure = afwImage.ExposureF(mi)
        for i in range(10):
            for j in range(10):
                exposure.getMaskedImage().getImage().set(i, j, i)
                exposure.getMaskedImage().getVariance().set(i, j, i)

        lookupTable = ipIsr.LookupTableFromPolicy(lookupPolicy)
        ipIsr.Linearization(exposure, self.policy, lookupTable=lookupTable)

        mi = exposure.getMaskedImage()
        for i in range(10):
            for j in range(10):
                if (i == 7):
                    self.assertAlmostEqual(mi.getImage().get(i,j),    i * 1.1,    5)
                    self.assertAlmostEqual(mi.getVariance().get(i,j), i * 1.1**2, 5)
                elif (i == 8):
                    self.assertAlmostEqual(mi.getImage().get(i,j),    i * 1.2,    5)
                    self.assertAlmostEqual(mi.getVariance().get(i,j), i * 1.2**2, 5)
                else:
                    self.assertAlmostEqual(mi.getImage().get(i,j),    i)
                    self.assertAlmostEqual(mi.getVariance().get(i,j), i)

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
