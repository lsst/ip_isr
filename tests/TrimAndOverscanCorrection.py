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

        # we need to grab a valid wcs because of ticket 686
        # http://dev.lsstcorp.org/trac/ticket/686
        self.wcs    = afwImage.ExposureF(InputExposure).getWcs()
        
    def tearDown(self):
        del self.policy

    def testOverscanCorrectionY(self):
        mi = afwImage.MaskedImageF(10,13)
        mi.getImage().set(10)

        # these should be functionally equivalent
        bbox     = afwImage.BBox(afwImage.PointI(0,10),
                                 afwImage.PointI(9,12))
        biassec  = '[1:10,11:13]'
        overscan = afwImage.MaskedImageF(mi, bbox)
        overscan.set(2)
        
        overscanKeyword = self.policy.getString('overscanPolicy.overscanKeyword')
        exposure = afwImage.ExposureF(mi, self.wcs)
        metadata = exposure.getMetadata()
        metadata.setString(overscanKeyword, biassec)

        ipIsr.OverscanCorrection(exposure, self.policy)

        height        = mi.getHeight()
        width         = mi.getWidth()
        for j in range(height):
            for i in range(width):
                if j >= 10:
                    self.assertEqual(mi.getImage().get(i,j), 0)
                else:
                    self.assertEqual(mi.getImage().get(i,j), 8)

    def testOverscanCorrectionX(self):
        mi = afwImage.MaskedImageF(13,10)
        mi.getImage().set(10)

        # these should be functionally equivalent
        bbox     = afwImage.BBox(afwImage.PointI(10,0),
                                 afwImage.PointI(12,9))
        biassec  = '[11:13,1:10]'
        overscan = afwImage.MaskedImageF(mi, bbox)
        overscan.set(2)
        
        overscanKeyword = self.policy.getString('overscanPolicy.overscanKeyword')
        exposure = afwImage.ExposureF(mi, self.wcs)
        metadata = exposure.getMetadata()
        metadata.setString(overscanKeyword, biassec)

        ipIsr.OverscanCorrection(exposure, self.policy)

        height        = mi.getHeight()
        width         = mi.getWidth()
        for j in range(height):
            for i in range(width):
                if i >= 10:
                    self.assertEqual(mi.getImage().get(i,j), 0)
                else:
                    self.assertEqual(mi.getImage().get(i,j), 8)

    def testTrimY(self):
        mi = afwImage.MaskedImageF(10,13)
        mi.getImage().set(10)

        # these should be functionally equivalent
        bbox     = afwImage.BBox(afwImage.PointI(0,10),
                                 afwImage.PointI(9,12))
        trimsec  = '[1:10,11:13]'
        trim     = afwImage.MaskedImageF(mi, bbox)
        trim.set(2)
        
        trimsecKeyword = self.policy.getString('trimPolicy.trimsecKeyword')
        exposure = afwImage.ExposureF(mi)
        metadata = exposure.getMetadata()
        metadata.setString(trimsecKeyword, trimsec)

        exposure2 = ipIsr.TrimNew(exposure, self.policy)
        mi2       = exposure2.getMaskedImage()

        height        = mi2.getHeight()
        width         = mi2.getWidth()
        self.assertEqual(height, 10)
        self.assertEqual(width,  10)
        for j in range(height):
            for i in range(width):
                self.assertEqual(mi2.getImage().get(i,j), 10)

    def testTrimX(self):
        mi = afwImage.MaskedImageF(13,10)
        mi.getImage().set(10)

        # these should be functionally equivalent
        bbox     = afwImage.BBox(afwImage.PointI(10,0),
                                 afwImage.PointI(12,9))
        trimsec  = '[11:13,1:10]'
        trim     = afwImage.MaskedImageF(mi, bbox)
        trim.set(2)
        
        trimsecKeyword = self.policy.getString('trimPolicy.trimsecKeyword')
        exposure = afwImage.ExposureF(mi)
        metadata = exposure.getMetadata()
        metadata.setString(trimsecKeyword, trimsec)

        exposure2 = ipIsr.TrimNew(exposure, self.policy)
        mi2       = exposure2.getMaskedImage()

        height        = mi2.getHeight()
        width         = mi2.getWidth()
        self.assertEqual(height, 10)
        self.assertEqual(width,  10)
        for j in range(height):
            for i in range(width):
                self.assertEqual(mi2.getImage().get(i,j), 10)
    
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
