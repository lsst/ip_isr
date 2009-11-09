#!/usr/bin/env python
"""
Run with:
   python IsrStageTest.py
"""

import sys, os, math

import eups
import pdb
import unittest

import lsst.utils.tests as utilsTests
import lsst.pex.harness.Queue as pexQueue
import lsst.pex.harness.Clipboard as pexClipboard
import lsst.pex.policy as pexPolicy
import lsst.pex.logging as pexLog
import lsst.ip.isr.IsrStages as isrStages
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.ip.isr.IsrStages as isrStages
import lsst.afw.display.ds9 as ds9

from lsst.ip.isr.MetadataStages import transformMetadata

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Verbosity = 4
pexLog.Trace_setVerbosity('lsst.ip.isr', Verbosity)

isrDataDir    = eups.productDir('isrdata')
if isrDataDir:
    inputImage    = os.path.join(isrDataDir, 'CFHT/D4/dc3a', 'raw-704893-e000-c000-a000.fits')
else:
    print >> sys.stderr, "Skipping tests as isrdata is not setup"

isrDir        = eups.productDir('ip_isr')

# For these tests, we need the policy files in ip_isr...
policyDir = isrDir+'/pipeline'
dc3MetadataPath  = os.path.join(policyDir, 'dc3MetadataPolicy.paf')
cfhtMetadataPath    = os.path.join(policyDir, 'cfhtDataTypePolicy.paf')
cfhtCalibrationPath = os.path.join(policyDir, 'cfhtCalibrationTypePolicy.paf')
ampPolicyPath = os.path.join(policyDir, 'cfhtAmpBBoxPolicy.paf')

dc3MetadataPolicy  = pexPolicy.Policy.createPolicy(dc3MetadataPath)
cfhtMetadataPolicy = pexPolicy.Policy.createPolicy(cfhtMetadataPath)
cfhtCalibrationPolicy = pexPolicy.Policy.createPolicy(cfhtCalibrationPath)
ampPolicy = pexPolicy.Policy.createPolicy(ampPolicyPath)

class IsrStageTestCase(unittest.TestCase):
    """A test case for IsrStages.py"""

    def setUp(self):
        # processing policy
        self.policy = pexPolicy.Policy()

        # PIPELINE INPUTS
        self.policy.add('inputImageKey',    'inputImage0')
        self.policy.add('inputMetadataKey', 'inputMetadata0')
        self.policy.add('calibDataKey',     'calibData') 
        # ISR PROCESSING
        self.policy.add('isrPolicy', pexPolicy.Policy.createPolicy(os.path.join(isrDir, 'pipeline', 'isrPolicy.paf')))
        # OUTPUTS
        self.policy.add('calibratedExposureKey', 'calibratedExposure0')
        self.policy.add('sdqaRatingSetKey',      'sdqaRatingSet0')
        
        
        # create clipboard and fill 'er up!
        clipboard = pexClipboard.Clipboard()

        # with : calibration information
        calibData = dafBase.PropertySet()
        biasPath  = os.path.join(isrDataDir, 'CFHT/D4/dc3a', 'bias-0-c000-a000')
        darkPath  = os.path.join(isrDataDir, 'CFHT/D4/dc3a', 'dark-300-c000-a000')
        flatPath  = os.path.join(isrDataDir, 'CFHT/D4/dc3a', 'flat-i-c000-a000')
        calibData.set('bias', biasPath)
        calibData.set('dark', darkPath)
        calibData.set('flat', flatPath)
        
        calibData.set('defectPath',    os.path.join(isrDataDir, 'CFHT/D4/dc3a', 'defect-c000-a000.paf'))
        calibData.set('linearizePath', os.path.join(isrDir, 'pipeline', 'linearizationLookupTable.paf'))
        clipboard.put('calibData', calibData)

        # with : an input image
        self.img = afwImage.ImageF(inputImage)
        clipboard.put('inputImage0', self.img)


        # with : input (transformed) metadata
        # 
        metadata = afwImage.readMetadata(inputImage)
        transformMetadata(metadata, cfhtMetadataPolicy, dc3MetadataPolicy, 'Keyword')
        clipboard.put('inputMetadata0', metadata)

        # with : calibration exposures
        bias         = afwImage.ExposureF(biasPath)
        biasMetadata = bias.getMetadata()
        transformMetadata(biasMetadata, cfhtCalibrationPolicy, dc3MetadataPolicy, 'Keyword')
        clipboard.put('biasExposure', bias)
        
        dark         = afwImage.ExposureF(darkPath)
        darkMetadata = dark.getMetadata()
        transformMetadata(darkMetadata, cfhtCalibrationPolicy, dc3MetadataPolicy, 'Keyword')
        clipboard.put('darkExposure', dark)

        flat         = afwImage.ExposureF(flatPath)
        flatMetadata = flat.getMetadata()
        transformMetadata(flatMetadata, cfhtCalibrationPolicy, dc3MetadataPolicy, 'Keyword')
        clipboard.put('flatExposure', flat)
        
        # with : an amp BBox
        ampstr  = "CcdBBox.Amp"+str(metadata.get("ampId"))
        ampheight  = ampPolicy.getInt(ampstr+".height")
        ampwidth   = ampPolicy.getInt(ampstr+".width")
        ampBBox = afwImage.BBox(afwImage.PointI(0,0),
                                ampwidth,
                                ampheight)
        clipboard.put('ampBBox', ampBBox)

        
        inQueue = pexQueue.Queue()
        inQueue.addDataset(clipboard)
        self.outQueue = pexQueue.Queue()
       
        self.stage = isrStages.IsrStage(1, self.policy)
        self.stage.initialize(self.outQueue, inQueue)
        self.stage.setUniverseSize(1)
        self.stage.setRun('SingleExposureTest')


    def tearDown(self):
        del self.stage
        del self.outQueue

    def testSingleInputExposure(self):
        self.stage.process()
        clipboard = self.outQueue.getNextDataset()
        assert(clipboard.contains(self.policy.getString('calibratedExposureKey')))
        assert(clipboard.contains(self.policy.getString('sdqaRatingSetKey')))

        if False:
            calibratedExposure = clipboard.get(self.policy.getString('calibratedExposureKey'))
            ds9.mtv(self.img, frame=1)
            ds9.mtv(calibratedExposure, frame=2)
        
def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    if isrDataDir:
        suites += unittest.makeSuite(IsrStageTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    if len(sys.argv) == 2:
        inputImage = sys.argv[1]
    run(True)
        
