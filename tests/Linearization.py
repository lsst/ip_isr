#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.pex.policy as pexPolicy
import lsst.ip.isr as ipIsr
import lsst.pex.logging as logging
import lsst.daf.base as dafBase

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.isr', Verbosity)

isrDir     = eups.productDir('ip_isr')

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
        mi       = afwImage.MaskedImageF(10,10)
        exposure = afwImage.ExposureF(mi, afwImage.Wcs())

        # needed for lookup table application
        metadata = exposure.getMetadata()
        metadata.set('gain', 1.0)
        
        for i in range(10):
            for j in range(10):
                exposure.getMaskedImage().getImage().set(i, j, i)
                exposure.getMaskedImage().getVariance().set(i, j, i)

        lookupTable = ipIsr.lookupTableFromPolicy(lookupPolicy)
        ipIsr.linearization(exposure, lookupTable=lookupTable)

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
        exposure = afwImage.ExposureF(mi, afwImage.Wcs())

        # needed for lookup table application
        metadata = exposure.getMetadata()
        metadata.set('gain', 1.0)
        
        for i in range(10):
            for j in range(10):
                exposure.getMaskedImage().getImage().set(i, j, i)
                exposure.getMaskedImage().getVariance().set(i, j, i)

        lookupTable = ipIsr.lookupTableFromPolicy(lookupPolicy)
        ipIsr.linearization(exposure, lookupTable=lookupTable)

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
