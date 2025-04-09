# This file is part of ip_isr.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
import unittest
import tempfile
import logging

import lsst.utils.tests

from lsst.ip.isr import IsrProvenance, IsrCalib


class FakeExposure():
    def __init__(self, **kwargs):
        self.metadata = kwargs

    def getMetadata(self):
        return self.metadata


class IsrCalibCases(lsst.utils.tests.TestCase):
    """Test unified calibration type.
    """
    def setUp(self):
        self.calib = IsrProvenance(detectorName='testCalibType Det00',
                                   detectorId=0,
                                   instrument="TestInst",
                                   calibType="Test Calib")
        self.calib.updateMetadata()
        self.calib.fromDataIds([{'exposure': 1234, 'detector': 0, 'filter': 'G'},
                                {'exposure': 1235, 'detector': 0, 'filter': 'G'},
                                {'exposure': 1234, 'detector': 1, 'filter': 'G'},
                                {'exposure': 1235, 'detector': 1, 'filter': 'G'}])
        self.exposure1 = FakeExposure(SEQNAME='sequencer A', SEQFILE='my.seq', SEQCKSUM='abcdef')
        self.exposure2 = FakeExposure(SEQNAME='sequencer B', SEQFILE='my.seq')

    def runText(self, textType):
        filename = tempfile.mktemp()
        usedFilename = self.calib.writeText(filename + textType)
        fromText = IsrProvenance.readText(usedFilename)
        self.assertEqual(self.calib, fromText)

        # Test generic interface:
        fromGeneric = IsrCalib.readText(usedFilename)
        self.assertEqual(self.calib, fromGeneric)

    def test_Text(self):
        self.runText('.yaml')
        self.runText('.ecsv')

    def test_Fits(self):
        filename = tempfile.mktemp()
        usedFilename = self.calib.writeFits(filename + '.fits')
        fromFits = IsrProvenance.readFits(usedFilename)
        self.assertEqual(self.calib, fromFits)

        with self.assertLogs(level=logging.WARNING) as cm:
            fromFits.updateMetadataFromExposures([self.exposure1, self.exposure2, None])
        self.assertIn("Metadata mismatch!", cm.output[0])
        fromFits.updateMetadata(setDate=True)
        self.assertNotEqual(self.calib, fromFits)

        # Test generic interface:
        fromGeneric = IsrCalib.readFits(usedFilename)
        self.assertEqual(self.calib, fromGeneric)

        self.assertEqual(fromFits.metadata["SEQNAME"], "sequencer A")
        self.assertEqual(fromFits.metadata["SEQFILE"], "my.seq")
        self.assertEqual(fromFits.metadata["SEQCKSUM"], "abcdef")
        self.assertEqual(fromFits.metadata["DET_NAME"], "testCalibType Det00")
        self.assertEqual(fromFits.metadata["DETECTOR"], 0)
        self.assertEqual(fromFits.metadata["INSTRUME"], "TestInst")


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
