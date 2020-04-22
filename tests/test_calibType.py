#
# LSST Data Management System
# Copyright 2008-2020 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import unittest
import tempfile

import lsst.utils.tests

from lsst.ip.isr import IsrProvenance


class IsrCalibCases(lsst.utils.tests.TestCase):
    """Test unified calibration type.
    """
    def setUp(self):
        self.calib = IsrProvenance(detectorName='test_calibType Det00',
                                   detectorSerial='Det00',
                                   calibType="Test Calib")
        self.calib.updateMetadata()
        self.calib.fromDataIds([{'exposure': 1234, 'detector': 0, 'filter': 'G'},
                                {'exposure': 1235, 'detector': 0, 'filter': 'G'},
                                {'exposure': 1234, 'detector': 1, 'filter': 'G'},
                                {'exposure': 1235, 'detector': 1, 'filter': 'G'}])

    def runText(self, textType):
        filename = tempfile.mktemp()
        self.calib.writeText(filename + textType)
        fromText = IsrProvenance.readText(filename + textType)
        self.assertEqual(self.calib, fromText)

    def test_Text(self):
        self.runText('.yaml')
        self.runText('.ecsv')

    def test_Fits(self):
        filename = tempfile.mktemp()
        self.calib.writeFits(filename + '.fits')
        fromFits = IsrProvenance.readFits(filename + '.fits')
        self.assertEqual(self.calib, fromFits)

        fromFits.updateMetadata(setDate=True)
        self.assertNotEqual(self.calib, fromFits)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
