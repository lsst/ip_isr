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
import os
import numpy as np
import tempfile

import lsst.geom
import lsst.utils.tests

from lsst.afw.image import TransmissionCurve
from lsst.ip.isr import IntermediateTransmissionCurve


TESTDIR = os.path.abspath(os.path.dirname(__file__))


class TransmissionCurveCases(lsst.utils.tests.TestCase):
    """Test intermediate transmission curve calibration type.
    """
    def setUp(self):
        rng = np.random.Generator(np.random.MT19937(1))
        self.points = [lsst.geom.Point2D(rng.random(), rng.random()) for i in range(5)]

        self.curve1 = IntermediateTransmissionCurve.readText(
            os.path.join(TESTDIR, "data", "test_curve1.ecsv"))
        self.curve2 = IntermediateTransmissionCurve.readText(
            os.path.join(TESTDIR, "data", "test_curve2.ecsv"))
        self.curve3 = IntermediateTransmissionCurve.readText(
            os.path.join(TESTDIR, "data", "test_curve3.ecsv"))

    def assertTransmissionCurvesEqual(self, a, b, rtol=1e-6, atol=0.0):
        """Test whether two TransimssionCurves are equivalent.
        From afw/tests/test_transmissionCurve.py
        """
        self.assertEqual(a.getWavelengthBounds(), b.getWavelengthBounds())
        self.assertEqual(a.getThroughputAtBounds(), b.getThroughputAtBounds())
        wavelengths = np.linspace(*(a.getWavelengthBounds() + (100,)))
        for point in self.points:
            self.assertFloatsAlmostEqual(
                a.sampleAt(point, wavelengths),
                b.sampleAt(point, wavelengths),
                rtol=rtol, atol=atol
            )

    def test_construction(self):
        self.assertTransmissionCurvesEqual(self.curve1.transmissionCurve,
                                           self.curve2.transmissionCurve)
        self.assertTransmissionCurvesEqual(self.curve1.transmissionCurve,
                                           self.curve3.transmissionCurve)

    def test_output(self):
        filename1 = tempfile.mktemp()
        self.curve1.writeFits(filename1)

        reread = TransmissionCurve.readFits(filename1)
        self.assertTransmissionCurvesEqual(reread, self.curve1.transmissionCurve)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
