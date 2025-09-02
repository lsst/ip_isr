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
from scipy.interpolate import Akima1DInterpolator

import numpy as np

import lsst.utils.tests

from lsst.ip.isr import FlatGradient


class FlatGradientTest(lsst.utils.tests.TestCase):
    """Test the FlatGradient dataset."""
    def setUp(self):

        radialNodes = np.array([0., 200., 250., 300., 310., 320., 330., 340., 350., 360., 368.])
        radialValues = np.array(
            [
                1.03285683,
                1.03457546,
                1.02751129,
                0.98331483,
                0.97569776,
                0.93527889,
                0.79704968,
                0.67758559,
                0.60344388,
                0.45119814,
                -0.22979413,
            ],
        )

        self.flatGradient = FlatGradient()
        self.flatGradient.setParameters(
            radialSplineNodes=radialNodes,
            radialSplineValues=radialValues,
            itlRatio=0.9,
            centroidX=0.0,
            centroidY=0.0,
            centroidDeltaX=0.5,
            centroidDeltaY=-0.5,
            gradientX=1e-4,
            gradientY=5e-5,
            normalizationFactor=1.1,
        )

    def _checkEqual(self, a, b):
        self.assertEqual(b.metadata, a.metadata)
        np.testing.assert_array_almost_equal(b.radialSplineNodes, a.radialSplineNodes)
        self.assertEqual(b.radialSplineNodes.dtype, a.radialSplineNodes.dtype)
        np.testing.assert_array_almost_equal(b.radialSplineValues, a.radialSplineValues)
        self.assertEqual(b.radialSplineValues.dtype, a.radialSplineValues.dtype)
        self.assertEqual(b.itlRatio, a.itlRatio)
        self.assertEqual(type(b.itlRatio), float)
        self.assertEqual(b.centroidX, a.centroidX)
        self.assertEqual(type(b.centroidX), float)
        self.assertEqual(b.centroidY, a.centroidY)
        self.assertEqual(type(b.centroidY), float)
        self.assertEqual(b.centroidDeltaX, a.centroidDeltaX)
        self.assertEqual(type(b.centroidDeltaX), float)
        self.assertEqual(b.centroidDeltaY, a.centroidDeltaY)
        self.assertEqual(type(b.centroidDeltaY), float)
        self.assertEqual(b.gradientX, a.gradientX)
        self.assertEqual(type(b.gradientX), float)
        self.assertEqual(b.gradientY, a.gradientY)
        self.assertEqual(type(b.gradientY), float)
        self.assertEqual(b.normalizationFactor, a.normalizationFactor)
        self.assertEqual(type(b.normalizationFactor), float)

    def testRoundTrip(self):
        """Test persistence round-tripping."""

        with tempfile.NamedTemporaryFile(suffix=".yaml") as f:
            usedFilename = self.flatGradient.writeText(f.name)
            fromText = FlatGradient.readText(usedFilename)
        self._checkEqual(fromText, self.flatGradient)

        with tempfile.NamedTemporaryFile(suffix=".fits") as f:
            usedFilename = self.flatGradient.writeFits(f.name)
            fromFits = FlatGradient.readFits(usedFilename)
        self._checkEqual(fromFits, self.flatGradient)

    def testComputeRadialSplineModelXY(self):
        x = np.arange(-300, 300, dtype=np.float64)
        y = np.arange(-300, 300, dtype=np.float64)

        spl = Akima1DInterpolator(
            self.flatGradient.radialSplineNodes,
            self.flatGradient.radialSplineValues,
        )
        centroidX = self.flatGradient.centroidX + self.flatGradient.centroidDeltaX
        centroidY = self.flatGradient.centroidY + self.flatGradient.centroidDeltaY
        radius = np.sqrt((x - centroidX)**2.
                         + (y - centroidY)**2.)
        radial = spl(np.clip(radius, 0.0, 368.0))

        model = self.flatGradient.computeRadialSplineModelXY(x, y)

        np.testing.assert_array_almost_equal(model, radial)

    def testComputeRadialSplineModel(self):
        radius = np.linspace(0, 368.0, 1000)

        spl = Akima1DInterpolator(
            self.flatGradient.radialSplineNodes,
            self.flatGradient.radialSplineValues,
        )
        radial = spl(np.clip(radius, 0.0, 368.0))

        model = self.flatGradient.computeRadialSplineModel(radius)

        np.testing.assert_array_almost_equal(model, radial)

    def testComputeGradientModel(self):
        x = np.arange(-300, 300, dtype=np.float64)
        y = np.arange(-300, 300, dtype=np.float64)

        gradient = (
            1
            + self.flatGradient.gradientX*(x - self.flatGradient.centroidX)
            + self.flatGradient.gradientY*(y - self.flatGradient.centroidY)
        )

        model = self.flatGradient.computeGradientModel(x, y)

        np.testing.assert_array_almost_equal(model, gradient)

    def test_computeFullModel(self):
        x = np.arange(-300, 300, dtype=np.float64)
        y = np.arange(-300, 300, dtype=np.float64)
        is_itl = np.zeros(len(x), dtype=np.bool_)
        is_itl[0: 50] = True

        # Compute the comparison model.
        spl = Akima1DInterpolator(
            self.flatGradient.radialSplineNodes,
            self.flatGradient.radialSplineValues,
        )
        centroidX = self.flatGradient.centroidX + self.flatGradient.centroidDeltaX
        centroidY = self.flatGradient.centroidY + self.flatGradient.centroidDeltaY
        radius = np.sqrt((x - centroidX)**2.
                         + (y - centroidY)**2.)
        radial = spl(np.clip(radius, 0.0, 368.0))

        gradient = (
            1
            + self.flatGradient.gradientX*(x - self.flatGradient.centroidX)
            + self.flatGradient.gradientY*(y - self.flatGradient.centroidY)
        )

        full_model = radial / gradient
        full_model[is_itl] *= self.flatGradient.itlRatio

        model = self.flatGradient.computeFullModel(x, y, is_itl)

        np.testing.assert_array_almost_equal(full_model, model)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
