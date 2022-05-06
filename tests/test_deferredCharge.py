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
import numpy as np
import tempfile
import lsst.utils.tests

from lsst.ip.isr import DeferredChargeCalib, DeferredChargeTask, SerialTrap


class SerialTrapTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.nx = 64
        self.ny = 64
        self.prescan = 8
        self.trapSize = 10.0
        self.trapDecay = 3.0
        self.trapPixel = self.nx + self.prescan - 3
        self.trapCoeffs = [5.0, 1e-3]

    def testTrapsLinear(self):
        trap = SerialTrap(self.trapSize, self.trapDecay,
                          self.trapPixel, 'linear', self.trapCoeffs)
        trap.initialize(self.ny, self.nx, self.prescan)

        signals = np.full((self.ny, self.nx + self.prescan), 100.0)
        signals[:, 0:self.prescan] = 0.0

        freeCharge = np.full((self.ny, self.nx + self.prescan), 100.0)
        freeCharge[:, 0:self.prescan] = 0.0

        # Captured charge is of size trapSize, stored in column
        # trapPixel, and introduced to all rows.
        capturedCharge = trap.trap_charge(signals)
        self.assertEqual(capturedCharge.shape, (self.ny, self.nx + self.prescan))
        self.assertEqual(np.sum(capturedCharge), 10.0 * (self.ny))

        # Released charge is the amount released after one decay.
        releasedCharge = trap.release_charge()
        self.assertEqual(releasedCharge.shape, (self.ny, self.nx + self.prescan))
        self.assertAlmostEqual(np.sum(releasedCharge), 181.41996, 4)

    def testTrapsLogistic(self):
        trap = SerialTrap(self.trapSize, self.trapDecay,
                          self.trapPixel, 'logistic', self.trapCoeffs)
        trap.initialize(self.ny, self.nx, self.prescan)

        signals = np.full((self.ny, self.nx + self.prescan), 100.0)
        freeCharge = np.full((self.ny, self.nx + self.prescan), 100.0)
        signals[:, 0:self.prescan] = 0.0
        freeCharge[:, 0:self.prescan] = 0.0

        # Captured charge is of size trapSize, stored in column
        # trapPixel, and introduced to all rows.
        capturedCharge = trap.trap_charge(signals)
        self.assertEqual(capturedCharge.shape, (self.ny, self.nx + self.prescan))
        self.assertAlmostEqual(np.sum(capturedCharge), 5.237321 * (self.ny), 4)

        # Released charge is the amount released after one decay.
        releasedCharge = trap.release_charge()
        self.assertEqual(releasedCharge.shape, (self.ny, self.nx + self.prescan))
        self.assertAlmostEqual(np.sum(releasedCharge), 95.015467, 4)


class DeferredChargeTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.nx = 64
        self.ny = 64
        self.prescan = 8
        self.overscan = 16
        self.trapSize = 10.0
        self.trapDecay = 3.0
        self.trapPixel = self.nx//2 + self.prescan
        self.trapCoeffs = [5.0, 1e-3]
        self.calib = DeferredChargeCalib()

        self.calib.driftScale['amp0'] = 1.8e-4
        self.calib.driftScale['amp1'] = 2.8e-4
        self.calib.decayTime['amp0'] = 3.1
        self.calib.decayTime['amp1'] = 3.2
        self.calib.globalCti['amp0'] = 1.4e-7
        self.calib.globalCti['amp1'] = 1.5e-7

        self.calib.serialTraps['amp0'] = SerialTrap(self.trapSize, self.trapDecay,
                                                    self.trapPixel, 'linear', self.trapCoeffs)
        self.calib.serialTraps['amp0'] = SerialTrap(self.trapSize, self.trapDecay,
                                                    self.trapPixel, 'logistic', self.trapCoeffs)

    def testFullyPopulated(self):
        # Do IO tests.
        outPath = tempfile.mktemp() + '.yaml'
        self.calib.writeText(outPath)
        newCalib = DeferredChargeCalib().readText(outPath)
        self.assertEqual(self.calib, newCalib)

        outPath = tempfile.mktemp() + '.fits'
        self.calib.writeFits(outPath)
        newCalib = DeferredChargeCalib().readFits(outPath)
        self.assertEqual(self.calib, newCalib)

    def testTaskMethods(self):
        task = DeferredChargeTask()

        image = np.full((self.ny, self.nx + self.prescan + self.overscan), 100.0)
        image[:, 0:self.prescan] = 0.0
        image[:, -self.overscan:] = 1.0

        corrected = task.local_offset_inverse(image, self.calib.driftScale['amp0'],
                                              self.calib.decayTime['amp0'],
                                              num_previous_pixels=15)
        # 64*64*100 + 16*64 * ~(1 - driftScale) ~= 77
        self.assertAlmostEqual(np.sum(corrected), 410547.255925, 5)

        corrected = task.local_trap_inverse(corrected, self.calib.serialTraps['amp0'],
                                            self.calib.globalCti['amp0'],
                                            num_previous_pixels=6)
        # As above + ~320 deposited in prescan
        self.assertAlmostEqual(np.sum(corrected), 410866.615852, 5)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
