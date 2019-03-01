#
# LSST Data Management System
# Copyright 2008-2017 AURA/LSST.
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

import lsst.afw.geom as afwGeom
import lsst.utils.tests
import lsst.ip.isr.straylight as straylight
import lsst.ip.isr.vignette as vignette
import lsst.ip.isr.masking as masking
import lsst.ip.isr.linearize as linearize

import lsst.ip.isr.isrMock as isrMock


class IsrMiscCases(lsst.utils.tests.TestCase):

    def setUp(self):
        self.inputExp = isrMock.TrimmedRawMock().run()
        self.mi = self.inputExp.getMaskedImage()

    def test_straylight(self):
        """Assert that the straylight task does not error when given an exposure.
        """
        task = straylight.StrayLightTask()
        task.run(self.inputExp)
        self.assertTrue(True)

    def test_vignette_noWrite(self):
        """Assert that the vignette task does not error when given an exposure
        """
        config = vignette.VignetteConfig()
        config.radius = 125.0
        config.xCenter = 100.0
        config.yCenter = 100.0

        config.doWriteVignettePolygon = False
        task = vignette.VignetteTask(config=config)
        result = task.run(self.inputExp)

        self.assertIsNone(result)

    def test_vignette_doWrite(self):
        """Assert that the vignette task does not error when given an exposure
        """
        config = vignette.VignetteConfig()
        config.radius = 125.0
        config.xCenter = 100.0
        config.yCenter = 100.0

        config.doWriteVignettePolygon = True
        task = vignette.VignetteTask(config=config)
        result = task.run(self.inputExp)

        self.assertIsInstance(result, afwGeom.Polygon)

    def test_masking(self):
        """Assert that the masking task does not error when given an exposure.
        """
        task = masking.MaskingTask()
        result = task.run(self.inputExp)

        self.assertIsNone(result)

    def test_linearize(self):
        """Assert that the linearize task does not error when a linearity is requested.
        """
        for linearityTypeName in ('LookupTable', 'Squared', 'Unknown'):
            result = linearize.getLinearityTypeByName(linearityTypeName)
            # These return the actual class to use, so use Equal instead of IsInstance.
            if linearityTypeName == 'LookupTable':
                self.assertEqual(result, linearize.LinearizeLookupTable, msg=f"{linearityTypeName}")
            elif linearityTypeName == 'Squared':
                self.assertEqual(result, linearize.LinearizeSquared, msg=f"{linearityTypeName}")
            else:
                self.assertIsNone(result)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
