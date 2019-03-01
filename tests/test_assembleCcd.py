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

import lsst.utils.tests
from lsst.ip.isr.assembleCcdTask import (AssembleCcdConfig, AssembleCcdTask)
import lsst.ip.isr.isrMock as isrMock


class AssembleCcdTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        self.outputExp = isrMock.TrimmedRawMock().run()
        self.outputUntrimmedExp = isrMock.RawMock().run()

    def runTest(self, inputExp=None, doTrim=False):
        """Test guts to avoid code duplication.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            Exposure to assemble.
        doTrim : `bool`
            To trim the input or not to trim.

        Returns
        -------
        assembleOutput : `lsst.afw.image.Exposure`
            Assembled exposure.

        Raises
        ------
        TypeError
            Raised if the inputExp is None.
        """
        self.config = AssembleCcdConfig(doTrim=doTrim,
                                        keysToRemove=['SHEEP', 'MONKEYS', 'ZSHEEP'])
        self.task = AssembleCcdTask(config=self.config)

        return self.task.assembleCcd(inputExp)

    def testAssembleCcdTask_exp_noTrim(self):
        self.assertEqual(self.runTest(inputExp=isrMock.RawMock().run(), doTrim=False).getBBox(),
                         self.outputUntrimmedExp.getBBox())

    def testAssembleCcdTask_exp_doTrim(self):
        self.assertEqual(self.runTest(inputExp=isrMock.RawMock().run(), doTrim=True).getBBox(),
                         self.outputExp.getBBox())

    def testAssembleCcdTask_expDict_noTrim(self):
        self.assertEqual(self.runTest(inputExp=isrMock.RawDictMock().run(), doTrim=False).getBBox(),
                         self.outputUntrimmedExp.getBBox())

    def testAssembleCcdTask_expDict_doTrim(self):
        self.assertEqual(self.runTest(inputExp=isrMock.RawDictMock().run(), doTrim=True).getBBox(),
                         self.outputExp.getBBox())

    def testAssembleCcdTask_fail(self):
        """Assembly should fail if no exposure is supplied.
        """
        with self.assertRaises(TypeError):
            self.runTest(inputExp=None, doTrim=False)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
