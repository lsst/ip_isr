#
# LSST Data Management System
# Copyright 2008-2017 LSST Corporation.
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

import unittest

import lsst.utils.tests
import lsst.geom
import lsst.afw.image as afwImage
import lsst.ip.isr as ipIsr


class IsrTestCases(unittest.TestCase):

    def setUp(self):
        self.pmin = lsst.geom.Point2I(1, 1)
        self.pmax = lsst.geom.Point2I(10, 10)
        self.meanCountsKeyword = "IMMODE"
        self.filenameKeyword = "filename"

    def tearDown(self):
        del self.pmin
        del self.pmax
        del self.meanCountsKeyword
        del self.filenameKeyword

    def testBias(self):
        maskedImage = afwImage.MaskedImageF(lsst.geom.Box2I(self.pmin, self.pmax))
        maskedImage.getImage().set(10)

        bias = afwImage.MaskedImageF(lsst.geom.Box2I(self.pmin, self.pmax))
        bias.getImage().set(1)
        biasexposure = afwImage.ExposureF(bias, None)
        bmetadata = biasexposure.getMetadata()
        bmetadata.setDouble(self.meanCountsKeyword, 1.)
        bmetadata.setString(self.filenameKeyword, 'Unittest Bias')

        ipIsr.biasCorrection(maskedImage, biasexposure.getMaskedImage())

        height = maskedImage.getHeight()
        width = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                self.assertEqual(maskedImage.image[i, j, afwImage.LOCAL], 9)

    def doDark(self, scaling):
        maskedImage = afwImage.MaskedImageF(lsst.geom.Box2I(self.pmin, self.pmax))
        maskedImage.getImage().set(10)

        dark = afwImage.MaskedImageF(lsst.geom.Box2I(self.pmin, self.pmax))
        dark.getImage().set(1)
        darkexposure = afwImage.ExposureF(dark, None)
        dmetadata = darkexposure.getMetadata()
        dmetadata.setString(self.filenameKeyword, 'Unittest Dark')

        ipIsr.darkCorrection(maskedImage, darkexposure.getMaskedImage(), 1., scaling)

        height = maskedImage.getHeight()
        width = maskedImage.getWidth()
        for j in range(height):
            for i in range(width):
                self.assertAlmostEqual(maskedImage.image[i, j, afwImage.LOCAL], 10 - 1./scaling, 5)

    def testDark1(self):
        self.doDark(scaling=10)

    def testDark2(self):
        self.doDark(scaling=0.1)

    def testDark3(self):
        self.doDark(scaling=3.7)

    def testDarkWithDarktime(self):
        darkTime = 128.0
        nan = float("NAN")

        exp = afwImage.ExposureF(1, 1)
        exp.getMaskedImage().getImage().set(1.0)
        exp.getMaskedImage().getMask().set(0)
        exp.getMaskedImage().getVariance().set(1.0)

        dark = afwImage.ExposureF(1, 1)
        dark.getMaskedImage().getImage().set(1.0/darkTime)
        dark.getMaskedImage().getMask().set(0)
        dark.getMaskedImage().getVariance().set(0.0)
        dark.getInfo().setVisitInfo(afwImage.VisitInfo())

        task = ipIsr.IsrTask()

        # No darktime set in at least one of the inputs
        exp.getInfo().setVisitInfo(afwImage.VisitInfo(darkTime=nan))
        dark.getInfo().setVisitInfo(afwImage.VisitInfo(darkTime=nan))
        with self.assertRaises(RuntimeError):
            task.darkCorrection(exp, dark)
        exp.getInfo().setVisitInfo(afwImage.VisitInfo(darkTime=nan))
        dark.getInfo().setVisitInfo(afwImage.VisitInfo(darkTime=1.0))
        with self.assertRaises(RuntimeError):
            task.darkCorrection(exp, dark)

        # With darktime set
        exp.getInfo().setVisitInfo(afwImage.VisitInfo(darkTime=darkTime))
        dark.getInfo().setVisitInfo(afwImage.VisitInfo(darkTime=1.0))
        task.darkCorrection(exp, dark)

        self.assertEqual(exp.image[0, 0, afwImage.LOCAL], 0.0)
        self.assertEqual(exp.mask[0, 0, afwImage.LOCAL], 0)
        self.assertEqual(exp.variance[0, 0, afwImage.LOCAL], 1.0)
        self.assertEqual(exp.getInfo().getVisitInfo().getDarkTime(), darkTime)  # Hasn't been modified

    def testDarkWithDarktimeNan(self):
        darkTime = 128.0
        nan = float("NAN")

        exp = afwImage.ExposureF(1, 1)
        exp.getMaskedImage().getImage().set(1.0)
        exp.getMaskedImage().getMask().set(0)
        exp.getMaskedImage().getVariance().set(1.0)

        dark = afwImage.ExposureF(1, 1)
        dark.getMaskedImage().getImage().set(1.0/darkTime)
        dark.getMaskedImage().getMask().set(0)
        dark.getMaskedImage().getVariance().set(0.0)
        dark.getInfo().setVisitInfo(afwImage.VisitInfo())

        task = ipIsr.IsrTask()

        # scale with darkScale=1 if the dark has darkTime=NaN.
        exp.getInfo().setVisitInfo(afwImage.VisitInfo(darkTime=darkTime))
        dark.getInfo().setVisitInfo(afwImage.VisitInfo(darkTime=nan))
        task.darkCorrection(exp, dark)

        self.assertEqual(exp.image[0, 0, afwImage.LOCAL], 0.0)
        self.assertEqual(exp.mask[0, 0, afwImage.LOCAL], 0)
        self.assertEqual(exp.variance[0, 0, afwImage.LOCAL], 1.0)
        self.assertEqual(exp.getInfo().getVisitInfo().getDarkTime(), darkTime)  # Hasn't been modified


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
