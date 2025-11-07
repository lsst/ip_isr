# This file is part of meas_algorithms.
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

import os
import unittest
import numpy
import copy

import lsst.geom
import lsst.afw.image as afwImage
import lsst.meas.algorithms as algorithms
import lsst.utils.tests
from lsst.daf.base import PropertyList
from lsst.ip.isr import Defects

try:
    type(display)
except NameError:
    display = False
else:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)

# Determine if we have afwdata
try:
    afwdataDir = lsst.utils.getPackageDir('afwdata')
except Exception:
    afwdataDir = None

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class DefectsTestCase(lsst.utils.tests.TestCase):
    """Tests for collections of Defect."""

    def assertMetadata(self, first, second):
        """Compare the metadata associated with Defects"""

        # Must strip out DATE metadata before comparison
        meta1 = first.getMetadata()
        meta2 = second.getMetadata()
        for d in (meta1, meta2):
            for k in ("DATE", "CALIB_CREATION_DATE", "CALIB_CREATION_TIME"):
                if k in d:
                    del d[k]

        self.assertEqual(meta1, meta2)
        meta1["NEW"] = "additional header"
        self.assertNotEqual(first.getMetadata(), second.getMetadata())
        del meta1["NEW"]

    def test_defectsReason(self):
        defects = Defects()

        defects.append( lsst.geom.Box2I(lsst.geom.Point2I(5, 6),
                                        lsst.geom.Point2I(41, 50)), reason='HOT_PIXEL')

        defects.append( lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                        lsst.geom.Point2I(4, 5)), reason='EDGE')

        self.assertEqual(len(defects), 2)

        for d in defects:
            self.assertIsInstance(d, algorithms.Defect)

        meta = PropertyList()
        meta["TESTHDR"] = "testing"
        defects.setMetadata(meta)

        # Test the defects are written and read via FITS correctly
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            defects.writeFits(tmpFile)
            defects2 = Defects.readFits(tmpFile)

        self.assertTrue(defects.__eq__(defects2))

        # Tests masking with reason is working properly
        ccdImage = afwImage.MaskedImageF(250, 225)
        # test 1. test mask according to a given reason
        reason = 'HOT_PIXEL'
        defects.maskPixelsReason(ccdImage.mask, reason)
        self.assertEqual(numpy.sum(ccdImage.mask.array), 3330)
        ccdImageMaskBefore = copy.copy(ccdImage.mask.array)
        # test that masking with a reason not in the defects will result in
        # the same mask
        reason = 'VAMPIRE_PIXEL'
        defects.maskPixelsReason(ccdImage.mask, reason)
        self.assertEqual(numpy.sum(ccdImage.mask.array), numpy.sum(ccdImageMaskBefore))
        # test that masking with an inexistent reason raises
        reason = 'TEST_PIXEL'
        with self.assertRaises(RuntimeError):
            defects.maskPixelsReason(ccdImage.mask, reason)

        # test 2. test mask plane with reasons is set properly
        defects.setMaskPlaneReason(ccdImage.mask)
        self.assertEqual(numpy.sum(ccdImage.mask.array), 5250)


        # add check on new methods if we pass _defects w/o defectsUnnormalized, it raises properly.


    def test_defects(self):
        defects = Defects()

        defects.append(algorithms.Defect(lsst.geom.Box2I(lsst.geom.Point2I(5, 6),
                                         lsst.geom.Point2I(41, 50))))

        defects.append(lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                       lsst.geom.Point2I(4, 5)))
        defects.append(lsst.geom.Point2I(50, 50))
        defects.append(afwImage.DefectBase(lsst.geom.Box2I(lsst.geom.Point2I(100, 200),
                                                           lsst.geom.Extent2I(5, 5))))
        self.assertEqual(len(defects), 4)

        for d in defects:
            self.assertIsInstance(d, algorithms.Defect)

        # Transposition
        transposed = defects.transpose()
        self.assertEqual(len(transposed), len(defects))

        # Check that an individual defect is found properly transposed within
        # the outputs.
        found = False
        for defect in transposed:
            if defect.getBBox() == lsst.geom.Box2I(lsst.geom.Point2I(6, 5), lsst.geom.Extent2I(45, 37)):
                found = True
                break
        self.assertTrue(found)

        # Serialization round trip
        meta = PropertyList()
        meta["TESTHDR"] = "testing"
        defects.setMetadata(meta)

        table = defects.toFitsRegionTable()

        defects2 = Defects.fromTable([table])

        self.assertEqual(defects2, defects)

        # via FITS
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            defects.writeFits(tmpFile)
            defects2 = Defects.readFits(tmpFile)

        # Equality tests the bounding boxes so metadata is tested separately.
        self.assertEqual(defects2, defects)
        self.assertMetadata(defects2, defects)

        # via text file
        with lsst.utils.tests.getTempFilePath(".ecsv") as tmpFile:
            defects.writeText(tmpFile)
            defects2 = Defects.readText(tmpFile)

        # Equality tests the bounding boxes so metadata is tested separately.
        self.assertEqual(defects2, defects)
        self.assertMetadata(defects2, defects)

        # Check bad values
        with self.assertRaises(ValueError):
            defects.append(lsst.geom.Box2D(lsst.geom.Point2D(0., 0.),
                                           lsst.geom.Point2D(3.1, 3.1)))
        with self.assertRaises(ValueError):
            defects.append("defect")

    def testAstropyRegion(self):
        """Read a FITS region file created by Astropy regions."""
        # The file contains three regions:
        #
        # - Point2I(340, 344)
        # - Point2I(340, 344)
        # - Box2I(minimum=Point2I(5, -5), dimensions=Extent2I(10, 20))
        #
        # The two coincident points are combined on read, so we end up with two
        # defects.

        with self.assertLogs():
            defects = Defects.readFits(os.path.join(TESTDIR, "data", "fits_region.fits"),
                                       normalize_on_init=True)

        self.assertEqual(len(defects), 2)

    def testLsstTextfile(self):
        """Read legacy LSST text file format"""
        with lsst.utils.tests.getTempFilePath(".txt") as tmpFile:
            with open(tmpFile, "w") as fh:
                print("""# X0  Y0  width height
     996        0       56       24
       0     4156     2048       20
       0        0       17     4176
    1998     4035       50      141
    1023        0        2     4176
    2027        0       21     4176
       0     4047       37      129
# Some rows without fixed column widths
14 20 2000 50
10 10 10 10
""", file=fh)

            defects = Defects.readLsstDefectsFile(tmpFile, normalize_on_init=True)

        # Although there are 9 defects listed above, we record 11 after
        # normalization. This is due to non-optimal behaviour in
        # Defects.fromMask; see DM-24781.
        self.assertEqual(len(defects), 11)

    def test_normalize_defects(self):
        """A test for the lsst.meas.algorithms.Defect.normalize() method.
        """
        defects = Defects()

        # First series of 1-pixel contiguous defects
        for yPix in range(1, 6):
            defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(15, yPix),
                           dimensions=lsst.geom.Extent2I(1, 1)))

        # Defects are normalized as they are added; check that the above have
        # been merged into a single bounding box.
        self.assertEqual(len(defects), 1)

        # Second series of 1-pixel contiguous defects in bulk mode
        with defects.bulk_update():
            for yPix in range(11, 16):
                defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(20, yPix),
                                               dimensions=lsst.geom.Extent2I(1, 1)))
            # In bulk mode, defects are not normalized.
            self.assertEqual(len(defects), 6)

        # Normalization applied on exiting bulk mode.
        self.assertEqual(len(defects), 2)

        boxesMeasured = []
        for defect in defects:
            boxesMeasured.append(defect.getBBox())

        # The normalizing function should have created the following two boxes
        # out of the individual 1-pixel defects from above
        expectedDefects = [lsst.geom.Box2I(corner=lsst.geom.Point2I(15, 1),
                                           dimensions=lsst.geom.Extent2I(1, 5)),
                           lsst.geom.Box2I(corner=lsst.geom.Point2I(20, 11),
                                           dimensions=lsst.geom.Extent2I(1, 5))]

        self.assertEqual(len(expectedDefects), len(boxesMeasured))
        for expDef, measDef in zip(expectedDefects, boxesMeasured):
            self.assertEqual(expDef, measDef)

        # Normalize two distinct sets of Defects and ensure they compare to the
        # same thing.
        defects = Defects()
        # Set 1
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 1), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 2), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 3), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 4), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 5), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 6), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 7), dimensions=lsst.geom.Extent2I(1, 1)))
        defects.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 8), dimensions=lsst.geom.Extent2I(1, 1)))

        # Set 2
        defects2 = Defects()
        defects2.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 1), dimensions=lsst.geom.Extent2I(1, 5)))
        defects2.append(lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 5), dimensions=lsst.geom.Extent2I(1, 4)))

        self.assertEqual(defects, defects2)

        boxesMeasured, boxesMeasured2 = [], []
        for defect, defect2 in zip(defects, defects2):
            boxesMeasured.append(defect.getBBox())
            boxesMeasured2.append(defect2.getBBox())

        expectedDefects = [lsst.geom.Box2I(corner=lsst.geom.Point2I(25, 1),
                                           dimensions=lsst.geom.Extent2I(1, 8))]

        self.assertEqual(len(expectedDefects), len(boxesMeasured))
        for expDef, measDef in zip(expectedDefects, boxesMeasured):
            self.assertEqual(expDef, measDef)

        self.assertEqual(len(expectedDefects), len(boxesMeasured2))
        for expDef, measDef in zip(expectedDefects, boxesMeasured2):
            self.assertEqual(expDef, measDef)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
