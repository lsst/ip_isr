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
import logging
import numpy as np

import lsst.geom as geom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.utils.tests
import lsst.ip.isr as ipIsr
import lsst.ip.isr.isrFunctions as ipIsrFunctions
import lsst.ip.isr.isrMock as isrMock


def computeImageMedianAndStd(image):
    """Function to calculate median and std of image data.

    Parameters
    ----------
    image : `lsst.afw.image.Image`
        Image to measure statistics on.

    Returns
    -------
    median : `float`
        Image median.
    std : `float`
        Image stddev.
    """
    median = np.nanmedian(image.getArray())
    std = np.nanstd(image.getArray())

    return (median, std)


def _setExposureSatCore(exposure, x, y, halfWidthX, halfWidthY,
                        satVal, satMaskBit):
    """Helper function to make a rectangular saturated core.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Input exposure.
    x : `int`
        X coordinate for the center of the rectangular saturated core.
    y : `int`
        Y coordinate for the center of the rectangular saturated core.
    halfWidthX : `int`
        Half the width of the rectangular saturated core along the x axis.
    halfWidthY : `int`
        Half the width of the rectangular saturated core along the y axis.
    satVal : `float`
        Value the saturated core is set to.
    satMaskBit : `int`
        Saturated mask bit the saturated core mask is set to.
    """
    xmax = exposure.image.array.shape[1]
    lowerRange = x-halfWidthX
    upperRange = x+halfWidthX
    if lowerRange < 0:
        lowerRange = 0
    if upperRange > xmax:
        upperRange = xmax
    exposure.image.array[y-halfWidthY:y+halfWidthY,
                         lowerRange:upperRange] = satVal
    exposure.mask.array[y-halfWidthY:y+halfWidthY,
                        lowerRange:upperRange] = satMaskBit


def _setExposureSatColumns(exposure, x, y, halfWidthX, limY,
                           satVal, satMaskBit, isTop=False):
    """Helper function to add saturated columns.
    We need to add these around a saturated core to test a typical edge bleed.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Input exposure.
    x : `int`
        X coordinate for the central saturated column.
    y : `int`
        Y coordinate for the central saturated column.
    halfWidthX : `int`
        Half the number of saturated columns.
    limY : `int`
        Y coordinate up to which the saturated columns go to.
    satVal : `float`
        Value the saturated columns are set to.
    satMaskBit : `int`
        Saturated mask bit the saturated columns mask is set to.
    isTop : `bool`
        True if the saturated columns go toward the top of the detector.
    """
    if isTop:
        exposure.image.array[y:limY,
                             x-halfWidthX:x+halfWidthX] = satVal
        exposure.mask.array[y:limY,
                            x-halfWidthX:x+halfWidthX] = satMaskBit
    else:
        exposure.image.array[limY:y,
                             x-halfWidthX:x+halfWidthX] = satVal
        exposure.mask.array[limY:y,
                            x-halfWidthX:x+halfWidthX] = satMaskBit


def _makeEdgeBleed(exposure, x, extentY, edgeBleedWidth,
                   edgeBleedConstant, edgeBleedWidthLimit,
                   saturationLevel, saturationFrac,
                   satVal, satMaskBit, isTop=False):
    """Helper function to add an edge bleed with a decaying exponential
    model.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Input exposure.
    x : `int`
        X coordinate for the center of the edge bleed (generally the same
        as the center of the saturated core and columns).
    extentY : `int`
        Y coordinate up to which the edge bleed goes to  (generally the same
        as limY of the saturated columns).
    edgeBleedWidth : `int`
        Width of the edge bleed at the edge of the detector.
    edgeBleedConstant : `float`
        Constant for the decaying exponential model.
    edgeBleedWidthLimit : `int`
        The width the edge bleed goes to away from the edge (generally the
        width of the saturated columns).
    saturationLevel : `float`
        Saturation level.
    saturationFrac : `float`
        The inside of the edge bleed is set to a value equal to this fraction
        of the saturation level.
    satVal : `float`
        The value the contour of the edge bleed is set to, greater or equal to
        the saturation level.
    satMaskBit : `int`
        Saturated mask bit the contour of the edge bleed mask is set to.
    isTop : `bool`
        True if the edge bleed is a the top edge of the detector.
    """
    xmax = exposure.image.array.shape[1]
    for y in range(extentY):
        edgeBleedWidthY = edgeBleedWidth*np.exp(-edgeBleedConstant*y) \
            + edgeBleedWidthLimit
        if isTop:
            # For edge bleed in top amplifier
            y = exposure.image.array.shape[0]-1-y
        lowerRange = x-int(edgeBleedWidthY/2.)
        upperRange = x+int(edgeBleedWidthY/2.)
        if lowerRange < 0:
            lowerRange = 0
        if upperRange >= xmax:
            upperRange = xmax-1
        exposure.image.array[y, lowerRange] = satVal
        exposure.mask.array[y, lowerRange] = satMaskBit
        exposure.image.array[y, upperRange] = satVal
        exposure.mask.array[y, upperRange] = satMaskBit
        exposure.image.array[y,
                             lowerRange+1:upperRange
                             ] = saturationFrac*saturationLevel


class MockITLAmp:
    def __init__(self, name, bbox):
        self._name = name
        self._bbox = bbox

    def getName(self):
        return self._name

    def getBBox(self):
        return self._bbox

    def __repr__(self):
        return f"MockITLAmp({self._name})"


class MockITLDetector(list):
    def __init__(self):
        amps = []
        for i in range(8):
            name = f"C1{i}"
            bbox = geom.Box2I(corner=geom.Point2I(i*509, 2000), dimensions=geom.Extent2I(509, 2000))
            amps.append(MockITLAmp(name, bbox))
        for i in reversed(range(8)):
            name = f"C0{i}"
            bbox = geom.Box2I(corner=geom.Point2I(i*509, 0), dimensions=geom.Extent2I(509, 2000))
            amps.append(MockITLAmp(name, bbox))

        super().__init__(amps)

    def getBBox(self):
        return geom.Box2I(corner=geom.Point2I(0, 0), dimensions=geom.Extent2I(4072, 4000))


class MockITLExposure(afwImage.ExposureF):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    # Monkey-patch our fake detector.
    def setDetector(self, detector):
        self._detector = detector

    def fillSaturationMetadata(self, saturationLevel):
        if not self._detector:
            raise RuntimeError("No detector set.")
        for amp in self.getDetector():
            self.metadata[f"LSST ISR SATURATION LEVEL {amp.getName()}"] = saturationLevel

    def getDetector(self):
        return self._detector


class IsrFunctionsCases(lsst.utils.tests.TestCase):
    """Test that functions for ISR produce expected outputs.
    """
    def setUp(self):
        self.inputExp = isrMock.TrimmedRawMock().run()
        self.mi = self.inputExp.getMaskedImage()

    def test_transposeMaskedImage(self):
        """Expect height and width to be exchanged.
        """
        transposed = ipIsr.transposeMaskedImage(self.mi)
        self.assertEqual(transposed.getImage().getBBox().getHeight(),
                         self.mi.getImage().getBBox().getWidth())
        self.assertEqual(transposed.getImage().getBBox().getWidth(),
                         self.mi.getImage().getBBox().getHeight())

    def test_interpolateDefectList(self):
        """Expect number of interpolated pixels to be non-zero.
        """
        defectList = isrMock.DefectMock().run()
        self.assertEqual(len(defectList), 1)

        for fallbackValue in (None, -999.0):
            for haveMask in (True, False):
                with self.subTest(fallbackValue=fallbackValue, haveMask=haveMask):
                    if haveMask is False:
                        if 'INTRP' in self.mi.getMask().getMaskPlaneDict():
                            self.mi.getMask().removeAndClearMaskPlane('INTRP')
                    else:
                        if 'INTRP' not in self.mi.getMask().getMaskPlaneDict():
                            self.mi.getMask().addMaskPlane('INTRP')
                    numBit = ipIsr.countMaskedPixels(self.mi, "INTRP")
                    self.assertEqual(numBit, 0)

    def test_transposeDefectList(self):
        """Expect bbox dimension values to flip.
        """
        defectList = isrMock.DefectMock().run()
        transposed = defectList.transpose()

        for d, t in zip(defectList, transposed):
            self.assertEqual(d.getBBox().getDimensions().getX(), t.getBBox().getDimensions().getY())
            self.assertEqual(d.getBBox().getDimensions().getY(), t.getBBox().getDimensions().getX())

    def test_makeThresholdMask(self):
        """Expect list of defects to have elements.
        """
        defectList = ipIsr.makeThresholdMask(self.mi, 200,
                                             growFootprints=2,
                                             maskName='SAT')

        self.assertEqual(len(defectList), 1)

    def test_ITLEdgeBleedMask(self):
        """Expect number of masked pixels according to edge bleed masking.
        """
        detector = MockITLDetector()
        exposure = MockITLExposure(detector.getBBox())
        exposure.setDetector(detector)
        exposure.mask.array[:, :] = 0
        satMaskBit = exposure.mask.getPlaneBitMask('SAT')

        saturationLevel = 12000.
        exposure.fillSaturationMetadata(saturationLevel)

        badAmpDict = {}
        for amp in exposure.getDetector():
            if amp.getName() == 'C10':
                badAmpDict[amp.getName()] = True
            else:
                badAmpDict[amp.getName()] = False

        # Case 1:
        # Create two edge bleeds that fall into the same footprint.

        # A saturated core in the bottom amplifier
        # that has a bottom edge bleed
        satVal = 15000.
        x = 60
        y = 1800
        halfWidthX = 40
        halfWidthY = 70
        halfWidthColX = 3
        limY = 100
        extentY = 100
        _setExposureSatCore(exposure, x=x, y=y,
                            halfWidthX=halfWidthX, halfWidthY=halfWidthY,
                            satVal=satVal, satMaskBit=satMaskBit)
        _setExposureSatColumns(exposure, x=x, y=y, halfWidthX=halfWidthColX,
                               limY=limY, satVal=satVal, satMaskBit=satMaskBit,
                               isTop=False)
        edgeBleedWidth = 60
        edgeBleedConstant = 0.04  # we make the edge bleed a bit smaller than the masking model
        edgeBleedWidthLimit = 2*halfWidthColX
        saturationFrac = 0.85  # higher than default in masking model
        _makeEdgeBleed(exposure, x=x, extentY=extentY,
                       edgeBleedWidth=edgeBleedWidth,
                       edgeBleedConstant=edgeBleedConstant,
                       edgeBleedWidthLimit=edgeBleedWidthLimit,
                       saturationLevel=saturationLevel,
                       saturationFrac=saturationFrac,
                       satVal=satVal, satMaskBit=satMaskBit, isTop=False)

        # A saturated core near the left edge of the detector, close to the
        # previous edge bleed
        x = 10
        y = 900
        halfWidthX = 14
        halfWidthY = 500
        halfWidthColX = 3
        _setExposureSatCore(exposure, x=x, y=y, halfWidthX=halfWidthX,
                            halfWidthY=halfWidthY, satVal=satVal,
                            satMaskBit=satMaskBit)
        # Bottom saturated columns to connect core and edge bleed
        limY = 100
        _setExposureSatColumns(exposure, x=x, y=y, halfWidthX=halfWidthColX,
                               limY=limY, satVal=satVal, satMaskBit=satMaskBit,
                               isTop=False)
        edgeBleedWidth = 70  # this would go beyond the edge of the detector
        edgeBleedWidthLimit = 2*halfWidthColX
        _makeEdgeBleed(exposure, x=x, extentY=extentY,
                       edgeBleedWidth=edgeBleedWidth,
                       edgeBleedConstant=edgeBleedConstant,
                       edgeBleedWidthLimit=edgeBleedWidthLimit,
                       saturationLevel=saturationLevel,
                       saturationFrac=saturationFrac,
                       satVal=satVal, satMaskBit=satMaskBit, isTop=False)

        thresh = afwDetection.Threshold(exposure.mask.getPlaneBitMask("SAT"),
                                        afwDetection.Threshold.BITMASK
                                        )
        fpList = afwDetection.FootprintSet(exposure.mask, thresh).getFootprints()

        satAreas = np.asarray([fp.getArea() for fp in fpList])
        largeAreas, = np.where((satAreas >= 10000)
                               & (satAreas < 100000))

        # We select the only footprint
        fpCore = fpList[largeAreas[0]]

        numPixSatBottomEdgeBeforeCase1 = len(np.where(exposure.mask.array[0, :] == satMaskBit)[0])
        ipIsr.maskITLEdgeBleed(exposure, badAmpDict,
                               fpCore, itlEdgeBleedThreshold=5000.,
                               itlEdgeBleedModelConstant=0.02,
                               saturatedMaskName='SAT')
        numPixSatBottomEdgeAfterCase1 = len(np.where(exposure.mask.array[0, :] == satMaskBit)[0])
        # Check the number of saturated pixels
        self.assertEqual(ipIsr.countMaskedPixels(exposure, 'SAT'), 56871)
        # Check there are more pixels along the bottom edge after masking
        self.assertGreater(numPixSatBottomEdgeAfterCase1, numPixSatBottomEdgeBeforeCase1)

        # Case 2:
        # A saturated core in the top amplifier
        # that has a top edge bleed
        x = 2100
        y = 3100
        halfWidthX = 50
        halfWidthY = 80
        halfWidthColX = 4
        limY = exposure.image.array.shape[0]-100
        _setExposureSatCore(exposure, x=x, y=y, halfWidthX=halfWidthX,
                            halfWidthY=halfWidthY, satVal=satVal,
                            satMaskBit=satMaskBit)
        _setExposureSatColumns(exposure, x=x, y=y, halfWidthX=halfWidthColX,
                               limY=limY, satVal=satVal, satMaskBit=satMaskBit,
                               isTop=True)
        edgeBleedWidth = 100
        edgeBleedWidthLimit = 2*halfWidthColX
        _makeEdgeBleed(exposure, x=x, extentY=extentY,
                       edgeBleedWidth=edgeBleedWidth,
                       edgeBleedConstant=edgeBleedConstant,
                       edgeBleedWidthLimit=edgeBleedWidthLimit,
                       saturationLevel=saturationLevel,
                       saturationFrac=saturationFrac,
                       satVal=satVal, satMaskBit=satMaskBit, isTop=True)

        # Add a saturation sag column around the saturated core
        xSatSag = x
        yminSatSag = y - halfWidthY - 100
        ymaxSatSag = y - halfWidthY
        exposure.image.array[yminSatSag:ymaxSatSag, xSatSag] = satVal*0.7

        # Re-measure the saturated footprints
        fpList = afwDetection.FootprintSet(exposure.mask, thresh).getFootprints()

        satAreas = np.asarray([fp.getArea() for fp in fpList])
        largeAreas, = np.where((satAreas >= 10000)
                               & (satAreas < 100000))

        for largeAreasIndex in largeAreas:
            # We select the new edge bleed footprint
            if fpList[largeAreasIndex].getBBox().contains(x, y):
                fpCore = fpList[largeAreasIndex]

                numPixSatTopEdgeBeforeCase2 = len(np.where(exposure.mask.array[-1, :] == satMaskBit)[0])
                ipIsr.maskITLEdgeBleed(exposure, badAmpDict,
                                       fpCore, itlEdgeBleedThreshold=5000.,
                                       itlEdgeBleedModelConstant=0.02,
                                       saturatedMaskName='SAT')
                numPixSatTopEdgeAfterCase2 = len(np.where(exposure.mask.array[-1, :] == satMaskBit)[0])
                # Check the number of saturated pixels
                self.assertEqual(ipIsr.countMaskedPixels(exposure, 'SAT'), 84649)
                # Check there are more pixels along the bottom edge
                # after masking
                self.assertGreater(numPixSatTopEdgeAfterCase2, numPixSatTopEdgeBeforeCase2)

                # This will mask the whole column containing the saturation sag
                ipIsrFunctions.maskITLSatSag(exposure, fpCore, saturatedMaskName='SAT')
                numPixColumnMask = len(np.where(exposure.mask.array[:, xSatSag] == satMaskBit)[0])
                # Test that the saturation sag column is masked
                self.assertEqual(numPixColumnMask, exposure.image.array.shape[0])

        # Case 3:
        # A saturated core with an edge bleed on both edges
        x = 3020
        y = 1000
        halfWidthX = 70
        halfWidthY = 150
        halfWidthColX = 7
        _setExposureSatCore(exposure, x=x, y=y, halfWidthX=halfWidthX,
                            halfWidthY=halfWidthY, satVal=satVal,
                            satMaskBit=satMaskBit)
        # Top saturated columns to connect core and edge bleed
        limY = exposure.image.array.shape[0]-100
        _setExposureSatColumns(exposure, x=x, y=y, halfWidthX=halfWidthColX,
                               limY=limY, satVal=satVal, satMaskBit=satMaskBit,
                               isTop=True)
        # Bottom saturated columns to connect core and edge bleed
        limY = 100
        _setExposureSatColumns(exposure, x=x, y=y, halfWidthX=halfWidthColX,
                               limY=limY, satVal=satVal, satMaskBit=satMaskBit,
                               isTop=False)
        edgeBleedWidth = 150
        edgeBleedWidthLimit = 2*halfWidthColX
        _makeEdgeBleed(exposure, x=x, extentY=extentY,
                       edgeBleedWidth=edgeBleedWidth,
                       edgeBleedConstant=edgeBleedConstant,
                       edgeBleedWidthLimit=edgeBleedWidthLimit,
                       saturationLevel=saturationLevel,
                       saturationFrac=saturationFrac,
                       satVal=satVal, satMaskBit=satMaskBit, isTop=True)
        _makeEdgeBleed(exposure, x=x, extentY=extentY,
                       edgeBleedWidth=edgeBleedWidth,
                       edgeBleedConstant=edgeBleedConstant,
                       edgeBleedWidthLimit=edgeBleedWidthLimit,
                       saturationLevel=saturationLevel,
                       saturationFrac=saturationFrac,
                       satVal=satVal, satMaskBit=satMaskBit, isTop=False)

        # Re-measure the saturated footprints
        fpList = afwDetection.FootprintSet(exposure.mask, thresh).getFootprints()

        satAreas = np.asarray([fp.getArea() for fp in fpList])
        largeAreas, = np.where((satAreas >= 10000)
                               & (satAreas < 100000))

        for largeAreasIndex in largeAreas:
            # We select the new edge bleed footprint
            if fpList[largeAreasIndex].getBBox().contains(x, y):
                fpCore = fpList[largeAreasIndex]

                # Number of saturated pixels at the bottom edge
                # before applying masking
                numPixSatBottomEdgeBefore = len(np.where(exposure.mask.array[0, :] == satMaskBit)[0])
                # Number of saturated pixels at the top edge
                # before applying masking
                numPixSatTopEdgeBefore = len(np.where(exposure.mask.array[-1, :] == satMaskBit)[0])
                # Apply edge bleed masking
                ipIsr.maskITLEdgeBleed(exposure, badAmpDict,
                                       fpCore, itlEdgeBleedThreshold=5000.,
                                       itlEdgeBleedModelConstant=0.02,
                                       saturatedMaskName='SAT')
                # Number of saturated pixels at the bottom edge
                # after applying edge bleed masking
                numPixSatBottomEdgeAfter = len(np.where(exposure.mask.array[0, :] == satMaskBit)[0])
                numPixSatTopEdgeAfter = len(np.where(exposure.mask.array[-1, :] == satMaskBit)[0])

                # Check the number of saturated pixels
                self.assertEqual(ipIsr.countMaskedPixels(exposure, 'SAT'), 228629)
                # Check there are more pixels along the bottom edge
                # after masking
                self.assertGreater(numPixSatBottomEdgeAfter, numPixSatBottomEdgeBefore)
                # Same with top edge
                self.assertGreater(numPixSatTopEdgeAfter, numPixSatTopEdgeBefore)

        # Case 4:
        # A saturated core near the right edge of the detector
        x = 4060
        y = 3000
        halfWidthX = 30
        halfWidthY = 500
        halfWidthColX = 3
        _setExposureSatCore(exposure, x=x, y=y, halfWidthX=halfWidthX,
                            halfWidthY=halfWidthY, satVal=satVal,
                            satMaskBit=satMaskBit)
        # Top saturated columns to connect core and edge bleed
        limY = exposure.image.array.shape[0]-100
        _setExposureSatColumns(exposure, x=x, y=y, halfWidthX=halfWidthColX,
                               limY=limY, satVal=satVal, satMaskBit=satMaskBit,
                               isTop=True)
        edgeBleedWidth = 100  # this would go beyond the edge of the detector
        edgeBleedWidthLimit = 2*halfWidthColX
        _makeEdgeBleed(exposure, x=x, extentY=extentY,
                       edgeBleedWidth=edgeBleedWidth,
                       edgeBleedConstant=edgeBleedConstant,
                       edgeBleedWidthLimit=edgeBleedWidthLimit,
                       saturationLevel=saturationLevel,
                       saturationFrac=saturationFrac,
                       satVal=satVal, satMaskBit=satMaskBit, isTop=True)

        # Re-measure the saturated footprints
        fpList = afwDetection.FootprintSet(exposure.mask, thresh).getFootprints()

        satAreas = np.asarray([fp.getArea() for fp in fpList])
        largeAreas, = np.where((satAreas >= 10000)
                               & (satAreas < 100000))

        for largeAreasIndex in largeAreas:
            # We select the new edge bleed footprint
            if fpList[largeAreasIndex].getBBox().contains(x, y):
                fpCore = fpList[largeAreasIndex]

                # Number of saturated pixels at the bottom edge
                # before applying masking
                numPixSatTopEdgeBefore = len(np.where(exposure.mask.array[-1, :] == satMaskBit)[0])
                # Apply edge bleed masking
                ipIsr.maskITLEdgeBleed(exposure, badAmpDict,
                                       fpCore, itlEdgeBleedThreshold=5000.,
                                       itlEdgeBleedModelConstant=0.02,
                                       saturatedMaskName='SAT')
                # Number of saturated pixels at the bottom edge
                # after applying edge bleed masking
                numPixSatTopEdgeAfter = len(np.where(exposure.mask.array[-1, :] == satMaskBit)[0])

                # Check the number of saturated pixels
                self.assertEqual(ipIsr.countMaskedPixels(exposure, 'SAT'), 289136)
                # Check there are more pixels along the bottom edge
                # after masking
                self.assertGreater(numPixSatTopEdgeAfter, numPixSatTopEdgeBefore)

    def test_itlDipMasking(self):
        """Test the ITL dip masking."""
        detector = MockITLDetector()
        exposure = MockITLExposure(detector.getBBox())
        exposure.setDetector(detector)
        exposure.mask.array[:, :] = 0
        exposure.mask.addMaskPlane("ITL_DIP")
        dipMaskValue = exposure.mask.getPlaneBitMask(["SUSPECT", "ITL_DIP"])

        detectorConfig = ipIsr.overscanAmpConfig.OverscanDetectorConfig()
        detectorConfig.itlDipBackgroundFraction = 0.0025

        # Set the background to 1000 electrons.
        exposure.image.array[:, :] = 1000.0

        # Add some saturated masks, below and above the trigger.
        satMaskValue = exposure.mask.getPlaneBitMask("SAT")

        # Above the threshold in width/height.
        exposure.mask.array[
            500: 500 + detectorConfig.itlDipMinHeight * 2,
            1000: 1000 + detectorConfig.itlDipMinWidth * 2
        ] |= satMaskValue

        # At the threshold in width/height.
        exposure.mask.array[
            1200: 1200 + detectorConfig.itlDipMinHeight + 1,
            1500: 1500 + detectorConfig.itlDipMinWidth + 1
        ] |= satMaskValue

        # Below the threshold in width.
        exposure.mask.array[
            1600: 1600 + detectorConfig.itlDipMinHeight + 1,
            2500: 2500 + detectorConfig.itlDipMinWidth // 2
        ] |= satMaskValue

        # Below the threshold in height.
        exposure.mask.array[
            2500: 2500 + detectorConfig.itlDipMinHeight // 2,
            3000: 3000 + detectorConfig.itlDipMinWidth + 1
        ] |= satMaskValue

        with self.assertLogs(level=logging.INFO) as cm:
            ipIsr.isrFunctions.maskITLDip(exposure, detectorConfig)
        self.assertEqual(len(cm[1]), 2)
        self.assertIn("Found ITL dip (width 30; bkg 1000.00); masking column 992 to 1037", cm[1][0])
        self.assertIn("Found ITL dip (width 16; bkg 1000.00); masking column 1495 to 1519", cm[1][1])

        # This includes the scaled edges
        np.testing.assert_array_equal(exposure.mask.array[:, 992: 1038] & dipMaskValue, dipMaskValue)
        np.testing.assert_array_equal(exposure.mask.array[:, 991] & dipMaskValue, 0)
        np.testing.assert_array_equal(exposure.mask.array[:, 1038] & dipMaskValue, 0)

        np.testing.assert_array_equal(exposure.mask.array[:, 1495: 1520] & dipMaskValue, dipMaskValue)
        np.testing.assert_array_equal(exposure.mask.array[:, 1494] & dipMaskValue, 0)
        np.testing.assert_array_equal(exposure.mask.array[:, 1520] & dipMaskValue, 0)

        # The other two should not be masked.
        np.testing.assert_array_equal(exposure.mask.array[:, 2500] & dipMaskValue, 0)
        np.testing.assert_array_equal(exposure.mask.array[:, 3000] & dipMaskValue, 0)

        # Change the background to a much lower value, only 1 should be masked.
        exposure.mask.array[:, :] &= ~dipMaskValue
        exposure.image.array[:, :] = 50.0

        with self.assertLogs(level=logging.INFO) as cm:
            ipIsr.isrFunctions.maskITLDip(exposure, detectorConfig)
        self.assertEqual(len(cm[1]), 1)
        self.assertIn("Found ITL dip (width 30; bkg 50.00); masking column 992 to 1037", cm[1][0])

        # This one should be masked.
        np.testing.assert_array_equal(exposure.mask.array[:, 992: 1038] & dipMaskValue, dipMaskValue)
        np.testing.assert_array_equal(exposure.mask.array[:, 991] & dipMaskValue, 0)
        np.testing.assert_array_equal(exposure.mask.array[:, 1038] & dipMaskValue, 0)

        # The other three should not be masked.
        np.testing.assert_array_equal(exposure.mask.array[:, 1500] & dipMaskValue, 0)
        np.testing.assert_array_equal(exposure.mask.array[:, 2500] & dipMaskValue, 0)
        np.testing.assert_array_equal(exposure.mask.array[:, 3000] & dipMaskValue, 0)

        # And blow things up so that masking is not applied.
        # This needs to be in several sections because of the max width.
        # We additionally avoid the edge. Note that there are rounding
        # offsets that occur due to the scaling such that this setting
        # will mask 42 columns per block or a total of 504.
        exposure.mask.array[:, :] = 0
        detectorConfig.itlDipWidthScale = 1.0
        for i in range(12):
            exposure.mask.array[:, (i + 1)*100: (i + 1)*100 + 41] |= satMaskValue

        maskBak = exposure.mask.array.copy()
        with self.assertLogs(level=logging.WARNING) as cm:
            ipIsr.isrFunctions.maskITLDip(exposure, detectorConfig)
        self.assertEqual(len(cm[1]), 1)
        self.assertIn("Too many (504) columns would be masked", cm[1][0])
        np.testing.assert_array_equal(exposure.mask.array, maskBak)

    def test_interpolateFromMask(self):
        """Expect number of interpolated pixels to be non-zero.
        """
        ipIsr.makeThresholdMask(self.mi, 200, growFootprints=2,
                                maskName='SAT')
        for growFootprints in range(0, 3):
            for useLegacyInterp in (False, True):
                interpMaskedImage = ipIsr.interpolateFromMask(self.mi, 2.0,
                                                              growSaturatedFootprints=growFootprints,
                                                              maskNameList=['SAT'],
                                                              useLegacyInterp=useLegacyInterp)
                numBit = ipIsr.countMaskedPixels(interpMaskedImage, "INTRP")
                if growFootprints == 0 and not useLegacyInterp:
                    # All pixel need to be interpolated over. There is
                    # no external information to interpolate over. In
                    # the GP code in this case, it is not doing
                    # anything.
                    self.assertEqual(numBit, 0,
                                     msg=f"interpolateFromMask with growFootprints={growFootprints}")
                else:
                    self.assertEqual(numBit, 40800,
                                     msg=f"interpolateFromMask with growFootprints={growFootprints}")

    def test_saturationCorrectionInterpolate(self):
        """Expect number of mask pixels with SAT marked to be non-zero.
        """
        corrMaskedImage = ipIsr.saturationCorrection(self.mi, 200, 2.0,
                                                     growFootprints=2, interpolate=True,
                                                     maskName='SAT')
        numBit = ipIsr.countMaskedPixels(corrMaskedImage, "SAT")
        self.assertEqual(numBit, 40800)

    def test_saturationCorrectionNoInterpolate(self):
        """Expect number of mask pixels with SAT marked to be non-zero.
        """
        corrMaskedImage = ipIsr.saturationCorrection(self.mi, 200, 2.0,
                                                     growFootprints=2, interpolate=False,
                                                     maskName='SAT')
        numBit = ipIsr.countMaskedPixels(corrMaskedImage, "SAT")
        self.assertEqual(numBit, 40800)

    def test_trimToMatchCalibBBox(self):
        """Expect bounding boxes to match.
        """
        darkExp = isrMock.DarkMock().run()
        darkMi = darkExp.getMaskedImage()

        nEdge = 2
        darkMi = darkMi[nEdge:-nEdge, nEdge:-nEdge, afwImage.LOCAL]
        newInput = ipIsr.trimToMatchCalibBBox(self.mi, darkMi)

        self.assertEqual(newInput.getImage().getBBox(), darkMi.getImage().getBBox())

    def test_darkCorrection(self):
        """Expect round-trip application to be equal.
        Expect RuntimeError if sizes are different.
        """
        darkExp = isrMock.DarkMock().run()
        darkMi = darkExp.getMaskedImage()

        mi = self.mi.clone()

        # The `invert` parameter controls the direction of the
        # application.  This will apply, and un-apply the dark.
        ipIsr.darkCorrection(self.mi, darkMi, 1.0, 1.0, trimToFit=True)
        ipIsr.darkCorrection(self.mi, darkMi, 1.0, 1.0, trimToFit=True, invert=True)

        self.assertMaskedImagesAlmostEqual(self.mi, mi, atol=1e-3)

        darkMi = darkMi[1:-1, 1:-1, afwImage.LOCAL]
        with self.assertRaises(RuntimeError):
            ipIsr.darkCorrection(self.mi, darkMi, 1.0, 1.0, trimToFit=False)

    def test_biasCorrection(self):
        """Expect smaller median image value after.
        Expect RuntimeError if sizes are different.
        """
        biasExp = isrMock.BiasMock().run()
        biasMi = biasExp.getMaskedImage()

        mi = self.mi.clone()
        ipIsr.biasCorrection(self.mi, biasMi, trimToFit=True)
        self.assertLess(computeImageMedianAndStd(self.mi.getImage())[0],
                        computeImageMedianAndStd(mi.getImage())[0])

        biasMi = biasMi[1:-1, 1:-1, afwImage.LOCAL]
        with self.assertRaises(RuntimeError):
            ipIsr.biasCorrection(self.mi, biasMi, trimToFit=False)

    def test_flatCorrection(self):
        """Expect round-trip application to be equal.
        Expect RuntimeError if sizes are different.
        """
        flatExp = isrMock.FlatMock().run()
        flatMi = flatExp.getMaskedImage()

        mi = self.mi.clone()
        for scaling in ('USER', 'MEAN', 'MEDIAN'):
            # The `invert` parameter controls the direction of the
            # application.  This will apply, and un-apply the flat.
            ipIsr.flatCorrection(self.mi, flatMi, scaling, userScale=1.0, trimToFit=True)
            ipIsr.flatCorrection(self.mi, flatMi, scaling, userScale=1.0,
                                 trimToFit=True, invert=True)

            self.assertMaskedImagesAlmostEqual(self.mi, mi, atol=1e-3,
                                               msg=f"flatCorrection with scaling {scaling}")

        flatMi = flatMi[1:-1, 1:-1, afwImage.LOCAL]
        with self.assertRaises(RuntimeError):
            ipIsr.flatCorrection(self.mi, flatMi, 'USER', userScale=1.0, trimToFit=False)

    def test_flatCorrectionUnknown(self):
        """Raise if an unknown scaling is used.

        The `scaling` parameter must be a known type.  If not, the
        flat correction will raise a RuntimeError.
        """
        flatExp = isrMock.FlatMock().run()
        flatMi = flatExp.getMaskedImage()

        with self.assertRaises(RuntimeError):
            ipIsr.flatCorrection(self.mi, flatMi, "UNKNOWN", userScale=1.0, trimToFit=True)

    def test_illumCorrection(self):
        """Expect larger median value after.
        Expect RuntimeError if sizes are different.
        """
        flatExp = isrMock.FlatMock().run()
        flatMi = flatExp.getMaskedImage()

        mi = self.mi.clone()
        ipIsr.illuminationCorrection(self.mi, flatMi, 1.0)
        self.assertGreater(computeImageMedianAndStd(self.mi.getImage())[0],
                           computeImageMedianAndStd(mi.getImage())[0])

        flatMi = flatMi[1:-1, 1:-1, afwImage.LOCAL]
        with self.assertRaises(RuntimeError):
            ipIsr.illuminationCorrection(self.mi, flatMi, 1.0, trimToFit=False)

    def test_brighterFatterCorrection(self):
        """Expect smoother image/smaller std before.
        """
        bfKern = isrMock.BfKernelMock().run()

        before = computeImageMedianAndStd(self.inputExp.getImage())
        ipIsr.brighterFatterCorrection(self.inputExp, bfKern, 10, 1e-2, False)
        after = computeImageMedianAndStd(self.inputExp.getImage())

        self.assertLess(before[1], after[1])

    def test_gainContext(self):
        """Expect image to be unmodified before and after
        """
        mi = self.inputExp.getMaskedImage().clone()
        with ipIsr.gainContext(self.inputExp, self.inputExp.getImage(), apply=True):
            pass

        self.assertIsNotNone(mi)
        self.assertMaskedImagesAlmostEqual(self.inputExp.getMaskedImage(), mi)

    def test_widenSaturationTrails(self):
        """Expect more mask pixels with SAT set after.
        """
        numBitBefore = ipIsr.countMaskedPixels(self.mi, "SAT")

        ipIsr.widenSaturationTrails(self.mi.getMask())
        numBitAfter = ipIsr.countMaskedPixels(self.mi, "SAT")

        self.assertGreaterEqual(numBitAfter, numBitBefore)

    def test_setBadRegions(self):
        """Expect RuntimeError if improper statistic given.
        Expect a float value otherwise.
        """
        for badStatistic in ('MEDIAN', 'MEANCLIP', 'UNKNOWN'):
            if badStatistic == 'UNKNOWN':
                with self.assertRaises(RuntimeError,
                                       msg=f"setBadRegions did not fail for stat {badStatistic}"):
                    nBad, value = ipIsr.setBadRegions(self.inputExp, badStatistic=badStatistic)
            else:
                nBad, value = ipIsr.setBadRegions(self.inputExp, badStatistic=badStatistic)
                self.assertGreaterEqual(abs(value), 0.0,
                                        msg=f"setBadRegions did not find valid value for stat {badStatistic}")

    def test_attachTransmissionCurve(self):
        """Expect no failure and non-None output from attachTransmissionCurve.
        """
        curve = isrMock.TransmissionMock().run()
        combined = ipIsr.attachTransmissionCurve(self.inputExp,
                                                 opticsTransmission=curve,
                                                 filterTransmission=curve,
                                                 sensorTransmission=curve,
                                                 atmosphereTransmission=curve)
        # DM-19707: ip_isr functionality not fully tested by unit tests
        self.assertIsNotNone(combined)

    def test_attachTransmissionCurve_None(self):
        """Expect no failure and non-None output from attachTransmissionCurve.
        """
        curve = None
        combined = ipIsr.attachTransmissionCurve(self.inputExp,
                                                 opticsTransmission=curve,
                                                 filterTransmission=curve,
                                                 sensorTransmission=curve,
                                                 atmosphereTransmission=curve)
        # DM-19707: ip_isr functionality not fully tested by unit tests
        self.assertIsNotNone(combined)

    def test_countMaskedPixels(self):
        mockImageConfig = isrMock.IsrMock.ConfigClass()

        # flatDrop is not really relevant as we replace the data
        # but good to note it in case we change how this image is made
        mockImageConfig.flatDrop = 0.99999
        mockImageConfig.isTrimmed = True

        flatExp = isrMock.FlatMock(config=mockImageConfig).run()
        (shapeY, shapeX) = flatExp.getDimensions()

        rng = np.random.RandomState(0)
        flatMean = 1000
        flatWidth = np.sqrt(flatMean)
        flatData = rng.normal(flatMean, flatWidth, (shapeX, shapeY))
        flatExp.image.array[:] = flatData

        exp = flatExp.clone()
        mi = exp.maskedImage
        self.assertEqual(ipIsr.countMaskedPixels(mi, 'NO_DATA'), 0)
        self.assertEqual(ipIsr.countMaskedPixels(mi, 'BAD'), 0)

        NODATABIT = mi.mask.getPlaneBitMask("NO_DATA")
        noDataBox = geom.Box2I(geom.Point2I(31, 49), geom.Extent2I(3, 6))
        mi.mask[noDataBox] |= NODATABIT

        self.assertEqual(ipIsr.countMaskedPixels(mi, 'NO_DATA'), noDataBox.getArea())
        self.assertEqual(ipIsr.countMaskedPixels(mi, 'BAD'), 0)

        mi.mask[noDataBox] ^= NODATABIT  # XOR to reset what we did
        self.assertEqual(ipIsr.countMaskedPixels(mi, 'NO_DATA'), 0)

    def test_getExposureGains(self):
        exposure = self.inputExp.clone()
        metadata = exposure.metadata
        amps = exposure.getDetector().getAmplifiers()

        # Check the default values
        values = ipIsr.getExposureGains(exposure).values()
        self.assertEqual(list(values), len(amps)*[1.0])

        # Set values using old ISR keys
        for amp in amps:
            metadata[f"LSST GAIN {amp.getName()}"] = 2.0
        values = ipIsr.getExposureGains(exposure).values()
        self.assertEqual(list(values), len(amps)*[2.0])

        # Set values using new ISR keys
        for amp in amps:
            metadata[f"LSST ISR GAIN {amp.getName()}"] = 3.0
        values = ipIsr.getExposureGains(exposure).values()
        self.assertEqual(list(values), len(amps)*[3.0])

    def test_getExposureReadNoises(self):
        exposure = self.inputExp.clone()
        metadata = exposure.metadata
        amps = exposure.getDetector().getAmplifiers()

        # Check the default values
        values = ipIsr.getExposureReadNoises(exposure).values()
        self.assertEqual(list(values), len(amps)*[10.0])

        # Set values using old ISR keys
        for amp in amps:
            metadata[f"LSST READNOISE {amp.getName()}"] = 2.0
        values = ipIsr.getExposureReadNoises(exposure).values()
        self.assertEqual(list(values), len(amps)*[2.0])

        # Set values using new ISR keys
        for amp in amps:
            metadata[f"LSST ISR READNOISE {amp.getName()}"] = 3.0
        values = ipIsr.getExposureReadNoises(exposure).values()
        self.assertEqual(list(values), len(amps)*[3.0])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
