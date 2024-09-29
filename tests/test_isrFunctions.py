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
import numpy as np

import lsst.geom as geom
import lsst.afw.image as afwImage
import lsst.utils.tests
import lsst.ip.isr as ipIsr
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
