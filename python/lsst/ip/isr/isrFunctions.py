#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
import math
import numpy
from deprecated.sphinx import deprecated

import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.pex.exceptions as pexExcept
import lsst.afw.cameraGeom as camGeom

from lsst.afw.geom.wcsUtils import makeDistortedTanWcs
from lsst.meas.algorithms.detection import SourceDetectionTask
from lsst.pipe.base import Struct

from contextlib import contextmanager


def createPsf(fwhm):
    """Make a double Gaussian PSF.

    Parameters
    ----------
    fwhm : scalar
        FWHM of double Gaussian smoothing kernel.

    Returns
    -------
    psf : `lsst.meas.algorithms.DoubleGaussianPsf`
        The created smoothing kernel.
    """
    ksize = 4*int(fwhm) + 1
    return measAlg.DoubleGaussianPsf(ksize, ksize, fwhm/(2*math.sqrt(2*math.log(2))))


def transposeMaskedImage(maskedImage):
    """Make a transposed copy of a masked image.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.

    Returns
    -------
    transposed : `lsst.afw.image.MaskedImage`
        The transposed copy of the input image.
    """
    transposed = maskedImage.Factory(lsst.geom.Extent2I(maskedImage.getHeight(), maskedImage.getWidth()))
    transposed.getImage().getArray()[:] = maskedImage.getImage().getArray().T
    transposed.getMask().getArray()[:] = maskedImage.getMask().getArray().T
    transposed.getVariance().getArray()[:] = maskedImage.getVariance().getArray().T
    return transposed


def interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=None):
    """Interpolate over defects specified in a defect list.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.
    defectList : `lsst.meas.algorithms.Defects`
        List of defects to interpolate over.
    fwhm : scalar
        FWHM of double Gaussian smoothing kernel.
    fallbackValue : scalar, optional
        Fallback value if an interpolated value cannot be determined.
        If None, then the clipped mean of the image is used.
    """
    psf = createPsf(fwhm)
    if fallbackValue is None:
        fallbackValue = afwMath.makeStatistics(maskedImage.getImage(), afwMath.MEANCLIP).getValue()
    if 'INTRP' not in maskedImage.getMask().getMaskPlaneDict():
        maskedImage.getMask().addMaskPlane('INTRP')
    measAlg.interpolateOverDefects(maskedImage, psf, defectList, fallbackValue, True)
    return maskedImage


@deprecated(reason="Replaced by Defects.fromFootPrintList() (will be removed after v18)",
            category=FutureWarning)
def defectListFromFootprintList(fpList):
    """Compute a defect list from a footprint list, optionally growing the footprints.

    Parameters
    ----------
    fpList : `list` of `lsst.afw.detection.Footprint`
        Footprint list to process.

    Returns
    -------
    defectList : `lsst.meas.algorithms.Defects`
        List of defects.
    """
    return measAlg.Defects.fromFootprintList(fpList)

@deprecated(reason="Replaced by Defects.transpose() (will be removed after v18)",
            category=FutureWarning)
def transposeDefectList(defectList):
    """Make a transposed copy of a defect list.

    Parameters
    ----------
    defectList : `lsst.meas.algorithms.Defects`
        Input list of defects.

    Returns
    -------
    retDefectList : `lsst.meas.algorithms.Defects`
        Transposed list of defects.
    """
    if isinstance(defectList, measAlg.Defects):
        return defectList.transpose()
    return measAlg.Defects(defectList).transpose()


@deprecated(reason="Replaced by Defects.maskPixels() (will be removed after v18)",
            category=FutureWarning)
def maskPixelsFromDefectList(maskedImage, defectList, maskName='BAD'):
    """Set mask plane based on a defect list.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.  Only the mask plane is updated.
    defectList : `lsst.meas.algorithms.Defects`
        Defect list to mask.
    maskName : str, optional
        Mask plane name to use.
    """
    return lsst.meas.algorithms.Defects(defectList).maskPixels(maskedImage, maskName=maskName)


@deprecated(reason="Replaced by Defects.fromMask() (will be removed after v18)",
            category=FutureWarning)
def getDefectListFromMask(maskedImage, maskName):
    """Compute a defect list from a specified mask plane.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.
    maskName : `str` or `list`
        Mask plane name, or list of names to convert.

    Returns
    -------
    defectList : `lsst.meas.algorithms.Defects`
        Defect list constructed from masked pixels.
    """
    return measAlg.Defects.fromMask(maskedImage, maskName)


def makeThresholdMask(maskedImage, threshold, growFootprints=1, maskName='SAT'):
    """Mask pixels based on threshold detection.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.  Only the mask plane is updated.
    threshold : scalar
        Detection threshold.
    growFootprints : scalar, optional
        Number of pixels to grow footprints of detected regions.
    maskName : str, optional
        Mask plane name, or list of names to convert

    Returns
    -------
    defectList : `lsst.meas.algorithms.Defects`
        Defect list constructed from pixels above the threshold.
    """
    # find saturated regions
    thresh = afwDetection.Threshold(threshold)
    fs = afwDetection.FootprintSet(maskedImage, thresh)

    if growFootprints > 0:
        fs = afwDetection.FootprintSet(fs, rGrow=growFootprints, isotropic=False)
    fpList = fs.getFootprints()

    # set mask
    mask = maskedImage.getMask()
    bitmask = mask.getPlaneBitMask(maskName)
    afwDetection.setMaskFromFootprintList(mask, fpList, bitmask)

    return measAlg.Defects.fromFootprintList(fpList)


def interpolateFromMask(maskedImage, fwhm, growSaturatedFootprints=1,
                        maskNameList=['SAT'], fallbackValue=None):
    """Interpolate over defects identified by a particular set of mask planes.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.
    fwhm : scalar
        FWHM of double Gaussian smoothing kernel.
    growSaturatedFootprints : scalar, optional
        Number of pixels to grow footprints for saturated pixels.
    maskNameList : `List` of `str`, optional
        Mask plane name.
    fallbackValue : scalar, optional
        Value of last resort for interpolation.
    """
    mask = maskedImage.getMask()

    if growSaturatedFootprints > 0 and "SAT" in maskNameList:
        thresh = afwDetection.Threshold(mask.getPlaneBitMask("SAT"), afwDetection.Threshold.BITMASK)
        fpSet = afwDetection.FootprintSet(mask, thresh)
        # If we are interpolating over an area larger than the original masked region, we need
        # to expand the original mask bit to the full area to explain why we interpolated there.
        fpSet = afwDetection.FootprintSet(fpSet, rGrow=growSaturatedFootprints, isotropic=False)
        fpSet.setMask(mask, "SAT")

    thresh = afwDetection.Threshold(mask.getPlaneBitMask(maskNameList), afwDetection.Threshold.BITMASK)
    fpSet = afwDetection.FootprintSet(mask, thresh)
    defectList = measAlg.Defects.fromFootprintList(fpSet.getFootprints())

    interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=fallbackValue)

    return maskedImage


def saturationCorrection(maskedImage, saturation, fwhm, growFootprints=1, interpolate=True, maskName='SAT',
                         fallbackValue=None):
    """Mark saturated pixels and optionally interpolate over them

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.
    saturation  : scalar
        Saturation level used as the detection threshold.
    fwhm : scalar
        FWHM of double Gaussian smoothing kernel.
    growFootprints : scalar, optional
        Number of pixels to grow footprints of detected regions.
    interpolate : Bool, optional
        If True, saturated pixels are interpolated over.
    maskName : str, optional
        Mask plane name.
    fallbackValue : scalar, optional
        Value of last resort for interpolation.
    """
    defectList = makeThresholdMask(
        maskedImage=maskedImage,
        threshold=saturation,
        growFootprints=growFootprints,
        maskName=maskName,
    )
    if interpolate:
        interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=fallbackValue)

    return maskedImage


def trimToMatchCalibBBox(rawMaskedImage, calibMaskedImage):
    """Compute number of edge trim pixels to match the calibration data.

    Use the dimension difference between the raw exposure and the
    calibration exposure to compute the edge trim pixels.  This trim
    is applied symmetrically, with the same number of pixels masked on
    each side.

    Parameters
    ----------
    rawMaskedImage : `lsst.afw.image.MaskedImage`
        Image to trim.
    calibMaskedImage : `lsst.afw.image.MaskedImage`
        Calibration image to draw new bounding box from.

    Returns
    -------
    replacementMaskedImage : `lsst.afw.image.MaskedImage`
        ``rawMaskedImage`` trimmed to the appropriate size
    Raises
    ------
    RuntimeError
       Rasied if ``rawMaskedImage`` cannot be symmetrically trimmed to
       match ``calibMaskedImage``.
    """
    nx, ny = rawMaskedImage.getBBox().getDimensions() - calibMaskedImage.getBBox().getDimensions()
    if nx != ny:
        raise RuntimeError("Raw and calib maskedImages are trimmed differently in X and Y.")
    if nx % 2 != 0:
        raise RuntimeError("Calibration maskedImage is trimmed unevenly in X.")
    if nx < 0:
        raise RuntimeError("Calibration maskedImage is larger than raw data.")

    nEdge = nx//2
    if nEdge > 0:
        replacementMaskedImage = rawMaskedImage[nEdge:-nEdge, nEdge:-nEdge, afwImage.LOCAL]
        SourceDetectionTask.setEdgeBits(
            rawMaskedImage,
            replacementMaskedImage.getBBox(),
            rawMaskedImage.getMask().getPlaneBitMask("EDGE")
        )
    else:
        replacementMaskedImage = rawMaskedImage

    return replacementMaskedImage


def biasCorrection(maskedImage, biasMaskedImage, trimToFit=False):
    """Apply bias correction in place.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
       Image to process.  The image is modified by this method.
    biasMaskedImage : `lsst.afw.image.MaskedImage`
        Bias image of the same size as ``maskedImage``
    trimToFit : `Bool`, optional
        If True, raw data is symmetrically trimmed to match
        calibration size.

    Raises
    ------
    RuntimeError
        Raised if ``maskedImage`` and ``biasMaskedImage`` do not have
        the same size.

    """
    if trimToFit:
        maskedImage = trimToMatchCalibBBox(maskedImage, biasMaskedImage)

    if maskedImage.getBBox(afwImage.LOCAL) != biasMaskedImage.getBBox(afwImage.LOCAL):
        raise RuntimeError("maskedImage bbox %s != biasMaskedImage bbox %s" %
                           (maskedImage.getBBox(afwImage.LOCAL), biasMaskedImage.getBBox(afwImage.LOCAL)))
    maskedImage -= biasMaskedImage


def darkCorrection(maskedImage, darkMaskedImage, expScale, darkScale, invert=False, trimToFit=False):
    """Apply dark correction in place.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
       Image to process.  The image is modified by this method.
    darkMaskedImage : `lsst.afw.image.MaskedImage`
        Dark image of the same size as ``maskedImage``.
    expScale : scalar
        Dark exposure time for ``maskedImage``.
    darkScale : scalar
        Dark exposure time for ``darkMaskedImage``.
    invert : `Bool`, optional
        If True, re-add the dark to an already corrected image.
    trimToFit : `Bool`, optional
        If True, raw data is symmetrically trimmed to match
        calibration size.

    Raises
    ------
    RuntimeError
        Raised if ``maskedImage`` and ``darkMaskedImage`` do not have
        the same size.

    Notes
    -----
    The dark correction is applied by calculating:
        maskedImage -= dark * expScaling / darkScaling
    """
    if trimToFit:
        maskedImage = trimToMatchCalibBBox(maskedImage, darkMaskedImage)

    if maskedImage.getBBox(afwImage.LOCAL) != darkMaskedImage.getBBox(afwImage.LOCAL):
        raise RuntimeError("maskedImage bbox %s != darkMaskedImage bbox %s" %
                           (maskedImage.getBBox(afwImage.LOCAL), darkMaskedImage.getBBox(afwImage.LOCAL)))

    scale = expScale / darkScale
    if not invert:
        maskedImage.scaledMinus(scale, darkMaskedImage)
    else:
        maskedImage.scaledPlus(scale, darkMaskedImage)


def updateVariance(maskedImage, gain, readNoise):
    """Set the variance plane based on the image plane.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.  The variance plane is modified.
    gain : scalar
        The amplifier gain in electrons/ADU.
    readNoise : scalar
        The amplifier read nmoise in ADU/pixel.
    """
    var = maskedImage.getVariance()
    var[:] = maskedImage.getImage()
    var /= gain
    var += readNoise**2


def flatCorrection(maskedImage, flatMaskedImage, scalingType, userScale=1.0, invert=False, trimToFit=False):
    """Apply flat correction in place.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.  The image is modified.
    flatMaskedImage : `lsst.afw.image.MaskedImage`
        Flat image of the same size as ``maskedImage``
    scalingType : str
        Flat scale computation method.  Allowed values are 'MEAN',
        'MEDIAN', or 'USER'.
    userScale : scalar, optional
        Scale to use if ``scalingType``='USER'.
    invert : `Bool`, optional
        If True, unflatten an already flattened image.
    trimToFit : `Bool`, optional
        If True, raw data is symmetrically trimmed to match
        calibration size.

    Raises
    ------
    RuntimeError
        Raised if ``maskedImage`` and ``flatMaskedImage`` do not have
        the same size or if ``scalingType`` is not an allowed value.
    """
    if trimToFit:
        maskedImage = trimToMatchCalibBBox(maskedImage, flatMaskedImage)

    if maskedImage.getBBox(afwImage.LOCAL) != flatMaskedImage.getBBox(afwImage.LOCAL):
        raise RuntimeError("maskedImage bbox %s != flatMaskedImage bbox %s" %
                           (maskedImage.getBBox(afwImage.LOCAL), flatMaskedImage.getBBox(afwImage.LOCAL)))

    # Figure out scale from the data
    # Ideally the flats are normalized by the calibration product pipeline, but this allows some flexibility
    # in the case that the flat is created by some other mechanism.
    if scalingType in ('MEAN', 'MEDIAN'):
        scalingType = afwMath.stringToStatisticsProperty(scalingType)
        flatScale = afwMath.makeStatistics(flatMaskedImage.image, scalingType).getValue()
    elif scalingType == 'USER':
        flatScale = userScale
    else:
        raise RuntimeError('%s : %s not implemented' % ("flatCorrection", scalingType))

    if not invert:
        maskedImage.scaledDivides(1.0/flatScale, flatMaskedImage)
    else:
        maskedImage.scaledMultiplies(1.0/flatScale, flatMaskedImage)


def illuminationCorrection(maskedImage, illumMaskedImage, illumScale, trimToFit=True):
    """Apply illumination correction in place.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.  The image is modified.
    illumMaskedImage : `lsst.afw.image.MaskedImage`
        Illumination correction image of the same size as ``maskedImage``.
    illumScale : scalar
        Scale factor for the illumination correction.
    trimToFit : `Bool`, optional
        If True, raw data is symmetrically trimmed to match
        calibration size.

    Raises
    ------
    RuntimeError
        Raised if ``maskedImage`` and ``illumMaskedImage`` do not have
        the same size.
    """
    if trimToFit:
        maskedImage = trimToMatchCalibBBox(maskedImage, illumMaskedImage)

    if maskedImage.getBBox(afwImage.LOCAL) != illumMaskedImage.getBBox(afwImage.LOCAL):
        raise RuntimeError("maskedImage bbox %s != illumMaskedImage bbox %s" %
                           (maskedImage.getBBox(afwImage.LOCAL), illumMaskedImage.getBBox(afwImage.LOCAL)))

    maskedImage.scaledDivides(1.0/illumScale, illumMaskedImage)


def overscanCorrection(ampMaskedImage, overscanImage, fitType='MEDIAN', order=1, collapseRej=3.0,
                       statControl=None, overscanIsInt=True):
    """Apply overscan correction in place.

    Parameters
    ----------
    ampMaskedImage : `lsst.afw.image.MaskedImage`
        Image of amplifier to correct; modified.
    overscanImage : `lsst.afw.image.Image` or `lsst.afw.image.MaskedImage`
        Image of overscan; modified.
    fitType : `str`
        Type of fit for overscan correction. May be one of:

        - ``MEAN``: use mean of overscan.
        - ``MEANCLIP``: use clipped mean of overscan.
        - ``MEDIAN``: use median of overscan.
        - ``POLY``: fit with ordinary polynomial.
        - ``CHEB``: fit with Chebyshev polynomial.
        - ``LEG``: fit with Legendre polynomial.
        - ``NATURAL_SPLINE``: fit with natural spline.
        - ``CUBIC_SPLINE``: fit with cubic spline.
        - ``AKIMA_SPLINE``: fit with Akima spline.

    order : `int`
        Polynomial order or number of spline knots; ignored unless
        ``fitType`` indicates a polynomial or spline.
    statControl : `lsst.afw.math.StatisticsControl`
        Statistics control object.  In particular, we pay attention to numSigmaClip
    overscanIsInt : `bool`
        Treat the overscan region as consisting of integers, even if it's been
        converted to float.  E.g. handle ties properly.

    Returns
    -------
    result : `lsst.pipe.base.Struct`
        Result struct with components:

        - ``imageFit``: Value(s) removed from image (scalar or
            `lsst.afw.image.Image`)
        - ``overscanFit``: Value(s) removed from overscan (scalar or
            `lsst.afw.image.Image`)
        - ``overscanImage``: Overscan corrected overscan region
            (`lsst.afw.image.Image`)
    Raises
    ------
    pexExcept.Exception
        Raised if ``fitType`` is not an allowed value.

    Notes
    -----
    The ``ampMaskedImage`` and ``overscanImage`` are modified, with the fit
    subtracted. Note that the ``overscanImage`` should not be a subimage of
    the ``ampMaskedImage``, to avoid being subtracted twice.

    Debug plots are available for the SPLINE fitTypes by setting the
    `debug.display` for `name` == "lsst.ip.isr.isrFunctions".  These
    plots show the scatter plot of the overscan data (collapsed along
    the perpendicular dimension) as a function of position on the CCD
    (normalized between +/-1).
    """
    ampImage = ampMaskedImage.getImage()
    if statControl is None:
        statControl = afwMath.StatisticsControl()

    numSigmaClip = statControl.getNumSigmaClip()

    if fitType in ('MEAN', 'MEANCLIP'):
        fitType = afwMath.stringToStatisticsProperty(fitType)
        offImage = afwMath.makeStatistics(overscanImage, fitType, statControl).getValue()
        overscanFit = offImage
    elif fitType in ('MEDIAN',):
        if overscanIsInt:
            # we need an image with integer pixels to handle ties properly
            if hasattr(overscanImage, "image"):
                imageI = overscanImage.image.convertI()
                overscanImageI = afwImage.MaskedImageI(imageI, overscanImage.mask, overscanImage.variance)
            else:
                overscanImageI = overscanImage.convertI()
        else:
            overscanImageI = overscanImage

        fitType = afwMath.stringToStatisticsProperty(fitType)
        offImage = afwMath.makeStatistics(overscanImageI, fitType, statControl).getValue()
        overscanFit = offImage

        if overscanIsInt:
            del overscanImageI
    elif fitType in ('POLY', 'CHEB', 'LEG', 'NATURAL_SPLINE', 'CUBIC_SPLINE', 'AKIMA_SPLINE'):
        if hasattr(overscanImage, "getImage"):
            biasArray = overscanImage.getImage().getArray()
            biasArray = numpy.ma.masked_where(overscanImage.getMask().getArray() & statControl.getAndMask(),
                                              biasArray)
        else:
            biasArray = overscanImage.getArray()
        # Fit along the long axis, so collapse along each short row and fit the resulting array
        shortInd = numpy.argmin(biasArray.shape)
        if shortInd == 0:
            # Convert to some 'standard' representation to make things easier
            biasArray = numpy.transpose(biasArray)

        # Do a single round of clipping to weed out CR hits and signal leaking into the overscan
        percentiles = numpy.percentile(biasArray, [25.0, 50.0, 75.0], axis=1)
        medianBiasArr = percentiles[1]
        stdevBiasArr = 0.74*(percentiles[2] - percentiles[0])  # robust stdev
        diff = numpy.abs(biasArray - medianBiasArr[:, numpy.newaxis])
        biasMaskedArr = numpy.ma.masked_where(diff > numSigmaClip*stdevBiasArr[:, numpy.newaxis], biasArray)
        collapsed = numpy.mean(biasMaskedArr, axis=1)
        if collapsed.mask.sum() > 0:
            collapsed.data[collapsed.mask] = numpy.mean(biasArray.data[collapsed.mask], axis=1)
        del biasArray, percentiles, stdevBiasArr, diff, biasMaskedArr

        if shortInd == 0:
            collapsed = numpy.transpose(collapsed)

        num = len(collapsed)
        indices = 2.0*numpy.arange(num)/float(num) - 1.0

        if fitType in ('POLY', 'CHEB', 'LEG'):
            # A numpy polynomial
            poly = numpy.polynomial
            fitter, evaler = {"POLY": (poly.polynomial.polyfit, poly.polynomial.polyval),
                              "CHEB": (poly.chebyshev.chebfit, poly.chebyshev.chebval),
                              "LEG": (poly.legendre.legfit, poly.legendre.legval),
                              }[fitType]

            coeffs = fitter(indices, collapsed, order)
            fitBiasArr = evaler(indices, coeffs)
        elif 'SPLINE' in fitType:
            # An afw interpolation
            numBins = order
            #
            # numpy.histogram needs a real array for the mask, but numpy.ma "optimises" the case
            # no-values-are-masked by replacing the mask array by a scalar, numpy.ma.nomask
            #
            # Issue DM-415
            #
            collapsedMask = collapsed.mask
            try:
                if collapsedMask == numpy.ma.nomask:
                    collapsedMask = numpy.array(len(collapsed)*[numpy.ma.nomask])
            except ValueError:      # If collapsedMask is an array the test fails [needs .all()]
                pass

            numPerBin, binEdges = numpy.histogram(indices, bins=numBins,
                                                  weights=1-collapsedMask.astype(int))
            # Binning is just a histogram, with weights equal to the values.
            # Use a similar trick to get the bin centers (this deals with different numbers per bin).
            with numpy.errstate(invalid="ignore"):  # suppress NAN warnings
                values = numpy.histogram(indices, bins=numBins,
                                         weights=collapsed.data*~collapsedMask)[0]/numPerBin
                binCenters = numpy.histogram(indices, bins=numBins,
                                             weights=indices*~collapsedMask)[0]/numPerBin
                interp = afwMath.makeInterpolate(binCenters.astype(float)[numPerBin > 0],
                                                 values.astype(float)[numPerBin > 0],
                                                 afwMath.stringToInterpStyle(fitType))
            fitBiasArr = numpy.array([interp.interpolate(i) for i in indices])

        import lsstDebug
        if lsstDebug.Info(__name__).display:
            import matplotlib.pyplot as plot
            figure = plot.figure(1)
            figure.clear()
            axes = figure.add_axes((0.1, 0.1, 0.8, 0.8))
            axes.plot(indices[~collapsedMask], collapsed[~collapsedMask], 'k+')
            if collapsedMask.sum() > 0:
                axes.plot(indices[collapsedMask], collapsed.data[collapsedMask], 'b+')
            axes.plot(indices, fitBiasArr, 'r-')
            plot.xlabel("centered/scaled position along overscan region")
            plot.ylabel("pixel value/fit value")
            figure.show()
            prompt = "Press Enter or c to continue [chp]... "
            while True:
                ans = input(prompt).lower()
                if ans in ("", "c",):
                    break
                if ans in ("p",):
                    import pdb
                    pdb.set_trace()
                elif ans in ("h", ):
                    print("h[elp] c[ontinue] p[db]")
            plot.close()

        offImage = ampImage.Factory(ampImage.getDimensions())
        offArray = offImage.getArray()
        overscanFit = afwImage.ImageF(overscanImage.getDimensions())
        overscanArray = overscanFit.getArray()
        if shortInd == 1:
            offArray[:, :] = fitBiasArr[:, numpy.newaxis]
            overscanArray[:, :] = fitBiasArr[:, numpy.newaxis]
        else:
            offArray[:, :] = fitBiasArr[numpy.newaxis, :]
            overscanArray[:, :] = fitBiasArr[numpy.newaxis, :]

        # We don't trust any extrapolation: mask those pixels as SUSPECT
        # This will occur when the top and or bottom edges of the overscan
        # contain saturated values. The values will be extrapolated from
        # the surrounding pixels, but we cannot entirely trust the value of
        # the extrapolation, and will mark the image mask plane to flag the
        # image as such.
        mask = ampMaskedImage.getMask()
        maskArray = mask.getArray() if shortInd == 1 else mask.getArray().transpose()
        suspect = mask.getPlaneBitMask("SUSPECT")
        try:
            if collapsed.mask == numpy.ma.nomask:
                # There is no mask, so the whole array is fine
                pass
        except ValueError:      # If collapsed.mask is an array the test fails [needs .all()]
            for low in range(num):
                if not collapsed.mask[low]:
                    break
            if low > 0:
                maskArray[:low, :] |= suspect
            for high in range(1, num):
                if not collapsed.mask[-high]:
                    break
            if high > 1:
                maskArray[-high:, :] |= suspect

    else:
        raise pexExcept.Exception('%s : %s an invalid overscan type' % ("overscanCorrection", fitType))
    ampImage -= offImage
    overscanImage -= overscanFit
    return Struct(imageFit=offImage, overscanFit=overscanFit, overscanImage=overscanImage)


def brighterFatterCorrection(exposure, kernel, maxIter, threshold, applyGain):
    """Apply brighter fatter correction in place for the image.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to have brighter-fatter correction applied.  Modified
        by this method.
    kernel : `numpy.ndarray`
        Brighter-fatter kernel to apply.
    maxIter : scalar
        Number of correction iterations to run.
    threshold : scalar
        Convergence threshold in terms of the sum of absolute
        deviations between an iteration and the previous one.
    applyGain : `Bool`
        If True, then the exposure values are scaled by the gain prior
        to correction.

    Returns
    -------
    diff : `float`
        Final difference between iterations achieved in correction.
    iteration : `int`
        Number of iterations used to calculate correction.

    Notes
    -----
    This correction takes a kernel that has been derived from flat
    field images to redistribute the charge.  The gradient of the
    kernel is the deflection field due to the accumulated charge.

    Given the original image I(x) and the kernel K(x) we can compute
    the corrected image Ic(x) using the following equation:

    Ic(x) = I(x) + 0.5*d/dx(I(x)*d/dx(int( dy*K(x-y)*I(y))))

    To evaluate the derivative term we expand it as follows:

    0.5 * ( d/dx(I(x))*d/dx(int(dy*K(x-y)*I(y))) + I(x)*d^2/dx^2(int(dy* K(x-y)*I(y))) )

    Because we use the measured counts instead of the incident counts
    we apply the correction iteratively to reconstruct the original
    counts and the correction.  We stop iterating when the summed
    difference between the current corrected image and the one from
    the previous iteration is below the threshold.  We do not require
    convergence because the number of iterations is too large a
    computational cost.  How we define the threshold still needs to be
    evaluated, the current default was shown to work reasonably well
    on a small set of images.  For more information on the method see
    DocuShare Document-19407.

    The edges as defined by the kernel are not corrected because they
    have spurious values due to the convolution.
    """
    image = exposure.getMaskedImage().getImage()

    # The image needs to be units of electrons/holes
    with gainContext(exposure, image, applyGain):

        kLx = numpy.shape(kernel)[0]
        kLy = numpy.shape(kernel)[1]
        kernelImage = afwImage.ImageD(kLx, kLy)
        kernelImage.getArray()[:, :] = kernel
        tempImage = image.clone()

        nanIndex = numpy.isnan(tempImage.getArray())
        tempImage.getArray()[nanIndex] = 0.

        outImage = afwImage.ImageF(image.getDimensions())
        corr = numpy.zeros_like(image.getArray())
        prev_image = numpy.zeros_like(image.getArray())
        convCntrl = afwMath.ConvolutionControl(False, True, 1)
        fixedKernel = afwMath.FixedKernel(kernelImage)

        # Define boundary by convolution region.  The region that the correction will be
        # calculated for is one fewer in each dimension because of the second derivative terms.
        # NOTE: these need to use integer math, as we're using start:end as numpy index ranges.
        startX = kLx//2
        endX = -kLx//2
        startY = kLy//2
        endY = -kLy//2

        for iteration in range(maxIter):

            afwMath.convolve(outImage, tempImage, fixedKernel, convCntrl)
            tmpArray = tempImage.getArray()
            outArray = outImage.getArray()

            with numpy.errstate(invalid="ignore", over="ignore"):
                # First derivative term
                gradTmp = numpy.gradient(tmpArray[startY:endY, startX:endX])
                gradOut = numpy.gradient(outArray[startY:endY, startX:endX])
                first = (gradTmp[0]*gradOut[0] + gradTmp[1]*gradOut[1])[1:-1, 1:-1]

                # Second derivative term
                diffOut20 = numpy.diff(outArray, 2, 0)[startY:endY, startX + 1:endX - 1]
                diffOut21 = numpy.diff(outArray, 2, 1)[startY + 1:endY - 1, startX:endX]
                second = tmpArray[startY + 1:endY - 1, startX + 1:endX - 1]*(diffOut20 + diffOut21)

                corr[startY + 1:endY - 1, startX + 1:endX - 1] = 0.5*(first + second)

                tmpArray[:, :] = image.getArray()[:, :]
                tmpArray[nanIndex] = 0.
                tmpArray[startY:endY, startX:endX] += corr[startY:endY, startX:endX]

            if iteration > 0:
                diff = numpy.sum(numpy.abs(prev_image - tmpArray))

                if diff < threshold:
                    break
                prev_image[:, :] = tmpArray[:, :]

        image.getArray()[startY + 1:endY - 1, startX + 1:endX - 1] += \
            corr[startY + 1:endY - 1, startX + 1:endX - 1]

    return diff, iteration


@contextmanager
def gainContext(exp, image, apply):
    """Context manager that applies and removes gain.

    Parameters
    ----------
    exp : `lsst.afw.image.Exposure`
        Exposure to apply/remove gain.
    image : `lsst.afw.image.Image`
        Image to apply/remove gain.
    apply : `Bool`
        If True, apply and remove the amplifier gain.

    Yields
    ------
    exp : `lsst.afw.image.Exposure`
        Exposure with the gain applied.
    """
    if apply:
        ccd = exp.getDetector()
        for amp in ccd:
            sim = image.Factory(image, amp.getBBox())
            sim *= amp.getGain()

    try:
        yield exp
    finally:
        if apply:
            ccd = exp.getDetector()
            for amp in ccd:
                sim = image.Factory(image, amp.getBBox())
                sim /= amp.getGain()


def attachTransmissionCurve(exposure, opticsTransmission=None, filterTransmission=None,
                            sensorTransmission=None, atmosphereTransmission=None):
    """Attach a TransmissionCurve to an Exposure, given separate curves for
    different components.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure object to modify by attaching the product of all given
        ``TransmissionCurves`` in post-assembly trimmed detector coordinates.
        Must have a valid ``Detector`` attached that matches the detector
        associated with sensorTransmission.
    opticsTransmission : `lsst.afw.image.TransmissionCurve`
        A ``TransmissionCurve`` that represents the throughput of the optics,
        to be evaluated in focal-plane coordinates.
    filterTransmission : `lsst.afw.image.TransmissionCurve`
        A ``TransmissionCurve`` that represents the throughput of the filter
        itself, to be evaluated in focal-plane coordinates.
    sensorTransmission : `lsst.afw.image.TransmissionCurve`
        A ``TransmissionCurve`` that represents the throughput of the sensor
        itself, to be evaluated in post-assembly trimmed detector coordinates.
    atmosphereTransmission : `lsst.afw.image.TransmissionCurve`
        A ``TransmissionCurve`` that represents the throughput of the
        atmosphere, assumed to be spatially constant.

    Returns
    -------
    combined : `lsst.afw.image.TransmissionCurve`
        The TransmissionCurve attached to the exposure.

    Notes
    -----
    All ``TransmissionCurve`` arguments are optional; if none are provided, the
    attached ``TransmissionCurve`` will have unit transmission everywhere.
    """
    combined = afwImage.TransmissionCurve.makeIdentity()
    if atmosphereTransmission is not None:
        combined *= atmosphereTransmission
    if opticsTransmission is not None:
        combined *= opticsTransmission
    if filterTransmission is not None:
        combined *= filterTransmission
    detector = exposure.getDetector()
    fpToPix = detector.getTransform(fromSys=camGeom.FOCAL_PLANE,
                                    toSys=camGeom.PIXELS)
    combined = combined.transformedBy(fpToPix)
    if sensorTransmission is not None:
        combined *= sensorTransmission
    exposure.getInfo().setTransmissionCurve(combined)
    return combined


@deprecated(reason="Camera geometry-based SkyWcs are now set when reading raws. To be removed after v19.",
            category=FutureWarning)
def addDistortionModel(exposure, camera):
    """!Update the WCS in exposure with a distortion model based on camera
    geometry.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to process.  Must contain a Detector and WCS.  The
        exposure is modified.
    camera : `lsst.afw.cameraGeom.Camera`
        Camera geometry.

    Raises
    ------
    RuntimeError
        Raised if ``exposure`` is lacking a Detector or WCS, or if
        ``camera`` is None.
    Notes
    -----
    Add a model for optical distortion based on geometry found in ``camera``
    and the ``exposure``'s detector. The raw input exposure is assumed
    have a TAN WCS that has no compensation for optical distortion.
    Two other possibilities are:
    - The raw input exposure already has a model for optical distortion,
    as is the case for raw DECam data.
    In that case you should set config.doAddDistortionModel False.
    - The raw input exposure has a model for distortion, but it has known
    deficiencies severe enough to be worth fixing (e.g. because they
    cause problems for fitting a better WCS). In that case you should
    override this method with a version suitable for your raw data.

    """
    wcs = exposure.getWcs()
    if wcs is None:
        raise RuntimeError("exposure has no WCS")
    if camera is None:
        raise RuntimeError("camera is None")
    detector = exposure.getDetector()
    if detector is None:
        raise RuntimeError("exposure has no Detector")
    pixelToFocalPlane = detector.getTransform(camGeom.PIXELS, camGeom.FOCAL_PLANE)
    focalPlaneToFieldAngle = camera.getTransformMap().getTransform(camGeom.FOCAL_PLANE,
                                                                   camGeom.FIELD_ANGLE)
    distortedWcs = makeDistortedTanWcs(wcs, pixelToFocalPlane, focalPlaneToFieldAngle)
    exposure.setWcs(distortedWcs)


def applyGains(exposure, normalizeGains=False):
    """Scale an exposure by the amplifier gains.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to process.  The image is modified.
    normalizeGains : `Bool`, optional
        If True, then amplifiers are scaled to force the median of
        each amplifier to equal the median of those medians.
    """
    ccd = exposure.getDetector()
    ccdImage = exposure.getMaskedImage()

    medians = []
    for amp in ccd:
        sim = ccdImage.Factory(ccdImage, amp.getBBox())
        sim *= amp.getGain()

        if normalizeGains:
            medians.append(numpy.median(sim.getImage().getArray()))

    if normalizeGains:
        median = numpy.median(numpy.array(medians))
        for index, amp in enumerate(ccd):
            sim = ccdImage.Factory(ccdImage, amp.getBBox())
            if medians[index] != 0.0:
                sim *= median/medians[index]


def widenSaturationTrails(mask):
    """Grow the saturation trails by an amount dependent on the width of the trail.

    Parameters
    ----------
    mask : `lsst.afw.image.Mask`
        Mask which will have the saturated areas grown.
    """

    extraGrowDict = {}
    for i in range(1, 6):
        extraGrowDict[i] = 0
    for i in range(6, 8):
        extraGrowDict[i] = 1
    for i in range(8, 10):
        extraGrowDict[i] = 3
    extraGrowMax = 4

    if extraGrowMax <= 0:
        return

    saturatedBit = mask.getPlaneBitMask("SAT")

    xmin, ymin = mask.getBBox().getMin()
    width = mask.getWidth()

    thresh = afwDetection.Threshold(saturatedBit, afwDetection.Threshold.BITMASK)
    fpList = afwDetection.FootprintSet(mask, thresh).getFootprints()

    for fp in fpList:
        for s in fp.getSpans():
            x0, x1 = s.getX0(), s.getX1()

            extraGrow = extraGrowDict.get(x1 - x0 + 1, extraGrowMax)
            if extraGrow > 0:
                y = s.getY() - ymin
                x0 -= xmin + extraGrow
                x1 -= xmin - extraGrow

                if x0 < 0:
                    x0 = 0
                if x1 >= width - 1:
                    x1 = width - 1

                mask.array[y, x0:x1+1] |= saturatedBit


def setBadRegions(exposure, badStatistic="MEDIAN"):
    """Set all BAD areas of the chip to the average of the rest of the exposure

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to mask.  The exposure mask is modified.
    badStatistic : `str`, optional
        Statistic to use to generate the replacement value from the
        image data.  Allowed values are 'MEDIAN' or 'MEANCLIP'.

    Returns
    -------
    badPixelCount : scalar
        Number of bad pixels masked.
    badPixelValue : scalar
        Value substituted for bad pixels.

    Raises
    ------
    RuntimeError
        Raised if `badStatistic` is not an allowed value.
    """
    if badStatistic == "MEDIAN":
        statistic = afwMath.MEDIAN
    elif badStatistic == "MEANCLIP":
        statistic = afwMath.MEANCLIP
    else:
        raise RuntimeError("Impossible method %s of bad region correction" % badStatistic)

    mi = exposure.getMaskedImage()
    mask = mi.getMask()
    BAD = mask.getPlaneBitMask("BAD")
    INTRP = mask.getPlaneBitMask("INTRP")

    sctrl = afwMath.StatisticsControl()
    sctrl.setAndMask(BAD)
    value = afwMath.makeStatistics(mi, statistic, sctrl).getValue()

    maskArray = mask.getArray()
    imageArray = mi.getImage().getArray()
    badPixels = numpy.logical_and((maskArray & BAD) > 0, (maskArray & INTRP) == 0)
    imageArray[:] = numpy.where(badPixels, value, imageArray)

    return badPixels.sum(), value
