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

__all__ = [
    "applyGains",
    "attachTransmissionCurve",
    "biasCorrection",
    "brighterFatterCorrection",
    "checkFilter",
    "countMaskedPixels",
    "createPsf",
    "darkCorrection",
    "flatCorrection",
    "fluxConservingBrighterFatterCorrection",
    "gainContext",
    "getPhysicalFilter",
    "growMasks",
    "illuminationCorrection",
    "interpolateDefectList",
    "interpolateFromMask",
    "makeThresholdMask",
    "saturationCorrection",
    "setBadRegions",
    "transferFlux",
    "transposeMaskedImage",
    "trimToMatchCalibBBox",
    "updateVariance",
    "widenSaturationTrails",
]

import math
import numpy

import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.afw.cameraGeom as camGeom

from lsst.meas.algorithms.detection import SourceDetectionTask

from contextlib import contextmanager

from .defects import Defects


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
    fwhm : `float`
        FWHM of double Gaussian smoothing kernel.
    fallbackValue : scalar, optional
        Fallback value if an interpolated value cannot be determined.
        If None, then the clipped mean of the image is used.

    Notes
    -----
    The ``fwhm`` parameter is used to create a PSF, but the underlying
    interpolation code (`lsst.meas.algorithms.interpolateOverDefects`) does
    not currently make use of this information.
    """
    psf = createPsf(fwhm)
    if fallbackValue is None:
        fallbackValue = afwMath.makeStatistics(maskedImage.getImage(), afwMath.MEANCLIP).getValue()
    if 'INTRP' not in maskedImage.getMask().getMaskPlaneDict():
        maskedImage.getMask().addMaskPlane('INTRP')
    measAlg.interpolateOverDefects(maskedImage, psf, defectList, fallbackValue, True)
    return maskedImage


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

    return Defects.fromFootprintList(fpList)


def growMasks(mask, radius=0, maskNameList=['BAD'], maskValue="BAD"):
    """Grow a mask by an amount and add to the requested plane.

    Parameters
    ----------
    mask : `lsst.afw.image.Mask`
        Mask image to process.
    radius : scalar
        Amount to grow the mask.
    maskNameList : `str` or `list` [`str`]
        Mask names that should be grown.
    maskValue : `str`
        Mask plane to assign the newly masked pixels to.
    """
    if radius > 0:
        thresh = afwDetection.Threshold(mask.getPlaneBitMask(maskNameList), afwDetection.Threshold.BITMASK)
        fpSet = afwDetection.FootprintSet(mask, thresh)
        fpSet = afwDetection.FootprintSet(fpSet, rGrow=radius, isotropic=False)
        fpSet.setMask(mask, maskValue)


def interpolateFromMask(maskedImage, fwhm, growSaturatedFootprints=1,
                        maskNameList=['SAT'], fallbackValue=None):
    """Interpolate over defects identified by a particular set of mask planes.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.
    fwhm : `float`
        FWHM of double Gaussian smoothing kernel.
    growSaturatedFootprints : scalar, optional
        Number of pixels to grow footprints for saturated pixels.
    maskNameList : `List` of `str`, optional
        Mask plane name.
    fallbackValue : scalar, optional
        Value of last resort for interpolation.

    Notes
    -----
    The ``fwhm`` parameter is used to create a PSF, but the underlying
    interpolation code (`lsst.meas.algorithms.interpolateOverDefects`) does
    not currently make use of this information.
    """
    mask = maskedImage.getMask()

    if growSaturatedFootprints > 0 and "SAT" in maskNameList:
        # If we are interpolating over an area larger than the original masked
        # region, we need to expand the original mask bit to the full area to
        # explain why we interpolated there.
        growMasks(mask, radius=growSaturatedFootprints, maskNameList=['SAT'], maskValue="SAT")

    thresh = afwDetection.Threshold(mask.getPlaneBitMask(maskNameList), afwDetection.Threshold.BITMASK)
    fpSet = afwDetection.FootprintSet(mask, thresh)
    defectList = Defects.fromFootprintList(fpSet.getFootprints())

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
    fwhm : `float`
        FWHM of double Gaussian smoothing kernel.
    growFootprints : scalar, optional
        Number of pixels to grow footprints of detected regions.
    interpolate : Bool, optional
        If True, saturated pixels are interpolated over.
    maskName : str, optional
        Mask plane name.
    fallbackValue : scalar, optional
        Value of last resort for interpolation.

    Notes
    -----
    The ``fwhm`` parameter is used to create a PSF, but the underlying
    interpolation code (`lsst.meas.algorithms.interpolateOverDefects`) does
    not currently make use of this information.
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
        ``rawMaskedImage`` trimmed to the appropriate size.

    Raises
    ------
    RuntimeError
       Raised if ``rawMaskedImage`` cannot be symmetrically trimmed to
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
        Scale to use if ``scalingType='USER'``.
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
    # Ideally the flats are normalized by the calibration product pipeline,
    # but this allows some flexibility in the case that the flat is created by
    # some other mechanism.
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


def brighterFatterCorrection(exposure, kernel, maxIter, threshold, applyGain, gains=None):
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
    gains : `dict` [`str`, `float`]
        A dictionary, keyed by amplifier name, of the gains to use.
        If gains is None, the nominal gains in the amplifier object are used.

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

    0.5 * ( d/dx(I(x))*d/dx(int(dy*K(x-y)*I(y)))
        + I(x)*d^2/dx^2(int(dy* K(x-y)*I(y))) )

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
    with gainContext(exposure, image, applyGain, gains):

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

        # Define boundary by convolution region.  The region that the
        # correction will be calculated for is one fewer in each dimension
        # because of the second derivative terms.
        # NOTE: these need to use integer math, as we're using start:end as
        # numpy index ranges.
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


def transferFlux(cFunc, fStep, correctionMode=True):
    """Take the input convolved deflection potential and the flux array
    to compute and apply the flux transfer into the correction array.

    Parameters
    ----------
    cFunc: `numpy.array`
        Deflection potential, being the convolution of the flux F with the
        kernel K.
    fStep: `numpy.array`
        The array of flux values which act as the source of the flux transfer.
    correctionMode: `bool`
        Defines if applying correction (True) or generating sims (False).

    Returns
    -------
    corr:
        BFE correction array
    """

    if cFunc.shape != fStep.shape:
        raise RuntimeError(f'transferFlux: array shapes do not match: {cFunc.shape}, {fStep.shape}')

    # set the sign of the correction and set its value for the
    # time averaged solution
    if correctionMode:
        # negative sign if applying BFE correction
        factor = -0.5
    else:
        # positive sign if generating BFE simulations
        factor = 0.5

    # initialise the BFE correction image to zero
    corr = numpy.zeros_like(cFunc)

    # Generate a 2D mesh of x,y coordinates
    yDim, xDim = cFunc.shape
    y = numpy.arange(yDim, dtype=int)
    x = numpy.arange(xDim, dtype=int)
    xc, yc = numpy.meshgrid(x, y)

    # process each axis in turn
    for ax in [0, 1]:

        # gradient of phi on right/upper edge of pixel
        diff = numpy.diff(cFunc, axis=ax)

        # expand array back to full size with zero gradient at the end
        gx = numpy.zeros_like(cFunc)
        yDiff, xDiff = diff.shape
        gx[:yDiff, :xDiff] += diff

        # select pixels with either positive gradients on the right edge,
        # flux flowing to the right/up
        # or negative gradients, flux flowing to the left/down
        for i, sel in enumerate([gx > 0, gx < 0]):
            xSelPixels = xc[sel]
            ySelPixels = yc[sel]
            # and add the flux into the pixel to the right or top
            # depending on which axis we are handling
            if ax == 0:
                xPix = xSelPixels
                yPix = ySelPixels+1
            else:
                xPix = xSelPixels+1
                yPix = ySelPixels
            # define flux as the either current pixel value or pixel
            # above/right
            # depending on whether positive or negative gradient
            if i == 0:
                # positive gradients, flux flowing to higher coordinate values
                flux = factor * fStep[sel]*gx[sel]
            else:
                # negative gradients, flux flowing to lower coordinate values
                flux = factor * fStep[yPix, xPix]*gx[sel]
            # change the fluxes of the donor and receiving pixels
            # such that flux is conserved
            corr[sel] -= flux
            corr[yPix, xPix] += flux

    # return correction array
    return corr


def fluxConservingBrighterFatterCorrection(exposure, kernel, maxIter, threshold, applyGain,
                                           gains=None, correctionMode=True):
    """Apply brighter fatter correction in place for the image.

    This version presents a modified version of the algorithm
    found in ``lsst.ip.isr.isrFunctions.brighterFatterCorrection``
    which conserves the image flux, resulting in improved
    correction of the cores of stars. The convolution has also been
    modified to mitigate edge effects.

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
    gains : `dict` [`str`, `float`]
        A dictionary, keyed by amplifier name, of the gains to use.
        If gains is None, the nominal gains in the amplifier object are used.
    correctionMode : `Bool`
        If True (default) the function applies correction for BFE.  If False,
        the code can instead be used to generate a simulation of BFE (sign
        change in the direction of the effect)

    Returns
    -------
    diff : `float`
        Final difference between iterations achieved in correction.
    iteration : `int`
        Number of iterations used to calculate correction.

    Notes
    -----
    Modified version of ``lsst.ip.isr.isrFunctions.brighterFatterCorrection``.

    This correction takes a kernel that has been derived from flat
    field images to redistribute the charge.  The gradient of the
    kernel is the deflection field due to the accumulated charge.

    Given the original image I(x) and the kernel K(x) we can compute
    the corrected image Ic(x) using the following equation:

    Ic(x) = I(x) + 0.5*d/dx(I(x)*d/dx(int( dy*K(x-y)*I(y))))

    Improved algorithm at this step applies the divergence theorem to
    obtain a pixelised correction.

    Because we use the measured counts instead of the incident counts
    we apply the correction iteratively to reconstruct the original
    counts and the correction.  We stop iterating when the summed
    difference between the current corrected image and the one from
    the previous iteration is below the threshold.  We do not require
    convergence because the number of iterations is too large a
    computational cost.  How we define the threshold still needs to be
    evaluated, the current default was shown to work reasonably well
    on a small set of images.

    Edges are handled in the convolution by padding.  This is still not
    a physical model for the edge, but avoids discontinuity in the correction.

    Author of modified version: Lance.Miller@physics.ox.ac.uk
    (see DM-38555).
    """
    image = exposure.getMaskedImage().getImage()

    # The image needs to be units of electrons/holes
    with gainContext(exposure, image, applyGain, gains):

        # get kernel and its shape
        kLy, kLx = kernel.shape
        kernelImage = afwImage.ImageD(kLx, kLy)
        kernelImage.getArray()[:, :] = kernel
        tempImage = image.clone()

        nanIndex = numpy.isnan(tempImage.getArray())
        tempImage.getArray()[nanIndex] = 0.

        outImage = afwImage.ImageF(image.getDimensions())
        corr = numpy.zeros_like(image.getArray())
        prevImage = numpy.zeros_like(image.getArray())
        convCntrl = afwMath.ConvolutionControl(False, False, 1)
        fixedKernel = afwMath.FixedKernel(kernelImage)

        # set the padding amount
        # ensure we pad by an even amount larger than the kernel
        kLy = 2 * ((1+kLy)//2)
        kLx = 2 * ((1+kLx)//2)

        # The deflection potential only depends on the gradient of
        # the convolution, so we can subtract the mean, which then
        # allows us to pad the image with zeros and avoid wrap-around effects
        # (although still not handling the image edges with a physical model)
        # This wouldn't be great if there were a strong image gradient.
        imYdimension, imXdimension = tempImage.array.shape
        imean = numpy.mean(tempImage.getArray()[~nanIndex])
        # subtract mean from image
        tempImage -= imean
        tempImage.array[nanIndex] = 0.
        padArray = numpy.pad(tempImage.getArray(), ((0, kLy), (0, kLx)))
        outImage = afwImage.ImageF(numpy.pad(outImage.getArray(), ((0, kLy), (0, kLx))))
        # Convert array to afw image so afwMath.convolve works
        padImage = afwImage.ImageF(padArray.shape[1], padArray.shape[0])
        padImage.array[:] = padArray

        for iteration in range(maxIter):

            # create deflection potential, convolution of flux with kernel
            # using padded counts array
            afwMath.convolve(outImage, padImage, fixedKernel, convCntrl)
            tmpArray = tempImage.getArray()
            outArray = outImage.getArray()

            # trim convolution output back to original shape
            outArray = outArray[:imYdimension, :imXdimension]

            # generate the correction array, with correctionMode set as input
            corr[...] = transferFlux(outArray, tmpArray, correctionMode=correctionMode)

            # update the arrays for the next iteration
            tmpArray[:, :] = image.getArray()[:, :]
            tmpArray += corr
            tmpArray[nanIndex] = 0.
            # update padded array
            # subtract mean
            tmpArray -= imean
            tempImage.array[nanIndex] = 0.
            padArray = numpy.pad(tempImage.getArray(), ((0, kLy), (0, kLx)))

            if iteration > 0:
                diff = numpy.sum(numpy.abs(prevImage - tmpArray))

                if diff < threshold:
                    break
                prevImage[:, :] = tmpArray[:, :]

        image.getArray()[:] += corr[:]

    return diff, iteration


@contextmanager
def gainContext(exp, image, apply, gains=None):
    """Context manager that applies and removes gain.

    Parameters
    ----------
    exp : `lsst.afw.image.Exposure`
        Exposure to apply/remove gain.
    image : `lsst.afw.image.Image`
        Image to apply/remove gain.
    apply : `Bool`
        If True, apply and remove the amplifier gain.
    gains : `dict` [`str`, `float`]
        A dictionary, keyed by amplifier name, of the gains to use.
        If gains is None, the nominal gains in the amplifier object are used.

    Yields
    ------
    exp : `lsst.afw.image.Exposure`
        Exposure with the gain applied.
    """
    # check we have all of them if provided because mixing and matching would
    # be a real mess
    if gains and apply is True:
        ampNames = [amp.getName() for amp in exp.getDetector()]
        for ampName in ampNames:
            if ampName not in gains.keys():
                raise RuntimeError(f"Gains provided to gain context, but no entry found for amp {ampName}")

    if apply:
        ccd = exp.getDetector()
        for amp in ccd:
            sim = image.Factory(image, amp.getBBox())
            if gains:
                gain = gains[amp.getName()]
            else:
                gain = amp.getGain()
            sim *= gain

    try:
        yield exp
    finally:
        if apply:
            ccd = exp.getDetector()
            for amp in ccd:
                sim = image.Factory(image, amp.getBBox())
                if gains:
                    gain = gains[amp.getName()]
                else:
                    gain = amp.getGain()
                sim /= gain


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


def applyGains(exposure, normalizeGains=False, ptcGains=None):
    """Scale an exposure by the amplifier gains.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to process.  The image is modified.
    normalizeGains : `Bool`, optional
        If True, then amplifiers are scaled to force the median of
        each amplifier to equal the median of those medians.
    ptcGains : `dict`[`str`], optional
        Dictionary keyed by amp name containing the PTC gains.
    """
    ccd = exposure.getDetector()
    ccdImage = exposure.getMaskedImage()

    medians = []
    for amp in ccd:
        sim = ccdImage.Factory(ccdImage, amp.getBBox())
        if ptcGains:
            sim *= ptcGains[amp.getName()]
        else:
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
    """Grow the saturation trails by an amount dependent on the width of the
    trail.

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


def checkFilter(exposure, filterList, log):
    """Check to see if an exposure is in a filter specified by a list.

    The goal of this is to provide a unified filter checking interface
    for all filter dependent stages.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to examine.
    filterList : `list` [`str`]
        List of physical_filter names to check.
    log : `logging.Logger`
        Logger to handle messages.

    Returns
    -------
    result : `bool`
        True if the exposure's filter is contained in the list.
    """
    if len(filterList) == 0:
        return False
    thisFilter = exposure.getFilter()
    if thisFilter is None:
        log.warning("No FilterLabel attached to this exposure!")
        return False

    thisPhysicalFilter = getPhysicalFilter(thisFilter, log)
    if thisPhysicalFilter in filterList:
        return True
    elif thisFilter.bandLabel in filterList:
        if log:
            log.warning("Physical filter (%s) should be used instead of band %s for filter configurations"
                        " (%s)", thisPhysicalFilter, thisFilter.bandLabel, filterList)
        return True
    else:
        return False


def getPhysicalFilter(filterLabel, log):
    """Get the physical filter label associated with the given filterLabel.

    If ``filterLabel`` is `None` or there is no physicalLabel attribute
    associated with the given ``filterLabel``, the returned label will be
    "Unknown".

    Parameters
    ----------
    filterLabel : `lsst.afw.image.FilterLabel`
        The `lsst.afw.image.FilterLabel` object from which to derive the
        physical filter label.
    log : `logging.Logger`
        Logger to handle messages.

    Returns
    -------
    physicalFilter : `str`
        The value returned by the physicalLabel attribute of ``filterLabel`` if
        it exists, otherwise set to \"Unknown\".
    """
    if filterLabel is None:
        physicalFilter = "Unknown"
        log.warning("filterLabel is None.  Setting physicalFilter to \"Unknown\".")
    else:
        try:
            physicalFilter = filterLabel.physicalLabel
        except RuntimeError:
            log.warning("filterLabel has no physicalLabel attribute.  Setting physicalFilter to \"Unknown\".")
            physicalFilter = "Unknown"
    return physicalFilter


def countMaskedPixels(maskedIm, maskPlane):
    """Count the number of pixels in a given mask plane.

    Parameters
    ----------
    maskedIm : `~lsst.afw.image.MaskedImage`
        Masked image to examine.
    maskPlane : `str`
        Name of the mask plane to examine.

    Returns
    -------
    nPix : `int`
        Number of pixels in the requested mask plane.
    """
    maskBit = maskedIm.mask.getPlaneBitMask(maskPlane)
    nPix = numpy.where(numpy.bitwise_and(maskedIm.mask.array, maskBit))[0].flatten().size
    return nPix
