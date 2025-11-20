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
    "compareCameraKeywords",
    "countMaskedPixels",
    "createPsf",
    "darkCorrection",
    "flatCorrection",
    "fluxConservingBrighterFatterCorrection",
    "gainContext",
    "getPhysicalFilter",
    "growMasks",
    "maskE2VEdgeBleed",
    "maskITLEdgeBleed",
    "maskITLSatSag",
    "maskITLDip",
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
    "getExposureGains",
    "getExposureReadNoises",
]

import logging
import math
import numpy
import scipy.signal
import pyfftw

import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.afw.cameraGeom as camGeom

from lsst.afw.geom import SpanSet, Stencil
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


def interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=None,
                          maskNameList=None, useLegacyInterp=True):
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
    maskNameList : `list [string]`
        List of the defects to interpolate over (used for GP interpolator).
    useLegacyInterp : `bool`
        Use the legacy interpolation (polynomial interpolation) if True. Use
        Gaussian Process interpolation if False.

    Notes
    -----
    The ``fwhm`` parameter is used to create a PSF, but the underlying
    interpolation code (`lsst.meas.algorithms.interpolateOverDefects`) does
    not currently make use of this information in legacy Interpolation, but use
    if for the Gaussian Process as an estimation of the correlation lenght.
    """
    psf = createPsf(fwhm)
    if fallbackValue is None:
        fallbackValue = afwMath.makeStatistics(maskedImage.getImage(), afwMath.MEANCLIP).getValue()
    if 'INTRP' not in maskedImage.getMask().getMaskPlaneDict():
        maskedImage.getMask().addMaskPlane('INTRP')

    # Hardcoded fwhm value. PSF estimated latter in step1,
    # not in ISR.
    if useLegacyInterp:
        kwargs = {}
        fwhm = fwhm
    else:
        # tested on a dozens of images and looks a good set of
        # hyperparameters, but cannot guarrenty this is optimal,
        # need further testing.
        kwargs = {"bin_spacing": 20,
                  "threshold_dynamic_binning": 2000,
                  "threshold_subdivide": 20000}
        fwhm = 15

    measAlg.interpolateOverDefects(maskedImage, psf, defectList,
                                   fallbackValue=fallbackValue,
                                   useFallbackValueAtEdge=True,
                                   fwhm=fwhm,
                                   useLegacyInterp=useLegacyInterp,
                                   maskNameList=maskNameList, **kwargs)
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
        spans = SpanSet.fromMask(mask, mask.getPlaneBitMask(maskNameList))
        # Use MANHATTAN for equivalence with 'isotropic=False` footprint grows,
        # but CIRCLE is probably better and might be just as fast.
        spans = spans.dilated(radius, Stencil.MANHATTAN)
        spans = spans.clippedTo(mask.getBBox())
        spans.setMask(mask, mask.getPlaneBitMask(maskValue))


def maskE2VEdgeBleed(exposure, e2vEdgeBleedSatMinArea=10000,
                     e2vEdgeBleedSatMaxArea=100000,
                     e2vEdgeBleedYMax=350,
                     saturatedMaskName="SAT", log=None):
    """Mask edge bleeds in E2V detectors.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to apply masking to.
    e2vEdgeBleedSatMinArea : `int`, optional
        Minimum limit of saturated cores footprint area.
    e2vEdgeBleedSatMaxArea : `int`, optional
        Maximum limit of saturated cores footprint area.
    e2vEdgeBleedYMax: `float`, optional
        Height of edge bleed masking.
    saturatedMaskName : `str`, optional
        Mask name for saturation.
    log : `logging.Logger`, optional
        Logger to handle messages.
    """

    log = log if log else logging.getLogger(__name__)

    maskedImage = exposure.maskedImage
    saturatedBit = maskedImage.mask.getPlaneBitMask(saturatedMaskName)

    thresh = afwDetection.Threshold(saturatedBit, afwDetection.Threshold.BITMASK)

    fpList = afwDetection.FootprintSet(exposure.mask, thresh).getFootprints()

    satAreas = numpy.asarray([fp.getArea() for fp in fpList])
    largeAreas, = numpy.where((satAreas >= e2vEdgeBleedSatMinArea)
                              & (satAreas < e2vEdgeBleedSatMaxArea))
    for largeAreasIndex in largeAreas:
        fpCore = fpList[largeAreasIndex]
        xCore, yCore = fpCore.getCentroid()
        xCore = int(xCore)
        yCore = int(yCore)

        for amp in exposure.getDetector():
            if amp.getBBox().contains(xCore, yCore):
                ampName = amp.getName()
                if ampName[:2] == 'C0':
                    # Check that the footprint reaches the bottom of the
                    # amplifier.
                    if fpCore.getBBox().getMinY() == 0:
                        # This is a large saturation footprint that hits the
                        # edge, and is thus classified as an edge bleed.

                        # TODO DM-50587: Optimize number of rows to mask by
                        # looking at the median signal level as a function of
                        # row number on the right side of the saturation trail.

                        log.info("Found E2V edge bleed in amp %s, column %d.", ampName, xCore)
                        maskedImage.mask[amp.getBBox()].array[:e2vEdgeBleedYMax, :] |= saturatedBit


def maskITLEdgeBleed(ccdExposure, badAmpDict,
                     fpCore, itlEdgeBleedSatMinArea=10000,
                     itlEdgeBleedSatMaxArea=100000,
                     itlEdgeBleedThreshold=5000.,
                     itlEdgeBleedModelConstant=0.02,
                     saturatedMaskName="SAT", log=None):
    """Mask edge bleeds in ITL detectors.

    Parameters
    ----------
    ccdExposure : `lsst.afw.image.Exposure`
        Exposure to apply masking to.
    badAmpDict : `dict` [`str`, `bool`]
        Dictionary of amplifiers, keyed by name, value is True if
        amplifier is fully masked.
    fpCore : `lsst.afw.detection._detection.Footprint`
        Footprint of saturated core.
    itlEdgeBleedThreshold : `float`, optional
        Threshold above median sky background for edge bleed detection
        (electron units).
    itlEdgeBleedModelConstant : `float`, optional
        Constant in the decaying exponential in the edge bleed masking.
    saturatedMaskName : `str`, optional
        Mask name for saturation.
    log : `logging.Logger`, optional
        Logger to handle messages.
    """

    log = log if log else logging.getLogger(__name__)

    # Get median of amplifier saturation level
    satLevel = numpy.nanmedian([ccdExposure.metadata[f"LSST ISR SATURATION LEVEL {amp.getName()}"]
                                for amp in ccdExposure.getDetector() if not badAmpDict[amp.getName()]])

    # 1. we check if there are several cores in the footprint:
    # Get centroid of saturated core
    xCore, yCore = fpCore.getCentroid()
    # Turn the Y detector coordinate into Y footprint coordinate
    yCoreFP = int(yCore) - fpCore.getBBox().getMinY()
    # Now test if there is one or more cores by checking if the slice at the
    # center is full of saturated pixels or has several segments of saturated
    # columns (i.e. several cores with trails)
    checkCoreNbRow = fpCore.getSpans().asArray()[yCoreFP, :]
    nbCore = 0
    indexSwitchTrue = []
    indexSwitchFalse = []
    if checkCoreNbRow[0]:
        # If the slice starts with saturated pixels
        inSatSegment = True
        nbCore = 1
        indexSwitchTrue.append(0)
    else:
        # If the slice starts with non saturated pixels
        inSatSegment = False

    for i, value in enumerate(checkCoreNbRow):
        if value:
            if not inSatSegment:
                indexSwitchTrue.append(i)
                # nbCore is the number of detected cores.
                nbCore += 1
                inSatSegment = True
        elif inSatSegment:
            indexSwitchFalse.append(i)
            inSatSegment = False

    # 1. we look for edge bleed in saturated cores in the footprint
    if nbCore == 2:
        # we now estimate the x coordinates of the edges of the subfootprint
        # for each core
        xEdgesCores = [0]
        xEdgesCores.append(int((indexSwitchTrue[1] + indexSwitchFalse[0])/2))
        xEdgesCores.append(fpCore.getSpans().asArray().shape[1])
        # Get the X and Y footprint coordinates of the cores
        for i in range(nbCore):
            subfp = fpCore.getSpans().asArray()[:, xEdgesCores[i]:xEdgesCores[i+1]]
            xCoreFP = int(xEdgesCores[i] + numpy.argmax(numpy.sum(subfp, axis=0)))
            # turn into X coordinate in detector space
            xCore = xCoreFP + fpCore.getBBox().getMinX()
            # get Y footprint coordinate of the core
            # by trimming the edges where edge bleeds are potentially dominant
            if subfp.shape[0] <= 200:
                yCoreFP = int(numpy.argmax(numpy.sum(subfp, axis=1)))
            else:
                yCoreFP = int(numpy.argmax(numpy.sum(subfp[100:-100, :],
                                                     axis=1)))
                yCoreFP = 100+yCoreFP

            # Estimate the width of the saturated core
            widthSat = numpy.sum(subfp[int(yCoreFP), :])

            subfpArea = numpy.sum(subfp)
            if subfpArea > itlEdgeBleedSatMinArea and subfpArea < itlEdgeBleedSatMaxArea:
                _applyMaskITLEdgeBleed(ccdExposure, xCore,
                                       satLevel, widthSat,
                                       itlEdgeBleedThreshold,
                                       itlEdgeBleedModelConstant,
                                       saturatedMaskName, log)
    elif nbCore > 2:
        # TODO DM-49736: support N cores in saturated footprint
        log.warning(
            "Too many (%d) cores in saturated footprint to mask edge bleeds.",
            nbCore,
        )
    else:
        # Get centroid of saturated core
        xCore, yCore = fpCore.getCentroid()
        # Turn the Y detector coordinate into Y footprint coordinate
        yCoreFP = yCore - fpCore.getBBox().getMinY()
        # Get the number of saturated columns around the centroid
        widthSat = numpy.sum(fpCore.getSpans().asArray()[int(yCoreFP), :])
        _applyMaskITLEdgeBleed(ccdExposure, xCore,
                               satLevel, widthSat, itlEdgeBleedThreshold,
                               itlEdgeBleedModelConstant, saturatedMaskName, log)


def _applyMaskITLEdgeBleed(ccdExposure, xCore,
                           satLevel, widthSat,
                           itlEdgeBleedThreshold=5000.,
                           itlEdgeBleedModelConstant=0.03,
                           saturatedMaskName="SAT", log=None):
    """Apply ITL edge bleed masking model.

    Parameters
    ----------
    ccdExposure : `lsst.afw.image.Exposure`
        Exposure to apply masking to.
    xCore: `int`
        X coordinate of the saturated core.
    satLevel: `float`
        Minimum saturation level of the detector.
    widthSat: `float`
        Width of the saturated core.
    itlEdgeBleedThreshold : `float`, optional
        Threshold above median sky background for edge bleed detection
        (electron units).
    itlEdgeBleedModelConstant : `float`, optional
        Constant in the decaying exponential in the edge bleed masking.
    saturatedMaskName : `str`, optional
        Mask name for saturation.
    log : `logging.Logger`, optional
        Logger to handle messages.
    """
    log = log if log else logging.getLogger(__name__)

    maskedImage = ccdExposure.maskedImage
    xmax = maskedImage.image.array.shape[1]
    saturatedBit = maskedImage.mask.getPlaneBitMask(saturatedMaskName)

    for amp in ccdExposure.getDetector():
        # Select the 2 top and bottom amplifiers around the saturated
        # core with a potential edge bleed by selecting the amplifiers
        # that have the same X coordinate as the saturated core.
        # As we don't care about the Y coordinate, we set it to the
        # center of the BBox.
        yBox = amp.getBBox().getCenter()[1]
        if amp.getBBox().contains(xCore, yBox):

            # Get the amp name
            ampName = amp.getName()

            # Because in ITLs the edge bleed happens on both edges
            # of the detector, we make a cutout around
            # both the top and bottom
            # edge bleed candidates around the saturated core.
            # We flip the cutout of the top amplifier
            # to then work with the same coordinates for both.
            # The way of selecting top vs bottom amp
            # is very specific to ITL.
            if ampName[:2] == 'C1':
                sliceImage = maskedImage.image.array[:200, :]
                sliceMask = maskedImage.mask.array[:200, :]
            elif ampName[:2] == 'C0':
                sliceImage = numpy.flipud(maskedImage.image.array[-200:, :])
                sliceMask = numpy.flipud(maskedImage.mask.array[-200:, :])

            # The middle columns of edge bleeds often have
            # high counts, so  we check there is an edge bleed
            # by looking at a small image up to 50 pixels from the edge
            # and around the saturated columns
            # of the saturated core, and checking its median is
            # above the sky background by itlEdgeBleedThreshold

            # If the centroid is too close to the edge of the detector
            # (within 5 pixels), we set the limit to the mean check
            # to the edge of the detector
            lowerRangeSmall = int(xCore)-5
            upperRangeSmall = int(xCore)+5
            if lowerRangeSmall < 0:
                lowerRangeSmall = 0
            if upperRangeSmall > xmax:
                upperRangeSmall = xmax
            ampImageBG = numpy.median(maskedImage[amp.getBBox()].image.array)
            edgeMedian = numpy.median(sliceImage[:50, lowerRangeSmall:upperRangeSmall])
            if edgeMedian > (ampImageBG + itlEdgeBleedThreshold):

                log.info("Found ITL edge bleed in amp %s, column %d.", ampName, xCore)

                # We need an estimate of the maximum width
                # of the edge bleed for our masking model
                # so we now estimate it by measuring the width of
                # areas above 60 percent of the saturation level
                # close to the edge,
                # in a cutout up to 100 pixels from the edge,
                # with a width of around the width of an amplifier.
                subImageXMin = int(xCore)-250
                subImageXMax = int(xCore)+250
                if subImageXMin < 0:
                    subImageXMin = 0
                elif subImageXMax > xmax:
                    subImageXMax = xmax

                subImage = sliceImage[:100, subImageXMin:subImageXMax]
                maxWidthEdgeBleed = numpy.max(numpy.sum(subImage > 0.45*satLevel,
                                                        axis=1))

                # Mask edge bleed with a decaying exponential model
                for y in range(200):
                    edgeBleedHalfWidth = \
                        int(((maxWidthEdgeBleed)*numpy.exp(-itlEdgeBleedModelConstant*y)
                             + widthSat)/2.)
                    lowerRange = int(xCore)-edgeBleedHalfWidth
                    upperRange = int(xCore)+edgeBleedHalfWidth
                    # If the edge bleed model goes outside the detector
                    # we set the limit for the masking
                    # to the edge of the detector
                    if lowerRange < 0:
                        lowerRange = 0
                    if upperRange > xmax:
                        upperRange = xmax
                    sliceMask[y, lowerRange:upperRange] |= saturatedBit


def maskITLSatSag(ccdExposure, fpCore, saturatedMaskName="SAT"):
    """Mask columns presenting saturation sag in saturated footprints in
    ITL detectors.

    Parameters
    ----------
    ccdExposure : `lsst.afw.image.Exposure`
        Exposure to apply masking to.
    fpCore : `lsst.afw.detection._detection.Footprint`
        Footprint of saturated core.
    saturatedMaskName : `str`, optional
        Mask name for saturation.
    """

    # TODO DM-49736: add a flux level check to apply masking

    maskedImage = ccdExposure.maskedImage
    saturatedBit = maskedImage.mask.getPlaneBitMask(saturatedMaskName)

    cc = numpy.sum(fpCore.getSpans().asArray(), axis=0)
    # Mask full columns that have 20 percent of the height of the footprint
    # saturated
    columnsToMaskFP = numpy.where(cc > fpCore.getSpans().asArray().shape[0]/5.)

    columnsToMask = [x + int(fpCore.getBBox().getMinX()) for x in columnsToMaskFP]
    maskedImage.mask.array[:, columnsToMask] |= saturatedBit


def maskITLDip(exposure, detectorConfig, maskPlaneNames=["SUSPECT", "ITL_DIP"], log=None):
    """Add mask bits according to the ITL dip model.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to do ITL dip masking.
    detectorConfig : `lsst.ip.isr.overscanAmpConfig.OverscanDetectorConfig`
        Configuration for this detector.
    maskPlaneNames : `list [`str`], optional
        Name of the ITL Dip mask planes.
    log : `logging.Logger`, optional
        If not set, a default logger will be used.
    """
    if detectorConfig.itlDipBackgroundFraction == 0.0:
        # Nothing to do.
        return

    if log is None:
        log = logging.getLogger(__name__)

    thresh = afwDetection.Threshold(
        exposure.mask.getPlaneBitMask("SAT"),
        afwDetection.Threshold.BITMASK,
    )
    fpList = afwDetection.FootprintSet(exposure.mask, thresh).getFootprints()

    heights = numpy.asarray([fp.getBBox().getHeight() for fp in fpList])

    largeHeights, = numpy.where(heights >= detectorConfig.itlDipMinHeight)

    if len(largeHeights) == 0:
        return

    # Get the approximate image background.
    approxBackground = numpy.median(exposure.image.array)
    maskValue = exposure.mask.getPlaneBitMask(maskPlaneNames)

    maskBak = exposure.mask.array.copy()
    nMaskedCols = 0

    for index in largeHeights:
        fp = fpList[index]
        center = fp.getCentroid()

        nSat = numpy.sum(fp.getSpans().asArray(), axis=0)
        width = numpy.sum(nSat > detectorConfig.itlDipMinHeight)

        if width < detectorConfig.itlDipMinWidth:
            continue

        width = numpy.clip(width, None, detectorConfig.itlDipMaxWidth)

        dipMax = detectorConfig.itlDipBackgroundFraction * approxBackground * width

        # Assume sky-noise dominated; we could add in read noise here.
        if dipMax < detectorConfig.itlDipMinBackgroundNoiseFraction * numpy.sqrt(approxBackground):
            continue

        minCol = int(center.getX() - (detectorConfig.itlDipWidthScale * width) / 2.)
        maxCol = int(center.getX() + (detectorConfig.itlDipWidthScale * width) / 2.)
        minCol = numpy.clip(minCol, 0, None)
        maxCol = numpy.clip(maxCol, None, exposure.mask.array.shape[1] - 1)

        log.info(
            "Found ITL dip (width %d; bkg %.2f); masking column %d to %d.",
            width,
            approxBackground,
            minCol,
            maxCol,
        )

        exposure.mask.array[:, minCol: maxCol + 1] |= maskValue

        nMaskedCols += (maxCol - minCol + 1)

    if nMaskedCols > detectorConfig.itlDipMaxColsPerImage:
        log.warning(
            "Too many (%d) columns would be masked on this image from dip masking; restoring original mask.",
            nMaskedCols,
        )
        exposure.mask.array[:, :] = maskBak


def interpolateFromMask(maskedImage, fwhm, growSaturatedFootprints=1,
                        maskNameList=['SAT'], fallbackValue=None, useLegacyInterp=True):
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

    interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=fallbackValue,
                          maskNameList=maskNameList, useLegacyInterp=useLegacyInterp)

    return maskedImage


def saturationCorrection(maskedImage, saturation, fwhm, growFootprints=1, interpolate=True, maskName='SAT',
                         fallbackValue=None, useLegacyInterp=True):
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
        interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=fallbackValue,
                              maskNameList=[maskName], useLegacyInterp=useLegacyInterp)

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


def updateVariance(maskedImage, gain, readNoise, replace=True):
    """Set the variance plane based on the image plane.

    The maskedImage must have units of `adu` (if gain != 1.0) or
    electron (if gain == 1.0). This routine will always produce a
    variance plane in the same units as the image.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image to process.  The variance plane is modified.
    gain : scalar
        The amplifier gain in electron/adu.
    readNoise : scalar
        The amplifier read noise in electron/pixel.
    replace : `bool`, optional
        Replace the current variance?  If False, the image
        variance will be added to the current variance plane.
    """
    var = maskedImage.variance
    if replace:
        var[:, :] = maskedImage.image
    else:
        var[:, :] += maskedImage.image
    var /= gain
    var += (readNoise/gain)**2


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
        tempImage = afwImage.ImageD(image, deep=True)

        nanIndex = numpy.isnan(tempImage.getArray())
        tempImage.getArray()[nanIndex] = 0.

        corr = numpy.zeros(image.array.shape, dtype=numpy.float64)
        prev_image = numpy.zeros(image.array.shape, dtype=numpy.float64)

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

            outArray = scipy.signal.convolve(
                tempImage.array,
                kernelImage.array,
                mode="same",
                method="fft",
            )
            tmpArray = tempImage.getArray()

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
                diff = numpy.sum(numpy.abs(prev_image - tmpArray), dtype=numpy.float64)

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
    corr = numpy.zeros(cFunc.shape, dtype=numpy.float64)

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
        gx = numpy.zeros(cFunc.shape, dtype=numpy.float64)
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
        tempImage = afwImage.ImageD(image, deep=True)

        nanIndex = numpy.isnan(tempImage.getArray())
        tempImage.getArray()[nanIndex] = 0.

        outImage = afwImage.ImageD(image.getDimensions())
        corr = numpy.zeros(image.array.shape, dtype=numpy.float64)
        prevImage = numpy.zeros(image.array.shape, dtype=numpy.float64)
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
        imean = numpy.mean(tempImage.getArray()[~nanIndex], dtype=numpy.float64)
        # subtract mean from image
        tempImage -= imean
        tempImage.array[nanIndex] = 0.0
        padArray = numpy.pad(tempImage.getArray(), ((0, kLy), (0, kLx)))
        outImage = afwImage.ImageD(numpy.pad(outImage.getArray(), ((0, kLy), (0, kLx))))
        # Convert array to afw image so afwMath.convolve works
        padImage = afwImage.ImageD(padArray.shape[1], padArray.shape[0])
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
                diff = numpy.sum(numpy.abs(prevImage - tmpArray), dtype=numpy.float64)

                if diff < threshold:
                    break
                prevImage[:, :] = tmpArray[:, :]

        image.getArray()[:] += corr[:]

    return diff, iteration


@contextmanager
def gainContext(exp, image, apply, gains=None, invert=False, isTrimmed=True):
    """Context manager that applies and removes gain.

    Parameters
    ----------
    exp : `lsst.afw.image.Exposure`
        Exposure to apply/remove gain.
    image : `lsst.afw.image.Image`
        Image to apply/remove gain.
    apply : `bool`
        If True, apply and remove the amplifier gain.
    gains : `dict` [`str`, `float`], optional
        A dictionary, keyed by amplifier name, of the gains to use.
        If gains is None, the nominal gains in the amplifier object are used.
    invert : `bool`, optional
        Invert the gains (e.g. convert electrons to adu temporarily)?
    isTrimmed : `bool`, optional
        Is this a trimmed exposure?

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
            sim = image.Factory(image, amp.getBBox() if isTrimmed else amp.getRawBBox())
            if gains:
                gain = gains[amp.getName()]
            else:
                gain = amp.getGain()
            if invert:
                sim /= gain
            else:
                sim *= gain

    try:
        yield exp
    finally:
        if apply:
            ccd = exp.getDetector()
            for amp in ccd:
                sim = image.Factory(image, amp.getBBox() if isTrimmed else amp.getRawBBox())
                if gains:
                    gain = gains[amp.getName()]
                else:
                    gain = amp.getGain()
                if invert:
                    sim *= gain
                else:
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


def applyGains(exposure, normalizeGains=False, ptcGains=None, isTrimmed=True):
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
    isTrimmed : `bool`, optional
        Is the input image trimmed?
    """
    ccd = exposure.getDetector()
    ccdImage = exposure.getMaskedImage()

    medians = []
    for amp in ccd:
        if isTrimmed:
            sim = ccdImage.Factory(ccdImage, amp.getBBox())
        else:
            sim = ccdImage.Factory(ccdImage, amp.getRawBBox())
        if ptcGains:
            sim *= ptcGains[amp.getName()]
        else:
            sim *= amp.getGain()

        if normalizeGains:
            medians.append(numpy.median(sim.getImage().getArray()))

    if normalizeGains:
        median = numpy.median(numpy.array(medians))
        for index, amp in enumerate(ccd):
            if isTrimmed:
                sim = ccdImage.Factory(ccdImage, amp.getBBox())
            else:
                sim = ccdImage.Factory(ccdImage, amp.getRawBBox())
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


def getExposureGains(exposure):
    """Get the per-amplifier gains used for this exposure.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The exposure to find gains for.

    Returns
    -------
    gains : `dict` [`str` `float`]
        Dictionary of gain values, keyed by amplifier name.
        Returns empty dict when detector is None.
    """
    det = exposure.getDetector()
    if det is None:
        return dict()

    metadata = exposure.getMetadata()
    gains = {}
    for amp in det:
        ampName = amp.getName()
        # The key may use the new LSST ISR or the old LSST prefix
        if (key1 := f"LSST ISR GAIN {ampName}") in metadata:
            gains[ampName] = metadata[key1]
        elif (key2 := f"LSST GAIN {ampName}") in metadata:
            gains[ampName] = metadata[key2]
        else:
            gains[ampName] = amp.getGain()
    return gains


def getExposureReadNoises(exposure):
    """Get the per-amplifier read noise used for this exposure.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The exposure to find read noise for.

    Returns
    -------
    readnoises : `dict` [`str` `float`]
        Dictionary of read noise values, keyed by amplifier name.
        Returns empty dict when detector is None.
    """
    det = exposure.getDetector()
    if det is None:
        return dict()

    metadata = exposure.getMetadata()
    readnoises = {}
    for amp in det:
        ampName = amp.getName()
        # The key may use the new LSST ISR or the old LSST prefix
        if (key1 := f"LSST ISR READNOISE {ampName}") in metadata:
            readnoises[ampName] = metadata[key1]
        elif (key2 := f"LSST READNOISE {ampName}") in metadata:
            readnoises[ampName] = metadata[key2]
        else:
            readnoises[ampName] = amp.getReadNoise()
    return readnoises


def isTrimmedExposure(exposure):
    """Check if the unused pixels (pre-/over-scan pixels) have
    been trimmed from an exposure.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The exposure to check.

    Returns
    -------
    result : `bool`
        True if the image is trimmed, else False.
    """
    return exposure.getDetector().getBBox() == exposure.getBBox()


def isTrimmedImage(image, detector):
    """Check if the unused pixels (pre-/over-scan pixels) have
    been trimmed from an image

    Parameters
    ----------
    image : `lsst.afw.image.Image`
        The image to check.
    detector : `lsst.afw.cameraGeom.Detector`
        The detector associated with the image.

    Returns
    -------
    result : `bool`
        True if the image is trimmed, else False.
    """
    return detector.getBBox() == image.getBBox()


def compareCameraKeywords(
    doRaiseOnCalibMismatch,
    cameraKeywordsToCompare,
    exposureMetadata,
    calib,
    calibName,
    log=None,
):
    """Compare header keywords to confirm camera states match.

    Parameters
    ----------
    doRaiseOnCalibMismatch : `bool`
        Raise on calibration mismatch?  Otherwise, log a warning.
    cameraKeywordsToCompare : `list` [`str`]
        List of camera keywords to compare.
    exposureMetadata : `lsst.daf.base.PropertyList`
        Header for the exposure being processed.
    calib : `lsst.afw.image.Exposure` or `lsst.ip.isr.IsrCalib`
        Calibration to be applied.
    calibName : `str`
        Calib type for log message.
    log : `logging.Logger`, optional
        Logger to handle messages.
    """
    try:
        calibMetadata = calib.metadata
    except AttributeError:
        return

    log = log if log else logging.getLogger(__name__)

    missingKeywords = []
    for keyword in cameraKeywordsToCompare:
        exposureValue = exposureMetadata.get(keyword, None)
        if exposureValue is None:
            log.debug("Sequencer keyword %s not found in exposure metadata.", keyword)
            continue

        calibValue = calibMetadata.get(keyword, None)

        # We don't log here if there is a missing keyword.
        if calibValue is None:
            missingKeywords.append(keyword)
            continue

        if exposureValue != calibValue:
            if doRaiseOnCalibMismatch:
                raise RuntimeError(
                    "Sequencer mismatch for %s [%s]: exposure: %s calib: %s",
                    calibName,
                    keyword,
                    exposureValue,
                    calibValue,
                )
            else:
                log.warning(
                    "Sequencer mismatch for %s [%s]: exposure: %s calib: %s",
                    calibName,
                    keyword,
                    exposureValue,
                    calibValue,
                )
                exposureMetadata[f"ISR {calibName.upper()} SEQUENCER MISMATCH"] = True

    if missingKeywords:
        log.info(
            "Calibration %s missing keywords %s, which were not checked.",
            calibName,
            ",".join(missingKeywords),
        )


def symmetrize(inputArray):
    """ Copy array over 4 quadrants prior to convolution.

    Parameters
    ----------
    inputarray : `numpy.array`
        Input array to symmetrize.

    Returns
    -------
    aSym : `numpy.array`
        Symmetrized array.
    """
    targetShape = list(inputArray.shape)
    r1, r2 = inputArray.shape[-1], inputArray.shape[-2]
    targetShape[-1] = 2*r1-1
    targetShape[-2] = 2*r2-1
    aSym = numpy.ndarray(tuple(targetShape))
    aSym[..., r2-1:, r1-1:] = inputArray
    aSym[..., r2-1:, r1-1::-1] = inputArray
    aSym[..., r2-1::-1, r1-1::-1] = inputArray
    aSym[..., r2-1::-1, r1-1:] = inputArray

    return aSym


class CustomFFTConvolution(object):
    """
    A class that performs image convolutions in Fourier space, using pyfftw.
    The constructor takes images as arguments and creates the FFTW plans.
    The convolutions are performed by the __call__ routine.
    This is faster than scipy.signal.fftconvolve, and it saves some transforms
    by allowing the same image to be convolved with several kernels.
    pyfftw does not accommodate float32 images, so everything
    should be double precision.

    Code adaped from :
    https://stackoverflow.com/questions/14786920/convolution-of-two-three-dimensional-arrays-with-padding-on-one-side-too-slow
    Code posted by Henry Gomersal
    """
    def __init__(self, im, kernel, threads=1):
        # Compute the minimum size of the convolution result.
        shape = (numpy.array(im.shape) + numpy.array(kernel.shape)) - 1

        # Find the next larger "fast size" for FFT computation.
        # This can provide a significant speedup.
        shape = numpy.array([pyfftw.next_fast_len(s) for s in shape])

        # Create FFTW building objects for the image and kernel.
        self.fftImageObj = pyfftw.builders.rfftn(im, s=shape, threads=threads)
        self.fftKernelObj = pyfftw.builders.rfftn(kernel, s=shape, threads=threads)
        self.ifftObj = pyfftw.builders.irfftn(
            self.fftImageObj.get_output_array(),
            s=shape,
            threads=threads,
        )

    def __call__(self, im, kernels):
        """
        Perform the convolution and trim the result to the
        size of the input image. If kernels is a list, then
        the routine returns a list of corresponding
        convolutions.
        """
        # Handle both a list of kernels and a single kernel.
        ks = [kernels] if type(kernels) is not list else kernels
        convolutions = []
        for k in ks:
            # Transform the image and the kernel using FFT.
            tim = self.fftImageObj(im)
            tk = self.fftKernelObj(k)

            # Multiply in Fourier space and perform the inverse FFT.
            convolution = self.ifftObj(tim * tk)
            # Trim the result to match the input image size, following
            # the 'same' policy of scipy.signal.fftconvolve.
            oy = k.shape[0] // 2
            ox = k.shape[1] // 2
            convolutions.append(convolution[oy:oy + im.shape[0], ox:ox + im.shape[1]].copy())
        # Return a single convolution if kernels was
        # not a list, otherwise return the list.
        return convolutions[0] if type(kernels) is not list else convolutions


def electrostaticBrighterFatterCorrection(exposure, electroBfDistortionMatrix, applyGain, gains=None):
    """
    Evaluates the correction of CCD images affected by the
    brighter-fatter effect, as described in
    https://arxiv.org/abs/2301.03274. Requires as input the result of
    an electrostatic fit to flat covariance data (or any other
    determination of pixel boundary shifts under the influence of a
    single electron).

    The filename refers to an input tuple that contains the
    boundary shifts for one electron. This file is produced by an
    electrostatic fit to data extracted from flat-field statistics,
    implemented in https://gitlab.in2p3.fr/astier/bfptc/tools/fit_cov.py.
    """

    # Use the symmetrize function to fill the four quadrants for each kernel
    r = electroBfDistortionMatrix.fitRange - 1
    aN = electroBfDistortionMatrix.aN
    aS = electroBfDistortionMatrix.aS
    aE = electroBfDistortionMatrix.aE
    aW = electroBfDistortionMatrix.aW

    # Initialize kN and kE arrays
    kN = numpy.zeros((2 * r + 1, 2 * r + 1))
    kE = numpy.zeros_like(kN)

    # Fill in the 4 quadrants for kN
    kN[r:, r:] = aN  # Quadrant 1 (bottom-right)
    kN[:r+1, r:] = numpy.flipud(aN)  # Quadrant 2 (top-right)
    kN[r:, :r+1] = numpy.fliplr(aS)  # Quadrant 3 (bottom-left)
    kN[:r+1, :r+1] = numpy.flipud(numpy.fliplr(aS))  # Quadrant 4 (top-left)

    # Fill in the 4 quadrants for kE
    kE[r:, r:] = aE  # Quadrant 1 (bottom-right)
    kE[:r+1, r:] = numpy.flipud(aW)  # Quadrant 2 (top-right)
    kE[r:, :r+1] = numpy.fliplr(aE)  # Quadrant 3 (bottom-left)
    kE[:r+1, :r+1] = numpy.flipud(numpy.fliplr(aW))  # Quadrant 4 (top-left)

    # Tweak the edges so that the sum rule applies.
    kN[:, 0] = -kN[:, -1]
    kE[0, :] = -kE[-1, :]

    # We use the normalization of Guyonnet et al. (2015),
    # which is compatible with the way the input file is produced.
    # The factor 1/2 is due to the fact that the charge distribution at the end
    # is twice the average, and the second 1/2 is due to
    # charge interpolation.
    kN *= 0.25
    kE *= 0.25

    # Indeed, i and j in the tuple refer to serial and parallel directions.
    # In most Python code, the image reads im[j, i], so we transpose:
    kN = kN.T
    kE = kE.T

    # Get the image to perform the correction
    image = exposure.getMaskedImage().getImage()

    # The image needs to be units of electrons/holes
    with gainContext(exposure, image, applyGain, gains):
        # Computes the correction and returns the "delta_image",
        # which should be subtracted from "im" in order to undo the BF effect.
        # The input image should be expressed in electrons
        im = image.getArray().copy()
        convolver = CustomFFTConvolution(im, kN)
        convolutions = convolver(im, [kN, kE])

        # The convolutions contain the boundary shifts (in pixel size units)
        # for [horizontal, vertical] boundaries.
        # We now compute the charge to move around.
        delta = numpy.zeros_like(im)
        boundaryCharge = numpy.zeros_like(im)

        # Horizontal boundaries (parallel direction).
        # We could use a more elaborate interpolator for estimating the
        # charge on the boundary.
        boundaryCharge[:-1, :] = im[1:, :] + im[:-1, :]
        # boundaryCharge[1:-2,:] = (9./8.)*(I[2:-1,:]+I[1:-2,:] -
        # (1./8.)*(I[0:-3,:]+I[3,:])

        # The charge to move around is the
        # product of the boundary shift (in pixel size units) times the
        # charge on the boundary (in charge per pixel unit).
        dq = boundaryCharge * convolutions[0]
        delta += dq

        # What is gained by a pixel is lost by its neighbor (the right one).
        delta[1:, :] -= dq[:-1, :]

        # Vertical boundaries.
        boundaryCharge = numpy.zeros_like(im)  # Reset to zero.
        boundaryCharge[:, :-1] = im[:, 1:] + im[:, :-1]
        dq = boundaryCharge * convolutions[1]
        delta += dq

        # What is gained by a pixel is lost by its neighbor.
        delta[:, 1:] -= dq[:, :-1]

        # TODO: One might check that delta.sum() ~ 0 (charge conservation).

        # Apply the correction to the original image
        exposure.image.array -= delta

    return exposure
