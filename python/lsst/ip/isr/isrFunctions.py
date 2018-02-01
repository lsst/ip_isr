from __future__ import division, print_function, absolute_import
from builtins import input
from builtins import range
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

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.pex.exceptions as pexExcept
import lsst.afw.cameraGeom as camGeom


def createPsf(fwhm):
    """Make a double Gaussian PSF

    @param[in] fwhm  FWHM of double Gaussian smoothing kernel
    @return measAlg.DoubleGaussianPsf
    """
    ksize = 4*int(fwhm) + 1
    return measAlg.DoubleGaussianPsf(ksize, ksize, fwhm/(2*math.sqrt(2*math.log(2))))

def transposeMaskedImage(maskedImage):
    """Make a transposed copy of a masked image

    @param[in] maskedImage  afw.image.MaskedImage to process
    @return transposed masked image
    """
    transposed = maskedImage.Factory(afwGeom.Extent2I(maskedImage.getHeight(), maskedImage.getWidth()))
    transposed.getImage().getArray()[:] = maskedImage.getImage().getArray().T
    transposed.getMask().getArray()[:] = maskedImage.getMask().getArray().T
    transposed.getVariance().getArray()[:] = maskedImage.getVariance().getArray().T
    return transposed


def interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=None):
    """Interpolate over defects specified in a defect list

    @param[in,out] maskedImage  masked image to process
    @param[in] defectList  defect list
    @param[in] fwhm  FWHM of double Gaussian smoothing kernel
    @param[in] fallbackValue  fallback value if an interpolated value cannot be determined;
                              if None then use clipped mean image value
    """
    psf = createPsf(fwhm)
    if fallbackValue is None:
        fallbackValue = afwMath.makeStatistics(maskedImage.getImage(), afwMath.MEANCLIP).getValue()
    if 'INTRP' not in maskedImage.getMask().getMaskPlaneDict():
        maskedImage.getMask.addMaskPlane('INTRP')
    measAlg.interpolateOverDefects(maskedImage, psf, defectList, fallbackValue, True)


def defectListFromFootprintList(fpList):
    """Compute a defect list from a footprint list, optionally growing the footprints

    @param[in] fpList  footprint list
    """
    defectList = []
    for fp in fpList:
        for bbox in afwDetection.footprintToBBoxList(fp):
            defect = measAlg.Defect(bbox)
            defectList.append(defect)
    return defectList


def transposeDefectList(defectList):
    """Make a transposed copy of a defect list

    @param[in] defectList  a list of defects (afw.meas.algorithms.Defect)
    @return a defect list with transposed defects
    """
    retDefectList = []
    for defect in defectList:
        bbox = defect.getBBox()
        nbbox = afwGeom.Box2I(afwGeom.Point2I(bbox.getMinY(), bbox.getMinX()),
                              afwGeom.Extent2I(bbox.getDimensions()[1], bbox.getDimensions()[0]))
        retDefectList.append(measAlg.Defect(nbbox))
    return retDefectList


def maskPixelsFromDefectList(maskedImage, defectList, maskName='BAD'):
    """Set mask plane based on a defect list

    @param[in,out] maskedImage  afw.image.MaskedImage to process; mask plane is updated
    @param[in] defectList  a list of defects (afw.meas.algorithms.Defect)
    @param[in] maskName  mask plane name
    """
    # mask bad pixels
    mask = maskedImage.getMask()
    bitmask = mask.getPlaneBitMask(maskName)
    for defect in defectList:
        bbox = defect.getBBox()
        afwGeom.SpanSet(bbox).clippedTo(mask.getBBox()).setMask(mask, bitmask)


def getDefectListFromMask(maskedImage, maskName):
    """Compute a defect list from a specified mask plane

    @param[in] maskedImage  masked image to process
    @param[in] maskName  mask plane name, or list of names
    """
    mask = maskedImage.getMask()
    thresh = afwDetection.Threshold(mask.getPlaneBitMask(maskName), afwDetection.Threshold.BITMASK)
    fpList = afwDetection.FootprintSet(mask, thresh).getFootprints()
    return defectListFromFootprintList(fpList)


def makeThresholdMask(maskedImage, threshold, growFootprints=1, maskName='SAT'):
    """Mask pixels based on threshold detection

    @param[in,out] maskedImage  afw.image.MaskedImage to process; the mask is altered
    @param[in] threshold  detection threshold
    @param[in] growFootprints  amount by which to grow footprints of detected regions
    @param[in] maskName  mask plane name
    @return a list of defects (meas.algrithms.Defect) of regions set in the mask.
    """
    # find saturated regions
    thresh = afwDetection.Threshold(threshold)
    fs = afwDetection.FootprintSet(maskedImage, thresh)

    if growFootprints > 0:
        fs = afwDetection.FootprintSet(fs, growFootprints)

    fpList = fs.getFootprints()
    # set mask
    mask = maskedImage.getMask()
    bitmask = mask.getPlaneBitMask(maskName)
    afwDetection.setMaskFromFootprintList(mask, fpList, bitmask)

    return defectListFromFootprintList(fpList)


def interpolateFromMask(maskedImage, fwhm, growFootprints=1, maskName='SAT', fallbackValue=None):
    """Interpolate over defects identified by a particular mask plane

    @param[in,out] maskedImage  afw.image.MaskedImage to process
    @param[in] fwhm  FWHM of double Gaussian smoothing kernel
    @param[in] growFootprints  amount by which to grow footprints of detected regions
    @param[in] maskName  mask plane name
    @param[in] fallbackValue  value of last resort for interpolation
    """
    mask = maskedImage.getMask()
    thresh = afwDetection.Threshold(mask.getPlaneBitMask(maskName), afwDetection.Threshold.BITMASK)
    fpSet = afwDetection.FootprintSet(mask, thresh)
    if growFootprints > 0:
        fpSet = afwDetection.FootprintSet(fpSet, rGrow=growFootprints, isotropic=False)
        # If we are interpolating over an area larger than the original masked region, we need
        # to expand the original mask bit to the full area to explain why we interpolated there.
        fpSet.setMask(mask, maskName)
    defectList = defectListFromFootprintList(fpSet.getFootprints())
    interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=fallbackValue)


def saturationCorrection(maskedImage, saturation, fwhm, growFootprints=1, interpolate=True, maskName='SAT',
                         fallbackValue=None):
    """Mark saturated pixels and optionally interpolate over them

    @param[in,out] maskedImage  afw.image.MaskedImage to process
    @param[in] saturation  saturation level (used as a detection threshold)
    @param[in] fwhm  FWHM of double Gaussian smoothing kernel
    @param[in] growFootprints  amount by which to grow footprints of detected regions
    @param[in] interpolate  interpolate over saturated pixels?
    @param[in] maskName  mask plane name
    @param[in] fallbackValue  value of last resort for interpolation
    """
    defectList = makeThresholdMask(
        maskedImage=maskedImage,
        threshold=saturation,
        growFootprints=growFootprints,
        maskName=maskName,
    )
    if interpolate:
        interpolateDefectList(maskedImage, defectList, fwhm, fallbackValue=fallbackValue)


def biasCorrection(maskedImage, biasMaskedImage):
    """Apply bias correction in place

    @param[in,out] maskedImage  masked image to correct
    @param[in] biasMaskedImage  bias, as a masked image
    """
    if maskedImage.getBBox(afwImage.LOCAL) != biasMaskedImage.getBBox(afwImage.LOCAL):
        raise RuntimeError("maskedImage bbox %s != biasMaskedImage bbox %s" %
                           (maskedImage.getBBox(afwImage.LOCAL), biasMaskedImage.getBBox(afwImage.LOCAL)))
    maskedImage -= biasMaskedImage


def darkCorrection(maskedImage, darkMaskedImage, expScale, darkScale, invert=False):
    """Apply dark correction in place

    maskedImage -= dark * expScaling / darkScaling

    @param[in,out] maskedImage  afw.image.MaskedImage to correct
    @param[in] darkMaskedImage  dark afw.image.MaskedImage
    @param[in] expScale  exposure scale
    @param[in] darkScale  dark scale
    @param[in] invert     if True, remove the dark from an already-corrected image
    """
    if maskedImage.getBBox(afwImage.LOCAL) != darkMaskedImage.getBBox(afwImage.LOCAL):
        raise RuntimeError("maskedImage bbox %s != darkMaskedImage bbox %s" %
                           (maskedImage.getBBox(afwImage.LOCAL), darkMaskedImage.getBBox(afwImage.LOCAL)))

    scale = expScale / darkScale
    if not invert:
        maskedImage.scaledMinus(scale, darkMaskedImage)
    else:
        maskedImage.scaledPlus(scale, darkMaskedImage)


def updateVariance(maskedImage, gain, readNoise):
    """Set the variance plane based on the image plane

    @param[in,out] maskedImage  afw.image.MaskedImage; image plane is read and variance plane is written
    @param[in] gain  amplifier gain (e-/ADU)
    @param[in] readNoise  amplifier read noise (ADU/pixel)
    """
    var = maskedImage.getVariance()
    var[:] = maskedImage.getImage()
    var /= gain
    var += readNoise**2


def flatCorrection(maskedImage, flatMaskedImage, scalingType, userScale=1.0, invert=False):
    """Apply flat correction in place

    @param[in,out] maskedImage  afw.image.MaskedImage to correct
    @param[in] flatMaskedImage  flat field afw.image.MaskedImage
    @param[in] scalingType  how to compute flat scale; one of 'MEAN', 'MEDIAN' or 'USER'
    @param[in] userScale  scale to use if scalingType is 'USER', else ignored
    @param[in] invert  if True, unflatten an already-flattened image instead.
    """
    if maskedImage.getBBox(afwImage.LOCAL) != flatMaskedImage.getBBox(afwImage.LOCAL):
        raise RuntimeError("maskedImage bbox %s != flatMaskedImage bbox %s" %
                           (maskedImage.getBBox(afwImage.LOCAL), flatMaskedImage.getBBox(afwImage.LOCAL)))

    # Figure out scale from the data
    # Ideally the flats are normalized by the calibration product pipelin, but this allows some flexibility
    # in the case that the flat is created by some other mechanism.
    if scalingType == 'MEAN':
        flatScale = afwMath.makeStatistics(flatMaskedImage.getImage(), afwMath.MEAN).getValue(afwMath.MEAN)
    elif scalingType == 'MEDIAN':
        flatScale = afwMath.makeStatistics(flatMaskedImage.getImage(),
                                           afwMath.MEDIAN).getValue(afwMath.MEDIAN)
    elif scalingType == 'USER':
        flatScale = userScale
    else:
        raise pexExcept.Exception('%s : %s not implemented' % ("flatCorrection", scalingType))

    if not invert:
        maskedImage.scaledDivides(1.0/flatScale, flatMaskedImage)
    else:
        maskedImage.scaledMultiplies(1.0/flatScale, flatMaskedImage)


def illuminationCorrection(maskedImage, illumMaskedImage, illumScale):
    """Apply illumination correction in place

    @param[in,out] maskedImage  afw.image.MaskedImage to correct
    @param[in] illumMaskedImage  illumination correction masked image
    @param[in] illumScale  scale value for illumination correction
    """
    if maskedImage.getBBox(afwImage.LOCAL) != illumMaskedImage.getBBox(afwImage.LOCAL):
        raise RuntimeError("maskedImage bbox %s != illumMaskedImage bbox %s" %
                           (maskedImage.getBBox(afwImage.LOCAL), illumMaskedImage.getBBox(afwImage.LOCAL)))

    maskedImage.scaledDivides(1./illumScale, illumMaskedImage)


def overscanCorrection(ampMaskedImage, overscanImage, fitType='MEDIAN', order=1, collapseRej=3.0,
                       statControl=None):
    """Apply overscan correction in place

    @param[in,out] ampMaskedImage  masked image to correct
    @param[in] overscanImage  overscan data as an afw.image.Image or afw.image.MaskedImage.
                              If a masked image is passed in the mask plane will be used
                              to constrain the fit of the bias level.
    @param[in] fitType  type of fit for overscan correction; one of:
                        - 'MEAN'
                        - 'MEDIAN'
                        - 'POLY' (ordinary polynomial)
                        - 'CHEB' (Chebyshev polynomial)
                        - 'LEG' (Legendre polynomial)
                        - 'NATURAL_SPLINE', 'CUBIC_SPLINE', 'AKIMA_SPLINE' (splines)
    @param[in] order  polynomial order or spline knots (ignored unless fitType
                      indicates a polynomial or spline)
    @param[in] collapseRej  Rejection threshold (sigma) for collapsing dimension of overscan
    @param[in] statControl  Statistics control object
    """
    ampImage = ampMaskedImage.getImage()
    if statControl is None:
        statControl = afwMath.StatisticsControl()
    if fitType == 'MEAN':
        offImage = afwMath.makeStatistics(overscanImage, afwMath.MEAN, statControl).getValue(afwMath.MEAN)
    elif fitType == 'MEDIAN':
        offImage = afwMath.makeStatistics(overscanImage, afwMath.MEDIAN, statControl).getValue(afwMath.MEDIAN)
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
        stdevBiasArr = 0.74*(percentiles[2] - percentiles[0]) # robust stdev
        diff = numpy.abs(biasArray - medianBiasArr[:, numpy.newaxis])
        biasMaskedArr = numpy.ma.masked_where(diff > collapseRej*stdevBiasArr[:, numpy.newaxis], biasArray)
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
                figure.close()

        offImage = ampImage.Factory(ampImage.getDimensions())
        offArray = offImage.getArray()
        if shortInd == 1:
            offArray[:, :] = fitBiasArr[:, numpy.newaxis]
        else:
            offArray[:, :] = fitBiasArr[numpy.newaxis, :]

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
        raise pexExcept.Exception('%s : %s an invalid overscan type' % \
            ("overscanCorrection", fitType))
    ampImage -= offImage


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

    All ``TransmissionCurve`` arguments are optional; if none are provided, the
    attached ``TransmissionCurve`` will have unit transmission everywhere.

    Returns
    -------
    combined : ``lsst.afw.image.TransmissionCurve``
        The TransmissionCurve attached to the exposure.
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
