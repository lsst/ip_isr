# This file is part of ip_isr.
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

__all__ = ["OverscanCorrectionTaskConfig", "OverscanCorrectionTask",
           "SerialOverscanCorrectionTaskConfig", "SerialOverscanCorrectionTask",
           "ParallelOverscanCorrectionTaskConfig", "ParallelOverscanCorrectionTask",
           ]

import numpy as np
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig

from .isr import fitOverscanImage
from .isrFunctions import makeThresholdMask


class OverscanCorrectionTaskConfigBase(pexConfig.Config):
    """Overscan correction options.
    """
    fitType = pexConfig.ChoiceField(
        dtype=str,
        doc="The method for fitting the overscan bias level.",
        default='MEDIAN',
        allowed={
            "POLY": "Fit ordinary polynomial to the longest axis of the overscan region",
            "CHEB": "Fit Chebyshev polynomial to the longest axis of the overscan region",
            "LEG": "Fit Legendre polynomial to the longest axis of the overscan region",
            "NATURAL_SPLINE": "Fit natural spline to the longest axis of the overscan region",
            "CUBIC_SPLINE": "Fit cubic spline to the longest axis of the overscan region",
            "AKIMA_SPLINE": "Fit Akima spline to the longest axis of the overscan region",
            "MEAN": "Correct using the mean of the overscan region",
            "MEANCLIP": "Correct using a clipped mean of the overscan region",
            "MEDIAN": "Correct using the median of the overscan region",
            "MEDIAN_PER_ROW": "Correct using the median per row of the overscan region",
        },
    )
    order = pexConfig.Field(
        dtype=int,
        doc=("Order of polynomial to fit if overscan fit type is a polynomial, "
             "or number of spline knots if overscan fit type is a spline."),
        default=1,
    )
    numSigmaClip = pexConfig.Field(
        dtype=float,
        doc="Rejection threshold (sigma) for collapsing overscan before fit",
        default=3.0,
    )
    maskPlanes = pexConfig.ListField(
        dtype=str,
        doc="Mask planes to reject when measuring overscan",
        default=['BAD', 'SAT'],
    )
    overscanIsInt = pexConfig.Field(
        dtype=bool,
        doc="Treat overscan as an integer image for purposes of fitType=MEDIAN"
            " and fitType=MEDIAN_PER_ROW.",
        default=True,
    )
    maxDeviation = pexConfig.Field(
        dtype=float,
        doc="Maximum deviation from median (in ADU) to mask in overscan correction.",
        default=1000.0, check=lambda x: x > 0,
    )


class OverscanCorrectionTaskConfig(OverscanCorrectionTaskConfigBase):
    doParallelOverscan = pexConfig.Field(
        dtype=bool,
        doc="Correct using parallel overscan after serial overscan correction?",
        default=False,
    )
    parallelOverscanMaskThreshold = pexConfig.Field(
        dtype=int,
        doc="Threshold above which pixels in the parallel overscan are masked as bleeds.",
        default=100000,
    )
    parallelOverscanMaskGrowSize = pexConfig.Field(
        dtype=int,
        doc="Masks created from saturated bleeds should be grown by this many "
            "pixels during construction of the parallel overscan mask. "
            "This value determined from the ITL chip in the LATISS camera",
        default=7,
    )
    leadingColumnsToSkip = pexConfig.Field(
        dtype=int,
        doc="Number of leading columns to skip in serial overscan correction.",
        default=0,
    )
    trailingColumnsToSkip = pexConfig.Field(
        dtype=int,
        doc="Number of trailing columns to skip in serial overscan correction.",
        default=0,
    )
    leadingRowsToSkip = pexConfig.Field(
        dtype=int,
        doc="Number of leading rows to skip in parallel overscan correction.",
        default=0,
    )
    trailingRowsToSkip = pexConfig.Field(
        dtype=int,
        doc="Number of trailing rows to skip in parallel overscan correction.",
        default=0,
    )


class OverscanCorrectionTaskBase(pipeBase.Task):
    """Base Correction task for overscan.

    This class contains a number of utilities that are easier to
    understand and use when they are not embedded in nested if/else
    loops.

    Parameters
    ----------
    statControl : `lsst.afw.math.StatisticsControl`, optional
        Statistics control object.
    """
    ConfigClass = OverscanCorrectionTaskConfigBase
    _DefaultName = "overscanBase"

    def __init__(self, statControl=None, **kwargs):
        super().__init__(**kwargs)
        self.allowDebug = True

        if statControl:
            self.statControl = statControl
        else:
            self.statControl = afwMath.StatisticsControl()
            self.statControl.setNumSigmaClip(self.config.numSigmaClip)
            self.statControl.setAndMask(afwImage.Mask.getPlaneBitMask(self.config.maskPlanes))

    def run(self, exposure, amp, isTransposed=False):
        raise NotImplementedError("run method is not defined for OverscanCorrectionTaskBase")

    def correctOverscan(self, exposure, amp, imageBBox, overscanBBox,
                        isTransposed=True, leadingToSkip=0, trailingToSkip=0):
        """Trim the exposure, fit the overscan, subtract the fit, and
        calculate statistics.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure containing the data.
        amp : `lsst.afw.cameraGeom.Amplifier`
            The amplifier that is to be corrected.
        imageBBox: `lsst.geom.Box2I`
            Bounding box of the image data that will have the overscan
            subtracted.  If parallel overscan will be performed, that
            area is added to the image bounding box during serial
            overscan correction.
        overscanBBox: `lsst.geom.Box2I`
            Bounding box for the overscan data.
        isTransposed: `bool`
            If true, then the data will be transposed before fitting
            the overscan.
        leadingToSkip : `int`, optional
            Leading rows/columns to skip.
        trailingToSkip : `int`, optional
            Leading rows/columns to skip.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            ``ampOverscanModel``
                Overscan model broadcast to the full image size.
                (`lsst.afw.image.Exposure`)
            ``overscanOverscanModel``
                Overscan model broadcast to the full overscan image
                size. (`lsst.afw.image.Exposure`)
            ``overscanImage``
                Overscan image with the overscan fit subtracted.
                (`lsst.afw.image.Exposure`)
            ``overscanValue``
                Overscan model. (`float` or `np.array`)
            ``overscanMean``
                Mean value of the overscan fit. (`float`)
            ``overscanMedian``
                Median value of the overscan fit. (`float`)
            ``overscanSigma``
                Standard deviation of the overscan fit. (`float`)
            ``overscanMeanResidual``
                Mean value of the overscan region after overscan
                subtraction. (`float`)
            ``overscanMedianResidual``
                Median value of the overscan region after overscan
                subtraction. (`float`)
            ``overscanSigmaResidual``
                Standard deviation of the overscan region after
                overscan subtraction. (`float`)
        """
        overscanBox = self.trimOverscan(exposure, amp, overscanBBox,
                                        leadingToSkip,
                                        trailingToSkip,
                                        transpose=isTransposed)
        overscanImage = exposure[overscanBox].getMaskedImage()
        overscanArray = overscanImage.image.array

        # Mask pixels.
        maskVal = overscanImage.mask.getPlaneBitMask(self.config.maskPlanes)
        overscanMask = ~((overscanImage.mask.array & maskVal) == 0)

        median = np.ma.median(np.ma.masked_where(overscanMask, overscanArray))
        bad = np.where(np.abs(overscanArray - median) > self.config.maxDeviation)
        overscanImage.mask.array[bad] = overscanImage.mask.getPlaneBitMask("SAT")

        # Do overscan fit.
        # CZW: Handle transposed correctly.
        overscanResults = self.fitOverscan(overscanImage, isTransposed=isTransposed)

        # Correct image region (and possibly parallel-overscan region).
        ampImage = exposure[imageBBox]
        ampOverscanModel = self.broadcastFitToImage(overscanResults.overscanValue,
                                                    ampImage.image.array,
                                                    transpose=isTransposed)
        ampImage.image.array -= ampOverscanModel

        # Correct overscan region (and possibly doubly-overscaned
        # region).
        overscanImage = exposure[overscanBBox]
        # CZW: Transposed?
        overscanOverscanModel = self.broadcastFitToImage(overscanResults.overscanValue,
                                                         overscanImage.image.array)
        self.debugView(overscanImage, overscanResults.overscanValue, amp, isTransposed=isTransposed)
        overscanImage.image.array -= overscanOverscanModel

        # Find residual fit statistics.
        stats = afwMath.makeStatistics(overscanImage.getMaskedImage(),
                                       afwMath.MEAN | afwMath.MEDIAN | afwMath.STDEVCLIP, self.statControl)
        residualMean = stats.getValue(afwMath.MEAN)
        residualMedian = stats.getValue(afwMath.MEDIAN)
        residualSigma = stats.getValue(afwMath.STDEVCLIP)

        return pipeBase.Struct(ampOverscanModel=ampOverscanModel,
                               overscanOverscanModel=overscanOverscanModel,
                               overscanImage=overscanImage,
                               overscanValue=overscanResults.overscanValue,

                               overscanMean=overscanResults.overscanMean,
                               overscanMedian=overscanResults.overscanMedian,
                               overscanSigma=overscanResults.overscanSigma,
                               overscanMeanResidual=residualMean,
                               overscanMedianResidual=residualMedian,
                               overscanSigmaResidual=residualSigma
                               )

    def broadcastFitToImage(self, overscanValue, imageArray, transpose=False):
        """Broadcast 0 or 1 dimension fit to appropriate shape.

        Parameters
        ----------
        overscanValue : `numpy.ndarray`, (Nrows, ) or scalar
            Overscan fit to broadcast.
        imageArray : `numpy.ndarray`, (Nrows, Ncols)
            Image array that we want to match.
        transpose : `bool`, optional
            Switch order to broadcast along the other axis.

        Returns
        -------
        overscanModel : `numpy.ndarray`, (Nrows, Ncols) or scalar
            Expanded overscan fit.

        Raises
        ------
        RuntimeError
            Raised if no axis has the appropriate dimension.
        """
        if isinstance(overscanValue, np.ndarray):
            overscanModel = np.zeros_like(imageArray)

            if transpose is False:
                if imageArray.shape[0] == overscanValue.shape[0]:
                    overscanModel[:, :] = overscanValue[:, np.newaxis]
                elif imageArray.shape[1] == overscanValue.shape[0]:
                    overscanModel[:, :] = overscanValue[np.newaxis, :]
                elif imageArray.shape[0] == overscanValue.shape[1]:
                    overscanModel[:, :] = overscanValue[np.newaxis, :]
                else:
                    raise RuntimeError(f"Could not broadcast {overscanValue.shape} to "
                                       f"match {imageArray.shape}")
            else:
                if imageArray.shape[1] == overscanValue.shape[0]:
                    overscanModel[:, :] = overscanValue[np.newaxis, :]
                elif imageArray.shape[0] == overscanValue.shape[0]:
                    overscanModel[:, :] = overscanValue[:, np.newaxis]
                elif imageArray.shape[1] == overscanValue.shape[1]:
                    overscanModel[:, :] = overscanValue[:, np.newaxis]
                else:
                    raise RuntimeError(f"Could not broadcast {overscanValue.shape} to "
                                       f"match {imageArray.shape}")
        else:
            overscanModel = overscanValue

        return overscanModel

    def trimOverscan(self, exposure, amp, bbox, skipLeading, skipTrailing, transpose=False):
        """Trim overscan region to remove edges.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure containing data.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Amplifier containing geometry information.
        bbox : `lsst.geom.Box2I`
            Bounding box of the overscan region.
        skipLeading : `int`
            Number of leading (towards data region) rows/columns to skip.
        skipTrailing : `int`
            Number of trailing (away from data region) rows/columns to skip.
        transpose : `bool`, optional
            Operate on the transposed array.

        Returns
        -------
        overscanArray : `numpy.array`, (N, M)
            Data array to fit.
        overscanMask : `numpy.array`, (N, M)
            Data mask.
        """
        dx0, dy0, dx1, dy1 = (0, 0, 0, 0)
        dataBBox = amp.getRawDataBBox()
        if transpose:
            if dataBBox.getBeginY() < bbox.getBeginY():
                dy0 += skipLeading
                dy1 -= skipTrailing
            else:
                dy0 += skipTrailing
                dy1 -= skipLeading
        else:
            if dataBBox.getBeginX() < bbox.getBeginX():
                dx0 += skipLeading
                dx1 -= skipTrailing
            else:
                dx0 += skipTrailing
                dx1 -= skipLeading

        overscanBBox = geom.Box2I(bbox.getBegin() + geom.Extent2I(dx0, dy0),
                                  geom.Extent2I(bbox.getWidth() - dx0 + dx1,
                                                bbox.getHeight() - dy0 + dy1))
        return overscanBBox

    def fitOverscan(self, overscanImage, isTransposed=False):
        if self.config.fitType in ('MEAN', 'MEANCLIP', 'MEDIAN'):
            # Transposition has no effect here.
            overscanResult = self.measureConstantOverscan(overscanImage)
            overscanValue = overscanResult.overscanValue
            overscanMean = overscanValue
            overscanMedian = overscanValue
            overscanSigma = 0.0
        elif self.config.fitType in ('MEDIAN_PER_ROW', 'POLY', 'CHEB', 'LEG',
                                     'NATURAL_SPLINE', 'CUBIC_SPLINE', 'AKIMA_SPLINE'):
            # Force transposes as needed
            overscanResult = self.measureVectorOverscan(overscanImage, isTransposed)
            overscanValue = overscanResult.overscanValue

            stats = afwMath.makeStatistics(overscanResult.overscanValue,
                                           afwMath.MEAN | afwMath.MEDIAN | afwMath.STDEVCLIP,
                                           self.statControl)
            overscanMean = stats.getValue(afwMath.MEAN)
            overscanMedian = stats.getValue(afwMath.MEDIAN)
            overscanSigma = stats.getValue(afwMath.STDEVCLIP)
        else:
            raise ValueError('%s : %s an invalid overscan type' %
                             ("overscanCorrection", self.config.fitType))

        return pipeBase.Struct(overscanValue=overscanValue,
                               overscanMean=overscanMean,
                               overscanMedian=overscanMedian,
                               overscanSigma=overscanSigma,
                               )

    @staticmethod
    def integerConvert(image):
        """Return an integer version of the input image.

        Parameters
        ----------
        image : `numpy.ndarray`, `lsst.afw.image.Image` or `MaskedImage`
            Image to convert to integers.

        Returns
        -------
        outI : `numpy.ndarray`, `lsst.afw.image.Image` or `MaskedImage`
            The integer converted image.

        Raises
        ------
        RuntimeError
            Raised if the input image could not be converted.
        """
        if hasattr(image, "image"):
            # Is a maskedImage:
            imageI = image.image.convertI()
            outI = afwImage.MaskedImageI(imageI, image.mask, image.variance)
        elif hasattr(image, "convertI"):
            # Is an Image:
            outI = image.convertI()
        elif hasattr(image, "astype"):
            # Is a numpy array:
            outI = image.astype(int)
        else:
            raise RuntimeError("Could not convert this to integers: %s %s %s",
                               image, type(image), dir(image))
        return outI

    def maskParallelOverscan(self, exposure, detector):
        """Mask the union of high values on all amplifiers in the parallel
        overscan.

        This operates on the image in-place.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            An untrimmed raw exposure.
        detector : `lsst.afw.cameraGeom.Detector`
            The detetor to use for amplifier geometry.
        """
        parallelMask = None

        for amp in detector:
            dataView = afwImage.MaskedImageF(exposure.getMaskedImage(),
                                             amp.getRawParallelOverscanBBox(),
                                             afwImage.PARENT)
            makeThresholdMask(
                maskedImage=dataView,
                threshold=self.config.parallelOverscanMaskThreshold,
                growFootprints=self.config.parallelOverscanMaskGrowSize,
                maskName="BAD"
            )
            if parallelMask is None:
                parallelMask = dataView.mask.array
            else:
                parallelMask |= dataView.mask.array
        for amp in detector:
            dataView = afwImage.MaskedImageF(exposure.getMaskedImage(),
                                             amp.getRawParallelOverscanBBox(),
                                             afwImage.PARENT)
            dataView.mask.array |= parallelMask

    # Constant methods
    def measureConstantOverscan(self, image):
        """Measure a constant overscan value.

        Parameters
        ----------
        image : `lsst.afw.image.Image` or `lsst.afw.image.MaskedImage`
            Image data to measure the overscan from.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Overscan result with entries:
            - ``overscanValue``: Overscan value to subtract (`float`)
            - ``isTransposed``: Orientation of the overscan (`bool`)
        """
        fitType = afwMath.stringToStatisticsProperty(self.config.fitType)
        overscanValue = afwMath.makeStatistics(image, fitType, self.statControl).getValue()

        return pipeBase.Struct(overscanValue=overscanValue,
                               isTransposed=False)

    # Vector correction utilities
    def getImageArray(self, image):
        """Extract the numpy array from the input image.

        Parameters
        ----------
        image : `lsst.afw.image.Image` or `lsst.afw.image.MaskedImage`
            Image data to pull array from.

        calcImage : `numpy.ndarray`
            Image data array for numpy operating.
        """
        if hasattr(image, "getImage"):
            calcImage = image.getImage().getArray()
            calcImage = np.ma.masked_where(image.getMask().getArray() & self.statControl.getAndMask(),
                                           calcImage)
        else:
            calcImage = image.getArray()
        return calcImage

    def maskOutliers(self, imageArray):
        """Mask  outliers in  a  row  of overscan  data  from  a robust  sigma
        clipping procedure.

        Parameters
        ----------
        imageArray : `numpy.ndarray`
            Image to filter along numpy axis=1.

        Returns
        -------
        maskedArray : `numpy.ma.masked_array`
            Masked image marking outliers.
        """
        lq, median, uq = np.percentile(np.ma.getdata(imageArray),
                                       [25.0, 50.0, 75.0], axis=1)
        axisMedians = median
        axisStdev = 0.74*(uq - lq)  # robust stdev

        # Replace pixels that have excessively large stdev values
        # with the median of stdev values.  A large stdev likely
        # indicates a bleed is spilling into the overscan.
        axisStdev = np.where(axisStdev > 2.0 * np.median(axisStdev),
                             np.median(axisStdev), axisStdev)

        # Mask pixels that are N-sigma away from their array medians.
        diff = np.abs(imageArray - axisMedians[:, np.newaxis])
        masked = np.ma.masked_where(diff > self.statControl.getNumSigmaClip()
                                    * axisStdev[:, np.newaxis], imageArray)

        return masked

    def fillMaskedPixels(self, overscanVector):
        """Fill masked/NaN pixels in the overscan.

        Parameters
        ----------
        overscanVector : `np.array` or `np.ma.masked_array`
            Overscan vector to fill.

        Returns
        -------
        overscanVector : `np.ma.masked_array`
            Filled vector.

        Notes
        -----
        Each maskSlice is a section of overscan with contiguous masks.
        Ideally this adds 5 pixels from the left and right of that
        mask slice, and takes the median of those values to fill the
        slice.  If this isn't possible, the median of all non-masked
        values is used.  The mask is removed for the pixels filled.
        """
        workingCopy = overscanVector
        if not isinstance(overscanVector, np.ma.MaskedArray):
            workingCopy = np.ma.masked_array(overscanVector,
                                             mask=~np.isfinite(overscanVector))

        defaultValue = np.median(workingCopy.data[~workingCopy.mask])
        for maskSlice in np.ma.clump_masked(workingCopy):
            neighborhood = []
            if maskSlice.start > 5:
                neighborhood.extend(workingCopy[maskSlice.start - 5:maskSlice.start].data)
            if maskSlice.stop < workingCopy.size - 5:
                neighborhood.extend(workingCopy[maskSlice.stop:maskSlice.stop+5].data)
            if len(neighborhood) > 0:
                workingCopy.data[maskSlice] = np.nanmedian(neighborhood)
                workingCopy.mask[maskSlice] = False
            else:
                workingCopy.data[maskSlice] = defaultValue
                workingCopy.mask[maskSlice] = False
        return workingCopy

    def collapseArray(self, maskedArray, fillMasked=True):
        """Collapse overscan array (and mask) to a 1-D vector of values.

        Parameters
        ----------
        maskedArray : `numpy.ma.masked_array`
            Masked array of input overscan data.
        fillMasked : `bool`, optional
            If true, fill any pixels that are masked with a median of
            neighbors.

        Returns
        -------
        collapsed : `numpy.ma.masked_array`
            Single dimensional overscan data, combined with the mean.

        """
        collapsed = np.mean(maskedArray, axis=1)
        if collapsed.mask.sum() > 0 and fillMasked:
            collapsed = self.fillMaskedPixels(collapsed)

        return collapsed

    def collapseArrayMedian(self, maskedArray):
        """Collapse overscan array (and mask) to a 1-D vector of using the
        correct integer median of row-values.

        Parameters
        ----------
        maskedArray : `numpy.ma.masked_array`
            Masked array of input overscan data.

        Returns
        -------
        collapsed : `numpy.ma.masked_array`
            Single dimensional overscan data, combined with the afwMath median.
        """
        integerMI = self.integerConvert(maskedArray)

        collapsed = []
        fitType = afwMath.stringToStatisticsProperty('MEDIAN')
        for row in integerMI:
            newRow = row.compressed()
            if len(newRow) > 0:
                rowMedian = afwMath.makeStatistics(newRow, fitType, self.statControl).getValue()
            else:
                rowMedian = np.nan
            collapsed.append(rowMedian)

        return np.array(collapsed)

    def splineFit(self, indices, collapsed, numBins):
        """Wrapper function to match spline fit API to polynomial fit API.

        Parameters
        ----------
        indices : `numpy.ndarray`
            Locations to evaluate the spline.
        collapsed : `numpy.ndarray`
            Collapsed overscan values corresponding to the spline
            evaluation points.
        numBins : `int`
            Number of bins to use in constructing the spline.

        Returns
        -------
        interp : `lsst.afw.math.Interpolate`
            Interpolation object for later evaluation.
        """
        if not np.ma.is_masked(collapsed):
            collapsed.mask = np.array(len(collapsed)*[np.ma.nomask])

        numPerBin, binEdges = np.histogram(indices, bins=numBins,
                                           weights=1 - collapsed.mask.astype(int))
        with np.errstate(invalid="ignore"):
            values = np.histogram(indices, bins=numBins,
                                  weights=collapsed.data*~collapsed.mask)[0]/numPerBin
            binCenters = np.histogram(indices, bins=numBins,
                                      weights=indices*~collapsed.mask)[0]/numPerBin

            if len(binCenters[numPerBin > 0]) < 5:
                self.log.warn("Cannot do spline fitting for overscan: %s valid points.",
                              len(binCenters[numPerBin > 0]))
                # Return a scalar value if we have one, otherwise
                # return zero.  This amplifier is hopefully already
                # masked.
                if len(values[numPerBin > 0]) != 0:
                    return float(values[numPerBin > 0][0])
                else:
                    return 0.0

            interp = afwMath.makeInterpolate(binCenters.astype(float)[numPerBin > 0],
                                             values.astype(float)[numPerBin > 0],
                                             afwMath.stringToInterpStyle(self.config.fitType))
        return interp

    @staticmethod
    def splineEval(indices, interp):
        """Wrapper function to match spline evaluation API to polynomial fit
        API.

        Parameters
        ----------
        indices : `numpy.ndarray`
            Locations to evaluate the spline.
        interp : `lsst.afw.math.interpolate`
            Interpolation object to use.

        Returns
        -------
        values : `numpy.ndarray`
            Evaluated spline values at each index.
        """

        return interp.interpolate(indices.astype(float))

    @staticmethod
    def maskExtrapolated(collapsed):
        """Create mask if edges are extrapolated.

        Parameters
        ----------
        collapsed : `numpy.ma.masked_array`
            Masked array to check the edges of.

        Returns
        -------
        maskArray : `numpy.ndarray`
            Boolean numpy array of pixels to mask.
        """
        maskArray = np.full_like(collapsed, False, dtype=bool)
        if np.ma.is_masked(collapsed):
            num = len(collapsed)
            for low in range(num):
                if not collapsed.mask[low]:
                    break
            if low > 0:
                maskArray[:low] = True
            for high in range(1, num):
                if not collapsed.mask[-high]:
                    break
            if high > 1:
                maskArray[-high:] = True
        return maskArray

    def measureVectorOverscan(self, image, isTransposed=False):
        """Calculate the 1-d vector overscan from the input overscan image.

        Parameters
        ----------
        image : `lsst.afw.image.MaskedImage`
            Image containing the overscan data.
        isTransposed : `bool`
            If true, the image has been transposed.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Overscan result with entries:

            ``overscanValue``
                Overscan value to subtract (`float`)
            ``maskArray``
                List of rows that should be masked as ``SUSPECT`` when the
                overscan solution is applied. (`list` [ `bool` ])
            ``isTransposed``
               Indicates if the overscan data was transposed during
               calcuation, noting along which axis the overscan should be
               subtracted. (`bool`)
        """
        calcImage = self.getImageArray(image)

        # operate on numpy-arrays from here
        if isTransposed:
            calcImage = np.transpose(calcImage)
        masked = self.maskOutliers(calcImage)

        if self.config.fitType == 'MEDIAN_PER_ROW':
            mi = afwImage.MaskedImageI(image.getBBox())
            masked = masked.astype(int)
            if isTransposed:
                masked = masked.transpose()

            mi.image.array[:, :] = masked.data[:, :]
            if bool(masked.mask.shape):
                mi.mask.array[:, :] = masked.mask[:, :]

            overscanVector = fitOverscanImage(mi, self.config.maskPlanes, isTransposed)
            overscanVector = self.fillMaskedPixels(overscanVector)
            maskArray = self.maskExtrapolated(overscanVector)
        else:
            collapsed = self.collapseArray(masked)

            num = len(collapsed)
            indices = 2.0*np.arange(num)/float(num) - 1.0

            poly = np.polynomial
            fitter, evaler = {
                'POLY': (poly.polynomial.polyfit, poly.polynomial.polyval),
                'CHEB': (poly.chebyshev.chebfit, poly.chebyshev.chebval),
                'LEG': (poly.legendre.legfit, poly.legendre.legval),
                'NATURAL_SPLINE': (self.splineFit, self.splineEval),
                'CUBIC_SPLINE': (self.splineFit, self.splineEval),
                'AKIMA_SPLINE': (self.splineFit, self.splineEval)
            }[self.config.fitType]

            # These are the polynomial coefficients, or an
            # interpolation object.
            coeffs = fitter(indices, collapsed, self.config.order)

            if isinstance(coeffs, float):
                self.log.warn("Using fallback value %f due to fitter failure. Amplifier will be masked.",
                              coeffs)
                overscanVector = np.full_like(indices, coeffs)
                maskArray = np.full_like(collapsed, True, dtype=bool)
            else:
                # Otherwise we can just use things as normal.
                overscanVector = evaler(indices, coeffs)
                maskArray = self.maskExtrapolated(collapsed)

        return pipeBase.Struct(overscanValue=np.array(overscanVector),
                               maskArray=maskArray,
                               isTransposed=isTransposed)

    def debugView(self, image, model, amp=None, isTransposed=True):
        """Debug display for the final overscan solution.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Input image the overscan solution was determined from.
        model : `numpy.ndarray` or `float`
            Overscan model determined for the image.
        amp : `lsst.afw.cameraGeom.Amplifier`, optional
            Amplifier to extract diagnostic information.
        isTransposed : `bool`, optional
            Does the data need to be transposed before display?
        """
        import lsstDebug
        if not lsstDebug.Info(__name__).display:
            return
        if not self.allowDebug:
            return

        calcImage = self.getImageArray(image)
        # CZW: Check that this is ok
        if isTransposed:
            calcImage = np.transpose(calcImage)
        masked = self.maskOutliers(calcImage)
        collapsed = self.collapseArray(masked, fillMasked=False)

        num = len(collapsed)
        indices = 2.0 * np.arange(num)/float(num) - 1.0
        indices = np.arange(num)

        if np.ma.is_masked(collapsed):
            collapsedMask = collapsed.mask
        else:
            collapsedMask = np.array(num*[np.ma.nomask])

        import matplotlib.pyplot as plot
        figure = plot.figure(1)
        figure.clear()
        axes = figure.add_axes((0.1, 0.1, 0.8, 0.8))
        axes.plot(indices[~collapsedMask], collapsed[~collapsedMask], 'k+')
        if collapsedMask.sum() > 0:
            axes.plot(indices[collapsedMask], collapsed.data[collapsedMask], 'b+')
        if isinstance(model, np.ndarray):
            plotModel = model
        else:
            plotModel = np.zeros_like(indices)
            plotModel += model

        axes.plot(indices, plotModel, 'r-')
        plot.xlabel("position along overscan region")
        plot.ylabel("pixel value/fit value")
        if amp:
            plot.title(f"{amp.getName()} DataX: "
                       f"[{amp.getRawDataBBox().getBeginX()}:{amp.getRawBBox().getEndX()}]"
                       f"OscanX: [{amp.getRawHorizontalOverscanBBox().getBeginX()}:"
                       f"{amp.getRawHorizontalOverscanBBox().getEndX()}] {self.config.fitType}")
        else:
            plot.title("No amp supplied.")
        figure.show()
        prompt = "Press Enter or c to continue [chp]..."
        while True:
            ans = input(prompt).lower()
            if ans in ("", " ", "c",):
                break
            elif ans in ("p", ):
                import pdb
                pdb.set_trace()
            elif ans in ('x', ):
                self.allowDebug = False
                break
            elif ans in ("h", ):
                print("[h]elp [c]ontinue [p]db e[x]itDebug")
        plot.close()


class OverscanCorrectionTask(OverscanCorrectionTaskBase):
    """Correction task for serial/parallel overscan.

    (Will be deprecated)

    This class contains a number of utilities that are easier to
    understand and use when they are not embedded in nested if/else
    loops.

    Parameters
    ----------
    statControl : `lsst.afw.math.StatisticsControl`, optional
        Statistics control object.
    """
    ConfigClass = OverscanCorrectionTaskConfig
    _DefaultName = "overscan"

    def run(self, exposure, amp, isTransposed=False):
        """Measure and remove serial/parallel overscan from an amplifier image.

        This will be deprecated.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Image data that will have the overscan corrections applied.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Amplifier to use for debugging purposes.
        isTransposed : `bool`, optional
            Is the image transposed, such that serial and parallel
            overscan regions are reversed?  Default is False.

        Returns
        -------
        overscanResults : `lsst.pipe.base.Struct`
            Result struct with components:

            ``imageFit``
                Value or fit subtracted from the amplifier image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanFit``
                Value or fit subtracted from the serial overscan image
                data (scalar or `lsst.afw.image.Image`).
            ``overscanImage``
                Image of the serial overscan region with the serial
                overscan correction applied
                (`lsst.afw.image.Image`). This quantity is used to
                estimate the amplifier read noise empirically.
            ``parallelOverscanFit``
                Value or fit subtracted from the parallel overscan
                image data (scalar, `lsst.afw.image.Image`, or None).
            ``parallelOverscanImage``
                Image of the parallel overscan region with the
                parallel overscan correction applied
                (`lsst.afw.image.Image` or None).
            ``overscanMean``
                Mean of the fit serial overscan region.
                This and the following values will be tuples of
                (serial, parallel) if doParallelOverscan=True.
            ``overscanMedian``
                Median of the fit serial overscan region.
            ``overscanSigma``
                Sigma of the fit serial overscan region.
            ``residualMean``
                Mean of the residual of the serial overscan region after
                correction.
            ``residualMedian``
                Median of the residual of the serial overscan region after
                correction.
            ``residualSigma``
                Mean of the residual of the serial overscan region after
                correction.


        Raises
        ------
        RuntimeError
            Raised if an invalid overscan type is set.
        """
        # Do Serial overscan first.
        serialOverscanBBox = amp.getRawSerialOverscanBBox()
        imageBBox = amp.getRawDataBBox()

        if self.config.doParallelOverscan:
            # We need to extend the serial overscan BBox to the full
            # size of the detector.
            parallelOverscanBBox = amp.getRawParallelOverscanBBox()
            imageBBox = imageBBox.expandedTo(parallelOverscanBBox)

            if isTransposed:
                serialOverscanBBox = geom.Box2I(
                    geom.Point2I(serialOverscanBBox.getMinX(), imageBBox.getEndY()),
                    geom.Extent2I(imageBBox.getWidth(), serialOverscanBBox.getHeight()),
                )
            else:
                serialOverscanBBox = geom.Box2I(
                    geom.Point2I(serialOverscanBBox.getMinX(),
                                 imageBBox.getMinY()),
                    geom.Extent2I(serialOverscanBBox.getWidth(),
                                  imageBBox.getHeight()),
                )

        serialResults = self.correctOverscan(
            exposure,
            amp,
            imageBBox,
            serialOverscanBBox,
            isTransposed=isTransposed,
            leadingToSkip=self.config.leadingColumnsToSkip,
            trailingToSkip=self.config.trailingColumnsToSkip,
        )
        overscanMean = serialResults.overscanMean
        overscanMedian = serialResults.overscanMedian
        overscanSigma = serialResults.overscanSigma
        residualMean = serialResults.overscanMeanResidual
        residualMedian = serialResults.overscanMedianResidual
        residualSigma = serialResults.overscanSigmaResidual

        # Do Parallel Overscan
        parallelResults = None
        if self.config.doParallelOverscan:
            # This does not need any extensions, as we'll only
            # subtract it from the data region.
            parallelOverscanBBox = amp.getRawParallelOverscanBBox()
            imageBBox = amp.getRawDataBBox()

            maskIm = exposure.getMaskedImage()
            maskIm = maskIm.Factory(maskIm, parallelOverscanBBox)

            # The serial overscan correction has removed some signal
            # from the parallel overscan region, but that is largely a
            # constant offset.  The collapseArray method now attempts
            # to fill fully masked columns with the median of
            # neighboring values, with a fallback to the median of the
            # correction in all other columns.  Filling with neighbor
            # values ensures that large variations in the parallel
            # overscan do not create new outlier points.  The
            # MEDIAN_PER_ROW method does this filling as a separate
            # operation, using the same method.
            parallelResults = self.correctOverscan(
                exposure,
                amp,
                imageBBox,
                parallelOverscanBBox,
                isTransposed=not isTransposed,
                leadingToSkip=self.config.leadingRowsToSkip,
                trailingToSkip=self.config.trailingRowsToSkip,
            )
            overscanMean = (overscanMean, parallelResults.overscanMean)
            overscanMedian = (overscanMedian, parallelResults.overscanMedian)
            overscanSigma = (overscanSigma, parallelResults.overscanSigma)
            residualMean = (residualMean, parallelResults.overscanMeanResidual)
            residualMedian = (residualMedian, parallelResults.overscanMedianResidual)
            residualSigma = (residualSigma, parallelResults.overscanSigmaResidual)

        parallelOverscanFit = parallelResults.overscanOverscanModel if parallelResults else None
        parallelOverscanImage = parallelResults.overscanImage if parallelResults else None

        return pipeBase.Struct(imageFit=serialResults.ampOverscanModel,
                               overscanFit=serialResults.overscanOverscanModel,
                               overscanImage=serialResults.overscanImage,

                               parallelOverscanFit=parallelOverscanFit,
                               parallelOverscanImage=parallelOverscanImage,
                               overscanMean=overscanMean,
                               overscanMedian=overscanMedian,
                               overscanSigma=overscanSigma,
                               residualMean=residualMean,
                               residualMedian=residualMedian,
                               residualSigma=residualSigma)


class SerialOverscanCorrectionTaskConfig(OverscanCorrectionTaskConfigBase):
    leadingToSkip = pexConfig.Field(
        dtype=int,
        doc="Number of leading values to skip in serial overscan correction.",
        default=0,
    )
    trailingToSkip = pexConfig.Field(
        dtype=int,
        doc="Number of trailing values to skip in serial overscan correction.",
        default=0,
    )


class SerialOverscanCorrectionTask(OverscanCorrectionTaskBase):
    """Correction task for serial overscan.

    Parameters
    ----------
    statControl : `lsst.afw.math.StatisticsControl`, optional
        Statistics control object.
    """
    ConfigClass = SerialOverscanCorrectionTaskConfig
    _DefaultName = "serialOverscan"

    def run(self, exposure, amp, isTransposed=False):
        """Measure and remove serial overscan from an amplifier image.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Image data that will have the overscan corrections applied.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Amplifier to use for debugging purposes.
        isTransposed : `bool`, optional
            Is the image transposed, such that serial and parallel
            overscan regions are reversed?  Default is False.

        Returns
        -------
        overscanResults : `lsst.pipe.base.Struct`
            Result struct with components:

            ``imageFit``
                Value or fit subtracted from the amplifier image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanFit``
                Value or fit subtracted from the serial overscan image
                data (scalar or `lsst.afw.image.Image`).
            ``overscanImage``
                Image of the serial overscan region with the serial
                overscan correction applied
                (`lsst.afw.image.Image`). This quantity is used to
                estimate the amplifier read noise empirically.
            ``overscanMean``
                Mean of the fit serial overscan region.
            ``overscanMedian``
                Median of the fit serial overscan region.
            ``overscanSigma``
                Sigma of the fit serial overscan region.
            ``residualMean``
                Mean of the residual of the serial overscan region after
                correction.
            ``residualMedian``
                Median of the residual of the serial overscan region after
                correction.
            ``residualSigma``
                Mean of the residual of the serial overscan region after
                correction.

        Raises
        ------
        RuntimeError
            Raised if an invalid overscan type is set.
        """
        serialOverscanBBox = amp.getRawSerialOverscanBBox()
        imageBBox = amp.getRawDataBBox()

        # We always want to extend the serial overscan bounding box to
        # the full size of the detector.
        parallelOverscanBBox = amp.getRawParallelOverscanBBox()
        imageBBox = imageBBox.expandedTo(parallelOverscanBBox)

        if isTransposed:
            serialOverscanBBox = geom.Box2I(
                geom.Point2I(serialOverscanBBox.getMinX(), imageBBox.getEndY()),
                geom.Extent2I(imageBBox.getWidth(), serialOverscanBBox.getHeight()),
            )
        else:
            serialOverscanBBox = geom.Box2I(
                geom.Point2I(serialOverscanBBox.getMinX(),
                             imageBBox.getMinY()),
                geom.Extent2I(serialOverscanBBox.getWidth(),
                              imageBBox.getHeight()),
            )

        results = self.correctOverscan(
            exposure,
            amp,
            imageBBox,
            serialOverscanBBox,
            isTransposed=isTransposed,
            leadingToSkip=self.config.leadingToSkip,
            trailingToSkip=self.config.trailingToSkip,
        )
        overscanMean = results.overscanMean
        overscanMedian = results.overscanMedian
        overscanSigma = results.overscanSigma
        residualMean = results.overscanMeanResidual
        residualMedian = results.overscanMedianResidual
        residualSigma = results.overscanSigmaResidual

        return pipeBase.Struct(
            imageFit=results.ampOverscanModel,
            overscanFit=results.overscanOverscanModel,
            overscanImage=results.overscanImage,
            overscanMean=overscanMean,
            overscanMedian=overscanMedian,
            overscanSigma=overscanSigma,
            residualMean=residualMean,
            residualMedian=residualMedian,
            residualSigma=residualSigma,
        )


class ParallelOverscanCorrectionTaskConfig(OverscanCorrectionTaskConfigBase):
    leadingToSkip = pexConfig.Field(
        dtype=int,
        doc="Number of leading values to skip in parallel overscan correction.",
        default=0,
    )
    trailingToSkip = pexConfig.Field(
        dtype=int,
        doc="Number of trailing values to skip in parallel overscan correction.",
        default=0,
    )


class ParallelOverscanCorrectionTask(OverscanCorrectionTaskBase):
    """Correction task for parallel overscan.

    Parameters
    ----------
    statControl : `lsst.afw.math.StatisticsControl`, optional
        Statistics control object.
    """
    ConfigClass = ParallelOverscanCorrectionTaskConfig
    _DefaultName = "parallelOverscan"

    def run(self, exposure, amp, isTransposed=False):
        """Measure and remove parallel overscan from an amplifier image.

        This method assumes that serial overscan has already been
        removed from the amplifier.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Image data that will have the overscan corrections applied.
        amp : `lsst.afw.cameraGeom.Amplifier`
            Amplifier to use for debugging purposes.
        isTransposed : `bool`, optional
            Is the image transposed, such that serial and parallel
            overscan regions are reversed?  Default is False.

        Returns
        -------
        overscanResults : `lsst.pipe.base.Struct`
            Result struct with components:

            ``imageFit``
                Value or fit subtracted from the amplifier image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanFit``
                Value or fit subtracted from the parallel overscan image
                data (scalar or `lsst.afw.image.Image`).
            ``overscanImage``
                Image of the parallel overscan region with the parallel
                overscan correction applied
                (`lsst.afw.image.Image`). This quantity is used to
                estimate the amplifier read noise empirically.
            ``overscanMean``
                Mean of the fit parallel overscan region.
            ``overscanMedian``
                Median of the fit parallel overscan region.
            ``overscanSigma``
                Sigma of the fit parallel overscan region.
            ``residualMean``
                Mean of the residual of the parallel overscan region after
                correction.
            ``residualMedian``
                Median of the residual of the parallel overscan region after
                correction.
            ``residualSigma``
                Mean of the residual of the parallel overscan region after
                correction.

        Raises
        ------
        RuntimeError
            Raised if an invalid overscan type is set.
        """
        # This does not need any extending, as we only subtract
        # from the data region.
        parallelOverscanBBox = amp.getRawParallelOverscanBBox()
        imageBBox = amp.getRawDataBBox()

        maskIm = exposure.getMaskedImage()
        maskIm = maskIm.Factory(maskIm, parallelOverscanBBox)

        # The serial overscan correction has removed some signal
        # from the parallel overscan region, but that is largely a
        # constant offset.  The collapseArray method now attempts
        # to fill fully masked columns with the median of
        # neighboring values, with a fallback to the median of the
        # correction in all other columns.  Filling with neighbor
        # values ensures that large variations in the parallel
        # overscan do not create new outlier points.  The
        # MEDIAN_PER_ROW method does this filling as a separate
        # operation, using the same method.
        results = self.correctOverscan(
            exposure,
            amp,
            imageBBox,
            parallelOverscanBBox,
            isTransposed=not isTransposed,
            leadingToSkip=self.config.leadingToSkip,
            trailingToSkip=self.config.trailingToSkip,
        )
        overscanMean = results.overscanMean
        overscanMedian = results.overscanMedian
        overscanSigma = results.overscanSigma
        residualMean = results.overscanMeanResidual
        residualMedian = results.overscanMedianResidual
        residualSigma = results.overscanSigmaResidual

        return pipeBase.Struct(
            imageFit=results.ampOverscanModel,
            overscanFit=results.overscanOverscanModel,
            overscanImage=results.overscanImage,
            overscanMean=overscanMean,
            overscanMedian=overscanMedian,
            overscanSigma=overscanSigma,
            residualMean=residualMean,
            residualMedian=residualMedian,
            residualSigma=residualSigma,
        )
