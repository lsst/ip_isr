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

import numpy as np
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig

__all__ = ["OverscanCorrectionTaskConfig", "OverscanCorrectionTask"]


class OverscanCorrectionTaskConfig(pexConfig.Config):
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
        default=['SAT'],
    )
    overscanIsInt = pexConfig.Field(
        dtype=bool,
        doc="Treat overscan as an integer image for purposes of fitType=MEDIAN"
            " and fitType=MEDIAN_PER_ROW.",
        default=True,
    )


class OverscanCorrectionTask(pipeBase.Task):
    """Correction task for overscan.

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

    def __init__(self, statControl=None, **kwargs):
        super().__init__(**kwargs)
        if statControl:
            self.statControl = statControl
        else:
            self.statControl = afwMath.StatisticsControl()
            self.statControl.setNumSigmaClip(self.config.numSigmaClip)
            self.statControl.setAndMask(afwImage.Mask.getPlaneBitMask(self.config.maskPlanes))

    def run(self, ampImage, overscanImage):
        """Measure and remove an overscan from an amplifier image.

        Parameters
        ----------
        ampImage : `lsst.afw.image.Image`
            Image data that will have the overscan removed.
        overscanImage : `lsst.afw.image.Image`
            Overscan data that the overscan is measured from.

        Returns
        -------
        overscanResults : `lsst.pipe.base.Struct`
            Result struct with components:

            ``imageFit``
                Value or fit subtracted from the amplifier image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanFit``
                Value or fit subtracted from the overscan image data
                (scalar or `lsst.afw.image.Image`).
            ``overscanImage``
                Image of the overscan region with the overscan
                correction applied (`lsst.afw.image.Image`). This
                quantity is used to estimate the amplifier read noise
                empirically.

        Raises
        ------
        RuntimeError
            Raised if an invalid overscan type is set.

        """
        if self.config.fitType in ('MEAN', 'MEANCLIP', 'MEDIAN'):
            overscanResult = self.measureConstantOverscan(overscanImage)
            overscanValue = overscanResult.overscanValue
            offImage = overscanValue
            overscanModel = overscanValue
            maskSuspect = None
        elif self.config.fitType in ('MEDIAN_PER_ROW', 'POLY', 'CHEB', 'LEG',
                                     'NATURAL_SPLINE', 'CUBIC_SPLINE', 'AKIMA_SPLINE'):
            overscanResult = self.measureVectorOverscan(overscanImage)
            overscanValue = overscanResult.overscanValue
            maskArray = overscanResult.maskArray
            isTransposed = overscanResult.isTransposed

            offImage = afwImage.ImageF(ampImage.getDimensions())
            offArray = offImage.getArray()
            overscanModel = afwImage.ImageF(overscanImage.getDimensions())
            overscanArray = overscanModel.getArray()

            if hasattr(ampImage, 'getMask'):
                maskSuspect = afwImage.Mask(ampImage.getDimensions())
            else:
                maskSuspect = None

            if isTransposed:
                offArray[:, :] = overscanValue[np.newaxis, :]
                overscanArray[:, :] = overscanValue[np.newaxis, :]
                if maskSuspect:
                    maskSuspect.getArray()[:, maskArray] |= ampImage.getMask().getPlaneBitMask("SUSPECT")
            else:
                offArray[:, :] = overscanValue[:, np.newaxis]
                overscanArray[:, :] = overscanValue[:, np.newaxis]
                if maskSuspect:
                    maskSuspect.getArray()[maskArray, :] |= ampImage.getMask().getPlaneBitMask("SUSPECT")
        else:
            raise RuntimeError('%s : %s an invalid overscan type' %
                               ("overscanCorrection", self.config.fitType))

        self.debugView(overscanImage, overscanValue)

        ampImage -= offImage
        if maskSuspect:
            ampImage.getMask().getArray()[:, :] |= maskSuspect.getArray()[:, :]
        overscanImage -= overscanModel
        return pipeBase.Struct(imageFit=offImage,
                               overscanFit=overscanModel,
                               overscanImage=overscanImage,
                               edgeMask=maskSuspect)

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
            - ``maskArray``: Placeholder for a mask array (`list`)
            - ``isTransposed``: Orientation of the overscan (`bool`)
        """
        if self.config.fitType == 'MEDIAN':
            calcImage = self.integerConvert(image)
        else:
            calcImage = image

        fitType = afwMath.stringToStatisticsProperty(self.config.fitType)
        overscanValue = afwMath.makeStatistics(calcImage, fitType, self.statControl).getValue()

        return pipeBase.Struct(overscanValue=overscanValue,
                               maskArray=None,
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

    @staticmethod
    def transpose(imageArray):
        """Transpose input numpy array if necessary.

        Parameters
        ----------
        imageArray : `numpy.ndarray`
            Image data to transpose.

        Returns
        -------
        imageArray : `numpy.ndarray`
            Transposed image data.
        isTransposed : `bool`
            Indicates whether the input data was transposed.
        """
        if np.argmin(imageArray.shape) == 0:
            return np.transpose(imageArray), True
        else:
            return imageArray, False

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
        lq, median, uq = np.percentile(imageArray, [25.0, 50.0, 75.0], axis=1)
        axisMedians = median
        axisStdev = 0.74*(uq - lq)  # robust stdev

        diff = np.abs(imageArray - axisMedians[:, np.newaxis])
        return np.ma.masked_where(diff > self.statControl.getNumSigmaClip()
                                  * axisStdev[:, np.newaxis], imageArray)

    @staticmethod
    def collapseArray(maskedArray):
        """Collapse overscan array (and mask) to a 1-D vector of values.

        Parameters
        ----------
        maskedArray : `numpy.ma.masked_array`
            Masked array of input overscan data.

        Returns
        -------
        collapsed : `numpy.ma.masked_array`
            Single dimensional overscan data, combined with the mean.
        """
        collapsed = np.mean(maskedArray, axis=1)
        if collapsed.mask.sum() > 0:
            collapsed.data[collapsed.mask] = np.mean(maskedArray.data[collapsed.mask], axis=1)
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
            rowMedian = afwMath.makeStatistics(row, fitType, self.statControl).getValue()
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
            interp = afwMath.makeInterpolate(binCenters.astype(float)[numPerBin > 0],
                                             values.astype(float)[numPerBin > 0],
                                             afwMath.stringToInterpStyle(self.config.fitType))
        return interp

    @staticmethod
    def splineEval(indices, interp):
        """Wrapper function to match spline evaluation API to polynomial fit API.

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

    def measureVectorOverscan(self, image):
        """Calculate the 1-d vector overscan from the input overscan image.

        Parameters
        ----------
        image : `lsst.afw.image.MaskedImage`
            Image containing the overscan data.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Overscan result with entries:
            - ``overscanValue``: Overscan value to subtract (`float`)
            - ``maskArray`` : `list` [ `bool` ]
                List of rows that should be masked as ``SUSPECT`` when the
                overscan solution is applied.
            - ``isTransposed`` : `bool`
               Indicates if the overscan data was transposed during
               calcuation, noting along which axis the overscan should be
               subtracted.
        """
        calcImage = self.getImageArray(image)

        # operate on numpy-arrays from here
        calcImage, isTransposed = self.transpose(calcImage)
        masked = self.maskOutliers(calcImage)

        if self.config.fitType == 'MEDIAN_PER_ROW':
            overscanVector = self.collapseArrayMedian(masked)
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

            coeffs = fitter(indices, collapsed, self.config.order)
            overscanVector = evaler(indices, coeffs)
            maskArray = self.maskExtrapolated(collapsed)
        return pipeBase.Struct(overscanValue=np.array(overscanVector),
                               maskArray=maskArray,
                               isTransposed=isTransposed)

    def debugView(self, image, model):
        """Debug display for the final overscan solution.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Input image the overscan solution was determined from.
        model : `numpy.ndarray` or `float`
            Overscan model determined for the image.
        """
        import lsstDebug
        if not lsstDebug.Info(__name__).display:
            return

        calcImage = self.getImageArray(image)
        calcImage, isTransposed = self.transpose(calcImage)
        masked = self.maskOutliers(calcImage)
        collapsed = self.collapseArray(masked)

        num = len(collapsed)
        indices = 2.0 * np.arange(num)/float(num) - 1.0

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
        plot.xlabel("centered/scaled position along overscan region")
        plot.ylabel("pixel value/fit value")
        figure.show()
        prompt = "Press Enter or c to continue [chp]..."
        while True:
            ans = input(prompt).lower()
            if ans in ("", " ", "c",):
                break
            elif ans in ("p", ):
                import pdb
                pdb.set_trace()
            elif ans in ("h", ):
                print("[h]elp [c]ontinue [p]db")
        plot.close()
