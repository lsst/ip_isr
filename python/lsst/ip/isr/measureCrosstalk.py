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
"""
Measure intra-detector crosstalk coefficients.
"""

__all__ = ["MeasureCrosstalkConfig", "MeasureCrosstalkTask"]


import itertools
import numpy as np

from lsstDebug import getDebugFrame
from lsst.afw.detection import FootprintSet, Threshold
from lsst.afw.display import getDisplay
from lsst.daf.persistence.butlerExceptions import NoResults
from lsst.pex.config import Config, Field, ListField, ConfigurableField
from lsst.pipe.base import CmdLineTask, Struct

from .crosstalk import calculateBackground, extractAmp, writeCrosstalkCoeffs
from .isrTask import IsrTask


class MeasureCrosstalkConfig(Config):
    """Configuration for MeasureCrosstalkTask."""
    isr = ConfigurableField(
        target=IsrTask,
        doc="Instrument signature removal task to use to process data."
    )
    threshold = Field(
        dtype=float,
        default=30000,
        doc="Minimum level of source pixels for which to measure crosstalk."
    )
    doRerunIsr = Field(
        dtype=bool,
        default=True,
        doc="Rerun the ISR, even if postISRCCD files are available?"
    )
    badMask = ListField(
        dtype=str,
        default=["SAT", "BAD", "INTRP"],
        doc="Mask planes to ignore when identifying source pixels."
    )
    rejIter = Field(
        dtype=int,
        default=3,
        doc="Number of rejection iterations for final coefficient calculation."
    )
    rejSigma = Field(
        dtype=float,
        default=2.0,
        doc="Rejection threshold (sigma) for final coefficient calculation."
    )

    def setDefaults(self):
        Config.setDefaults(self)
        # Set ISR processing to run up until we would be applying the CT
        # correction.  Applying subsequent stages may corrupt the signal.
        self.isr.doWrite = False
        self.isr.doOverscan = True
        self.isr.doAssembleCcd = True
        self.isr.doBias = True
        self.isr.doVariance = False  # This isn't used in the calculation below.
        self.isr.doLinearize = True  # This is the last ISR step we need.
        self.isr.doCrosstalk = False
        self.isr.doBrighterFatter = False
        self.isr.doDark = False
        self.isr.doStrayLight = False
        self.isr.doFlat = False
        self.isr.doFringe = False
        self.isr.doApplyGains = False
        self.isr.doDefect = True  # Masking helps remove spurious pixels.
        self.isr.doSaturationInterpolation = False
        self.isr.growSaturationFootprintSize = 0  # We want the saturation spillover: it's good signal.


class MeasureCrosstalkTask(CmdLineTask):
    """Measure intra-detector crosstalk.

    Notes
    -----
    The crosstalk this method measures assumes that when a bright
    pixel is found in one detector amplifier, all other detector
    amplifiers may see an increase in the same pixel location
    (relative to the readout amplifier) as these other pixels are read
    out at the same time.

    After processing each input exposure through a limited set of ISR
    stages, bright unmasked pixels above the threshold are identified.
    The potential CT signal is found by taking the ratio of the
    appropriate background-subtracted pixel value on the other
    amplifiers to the input value on the source amplifier.  If the
    source amplifier has a large number of bright pixels as well, the
    background level may be elevated, leading to poor ratio
    measurements.

    The set of ratios found between each pair of amplifiers across all
    input exposures is then gathered to produce the final CT
    coefficients.  The sigma-clipped mean and sigma are returned from
    these sets of ratios, with the coefficient to supply to the ISR
    CrosstalkTask() being the multiplicative inverse of these values.
    """
    ConfigClass = MeasureCrosstalkConfig
    _DefaultName = "measureCrosstalk"

    def __init__(self, *args, **kwargs):
        CmdLineTask.__init__(self, *args, **kwargs)
        self.makeSubtask("isr")

    @classmethod
    def _makeArgumentParser(cls):
        parser = super(MeasureCrosstalkTask, cls)._makeArgumentParser()
        parser.add_argument("--crosstalkName",
                            help="Name for this set of crosstalk coefficients", default="Unknown")
        parser.add_argument("--outputFileName",
                            help="Name of yaml file to which to write crosstalk coefficients")
        parser.add_argument("--dump-ratios", dest="dumpRatios",
                            help="Name of pickle file to which to write crosstalk ratios")
        return parser

    @classmethod
    def parseAndRun(cls, *args, **kwargs):
        """Implement scatter/gather

        Returns
        -------
        coeff : `numpy.ndarray`
            Crosstalk coefficients.
        coeffErr : `numpy.ndarray`
            Crosstalk coefficient errors.
        coeffNum : `numpy.ndarray`
            Number of pixels used for crosstalk measurement.
        """
        kwargs["doReturnResults"] = True
        results = super(MeasureCrosstalkTask, cls).parseAndRun(*args, **kwargs)
        task = cls(config=results.parsedCmd.config, log=results.parsedCmd.log)
        resultList = [rr.result for rr in results.resultList]
        if results.parsedCmd.dumpRatios:
            import pickle
            pickle.dump(resultList, open(results.parsedCmd.dumpRatios, "wb"))
        coeff, coeffErr, coeffNum = task.reduce(resultList)

        outputFileName = results.parsedCmd.outputFileName
        if outputFileName is not None:
            butler = results.parsedCmd.butler
            dataId = results.parsedCmd.id.idList[0]
            dataId["detector"] = butler.queryMetadata("raw", ["detector"], dataId)[0]

            det = butler.get('raw', dataId).getDetector()
            writeCrosstalkCoeffs(outputFileName, coeff, det=det,
                                 crosstalkName=results.parsedCmd.crosstalkName, indent=2)

        return Struct(
            coeff=coeff,
            coeffErr=coeffErr,
            coeffNum=coeffNum
        )

    def _getConfigName(self):
        """Disable config output."""
        return None

    def _getMetadataName(self):
        """Disable metdata output."""
        return None

    def runDataRef(self, dataRef):
        """Get crosstalk ratios for detector.

        Parameters
        ----------
        dataRef : `lsst.daf.peristence.ButlerDataRef`
            Data references for detectors to process.

        Returns
        -------
        ratios : `list` of `list` of `numpy.ndarray`
            A matrix of pixel arrays.
        """
        exposure = None
        if not self.config.doRerunIsr:
            try:
                exposure = dataRef.get("postISRCCD")
            except NoResults:
                pass

        if exposure is None:
            exposure = self.isr.runDataRef(dataRef).exposure

        dataId = dataRef.dataId
        return self.run(exposure, dataId=dataId)

    def run(self, exposure, dataId=None):
        """Extract and return cross talk ratios for an exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Image data to measure crosstalk ratios from.
        dataId   :
            Optional data ID for the exposure to process; used for logging.

        Returns
        -------
        ratios : `list` of `list` of `numpy.ndarray`
            A matrix of pixel arrays.
        """
        ratios = self.extractCrosstalkRatios(exposure)
        self.log.info("Extracted %d pixels from %s",
                      sum(len(jj) for ii in ratios for jj in ii if jj is not None), dataId)
        return ratios

    def extractCrosstalkRatios(self, exposure, threshold=None, badPixels=None):
        """Extract crosstalk ratios between different amplifiers.

        For pixels above ``threshold``, we calculate the ratio between
        each background-subtracted target amp and the source amp. We
        return a list of ratios for each pixel for each target/source
        combination, as a matrix of lists.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure for which to measure crosstalk.
        threshold : `float`, optional
            Lower limit on pixels for which we measure crosstalk.
        badPixels : `list` of `str`, optional
            Mask planes indicating a pixel is bad.

        Returns
        -------
        ratios : `list` of `list` of `numpy.ndarray`
           A matrix of pixel arrays. ``ratios[i][j]`` is an array of
           the fraction of the ``j``-th amp present on the ``i``-th amp.
           The value is `None` for the diagonal elements.

        Notes
        -----
        This has been moved into MeasureCrosstalkTask to allow for easier
        debugging.

        The lsstDebug.Info() method can be rewritten for __name__ =
        `lsst.ip.isr.measureCrosstalk`, and supports the parameters:

        debug.display['extract'] : `bool`
            Display the exposure under consideration, with the pixels used
            for crosstalk measurement indicated by the DETECTED mask plane.
        debug.display['pixels'] : `bool`
            Display a plot of the ratio calculated for each pixel used in this
            exposure, split by amplifier pairs.  The median value is listed
            for reference.
        """
        if threshold is None:
            threshold = self.config.threshold
        if badPixels is None:
            badPixels = list(self.config.badMask)

        mi = exposure.getMaskedImage()
        FootprintSet(mi, Threshold(threshold), "DETECTED")
        detected = mi.getMask().getPlaneBitMask("DETECTED")
        bad = mi.getMask().getPlaneBitMask(badPixels)
        bg = calculateBackground(mi, badPixels + ["DETECTED"])

        self.debugView('extract', exposure)

        ccd = exposure.getDetector()
        ratios = [[None for iAmp in ccd] for jAmp in ccd]

        for ii, iAmp in enumerate(ccd):
            iImage = mi[iAmp.getBBox()]
            iMask = iImage.mask.array
            select = (iMask & detected > 0) & (iMask & bad == 0) & np.isfinite(iImage.image.array)
            for jj, jAmp in enumerate(ccd):
                if ii == jj:
                    continue
                jImage = extractAmp(mi.image, jAmp, iAmp.getReadoutCorner(), isTrimmed=True)
                ratios[jj][ii] = (jImage.array[select] - bg)/iImage.image.array[select]
                self.debugPixels('pixels', iImage.image.array[select], jImage.array[select] - bg, ii, jj)
        return ratios

    def reduce(self, ratioList):
        """Combine ratios to produce crosstalk coefficients.

        Parameters
        ----------
        ratioList : `list` of `list` of `list` of `numpy.ndarray`
            A list of matrices of arrays; a list of results from
            `extractCrosstalkRatios`.

        Returns
        -------
        coeff : `numpy.ndarray`
            Crosstalk coefficients.
        coeffErr : `numpy.ndarray`
            Crosstalk coefficient errors.
        coeffNum : `numpy.ndarray`
            Number of pixels used for crosstalk measurement.

        Raises
        ------
        RuntimeError
            Raised if there is no crosstalk data available.

        Notes
        -----
        The lsstDebug.Info() method can be rewritten for __name__ =
        `lsst.ip.isr.measureCrosstalk`, and supports the parameters:

        debug.display['reduce'] : `bool`
            Display a histogram of the combined ratio measurements for
            a pair of source/target amplifiers from all input
            exposures/detectors.
        """
        numAmps = None
        for rr in ratioList:
            if rr is None:
                continue

            if numAmps is None:
                numAmps = len(rr)

            assert len(rr) == numAmps
            assert all(len(xx) == numAmps for xx in rr)

        if numAmps is None:
            raise RuntimeError("Unable to measure crosstalk signal for any amplifier")

        ratios = [[None for jj in range(numAmps)] for ii in range(numAmps)]
        for ii, jj in itertools.product(range(numAmps), range(numAmps)):
            if ii == jj:
                result = []
            else:
                values = [rr[ii][jj] for rr in ratioList]
                num = sum(len(vv) for vv in values)
                if num == 0:
                    self.log.warn("No values for matrix element %d,%d" % (ii, jj))
                    result = np.nan
                else:
                    result = np.concatenate([vv for vv in values if len(vv) > 0])
            ratios[ii][jj] = result
            self.debugRatios('reduce', ratios, ii, jj)
        coeff, coeffErr, coeffNum = self.measureCrosstalkCoefficients(ratios, self.config.rejIter,
                                                                      self.config.rejSigma)
        self.log.info("Coefficients:\n%s\n", coeff)
        self.log.info("Errors:\n%s\n", coeffErr)
        self.log.info("Numbers:\n%s\n", coeffNum)
        return coeff, coeffErr, coeffNum

    def measureCrosstalkCoefficients(self, ratios, rejIter=3, rejSigma=2.0):
        """Measure crosstalk coefficients from the ratios.

        Given a list of ratios for each target/source amp combination,
        we measure a sigma clipped mean and error.

        The coefficient errors returned are the standard deviation of
        the final set of clipped input ratios.

        Parameters
        ----------
        ratios : `list` of `list` of `numpy.ndarray`
           Matrix of arrays of ratios.
        rejIter : `int`
           Number of rejection iterations.
        rejSigma : `float`
           Rejection threshold (sigma).

        Returns
        -------
        coeff : `numpy.ndarray`
            Crosstalk coefficients.
        coeffErr : `numpy.ndarray`
            Crosstalk coefficient errors.
        coeffNum : `numpy.ndarray`
            Number of pixels for each measurement.

        Notes
        -----
        This has been moved into MeasureCrosstalkTask to allow for easier
        debugging.

        The lsstDebug.Info() method can be rewritten for __name__ =
        `lsst.ip.isr.measureCrosstalk`, and supports the parameters:

        debug.display['measure'] : `bool`
            Display a histogram of the combined ratio measurements for
            a pair of source/target amplifiers from the final set of
            clipped input ratios.
        """
        if rejIter is None:
            rejIter = self.config.rejIter
        if rejSigma is None:
            rejSigma = self.config.rejSigma

        numAmps = len(ratios)
        assert all(len(rr) == numAmps for rr in ratios)

        coeff = np.zeros((numAmps, numAmps))
        coeffErr = np.zeros((numAmps, numAmps))
        coeffNum = np.zeros((numAmps, numAmps), dtype=int)

        for ii, jj in itertools.product(range(numAmps), range(numAmps)):
            if ii == jj:
                values = [0.0]
            else:
                values = np.array(ratios[ii][jj])
                values = values[np.abs(values) < 1.0]  # Discard unreasonable values

            coeffNum[ii][jj] = len(values)

            if len(values) == 0:
                self.log.warn("No values for matrix element %d,%d" % (ii, jj))
                coeff[ii][jj] = np.nan
                coeffErr[ii][jj] = np.nan
            else:
                if ii != jj:
                    for rej in range(rejIter):
                        lo, med, hi = np.percentile(values, [25.0, 50.0, 75.0])
                        sigma = 0.741*(hi - lo)
                        good = np.abs(values - med) < rejSigma*sigma
                        if good.sum() == len(good):
                            break
                        values = values[good]

            coeff[ii][jj] = np.mean(values)
            coeffErr[ii][jj] = np.nan if coeffNum[ii][jj] == 1 else np.std(values)
            self.debugRatios('measure', ratios, ii, jj)

        return coeff, coeffErr, coeffNum

    def debugView(self, stepname, exposure):
        """Utility function to examine the image being processed.

        Parameters
        ----------
        stepname : `str`
            State of processing to view.
        exposure : `lsst.afw.image.Exposure`
            Exposure to view.
        """
        frame = getDebugFrame(self._display, stepname)
        if frame:
            display = getDisplay(frame)
            display.scale('asinh', 'zscale')
            display.mtv(exposure)

            prompt = "Press Enter to continue: "
            while True:
                ans = input(prompt).lower()
                if ans in ("", "c",):
                    break

    def debugPixels(self, stepname, pixelsIn, pixelsOut, i, j):
        """Utility function to examine the CT ratio pixel values.

        Parameters
        ----------
        stepname : `str`
            State of processing to view.
        pixelsIn : `np.ndarray`
            Pixel values from the potential crosstalk "source".
        pixelsOut : `np.ndarray`
            Pixel values from the potential crosstalk "victim".
        i : `int`
            Index of the source amplifier.
        j : `int`
            Index of the target amplifier.
        """
        frame = getDebugFrame(self._display, stepname)
        if frame:
            if i == j or len(pixelsIn) == 0 or len(pixelsOut) < 1:
                pass
            import matplotlib.pyplot as plot
            figure = plot.figure(1)
            figure.clear()

            axes = figure.add_axes((0.1, 0.1, 0.8, 0.8))
            axes.plot(pixelsIn, pixelsOut / pixelsIn, 'k+')
            plot.xlabel("Source amplifier pixel value")
            plot.ylabel("Measured pixel ratio")
            plot.title("(Source %d -> Victim %d) median ratio: %f" %
                       (i, j, np.median(pixelsOut / pixelsIn)))
            figure.show()

            prompt = "Press Enter to continue: "
            while True:
                ans = input(prompt).lower()
                if ans in ("", "c",):
                    break
            plot.close()

    def debugRatios(self, stepname, ratios, i, j):
        """Utility function to examine the final CT ratio set.

        Parameters
        ----------
        stepname : `str`
            State of processing to view.
        ratios : `List` of `List` of `np.ndarray`
            Array of measured CT ratios, indexed by source/victim
            amplifier.
        i : `int`
            Index of the source amplifier.
        j : `int`
            Index of the target amplifier.
        """
        frame = getDebugFrame(self._display, stepname)
        if frame:
            if i == j or ratios is None or len(ratios) < 1:
                pass

            RR = ratios[i][j]
            if RR is None or len(RR) < 1:
                pass

            value = np.mean(RR)

            import matplotlib.pyplot as plot
            figure = plot.figure(1)
            figure.clear()
            plot.hist(x=RR, bins='auto', color='b', rwidth=0.9)
            plot.xlabel("Measured pixel ratio")
            plot.axvline(x=value, color="k")
            plot.title("(Source %d -> Victim %d) clipped mean ratio: %f" % (i, j, value))
            figure.show()

            prompt = "Press Enter to continue: "
            while True:
                ans = input(prompt).lower()
                if ans in ("", "c",):
                    break
            plot.close()
