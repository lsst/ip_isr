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

import numpy

import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display as afwDisplay

from lsst.pipe.base import Task, Struct, timeMethod
from lsst.pex.config import Config, Field, ListField, ConfigField

afwDisplay.setDefaultMaskTransparency(75)


def getFrame():
    """Produce a new frame number each time"""
    getFrame.frame += 1
    return getFrame.frame


getFrame.frame = 0


class FringeStatisticsConfig(Config):
    """Options for measuring fringes on an exposure"""
    badMaskPlanes = ListField(dtype=str, default=["SAT"], doc="Ignore pixels with these masks")
    stat = Field(dtype=int, default=int(afwMath.MEDIAN), doc="Statistic to use")
    clip = Field(dtype=float, default=3.0, doc="Sigma clip threshold")
    iterations = Field(dtype=int, default=3, doc="Number of fitting iterations")
    rngSeedOffset = Field(dtype=int, default=0,
                          doc="Offset to the random number generator seed (full seed includes exposure ID)")


class FringeConfig(Config):
    """Fringe subtraction options"""
    filters = ListField(dtype=str, default=[], doc="Only fringe-subtract these filters")
    num = Field(dtype=int, default=30000, doc="Number of fringe measurements")
    small = Field(dtype=int, default=3, doc="Half-size of small (fringe) measurements (pixels)")
    large = Field(dtype=int, default=30, doc="Half-size of large (background) measurements (pixels)")
    iterations = Field(dtype=int, default=20, doc="Number of fitting iterations")
    clip = Field(dtype=float, default=3.0, doc="Sigma clip threshold")
    stats = ConfigField(dtype=FringeStatisticsConfig, doc="Statistics for measuring fringes")
    pedestal = Field(dtype=bool, default=False, doc="Remove fringe pedestal?")


class FringeTask(Task):
    """Task to remove fringes from a science exposure

    We measure fringe amplitudes at random positions on the science exposure
    and at the same positions on the (potentially multiple) fringe frames
    and solve for the scales simultaneously.
    """
    ConfigClass = FringeConfig
    _DefaultName = 'isrFringe'

    def readFringes(self, dataRef, assembler=None):
        """Read the fringe frame(s)

        The current implementation assumes only a single fringe frame and
        will have to be updated to support multi-mode fringe subtraction.

        This implementation could be optimised by persisting the fringe
        positions and fluxes.

        @param dataRef     Data reference for the science exposure
        @param assembler   An instance of AssembleCcdTask (for assembling fringe frames)
        @return Struct(fringes: fringe exposure or list of fringe exposures;
                       seed: 32-bit uint derived from ccdExposureId for random number generator
        """
        try:
            fringe = dataRef.get("fringe", immediate=True)
        except Exception as e:
            raise RuntimeError("Unable to retrieve fringe for %s: %s" % (dataRef.dataId, e))
        if assembler is not None:
            fringe = assembler.assembleCcd(fringe)

        seed = self.config.stats.rngSeedOffset + dataRef.get("ccdExposureId", immediate=True)
        # Seed for numpy.random.RandomState must be convertable to a 32 bit unsigned integer
        seed %= 2**32

        return Struct(fringes=fringe,
                      seed=seed)

    @timeMethod
    def run(self, exposure, fringes, seed=None):
        """Remove fringes from the provided science exposure.

        Primary method of FringeTask.  Fringes are only subtracted if the
        science exposure has a filter listed in the configuration.

        @param exposure    Science exposure from which to remove fringes
        @param fringes     Exposure or list of Exposures
        @param seed        32-bit unsigned integer for random number generator
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display

        if not self.checkFilter(exposure):
            return

        if seed is None:
            seed = self.config.stats.rngSeedOffset
        rng = numpy.random.RandomState(seed=seed)

        if not hasattr(fringes, '__iter__'):
            fringes = [fringes]

        mask = exposure.getMaskedImage().getMask()
        for fringe in fringes:
            fringe.getMaskedImage().getMask().__ior__(mask)
            if self.config.pedestal:
                self.removePedestal(fringe)

        positions = self.generatePositions(fringes[0], rng)
        fluxes = numpy.ndarray([self.config.num, len(fringes)])
        for i, f in enumerate(fringes):
            fluxes[:, i] = self.measureExposure(f, positions, title="Fringe frame")

        expFringes = self.measureExposure(exposure, positions, title="Science")
        solution, rms = self.solve(expFringes, fluxes)
        self.subtract(exposure, fringes, solution)
        if display:
            afwDisplay.Display(frame=getFrame()).mtv(exposure, title="Fringe subtracted")
        return solution, rms

    @timeMethod
    def runDataRef(self, exposure, dataRef, assembler=None):
        """Remove fringes from the provided science exposure.

        Retrieve fringes from butler dataRef provided and remove from
        provided science exposure.
        Fringes are only subtracted if the science exposure has a filter
        listed in the configuration.

        @param exposure    Science exposure from which to remove fringes
        @param dataRef     Data reference for the science exposure
        @param assembler   An instance of AssembleCcdTask (for assembling fringe frames)
        """
        if not self.checkFilter(exposure):
            return
        fringeStruct = self.readFringes(dataRef, assembler=assembler)
        return self.run(exposure, **fringeStruct.getDict())

    def checkFilter(self, exposure):
        """Check whether we should fringe-subtract the science exposure"""
        return exposure.getFilter().getName() in self.config.filters

    def removePedestal(self, fringe):
        """Remove pedestal from fringe exposure"""
        stats = afwMath.StatisticsControl()
        stats.setNumSigmaClip(self.config.stats.clip)
        stats.setNumIter(self.config.stats.iterations)
        mi = fringe.getMaskedImage()
        pedestal = afwMath.makeStatistics(mi, afwMath.MEDIAN, stats).getValue()
        self.log.info("Removing fringe pedestal: %f", pedestal)
        mi -= pedestal

    def generatePositions(self, exposure, rng):
        """Generate a random distribution of positions for measuring fringe amplitudes"""
        start = self.config.large
        num = self.config.num
        width = exposure.getWidth() - self.config.large
        height = exposure.getHeight() - self.config.large
        return numpy.array([rng.randint(start, width, size=num),
                            rng.randint(start, height, size=num)]).swapaxes(0, 1)

    @timeMethod
    def measureExposure(self, exposure, positions, title="Fringe"):
        """Measure fringe amplitudes for an exposure

        The fringe amplitudes are measured as the statistic within a square
        aperture.  The statistic within a larger aperture are subtracted so
        as to remove the background.

        @param exposure    Exposure to measure
        @param positions   Array of (x,y) for fringe measurement
        @param title       Title for display
        @return Array of fringe measurements
        """
        stats = afwMath.StatisticsControl()
        stats.setNumSigmaClip(self.config.stats.clip)
        stats.setNumIter(self.config.stats.iterations)
        stats.setAndMask(exposure.getMaskedImage().getMask().getPlaneBitMask(self.config.stats.badMaskPlanes))

        num = self.config.num
        fringes = numpy.ndarray(num)

        for i in range(num):
            x, y = positions[i]
            small = measure(exposure.getMaskedImage(), x, y, self.config.small, self.config.stats.stat, stats)
            large = measure(exposure.getMaskedImage(), x, y, self.config.large, self.config.stats.stat, stats)
            fringes[i] = small - large

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        if display:
            disp = afwDisplay.Display(frame=getFrame())
            disp.mtv(exposure, title=title)
            if False:
                with disp.Buffering():
                    for x, y in positions:
                        corners = numpy.array([[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, -1]]) + [[x, y]]
                        disp.line(corners*self.config.small, ctype=afwDisplay.GREEN)
                        disp.line(corners*self.config.large, ctype=afwDisplay.BLUE)

        return fringes

    @timeMethod
    def solve(self, science, fringes):
        """Solve (with iterative clipping) for the scale factors

        @param science     Array of science exposure fringe amplitudes
        @param fringes     Array of arrays of fringe frame fringe amplitudes
        @return Array of scale factors for the fringe frames
        """
        import lsstDebug
        doPlot = lsstDebug.Info(__name__).plot

        origNum = len(science)

        def emptyResult(msg=""):
            """Generate an empty result for return to the user

            There are no good pixels; doesn't matter what we return.
            """
            self.log.warn("Unable to solve for fringes: no good pixels%s", msg)
            out = [0]
            if len(fringes) > 1:
                out = out*len(fringes)
            return numpy.array(out), numpy.nan

        good = numpy.where(numpy.logical_and(numpy.isfinite(science), numpy.any(numpy.isfinite(fringes), 1)))
        science = science[good]
        fringes = fringes[good]
        oldNum = len(science)
        if oldNum == 0:
            return emptyResult()

        # Up-front rejection to get rid of extreme, potentially troublesome values
        # (e.g., fringe apertures that fall on objects).
        good = select(science, self.config.clip)
        for ff in range(fringes.shape[1]):
            good &= select(fringes[:, ff], self.config.clip)
        science = science[good]
        fringes = fringes[good]
        oldNum = len(science)
        if oldNum == 0:
            return emptyResult(" after initial rejection")

        for i in range(self.config.iterations):
            solution = self._solve(science, fringes)
            resid = science - numpy.sum(solution*fringes, 1)
            rms = stdev(resid)
            good = numpy.logical_not(abs(resid) > self.config.clip*rms)
            self.log.debug("Iteration %d: RMS=%f numGood=%d", i, rms, good.sum())
            self.log.debug("Solution %d: %s", i, solution)
            newNum = good.sum()
            if newNum == 0:
                return emptyResult(" after %d rejection iterations" % i)

            if doPlot:
                import matplotlib.pyplot as plot
                for j in range(fringes.shape[1]):
                    fig = plot.figure(j)
                    fig.clf()
                    try:
                        fig.canvas._tkcanvas._root().lift()  # == Tk's raise
                    except Exception:
                        pass
                    ax = fig.add_subplot(1, 1, 1)
                    adjust = science.copy()
                    others = set(range(fringes.shape[1]))
                    others.discard(j)
                    for k in others:
                        adjust -= solution[k]*fringes[:, k]
                    ax.plot(fringes[:, j], adjust, 'r.')
                    xmin = fringes[:, j].min()
                    xmax = fringes[:, j].max()
                    ymin = solution[j]*xmin
                    ymax = solution[j]*xmax
                    ax.plot([xmin, xmax], [ymin, ymax], 'b-')
                    ax.set_title("Fringe %d: %f" % (j, solution[j]))
                    ax.set_xlabel("Fringe amplitude")
                    ax.set_ylabel("Science amplitude")
                    ax.set_autoscale_on(False)
                    ax.set_xbound(lower=xmin, upper=xmax)
                    ax.set_ybound(lower=ymin, upper=ymax)
                    fig.show()
                while True:
                    ans = input("Enter or c to continue [chp]").lower()
                    if ans in ("", "c",):
                        break
                    if ans in ("p",):
                        import pdb
                        pdb.set_trace()
                    elif ans in ("h", ):
                        print("h[elp] c[ontinue] p[db]")

            if newNum == oldNum:
                # Not gaining
                break
            oldNum = newNum
            good = numpy.where(good)
            science = science[good]
            fringes = fringes[good]

        # Final solution without rejection
        solution = self._solve(science, fringes)
        self.log.info("Fringe solution: %s RMS: %f Good: %d/%d", solution, rms, len(science), origNum)
        return solution, rms

    def _solve(self, science, fringes):
        """Solve for the scale factors

        @param science     Array of science exposure fringe amplitudes
        @param fringes     Array of arrays of fringe frame fringe amplitudes
        @return Array of scale factors for the fringe frames
        """
        return afwMath.LeastSquares.fromDesignMatrix(fringes, science,
                                                     afwMath.LeastSquares.DIRECT_SVD).getSolution()

    def subtract(self, science, fringes, solution):
        """Subtract the fringes

        @param science     Science exposure
        @param fringes     List of fringe frames
        @param solution    Array of scale factors for the fringe frames
        """
        if len(solution) != len(fringes):
            raise RuntimeError("Number of fringe frames (%s) != number of scale factors (%s)" %
                               (len(fringes), len(solution)))

        for s, f in zip(solution, fringes):
            science.getMaskedImage().scaledMinus(s, f.getMaskedImage())


def measure(mi, x, y, size, statistic, stats):
    """Measure a statistic within an aperture

    @param mi          MaskedImage to measure
    @param x, y        Center for aperture
    @param size        Size of aperture
    @param statistic   Statistic to measure
    @param stats       StatisticsControl object
    @return Value of statistic within aperture
    """
    bbox = lsst.geom.Box2I(lsst.geom.Point2I(int(x) - size, int(y - size)),
                           lsst.geom.Extent2I(2*size, 2*size))
    subImage = mi.Factory(mi, bbox, afwImage.LOCAL)
    return afwMath.makeStatistics(subImage, statistic, stats).getValue()


def stdev(vector):
    """Calculate a robust standard deviation of an array of values

    @param vector  Array of values
    @return Standard deviation
    """
    q1, q3 = numpy.percentile(vector, (25, 75))
    return 0.74*(q3 - q1)


def select(vector, clip):
    """Select values within 'clip' standard deviations of the median

    Returns a boolean array.
    """
    q1, q2, q3 = numpy.percentile(vector, (25, 50, 75))
    return numpy.abs(vector - q2) < clip*0.74*(q3 - q1)
