#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2012 LSST Corporation.
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

import numpy

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9

from lsst.pipe.base import Task, Struct, timeMethod
from lsst.pex.config import Config, Field, ListField, ConfigField

def getFrame():
    """Produce a new frame number each time"""
    getFrame.frame += 1
    return getFrame.frame
getFrame.frame = 0

class FringeStatisticsConfig(Config):
    """Options for measuring fringes on an exposure"""
    badMaskPlanes = ListField(dtype=str, default=["SAT"], doc="Ignore pixels with these masks")
    stat = Field(dtype=int, default=afwMath.MEDIAN, doc="Statistic to use")
    clip = Field(dtype=float, default=3.0, doc="Sigma clip threshold")
    iterations = Field(dtype=int, default=3, doc="Number of fitting iterations")


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

    @timeMethod
    def run(self, exposure, dataRef, assembler=None):
        """Remove fringes from the provided science exposure

        Fringes are only subtracted if the science exposure has a filter
        listed in the configuration.

        @param exposure    Science exposure from which to remove fringes
        @param dataRef     Data reference for the science exposure
        @param assembler   An instance of AssembleCcdTask (for assembling fringe frames)
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display

        if not self.checkFilter(exposure):
            return
        fringes = self.readFringes(dataRef, assembler=assembler)
        expFringes = self.measureExposure(exposure, fringes.positions, title="Science")
        solution = self.solve(expFringes, fringes.fluxes)
        self.subtract(exposure, fringes.fringes, solution)
        if display:
            ds9.mtv(exposure, title="Fringe subtracted", frame=getFrame())

    def checkFilter(self, exposure):
        """Check whether we should fringe-subtract the science exposure"""
        return exposure.getFilter().getName() in self.config.filters

    def readFringes(self, dataRef, assembler=None):
        """Read the fringe frame(s) and measure fringe amplitudes.

        The current implementation assumes only a single fringe frame and
        will have to be updated to support multi-mode fringe subtraction.

        This implementation could be optimised by persisting the fringe
        positions and fluxes.

        @param dataRef     Data reference for the science exposure
        @param assembler   An instance of AssembleCcdTask (for assembling fringe frames)
        @return Struct(fringes: list of fringe frames;
                       fluxes: fringe amplitues;
                       positions: array of (x,y) for fringe amplitude measurements)
        """
        try:
            fringe = dataRef.get("fringe", immediate=True)
        except Exception as e:
            raise RuntimeError("Unable to retrieve fringe for %s: %s" % (dataRef.dataId, e))
        if assembler is not None:
            fringe = assembler.assembleCcd(fringe)

        if self.config.pedestal:
            stats = afwMath.StatisticsControl()
            stats.setNumSigmaClip(self.config.stats.clip)
            stats.setNumIter(self.config.stats.iterations)
            mi = fringe.getMaskedImage()
            pedestal = afwMath.makeStatistics(mi, afwMath.MEDIAN, stats).getValue()
            self.log.info("Removing fringe pedestal: %f" % pedestal)
            mi -= pedestal

        positions = self.generatePositions(fringe)
        fluxes = self.measureExposure(fringe, positions, title="Fringe frame")

        return Struct(fringes=[fringe],
                      fluxes=fluxes.reshape([len(positions), 1]),
                      positions=positions
                      )

    def generatePositions(self, exposure):
        """Generate a random distribution of positions for measuring fringe amplitudes"""
        start = self.config.large
        num = self.config.num
        width = exposure.getWidth() - self.config.large
        height = exposure.getHeight() - self.config.large
        return numpy.array([numpy.random.randint(start, width, size=num),
                            numpy.random.randint(start, height, size=num)]).swapaxes(0, 1)

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
            frame = getFrame()
            ds9.mtv(exposure, frame=frame, title=title)
            if False:
                with ds9.Buffering():
                    for x,y in positions:
                        corners = numpy.array([[-1, -1], [ 1, -1], [ 1,  1], [-1,  1], [-1, -1]]) + [[x, y]]
                        ds9.line(corners * self.config.small, frame=frame, ctype="green")
                        ds9.line(corners * self.config.large, frame=frame, ctype="blue")

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

        good = numpy.where(numpy.logical_and(numpy.isfinite(science), numpy.any(numpy.isfinite(fringes), 1)))
        science = science[good]
        fringes = fringes[good]
        oldNum = len(science)

        for i in range(self.config.iterations):
            solution = self._solve(science, fringes)
            resid = science - numpy.sum(solution * fringes, 1)
            rms = stdev(resid)
            good = numpy.logical_not(abs(resid) > self.config.clip * rms)
            self.log.logdebug("Iteration %d: RMS=%f numGood=%d" % (i, rms, good.sum()))
            self.log.logdebug("Solution %d: %s" % (i, solution))
            newNum = good.sum()

            if doPlot:
                import matplotlib.pyplot as plot
                for j in range(fringes.shape[1]):
                    fig = plot.figure(j)
                    fig.clf()
                    try:
                        fig.canvas._tkcanvas._root().lift() # == Tk's raise
                    except:
                        pass
                    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
                    adjust = science.copy()
                    others = set(range(fringes.shape[1]))
                    others.discard(j)
                    for k in others:
                        adjust -= solution[k] * fringes[:,k]
                    ax.plot(fringes[:,j], adjust, 'r.')
                    xmin = fringes[:,j].min()
                    xmax = fringes[:,j].max()
                    ymin = solution[j] * xmin
                    ymax = solution[j] * xmax
                    ax.plot([xmin, xmax], [ymin, ymax], 'b-')
                    ax.set_title("Fringe %d: %f" % (j, solution[j]))
                    ax.set_xlabel("Fringe amplitude")
                    ax.set_ylabel("Science amplitude")
                    ax.set_autoscale_on(False)
                    ax.set_xbound(lower=xmin, upper=xmax)
                    ax.set_ybound(lower=ymin, upper=ymax)
                    fig.show()
                while True:
                    ans = raw_input("Enter or c to continue [chp]").lower()
                    if ans in ("", "c",):
                        break
                    if ans in ("p",):
                        import pdb; pdb.set_trace()
                    elif ans in ("h", ):
                        print "h[elp] c[ontinue] p[db]"

            if newNum == oldNum:
                # Not gaining
                break
            oldNum = newNum
            good = numpy.where(good)
            science = science[good]
            fringes = fringes[good]

        # Final solution without rejection
        solution = self._solve(science, fringes)
        self.log.info("Fringe solution: %s RMS: %f Good: %d/%d" % (solution, rms, len(science), origNum))
        return solution

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
    bbox = afwGeom.Box2I(afwGeom.Point2I(int(x) - size, int(y - size)), afwGeom.Extent2I(2*size, 2*size))
    subImage = mi.Factory(mi, bbox, afwImage.LOCAL)
    return afwMath.makeStatistics(subImage, statistic, stats).getValue()

def stdev(vector):
    """Calculate a robust standard deviation of an array of values

    @param vector  Array of values
    @return Standard deviation
    """
    num = len(vector)
    vector = vector.copy()
    vector.sort()
    return 0.74 * (vector[int(0.75 * num)] - vector[int(0.25 * num)])

