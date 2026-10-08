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
"""Per-pixel mid-exposure times of LSSTCam exposures from the shutter motion.

The two LSSTCam shutter blades cross the focal plane in ~0.9 s, so the
mid-exposure time of a pixel depends on its position (by up to ~0.43 s across
the focal plane).  For one detector, `computeShutterTiming` gives the
flux-weighted mid-exposure time at the detector centre and a quadratic in
pixel coordinates for the rest of it.

Method:

1. *Blade trajectories.*  The camera fits each blade motion with the
   ``ThreeJerksModelv1`` model and writes the fit to the ``SHUTTER {OPEN,CLOSE}
   ...`` metadata cards (present since 2025-10-29).  Inverting it gives the
   time at which a blade edge reaches a given position.
2. *Beam-weighted crossing times.*  A pixel's converging ray bundle is
   ~60-70 mm wide at the shutter plane, so an edge uncovers (or covers) it
   gradually.  The flux-weighted crossing time of one edge is
   ``<T> = integral T(s) dF(s)``, with ``F`` the cumulative flux profile of
   the bundle along the blade travel (a raytrace beam table, LCA-20578,
   shipped by the obs package).  The mid-exposure time is
   ``(<T>_open + <T>_close) / 2``.
3. *Per-detector quadratic.*  These times are computed on a grid spanning
   the detector, and a quadratic in the pixel offsets along and across the
   blade direction is fit to them.

Times are TAI MJD.  The shutter motion and its Hall-sensor fit are described
in CTN-002, *Camera Shutter Motion Analysis* (https://ctn-002.lsst.io).
`lsst.ip.isr.ShutterMotionProfile` reads the same cards from an exposure.
"""

__all__ = [
    "ShutterTiming",
    "ShutterTimingConfig",
    "ShutterTimingFlag",
    "ShutterTimingStatus",
    "computeShutterTiming",
]

import dataclasses
import enum
import functools
import io
import math

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import ConvexHull

import lsst.geom
import lsst.pex.config as pexConfig
from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS, Detector
from lsst.daf.base import PropertyList
from lsst.resources import ResourcePath


_SECONDS_PER_DAY = 86400.0

#: Cumulative flux levels q of the beam table.
_BEAM_LEVELS = np.array([0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99])

#: Extrapolation margin (mm) outside the beam table's convex hull: the worst
#: science pixel (28.64 mm outside) plus 10 mm; NaN beyond.
_BEAM_MARGIN_MM = 28.64 + 10.0

#: 2-point Gauss-Legendre abscissae on [0, 1].
_GL2 = np.array([0.5 - 0.5 / np.sqrt(3.0), 0.5 + 0.5 / np.sqrt(3.0)])

_THREE_JERKS_V1 = "ThreeJerksModelv1"

#: Card suffixes of the _ThreeJerksParams fields, in order.
_FIT_CARDS = ("MODELSTARTTIME", "PIVOTPOINT1", "PIVOTPOINT2", "JERK0", "JERK1", "JERK2")

#: Encoder position (mm) of the focal-plane centre: x_ccs = c - e.
_ENCODER_CENTER = 375.0

#: Start position (mm) of a 750 -> 0 move (the cards carry none).
_NOMINAL_STROKE = 750.76

#: Start position (mm) of a 0 -> 750 move.
_NOMINAL_START_INCREASING = -0.05

#: Nominal Hall-fit ranges: pivot points (s), |jerks| (mm/s^3),
#: |ModelStartTime| (s) and the model displacement at t = 0.9 s (mm).
_PIVOT1_RANGE = (0.20, 0.25)
_PIVOT2_RANGE = (0.65, 0.70)
_ABS_JERK_RANGE = (30000.0, 38000.0)
_ABS_MODEL_START_TIME_RANGE = (0.0, 0.003)
_DISPLACEMENT_AT_0P9S_RANGE = (745.0, 757.0)

#: Range (ms) of close - open STARTTIME - requested exposure time; outside it
#: the shutter clocks disagree and both start times are re-anchored to the
#: header.
_CLOSE_MINUS_OPEN_MINUS_EXPOSURE_TIME_RANGE = (-2.0, 4.0)

#: Re-anchoring offsets (s): open STARTTIME - MJD-BEG, and close - open
#: STARTTIME - requested exposure time.
_HEADER_ANCHOR_OFFSETS = (0.00809, 0.00052)

#: Fit grid points along and across the blade axis.
_GRID = 9


class ShutterTimingStatus(enum.IntEnum):
    """Whether a detector has a corrected time."""

    OK = 0
    """Corrected.  The quadratic is within 0.1 ms of the exact times on most
    science detectors, and within ~2 ms at the outer corners
    (``maxAbsResidual``).
    """
    UNAVAILABLE = 2
    """No corrected time: callers keep the header midpoint."""


class ShutterTimingFlag(enum.IntFlag):
    """Details of a result (fixed bit values)."""

    NONE = 0
    NO_PROFILE = 1
    """No shutter cards (exposures before 2025-10-29): UNAVAILABLE."""
    PARAM_RANGE = 2
    """A Hall fit outside the nominal ranges: UNAVAILABLE."""
    CLOCK_CLOSE_VS_OPEN = 16
    """The two shutter clocks disagree: both start times are re-anchored to
    MJD-BEG and the requested exposure time.
    """
    BEAM_EXTRAPOLATED = 2048
    """Part of the detector lies outside the beam table's coverage
    (informational).
    """


class ShutterTimingConfig(pexConfig.Config):
    """Configuration of `computeShutterTiming` (the obs package sets
    ``beamFile``).
    """

    beamFile = pexConfig.Field(
        dtype=str, default="",
        doc="Path or URI of the shutter-plane beam table (raytrace flux quantiles, LCA-20578 format).",
    )


@dataclasses.dataclass(frozen=True)
class _ThreeJerksParams:
    """One ThreeJerksModelv1 fit, in the order of the header cards."""

    modelStartTime: float
    """Delay (s) from the motion start time to model time zero."""

    pivot1: float
    """Model time (s) at which the jerk changes from ``jerk0`` to ``jerk1``."""

    pivot2: float
    """Model time (s) at which the jerk changes from ``jerk1`` to ``jerk2``."""

    jerk0: float
    """Jerk (mm/s^3) before ``pivot1``."""

    jerk1: float
    """Jerk (mm/s^3) between ``pivot1`` and ``pivot2``."""

    jerk2: float
    """Jerk (mm/s^3) after ``pivot2``."""


class _ThreeJerksTrajectory:
    """ThreeJerksModelv1 trajectory of one blade motion.

    The camera fits each blade motion, from its Hall-sensor positions, with
    this model (CTN-002): piecewise-constant jerk ``j0`` on ``[0, t1)``,
    ``j1`` on ``[t1, t2)``, ``j2`` on ``[t2, inf)``, from rest at model time
    0 (``(T - startTime) - ModelStartTime``, s).  It gives the displacement
    ``d(t) >= 0`` along the travel direction; the encoder position is
    ``startPosition + travelSign * d(t)``.

    Parameters
    ----------
    params : ``_ThreeJerksParams``
        The fit.
    startPosition : `float`
        Encoder position (mm) at model time zero.
    travelSign : `int`
        +1 if the encoder position increases during the move, else -1.

    Raises
    ------
    ValueError
        Raised if the pivots do not satisfy ``0 <= pivot1 <= pivot2``.
    """

    def __init__(self, params: _ThreeJerksParams, startPosition: float, travelSign: int):
        t1, t2 = params.pivot1, params.pivot2
        if not 0.0 <= t1 <= t2:
            raise ValueError(f"invalid pivots (need 0 <= pivot1 <= pivot2): {t1}, {t2}")

        self.startPosition = startPosition
        self.travelSign = travelSign

        # Segment start times and jerks.
        self._tb = np.array([0.0, t1, t2])
        self._j = np.array([params.jerk0, params.jerk1, params.jerk2])

        # Displacement, velocity and acceleration at the start of each
        # segment, integrating from rest.
        state = np.zeros((3, 3))
        s = v = a = 0.0
        for k in range(3):
            state[k] = (s, v, a)
            if k < 2:
                h, j = self._tb[k + 1] - self._tb[k], self._j[k]
                s, v, a = (
                    s + v*h + a*h*h/2 + j*h**3/6,
                    v + a*h + j*h*h/2,
                    a + j*h,
                )
        self._state = state

        # End of the monotonic branch, where the inversion is defined.
        self.tStop = self._branchEnd()
        self.dMax = float(self.sva(self.tStop)[0]) if np.isfinite(self.tStop) else np.inf

    def _branchEnd(self) -> float:
        """End (model time, s) of the monotonic branch: the first velocity
        zero or velocity minimum after 0; 0 if ``j0 <= 0``; inf if none.
        """
        if self._j[0] <= 0:
            return 0.0

        ends = np.append(self._tb[1:], np.inf)
        for k in (1, 2):
            _, v0, a0 = self._state[k]
            j = self._j[k]
            hMax = ends[k] - self._tb[k]

            # Candidate times within this segment: velocity zeros, and the
            # velocity minimum.
            cands = []
            if j != 0:
                disc = a0*a0 - 2.0*j*v0
                if disc >= 0:
                    sq = np.sqrt(disc)
                    cands += [(-a0 - sq)/j, (-a0 + sq)/j]
            elif a0 != 0:
                cands.append(-v0/a0)
            if a0 < 0 < j:
                cands.append(-a0/j)

            cands = [h for h in cands if 0 < h <= hMax]
            if cands:
                return float(self._tb[k] + min(cands))

        return np.inf

    def sva(self, dt: np.ndarray | float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate the trajectory.

        Parameters
        ----------
        dt : `numpy.ndarray` or `float`
            Model time (s).

        Returns
        -------
        s, v, a : `tuple` [`numpy.ndarray`]
            Displacement (mm), velocity (mm/s) and acceleration (mm/s^2);
            zero for ``dt < 0``.
        """
        # The segment of each time, and the time since its start.
        tc = np.maximum(dt, 0.0)
        k = np.clip(np.searchsorted(self._tb, tc, side="right") - 1, 0, 2)
        h = tc - self._tb[k]

        # Constant-jerk motion from that segment's starting state.
        s0, v0, a0 = (self._state[k, i] for i in range(3))
        j = self._j[k]
        s = s0 + h*(v0 + h*(a0/2 + h*j/6))
        v = v0 + h*(a0 + h*j/2)
        a = a0 + h*j

        return s, v, a

    def timeOfDisplacement(self, d: np.ndarray) -> np.ndarray:
        """Invert the trajectory on its monotonic branch.

        Parameters
        ----------
        d : `numpy.ndarray`
            Displacement (mm).

        Returns
        -------
        t : `numpy.ndarray`
            Model time (s) at which the displacement is ``d``, by bisection
            within its jerk segment; NaN outside ``[0, dMax]``.
        """
        out = np.full(d.shape, np.nan)
        valid = (d >= 0.0) & (d <= self.dMax)
        if self.tStop <= 0.0 or not np.any(valid):
            return out
        dv = d[valid]

        # End of the branch; if the motion never stops, far enough to reach
        # the largest d (s(t) may never reach it).
        tEnd = self.tStop
        if not np.isfinite(tEnd):
            tEnd = 1.0
            while self.sva(tEnd)[0] < dv.max() and tEnd < 1024.0:
                tEnd *= 2.0

        # The segment containing each d, with breakpoints cut at the end of
        # the branch.
        tb = np.minimum(np.append(self._tb, tEnd), tEnd)
        k = np.clip(np.searchsorted(self.sva(tb[:3])[0], dv, side="right") - 1, 0, 2)
        s0, v0, a0 = self._state[k].T
        a0, j = a0/2, self._j[k]/6

        # Bisect within the segment; bracket width / 2**45 < 1e-13 s for the
        # < 1 s shutter segments.
        lo = np.zeros_like(dv)
        hi = tb[k + 1] - tb[k]
        for _ in range(45 + max(0, int(np.log2(max(hi.max(), 1.0))))):
            h = 0.5*(lo + hi)
            below = s0 + h*(v0 + h*(a0 + h*j)) < dv
            lo = np.where(below, h, lo)
            hi = np.where(below, hi, h)
        t = tb[k] + 0.5*(lo + hi)
        t[dv > self.sva(tEnd)[0]] = np.nan

        out[valid] = t
        return out


@dataclasses.dataclass(frozen=True)
class _DetectorGeometry:
    """Affine pixel -> DVCS (afw FOCAL_PLANE, mm) map of one detector:
    ``fp = centerMm + jacobian @ (p - centerPixel)``, with pixel centres at
    integer coordinates.
    """

    detectorId: int
    """The detector ID."""

    nx: int
    """Width (pixels)."""

    ny: int
    """Height (pixels)."""

    centerPixel: tuple[float, float]
    """Pixel coordinates of the detector centre."""

    centerMm: tuple[float, float]
    """DVCS coordinates (mm) of the detector centre."""

    jacobian: tuple[tuple[float, float], tuple[float, float]]
    """[[dX/dx, dX/dy], [dY/dx, dY/dy]] in mm per pixel."""

    @classmethod
    def fromDetector(cls, detector: Detector) -> "_DetectorGeometry":
        """Linearize a detector's pixel -> focal-plane map at its centre.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            The detector; for LSSTCam the map is affine to < 1e-3 mm.

        Returns
        -------
        geometry : ``_DetectorGeometry``
            The linearized map.
        """
        bbox = detector.getBBox()
        center = lsst.geom.Box2D(bbox).getCenter()
        transform = detector.getTransform(PIXELS, FOCAL_PLANE)
        fp = transform.applyForward(center)
        jac = transform.getJacobian(center)

        return cls(
            detectorId=detector.getId(),
            nx=bbox.getWidth(),
            ny=bbox.getHeight(),
            centerPixel=(center.getX(), center.getY()),
            centerMm=(fp.getX(), fp.getY()),
            jacobian=((jac[0, 0], jac[0, 1]), (jac[1, 0], jac[1, 1])),
        )

    def pixelToCcs(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Pixel -> CCS (mm): ``x_ccs = Y_dvcs``, ``y_ccs = X_dvcs``."""
        dx = x - self.centerPixel[0]
        dy = y - self.centerPixel[1]

        j = self.jacobian
        xCcs = self.centerMm[1] + j[1][0]*dx + j[1][1]*dy
        yCcs = self.centerMm[0] + j[0][0]*dx + j[0][1]*dy

        return xCcs, yCcs

    def bladeAxis(self) -> str:
        """The pixel axis (anti)parallel to DVCS y (the blade motion)."""
        j = self.jacobian
        return "x" if abs(j[1][0]) > abs(j[1][1]) else "y"

    def pixelCorners(self) -> np.ndarray:
        """The four pixel-edge corners, shape (4, 2)."""
        lo, hx, hy = -0.5, self.nx - 0.5, self.ny - 0.5
        return np.array([[lo, lo], [hx, lo], [hx, hy], [lo, hy]])


class _ShutterBeamModel:
    """Shutter-plane beam model: for each field position (CCS mm), the blade
    coordinate s_q (CCS-x mm) below which a fraction q of the pixel's beam
    flux lies, at the levels ``_BEAM_LEVELS``.

    The cumulative profile F(s) is piecewise linear through (s_q, q), with
    linear tails reaching F = 0 at ``s_0.01 - (s_0.05 - s_0.01)`` and F = 1 at
    ``s_0.99 + (s_0.99 - s_0.95)``.  ``R_q = s_q - x_ccs`` is interpolated
    bilinearly in grid cells with all four nodes tabulated, barycentrically on
    the Delaunay triangulation elsewhere inside the convex hull of the nodes
    (32 corner nodes are missing), and held at the nearest hull point outside
    it, up to ``_BEAM_MARGIN_MM`` (NaN beyond).

    Parameters
    ----------
    x, y : `numpy.ndarray`
        CCS coordinates (mm) of the table nodes, shape (n,).
    s : `numpy.ndarray`
        Blade coordinates s_q (mm) at the levels ``_BEAM_LEVELS``, shape
        (n, 9).

    Raises
    ------
    ValueError
        Raised if the quantiles are not strictly increasing in level.
    """

    def __init__(self, x: np.ndarray, y: np.ndarray, s: np.ndarray):
        if not np.all(np.diff(s, axis=1) > 0):
            raise ValueError("beam quantiles are not strictly increasing in level")

        # The nodes, and R_q on the regular grid they span (NaN where a node
        # is missing).
        self._points = np.column_stack([x, y])
        self._r = s - x[:, None]
        self._gx = np.unique(x)
        self._gy = np.unique(y)
        grid = np.full((self._gx.size, self._gy.size, _BEAM_LEVELS.size), np.nan)
        grid[np.searchsorted(self._gx, x), np.searchsorted(self._gy, y)] = self._r
        self._grid = grid

        # Grid cells with all four nodes, for bilinear interpolation.
        have = np.isfinite(grid[..., 0])
        self._cellOk = have[:-1, :-1] & have[1:, :-1] & have[:-1, 1:] & have[1:, 1:]

        # The triangulation and convex hull, for everywhere else.
        self._tri = LinearNDInterpolator(self._points, self._r)
        hull = ConvexHull(self._points)
        self._hullVertices = hull.vertices
        self._hullEq = hull.equations

        # Gauss-Legendre nodes and weights in flux fraction: two per interval
        # between levels, including the tails.
        knots = np.concatenate([[0.0], _BEAM_LEVELS, [1.0]])
        du = np.diff(knots)
        self._knots = knots
        self._nodes = (knots[:-1, None] + du[:, None]*_GL2[None, :]).ravel()
        self._weights = np.repeat(du/2.0, 2)

    @classmethod
    def fromFile(cls, path: str) -> "_ShutterBeamModel":
        """Read a beam table in the raytrace ``.tnt`` format.

        Parameters
        ----------
        path : `str`
            Path or URI of the table: a ``DATA`` line, a column-header line,
            then whitespace-separated columns x_ccs, y_ccs, level, (unused),
            s_q.

        Returns
        -------
        beam : ``_ShutterBeamModel``
            The model.

        Raises
        ------
        ValueError
            Raised if the file is not a beam table, its levels are not
            ``_BEAM_LEVELS``, or a (position, level) entry is missing or
            duplicated.
        """
        f = io.StringIO(ResourcePath(path).read().decode())
        if f.readline().strip().upper() != "DATA":
            raise ValueError(f"{path}: not a beam table (no DATA line)")
        f.readline()
        rows = np.loadtxt(f, ndmin=2)[:, [0, 1, 2, 4]]

        levels = np.unique(rows[:, 2])
        if levels.shape != _BEAM_LEVELS.shape or not np.allclose(levels, _BEAM_LEVELS):
            raise ValueError(f"{path}: unexpected levels {levels}")

        # One row per position, one column per level.
        pos, inv = np.unique(rows[:, :2], axis=0, return_inverse=True)
        s = np.full((pos.shape[0], levels.size), np.nan)
        s[inv.ravel(), np.searchsorted(levels, rows[:, 2])] = rows[:, 3]
        if np.isnan(s).any() or len(rows) != s.size:
            raise ValueError(f"{path}: incomplete or duplicated (position, level) table")

        return cls(pos[:, 0], pos[:, 1], s)

    def isOutside(self, p: np.ndarray) -> np.ndarray:
        """Whether each row of ``p`` (CCS mm) lies outside the convex hull of
        the table nodes.
        """
        return (p @ self._hullEq[:, :2].T + self._hullEq[:, 2]).max(axis=1) > 1e-9

    def _hullProject(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Nearest hull-boundary point of each row of p, and its distance."""
        # Each hull edge runs from a to a + ab.
        a = self._points[self._hullVertices]
        ab = np.roll(a, -1, axis=0) - a

        # The nearest point of every edge, and the nearest edge.
        ap = p[:, None, :] - a[None]
        t = np.clip((ap*ab).sum(-1)/(ab*ab).sum(-1), 0.0, 1.0)
        d2 = ((ap - t[..., None]*ab)**2).sum(-1)
        k = np.argmin(d2, axis=1)

        n = np.arange(p.shape[0])
        return a[k] + t[n, k, None]*ab[k], np.sqrt(d2[n, k])

    def _residual(self, x: np.ndarray, y: np.ndarray, extrapolate: bool = True) -> np.ndarray:
        """R_q = s_q - x_ccs at flat arrays x, y; shape (n, nlevels)."""
        out = np.full((x.size, _BEAM_LEVELS.size), np.nan)

        # Bilinear in grid cells with all four nodes.
        gx, gy = self._gx, self._gy
        ix = np.clip(np.searchsorted(gx, x, side="right") - 1, 0, gx.size - 2)
        iy = np.clip(np.searchsorted(gy, y, side="right") - 1, 0, gy.size - 2)
        inbox = (x >= gx[0]) & (x <= gx[-1]) & (y >= gy[0]) & (y <= gy[-1])
        bil = inbox & self._cellOk[ix, iy]
        if bil.any():
            i, j = ix[bil], iy[bil]
            tx = ((x[bil] - gx[i])/(gx[i + 1] - gx[i]))[:, None]
            ty = ((y[bil] - gy[j])/(gy[j + 1] - gy[j]))[:, None]
            g = self._grid
            out[bil] = (
                (1 - tx)*(1 - ty)*g[i, j]
                + tx*(1 - ty)*g[i + 1, j]
                + (1 - tx)*ty*g[i, j + 1]
                + tx*ty*g[i + 1, j + 1]
            )

        rest = ~bil
        if rest.any():
            # Barycentric on the triangulation inside the hull.
            p = np.column_stack([x[rest], y[rest]])
            inside = ~self.isOutside(p)
            r = np.full((p.shape[0], _BEAM_LEVELS.size), np.nan)
            if inside.any():
                r[inside] = self._tri(p[inside])

            # Outside the hull: the value just inside the nearest hull point,
            # up to the margin.
            need = ~np.isfinite(r[:, 0])
            if extrapolate and need.any():
                q, dist = self._hullProject(p[need])
                c = self._points.mean(axis=0)
                q += 1e-6*(c - q)/np.linalg.norm(c - q, axis=1, keepdims=True)
                rr = self._residual(q[:, 0], q[:, 1], extrapolate=False)
                rr[dist > _BEAM_MARGIN_MM] = np.nan
                r[need] = rr

            out[rest] = r

        return out

    def quadratureNodes(self, xCcs: np.ndarray, yCcs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Quadrature of the beam's cumulative flux profile.

        Parameters
        ----------
        xCcs, yCcs : `numpy.ndarray`
            Field positions (CCS mm), shape (n,).

        Returns
        -------
        s : `numpy.ndarray`
            Blade coordinates (mm) at the 20 Gauss nodes u_k of each
            position's profile, shape (n, 20).
        weights : `numpy.ndarray`
            The weights of the nodes, shape (20,).
        """
        # s_q at the tabulated levels, plus the linear tails at F = 0 and 1.
        s = self._residual(xCcs, yCcs) + xCcs[:, None]
        lo = s[:, :1] - (s[:, 1:2] - s[:, :1])
        hi = s[:, -1:] + (s[:, -1:] - s[:, -2:-1])
        ext = np.concatenate([lo, s, hi], axis=-1)

        # Piecewise-linear s(u) at the nodes.
        knots, u = self._knots, self._nodes
        k = np.clip(np.searchsorted(knots, u, side="right") - 1, 0, knots.size - 2)
        f = (u - knots[k])/(knots[k + 1] - knots[k])

        return ext[..., k]*(1 - f) + ext[..., k + 1]*f, self._weights


@functools.cache
def _loadShutterBeam(path: str) -> _ShutterBeamModel:
    """Read a shutter-plane beam table, cached per process.

    Parameters
    ----------
    path : `str`
        Path or URI of the table; the cache is keyed by this string.

    Returns
    -------
    beam : ``_ShutterBeamModel``
        The model.  It is shared between callers: do not modify it.
    """
    return _ShutterBeamModel.fromFile(path)


@dataclasses.dataclass(frozen=True)
class ShutterTiming:
    """Shutter-corrected mid-exposure times of one detector of one exposure.

    Times are MJD TAI, durations seconds.  Per source,
    ``t(x, y) = centerMidpointMjdTai + (c_u u + c_uu u^2 + c_v v + c_uv u v
    + c_vv v^2) / 86400``, with ``(u, v)`` the pixel offsets from the detector
    centre along and across the blade axis (``axis``: ``"x"`` -> u = x - cx,
    v = y - cy; ``"y"`` -> u = y - cy, v = x - cx).
    """

    status: ShutterTimingStatus = ShutterTimingStatus.UNAVAILABLE
    """OK, or UNAVAILABLE (callers keep the header midpoint)."""

    flags: ShutterTimingFlag = ShutterTimingFlag.NONE
    """Why the result is UNAVAILABLE, or notes on an OK result."""

    message: str = ""
    """Why the result is UNAVAILABLE (for logs and task metadata); "" if OK."""

    detectorId: int = -1
    """The detector ID."""

    centerMidpointMjdTai: float = math.nan
    """Flux-weighted mid-exposure time at the detector centre."""

    focalPlaneMidpointMjdTai: float = math.nan
    """Mid-exposure time at the focal-plane centre: the visit epoch."""

    headerMidpointMjdTai: float = math.nan
    """(MJD-BEG + MJD-END) / 2 (diagnostic only)."""

    effectiveExposureTime: float = math.nan
    """Flux-weighted open time at the detector centre (s)."""

    maxAbsResidual: float = math.nan
    """Largest |quadratic - exact| over the fit grid, its cell centres and the
    detector centre (s).
    """

    coefficients: tuple[float, float, float, float, float] = (math.nan,) * 5
    """(c_u, c_uu, c_v, c_uv, c_vv) in s / pixel^n."""

    axis: str = ""
    """The pixel axis along the blade motion, "x" or "y"; "" if UNAVAILABLE."""

    geometry: _DetectorGeometry | None = dataclasses.field(default=None, repr=False)
    """The detector's pixel-to-focal-plane map, used by `midpointMjdTai`."""

    def midpointMjdTai(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Per-source mid-exposure times at pixel positions.

        Parameters
        ----------
        x, y : `numpy.ndarray`
            Pixel coordinates; broadcast against each other.

        Returns
        -------
        tMid : `numpy.ndarray`
            Mid-exposure times (MJD TAI), shape ``np.broadcast(x, y).shape``.
            NaN everywhere if the detector is UNAVAILABLE, and NaN where x or
            y is NaN.
        """
        if self.status == ShutterTimingStatus.UNAVAILABLE:
            return np.full(np.broadcast(x, y).shape, np.nan)

        # Pixel offsets from the detector centre, along (u) and across (v)
        # the blade motion.
        cx, cy = self.geometry.centerPixel
        if self.axis == "x":
            u, v = x - cx, y - cy
        else:
            u, v = y - cy, x - cx

        # The quadratic, in seconds from the centre time.
        cU, cUU, cV, cUV, cVV = self.coefficients
        dt = cU*u + cUU*u*u + cV*v + cUV*u*v + cVV*v*v

        return self.centerMidpointMjdTai + dt/_SECONDS_PER_DAY


def computeShutterTiming(
    metadata: PropertyList | dict,
    detector: Detector,
    exposureTime: float,
    config: ShutterTimingConfig | None = None,
) -> ShutterTiming:
    """Shutter-corrected times of one detector of one exposure.

    Parameters
    ----------
    metadata : `lsst.daf.base.PropertyList`
        The exposure metadata (a `dict` works too).
    detector : `lsst.afw.cameraGeom.Detector`
        The detector.
    exposureTime : `float`
        The requested (nominal) exposure time (s).  Callers take it from the
        butler ``visit`` dimension record (``exposure_time``): ``EXPTIME`` is
        stripped from exposure metadata at ingest, and
        ``visitInfo.exposureTime`` is the measured shutter time, not the
        requested one.
    config : `ShutterTimingConfig`, optional
        Defaults to ``ShutterTimingConfig()``; ``beamFile`` must be set.

    Returns
    -------
    timing : `ShutterTiming`
        UNAVAILABLE (callers keep the header midpoint) if:

        - there are no shutter cards (``SHUTTER OPEN STARTTIME TAI MJD``
          absent; exposures before 2025-10-29): flag NO_PROFILE;
        - either Hall fit is outside the nominal ranges: flag PARAM_RANGE;
        - there is no time at the focal-plane centre, the detector centre or
          any fit-grid point (overlapping blades in very short exposures, or
          more than 38.64 mm outside the beam table).

        Otherwise OK.

    Raises
    ------
    ValueError
        If ``config.beamFile`` is empty, a SIDE card is not PLUSX or MINUSX,
        or a MODEL card is not ThreeJerksModelv1.
    KeyError, ValueError, TypeError
        If the shutter cards are present but one of them, or MJD-BEG or
        MJD-END, is missing or not a number.

    Notes
    -----
    Rules beyond the method in the module docstring:

    - Blade start positions are nominal: 750.76 mm for 750 -> 0 moves, -0.05
      mm for 0 -> 750 moves.
    - If close - open STARTTIME - ``exposureTime`` is outside -2 to 4 ms, the
      shutter clocks disagree: both start times are re-anchored to MJD-BEG
      and ``exposureTime`` (flag CLOCK_CLOSE_VS_OPEN).
    - Each edge's crossing time is integrated by 2-point Gauss-Legendre
      quadrature in each interval between tabulated flux levels.  The
      quadratic is fit on a 9 x 9 grid; ``maxAbsResidual`` is over that grid,
      its cell centres and the detector centre.
    """
    if config is None:
        config = ShutterTimingConfig()
    if not config.beamFile:
        raise ValueError("ShutterTimingConfig.beamFile is empty")
    beam = _loadShutterBeam(config.beamFile)

    # The tests pass a ``_DetectorGeometry`` directly.
    geometry = detector
    if not isinstance(detector, _DetectorGeometry):
        geometry = _DetectorGeometry.fromDetector(detector)
    common = dict(detectorId=geometry.detectorId, axis=geometry.bladeAxis(), geometry=geometry)

    # Exposures before the shutter cards were added.
    if "SHUTTER OPEN STARTTIME TAI MJD" not in metadata:
        return ShutterTiming(flags=ShutterTimingFlag.NO_PROFILE, message="no shutter motion cards", **common)

    common["headerMidpointMjdTai"] = 0.5*(float(metadata["MJD-BEG"]) + float(metadata["MJD-END"]))
    try:
        values = _compute(metadata, geometry, beam, exposureTime)
    except _Unavailable as e:
        return ShutterTiming(flags=e.flags, message=e.message, **common)

    return ShutterTiming(status=ShutterTimingStatus.OK, **values, **common)


class _Unavailable(Exception):
    """No corrected time for this detector.

    Parameters
    ----------
    flags : `ShutterTimingFlag`
        Why.
    message : `str`
        Why, for logs and task metadata.
    """

    def __init__(self, flags: ShutterTimingFlag, message: str):
        super().__init__(message)
        self.flags = flags
        self.message = message


def _inRange(value: float, bounds: tuple[float, float]) -> bool:
    """lo <= value <= hi; False for NaN."""
    return bounds[0] <= value <= bounds[1]


def _fitInRange(fit: _ThreeJerksParams) -> bool:
    """True if a Hall fit is within the nominal ranges."""
    displacementAt0p9s = _ThreeJerksTrajectory(fit, 0.0, -1).sva(0.9)[0]

    return (
        0.0 <= fit.pivot1 <= fit.pivot2
        and _inRange(fit.pivot1, _PIVOT1_RANGE)
        and _inRange(fit.pivot2, _PIVOT2_RANGE)
        and all(_inRange(abs(j), _ABS_JERK_RANGE) for j in (fit.jerk0, fit.jerk1, fit.jerk2))
        and _inRange(abs(fit.modelStartTime), _ABS_MODEL_START_TIME_RANGE)
        and _inRange(displacementAt0p9s, _DISPLACEMENT_AT_0P9S_RANGE)
    )


def _readMotion(metadata: PropertyList | dict, which: str) -> tuple[float, int, _ThreeJerksParams]:
    """Read the cards of one blade motion.

    Parameters
    ----------
    metadata : `lsst.daf.base.PropertyList` or `dict`
        The exposure metadata.
    which : `str`
        ``"OPEN"`` or ``"CLOSE"``.

    Returns
    -------
    startTime : `float`
        When the blade started moving (MJD TAI).
    travelSign : `int`
        -1 for a 750 -> 0 mm move (PLUSX-open, MINUSX-close), else +1.
    fit : ``_ThreeJerksParams``
        The Hall-sensor fit.

    Raises
    ------
    KeyError, ValueError, TypeError
        Raised if a card is missing or malformed.
    """
    pre = f"SHUTTER {which}"

    side = metadata[f"{pre} SIDE"]
    if side not in ("PLUSX", "MINUSX"):
        raise ValueError(f"'{pre} SIDE' is {side!r}, not PLUSX or MINUSX")
    if f"{pre} MODEL" in metadata and metadata[f"{pre} MODEL"] != _THREE_JERKS_V1:
        raise ValueError(f"'{pre} MODEL' is {metadata[f'{pre} MODEL']!r}, not {_THREE_JERKS_V1}")

    startTime = float(metadata[f"{pre} STARTTIME TAI MJD"])
    decreasing = (side == "PLUSX") == (which == "OPEN")
    fit = _ThreeJerksParams(*(float(metadata[f"{pre} HALLSENSORFIT {c}"]) for c in _FIT_CARDS))

    return startTime, -1 if decreasing else 1, fit


def _compute(
    metadata: PropertyList | dict,
    geometry: _DetectorGeometry,
    beam: _ShutterBeamModel,
    exposureTime: float,
) -> dict:
    """Compute the data-dependent fields of a `ShutterTiming`.

    Parameters
    ----------
    metadata : `lsst.daf.base.PropertyList` or `dict`
        The exposure metadata, with shutter cards.
    geometry : ``_DetectorGeometry``
        The detector.
    beam : ``_ShutterBeamModel``
        The beam model.
    exposureTime : `float`
        The requested exposure time (s).

    Returns
    -------
    values : `dict`
        ``flags``, ``centerMidpointMjdTai``, ``focalPlaneMidpointMjdTai``,
        ``effectiveExposureTime``, ``maxAbsResidual`` and ``coefficients``.

    Raises
    ------
    _Unavailable
        Raised if the detector has no corrected time.
    """
    motions = (_readMotion(metadata, "OPEN"), _readMotion(metadata, "CLOSE"))
    (openStart, _, _), (closeStart, _, _) = motions
    beg = float(metadata["MJD-BEG"])

    if not all(_fitInRange(fit) for _, _, fit in motions):
        raise _Unavailable(ShutterTimingFlag.PARAM_RANGE, "a Hall fit is outside the nominal ranges")

    # Zero point: re-anchor to the header if the two shutter clocks disagree.
    flags = ShutterTimingFlag.NONE
    closeRel = (closeStart - openStart)*_SECONDS_PER_DAY
    if not _inRange((closeRel - exposureTime)*1e3, _CLOSE_MINUS_OPEN_MINUS_EXPOSURE_TIME_RANGE):
        flags |= ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN
        openStart = beg + _HEADER_ANCHOR_OFFSETS[0]/_SECONDS_PER_DAY
        closeRel = exposureTime + _HEADER_ANCHOR_OFFSETS[1]

    # Beam coverage of the detector (its corners, by convexity).
    corners = geometry.pixelCorners()
    cornersCcs = np.column_stack(geometry.pixelToCcs(corners[:, 0], corners[:, 1]))
    if beam.isOutside(cornersCcs).any():
        flags |= ShutterTimingFlag.BEAM_EXTRAPOLATED

    # Evaluation points, along (pa) and across (pc) the blade axis: the
    # detector centre, the fit grid, and its cell centres.
    ia = 0 if geometry.bladeAxis() == "x" else 1
    c = geometry.centerPixel
    n = (geometry.nx, geometry.ny)
    along = np.linspace(0.0, n[ia] - 1.0, _GRID)
    across = np.linspace(0.0, n[1 - ia] - 1.0, _GRID)
    ga, gc = np.meshgrid(along, across, indexing="ij")
    ha, hc = np.meshgrid(0.5*(along[1:] + along[:-1]), 0.5*(across[1:] + across[:-1]), indexing="ij")
    pa = np.concatenate([[c[ia]], ga.ravel(), ha.ravel()])
    pc = np.concatenate([[c[1 - ia]], gc.ravel(), hc.ravel()])

    # The same points in CCS, plus the focal-plane centre.
    xc, yc = geometry.pixelToCcs(*((pa, pc) if ia == 0 else (pc, pa)))
    xc = np.append(xc, 0.0)
    yc = np.append(yc, 0.0)

    # Crossing times of each blade edge at the quadrature nodes (s after the
    # open start time).
    s, w = beam.quadratureNodes(xc, yc)
    tau = []
    for _, travelSign, fit in motions:
        start = _NOMINAL_STROKE if travelSign < 0 else _NOMINAL_START_INCREASING
        traj = _ThreeJerksTrajectory(fit, start, travelSign)
        d = travelSign*((_ENCODER_CENTER - start) - s)
        tau.append(traj.timeOfDisplacement(d) + fit.modelStartTime)
    tO, tauC = tau
    tC = tauC + closeRel

    # Flux-weighted times; none where the blades overlap in the beam.
    with np.errstate(invalid="ignore"):
        overlap = np.min(tC, axis=-1) - np.max(tO, axis=-1) <= 0
    eO = np.where(overlap, np.nan, tO @ w)
    eC = np.where(overlap, np.nan, tC @ w)
    tMid = 0.5*(eO + eC)
    if not np.all(np.isfinite(tMid)):
        raise _Unavailable(flags, f"no shutter time at {int((~np.isfinite(tMid)).sum())} of {tMid.size} "
                           "evaluation points (overlapping blades or outside the beam table's margin)")

    # Quadratic fit on the grid, in offsets from the centre time, with the
    # basis normalized to +-1 for conditioning.
    focalPlane, tMid = tMid[-1], tMid[:-1]
    dt = tMid - tMid[0]
    su = max(0.5*(n[ia] - 1.0), 1.0)
    sv = max(0.5*(n[1 - ia] - 1.0), 1.0)
    z = (pa - c[ia])/su
    wv = (pc - c[1 - ia])/sv
    basis = np.stack([z, z*z, wv, z*wv, wv*wv], axis=-1)
    unit = np.array([1/su, su**-2.0, 1/sv, 1/(su*sv), sv**-2.0])
    cz = np.linalg.lstsq(basis[1:1 + _GRID**2], dt[1:1 + _GRID**2], rcond=None)[0]

    return dict(
        flags=flags,
        centerMidpointMjdTai=float(openStart + tMid[0]/_SECONDS_PER_DAY),
        focalPlaneMidpointMjdTai=float(openStart + focalPlane/_SECONDS_PER_DAY),
        effectiveExposureTime=float(eC[0] - eO[0]),
        maxAbsResidual=float(np.max(np.abs(basis @ cz - dt))),
        coefficients=tuple(float(v) for v in cz*unit),
    )
