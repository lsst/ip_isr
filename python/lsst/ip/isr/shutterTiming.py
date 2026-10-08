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

from __future__ import annotations

__all__ = [
    "ShutterTiming",
    "ShutterTimingConfig",
    "ShutterTimingFlag",
    "ShutterTimingStatus",
    "computeShutterTiming",
    "loadShutterBeam",
]

import dataclasses
import enum
import functools
import io
import math
import numbers

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import ConvexHull

import lsst.pex.config as pexConfig

from ._shutterTrajectory import ThreeJerksParams, ThreeJerksTrajectory

_SECONDS_PER_DAY = 86400.0

#: Cumulative flux levels q of the beam table.
_BEAM_LEVELS = np.array([0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99])

#: Extrapolation margin (mm) outside the beam table's convex hull: the worst
#: science pixel (28.64 mm outside) plus 10 mm; NaN beyond.
_BEAM_MARGIN_MM = 28.64 + 10.0

#: 2-point Gauss-Legendre abscissae on [0, 1].
_GL2 = np.array([0.5 - 0.5 / np.sqrt(3.0), 0.5 + 0.5 / np.sqrt(3.0)])

_THREE_JERKS_V1 = "ThreeJerksModelv1"

#: Card suffixes of the ThreeJerksParams fields, in order.
_FIT_CARDS = ("MODELSTARTTIME", "PIVOTPOINT1", "PIVOTPOINT2", "JERK0", "JERK1", "JERK2")

#: Encoder position (mm) of the focal-plane centre: x_ccs = c - e.
_ENCODER_CENTER = 375.0

#: Start position (mm) of a 750 -> 0 move (the cards carry none).
_NOMINAL_STROKE = 750.76

#: Start position (mm) of a 0 -> 750 move.
_NOMINAL_START_INCREASING = -0.05

#: Mean Hall fits per travel sign, used when a fit is missing or out of range.
_MEAN_PROFILES = {
    -1: ThreeJerksParams(0.0007497095570137927, 0.22426363506685815, 0.6776782303465883,
                         33100.38075466996, -32888.35381319897, 34383.386427176374),
    1: ThreeJerksParams(0.0009719394948759783, 0.2237721807772504, 0.6788008740104505,
                        33178.9312914494, -32841.22646579422, 35458.25981717702),
}

#: Re-anchoring offsets (s): open STARTTIME - MJD-BEG, and close - open
#: STARTTIME - EXPTIME.
_HEADER_ANCHOR_OFFSETS = (0.00809, 0.00052)

#: Fit grid points along and across the blade axis.
_GRID = 9


class ShutterTimingStatus(enum.IntEnum):
    """Quality of a detector's corrected time."""

    OK = 0
    """Corrected; accurate to ~1 ms."""
    DEGRADED = 1
    """Corrected with reduced accuracy (see ``flags``); still better than the
    header midpoint.
    """
    UNAVAILABLE = 2
    """No corrected time: callers keep the header midpoint."""


class ShutterTimingFlag(enum.IntFlag):
    """Reasons for a DEGRADED or UNAVAILABLE result (fixed bit values)."""

    NONE = 0
    NO_PROFILE = 1
    """Missing or unusable Hall-fit cards: the mean profile is used.
    UNAVAILABLE if the start time or side cards are unusable too.
    """
    PARAM_RANGE = 2
    """A Hall fit outside the nominal ranges (the mean profile is used)."""
    CLOCK_OPEN_VS_BEG = 4
    """Open STARTTIME - MJD-BEG outside ``openMinusBegRange`` (flag only)."""
    CLOCK_END_VS_CLOSE = 8
    """MJD-END - close STARTTIME outside ``endMinusCloseRange`` (a late
    readout; flag only).
    """
    CLOCK_CLOSE_VS_OPEN = 16
    """Close - open STARTTIME - EXPTIME outside
    ``closeMinusOpenMinusExptimeRange``: both start times are re-anchored to
    MJD-BEG and EXPTIME.
    """
    PRE_CLOCK_EPOCH = 32
    """Zero point re-anchored to the header (with CLOCK_CLOSE_VS_OPEN)."""
    MEAN_PROFILE = 64
    """A mean profile replaced a missing or out-of-range Hall fit."""
    SHUTTIME_MISMATCH = 1024
    """|T_eff(focal-plane centre) - SHUTTIME| > ``shuttimeTolerance``
    (diagnostic only).
    """
    BEAM_EXTRAPOLATED = 2048
    """Part of the detector lies outside the beam table's coverage
    (diagnostic only).
    """


#: Flags that make a corrected time DEGRADED.
_DEGRADED_FLAGS = (
    ShutterTimingFlag.NO_PROFILE | ShutterTimingFlag.PARAM_RANGE | ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN
    | ShutterTimingFlag.PRE_CLOCK_EPOCH | ShutterTimingFlag.MEAN_PROFILE
)


class ShutterTimingConfig(pexConfig.Config):
    """Configuration of `computeShutterTiming` (LSSTCam defaults; the obs
    package sets ``beamFile``).
    """

    beamFile = pexConfig.Field(
        dtype=str, default="",
        doc="Path or URI of the shutter-plane beam table (raytrace flux quantiles, LCA-20578 format).",
    )
    pivot1Range = pexConfig.ListField(
        dtype=float, default=[0.20, 0.25], length=2, doc="PivotPoint1 range (s).")
    pivot2Range = pexConfig.ListField(
        dtype=float, default=[0.65, 0.70], length=2, doc="PivotPoint2 range (s).")
    absJerkRange = pexConfig.ListField(
        dtype=float, default=[30000.0, 38000.0], length=2, doc="|Jerk0,1,2| range (mm/s^3).")
    absModelStartTimeRange = pexConfig.ListField(
        dtype=float, default=[0.0, 0.003], length=2, doc="|ModelStartTime| range (s).")
    displacementAt0p9sRange = pexConfig.ListField(
        dtype=float, default=[745.0, 757.0], length=2, doc="Model displacement at t = 0.9 s range (mm).")
    openMinusBegRange = pexConfig.ListField(
        dtype=float, default=[5.0, 15.0], length=2, doc="Open STARTTIME - MJD-BEG range (ms).")
    endMinusCloseRange = pexConfig.ListField(
        dtype=float, default=[0.89, 0.96], length=2, doc="MJD-END - close STARTTIME range (s).")
    closeMinusOpenMinusExptimeRange = pexConfig.ListField(
        dtype=float, default=[-2.0, 4.0], length=2, doc="Close - open STARTTIME - EXPTIME range (ms).")
    shuttimeTolerance = pexConfig.Field(
        dtype=float, default=1.0, doc="|T_eff(centre) - SHUTTIME| tolerance (ms).")
    degradedResidual = pexConfig.Field(
        dtype=float, default=1e-3, doc="Quadratic fit residual (s) above which the detector is DEGRADED.")
    offDetectorLimit = pexConfig.Field(
        dtype=float, default=100.0,
        doc="Per-source times are given up to this far (pixels) outside the detector; NaN beyond.",
    )


@dataclasses.dataclass(frozen=True)
class _DetectorGeometry:
    """Affine pixel -> DVCS (afw FOCAL_PLANE, mm) map of one detector:
    ``fp = centerMm + jacobian @ (p - centerPixel)``, with pixel centres at
    integer coordinates.
    """

    detectorId: int
    nx: int
    ny: int
    centerPixel: tuple[float, float]
    centerMm: tuple[float, float]
    jacobian: tuple[tuple[float, float], tuple[float, float]]
    """[[dX/dx, dX/dy], [dY/dx, dY/dy]] in mm per pixel."""
    isScience: bool = True

    @classmethod
    def fromDetector(cls, detector) -> _DetectorGeometry:
        """From an `lsst.afw.cameraGeom.Detector`, linearized at its centre
        (the LSSTCam map is affine to < 1e-3 mm).
        """
        import lsst.geom
        from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS, DetectorType

        bbox = detector.getBBox()
        center = lsst.geom.Box2D(bbox).getCenter()
        transform = detector.getTransform(PIXELS, FOCAL_PLANE)
        fp = transform.applyForward(center)
        jac = np.asarray(transform.getJacobian(center), dtype=float)
        return cls(
            detectorId=int(detector.getId()),
            nx=int(bbox.getWidth()),
            ny=int(bbox.getHeight()),
            centerPixel=(float(center.getX()), float(center.getY())),
            centerMm=(float(fp.getX()), float(fp.getY())),
            jacobian=((float(jac[0, 0]), float(jac[0, 1])), (float(jac[1, 0]), float(jac[1, 1]))),
            isScience=detector.getType() == DetectorType.SCIENCE,
        )

    def pixelToCcs(self, x, y):
        """Pixel -> CCS (mm): ``x_ccs = Y_dvcs``, ``y_ccs = X_dvcs``."""
        dx = np.asarray(x, dtype=float) - self.centerPixel[0]
        dy = np.asarray(y, dtype=float) - self.centerPixel[1]
        j = self.jacobian
        return (self.centerMm[1] + j[1][0] * dx + j[1][1] * dy,
                self.centerMm[0] + j[0][0] * dx + j[0][1] * dy)

    def offDetector(self, x, y):
        """Distance (pixels, the larger axis) outside the pixel bounds; 0
        inside.
        """
        x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
        return np.maximum.reduce([np.zeros(x.shape), -0.5 - x, x - (self.nx - 0.5),
                                  -0.5 - y, y - (self.ny - 0.5)])

    def bladeAxis(self):
        """The pixel axis (anti)parallel to DVCS y (the blade motion)."""
        j = self.jacobian
        return "x" if abs(j[1][0]) > abs(j[1][1]) else "y"

    def pixelCorners(self):
        """The four pixel-edge corners (4, 2)."""
        lo, hx, hy = -0.5, self.nx - 0.5, self.ny - 0.5
        return np.array([[lo, lo], [hx, lo], [hx, hy], [lo, hy]])


class _ShutterBeamModel:
    """Shutter-plane beam model: for each field position (CCS mm), the blade
    coordinate s_q (CCS-x mm) below which a fraction q of the pixel's beam
    flux lies, at the levels `_BEAM_LEVELS`.

    The cumulative profile F(s) is piecewise linear through (s_q, q), with
    linear tails reaching F = 0 at ``s_0.01 - (s_0.05 - s_0.01)`` and F = 1 at
    ``s_0.99 + (s_0.99 - s_0.95)``.  ``R_q = s_q - x_ccs`` is interpolated
    bilinearly in grid cells with all four nodes tabulated, barycentrically on
    the Delaunay triangulation elsewhere inside the convex hull of the nodes
    (32 corner nodes are missing), and held at the nearest hull point outside
    it, up to `_BEAM_MARGIN_MM` (NaN beyond).
    """

    def __init__(self, xCcs, yCcs, sQ):
        x = np.asarray(xCcs, dtype=float)
        y = np.asarray(yCcs, dtype=float)
        s = np.asarray(sQ, dtype=float)
        if not np.all(np.diff(s, axis=1) > 0):
            raise ValueError("beam quantiles are not strictly increasing in level")
        self._points = np.column_stack([x, y])
        self._r = s - x[:, None]
        self._gx = np.unique(x)
        self._gy = np.unique(y)
        grid = np.full((self._gx.size, self._gy.size, _BEAM_LEVELS.size), np.nan)
        grid[np.searchsorted(self._gx, x), np.searchsorted(self._gy, y)] = self._r
        self._grid = grid
        have = np.isfinite(grid[..., 0])
        self._cellOk = have[:-1, :-1] & have[1:, :-1] & have[:-1, 1:] & have[1:, 1:]
        self._tri = LinearNDInterpolator(self._points, self._r)
        hull = ConvexHull(self._points)
        self._hullVertices = hull.vertices
        self._hullEq = hull.equations
        knots = np.concatenate([[0.0], _BEAM_LEVELS, [1.0]])
        du = np.diff(knots)
        self._knots = knots
        self._nodes = (knots[:-1, None] + du[:, None] * _GL2[None, :]).ravel()
        self._weights = np.repeat(du / 2.0, 2)

    @classmethod
    def fromFile(cls, path):
        """Read a beam table in the raytrace ``.tnt`` format (a ``DATA``
        line, a column-header line, then whitespace-separated columns x_ccs,
        y_ccs, level, (unused), s_q).
        """
        from lsst.resources import ResourcePath

        f = io.StringIO(ResourcePath(path).read().decode())
        if f.readline().strip().upper() != "DATA":
            raise ValueError(f"{path}: not a beam table (no DATA line)")
        f.readline()
        rows = np.loadtxt(f, ndmin=2)[:, [0, 1, 2, 4]]
        levels = np.unique(rows[:, 2])
        if levels.shape != _BEAM_LEVELS.shape or not np.allclose(levels, _BEAM_LEVELS):
            raise ValueError(f"{path}: unexpected levels {levels}")
        pos, inv = np.unique(rows[:, :2], axis=0, return_inverse=True)
        s = np.full((pos.shape[0], levels.size), np.nan)
        s[inv.ravel(), np.searchsorted(levels, rows[:, 2])] = rows[:, 3]
        if np.isnan(s).any() or len(rows) != s.size:
            raise ValueError(f"{path}: incomplete or duplicated (position, level) table")
        return cls(pos[:, 0], pos[:, 1], s)

    def isOutside(self, p):
        """Whether each row of ``p`` (CCS mm) lies outside the convex hull of
        the table nodes.
        """
        return (p @ self._hullEq[:, :2].T + self._hullEq[:, 2]).max(axis=1) > 1e-9

    def _hullProject(self, p):
        """Nearest hull-boundary point of each row of p, and its distance."""
        a = self._points[self._hullVertices]
        ab = np.roll(a, -1, axis=0) - a
        ap = p[:, None, :] - a[None]
        t = np.clip((ap * ab).sum(-1) / (ab * ab).sum(-1), 0.0, 1.0)
        d2 = ((ap - t[..., None] * ab) ** 2).sum(-1)
        k = np.argmin(d2, axis=1)
        n = np.arange(p.shape[0])
        return a[k] + t[n, k, None] * ab[k], np.sqrt(d2[n, k])

    def _residual(self, x, y, extrapolate=True):
        """R_q = s_q - x_ccs at flat arrays x, y; shape (n, nlevels)."""
        out = np.full((x.size, _BEAM_LEVELS.size), np.nan)
        gx, gy = self._gx, self._gy
        ix = np.clip(np.searchsorted(gx, x, side="right") - 1, 0, gx.size - 2)
        iy = np.clip(np.searchsorted(gy, y, side="right") - 1, 0, gy.size - 2)
        inbox = (x >= gx[0]) & (x <= gx[-1]) & (y >= gy[0]) & (y <= gy[-1])
        bil = inbox & self._cellOk[ix, iy]
        if bil.any():
            i, j = ix[bil], iy[bil]
            tx = ((x[bil] - gx[i]) / (gx[i + 1] - gx[i]))[:, None]
            ty = ((y[bil] - gy[j]) / (gy[j + 1] - gy[j]))[:, None]
            g = self._grid
            out[bil] = (
                (1 - tx) * (1 - ty) * g[i, j]
                + tx * (1 - ty) * g[i + 1, j]
                + (1 - tx) * ty * g[i, j + 1]
                + tx * ty * g[i + 1, j + 1]
            )
        rest = ~bil & np.isfinite(x) & np.isfinite(y)
        if rest.any():
            p = np.column_stack([x[rest], y[rest]])
            inside = ~self.isOutside(p)
            r = np.full((p.shape[0], _BEAM_LEVELS.size), np.nan)
            if inside.any():
                r[inside] = self._tri(p[inside])
            need = ~np.isfinite(r[:, 0])
            if extrapolate and need.any():
                q, dist = self._hullProject(p[need])
                c = self._points.mean(axis=0)
                q += 1e-6 * (c - q) / np.linalg.norm(c - q, axis=1, keepdims=True)
                rr = self._residual(q[:, 0], q[:, 1], extrapolate=False)
                rr[dist > _BEAM_MARGIN_MM] = np.nan
                r[need] = rr
            out[rest] = r
        return out

    def quadratureNodes(self, xCcs, yCcs):
        """Blade coordinates s(u_k) at the 20 Gauss nodes u_k of the beam CDF,
        shape ``(n, 20)`` for flat ``xCcs``, ``yCcs``; and the weights.
        """
        s = self._residual(xCcs, yCcs) + xCcs[:, None]
        lo = s[:, :1] - (s[:, 1:2] - s[:, :1])
        hi = s[:, -1:] + (s[:, -1:] - s[:, -2:-1])
        ext = np.concatenate([lo, s, hi], axis=-1)
        knots, u = self._knots, self._nodes
        k = np.clip(np.searchsorted(knots, u, side="right") - 1, 0, knots.size - 2)
        f = (u - knots[k]) / (knots[k + 1] - knots[k])
        return ext[..., k] * (1 - f) + ext[..., k + 1] * f, self._weights


@functools.cache
def loadShutterBeam(path: str) -> _ShutterBeamModel:
    """Read the beam table at ``path`` (a path or URI), cached per process by
    the ``path`` string.  The returned model is shared: do not modify it.
    """
    return _ShutterBeamModel.fromFile(path)


@dataclasses.dataclass(frozen=True)
class ShutterTiming:
    """Shutter-corrected mid-exposure times of one detector of one exposure.

    Times are MJD TAI, durations seconds.  Per source,
    ``t(x, y) = centerMjdTai + (c_u u + c_uu u^2 + c_v v + c_uv u v + c_vv
    v^2) / 86400``, with ``(u, v)`` the pixel offsets from the detector
    centre along and across the blade axis (``axis``: ``"x"`` -> u = x - cx,
    v = y - cy; ``"y"`` -> u = y - cy, v = x - cx).

    A DEGRADED detector may have a centre time but no quadratic (NaN
    ``coefficients`` and ``maxAbsResidual``) when no time exists at some grid
    points, e.g. where the blades overlap in exposures shorter than ~0.2 s.
    """

    status: ShutterTimingStatus = ShutterTimingStatus.UNAVAILABLE
    """UNAVAILABLE if there is no time at the detector or focal-plane centre;
    DEGRADED if ``flags`` has a degrading flag or not ``maxAbsResidual <=
    degradedResidual``; else OK.
    """
    flags: ShutterTimingFlag = ShutterTimingFlag.NONE
    message: str = ""
    """Why the result is not OK (for logs and task metadata); "" if OK."""
    detectorId: int = -1
    centerMjdTai: float = math.nan
    """Flux-weighted mid-exposure time at the detector centre."""
    focalPlaneMjdTai: float = math.nan
    """Mid-exposure time at the focal-plane centre: the visit epoch."""
    headerMidMjdTai: float = math.nan
    """(MJD-BEG + MJD-END) / 2, NaN if missing (diagnostic only)."""
    effectiveExposureTime: float = math.nan
    """Flux-weighted open time at the detector centre (s)."""
    maxAbsResidual: float = math.nan
    """Largest |quadratic - exact| over the fit grid (s)."""
    coefficients: tuple[float, float, float, float, float] = (math.nan,) * 5
    """(c_u, c_uu, c_v, c_uv, c_vv) in s / pixel^n."""
    axis: str = ""
    geometry: _DetectorGeometry | None = dataclasses.field(default=None, repr=False)
    offDetectorLimit: float = dataclasses.field(default=100.0, repr=False)

    def tMidMjdTai(self, x, y) -> np.ndarray:
        """Per-source mid-exposure times (MJD TAI) at pixel positions, shape
        ``np.broadcast(x, y)``.  NaN where the detector is UNAVAILABLE or has
        no quadratic, x or y is not finite, or the position is more than
        ``offDetectorLimit`` pixels outside the detector.
        """
        x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
        t = np.full(x.shape, np.nan)
        if self.status == ShutterTimingStatus.UNAVAILABLE or self.geometry is None:
            return t
        good = self.geometry.offDetector(x, y) <= self.offDetectorLimit
        cx, cy = self.geometry.centerPixel
        xs, ys = x[good] - cx, y[good] - cy
        u, v = (xs, ys) if self.axis == "x" else (ys, xs)
        cU, cUU, cV, cUV, cVV = self.coefficients
        ds = cU * u + cUU * u * u + cV * v + cUV * u * v + cVV * v * v
        t[good] = self.centerMjdTai + ds / _SECONDS_PER_DAY
        return t

    def summary(self) -> dict:
        """Scalars for task metadata."""
        return dict(
            status=self.status.name,
            flags=int(self.flags),
            message=self.message,
            centerMjdTai=float(self.centerMjdTai),
            focalPlaneMjdTai=float(self.focalPlaneMjdTai),
            centerMinusHeaderMid=float((self.centerMjdTai - self.headerMidMjdTai) * _SECONDS_PER_DAY),
            maxAbsResidual=float(self.maxAbsResidual),
        )


def computeShutterTiming(metadata, detector, config: ShutterTimingConfig | None = None) -> ShutterTiming:
    """Shutter-corrected times of one detector of one exposure.

    Parameters
    ----------
    metadata : `lsst.daf.base.PropertyList`
        The exposure metadata (a `dict` works too).
    detector : `lsst.afw.cameraGeom.Detector`
        The detector.
    config : `ShutterTimingConfig`, optional
        Defaults to ``ShutterTimingConfig()``; ``beamFile`` must be set.

    Returns
    -------
    timing : `ShutterTiming`
        Problems in the data (missing or non-numeric cards, failed fits,
        overlapping blades, a non-science detector) give UNAVAILABLE or
        DEGRADED with ``flags`` and ``message``.

    Raises
    ------
    ValueError
        If ``config.beamFile`` is empty.

    Notes
    -----
    Rules beyond the method in the module docstring:

    - A card value that is not a finite number counts as missing.  Missing
      start time or side cards give UNAVAILABLE (NO_PROFILE).  A missing Hall
      fit, an unsupported ``MODEL`` or a fit outside the configured ranges is
      replaced by the mean profile of its travel direction (MEAN_PROFILE).
    - Blade start positions are nominal: 750.76 mm for 750 -> 0 moves, -0.05
      mm for 0 -> 750 moves.
    - CLOCK_CLOSE_VS_OPEN re-anchors both start times to MJD-BEG and EXPTIME
      (PRE_CLOCK_EPOCH); without them the result is UNAVAILABLE.
    - Each edge's crossing time is integrated by 2-point Gauss-Legendre
      quadrature in each interval between tabulated flux levels.  The
      quadratic is fit on a 9 x 9 grid; ``maxAbsResidual`` is over that grid,
      its cell centres and the detector centre.
    """
    if config is None:
        config = ShutterTimingConfig()
    if not config.beamFile:
        raise ValueError("ShutterTimingConfig.beamFile is empty")
    beam = loadShutterBeam(config.beamFile)
    geometry = detector
    if not isinstance(detector, _DetectorGeometry):
        geometry = _DetectorGeometry.fromDetector(detector)
    beg, end = _number(metadata, "MJD-BEG"), _number(metadata, "MJD-END")
    common = dict(
        detectorId=geometry.detectorId, axis=geometry.bladeAxis(), geometry=geometry,
        headerMidMjdTai=0.5 * (beg + end) if beg is not None and end is not None else math.nan,
        offDetectorLimit=config.offDetectorLimit,
    )
    try:
        values, noQuadratic = _compute(metadata, geometry, config, beam)
    except _Unavailable as e:
        return ShutterTiming(flags=e.flags, message=e.message, **common)
    flags = values["flags"]
    reasons = []
    if flags & _DEGRADED_FLAGS:
        reasons.append("|".join(f.name for f in ShutterTimingFlag if f & flags & _DEGRADED_FLAGS))
    if noQuadratic:
        reasons.append(noQuadratic)
    elif not values["maxAbsResidual"] <= config.degradedResidual:
        reasons.append(f"quadratic residual {values['maxAbsResidual'] * 1e3:.3f} ms > "
                       f"{config.degradedResidual * 1e3:g} ms")
    status = ShutterTimingStatus.DEGRADED if reasons else ShutterTimingStatus.OK
    return ShutterTiming(status=status, message="; ".join(reasons), **values, **common)


class _Unavailable(Exception):
    """No corrected time for this detector."""

    def __init__(self, flags, message):
        super().__init__(message)
        self.flags = flags
        self.message = message


def _number(metadata, key):
    """The card's value as a float, or None unless it is a finite real
    number (not a bool).
    """
    value = metadata.get(key)
    if not isinstance(value, numbers.Real) or isinstance(value, (bool, np.bool_)):
        return None
    try:
        value = float(value)
    except OverflowError:
        return None
    return value if math.isfinite(value) else None


def _inRange(value, bounds):
    lo, hi = bounds
    return bool(np.isfinite(value) and lo <= value <= hi)


def _fitInRange(fit, config):
    """True if a Hall fit is within the configured nominal ranges."""
    return (
        0.0 <= fit.pivot1 <= fit.pivot2
        and _inRange(fit.pivot1, config.pivot1Range)
        and _inRange(fit.pivot2, config.pivot2Range)
        and all(_inRange(abs(j), config.absJerkRange) for j in (fit.jerk0, fit.jerk1, fit.jerk2))
        and _inRange(abs(fit.modelStartTime), config.absModelStartTimeRange)
        and _inRange(float(ThreeJerksTrajectory(fit, 0.0, -1).sva(0.9)[0]), config.displacementAt0p9sRange)
    )


@dataclasses.dataclass
class _Motion:
    """One blade motion from the cards."""

    startMjdTai: float
    travelSign: int
    fit: ThreeJerksParams | None
    """The Hall fit; None if a card is missing."""
    usable: bool
    """A Hall fit of the supported model is present."""


def _readMotion(metadata, which):
    pre = f"SHUTTER {which}"
    start = _number(metadata, f"{pre} STARTTIME TAI MJD")
    side = metadata.get(f"{pre} SIDE")
    if start is None or side not in ("PLUSX", "MINUSX"):
        raise _Unavailable(ShutterTimingFlag.NO_PROFILE,
                           f"no usable '{pre} STARTTIME TAI MJD' / '{pre} SIDE' cards")
    values = [_number(metadata, f"{pre} HALLSENSORFIT {c}") for c in _FIT_CARDS]
    fit = None if None in values else ThreeJerksParams(*values)
    model = metadata.get(f"{pre} MODEL")
    decreasing = (side == "PLUSX") == (which == "OPEN")  # PLUSX-open, MINUSX-close: 750 -> 0
    return _Motion(startMjdTai=start, travelSign=-1 if decreasing else 1, fit=fit,
                   usable=fit is not None and (model is None or model == _THREE_JERKS_V1))


def _compute(metadata, geometry, config, beam):
    """The `ShutterTiming` fields that depend on the data, and why there is
    no quadratic ("" if there is one); raises `_Unavailable`.
    """
    if not geometry.isScience:
        raise _Unavailable(ShutterTimingFlag.NONE,
                           f"detector {geometry.detectorId} is not SCIENCE: the beam table does not cover it")
    motions = (_readMotion(metadata, "OPEN"), _readMotion(metadata, "CLOSE"))
    mOpen, mClose = motions
    exptime, shuttime, beg, end = (_number(metadata, k)
                                   for k in ("EXPTIME", "SHUTTIME", "MJD-BEG", "MJD-END"))

    # Exposure QC.
    flags = ShutterTimingFlag.NONE
    if not (mOpen.usable and mClose.usable):
        flags |= ShutterTimingFlag.NO_PROFILE
    fitOk = []
    for m in motions:
        inRange = m.fit is not None and _fitInRange(m.fit, config)
        if m.fit is not None and not inRange:
            flags |= ShutterTimingFlag.PARAM_RANGE
        fitOk.append(m.usable and inRange)
    if beg is not None:
        v = (mOpen.startMjdTai - beg) * _SECONDS_PER_DAY * 1e3
        if not _inRange(v, config.openMinusBegRange):
            flags |= ShutterTimingFlag.CLOCK_OPEN_VS_BEG
    if end is not None:
        v = (end - mClose.startMjdTai) * _SECONDS_PER_DAY
        if not _inRange(v, config.endMinusCloseRange):
            flags |= ShutterTimingFlag.CLOCK_END_VS_CLOSE
    if exptime is not None:
        v = ((mClose.startMjdTai - mOpen.startMjdTai) * _SECONDS_PER_DAY - exptime) * 1e3
        if not _inRange(v, config.closeMinusOpenMinusExptimeRange):
            flags |= ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN

    # Missing, unusable or out-of-range Hall fits -> the mean profile.
    for m, ok in zip(motions, fitOk):
        if not ok:
            m.fit = _MEAN_PROFILES[m.travelSign]
            flags |= ShutterTimingFlag.MEAN_PROFILE

    # Zero point: re-anchor to the header if the two shutter clocks disagree.
    closeRelProfile = (mClose.startMjdTai - mOpen.startMjdTai) * _SECONDS_PER_DAY
    if flags & ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN:
        if beg is None or exptime is None:
            raise _Unavailable(flags, "CLOCK_CLOSE_VS_OPEN, and re-anchoring needs MJD-BEG and EXPTIME")
        openStart = beg + _HEADER_ANCHOR_OFFSETS[0] / _SECONDS_PER_DAY
        closeRel = exptime + _HEADER_ANCHOR_OFFSETS[1]
        flags |= ShutterTimingFlag.PRE_CLOCK_EPOCH
    else:
        openStart = mOpen.startMjdTai
        closeRel = closeRelProfile

    # Evaluation points: detector centre, fit grid, its cell centres, and the
    # focal-plane centre.
    ia = 0 if geometry.bladeAxis() == "x" else 1
    c = geometry.centerPixel
    n = (geometry.nx, geometry.ny)
    along = np.linspace(0.0, n[ia] - 1.0, _GRID)
    across = np.linspace(0.0, n[1 - ia] - 1.0, _GRID)
    ga, gc = np.meshgrid(along, across, indexing="ij")
    ha, hc = np.meshgrid(0.5 * (along[1:] + along[:-1]), 0.5 * (across[1:] + across[:-1]), indexing="ij")
    pa = np.concatenate([[c[ia]], ga.ravel(), ha.ravel()])
    pc = np.concatenate([[c[1 - ia]], gc.ravel(), hc.ravel()])
    xc, yc = geometry.pixelToCcs(*((pa, pc) if ia == 0 else (pc, pa)))
    xc = np.append(xc, 0.0)
    yc = np.append(yc, 0.0)

    # Flux-weighted crossing times (s after the open start time).
    s, w = beam.quadratureNodes(xc, yc)
    tau = []
    for m in motions:
        start = _NOMINAL_STROKE if m.travelSign < 0 else _NOMINAL_START_INCREASING
        traj = ThreeJerksTrajectory(m.fit, start, m.travelSign)
        d = m.travelSign * ((_ENCODER_CENTER - start) - s)
        tau.append(traj.timeOfDisplacement(d) + m.fit.modelStartTime)
    tO, tauC = tau
    tC = tauC + closeRel
    with np.errstate(invalid="ignore"):
        overlap = np.min(tC, axis=-1) - np.max(tO, axis=-1) <= 0
    eO = np.where(overlap, np.nan, tO @ w)
    eC = np.where(overlap, np.nan, tC @ w)

    # SHUTTIME diagnostic (profile zero point, focal-plane centre).
    if shuttime is not None and not (flags & (ShutterTimingFlag.NO_PROFILE | ShutterTimingFlag.PARAM_RANGE)):
        tCp = tauC[-1] + closeRelProfile
        tEff = (tCp @ w) - (tO[-1] @ w) if np.min(tCp) > np.max(tO[-1]) else math.nan
        if not abs((tEff - shuttime) * 1e3) <= config.shuttimeTolerance:
            flags |= ShutterTimingFlag.SHUTTIME_MISMATCH

    # Beam coverage of the detector (its corners, by convexity).
    corners = geometry.pixelCorners()
    if beam.isOutside(np.column_stack(geometry.pixelToCcs(corners[:, 0], corners[:, 1]))).any():
        flags |= ShutterTimingFlag.BEAM_EXTRAPOLATED

    focalPlane = float(openStart + 0.5 * (eO[-1] + eC[-1]) / _SECONDS_PER_DAY)
    eO, eC = eO[:-1], eC[:-1]
    tMid = 0.5 * (eO + eC)
    why = "(outside the beam table's margin, off the trajectory branch, or overlapping blades)"
    if not math.isfinite(focalPlane):
        raise _Unavailable(flags, f"no shutter time at the focal-plane centre {why}")
    if not math.isfinite(tMid[0]):
        raise _Unavailable(flags, f"no shutter time at the detector centre {why}")
    values = dict(flags=flags, centerMjdTai=float(openStart + tMid[0] / _SECONDS_PER_DAY),
                  focalPlaneMjdTai=focalPlane, effectiveExposureTime=float(eC[0] - eO[0]))

    # Quadratic fit, only with a time at every grid point.
    if not np.all(np.isfinite(tMid)):
        nBad = int((~np.isfinite(tMid)).sum())
        return values, (f"no shutter time at {nBad} of {tMid.size} detector grid points {why}: "
                        "no per-source times")
    dt = tMid - tMid[0]
    su = max(0.5 * (n[ia] - 1.0), 1.0)
    sv = max(0.5 * (n[1 - ia] - 1.0), 1.0)
    z = (pa - c[ia]) / su
    wv = (pc - c[1 - ia]) / sv
    basis = np.stack([z, z * z, wv, z * wv, wv * wv], axis=-1)
    unit = np.array([1 / su, su**-2.0, 1 / sv, 1 / (su * sv), sv**-2.0])
    cz = np.linalg.lstsq(basis[1:1 + _GRID**2], dt[1:1 + _GRID**2], rcond=None)[0]
    values["maxAbsResidual"] = float(np.max(np.abs(basis @ cz - dt)))
    values["coefficients"] = tuple(float(v) for v in cz * unit)
    if not (np.all(np.isfinite(values["coefficients"])) and math.isfinite(values["maxAbsResidual"])):
        raise _Unavailable(flags, "non-finite quadratic fit")
    return values, ""
