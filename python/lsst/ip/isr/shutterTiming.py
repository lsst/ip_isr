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

The LSSTCam shutter has two blades that cross the focal plane in ~0.9 s, so the
mid-exposure time of a pixel depends on its position: across the focal plane it
spans ~0.43 s.  This module computes, per detector, the flux-weighted
mid-exposure time at the detector centre and a quadratic in pixel coordinates
for the rest of the detector, from:

- the shutter Hall-sensor fit cards in the exposure metadata (``SHUTTER
  {OPEN,CLOSE} STARTTIME TAI MJD``, ``... SIDE``, ``... MODEL``, ``...
  HALLSENSORFIT {MODELSTARTTIME, PIVOTPOINT1, PIVOTPOINT2, JERK0, JERK1,
  JERK2}``; present since 2025-10-29), plus ``MJD-BEG``, ``MJD-END``,
  ``EXPTIME`` and ``SHUTTIME`` for the clock cross-checks;
- a beam model: the ray bundle's flux quantiles at the shutter plane as a
  function of field position (raytrace table, LCA-20578), shipped by the obs
  package;
- the detector geometry (``PIXELS -> FOCAL_PLANE``).

It is a port of the ``shutter_timing`` header-card path
(github.com/mjuric/shutter-timing, docs/design/stack-port.md).  Times are TAI
MJD.

WAVE-0 CONTRACT: the API below is fixed by the integrator; work packages
implement the bodies and must not change signatures, field names, enum values
or semantics.
"""

from __future__ import annotations

__all__ = [
    "DEGRADED_FLAGS",
    "DetectorGeometry",
    "ShutterBeamModel",
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

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import ConvexHull

import lsst.pex.config as pexConfig

from ._shutterTrajectory import ThreeJerksParams, ThreeJerksTrajectory, checkFitParams

_SECONDS_PER_DAY = 86400.0

#: Cumulative flux levels q of the beam table.
_BEAM_LEVELS = np.array([0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99])

#: Extrapolation margin (mm) outside the beam table's convex hull: the worst
#: science pixel (28.64 mm outside) plus 10 mm; NaN beyond.
_BEAM_MARGIN_MM = 28.64 + 10.0

#: 2-point Gauss-Legendre abscissae on [0, 1].
_GL2 = np.array([0.5 - 0.5 / np.sqrt(3.0), 0.5 + 0.5 / np.sqrt(3.0)])

_THREE_JERKS_V1 = "ThreeJerksModelv1"

#: ThreeJerksParams field -> card suffix.
_FIT_CARDS = (
    ("modelStartTime", "MODELSTARTTIME"),
    ("pivot1", "PIVOTPOINT1"),
    ("pivot2", "PIVOTPOINT2"),
    ("jerk0", "JERK0"),
    ("jerk1", "JERK1"),
    ("jerk2", "JERK2"),
)


class ShutterTimingStatus(enum.IntEnum):
    """Quality of a corrected time (detector level, or per source)."""

    OK = 0
    """Corrected; accurate to ~1 ms."""
    DEGRADED = 1
    """Corrected, reduced accuracy (see `ShutterTimingFlag` and
    `DEGRADED_FLAGS`). Still better than the header midpoint: callers use it.
    """
    UNAVAILABLE = 2
    """No corrected time: callers keep the header midpoint
    (``visitInfo.date``).
    """


class ShutterTimingFlag(enum.IntFlag):
    """Reasons, as a bit mask.  Bit values equal
    ``shutter_timing._contract.QC`` so the stack and the standalone package can
    be compared directly.
    """

    NONE = 0
    NO_PROFILE = 1
    """Missing or unusable Hall-fit cards (no fit, or a model other than
    ``ThreeJerksModelv1``): the mean profile is used (with MEAN_PROFILE). The
    result is UNAVAILABLE only if the start time or side cards are unusable
    too.
    """
    PARAM_RANGE = 2
    """A Hall fit outside the nominal parameter ranges (mean profile used
    instead).
    """
    CLOCK_OPEN_VS_BEG = 4
    """Open STARTTIME - MJD-BEG outside ``openMinusBegRange`` (flag only)."""
    CLOCK_END_VS_CLOSE = 8
    """MJD-END - close STARTTIME outside ``endMinusCloseRange`` (a late
    readout; flag only: MJD-END is not used).
    """
    CLOCK_CLOSE_VS_OPEN = 16
    """Close - open STARTTIME - EXPTIME outside
    ``closeMinusOpenMinusExptimeRange``: the two shutter clocks disagree, so
    both start times are re-anchored to MJD-BEG and EXPTIME with
    ``headerAnchorOffsets``.
    """
    PRE_CLOCK_EPOCH = 32
    """Zero point re-anchored to the header (set together with
    CLOCK_CLOSE_VS_OPEN).
    """
    MEAN_PROFILE = 64
    """A per-direction mean profile replaced a missing or out-of-range Hall
    fit.
    """
    FIT_RESIDUAL = 128
    """Not computed from header cards (reserved; JSON-profile path)."""
    HALL_VS_ENCODER = 256
    """Not computed from header cards (reserved; JSON-profile path)."""
    ACTION_DURATION = 512
    """Not computed from header cards (reserved; JSON-profile path)."""
    SHUTTIME_MISMATCH = 1024
    """|T_eff(focal-plane centre) - SHUTTIME| > ``shuttimeTolerance``
    (diagnostic only).
    """
    BEAM_EXTRAPOLATED = 2048
    """Part of the detector lies outside the beam table's coverage (detector
    level; per source, the hull test decides).
    """
    NO_HEADER = 4096
    """Reserved (profiles without header times; not reachable from exposure
    metadata).
    """


#: Detector-level flags that make a corrected time DEGRADED (``shutter_timing``
#: ``DEGRADED_QC``).  Not included: CLOCK_OPEN_VS_BEG (the shutter clock is
#: kept),
#: CLOCK_END_VS_CLOSE (MJD-END unused), SHUTTIME_MISMATCH (diagnostic) and
#: BEAM_EXTRAPOLATED (replaced by the per-source hull test).
DEGRADED_FLAGS = (
    ShutterTimingFlag.NO_PROFILE | ShutterTimingFlag.PARAM_RANGE | ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN
    | ShutterTimingFlag.PRE_CLOCK_EPOCH | ShutterTimingFlag.MEAN_PROFILE | ShutterTimingFlag.FIT_RESIDUAL
    | ShutterTimingFlag.HALL_VS_ENCODER | ShutterTimingFlag.ACTION_DURATION | ShutterTimingFlag.NO_HEADER
)


class ShutterTimingConfig(pexConfig.Config):
    """Configuration of `computeShutterTiming`.  Defaults are the values
    validated in ``shutter_timing`` v0.4.0 (calibration b8554546ab84);
    ``beamFile`` is set by the obs package's config overrides.
    """

    beamFile = pexConfig.Field(
        dtype=str, default="",
        doc="Path of the shutter-plane beam table (raytrace flux quantiles, LCA-20578 format). "
            "Empty: no beam model, so every result is UNAVAILABLE.",
    )
    # --- Hall-fit parameter ranges (PARAM_RANGE)
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
    # --- clock cross-checks
    openMinusBegRange = pexConfig.ListField(
        dtype=float, default=[5.0, 15.0], length=2, doc="Open STARTTIME - MJD-BEG range (ms).")
    endMinusCloseRange = pexConfig.ListField(
        dtype=float, default=[0.89, 0.96], length=2, doc="MJD-END - close STARTTIME range (s).")
    closeMinusOpenMinusExptimeRange = pexConfig.ListField(
        dtype=float, default=[-2.0, 4.0], length=2, doc="Close - open STARTTIME - EXPTIME range (ms).")
    shuttimeTolerance = pexConfig.Field(
        dtype=float, default=1.0, doc="|T_eff(centre) - SHUTTIME| tolerance (ms).")
    headerAnchorOffsets = pexConfig.ListField(
        dtype=float, default=[0.00809, 0.00052], length=2,
        doc="Re-anchoring offsets (s): open STARTTIME - MJD-BEG, and close - open STARTTIME - EXPTIME.",
    )
    # --- mechanics and coordinates
    a1Decreasing = pexConfig.Field(
        dtype=float, default=0.0, doc="Blade-edge offset A1 (mm) of motions with travel sign -1.")
    a1Increasing = pexConfig.Field(
        dtype=float, default=0.0, doc="Blade-edge offset A1 (mm) of motions with travel sign +1.")
    encoderCenter = pexConfig.Field(
        dtype=float, default=375.0,
        doc="Encoder position (mm) of the focal-plane centre: x_ccs = c - e + A1.")
    nominalStroke = pexConfig.Field(
        dtype=float, default=750.76, doc="Start position (mm) of a 750 -> 0 move (header cards carry none).")
    nominalStartIncreasing = pexConfig.Field(
        dtype=float, default=-0.05, doc="Start position (mm) of a 0 -> 750 move.")
    ccsFromDvcsSign = pexConfig.ChoiceField(
        dtype=int, default=1, allowed={1: "x_ccs = y_dvcs", -1: "x_ccs = -y_dvcs"},
        doc="CCS from DVCS (afw FOCAL_PLANE): x_ccs = s * y_dvcs, y_ccs = s * x_dvcs.",
    )
    meanProfileDecreasing = pexConfig.ListField(
        dtype=float, length=6,
        default=[0.0007497095570137927, 0.22426363506685815, 0.6776782303465883,
                 33100.38075466996, -32888.35381319897, 34383.386427176374],
        doc="Mean ThreeJerksModelv1 fit (ModelStartTime, PivotPoint1, PivotPoint2, Jerk0, Jerk1, Jerk2) "
            "for travel sign -1, used when a Hall fit is missing or out of range.",
    )
    meanProfileIncreasing = pexConfig.ListField(
        dtype=float, length=6,
        default=[0.0009719394948759783, 0.2237721807772504, 0.6788008740104505,
                 33178.9312914494, -32841.22646579422, 35458.25981717702],
        doc="Mean ThreeJerksModelv1 fit for travel sign +1.",
    )
    # --- per-detector representation and per-source rules
    gridAlong = pexConfig.Field(dtype=int, default=9, doc="Grid points along the blade axis for the fit.")
    gridAcross = pexConfig.Field(dtype=int, default=9, doc="Grid points across the blade axis for the fit.")
    degradedResidual = pexConfig.Field(
        dtype=float, default=1e-3, doc="Quadratic fit residual (s) above which the detector is DEGRADED.")
    offDetectorLimit = pexConfig.Field(
        dtype=float, default=100.0,
        doc="Sources up to this far (pixels) outside the detector are DEGRADED; beyond, UNAVAILABLE.",
    )


@dataclasses.dataclass(frozen=True)
class DetectorGeometry:
    """Affine pixel -> DVCS (afw FOCAL_PLANE, mm) map of one detector.

    Pixel convention: LSST/afw, integer coordinates at pixel centres; the
    detector spans ``[-0.5, nx - 0.5] x [-0.5, ny - 0.5]``.  ``fp = centerMm +
    jacobian @ (p - centerPixel)``.
    """

    detectorId: int
    nx: int
    ny: int
    centerPixel: tuple[float, float]
    centerMm: tuple[float, float]
    jacobian: tuple[tuple[float, float], tuple[float, float]]
    """[[dX/dx, dX/dy], [dY/dx, dY/dy]] in mm per pixel."""

    @classmethod
    def fromDetector(cls, detector) -> DetectorGeometry:
        """From an `lsst.afw.cameraGeom.Detector` (``PIXELS -> FOCAL_PLANE``),
        linearized at the detector centre (the LSSTCam map is affine to < 1e-3
        mm).
        """
        import lsst.geom
        from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS

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
        )

    def pixelToDvcs(self, x, y) -> tuple[np.ndarray, np.ndarray]:
        """Pixel -> DVCS (afw FOCAL_PLANE, mm), broadcast."""
        dx = np.asarray(x, dtype=float) - self.centerPixel[0]
        dy = np.asarray(y, dtype=float) - self.centerPixel[1]
        j = self.jacobian
        return (self.centerMm[0] + j[0][0] * dx + j[0][1] * dy,
                self.centerMm[1] + j[1][0] * dx + j[1][1] * dy)

    def pixelToCcs(self, x, y, ccsFromDvcsSign: int = 1) -> tuple[np.ndarray, np.ndarray]:
        """Pixel -> CCS (mm): ``x_ccs = s * Y_dvcs``, ``y_ccs = s * X_dvcs``.
        """
        xd, yd = self.pixelToDvcs(x, y)
        return ccsFromDvcsSign * yd, ccsFromDvcsSign * xd

    def offDetector(self, x, y) -> np.ndarray:
        """How far (pixels; the larger axis) each position lies outside the
        pixel bounds; 0 inside.
        """
        x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
        return np.maximum.reduce([np.zeros(x.shape), -0.5 - x, x - (self.nx - 0.5),
                                  -0.5 - y, y - (self.ny - 0.5)])

    def _bladeAxis(self) -> str:
        """The pixel axis (anti)parallel to DVCS y (the blade motion)."""
        j = self.jacobian
        return "x" if abs(j[1][0]) > abs(j[1][1]) else "y"

    def _pixelCorners(self) -> np.ndarray:
        """The four pixel-edge corners (4, 2)."""
        lo, hx, hy = -0.5, self.nx - 0.5, self.ny - 0.5
        return np.array([[lo, lo], [hx, lo], [hx, hy], [lo, hy]])


class ShutterBeamModel:
    """Shutter-plane beam model: for each field position (CCS mm), the blade
    coordinate s_q (CCS-x mm) below which a fraction q of the pixel's beam flux
    lies, at the tabulated levels q (0.01 ... 0.99); interpolated in the field
    position, and extrapolated (with a hull-distance test) outside the table's
    coverage.

    Parameters
    ----------
    xCcs, yCcs : `numpy.ndarray`
        The tabulated field positions (CCS mm), shape ``(npos,)``.
    sQ : `numpy.ndarray`
        s_q at each position and level, shape ``(npos, nlevels)``, strictly
        increasing in q.
    levels : `numpy.ndarray`
        The levels q.
    marginMm : `float`, optional
        Extrapolation margin outside the hull (mm); NaN beyond.

    Notes
    -----
    The cumulative profile F(s) is piecewise linear through (s_q, q), with
    linear tails reaching F = 0 at ``s_0.01 - (s_0.05 - s_0.01)`` and F = 1 at
    ``s_0.99 + (s_0.99 - s_0.95)``.  ``R_q = s_q - x_ccs`` is interpolated
    bilinearly in grid cells whose four nodes are tabulated, barycentrically
    on the Delaunay triangulation of the nodes elsewhere inside their convex
    hull, and held at the nearest hull point outside it (up to ``marginMm``).
    Port of ``shutter_timing.geometry.BeamModel``.
    """

    def __init__(self, xCcs, yCcs, sQ, levels=_BEAM_LEVELS, *, marginMm: float = _BEAM_MARGIN_MM):
        x = np.asarray(xCcs, dtype=float)
        y = np.asarray(yCcs, dtype=float)
        s = np.asarray(sQ, dtype=float)
        levels = np.asarray(levels, dtype=float)
        if s.shape != (x.size, levels.size):
            raise ValueError("sQ must have shape (npos, nlevels)")
        if not np.all(np.diff(s, axis=1) > 0):
            raise ValueError("beam quantiles are not strictly increasing in level")
        self._levels = levels
        self.marginMm = float(marginMm)
        self._points = np.column_stack([x, y])
        self._r = s - x[:, None]
        self._gx = np.unique(x)
        self._gy = np.unique(y)
        ix = np.searchsorted(self._gx, x)
        iy = np.searchsorted(self._gy, y)
        grid = np.full((self._gx.size, self._gy.size, levels.size), np.nan)
        grid[ix, iy] = self._r
        self._grid = grid
        have = np.isfinite(grid[..., 0])
        self._cellOk = have[:-1, :-1] & have[1:, :-1] & have[:-1, 1:] & have[1:, 1:]
        self._tri = LinearNDInterpolator(self._points, self._r)
        hull = ConvexHull(self._points)
        self._hullVertices = hull.vertices
        self._hullEq = hull.equations
        knots = np.concatenate([[0.0], levels, [1.0]])
        du = np.diff(knots)
        self._knots = knots
        self._nodes = (knots[:-1, None] + du[:, None] * _GL2[None, :]).ravel()
        self._weights = np.repeat(du / 2.0, 2)

    @classmethod
    def fromFile(cls, path: str) -> ShutterBeamModel:
        """Read a beam table in the raytrace ``.tnt`` format."""
        from lsst.resources import ResourcePath

        text = ResourcePath(path).read().decode()
        f = io.StringIO(text)
        first = f.readline().strip()
        if first.upper() == "DATA":
            f.readline()
            data = np.loadtxt(f, ndmin=2)
        else:
            data = np.loadtxt(f, delimiter=",", ndmin=2)
        if data.shape[1] != 5:
            raise ValueError(f"{path}: expected 5 columns, got {data.shape[1]}")
        rows = data[:, [0, 1, 2, 4]]
        levels = np.unique(rows[:, 2])
        if levels.shape != _BEAM_LEVELS.shape or not np.allclose(levels, _BEAM_LEVELS):
            raise ValueError(f"{path}: unexpected levels {levels}")
        pos, inv = np.unique(rows[:, :2], axis=0, return_inverse=True)
        inv = inv.ravel()
        il = np.searchsorted(levels, rows[:, 2])
        s = np.full((pos.shape[0], levels.size), np.nan)
        s[inv, il] = rows[:, 3]
        if np.isnan(s).any() or len(rows) != s.size:
            raise ValueError(f"{path}: incomplete or duplicated (position, level) table")
        return cls(pos[:, 0], pos[:, 1], s, _BEAM_LEVELS.copy())

    @property
    def levels(self) -> np.ndarray:
        """The flux levels q."""
        return self._levels.copy()

    def _isOutside(self, p):
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

    def hullDistance(self, xCcs, yCcs) -> np.ndarray:
        """Distance (mm) outside the convex hull of the table nodes; 0 inside.
        """
        x, y = np.broadcast_arrays(np.asarray(xCcs, dtype=float), np.asarray(yCcs, dtype=float))
        p = np.column_stack([x.ravel(), y.ravel()])
        out = np.zeros(p.shape[0])
        with np.errstate(invalid="ignore"):
            outside = self._isOutside(p)
        if outside.any():
            out[outside] = self._hullProject(p[outside])[1]
        out[~np.all(np.isfinite(p), axis=1)] = np.nan
        return out.reshape(x.shape)

    def _residual(self, x, y, extrapolate=True):
        """R_q = s_q - x_ccs at flat arrays x, y; shape (n, nlevels)."""
        n = x.size
        out = np.full((n, self._levels.size), np.nan)
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
            inside = ~self._isOutside(p)
            r = np.full((p.shape[0], self._levels.size), np.nan)
            if inside.any():
                r[inside] = self._tri(p[inside])
            need = ~np.isfinite(r[:, 0])
            if extrapolate and need.any():
                q, dist = self._hullProject(p[need])
                c = self._points.mean(axis=0)
                q += 1e-6 * (c - q) / np.linalg.norm(c - q, axis=1, keepdims=True)
                rr = self._residual(q[:, 0], q[:, 1], extrapolate=False)
                rr[dist > self.marginMm] = np.nan
                r[need] = rr
            out[rest] = r
        return out

    def _quadratureNodes(self, xCcs, yCcs):
        """Blade coordinates s(u_k) at the 20 Gauss nodes of the beam CDF,
        shape ``(n, 20)`` for flat ``xCcs``, ``yCcs``; and the weights.
        """
        xf = np.asarray(xCcs, dtype=float).ravel()
        yf = np.asarray(yCcs, dtype=float).ravel()
        s = self._residual(xf, yf) + xf[:, None]
        lo = s[:, :1] - (s[:, 1:2] - s[:, :1])
        hi = s[:, -1:] + (s[:, -1:] - s[:, -2:-1])
        ext = np.concatenate([lo, s, hi], axis=-1)
        knots, u = self._knots, self._nodes
        k = np.clip(np.searchsorted(knots, u, side="right") - 1, 0, knots.size - 2)
        f = (u - knots[k]) / (knots[k + 1] - knots[k])
        return ext[..., k] * (1 - f) + ext[..., k + 1] * f, self._weights


@functools.lru_cache(maxsize=4)
def loadShutterBeam(path: str) -> ShutterBeamModel:
    """`ShutterBeamModel.fromFile`, cached per process (call sites load it per
    quantum).
    """
    return ShutterBeamModel.fromFile(path)


@dataclasses.dataclass(frozen=True)
class ShutterTiming:
    """Shutter-corrected mid-exposure times of one detector of one exposure.

    Times are MJD TAI; durations in seconds.  The per-detector representation
    is ``t(x, y) = centerMjdTai + (c_u u + c_uu u^2 + c_v v + c_uv u v + c_vv
    v^2) / 86400`` with ``(u, v)`` the pixel offsets from
    ``geometry.centerPixel`` along and across the blade axis (``axis`` names
    the pixel axis along the blades: ``"x"`` -> u = x - cx, v = y - cy; ``"y"``
    -> u = y - cy, v = x - cx).

    When ``status`` is UNAVAILABLE the numeric fields are NaN and the per-
    source methods return NaN / UNAVAILABLE everywhere.
    """

    status: ShutterTimingStatus
    """Detector level: UNAVAILABLE if no timing; DEGRADED if ``flags &
    DEGRADED_FLAGS`` or ``maxAbsResidual > degradedResidual``; else OK.
    """
    flags: ShutterTimingFlag
    """Exposure-level flags | detector-level flags (``shutter_timing`` row
    ``qc_flags``).
    """
    message: str
    """Human-readable reason when not OK (for logs and task metadata); "" when
    OK.
    """
    detectorId: int
    axis: str
    centerMjdTai: float
    """Flux-weighted mid-exposure time at the detector centre."""
    coefficients: tuple[float, float, float, float, float]
    """(c_u, c_uu, c_v, c_uv, c_vv) in s / pixel^n."""
    maxAbsResidual: float
    """Largest |quadratic - exact| over the fit grid (s)."""
    effectiveExposureTime: float
    """Flux-weighted open time at the detector centre (s)."""
    focalPlaneMjdTai: float
    """Mid-exposure time at the focal-plane centre (DVCS 0, 0): the visit
    epoch.
    """
    headerMidMjdTai: float
    """(MJD-BEG + MJD-END) / 2 from the metadata, NaN if missing (what
    ``visitInfo.date`` holds today; for diagnostics only).
    """
    policy: str
    """Zero point: "profile" (shutter clock) or "header_anchor" (re-anchored to
    MJD-BEG).
    """
    geometry: DetectorGeometry | None
    beam: ShutterBeamModel | None = dataclasses.field(repr=False, compare=False)
    config: ShutterTimingConfig | None = dataclasses.field(repr=False, compare=False)

    def tMidMjdTai(self, x, y) -> np.ndarray:
        """Per-source mid-exposure times (MJD TAI) at pixel positions; NaN
        where the per-source status is UNAVAILABLE.  Vectorized; shape of
        ``np.broadcast(x, y)``.
        """
        x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
        status = self.sourceStatus(x, y)
        out = np.full(x.shape, np.nan)
        good = status != ShutterTimingStatus.UNAVAILABLE
        if good.any():
            cx, cy = self.geometry.centerPixel
            xs, ys = x[good] - cx, y[good] - cy
            u, v = (xs, ys) if self.axis == "x" else (ys, xs)
            cU, cUU, cV, cUV, cVV = self.coefficients
            ds = cU * u + cUU * u * u + cV * v + cUV * u * v + cVV * v * v
            out[good] = self.centerMjdTai + ds / _SECONDS_PER_DAY
        return out

    def sourceStatus(self, x, y) -> np.ndarray:
        """Per-source `ShutterTimingStatus` values (uint8), as
        ``shutter_timing`` ``corrected_midpoints``:

        - UNAVAILABLE: detector UNAVAILABLE, non-finite x or y, or more than
          ``offDetectorLimit`` pixels outside the detector;
        - DEGRADED: detector DEGRADED, or outside the detector (up to the
          limit), or the position's CCS coordinates outside the beam table's
          hull (``hullDistance > 0``);
        - OK otherwise.
        """
        x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
        if self.status == ShutterTimingStatus.UNAVAILABLE or self.geometry is None:
            return np.full(x.shape, ShutterTimingStatus.UNAVAILABLE, dtype=np.uint8)
        out = np.full(x.shape, self.status, dtype=np.uint8)
        finite = np.isfinite(x) & np.isfinite(y)
        off = np.where(finite, self.geometry.offDetector(np.where(finite, x, 0.0),
                                                         np.where(finite, y, 0.0)), np.inf)
        limit = self.config.offDetectorLimit if self.config is not None else 100.0
        sign = self.config.ccsFromDvcsSign if self.config is not None else 1
        unavailable = ~(off <= limit)
        degraded = off > 0
        if self.beam is not None and (~unavailable).any():
            xc, yc = self.geometry.pixelToCcs(x[~unavailable], y[~unavailable], sign)
            hull = np.zeros(x.shape)
            hull[~unavailable] = self.beam.hullDistance(xc, yc)
            degraded |= hull > 0
        out[degraded] = np.maximum(out[degraded], ShutterTimingStatus.DEGRADED)
        out[unavailable] = ShutterTimingStatus.UNAVAILABLE
        return out

    def summary(self) -> dict:
        """Scalars for task metadata: status (name), flags (int), message,
        centerMjdTai, focalPlaneMjdTai, centerMinusHeaderMid (s),
        maxAbsResidual, policy.
        """
        return dict(
            status=self.status.name,
            flags=int(self.flags),
            message=self.message,
            centerMjdTai=float(self.centerMjdTai),
            focalPlaneMjdTai=float(self.focalPlaneMjdTai),
            centerMinusHeaderMid=float((self.centerMjdTai - self.headerMidMjdTai) * _SECONDS_PER_DAY),
            maxAbsResidual=float(self.maxAbsResidual),
            policy=self.policy,
        )


def computeShutterTiming(metadata, detector, config: ShutterTimingConfig | None = None, *,
                         beam: ShutterBeamModel | None = None) -> ShutterTiming:
    """Shutter-corrected times of one detector of one exposure.

    Parameters
    ----------
    metadata : `lsst.daf.base.PropertyList` or `collections.abc.Mapping`
        The exposure metadata (``exposure.metadata``); keys without the
        ``HIERARCH`` prefix (``"SHUTTER OPEN STARTTIME TAI MJD"``).
    detector : `lsst.afw.cameraGeom.Detector` or `DetectorGeometry`
        The detector (an afw Detector is converted with
        `DetectorGeometry.fromDetector`).
    config : `ShutterTimingConfig`, optional
        Defaults to ``ShutterTimingConfig()``.
    beam : `ShutterBeamModel`, optional
        Defaults to ``loadShutterBeam(config.beamFile)``; if that is empty, the
        result is UNAVAILABLE with a message.

    Returns
    -------
    timing : `ShutterTiming`
        Never raises for problems in the data (missing or malformed cards,
        failed fits, geometry outside the beam table): those give UNAVAILABLE
        or DEGRADED with ``flags`` and ``message``.  Raises only for invalid
        ``config`` or arguments.

    Notes
    -----
    Algorithm (``shutter_timing`` header-card path, ``table.compute_rows`` with
    ``TABLE_SETTINGS = dict(kind="hall_fit", quadrature="weighted", order=2,
    n_along=9, n_across=9)`` and policy "auto"):

    1. Parse the open and close motions; a missing or malformed start time
       or side card -> NO_PROFILE, UNAVAILABLE; a missing or malformed
       Hall-fit card, or a model other than ``ThreeJerksModelv1`` ->
       NO_PROFILE and the mean profile of its travel direction
       (MEAN_PROFILE, DEGRADED).  Start positions: ``nominalStroke`` for 750
       -> 0 moves, ``nominalStartIncreasing`` for 0 -> 750 moves.
    2. Exposure QC: PARAM_RANGE, the clock checks, SHUTTIME_MISMATCH.
    3. A Hall fit that is out of range -> the mean profile of its travel
       direction (MEAN_PROFILE).
    4. CLOCK_CLOSE_VS_OPEN -> re-anchor both start times to MJD-BEG / EXPTIME
       with ``headerAnchorOffsets`` (PRE_CLOCK_EPOCH; policy "header_anchor").
    5. ThreeJerks trajectories; flux-weighted crossing times <T> of each blade
       at each grid point by Gauss quadrature over the beam's flux quantiles;
       t_mid = (<T>o + <T>c) / 2.
    6. Least-squares quadratic in (u, v) over a ``gridAlong x gridAcross`` grid
       spanning the detector; residual; BEAM_EXTRAPOLATED where the detector
       leaves the hull.

    Agreement with ``shutter_timing`` (wave-0 fixtures,
    tests/data/shutterTiming): centre times and the quadratic terms at a
    2000-pixel lever arm to <= 10 us, identical ``axis``, ``flags`` and
    statuses.
    """
    if config is None:
        config = ShutterTimingConfig()
    elif not isinstance(config, ShutterTimingConfig):
        raise TypeError(f"config must be a ShutterTimingConfig, not {type(config).__name__}")
    config.validate()
    if isinstance(detector, DetectorGeometry):
        geometry = detector
        notScience = None
    elif hasattr(detector, "getTransform") and hasattr(detector, "getBBox"):
        geometry = DetectorGeometry.fromDetector(detector)
        notScience = _notScienceReason(detector)
    else:
        raise TypeError("detector must be an lsst.afw.cameraGeom.Detector or a DetectorGeometry, "
                        f"not {type(detector).__name__}")
    if not (hasattr(metadata, "get") and hasattr(metadata, "__contains__")):
        raise TypeError(f"metadata must be a PropertyList or a Mapping, not {type(metadata).__name__}")
    if beam is None and config.beamFile:
        beam = loadShutterBeam(config.beamFile)

    ctx = _Result(geometry=geometry, config=config, beam=beam)
    try:
        ctx.headerMid = _headerMid(metadata)
        if beam is None:
            raise _Unavailable(ShutterTimingFlag.NONE, "no beam model (config.beamFile is empty)")
        if notScience:
            raise _Unavailable(ShutterTimingFlag.NONE, notScience)
        _compute(ctx, metadata)
    except _Unavailable as e:
        return ctx.unavailable(e.flags, e.message)
    except Exception as e:  # noqa: BLE001 -- bad data must never raise (contract)
        return ctx.unavailable(ctx.flags, f"shutter timing failed: {type(e).__name__}: {e}")
    return ctx.result()


# --------------------------------------------------------------------------
# Private implementation (port of the shutter_timing header-card path).
# --------------------------------------------------------------------------


class _Unavailable(Exception):
    """No corrected time for this detector (flags, message)."""

    def __init__(self, flags, message):
        super().__init__(message)
        self.flags = ShutterTimingFlag(int(flags))
        self.message = message


def _notScienceReason(detector):
    """A message if an afw detector is not a science detector (the beam table
    describes only the science beams), else None.
    """
    try:
        from lsst.afw.cameraGeom import DetectorType

        if detector.getType() != DetectorType.SCIENCE:
            return (f"detector {detector.getId()} ({detector.getName()}) is "
                    f"{detector.getType().name}, not SCIENCE: the beam table does not cover it")
    except Exception:  # noqa: BLE001 -- detectors without a type are treated as science
        return None
    return None


def _card(metadata, key):
    """The value of ``key`` (or ``HIERARCH key``) in ``metadata``, or None.
    """
    for k in (key, "HIERARCH " + key):
        try:
            if k in metadata:
                return metadata.get(k)
        except Exception:  # noqa: BLE001 -- unreadable card: treat as missing
            return None
    return None


def _float(value):
    """float(value), or None for missing, empty, non-numeric or non-finite.
    """
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, str):
        value = value.strip()
        if not value:
            return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _str(value):
    if value is None:
        return None
    value = str(value).strip()
    return value or None


def _headerMid(metadata):
    beg = _float(_card(metadata, "MJD-BEG"))
    end = _float(_card(metadata, "MJD-END"))
    return 0.5 * (beg + end) if beg is not None and end is not None else math.nan


@dataclasses.dataclass
class _Motion:
    """One blade motion from the header cards."""

    which: str
    side: str
    startMjdTai: float
    fit: ThreeJerksParams | None
    isOpen: bool
    model: str | None = None

    @property
    def usable(self):
        """A Hall fit of the supported model is present."""
        return self.fit is not None and (self.model is None or self.model == _THREE_JERKS_V1)

    @property
    def travelSign(self):
        decreasing = (self.side == "PLUSX") == self.isOpen  # PLUSX-open, MINUSX-close: 750 -> 0
        return -1 if decreasing else 1


def _readMotion(metadata, which):
    pre = f"SHUTTER {which}"
    start = _float(_card(metadata, f"{pre} STARTTIME TAI MJD"))
    side = _str(_card(metadata, f"{pre} SIDE"))
    if start is None or side is None:
        raise _Unavailable(ShutterTimingFlag.NO_PROFILE,
                           f"no usable '{pre} STARTTIME TAI MJD' / '{pre} SIDE' cards")
    if side.upper() not in ("PLUSX", "MINUSX"):
        raise _Unavailable(ShutterTimingFlag.NO_PROFILE, f"unknown shutter side {pre} SIDE = {side!r}")
    model = _str(_card(metadata, f"{pre} MODEL"))
    vals = {f: _float(_card(metadata, f"{pre} HALLSENSORFIT {c}")) for f, c in _FIT_CARDS}
    fit = None if any(v is None for v in vals.values()) else ThreeJerksParams(**vals)
    return _Motion(which=which, side=side.upper(), startMjdTai=start, fit=fit,
                   isOpen=(which == "OPEN"), model=model)


def _inRange(value, bounds):
    lo, hi = bounds
    return bool(np.isfinite(value) and lo <= value <= hi)


def _flagNames(flags):
    return "|".join(f.name for f in ShutterTimingFlag if f and (flags & f)) or "NONE"


@dataclasses.dataclass
class _Result:
    """Mutable accumulator of `computeShutterTiming`."""

    geometry: DetectorGeometry
    config: ShutterTimingConfig
    beam: ShutterBeamModel | None
    flags: ShutterTimingFlag = ShutterTimingFlag.NONE
    headerMid: float = math.nan
    policy: str = ""
    centerMjdTai: float = math.nan
    coefficients: tuple = (math.nan,) * 5
    maxAbsResidual: float = math.nan
    effectiveExposureTime: float = math.nan
    focalPlaneMjdTai: float = math.nan

    def unavailable(self, flags, message):
        return ShutterTiming(
            status=ShutterTimingStatus.UNAVAILABLE, flags=ShutterTimingFlag(int(flags)),
            message=message, detectorId=self.geometry.detectorId, axis=self.geometry._bladeAxis(),
            centerMjdTai=math.nan, coefficients=(math.nan,) * 5, maxAbsResidual=math.nan,
            effectiveExposureTime=math.nan, focalPlaneMjdTai=math.nan, headerMidMjdTai=self.headerMid,
            policy=self.policy, geometry=self.geometry, beam=self.beam, config=self.config,
        )

    def result(self):
        flags = ShutterTimingFlag(int(self.flags))
        reasons = []
        if flags & DEGRADED_FLAGS:
            reasons.append(_flagNames(flags & DEGRADED_FLAGS))
        if not self.maxAbsResidual <= self.config.degradedResidual:
            reasons.append(f"quadratic residual {self.maxAbsResidual * 1e3:.3f} ms > "
                           f"{self.config.degradedResidual * 1e3:g} ms")
        status = ShutterTimingStatus.DEGRADED if reasons else ShutterTimingStatus.OK
        return ShutterTiming(
            status=status, flags=flags, message="; ".join(reasons),
            detectorId=self.geometry.detectorId, axis=self.geometry._bladeAxis(),
            centerMjdTai=self.centerMjdTai, coefficients=tuple(self.coefficients),
            maxAbsResidual=self.maxAbsResidual, effectiveExposureTime=self.effectiveExposureTime,
            focalPlaneMjdTai=self.focalPlaneMjdTai, headerMidMjdTai=self.headerMid, policy=self.policy,
            geometry=self.geometry, beam=self.beam, config=self.config,
        )


def _compute(ctx, metadata):
    """Fill ``ctx`` (`_Result`) from the metadata; raises `_Unavailable`."""
    config, geometry, beam = ctx.config, ctx.geometry, ctx.beam
    mOpen = _readMotion(metadata, "OPEN")
    mClose = _readMotion(metadata, "CLOSE")
    exptime = _float(_card(metadata, "EXPTIME"))
    shuttime = _float(_card(metadata, "SHUTTIME"))
    beg = _float(_card(metadata, "MJD-BEG"))
    end = _float(_card(metadata, "MJD-END"))

    # ---- exposure QC (shutter_timing core.qc_exposure, header path)
    flags = ShutterTimingFlag.NONE
    if not (mOpen.usable and mClose.usable):
        flags |= ShutterTimingFlag.NO_PROFILE
    inRange = [checkFitParams(m.fit, config) for m in (mOpen, mClose)]
    if any(m.fit is not None and not ok for m, ok in zip((mOpen, mClose), inRange)):
        flags |= ShutterTimingFlag.PARAM_RANGE
    fitOk = [m.usable and ok for m, ok in zip((mOpen, mClose), inRange)]
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
    ctx.flags = flags

    # ---- shape: missing, unusable or out-of-range Hall fits -> the
    # per-direction mean profile
    for m, ok in zip((mOpen, mClose), fitOk):
        if not ok:
            mean = config.meanProfileDecreasing if m.travelSign < 0 else config.meanProfileIncreasing
            m.fit = ThreeJerksParams(*(float(v) for v in mean))
            flags |= ShutterTimingFlag.MEAN_PROFILE
    ctx.flags = flags

    # ---- zero point: re-anchor to the header if the two shutter clocks
    # disagree
    closeRelProfile = (mClose.startMjdTai - mOpen.startMjdTai) * _SECONDS_PER_DAY
    if flags & ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN:
        if beg is None or exptime is None:
            raise _Unavailable(flags, "CLOCK_CLOSE_VS_OPEN, and re-anchoring needs MJD-BEG and EXPTIME")
        offBeg, offCloseOpen = (float(v) for v in config.headerAnchorOffsets)
        openStart = beg + offBeg / _SECONDS_PER_DAY
        closeRel = float(exptime) + offCloseOpen
        flags |= ShutterTimingFlag.PRE_CLOCK_EPOCH
        ctx.policy = "header_anchor"
    else:
        openStart = mOpen.startMjdTai
        closeRel = closeRelProfile
        ctx.policy = "profile"
    ctx.flags = flags

    # ---- trajectories
    trajs, offsets, a1s = [], [], []
    for m in (mOpen, mClose):
        sign = m.travelSign
        start = config.nominalStroke if sign < 0 else config.nominalStartIncreasing
        trajs.append(ThreeJerksTrajectory(m.fit, start, sign))
        offsets.append(float(m.fit.modelStartTime))
        a1s.append(config.a1Decreasing if sign < 0 else config.a1Increasing)

    # ---- evaluation points: centre, fit grid, staggered grid, focal-plane
    # centre
    na, nc = int(config.gridAlong), int(config.gridAcross)
    if na < 3 or nc < 3:
        raise ValueError("gridAlong and gridAcross must be >= 3")
    axis = geometry._bladeAxis()
    ia = 0 if axis == "x" else 1
    c = geometry.centerPixel
    n = (geometry.nx, geometry.ny)
    along = np.linspace(0.0, n[ia] - 1.0, na)
    across = np.linspace(0.0, n[1 - ia] - 1.0, nc)
    ga, gc = np.meshgrid(along, across, indexing="ij")
    ha, hc = np.meshgrid(0.5 * (along[1:] + along[:-1]), 0.5 * (across[1:] + across[:-1]), indexing="ij")
    pa = np.concatenate([[c[ia]], ga.ravel(), ha.ravel()])
    pc = np.concatenate([[c[1 - ia]], gc.ravel(), hc.ravel()])
    px, py = (pa, pc) if ia == 0 else (pc, pa)
    sign = config.ccsFromDvcsSign
    xc, yc = geometry.pixelToCcs(px, py, sign)
    xc = np.append(xc, 0.0)
    yc = np.append(yc, 0.0)

    # ---- flux-weighted crossing times (s after the open start time)
    s, w = beam._quadratureNodes(xc, yc)
    rel = []
    for traj, off, a1 in zip(trajs, offsets, a1s):
        d = traj.travelSign * ((config.encoderCenter + a1 - traj.startPosition) - s)
        rel.append(np.asarray(traj.timeOfDisplacement(d)) + off)
    tauO, tauC = rel
    tO = tauO
    tC = tauC + closeRel
    with np.errstate(invalid="ignore"):
        gap = np.min(tC, axis=-1) - np.max(tO, axis=-1)
    overlap = gap <= 0
    eO, eC = tO @ w, tC @ w
    eO = np.where(overlap, np.nan, eO)
    eC = np.where(overlap, np.nan, eC)

    # ---- SHUTTIME diagnostic (profile zero point, focal-plane centre)
    if shuttime is not None and not (flags & (ShutterTimingFlag.NO_PROFILE | ShutterTimingFlag.PARAM_RANGE)):
        tCp = tauC[-1] + closeRelProfile
        tEff = (tCp @ w) - (tauO[-1] @ w) if np.min(tCp) > np.max(tauO[-1]) else math.nan
        if not abs((tEff - shuttime) * 1e3) <= config.shuttimeTolerance:
            flags |= ShutterTimingFlag.SHUTTIME_MISMATCH

    # ---- beam coverage of the detector (its corners, by convexity)
    corners = geometry._pixelCorners()
    if np.max(beam.hullDistance(*geometry.pixelToCcs(corners[:, 0], corners[:, 1], sign))) > 0:
        flags |= ShutterTimingFlag.BEAM_EXTRAPOLATED
    ctx.flags = flags

    # ---- visit epoch
    ctx.focalPlaneMjdTai = float(openStart + 0.5 * (eO[-1] + eC[-1]) / _SECONDS_PER_DAY)

    # ---- quadratic fit
    eO, eC = eO[:-1], eC[:-1]
    tMid = 0.5 * (eO + eC)
    if not np.all(np.isfinite(tMid)):
        nBad = int((~np.isfinite(tMid)).sum())
        raise _Unavailable(flags, f"no shutter time at {nBad} of {tMid.size} detector grid points "
                           "(outside the beam table's margin, off the trajectory branch, or "
                           "overlapping blades)")
    dt = tMid - tMid[0]
    su = max(0.5 * (n[ia] - 1.0), 1.0)
    sv = max(0.5 * (n[1 - ia] - 1.0), 1.0)
    z = (pa - c[ia]) / su
    wv = (pc - c[1 - ia]) / sv
    basis = np.stack([z, z * z, wv, z * wv, wv * wv], axis=-1)
    unit = np.array([1 / su, su**-2.0, 1 / sv, 1 / (su * sv), sv**-2.0])
    fit = slice(1, 1 + na * nc)
    cz = np.linalg.lstsq(basis[fit], dt[fit], rcond=None)[0]
    ctx.maxAbsResidual = float(np.max(np.abs(basis @ cz - dt)))
    ctx.coefficients = tuple(float(v) for v in cz * unit)
    ctx.centerMjdTai = float(openStart + tMid[0] / _SECONDS_PER_DAY)
    ctx.effectiveExposureTime = float(eC[0] - eO[0])
    if not (np.all(np.isfinite(ctx.coefficients)) and math.isfinite(ctx.focalPlaneMjdTai)
            and math.isfinite(ctx.maxAbsResidual)):
        raise _Unavailable(flags, "non-finite quadratic fit or visit epoch")
