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
   ...`` metadata cards (present since 2025-10-29), which
   `lsst.ip.isr.ShutterMotionProfileFull` reads.  Inverting the fit gives the
   time at which a blade edge reaches a given position.
2. *Beam-weighted crossing times.*  A pixel's converging ray bundle is
   ~67 mm wide at the shutter plane, so an edge uncovers (or covers) it
   gradually.  The flux-weighted crossing time of one edge is
   ``<T> = integral T(s) dF(s)``, with ``F`` the cumulative flux profile of
   the bundle along the blade travel.  The profile has the same shape
   everywhere on the focal plane and moves with the pixel: its levels lie
   at ``beamScale * x + beamOffsets`` (`ShutterTimingConfig`), a linear fit
   to the LSST Camera raytrace of the beam (LCA-20578).  The mid-exposure
   time is ``(<T>_open + <T>_close) / 2``.
3. *Per-detector quadratic.*  These times are computed on a grid spanning
   the detector, and a quadratic in the pixel offsets along and across the
   blade direction is fit to them.

Times are TAI MJD.  The shutter motion and its Hall-sensor fit are described
in CTN-002, *Camera Shutter Motion Analysis* (https://ctn-002.lsst.io).
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
import math

import numpy as np

import lsst.geom
import lsst.pex.config as pexConfig
from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS, Detector
from lsst.daf.base import PropertyList

from .shutterMotion import ShutterMotionProfile, ShutterMotionProfileFull


_SECONDS_PER_DAY = 86400.0

#: Cumulative flux levels q of ``ShutterTimingConfig.beamOffsets``.
_BEAM_LEVELS = np.array([0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99])

#: 2-point Gauss-Legendre abscissae on [0, 1].
_GL2 = np.array([0.5 - 0.5 / np.sqrt(3.0), 0.5 + 0.5 / np.sqrt(3.0)])

_THREE_JERKS_V1 = "ThreeJerksModelv1"

#: Name of the Hall-sensor fit in a ShutterMotionProfile.
_HALL_SENSOR_FIT = "hallSensorFit"

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


class ShutterTimingConfig(pexConfig.Config):
    """Configuration of `computeShutterTiming`.

    The defaults are for LSSTCam: a least-squares fit to the LSST Camera
    raytrace of the beam at the shutter plane (LCA-20578, rotator angle 0),
    on its 137 field positions within 280 mm of the focal-plane centre.
    Against on-sky satellite edge times it does as well as the raytrace
    table within 280 mm, and better beyond, where the table is coarse or
    extrapolated.
    """

    beamScale = pexConfig.Field(
        dtype=float, default=1.03014,
        doc="Position of the centre of a pixel's beam at the shutter plane along the blade travel, "
            "per mm of the pixel's position along it (CCS x).",
    )
    beamOffsets = pexConfig.ListField(
        dtype=float,
        default=[-33.618, -30.133, -27.032, -19.493, 0.001, 19.494, 27.031, 30.130, 33.616],
        length=len(_BEAM_LEVELS),
        listCheck=lambda offsets: all(a < b for a, b in zip(offsets[:-1], offsets[1:])),
        doc="Offsets (mm, increasing) from beamScale * (CCS x) of the blade positions below which "
            "1, 5, 10, 25, 50, 75, 90, 95 and 99% of a pixel's beam flux crosses the shutter plane.",
    )


@dataclasses.dataclass(frozen=True)
class _ThreeJerksParams:
    """One ThreeJerksModelv1 fit."""

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
        """Return the end (model time, s) of the monotonic branch: the first
        velocity zero or velocity minimum after 0; 0 if ``j0 <= 0``; inf if
        none.
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
        """Convert pixel to CCS coordinates (mm): ``x_ccs = Y_dvcs``,
        ``y_ccs = X_dvcs``.
        """
        dx = x - self.centerPixel[0]
        dy = y - self.centerPixel[1]

        j = self.jacobian
        xCcs = self.centerMm[1] + j[1][0]*dx + j[1][1]*dy
        yCcs = self.centerMm[0] + j[0][0]*dx + j[0][1]*dy

        return xCcs, yCcs

    def bladeAxis(self) -> str:
        """Return the pixel axis (anti)parallel to DVCS y (the blade
        motion).
        """
        j = self.jacobian
        return "x" if abs(j[1][0]) > abs(j[1][1]) else "y"


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
        """Compute per-source mid-exposure times at pixel positions.

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
    """Compute the shutter-corrected times of one detector of one exposure.

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
        Defaults to ``ShutterTimingConfig()``.

    Returns
    -------
    timing : `ShutterTiming`
        UNAVAILABLE (callers keep the header midpoint) if:

        - there are no shutter cards (``SHUTTER OPEN STARTTIME TAI MJD``
          absent; exposures before 2025-10-29): flag NO_PROFILE;
        - either Hall fit is outside the nominal ranges: flag PARAM_RANGE;
        - there is no time at the focal-plane centre, the detector centre or
          any fit-grid point (overlapping blades in very short exposures).

        Otherwise OK.

    Raises
    ------
    ValueError
        Raised if a SIDE card is not PLUSX or MINUSX, or a MODEL card is not
        ThreeJerksModelv1.
    RuntimeError
        Raised if the shutter cards are present but one of them is missing.
    KeyError, ValueError, TypeError
        Raised if MJD-BEG or MJD-END is missing, or a card is not a number.

    Notes
    -----
    Rules beyond the method in the module docstring:

    - Blade start positions are nominal: 750.76 mm for 750 -> 0 moves, -0.05
      mm for 0 -> 750 moves.
    - If close - open STARTTIME - ``exposureTime`` is outside -2 to 4 ms, the
      shutter clocks disagree: both start times are re-anchored to MJD-BEG
      and ``exposureTime`` (flag CLOCK_CLOSE_VS_OPEN).
    - The beam's cumulative flux profile is piecewise linear through the
      levels of ``config.beamOffsets``, with linear tails reaching 0 and 1
      one interval beyond the first and last level.  Each edge's crossing
      time is integrated over it by 2-point Gauss-Legendre quadrature in
      each interval between levels.  The quadratic is fit on a 9 x 9 grid;
      ``maxAbsResidual`` is over that grid, its cell centres and the
      detector centre.
    """
    if config is None:
        config = ShutterTimingConfig()

    geometry = _DetectorGeometry.fromDetector(detector)
    common = dict(detectorId=geometry.detectorId, axis=geometry.bladeAxis(), geometry=geometry)

    # Exposures before the shutter cards were added.
    if "SHUTTER OPEN STARTTIME TAI MJD" not in metadata:
        return ShutterTiming(flags=ShutterTimingFlag.NO_PROFILE, message="no shutter motion cards", **common)

    profile = ShutterMotionProfileFull.fromMetadata(metadata)
    begin = float(metadata["MJD-BEG"])
    common["headerMidpointMjdTai"] = 0.5*(begin + float(metadata["MJD-END"]))
    try:
        values = _compute(profile, begin, geometry, config, exposureTime)
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
    """Return whether lo <= value <= hi (False for NaN)."""
    return bounds[0] <= value <= bounds[1]


def _fitInRange(fit: _ThreeJerksParams) -> bool:
    """Return whether a Hall fit is within the nominal ranges."""
    displacementAt0p9s = _ThreeJerksTrajectory(fit, 0.0, -1).sva(0.9)[0]

    return (
        0.0 <= fit.pivot1 <= fit.pivot2
        and _inRange(fit.pivot1, _PIVOT1_RANGE)
        and _inRange(fit.pivot2, _PIVOT2_RANGE)
        and all(_inRange(abs(j), _ABS_JERK_RANGE) for j in (fit.jerk0, fit.jerk1, fit.jerk2))
        and _inRange(abs(fit.modelStartTime), _ABS_MODEL_START_TIME_RANGE)
        and _inRange(displacementAt0p9s, _DISPLACEMENT_AT_0P9S_RANGE)
    )


def _readMotion(profile: ShutterMotionProfile, isOpen: bool) -> tuple[float, int, _ThreeJerksParams]:
    """Read one blade motion from its profile.

    Parameters
    ----------
    profile : `lsst.ip.isr.ShutterMotionProfile`
        The profile of the motion.
    isOpen : `bool`
        Whether the motion opens the shutter.

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
    ValueError, TypeError
        Raised if the side is not PLUSX or MINUSX, the model is not
        ThreeJerksModelv1, or a value is not a number.
    """
    side = profile.metadata["side"]
    if side not in ("PLUSX", "MINUSX"):
        raise ValueError(f"shutter side is {side!r}, not PLUSX or MINUSX")
    if profile.fit_model != _THREE_JERKS_V1:
        raise ValueError(f"shutter motion model is {profile.fit_model!r}, not {_THREE_JERKS_V1}")

    # The Hall-sensor fit.
    i = profile.fit_name.index(_HALL_SENSOR_FIT)
    fit = _ThreeJerksParams(
        modelStartTime=float(profile.fit_start_time[i]),
        pivot1=float(profile.fit_pivot1[i]),
        pivot2=float(profile.fit_pivot2[i]),
        jerk0=float(profile.fit_jerk0[i]),
        jerk1=float(profile.fit_jerk1[i]),
        jerk2=float(profile.fit_jerk2[i]),
    )

    startTime = float(profile.metadata["startTime_mjd"])
    decreasing = (side == "PLUSX") == isOpen

    return startTime, -1 if decreasing else 1, fit


def _beamQuadrature(offsets: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute the quadrature nodes of a beam's cumulative flux profile.

    Parameters
    ----------
    offsets : `numpy.ndarray`
        Blade positions (mm, from the beam centre) at the levels
        ``_BEAM_LEVELS``.

    Returns
    -------
    nodes : `numpy.ndarray`
        Blade positions (mm, from the beam centre) of the 2-point
        Gauss-Legendre nodes in each interval of flux fraction, including the
        tails, shape (20,).
    weights : `numpy.ndarray`
        Their weights, shape (20,).
    """
    # Blade positions at flux fractions 0, the levels, and 1 (linear tails).
    s = np.concatenate([[2*offsets[0] - offsets[1]], offsets, [2*offsets[-1] - offsets[-2]]])
    knots = np.concatenate([[0.0], _BEAM_LEVELS, [1.0]])
    du = np.diff(knots)

    # Two nodes per interval, interpolated linearly in flux fraction.
    k = np.repeat(np.arange(du.size), 2)
    f = np.tile(_GL2, du.size)
    nodes = s[k]*(1 - f) + s[k + 1]*f

    return nodes, np.repeat(du/2, 2)


def _compute(
    profile: ShutterMotionProfileFull,
    begin: float,
    geometry: _DetectorGeometry,
    config: ShutterTimingConfig,
    exposureTime: float,
) -> dict:
    """Compute the data-dependent fields of a `ShutterTiming`.

    Parameters
    ----------
    profile : `lsst.ip.isr.ShutterMotionProfileFull`
        The two blade motions.
    begin : `float`
        The exposure start time in the header, MJD-BEG (MJD TAI).
    geometry : ``_DetectorGeometry``
        The detector.
    config : `ShutterTimingConfig`
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
    motions = (_readMotion(profile.profile_open, True), _readMotion(profile.profile_close, False))
    (openStart, _, _), (closeStart, _, _) = motions

    if not all(_fitInRange(fit) for _, _, fit in motions):
        raise _Unavailable(ShutterTimingFlag.PARAM_RANGE, "a Hall fit is outside the nominal ranges")

    # Zero point: re-anchor to the header if the two shutter clocks disagree.
    flags = ShutterTimingFlag.NONE
    closeRel = (closeStart - openStart)*_SECONDS_PER_DAY
    if not _inRange((closeRel - exposureTime)*1e3, _CLOSE_MINUS_OPEN_MINUS_EXPOSURE_TIME_RANGE):
        flags |= ShutterTimingFlag.CLOCK_CLOSE_VS_OPEN
        openStart = begin + _HEADER_ANCHOR_OFFSETS[0]/_SECONDS_PER_DAY
        closeRel = exposureTime + _HEADER_ANCHOR_OFFSETS[1]

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

    # Their position along the blade travel (CCS x), plus the focal-plane
    # centre.
    xc, _ = geometry.pixelToCcs(*((pa, pc) if ia == 0 else (pc, pa)))
    xc = np.append(xc, 0.0)

    # Blade positions of each point's beam quadrature nodes.
    nodes, w = _beamQuadrature(np.array(config.beamOffsets))
    s = config.beamScale*xc[:, None] + nodes

    # Crossing times of each blade edge at the quadrature nodes (s after the
    # open start time).
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
                           "evaluation points (overlapping blades)")

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
