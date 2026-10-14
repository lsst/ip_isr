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
flux-weighted mid-exposure time at any of its pixels.

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
3. *Interpolation.*  The time therefore depends only on the position along
   the blade travel, focal-plane y.  It is computed every
   ``_FOCAL_PLANE_Y_STEP`` mm over the detector and interpolated linearly.

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

import lsst.pex.config as pexConfig
from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS, Detector
from lsst.afw.geom import TransformPoint2ToPoint2
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

#: Spacing (mm) of the focal-plane y at which times are computed; linear
#: interpolation between them is good to ~1 us.
_FOCAL_PLANE_Y_STEP = 0.5

#: Margin (mm) beyond the detector edges over which times are computed; NaN
#: beyond.
_FOCAL_PLANE_Y_MARGIN = 5.0


class ShutterTimingStatus(enum.IntEnum):
    """Whether a detector has a corrected time."""

    OK = 0
    """Corrected."""
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


@dataclasses.dataclass(frozen=True, eq=False)
class ShutterTiming:
    """Shutter-corrected mid-exposure times of one detector of one exposure.

    Times are MJD TAI, durations seconds.  The mid-exposure time depends only
    on focal-plane y; it is tabulated over the detector and interpolated
    linearly.
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

    focalPlaneY: np.ndarray | None = dataclasses.field(default=None, repr=False)
    """Focal-plane y (mm, increasing) at which the times are tabulated."""

    midpointOffsets: np.ndarray | None = dataclasses.field(default=None, repr=False)
    """Mid-exposure times at ``focalPlaneY``, s after
    ``centerMidpointMjdTai``.
    """

    pixelToFocalPlane: TransformPoint2ToPoint2 | None = dataclasses.field(default=None, repr=False)
    """The detector's pixel to focal-plane (mm) transform."""

    def midpointMjdTai(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Compute per-source mid-exposure times at pixel positions.

        Parameters
        ----------
        x, y : `numpy.ndarray`
            Pixel coordinates, of the same shape.

        Returns
        -------
        tMid : `numpy.ndarray`
            Mid-exposure times (MJD TAI), of the shape of ``x``.  NaN
            everywhere if the detector is UNAVAILABLE, and NaN where x or y is
            NaN or more than 5 mm (500 pixels) beyond the detector.
        """
        if self.status == ShutterTimingStatus.UNAVAILABLE:
            return np.full(np.shape(x), np.nan)

        # Focal-plane y of each position.
        pixels = np.array([np.ravel(x), np.ravel(y)], dtype=float)
        focalPlaneY = self.pixelToFocalPlane.getMapping().applyForward(pixels)[1]

        # Interpolated time, in seconds from the centre time.
        dt = np.interp(focalPlaneY, self.focalPlaneY, self.midpointOffsets, left=np.nan, right=np.nan)

        return np.reshape(self.centerMidpointMjdTai + dt/_SECONDS_PER_DAY, np.shape(x))


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
        - there is no time at the focal-plane centre or somewhere on the
          detector (overlapping blades in very short exposures).

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
      each interval between levels.
    """
    if config is None:
        config = ShutterTimingConfig()

    common = dict(detectorId=detector.getId())

    # Exposures before the shutter cards were added.
    if "SHUTTER OPEN STARTTIME TAI MJD" not in metadata:
        return ShutterTiming(flags=ShutterTimingFlag.NO_PROFILE, message="no shutter motion cards", **common)

    profile = ShutterMotionProfileFull.fromMetadata(metadata)
    begin = float(metadata["MJD-BEG"])
    common["headerMidpointMjdTai"] = 0.5*(begin + float(metadata["MJD-END"]))

    # Focal-plane y of the detector centre and of the table, which spans the
    # detector plus a margin.
    pixelToFocalPlane = detector.getTransform(PIXELS, FOCAL_PLANE)
    centerY = detector.getCenter(FOCAL_PLANE).getY()
    cornersY = [p.getY() for p in detector.getCorners(FOCAL_PLANE)]
    lo = min(cornersY) - _FOCAL_PLANE_Y_MARGIN
    hi = max(cornersY) + _FOCAL_PLANE_Y_MARGIN
    focalPlaneY = np.linspace(lo, hi, int(np.ceil((hi - lo)/_FOCAL_PLANE_Y_STEP)) + 1)

    # Times at the detector centre, the focal-plane centre and the table.
    try:
        flags, openStart, tMid, tEff = _compute(profile, begin, config, exposureTime,
                                                np.concatenate([[centerY, 0.0], focalPlaneY]))
    except _Unavailable as e:
        return ShutterTiming(flags=e.flags, message=e.message, **common)

    return ShutterTiming(
        status=ShutterTimingStatus.OK,
        flags=flags,
        centerMidpointMjdTai=float(openStart + tMid[0]/_SECONDS_PER_DAY),
        focalPlaneMidpointMjdTai=float(openStart + tMid[1]/_SECONDS_PER_DAY),
        effectiveExposureTime=float(tEff[0]),
        focalPlaneY=focalPlaneY,
        midpointOffsets=tMid[2:] - tMid[0],
        pixelToFocalPlane=pixelToFocalPlane,
        **common,
    )


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
    config: ShutterTimingConfig,
    exposureTime: float,
    focalPlaneY: np.ndarray,
) -> tuple[ShutterTimingFlag, float, np.ndarray, np.ndarray]:
    """Compute mid-exposure and effective exposure times along the blade
    travel.

    Parameters
    ----------
    profile : `lsst.ip.isr.ShutterMotionProfileFull`
        The two blade motions.
    begin : `float`
        The exposure start time in the header, MJD-BEG (MJD TAI).
    config : `ShutterTimingConfig`
        The beam model.
    exposureTime : `float`
        The requested exposure time (s).
    focalPlaneY : `numpy.ndarray`
        Focal-plane y (mm) of the points.

    Returns
    -------
    flags : `ShutterTimingFlag`
        CLOCK_CLOSE_VS_OPEN if the start times were re-anchored.
    openStart : `float`
        Start time of the opening motion (MJD TAI).
    tMid : `numpy.ndarray`
        Flux-weighted mid-exposure times at the points (s after
        ``openStart``).
    tEff : `numpy.ndarray`
        Flux-weighted open times at the points (s).

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

    # Blade positions of each point's beam quadrature nodes (focal-plane y is
    # CCS x, the blade coordinate).
    nodes, w = _beamQuadrature(np.array(config.beamOffsets))
    s = config.beamScale*focalPlaneY[:, None] + nodes

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

    return flags, openStart, tMid, eC - eO
