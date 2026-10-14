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
   the blade travel, focal-plane y.  It is computed every 0.5 mm over the
   detector and interpolated linearly.

Times are TAI MJD.  The shutter motion and its Hall-sensor fit are described
in CTN-002, *Camera Shutter Motion Analysis* (https://ctn-002.lsst.io).
"""

__all__ = [
    "ShutterTiming",
    "ShutterTimingConfig",
    "computeShutterTiming",
]

import dataclasses
import logging

import numpy as np

import lsst.pex.config as pexConfig
from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS
from lsst.afw.geom import TransformPoint2ToPoint2
from lsst.afw.image import Exposure

from .shutterMotion import ShutterMotionProfile, ShutterMotionProfileFull

log = logging.getLogger(__name__)

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
#: the shutter timestamps disagree and both start times are re-anchored to
#: the header.
_CLOSE_MINUS_OPEN_MINUS_EXPOSURE_TIME_RANGE = (-2.0, 4.0)

#: Time step (s) and span (s) of the model times at which blade motions are
#: evaluated for inversion; a motion lasts ~0.9 s.
_MOTION_TIME_STEP = 1e-4
_MAX_MOTION_TIME = 1.5

#: Re-anchoring offsets (s): open STARTTIME - MJD-BEG, and close - open
#: STARTTIME - requested exposure time.
_HEADER_ANCHOR_OFFSETS = (0.00809, 0.00052)

#: Spacing (mm) of the focal-plane y at which times are computed; linear
#: interpolation between them is good to ~1 us.
_FOCAL_PLANE_Y_STEP = 0.5

#: Margin (mm) beyond the detector edges over which times are computed; NaN
#: beyond.
_FOCAL_PLANE_Y_MARGIN = 5.0


class ShutterTimingConfig(pexConfig.Config):
    """Configuration of `computeShutterTiming`.

    The defaults are for LSSTCam: a least-squares fit to the LSST Camera
    raytrace of the beam at the shutter plane (LCA-20578, rotator angle 0),
    on its 137 field positions within 280 mm of the focal-plane centre.
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


def _displacement(fit: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Compute the displacement of a ThreeJerksModelv1 blade motion.

    The camera fits each blade motion, from its Hall-sensor positions, with
    this model (CTN-002): piecewise-constant jerk ``jerk0`` before
    ``pivot1``, ``jerk1`` to ``pivot2`` and ``jerk2`` after, from rest at
    model time 0 (the motion start time plus ``modelStartTime``).

    Parameters
    ----------
    fit : `numpy.ndarray`
        The fit: modelStartTime (s), pivot1 (s), pivot2 (s), jerk0, jerk1 and
        jerk2 (mm/s^3).
    t : `numpy.ndarray`
        Model times (s), >= 0.

    Returns
    -------
    d : `numpy.ndarray`
        Displacement (mm) along the travel direction.
    """
    # Start time and jerk of each segment.
    starts = np.array([0.0, fit[1], fit[2]])
    jerks = fit[3:]

    # Displacement, velocity and acceleration at the start of each segment.
    state = np.zeros((3, 3))
    for k in (1, 2):
        s0, v0, a0 = state[k - 1]
        h, j = starts[k] - starts[k - 1], jerks[k - 1]
        state[k] = (s0 + h*(v0 + h*(a0/2 + h*j/6)), v0 + h*(a0 + h*j/2), a0 + h*j)

    # Constant-jerk motion within each time's segment.
    k = np.searchsorted(starts, t, side="right") - 1
    h = t - starts[k]
    s0, v0, a0 = state[k].T

    return s0 + h*(v0 + h*(a0/2 + h*jerks[k]/6))


def _timeOfDisplacement(fit: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Invert a ThreeJerksModelv1 blade motion.

    Parameters
    ----------
    fit : `numpy.ndarray`
        The fit, as for ``_displacement``.
    d : `numpy.ndarray`
        Displacements (mm) along the travel direction.

    Returns
    -------
    t : `numpy.ndarray`
        Model times (s) at which the blade reaches ``d``, interpolated
        linearly between times ``_MOTION_TIME_STEP`` apart (to < 0.1 us); NaN
        before the start and beyond the end of the forward motion.
    """
    # The displacement on a time grid, up to where the motion stops.
    tGrid = np.arange(0.0, _MAX_MOTION_TIME + _MOTION_TIME_STEP/2, _MOTION_TIME_STEP)
    dGrid = _displacement(fit, tGrid)
    forward = np.flatnonzero(np.diff(dGrid) <= 0)
    end = forward[0] + 1 if forward.size else dGrid.size

    return np.interp(d, dGrid[:end], tGrid[:end], left=np.nan, right=np.nan)


@dataclasses.dataclass(frozen=True, eq=False)
class ShutterTiming:
    """Shutter-corrected mid-exposure times of one detector of one exposure.

    Times are MJD TAI, durations seconds.  The mid-exposure time depends only
    on focal-plane y; it is tabulated over the detector and interpolated
    linearly.
    """

    centerMidpointMjdTai: float
    """Flux-weighted mid-exposure time at the detector centre."""

    focalPlaneY: np.ndarray = dataclasses.field(repr=False)
    """Focal-plane y (mm, increasing) at which the times are tabulated."""

    midpointOffsets: np.ndarray = dataclasses.field(repr=False)
    """Mid-exposure times at ``focalPlaneY``, s after
    ``centerMidpointMjdTai``.
    """

    pixelToFocalPlane: TransformPoint2ToPoint2 = dataclasses.field(repr=False)
    """The detector's pixel to focal-plane (mm) transform."""

    reanchored: bool = False
    """Whether the shutter timestamps disagreed with the requested exposure
    time, so that the start times were re-anchored to ``MJD-BEG``.
    """

    def midpointMjdTai(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Compute per-source mid-exposure times at pixel positions.

        Parameters
        ----------
        x, y : `numpy.ndarray`
            Pixel coordinates, of the same shape.

        Returns
        -------
        tMid : `numpy.ndarray`
            Mid-exposure times (MJD TAI), of the shape of ``x``; NaN where x
            or y is NaN or more than 5 mm (500 pixels) beyond the detector.
        """
        # Focal-plane y of each position.
        pixels = np.array([np.ravel(x), np.ravel(y)], dtype=float)
        focalPlaneY = self.pixelToFocalPlane.getMapping().applyForward(pixels)[1]

        # Interpolated time, in seconds from the centre time.
        dt = np.interp(focalPlaneY, self.focalPlaneY, self.midpointOffsets, left=np.nan, right=np.nan)

        return np.reshape(self.centerMidpointMjdTai + dt/_SECONDS_PER_DAY, np.shape(x))


def computeShutterTiming(
    exposure: Exposure,
    requestedExposureTime: float | None = None,
    config: ShutterTimingConfig | None = None,
) -> ShutterTiming | None:
    """Compute the shutter-corrected times of one detector of one exposure.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The exposure; its metadata carries the shutter cards and ``MJD-BEG``.
    requestedExposureTime : `float`, optional
        The exposure time requested for this exposure (s), from the butler
        ``visit`` dimension record (``exposure_time``).  It is an independent
        check on the two shutter timestamps in the header: if they are not
        this far apart, one of them is wrong, and both start times are
        re-anchored to ``MJD-BEG``.  (``EXPTIME`` is not in the exposure
        metadata, and ``visitInfo.exposureTime`` is measured from the same
        timestamps, so neither can serve.)  If `None`, the timestamps are
        used unchecked.
    config : `ShutterTimingConfig`, optional
        Defaults to ``ShutterTimingConfig()``.

    Returns
    -------
    timing : `ShutterTiming` or `None`
        The times, or `None` (callers keep the header midpoint; the reason is
        logged) if:

        - there are no shutter cards (``SHUTTER OPEN STARTTIME TAI MJD``
          absent; exposures before 2025-10-29);
        - either Hall fit is outside the nominal ranges;
        - a blade does not reach part of the detector's beam (not seen in
          data that passes the fit checks).

    Raises
    ------
    ValueError
        Raised if a SIDE card is not PLUSX or MINUSX, or a MODEL card is not
        ThreeJerksModelv1.
    RuntimeError
        Raised if the shutter cards are present but one of them is missing.
    KeyError, ValueError, TypeError
        Raised if MJD-BEG is missing, or a card is not a number.

    Notes
    -----
    Rules beyond the method in the module docstring:

    - Blade start positions are nominal: 750.76 mm for 750 -> 0 moves, -0.05
      mm for 0 -> 750 moves.
    - If ``requestedExposureTime`` is given and close - open STARTTIME -
      ``requestedExposureTime`` is outside -2 to 4 ms, the shutter timestamps
      disagree: the open start time is set from MJD-BEG and the close start
      time from ``requestedExposureTime``.
    """
    if config is None:
        config = ShutterTimingConfig()

    metadata = exposure.metadata
    detector = exposure.detector

    # Exposures before the shutter cards were added.
    if "SHUTTER OPEN STARTTIME TAI MJD" not in metadata:
        log.info("No shutter motion cards (exposures before 2025-10-29): no shutter-corrected times.")
        return None

    # The two blade motions; none if a Hall fit is out of range.
    profile = ShutterMotionProfileFull.fromExposure(exposure)
    motions = (_readMotion(profile.profile_open, True), _readMotion(profile.profile_close, False))
    if not all(_fitInRange(fit) for _, _, fit in motions):
        log.warning("A shutter Hall-sensor fit is outside the nominal ranges: no shutter-corrected times.")
        return None

    # Zero point: re-anchor to the header if the shutter timestamps disagree
    # with the requested exposure time.
    (openStart, _, _), (closeStart, _, _) = motions
    closeRel = (closeStart - openStart)*_SECONDS_PER_DAY
    reanchored = requestedExposureTime is not None and not _inRange(
        (closeRel - requestedExposureTime)*1e3, _CLOSE_MINUS_OPEN_MINUS_EXPOSURE_TIME_RANGE
    )
    if reanchored:
        log.warning("Shutter timestamps %.4f s apart for a requested exposure time of %.4f s: "
                    "start times re-anchored to MJD-BEG.", closeRel, requestedExposureTime)
        openStart = float(metadata["MJD-BEG"]) + _HEADER_ANCHOR_OFFSETS[0]/_SECONDS_PER_DAY
        closeRel = requestedExposureTime + _HEADER_ANCHOR_OFFSETS[1]

    # Focal-plane y of the detector centre and of the table, which spans the
    # detector plus a margin.
    pixelToFocalPlane = detector.getTransform(PIXELS, FOCAL_PLANE)
    centerY = detector.getCenter(FOCAL_PLANE).getY()
    cornersY = [p.getY() for p in detector.getCorners(FOCAL_PLANE)]
    lo = min(cornersY) - _FOCAL_PLANE_Y_MARGIN
    hi = max(cornersY) + _FOCAL_PLANE_Y_MARGIN
    focalPlaneY = np.linspace(lo, hi, int(np.ceil((hi - lo)/_FOCAL_PLANE_Y_STEP)) + 1)

    # Times at the detector centre and on the table; none where a blade
    # does not reach the beam.
    tMid = _midpointTimes(motions, closeRel, config, np.concatenate([[centerY], focalPlaneY]))
    if not np.all(np.isfinite(tMid)):
        log.warning("A shutter blade does not reach part of the detector's beam: no shutter-corrected times.")
        return None

    return ShutterTiming(
        centerMidpointMjdTai=float(openStart + tMid[0]/_SECONDS_PER_DAY),
        focalPlaneY=focalPlaneY,
        midpointOffsets=tMid[1:] - tMid[0],
        pixelToFocalPlane=pixelToFocalPlane,
        reanchored=reanchored,
    )


def _inRange(value: float, bounds: tuple[float, float]) -> bool:
    """Return whether lo <= value <= hi (False for NaN)."""
    return bounds[0] <= value <= bounds[1]


def _fitInRange(fit: np.ndarray) -> bool:
    """Return whether a Hall fit (as for ``_displacement``) is within the
    nominal ranges.
    """
    modelStartTime, pivot1, pivot2 = fit[:3]
    displacementAt0p9s = _displacement(fit, np.array([0.9]))[0]

    return (
        0.0 <= pivot1 <= pivot2
        and _inRange(pivot1, _PIVOT1_RANGE)
        and _inRange(pivot2, _PIVOT2_RANGE)
        and all(_inRange(abs(j), _ABS_JERK_RANGE) for j in fit[3:])
        and _inRange(abs(modelStartTime), _ABS_MODEL_START_TIME_RANGE)
        and _inRange(displacementAt0p9s, _DISPLACEMENT_AT_0P9S_RANGE)
    )


def _readMotion(profile: ShutterMotionProfile, isOpen: bool) -> tuple[float, int, np.ndarray]:
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
    fit : `numpy.ndarray`
        The Hall-sensor fit, as for ``_displacement``.

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
    fit = np.array([profile.fit_start_time[i], profile.fit_pivot1[i], profile.fit_pivot2[i],
                    profile.fit_jerk0[i], profile.fit_jerk1[i], profile.fit_jerk2[i]], dtype=float)

    startTime = float(profile.metadata["startTime_mjd"])
    decreasing = (side == "PLUSX") == isOpen

    return startTime, -1 if decreasing else 1, fit


def _beamQuadrature(offsets: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute the quadrature nodes of a beam's cumulative flux profile.

    The profile is piecewise linear through the levels ``_BEAM_LEVELS``, with
    linear tails reaching 0 and 1 one interval beyond the first and last
    level; each interval gets 2-point Gauss-Legendre nodes.

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


def _midpointTimes(
    motions: tuple[tuple[float, int, np.ndarray], tuple[float, int, np.ndarray]],
    closeRel: float,
    config: ShutterTimingConfig,
    focalPlaneY: np.ndarray,
) -> np.ndarray:
    """Compute flux-weighted mid-exposure times along the blade travel.

    Parameters
    ----------
    motions : `tuple`
        The opening and closing motions, as returned by ``_readMotion``.
    closeRel : `float`
        Start time of the closing motion, s after the opening one.
    config : `ShutterTimingConfig`
        The beam model.
    focalPlaneY : `numpy.ndarray`
        Focal-plane y (mm) of the points.

    Returns
    -------
    tMid : `numpy.ndarray`
        Mid-exposure times at the points, s after the start of the opening
        motion; NaN where a blade does not reach the beam.
    """
    # Blade positions of each point's beam quadrature nodes (focal-plane y is
    # CCS x, the blade coordinate).
    nodes, w = _beamQuadrature(np.array(config.beamOffsets))
    s = config.beamScale*focalPlaneY[:, None] + nodes

    # Crossing times of each blade edge at the quadrature nodes, s after the
    # start of its motion.
    tau = []
    for _, travelSign, fit in motions:
        start = _NOMINAL_STROKE if travelSign < 0 else _NOMINAL_START_INCREASING
        d = travelSign*((_ENCODER_CENTER - start) - s)
        tau.append(_timeOfDisplacement(fit, d) + fit[0])
    tO, tauC = tau
    tC = tauC + closeRel

    # Flux-weighted times of the two edges, and their mean.
    return 0.5*(tO @ w + tC @ w)
