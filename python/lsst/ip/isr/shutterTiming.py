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

import numpy as np

import lsst.pex.config as pexConfig


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
    """Missing or unusable shutter cards (no Hall fit, or a model other than
    ``ThreeJerksModelv1``).  Implies `ShutterTimingStatus.UNAVAILABLE`.
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
        raise NotImplementedError("WP-H")

    def pixelToCcs(self, x, y, ccsFromDvcsSign: int = 1) -> tuple[np.ndarray, np.ndarray]:
        """Pixel -> CCS (mm): ``x_ccs = s * Y_dvcs``, ``y_ccs = s * X_dvcs``.
        """
        raise NotImplementedError("WP-H")

    def offDetector(self, x, y) -> np.ndarray:
        """How far (pixels; the larger axis) each position lies outside the
        pixel bounds; 0 inside.
        """
        raise NotImplementedError("WP-H")


class ShutterBeamModel:
    """Shutter-plane beam model: for each field position (CCS mm), the blade
    coordinate s_q (CCS-x mm) below which a fraction q of the pixel's beam flux
    lies, at the tabulated levels q (0.01 ... 0.99); interpolated in the field
    position, and extrapolated (with a hull-distance test) outside the table's
    coverage.
    """

    @classmethod
    def fromFile(cls, path: str) -> ShutterBeamModel:
        """Read a beam table in the raytrace ``.tnt`` format."""
        raise NotImplementedError("WP-H")

    @property
    def levels(self) -> np.ndarray:
        """The flux levels q."""
        raise NotImplementedError("WP-H")

    def hullDistance(self, xCcs, yCcs) -> np.ndarray:
        """Distance (mm) outside the convex hull of the table nodes; 0 inside.
        """
        raise NotImplementedError("WP-H")


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
        raise NotImplementedError("WP-H")

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
        raise NotImplementedError("WP-H")

    def summary(self) -> dict:
        """Scalars for task metadata: status (name), flags (int), message,
        centerMjdTai, focalPlaneMjdTai, centerMinusHeaderMid (s),
        maxAbsResidual, policy.
        """
        raise NotImplementedError("WP-H")


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

    1. Parse the open and close motions; a missing or malformed card ->
       NO_PROFILE, UNAVAILABLE.  Start positions: ``nominalStroke`` for 750 ->
       0 moves, ``nominalStartIncreasing`` for 0 -> 750 moves.
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
    raise NotImplementedError("WP-H")
