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
the focal plane).  For one detector, `computeShutterTiming` gives a compact
model of the flux-weighted mid-exposure time at any of its pixels,
`ShutterTiming`, which ISR records in the image metadata and
`getShutterTiming` returns.

Method:

1. *Blade trajectories.*  The camera fits each blade motion with the
   ``ThreeJerksModelv1`` model; the fit is in its motion profiles, and in the
   ``SHUTTER {OPEN,CLOSE} ...`` metadata cards (since 2025-10-29).
   `lsst.ip.isr.ShutterMotionProfile` gives the time at which a blade edge
   reaches a given position.
2. *Beam-weighted crossing times.*  A pixel's converging ray bundle is
   ~67 mm wide at the shutter plane, so an edge uncovers (or covers) it
   gradually.  The flux-weighted crossing time of one edge is
   ``<T> = integral T(s) dF(s)``, with ``F`` the cumulative flux profile of
   the bundle along the blade travel.  The profile has the same shape
   everywhere on the focal plane and moves with the pixel: its levels lie
   at ``beamScale * x + beamOffsets`` (`ShutterTimingConfig`), a linear fit
   to the LSST Camera raytrace of the beam (LCA-20578).  The mid-exposure
   time is ``(<T>_open + <T>_close) / 2``.
3. *Chebyshev series.*  The time therefore depends only on the position
   along the blade travel, focal-plane y, which is pixel y on every LSSTCam
   detector.  Over one detector it is a degree-5 Chebyshev series in pixel
   y, fitted to the exact times (to ~2 us).

Times are TAI MJD.  The shutter motion and its Hall-sensor fit are described
in CTN-002, *Camera Shutter Motion Analysis* (https://ctn-002.lsst.io).
"""

__all__ = [
    "ShutterTiming",
    "ShutterTimingConfig",
    "computeShutterTiming",
    "getShutterTiming",
]

import dataclasses
import logging

import numpy as np

import lsst.pex.config as pexConfig
from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS, Detector
from lsst.afw.image import Exposure, ExposureSummaryStats
from lsst.daf.base import PropertyList

from .shutterMotion import ShutterMotionProfile, ShutterMotionProfileFull

log = logging.getLogger(__name__)

_SECONDS_PER_DAY = 86400.0

#: Cumulative flux levels q of ``ShutterTimingConfig.beamOffsets``.
_BEAM_LEVELS = np.array([0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99])

#: 2-point Gauss-Legendre abscissae on [0, 1].
_GL2 = np.array([0.5 - 0.5 / np.sqrt(3.0), 0.5 + 0.5 / np.sqrt(3.0)])

#: Encoder position (mm) of the focal-plane centre: x_ccs = c - e.
_ENCODER_CENTER = 375.0

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

#: Metadata card with the requested exposure time (s), added by the obs_lsst
#: raw formatter as a workaround until VisitInfo records it.
_REQUESTED_EXPOSURE_TIME_CARD = "EXPTIME REQUESTED"

#: Requested exposure times (s) of the exposures whose shutter timestamps are
#: known to disagree, for images made before the card above was added: all
#: such LSSTCam exposures from 2025-10-29 to 2026-07-14.  Remove once those
#: images are no longer processed.
_KNOWN_CLOCK_FAULT_EXPOSURES = {2026031900172: 30.0, 2026031900173: 30.0, 2026052200119: 30.0}

#: Re-anchoring offsets (s): open STARTTIME - MJD-BEG, and close - open
#: STARTTIME - requested exposure time.
_HEADER_ANCHOR_OFFSETS = (0.00809, 0.00052)

#: Degree of the Chebyshev series in pixel y, and the number of Chebyshev
#: points over the detector to which it is fitted.
_CHEBYSHEV_DEGREE = 5
_CHEBYSHEV_FIT_POINTS = 16

#: Metadata cards of a `ShutterTiming`.
_MIDPOINT_CARD = "SHUTTER MIDPT MJD"
_CHEBYSHEV_CARD = "SHUTTER CHEB{}"
_CHEBYSHEV_Y_CENTER_CARD = "SHUTTER CHEB YCENTER"
_CHEBYSHEV_Y_HALF_RANGE_CARD = "SHUTTER CHEB YHALF"
_EXPOSURE_TIME_CARD = "SHUTTER EXPOSURE TIME"
_REANCHORED_CARD = "SHUTTER REANCHORED"
_SOURCE_CARD = "SHUTTER TIMING SOURCE"


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


@dataclasses.dataclass(frozen=True, eq=False)
class ShutterTiming:
    """Shutter-corrected mid-exposure times of one detector of one exposure.

    The mid-exposure time at pixel row ``y`` is
    ``centerMidpointMjdTai + sum(c[k] T_k(u))/86400``, with ``T_k`` the
    Chebyshev polynomials, ``c`` the ``chebyshevCoefficients`` and
    ``u = (y - chebyshevYCenter)/chebyshevYHalfRange``, which runs from -1 to
    1 across the detector.
    """

    centerMidpointMjdTai: float
    """Flux-weighted mid-exposure time at the detector centre (MJD TAI)."""

    chebyshevCoefficients: np.ndarray
    """Chebyshev coefficients (s) of the mid-exposure time relative to
    ``centerMidpointMjdTai``, in increasing degree.
    """

    chebyshevYCenter: float
    """Pixel y at which the Chebyshev variable is 0: the detector's central
    row.
    """

    chebyshevYHalfRange: float
    """Pixels per unit of the Chebyshev variable: half the detector height.
    """

    centerExposureTime: float
    """Flux-weighted time (s) between the opening and closing of the shutter
    at the detector centre.
    """

    reanchored: bool | None
    """Whether the shutter timestamps disagreed with the requested exposure
    time, so that the start times were re-anchored to ``MJD-BEG``; `None` if
    read from summary statistics, which do not record it.
    """

    source: str | None
    """Source of the blade motions: ``"profile"`` (motion profiles with the
    measured start positions) or ``"header"`` (exposure header cards); `None`
    if read from summary statistics, which do not record it.
    """

    def midpointMjdTai(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Compute per-source mid-exposure times at pixel positions.

        Parameters
        ----------
        x, y : `numpy.ndarray`
            Pixel coordinates, of the same shape; the times depend only on
            ``y``.

        Returns
        -------
        tMid : `numpy.ndarray`
            Mid-exposure times (MJD TAI), of the shape of ``y``.
        """
        u = (y - self.chebyshevYCenter)/self.chebyshevYHalfRange
        dt = np.polynomial.chebyshev.chebval(u, self.chebyshevCoefficients)

        return self.centerMidpointMjdTai + dt/_SECONDS_PER_DAY

    def toMetadata(self, metadata: PropertyList) -> None:
        """Write the times to metadata cards.

        Parameters
        ----------
        metadata : `lsst.daf.base.PropertyList`
            Metadata to update.
        """
        metadata.set(_MIDPOINT_CARD, self.centerMidpointMjdTai,
                     "Shutter-corrected mid-exposure time at detector centre (MJD TAI)")
        for k, c in enumerate(self.chebyshevCoefficients):
            metadata.set(_CHEBYSHEV_CARD.format(k), c,
                         f"Chebyshev coefficient {k} (s) of mid-exposure time in pixel y")
        metadata.set(_CHEBYSHEV_Y_CENTER_CARD, self.chebyshevYCenter,
                     "Pixel y at which the Chebyshev variable is 0")
        metadata.set(_CHEBYSHEV_Y_HALF_RANGE_CARD, self.chebyshevYHalfRange,
                     "Pixels per unit of the Chebyshev variable")
        metadata.set(_EXPOSURE_TIME_CARD, self.centerExposureTime,
                     "Shutter open time (s) at detector centre")
        metadata.set(_REANCHORED_CARD, self.reanchored, "Shutter start times re-anchored to MJD-BEG")
        metadata.set(_SOURCE_CARD, self.source, "Source of the shutter motions: profile or header")

    @classmethod
    def fromMetadata(cls, metadata: PropertyList) -> "ShutterTiming | None":
        """Read the times from metadata cards.

        Parameters
        ----------
        metadata : `lsst.daf.base.PropertyList`
            Metadata written by `toMetadata`.

        Returns
        -------
        timing : `ShutterTiming` or `None`
            The times, or `None` if the metadata has none (images that ISR
            processed without shutter timing).
        """
        if _MIDPOINT_CARD not in metadata:
            return None

        return cls(
            centerMidpointMjdTai=metadata[_MIDPOINT_CARD],
            chebyshevCoefficients=np.array(
                [metadata[_CHEBYSHEV_CARD.format(k)] for k in range(_CHEBYSHEV_DEGREE + 1)]
            ),
            chebyshevYCenter=metadata[_CHEBYSHEV_Y_CENTER_CARD],
            chebyshevYHalfRange=metadata[_CHEBYSHEV_Y_HALF_RANGE_CARD],
            centerExposureTime=metadata[_EXPOSURE_TIME_CARD],
            reanchored=metadata[_REANCHORED_CARD],
            source=metadata[_SOURCE_CARD],
        )

    def updateSummaryStats(self, summary: ExposureSummaryStats) -> None:
        """Write the times to summary statistics.

        Parameters
        ----------
        summary : `lsst.afw.image.ExposureSummaryStats`
            Summary statistics to update.
        """
        summary.shutterMidpointMjdTai = self.centerMidpointMjdTai
        summary.shutterChebyshevCoefficients = self.chebyshevCoefficients.tolist()
        summary.shutterChebyshevYCenter = self.chebyshevYCenter
        summary.shutterChebyshevYHalfRange = self.chebyshevYHalfRange
        summary.shutterExposureTime = self.centerExposureTime

    @classmethod
    def fromSummaryStats(cls, summary: ExposureSummaryStats) -> "ShutterTiming | None":
        """Read the times from summary statistics.

        Parameters
        ----------
        summary : `lsst.afw.image.ExposureSummaryStats`
            Summary statistics written by `updateSummaryStats`.

        Returns
        -------
        timing : `ShutterTiming` or `None`
            The times, or `None` if the summary has none.
        """
        if np.isnan(summary.shutterMidpointMjdTai):
            return None

        return cls(
            centerMidpointMjdTai=summary.shutterMidpointMjdTai,
            chebyshevCoefficients=np.array(summary.shutterChebyshevCoefficients),
            chebyshevYCenter=summary.shutterChebyshevYCenter,
            chebyshevYHalfRange=summary.shutterChebyshevYHalfRange,
            centerExposureTime=summary.shutterExposureTime,
            reanchored=None,
            source=None,
        )


def computeShutterTiming(
    profileOpen: ShutterMotionProfile,
    profileClose: ShutterMotionProfile,
    detector: Detector,
    mjdBeg: float,
    requestedExposureTime: float | None,
    config: ShutterTimingConfig | None = None,
) -> ShutterTiming | None:
    """Compute the shutter-corrected times of one detector of one exposure.

    Parameters
    ----------
    profileOpen, profileClose : `lsst.ip.isr.ShutterMotionProfile`
        The opening and closing blade motions, from the camera's motion
        profiles or from the exposure header.
    detector : `lsst.afw.cameraGeom.Detector`
        The detector.
    mjdBeg : `float`
        The exposure's ``MJD-BEG`` (MJD TAI).
    requestedExposureTime : `float` or `None`
        The requested exposure time (s), used only to check the shutter
        timestamps; `None` to use them unchecked.
    config : `ShutterTimingConfig`, optional
        Defaults to ``ShutterTimingConfig()``.

    Returns
    -------
    timing : `ShutterTiming` or `None`
        The times, or `None` (callers keep the header midpoint; the reason is
        logged) if:

        - either Hall fit is outside the nominal ranges;
        - a blade does not reach part of the detector's beam (not seen in
          data that passes the fit checks).

    Raises
    ------
    ValueError
        Raised if a side is not PLUSX or MINUSX, or a model is not
        ThreeJerksModelv1.

    Notes
    -----
    Rules beyond the method in the module docstring:

    - Blade start positions are the measured ones of motion profiles; the
      header has none, and they are nominal: 750.76 mm for 750 -> 0 moves,
      -0.05 mm for 0 -> 750 moves.
    - The requested exposure time is an independent check on the two shutter
      timestamps: if close - open start time differs from it by more than -2
      to 4 ms, one timestamp is wrong, and the open start time is set from
      ``MJD-BEG`` and the close start time from the requested exposure time.
      (``visitInfo.exposureTime`` is measured from the same timestamps, so
      it cannot serve.)
    """
    if config is None:
        config = ShutterTimingConfig()
    motions = (profileOpen, profileClose)

    # No times if a Hall fit is out of range.
    if not all(_fitInRange(motion) for motion in motions):
        log.warning("A shutter Hall-sensor fit is outside the nominal ranges: no shutter-corrected times.")
        return None

    # Zero point: re-anchor to the header if the shutter timestamps disagree
    # with the requested exposure time.
    openStart, closeStart = (float(motion.metadata["startTime_mjd"]) for motion in motions)
    closeRel = (closeStart - openStart)*_SECONDS_PER_DAY
    reanchored = requestedExposureTime is not None and not _inRange(
        (closeRel - requestedExposureTime)*1e3, _CLOSE_MINUS_OPEN_MINUS_EXPOSURE_TIME_RANGE
    )
    if reanchored:
        log.warning("Shutter timestamps %.4f s apart for a requested exposure time of %.4f s: "
                    "start times re-anchored to MJD-BEG.", closeRel, requestedExposureTime)
        openStart = mjdBeg + _HEADER_ANCHOR_OFFSETS[0]/_SECONDS_PER_DAY
        closeRel = requestedExposureTime + _HEADER_ANCHOR_OFFSETS[1]

    # The detector centre and Chebyshev points over the detector, along its
    # central column.
    bbox = detector.getBBox()
    xCenter = bbox.getMinX() + 0.5*(bbox.getWidth() - 1)
    yCenter = bbox.getMinY() + 0.5*(bbox.getHeight() - 1)
    yHalfRange = 0.5*bbox.getHeight()
    u = np.concatenate([[0.0], np.polynomial.chebyshev.chebpts1(_CHEBYSHEV_FIT_POINTS)])
    pixels = np.array([np.full(u.size, xCenter), yCenter + yHalfRange*u])
    focalPlaneY = detector.getTransform(PIXELS, FOCAL_PLANE).getMapping().applyForward(pixels)[1]

    # Opening and closing times there; none where a blade does not reach the
    # beam.
    tOpen, tClose = _edgeTimes(motions, closeRel, config, focalPlaneY)
    if not np.all(np.isfinite(tOpen + tClose)):
        log.warning("A shutter blade does not reach part of the detector's beam: no shutter-corrected times.")
        return None
    tMid = 0.5*(tOpen + tClose)

    # Mid-exposure time at the centre, and the series relative to it.
    coefficients = np.polynomial.chebyshev.chebfit(u[1:], tMid[1:] - tMid[0], _CHEBYSHEV_DEGREE)

    return ShutterTiming(
        centerMidpointMjdTai=openStart + tMid[0]/_SECONDS_PER_DAY,
        chebyshevCoefficients=coefficients,
        chebyshevYCenter=yCenter,
        chebyshevYHalfRange=yHalfRange,
        centerExposureTime=float(tClose[0] - tOpen[0]),
        reanchored=reanchored,
        source="profile" if "startPosition" in profileOpen.metadata else "header",
    )


def getShutterTiming(
    exposure: Exposure,
    config: ShutterTimingConfig | None = None,
) -> ShutterTiming | None:
    """Return the shutter-corrected times of an image.

    These are the times ISR recorded in the metadata; images that ISR
    processed without shutter timing get them computed from the header (remove
    this once those are no longer processed).

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The exposure; without recorded times, its metadata carries the
        shutter cards, ``MJD-BEG`` and ``EXPTIME REQUESTED``.
    config : `ShutterTimingConfig`, optional
        Configuration for computing the times; defaults to
        ``ShutterTimingConfig()``.

    Returns
    -------
    timing : `ShutterTiming` or `None`
        The times, or `None` if there are none recorded and no shutter cards
        (exposures before 2025-10-29), or as for `computeShutterTiming`.

    Raises
    ------
    RuntimeError
        Raised if the shutter cards are present but one of them is missing.

    Notes
    -----
    Without recorded times, the requested exposure time is read from the
    ``EXPTIME REQUESTED`` card, which the obs_lsst raw formatter added from
    DM-56570 to DM-56693.  Images made before that card existed fall back to
    a list of the exposures known to need re-anchoring; other timestamps in
    them are used unchecked.
    """
    metadata = exposure.metadata

    # The times recorded by ISR.
    timing = ShutterTiming.fromMetadata(metadata)
    if timing is not None:
        return timing

    # Exposures before the shutter cards were added.
    if "SHUTTER OPEN STARTTIME TAI MJD" not in metadata:
        log.info("No shutter motion cards (exposures before 2025-10-29): no shutter-corrected times.")
        return None

    # The requested exposure time: from its card, or for images made before
    # the card existed, from the list of exposures known to need it.
    requestedExposureTime = metadata.get(_REQUESTED_EXPOSURE_TIME_CARD)
    if requestedExposureTime is None:
        requestedExposureTime = _KNOWN_CLOCK_FAULT_EXPOSURES.get(exposure.visitInfo.id)

    profiles = ShutterMotionProfileFull.fromExposure(exposure)

    return computeShutterTiming(profiles.profile_open, profiles.profile_close, exposure.detector,
                                float(metadata["MJD-BEG"]), requestedExposureTime, config)


def _inRange(value: float, bounds: tuple[float, float]) -> bool:
    """Return whether lo <= value <= hi (False for NaN)."""
    return bounds[0] <= value <= bounds[1]


def _fitInRange(motion: ShutterMotionProfile) -> bool:
    """Return whether the Hall fit of a blade motion is within the nominal
    ranges.
    """
    modelStartTime, pivot1, pivot2, *jerks = motion._fitParameters("hallSensorFit")
    displacementAt0p9s = motion.displacement(np.array([modelStartTime + 0.9]))[0]

    return (
        0.0 <= pivot1 <= pivot2
        and _inRange(pivot1, _PIVOT1_RANGE)
        and _inRange(pivot2, _PIVOT2_RANGE)
        and all(_inRange(abs(j), _ABS_JERK_RANGE) for j in jerks)
        and _inRange(abs(modelStartTime), _ABS_MODEL_START_TIME_RANGE)
        and _inRange(displacementAt0p9s, _DISPLACEMENT_AT_0P9S_RANGE)
    )


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


def _edgeTimes(
    motions: tuple[ShutterMotionProfile, ShutterMotionProfile],
    closeRel: float,
    config: ShutterTimingConfig,
    focalPlaneY: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute flux-weighted opening and closing times along the blade travel.

    Parameters
    ----------
    motions : `tuple` [`lsst.ip.isr.ShutterMotionProfile`]
        The opening and closing motions.
    closeRel : `float`
        Start time of the closing motion, s after the opening one.
    config : `ShutterTimingConfig`
        The beam model.
    focalPlaneY : `numpy.ndarray`
        Focal-plane y (mm) of the points.

    Returns
    -------
    tOpen, tClose : `numpy.ndarray`
        Flux-weighted times at which the shutter opens and closes at the
        points, s after the start of the opening motion; NaN where a blade
        does not reach the beam.
    """
    # Blade positions of each point's beam quadrature nodes (focal-plane y is
    # CCS x, the blade coordinate).
    nodes, w = _beamQuadrature(np.array(config.beamOffsets))
    s = config.beamScale*focalPlaneY[:, None] + nodes

    # Crossing times of each blade edge at the quadrature nodes, s after the
    # start of the opening motion.
    tO, tauC = (motion.timeAtPosition(_ENCODER_CENTER - s) for motion in motions)
    tC = tauC + closeRel

    # Flux-weighted times of the two edges.
    return tO @ w, tC @ w
