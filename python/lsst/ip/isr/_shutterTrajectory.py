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
"""ThreeJerksModelv1 blade trajectory (private helper of `shutterTiming`).

The camera fits each blade motion, from its Hall-sensor positions, with this
model (CTN-002, https://ctn-002.lsst.io): piecewise-constant jerk ``j0`` on
``[0, t1)``, ``j1`` on ``[t1, t2)``, ``j2`` on ``[t2, inf)``, from rest at
model time 0 (``(T - startTime) - ModelStartTime``, s).  It gives the
displacement ``d(t) >= 0`` along the travel direction; the encoder position
is ``startPosition + travelSign * d(t)``.
"""

from __future__ import annotations

__all__ = ["ThreeJerksParams", "ThreeJerksTrajectory"]

import dataclasses

import numpy as np


@dataclasses.dataclass(frozen=True)
class ThreeJerksParams:
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


class ThreeJerksTrajectory:
    """ThreeJerksModelv1 trajectory of one blade motion.

    Parameters
    ----------
    params : `ThreeJerksParams`
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

    def __init__(self, params: ThreeJerksParams, startPosition: float, travelSign: int):
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

    def _branchEnd(self):
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
