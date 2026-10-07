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

The camera control system fits each shutter blade motion, from its
Hall-sensor positions, with this model and writes the fit to the exposure
metadata.  Conventions (CTN-002, https://ctn-002.lsst.io):

- piecewise-constant jerk ``j0`` on ``[0, t1)``, ``j1`` on ``[t1, t2)``,
  ``j2`` on ``[t2, inf)``, from rest at model time 0;
- model time is ``(T - startTime) - ModelStartTime`` (s); pivots are in model
  time;
- the model gives the displacement magnitude ``d(t) >= 0`` along the travel
  direction, and the encoder position is
  ``startPosition + travelSign * d(t)``.

The inverse is restricted to the monotonic branch (up to the first velocity
zero or velocity minimum) and is NaN outside it.
"""

from __future__ import annotations

__all__ = ["ThreeJerksParams", "ThreeJerksTrajectory", "checkFitParams"]

import dataclasses

import numpy as np


@dataclasses.dataclass(frozen=True)
class ThreeJerksParams:
    """One ThreeJerksModelv1 fit (s, s, s, mm/s^3, mm/s^3, mm/s^3)."""

    modelStartTime: float
    pivot1: float
    pivot2: float
    jerk0: float
    jerk1: float
    jerk2: float

    def values(self):
        return (self.modelStartTime, self.pivot1, self.pivot2, self.jerk0, self.jerk1, self.jerk2)


def _cubicRoots(a, b, c, d):
    """All three complex roots of ``a x^3 + b x^2 + c x + d = 0``
    (vectorised Cardano; shape ``(3,) + broadcast shape``).
    """
    a, b, c, d = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in (a, b, c, d)))
    with np.errstate(all="ignore"):
        d0 = b * b - 3.0 * a * c
        d1 = 2.0 * b**3 - 9.0 * a * b * c + 27.0 * a * a * d
        disc = np.sqrt((d1 * d1 - 4.0 * d0**3).astype(complex))
        cp = (d1 + disc) / 2.0
        cm = (d1 - disc) / 2.0
        cc = np.where(np.abs(cp) >= np.abs(cm), cp, cm)
        bigC = cc ** (1.0 / 3.0)
        xi = np.exp(2j * np.pi / 3.0)
        roots = []
        for k in range(3):
            ck = bigC * xi**k
            roots.append(-(b + ck + np.where(ck != 0, d0 / ck, 0.0)) / (3.0 * a))
    return np.stack(roots)


class ThreeJerksTrajectory:
    """ThreeJerksModelv1 trajectory of one blade motion.

    Parameters
    ----------
    params : `ThreeJerksParams`
        The fit; all values finite and ``0 <= pivot1 <= pivot2``
        (`ValueError` otherwise).
    startPosition : `float`
        Encoder position (mm) at model time zero.
    travelSign : `int`
        +1 if the encoder position increases during the move, else -1.
    """

    def __init__(self, params, startPosition, travelSign):
        if travelSign not in (-1, 1):
            raise ValueError(f"travelSign must be +1 or -1, got {travelSign!r}")
        if not np.all(np.isfinite(params.values())):
            raise ValueError(f"non-finite ThreeJerks parameters: {params}")
        t1, t2 = float(params.pivot1), float(params.pivot2)
        if not 0.0 <= t1 <= t2:
            raise ValueError(f"invalid pivots (need 0 <= pivot1 <= pivot2): {t1}, {t2}")
        self.params = params
        self.startPosition = float(startPosition)
        self.travelSign = int(travelSign)
        self._tb = np.array([0.0, t1, t2])
        self._j = np.array([params.jerk0, params.jerk1, params.jerk2], dtype=float)
        state = np.zeros((3, 3))
        s = v = a = 0.0
        for k in range(3):
            state[k] = (s, v, a)
            if k < 2:
                h, j = self._tb[k + 1] - self._tb[k], self._j[k]
                s, v, a = (
                    s + v * h + a * h * h / 2 + j * h**3 / 6,
                    v + a * h + j * h * h / 2,
                    a + j * h,
                )
        self._state = state
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
            cands = []
            if j != 0:
                disc = a0 * a0 - 2.0 * j * v0
                if disc >= 0:
                    sq = np.sqrt(disc)
                    cands += [(-a0 - sq) / j, (-a0 + sq) / j]
            elif a0 != 0:
                cands.append(-v0 / a0)
            if a0 < 0 < j:
                cands.append(-a0 / j)
            cands = [h for h in cands if 0 < h <= hMax]
            if cands:
                return float(self._tb[k] + min(cands))
        return np.inf

    def sva(self, dt):
        """Displacement (mm), velocity (mm/s), acceleration (mm/s^2) at model
        time ``dt`` (zeros for ``dt < 0``).
        """
        t = np.asarray(dt, dtype=float)
        tc = np.maximum(t, 0.0)
        k = np.clip(np.searchsorted(self._tb, tc, side="right") - 1, 0, 2)
        s0, v0, a0 = (self._state[k, i] for i in range(3))
        j = self._j[k]
        h = tc - self._tb[k]
        s = s0 + h * (v0 + h * (a0 / 2 + h * j / 6))
        v = v0 + h * (a0 + h * j / 2)
        a = a0 + h * j
        return s, v, a

    def displacement(self, dt):
        """Displacement magnitude (mm) at model time ``dt`` (s)."""
        return self.sva(dt)[0]

    def timeOfDisplacement(self, d):
        """Model time (s) at which the displacement is ``d`` (mm); NaN
        outside ``[0, dMax]``.

        Segment 0 is solved exactly; segments 1-2 by Cardano as a first guess
        finished by a safeguarded Newton-bisection on the segment's bracket.
        """
        d = np.asarray(d, dtype=float)
        out = np.full(d.shape, np.nan)
        valid = np.isfinite(d) & (d >= 0.0) & (d <= self.dMax)
        if self.tStop <= 0.0 or not np.any(valid):
            return out
        tb = np.minimum(np.append(self._tb, np.inf), self.tStop)
        sb = np.append(self.sva(tb[:3])[0], self.dMax)
        dv = d[valid]
        k = np.clip(np.searchsorted(sb, dv, side="right") - 1, 0, 2)
        for _ in range(2):
            k = np.where((tb[k + 1] <= tb[k]) & (k > 0), k - 1, k)
        s0, v0, a0 = (self._state[k, i] for i in range(3))
        j = self._j[k]
        lo = np.zeros_like(dv)
        hi = tb[k + 1] - tb[k]

        def seg(h):
            return s0 + h * (v0 + h * (a0 / 2 + h * j / 6))

        unb = ~np.isfinite(hi)
        if np.any(unb):
            hb = np.where(unb, 1.0, 0.0)
            for _ in range(64):
                short = unb & (seg(hb) < dv)
                if not np.any(short):
                    break
                hb = np.where(short, 2.0 * hb, hb)
            hi = np.where(unb, hb, hi)

        h = 0.5 * (lo + hi)
        seg0 = k == 0
        h[seg0] = np.cbrt(6.0 * dv[seg0] / self._j[0])
        other = ~seg0
        if np.any(other):
            roots = _cubicRoots(j[other] / 6, a0[other] / 2, v0[other], s0[other] - dv[other])
            hm = hi[other]
            re = roots.real
            tol = 1e-9 + 1e-6 * np.minimum(hm, 1.0)
            ok = np.isfinite(re) & (np.abs(roots.imag) <= 1e-5) & (re >= -tol) & (re <= hm + tol)
            cand = np.where(ok, re, np.inf).min(axis=0)
            h[other] = np.where(np.isfinite(cand), np.clip(cand, 0.0, hm), h[other])

        rtol = 1e-13 * (1.0 + np.abs(dv))
        idx = np.flatnonzero(~seg0)
        s0a, v0a, a0a, ja, da, loA, hiA, hA, rtA = (x[idx] for x in (s0, v0, a0, j, dv, lo, hi, h, rtol))
        for _ in range(200):
            if idx.size == 0:
                break
            s = s0a + hA * (v0a + hA * (a0a / 2 + hA * ja / 6))
            v = v0a + hA * (a0a + hA * ja / 2)
            r = s - da
            loA = np.where(r <= 0, hA, loA)
            hiA = np.where(r >= 0, hA, hiA)
            with np.errstate(divide="ignore", invalid="ignore"):
                hn = hA - r / v
            inside = np.isfinite(hn) & (hn > loA) & (hn < hiA)
            hn = np.where(inside, hn, 0.5 * (loA + hiA))
            smallStep = inside & (np.abs(hn - hA) <= 1e-14)
            done = (np.abs(r) <= rtA) | (hiA - loA <= 1e-15 * (1.0 + hiA)) | smallStep
            hA = np.where(smallStep, hn, hA)
            h[idx[done]] = hA[done]
            keep = ~done
            hA = np.where(keep, hn, hA)
            idx, s0a, v0a, a0a, ja, da, loA, hiA, hA, rtA = (
                x[keep] for x in (idx, s0a, v0a, a0a, ja, da, loA, hiA, hA, rtA)
            )
        h[idx] = hA
        out[valid] = tb[k] + h
        return out


def _inRange(value, bounds):
    lo, hi = bounds
    return bool(np.isfinite(value) and lo <= value <= hi)


def checkFitParams(params, config):
    """True if a fit lies within the nominal ranges of ``config``
    (`ShutterTimingConfig`): pivots, every |jerk|, |ModelStartTime| and the
    model displacement at model time 0.9 s.  None or invalid fits fail.
    """
    if params is None:
        return False
    ok = (
        _inRange(params.pivot1, config.pivot1Range)
        and _inRange(params.pivot2, config.pivot2Range)
        and all(_inRange(abs(j), config.absJerkRange) for j in (params.jerk0, params.jerk1, params.jerk2))
        and _inRange(abs(params.modelStartTime), config.absModelStartTimeRange)
    )
    if not ok:
        return False
    try:
        traj = ThreeJerksTrajectory(params, 0.0, -1)
    except ValueError:
        return False
    return _inRange(float(traj.displacement(0.9)), config.displacementAt0p9sRange)
