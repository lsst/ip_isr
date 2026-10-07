#!/usr/bin/env python3
"""Expected `shutter_timing` results for mutated fixture headers.

Run with the venv of github.com/mjuric/shutter-timing (not the LSST stack),
from any directory:

    <shutter-timing>/.venv/bin/python make_mutations.py [SHUTTER_TIMING_REPO]

For each mutation of MC_O_20260712_000100.json's metadata, read the header
through the reference header-card path (``io.read_header``), compute the
correction-table rows (``table.compute_rows``) for a few detectors and the
per-source samples (same rules as ``corrections.corrected_midpoints``), and
write ``mutations.json`` next to this script.  A header the reference cannot
read is recorded with ``status`` 2 (UNAVAILABLE) and the reason.
"""

import importlib.util
import json
import pathlib
import sys

from astropy.io import fits

HERE = pathlib.Path(__file__).resolve().parent
BASE = "MC_O_20260712_000100"
DETECTORS = [0, 4, 30, 94, 120, 168, 188]
MJD_PER_MS = 1e-3 / 86400.0
COEFFS = ("c_u_s", "c_uu_s", "c_v_s", "c_uv_s", "c_vv_s")

#: name -> (card, new value or callable(old) -> new; None deletes the card);
#: card None: the header is unchanged (see A1_MM).
MUTATIONS = {
    "close_start_plus_10ms": ("SHUTTER CLOSE STARTTIME TAI MJD", lambda v: v + 10 * MJD_PER_MS),
    "close_start_minus_5ms": ("SHUTTER CLOSE STARTTIME TAI MJD", lambda v: v - 5 * MJD_PER_MS),
    "open_pivot1_out_of_range": ("SHUTTER OPEN HALLSENSORFIT PIVOTPOINT1", 0.30),
    "close_jerk2_out_of_range": ("SHUTTER CLOSE HALLSENSORFIT JERK2", 50000.0),
    "close_jerk1_string": ("SHUTTER CLOSE HALLSENSORFIT JERK1", "abc"),
    "open_jerk0_missing": ("SHUTTER OPEN HALLSENSORFIT JERK0", None),
    "open_start_string": ("SHUTTER OPEN STARTTIME TAI MJD", "not-a-number"),
    "open_beg_minus_20ms": ("MJD-BEG", lambda v: v - 20 * MJD_PER_MS),
    "open_jerk0_negative": ("SHUTTER OPEN HALLSENSORFIT JERK0", lambda v: -v),
    "exptime_zero": ("EXPTIME", 0.0),
    "a1_decreasing": (None, None),
    "a1_increasing": (None, None),
}

#: name -> reference a1_mm (travel sign -> mm) for the A1 cases; the stack's
#: a1Decreasing / a1Increasing are the -1 / +1 values.
A1_MM = {
    "a1_decreasing": {-1: 0.4, 1: -0.3},
    "a1_increasing": {-1: 0.4, 1: -0.3},
}

#: Mutations of another fixture than BASE.  MC_O_20260105_000250 opens with
#: MINUSX: both motions move 0 -> 750 (travel sign +1), whereas BASE opens
#: with PLUSX (both -1).
OTHER_BASE = {"a1_increasing": "MC_O_20260105_000250"}


def _header(md):
    hdr = fits.Header()
    for k, v in md.items():
        hdr[("HIERARCH " + k) if (len(k) > 8 or " " in k) else k] = v
    return hdr


def main(argv):
    repo = pathlib.Path(argv[1]) if len(argv) > 1 else pathlib.Path(
        "/sdf/data/rubin/user/mjuric/projects/github.com/mjuric/shutter-timing")
    sys.path.insert(0, str(repo))
    from shutter_timing import io, table

    spec = importlib.util.spec_from_file_location("make_fixtures", repo / "stackport" / "make_fixtures.py")
    mf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mf)

    bases = {b: json.loads((HERE.parent / f"{b}.json").read_text())
             for b in {BASE, *OTHER_BASE.values()}}
    base = bases[BASE]
    ctx = table.TableContext.load(raw_root=None, detectors=DETECTORS)
    out = dict(base=BASE, visit=base["visit"], detectors=DETECTORS, generator=mf._describe(), mutations={})
    for name, (card, change) in MUTATIONS.items():
        mbase = OTHER_BASE.get(name, BASE)
        md = dict(bases[mbase]["metadata"])
        visit = bases[mbase]["visit"]
        if card is None:
            pass
        elif change is None:
            del md[card]
        else:
            md[card] = change(md[card]) if callable(change) else change
        a1 = A1_MM.get(name)
        doc = dict(base=mbase, visit=visit, card=card, metadata=md)
        if a1 is not None:
            doc["a1_mm"] = {str(k): v for k, v in a1.items()}
        ctx.a1_mm = a1
        try:
            exp = io.read_header(_header(md), obs_id=BASE)
            rows, info = table.compute_rows(exp, visit, ctx)
        except Exception as e:  # noqa: BLE001 -- recorded as UNAVAILABLE
            doc.update(status=2, reason=f"{type(e).__name__}: {e}")
        else:
            doc.update(
                status=None, policy=info["policy"], exposure_qc_flags=int(info["qc_flags"]),
                visit_mjd_tai=float(info["t_mid_visit_mjd_tai"]),
                detectors={
                    str(int(r["detector"])): dict(
                        axis=str(r["axis"]), center_mjd_tai=float(r["t_mid_center_mjd_tai"]),
                        coefficients_s=[float(r[c]) for c in COEFFS],
                        max_abs_residual_s=float(r["max_abs_residual_s"]),
                        effective_exposure_time_s=float(r["t_eff_s"]), qc_flags=int(r["qc_flags"]),
                    ) for r in rows
                },
                samples=mf.per_source(rows, ctx, visit),
            )
        out["mutations"][name] = doc
        print(name, doc.get("policy"), doc.get("exposure_qc_flags"), doc.get("reason", ""))
    (HERE / "mutations.json").write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
