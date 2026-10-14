# Shutter-timing test data

Expected results for `lsst.ip.isr.shutterTiming`, computed by an independent
reference implementation of the same algorithm (`shutter_timing`, header-card
path; its version is in the `generator` field of each JSON file), for seven
science detectors from the centre to the corners of the focal plane (0, 4, 30,
94, 120, 168, 188). Only header cards are included, no pixel data. Read-only.

## Files

- `MC_O_<obs>.json`, one per exposure (`MC_O_20260105_000250`, open blade
  MINUSX, 10 s; `MC_O_20260712_000100`, open blade PLUSX, 30 s):
  - `metadata`: the primary-header cards the code reads (all `SHUTTER*` cards,
    `MJD-BEG`, `MJD-END`, `EXPTIME`, `SHUTTIME`, ...), keys without
    `HIERARCH`.
  - `detectors`: per detector, its `geometry` (`center_mm`, the focal-plane
    position of its centre; `center_pix`, the centre's pixel coordinates;
    `size` in pixels; 10 um pixels and yaw 0, as for every LSSTCam detector),
    from which the tests build it with `DetectorWrapper`, and the expected
    results (`axis`, `center_mjd_tai`, `coefficients_s` = c_u, c_uu, c_v,
    c_uv, c_vv, `max_abs_residual_s`, `effective_exposure_time_s`,
    `qc_flags`).
  - `visit_mjd_tai`: the focal-plane-centre time; `header_mid_mjd_tai`;
    `exposure_qc_flags`.
  - `samples`: `[detector, x, y, t_mid_mjd_tai or null, status]` at 8
    positions per detector (centre, 4 pixel-edge corners, one interior point,
    50 px off the detector, and 150 px off); status as the reference
    implementation reports it, 0 OK, 1 DEGRADED, 2 UNAVAILABLE. The stack code
    has no DEGRADED status: a per-source time exists unless the status is 2,
    and the tests only distinguish 2 from not 2.
- `mutations.json`: single-card changes to the metadata of `MC_O_20260712_000100`
  (`base`), each with the cards it sets (`set`) and the reference results for
  the same detectors in the same format as above. The mutations: the close
  start time shifted by +10 ms and -5 ms (clock re-anchoring); PIVOTPOINT1 and
  JERK2 out of range; and a negative JERK0.
- `beam_at_L3S1_z9.618_rot0_evaluated.tnt`: the LSSTCam shutter-plane beam
  table (an LSST Camera raytrace, LCA-20578), identical to the copy shipped
  by obs_lsst (see its `resources/shutter/README.md` for credit). A byte-identical copy is public in
  lsst-dm/ap_pipe-notebooks (branch `tickets/DM-50985`, commit 3e6fc39,
  `data/beam_at_L3S1_-z9.618_rot0_evaluated.tnt`).

Agreement target: centre times and quadratic terms at a 2000-pixel lever arm
to <= 10 us; identical `axis`, flags and detector statuses; a per-source time
exactly where the expected status is not 2.

The JSON files are written compactly (one line each).
