# Shutter-timing test data

Expected results for `lsst.ip.isr.shutterTiming`, computed by an independent
reference implementation of the same algorithm (`shutter_timing`, header-card
path; its version is in the `generator` field of each JSON file). Only header
cards are included, no pixel data. Read-only.

## Files

- `MC_O_<obs>.json`, one per exposure:
  - `metadata`: the primary-header cards the code reads (all `SHUTTER*` cards,
    `MJD-BEG`, `MJD-END`, `EXPTIME`, `SHUTTIME`, ...), keys without
    `HIERARCH`.
  - `detectors`: the expected per-detector results for the 189 science
    detectors (`axis`, `center_mjd_tai`, `coefficients_s` = c_u, c_uu, c_v,
    c_uv, c_vv, `max_abs_residual_s`, `effective_exposure_time_s`,
    `qc_flags`).
  - `visit_mjd_tai`: the focal-plane-centre time; `header_mid_mjd_tai`;
    `policy`; `exposure_qc_flags`.
  - `samples`: `[detector, x, y, t_mid_mjd_tai or null, status]` at 8
    positions per detector (centre, 4 pixel-edge corners, one interior point,
    50 px off the detector, and 150 px off); status 0 OK, 1 DEGRADED,
    2 UNAVAILABLE.
  - `constants`: the parameter ranges, clock thresholds and other constants
    the values were computed with (the `ShutterTimingConfig` defaults).
  - `MC_O_20250810_000030.json` predates the shutter cards: the expected
    result is UNAVAILABLE (NO_PROFILE).
- `mutations.json`: the metadata of `MC_O_20260712_000100` (for
  `a1_increasing`, of `MC_O_20260105_000250`) with one card changed per entry
  (`card`; the full mutated `metadata` is included), and the reference results
  for detectors 0, 4, 30, 94, 120, 168 and 188 in the same format as above.
  Entries the reference rejects have `status` 2 (UNAVAILABLE) and a `reason`.
  The mutations: the close start time shifted by +10 ms and -5 ms (clock
  re-anchoring); PIVOTPOINT1 and JERK2 out of range; a non-numeric JERK1; a
  missing JERK0; a non-numeric open start time; MJD-BEG 20 ms early; a
  negative JERK0; EXPTIME 0 (overlapping blades); and two unchanged headers
  computed with a nonzero blade-edge offset A1 (`a1_mm`: travel sign ->
  offset in mm), one per travel direction.
- `beam_at_L3S1_z9.618_rot0_evaluated.tnt`: the LSSTCam shutter-plane beam
  table (raytrace by A. Rasmussen, LSST Camera, LCA-20578), identical to the
  copy shipped by obs_lsst. A byte-identical copy is public in
  lsst-dm/ap_pipe-notebooks (branch `tickets/DM-50985`, commit 3e6fc39,
  `data/beam_at_L3S1_-z9.618_rot0_evaluated.tnt`).
- `lsstcam_detector_geometry.csv`: the LSSTCam detector geometry exported from
  obs_lsst w_2026_10 (`PIXELS -> FOCAL_PLANE`, affine per detector, and the
  detector type), so that the ip_isr tests need not import obs_lsst.

Agreement target: centre times and quadratic terms at a 2000-pixel lever arm
to <= 10 us; identical `axis`, flags and statuses.
