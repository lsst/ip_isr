# Shutter-timing test data

Expected results for `lsst.ip.isr.shutterTiming`, computed by an independent
reference implementation of the same algorithm (`shutter_timing`, header-card
path, with the linear beam model of `ShutterTimingConfig`'s defaults; its
version is in the `generator` field of each JSON file), for seven
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
    results (`center_mjd_tai`, `qc_flags`).
  - `samples`: `[detector, x, y, t_mid_mjd_tai]` at 8 positions per detector
    (centre, 4 pixel-edge corners, one interior point, 50 px off the detector,
    and 150 px off), computed exactly at each position.
- `MC_O_20260712_000100.json` also has `reanchored`: the expected results
  (same format) when its close start time is shifted by +10 ms, so that the
  shutter clocks disagree and the start times are re-anchored to `MJD-BEG`
  (a -5 ms shift gives identical results).

Agreement target: centre and per-source times to <= 10 us; identical flags.

The JSON files are written compactly (one line each).
