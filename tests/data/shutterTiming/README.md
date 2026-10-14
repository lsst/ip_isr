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
    results (`center_mjd_tai`, `effective_exposure_time_s`, `qc_flags`).
  - `visit_mjd_tai`: the focal-plane-centre time; `header_mid_mjd_tai`;
    `exposure_qc_flags`.
  - `samples`: `[detector, x, y, t_mid_mjd_tai]` at 8 positions per detector
    (centre, 4 pixel-edge corners, one interior point, 50 px off the detector,
    and 150 px off), computed exactly at each position.
- `mutations.json`: single-card changes to the metadata of `MC_O_20260712_000100`
  (`base`), each with the cards it sets (`set`) and the reference results for
  the same detectors in the same format as above. The mutations: the close
  start time shifted by +10 ms and -5 ms (clock re-anchoring); PIVOTPOINT1 and
  JERK2 out of range; and a negative JERK0.

Agreement target: centre, focal-plane-centre and per-source times and the
effective exposure time to <= 10 us; identical flags and detector statuses.

The JSON files are written compactly (one line each).
