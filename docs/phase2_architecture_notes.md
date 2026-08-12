# Phase 2 shared-core design note

## Current duplication

The Phase 1 code has several correct but repeated implementations:

- assay scripts and `Behavioral_Metrics_Phase1.R` read `Tracking$data`,
  `Tracking$frames`, and `Tracking$fps` directly;
- distance, vector-angle, polygon, and arena checks are split between the
  upstream-derived DLCAnalyzer file, the Phase 1 helper, and SocInt;
- event normalization, minimum-duration filtering, entry counting, latency,
  bout duration, and interbout intervals are implemented independently in the
  Phase 1 helper, OFT, and SocInt;
- frame/second and distance/speed conversions are embedded in calculations;
- likelihood, missing-coordinate, geometry, and calibration checks are local
  to individual workflows.

This makes it easy for fixes to diverge while preserving identical biological
definitions on paper.

## Proposed shared modules

Add dependency-light files under `02_SLEAPanalzyer/core/`:

- `tracking_data.R`: validate and access the existing TrackingData list shape;
- `io.R`: deterministic core loading and thin import wrappers only;
- `validation.R`: reusable scalar, fps, landmark, and vector checks;
- `units.R`: explicit coordinate/threshold units and time/speed conversion;
- `geometry.R`: distances, angles, boxes, polygons, and arena diagnostics;
- `calibration.R`: calibration validation and scale calculations, without I/O;
- `zones.R`: point-in-polygon and zone-membership helpers;
- `events.R`: event filtering and canonical bout tables;
- `metrics.R`: generic duration/frequency/latency summaries built from bouts;
- `qc.R`: report-only missingness, likelihood, coordinate, calibration, and
  geometry summaries.

These are sourced as plain R files. Phase 2 does not introduce S3 classes or
convert the repository into a package.

## Functions to extract

- TrackingData: `validate_tracking_data()`, `get_tracking_frames()`,
  `get_tracking_fps()`, `get_tracking_duration()`,
  `get_point_coordinates()`, and `has_landmarks()`.
- Units: `validate_coordinate_unit()`, `validate_threshold_unit()`,
  `frames_to_seconds()`, and `distance_to_speed()`.
- Geometry: point/vector distance, vector angle, tracking-point distance and
  target angle, polygon validation/area/membership, axis-aligned boxes, and
  rectangular arena diagnostics.
- Events: logical normalization, minimum-duration filtering, canonical bout
  table generation, event duration, latency, entry count, and interbout gaps.
- QC: coordinate missingness/validity and likelihood summaries plus reporting
  wrappers for calibration and geometry validation.

`Behavioral_Metrics_Phase1.R` remains as a compatibility and assay-definition
layer. Its existing public helper names will delegate to the shared core so
Phase 1 callers continue to work.

## Intentionally assay-specific

The following remain outside the generic core:

- NOR contact geometry asymmetry, orientation windows, proximity thresholds,
  body-exclusion rule, and legacy `NovelLoc` remapping;
- SocP contact/proximity thresholds and novel/familiar remapping;
- EPM zone definitions and 60 cm calibration assumption;
- SocInt contact, following, avoidance, approach, and retreat definitions;
- OFT biological zones, immobility/movement thresholds, plotting, batch paths,
  metadata lookup, and all export schemas.

NOR is the only assay script migrated in this phase. Other scripts and the
upstream-derived DLCAnalyzer functions remain in place for compatibility.
