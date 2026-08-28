# v2 prioritized issue list

Every entry below was verified by reading the code on `v2`, not inferred from
file names or comments. Items are removed only when a fix plus a regression
test is committed.

## P0 - result-changing defects

| # | Location | Defect |
|---|----------|--------|
| P0-1 FIXED | `Behavioral_Metrics_Phase1.R` `compute_nor_metrics()` | The novel object is always scored with a 9x7 cm axis-aligned box while the familiar object is always scored with a 4 cm radius. The two detectors cover different areas (63 vs 50.3 cm^2), so novel contact time is inflated by construction and the discrimination index is biased. |
| P0-2 FIXED | `DLCA_SocInt v.0.0.2.r` | `a1_avoidance_from_a2` / `a2_avoidance_from_a1` share an identical distance condition and differ only in which animal's *scalar* speed exceeds a cutoff. Scalar speed carries no direction, so both flags fire whenever both animals move. `approach_event` / `retreat_event` are fully symmetric and cannot attribute at all. |
| P0-3 FIXED | core `events.R`, all assays | `normalize_event_vector()` maps `NA` to `FALSE`, so untracked frames are scored as confident absence of behavior. No validity mask exists anywhere and denominators use every frame. |
| P0-4 FIXED | `CleanTrackingData()`, NOR `fill_edges_and_gaps()` | Interpolation is unbounded. `imputeTS::na_interpolation()` fills gaps of any length and also overwrites the `likelihood` column; NOR's `zoo::na.locf` fills leading and trailing runs of any length. Fabricated coordinates are indistinguishable from observed ones downstream. |
| P0-5 | `compute_socp_metrics()` | No unit, fps or TrackingData validation, unlike `compute_nor_metrics()`. Centimetre thresholds can be applied silently to pixel coordinates. |
| P0-6 | `CalculateTransitions()` | Counts onsets *and* offsets, so zone `transitions` is roughly twice the entry count and is biased by whether the animal starts or ends inside the zone. |

## P1 - tracking, geometry and QC core

| # | Location | Defect |
|---|----------|--------|
| P1-1 FIXED | core `events.R` | No validity mask, no seconds-based minimum bout, no maximum merge gap, no valid analyzed time, no median bout, no percent-valid-time. |
| P1-2 | core `calibration.R` | Only scalar x/y scaling. No projective homography, no canonical arena frame, no reprojection error. |
| P1-3 | core `qc.R` | Report-only. No observed/interpolated/invalid status, no gap limits in seconds, no implausible-displacement check, no frame-continuity check. |
| P1-4 | `ReadDLCDataFromCSV()` | Accepts non-monotonic and non-contiguous frame numbering; `integratevector()` then treats non-adjacent frames as adjacent and inflates speed. |
| P1-5 FIXED | `IsInZone()`, `ZoneReport()` | `sp::point.in.polygon(...) == 1` excludes boundary points, so a point on a shared zone edge belongs to no zone. Inconsistent with core `points_in_polygon()`, which includes the boundary. |

## P1 - tests and CI

| # | Defect |
|---|--------|
| P1-6 FIXED | No CI configuration of any kind. |
| P1-7 | No synthetic assay tests for OFT or EPM; EPM is not executable under test. |

## P3 - reproducibility and architecture

| # | Defect |
|---|--------|
| P3-1 | Machine-specific `S:/` and `C:/Users/...` paths in every production script. |
| P3-2 | No configuration files; arena size, thresholds and fps are literals inside loops. |
| P3-3 | No run manifest, commit SHA, or dependency versions in outputs. |
| P3-4 | `DLCA_SocP` re-implements `read_metadata_table()` / `metadata_lookup()` locally. |
