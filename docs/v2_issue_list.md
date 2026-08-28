# v2 issue list

Every entry was verified by reading the code on `v2`, not inferred from file
names or comments. Items move to **Resolved** only when a fix plus a regression
test is committed.

---

## Resolved

### P0 — result-changing defects

| # | Location | Defect | Fixed in |
|---|---|---|---|
| P0-1 | `compute_nor_metrics()` | The novel object was always scored with a 9x7 cm box and the familiar object with a 4 cm radius. The regions have different areas (63 vs 50.3 cm²), so novel contact time was inflated by construction and the discrimination index was biased independently of behaviour. | `8ea9013` |
| P0-2 | `DLCA_SocInt` | Avoidance and following thresholded each animal's *scalar* speed, which carries no direction. Both avoidance flags fired whenever both animals moved, and parallel locomotion satisfied every rule. Approach/retreat were symmetric and could not attribute at all. | `0f46884` |
| P0-3 | core `events.R`, all assays | `normalize_event_vector()` mapped `NA` to `FALSE`, so untracked frames were scored as confident absence of behaviour, and denominators counted every frame. | `f491dad` |
| P0-4 | `CleanTrackingData()`, NOR | Interpolation was unbounded: `na_interpolation()` filled gaps of any length, overwrote the `likelihood` column, and `zoo::na.locf` forward-filled leading and trailing runs. | `ec0f4ef` |
| P0-5 | `compute_socp_metrics()` | No unit, fps or TrackingData validation. Centimetre thresholds could be applied silently to pixel coordinates. | `7e1ee18` |
| P0-6 | `ZoneReport()`, `EPMAnalysis()` | Zone `transitions` counted onsets *and* offsets. Unobserved frames were credited to inverted zones (OFT periphery, EPM arms). Nose dips were counted as `transitions / 2`, giving half-integers. | `927d915` |
| P0-7 | `integratevector()` | Returned `0` for the first frame, averaging n frames over n-1 real intervals and biasing every mean speed by (n-1)/n — 3.3% over 30-frame bins. | `927d915` |
| P0-8 | `MultiFileReport()` | Called `rbindlist()` with `data.table` declared nowhere, so **every EPM run failed** unless the package happened to be attached. Invisible until the assay was actually executed. | `519b9bd` |
| P0-9 | `DLCA_OFT_advanced_metrics` | Model chosen by `requireNamespace("lmerTest")`. With one session per animal, `lme4` refuses `(1|ID)` outright and the error was swallowed, so installing lmerTest changed the output from p-values to nothing at all. | `8d501b1` |
| P0-10 | autoencoder motifs | Train/validation split was random over *windows*, so the same animal appeared on both sides and `val_loss` measured memorisation rather than generalisation. | `ab46426` |

### P1 — tracking, geometry, QC and tests

| # | Defect | Fixed in |
|---|---|---|
| P1-1 | Event engine had no validity mask, seconds-valued minimum bout, merge gap, valid analysed time, median bout or percent-valid-time. | `f491dad` |
| P1-2 | Calibration was scalar x/y scaling only; no homography, canonical frame or reprojection error. | `789fdb3` |
| P1-3 | QC was report-only and partial: no gap length, interpolation fraction, implausible displacement, skeleton check, identity swaps or arena violations. | `47e6624` |
| P1-4 | `ReadDLCDataFromCSV()` accepted non-monotonic and non-contiguous frame numbering, making non-adjacent frames look adjacent. | `789fdb3` |
| P1-5 | `point.in.polygon(...) == 1` excluded boundary points, so a frame on a shared zone edge belonged to no zone. | `ec0f4ef` |
| P1-6 | No CI of any kind. | `ec0f4ef` |
| P1-7 | EPM was not executable under test. | `927d915` |
| P1-8 | `AddZones()` accepted self-intersecting or zero-area zone polygons, which silently under-count occupancy. | `519b9bd` |

### P3 — reproducibility and architecture

| # | Defect | Fixed in |
|---|---|---|
| P3-2 | No configuration files; arena size, fps and thresholds were literals inside loops. | `2992237` |
| P3-3 | No run manifest, commit SHA or dependency versions in outputs. | `2992237` |
| P3-4 | `DLCA_SocP` re-implemented the metadata helpers locally. | `7e1ee18` |
| P3-5 | `oft_center_exploration_score` was cohort-relative but named as a subject-level measurement; habituation p-values were within-animal and nominal but named as ordinary tests. | `f35a400` |

---

## Outstanding

| # | Item | Why it is not done |
|---|---|---|
| **O-1** | Assay scripts still use `CalibrateTrackingData()` scalar scaling rather than the tested `arena_calibration()` homography. | Migrating changes every calibrated coordinate and therefore every distance, speed and zone boundary. It needs a validation run against existing results on real data, which this repository does not contain. |
| **O-2** | OFT and SocInt accept a configuration *overlay* but are not on the full `load_assay_config()` schema, and still carry machine-specific defaults (8 remaining absolute paths, all overridable). | Both are large evolved scripts; a full schema migration risks breaking working analyses without real data to verify against. |
| **O-3** | No native SLEAP HDF5 import. Conversion still goes through the `01_SLEAPcoords/` notebooks to DLC-style CSV. | Requires fixture data in the actual SLEAP analysis HDF5 layout to verify axis order, track identity and node names. Inventing the schema would be worse than not implementing it. |
| **O-4** | SocInt thresholds remain in pixels. | Converting them requires arena calibration *and* re-validation of every distance against manually scored video. `validate_socint_units()` blocks centimetre mode until that is done deliberately. |
| **O-5** | Motif analyses have stability tooling (`adjusted_rand_index()`, `cluster_seed_stability()`, `group_kfold()`) but no stability *report* is generated by the motif scripts. | Needs real data to be meaningful; the tooling is tested and ready to wire in. |
| **O-6** | `LabelReport()` still uses the legacy `CalculateTransitions(...) / 2`. | Only affects the supervised-classification path, which no production assay script calls. |
| **O-7** | No dependency lock (`renv`). | CI pins the package set explicitly; a lockfile would need generating against a verified working environment. |

---

## Validation debt

Nothing in this repository has been validated against manually scored video.
The full list is at the end of [assay_definitions.md](assay_definitions.md).
The highest-priority items are the NOR contact distance and orientation angle,
the SocP contact distance, the SocInt distance thresholds, and the inherited
novel/familiar metadata mapping.
