# SLEAPanalyzer

Behavioural analysis of SLEAP and DeepLabCut tracking data for rodent assays:
Open Field (OFT), Elevated Plus Maze (EPM), Novel Object Recognition (NOR),
Social Preference (SocP) and Social Interaction (SocInt).

Built on [DLCAnalyzer](https://github.com/ETHZ-INS/DLCAnalyzer), restructured
in v2 around a shared, tested measurement core.

> **v2 status.** The five assays run through shared import, interpolation, QC,
> geometry and event segmentation. NOR, SocP and EPM are fully configuration
> driven and emit provenance manifests. OFT and SocInt accept a configuration
> overlay but have not yet been migrated onto the full schema. See
> [docs/v2_issue_list.md](docs/v2_issue_list.md) for what remains.

---

## What v2 changed, and why it matters for your results

If you have results from a previous version, these changes alter numbers.

| Change | Effect on outputs |
|---|---|
| NOR scored the novel object with a 9x7 cm box and the familiar object with a 4 cm radius. Those regions have different areas, so novel contact was inflated by construction. | Discrimination indices shrink toward zero. **Re-run any NOR analysis.** |
| SocInt attributed avoidance by scalar speed, which carries no direction, so both animals were flagged whenever both moved. | Avoidance durations fall substantially. |
| Untracked frames were scored as confident absence of behaviour and still counted in denominators. | Durations and percentages change wherever tracking dropped out. |
| Interpolation was unbounded and forward-filled leading and trailing gaps. | Long dropouts are no longer fabricated; distance and speed change. |
| Unobserved frames were credited to inverted zones (OFT periphery, EPM arms). | Periphery and arm times fall. |
| Mean speed averaged n frames over n-1 real intervals. | Mean speeds rise by n/(n-1); distances unchanged. |
| Zone `transitions` counted onsets **and** offsets. | New `entries` column is the entry count; `transitions` is retained. |

Every change is documented with old behaviour, why it was wrong, new behaviour
and expected effect in the commit that made it.

---

## Requirements

R >= 4.4, plus:

```r
install.packages(c(
  "sp", "ggplot2", "cowplot", "stringr", "yaml",     # core assays
  "dplyr", "tidyr", "purrr", "readr", "tibble", "zoo", "fs", "rlang",  # OFT, SocInt
  "testthat"                                          # tests
))
```

Optional, for the statistical layer and motif analyses: `lmerTest`, `emmeans`,
`mclust`, `factoextra`, `keras`, `tensorflow`.

Dependencies are declared, never installed at runtime.

---

## Input format

A DLC/SLEAP-style CSV with three header rows and one row per frame:

```
scorer,scorer,scorer,scorer,...
bodyparts,nose,nose,nose,...
coords,x,y,likelihood,...
0,251.3,180.7,0.99,...
```

**Frame numbering must be monotonic and contiguous.** Every displacement is
computed between adjacent rows, so a gap in the numbering would make two frames
recorded seconds apart look adjacent. Non-contiguous files are rejected with
the number of skipped frames; re-export with every frame present, or insert the
missing frames as `NA` rows so they are treated as untracked.

### Required landmarks

| Assay | Landmarks |
|---|---|
| all | four arena corners, by default `tl`, `tr`, `br`, `bl` |
| OFT | `bodycentre` |
| EPM | `bodycentre`, `headcentre`, `neck`, plus every landmark named in the zone file |
| NOR | `nose`, `bodycentre`, `objL`, `objR` |
| SocP | `nose`, `bodycentre`, `socl`, `socr` |
| SocInt | `nose_N`, `bodycentre_N`, `tailBase_N`, `leftEar_N`, `rightEar_N`, `leftSide_N`, `rightSide_N`, `tailEnd_N` for N in 1, 2 |

---

## Running an analysis

Nothing requires editing source code. Copy an example configuration, edit it,
and point the pipeline at it:

```bash
cp config/nor.example.yaml my_nor.yaml
# edit paths, fps, arena size and thresholds

SLEAP_ANALYZER_CONFIG=my_nor.yaml Rscript "02_SLEAPanalzyer/DLCA_NOR v1.2.1.R"
```

The same pattern works for `DLCA_SocP v.1.1.0.R`, `DLCA_EPM v1.0.0.R`,
`DLCA_OFT v1.2.0.R` and `DLCA_SocInt v.0.0.2.r`.

Configurations are validated before anything runs, and every problem is
reported at once rather than one per attempt. Relative paths resolve against
the directory holding the configuration file.

### Outputs

| File | Contents |
|---|---|
| `<file>_output.csv` | Per-animal measurements |
| `combined_output.csv` | All animals in the batch |
| `tracking_qc.csv` | Per-animal tracking quality |
| `run_manifest.yaml` | Commit SHA, timestamp, full configuration, package versions, input file hashes |
| `plots/` | Density paths and overview plots |

The manifest is what makes a run reproducible later. It records whether the
working tree was dirty, because a run made from uncommitted changes cannot be
reproduced from its commit alone.

---

## Calibration

Two mechanisms exist.

**`CalibrateTrackingData()`** applies scalar x/y scaling from a known distance
or arena area. This is what the assay scripts currently use. It is correct when
the camera looks straight down at the arena centre.

**`arena_calibration()`** solves a projective homography from four known arena
corners and maps coordinates into a canonical frame where a rectangular arena
spans `(0, 0)` to `(width_cm, height_cm)`. Under camera tilt the
pixels-per-centimetre ratio varies across the image, and a scalar scale is then
wrong by an amount that depends on where the animal is; the test suite measures
this at over 10% for a modest perspective. It stores the transform, source and
target corners, per-corner reprojection error and a perspective index.

The homography module is tested and available. Migrating the assay scripts onto
it is outstanding work, because it changes calibrated coordinates and so
requires re-validation against existing results.

---

## Quality control

`tracking_qc_report()` reports missing fraction, low-confidence fraction,
longest invalid gap, observed versus interpolated time, implausible
displacement, body-length abnormalities, suspected identity swaps, arena
violations and landmark frame-count consistency.

`qc_flags()` applies thresholds and returns every failing reason. **It never
drops anything.** Excluding a recording is an explicit decision for the
analyst, so a failure warns, records `qcPass = FALSE` with reasons, and still
writes the result.

Every output reports `validTime` alongside `totalTime`. A duration is only
interpretable together with the valid time it came from.

---

## Interpreting the numbers

Metrics fall into four tiers, documented per assay in
[docs/assay_definitions.md](docs/assay_definitions.md):

| Tier | Meaning |
|---|---|
| **Core** | Follows from calibrated coordinates and the event engine |
| **Configurable** | Machinery is tested, but a threshold determines the result and needs assay-specific validation |
| **Experimental** | Heuristic, not validated against manual scoring, suffixed `_experimental` |
| **Exploratory** | Unsupervised or composite, cohort-relative or model-dependent |

Thresholds shipped in the example configurations came from one laboratory's
setup. **They are starting points, not constants.** Validate them against
manually scored video for your apparatus before publishing. The list of what
still needs validation is at the end of
[docs/assay_definitions.md](docs/assay_definitions.md).

Two specific cautions:

* **OFT is a profile, not an anxiety score.** Reduced centre occupancy must be
  read together with locomotion, immobility and tracking quality.
* **SocInt following and avoidance are experimental.** They are heuristics that
  have not been validated as behavioural categories here.

---

## Statistical analysis

Measurement extraction and inference are separate. Extraction produces
measurements and QC; `03_statistics/` decides how to model them, from the
**design of the experiment**, never from which packages are installed.

```r
source("03_statistics/design.R")
source("03_statistics/models.R")

results <- analyze_responses(
  measurements,
  responses = c("contactNov", "contactFam", "latency"),
  group = "treatment", animal_id = "ID", sex = "sex", batch = "batch",
  adjust_method = "holm"
)
```

A random intercept for animal is included only when at least one animal
contributes more than one observation, because `(1|ID)` with one observation
per animal is unidentifiable and `lme4` refuses it. Time bins declared with
`within =` are not treated as independent animals. Multiplicity adjustment
requires an explicit family label and records the family size.

Nothing in this layer interprets a p-value. See
[docs/statistical_layer.md](docs/statistical_layer.md).

---

## Repository layout

```
01_SLEAPcoords/        Notebooks converting SLEAP HDF5 to CSV
02_SLEAPanalzyer/      Assay scripts and the shared core
  core/                Measurement core (see below)
  DLCA_*.R             Per-assay batch workflows
03_statistics/         Experiment-level inference, separate from extraction
config/                Documented example configurations
docs/                  Architecture and definitions
tests/testthat/        Automated test suite
tools/                 Parse check and static audit, run in CI
90_Testing/            Legacy exploratory scripts, not maintained
99_deprecated/         Retained for provenance only
```

### The measurement core

| Module | Responsibility |
|---|---|
| `validation.R`, `units.R` | Argument checks, unit conversion, seconds/frames |
| `tracking_data.R` | TrackingData accessors and per-frame validity |
| `geometry.R` | Distances, angles, polygons, self-intersection |
| `interpolation.R` | Bounded gap filling with observed/interpolated/invalid status |
| `homography.R` | Projective arena rectification |
| `dyadic.R` | Inter-animal geometry and directional relative motion |
| `events.R` | The single event/bout segmentation engine |
| `qc.R` | Canonical quality control |
| `config.R`, `assay_config.R` | Declarative configuration |
| `provenance.R` | Run manifests |

---

## Tests

```bash
Rscript tests/testthat.R      # full suite
Rscript tools/check_parse.R   # every R file parses
Rscript tools/audit_repo.R    # static audit ratchet
```

The audit is a ratchet: it records a budget per problem class (absolute paths,
unbounded interpolation, legacy transition counting, boundary-excluding zone
tests, runtime installation) and fails when a count rises. Budgets are lowered
as problems are fixed.

Integration tests run the NOR, SocP and EPM batch scripts as subprocesses
against synthetic tracking with analytically known answers, so the tests cover
the real entry points rather than the library alone.

CI runs the parse check, the suite and the audit on every push and pull request
to `main` and `v2`. Optional heavy dependencies are excluded from the baseline
so it stays fast and does not break on unrelated upstream changes.

---

## Limitations

* Assay scripts still use scalar calibration; the homography is available but
  not yet wired in.
* SocInt thresholds are in pixels and specific to one camera setup.
* No metric in this repository has been validated against manually scored video.
* Motif analyses (PCA/GMM, autoencoder) are exploratory and have not been
  assessed for seed sensitivity, bootstrap stability or held-out-animal
  generalisation.
* SLEAP HDF5 is read via the conversion notebooks in `01_SLEAPcoords/`; there
  is no native HDF5 import, so SLEAP track identity and confidence beyond the
  CSV columns are not carried through.

## License and contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) and
[CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).
