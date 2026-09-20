# Validation: SLEAPanalyzer v2 against manual BORIS scoring (Exp9)

The evidence behind several parameter choices that ship in `config/`. Run against the
Exp9 chronic social-defeat cohort — 117 animals, six cohorts, four assays (EPM, OFT,
NOR, social preference) — with manual BORIS scoring as the comparator.

> **These scripts are not runnable from this repository.** Their inputs are ~4.8 GB of
> SLEAP tracking output (`sleap_input/`, `sleap_output/`, `sleap_output_all/`) which is
> deliberately not versioned, and their paths still point at the working directory they
> were written in. What is preserved here is the **code, the derived tables and the
> figures** — the audit trail for the numbers quoted below. Treat them as a record, not
> as a pipeline.

---

## What this established

| Finding | Where |
|---|---|
| SLEAP vs manual agreement: **r = 0.95–0.99, Lin's CCC 0.93–0.97** on durations | `05_correlate.R`, `results/method_agreement.csv` |
| NOR contact distance **4 cm** chosen by sweep — ratio 0.96, CCC 0.939 | `12_calibrate_nor_contact.R`, `results/nor_contact_sweep.csv` |
| Nose-dip integration window **20 frames**, not 5 — the 5-frame default over-counted 1.96× (CCC 0.42) vs 1.13× (CCC 0.82) | `11_recalibrate_nosedips.R`, `results/nosedip_window_sweep.csv` |
| Shape-matched NOR detector (round vs square objects) | `13_compare_nor_detectors.R` |
| NOR spreadsheet header swap found by reconciling against raw BORIS exports | `06_validate_nor_against_raw.R` |
| 574 files staged across six cohorts, 0 unresolved identifiers | `08_stage_all_batches.R`, `results/staging_report.csv` |

The corresponding fixes are committed to the analysis code itself, not here — see the
EPM calibration, OFT, nose-dip and NOR detector commits on `v2`.

## Scripts

**`01`–`13` — tool validation and calibration.** This is why the study lives in this
repository: it is how the shipped parameters were chosen and checked.

| | |
|---|---|
| `01`–`04` | Build animal metadata, enrich the manual assay exports, prepare SLEAP inputs, assemble |
| `05` | **Core comparison**: SLEAP vs manual, per assay and metric (Pearson, Spearman, Lin's CCC, Bland–Altman) |
| `06` | Reconcile the NOR summary sheet against raw BORIS exports |
| `07`–`09` | Pre-flight audit, staging and assembly across all six cohorts |
| `10` | Cross-cohort group analysis |
| `11`–`13` | Parameter sweeps: nose-dip window, NOR contact distance, detector shape matching |

**`14`–`22` — Exp9 phenotype analyses.** Retained for provenance only. They are
**superseded by [`SISanalyzer/exp9_publication/`](../../../SISanalyzer/exp9_publication)**,
which reproduces their conclusions from source with assertions and a single entry point.
They are kept here rather than deleted because several intermediate findings are quoted
in the analysis record: leave-one-out phenotyping (`14`), labelling-scheme comparison and
bootstrap stability (`15`, `21`), classifier reconstruction (`16`, `19`), CON/RES/SUS
characterisation (`17`), within-sex composites (`18`), cohort handling (`20`), and the
effect of the improved NOR pipeline on the classification (`22`).

**Do not run `14`–`22` for new work.** They predate the 20 September 2026 workbook
correction and use the older 39-animal susceptible list.

## Layout

```
scripts/    22 R scripts, as described above
config/     assay YAMLs used for the Exp9 runs (EPM, NOR, SocP, OFT)
metadata/   animal metadata, identifier crosswalks, conflict log
enriched/   assembled wide/long tables, including sleap_all_batches_wide.tsv
results/    33 derived tables — sweeps, agreement statistics, audits
figures/    12 validation figures
```

`enriched/sleap_all_batches_wide.tsv` is the handoff to downstream analysis; a copy is
vendored in `SISanalyzer/exp9_publication/data/raw/`.

## One caveat worth repeating

The SLEAP–manual agreement is **not** two independent methods agreeing. The contact
detector was *calibrated* to match manual scoring (that is what `12` does), so
r = 0.995 on NOR D2 demonstrates that one fixed parameter set reproduces manual scoring
consistently across all six cohorts. That is a reproducibility claim, and a useful one —
but it is not independent validation, and should not be reported as such.
