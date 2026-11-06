# Validation: SLEAPanalyzer v2 against manual BORIS scoring (Exp9)

The evidence behind several parameter choices that ship in `config/`. Run against the
Exp9 social instability stress cohort — 117 animals, six batches, four assays (EPM, OFT,
NOR, social preference) — with manual BORIS scoring as the comparator.

> **Runnable, with one exception.** Paths resolve from each script's own location, and
> the derived tables in `enriched/` are versioned, so every figure regenerates from this
> checkout:
>
> ```bash
> Rscript scripts/05_correlate.R        # and 10-18
> ```
>
> The ~4.8 GB of SLEAP tracking output those tables were built from (`sleap_input/`,
> `sleap_output/`, `sleap_output_all/`) is deliberately not versioned. Scripts `01`-`04`
> therefore cannot be re-run from a clone alone. Scripts `07`-`10` can be run against
> the external Exp9 inputs through the release builder below.

## How the Exp9 scores are made

Each test was filmed from above (GoPro, 29.97 frames per second; the scripts count one
tracking row as 1/30 s), tracked with SLEAP (separate single-animal models for the mouse and
for the fixed parts of the apparatus) and scored by its assay script with the `*_all.yaml`
configuration in `config/`. Coordinates are scaled to cm per video from the tracked
arena, and zones and contact regions are placed from tracked landmarks, never drawn by
hand. Definitions are in [docs/assay_definitions.md](../../docs/assay_definitions.md).

| Test | Script | Scale | Scored | Gaps |
|---|---|---|---|---|
| Open field, 10 min | `DLCA_OFT v1.2.0.R` | 49 x 49 cm floor | body centre in the centre zone (central square, half the side); moving = faster than 5 cm/s by majority over 11 frames | points outside the arena or moving more than 15 cm between frames removed; up to 0.5 s interpolated |
| Plus maze, 5 min | `DLCA_EPM v1.0.0.R` | 60 cm, end to end of the open arms | body centre on open arms, closed arms, centre; head dips = head-centre point beyond the edge, body centre on the maze, neck not in a closed arm, majority over 41 frames | up to 0.2 s interpolated |
| Novel object, 5 min | `DLCA_NOR v1.2.1.R` | 49 x 49 cm | nose in a 9 x 7 cm box at the rectangular novel object or within 4 cm of the round familiar one, body centre more than 1 cm from it, head not turned away (70 degrees or more); objects at fixed positions | up to 0.2 s |
| Social preference, 3 x 10 min | `DLCA_SocP v.1.1.0.R` | 44 x 24 cm | nose within 6 cm of the cup centre; time on top of a cup is not contact | landmark repair first; up to 0.2 s |

Quality control, in order: session QC in every assay (valid time, longest gap,
interpolated time; a failure is recorded, never dropped), the frame-level plausibility
screen (`23`), the light checks of every recording (`31`) and, for the social preference
test, the reviews and validations behind its two corrections (`24`-`30`).

The head-length restoration that `socp_all.yaml` switches on (2026-11-06, see "A head
hidden under the cup lid") is not yet in the released scores; the release of 2026-11-05
(landmark repair) is the one the paper uses.

## Canonical all-batch release bundle

From the repository root in PowerShell, build a new dated, non-overwriting release with:

```powershell
Set-Location 'C:\Users\topohl\Documents\GitHub\SLEAPanalyzer'
.\tools\build_exp9_behavior_bundle.ps1 `
  -ReleaseName '2026-09-22_exp9_all-batches_canonical'
```

The default destination is:

```text
S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Analysis\Behavior\SLEAPanalyzer_v2\releases\<release-name>
```

Each release contains canonical all-animal, male (B1/B2/B5), and female
(B3/B4/B6) tables; pooled and sex-stratified batch-adjusted statistics; eight
regenerated figures as PDF, SVG and PNG (Figures 1-4, 6, and 7-9); compact
assay summaries; QC; exact scripts/configs; stage logs; source and bundle
SHA-256 manifests; and Git/R provenance. Figure 5 is excluded because it has no
generating script.
Figures 10-12 are excluded because they are produced by the superseded phenotype
scripts 14-22. Bulk coordinate files, per-animal plots, and EPM TIFFs remain in
their recorded source locations. The builder refuses to overwrite an existing
release name.

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

**`23` — plausibility screen.** Every check of `core/plausibility.R` on every staged
recording of the four assays (report only): `results/plausibility_episodes.csv`,
`results/plausibility_summary.csv`, with what the landmark repair would remove and
rebuild in each recording and, for the novel-object test, contact before and after it.
Review any episode with `tools/qc_review.py`. See "Plausibility checks" in
`docs/assay_definitions.md`.

**`24` — QC report.** The tracking QC of the social preference test as one workbook,
the same tables as CSV and figures: per recording what was removed, rebuilt, left
untracked or excluded on a cup, contact before and after against a reference release,
the anchor models, the reviewed decisions, the batch 1 manual check and the screen of
`23` (report only). `Rscript 24_qc_report.R <socp_output_root> <reference_release> <out_dir>`.

**`25` — repair review.** Picks repaired episodes of the social preference test (the
longest rebuilt noses, the named cases, every collapsed and the longest split
skeletons), writes those recordings' tracking after the repair and the lists for
`tools/qc_review.py`, which draws the original and the repaired nose and body centre
on stills and clips.

**`26` — frame agreement.** The social-preference contact against the raw manual events of
batch 1 frame by frame: 19 animals (camera 2, where the head length restoration was developed;
camera 1, held out), habituation and both phases, each record aligned on the ear midpoint and
compared over the frames it covers (four records cover a shifted window). With the nose as
tracked and with the head length restored (`restore_head_length()`): agreement by camera and cup
with 95% intervals over animals, where contact is missed or added, along the cup lid, against the
contact distance; the one trial scored twice by hand; bout counts, mean bout and latency against
the automated merge gap and minimum bout. `Rscript 26_frame_agreement.R <reference_release>`;
results in `results/socp_frame_*.csv` and `results/socp_bout_*.csv`. See "A head hidden under
the cup lid" in `docs/assay_definitions.md`.

**`27` — figure stills.** One video frame per kind of misplaced nose (a fixed spot at the
image edge, an object outside the arena, a head under the cup lid), with the nose as tracked
and as the analysis uses it: this script finds the coordinates, `tools/figure_stills.py`
draws them. `Rscript 27_figure_stills.R <out_dir>`, then
`python -I tools/figure_stills.py <out_dir>/stills.csv <out_dir>`.

**`28` — cup geometry.** Calibration and cup geometry of every social-preference recording
as the tracking sees them: pixels per cm, the arena's aspect ratio, the distance between the
cup points and their jitter, head and body length, and the apparent cup radius (where observed
noses stop on the open side of the cup). The cups are the same in every batch, so their apparent
radius checks that one contact distance means the same everywhere (2026-11-05: 3.8 to 4.15 cm by
batch, no difference between the sexes; a per-recording correction for it agreed worse with the
manual scoring). `Rscript 28_cup_geometry.R`; results in `results/socp_cup_geometry.csv`.

**`29` — head restoration impact.** Contact with each cup in every social-preference recording
without and with the head length restoration as configured, the median head length and the time
restored (report only). `Rscript 29_head_restore_impact.R`; results in
`results/socp_head_restore_impact.csv`.

**`30` — nose rebuild benchmark.** How well the landmark repair rebuilds a missing nose where the
truth is known: stretches of observed nose (0.3 to 10 s) are hidden and rebuilt with several
time constants for the pose of the moment (`landmark_repair_pose_decay_s`). `Rscript
30_nose_rebuild_benchmark.R`; results in `results/socp_nose_rebuild_benchmark.csv`.

**`31` — light checks.** One row per recording of every test: its length against the
protocol, an open field last minute too short to count, the same trial twice (identical
files, or two body-centre tracks that coincide at some offset), landmarks the test uses
that are missing, arena shape, corner drift, corner loss and camera shift, the length and
nose-tail order of the moving mouse, poses taken from another video frame, the novel object
positions, plus maze time off every zone, open field time the scoring dropped, novel object
total contact, and readouts more than 3.5 robust SD from their test and sex. Report only.
`Rscript 31_light_checks.R` (about 20 min); results in `results/light_checks_*.csv`. The
2026-11-06 run and its review on video found three input errors (two social preference
trials from the wrong part of a video or from another animal, one plus maze video cut 39 s
short), one open field recording whose poses run out of step with its video, and the
partial last bin in the open field habituation slope (fixed in `DLCA_OFT v1.2.0.R`); see
`Behavior/SLEAP/qc/light_checks` on the lab share. The outcome of each flag looked at is in
`config/light_check_decisions.csv` (keep, correct, code_fixed or open, with the finding, the
evidence and the md5 of the tracking file); `31` shows it under `reviewed` and lists only the
flags without one, or on a tracking file changed since, under `to_review`.

**`32` — video QC manifest.** The manifest for `tools/video_qc` (see the main README): every
recording `31` checks, with the clip it was tracked from (matched by name; eight assigned by
hand, with the evidence in the script), its test folder for the raw camera files, protocol
length, arena size and objects. `31` then lists the tool's flags as `video:<check>:<flag>`.
`frames` and `objects` take about 25 min, `sync` about 1.5 h and `clip` several hours; each
can run on its own into the same folder.

**Running the checks.** `tools/run_exp9_qc.ps1` runs `32`, the video checks and `31`, copies
the results for review to the lab share and ends with the flags to review (about 45 min;
`-SkipVideo` for the light checks only, `-Checks sync` or `-Checks clip` to add a slower check).

## Layout

```
scripts/    34 R scripts, as described above
config/     assay YAMLs used for the Exp9 runs (EPM, NOR, SocP, OFT) and the QC review decisions
metadata/   animal metadata, identifier crosswalks, conflict log
enriched/   assembled wide/long tables, including sleap_all_batches_wide.tsv
results/    derived tables — sweeps, agreement statistics, audits; video_qc/ from tools/video_qc
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

## Figure style

All figures use `scripts/00_theme.R`, which is deliberately identical to `theme_exp9()`
in `SISanalyzer/exp9_publication/R/00_setup.R` — same palette, type size, gridline
convention and output widths — so the validation and publication sets read as one.

| | |
|---|---|
| Palette | navy `#3F4576` / grey `#C2C2C2` / coral `#F4636E`; series extend with teal and mauve |
| Type | 7 pt sans-serif |
| Widths | 183 mm double-column, 120 mm, 89 mm single |
| Output | vector PDF for submission + editable SVG (svglite, live text) + 600 dpi PNG preview, from one `save_fig()` call |
| Non-data ink | no panel borders, no axis lines, no tick marks; gridlines on the value axis only |

Navy and coral differ in hue *and* lightness, so they stay separable under deuteranopia
and protanopia (where the coral reads as tan), with grey between them.

### `fig5_all_batches_matrix.png` is an orphan

No script in `scripts/` produces it — it came from a figure block that was removed
from the code before this study was archived. It is the only figure here that cannot be
regenerated, and it is the only one still at the old style and size. Kept because the
cross-cohort correlation structure it shows is quoted in the analysis record; treat its
provenance as unverified.

