# Assay definitions

What each measurement actually is, in what units, and how far it can be
trusted. Every threshold below is a **configuration value**, not a constant of
nature: the defaults came from this laboratory's setup and must be re-validated
against manually scored video for any new apparatus, camera or cohort.

## Metric confidence tiers

Every reported quantity falls into one of four tiers. The tier is the first
thing to check before using a number.

| Tier | Meaning |
|---|---|
| **Core** | Geometric or temporal quantity that follows directly from calibrated coordinates and the event engine. Correct if the tracking and calibration are correct. |
| **Configurable** | A behavioural definition with a threshold. The machinery is tested, but the threshold determines the result and needs assay-specific validation. |
| **Experimental** | A heuristic that has not been validated against manual scoring in this repository. Named with an `_experimental` suffix. Do not present as an established behaviour. |
| **Exploratory** | Unsupervised or composite output with no ground truth. Cohort-relative or model-dependent. |

## Units

| Quantity | Unit |
|---|---|
| coordinates after calibration | cm |
| coordinates before calibration | px |
| distance | cm |
| speed | cm/s |
| duration, latency, bout length | s |
| angle | degrees |
| `fps` | frames per second of the **exported tracking file** |

`TrackingData$distance.units` records the coordinate unit. Assay functions
compare it against the unit their thresholds are expressed in and refuse to run
on a mismatch, so centimetre thresholds cannot be applied to pixel coordinates.

The per-frame `speed` column produced by `CalculateMovement()` holds
**displacement per frame**, not speed. Summing it gives a distance; multiplying
its mean by `fps` gives a speed. Its first element is `NA`, because a frame
with no predecessor has no displacement.

## Shared machinery

All assays use the same event engine, interpolation, QC and geometry. See
[event_engine.md](event_engine.md) for bout, latency and entry semantics.

### Validity

A frame is **valid** when every landmark a measurement depends on has finite
coordinates and, if a cutoff is configured, adequate confidence. Invalid frames
are never counted as behaviour and never counted as analysed time. Every assay
reports `validTime` alongside `totalTime`.

**A duration is only interpretable together with the valid time it came from.**
40 s of investigation out of 300 s valid is not the same result as 40 s out of
600 s.

### Interpolation

Interior gaps up to `max_interpolation_gap_s` are linearly interpolated and
labelled `interpolated`. Longer gaps and *all* leading and trailing gaps stay
missing. Nothing forward-fills, so a dropout never becomes a stationary animal.

### Entries

`<zone>.entries` counts observed onsets. `<zone>.transitions` applies the
legacy onsets-plus-offsets formula and is roughly double; it is retained for
continuity only.

---

## Open Field Test (OFT)

| Metric | Tier | Definition |
|---|---|---|
| `distance` | Core | Sum of per-frame displacement, cm |
| `raw.speed` | Core | Mean per-frame displacement x fps, cm/s |
| `center/periphery/corner .total.time` | Core | Occupancy of the zone, s |
| `*.entries` | Core | Observed onsets |
| `mean_wall_distance_cm` | Core | Mean distance to nearest wall |
| `occupancy_entropy` | Configurable | Shannon entropy over an N x N grid; depends on `occupancy_grid_n` |
| habituation `slope_per_min_*` | Configurable | Within-animal regression on time bins |
| `p_value_within_animal_nominal_*` | Exploratory | See below |
| `oft_center_exploration_score_cohort_z_experimental` | Exploratory | See below |

Zones are derived from the four tracked corners by scaling: `center_scale`
0.5, `periphery_scale` 0.8, `corner_scale` 0.4. These are geometric
conventions, not biological ones.

**Interpretation.** OFT metrics form a behavioural profile. Reduced centre
occupancy is not a measurement of anxiety. It must be read together with
locomotion, immobility, tracking quality and the time course of exploration; an
animal that moves little occupies the centre little for reasons that have
nothing to do with avoidance.

**The composite score is cohort-relative.** Each component is z-scored against
the other animals in the same run, so the same recording scores differently
depending on which animals it was analysed with, and scores cannot be compared
across batches. It is a within-cohort ranking aid. Model the raw components
instead.

**Habituation p-values are within-animal and nominal.** They come from
regressing one animal's metric on its own time bins, treating consecutive bins
as independent. Locomotion is strongly autocorrelated across adjacent bins, so
the standard error is underestimated. Use them as fit diagnostics. Group
differences in habituation belong in the statistical layer, with bins nested
within animal.

---

## Elevated Plus Maze (EPM)

| Metric | Tier | Definition |
|---|---|---|
| `bodycentre.open/closed/center .total.time` | Core | Arm occupancy, s |
| `*.entries` | Core | Observed onsets |
| `bodycentre.raw.distance` | Core | Path length, cm |
| `nose.dip` | Configurable | Head outside the maze outline, body on the maze, neck not in a closed arm |

Zones come from `EPM_zoneinfo.csv`, which lists the tracked landmarks forming
each zone's outline **in order around the perimeter**. Listing them diagonally
produces a self-intersecting polygon whose point-in-polygon test silently
reports interior frames as outside; `AddZones()` now refuses such a zone.

Unobserved frames belong to no zone, including inverted zones. This matters
because the periphery-style inverted test would otherwise credit every dropout
to the complementary zone.

`nose.dip` counts whole onsets. It is a geometric proxy for head-dipping and
should be validated against manual scoring before use as a primary outcome.

### Nose dips over-count at the default smoothing window

That validation has now been done once, against the manually scored Exp9
Batch 1 (n = 20 — the only batch with manual nose-dip scoring).

The detector has no threshold: a frame is a dip when the head is outside the
arena outline, the body is inside it, and the neck is not in a closed arm. The
only lever is the half-width of the majority filter applied before onsets are
counted, `nosedip_integration_period`, which decides how long a head must stay
over the edge to count as one dip rather than several.

| window (frames) | SLEAP / manual | Lin's CCC | Spearman rho |
|---|---|---|---|
| 1 | 2.46 | 0.27 | 0.82 |
| **5 (default)** | **1.96** | **0.42** | **0.84** |
| 10 | 1.71 | 0.54 | 0.83 |
| 15 | 1.43 | 0.70 | 0.82 |
| **20** | **1.13** | **0.82** | **0.86** |
| 30 | 0.76 | 0.68 | 0.82 |
| 45 | 0.49 | 0.40 | 0.76 |

At the shipped default the detector reports **twice** as many dips as a human
scorer. A window of 20 frames (±0.67 s at 30 fps) brings it to within 13% and
roughly doubles concordance. Rank agreement is flat across the whole sweep, so
the window trades count magnitude, not ordering: `nose.dip` is usable as a
ranked measure at any setting, and as an absolute count only near 20.

**The default is deliberately unchanged.** One cohort, one apparatus and one
scorer is thin ground for changing a shipped default, and the right window
depends on frame rate and on how the scorer segments dips. Set it explicitly:

```yaml
integration_period_frames: 5           # locomotion: moving/stationary split
nosedip_integration_period_frames: 20  # dips: calibrated against manual scoring
```

The two were a single parameter before, which made this untunable — widening
the dip window also blurred the moving/stationary split. They are now separate,
and `nosedip_integration_period_frames` defaults to `integration_period_frames`
so existing configurations are unaffected.

### Calibration: distance, not arena area

EPM is calibrated from a **measured distance between two landmarks**
(`calibration_method: distance`), and that is the default:

```yaml
calibration_method: distance
calibration_points: [tl, bl]
calibration_distance_cm: 60
```

Pick the two landmarks on the **same edge of one arm**, so the value you
declare is the arm-axis span you measured. `tl` and `bl` are the left corners
of the top and bottom arms, so `tl`-`bl` runs tip to tip along the left edge of
the vertical arm.

Pairing diagonally opposite corners such as `tl`-`br` is a subtler version of
the same mistake as the area trap below: that distance is
`sqrt(span² + arm_width²)`, so calling it the span under-scales the maze by
roughly `arm_width² / (2 · span²)` — about 0.3% for a 5 cm arm over 60 cm, and
larger on a wider arm. It is small enough to survive review and large enough to
bias a reported path length.

The other assays calibrate an area against `arena_corner_names`, which is
sound for them because `tl, tr, br, bl` really are the corners of a
rectangular arena — on OFT, NOR and SocP the polygon through those four
points matches their bounding box to within about 1%.

**The plus maze breaks that assumption.** There, `tl`/`tr` and `bl`/`br` are
the outer corners of the two *opposing arms*, so the quadrilateral through
them is a narrow corridor along one arm axis, not the maze. On a maze whose
arms are 5 cm wide and span 60 cm, that polygon covers roughly a *ninth* of
the plus outline's bounding box. Equating it to `arena_width_cm *
arena_height_cm` therefore understates the pixel area, and since
`px.to.cm = sqrt(metric_area / pixel_area)`, it inflates the scale — by about
3.4x on the reference Batch-1 recordings. Every distance, speed and the
`movement_cutoff_cm_s` threshold is wrong by that factor, and nothing in the
output reveals it: occupancy times and entry counts are unaffected, so the
run looks healthy.

If you do calibrate an area on EPM, pass the full maze outline via
`calibration_points` **in perimeter order** (the `arena` column of
`EPM_zoneinfo.csv` is exactly that order) and set
`arena_width_cm * arena_height_cm` to the true area of the plus — which is
`2 * arm_width * tip_to_tip - arm_width^2`, not the bounding square.

A useful sanity check on any calibration: multiply a known short landmark
distance, such as the arm width `tl`-`tr`, by the resulting `px.to.cm` and
confirm it returns the measured width.

---

## Novel Object Recognition (NOR)

| Metric | Tier | Definition |
|---|---|---|
| `contactLeft` / `contactRight` | Configurable | Nose in the contact region, body beyond `body_exclusion_distance_cm`, head oriented toward the object |
| `contactNov` / `contactFam` | Configurable | The above, mapped through novel-location metadata |
| `proxLeft` / `proxRight` | Configurable | Nose distance within `proximity_range_cm` |
| `latency*` | Core | Time to first contact bout, `NA` if none |
| `frequency*`, `entries*`, `meanBout*` | Core | From the shared event engine |
| `frequencyRear_experimental` | Experimental | Spine compression proxy for rearing |

### The contact detector is symmetric

`contact_geometry` selects **one** detector applied identically to both
objects:

* `radial` — nose within `contact_distance_cm` of the object point (default)
* `box` — nose inside a `object_box_width_cm` x `object_box_height_cm` footprint
* `legacy_asymmetric` — reproduces the pre-v2 behaviour and **warns**

Before v2, the novel object was always scored with the 9x7 cm box and the
familiar object with the 4 cm radius. Those regions have different areas
(63 vs 50.3 cm²), so on an apparatus where the two objects are physically
**identical** the asymmetric detector inflates novel contact time by
construction and biases the discrimination index independently of behaviour.

### But validate it: on Exp9 the shape-matched detector wins

That argument assumes the objects are interchangeable. On the Exp9 apparatus
they are not — one is round and small, the other rectangular — and the
pre-v2 dimensions describe the real objects. Scored against the manually
scored Batch 1 (n = 20, the only batch with manual NOR scoring):

| detector | novel ratio | familiar ratio | D2 r | D2 CCC | D2 bias |
|---|---|---|---|---|---|
| `radial`, 4 cm | 0.93 | 1.01 | 0.740 | 0.715 | **−0.058** |
| `radial`, 5 cm | 1.11 | 1.18 | 0.800 | 0.791 | −0.022 |
| **`legacy_asymmetric` 9x7 / r4** | **1.02** | **1.01** | 0.845 | 0.842 | **+0.009** |
| 9x7 box / area-matched r4.48 | 1.02 | 1.10 | **0.860** | **0.848** | −0.032 |

A single radius applied to both **under-detects the larger object** (0.93
against 1.01), and that asymmetry is itself a D2 bias — a bigger one
(−0.058) than the area difference it was meant to avoid (+0.009). Matching
the detector to each object recovers contact time for both at ~1.0 and nearly
doubles the reduction in D2 bias.

So `legacy_asymmetric` is the validated choice **for this apparatus**, and the
name is misleading: it is shape-matching, not a legacy quirk. The warning it
emits should be read as "confirm this matches your objects", not "do not use".

The general rule stands: **the detector must match the apparatus, and which
detector that is can only be settled against manually scored video.** Identical
objects → one detector. Physically different objects → match each, and check
that both recover contact time at a ratio near 1, because an uneven recovery
is what actually biases D2.

### Orientation convention

`contact_angle_deg` is the angle between the **body-centre-to-nose** vector and
the **object-to-nose** vector. Under this convention an animal looking straight
at the object gives ~180°, and the default `[70, 290]` accepts anything not
pointing away. Because `vector_angle_degrees()` returns 0–180, the upper bound
is inert; the effective criterion is `>= 70`. This is permissive and is
retained from the legacy definition; it needs validation.

### Novel/familiar mapping

Metadata `"R"` maps the **left** side to novel and `"L"` maps the **right**
side to novel. This inverted-looking convention is inherited and preserved.
**Confirm it against your metadata before interpreting `contactNov`.**

Contact no longer depends on this metadata: a missing row leaves `contactNov`
and `contactFam` as `NA` while `contactLeft` and `contactRight` remain valid.

`frequencyRear_experimental` treats simultaneous spine1 and spine2 compression
as rearing. It has not been validated against manual scoring.

---

## Social Preference (SocP)

| Metric | Tier | Definition |
|---|---|---|
| `contactLeft` / `contactRight` | Configurable | Nose within `contact_distance_cm` of the stimulus, body beyond the exclusion distance |
| `contactNovel` / `contactFamiliar` | Configurable | The above, mapped through metadata |
| `proxNovel` / `proxFamiliar` | Configurable | Nose distance within `proximity_range_cm` |
| `latency*`, `frequency*`, `entries*` | Core | Shared event engine |

The same detector is applied to both chambers. Unlike NOR, no head-orientation
criterion is applied by default; this preserves the established SocP
definition. `require_orientation = TRUE` adds one, which **changes the
definition** and needs its own validation.

The novel/familiar mapping follows the same inherited convention as NOR.

---

## Social Interaction (SocInt)

SocInt is dyadic. Symmetric measures describing the pair are reported
separately from directed measures attributed to an individual.

### Units

SocInt runs in **pixels** by default. `validate_socint_units()` refuses to run
if `threshold_unit` disagrees with the analysed coordinate unit, and refuses
centimetre mode until `scientific_cm_thresholds_confirmed` is set, because
every distance threshold would need converting and re-validating first.

### Directional relative motion

Each animal's velocity is projected onto the inter-animal axis, splitting the
observed closing speed into the part each animal contributed:

```
closing = a1_closing_speed + a2_closing_speed
```

`a1_closing_speed` is positive when animal 1 moves toward animal 2.

| Metric | Tier | Definition |
|---|---|---|
| `a1_approaches_a2`, `a2_approaches_a1` | Configurable | That animal's closing contribution exceeds `movement_cutoff` |
| `both_approach` | Configurable | Both contributions positive |
| `a1_retreats_from_a2`, `a2_retreats_from_a1` | Configurable | Contribution below `-movement_cutoff` |
| `separation_led_by_a1` / `_a2` | Configurable | Pair separating, attributed to the larger contributor |
| `approach_event`, `retreat_event` | Configurable | Symmetric, pair-level |
| `a1_nose_to_nose2` and similar | Configurable | Directed investigation: distance plus facing angle |
| `*_following_*_experimental` | Experimental | Not validated |
| `*_avoidance_*_experimental` | Experimental | Not validated |

Before v2 the avoidance rules thresholded each animal's **scalar** speed, which
carries no direction: whenever both animals moved, both avoidance flags fired
regardless of who produced the separation, and two animals walking in parallel
with no relative motion satisfied every rule. Following and avoidance now
require the correct directional term but remain heuristics, hence the
`_experimental` suffix.

---

## What needs empirical validation

Nothing below has been validated against manually scored video **in this
repository**. Each requires it before use as a primary outcome.

1. **NOR contact distance and orientation angle** — against your object size,
   camera height and species.
2. **NOR proximity band** (`proximity_range_cm`).
3. **SocP contact distance and proximity band** — against your apparatus.
4. **SocInt contact, proximity, follow and avoidance distances** — currently in
   pixels and specific to one camera setup.
5. **SocInt following and avoidance** as behavioural categories.
6. **EPM nose-dip geometry.**
7. **NOR rearing proxy** (`frequencyRear_experimental`).
8. **The inherited novel/familiar mapping**, which is not self-evidently correct.
9. **`max_plausible_speed_cm_s`** — should reflect the species' actual maximum.
10. **`movement_cutoff_cm_s`** — separates locomotion from tracking jitter and
    depends on tracking noise.
