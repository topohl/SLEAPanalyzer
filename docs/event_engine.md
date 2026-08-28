# Shared event engine

`02_SLEAPanalzyer/core/events.R` is the single behavioral event/bout
segmentation implementation. Every assay must route bout, duration, latency and
entry calculations through it so the definitions cannot drift apart.

## Entry point

```r
segment_events(event, fps, valid = NULL, min_bout_s = 0, max_gap_s = 0)
```

| Argument | Meaning |
|---|---|
| `event` | per-frame logical behavior vector |
| `fps` | frames per second, used for every seconds/frames conversion |
| `valid` | per-frame validity mask; `NULL` means every frame was observed |
| `min_bout_s` | minimum bout duration, **in seconds** |
| `max_gap_s` | maximum interior gap bridged before filtering, **in seconds** |

Thresholds are specified in seconds and converted with `fps`, so one
configuration is valid at any acquisition rate. `min_bout_frames` /
`max_gap_frames` are available when an exact frame count is required.

## Semantics

### Validity is not absence

`NA` in an event vector means *unknown*, not *no behavior*. Frames marked
invalid are:

* never counted as behavior,
* never counted as analyzed time.

This is the difference between "the animal did not investigate the object" and
"we could not see the animal". The old code collapsed both to `FALSE`.

### Two vectors, deliberately

`segment_events()` returns both:

* `event` — the **segmented** vector, which defines bout *identity*. A behavior
  briefly interrupted by a tracking dropout is one episode, not two.
* `observed_event` — the **observed** vector, which never includes a frame that
  was not observed.

Durations, percentages and bout statistics are computed from `observed_event`,
so bridging a gap can join two bouts but can never inflate a duration. The bout
table reports both `duration_seconds` (the span of the episode) and
`observed_seconds` (how much of it was actually seen).

The sum of `observed_seconds` over bouts always equals `duration_s`.

### Order of operations

1. apply the validity mask,
2. bridge interior gaps up to `max_gap_s`,
3. drop bouts shorter than `min_bout_s`,
4. tabulate.

Bridging happens *before* filtering so that one episode split by a dropout is
not discarded as two sub-threshold fragments. Leading and trailing gaps are
never bridged: the engine will not invent behavior before the first or after
the last observation.

### No-event behavior

When no bout occurs, `latency_s` is `NA_real_`. It is never `Inf` and never
silently set to the trial duration. This is a **right-censored** observation
and must be handled as such in the statistical layer — a survival model, or an
explicit censoring convention, not a mean over `NA`-dropped values.

`percent_valid_time` is `NA_real_` when there is no valid time at all, and `0`
when there is valid time but no behavior. These are different situations.

## Entries versus transitions

| Function | Counts |
|---|---|
| `count_entries()` | observed onsets only |
| `count_state_changes()` | onsets **and** offsets (legacy `CalculateTransitions()`) |

`count_state_changes()` is roughly twice the entry count and is biased by
whether the animal starts or ends inside the zone. It is retained only so that
summaries produced with the legacy definition can be reproduced; new work
should use `count_entries()`.

`count_entries()` compares consecutive **observed** frames, so a dropout in the
middle of a zone visit does not fabricate an extra entry. A visit already in
progress at the first observed frame is not counted as an entry, because its
onset was never observed; pass `count_initial = TRUE` to change that.

Note that `count_entries()` and the `n_bouts` reported by `segment_events()`
can legitimately differ by one: an episode in progress at the first frame is a
real observed bout, it simply has no observed onset.

## Reported quantities

`segment_events()` returns `duration_s`, `percent_valid_time`, `n_bouts`,
`latency_s`, `mean_bout_s`, `median_bout_s`, `max_bout_s`,
`interbout_intervals_s`, `valid_time_s`, `total_time_s`, `valid_frames` and
`total_frames`.

Always report `valid_time_s` alongside any duration. A 40 s investigation out
of 300 s of valid tracking is not the same result as 40 s out of 600 s.

`event_summary_row(segmentation, prefix)` flattens the summary into a one-row
data frame with prefixed column names for assay output tables.
