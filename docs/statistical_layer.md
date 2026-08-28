# Statistical layer

Measurement extraction and inferential statistics are separate. Extraction
(`02_SLEAPanalzyer/`) produces per-animal behavioral measurements and tracking
QC. Inference (`03_statistics/`) decides how those measurements should be
modelled, and decides it from the **design of the experiment**, never from
which packages happen to be installed.

## Why the separation exists

The pre-v2 statistics block inside
`DLCA_OFT_advanced_metrics_extension_v1.0.0.R` chose its model like this:

```r
use_lmm <- requireNamespace("lmerTest", quietly = TRUE) && "ID" %in% names(full_tbl)
```

Three problems follow from that one line.

1. **The inference depends on the library state.** Two researchers analysing
   identical data with different installed packages get different p-values.
   That is not a property any analysis should have.

2. **`(1|ID)` was fitted with one observation per animal.** OFT is one session
   per animal, so every ID appeared exactly once. A random intercept per animal
   is then unidentifiable — its variance cannot be separated from the residual
   variance — and `lme4` refuses it outright:

   > number of levels of each grouping factor must be < number of observations

3. **That error was swallowed.** The `tryCatch(..., error = function(e) NULL)`
   turned the refusal into `"Model failed."` for every metric. So *installing*
   `lmerTest` changed the output from a table of p-values to no results at all,
   with no message explaining why.

The block is retained but deprecated and now warns on use.

## How the model is chosen

`03_statistics/design.R` derives the structure from the data alone:

| Condition | Model |
|---|---|
| every animal contributes exactly one observation | linear model, no random intercept |
| at least one animal contributes more than one | linear mixed model with `(1|ID)` |

Blocking factors (`sex`, `batch`, `cage`, `litter`) enter as fixed effects when
they have more than one level. A declared `within` factor (assay phase, time
bin) enters as a fixed effect when the design is genuinely repeated.

`model_structure()` returns the formula together with a `rationale` string
stating why that structure was chosen. `describe_analysis_plan()` reports the
whole plan without fitting anything, which makes it usable as a
pre-registration artefact.

If a design requires a mixed model and `lmerTest` is absent, `fit_group_model()`
**errors**. It does not quietly fit a model that ignores the repeated-measures
structure.

## Design validation

`validate_experiment_design()` reports, rather than throws, so every problem
with a dataset surfaces at once:

**Errors**

* fewer than three complete observations,
* the grouping column has fewer than two levels,
* an animal appears in more than one treatment group with no declared
  within-animal factor,
* a blocking factor is completely confounded with the group — every batch,
  cage or litter contains only one group, so their effects cannot be
  separated.

**Warnings**

* some animals contribute multiple observations but no `within` factor was
  declared. If these are time bins or repeated sessions, name them; otherwise
  they are treated as exchangeable repeated measures.

### Time bins are not animals

Splitting a 10-minute session into six 1-minute bins produces six rows per
animal. Treating those as independent observations inflates the apparent
sample size sixfold. Declare the bin column with `within = "time_bin"`; the
reported `n_observations` and `n_animals` then stay distinct, and the model
carries `(1|ID)`.

## Multiplicity

`adjust_multiplicity()` requires an explicit `family_label` and records it,
along with `adjustment_method` and `adjustment_family_size`, on every row. An
adjusted p-value without a stated family is not interpretable.

Models that failed are excluded from the family size rather than counted.
Passing `NA` values straight to `stats::p.adjust()` inflates the correction by
the number of failures, which makes the surviving results more conservative for
a reason that has nothing to do with the hypotheses.

## What this layer does not do

It does not interpret p-values. Results carry the statistic, degrees of
freedom, raw and adjusted p-value, the adjustment family, the number of
observations and the number of animals. Whether that constitutes evidence is a
scientific judgement that depends on the pre-registered hypothesis, the effect
size, and the rest of the experiment — not on whether a number fell below 0.05.

## Example

```r
source("03_statistics/design.R")
source("03_statistics/models.R")

measurements <- read.csv("output/combined_output.csv")
measurements <- merge(measurements, read.csv("metadata/animals.csv"), by = "ID")

# State the plan before looking at any p-value.
plan <- describe_analysis_plan(experiment_design(
  measurements, response = "contactNov", group = "treatment",
  animal_id = "ID", sex = "sex", batch = "batch"
))
str(plan)

results <- analyze_responses(
  measurements,
  responses = c("contactNov", "contactFam", "latency"),
  group = "treatment", animal_id = "ID", sex = "sex", batch = "batch",
  adjust_method = "holm"
)
```
