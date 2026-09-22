# ============================================================================
# 10_analyse_all_batches.R
#
# Group comparisons for the full B1-B6 cohort. Effect sizes with confidence
# intervals are the primary output; p-values are secondary.
#
# WHY THIS IS STRUCTURED THE WAY IT IS
#
# An earlier version screened 11 metrics x 3 groupings and corrected with Holm
# across all 33, reported only p-values, and concluded "nothing survives". That
# was wrong in three ways, and each is fixed here:
#
#   1. Holm controls the probability of ANY false positive -- a confirmatory
#      standard applied to an exploratory screen. Benjamini-Hochberg is used
#      instead, with Holm still reported for the small pre-specified set.
#      Uncorrected, 5 of 11 metrics separate CON/RES/SUS; Holm across 33
#      removed all of them.
#
#   2. Reporting p alone hid five medium effects (|d| 0.44-0.67, all pointing
#      the same way). Every contrast now carries a standardised effect and a
#      95% CI, so "not significant" and "no effect" stay distinguishable.
#
#   3. Pooling the sexes cancels opposite-signed effects for some metrics --
#      OFT distance is male p = 0.009, female p = 0.90, pooled p = 0.13. Every
#      contrast is therefore reported pooled AND stratified by sex.
#
# DESIGN CONSTRAINTS, unchanged:
#
#   Sex is perfectly nested in batch (B1/B2/B5 male, B3/B4/B6 female), so a sex
#   MAIN effect is not estimable and none is reported. Stratifying by sex is
#   legitimate -- each stratum still holds 3 batches to block on -- and the
#   phenotype x sex interaction is estimable from the batch-level contrasts,
#   though with only 3 batches per sex it is badly underpowered.
#
#   Condition is estimable: 4 CON in every batch. Batch is blocked throughout,
#   which buys power rather than costing it (it lowered p for 9 of 11 metrics).
# ============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(writexl)
})

# Resolve code from this script's own location. Data and outputs default to the
# validation study beside it, but a release builder can redirect them together
# with EXP9_SLEAP_RUN_ROOT and the more specific directory variables below.
SCRIPTS <- local({
  a <- commandArgs(trailingOnly = FALSE)
  p <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/")
})
DEFAULT_ROOT <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")
RUN_ROOT <- Sys.getenv("EXP9_SLEAP_RUN_ROOT", unset = DEFAULT_ROOT)

# Shared Nature-style theme, identical to the publication figure set.
source(file.path(SCRIPTS, "00_theme.R"))
ENRICH <- Sys.getenv("EXP9_SLEAP_DATA_DIR", unset = file.path(RUN_ROOT, "enriched"))
FIG <- Sys.getenv("EXP9_SLEAP_FIGURES_DIR", unset = file.path(RUN_ROOT, "figures"))
RES <- Sys.getenv("EXP9_SLEAP_RESULTS_DIR", unset = file.path(RUN_ROOT, "results"))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

dat <- read.delim(file.path(ENRICH, "sleap_all_batches_wide.tsv"),
                  stringsAsFactors = FALSE, na.strings = "")
cat(sprintf("cohort: %d animals\n", nrow(dat)))


METRICS <- c(
  "EPM open-arm fraction" = "sleap_EPM_open_frac",
  "EPM centre time"       = "sleap_EPM_center_s",
  "EPM nose dips"         = "sleap_EPM_nosedip_n",
  "EPM distance"          = "sleap_EPM_distance_cm",
  "NOR D2"                = "sleap_NOR_D2",
  "NOR interaction"       = "sleap_NOR_total_s",
  "SocP S1 preference"    = "sleap_SocP_S1_pref_index",
  "SocP S2 preference"    = "sleap_SocP_S2_pref_index",
  "OFT centre %"          = "sleap_OFT_center_pct",
  "OFT corner %"          = "sleap_OFT_corner_pct",
  "OFT distance"          = "sleap_OFT_distance_cm"
)

# One conventional primary readout per assay. These are the standard primary
# measures for each task, so the choice is defensible on prior grounds -- but
# it was made AFTER seeing this data, in the same session, and that is not the
# same as pre-registration. EPM centre time (p = 0.006 for Condition) is not in
# the set and would change the picture if it were. Read the primary tier as
# "the conventional readouts", not as a confirmatory test.
PRIMARY <- c("EPM open-arm fraction", "NOR D2", "SocP S1 preference", "OFT centre %")

# Two-level contrasts, each with a signed direction the column name states.
CONTRASTS <- list(
  list(label = "SIS - CON",  col = "Condition", ref = "CON", alt = "SIS",  term = "ConditionSIS"),
  list(label = "SUS - RES",  col = "Phenotype", ref = "RES", alt = "SUS",  term = "PhenotypeSUS"),
  list(label = "SUS - CON",  col = "Phenotype", ref = "CON", alt = "SUS",  term = "PhenotypeSUS"),
  list(label = "RES - CON",  col = "Phenotype", ref = "CON", alt = "RES",  term = "PhenotypeRES")
)
STRATA <- c("pooled", "Male", "Female")

# --- Estimation -------------------------------------------------------------
# Standardised effect = batch-adjusted group coefficient / residual SD, i.e. a
# partial Cohen's d. The CI is the coefficient CI on the same scale, so it is
# directly comparable across metrics with different units.
estimate <- function(metric_label, metric, ct, stratum) {
  d <- dat %>%
    filter(!is.na(.data[[metric]]), .data[[ct$col]] %in% c(ct$ref, ct$alt))
  if (stratum != "pooled") d <- d %>% filter(Sex == stratum)
  d[[ct$col]] <- factor(d[[ct$col]], levels = c(ct$ref, ct$alt))
  if (nrow(d) < 10 || dplyr::n_distinct(d[[ct$col]]) < 2) return(NULL)
  # Batch must still vary within the stratum for it to be a usable block.
  terms <- if (dplyr::n_distinct(d$Batch) > 1) c(ct$col, "Batch") else ct$col
  fit <- try(lm(reformulate(terms, metric), data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  co <- summary(fit)$coefficients
  if (!ct$term %in% rownames(co)) return(NULL)
  s <- summary(fit)$sigma
  ci <- confint(fit)[ct$term, ]
  n_ref <- sum(d[[ct$col]] == ct$ref); n_alt <- sum(d[[ct$col]] == ct$alt)
  data.frame(
    metric = metric_label, contrast = ct$label, stratum = stratum,
    n = nrow(d), n_ref = n_ref, n_alt = n_alt,
    raw_diff = co[ct$term, 1], d = co[ct$term, 1] / s,
    d_lo = ci[1] / s, d_hi = ci[2] / s,
    p_raw = co[ct$term, 4],
    tier = if (metric_label %in% PRIMARY) "primary" else "exploratory",
    stringsAsFactors = FALSE)
}

est <- bind_rows(lapply(STRATA, function(s)
  bind_rows(lapply(CONTRASTS, function(ct)
    bind_rows(lapply(names(METRICS), function(nm)
      estimate(nm, METRICS[[nm]], ct, s)))))))

# Multiplicity family = the metrics tested within one contrast and stratum.
# That is the set a reader scans together; it does not pool across contrasts,
# which ask different questions of overlapping animals.
est <- est %>%
  group_by(contrast, stratum) %>%
  mutate(p_BH = p.adjust(p_raw, "BH"), family_n = n()) %>%
  ungroup()

primary <- est %>%
  filter(tier == "primary") %>%
  group_by(contrast, stratum) %>%
  mutate(p_holm_primary = p.adjust(p_raw, "holm"),
         p_BH_primary   = p.adjust(p_raw, "BH")) %>%
  ungroup()

write.csv(est, file.path(RES, "all_batches_effects.csv"), row.names = FALSE)
write.csv(primary, file.path(RES, "all_batches_primary.csv"), row.names = FALSE)
write.csv(
  est %>% filter(stratum == "pooled"),
  file.path(RES, "exp9_sleap_effects_all_animals_batch_adjusted.csv"),
  row.names = FALSE
)
write.csv(
  est %>% filter(stratum == "Male"),
  file.path(RES, "exp9_sleap_effects_males_B1-B2-B5_batch_adjusted.csv"),
  row.names = FALSE
)
write.csv(
  est %>% filter(stratum == "Female"),
  file.path(RES, "exp9_sleap_effects_females_B3-B4-B6_batch_adjusted.csv"),
  row.names = FALSE
)

# --- Reporting --------------------------------------------------------------
fmt <- function(x) x %>%
  mutate(across(c(d, d_lo, d_hi), ~round(., 2)),
         across(starts_with("p_"), ~signif(., 3))) %>%
  as.data.frame()

cat("\n=== PRIMARY TIER: the conventional readout per assay (4 metrics) ===\n")
cat("Chosen after seeing the data -- conventional, but not pre-registered.\n")
print(fmt(primary %>% filter(stratum == "pooled") %>%
  select(contrast, metric, n_ref, n_alt, d, d_lo, d_hi,
         p_raw, p_BH_primary, p_holm_primary) %>%
  arrange(contrast, p_raw)), row.names = FALSE)

cat("\n=== EXPLORATORY: all 11 metrics, BH within contrast x stratum ===\n")
for (ct in vapply(CONTRASTS, `[[`, character(1), "label")) {
  sub <- est %>% filter(contrast == ct, stratum == "pooled") %>% arrange(p_raw)
  hits <- sub %>% filter(p_BH < 0.05)
  cat(sprintf("\n-- %s (pooled, n = %d vs %d) : %d of %d survive BH\n",
              ct, sub$n_ref[1], sub$n_alt[1], nrow(hits), nrow(sub)))
  print(fmt(sub %>% select(metric, d, d_lo, d_hi, p_raw, p_BH) %>% head(5)),
        row.names = FALSE)
}

cat("\n=== SEX-STRATIFIED: where pooling changes the answer ===\n")
wide <- est %>%
  select(contrast, metric, stratum, d, p_raw) %>%
  pivot_wider(names_from = stratum, values_from = c(d, p_raw)) %>%
  mutate(opposite_sign = sign(d_Male) != sign(d_Female),
         pooled_weaker = abs(d_pooled) < pmin(abs(d_Male), abs(d_Female)))
flip <- wide %>% filter(opposite_sign | pooled_weaker) %>%
  arrange(contrast, desc(abs(d_Male - d_Female)))
print(flip %>% mutate(across(starts_with("d_"), ~round(., 2)),
                      across(starts_with("p_raw"), ~signif(., 3))) %>%
        select(contrast, metric, d_Male, d_Female, d_pooled,
               p_raw_Male, p_raw_Female, p_raw_pooled) %>%
        as.data.frame(), row.names = FALSE)
cat(sprintf("\n%d of %d contrast x metric combinations flip sign or weaken when pooled\n",
            nrow(flip), nrow(wide)))

# Interaction, reported with its own caveat rather than buried.
cat("\n=== PHENOTYPE x SEX INTERACTION (SUS vs RES) ===\n")
inter <- bind_rows(lapply(names(METRICS), function(nm) {
  v <- METRICS[[nm]]
  d <- dat %>% filter(!is.na(.data[[v]]), Phenotype %in% c("RES", "SUS"))
  an <- try(anova(lm(reformulate(c("Phenotype*Sex", "Sex/Batch"), v), data = d)),
            silent = TRUE)
  if (inherits(an, "try-error")) return(NULL)
  r <- grep("Phenotype:Sex", rownames(an), value = TRUE)
  if (!length(r)) return(NULL)
  data.frame(metric = nm, p_interaction = an[r, "Pr(>F)"], stringsAsFactors = FALSE)
}))
inter$p_BH <- p.adjust(inter$p_interaction, "BH")
print(inter %>% arrange(p_interaction) %>%
        mutate(across(starts_with("p_"), ~signif(., 3))) %>% as.data.frame(),
      row.names = FALSE)
cat("Underpowered: 3 batches per sex, so the batch-level error term is thin.\n",
    "Treat sign flips as a signal to investigate, not an established finding.\n", sep = "")
write.csv(inter, file.path(RES, "all_batches_sex_interaction.csv"), row.names = FALSE)

# --- Figure: forest plot of effect sizes -----------------------------------
fp <- est %>%
  filter(contrast %in% c("SIS - CON", "SUS - RES")) %>%
  mutate(metric = factor(metric, levels = rev(names(METRICS))),
         stratum = factor(stratum, levels = c("pooled", "Male", "Female")))

p1 <- ggplot(fp, aes(d, metric, colour = stratum)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.5) +
  geom_linerange(aes(xmin = d_lo, xmax = d_hi),
                 position = position_dodge(width = 0.65), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.65), size = 1.9) +
  scale_colour_manual(values = c(pooled = SER1, Male = PAL[["Male"]], Female = PAL[["Female"]]),
                      name = NULL) +
  facet_wrap(~ contrast, scales = "free_x") +
  labs(title = "Effect sizes, not just p-values",
       subtitle = "Batch-adjusted standardised difference with 95% CI. Positive = the first-named group scores higher.",
       x = "standardised effect (partial Cohen's d)", y = NULL,
       caption = paste(
         "Sex is nested in batch, so no sex main effect is estimable; the strata are shown because pooling cancels opposite-signed effects.",
         "A CI crossing zero is an inconclusive result, not a demonstrated absence of effect.", sep = "\n")) +
  theme_exp9() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())
save_fig(p1, "fig6_effect_sizes", W2, MM(100))

write_xlsx(list(effects = est, primary = primary, sex_interaction = inter),
           file.path(RES, "all_batches_results.xlsx"))

cat("\nfigures ->", FIG, "\nresults  ->", RES, "\n")
