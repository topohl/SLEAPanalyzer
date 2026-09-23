# ============================================================================
# 05_correlate.R
#
# Three questions, kept separate because they need different statistics:
#
#   1. METHOD AGREEMENT. Does SLEAPanalyzer measure what the human scorer
#      measured? Correlation alone cannot answer this -- a metric can correlate
#      perfectly and still be biased -- so each matched pair also gets Lin's
#      concordance coefficient and Bland-Altman bias with limits of agreement.
#
#   2. CROSS-ASSAY STRUCTURE. How do the behavioural measures relate to each
#      other? Spearman, because n = 20 and several metrics are skewed, with
#      Benjamini-Hochberg correction across the reported matrix.
#
#   3. GROUP DIFFERENCES, reported for Condition and the current canonical
#      phenotype. Deprecated batch-corrected calls are neither read nor used.
#
# n = 20 (4 CON / 16 SIS). Every correlation here has a 95% CI roughly +/- 0.4
# wide, so these are estimates with wide uncertainty, not confirmations.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(writexl)
})

# Resolve paths from this script's own location, so the study runs from any
# checkout of the repository.
SCRIPTS <- local({
  a <- commandArgs(trailingOnly = FALSE)
  p <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/")
})
PROJ   <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")
ENRICH <- file.path(PROJ, "enriched")
RUN_ROOT <- Sys.getenv("EXP9_SLEAP_RUN_ROOT", unset = PROJ)
FIG <- Sys.getenv("EXP9_SLEAP_FIGURES_DIR", unset = file.path(RUN_ROOT, "figures"))
RES <- Sys.getenv("EXP9_SLEAP_RESULTS_DIR", unset = file.path(RUN_ROOT, "results"))
SOURCE_DATA <- Sys.getenv(
  "EXP9_SLEAP_SOURCE_DATA_DIR",
  unset = file.path(RUN_ROOT, "source_data", "validation")
)
EXP9 <- Sys.getenv(
  "EXP9_ROOT",
  unset = "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress"
)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
dir.create(SOURCE_DATA, showWarnings = FALSE, recursive = TRUE)

# Shared Nature-style theme, identical to the publication figure set.
source(file.path(SCRIPTS, "00_theme.R"))

boris <- read.delim(file.path(ENRICH, "analysis_ready_wide.tsv"),
                    stringsAsFactors = FALSE, na.strings = "")
forbidden_phenotype_columns <- c(
  "Phenotype_batchCorrected", "Phenotype_bc_complement", "Phenotype_conflict"
)
sleap <- read.delim(file.path(ENRICH, "sleap_wide.tsv"),
                    stringsAsFactors = FALSE, na.strings = "") %>%
  select(-any_of(forbidden_phenotype_columns))
dat <- inner_join(sleap, boris %>% select(-c(ID, Batch, Sex, Condition, Phenotype,
                                             Phenotype_batchCorrected,
                                             Phenotype_bc_complement,
                                             Phenotype_conflict)),
                  by = "Code")
stopifnot(nrow(dat) == 20)

canon <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(grepl("^[0-9]+$", x), sub("^0+(?=[0-9])", "", x, perl = TRUE), x)
}
read_list <- function(path) {
  values <- trimws(readLines(path, warn = FALSE))
  canon(values[nzchar(values)])
}
con_list <- read_list(file.path(EXP9, "Analysis/con_animals.csv"))
sus_list <- read_list(file.path(EXP9, "Analysis/sus_animals.txt"))
dat <- dat %>%
  mutate(
    key = canon(ID),
    Condition = ifelse(key %in% con_list, "CON", "SIS"),
    Phenotype = case_when(
      Condition == "CON" ~ "CON",
      key %in% sus_list ~ "SUS",
      TRUE ~ "RES"
    ),
    Phenotype_source = "Analysis/{con_animals.csv, sus_animals.txt}, RES by complement"
  ) %>%
  select(-key)
write.table(
  dat,
  file.path(SOURCE_DATA, "method_validation_matched_data.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = ""
)


# --- Agreement statistics ---------------------------------------------------
# Lin's concordance correlation coefficient: Pearson r penalised by any shift
# in mean or scale, so it separates "tracks the same ranking" from "reports
# the same number".
ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 3) return(NA_real_)
  vx <- var(x) * (n - 1) / n; vy <- var(y) * (n - 1) / n
  cxy <- cov(x, y) * (n - 1) / n
  2 * cxy / (vx + vy + (mean(x) - mean(y))^2)
}

agreement <- function(sleap_v, boris_v, label, assay, unit, same_quantity = TRUE) {
  ok <- is.finite(sleap_v) & is.finite(boris_v)
  s <- sleap_v[ok]; b <- boris_v[ok]
  pe <- suppressWarnings(cor.test(s, b, method = "pearson"))
  sp <- suppressWarnings(cor.test(s, b, method = "spearman", exact = FALSE))
  d  <- s - b
  data.frame(
    assay = assay, metric = label, unit = unit, same_quantity = same_quantity,
    n = length(s),
    pearson_r = unname(pe$estimate),
    ci_low = pe$conf.int[1], ci_high = pe$conf.int[2],
    p_pearson = pe$p.value,
    spearman_rho = unname(sp$estimate), p_spearman = sp$p.value,
    ccc = ccc(s, b),
    bias_sleap_minus_boris = mean(d),
    bias_pct = 100 * mean(d) / mean(b),
    loa_low = mean(d) - 1.96 * sd(d), loa_high = mean(d) + 1.96 * sd(d),
    sleap_median = median(s), boris_median = median(b),
    stringsAsFactors = FALSE
  )
}

# same_quantity = FALSE marks a pair that is related but is NOT the same
# measurement, so its bias and CCC must not be read as method error. BORIS's
# grouped "Total number" for OpenTime counts scored behaviour bouts in the open
# arms (OpenExplor + OpenSAP + OpenGroom + OpenND), whereas SLEAP's
# open.entries counts zone crossings. There is no manual zone-entry count in
# this dataset, so the two cannot be reconciled; the correlation is reported
# because it is informative, the agreement statistics are not.
pairs_spec <- list(
  list("EPM", "Open-arm time",     "s",     "sleap_EPM_open_s",      "EPM_OpenTime_dur", TRUE),
  list("EPM", "Closed-arm time",   "s",     "sleap_EPM_closed_s",    "EPM_ClosedTime_dur", TRUE),
  list("EPM", "Centre time",       "s",     "sleap_EPM_center_s",    "EPM_CenterTime_dur", TRUE),
  list("EPM", "Open entries vs bouts", "count", "sleap_EPM_open_entries","EPM_OpenTime_n", FALSE),
  list("EPM", "Nose dips",         "count", "sleap_EPM_nosedip_n",   "EPM_ND_n", TRUE),
  list("NOR", "Novel object",      "s",     "sleap_NOR_nov_s",       "NOR_nov", TRUE),
  list("NOR", "Familiar object",   "s",     "sleap_NOR_fam_s",       "NOR_fam", TRUE),
  list("NOR", "Total interaction", "s",     "sleap_NOR_total_s",     "NOR_total_interaction", TRUE),
  list("NOR", "Discrimination D2", "index", "sleap_NOR_D2",          "NOR_D2", TRUE),
  list("SocP", "S1 novel",         "s",     "sleap_SocP_S1_novel",    "SocP_S1_novel", TRUE),
  list("SocP", "S1 familiar",      "s",     "sleap_SocP_S1_familiar", "SocP_S1_familiar", TRUE),
  list("SocP", "S1 preference",    "index", "sleap_SocP_S1_pref_index","SocP_S1_pref_index", TRUE),
  list("SocP", "S2 novel",         "s",     "sleap_SocP_S2_novel",    "SocP_S2_novel", TRUE),
  list("SocP", "S2 familiar",      "s",     "sleap_SocP_S2_familiar", "SocP_S2_familiar", TRUE),
  list("SocP", "S2 preference",    "index", "sleap_SocP_S2_pref_index","SocP_S2_pref_index", TRUE)
)

agree <- do.call(rbind, lapply(pairs_spec, function(p) {
  agreement(dat[[p[[4]]]], dat[[p[[5]]]], p[[2]], p[[1]], p[[3]], p[[6]])
}))
agree$p_pearson_BH <- p.adjust(agree$p_pearson, method = "BH")

write.csv(agree, file.path(RES, "method_agreement.csv"), row.names = FALSE)

cat("=== 1. METHOD AGREEMENT: SLEAP vs manual BORIS ===\n")
print(agree %>%
  mutate(across(c(pearson_r, ci_low, ci_high, spearman_rho, ccc), ~round(., 3)),
         across(c(bias_sleap_minus_boris, bias_pct), ~round(., 2)),
         p_BH = signif(p_pearson_BH, 3)) %>%
  select(assay, metric, unit, same_qty = same_quantity, n,
         pearson_r, ci_low, ci_high, p_BH,
         spearman_rho, ccc, bias = bias_sleap_minus_boris, bias_pct),
  row.names = FALSE)
cat("\nsame_qty = FALSE means the two columns are related but are not the same\n",
    "measurement, so CCC and bias there are not method error.\n", sep = "")

# --- Figure 1: agreement scatter, small multiples ---------------------------
long <- do.call(rbind, lapply(pairs_spec, function(p) {
  data.frame(assay = p[[1]],
             metric = factor(p[[2]], levels = vapply(pairs_spec, `[[`, character(1), 2)),
             unit = p[[3]], Code = dat$Code,
             sleap = dat[[p[[4]]]], boris = dat[[p[[5]]]],
             stringsAsFactors = FALSE)
}))
lab <- agree %>%
  mutate(txt = sprintf("r = %.2f\nCCC = %.2f", pearson_r, ccc),
         metric = factor(metric, levels = levels(long$metric)))

p1 <- ggplot(long, aes(boris, sleap)) +
  # Identity line first, so it reads as the reference the points are judged
  # against rather than as a fitted trend.
  geom_abline(slope = 1, intercept = 0, colour = RULE, linetype = "22",
              linewidth = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = SER2, linewidth = 0.5) +
  geom_point(colour = SER1, size = 1.5, alpha = 0.7, stroke = 0, shape = 16) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt),
            hjust = -0.14, vjust = 1.2, size = 1.9, colour = INK,
            lineheight = 1.05, inherit.aes = FALSE) +
  facet_wrap(~ metric, scales = "free", ncol = 5) +
  labs(
    title = "SLEAPanalyzer v2 against manual BORIS scoring",
    subtitle = paste("Exp9 cohort 1, n = 20. Dashed grey is identity;",
                     "<span style='color:#F4636E'>**coral**</span> is the least-squares fit."),
    x = "Manual BORIS", y = "SLEAP",
    caption = paste(
      "CCC is Lin's concordance: it penalises bias, so a high r with a low CCC means the same ranking at a different value.",
      "NOR novel/familiar use contactNov/contactFam, which compensate for the left/right mirror between the two sources.",
      sep = "\n")
  ) +
  theme_exp9(grid = "both")
save_fig(p1, "fig1_method_agreement", W2, MM(105))

# --- Figure 2: Bland-Altman for the same-unit metrics -----------------------
ba <- long %>%
  filter(unit == "s") %>%
  mutate(mean_val = (sleap + boris) / 2, diff_val = sleap - boris)
ba_lines <- ba %>%
  group_by(metric) %>%
  summarise(bias = mean(diff_val),
            lo = mean(diff_val) - 1.96 * sd(diff_val),
            hi = mean(diff_val) + 1.96 * sd(diff_val), .groups = "drop")

p2 <- ggplot(ba, aes(mean_val, diff_val)) +
  geom_hline(yintercept = 0, colour = RULE, linewidth = 0.35) +
  geom_hline(data = ba_lines, aes(yintercept = bias), colour = SER2, linewidth = 0.5) +
  geom_hline(data = ba_lines, aes(yintercept = lo), colour = RULE,
             linetype = "22", linewidth = 0.35) +
  geom_hline(data = ba_lines, aes(yintercept = hi), colour = RULE,
             linetype = "22", linewidth = 0.35) +
  geom_point(colour = SER1, size = 1.5, alpha = 0.7, stroke = 0, shape = 16) +
  facet_wrap(~ metric, scales = "free", ncol = 5) +
  labs(title = "Bland-Altman: where the two methods disagree",
       subtitle = paste("Seconds-valued metrics only.",
                        "<span style='color:#F4636E'>**Coral**</span> is the mean difference;",
                        "dashed are the 95% limits of agreement."),
       x = "Mean of the two methods (s)", y = "SLEAP - BORIS (s)",
       caption = "A sloped cloud means the disagreement grows with the value; an offset orange line means constant bias.") +
  theme_exp9(grid = "both")
save_fig(p2, "fig2_bland_altman", W2, MM(82))

# --- 2. Cross-assay correlation matrix --------------------------------------
# Manual measures only, so the structure reported is not an artefact of the
# tracking pipeline.
matrix_vars <- c(
  "EPM open time"      = "EPM_OpenTime_dur",
  "EPM centre time"    = "EPM_CenterTime_dur",
  "EPM open entries"   = "EPM_OpenTime_n",
  "EPM SAP"            = "EPM_SAP_dur",
  "EPM nose dips"      = "EPM_ND_n",
  "EPM grooming"       = "EPM_groom_dur",
  "NOR D2"             = "NOR_D2",
  "NOR interaction"    = "NOR_total_interaction",
  "SocP S1 preference" = "SocP_S1_pref_index",
  "SocP S2 preference" = "SocP_S2_pref_index",
  "SocP S1 total"      = "SocP_S1_total"
)
M <- dat[, matrix_vars]; colnames(M) <- names(matrix_vars)

cm <- expand.grid(x = names(matrix_vars), y = names(matrix_vars),
                  stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    rho = suppressWarnings(cor(M[[x]], M[[y]], method = "spearman",
                               use = "pairwise.complete.obs")),
    p   = if (x == y) NA_real_ else suppressWarnings(
            cor.test(M[[x]], M[[y]], method = "spearman", exact = FALSE)$p.value)
  ) %>%
  ungroup()
# Correct over the unique off-diagonal pairs only, not the mirrored matrix.
upper <- cm %>%
  filter(match(x, names(matrix_vars)) < match(y, names(matrix_vars))) %>%
  mutate(p_BH = p.adjust(p, method = "BH"))
cm <- cm %>%
  left_join(upper %>% select(x, y, p_BH), by = c("x", "y")) %>%
  left_join(upper %>% select(x = y, y = x, p_BH2 = p_BH), by = c("x", "y")) %>%
  mutate(p_BH = coalesce(p_BH, p_BH2)) %>%
  select(-p_BH2)

write.csv(upper %>% arrange(p_BH), file.path(RES, "cross_assay_spearman.csv"),
          row.names = FALSE)

cat("\n=== 2. CROSS-ASSAY correlations (manual measures, Spearman, BH-corrected) ===\n")
cat("strongest 12 of", nrow(upper), "unique pairs:\n")
print(upper %>% arrange(desc(abs(rho))) %>% head(12) %>%
  mutate(rho = round(rho, 3), p = signif(p, 3), p_BH = signif(p_BH, 3)) %>%
  as.data.frame(), row.names = FALSE)
cat(sprintf("\npairs surviving BH < 0.05: %d of %d\n",
            sum(upper$p_BH < 0.05, na.rm = TRUE), nrow(upper)))

lvl <- names(matrix_vars)
p3 <- ggplot(cm %>% filter(x != y) %>%
               mutate(x = factor(x, lvl), y = factor(y, rev(lvl))),
             aes(x, y, fill = rho)) +
  # The diagonal is dropped: rho = 1 by construction carries no information,
  # and rendering it at full saturation pulls the eye to the one part of the
  # matrix that cannot say anything.
  geom_tile(colour = "white", linewidth = 0.9) +
  geom_text(aes(label = ifelse(is.na(p_BH) | p_BH >= 0.05, "",
                               sprintf("%.2f", rho))),
            size = 1.9, colour = INK) +
  scale_fill_gradient2(low = DIV_NEG, mid = DIV_MID, high = DIV_POS,
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman\nrho") +
  coord_fixed() +
  labs(title = "Cross-assay correlations among the manual measures",
       subtitle = "Exp9 Batch 1, n = 20. Values shown only where Benjamini-Hochberg q < 0.05.",
       x = NULL, y = NULL,
       caption = paste("Blank cells are not null results: at n = 20 only very large<br>",
                       "effects clear correction. Every estimate is in<br>",
                       "results/cross_assay_spearman.csv.")) +
  theme_exp9(grid = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right", legend.key.width = unit(5, "pt"),
        legend.key.height = unit(22, "pt"))
save_fig(p3, "fig3_cross_assay_matrix", W15, MM(112))

# --- 3. Group comparisons using the canonical phenotype ---------------------
group_vars <- c("EPM open time" = "EPM_OpenTime_dur",
                "EPM centre time" = "EPM_CenterTime_dur",
                "EPM SAP" = "EPM_SAP_dur",
                "NOR D2" = "NOR_D2",
                "SocP S1 preference" = "SocP_S1_pref_index",
                "SocP S2 preference" = "SocP_S2_pref_index")

group_long <- bind_rows(lapply(names(group_vars), function(nm) {
  bind_rows(
    data.frame(metric = nm, value = dat[[group_vars[[nm]]]],
               grouping = "Condition", group = dat$Condition,
               Code = dat$Code, stringsAsFactors = FALSE),
    data.frame(metric = nm, value = dat[[group_vars[[nm]]]],
               grouping = "Phenotype (canonical)", group = dat$Phenotype,
               Code = dat$Code, stringsAsFactors = FALSE)
  )
})) %>%
  mutate(group = ifelse(is.na(group), "unknown", group),
         metric = factor(metric, levels = names(group_vars)),
         grouping = factor(grouping, levels = c("Condition", "Phenotype (canonical)")),
         group = factor(group, levels = c("CON", "SIS", "RES", "SUS", "unknown")))

group_stats <- group_long %>%
  group_by(grouping, metric, group) %>%
  summarise(n = sum(is.finite(value)), median = median(value, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
            .groups = "drop")
write.csv(group_stats, file.path(RES, "group_summaries.csv"), row.names = FALSE)

# Kruskal-Wallis across the known groups; "unknown" is never a group.
group_tests <- group_long %>%
  filter(group != "unknown") %>%
  group_by(grouping, metric) %>%
  filter(n_distinct(group) > 1) %>%
  summarise(
    groups = paste(sprintf("%s(n=%d)", levels(droplevels(group)),
                           as.integer(table(droplevels(group)))), collapse = " "),
    p = suppressWarnings(kruskal.test(value ~ droplevels(group))$p.value),
    .groups = "drop") %>%
  group_by(grouping) %>%
  mutate(p_BH = p.adjust(p, method = "BH")) %>%
  ungroup()
write.csv(group_tests, file.path(RES, "group_tests.csv"), row.names = FALSE)

cat("\n=== 3. GROUP DIFFERENCES using the canonical phenotype ===\n")
print(group_tests %>% mutate(p = signif(p, 3), p_BH = signif(p_BH, 3)) %>%
        as.data.frame(), row.names = FALSE)

# Identity comes from x position, not colour, so the aqua contrast caveat in
# the palette reference does not apply.
p4 <- ggplot(group_long %>% filter(group != "unknown") %>% mutate(group = droplevels(group)),
             aes(group, value, colour = group)) +
  # A bare median line, not a crossbar: the box of a crossbar reads as an
  # interval it does not represent.
  geom_jitter(width = 0.18, height = 0, size = 1.4, alpha = 0.62,
              stroke = 0, shape = 16) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "errorbar", width = 0.55, linewidth = 0.6,
               colour = INK, show.legend = FALSE) +
  scale_colour_manual(values = c(CON = PAL[["CON"]], SIS = PAL[["SIS"]],
                                 RES = PAL[["RES"]], SUS = PAL[["SUS"]]),
                      guide = "none") +
  # scales = "free" in facet_grid frees x per column and y per row, so each
  # grouping shows only its own groups while a metric keeps one y scale
  # across both groupings.
  facet_grid(metric ~ grouping, scales = "free", switch = "y") +
  labs(title = "Group differences use the canonical phenotype",
       subtitle = paste0("Bar is the group median. ",
                         colour_key(c("CON", "RES", "SUS")),
                         "<br>RES is the complement among SIS animals after applying ",
                         "the current<br>canonical susceptible roster."),
       x = NULL, y = NULL,
       caption = paste0("Groups are identified by position, not colour alone.<br>",
                        "With 4 CON animals these comparisons are severely underpowered; ",
                        "see statistics/group_tests.csv.")) +
  theme_exp9() +
  theme(strip.placement = "outside",
        strip.text.y.left = element_markdown(angle = 0, hjust = 1, size = 6,
                                             face = "plain"),
        strip.text.x = element_markdown(hjust = 0.5, size = 7.5))
save_fig(p4, "fig4_group_differences", W15, MM(175))

# --- Per-animal discrepancies -----------------------------------------------
# Which animal x metric combinations disagree most between the two methods,
# scaled by that metric's own spread so the assays are comparable. These are
# the recordings to re-check against video.
discrep <- long %>%
  filter(unit %in% c("s", "index")) %>%
  group_by(metric) %>%
  mutate(diff = sleap - boris,
         z_diff = (diff - mean(diff)) / sd(diff)) %>%
  ungroup() %>%
  arrange(desc(abs(z_diff))) %>%
  select(assay, metric, Code, boris, sleap, diff, z_diff)

write.csv(discrep, file.path(RES, "per_animal_discrepancies.csv"), row.names = FALSE)

cat("\n=== 4. LARGEST SLEAP-vs-BORIS DISAGREEMENTS (re-check these recordings) ===\n")
print(discrep %>% filter(abs(z_diff) > 2.2) %>%
        mutate(across(c(boris, sleap, diff, z_diff), ~round(., 3))) %>%
        as.data.frame(), row.names = FALSE)

flagged <- discrep %>% filter(abs(z_diff) > 2.2) %>% count(Code, sort = TRUE)
cat("\nanimals flagged on more than one metric:\n")
print(as.data.frame(flagged %>% filter(n > 1)), row.names = FALSE)

write_xlsx(
  list(method_agreement = agree,
       cross_assay = upper %>% arrange(p_BH),
       group_summaries = group_stats,
       group_tests = group_tests,
       discrepancies = discrep),
  file.path(RES, "correlation_results.xlsx")
)

cat("\nfigures ->", FIG, "\nresults  ->", RES, "\n")
