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
#   3. GROUP DIFFERENCES, reported under BOTH phenotype definitions, because
#      the curated and batch-corrected calls disagree for T2H7 and the curated
#      call is missing for five animals.
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

PROJ   <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
ENRICH <- file.path(PROJ, "enriched")
FIG    <- file.path(PROJ, "figures")
RES    <- file.path(PROJ, "results")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

boris <- read.delim(file.path(ENRICH, "analysis_ready_wide.tsv"),
                    stringsAsFactors = FALSE, na.strings = "")
sleap <- read.delim(file.path(ENRICH, "sleap_wide.tsv"),
                    stringsAsFactors = FALSE, na.strings = "")
dat <- inner_join(sleap, boris %>% select(-c(ID, Batch, Sex, Condition, Phenotype,
                                             Phenotype_batchCorrected,
                                             Phenotype_bc_complement,
                                             Phenotype_conflict)),
                  by = "Code")
stopifnot(nrow(dat) == 20)

# --- Theme ------------------------------------------------------------------
# Values from the data-viz reference palette, used unchanged.
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
SER1 <- "#2a78d6"; SER2 <- "#eb6834"; SER3 <- "#1baf7a"
DIV_NEG <- "#2a78d6"; DIV_MID <- "#f0efec"; DIV_POS <- "#d03b3b"

theme_viz <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background   = element_rect(fill = SURFACE, colour = NA),
      panel.background  = element_rect(fill = SURFACE, colour = NA),
      panel.grid.major  = element_line(colour = GRID, linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(colour = AXIS, linewidth = 0.4),
      axis.text         = element_text(colour = MUTED, size = base_size - 2),
      axis.title        = element_text(colour = INK2, size = base_size - 1),
      plot.title        = element_text(colour = INK, face = "bold", size = base_size + 2),
      plot.subtitle     = element_text(colour = INK2, size = base_size - 1),
      plot.caption      = element_text(colour = MUTED, size = base_size - 2, hjust = 0),
      strip.text        = element_text(colour = INK, face = "bold", size = base_size - 1),
      strip.background  = element_blank(),
      legend.text       = element_text(colour = INK2, size = base_size - 2),
      legend.title      = element_text(colour = INK2, size = base_size - 1),
      plot.margin       = margin(10, 14, 10, 10)
    )
}

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
  geom_abline(slope = 1, intercept = 0, colour = AXIS, linetype = "22", linewidth = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = SER2, linewidth = 0.7) +
  geom_point(colour = SER1, fill = SURFACE, shape = 21, stroke = 0.6, size = 2.4) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt),
            hjust = -0.12, vjust = 1.15, size = 2.7, colour = INK2,
            lineheight = 0.95, inherit.aes = FALSE) +
  facet_wrap(~ metric, scales = "free", ncol = 5) +
  labs(
    title = "SLEAPanalyzer v2 against manual BORIS scoring",
    subtitle = "Exp9 Batch 1, n = 20. Dashed line is identity; orange is the least-squares fit.",
    x = "Manual BORIS", y = "SLEAP",
    caption = paste(
      "CCC is Lin's concordance: it penalises bias, so a high r with a low CCC means the same ranking at a different value.",
      "NOR novel/familiar use contactNov/contactFam, which compensate for the left/right mirror between the two sources.",
      sep = "\n")
  ) +
  theme_viz()
ggsave(file.path(FIG, "fig1_method_agreement.png"), p1,
       width = 13, height = 7.2, dpi = 200, bg = SURFACE)
ggsave(file.path(FIG, "fig1_method_agreement.pdf"), p1, width = 13, height = 7.2)

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
  geom_hline(yintercept = 0, colour = AXIS, linewidth = 0.4) +
  geom_hline(data = ba_lines, aes(yintercept = bias), colour = SER2, linewidth = 0.6) +
  geom_hline(data = ba_lines, aes(yintercept = lo), colour = MUTED,
             linetype = "22", linewidth = 0.4) +
  geom_hline(data = ba_lines, aes(yintercept = hi), colour = MUTED,
             linetype = "22", linewidth = 0.4) +
  geom_point(colour = SER1, fill = SURFACE, shape = 21, stroke = 0.6, size = 2.4) +
  facet_wrap(~ metric, scales = "free", ncol = 5) +
  labs(title = "Bland-Altman: where the two methods disagree",
       subtitle = "Seconds-valued metrics only. Orange is mean difference, dashed are 95% limits of agreement.",
       x = "Mean of the two methods (s)", y = "SLEAP - BORIS (s)",
       caption = "A sloped cloud means the disagreement grows with the value; an offset orange line means constant bias.") +
  theme_viz()
ggsave(file.path(FIG, "fig2_bland_altman.png"), p2,
       width = 13, height = 5.6, dpi = 200, bg = SURFACE)

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
  geom_tile(colour = SURFACE, linewidth = 1.2) +
  geom_text(aes(label = ifelse(is.na(p_BH) | p_BH >= 0.05, "",
                               sprintf("%.2f", rho))),
            size = 2.6, colour = INK) +
  scale_fill_gradient2(low = DIV_NEG, mid = DIV_MID, high = DIV_POS,
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman\nrho") +
  coord_fixed() +
  labs(title = "Cross-assay correlations among the manual measures",
       subtitle = "Exp9 Batch 1, n = 20. Values shown only where Benjamini-Hochberg q < 0.05.",
       x = NULL, y = NULL,
       caption = "Blank cells are not null results: at n = 20 only very large effects clear correction. See results/cross_assay_spearman.csv for every estimate.") +
  theme_viz() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), axis.line = element_blank())
ggsave(file.path(FIG, "fig3_cross_assay_matrix.png"), p3,
       width = 8.4, height = 7.6, dpi = 200, bg = SURFACE)

# --- 3. Group comparisons under both phenotype definitions ------------------
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
               grouping = "Phenotype (curated)", group = dat$Phenotype,
               Code = dat$Code, stringsAsFactors = FALSE),
    data.frame(metric = nm, value = dat[[group_vars[[nm]]]],
               grouping = "Phenotype (batch-corrected)",
               group = dat$Phenotype_batchCorrected,
               Code = dat$Code, stringsAsFactors = FALSE)
  )
})) %>%
  mutate(group = ifelse(is.na(group), "unknown", group),
         metric = factor(metric, levels = names(group_vars)),
         grouping = factor(grouping, levels = c("Condition", "Phenotype (curated)",
                                                "Phenotype (batch-corrected)")),
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

cat("\n=== 3. GROUP DIFFERENCES under both phenotype definitions ===\n")
print(group_tests %>% mutate(p = signif(p, 3), p_BH = signif(p_BH, 3)) %>%
        as.data.frame(), row.names = FALSE)

# Identity comes from x position, not colour, so the aqua contrast caveat in
# the palette reference does not apply.
p4 <- ggplot(group_long %>% filter(group != "unknown") %>% mutate(group = droplevels(group)),
             aes(group, value, colour = group)) +
  # A bare median line, not a crossbar: the box of a crossbar reads as an
  # interval it does not represent.
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "errorbar", width = 0.5, linewidth = 0.6, colour = INK2) +
  geom_jitter(width = 0.14, height = 0, size = 2.2, alpha = 0.9) +
  scale_colour_manual(values = c(CON = SER1, SIS = SER2, RES = SER3, SUS = "#4a3aa7"),
                      guide = "none") +
  # scales = "free" in facet_grid frees x per column and y per row, so each
  # grouping shows only its own groups while a metric keeps one y scale
  # across all three groupings.
  facet_grid(metric ~ grouping, scales = "free", switch = "y") +
  labs(title = "Group differences are reported under both phenotype definitions",
       subtitle = "Bar is the group median. The curated and batch-corrected calls disagree for T2H7, and the curated call is missing for 5 animals.",
       x = NULL, y = NULL,
       caption = "Groups are identified by position, not colour alone. With 4 CON animals these comparisons are severely underpowered; see results/group_tests.csv.") +
  theme_viz() +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1))
ggsave(file.path(FIG, "fig4_group_differences.png"), p4,
       width = 9.5, height = 11, dpi = 200, bg = SURFACE)

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
