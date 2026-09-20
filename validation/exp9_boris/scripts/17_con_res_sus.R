# ============================================================================
# 17_con_res_sus.R
#
# Characterising CON / RES / SUS across every assay.
#
# The design, stated plainly so the output can be read correctly:
#
#   CON animals are unstressed controls and are NOT classified. SIS animals
#   are split into RES and SUS by thresholding a six-component composite
#   z-scored against the controls: NOR D2, sucrose preference, delta CORT,
#   body-weight development, adrenal weight, spleen weight.
#
# That makes two kinds of outcome, and they support different claims:
#
#   DEFINING    the six components. RES and SUS differ on these BY
#               CONSTRUCTION. Reporting them is right -- they describe what
#               the groups are -- but a difference here is not evidence for
#               anything, it is a restatement of the classification. No
#               multiplicity correction is applied and no p-value is starred.
#
#   INDEPENDENT everything else (EPM, OFT, social preference). These played no
#               part in defining the groups, so a difference here IS evidence
#               about what susceptibility means behaviourally. These carry the
#               inferential weight and are BH-corrected as a family.
#
# CON vs SIS is separate again: Condition is the experimental assignment, so
# it is independent of all of this.
#
# Batch is in every model. It is justified twice over: controls themselves
# differ between batches (F = 4.91, p = 0.0052), and batch explains ~15% of
# the variance in the SLEAP outcomes.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(writexl)
})

PROJ <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
RES  <- file.path(PROJ, "results"); FIG <- file.path(PROJ, "figures")

dat <- read.delim(file.path(PROJ, "enriched", "sleap_all_batches_wide.tsv"),
                  stringsAsFactors = FALSE, na.strings = "")

# Which SLEAP outcomes correspond to a classifier component.
METRICS <- tibble::tribble(
  ~label,                   ~col,                        ~kind,
  "NOR D2",                 "sleap_NOR_D2",              "defining",
  "NOR interaction time",   "sleap_NOR_total_s",         "shares measurement",
  "EPM open-arm fraction",  "sleap_EPM_open_frac",       "independent",
  "EPM centre time",        "sleap_EPM_center_s",        "independent",
  "EPM nose dips",          "sleap_EPM_nosedip_n",       "independent",
  "EPM distance",           "sleap_EPM_distance_cm",     "independent",
  "OFT centre %",           "sleap_OFT_center_pct",      "independent",
  "OFT corner %",           "sleap_OFT_corner_pct",      "independent",
  "OFT distance",           "sleap_OFT_distance_cm",     "independent",
  "SocP S1 preference",     "sleap_SocP_S1_pref_index",  "independent",
  "SocP S2 preference",     "sleap_SocP_S2_pref_index",  "independent"
)

fit_group <- function(metric, phenotype_col) {
  d <- dat %>% filter(!is.na(.data[[metric]]), !is.na(.data[[phenotype_col]]))
  d$g <- factor(d[[phenotype_col]], levels = c("CON", "RES", "SUS"))
  if (dplyr::n_distinct(d$g) < 3) return(NULL)
  fit <- lm(reformulate(c("g", "Batch"), metric), data = d)
  s <- summary(fit)$sigma
  om <- anova(fit)["g", "Pr(>F)"]

  # Pairwise contrasts from the same batch-adjusted fit, so every comparison
  # shares one error term and one batch adjustment.
  em <- c(CON = 0, RES = unname(coef(fit)["gRES"]), SUS = unname(coef(fit)["gSUS"]))
  V <- vcov(fit)
  pair <- function(a, b) {
    cf <- rep(0, length(coef(fit))); names(cf) <- names(coef(fit))
    if (a != "CON") cf[paste0("g", a)] <- -1
    if (b != "CON") cf[paste0("g", b)] <- 1
    est <- sum(cf * coef(fit)); se <- sqrt(drop(t(cf) %*% V %*% cf))
    tv <- est / se; df <- fit$df.residual
    data.frame(contrast = paste(b, "-", a), d = est / s,
               lo = (est - qt(.975, df) * se) / s, hi = (est + qt(.975, df) * se) / s,
               p = 2 * pt(-abs(tv), df))
  }
  # Captured before the pipe: inside mutate(), `d` is the effect-size column.
  sizes <- c(n = nrow(d), CON = sum(d$g == "CON"),
             RES = sum(d$g == "RES"), SUS = sum(d$g == "SUS"))
  bind_rows(pair("CON", "RES"), pair("CON", "SUS"), pair("RES", "SUS")) %>%
    mutate(metric = metric, n = sizes[["n"]], omnibus_p = om,
           n_CON = sizes[["CON"]], n_RES = sizes[["RES"]], n_SUS = sizes[["SUS"]])
}

run_all <- function(phenotype_col, tag) {
  bind_rows(lapply(METRICS$col, fit_group, phenotype_col = phenotype_col)) %>%
    left_join(METRICS, by = c("metric" = "col")) %>%
    mutate(scheme = tag) %>%
    # BH across the INDEPENDENT metrics only, within each contrast. The
    # defining metrics are not inferential, so including them in the family
    # would both dilute the correction and imply they are being tested.
    group_by(scheme, contrast) %>%
    mutate(p_BH = ifelse(kind == "independent",
                         p.adjust(ifelse(kind == "independent", p, NA), "BH"),
                         NA_real_)) %>%
    ungroup()
}

res <- bind_rows(
  run_all("Phenotype", "your list (sus_animals.txt)"),
  run_all("Phenotype_batchCorrected", "batch z-scored list")
)
write.csv(res, file.path(RES, "con_res_sus_comparisons.csv"), row.names = FALSE)

show <- function(sc, kd, ttl) {
  cat(sprintf("\n=== %s  [%s] ===\n", ttl, sc))
  x <- res %>% filter(scheme == sc, kind == kd) %>%
    mutate(d = round(d, 2), CI = sprintf("[%.2f, %.2f]", lo, hi),
           p = signif(p, 3), q = ifelse(is.na(p_BH), "-", signif(p_BH, 3))) %>%
    select(label, contrast, n_CON, n_RES, n_SUS, d, CI, p, q) %>%
    arrange(label, contrast)
  print(as.data.frame(x), row.names = FALSE)
}

PRIMARY <- "your list (sus_animals.txt)"
cat("group sizes and design:\n")
cat(sprintf("  CON %d | RES %d | SUS %d\n",
            sum(dat$Phenotype == "CON", na.rm = TRUE),
            sum(dat$Phenotype == "RES", na.rm = TRUE),
            sum(dat$Phenotype == "SUS", na.rm = TRUE)))

show(PRIMARY, "independent", "INDEPENDENT outcomes -- these carry evidential weight (BH-corrected)")
show(PRIMARY, "defining", "DEFINING outcome -- differs by construction, descriptive only")
show(PRIMARY, "shares measurement", "SHARES MEASUREMENT with a defining variable -- read with care")

sig <- res %>% filter(scheme == PRIMARY, kind == "independent", p_BH < 0.05)
cat(sprintf("\nindependent outcomes surviving BH within contrast: %d of %d\n",
            nrow(sig), sum(res$scheme == PRIMARY & res$kind == "independent")))
if (nrow(sig)) print(sig %>% mutate(d = round(d, 2), q = signif(p_BH, 3)) %>%
                       select(label, contrast, d, q) %>% as.data.frame(), row.names = FALSE)

cat("\n=== does the conclusion depend on which list is used? ===\n")
cmp <- res %>% filter(kind == "independent") %>%
  select(scheme, label, contrast, d, p) %>%
  pivot_wider(names_from = scheme, values_from = c(d, p))
names(cmp) <- gsub("your list \\(sus_animals.txt\\)", "yours", names(cmp))
names(cmp) <- gsub("batch z-scored list", "zscored", names(cmp))
print(cmp %>% mutate(across(starts_with("d_"), ~round(., 2)),
                     across(starts_with("p_"), ~signif(., 3))) %>%
        arrange(`p_yours`) %>% head(8) %>% as.data.frame(), row.names = FALSE)

# --- Figure -----------------------------------------------------------------
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
pl <- res %>% filter(scheme == PRIMARY) %>%
  mutate(kind = factor(kind, levels = c("independent", "shares measurement", "defining")),
         label = factor(label, levels = rev(METRICS$label)),
         contrast = factor(contrast, levels = c("RES - CON", "SUS - CON", "SUS - RES")))
p <- ggplot(pl, aes(d, label, colour = kind)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.5) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.5) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(independent = "#2a78d6",
                                 `shares measurement` = "#eda100",
                                 defining = "#e34948"), name = NULL) +
  facet_wrap(~ contrast) +
  labs(title = "CON, RES and SUS across every assay",
       subtitle = "Batch-adjusted standardised difference with 95% CI. Positive = the first-named group scores higher.",
       x = "standardised effect (partial Cohen's d)", y = NULL,
       caption = paste(
         "Red: a component of the classifier -- RES and SUS differ on it by construction, so it describes the groups rather than testing them.",
         "Amber: a different measure from the same recordings as a component. Blue: played no part in the classification, so these carry the evidence.",
         sep = "\n")) +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA),
        panel.background = element_rect(fill = SURFACE, colour = NA),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = AXIS, linewidth = 0.4),
        axis.text = element_text(colour = MUTED, size = 8),
        axis.title = element_text(colour = INK2, size = 9),
        plot.title = element_text(colour = INK, face = "bold", size = 12),
        plot.subtitle = element_text(colour = INK2, size = 9),
        plot.caption = element_text(colour = MUTED, size = 8, hjust = 0),
        strip.text = element_text(colour = INK, face = "bold", size = 9),
        legend.position = "top")
ggsave(file.path(FIG, "fig11_con_res_sus.png"), p,
       width = 11, height = 6, dpi = 200, bg = SURFACE)

write_xlsx(list(comparisons = res), file.path(RES, "con_res_sus.xlsx"))
cat("\nwrote:", file.path(RES, "con_res_sus_comparisons.csv"), "\n")
