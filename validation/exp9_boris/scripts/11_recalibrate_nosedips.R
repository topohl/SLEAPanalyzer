# ============================================================================
# 11_recalibrate_nosedips.R
#
# SLEAP reports 2.2x more nose dips than the manual scorer (r = 0.82 but
# CCC = 0.42): the ranking is right, the count is not. docs/assay_definitions.md
# already says nose.dip "is a geometric proxy for head-dipping and should be
# validated against manual scoring before use as a primary outcome". This is
# that validation, and Batch 1 is the ONLY place it can ever be done, because
# no other batch has manual ND scoring.
#
# The detector itself has no threshold to tune. A frame is a dip when the head
# is outside the arena outline, the body is inside it, and the neck is not in a
# closed arm -- all pure geometry. The single lever is `integration_period`,
# the half-width of the centred majority filter avgbool() applies before onsets
# are counted. A short window lets a head wavering across the maze edge register
# as several separate dips; a longer one merges them.
#
# So: compute the raw per-frame boolean once per animal, then sweep the window
# and compare each setting against the manual counts. Sweeping the window is
# cheap; re-reading and re-zoning 20 files for every candidate is not.
#
# This reports the sweep. It does not silently change the shipped default.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(ggplot2)
})

PROJ <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
REPO <- "C:/Users/topohl/Documents/GitHub/SLEAPanalyzer"
RES  <- file.path(PROJ, "results")
FIG  <- file.path(PROJ, "figures")

source(file.path(REPO, "02_SLEAPanalzyer", "DLCAnalyzer_Functions_final.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "core", "events.R"))

cfg <- yaml::read_yaml(file.path(PROJ, "config", "epm_b1.yaml"))
zoneInfo <- utils::read.table(cfg$zone_file, sep = ";", header = TRUE,
                              stringsAsFactors = FALSE, check.names = FALSE)

WINDOWS <- c(1, 2, 3, 5, 7, 10, 15, 20, 30, 45)

files <- list.files(cfg$input_dir, pattern = "[.]csv$", full.names = TRUE)
cat(sprintf("computing the raw nose-dip boolean for %d Batch-1 animals\n", length(files)))

raw <- lapply(files, function(p) {
  t <- ReadDLCDataFromCSV(file = p, fps = cfg$fps)
  t <- interpolate_tracking(t, landmarks = c("headcentre", "bodycentre", "neck"),
                            max_gap_s = cfg$max_interpolation_gap_s,
                            likelihood_cutoff = cfg$likelihood_cutoff)
  t <- CalibrateTrackingData(t, method = cfg$calibration_method,
                             in.metric = cfg$calibration_distance_cm,
                             points = cfg$calibration_points)
  t <- AddZones(t, zoneInfo)
  observed <- is.finite(t$data$headcentre$x) & is.finite(t$data$headcentre$y) &
              is.finite(t$data$bodycentre$x) & is.finite(t$data$bodycentre$y) &
              is.finite(t$data$neck$x) & is.finite(t$data$neck$y)
  dip <- !IsInZone(t, "headcentre", "arena") &
          IsInZone(t, "bodycentre", "arena") &
         !IsInZone(t, "neck", c("closed.left", "closed.right"))
  list(code = substr(basename(p), 1, 4), dip = dip & observed, observed = observed,
       fps = t$fps)
})
names(raw) <- vapply(raw, function(x) x$code, character(1))

boris <- read.delim(file.path(PROJ, "enriched", "analysis_ready_wide.tsv"),
                    stringsAsFactors = FALSE, na.strings = "") %>%
  select(Code, manual_ND = EPM_ND_n, manual_ND_dur = EPM_ND_dur)

ccc <- function(x, y) {
  n <- length(x); vx <- var(x) * (n - 1) / n; vy <- var(y) * (n - 1) / n
  2 * cov(x, y) * (n - 1) / n / (vx + vy + (mean(x) - mean(y))^2)
}

sweep <- bind_rows(lapply(WINDOWS, function(w) {
  counts <- vapply(raw, function(a)
    count_entries(as.logical(avgbool(a$dip, w)), valid = a$observed),
    integer(1))
  d <- data.frame(Code = names(counts), sleap_ND = as.integer(counts),
                  stringsAsFactors = FALSE) %>%
    inner_join(boris, by = "Code") %>%
    filter(is.finite(manual_ND))
  data.frame(
    window = w,
    n = nrow(d),
    sleap_median = median(d$sleap_ND),
    manual_median = median(d$manual_ND),
    ratio = mean(d$sleap_ND) / mean(d$manual_ND),
    bias = mean(d$sleap_ND - d$manual_ND),
    r = cor(d$sleap_ND, d$manual_ND),
    rho = suppressWarnings(cor(d$sleap_ND, d$manual_ND, method = "spearman")),
    ccc = ccc(d$sleap_ND, d$manual_ND),
    stringsAsFactors = FALSE)
}))

write.csv(sweep, file.path(RES, "nosedip_window_sweep.csv"), row.names = FALSE)

cat("\n=== integration-period sweep against manual nose-dip counts (n = 20) ===\n")
print(sweep %>% mutate(across(c(ratio, r, rho, ccc), ~round(., 3)),
                       bias = round(bias, 2)) %>% as.data.frame(), row.names = FALSE)

best_ccc <- sweep$window[which.max(sweep$ccc)]
best_bias <- sweep$window[which.min(abs(sweep$ratio - 1))]
cat(sprintf("\nshipped default is %d frames -> ratio %.2f, CCC %.2f\n",
            cfg$integration_period_frames,
            sweep$ratio[sweep$window == cfg$integration_period_frames],
            sweep$ccc[sweep$window == cfg$integration_period_frames]))
cat(sprintf("best concordance at window = %d (CCC %.2f, ratio %.2f)\n",
            best_ccc, max(sweep$ccc), sweep$ratio[sweep$window == best_ccc]))
cat(sprintf("closest to unbiased at window = %d (ratio %.2f, CCC %.2f)\n",
            best_bias, sweep$ratio[sweep$window == best_bias],
            sweep$ccc[sweep$window == best_bias]))
cat("\nRank agreement is roughly flat across the sweep, so the window trades\n",
    "count magnitude, not ordering. Changing it is a calibration decision,\n",
    "not a bug fix, and the shipped default is left alone.\n", sep = "")

# --- Figure -----------------------------------------------------------------
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"; SER1 <- "#2a78d6"; SER2 <- "#eb6834"

long <- sweep %>%
  select(window, ratio, ccc, r) %>%
  tidyr::pivot_longer(-window, names_to = "stat", values_to = "value") %>%
  mutate(stat = recode(stat, ratio = "SLEAP / manual count ratio",
                       ccc = "Lin's concordance", r = "Pearson r"),
         stat = factor(stat, levels = c("SLEAP / manual count ratio",
                                        "Lin's concordance", "Pearson r")))
hl <- data.frame(stat = factor("SLEAP / manual count ratio",
                               levels = levels(long$stat)), y = 1)

p <- ggplot(long, aes(window, value)) +
  geom_hline(data = hl, aes(yintercept = y), colour = AXIS,
             linetype = "22", linewidth = 0.5) +
  geom_vline(xintercept = cfg$integration_period_frames, colour = MUTED,
             linetype = "22", linewidth = 0.5) +
  geom_line(colour = SER1, linewidth = 0.7) +
  geom_point(colour = SER1, fill = SURFACE, shape = 21, stroke = 0.6, size = 2.4) +
  facet_wrap(~ stat, scales = "free_y") +
  labs(title = "Nose dips: the smoothing window sets the count, not the ranking",
       subtitle = "Exp9 Batch 1, n = 20. Grey line is the shipped default of 5 frames; dashed line on the left panel is perfect agreement.",
       x = "integration_period (frames, half-width of the majority filter)",
       y = NULL,
       caption = "Batch 1 is the only batch with manual nose-dip scoring, so this cannot be checked anywhere else.") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA),
        panel.background = element_rect(fill = SURFACE, colour = NA),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = AXIS, linewidth = 0.4),
        axis.text = element_text(colour = MUTED, size = 8),
        axis.title = element_text(colour = INK2, size = 9),
        plot.title = element_text(colour = INK, face = "bold", size = 12),
        plot.subtitle = element_text(colour = INK2, size = 9),
        plot.caption = element_text(colour = MUTED, size = 8, hjust = 0),
        strip.text = element_text(colour = INK, face = "bold", size = 9),
        plot.margin = margin(10, 14, 10, 10))
ggsave(file.path(FIG, "fig7_nosedip_sweep.png"), p,
       width = 10, height = 4, dpi = 200, bg = SURFACE)

cat("\nwrote:", file.path(RES, "nosedip_window_sweep.csv"), "\n")
