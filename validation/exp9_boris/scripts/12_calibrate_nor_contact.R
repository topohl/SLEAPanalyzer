# ============================================================================
# 12_calibrate_nor_contact.R
#
# contact_distance_cm and contact_angle_deg are the last unvalidated parameters
# in the NOR pipeline. docs/assay_definitions.md says so explicitly:
# "REQUIRES VALIDATION. These thresholds must be checked against manually
# scored video for your object size and camera before results are published."
#
# The scaled analysis used the shipped example values (4 cm, [70, 290]). Prior
# work in Raw Data/.../NOR/SLEAP/output_angle/ swept 5 cm, 6 cm, 70-270 and
# 80-280, so the choice was already suspected to matter.
#
# Batch 1 is the only batch with manual NOR scoring, so this is the only place
# the detector can be calibrated. Ground truth is the RAW BORIS export
# (Interaction left / Interaction Right), not NOR.xlsx, whose side headers are
# transposed -- see 06_validate_nor_against_raw.R.
#
# The expensive part is reading, interpolating and calibrating each file. The
# per-frame nose-object distances, body-object distances and orientation angles
# are computed ONCE per animal; sweeping thresholds over them is cheap.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(ggplot2)
})

# Resolve paths from this script's own location, so the study runs from any
# checkout of the repository.
SCRIPTS <- local({
  a <- commandArgs(trailingOnly = FALSE)
  p <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/")
})
PROJ <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")

# Shared Nature-style theme, identical to the publication figure set.
source(file.path(SCRIPTS, "00_theme.R"))
REPO  <- "C:/Users/topohl/Documents/GitHub/SLEAPanalyzer"
BORIS <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B1/NOR/BORIS"
RES   <- file.path(PROJ, "results")
FIG   <- file.path(PROJ, "figures")

source(file.path(REPO, "02_SLEAPanalzyer", "DLCAnalyzer_Functions_final.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "Behavioral_Metrics_Phase1.R"))
for (f in c("io", "events", "interpolation", "geometry", "validation"))
  try(source(file.path(REPO, "02_SLEAPanalzyer", "core", paste0(f, ".R"))), silent = TRUE)

cfg <- yaml::read_yaml(file.path(PROJ, "config", "nor_b1.yaml"))

# --- Ground truth: the raw per-animal BORIS exports -------------------------
# Keyed on the Subject column, never the file name: three exports are misnamed
# with lookalike characters (F2l1, S7l3, R803).
raw <- bind_rows(lapply(list.files(BORIS, pattern = "_nov[.]tsv$", full.names = TRUE), function(p) {
  d <- read.delim(p, stringsAsFactors = FALSE, check.names = FALSE)
  pick <- function(cat) {
    v <- d[["Total duration (s)"]][trimws(d$Category) == cat]
    if (length(v) == 0) NA_real_ else v[1]
  }
  subject <- unique(trimws(d$Subject)); subject <- subject[nzchar(subject)]
  data.frame(Code = subject[1],
             manual_left = pick("Interaction left"),
             manual_right = pick("Interaction Right"),
             stringsAsFactors = FALSE)
}))
cat(sprintf("manual ground truth: %d animals\n", nrow(raw)))

# --- Per-frame geometry, computed once per animal ---------------------------
files <- list.files(cfg$input_dir, pattern = "[.]csv$", full.names = TRUE)
cat(sprintf("computing per-frame geometry for %d tracking files\n", length(files)))

geo <- lapply(files, function(p) {
  t <- ReadDLCDataFromCSV(file = p, fps = cfg$fps)
  t <- interpolate_tracking(t, landmarks = c("nose", "bodycentre", "objL", "objR"),
                            max_gap_s = cfg$max_interpolation_gap_s,
                            likelihood_cutoff = cfg$likelihood_cutoff)
  t <- CalibrateTrackingData(t, method = "area",
                             in.metric = cfg$arena_width_cm * cfg$arena_height_cm,
                             points = cfg$arena_corner_names)
  valid <- landmark_validity(t, c("nose", "bodycentre", "objL", "objR"))
  list(code = substr(basename(p), 1, 4), fps = t$fps, valid = valid,
       dL = tracking_point_distance(t, "objL", "nose"),
       dR = tracking_point_distance(t, "objR", "nose"),
       bL = tracking_point_distance(t, "objL", "bodycentre"),
       bR = tracking_point_distance(t, "objR", "bodycentre"),
       aL = tracking_target_angle(t, "objL"),
       aR = tracking_target_angle(t, "objR"))
})
names(geo) <- vapply(geo, function(x) x$code, character(1))

# Contact exactly as Behavioral_Metrics_Phase1.R assembles it for "radial":
# nose inside the radius, body beyond the exclusion distance, head oriented.
duration <- function(g, side, radius, angle, body_excl) {
  d <- if (side == "L") g$dL else g$dR
  b <- if (side == "L") g$bL else g$bR
  a <- if (side == "L") g$aL else g$aR
  hit <- (d <= radius) & (b > body_excl) & (abs(a) >= angle[1]) & (abs(a) <= angle[2])
  hit[is.na(hit)] <- FALSE
  sum(hit & g$valid, na.rm = TRUE) / g$fps
}

ccc <- function(x, y) {
  n <- length(x); vx <- var(x) * (n - 1) / n; vy <- var(y) * (n - 1) / n
  2 * cov(x, y) * (n - 1) / n / (vx + vy + (mean(x) - mean(y))^2)
}

RADII  <- c(2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8)
ANGLES <- list("70-290 (shipped)" = c(70, 290), "80-280" = c(80, 280),
               "90-270" = c(90, 270), "none (0-360)" = c(0, 360))

sweep <- bind_rows(lapply(names(ANGLES), function(an) {
  bind_rows(lapply(RADII, function(r) {
    est <- bind_rows(lapply(names(geo), function(cd) {
      g <- geo[[cd]]
      data.frame(Code = cd,
                 sleap_left  = duration(g, "L", r, ANGLES[[an]], cfg$body_exclusion_distance_cm),
                 sleap_right = duration(g, "R", r, ANGLES[[an]], cfg$body_exclusion_distance_cm),
                 stringsAsFactors = FALSE)
    }))
    j <- inner_join(est, raw, by = "Code") %>%
      filter(is.finite(manual_left), is.finite(manual_right))
    # Both sides stacked: the detector is symmetric, so one agreement figure.
    s <- c(j$sleap_left, j$sleap_right); b <- c(j$manual_left, j$manual_right)
    tot_s <- j$sleap_left + j$sleap_right; tot_b <- j$manual_left + j$manual_right
    data.frame(angle = an, radius_cm = r, n = nrow(j),
               ratio = mean(s) / mean(b), bias = mean(s - b),
               r = cor(s, b), ccc = ccc(s, b),
               r_total = cor(tot_s, tot_b), ccc_total = ccc(tot_s, tot_b),
               stringsAsFactors = FALSE)
  }))
}))

write.csv(sweep, file.path(RES, "nor_contact_sweep.csv"), row.names = FALSE)

cat("\n=== per-side agreement against the raw BORIS export (n = 20 animals, 40 sides) ===\n")
print(sweep %>%
        mutate(across(c(ratio, r, ccc, r_total, ccc_total), ~round(., 3)),
               bias = round(bias, 2)) %>%
        select(angle, radius_cm, ratio, bias, r, ccc) %>%
        as.data.frame(), row.names = FALSE)

best <- sweep %>% slice_max(ccc, n = 1)
ship <- sweep %>% filter(angle == "70-290 (shipped)", radius_cm == cfg$contact_distance_cm)
cat(sprintf("\nshipped setting (%.1f cm, %s): ratio %.2f, CCC %.3f, r %.3f\n",
            cfg$contact_distance_cm, "70-290", ship$ratio, ship$ccc, ship$r))
cat(sprintf("best concordance   (%.1f cm, %s): ratio %.2f, CCC %.3f, r %.3f\n",
            best$radius_cm, best$angle, best$ratio, best$ccc, best$r))
unb <- sweep %>% slice_min(abs(ratio - 1), n = 1)
cat(sprintf("closest to unbiased (%.1f cm, %s): ratio %.2f, CCC %.3f\n",
            unb$radius_cm, unb$angle, unb$ratio, unb$ccc))

# --- Figure -----------------------------------------------------------------
p <- sweep %>%
  tidyr::pivot_longer(c(ratio, ccc), names_to = "stat", values_to = "value") %>%
  mutate(stat = recode(stat, ratio = "SLEAP / manual contact time",
                       ccc = "Lin's concordance"),
         angle = factor(angle, levels = names(ANGLES))) %>%
  ggplot(aes(radius_cm, value, colour = angle)) +
  geom_hline(data = data.frame(stat = "SLEAP / manual contact time", y = 1),
             aes(yintercept = y), colour = AXIS, linetype = "22", linewidth = 0.5) +
  geom_vline(xintercept = cfg$contact_distance_cm, colour = MUTED,
             linetype = "22", linewidth = 0.5) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_colour_manual(values = c(SER1, SER2, SER3, SER4), name = NULL) +
  facet_wrap(~ stat, scales = "free_y") +
  labs(title = "NOR contact detector, calibrated against manual scoring",
       subtitle = "Exp9 Batch 1, 20 animals / 40 object sides. Grey line is the shipped 4 cm; dashed is perfect agreement.",
       x = "contact_distance_cm (radial)", y = NULL,
       caption = "Ground truth is the raw BORIS export, not NOR.xlsx. Batch 1 is the only batch with manual NOR scoring.") +
  theme_exp9() +
  theme(legend.position = "top")
save_fig(p, "fig8_nor_contact_sweep", W2, MM(81))

cat("\nwrote:", file.path(RES, "nor_contact_sweep.csv"), "\n")
