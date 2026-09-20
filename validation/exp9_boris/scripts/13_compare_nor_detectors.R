# ============================================================================
# 13_compare_nor_detectors.R
#
# The two NOR objects are physically different: one round (small, a radius),
# one rectangular (distinct width and height). A single radial detector cannot
# be geometrically correct for both, but shape-matching has its own problem.
#
# The dimensions were defined in the pre-v2 code (commit 8ea9013^):
#
#     object_box_width  = 9      # the NOVEL object got a 9 x 7 cm box
#     object_box_height = 7
#     contact_distance  = 4      # the FAMILIAR object got a 4 cm radius
#
# v2 replaced that with one detector for both, because the two regions have
# different areas (63 vs 50.3 cm2, a factor of 1.25) and the assignment
# follows the novel/familiar ROLE, which rotates between animals. That makes
# the detection area a function of the very quantity D2 measures.
#
# So there are three candidate designs and no purely technical winner:
#
#   A  symmetric      one radius for both. D2 unbiased by construction,
#                     geometrically wrong for at least one object.
#   B  legacy shape   9x7 box on novel, 4 cm circle on familiar. Geometrically
#                     right, but D2 carries a 1.25x area term.
#   C  area-equalised 9x7 box on novel, circle enlarged to the same 63 cm2
#                     (r = 4.48 cm) on familiar. Shape-aware AND area-neutral.
#
# This does not have to be settled by argument. The human scorer watched the
# real objects, so whichever detector best reproduces their scoring is the
# empirically correct one. Batch 1 is the only batch with manual NOR scoring.
#
# Ground truth is the RAW BORIS export, not NOR.xlsx (transposed side headers).
# ============================================================================

suppressMessages({
  library(dplyr)
  library(ggplot2)
})

PROJ  <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
REPO  <- "C:/Users/topohl/Documents/GitHub/SLEAPanalyzer"
BORIS <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B1/NOR/BORIS"
RES   <- file.path(PROJ, "results")
FIG   <- file.path(PROJ, "figures")

source(file.path(REPO, "02_SLEAPanalzyer", "DLCAnalyzer_Functions_final.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "Behavioral_Metrics_Phase1.R"))
for (f in c("io", "events", "interpolation", "geometry", "validation"))
  try(source(file.path(REPO, "02_SLEAPanalzyer", "core", paste0(f, ".R"))), silent = TRUE)

cfg <- yaml::read_yaml(file.path(PROJ, "config", "nor_b1.yaml"))
novelLoc <- read.delim(file.path(PROJ, "metadata", "novelLoc.txt"), stringsAsFactors = FALSE)

BOX_W <- 9; BOX_H <- 7            # novel object, from the pre-v2 definition
CIRC_R <- 4                        # familiar object, from the pre-v2 definition
EQ_R <- sqrt(BOX_W * BOX_H / pi)   # circle with the same area as the box

# --- Ground truth -----------------------------------------------------------
raw <- bind_rows(lapply(list.files(BORIS, pattern = "_nov[.]tsv$", full.names = TRUE), function(p) {
  d <- read.delim(p, stringsAsFactors = FALSE, check.names = FALSE)
  pick <- function(cat) {
    v <- d[["Total duration (s)"]][trimws(d$Category) == cat]
    if (length(v) == 0) NA_real_ else v[1]
  }
  s <- unique(trimws(d$Subject)); s <- s[nzchar(s)]
  data.frame(Code = s[1], manual_left = pick("Interaction left"),
             manual_right = pick("Interaction Right"), stringsAsFactors = FALSE)
}))

# --- Per-frame geometry, once per animal ------------------------------------
files <- list.files(cfg$input_dir, pattern = "[.]csv$", full.names = TRUE)
cat(sprintf("reading %d Batch-1 NOR files\n", length(files)))
geo <- lapply(files, function(p) {
  t <- ReadDLCDataFromCSV(file = p, fps = cfg$fps)
  t <- interpolate_tracking(t, landmarks = c("nose", "bodycentre", "objL", "objR"),
                            max_gap_s = cfg$max_interpolation_gap_s,
                            likelihood_cutoff = cfg$likelihood_cutoff)
  t <- CalibrateTrackingData(t, method = "area",
                             in.metric = cfg$arena_width_cm * cfg$arena_height_cm,
                             points = cfg$arena_corner_names)
  list(code = substr(basename(p), 1, 4), fps = t$fps,
       valid = landmark_validity(t, c("nose", "bodycentre", "objL", "objR")),
       nose = get_point_coordinates(t, "nose"),
       oL = get_point_coordinates(t, "objL"), oR = get_point_coordinates(t, "objR"),
       dL = tracking_point_distance(t, "objL", "nose"),
       dR = tracking_point_distance(t, "objR", "nose"),
       bL = tracking_point_distance(t, "objL", "bodycentre"),
       bR = tracking_point_distance(t, "objR", "bodycentre"),
       aL = tracking_target_angle(t, "objL"), aR = tracking_target_angle(t, "objR"))
})
names(geo) <- vapply(geo, function(x) x$code, character(1))

ANGLE <- cfg$contact_angle_deg
BODY  <- cfg$body_exclusion_distance_cm

# region = TRUE where the nose is inside the detection shape for that side
dur <- function(g, side, region) {
  b <- if (side == "L") g$bL else g$bR
  a <- if (side == "L") g$aL else g$aR
  hit <- region & (b > BODY) & (abs(a) >= ANGLE[1]) & (abs(a) <= ANGLE[2])
  hit[is.na(hit)] <- FALSE
  sum(hit & g$valid, na.rm = TRUE) / g$fps
}
circ <- function(g, side, r) (if (side == "L") g$dL else g$dR) <= r
boxr <- function(g, side, w, h) {
  o <- if (side == "L") g$oL else g$oR
  abs(g$nose$x - o$x) <= w / 2 & abs(g$nose$y - o$y) <= h / 2
}

# Which side holds the NOVEL object. The pipeline's inherited convention --
# metadata "R" means novel is on the LEFT -- is correct for image-consistent
# sides, verified against the raw export in 06_validate_nor_against_raw.R.
novel_side <- function(code) {
  loc <- novelLoc$NovelLoc[novelLoc$Code == code]
  if (length(loc) == 0 || is.na(loc)) return(NA_character_)
  if (loc == "R") "L" else "R"
}

design <- function(name, fun) {
  bind_rows(lapply(names(geo), function(cd) {
    g <- geo[[cd]]; ns <- novel_side(cd)
    if (is.na(ns)) return(NULL)
    fs <- if (ns == "L") "R" else "L"
    data.frame(design = name, Code = cd,
               sleap_nov = fun(g, ns, "novel"), sleap_fam = fun(g, fs, "familiar"),
               stringsAsFactors = FALSE)
  }))
}

designs <- bind_rows(
  design("A symmetric 4 cm",      function(g, s, role) dur(g, s, circ(g, s, 4))),
  design("A symmetric 5 cm",      function(g, s, role) dur(g, s, circ(g, s, 5))),
  design("B legacy shape 9x7 / r4",
         function(g, s, role) if (role == "novel") dur(g, s, boxr(g, s, BOX_W, BOX_H))
                              else dur(g, s, circ(g, s, CIRC_R))),
  design("C area-equalised 9x7 / r4.48",
         function(g, s, role) if (role == "novel") dur(g, s, boxr(g, s, BOX_W, BOX_H))
                              else dur(g, s, circ(g, s, EQ_R)))
)

# --- Score every design against the manual scoring --------------------------
truth <- raw %>%
  rowwise() %>%
  mutate(ns = novel_side(Code),
         manual_nov = if (is.na(ns)) NA_real_ else if (ns == "L") manual_left else manual_right,
         manual_fam = if (is.na(ns)) NA_real_ else if (ns == "L") manual_right else manual_left) %>%
  ungroup() %>%
  mutate(manual_D2 = (manual_nov - manual_fam) / (manual_nov + manual_fam))

ccc <- function(x, y) {
  n <- length(x); vx <- var(x) * (n - 1) / n; vy <- var(y) * (n - 1) / n
  2 * cov(x, y) * (n - 1) / n / (vx + vy + (mean(x) - mean(y))^2)
}

score <- designs %>%
  inner_join(truth %>% select(Code, manual_nov, manual_fam, manual_D2), by = "Code") %>%
  filter(is.finite(manual_nov), is.finite(manual_fam)) %>%
  mutate(sleap_D2 = (sleap_nov - sleap_fam) / (sleap_nov + sleap_fam)) %>%
  group_by(design) %>%
  summarise(
    n = n(),
    nov_ratio = mean(sleap_nov) / mean(manual_nov),
    fam_ratio = mean(sleap_fam) / mean(manual_fam),
    r_contact = cor(c(sleap_nov, sleap_fam), c(manual_nov, manual_fam)),
    ccc_contact = ccc(c(sleap_nov, sleap_fam), c(manual_nov, manual_fam)),
    r_D2 = cor(sleap_D2, manual_D2),
    ccc_D2 = ccc(sleap_D2, manual_D2),
    D2_bias = mean(sleap_D2 - manual_D2),
    .groups = "drop")

write.csv(score, file.path(RES, "nor_detector_comparison.csv"), row.names = FALSE)

cat(sprintf("\nnovel object: %g x %g cm box (%.1f cm2)\n", BOX_W, BOX_H, BOX_W * BOX_H))
cat(sprintf("familiar object: r = %g cm circle (%.1f cm2); area-matched r = %.2f cm\n\n",
            CIRC_R, pi * CIRC_R^2, EQ_R))
cat("=== detector designs scored against the manual scoring (n = 20) ===\n")
print(score %>% mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(),
      row.names = FALSE)

cat("\nnov_ratio / fam_ratio near 1 means the detector recovers the right amount\n")
cat("of contact time for that object. A design that inflates one and not the\n")
cat("other is exactly the D2 bias v2 was written to remove.\n")

best_c <- score$design[which.max(score$ccc_contact)]
best_d <- score$design[which.max(score$ccc_D2)]
cat(sprintf("\nbest contact-time concordance: %s (CCC %.3f)\n",
            best_c, max(score$ccc_contact)))
cat(sprintf("best D2 concordance:           %s (CCC %.3f)\n",
            best_d, max(score$ccc_D2)))

# --- Figure -----------------------------------------------------------------
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
pl <- score %>%
  select(design, `novel object` = nov_ratio, `familiar object` = fam_ratio) %>%
  tidyr::pivot_longer(-design, names_to = "object", values_to = "ratio")
p <- ggplot(pl, aes(ratio, design, colour = object)) +
  geom_vline(xintercept = 1, colour = AXIS, linetype = "22", linewidth = 0.5) +
  geom_point(size = 2.6, position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("#2a78d6", "#eb6834"), name = NULL) +
  labs(title = "Does shape-matching the NOR detector help?",
       subtitle = "SLEAP / manual contact-time ratio per object, Exp9 Batch 1, n = 20. 1.0 is perfect.",
       x = "SLEAP contact time / manual contact time", y = NULL,
       caption = "A design that lands one object on 1.0 and the other off it recovers contact time unevenly, which biases D2.") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA),
        panel.background = element_rect(fill = SURFACE, colour = NA),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(colour = AXIS, linewidth = 0.4),
        axis.text = element_text(colour = MUTED, size = 8),
        axis.title = element_text(colour = INK2, size = 9),
        plot.title = element_text(colour = INK, face = "bold", size = 12),
        plot.subtitle = element_text(colour = INK2, size = 9),
        plot.caption = element_text(colour = MUTED, size = 8, hjust = 0),
        legend.position = "top")
ggsave(file.path(FIG, "fig9_nor_detector_designs.png"), p,
       width = 9, height = 4, dpi = 200, bg = SURFACE)

cat("\nwrote:", file.path(RES, "nor_detector_comparison.csv"), "\n")
