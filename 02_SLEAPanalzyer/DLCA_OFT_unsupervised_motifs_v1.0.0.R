# ================================================================
# SLEAP/DLC OFT Unsupervised Behavioral Motif Discovery
# Author: Tobias Pohl / ChatGPT revision
# Version: 1.0.0
# Date: 2026-05-11
# ================================================================
#
# Run after DLCA_OFT v1.2.0.R.
#
# Purpose:
# Discover recurring OFT behavioral motifs from short time windows without
# using Group, Treatment, Sex, or stress labels. Biological interpretation is
# assigned only after motif discovery from motif feature profiles.
#
# Pipeline:
# 1. Load frame-level OFT outputs.
# 2. Split each animal trajectory into short windows.
# 3. Extract low-level kinematic/spatial descriptors per window.
# 4. Scale features and perform PCA.
# 5. Cluster windows using Gaussian mixture models if mclust is available,
#    otherwise k-means with silhouette-based k selection.
# 6. Summarise motif occupancy, transitions, entropy, and representative windows.
# 7. Generate post hoc motif annotations from quantitative profiles.
#
# Important:
# The discovery step is unsupervised. Motif labels are post hoc biological
# annotations, analogous to annotating unsupervised scRNA-seq clusters.
# ================================================================

required_packages <- c(
  "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble", "ggplot2", "cluster"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

optional_packages <- c("mclust")
for (pkg in optional_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Optional package not installed: ", pkg, ". Falling back to k-means motif discovery.")
  }
}

config <- list(
  behavior_root = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior",
  batch_output_subdir = file.path("OFT", "SLEAP", "output_v1.2.0"),
  global_output_folder = "OFT_SLEAP_output_v1.2.0",
  batches = c("B1"),
  fps = 30,

  # Windowing.
  window_size_s = 5,
  window_step_s = 5,
  min_valid_fraction_per_window = 0.70,

  # Unsupervised model.
  pca_variance_explained_target = 0.85,
  min_pca_components = 2,
  max_pca_components = 10,
  k_min = 2,
  k_max = 8,
  kmeans_nstart = 100,
  kmeans_iter_max = 100,
  random_seed = 123,

  # Biological post hoc annotation thresholds.
  high_quantile = 0.67,
  low_quantile = 0.33,

  # Plots.
  save_plots = TRUE,
  plot_width_single_col = 85 / 25.4,
  plot_height_single_col = 65 / 25.4,
  plot_width_double_col = 180 / 25.4,
  plot_height_double_col = 100 / 25.4
)

safe_divide <- function(x, y) {
  ifelse(is.na(y) | y == 0, NA_real_, x / y)
}

entropy <- function(x) {
  tab <- table(x, useNA = "no")
  if (length(tab) == 0) return(NA_real_)
  p <- as.numeric(tab) / sum(tab)
  -sum(p * log2(p))
}

entry_count <- function(x) {
  x <- tidyr::replace_na(as.logical(x), FALSE)
  sum(x & !dplyr::lag(x, default = FALSE), na.rm = TRUE)
}

make_theme_nature <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.text = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0)
    )
}

get_frame_files <- function(config) {
  purrr::map(config$batches, function(batch) {
    frame_dir <- file.path(config$behavior_root, batch, config$batch_output_subdir, "frames")
    if (!dir.exists(frame_dir)) return(character())
    list.files(frame_dir, pattern = "_frames\\.csv$", full.names = TRUE)
  }) %>% unlist()
}

assign_zone_state <- function(df) {
  df %>%
    mutate(
      zone_state = case_when(
        in_center ~ "center",
        in_corners ~ "corner",
        in_periphery ~ "periphery",
        TRUE ~ "other"
      )
    )
}

make_window_features_one_file <- function(frame_file, config) {
  df <- readr::read_csv(frame_file, show_col_types = FALSE) %>%
    arrange(frame) %>%
    assign_zone_state() %>%
    mutate(
      dx = x - dplyr::lag(x),
      dy = y - dplyr::lag(y),
      step_cm = sqrt(dx^2 + dy^2),
      heading_rad = atan2(dy, dx),
      turn_angle_rad = atan2(sin(heading_rad - dplyr::lag(heading_rad)), cos(heading_rad - dplyr::lag(heading_rad))),
      turn_angle_deg = abs(turn_angle_rad) * 180 / pi,
      acceleration_abs_cm_s2 = abs(acceleration_cm_s2),
      window_id = floor(time_s / config$window_step_s),
      window_start_s = window_id * config$window_step_s,
      window_end_s = window_start_s + config$window_size_s
    )

  df %>%
    filter(time_s >= window_start_s, time_s < window_end_s) %>%
    group_by(file, Batch, Code, window_id, window_start_s, window_end_s) %>%
    summarise(
      n_frames = n(),
      valid_fraction = mean(!is.na(x) & !is.na(y) & !is.na(speed_cm_s)),
      mean_speed_cm_s = mean(speed_cm_s, na.rm = TRUE),
      median_speed_cm_s = median(speed_cm_s, na.rm = TRUE),
      sd_speed_cm_s = sd(speed_cm_s, na.rm = TRUE),
      max_speed_cm_s = max(speed_cm_s, na.rm = TRUE),
      p95_speed_cm_s = quantile(speed_cm_s, 0.95, na.rm = TRUE, names = FALSE),
      distance_cm = sum(speed_cm_s / config$fps, na.rm = TRUE),
      mean_acceleration_abs_cm_s2 = mean(acceleration_abs_cm_s2, na.rm = TRUE),
      sd_acceleration_abs_cm_s2 = sd(acceleration_abs_cm_s2, na.rm = TRUE),
      mean_turn_angle_deg = mean(turn_angle_deg, na.rm = TRUE),
      sd_turn_angle_deg = sd(turn_angle_deg, na.rm = TRUE),
      high_turn_fraction = mean(turn_angle_deg >= 90, na.rm = TRUE),
      x_sd = sd(x, na.rm = TRUE),
      y_sd = sd(y, na.rm = TRUE),
      spatial_spread_cm = sqrt(x_sd^2 + y_sd^2),
      mean_wall_distance_cm = mean(wall_distance_cm, na.rm = TRUE),
      sd_wall_distance_cm = sd(wall_distance_cm, na.rm = TRUE),
      min_wall_distance_cm = min(wall_distance_cm, na.rm = TRUE),
      center_fraction = mean(in_center, na.rm = TRUE),
      periphery_fraction = mean(in_periphery, na.rm = TRUE),
      corner_fraction = mean(in_corners, na.rm = TRUE),
      immobile_fraction = mean(is_immobile, na.rm = TRUE),
      moving_fraction = mean(is_moving, na.rm = TRUE),
      center_entries = entry_count(in_center),
      zone_entropy = entropy(zone_state),
      .groups = "drop"
    ) %>%
    filter(valid_fraction >= config$min_valid_fraction_per_window) %>%
    mutate(
      window_uid = paste(file, window_id, sep = "__")
    )
}

select_motif_features <- function(window_tbl) {
  candidate_features <- c(
    "mean_speed_cm_s", "median_speed_cm_s", "sd_speed_cm_s", "max_speed_cm_s",
    "p95_speed_cm_s", "distance_cm", "mean_acceleration_abs_cm_s2",
    "sd_acceleration_abs_cm_s2", "mean_turn_angle_deg", "sd_turn_angle_deg",
    "high_turn_fraction", "spatial_spread_cm", "mean_wall_distance_cm",
    "sd_wall_distance_cm", "min_wall_distance_cm", "center_fraction",
    "periphery_fraction", "corner_fraction", "immobile_fraction", "moving_fraction",
    "center_entries", "zone_entropy"
  )

  present <- intersect(candidate_features, names(window_tbl))
  present[vapply(window_tbl[present], is.numeric, logical(1))]
}

prepare_scaled_matrix <- function(window_tbl, features) {
  mat <- window_tbl %>%
    select(all_of(features)) %>%
    mutate(across(everything(), ~ ifelse(is.infinite(.x), NA_real_, .x)))

  # Median imputation is unsupervised and only uses feature distributions.
  mat <- mat %>%
    mutate(across(everything(), ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)))

  keep <- vapply(mat, function(x) sd(x, na.rm = TRUE) > 0, logical(1))
  mat <- mat[, keep, drop = FALSE]

  scaled <- scale(as.matrix(mat))
  list(matrix = scaled, features = colnames(mat))
}

choose_pca_components <- function(pca, config) {
  var_exp <- pca$sdev^2 / sum(pca$sdev^2)
  cum <- cumsum(var_exp)
  n <- which(cum >= config$pca_variance_explained_target)[1]
  n <- max(config$min_pca_components, n, na.rm = TRUE)
  n <- min(config$max_pca_components, n, ncol(pca$x), na.rm = TRUE)
  n
}

cluster_windows <- function(score_mat, config) {
  set.seed(config$random_seed)

  if (requireNamespace("mclust", quietly = TRUE)) {
    fit <- mclust::Mclust(score_mat, G = config$k_min:config$k_max, verbose = FALSE)
    return(list(
      method = "Gaussian mixture model (mclust)",
      cluster = as.integer(fit$classification),
      k = length(unique(fit$classification)),
      model = fit$modelName,
      bic = fit$bic
    ))
  }

  k_values <- config$k_min:min(config$k_max, nrow(score_mat) - 1)
  sil_tbl <- purrr::map_dfr(k_values, function(k) {
    km <- stats::kmeans(score_mat, centers = k, nstart = config$kmeans_nstart, iter.max = config$kmeans_iter_max)
    sil <- cluster::silhouette(km$cluster, stats::dist(score_mat))
    tibble(k = k, mean_silhouette = mean(sil[, "sil_width"]))
  })

  best_k <- sil_tbl$k[which.max(sil_tbl$mean_silhouette)]
  km <- stats::kmeans(score_mat, centers = best_k, nstart = config$kmeans_nstart, iter.max = config$kmeans_iter_max)

  list(
    method = "k-means on PCA scores",
    cluster = as.integer(km$cluster),
    k = best_k,
    model = NA_character_,
    silhouette = sil_tbl
  )
}

annotate_motifs <- function(profile_tbl, config) {
  q_hi <- function(x) quantile(x, config$high_quantile, na.rm = TRUE, names = FALSE)
  q_lo <- function(x) quantile(x, config$low_quantile, na.rm = TRUE, names = FALSE)

  profile_tbl %>%
    mutate(
      high_speed = mean_speed_cm_s >= q_hi(mean_speed_cm_s),
      low_speed = mean_speed_cm_s <= q_lo(mean_speed_cm_s),
      high_immobility = immobile_fraction >= q_hi(immobile_fraction),
      high_center = center_fraction >= q_hi(center_fraction),
      high_wall_distance = mean_wall_distance_cm >= q_hi(mean_wall_distance_cm),
      low_wall_distance = mean_wall_distance_cm <= q_lo(mean_wall_distance_cm),
      high_corner = corner_fraction >= q_hi(corner_fraction),
      high_turning = mean_turn_angle_deg >= q_hi(mean_turn_angle_deg),
      high_burstiness = sd_speed_cm_s >= q_hi(sd_speed_cm_s) | p95_speed_cm_s >= q_hi(p95_speed_cm_s),
      high_spread = spatial_spread_cm >= q_hi(spatial_spread_cm),
      motif_annotation = case_when(
        high_immobility & (low_wall_distance | high_corner) ~ "wall/corner immobility",
        high_speed & low_wall_distance & !high_center ~ "thigmotactic locomotion",
        high_center & high_wall_distance & high_spread ~ "center exploration",
        high_burstiness & high_turning ~ "fragmented burst exploration",
        low_speed & !high_immobility & low_wall_distance ~ "cautious wall-local exploration",
        low_speed & !high_immobility & high_wall_distance ~ "slow centerward scanning",
        TRUE ~ "mixed exploratory motif"
      ),
      annotation_basis = paste0(
        "speed=", ifelse(high_speed, "high", ifelse(low_speed, "low", "mid")),
        "; immobility=", ifelse(high_immobility, "high", "not-high"),
        "; center=", ifelse(high_center, "high", "not-high"),
        "; wall_distance=", ifelse(high_wall_distance, "high", ifelse(low_wall_distance, "low", "mid")),
        "; turning=", ifelse(high_turning, "high", "not-high"),
        "; burstiness=", ifelse(high_burstiness, "high", "not-high")
      )
    )
}

make_animal_motif_summary <- function(window_scores) {
  motif_counts <- window_scores %>%
    count(file, Batch, Code, motif, name = "n_windows") %>%
    group_by(file, Batch, Code) %>%
    mutate(
      motif_fraction = n_windows / sum(n_windows),
      motif_entropy = entropy(motif)
    ) %>%
    ungroup()

  wide <- motif_counts %>%
    select(file, Batch, Code, motif, motif_fraction) %>%
    pivot_wider(
      names_from = motif,
      values_from = motif_fraction,
      names_prefix = "motif_fraction_",
      values_fill = 0
    )

  entropy_tbl <- motif_counts %>%
    group_by(file, Batch, Code) %>%
    summarise(
      motif_entropy = entropy(rep(motif, n_windows)),
      dominant_motif = motif[which.max(n_windows)],
      dominant_motif_fraction = max(motif_fraction),
      .groups = "drop"
    )

  wide %>% left_join(entropy_tbl, by = c("file", "Batch", "Code"))
}

make_motif_transition_summary <- function(window_scores) {
  window_scores %>%
    arrange(file, window_id) %>%
    group_by(file, Batch, Code) %>%
    mutate(next_motif = lead(motif)) %>%
    ungroup() %>%
    filter(!is.na(next_motif)) %>%
    count(file, Batch, Code, motif, next_motif, name = "n") %>%
    group_by(file, Batch, Code, motif) %>%
    mutate(transition_probability = n / sum(n)) %>%
    ungroup()
}

save_plots <- function(window_scores, loadings_tbl, motif_profiles, output_dir, config) {
  if (!isTRUE(config$save_plots)) return(invisible(NULL))
  plot_dir <- file.path(output_dir, "unsupervised_motif_plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  if (all(c("PC1", "PC2") %in% names(window_scores))) {
    p <- ggplot(window_scores, aes(x = PC1, y = PC2, factor(motif))) +
      geom_point(size = 0.7, alpha = 0.55) +
      labs(title = "Unsupervised OFT motif space", x = "PC1", y = "PC2") +
      make_theme_nature()
    ggsave(file.path(plot_dir, "oft_unsupervised_motif_pca_space.svg"), p,
           width = config$plot_width_single_col, height = config$plot_height_single_col)
  }

  top_loadings <- loadings_tbl %>%
    filter(PC %in% c("PC1", "PC2", "PC3")) %>%
    group_by(PC) %>%
    slice_max(abs_loading, n = 8, with_ties = FALSE) %>%
    ungroup()

  if (nrow(top_loadings) > 0) {
    p <- ggplot(top_loadings, aes(x = reorder(feature, loading), y = loading)) +
      geom_col(width = 0.75) +
      coord_flip() +
      facet_wrap(~ PC, scales = "free_y") +
      labs(title = "Top PCA loadings for unsupervised OFT motifs", x = NULL, y = "Loading") +
      make_theme_nature()
    ggsave(file.path(plot_dir, "oft_unsupervised_pca_top_loadings.svg"), p,
           width = config$plot_width_double_col, height = config$plot_height_double_col)
  }

  profile_long <- motif_profiles %>%
    select(motif, motif_annotation, mean_speed_cm_s, center_fraction, immobile_fraction,
           mean_wall_distance_cm, mean_turn_angle_deg, spatial_spread_cm, sd_speed_cm_s) %>%
    pivot_longer(-c(motif, motif_annotation), names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    mutate(z_value = as.numeric(scale(value))) %>%
    ungroup()

  p <- ggplot(profile_long, aes(x = feature, y = factor(motif), fill = z_value)) +
    geom_tile() +
    coord_flip() +
    labs(title = "Post hoc motif feature profiles", x = NULL, y = "Motif") +
    make_theme_nature()
  ggsave(file.path(plot_dir, "oft_unsupervised_motif_profile_heatmap.svg"), p,
         width = config$plot_width_double_col, height = config$plot_height_single_col)

  invisible(NULL)
}

# -------------------------------
# Main execution
# -------------------------------

set.seed(config$random_seed)

frame_files <- get_frame_files(config)
if (length(frame_files) == 0) stop("No frame-level OFT files found. Run DLCA_OFT v1.2.0.R first.")

output_dir <- file.path(config$behavior_root, config$global_output_folder, "unsupervised_motifs_v1.0.0")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Creating time-window feature matrix...")
window_features <- purrr::map_dfr(frame_files, make_window_features_one_file, config = config)
if (nrow(window_features) < 10) stop("Too few valid windows for unsupervised motif discovery.")

features <- select_motif_features(window_features)
scaled <- prepare_scaled_matrix(window_features, features)

message("Running PCA...")
pca <- stats::prcomp(scaled$matrix, center = FALSE, scale. = FALSE)
n_pc <- choose_pca_components(pca, config)
pc_scores <- as_tibble(pca$x[, seq_len(n_pc), drop = FALSE])
colnames(pc_scores) <- paste0("PC", seq_len(n_pc))

var_exp <- pca$sdev^2 / sum(pca$sdev^2)
pca_variance <- tibble(
  PC = paste0("PC", seq_along(var_exp)),
  variance_explained = var_exp,
  cumulative_variance_explained = cumsum(var_exp)
)

loadings_tbl <- as.data.frame(pca$rotation[, seq_len(n_pc), drop = FALSE]) %>%
  rownames_to_column("feature") %>%
  as_tibble() %>%
  pivot_longer(-feature, names_to = "PC", values_to = "loading") %>%
  mutate(abs_loading = abs(loading))

message("Discovering unsupervised motifs...")
clustering <- cluster_windows(as.matrix(pc_scores), config)

window_scores <- bind_cols(window_features, pc_scores) %>%
  mutate(
    motif = paste0("M", clustering$cluster),
    motif_numeric = clustering$cluster,
    clustering_method = clustering$method,
    n_motifs = clustering$k
  )

motif_profiles <- window_scores %>%
  group_by(motif) %>%
  summarise(
    n_windows = n(),
    n_animals = n_distinct(file),
    mean_speed_cm_s = mean(mean_speed_cm_s, na.rm = TRUE),
    sd_speed_cm_s = mean(sd_speed_cm_s, na.rm = TRUE),
    p95_speed_cm_s = mean(p95_speed_cm_s, na.rm = TRUE),
    distance_cm = mean(distance_cm, na.rm = TRUE),
    mean_acceleration_abs_cm_s2 = mean(mean_acceleration_abs_cm_s2, na.rm = TRUE),
    mean_turn_angle_deg = mean(mean_turn_angle_deg, na.rm = TRUE),
    high_turn_fraction = mean(high_turn_fraction, na.rm = TRUE),
    spatial_spread_cm = mean(spatial_spread_cm, na.rm = TRUE),
    mean_wall_distance_cm = mean(mean_wall_distance_cm, na.rm = TRUE),
    center_fraction = mean(center_fraction, na.rm = TRUE),
    periphery_fraction = mean(periphery_fraction, na.rm = TRUE),
    corner_fraction = mean(corner_fraction, na.rm = TRUE),
    immobile_fraction = mean(immobile_fraction, na.rm = TRUE),
    moving_fraction = mean(moving_fraction, na.rm = TRUE),
    center_entries = mean(center_entries, na.rm = TRUE),
    zone_entropy = mean(zone_entropy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  annotate_motifs(config)

animal_motif_summary <- make_animal_motif_summary(window_scores)
motif_transitions <- make_motif_transition_summary(window_scores)

representative_windows <- window_scores %>%
  group_by(motif) %>%
  mutate(
    motif_PC1_center = mean(PC1, na.rm = TRUE),
    motif_PC2_center = mean(PC2, na.rm = TRUE),
    distance_to_motif_center = sqrt((PC1 - motif_PC1_center)^2 + (PC2 - motif_PC2_center)^2)
  ) %>%
  slice_min(distance_to_motif_center, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  select(file, Batch, Code, window_id, window_start_s, window_end_s, motif, distance_to_motif_center)

model_info <- tibble(
  clustering_method = clustering$method,
  n_motifs = clustering$k,
  model = ifelse(is.null(clustering$model), NA_character_, clustering$model),
  n_windows = nrow(window_scores),
  n_animals = n_distinct(window_scores$file),
  n_pca_components = n_pc,
  pca_variance_target = config$pca_variance_explained_target,
  pca_variance_explained_used = sum(var_exp[seq_len(n_pc)])
)

readr::write_csv(window_features, file.path(output_dir, "oft_unsupervised_window_feature_matrix.csv"))
readr::write_csv(window_scores, file.path(output_dir, "oft_unsupervised_window_scores_and_motifs.csv"))
readr::write_csv(pca_variance, file.path(output_dir, "oft_unsupervised_pca_variance.csv"))
readr::write_csv(loadings_tbl, file.path(output_dir, "oft_unsupervised_pca_loadings.csv"))
readr::write_csv(motif_profiles, file.path(output_dir, "oft_unsupervised_motif_profiles_posthoc_annotations.csv"))
readr::write_csv(animal_motif_summary, file.path(output_dir, "oft_unsupervised_animal_motif_summary.csv"))
readr::write_csv(motif_transitions, file.path(output_dir, "oft_unsupervised_motif_transitions.csv"))
readr::write_csv(representative_windows, file.path(output_dir, "oft_unsupervised_representative_windows.csv"))
readr::write_csv(model_info, file.path(output_dir, "oft_unsupervised_model_info.csv"))

if (!is.null(clustering$silhouette)) {
  readr::write_csv(clustering$silhouette, file.path(output_dir, "oft_unsupervised_kmeans_silhouette_selection.csv"))
}

save_plots(window_scores, loadings_tbl, motif_profiles, output_dir, config)

message("Unsupervised OFT motif discovery complete. Output written to: ", output_dir)
