# ================================================================
# SLEAP/DLC OFT Optional Autoencoder Motif Discovery
# Author: Tobias Pohl / ChatGPT revision
# Version: 1.0.0
# Date: 2026-05-11
# ================================================================
#
# Run after DLCA_OFT_unsupervised_motifs_v1.0.0.R.
#
# Purpose:
# Exploratory nonlinear unsupervised motif discovery from OFT time-window
# feature matrices using a shallow autoencoder.
#
# Important interpretation note:
# This is secondary/exploratory. The PCA/GMM motif pipeline should remain the
# primary interpretable analysis. This script is intended to test whether a
# nonlinear latent embedding recovers similar or additional behavioral motifs.
#
# No Group, Treatment, Sex, CON/RES/SUS, or stress labels are used for training.
# Biological interpretation is performed post hoc from motif profiles.
# ================================================================

required_packages <- c(
  "dplyr", "tidyr", "purrr", "readr", "tibble", "ggplot2", "cluster"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

if (!requireNamespace("keras", quietly = TRUE)) {
  stop(
    "Package 'keras' is required for this optional autoencoder module.\n",
    "Install it with install.packages('keras') and configure TensorFlow before running.\n",
    "Use DLCA_OFT_unsupervised_motifs_v1.0.0.R as the primary non-deep unsupervised analysis."
  )
}

suppressPackageStartupMessages(library(keras))

config <- list(
  behavior_root = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior",
  global_output_folder = "OFT_SLEAP_output_v1.2.0",
  unsupervised_folder = "unsupervised_motifs_v1.0.0",

  # Autoencoder architecture.
  latent_dim = 2,
  hidden_dim_1 = 16,
  hidden_dim_2 = 8,
  activation = "relu",
  output_activation = "linear",
  dropout_rate = 0.05,

  # Training.
  validation_split = 0.20,
  epochs = 300,
  batch_size = 32,
  learning_rate = 0.001,
  patience = 25,
  random_seed = 123,

  # Clustering latent space.
  k_min = 2,
  k_max = 8,
  kmeans_nstart = 100,
  kmeans_iter_max = 100,

  # Post hoc biological annotation thresholds.
  high_quantile = 0.67,
  low_quantile = 0.33,

  # Plots.
  save_plots = TRUE,
  plot_width_single_col = 85 / 25.4,
  plot_height_single_col = 65 / 25.4,
  plot_width_double_col = 180 / 25.4,
  plot_height_double_col = 100 / 25.4
)

set.seed(config$random_seed)
tensorflow::set_random_seed(config$random_seed)

safe_divide <- function(x, y) {
  ifelse(is.na(y) | y == 0, NA_real_, x / y)
}

entropy <- function(x) {
  tab <- table(x, useNA = "no")
  if (length(tab) == 0) return(NA_real_)
  p <- as.numeric(tab) / sum(tab)
  -sum(p * log2(p))
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

select_autoencoder_features <- function(window_tbl) {
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
    mutate(across(everything(), ~ ifelse(is.infinite(.x), NA_real_, .x))) %>%
    mutate(across(everything(), ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)))

  keep <- vapply(mat, function(x) sd(x, na.rm = TRUE) > 0, logical(1))
  mat <- mat[, keep, drop = FALSE]

  scaled <- scale(as.matrix(mat))
  list(
    matrix = scaled,
    features = colnames(mat),
    center = attr(scaled, "scaled:center"),
    scale = attr(scaled, "scaled:scale")
  )
}

build_autoencoder <- function(input_dim, config) {
  input <- layer_input(shape = input_dim, name = "window_features")

  encoded <- input %>%
    layer_dense(units = config$hidden_dim_1, activation = config$activation, name = "encoder_dense_1") %>%
    layer_dropout(rate = config$dropout_rate) %>%
    layer_dense(units = config$hidden_dim_2, activation = config$activation, name = "encoder_dense_2") %>%
    layer_dense(units = config$latent_dim, activation = "linear", name = "latent")

  decoded <- encoded %>%
    layer_dense(units = config$hidden_dim_2, activation = config$activation, name = "decoder_dense_1") %>%
    layer_dense(units = config$hidden_dim_1, activation = config$activation, name = "decoder_dense_2") %>%
    layer_dense(units = input_dim, activation = config$output_activation, name = "reconstruction")

  autoencoder <- keras_model(inputs = input, outputs = decoded)
  encoder <- keras_model(inputs = input, outputs = encoded)

  autoencoder %>% compile(
    optimizer = optimizer_adam(learning_rate = config$learning_rate),
    loss = "mse",
    metrics = c("mae")
  )

  list(autoencoder = autoencoder, encoder = encoder)
}

choose_k_latent <- function(latent_mat, config) {
  k_values <- config$k_min:min(config$k_max, nrow(latent_mat) - 1)

  sil_tbl <- purrr::map_dfr(k_values, function(k) {
    km <- stats::kmeans(
      latent_mat,
      centers = k,
      nstart = config$kmeans_nstart,
      iter.max = config$kmeans_iter_max
    )
    sil <- cluster::silhouette(km$cluster, stats::dist(latent_mat))
    tibble(k = k, mean_silhouette = mean(sil[, "sil_width"]))
  })

  best_k <- sil_tbl$k[which.max(sil_tbl$mean_silhouette)]
  km <- stats::kmeans(
    latent_mat,
    centers = best_k,
    nstart = config$kmeans_nstart,
    iter.max = config$kmeans_iter_max
  )

  list(k = best_k, cluster = as.integer(km$cluster), silhouette = sil_tbl, model = km)
}

annotate_autoencoder_motifs <- function(profile_tbl, config) {
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
      autoencoder_motif_annotation = case_when(
        high_immobility & (low_wall_distance | high_corner) ~ "wall/corner immobility",
        high_speed & low_wall_distance & !high_center ~ "thigmotactic locomotion",
        high_center & high_wall_distance & high_spread ~ "center exploration",
        high_burstiness & high_turning ~ "fragmented burst exploration",
        low_speed & !high_immobility & low_wall_distance ~ "cautious wall-local exploration",
        low_speed & !high_immobility & high_wall_distance ~ "slow centerward scanning",
        TRUE ~ "mixed nonlinear motif"
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

make_animal_summary <- function(window_scores) {
  motif_counts <- window_scores %>%
    count(file, Batch, Code, autoencoder_motif, name = "n_windows") %>%
    group_by(file, Batch, Code) %>%
    mutate(motif_fraction = n_windows / sum(n_windows)) %>%
    ungroup()

  wide <- motif_counts %>%
    select(file, Batch, Code, autoencoder_motif, motif_fraction) %>%
    pivot_wider(
      names_from = autoencoder_motif,
      values_from = motif_fraction,
      names_prefix = "autoencoder_motif_fraction_",
      values_fill = 0
    )

  entropy_tbl <- motif_counts %>%
    group_by(file, Batch, Code) %>%
    summarise(
      autoencoder_motif_entropy = entropy(rep(autoencoder_motif, n_windows)),
      dominant_autoencoder_motif = autoencoder_motif[which.max(n_windows)],
      dominant_autoencoder_motif_fraction = max(motif_fraction),
      .groups = "drop"
    )

  wide %>% left_join(entropy_tbl, by = c("file", "Batch", "Code"))
}

make_transition_summary <- function(window_scores) {
  window_scores %>%
    arrange(file, window_id) %>%
    group_by(file, Batch, Code) %>%
    mutate(next_autoencoder_motif = lead(autoencoder_motif)) %>%
    ungroup() %>%
    filter(!is.na(next_autoencoder_motif)) %>%
    count(file, Batch, Code, autoencoder_motif, next_autoencoder_motif, name = "n") %>%
    group_by(file, Batch, Code, autoencoder_motif) %>%
    mutate(transition_probability = n / sum(n)) %>%
    ungroup()
}

compare_with_pca_gmm_motifs <- function(autoencoder_scores, unsupervised_dir) {
  pca_file <- file.path(unsupervised_dir, "oft_unsupervised_window_scores_and_motifs.csv")
  if (!file.exists(pca_file)) {
    return(tibble(note = "PCA/GMM motif file not found; comparison skipped."))
  }

  pca_scores <- readr::read_csv(pca_file, show_col_types = FALSE) %>%
    select(window_uid, pca_gmm_motif = motif)

  joined <- autoencoder_scores %>%
    select(window_uid, autoencoder_motif) %>%
    inner_join(pca_scores, by = "window_uid")

  if (nrow(joined) == 0) {
    return(tibble(note = "No overlapping window_uid values; comparison skipped."))
  }

  joined %>%
    count(autoencoder_motif, pca_gmm_motif, name = "n") %>%
    group_by(autoencoder_motif) %>%
    mutate(row_fraction = n / sum(n)) %>%
    ungroup()
}

save_autoencoder_plots <- function(window_scores, motif_profiles, training_history, output_dir, config) {
  if (!isTRUE(config$save_plots)) return(invisible(NULL))

  plot_dir <- file.path(output_dir, "autoencoder_plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  if (all(c("AE1", "AE2") %in% names(window_scores))) {
    p <- ggplot(window_scores, aes(x = AE1, y = AE2, color = factor(autoencoder_motif))) +
      geom_point(size = 0.7, alpha = 0.55) +
      labs(title = "Autoencoder latent OFT motif space", x = "AE latent 1", y = "AE latent 2") +
      make_theme_nature()
    ggsave(file.path(plot_dir, "oft_autoencoder_latent_motif_space.svg"), p,
           width = config$plot_width_single_col, height = config$plot_height_single_col)
  }

  hist_tbl <- as_tibble(training_history$metrics) %>%
    mutate(epoch = seq_len(n())) %>%
    select(epoch, everything()) %>%
    pivot_longer(-epoch, names_to = "metric", values_to = "value")

  p <- ggplot(hist_tbl, aes(x = epoch, y = value, group = metric)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(title = "Autoencoder training history", x = "Epoch", y = "Value") +
    make_theme_nature()
  ggsave(file.path(plot_dir, "oft_autoencoder_training_history.svg"), p,
         width = config$plot_width_double_col, height = config$plot_height_single_col)

  profile_long <- motif_profiles %>%
    select(autoencoder_motif, autoencoder_motif_annotation, mean_speed_cm_s, center_fraction,
           immobile_fraction, mean_wall_distance_cm, mean_turn_angle_deg,
           spatial_spread_cm, sd_speed_cm_s) %>%
    pivot_longer(-c(autoencoder_motif, autoencoder_motif_annotation), names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    mutate(z_value = as.numeric(scale(value))) %>%
    ungroup()

  p <- ggplot(profile_long, aes(x = feature, y = factor(autoencoder_motif), fill = z_value)) +
    geom_tile() +
    coord_flip() +
    labs(title = "Autoencoder motif profiles", x = NULL, y = "Motif") +
    make_theme_nature()
  ggsave(file.path(plot_dir, "oft_autoencoder_motif_profile_heatmap.svg"), p,
         width = config$plot_width_double_col, height = config$plot_height_single_col)

  invisible(NULL)
}

# -------------------------------
# Main execution
# -------------------------------

unsupervised_dir <- file.path(config$behavior_root, config$global_output_folder, config$unsupervised_folder)
window_file <- file.path(unsupervised_dir, "oft_unsupervised_window_feature_matrix.csv")

if (!file.exists(window_file)) {
  stop(
    "Missing unsupervised window feature matrix. Run DLCA_OFT_unsupervised_motifs_v1.0.0.R first.\n",
    "Expected: ", window_file
  )
}

output_dir <- file.path(config$behavior_root, config$global_output_folder, "autoencoder_motifs_v1.0.0")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

window_tbl <- readr::read_csv(window_file, show_col_types = FALSE)
features <- select_autoencoder_features(window_tbl)
scaled <- prepare_scaled_matrix(window_tbl, features)
x <- scaled$matrix

if (nrow(x) < 20) {
  stop("Too few windows for autoencoder training. Use PCA/GMM motif discovery instead.")
}

message("Training autoencoder on ", nrow(x), " windows and ", ncol(x), " features...")
models <- build_autoencoder(ncol(x), config)

early_stop <- callback_early_stopping(
  monitor = "val_loss",
  patience = config$patience,
  restore_best_weights = TRUE
)

history <- models$autoencoder %>% fit(
  x = x,
  y = x,
  epochs = config$epochs,
  batch_size = config$batch_size,
  validation_split = config$validation_split,
  callbacks = list(early_stop),
  verbose = 0
)

latent <- models$encoder %>% predict(x, verbose = 0)
latent_tbl <- as_tibble(latent)
colnames(latent_tbl) <- paste0("AE", seq_len(ncol(latent_tbl)))

clustering <- choose_k_latent(latent, config)

reconstruction <- models$autoencoder %>% predict(x, verbose = 0)
reconstruction_mse <- rowMeans((x - reconstruction)^2)

window_scores <- bind_cols(window_tbl, latent_tbl) %>%
  mutate(
    autoencoder_motif = paste0("AE_M", clustering$cluster),
    autoencoder_motif_numeric = clustering$cluster,
    reconstruction_mse = reconstruction_mse,
    autoencoder_latent_dim = config$latent_dim
  )

motif_profiles <- window_scores %>%
  group_by(autoencoder_motif) %>%
  summarise(
    n_windows = n(),
    n_animals = n_distinct(file),
    reconstruction_mse_mean = mean(reconstruction_mse, na.rm = TRUE),
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
  annotate_autoencoder_motifs(config)

animal_summary <- make_animal_summary(window_scores)
transition_summary <- make_transition_summary(window_scores)
comparison_tbl <- compare_with_pca_gmm_motifs(window_scores, unsupervised_dir)

training_history_tbl <- as_tibble(history$metrics) %>%
  mutate(epoch = seq_len(n())) %>%
  select(epoch, everything())

model_info <- tibble(
  model = "shallow dense autoencoder",
  n_windows = nrow(x),
  n_features = ncol(x),
  latent_dim = config$latent_dim,
  hidden_dim_1 = config$hidden_dim_1,
  hidden_dim_2 = config$hidden_dim_2,
  dropout_rate = config$dropout_rate,
  validation_split = config$validation_split,
  epochs_requested = config$epochs,
  epochs_run = nrow(training_history_tbl),
  final_loss = tail(training_history_tbl$loss, 1),
  final_val_loss = ifelse("val_loss" %in% names(training_history_tbl), tail(training_history_tbl$val_loss, 1), NA_real_),
  cluster_method = "k-means on autoencoder latent space",
  n_autoencoder_motifs = clustering$k,
  random_seed = config$random_seed
)

readr::write_csv(window_scores, file.path(output_dir, "oft_autoencoder_window_scores_and_motifs.csv"))
readr::write_csv(motif_profiles, file.path(output_dir, "oft_autoencoder_motif_profiles_posthoc_annotations.csv"))
readr::write_csv(animal_summary, file.path(output_dir, "oft_autoencoder_animal_motif_summary.csv"))
readr::write_csv(transition_summary, file.path(output_dir, "oft_autoencoder_motif_transitions.csv"))
readr::write_csv(comparison_tbl, file.path(output_dir, "oft_autoencoder_vs_pca_gmm_motif_overlap.csv"))
readr::write_csv(training_history_tbl, file.path(output_dir, "oft_autoencoder_training_history.csv"))
readr::write_csv(clustering$silhouette, file.path(output_dir, "oft_autoencoder_kmeans_silhouette_selection.csv"))
readr::write_csv(model_info, file.path(output_dir, "oft_autoencoder_model_info.csv"))

save_autoencoder_plots(window_scores, motif_profiles, history, output_dir, config)

message("Autoencoder OFT motif discovery complete. Output written to: ", output_dir)
