# ================================================================
# SLEAP/DLC OFT Advanced Metrics Extension
# Author: Tobias Pohl / ChatGPT revision
# Version: 1.0.0
# Date: 2026-05-11
# ================================================================
#
# Run this after DLCA_OFT v1.2.0.R.
# It consumes the frame-, summary-, and enhanced-metrics outputs produced by
# DLCA_OFT v1.2.0.R and adds:
# - zone-transition matrices
# - fragmentation / state-switching metrics
# - sustained center-entry metrics
# - locomotion-normalized OFT metrics
# - burst locomotion metrics
# - QC sensitivity tables excluding flagged animals
# - optional group statistics if metadata columns are present
#
# Interpretation:
# These outputs should be interpreted as an OFT behavioral profile, not as a
# pure anxiety score. Locomotor confounding is explicitly quantified via
# distance-normalized and speed-normalized metrics.
# ================================================================

required_packages <- c(
  "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble", "ggplot2"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

optional_packages <- c("lmerTest", "emmeans", "effectsize")
for (pkg in optional_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Optional package not installed: ", pkg, ". Statistics module will skip features requiring it.")
  }
}

config <- list(
  behavior_root = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior",
  global_output_folder = "OFT_SLEAP_output_v1.2.0",
  batch_output_subdir = file.path("OFT", "SLEAP", "output_v1.2.0"),
  batches = c("B1"),
  fps = 30,

  sustained_center_min_s = c(1, 3),
  burst_speed_threshold_cm_s = 20,
  burst_min_duration_s = 0.20,
  immobility_min_duration_s = 1.0,
  qc_exclusion_column = "enhanced_qc_flag",

  group_columns_priority = c("Group", "Phenotype", "Treatment"),
  covariate_columns_priority = c("Sex", "Batch"),
  statistics_metrics = c(
    "center_time_percent", "center_entries", "center_latency_s",
    "distance_cm", "mean_speed_cm_s", "mean_wall_distance_cm",
    "oft_center_exploration_score", "zone_fragmentation_rate_per_min",
    "movement_fragmentation_rate_per_min", "sustained_center_entries_3s",
    "center_time_per_100cm", "wall_5cm_percent_per_100cm",
    "burst_count", "burst_rate_per_min"
  ),

  save_plots = TRUE,
  plot_width_single_col = 85 / 25.4,
  plot_height_single_col = 65 / 25.4,
  plot_width_double_col = 180 / 25.4,
  plot_height_double_col = 85 / 25.4
)

safe_divide <- function(x, y) {
  ifelse(is.na(y) | y == 0, NA_real_, x / y)
}

entry_count <- function(x) {
  x <- tidyr::replace_na(as.logical(x), FALSE)
  sum(x & !dplyr::lag(x, default = FALSE), na.rm = TRUE)
}

bout_durations <- function(x, fps, min_duration_s = 0) {
  x <- tidyr::replace_na(as.logical(x), FALSE)
  if (length(x) == 0) return(numeric())
  r <- rle(x)
  d <- r$lengths[r$values] / fps
  d[d >= min_duration_s]
}

count_bouts <- function(x, fps, min_duration_s = 0) {
  length(bout_durations(x, fps, min_duration_s))
}

mean_bout <- function(x, fps, min_duration_s = 0) {
  d <- bout_durations(x, fps, min_duration_s)
  if (length(d) == 0) NA_real_ else mean(d)
}

max_bout <- function(x, fps, min_duration_s = 0) {
  d <- bout_durations(x, fps, min_duration_s)
  if (length(d) == 0) NA_real_ else max(d)
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

classify_zone <- function(frame_tbl) {
  frame_tbl %>%
    mutate(
      oft_zone_state = case_when(
        isTRUE(in_center) | in_center ~ "center",
        isTRUE(in_corners) | in_corners ~ "corner",
        isTRUE(in_periphery) | in_periphery ~ "periphery",
        TRUE ~ "other"
      ),
      movement_state = case_when(
        isTRUE(is_moving) | is_moving ~ "moving",
        isTRUE(is_immobile) | is_immobile ~ "immobile",
        TRUE ~ "unknown"
      )
    )
}

make_transition_matrix <- function(frame_tbl, state_col, fps) {
  state_sym <- rlang::sym(state_col)
  dat <- frame_tbl %>%
    arrange(frame) %>%
    transmute(
      from = as.character(!!state_sym),
      to = dplyr::lead(as.character(!!state_sym))
    ) %>%
    filter(!is.na(from), !is.na(to), from != "unknown", to != "unknown")

  if (nrow(dat) == 0) {
    return(tibble(from = character(), to = character(), n = integer(), probability = numeric()))
  }

  dat %>%
    count(from, to, name = "n") %>%
    group_by(from) %>%
    mutate(probability = n / sum(n)) %>%
    ungroup()
}

make_transition_wide <- function(transition_tbl, prefix) {
  if (nrow(transition_tbl) == 0) return(tibble())
  transition_tbl %>%
    mutate(metric = paste0(prefix, "_", from, "_to_", to)) %>%
    select(metric, probability) %>%
    pivot_wider(names_from = metric, values_from = probability)
}

make_fragmentation_metrics <- function(frame_tbl, config) {
  dat <- classify_zone(frame_tbl) %>% arrange(frame)
  total_min <- nrow(dat) / config$fps / 60

  zone_switches <- sum(dat$oft_zone_state != dplyr::lag(dat$oft_zone_state), na.rm = TRUE)
  movement_switches <- sum(dat$movement_state != dplyr::lag(dat$movement_state), na.rm = TRUE)

  tibble(
    zone_switch_count = zone_switches,
    zone_fragmentation_rate_per_min = safe_divide(zone_switches, total_min),
    movement_switch_count = movement_switches,
    movement_fragmentation_rate_per_min = safe_divide(movement_switches, total_min),
    moving_to_immobile_switches = sum(dat$movement_state == "immobile" & dplyr::lag(dat$movement_state) == "moving", na.rm = TRUE),
    immobile_to_moving_switches = sum(dat$movement_state == "moving" & dplyr::lag(dat$movement_state) == "immobile", na.rm = TRUE)
  )
}

make_sustained_center_metrics <- function(frame_tbl, config) {
  purrr::map_dfr(config$sustained_center_min_s, function(min_s) {
    d <- bout_durations(frame_tbl$in_center, config$fps, min_s)
    tibble(
      min_s = min_s,
      sustained_center_entries = length(d),
      sustained_center_time_s = sum(d),
      sustained_center_mean_bout_s = ifelse(length(d) == 0, NA_real_, mean(d)),
      sustained_center_max_bout_s = ifelse(length(d) == 0, NA_real_, max(d))
    )
  }) %>%
    pivot_wider(
      names_from = min_s,
      values_from = c(
        sustained_center_entries,
        sustained_center_time_s,
        sustained_center_mean_bout_s,
        sustained_center_max_bout_s
      ),
      names_glue = "{.value}_{min_s}s"
    )
}

make_locomotion_normalized_metrics <- function(summary_tbl, enhanced_tbl) {
  joined <- summary_tbl %>%
    left_join(enhanced_tbl, by = c("file", "Batch", "Code"), suffix = c("", ".enh"))

  joined %>%
    transmute(
      file, Batch, Code,
      center_time_per_100cm = safe_divide(center_time_s, distance_cm) * 100,
      center_entries_per_100cm = safe_divide(center_entries, distance_cm) * 100,
      corner_time_per_100cm = safe_divide(corner_time_s, distance_cm) * 100,
      mean_wall_distance_per_100cm = safe_divide(mean_wall_distance_cm, distance_cm) * 100,
      wall_5cm_percent_per_100cm = if ("wall_time_percent_5cm" %in% names(joined)) {
        safe_divide(wall_time_percent_5cm, distance_cm) * 100
      } else NA_real_,
      wall_10cm_percent_per_100cm = if ("wall_time_percent_10cm" %in% names(joined)) {
        safe_divide(wall_time_percent_10cm, distance_cm) * 100
      } else NA_real_,
      center_time_per_min_moving = safe_divide(center_time_s, moving_time_s / 60),
      wall_distance_per_mean_speed = safe_divide(mean_wall_distance_cm, mean_speed_cm_s)
    )
}

make_burst_metrics <- function(frame_tbl, config) {
  is_burst <- tidyr::replace_na(frame_tbl$speed_cm_s >= config$burst_speed_threshold_cm_s, FALSE)
  burst_d <- bout_durations(is_burst, config$fps, config$burst_min_duration_s)
  total_min <- nrow(frame_tbl) / config$fps / 60

  tibble(
    burst_speed_threshold_cm_s = config$burst_speed_threshold_cm_s,
    burst_min_duration_s = config$burst_min_duration_s,
    burst_count = length(burst_d),
    burst_rate_per_min = safe_divide(length(burst_d), total_min),
    burst_total_time_s = sum(burst_d),
    burst_mean_duration_s = ifelse(length(burst_d) == 0, NA_real_, mean(burst_d)),
    burst_max_duration_s = ifelse(length(burst_d) == 0, NA_real_, max(burst_d)),
    high_speed_time_percent = mean(is_burst, na.rm = TRUE) * 100,
    speed_p95_cm_s = stats::quantile(frame_tbl$speed_cm_s, 0.95, na.rm = TRUE, names = FALSE),
    speed_p99_cm_s = stats::quantile(frame_tbl$speed_cm_s, 0.99, na.rm = TRUE, names = FALSE)
  )
}

make_advanced_row <- function(frame_tbl, summary_tbl, enhanced_tbl, config) {
  dat <- classify_zone(frame_tbl)
  zone_trans <- make_transition_matrix(dat, "oft_zone_state", config$fps)
  move_trans <- make_transition_matrix(dat, "movement_state", config$fps)

  bind_cols(
    summary_tbl %>% select(file, Batch, Code),
    make_fragmentation_metrics(dat, config),
    make_sustained_center_metrics(dat, config),
    make_burst_metrics(dat, config),
    make_transition_wide(zone_trans, "zone_transition_prob"),
    make_transition_wide(move_trans, "movement_transition_prob")
  )
}

get_frame_files <- function(config) {
  purrr::map(config$batches, function(batch) {
    frame_dir <- file.path(config$behavior_root, batch, config$batch_output_subdir, "frames")
    if (!dir.exists(frame_dir)) return(character())
    list.files(frame_dir, pattern = "_frames\\.csv$", full.names = TRUE)
  }) %>% unlist()
}

read_existing_outputs <- function(config) {
  global_output <- file.path(config$behavior_root, config$global_output_folder)
  summary_file <- file.path(global_output, "all_batches_summary.csv")
  enhanced_file <- file.path(global_output, "all_batches_enhanced_metrics.csv")

  if (!file.exists(summary_file)) stop("Missing all_batches_summary.csv. Run DLCA_OFT v1.2.0.R first.")
  if (!file.exists(enhanced_file)) stop("Missing all_batches_enhanced_metrics.csv. Run DLCA_OFT v1.2.0.R first.")

  list(
    global_output = global_output,
    summary = readr::read_csv(summary_file, show_col_types = FALSE),
    enhanced = readr::read_csv(enhanced_file, show_col_types = FALSE)
  )
}

make_qc_sensitivity_table <- function(full_tbl, config) {
  qc_col <- config$qc_exclusion_column
  metric_cols <- names(full_tbl)[vapply(full_tbl, is.numeric, logical(1))]

  if (!qc_col %in% names(full_tbl)) {
    return(tibble(note = paste("QC column not found:", qc_col)))
  }

  purrr::map_dfr(metric_cols, function(metric) {
    all_x <- full_tbl[[metric]]
    keep_x <- full_tbl %>% filter(!.data[[qc_col]]) %>% pull(.data[[metric]])

    tibble(
      metric = metric,
      n_all = sum(!is.na(all_x)),
      n_qc_pass = sum(!is.na(keep_x)),
      mean_all = mean(all_x, na.rm = TRUE),
      mean_qc_pass = mean(keep_x, na.rm = TRUE),
      delta_qc_pass_minus_all = mean_qc_pass - mean_all,
      sd_all = stats::sd(all_x, na.rm = TRUE),
      sd_qc_pass = stats::sd(keep_x, na.rm = TRUE)
    )
  })
}

run_optional_group_statistics <- function(full_tbl, config) {
  group_col <- config$group_columns_priority[config$group_columns_priority %in% names(full_tbl)][1]
  if (is.na(group_col)) {
    return(tibble(note = "No group column found. Add Group, Phenotype, or Treatment to outputs/metadata before running group statistics."))
  }

  available_metrics <- intersect(config$statistics_metrics, names(full_tbl))
  available_metrics <- available_metrics[vapply(full_tbl[available_metrics], is.numeric, logical(1))]

  if (length(available_metrics) == 0) {
    return(tibble(note = "No configured numeric statistics metrics found."))
  }

  covariates <- setdiff(intersect(config$covariate_columns_priority, names(full_tbl)), group_col)
  use_lmm <- requireNamespace("lmerTest", quietly = TRUE) && "ID" %in% names(full_tbl)
  use_emm <- requireNamespace("emmeans", quietly = TRUE)

  purrr::map_dfr(available_metrics, function(metric) {
    dat <- full_tbl %>%
      select(all_of(c(metric, group_col, covariates, "ID"["ID" %in% names(full_tbl)]))) %>%
      filter(!is.na(.data[[metric]]), !is.na(.data[[group_col]]))

    if (nrow(dat) < 4 || dplyr::n_distinct(dat[[group_col]]) < 2) {
      return(tibble(metric = metric, note = "Insufficient data or fewer than two groups."))
    }

    rhs <- paste(c(group_col, covariates), collapse = " + ")

    if (use_lmm) {
      form <- stats::as.formula(paste(metric, "~", rhs, "+ (1|ID)"))
      fit <- tryCatch(lmerTest::lmer(form, data = dat), error = function(e) NULL)
    } else {
      form <- stats::as.formula(paste(metric, "~", rhs))
      fit <- tryCatch(stats::lm(form, data = dat), error = function(e) NULL)
    }

    if (is.null(fit)) return(tibble(metric = metric, note = "Model failed."))

    an <- tryCatch(as.data.frame(stats::anova(fit)), error = function(e) NULL)
    an_tbl <- if (!is.null(an)) {
      an %>%
        tibble::rownames_to_column("term") %>%
        filter(term == group_col) %>%
        transmute(
          metric = metric,
          model = ifelse(use_lmm, "lmer", "lm"),
          term = term,
          statistic = dplyr::coalesce(.data[[intersect(c("F value", "F.value"), names(.))[1]]], NA_real_),
          p_raw = dplyr::coalesce(.data[[intersect(c("Pr(>F)", "Pr..F."), names(.))[1]]], NA_real_),
          n_used = nrow(dat),
          note = NA_character_
        )
    } else {
      tibble(metric = metric, model = ifelse(use_lmm, "lmer", "lm"), term = group_col, statistic = NA_real_, p_raw = NA_real_, n_used = nrow(dat), note = "ANOVA extraction failed.")
    }

    if (use_emm) {
      emm_tbl <- tryCatch({
        emmeans::emmeans(fit, specs = stats::as.formula(paste("~", group_col))) %>%
          emmeans::contrast(method = "pairwise", adjust = "holm") %>%
          as.data.frame() %>%
          as_tibble() %>%
          transmute(
            metric = metric,
            model = ifelse(use_lmm, "lmer", "lm"),
            term = paste0(group_col, " pairwise"),
            contrast = contrast,
            estimate = estimate,
            SE = SE,
            df = df,
            statistic = t.ratio,
            p_raw = p.value,
            p_holm = p.value,
            n_used = nrow(dat),
            note = NA_character_
          )
      }, error = function(e) tibble(metric = metric, note = paste("emmeans failed:", conditionMessage(e))))

      bind_rows(an_tbl, emm_tbl)
    } else {
      an_tbl
    }
  }) %>%
    group_by(term) %>%
    mutate(p_holm_global = ifelse(!is.na(p_raw), p.adjust(p_raw, method = "holm"), NA_real_)) %>%
    ungroup()
}

save_advanced_plots <- function(full_tbl, output_dir, config) {
  if (!isTRUE(config$save_plots)) return(invisible(NULL))
  plot_dir <- file.path(output_dir, "advanced_plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  if (all(c("zone_fragmentation_rate_per_min", "center_time_percent") %in% names(full_tbl))) {
    p <- ggplot(full_tbl, aes(x = zone_fragmentation_rate_per_min, y = center_time_percent)) +
      geom_point(size = 1.8, alpha = 0.85) +
      labs(
        title = "Center time versus zone fragmentation",
        x = "Zone switches per min",
        y = "Center time (%)"
      ) +
      make_theme_nature()
    ggsave(file.path(plot_dir, "advanced_center_vs_zone_fragmentation.svg"), p,
           width = config$plot_width_single_col, height = config$plot_height_single_col)
  }

  if (all(c("burst_rate_per_min", "distance_cm") %in% names(full_tbl))) {
    p <- ggplot(full_tbl, aes(x = burst_rate_per_min, y = distance_cm)) +
      geom_point(size = 1.8, alpha = 0.85) +
      labs(
        title = "Burst locomotion versus total distance",
        x = "Burst rate per min",
        y = "Distance (cm)"
      ) +
      make_theme_nature()
    ggsave(file.path(plot_dir, "advanced_burst_rate_vs_distance.svg"), p,
           width = config$plot_width_single_col, height = config$plot_height_single_col)
  }

  if (all(c("center_time_per_100cm", "wall_5cm_percent_per_100cm") %in% names(full_tbl))) {
    p <- ggplot(full_tbl, aes(x = wall_5cm_percent_per_100cm, y = center_time_per_100cm)) +
      geom_point(size = 1.8, alpha = 0.85) +
      labs(
        title = "Locomotion-normalized OFT profile",
        x = "Wall 5 cm % per 100 cm",
        y = "Center time per 100 cm"
      ) +
      make_theme_nature()
    ggsave(file.path(plot_dir, "advanced_locomotion_normalized_profile.svg"), p,
           width = config$plot_width_single_col, height = config$plot_height_single_col)
  }

  invisible(NULL)
}

# -------------------------------
# Main execution
# -------------------------------

existing <- read_existing_outputs(config)
summary_tbl <- existing$summary
enhanced_tbl <- existing$enhanced
output_dir <- file.path(existing$global_output, "advanced_metrics_v1.0.0")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

frame_files <- get_frame_files(config)
if (length(frame_files) == 0) stop("No *_frames.csv files found. Run DLCA_OFT v1.2.0.R first.")

advanced_rows <- purrr::map_dfr(frame_files, function(frame_file) {
  frame_tbl <- readr::read_csv(frame_file, show_col_types = FALSE)
  key <- frame_tbl %>% distinct(file, Batch, Code) %>% slice(1)

  summary_one <- summary_tbl %>% semi_join(key, by = c("file", "Batch", "Code"))
  enhanced_one <- enhanced_tbl %>% semi_join(key, by = c("file", "Batch", "Code"))

  if (nrow(summary_one) == 0) stop("No summary row for frame file: ", frame_file)
  if (nrow(enhanced_one) == 0) enhanced_one <- key

  make_advanced_row(frame_tbl, summary_one, enhanced_one, config)
})

loc_norm <- make_locomotion_normalized_metrics(summary_tbl, enhanced_tbl)

advanced_full <- summary_tbl %>%
  left_join(enhanced_tbl, by = c("file", "Batch", "Code"), suffix = c("", ".enh")) %>%
  left_join(advanced_rows, by = c("file", "Batch", "Code")) %>%
  left_join(loc_norm, by = c("file", "Batch", "Code"))

qc_sensitivity <- make_qc_sensitivity_table(advanced_full, config)
statistics_tbl <- run_optional_group_statistics(advanced_full, config)

readr::write_csv(advanced_rows, file.path(output_dir, "all_batches_advanced_behavior_metrics.csv"))
readr::write_csv(loc_norm, file.path(output_dir, "all_batches_locomotion_normalized_metrics.csv"))
readr::write_csv(advanced_full, file.path(output_dir, "all_batches_oft_full_advanced_profile.csv"))
readr::write_csv(qc_sensitivity, file.path(output_dir, "all_batches_qc_sensitivity_metrics.csv"))
readr::write_csv(statistics_tbl, file.path(output_dir, "all_batches_oft_statistics.csv"))

save_advanced_plots(advanced_full, output_dir, config)

message("Advanced OFT metrics complete. Output written to: ", output_dir)
