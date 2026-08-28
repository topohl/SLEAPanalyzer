# Shared configuration schema and defaults for the batch assay workflows.
#
# Every assay script resolves its configuration through here, so no production
# script contains a machine-specific path, a fixed acquisition rate or an
# undocumented threshold.

#' Fields common to every assay configuration.
common_assay_schema <- function() {
  list(
    input_dir = config_field("path"),
    output_dir = config_field("path"),
    fps = config_field("number", positive = TRUE),
    batches = config_field("string_vector", required = FALSE),

    arena_width_cm = config_field("number", positive = TRUE),
    arena_height_cm = config_field("number", positive = TRUE),
    arena_corner_names = config_field("string_vector", length = 4),

    likelihood_cutoff = config_field("number", required = FALSE),
    max_interpolation_gap_s = config_field("number", positive = TRUE, allow_zero = TRUE),
    max_plausible_speed_cm_s = config_field("number", positive = TRUE, required = FALSE),

    qc_min_valid_fraction = config_field("number", positive = TRUE, allow_zero = TRUE),
    qc_max_longest_gap_s = config_field("number", positive = TRUE, allow_zero = TRUE),
    qc_max_interpolated_fraction = config_field("number", positive = TRUE, allow_zero = TRUE),

    movement_cutoff_cm_s = config_field("number", positive = TRUE),
    integration_period_frames = config_field("integer", required = FALSE),

    min_bout_s = config_field("number", positive = TRUE, allow_zero = TRUE, required = FALSE),
    max_gap_s = config_field("number", positive = TRUE, allow_zero = TRUE, required = FALSE),

    metadata_dir = config_field("path", required = FALSE),
    animal_id_code_file = config_field("path", required = FALSE),
    write_manifest = config_field("logical", required = FALSE)
  )
}

#' Defaults shared by every assay.
#'
#' Thresholds that carry a biological meaning are deliberately absent here and
#' must be stated in the assay configuration, so that a value is never applied
#' silently just because it happened to be the default.
common_assay_defaults <- function() {
  list(
    fps = NULL,
    batches = NULL,
    arena_corner_names = c("tl", "tr", "br", "bl"),
    likelihood_cutoff = NULL,
    max_interpolation_gap_s = 0.2,
    max_plausible_speed_cm_s = NULL,
    qc_min_valid_fraction = 0.8,
    qc_max_longest_gap_s = 5,
    qc_max_interpolated_fraction = 0.2,
    movement_cutoff_cm_s = 5,
    integration_period_frames = 5,
    min_bout_s = 0,
    max_gap_s = 0,
    write_manifest = TRUE
  )
}

#' Assay-specific schema additions.
assay_schema <- function(assay) {
  extra <- switch(
    assay,
    OFT = list(),
    EPM = list(
      zone_file = config_field("path"),
      nose_dips = config_field("logical", required = FALSE)
    ),
    NOR = list(
      novel_location_file = config_field("string"),
      contact_geometry = config_field(
        "string", choices = c("radial", "box", "legacy_asymmetric")
      ),
      contact_distance_cm = config_field("number", positive = TRUE),
      body_exclusion_distance_cm = config_field("number", positive = TRUE),
      object_box_width_cm = config_field("number", positive = TRUE, required = FALSE),
      object_box_height_cm = config_field("number", positive = TRUE, required = FALSE),
      contact_angle_deg = config_field("number_vector", length = 2),
      proximity_range_cm = config_field("number_vector", length = 2),
      proximity_angle_deg = config_field("number_vector", length = 2),
      report_rearing = config_field("logical", required = FALSE),
      rearing_spine_distance_cm = config_field("number", positive = TRUE, required = FALSE)
    ),
    SocP = list(
      phases = config_field("string_vector"),
      novel_location_file_prefix = config_field("string"),
      contact_distance_cm = config_field("number", positive = TRUE),
      body_exclusion_distance_cm = config_field("number", positive = TRUE),
      proximity_range_cm = config_field("number_vector", length = 2),
      require_orientation = config_field("logical", required = FALSE)
    ),
    stop("Unknown assay: ", assay)
  )
  c(common_assay_schema(), extra)
}

#' Load and validate an assay configuration file.
#'
#' @param path path to a YAML configuration file
#' @param assay one of "OFT", "EPM", "NOR", "SocP"
load_assay_config <- function(path, assay) {
  load_analysis_config(
    path,
    defaults = common_assay_defaults(),
    schema = assay_schema(assay),
    path_fields = c(
      "input_dir", "output_dir", "metadata_dir", "animal_id_code_file", "zone_file"
    )
  )
}

#' Resolve the configuration path for a run.
#'
#' Order of precedence: an explicit argument, the SLEAP_ANALYZER_CONFIG
#' environment variable, then the bundled example for the assay. Making this
#' explicit is what removes machine-specific paths from the scripts.
resolve_config_path <- function(path = NULL, assay, script_dir = getwd()) {
  if (!is.null(path) && nzchar(path)) return(path)
  from_env <- Sys.getenv("SLEAP_ANALYZER_CONFIG")
  if (nzchar(from_env)) return(from_env)
  candidate <- file.path(
    dirname(script_dir), "config", paste0(tolower(assay), ".example.yaml")
  )
  if (file.exists(candidate)) {
    message(
      "No configuration supplied; using the bundled example at ", candidate,
      ". Copy it and set SLEAP_ANALYZER_CONFIG to your own file."
    )
    return(candidate)
  }
  stop(
    "No configuration file found. Set SLEAP_ANALYZER_CONFIG to a YAML file, ",
    "or copy config/", tolower(assay), ".example.yaml and edit it."
  )
}
