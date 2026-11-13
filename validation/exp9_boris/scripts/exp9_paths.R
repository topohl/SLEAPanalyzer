# Where the Exp9 tracking and the behaviour records are.
#
#   TRACKING_ROOT  <EXP>/Analysis/Behavior/SLEAP: the SLEAP tracking, one folder per test,
#                  <batch>/<assay>/ with projects/, predictions/, models/, geom/, animal/,
#                  merged/ (<CODE>[_<PHASE>].merged_locs.csv), formatted/
#                  (<CODE>[_<PHASE>].merged_locs_formatted.csv; SocP in formatted/<HAB|S1|S2>/)
#                  and old/. Its qc/ and old/ folders are not batches.
#   RECORDS_ROOT   <EXP>/Raw Data/Behavior: the records of each test, <batch>/<assay>/
#                  (novelLoc.txt for NOR, novelLocHAB/S1/S2.txt for SocP, run sheets, BORIS/, TSE/).
#
# EXP9_TRACKING_ROOT and EXP9_RECORDS_ROOT point them elsewhere. The Python tools use the
# same layout (01_SLEAPcoords/exp9_paths.py).

exp9_root <- function(name, default) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) value <- default
  sub("/+$", "", chartr("\\", "/", value))
}

TRACKING_ROOT <- exp9_root("EXP9_TRACKING_ROOT",
                           "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/SLEAP")
RECORDS_ROOT  <- exp9_root("EXP9_RECORDS_ROOT",
                           "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior")

# The tracking of one test (formatted/ and the rest) and its records (novelLoc sheets and the rest).
tracking_dir <- function(b, assay) file.path(TRACKING_ROOT, b, assay)
records_dir  <- function(b, assay) file.path(RECORDS_ROOT, b, assay)
