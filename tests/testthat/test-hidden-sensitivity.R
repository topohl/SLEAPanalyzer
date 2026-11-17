# The hidden-contact sensitivity (validation/exp9_boris/scripts/
# 35_hidden_sensitivity.R) on two small SocP output trees.

sensitivity_script <- file.path(repo_root, "validation", "exp9_boris", "scripts", "35_hidden_sensitivity.R")

# Recordings scored with hidden time as contact. A002_S1 is hidden 4 s at the
# right cup (novel), A003_S1 2 s at the left cup (novel); A005_S1 has no
# contact, A006_S1 no group.
hidden_as_contact <- function() {
  data.frame(
    batch = c("B1", "B1", "B2", "B2", "B2", "B2", "B1", "B1"),
    phase = c("S1", "S1", "S1", "S1", "S1", "S1", "HAB", "HAB"),
    Code = c("A001", "A002", "A003", "A004", "A005", "A006", "A001", "A002"),
    novelLoc = c("R", "L", "R", "L", "R", "R", "R", "L"),
    contactLeft = c(30, 8, 12, 5, 0, 10, 5, 6),
    contactRight = c(10, 24, 12, 15, 0, 10, 5, 2),
    validTime = c(600, 600, 590, 600, 500, 600, 600, 600),
    hiddenLeft = c(0, 0, 2, 0, 0, 0, 0, 0),
    hiddenRight = c(0, 4, 0, 0, 0, 0, 0, 0),
    stringsAsFactors = FALSE
  )
}

hidden_as_unobserved <- function() {
  d <- hidden_as_contact()
  d$contactRight[2] <- 20; d$validTime[2] <- 596
  d$contactLeft[3] <- 10; d$validTime[3] <- 588
  d
}

write_socp_tree <- function(root, d) {
  novel_left <- d$novelLoc == "R"                        # "R" maps the novel side to the left cup
  d$file <- paste0(d$Code, "_", d$phase)
  d$ID <- paste0("mouse-", d$Code)
  d$contactNovel <- ifelse(novel_left, d$contactLeft, d$contactRight)
  d$contactFamiliar <- ifelse(novel_left, d$contactRight, d$contactLeft)
  d$hiddenNovel <- ifelse(novel_left, d$hiddenLeft, d$hiddenRight)
  d$hiddenFamiliar <- ifelse(novel_left, d$hiddenRight, d$hiddenLeft)
  columns <- c("file", "ID", "Code", "batch", "phase", "contactNovel", "contactFamiliar", "contactLeft", "contactRight",
               "validTime", "novelLoc", "hiddenLeft", "hiddenRight", "hiddenNovel", "hiddenFamiliar")
  for (part in split(d, paste(d$batch, d$phase))) {
    folder <- file.path(root, part$batch[1], part$phase[1])
    dir.create(file.path(folder, "plots"), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(part[, columns], file.path(folder, "combined_output.csv"), row.names = FALSE)
  }
  root
}

sensitivity_fixture <- function(contact = hidden_as_contact(), unobserved = hidden_as_unobserved()) {
  root <- temp_run_dir("hidden_sensitivity")
  groups <- file.path(root, "groups.csv")
  dir.create(root, recursive = TRUE)
  writeLines(c("Code,group", "A001,RES", "A002,SUS", "A003,SUS", "A004,RES", "A005,SUS", "A999,CON"), groups)
  list(contact = write_socp_tree(file.path(root, "contact"), contact),
       unobserved = write_socp_tree(file.path(root, "unobserved"), unobserved),
       groups = groups, out = file.path(root, "out"))
}

test_that("the sensitivity compares hidden time as contact and as unobserved, recording by recording and by group", {
  f <- sensitivity_fixture()
  output <- run_r_script(sensitivity_script, c(f$contact, f$unobserved, f$out, "--groups", f$groups))
  expect_script_succeeded(output, "35_hidden_sensitivity.R")
  expect_true(any(grepl("left out of the group tables\\): A006", output)))

  rec <- utils::read.csv(file.path(f$out, "recordings.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(rec), 8 * 8)
  value <- function(d, file, measure, column) d[[column]][d$file == file & d$measure == measure]
  expect_equal(value(rec, "A002_S1", "contactRight", "change"), -4)
  expect_equal(value(rec, "A002_S1", "contactNovel", "as_unobserved"), 20)
  expect_equal(value(rec, "A002_S1", "validTime", "change"), -4)
  expect_equal(value(rec, "A002_S1", "hiddenRight", "as_unobserved"), 4)
  expect_equal(value(rec, "A002_S1", "pref_index", "as_contact"), 0.5)
  expect_equal(value(rec, "A002_S1", "pref_index", "as_unobserved"), 12 / 28)
  expect_equal(value(rec, "A003_S1", "pref_index", "as_unobserved"), -2 / 22)
  expect_true(is.na(value(rec, "A005_S1", "pref_index", "as_contact")))   # no contact: no index
  expect_true(all(rec$change[!rec$file %in% c("A002_S1", "A003_S1")] %in% c(0, NA)))

  groups <- utils::read.csv(file.path(f$out, "groups.csv"), stringsAsFactors = FALSE)
  g <- function(phase, measure, group) groups[groups$phase == phase & groups$measure == measure & groups$group == group, ]
  sus <- g("S1", "pref_index", "SUS")
  expect_equal(c(sus$n_as_contact, sus$median_as_contact, sus$median_as_unobserved), c(2, 0.25, (12 / 28 - 2 / 22) / 2))
  expect_equal(g("S1", "validTime", "SUS")$median_as_unobserved, 588)
  expect_equal(g("S1", "contactRight", "RES")$median_as_contact, 12.5)
  expect_equal(g("HAB", "pref_index", "SUS")$n_as_contact, 1)
  expect_false("CON" %in% groups$group)                  # no CON animal is in the trees
  expect_false("A006" %in% groups$group)

  diffs <- utils::read.csv(file.path(f$out, "group_differences.csv"), stringsAsFactors = FALSE)
  expect_equal(unique(diffs$comparison), "RES - SUS")
  pref <- diffs[diffs$phase == "S1" & diffs$measure == "pref_index", ]
  expect_equal(c(pref$as_contact, pref$as_unobserved), c(0.5 - 0.25, 0.5 - (12 / 28 - 2 / 22) / 2))
  expect_equal(diffs$as_contact[diffs$phase == "HAB" & diffs$measure == "pref_index"], 0 - (-0.5))
})

test_that("the sensitivity writes only the recordings without --groups", {
  f <- sensitivity_fixture()
  expect_script_succeeded(run_r_script(sensitivity_script, c(f$contact, f$unobserved, f$out)), "35_hidden_sensitivity.R")
  expect_true(file.exists(file.path(f$out, "recordings.csv")))
  expect_false(file.exists(file.path(f$out, "groups.csv")))
})

test_that("the sensitivity stops on trees that do not match", {
  fails_with <- function(output, pattern) {
    expect_false(is.null(attr(output, "status")))
    expect_true(any(grepl(pattern, output)), info = paste(output, collapse = "\n"))
  }
  missing <- hidden_as_unobserved()
  f <- sensitivity_fixture(unobserved = missing[missing$Code != "A004", ])
  fails_with(run_r_script(sensitivity_script, c(f$contact, f$unobserved, f$out)), "only in the contact tree: B2/S1/A004_S1")

  other <- hidden_as_unobserved()
  other$hiddenRight[2] <- 3
  f <- sensitivity_fixture(unobserved = other)
  fails_with(run_r_script(sensitivity_script, c(f$contact, f$unobserved, f$out)), "hidden time differs.*B1/S1/A002_S1")

  f <- sensitivity_fixture()
  fails_with(run_r_script(sensitivity_script, c(f$unobserved, f$contact, f$out)), "trees swapped")
  fails_with(run_r_script(sensitivity_script, c(f$contact, f$unobserved, f$out, "--groups")), "usage")
  expect_false(dir.exists(f$out))
})
