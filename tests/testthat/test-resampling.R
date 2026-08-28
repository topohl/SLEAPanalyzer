# Grouped splitting and cluster-agreement measures for motif analyses.

windows_from <- function(animals, per_animal) {
  rep(animals, each = per_animal)
}

test_that("a holdout split never puts one animal on both sides", {
  groups <- windows_from(paste0("m", 1:10), 30)
  split <- group_holdout_split(groups, validation_fraction = 0.2, seed = 42)

  train_animals <- unique(groups[split$train_index])
  validation_animals <- unique(groups[split$validation_index])

  # This is the property a random split over windows does not have.
  expect_length(intersect(train_animals, validation_animals), 0)
  expect_equal(
    sort(c(train_animals, validation_animals)), sort(unique(groups))
  )
  expect_equal(
    length(split$train_index) + length(split$validation_index), length(groups)
  )
})

test_that("a random split over windows does leak animals, unlike the grouped split", {
  set.seed(1)
  groups <- windows_from(paste0("m", 1:10), 30)
  # Demonstrate the failure mode the grouped split exists to prevent.
  naive_validation <- sample.int(length(groups), size = 60)
  leaked <- intersect(
    unique(groups[naive_validation]), unique(groups[-naive_validation])
  )
  expect_gt(length(leaked), 0)

  split <- group_holdout_split(groups, seed = 1)
  expect_length(
    intersect(unique(groups[split$validation_index]), unique(groups[split$train_index])),
    0
  )
})

test_that("the split holds out roughly the requested fraction of animals", {
  groups <- windows_from(paste0("m", 1:20), 10)
  split <- group_holdout_split(groups, validation_fraction = 0.25, seed = 7)
  expect_equal(split$n_groups, 20L)
  expect_equal(split$n_validation_groups, 5L)
  expect_length(split$validation_groups, 5L)
})

test_that("the split is reproducible from its seed and does not disturb the RNG", {
  groups <- windows_from(paste0("m", 1:12), 5)
  first <- group_holdout_split(groups, seed = 99)
  second <- group_holdout_split(groups, seed = 99)
  expect_equal(first$validation_groups, second$validation_groups)
  expect_false(identical(
    sort(group_holdout_split(groups, seed = 1)$validation_groups),
    sort(group_holdout_split(groups, seed = 2)$validation_groups)
  ))

  # A caller's random stream must not be silently reset.
  set.seed(123)
  before <- runif(1)
  set.seed(123)
  invisible(group_holdout_split(groups, seed = 5))
  expect_equal(runif(1), before)
})

test_that("a single animal cannot be split", {
  expect_error(
    group_holdout_split(rep("m1", 50)),
    "at least two groups"
  )
})

test_that("group k-fold keeps each animal within one fold", {
  groups <- windows_from(paste0("m", 1:10), 6)
  folds <- group_kfold(groups, k = 5, seed = 3)
  expect_length(folds, 5L)
  expect_equal(sum(lengths(folds)), length(groups))

  fold_animals <- lapply(folds, function(idx) unique(groups[idx]))
  for (i in seq_along(fold_animals)) {
    for (j in seq_along(fold_animals)) {
      if (i >= j) next
      expect_length(intersect(fold_animals[[i]], fold_animals[[j]]), 0)
    }
  }
})

test_that("group k-fold refuses more folds than groups", {
  expect_error(group_kfold(windows_from(paste0("m", 1:3), 5), k = 5), "cannot build 5 folds")
})

test_that("the adjusted Rand index is 1 for identical clusterings", {
  labels <- c(1, 1, 2, 2, 3, 3)
  expect_equal(adjusted_rand_index(labels, labels), 1)
  # Relabelling must not matter.
  expect_equal(adjusted_rand_index(labels, c("a", "a", "b", "b", "c", "c")), 1)
})

test_that("the adjusted Rand index is near zero for unrelated clusterings", {
  set.seed(11)
  a <- sample(1:4, 400, replace = TRUE)
  b <- sample(1:4, 400, replace = TRUE)
  expect_lt(abs(adjusted_rand_index(a, b)), 0.05)
})

test_that("the adjusted Rand index falls as clusterings diverge", {
  base <- rep(1:4, each = 25)
  perturbed <- base
  perturbed[sample.int(100, 10)] <- 1L
  strong <- adjusted_rand_index(base, perturbed)

  scrambled <- base
  set.seed(5)
  scrambled[sample.int(100, 70)] <- sample(1:4, 70, replace = TRUE)
  weak <- adjusted_rand_index(base, scrambled)

  expect_gt(strong, weak)
})

test_that("normalized mutual information behaves like ARI at the extremes", {
  labels <- rep(1:3, each = 20)
  expect_equal(normalized_mutual_information(labels, labels), 1)
  set.seed(12)
  random <- sample(1:3, 60, replace = TRUE)
  expect_lt(normalized_mutual_information(labels, random), 0.2)
})

test_that("seed stability reports agreement across runs", {
  # A clustering that ignores its seed is perfectly stable.
  stable <- cluster_seed_stability(function(seed) rep(1:3, each = 10), seeds = 1:4)
  expect_equal(stable$mean_ari, 1)
  expect_equal(stable$min_ari, 1)
  expect_equal(nrow(stable$pairwise), 6L)

  # A clustering that is essentially random is not.
  unstable <- cluster_seed_stability(function(seed) {
    set.seed(seed)
    sample(1:3, 300, replace = TRUE)
  }, seeds = 1:4)
  expect_lt(unstable$mean_ari, 0.1)
})

test_that("cluster_seed_stability validates its arguments", {
  expect_error(cluster_seed_stability("not a function"), "must be a function")
  expect_error(cluster_seed_stability(function(s) 1:5, seeds = 1), "at least two seeds")
})

test_that("window labels are summarised to one row per animal", {
  animal <- c(rep("m1", 10), rep("m2", 10))
  labels <- c(rep("rest", 8), rep("run", 2), rep("rest", 3), rep("run", 7))
  summary_tbl <- animal_motif_fractions(animal, labels)

  # Two animals, not twenty windows: this is what inference must use.
  expect_equal(nrow(summary_tbl), 2L)
  expect_equal(summary_tbl$n_windows, c(10L, 10L))
  expect_equal(summary_tbl$motif_fraction_rest, c(0.8, 0.3))
  expect_equal(summary_tbl$motif_fraction_run, c(0.2, 0.7))
  expect_equal(
    unname(rowSums(summary_tbl[, grep("^motif_fraction_", names(summary_tbl))])),
    c(1, 1)
  )
})

test_that("animal_motif_fractions rejects mismatched lengths", {
  expect_error(animal_motif_fractions(c("a", "b"), "x"), "same length")
})
