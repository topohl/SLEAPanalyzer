# Grouped resampling and cluster-agreement measures for motif analyses.
#
# The recurring mistake in window-based motif analysis is treating windows as
# independent. Windows from one animal share its body size, its tracking
# idiosyncrasies and its behavioural style, so a model validated on windows
# from animals it trained on is measuring memorisation, not generalisation.
# Everything here splits and summarises at the level of the animal.

#' Split rows into training and validation sets by group.
#'
#' Whole groups go to one side or the other, so no animal appears in both. A
#' random split over rows would put windows from the same animal in training
#' and validation, making the validation loss optimistically biased and causing
#' early stopping to select an over-fitted model.
#'
#' @param groups a vector identifying the animal (or batch) of each row
#' @param validation_fraction target fraction of *groups* held out
#' @param seed random seed, recorded so the split is reproducible
#' @return a list of row indices and the group assignment
group_holdout_split <- function(groups, validation_fraction = 0.2, seed = 1) {
  if (length(groups) == 0) stop("groups must not be empty")
  validate_scalar_number(validation_fraction, "validation_fraction", positive = TRUE)
  if (validation_fraction >= 1) stop("validation_fraction must be below 1")
  validate_scalar_number(seed, "seed")

  unique_groups <- unique(groups)
  n_groups <- length(unique_groups)
  if (n_groups < 2) {
    stop(
      "A held-out validation set needs at least two groups; got ", n_groups,
      ". With one animal there is nothing to generalise to."
    )
  }
  n_validation <- max(1L, min(n_groups - 1L, as.integer(round(n_groups * validation_fraction))))

  original_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else NULL
  on.exit({
    if (!is.null(original_seed)) assign(".Random.seed", original_seed, envir = globalenv())
  }, add = TRUE)

  set.seed(seed)
  validation_groups <- sample(unique_groups, n_validation)
  in_validation <- groups %in% validation_groups

  list(
    train_index = which(!in_validation),
    validation_index = which(in_validation),
    train_groups = setdiff(unique_groups, validation_groups),
    validation_groups = validation_groups,
    n_groups = n_groups,
    n_validation_groups = n_validation,
    seed = seed
  )
}

#' Assign groups to k folds, keeping each group whole.
group_kfold <- function(groups, k = 5, seed = 1) {
  validate_positive_integer(k, "k")
  unique_groups <- unique(groups)
  if (length(unique_groups) < k) {
    stop("cannot build ", k, " folds from ", length(unique_groups), " group(s)")
  }
  original_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else NULL
  on.exit({
    if (!is.null(original_seed)) assign(".Random.seed", original_seed, envir = globalenv())
  }, add = TRUE)

  set.seed(seed)
  shuffled <- sample(unique_groups)
  assignment <- rep(seq_len(k), length.out = length(shuffled))
  fold_of_group <- stats::setNames(assignment, shuffled)
  lapply(seq_len(k), function(fold) which(fold_of_group[as.character(groups)] == fold))
}

# ---------------------------------------------------------------------------
# Cluster agreement
# ---------------------------------------------------------------------------

#' Adjusted Rand Index between two clusterings.
#'
#' Measures agreement corrected for chance: 1 is identical, 0 is what random
#' labelling would give. Use it to check that a motif solution survives a
#' change of seed, window length or subsample. A solution that does not is a
#' property of the run, not of the behaviour.
adjusted_rand_index <- function(a, b) {
  if (length(a) != length(b)) stop("clusterings must have the same length")
  if (length(a) < 2) stop("need at least two observations")
  contingency <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2

  sum_cells <- sum(choose2(contingency))
  sum_rows <- sum(choose2(rowSums(contingency)))
  sum_cols <- sum(choose2(colSums(contingency)))
  total <- choose2(length(a))

  expected <- sum_rows * sum_cols / total
  maximum <- (sum_rows + sum_cols) / 2
  if (isTRUE(all.equal(maximum, expected))) return(1)
  unname((sum_cells - expected) / (maximum - expected))
}

#' Normalised mutual information between two clusterings.
normalized_mutual_information <- function(a, b) {
  if (length(a) != length(b)) stop("clusterings must have the same length")
  n <- length(a)
  if (n < 2) stop("need at least two observations")
  joint <- table(a, b) / n
  pa <- rowSums(joint)
  pb <- colSums(joint)

  entropy <- function(p) {
    p <- p[p > 0]
    -sum(p * log(p))
  }
  mutual <- 0
  for (i in seq_along(pa)) {
    for (j in seq_along(pb)) {
      if (joint[i, j] > 0) {
        mutual <- mutual + joint[i, j] * log(joint[i, j] / (pa[i] * pb[j]))
      }
    }
  }
  denominator <- sqrt(entropy(pa) * entropy(pb))
  if (denominator == 0) return(1)
  unname(mutual / denominator)
}

#' Quantify how stable a clustering is across random seeds.
#'
#' @param cluster_fn a function taking a seed and returning cluster labels
#' @param seeds the seeds to compare
#' @return a data frame of pairwise agreement, plus the summary statistics
cluster_seed_stability <- function(cluster_fn, seeds = 1:5) {
  if (!is.function(cluster_fn)) stop("cluster_fn must be a function")
  if (length(seeds) < 2) stop("need at least two seeds to compare")
  labelings <- lapply(seeds, cluster_fn)

  pairs <- utils::combn(seq_along(seeds), 2)
  rows <- lapply(seq_len(ncol(pairs)), function(i) {
    first <- pairs[1, i]
    second <- pairs[2, i]
    data.frame(
      seed_a = seeds[first],
      seed_b = seeds[second],
      adjusted_rand_index = adjusted_rand_index(labelings[[first]], labelings[[second]]),
      normalized_mutual_information = normalized_mutual_information(
        labelings[[first]], labelings[[second]]
      ),
      stringsAsFactors = FALSE
    )
  })
  agreement <- do.call(rbind, rows)
  list(
    pairwise = agreement,
    mean_ari = mean(agreement$adjusted_rand_index),
    min_ari = min(agreement$adjusted_rand_index),
    mean_nmi = mean(agreement$normalized_mutual_information),
    n_seeds = length(seeds)
  )
}

#' Summarise window-level labels to one row per animal.
#'
#' Downstream inference must use these, not the windows. A 10-minute session
#' split into 2-second windows yields 300 rows per animal; treating them as
#' independent observations inflates the apparent sample size 300-fold.
#'
#' @param animal a vector identifying the animal of each window
#' @param labels the cluster label of each window
#' @return a data frame with one row per animal and one column per motif
animal_motif_fractions <- function(animal, labels) {
  if (length(animal) != length(labels)) {
    stop("animal and labels must have the same length")
  }
  counts <- table(animal, labels)
  fractions <- counts / rowSums(counts)
  out <- as.data.frame.matrix(fractions)
  names(out) <- paste0("motif_fraction_", names(out))
  cbind(
    data.frame(
      animal = rownames(fractions),
      n_windows = as.integer(rowSums(counts)),
      stringsAsFactors = FALSE
    ),
    out
  )
}
