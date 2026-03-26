# test_tox_jensen_shannon_test.R
source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

test_compute_gene_means_basic <- function() {
  n_genes <- 4
  n_reps  <- 3

  expr <- matrix(c(
    1, 2, 3,      # gene 1 → mean 2
    4, 5, 6,      # gene 2 → mean 5
    10, 20, 30,   # gene 3 → mean 20
    0, 0, 0       # gene 4 → mean 0
  ), nrow = n_reps, ncol = n_genes, byrow = FALSE)

  expected_means <- c(2, 5, 20, 0)

  out <- compute_gene_means(list(expr, expr))

  assert_equal_int(out$ierr, 0L,
                   "test_compute_gene_means_basic: should succeed")

  assert_equal_numeric(out$means[,1], expected_means,
                       msg = "test_compute_gene_means_basic: S1 means")

  assert_equal_numeric(out$means[,2], expected_means,
                       msg = "test_compute_gene_means_basic: S2 means")
}

test_compute_gene_means_with_nan <- function() {
  n_genes <- 3
  n_reps  <- 4

  expr <- matrix(NA_real_, nrow = n_reps, ncol = n_genes)

  expr[,1] <- c(1, 2, NaN, 3)      # mean = (1+2+3)/3 = 2
  expr[,2] <- c(NaN, 5, 7, 9)      # mean = (5+7+9)/3 = 7
  expr[,3] <- c(10, 20, 30, 40)    # mean = 25

  out <- compute_gene_means(list(expr, expr))

  assert_equal_int(out$ierr, 0L,
                   "test_compute_gene_means_with_nan: should succeed")
  assert_equal_numeric(out$means[1,1], 2,
                       msg = "test_compute_gene_means_with_nan: S1 gene 1 mean")
  assert_equal_numeric(out$means[2,1], 7,
                       msg = "test_compute_gene_means_with_nan: S1 gene 2 mean")
  assert_equal_numeric(out$means[3,1], 25,
                       msg = "test_compute_gene_means_with_nan: S1 gene 3 mean")

  assert_equal_numeric(out$means[1,2], 2,
                       msg = "test_compute_gene_means_with_nan: S2 gene 1 mean")
  assert_equal_numeric(out$means[2,2], 7,
                       msg = "test_compute_gene_means_with_nan: S2 gene 2 mean")
  assert_equal_numeric(out$means[3,2], 25,
                       msg = "test_compute_gene_means_with_nan: S2 gene 3 mean")
}

test_compute_gene_means_all_nan <- function() {
  n_genes <- 2
  n_reps  <- 3

  expr <- matrix(NA_real_, nrow = n_reps, ncol = n_genes)

  expr[,1] <- c(1, 2, 3)   # mean = 2
  expr[,2] <- c(NaN, NaN, NaN)

  out <- compute_gene_means(list(expr))

  assert_equal_int(out$ierr, 0L,
                   "test_compute_gene_means_all_nan: should succeed")

  assert_equal_numeric(out$means[1,1], 2,
                       msg = "test_compute_gene_means_all_nan: gene 1 mean")

  assert_true(is.nan(out$means[2,1]),
              "test_compute_gene_means_all_nan: gene 2 should be NaN")
}

test_compute_residuals_basic <- function() {
  n_genes <- 4
  n_reps  <- 3

  expr <- matrix(c(
    1, 2, 3,      # gene 1
    4, 5, 6,      # gene 2
    10, 20, 30,   # gene 3
    0, 0, 0       # gene 4
  ), nrow = n_reps, ncol = n_genes, byrow = FALSE)

  means <- c(2, 5, 20, 0, 2, 5, 20, 0)

  expected_resid <- matrix(c(
    -1, 0, 1,       # gene 1
    -1, 0, 1,       # gene 2
    -10, 0, 10,     # gene 3
    0, 0, 0         # gene 4
  ), nrow = n_reps, ncol = n_genes, byrow = FALSE)

  out <- compute_residuals(
    expr_list = list(expr, expr),
    means     = matrix(means, nrow = n_genes, ncol = 2),
    max_n_reps_all_studies = n_reps
  )

  assert_equal_int(out$ierr, 0L,
                   "test_compute_residuals_basic: should succeed")

  resid <- out$residuals[,,1]   # extract study 1

  assert_equal_numeric(as.vector(out$residuals[,,1]),
                       as.vector(expected_resid),
                       msg = "test_compute_residuals_basic: residuals 1")

  assert_equal_numeric(as.vector(out$residuals[,,2]),
                       as.vector(expected_resid),
                       msg = "test_compute_residuals_basic: residuals 2")
}

test_compute_residuals_with_nan <- function() {
  n_genes <- 2
  n_reps  <- 4

  expr <- matrix(NA_real_, nrow = n_reps, ncol = n_genes)

  expr[,1] <- c(1, 2, NaN, 3)
  expr[,2] <- c(NaN, 5, 7, 9)

  means <- c(2, 7)

  out <- compute_residuals(
    expr_list = list(expr),
    means     = matrix(means, nrow = n_genes, ncol = 1),
    max_n_reps_all_studies = n_reps
  )

  assert_equal_int(out$ierr, 0L,
                   "test_compute_residuals_with_nan: should succeed")

  resid <- out$residuals[,,1]

  # gene 1
  assert_equal_numeric(resid[1,1], -1, msg = "resid(1,1)")
  assert_equal_numeric(resid[2,1],  0, msg = "resid(2,1)")
  assert_true(is.nan(resid[3,1]),        "resid(3,1) should be NaN")
  assert_equal_numeric(resid[4,1],  1, msg = "resid(4,1)")

  # gene 2
  assert_true(is.nan(resid[1,2]),        "resid(1,2) should be NaN")
  assert_equal_numeric(resid[2,2], -2, msg = "resid(2,2)")
  assert_equal_numeric(resid[3,2],  0, msg = "resid(3,2)")
  assert_equal_numeric(resid[4,2],  2, msg = "resid(4,2)")
}

test_compute_residuals_all_nan <- function() {
  n_genes <- 2
  n_reps  <- 3

  expr <- matrix(NA_real_, nrow = n_reps, ncol = n_genes)

  expr[,1] <- c(1, 2, 3)
  expr[,2] <- c(NaN, NaN, NaN)

  means <- c(2, 0)   # second mean irrelevant

  out <- compute_residuals(
    expr_list = list(expr),
    means     = matrix(means, nrow = n_genes, ncol = 1),
    max_n_reps_all_studies = n_reps
  )

  assert_equal_int(out$ierr, 0L,
                   "test_compute_residuals_all_nan: should succeed")

  resid <- out$residuals[,,1]

  assert_true(all(is.nan(resid[,2])),
              "all residuals for NaN gene should be NaN")
}

test_compute_residuals_invalid_input <- function() {
  expr  <- matrix(1, nrow = 3, ncol = 1)
  means <- matrix(1, nrow = 1, ncol = 1)

  out <- compute_residuals(
    expr_list = list(expr),
    means     = means,
    max_n_reps_all_studies = 3
  )

  # zero genes → simulate by passing empty means
  assert_fails(function () compute_residuals(
    expr_list = list(expr),
    means     = matrix(numeric(0), nrow = 0, ncol = 1),
    max_n_reps_all_studies = 3
  ), "test_compute_residuals_invalid_input: empty input should fail")
}

test_determine_js_comp_test_n_points_n_neighbors <- function() {

  n_studies <- 3
  max_n_genes <- 5
  max_n_reps  <- 4

  join_methods <- c("min", "max", "median")

  # constant data
  R <- 3.0
  residuals <- array(R, dim = c(max_n_reps, max_n_genes, n_studies))
  gene_means <- matrix(R, nrow = max_n_genes, ncol = n_studies)

  expected_n_points    <- 200L # is minimum value for n_points for such few genes
  expected_n_neighbors <- 1L # minimum value for n_neighbors for such few genes

  for (join_method in join_methods) {
    out <- determine_js_comp_test_n_points_n_neighbors(
      residuals = residuals,
      gene_means = gene_means,
      n_bootstraps = 10,
      join_method = join_method,
      min_count_per_mean_bin = 0,
      min_neighbor_overlap = 0.0,
      succeeding_ci_overlap = 0.0,
      random_seed = 42,
      residual_range_quantile = 95
    )

    assert_equal_int(out$n_points, expected_n_points,
      "Case 1: no fallback: expected n_points")

    assert_equal_int(out$n_neighbors, expected_n_neighbors,
      "Case 1: no fallback: expected n_neighbors")

    # best CI = [1,1]    only one candidate, but with succeeding overlap zero -> plateau
    ci <- out$best_candidate_pair_confidence_interval
    assert_equal_numeric(ci[1,], rep(0, n_studies),
      msg = "Case 1: no fallback: CI lower mismatch")
    assert_equal_numeric(ci[2,], rep(0, n_studies),
      msg = "Case 1: no fallback: CI upper mismatch")

    out <- determine_js_comp_test_n_points_n_neighbors(
      residuals = residuals,
      gene_means = gene_means,
      n_bootstraps = 10,
      join_method = join_method,
      min_count_per_mean_bin = 0,
      min_neighbor_overlap = 0.0,
      succeeding_ci_overlap = 1.0,   # not zero, should fail, but only one candidate -> plateau anyway
      random_seed = 42,
      residual_range_quantile = 95
    )

    assert_equal_int(out$n_points, expected_n_points,
      "Case 2: high CI overlap: n_points")

    assert_equal_int(out$n_neighbors, expected_n_neighbors,
      "Case 2: high CI overlap: n_neighbors")

    ci <- out$best_candidate_pair_confidence_interval
    assert_equal_numeric(ci[1,], rep(0, n_studies),
      msg = "Case 2: high CI overlap: CI lower")
    assert_equal_numeric(ci[2,], rep(0, n_studies),
      msg = "Case 2: high CI overlap: CI upper")

    huge_min_count <- max_n_reps * max_n_genes * n_studies + 1

    out <- determine_js_comp_test_n_points_n_neighbors(
      residuals = residuals,
      gene_means = gene_means,
      n_bootstraps = 10,
      join_method = join_method,
      min_count_per_mean_bin = huge_min_count, # should fail, but only one candidate -> plateau anyway
      min_neighbor_overlap = 0.0,
      succeeding_ci_overlap = 0.0,
      random_seed = 42,
      residual_range_quantile = 95
    )

    assert_equal_int(out$n_points, expected_n_points,
      "Case 3: high pmf min count: n_points")

    assert_equal_int(out$n_neighbors, expected_n_neighbors,
      "Case 3: high pmf min count: n_neighbors")

    # CI is -1, as the bootstrapping isn't performed if the min count is not sufficient
    ci <- out$best_candidate_pair_confidence_interval
    assert_equal_numeric(ci[1,], rep(-1, n_studies),
      msg = "Case 3: high pmf min count: CI lower")
    assert_equal_numeric(ci[2,], rep(-1, n_studies),
      msg = "Case 3: high pmf min count: CI upper")

    out <- determine_js_comp_test_n_points_n_neighbors(
      residuals = residuals,
      gene_means = gene_means,
      n_bootstraps = 10,
      join_method = join_method,
      min_count_per_mean_bin = 0,
      min_neighbor_overlap = 0.1, # too high for only one neighbor, but as all gene means are equal, it succeeds anyway
      succeeding_ci_overlap = 0.9,
      random_seed = 42,
      residual_range_quantile = 95
    )

    assert_equal_int(out$n_points, expected_n_points,
      "Case 4: high neighbor overlap: n_points")

    assert_equal_int(out$n_neighbors, expected_n_neighbors,
      "Case 4: high neighbor overlap: n_neighbors")

    ci <- out$best_candidate_pair_confidence_interval
    assert_equal_numeric(ci[1,], rep(0, n_studies),
      msg = "Case 4: high neighbor overlap: CI lower")
    assert_equal_numeric(ci[2,], rep(0, n_studies),
      msg = "Case 4: high neighbor overlap: CI upper")
  }
}

test_js_comp_test <- function() {

  n_studies   <- 2L
  max_n_genes <- 4L
  max_n_reps  <- 3L
  n_points    <- 2L
  n_neighbors <- 2L
  n_bins      <- 3L
  n_perms     <- 5L

  # ------------------------------------------------------------
  # Synthetic deterministic data (same as Fortran)
  # ------------------------------------------------------------
  gene_means <- matrix(c(
    1.0, 2.0, 3.0, 4.0,
    1.5, 2.5, 3.5, 4.5
  ), nrow = max_n_genes, ncol = n_studies, byrow = FALSE)

  residuals <- array(0.1, dim = c(max_n_reps, max_n_genes, n_studies))

  x_star <- c(1.5, 3.5)
  shared_residual_range <- 1.0
  seed <- 314159256L

  # ------------------------------------------------------------
  # Run the R wrapper
  # ------------------------------------------------------------
  out <- js_comp_test(
    residuals = residuals,
    gene_means = gene_means,
    shared_residual_range = shared_residual_range,
    n_bins = n_bins,
    n_points = n_points,
    n_neighbors = n_neighbors,
    n_permutations = n_perms,
    random_seed = seed
  )

  # ------------------------------------------------------------
  # Assertions
  # ------------------------------------------------------------

  # 1. p-values must be exactly 1.0, as all residuals are same
  assert_equal_numeric(
    out$p_values,
    rep(1.0, n_studies),
    msg = "test_js_comp_test: p_values must be 1.0"
  )

  # 2. ierr must be zero
  assert_equal_int(
    out$ierr,
    0L,
    msg = "test_js_comp_test: ierr must be zero"
  )

  # 3. Optional: check shapes (not values)
  assert_equal_int(length(out$x_star), n_points,
    "test_js_comp_test: x_star length mismatch")

  assert_equal_int(length(out$global_js_divergence), n_studies,
    "test_js_comp_test: global_js_divergence length mismatch")

  assert_equal_int(length(out$weights), n_points * n_studies,
    "test_js_comp_test: weights length mismatch")

  invisible(TRUE)
}

run_all_tests()