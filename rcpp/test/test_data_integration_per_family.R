source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

TOL <- 1e-12

test_fjct <- function() {
  # ============================================================
  # Test fjct_compute_jsd + contribution scores
  # ============================================================

  n_reps_S1 <- 3L
  n_reps_S2 <- 4L
  n_neighbors <- 5L
  n_points <- 3L
  k_families <- 2L

  n_bins <- 4L
  n_genes_S1 <- 10L
  n_genes_S2 <- 10L

  gene_to_family_S1 <- c(1L,0L,0L,0L,0L,0L,0L,0L,0L,2L)
  gene_to_family_S2 <- c(1L,0L,0L,0L,0L,0L,0L,0L,3L,0L)

  neighborhood_genes_S1 <- matrix(0L, n_neighbors, n_points)
  neighborhood_genes_S2 <- matrix(0L, n_neighbors, n_points)

  neighborhood_residuals_S1 <- array(0, dim = c(n_reps_S1, n_neighbors, n_points))
  neighborhood_residuals_S2 <- array(0, dim = c(n_reps_S2, n_neighbors, n_points))

  # ============================================================
  # Test 1 — Simple symmetric case
  # ============================================================

  shared_residual_range <- 2.0

  # S1 residuals
  neighborhood_residuals_S1[,1,1] <- c(-2.0, -0.5, 0.2)
  neighborhood_residuals_S1[,1,3] <- c(2.5, -3.0, 1.2)

  neighborhood_genes_S1[,1] <- c(1L,2L,3L,4L,5L)
  neighborhood_genes_S1[,2] <- c(2L,3L,4L,5L,6L)
  neighborhood_genes_S1[,3] <- c(10L,2L,3L,4L,5L)

  # S2 residuals
  neighborhood_residuals_S2[,1,1] <- c(-2.0, NaN, -0.5, 0.2)
  neighborhood_residuals_S2[,3,3] <- c(0.4, -0.1, 0.0, NaN)

  neighborhood_genes_S2[,1] <- c(1L,2L,3L,4L,5L)
  neighborhood_genes_S2[,2] <- c(2L,3L,4L,5L,6L)
  neighborhood_genes_S2[,3] <- c(9L,2L,3L,4L,5L)

  # ============================================================
  # Family 1
  # ============================================================

  family_idx <- 1L

  mask_S1 <- matrix(FALSE, n_neighbors, n_points)
  mask_S2 <- matrix(FALSE, n_neighbors, n_points)
  mask_S1[1,1] <- TRUE
  mask_S2[1,1] <- TRUE

  expected_included_n_reps_S1 <- c(3L,0L,0L)
  expected_included_n_reps_S2 <- c(3L,0L,0L)
  expected_total_included_n_reps <- 6L
  expected_js_divergences <- c(0,0,0)
  expected_weights <- c(1,0,0)

  # ---- alloc variant ----
  res_alloc <- tox_fjct_compute_jsd(
    family_idx,
    gene_to_family_S1,
    gene_to_family_S2,
    neighborhood_residuals_S1,
    neighborhood_residuals_S2,
    neighborhood_genes_S1,
    neighborhood_genes_S2,
    n_bins,
    shared_residual_range
  )

  res_expert <- tox_fjct_compute_jsd_expert(
    neighborhood_residuals_S1,
    neighborhood_residuals_S2,
    mask_S1,
    mask_S2,
    n_bins,
    shared_residual_range
  )

  # ---- compare alloc vs expert ----
  assert_equal_int(res_expert$included_n_reps_S1, res_alloc$included_n_reps_S1,
                         "Family 1: expert vs alloc included_n_reps_S1 mismatch")
  assert_equal_int(res_expert$included_n_reps_S2, res_alloc$included_n_reps_S2,
                         "Family 1: expert vs alloc included_n_reps_S2 mismatch")
  assert_equal_numeric(res_expert$js_divergences, res_alloc$js_divergences, TOL,
                          "Family 1: expert vs alloc js_divergences mismatch")
  assert_equal_numeric(res_expert$weights, res_alloc$weights, TOL,
                          "Family 1: expert vs alloc weights mismatch")
  assert_equal_numeric(res_expert$global_js_divergence, res_alloc$global_js_divergence, TOL,
                    "Family 1: expert vs alloc global_js_divergence mismatch")

  # ---- compare alloc vs expected ----
  assert_equal_int(res_alloc$total_included_n_reps, expected_total_included_n_reps,
                   "Family 1: total_included_n_reps mismatch")
  assert_equal_int(res_alloc$included_n_reps_S1, expected_included_n_reps_S1,
                         "Family 1: included_n_reps_S1 mismatch")
  assert_equal_int(res_alloc$included_n_reps_S2, expected_included_n_reps_S2,
                         "Family 1: included_n_reps_S2 mismatch")
  assert_equal_numeric(res_alloc$js_divergences, expected_js_divergences, TOL,
                          "Family 1: js_divergences mismatch")
  assert_equal_numeric(res_alloc$weights, expected_weights, TOL,
                          "Family 1: weights mismatch")

  expected_global_jsd <- expected_weights[3] * expected_js_divergences[3]
  assert_equal_numeric(res_alloc$global_js_divergence, expected_global_jsd, TOL,
                    "Family 1: global_js_divergence mismatch")
  # ============================================================
  # Family 2
  # ============================================================

  family_idx <- 2L

  mask_S1 <- matrix(FALSE, n_neighbors, n_points)
  mask_S2 <- matrix(FALSE, n_neighbors, n_points)
  mask_S1[1,3] <- TRUE

  expected_included_n_reps_S1 <- c(0L,0L,3L)
  expected_included_n_reps_S2 <- c(0L,0L,0L)
  expected_total_included_n_reps <- 3L

  expected_js_divergences <- c(
    0.0,
    0.0,
    0.5 * log(2.0)
  )

  expected_weights <- c(0,0,1)

  # ---- alloc variant ----
  res_alloc <- tox_fjct_compute_jsd(
    family_idx,
    gene_to_family_S1,
    gene_to_family_S2,
    neighborhood_residuals_S1,
    neighborhood_residuals_S2,
    neighborhood_genes_S1,
    neighborhood_genes_S2,
    n_bins,
    shared_residual_range
  )

  res_expert <- tox_fjct_compute_jsd_expert(
    neighborhood_residuals_S1,
    neighborhood_residuals_S2,
    mask_S1,
    mask_S2,
    n_bins,
    shared_residual_range
  )

  # ---- compare alloc vs expert ----
  assert_equal_int(res_expert$included_n_reps_S1, res_alloc$included_n_reps_S1,
                         "Family 2: expert vs alloc included_n_reps_S1 mismatch")
  assert_equal_int(res_expert$included_n_reps_S2, res_alloc$included_n_reps_S2,
                         "Family 2: expert vs alloc included_n_reps_S2 mismatch")
  assert_equal_numeric(res_expert$js_divergences, res_alloc$js_divergences, TOL,
                          "Family 2: expert vs alloc js_divergences mismatch")
  assert_equal_numeric(res_expert$weights, res_alloc$weights, TOL,
                          "Family 2: expert vs alloc weights mismatch")
  assert_equal_numeric(res_expert$global_js_divergence, res_alloc$global_js_divergence, TOL,
                    "Family 2: expert vs alloc global_js_divergence mismatch")

  # ---- compare alloc vs expected ----
  assert_equal_int(res_alloc$total_included_n_reps, expected_total_included_n_reps,
                   "Family 2: total_included_n_reps mismatch")
  assert_equal_int(res_alloc$included_n_reps_S1, expected_included_n_reps_S1,
                         "Family 2: included_n_reps_S1 mismatch")
  assert_equal_int(res_alloc$included_n_reps_S2, expected_included_n_reps_S2,
                         "Family 2: included_n_reps_S2 mismatch")
  assert_equal_numeric(res_alloc$js_divergences, expected_js_divergences, TOL,
                          "Family 2: js_divergences mismatch")
  assert_equal_numeric(res_alloc$weights, expected_weights, TOL,
                          "Family 2: weights mismatch")

  expected_global_jsd <- expected_weights[3] * expected_js_divergences[3]
  assert_equal_numeric(res_alloc$global_js_divergence, expected_global_jsd, TOL,
                    "Family 2: global_js_divergence mismatch")

  # ============================================================
  # Contribution scores
  # ============================================================

  total_included_n_reps <- c(6L, 3L)
  global_js_divergence <- c(0.0, 0.5 * log(2.0))

  expected_support_weights <- c(6/9, 3/9)
  expected_contribution_scores <- expected_support_weights * global_js_divergence

  res2 <- tox_fjct_compute_contribution_scores(
    global_js_divergence,
    total_included_n_reps
  )

  assert_equal_numeric(res2$support_weights, expected_support_weights, TOL,
                          "Contribution scores: support_weights mismatch")
  assert_equal_numeric(res2$contribution_scores, expected_contribution_scores, TOL,
                          "Contribution scores: contribution_scores mismatch")
}

run_all_tests()
