# Comprehensive R test suite for the gJCT/fJCT JSD-Comp-Test parameter-search building blocks
# in tensoromics (Issue #126). Fixtures below are translated from the already
# hand-computed/verified numbers in the Fortran test suite
# mod_test_data_integration_js_comp_test.F90, rather than re-derived.
source("r/load_tensor_omics.R")
source("r/test_helpers.R")

TOL <- 1e-12

test_estimate_bin_count <- function() {
  # ============================================================
  # Test 1 -- basic hand-computed (test_estimate_bin_count_basic_hand_computed):
  # Sturges gives 1 + nint(log(1)/LOG_2) = 1. 25th/75th percentiles (rank 3.25 and 7.75)
  # interpolate to -1.75 and 2.75, so the Freedman-Diaconis half-width is
  # (2.75 - (-1.75)) / 1^(1/3) = 4.5; with shared_residual_range=9.0, nint(9.0/4.5)=2.
  # ============================================================
  residuals <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5)

  n_bins <- estimate_bin_count(residuals, max_n_reps_all_studies = 1, n_neighbors = 1,
                                shared_residual_range = 9.0)
  assert_equal_int(n_bins, 2L, "Test 1 failed: expected n_bins=2")

  residuals_perm <- order(residuals)
  n_bins_expert <- estimate_bin_count_expert(residuals, residuals_perm, max_n_reps_all_studies = 1,
                                              n_neighbors = 1, shared_residual_range = 9.0)
  assert_equal_int(n_bins_expert, 2L, "Test 1 (expert) failed: expected n_bins=2")

  # ============================================================
  # Test 2 -- all-NaN pool falls back to a single bin (test_estimate_bin_count_all_nan_gives_one_bin)
  # ============================================================
  residuals_nan <- rep(NaN, 4)
  n_bins_nan <- estimate_bin_count(residuals_nan, max_n_reps_all_studies = 1, n_neighbors = 1,
                                    shared_residual_range = 9.0)
  assert_equal_int(n_bins_nan, 1L, "Test 2 failed: expected n_bins=1")

  # ============================================================
  # Test 3 -- clamped to MAX_N_BINS=256 (test_estimate_bin_count_clamped_to_max_n_bins)
  # ============================================================
  n_bins_clamped <- estimate_bin_count(residuals, max_n_reps_all_studies = 1, n_neighbors = 1,
                                        shared_residual_range = 5000.0)
  assert_equal_int(n_bins_clamped, 256L, "Test 3 failed: expected n_bins=256")
}

test_generate_js_comp_test_candidates <- function() {
  residuals <- c(1, 2, 3, 4, 5)

  # ============================================================
  # Test 1 -- small-N candidate-grid collapse at its exact threshold
  # (test_generate_js_comp_test_candidates_collapses_at_8742): every candidate shares n_points=374.
  # ============================================================
  out <- generate_js_comp_test_candidates(8742L, residuals, max_n_reps_all_studies = 1,
                                           shared_residual_range = 1.0)
  # NOTE: the generated R wrapper truncates `n_bins_candidates` to `n_candidates` correctly (it
  # is a vector), but truncates `candidates_n_points_n_neighbors` via `utils::head(matrix,
  # n_candidates)`, which trims ROWS, not columns -- a no-op here since the matrix always has
  # exactly 2 rows. The columns past `n_candidates` therefore hold uninitialized garbage; slice
  # manually using the (correctly truncated) length of `n_bins_candidates` instead.
  n_candidates <- length(out$n_bins_candidates)
  candidates <- out$candidates_n_points_n_neighbors[, 1:n_candidates, drop = FALSE]
  assert_true(n_candidates >= 1, "Test 1 failed: expected at least one candidate")
  assert_true(all(candidates[1, ] == candidates[1, 1]),
              "Test 1 failed: every candidate must share the same n_points")
  assert_equal_int(as.integer(candidates[1, 1]), 374L, "Test 1 failed: n_points should be 374")

  # ============================================================
  # Test 2 -- the other side of the bracket (test_generate_js_comp_test_candidates_has_two_distinct_at_8743):
  # at least two distinct n_points values, first=375, last=300.
  # ============================================================
  out2 <- generate_js_comp_test_candidates(8743L, residuals, max_n_reps_all_studies = 1,
                                            shared_residual_range = 1.0)
  # See the NOTE above test 1: slice manually using n_bins_candidates' length, not head().
  n_candidates2 <- length(out2$n_bins_candidates)
  candidates2 <- out2$candidates_n_points_n_neighbors[, 1:n_candidates2, drop = FALSE]
  n_distinct <- length(unique(candidates2[1, ]))
  assert_true(n_distinct >= 2, "Test 2 failed: expected at least two distinct n_points values")
  assert_equal_int(as.integer(candidates2[1, 1]), 375L, "Test 2 failed: first n_points should be 375")
  assert_equal_int(as.integer(candidates2[1, n_candidates2]), 300L,
                    "Test 2 failed: last n_points should collapse to 300")

  # ============================================================
  # Test 3 -- validation (test_generate_js_comp_test_candidates_validation)
  # ============================================================
  assert_error(
    generate_js_comp_test_candidates(0L, residuals, max_n_reps_all_studies = 1, shared_residual_range = 1.0),
    "Test 3 failed: expected ERR_INVALID_INPUT for max_n_genes_all_studies=0", ERR_INVALID_INPUT
  )

  residuals_perm <- c(1L, 2L, 3L, 4L, 5L)
  out3 <- generate_js_comp_test_candidates_expert(100L, residuals, residuals_perm, max_n_reps_all_studies = 1,
                                                   shared_residual_range = 1.0)
  assert_true(ncol(out3$candidates_n_points_n_neighbors) >= 1,
              "Test 3 failed: expert tier should accept a valid permutation")
}

test_check_neighborhood_overlaps <- function() {
  # ============================================================
  # Test 1 -- identical ranges always pass (test_neighborhood_overlaps_all_pass)
  # ============================================================
  neighborhood_range <- matrix(c(1, 10, 1, 10, 1, 10), nrow = 2)
  assert_true(check_neighborhood_overlaps(neighborhood_range, min_neighbor_overlap = 1.0),
              "Test 1 failed: identical ranges should always pass")

  # ============================================================
  # Test 2 -- overlap exactly at threshold passes (test_neighborhood_overlaps_exactly_at_threshold_passes):
  # ranges [1,10] and [6,10]: overlap = (10-6)/(10-1) = 4/9.
  # ============================================================
  neighborhood_range_2 <- matrix(c(1, 10, 6, 10), nrow = 2)
  assert_true(check_neighborhood_overlaps(neighborhood_range_2, min_neighbor_overlap = 4.0 / 9.0),
              "Test 2 failed: overlap exactly at threshold must pass")

  # ============================================================
  # Test 3 -- same 4/9 overlap, threshold 0.5 > 4/9, must fail
  # (test_neighborhood_overlaps_below_threshold_fails)
  # ============================================================
  assert_false(check_neighborhood_overlaps(neighborhood_range_2, min_neighbor_overlap = 0.5),
               "Test 3 failed: overlap below threshold must fail")
}

test_check_mean_pmf_min_counts <- function() {
  # Test 1 -- all pass (test_mean_pmf_min_counts_all_pass)
  mean_pmf_counts <- matrix(c(10, 10, 10, 10), nrow = 2)
  assert_true(check_mean_pmf_min_counts(mean_pmf_counts, min_count = 5L),
              "Test 1 failed: all bins above minimum should pass")

  # Test 2 -- exactly at threshold passes (test_mean_pmf_min_counts_exactly_at_threshold_passes)
  mean_pmf_counts_2 <- matrix(c(3, 3, 3, 3), nrow = 2)
  assert_true(check_mean_pmf_min_counts(mean_pmf_counts_2, min_count = 3L),
              "Test 2 failed: count exactly at minimum must pass")

  # Test 3 -- one bin below threshold fails (test_mean_pmf_min_counts_below_threshold_fails)
  mean_pmf_counts_3 <- matrix(c(3, 3, 3, 2), nrow = 2)
  assert_false(check_mean_pmf_min_counts(mean_pmf_counts_3, min_count = 3L),
               "Test 3 failed: one bin below minimum must fail the whole gate")
}

test_check_plateau_condition <- function() {
  # ============================================================
  # Test 1 -- METHOD_JOIN_MIN requires every study to pass (test_check_plateau_condition_join_min_requires_all).
  # 2-of-3 studies passing must fail; 3-of-3 must succeed.
  # ============================================================
  confidence_interval <- matrix(c(0.4, 0.6, 0.4, 0.6, -2.0, -1.0), nrow = 2)
  best_ci <- matrix(c(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), nrow = 2)

  result <- check_plateau_condition(confidence_interval, best_ci, best_candidate_index = 0L,
                                     best_exceeded_ci_overlap_count = 0L, candidate_index = 1L,
                                     join_method = "join_min", succeeding_ci_overlap = 0.9)
  assert_false(result$plateau_found, "Test 1 failed: 2-of-3 studies must not plateau")

  confidence_interval[, 3] <- c(0.4, 0.6)
  best_ci_2 <- matrix(c(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), nrow = 2)
  result2 <- check_plateau_condition(confidence_interval, best_ci_2, best_candidate_index = 0L,
                                      best_exceeded_ci_overlap_count = 0L, candidate_index = 1L,
                                      join_method = "join_min", succeeding_ci_overlap = 0.9)
  assert_true(result2$plateau_found, "Test 1 failed: 3-of-3 studies must plateau")

  # ============================================================
  # Test 2 -- METHOD_JOIN_MAX succeeds once any one study passes
  # (test_check_plateau_condition_join_max_requires_any). 0-of-3 must fail, 1-of-3 must succeed.
  # ============================================================
  confidence_interval_max <- matrix(c(-2.0, -1.0, -2.0, -1.0, -2.0, -1.0), nrow = 2)
  best_ci_max <- matrix(c(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), nrow = 2)
  result_max_0 <- check_plateau_condition(confidence_interval_max, best_ci_max, best_candidate_index = 0L,
                                           best_exceeded_ci_overlap_count = 0L, candidate_index = 1L,
                                           join_method = "join_max", succeeding_ci_overlap = 0.9)
  assert_false(result_max_0$plateau_found, "Test 2 failed: 0-of-3 studies must not plateau")

  confidence_interval_max[, 1] <- c(0.4, 0.6)
  best_ci_max_2 <- matrix(c(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), nrow = 2)
  result_max_1 <- check_plateau_condition(confidence_interval_max, best_ci_max_2, best_candidate_index = 0L,
                                           best_exceeded_ci_overlap_count = 0L, candidate_index = 1L,
                                           join_method = "join_max", succeeding_ci_overlap = 0.9)
  assert_true(result_max_1$plateau_found, "Test 2 failed: 1-of-3 studies must plateau")

  # ============================================================
  # Test 3 -- METHOD_JOIN_MEDIAN succeeds once a majority pass (n_studies=3 -> threshold 2)
  # (test_check_plateau_condition_join_median_requires_majority). 1-of-3 fails, 2-of-3 succeeds.
  # ============================================================
  confidence_interval_med <- matrix(c(0.4, 0.6, -2.0, -1.0, -2.0, -1.0), nrow = 2)
  best_ci_med <- matrix(c(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), nrow = 2)
  result_med_1 <- check_plateau_condition(confidence_interval_med, best_ci_med, best_candidate_index = 0L,
                                           best_exceeded_ci_overlap_count = 0L, candidate_index = 1L,
                                           join_method = "join_median", succeeding_ci_overlap = 0.9)
  assert_false(result_med_1$plateau_found, "Test 3 failed: 1-of-3 must not plateau")

  confidence_interval_med[, 2] <- c(0.4, 0.6)
  best_ci_med_2 <- matrix(c(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), nrow = 2)
  result_med_2 <- check_plateau_condition(confidence_interval_med, best_ci_med_2, best_candidate_index = 0L,
                                           best_exceeded_ci_overlap_count = 0L, candidate_index = 1L,
                                           join_method = "join_median", succeeding_ci_overlap = 0.9)
  assert_true(result_med_2$plateau_found, "Test 3 failed: 2-of-3 must plateau")

  # ============================================================
  # Test 4 -- a candidate worse than the previous best plateaus immediately, and none of the
  # best_* state is overwritten (test_check_plateau_condition_worse_than_previous_short_circuits).
  # ============================================================
  confidence_interval_worse <- matrix(c(-5.0, -4.0, -5.0, -4.0), nrow = 2)
  best_ci_worse <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
  expected_best_ci_worse <- best_ci_worse

  result_worse <- check_plateau_condition(confidence_interval_worse, best_ci_worse, best_candidate_index = 5L,
                                           best_exceeded_ci_overlap_count = 2L, candidate_index = 9L,
                                           join_method = "join_min", succeeding_ci_overlap = 0.5)
  assert_true(result_worse$plateau_found, "Test 4 failed: a worse candidate must plateau immediately")
  assert_equal_int(as.integer(result_worse$best_candidate_index), 5L,
                    "Test 4 failed: best_candidate_index must not be overwritten")
  assert_equal_int(as.integer(result_worse$best_exceeded_ci_overlap_count), 2L,
                    "Test 4 failed: best_exceeded_ci_overlap_count must not be overwritten")
  assert_equal_numeric(best_ci_worse, expected_best_ci_worse, TOL,
                        "Test 4 failed: best_candidate_pair_confidence_interval must not be overwritten")
}

test_create_mean_pmf <- function() {
  # ============================================================
  # Test 1 -- two studies, hand-computed (test_create_mean_pmf_two_studies_hand_computed):
  # mean_pmf = ([0.3,0.7]+[0.5,0.5])/2 = [0.4,0.6], mean_pmf_counts=[8,12], included=20.
  # ============================================================
  pmfs <- array(0, dim = c(2, 1, 2))
  pmfs[, 1, 1] <- c(0.3, 0.7)
  pmfs[, 1, 2] <- c(0.5, 0.5)
  counts <- array(0L, dim = c(2, 1, 2))
  counts[, 1, 1] <- c(3L, 7L)
  counts[, 1, 2] <- c(5L, 5L)
  included_n_reps <- matrix(c(10L, 10L), nrow = 1)

  result <- create_mean_pmf(pmfs, counts, included_n_reps)
  assert_equal_numeric(result$mean_pmf[, 1], c(0.4, 0.6), TOL, "Test 1 failed: mean_pmf mismatch")
  assert_equal_int(as.integer(result$mean_pmf_counts[, 1]), c(8L, 12L), "Test 1 failed: mean_pmf_counts mismatch")
  assert_equal_int(as.integer(result$mean_pmf_included_n_reps[1]), 20L,
                    "Test 1 failed: mean_pmf_included_n_reps mismatch")

  # ============================================================
  # Test 2 -- three studies, DOCUMENTS THE KNOWN LIMITATION: averages over ALL studies including
  # self, not a true leave-one-out background (test_create_mean_pmf_three_studies_includes_self).
  # ============================================================
  pmfs3 <- array(0, dim = c(2, 1, 3))
  pmfs3[, 1, 1] <- c(0.2, 0.8)
  pmfs3[, 1, 2] <- c(0.4, 0.6)
  pmfs3[, 1, 3] <- c(0.9, 0.1)
  counts3 <- array(0L, dim = c(2, 1, 3))
  counts3[, 1, 1] <- c(2L, 8L)
  counts3[, 1, 2] <- c(4L, 6L)
  counts3[, 1, 3] <- c(9L, 1L)
  included_n_reps3 <- matrix(c(10L, 10L, 10L), nrow = 1)

  result3 <- create_mean_pmf(pmfs3, counts3, included_n_reps3)
  assert_equal_numeric(result3$mean_pmf[, 1], c(0.5, 0.5), TOL, "Test 2 failed: mean_pmf mismatch")
  assert_equal_int(as.integer(result3$mean_pmf_counts[, 1]), c(15L, 15L), "Test 2 failed: mean_pmf_counts mismatch")
  assert_equal_int(as.integer(result3$mean_pmf_included_n_reps[1]), 30L,
                    "Test 2 failed: mean_pmf_included_n_reps mismatch")

  leave_one_out_study1 <- (c(0.4, 0.6) + c(0.9, 0.1)) / 2.0
  assert_true(abs(result3$mean_pmf[1, 1] - leave_one_out_study1[1]) > TOL,
              "Test 2 failed: mean_pmf must differ from the leave-one-out background -- known limitation was fixed "
              %+% "without updating this test")

  # ============================================================
  # Test 3 -- validation: pmfs must be in [0,1], counts must be non-negative
  # (test_create_mean_pmf_validation)
  # ============================================================
  pmfs_bad <- array(1.5, dim = c(1, 1, 1))
  counts_bad <- array(1L, dim = c(1, 1, 1))
  included_n_reps_bad <- matrix(1L, nrow = 1)
  assert_error(create_mean_pmf(pmfs_bad, counts_bad, included_n_reps_bad),
               "Test 3 failed: expected ERR_INVALID_INPUT for pmf above 1.0", ERR_INVALID_INPUT)

  pmfs_ok <- array(0.5, dim = c(1, 1, 1))
  counts_neg <- array(-1L, dim = c(1, 1, 1))
  assert_error(create_mean_pmf(pmfs_ok, counts_neg, included_n_reps_bad),
               "Test 3 failed: expected ERR_INVALID_INPUT for negative count", ERR_INVALID_INPUT)
}

# small helper: base R has no built-in string concatenation infix
`%+%` <- function(a, b) paste0(a, b)

test_create_mean_pmf_only <- function() {
  # `create_mean_pmf_only` must produce exactly the same mean_pmf as `create_mean_pmf`
  # (test_create_mean_pmf_only_matches_create_mean_pmf).
  pmfs <- array(0, dim = c(2, 1, 2))
  pmfs[, 1, 1] <- c(0.3, 0.7)
  pmfs[, 1, 2] <- c(0.5, 0.5)
  counts <- array(0L, dim = c(2, 1, 2))
  counts[, 1, 1] <- c(3L, 7L)
  counts[, 1, 2] <- c(5L, 5L)
  included_n_reps <- matrix(c(10L, 10L), nrow = 1)

  full <- create_mean_pmf(pmfs, counts, included_n_reps)
  only <- create_mean_pmf_only(pmfs)

  assert_equal_numeric(only, full$mean_pmf, TOL, "mean_pmf mismatch between create_mean_pmf and create_mean_pmf_only")
}

test_bootstrap_histogram <- function() {
  # ============================================================
  # Test 1 -- same random_seed twice must give bit-for-bit identical confidence_interval
  # (test_bootstrap_histogram_seeded_reproducibility).
  # ============================================================
  mean_pmf_counts <- matrix(c(5, 3, 2, 2, 2, 6), nrow = 3)
  mean_pmf_included_n_reps <- c(10L, 10L)
  included_n_reps <- matrix(c(4L, 4L, 6L, 6L), nrow = 2)

  ci_first <- matrix(c(0.2, 0.2, 0.3, 0.3), nrow = 2)
  ci_first <- bootstrap_histogram(15L, mean_pmf_counts, mean_pmf_included_n_reps, included_n_reps, ci_first,
                                   random_seed = 123L)

  ci_second <- matrix(c(0.2, 0.2, 0.3, 0.3), nrow = 2)
  ci_second <- bootstrap_histogram(15L, mean_pmf_counts, mean_pmf_included_n_reps, included_n_reps, ci_second,
                                    random_seed = 123L)

  assert_equal_numeric(ci_first, ci_second, TOL, "same random_seed twice must give identical confidence_interval")

  # ============================================================
  # Test 2 -- degenerate single-bin histogram: JSD is always exactly 0.0, so the CI must
  # collapse to [0.0, 0.0] (test_bootstrap_histogram_degenerate_single_bin_zero_ci).
  # ============================================================
  mean_pmf_counts_deg <- matrix(10L, nrow = 1, ncol = 1)
  mean_pmf_included_n_reps_deg <- c(10L)
  included_n_reps_deg <- matrix(c(5L, 5L), nrow = 1)
  confidence_interval_deg <- matrix(0.0, nrow = 2, ncol = 2)

  confidence_interval_deg <- bootstrap_histogram(5L, mean_pmf_counts_deg, mean_pmf_included_n_reps_deg,
                                                  included_n_reps_deg, confidence_interval_deg, random_seed = 7L)

  assert_equal_numeric(confidence_interval_deg, matrix(0.0, nrow = 2, ncol = 2), TOL,
                        "a single-bin histogram always has JSD=0, so the CI must collapse to [0,0]")
}

test_calc_js_comp_test_candidate_bounds <- function() {
  # Cross-check against the already-verified candidate-grid collapse fixture: with
  # max_n_genes_all_studies=8742 the grid's first (largest) n_points candidate is 374, and
  # `max_n_points_candidate` is documented as being exactly that value.
  bounds <- calc_js_comp_test_candidate_bounds(8742L)
  assert_equal_int(as.integer(bounds$max_n_points_candidate), 374L, "expected max_n_points_candidate=374")

  out <- generate_js_comp_test_candidates(8742L, c(1, 2, 3, 4, 5), max_n_reps_all_studies = 1,
                                           shared_residual_range = 1.0)
  # See the NOTE in test_generate_js_comp_test_candidates: slice manually using
  # n_bins_candidates' length, since the generated wrapper's own truncation of the matrix is a
  # no-op (utils::head() on a matrix trims rows, not columns).
  n_candidates <- length(out$n_bins_candidates)
  candidates <- out$candidates_n_points_n_neighbors[, 1:n_candidates, drop = FALSE]
  assert_true(bounds$max_n_neighbors_candidate >= max(candidates[2, ]),
              "max_n_neighbors_candidate must be a safe upper bound on the grid's n_neighbors")
}

test_calc_js_comp_test_n_top_k_jsds <- function() {
  # Test 1 -- default 2.5% (test_calc_js_comp_test_n_top_k_jsds_default_hand_computed): n_bootstraps=1000 -> 25.
  assert_equal_int(as.integer(calc_js_comp_test_n_top_k_jsds(1000L)), 25L, "Test 1 failed: expected n_top_k=25")

  # Test 2 -- explicit level (test_calc_js_comp_test_n_top_k_jsds_explicit_significance_level): 200, 5.0% -> 10.
  assert_equal_int(as.integer(calc_js_comp_test_n_top_k_jsds(200L, two_sided_bootstrapping_significance_level = 5.0)),
                    10L, "Test 2 failed: expected n_top_k=10")

  # Test 3 -- clamped to 1 (test_calc_js_comp_test_n_top_k_jsds_clamped_to_at_least_one): 10 -> floor(0.25)=0 -> 1.
  assert_equal_int(as.integer(calc_js_comp_test_n_top_k_jsds(10L)), 1L, "Test 3 failed: expected n_top_k=1")
}

test_run_js_comp_test_two_studies_hand_traceable <- function() {
  # End-to-end, fully closed-form 2-study case (test_run_js_comp_test_two_studies_hand_traceable),
  # with n_permutations=0 so the pipeline is deterministic and hand-traceable.
  log2_3 <- log(3.0) / log(2.0)
  expected_jsd <- 1.5 - 0.75 * log2_3

  gene_means <- matrix(c(1.0, 5.0, 1.0, 5.0), nrow = 2)
  gene_means_perms <- matrix(c(1L, 2L, 1L, 2L), nrow = 2)

  residuals <- array(NaN, dim = c(2, 2, 2))
  residuals[, 1, 1] <- c(-1.0, 1.0)
  residuals[, 1, 2] <- c(-3.0, 3.0)

  x_star <- c(1.0)

  result <- run_js_comp_test(n_neighbors = 1L, n_bins = 4L, shared_residual_range = 4.0, gene_means = gene_means,
                              gene_means_perms = gene_means_perms, residuals = residuals, x_star = x_star,
                              n_permutations = 0L, random_seed = 1L)

  assert_equal_int(as.integer(result$neighborhood_indices[, 1, 1]), 1L, "study 1 neighbor should be gene 1")
  assert_equal_int(as.integer(result$neighborhood_indices[, 1, 2]), 1L, "study 2 neighbor should be gene 1")

  assert_equal_int(as.integer(result$counts[, 1, 1]), c(0L, 1L, 1L, 0L), "study 1 counts mismatch")
  assert_equal_int(as.integer(result$counts[, 1, 2]), c(1L, 0L, 0L, 1L), "study 2 counts mismatch")
  assert_equal_int(as.integer(result$included_n_reps[1, 1]), 2L, "study 1 included_n_reps mismatch")
  assert_equal_int(as.integer(result$included_n_reps[1, 2]), 2L, "study 2 included_n_reps mismatch")

  assert_equal_numeric(result$mean_pmf[, 1], c(0.25, 0.25, 0.25, 0.25), TOL, "mean_pmf should be uniform")
  assert_equal_int(as.integer(result$mean_pmf_counts[, 1]), c(1L, 1L, 1L, 1L), "mean_pmf_counts mismatch")
  assert_equal_int(as.integer(result$mean_pmf_included_n_reps[1]), 4L, "mean_pmf_included_n_reps mismatch")

  assert_true(abs(result$global_js_divergence[1] - expected_jsd) < 1e-9, "study 1 global JSD, closed form")
  assert_true(abs(result$global_js_divergence[2] - expected_jsd) < 1e-9, "study 2 global JSD, closed form")

  assert_true(abs(result$weights[1, 1] - 1.0) < TOL, "single reference point weight is 1.0")
  assert_true(abs(result$weights[1, 2] - 1.0) < TOL, "single reference point weight is 1.0")

  assert_true(abs(result$p_values[1] - 0.0) < TOL, "n_permutations=0 -> p_values stay 0.0")
  assert_true(abs(result$p_values[2] - 0.0) < TOL, "n_permutations=0 -> p_values stay 0.0")
}

test_run_js_comp_test_parameter_search <- function() {
  max_n_genes_all_studies <- 2000L
  max_n_reps_all_studies <- 3L
  n_studies <- 2L

  gene_means <- matrix(0.0, nrow = max_n_genes_all_studies, ncol = n_studies)
  residuals <- array(0.0, dim = c(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
  for (i_study in 1:n_studies) {
    gene_means[, i_study] <- seq_len(max_n_genes_all_studies)
    for (i_gene in 1:max_n_genes_all_studies) {
      residuals[, i_gene, i_study] <- c(-0.5, 0.0, 0.5)
    }
  }

  # ============================================================
  # Test 1 -- no candidate ever plateaus (min_count_per_mean_bin impossibly high): falls back to
  # the FIRST (finest-resolution) candidate, (n_points, n_neighbors)=(300, 3), and resets
  # best_candidate_pair_confidence_interval to -1.0 throughout
  # (test_param_search_no_plateau_falls_back_to_finest).
  # ============================================================
  result <- run_js_comp_test_parameter_search(gene_means, residuals, shared_residual_range = 1.0, n_bootstraps = 10L,
                                               join_method = "join_min", min_count_per_mean_bin = 1000000L,
                                               random_seed = 1L)

  assert_equal_int(as.integer(result$n_points), 300L, "expected n_points=300")
  assert_equal_int(as.integer(result$n_neighbors), 3L, "expected n_neighbors=3")
  assert_equal_numeric(result$best_candidate_pair_confidence_interval[, 1], c(-1.0, -1.0), TOL, "study 1 CI reset")
  assert_equal_numeric(result$best_candidate_pair_confidence_interval[, 2], c(-1.0, -1.0), TOL, "study 2 CI reset")

  # ============================================================
  # Test 2 -- the GAMMA-decay grid collapses to exactly ONE candidate at
  # max_n_genes_all_studies=1000 ((n_points, n_neighbors)=(300, 1)); with both admissibility
  # gates relaxed, the sole candidate is returned regardless of plateau, and the bootstrap
  # actually runs (test_param_search_single_candidate_bypasses_plateau).
  # ============================================================
  max_n_genes_2 <- 1000L
  gene_means_2 <- matrix(0.0, nrow = max_n_genes_2, ncol = n_studies)
  residuals_2 <- array(0.0, dim = c(max_n_reps_all_studies, max_n_genes_2, n_studies))
  for (i_study in 1:n_studies) {
    gene_means_2[, i_study] <- seq_len(max_n_genes_2)
    for (i_gene in 1:max_n_genes_2) {
      residuals_2[, i_gene, i_study] <- c(-0.5, 0.0, 0.5)
    }
  }

  result2 <- run_js_comp_test_parameter_search(gene_means_2, residuals_2, shared_residual_range = 1.0,
                                                n_bootstraps = 5L, join_method = "join_min",
                                                min_count_per_mean_bin = 0L, min_neighbor_overlap = 0.0,
                                                random_seed = 1L)

  assert_equal_int(as.integer(result2$n_points), 300L, "expected n_points=300")
  assert_equal_int(as.integer(result2$n_neighbors), 1L, "expected n_neighbors=1")
  ci <- result2$best_candidate_pair_confidence_interval
  assert_true(ci[1, 1] >= 0.0 && ci[2, 1] <= 1.0,
              "relaxed gates should let bootstrap actually run, giving a real (not -1.0) CI")

  # ============================================================
  # Test 3 -- callability/shape coverage only (numerical correctness of the effect-size criterion
  # itself is Fortran-only, per this project's testing convention) for the new Issue #178
  # diagnostics: plateau_mode as a mode string, the new threshold optionals, and the trace_*
  # outputs' presence and shape. Reuses Test 2's single-candidate setup.
  # ============================================================
  result3 <- run_js_comp_test_parameter_search(gene_means_2, residuals_2, shared_residual_range = 1.0,
                                                n_bootstraps = 5L, join_method = "join_min",
                                                min_count_per_mean_bin = 0L, min_neighbor_overlap = 0.0,
                                                plateau_mode = "plateau_effect_size", delta_median_threshold = 0.05,
                                                delta_max_threshold = 0.10, delta_epsilon = 1e-10,
                                                delta_min_consecutive_transitions = 2L, random_seed = 1L)

  # n_admissible_evaluated is DM_RESULT_SIZE_IS's own count argument, dropped from the R return
  # since every trace_* output already comes back trimmed to exactly that length.
  for (key in c("trace_n_points", "trace_n_neighbors", "trace_global_js_divergence", "trace_ci_lower",
                "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_delta",
                "trace_delta_median", "trace_delta_max")) {
    assert_true(!is.null(result3[[key]]), paste0("missing expected output '", key, "'"))
  }
  n_admissible <- length(result3$trace_n_points)
  assert_true(n_admissible >= 1L, "expected at least one admissible candidate")
  assert_true(ncol(result3$trace_ci_width) == n_admissible, "trace_ci_width has n_admissible columns")
  assert_true(all(dim(result3$trace_ci_width_relative) == dim(result3$trace_ci_width)),
              "trace_ci_width_relative shape matches trace_ci_width")
  # The first (and here, only) admissible candidate has no predecessor to diff against.
  assert_equal_numeric(result3$trace_delta[, 1], c(-1.0, -1.0), TOL, "first candidate's delta is the -1.0 sentinel")
}

run_all_tests()
