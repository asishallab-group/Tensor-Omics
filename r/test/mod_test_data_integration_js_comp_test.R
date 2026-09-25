# Comprehensive R test suite for the gJCT/fJCT JSD-Comp-Test parameter-search building blocks
# in tensoromics (Issue #126). Fixtures below are translated from the already
# hand-computed/verified numbers in the Fortran test suite
# mod_test_data_integration_js_comp_test.F90, rather than re-derived.
source("r/load_tensor_omics.R")
source("r/test_helpers.R")

TOL <- 1e-12

test_estimate_bin_count <- function() {
  # NOTE: estimate_bin_count/estimate_bin_count_expert now return a named list with n_bins,
  # sturges_bins and fd_bins (Issue #187 checkpoint 1). Per this project's testing philosophy
  # (Fortran_Coding_Guides.pdf Sec 17.1), the R suite only checks call-ability, return
  # shape/type, and the n_bins == max(sturges_bins, fd_bins) structural invariant -- individual
  # numerical correctness of sturges_bins/fd_bins is the Fortran suite's job
  # (mod_test_data_integration_js_comp_test.F90).

  # ============================================================
  # Test 1 -- basic hand-computed (test_estimate_bin_count_basic_hand_computed):
  # Sturges gives 1 + nint(log(1)/LOG_2) = 1. 25th/75th percentiles (rank 3.25 and 7.75)
  # interpolate to -1.75 and 2.75, so the Freedman-Diaconis half-width is
  # (2.75 - (-1.75)) / 1^(1/3) = 4.5; with shared_residual_range=9.0, nint(9.0/4.5)=2.
  # ============================================================
  residuals <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5)

  out <- estimate_bin_count(residuals, max_n_reps_all_studies = 1, n_neighbors = 1,
                             shared_residual_range = 9.0)
  assert_equal_int(out$n_bins, 2L, "Test 1 failed: expected n_bins=2")
  assert_true(is.numeric(out$sturges_bins) && is.numeric(out$fd_bins),
              "Test 1 failed: sturges_bins/fd_bins must be numeric")
  assert_equal_int(out$n_bins, max(as.integer(out$sturges_bins), as.integer(out$fd_bins)),
                    "Test 1 failed: n_bins must equal max(sturges_bins, fd_bins)")

  residuals_perm <- order(residuals)
  out_expert <- estimate_bin_count_expert(residuals, residuals_perm, max_n_reps_all_studies = 1,
                                           n_neighbors = 1, shared_residual_range = 9.0)
  assert_equal_int(out_expert$n_bins, 2L, "Test 1 (expert) failed: expected n_bins=2")
  assert_equal_int(out_expert$n_bins, max(as.integer(out_expert$sturges_bins), as.integer(out_expert$fd_bins)),
                    "Test 1 (expert) failed: n_bins must equal max(sturges_bins, fd_bins)")

  # ============================================================
  # Test 2 -- all-NaN pool falls back to a single bin (test_estimate_bin_count_all_nan_gives_one_bin)
  # ============================================================
  residuals_nan <- rep(NaN, 4)
  out_nan <- estimate_bin_count(residuals_nan, max_n_reps_all_studies = 1, n_neighbors = 1,
                                 shared_residual_range = 9.0)
  assert_equal_int(out_nan$n_bins, 1L, "Test 2 failed: expected n_bins=1")
  assert_equal_int(as.integer(out_nan$sturges_bins), 1L, "Test 2 failed: expected sturges_bins=1")
  assert_equal_int(as.integer(out_nan$fd_bins), 1L, "Test 2 failed: expected fd_bins=1")

  # ============================================================
  # Test 3 -- clamped to MAX_N_BINS=256 (test_estimate_bin_count_clamped_to_max_n_bins)
  # ============================================================
  out_clamped <- estimate_bin_count(residuals, max_n_reps_all_studies = 1, n_neighbors = 1,
                                     shared_residual_range = 5000.0)
  assert_equal_int(out_clamped$n_bins, 256L, "Test 3 failed: expected n_bins=256")
  assert_equal_int(out_clamped$n_bins, max(as.integer(out_clamped$sturges_bins), as.integer(out_clamped$fd_bins)),
                    "Test 3 failed: n_bins must equal max(sturges_bins, fd_bins)")
}

test_determine_bin_count_occupancy <- function() {
  # Call-ability and return type/shape only, per this project's testing philosophy
  # (Fortran_Coding_Guides.pdf Sec 17.1) -- numerical branch coverage of Issue #187's
  # geometric-search-then-refinement algorithm is the Fortran suite's job
  # (mod_test_data_integration_js_comp_test.F90's own 9 enumerated
  # test_determine_bin_count_occupancy_* tests).
  residuals <- seq(-60, 59, by = 1) # 120 values, matches a Fortran fixture

  out <- determine_bin_count_occupancy(residuals, max_n_reps_all_studies = 1, n_neighbors = 1)
  for (key in c("selected_n_bins", "occupancy_failed", "shared_residual_range_low", "shared_residual_range_high",
                "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy",
                "sturges_bins", "fd_bins")) {
    assert_true(key %in% names(out), paste0("missing expected output key '", key, "'"))
  }
  assert_true(is.logical(out$occupancy_failed), "occupancy_failed should be logical")
  assert_equal_int(as.integer(out$n_pooled_residuals), 120L, "expected n_pooled_residuals=120")
  # Call-ability/type/shape only, per this project's testing philosophy (Fortran_Coding_Guides.pdf
  # Sec 17.1) -- shared_residual_range_low/high are now this neighborhood's own 5th/95th
  # percentiles (Step 3), and selected_n_bins is Issue #187's own occupancy-search result; the
  # exact numerical values for this fixture are hand-derived and asserted in the Fortran suite's
  # own test_occupancy_finds_valid_below_m_max, not here.
  assert_true(is.numeric(out$shared_residual_range_low), "shared_residual_range_low should be numeric")
  assert_true(is.numeric(out$shared_residual_range_high), "shared_residual_range_high should be numeric")
  assert_true(is.finite(out$shared_residual_range_low), "shared_residual_range_low must be finite")
  assert_true(is.finite(out$shared_residual_range_high), "shared_residual_range_high must be finite")
  assert_true(out$shared_residual_range_low < out$shared_residual_range_high,
              "shared_residual_range_low must be strictly below shared_residual_range_high on this non-degenerate fixture")
  assert_true(is.numeric(out$selected_n_bins), "selected_n_bins should be numeric")
  assert_true(as.integer(out$selected_n_bins) > 0L, "selected_n_bins must be positive")
  assert_true(!out$occupancy_failed, "occupancy should not fail on this dense fixture")

  # The expert entry point, given the same sorting permutation, must agree.
  residuals_perm <- order(residuals)
  out_expert <- determine_bin_count_occupancy_expert(residuals, residuals_perm, max_n_reps_all_studies = 1,
                                                       n_neighbors = 1)
  assert_equal_int(as.integer(out_expert$selected_n_bins), as.integer(out$selected_n_bins),
                    "expert entry point should agree with the plain one given the same sorted permutation")

  # All-NaN pool: occupancy_failed must be TRUE, and every occupancy diagnostic must be 0.
  residuals_nan <- rep(NaN, 4)
  out_nan <- determine_bin_count_occupancy(residuals_nan, max_n_reps_all_studies = 1, n_neighbors = 1)
  assert_true(out_nan$occupancy_failed, "all-NaN pool must FAIL")
  assert_equal_int(as.integer(out_nan$n_pooled_residuals), 0L, "expected n_pooled_residuals=0")
  assert_equal_int(as.integer(out_nan$min_bin_occupancy), 0L, "expected min_bin_occupancy=0")
  assert_equal_int(as.integer(out_nan$max_bin_occupancy), 0L, "expected max_bin_occupancy=0")
  assert_equal_numeric(out_nan$shared_residual_range_low, 0.0, msg = "expected shared_residual_range_low=0")
  assert_equal_numeric(out_nan$shared_residual_range_high, 0.0, msg = "expected shared_residual_range_high=0")
}

test_determine_bin_count_occupancy_exhaustive <- function() {
  # Call-ability and return type/shape only, per this project's testing philosophy
  # (Fortran_Coding_Guides.pdf Sec 17.1) -- numerical correctness (including the adversarial
  # fixture proving this routine's own reason for existing) is the Fortran suite's job
  # (mod_test_data_integration_js_comp_test.F90's own test_occupancy_exhaustive_* tests).
  # Unlike determine_bin_count_occupancy, this routine takes shared_residual_range_low/high and
  # n_pooled_residuals as direct inputs rather than deriving them, so literal values are enough.
  residuals <- seq(-60, 59, by = 1) # 120 values, matches a Fortran fixture

  out <- determine_bin_count_occupancy_exhaustive(residuals, n_pooled_residuals = 120L,
                                                    shared_residual_range_low = -60.0,
                                                    shared_residual_range_high = 60.0)
  for (key in c("selected_n_bins", "occupancy_failed", "min_bin_occupancy", "mean_bin_occupancy",
                "max_bin_occupancy")) {
    assert_true(key %in% names(out), paste0("missing expected output key '", key, "'"))
  }
  assert_true(is.logical(out$occupancy_failed), "occupancy_failed should be logical")
  assert_true(is.numeric(out$selected_n_bins), "selected_n_bins should be numeric")
  assert_true(as.integer(out$selected_n_bins) > 0L, "selected_n_bins must be positive")
  assert_true(!out$occupancy_failed, "occupancy should not fail on this dense fixture")

  # The expert entry point, given the same sorting permutation, must agree.
  residuals_perm <- order(residuals)
  out_expert <- determine_bin_count_occupancy_exhaustive_expert(residuals, residuals_perm,
                                                                  n_pooled_residuals = 120L,
                                                                  shared_residual_range_low = -60.0,
                                                                  shared_residual_range_high = 60.0)
  assert_equal_int(as.integer(out_expert$selected_n_bins), as.integer(out$selected_n_bins),
                    "expert entry point should agree with the plain one given the same sorted permutation")

  # All-NaN pool: occupancy_failed must be TRUE (n_pooled_residuals=0, default min_residuals_per_bin=10).
  residuals_nan <- rep(NaN, 4)
  out_nan <- determine_bin_count_occupancy_exhaustive(residuals_nan, n_pooled_residuals = 0L,
                                                        shared_residual_range_low = 0.0,
                                                        shared_residual_range_high = 0.0)
  assert_true(out_nan$occupancy_failed, "all-NaN pool must FAIL")
  assert_equal_int(as.integer(out_nan$min_bin_occupancy), 0L, "expected min_bin_occupancy=0")
  assert_equal_int(as.integer(out_nan$max_bin_occupancy), 0L, "expected max_bin_occupancy=0")
}

test_generate_js_comp_test_candidates <- function() {
  # ============================================================
  # Test 1 -- small-N candidate-grid collapse at its exact threshold
  # (test_generate_js_comp_test_candidates_collapses_at_8742): every candidate shares n_points=374.
  # ============================================================
  candidates <- generate_js_comp_test_candidates(8742L)
  # Issue #187, Step 2.7: `candidates_n_points_n_neighbors` is now the routine's only
  # DM_RESULT_SIZE_IS output, so the generated R wrapper's own bracket slice
  # (`[, seq_len(n_candidates), drop = FALSE]`) already returns it correctly trimmed by column --
  # no manual re-slicing needed (contrast the old n_bins_candidates-vs-head() workaround this test
  # used to need, now gone along with n_bins_candidates itself).
  n_candidates <- ncol(candidates)
  assert_true(n_candidates >= 1, "Test 1 failed: expected at least one candidate")
  assert_true(all(candidates[1, ] == candidates[1, 1]),
              "Test 1 failed: every candidate must share the same n_points")
  assert_equal_int(as.integer(candidates[1, 1]), 374L, "Test 1 failed: n_points should be 374")

  # ============================================================
  # Test 2 -- the other side of the bracket (test_generate_js_comp_test_candidates_has_two_distinct_at_8743):
  # at least two distinct n_points values, first=375, last=300.
  # ============================================================
  candidates2 <- generate_js_comp_test_candidates(8743L)
  n_candidates2 <- ncol(candidates2)
  n_distinct <- length(unique(candidates2[1, ]))
  assert_true(n_distinct >= 2, "Test 2 failed: expected at least two distinct n_points values")
  assert_equal_int(as.integer(candidates2[1, 1]), 375L, "Test 2 failed: first n_points should be 375")
  assert_equal_int(as.integer(candidates2[1, n_candidates2]), 300L,
                    "Test 2 failed: last n_points should collapse to 300")

  # ============================================================
  # Test 3 -- validation (test_generate_js_comp_test_candidates_validation). Issue #187, Step 2.7:
  # this routine's `_expert` tier no longer exists (no tmp_/work-array/permutation left in its
  # signature once the bin-estimate side effect was removed), so this test no longer exercises it.
  # ============================================================
  assert_error(
    generate_js_comp_test_candidates(0L),
    "Test 3 failed: expected ERR_INVALID_INPUT for max_n_genes_all_studies=0", ERR_INVALID_INPUT
  )
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
  # n_bins_per_point uniformly equals n_bins (call-ability/shape only, per this file's convention)
  mean_pmf_counts <- matrix(c(10, 10, 10, 10), nrow = 2)
  assert_true(check_mean_pmf_min_counts(mean_pmf_counts, n_bins_per_point = c(2L, 2L), min_count = 5L),
              "Test 1 failed: all bins above minimum should pass")

  # Test 2 -- exactly at threshold passes (test_mean_pmf_min_counts_exactly_at_threshold_passes)
  mean_pmf_counts_2 <- matrix(c(3, 3, 3, 3), nrow = 2)
  assert_true(check_mean_pmf_min_counts(mean_pmf_counts_2, n_bins_per_point = c(2L, 2L), min_count = 3L),
              "Test 2 failed: count exactly at minimum must pass")

  # Test 3 -- one bin below threshold fails (test_mean_pmf_min_counts_below_threshold_fails)
  mean_pmf_counts_3 <- matrix(c(3, 3, 3, 2), nrow = 2)
  assert_false(check_mean_pmf_min_counts(mean_pmf_counts_3, n_bins_per_point = c(2L, 2L), min_count = 3L),
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

  # Issue #187, Step 2.7: generate_js_comp_test_candidates now returns the already-correctly-
  # trimmed matrix directly (see the NOTE in test_generate_js_comp_test_candidates).
  candidates <- generate_js_comp_test_candidates(8742L)
  assert_true(bounds$max_n_neighbors_candidate >= max(candidates[2, ]),
              "max_n_neighbors_candidate must be a safe upper bound on the grid's n_neighbors")
}

test_gather_pooled_neighborhood_residuals <- function() {
  # Call-ability and return type/shape only, per this project's testing philosophy
  # (Fortran_Coding_Guides.pdf Sec 17.1) -- numerical correctness is the Fortran suite's job
  # (mod_test_data_integration_js_comp_test.F90's own test_gather_pooled_residuals_* tests).
  # max_n_reps_all_studies/max_n_genes_all_studies/n_neighbors/n_studies must all auto-derive
  # from residuals'/neighborhood_indices_point's own shapes -- not asked of the caller.
  max_n_reps_all_studies <- 2L
  max_n_genes_all_studies <- 3L
  n_neighbors <- 2L
  n_studies <- 2L
  residuals <- array(as.double(seq_len(max_n_reps_all_studies*max_n_genes_all_studies*n_studies)),
                      dim = c(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
  neighborhood_indices_point <- matrix(c(1L, 2L, 3L, 1L), nrow = n_neighbors, ncol = n_studies)

  pooled_residuals <- gather_pooled_neighborhood_residuals(residuals, neighborhood_indices_point)

  assert_true(is.numeric(pooled_residuals), "expected a numeric vector")
  assert_equal_int(length(pooled_residuals), max_n_reps_all_studies*n_neighbors*n_studies,
                    "expected pooled_residuals of length max_n_reps_all_studies*n_neighbors*n_studies")

  # An out-of-range gene index must raise, confirming the new bounds check is wired through.
  bad_indices <- matrix(c(1L, 2L, max_n_genes_all_studies + 1L, 1L), nrow = n_neighbors, ncol = n_studies)
  assert_error(gather_pooled_neighborhood_residuals(residuals, bad_indices),
               "expected ERR_INVALID_INPUT for an out-of-range gene index", ERR_INVALID_INPUT)
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
  # Call-ability/shape check only (Fortran_Coding_Guides.pdf Sec 17.1) -- the same fixture as the
  # Fortran suite's own test_run_js_comp_test_two_studies_hand_traceable, which is where the exact
  # hand-derived closed-form JSD assertions live. Issue #187 removed the mandatory scalar `n_bins`
  # input (the bin count is now decided per reference point by an internal occupancy search) and
  # added `n_bins_per_point`/`max_n_bins_per_point` plus several per-point diagnostics;
  # `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts` are now always sized to the fixed 256-bin ceiling
  # (MAX_N_BINS), of which only the leading `max_n_bins_per_point` rows are meaningful -- a caller
  # must slice `[1:max_n_bins_per_point, ...]` themselves.
  gene_means <- matrix(c(1.0, 5.0, 1.0, 5.0), nrow = 2)
  gene_means_perms <- matrix(c(1L, 2L, 1L, 2L), nrow = 2)

  residuals <- array(NaN, dim = c(2, 2, 2))
  residuals[, 1, 1] <- c(-1.0, 1.0)
  residuals[, 1, 2] <- c(-3.0, 3.0)

  x_star <- c(1.0)

  # min_residuals_per_bin=1 (its default is 10, which this tiny 4-residual fixture could never
  # satisfy) so the occupancy search actually succeeds, matching the Fortran fixture's own
  # override -- exercised here for call-ability, not to reproduce that test's exact numbers.
  result <- run_js_comp_test(n_neighbors = 1L, gene_means = gene_means,
                              gene_means_perms = gene_means_perms, residuals = residuals, x_star = x_star,
                              n_permutations = 0L, random_seed = 1L, min_residuals_per_bin = 1L)

  n_points <- 1L
  n_studies <- 2L
  assert_equal_int(dim(result$neighborhood_indices), c(1L, n_points, n_studies), "neighborhood_indices shape")
  assert_equal_int(dim(result$neighborhood_range), c(2L, n_points, n_studies), "neighborhood_range shape")
  assert_equal_int(length(result$n_bins_per_point), n_points, "n_bins_per_point shape")
  assert_equal_int(length(result$shared_residual_range_low), n_points, "shared_residual_range_low shape")
  assert_equal_int(length(result$shared_residual_range_high), n_points, "shared_residual_range_high shape")
  assert_true(is.numeric(result$max_n_bins_per_point), "max_n_bins_per_point must be numeric")
  assert_equal_int(length(result$occupancy_failed), n_points, "occupancy_failed shape")
  assert_equal_int(length(result$n_pooled_residuals), n_points, "n_pooled_residuals shape")
  assert_equal_int(length(result$min_bin_occupancy), n_points, "min_bin_occupancy shape")
  assert_equal_int(length(result$mean_bin_occupancy), n_points, "mean_bin_occupancy shape")
  assert_equal_int(length(result$max_bin_occupancy), n_points, "max_bin_occupancy shape")
  assert_equal_int(length(result$sturges_bins), n_points, "sturges_bins shape")
  assert_equal_int(length(result$fd_bins), n_points, "fd_bins shape")
  assert_equal_int(dim(result$pmfs), c(256L, n_points, n_studies), "pmfs shape")
  assert_equal_int(dim(result$counts), c(256L, n_points, n_studies), "counts shape")
  assert_equal_int(dim(result$included_n_reps), c(n_points, n_studies), "included_n_reps shape")
  assert_equal_int(dim(result$mean_pmf), c(256L, n_points), "mean_pmf shape")
  assert_equal_int(dim(result$mean_pmf_counts), c(256L, n_points), "mean_pmf_counts shape")
  assert_equal_int(length(result$mean_pmf_included_n_reps), n_points, "mean_pmf_included_n_reps shape")
  assert_equal_int(dim(result$js_divergences), c(n_points, n_studies), "js_divergences shape")
  assert_equal_int(dim(result$weights), c(n_points, n_studies), "weights shape")
  assert_equal_int(length(result$global_js_divergence), n_studies, "global_js_divergence shape")
  assert_equal_int(length(result$p_values), n_studies, "p_values shape")

  # A caller must be able to slice down to the meaningful leading bins.
  max_n_bins <- as.integer(result$max_n_bins_per_point)
  assert_equal_int(dim(result$pmfs[1:max_n_bins, , , drop = FALSE]), c(max_n_bins, n_points, n_studies),
                    "pmfs slices down to max_n_bins_per_point")

  assert_equal_int(as.integer(result$neighborhood_indices[, 1, 1]), 1L, "study 1 neighbor should be gene 1")
  assert_equal_int(as.integer(result$neighborhood_indices[, 1, 2]), 1L, "study 2 neighbor should be gene 1")

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
  # Test 1 -- no candidate ever plateaus (min_residuals_per_bin impossibly high): falls back to
  # the FIRST (finest-resolution) candidate, (n_points, n_neighbors)=(300, 26), and resets
  # best_candidate_pair_confidence_interval to -1.0 throughout
  # (test_param_search_no_plateau_falls_back_to_finest).
  # ============================================================
  bounds <- calc_js_comp_test_candidate_bounds(max_n_genes_all_studies)
  result <- run_js_comp_test_parameter_search(gene_means, residuals, n_bootstraps = 10L,
                                               join_method = "join_min",
                                               max_n_points_candidate = bounds$max_n_points_candidate,
                                               max_n_neighbors_candidate = bounds$max_n_neighbors_candidate,
                                               min_residuals_per_bin = 1000000L,
                                               random_seed = 1L)

  assert_equal_int(as.integer(result$n_points), 300L, "expected n_points=300")
  assert_equal_int(as.integer(result$n_neighbors), 26L, "expected n_neighbors=26")
  assert_false(result$plateau_established, "no candidate ever plateaus here")
  assert_equal_numeric(result$best_candidate_pair_confidence_interval[, 1], c(-1.0, -1.0), TOL, "study 1 CI reset")
  assert_equal_numeric(result$best_candidate_pair_confidence_interval[, 2], c(-1.0, -1.0), TOL, "study 2 CI reset")

  # ============================================================
  # Test 2 -- the GAMMA-decay grid collapses to exactly ONE candidate at
  # max_n_genes_all_studies=100 ((n_points, n_neighbors)=(300, 1)); with both admissibility
  # gates relaxed, the sole candidate is returned regardless of plateau, and the bootstrap
  # actually runs (test_param_search_single_candidate_bypasses_plateau).
  # ============================================================
  max_n_genes_2 <- 100L
  gene_means_2 <- matrix(0.0, nrow = max_n_genes_2, ncol = n_studies)
  residuals_2 <- array(0.0, dim = c(max_n_reps_all_studies, max_n_genes_2, n_studies))
  for (i_study in 1:n_studies) {
    gene_means_2[, i_study] <- seq_len(max_n_genes_2)
    for (i_gene in 1:max_n_genes_2) {
      residuals_2[, i_gene, i_study] <- c(-0.5, 0.0, 0.5)
    }
  }

  bounds_2 <- calc_js_comp_test_candidate_bounds(max_n_genes_2)
  result2 <- run_js_comp_test_parameter_search(gene_means_2, residuals_2,
                                                n_bootstraps = 5L, join_method = "join_min",
                                                max_n_points_candidate = bounds_2$max_n_points_candidate,
                                                max_n_neighbors_candidate = bounds_2$max_n_neighbors_candidate,
                                                min_residuals_per_bin = 0L, min_neighbor_overlap = 0.0,
                                                random_seed = 1L)

  assert_equal_int(as.integer(result2$n_points), 300L, "expected n_points=300")
  assert_equal_int(as.integer(result2$n_neighbors), 1L, "expected n_neighbors=1")
  assert_true(result2$plateau_established, "the single-candidate bypass counts as established")
  ci <- result2$best_candidate_pair_confidence_interval
  assert_true(ci[1, 1] >= 0.0 && ci[2, 1] <= 1.0,
              "relaxed gates should let bootstrap actually run, giving a real (not -1.0) CI")

  # ============================================================
  # Test 3 -- callability/shape coverage only (numerical correctness of the effect-size criterion
  # itself is Fortran-only, per this project's testing convention) for the new Issue #178
  # diagnostics: plateau_mode as a mode string, the new threshold optionals, and the trace_*
  # outputs' presence and shape. Reuses Test 2's single-candidate setup.
  # ============================================================
  result3 <- run_js_comp_test_parameter_search(gene_means_2, residuals_2,
                                                n_bootstraps = 5L, join_method = "join_min",
                                                max_n_points_candidate = bounds_2$max_n_points_candidate,
                                                max_n_neighbors_candidate = bounds_2$max_n_neighbors_candidate,
                                                min_residuals_per_bin = 0L, min_neighbor_overlap = 0.0,
                                                plateau_mode = "plateau_effect_size", delta_median_threshold = 0.05,
                                                delta_max_threshold = 0.10, delta_epsilon = 1e-10,
                                                delta_min_consecutive_transitions = 2L, random_seed = 1L)

  # n_admissible_evaluated is DM_RESULT_SIZE_IS's own count argument, dropped from the R return
  # since every trace_* output already comes back trimmed to exactly that length.
  for (key in c("trace_n_points", "trace_n_neighbors", "trace_global_js_divergence", "trace_ci_lower",
                "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_delta",
                "trace_delta_median", "trace_delta_max", "plateau_established")) {
    assert_true(!is.null(result3[[key]]), paste0("missing expected output '", key, "'"))
  }
  assert_true(is.logical(result3$plateau_established), "plateau_established should be logical")
  n_admissible <- length(result3$trace_n_points)
  assert_true(n_admissible >= 1L, "expected at least one admissible candidate")
  assert_true(ncol(result3$trace_ci_width) == n_admissible, "trace_ci_width has n_admissible columns")
  assert_true(all(dim(result3$trace_ci_width_relative) == dim(result3$trace_ci_width)),
              "trace_ci_width_relative shape matches trace_ci_width")
  # The first (and here, only) admissible candidate has no predecessor to diff against.
  assert_equal_numeric(result3$trace_delta[, 1], c(-1.0, -1.0), TOL, "first candidate's delta is the -1.0 sentinel")

  # ============================================================
  # Test 4 -- callability/shape coverage for the Issue #187 per-point outputs added by this
  # signature change: n_bins_per_point (replacing the old scalar n_bins) and the 8 new trace_*
  # diagnostics from determine_bin_count_occupancy. Numerical correctness of the occupancy search
  # itself is Fortran-only (mod_test_data_integration_js_comp_test.F90's
  # test_param_search_occupancy_failure_rejects_candidate/
  # test_param_search_different_neighborhoods_different_m_j/
  # test_param_search_final_n_bins_matches_selected_trace_column).
  # ============================================================
  for (key in c("n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high",
                "trace_selected_n_bins", "trace_occupancy_failed",
                "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy",
                "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins",
                "trace_shared_residual_range_low", "trace_shared_residual_range_high")) {
    assert_true(!is.null(result2[[key]]), paste0("missing expected output '", key, "'"))
  }
  assert_true(is.integer(result2$n_bins_per_point), "n_bins_per_point should be an integer vector")
  assert_true(length(result2$n_bins_per_point) == bounds_2$max_n_points_candidate,
              "n_bins_per_point should be sized to max_n_points_candidate")
  assert_true(is.numeric(result2$shared_residual_range_low), "shared_residual_range_low should be numeric")
  assert_true(length(result2$shared_residual_range_low) == bounds_2$max_n_points_candidate,
              "shared_residual_range_low should be sized to max_n_points_candidate")
  assert_true(length(result2$shared_residual_range_high) == length(result2$shared_residual_range_low),
              "shared_residual_range_high shape matches shared_residual_range_low")
  assert_true(is.integer(result2$trace_selected_n_bins), "trace_selected_n_bins should be an integer matrix")
  assert_true(nrow(result2$trace_selected_n_bins) == bounds_2$max_n_points_candidate,
              "trace_selected_n_bins should have max_n_points_candidate rows")
  assert_true(is.logical(result2$trace_occupancy_failed), "trace_occupancy_failed should be logical")
  assert_true(all(dim(result2$trace_occupancy_failed) == dim(result2$trace_selected_n_bins)),
              "trace_occupancy_failed shape matches trace_selected_n_bins")
  for (key in c("trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_max_bin_occupancy",
                "trace_sturges_bins", "trace_fd_bins",
                "trace_shared_residual_range_low", "trace_shared_residual_range_high")) {
    assert_true(all(dim(result2[[key]]) == dim(result2$trace_selected_n_bins)),
                paste0(key, " shape mismatch"))
  }
  # n_points=300 leading entries of n_bins_per_point must match the sole admissible candidate's
  # own trace_selected_n_bins column (test_param_search_final_n_bins_matches_selected_trace_column
  # is the Fortran suite's rigorous version of this same wiring check).
  n_points_2 <- result2$n_points
  assert_true(all(result2$n_bins_per_point[1:n_points_2] == result2$trace_selected_n_bins[1:n_points_2, 1]),
              "n_bins_per_point matches the sole admissible candidate's trace column")
}

run_all_tests()
