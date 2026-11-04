# Comprehensive R test suite for the gJCT/fJCT JSD-Comp-Test parameter-search building blocks
# in tensoromics (Issue #126). Fixtures below are translated from the already
# hand-computed/verified numbers in the Fortran test suite
# mod_test_data_integration_js_comp_test.F90, rather than re-derived.
source("r/load_tensor_omics.R")
source("r/test_helpers.R")

TOL <- 1e-12

# Every `bin_count_mode` string of run_js_comp_test, run_js_comp_test_adaptive and both
# parameter searches.
BIN_COUNT_MODES <- c("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted",
                     "bin_count_stone_pooled")

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

test_determine_bin_count_stone <- function() {
  # Call-ability and return type/shape only, per this project's testing philosophy
  # (Fortran_Coding_Guides.pdf Sec 17.1) -- the numerics of Stone's rule are the Fortran suite's
  # job (mod_test_data_integration_js_comp_test.F90's own test_stone_* tests).
  set.seed(7)
  residuals <- rnorm(300)
  residual_to_study <- sample(rep(1:3, times = c(120, 60, 120)))

  for (mode in c("bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled")) {
    out <- determine_bin_count_stone(residuals, residual_to_study, n_studies = 3L, max_n_reps_all_studies = 1L,
                                     n_neighbors = 1L, bin_count_mode = mode)
    for (key in c("selected_n_bins", "stone_failed", "shared_residual_range_low", "shared_residual_range_high",
                  "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy",
                  "sturges_bins", "fd_bins")) {
      assert_true(key %in% names(out), paste0("missing expected output key '", key, "' (", mode, ")"))
    }
    assert_true(is.logical(out$stone_failed), paste("stone_failed should be logical", mode))
    assert_true(is.numeric(out$selected_n_bins), paste("selected_n_bins should be numeric", mode))
    assert_true(is.numeric(out$mean_bin_occupancy), paste("mean_bin_occupancy should be numeric", mode))
    assert_true(!out$stone_failed, paste("Stone should not fail on this dense fixture", mode))
    assert_true(out$selected_n_bins >= 3 && out$selected_n_bins <= 120,
                paste("selected_n_bins within [m_min, m_max]", mode))
    assert_true(out$shared_residual_range_low < out$shared_residual_range_high,
                paste("non-degenerate range expected", mode))
  }

  # The default mode is accepted without being named.
  out <- determine_bin_count_stone(residuals, residual_to_study, n_studies = 3L, max_n_reps_all_studies = 1L,
                                   n_neighbors = 1L)
  assert_true(!out$stone_failed, "default mode should not fail")

  # The expert entry point, given the same sorting permutation, must agree.
  out_expert <- determine_bin_count_stone_expert(residuals, order(residuals), residual_to_study, n_studies = 3L,
                                                 max_n_reps_all_studies = 1L, n_neighbors = 1L)
  assert_equal_int(as.integer(out_expert$selected_n_bins), as.integer(out$selected_n_bins),
                   "expert entry point should agree with the plain one given the same sorted permutation")

  # An unknown mode string, and the occupancy mode (not a Stone mode), are rejected before Fortran.
  assert_error(determine_bin_count_stone(residuals, residual_to_study, n_studies = 3L, max_n_reps_all_studies = 1L,
                                         n_neighbors = 1L, bin_count_mode = "bin_count_occupancy"),
               "expected the occupancy mode to be rejected")
  assert_error(determine_bin_count_stone(residuals, residual_to_study, n_studies = 3L, max_n_reps_all_studies = 1L,
                                         n_neighbors = 1L, bin_count_mode = "nonsense"),
               "expected an unknown mode string to be rejected")
  # A study label outside [1, n_studies] is invalid input.
  assert_error(determine_bin_count_stone(residuals, residual_to_study, n_studies = 2L, max_n_reps_all_studies = 1L,
                                         n_neighbors = 1L),
               "expected ERR_INVALID_INPUT for a label above n_studies", ERR_INVALID_INPUT)

  # All-NaN pool: the search fails.
  out_nan <- determine_bin_count_stone(rep(NaN, 4), c(1L, 2L, 1L, 2L), n_studies = 2L,
                                       max_n_reps_all_studies = 1L, n_neighbors = 1L)
  assert_true(out_nan$stone_failed, "all-NaN pool must fail")
}

test_generate_js_comp_test_candidates <- function() {
  # ============================================================
  # Test 1 -- the default grid by hand (grid_points_factor 4, grid_min_points 200), as in the
  # Python suite: 8742 genes give n_points 374, 299, 239 and the floor 200 as the last candidate,
  # with n_neighbors floor(8742/(0.25*n_points)) = 93, 116, 146, 174. The routine's only
  # DM_RESULT_SIZE_IS output, so the wrapper returns the trimmed matrix.
  # ============================================================
  candidates <- generate_js_comp_test_candidates(8742L)
  assert_true(is.matrix(candidates) && is.integer(candidates), "Test 1 failed: an integer matrix")
  assert_equal_int(dim(candidates), c(2L, 4L), "Test 1 failed: a 2 x 4 grid")
  assert_equal_int(as.integer(candidates[1, ]), c(374L, 299L, 239L, 200L), "Test 1 failed: n_points")
  assert_equal_int(as.integer(candidates[2, ]), c(93L, 116L, 146L, 174L), "Test 1 failed: n_neighbors")

  # ============================================================
  # Test 2 -- the two grid options: grid_min_points = 300 gives (374, 93), (300, 116);
  # grid_points_factor = 2 gives the single candidate (200, 174); bad values are rejected.
  # ============================================================
  floor_300 <- generate_js_comp_test_candidates(8742L, grid_min_points = 300L)
  assert_equal_int(as.integer(floor_300), c(374L, 93L, 300L, 116L), "Test 2 failed: grid_min_points = 300")
  factor_2 <- generate_js_comp_test_candidates(8742L, grid_points_factor = 2.0)
  assert_equal_int(as.integer(factor_2), c(200L, 174L), "Test 2 failed: grid_points_factor = 2")
  assert_error(generate_js_comp_test_candidates(8742L, grid_points_factor = 0.0),
               "Test 2 failed: expected ERR_INVALID_INPUT for grid_points_factor = 0", ERR_INVALID_INPUT)
  assert_error(generate_js_comp_test_candidates(8742L, grid_min_points = 1501L),
               "Test 2 failed: expected ERR_INVALID_INPUT for grid_min_points above 1500", ERR_INVALID_INPUT)

  # ============================================================
  # Test 3 -- validation (test_generate_js_comp_test_candidates_validation).
  # ============================================================
  assert_error(
    generate_js_comp_test_candidates(0L),
    "Test 3 failed: expected ERR_INVALID_INPUT for max_n_genes_all_studies=0", ERR_INVALID_INPUT
  )

  # ============================================================
  # Test 4 -- neighborhood_size_mode as a mode string (call, type, shape): both modes share the
  # n_points, the small one has fewer neighbors; an unknown string (also the removed "both") is
  # rejected before Fortran.
  # ============================================================
  large <- generate_js_comp_test_candidates(8742L, neighborhood_size_mode = "neighborhood_size_large")
  small <- generate_js_comp_test_candidates(8742L, neighborhood_size_mode = "neighborhood_size_small")
  for (grid in list(large, small)) {
    assert_true(is.matrix(grid) && nrow(grid) == 2L && is.integer(grid) && all(grid > 0L), "Test 4 failed: shape")
  }
  assert_equal_int(as.integer(small[1, ]), as.integer(large[1, ]), "Test 4 failed: both modes share the n_points")
  assert_true(all(small[2, ] < large[2, ]), "Test 4 failed: the small neighborhoods have fewer neighbors")
  assert_error(generate_js_comp_test_candidates(2000L, neighborhood_size_mode = "nonsense"),
               "Test 4 failed: expected an unknown neighborhood_size_mode to be rejected before Fortran")
  assert_error(generate_js_comp_test_candidates(2000L, neighborhood_size_mode = "neighborhood_size_both"),
               "Test 4 failed: expected the removed mode 'both' to be rejected before Fortran")
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
  # inclusive position ranges [1,10] and [6,10] share the 5 positions 6..10 of the first range's
  # 10, so overlap = 5/10 = 0.5.
  # ============================================================
  neighborhood_range_2 <- matrix(c(1, 10, 6, 10), nrow = 2)
  assert_true(check_neighborhood_overlaps(neighborhood_range_2, min_neighbor_overlap = 0.5),
              "Test 2 failed: overlap exactly at threshold must pass")

  # ============================================================
  # Test 3 -- same 5/10 overlap, threshold 0.51 > 0.5, must fail
  # (test_neighborhood_overlaps_below_threshold_fails)
  # ============================================================
  assert_false(check_neighborhood_overlaps(neighborhood_range_2, min_neighbor_overlap = 0.51),
               "Test 3 failed: overlap below threshold must fail")

  # ============================================================
  # Test 4 -- the overlap counts shared positions of inclusive ranges
  # (test_neighborhood_overlaps_counts_inclusive_positions): [5,6] and [6,7] share 1 of 2
  # positions (0.5), [1,100] and [91,190] share 10 of 100 (exactly 0.1); both pass at those
  # exact thresholds, and the second still fails just above, at 0.11.
  # ============================================================
  neighborhood_range_3 <- matrix(c(5, 6, 6, 7), nrow = 2)
  assert_true(check_neighborhood_overlaps(neighborhood_range_3, min_neighbor_overlap = 0.5),
              "Test 4 failed: [5,6]/[6,7] share 1 of 2 positions, must pass at 0.5")
  neighborhood_range_4 <- matrix(c(1, 100, 91, 190), nrow = 2)
  assert_true(check_neighborhood_overlaps(neighborhood_range_4, min_neighbor_overlap = 0.1),
              "Test 4 failed: [1,100]/[91,190] share 10 of 100 positions, must pass at 0.1")
  assert_false(check_neighborhood_overlaps(neighborhood_range_4, min_neighbor_overlap = 0.11),
               "Test 4 failed: [1,100]/[91,190] share only 10 of 100 positions, must fail at 0.11")
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
  # Studies 1,2 sit inside the best interval (overlap 1.0), study 3 is far outside (0.0).
  assert_equal_numeric(as.numeric(result$ci_overlap), c(1.0, 1.0, 0.0), TOL, "Test 1 failed: ci_overlap of 2-of-3")

  confidence_interval[, 3] <- c(0.4, 0.6)
  best_ci_2 <- matrix(c(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), nrow = 2)
  result2 <- check_plateau_condition(confidence_interval, best_ci_2, best_candidate_index = 0L,
                                      best_exceeded_ci_overlap_count = 0L, candidate_index = 1L,
                                      join_method = "join_min", succeeding_ci_overlap = 0.9)
  assert_true(result2$plateau_found, "Test 1 failed: 3-of-3 studies must plateau")
  assert_equal_numeric(as.numeric(result2$ci_overlap), c(1.0, 1.0, 1.0), TOL, "Test 1 failed: ci_overlap of 3-of-3")

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
  assert_equal_numeric(as.numeric(result_max_0$ci_overlap), c(0.0, 0.0, 0.0), TOL, "Test 2 failed: ci_overlap of 0-of-3")

  confidence_interval_max[, 1] <- c(0.4, 0.6)
  best_ci_max_2 <- matrix(c(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), nrow = 2)
  result_max_1 <- check_plateau_condition(confidence_interval_max, best_ci_max_2, best_candidate_index = 0L,
                                           best_exceeded_ci_overlap_count = 0L, candidate_index = 1L,
                                           join_method = "join_max", succeeding_ci_overlap = 0.9)
  assert_true(result_max_1$plateau_found, "Test 2 failed: 1-of-3 studies must plateau")
  assert_equal_numeric(as.numeric(result_max_1$ci_overlap), c(1.0, 0.0, 0.0), TOL, "Test 2 failed: ci_overlap of 1-of-3")

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
  assert_equal_numeric(as.numeric(result_med_1$ci_overlap), c(1.0, 0.0, 0.0), TOL, "Test 3 failed: ci_overlap of 1-of-3")

  confidence_interval_med[, 2] <- c(0.4, 0.6)
  best_ci_med_2 <- matrix(c(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), nrow = 2)
  result_med_2 <- check_plateau_condition(confidence_interval_med, best_ci_med_2, best_candidate_index = 0L,
                                           best_exceeded_ci_overlap_count = 0L, candidate_index = 1L,
                                           join_method = "join_median", succeeding_ci_overlap = 0.9)
  assert_true(result_med_2$plateau_found, "Test 3 failed: 2-of-3 must plateau")
  assert_equal_numeric(as.numeric(result_med_2$ci_overlap), c(1.0, 1.0, 0.0), TOL, "Test 3 failed: ci_overlap of 2-of-3")

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
  # The overlaps are still reported on the worse-than-best shortcut.
  assert_equal_numeric(as.numeric(result_worse$ci_overlap), c(0.0, 0.0), TOL, "Test 4 failed: ci_overlap still reported")

  # ============================================================
  # Test 5 -- fractional overlaps (test_check_plateau_condition_ci_overlap_fractional_values):
  # against the best interval [0, 1], [0.4, 0.6] lies inside (1.0), [-0.5, 0.5] shares half of its
  # width (0.5), [-2, -1] shares nothing (0.0), and [0, 2] is twice as wide as the best interval,
  # so only half of the candidate is covered (0.5: the candidate's width is the denominator).
  # ============================================================
  confidence_interval_frac <- matrix(c(0.4, 0.6, -0.5, 0.5, -2.0, -1.0, 0.0, 2.0), nrow = 2)
  best_ci_frac <- matrix(c(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0), nrow = 2)
  result_frac <- check_plateau_condition(confidence_interval_frac, best_ci_frac, best_candidate_index = 0L,
                                          best_exceeded_ci_overlap_count = 0L, candidate_index = 1L,
                                          join_method = "join_min", succeeding_ci_overlap = 0.9)
  assert_equal_numeric(as.numeric(result_frac$ci_overlap), c(1.0, 0.5, 0.0, 0.5), TOL,
                        "Test 5 failed: fractional ci_overlap values")
}

test_check_effect_size_plateau_condition <- function() {
  # Call-ability and return type/shape (Fortran_Coding_Guides.pdf Sec 17.1); the counting rules
  # themselves are the Fortran suite's test_effect_size_plateau_* tests.
  # ============================================================
  # Test 1 -- first admissible candidate (test_effect_size_plateau_first_candidate_no_delta):
  # no transition, so delta/delta_median/delta_max are the documented -1.0 and the count resets.
  # ============================================================
  result <- check_effect_size_plateau_condition(c(0.1, 0.2), c(0.0, 0.0), has_previous = FALSE,
                                                delta_median_threshold = 0.05, delta_max_threshold = 0.10,
                                                delta_epsilon = 1.0e-10, delta_min_consecutive_transitions = 2L,
                                                n_consecutive_ok = 5L)
  assert_true(is.list(result), "expected a named list")
  assert_true(is.numeric(result$delta), "expected delta to be numeric")
  assert_equal_int(length(result$delta), 2L, "expected delta of length n_studies")
  assert_equal_numeric(as.numeric(result$delta), c(-1.0, -1.0), TOL, "Test 1 failed: delta is the -1.0 sentinel")
  assert_equal_numeric(c(result$delta_median, result$delta_max), c(-1.0, -1.0), TOL,
                       "Test 1 failed: delta_median and delta_max are the -1.0 sentinel")
  assert_equal_int(as.integer(result$n_consecutive_ok), 0L, "Test 1 failed: n_consecutive_ok must reset to 0")
  assert_false(result$plateau_found, "Test 1 failed: no plateau on the first admissible candidate")

  # ============================================================
  # Test 2 -- a qualifying transition (|0.101 - 0.100| / 0.100 = 0.01 for both studies): one alone
  # is not enough (test_effect_size_plateau_single_transition_insufficient), the second in a row
  # is (test_effect_size_plateau_two_consecutive_transitions).
  # ============================================================
  for (n_ok_before in 0:1) {
    result <- check_effect_size_plateau_condition(c(0.101, 0.101), c(0.100, 0.100), has_previous = TRUE,
                                                  delta_median_threshold = 0.05, delta_max_threshold = 0.10,
                                                  delta_epsilon = 1.0e-10, delta_min_consecutive_transitions = 2L,
                                                  n_consecutive_ok = n_ok_before)
    assert_equal_int(as.integer(result$n_consecutive_ok), n_ok_before + 1L,
                     paste("Test 2 failed: n_consecutive_ok after", n_ok_before + 1L, "transition(s)"))
    assert_true(is.logical(result$plateau_found) && length(result$plateau_found) == 1,
                "Test 2 failed: plateau_found must be a single logical")
    assert_true(identical(result$plateau_found, n_ok_before == 1L),
                paste("Test 2 failed: plateau_found after", n_ok_before + 1L, "transition(s)"))
    assert_true(is.numeric(result$delta_median) && is.numeric(result$delta_max),
                "Test 2 failed: delta_median and delta_max must be numeric")
  }

  # ============================================================
  # Test 3 -- delta_abs_tolerance (call, type, shape; the values themselves are
  # test_effect_size_abs_tolerance_*), a negative one rejected.
  # ============================================================
  result_tol <- check_effect_size_plateau_condition(c(0.0012, 0.0012), c(0.0010, 0.0010), has_previous = TRUE,
                                                    delta_median_threshold = 0.05, delta_max_threshold = 0.10,
                                                    delta_epsilon = 1.0e-10, delta_min_consecutive_transitions = 2L,
                                                    n_consecutive_ok = 0L, delta_abs_tolerance = 0.0002)
  assert_equal_int(length(result_tol$delta), 2L, "Test 3 failed: delta of length n_studies")
  assert_true(is.numeric(result_tol$delta) && all(is.finite(result_tol$delta)) && all(result_tol$delta >= 0),
              "Test 3 failed: delta finite and >= 0")
  result_no_tol <- check_effect_size_plateau_condition(c(0.0012, 0.0012), c(0.0010, 0.0010), has_previous = TRUE,
                                                       delta_median_threshold = 0.05, delta_max_threshold = 0.10,
                                                       delta_epsilon = 1.0e-10, delta_min_consecutive_transitions = 2L,
                                                       n_consecutive_ok = 0L)
  assert_equal_int(length(result_no_tol$delta), 2L, "Test 3 failed: delta of length n_studies without the tolerance")
  assert_error(check_effect_size_plateau_condition(c(0.1), c(0.1), has_previous = TRUE,
                                                   delta_median_threshold = 0.05, delta_max_threshold = 0.10,
                                                   delta_epsilon = 1.0e-10, delta_min_consecutive_transitions = 2L,
                                                   n_consecutive_ok = 0L, delta_abs_tolerance = -1.0),
               "Test 3 failed: expected a negative delta_abs_tolerance to be rejected", ERR_INVALID_INPUT)
}

test_create_mean_pmf <- function() {
  # ============================================================
  # Test 1 -- two studies, hand-computed (test_create_mean_pmf_two_studies_hand_computed):
  # mean_pmf = ([0.3,0.7]+[0.5,0.5])/2 = [0.4,0.6], mean_pmf_counts=[8,12], included=20.
  # ============================================================
  pmfs <- array(0, dim = c(1, 2, 2))
  pmfs[1, , 1] <- c(0.3, 0.7)
  pmfs[1, , 2] <- c(0.5, 0.5)
  counts <- array(0L, dim = c(1, 2, 2))
  counts[1, , 1] <- c(3L, 7L)
  counts[1, , 2] <- c(5L, 5L)
  included_n_reps <- matrix(c(10L, 10L), nrow = 1)

  result <- create_mean_pmf(pmfs, counts, included_n_reps)
  assert_equal_numeric(result$mean_pmf[1, ], c(0.4, 0.6), TOL, "Test 1 failed: mean_pmf mismatch")
  assert_equal_int(as.integer(result$mean_pmf_counts[1, ]), c(8L, 12L), "Test 1 failed: mean_pmf_counts mismatch")
  assert_equal_int(as.integer(result$mean_pmf_included_n_reps[1]), 20L,
                    "Test 1 failed: mean_pmf_included_n_reps mismatch")

  # ============================================================
  # Test 2 -- three studies, DOCUMENTS THE KNOWN LIMITATION: averages over ALL studies including
  # self, not a true leave-one-out background (test_create_mean_pmf_three_studies_includes_self).
  # ============================================================
  pmfs3 <- array(0, dim = c(1, 2, 3))
  pmfs3[1, , 1] <- c(0.2, 0.8)
  pmfs3[1, , 2] <- c(0.4, 0.6)
  pmfs3[1, , 3] <- c(0.9, 0.1)
  counts3 <- array(0L, dim = c(1, 2, 3))
  counts3[1, , 1] <- c(2L, 8L)
  counts3[1, , 2] <- c(4L, 6L)
  counts3[1, , 3] <- c(9L, 1L)
  included_n_reps3 <- matrix(c(10L, 10L, 10L), nrow = 1)

  result3 <- create_mean_pmf(pmfs3, counts3, included_n_reps3)
  assert_equal_numeric(result3$mean_pmf[1, ], c(0.5, 0.5), TOL, "Test 2 failed: mean_pmf mismatch")
  assert_equal_int(as.integer(result3$mean_pmf_counts[1, ]), c(15L, 15L), "Test 2 failed: mean_pmf_counts mismatch")
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
  pmfs <- array(0, dim = c(1, 2, 2))
  pmfs[1, , 1] <- c(0.3, 0.7)
  pmfs[1, , 2] <- c(0.5, 0.5)
  counts <- array(0L, dim = c(1, 2, 2))
  counts[1, , 1] <- c(3L, 7L)
  counts[1, , 2] <- c(5L, 5L)
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
  counts <- bootstrap_test_counts()
  mean_pmf_included_n_reps <- c(10L, 10L)
  included_n_reps <- matrix(c(4L, 4L, 6L, 6L), nrow = 2)

  ci_first <- matrix(c(0.2, 0.2, 0.3, 0.3), nrow = 2)
  ci_first <- bootstrap_histogram(15L, counts, mean_pmf_included_n_reps, included_n_reps, ci_first,
                                   random_seed = 123L)

  ci_second <- matrix(c(0.2, 0.2, 0.3, 0.3), nrow = 2)
  ci_second <- bootstrap_histogram(15L, counts, mean_pmf_included_n_reps, included_n_reps, ci_second,
                                    random_seed = 123L)

  assert_true(identical(ci_first, ci_second), "same random_seed twice must give identical confidence_interval")

  # ============================================================
  # Test 2 -- degenerate single-bin histogram: JSD is always exactly 0.0, so the CI must
  # collapse to [0.0, 0.0] (test_bootstrap_histogram_degenerate_single_bin_zero_ci).
  # ============================================================
  mean_pmf_included_n_reps_deg <- c(10L)
  included_n_reps_deg <- matrix(c(5L, 5L), nrow = 1)
  confidence_interval_deg <- matrix(0.0, nrow = 2, ncol = 2)

  counts_deg <- array(c(5L, 5L), dim = c(1, 1, 2))
  confidence_interval_deg <- bootstrap_histogram(5L, counts_deg, mean_pmf_included_n_reps_deg,
                                                  included_n_reps_deg, confidence_interval_deg, random_seed = 7L)

  assert_true(identical(confidence_interval_deg, matrix(0.0, nrow = 2, ncol = 2)),
              "a single-bin histogram always has JSD=0, so the CI must collapse to [0,0]")
}

test_bootstrap_histogram_rejects_draws_from_empty_pool <- function() {
  # Each study draws 3 residuals from its own counts, which are empty (c(0, 0)) at the point, so
  # the call fails with invalid input, not tied to any one argument
  # (test_bootstrap_histogram_rejects_draws_from_empty_pool in the Fortran suite).
  counts <- array(0L, dim = c(1, 2, 2))
  mean_pmf_included_n_reps <- c(6L)
  included_n_reps <- matrix(c(3L, 3L), nrow = 1)

  call <- function() {
    bootstrap_histogram(3L, counts, mean_pmf_included_n_reps, included_n_reps,
                        matrix(0.0, nrow = 2, ncol = 2), random_seed = 7L)
  }
  assert_error(call(), "expected ERR_INVALID_INPUT for draws from an empty pool", ERR_INVALID_INPUT)
  err <- tryCatch({ call(); NULL }, error = function(e) e)
  assert_true(is.na(err$argument), "expected no argument named")
}

test_create_mean_pmf_rejects_zero_included_n_reps <- function() {
  # Every study needs at least one included residual at every reference point, so a single 0 in
  # included_n_reps is rejected on that argument
  # (test_create_mean_pmf_rejects_zero_included_n_reps in the Fortran suite).
  pmfs <- array(0.0, dim = c(2, 2, 2))
  pmfs[1, , 1] <- c(0.5, 0.5)
  pmfs[2, , 1] <- c(0.25, 0.75)
  pmfs[1, , 2] <- c(1.0, 0.0)
  counts <- array(0L, dim = c(2, 2, 2))
  counts[1, , 1] <- c(1L, 1L)
  counts[2, , 1] <- c(1L, 3L)
  counts[1, , 2] <- c(2L, 0L)
  included_n_reps <- matrix(c(2L, 4L, 2L, 0L), nrow = 2)

  call <- function() {
    create_mean_pmf(pmfs, counts, included_n_reps)
  }
  assert_error(call(), "expected ERR_INVALID_INPUT for a study without a residual at a point", ERR_INVALID_INPUT)
  err <- tryCatch({ call(); NULL }, error = function(e) e)
  assert_true(identical(err$argument, "included_n_reps"), "expected the error on included_n_reps")
}

test_bootstrap_histogram_rejects_zero_included_n_reps <- function() {
  # Every study needs at least one included residual at every reference point, so a single 0 in
  # included_n_reps is rejected on that argument, even where the counts otherwise add up
  # (test_bootstrap_histogram_rejects_zero_included_n_reps in the Fortran suite).
  counts <- array(0L, dim = c(1, 2, 2))
  counts[1, , 2] <- c(2L, 2L)
  mean_pmf_included_n_reps <- c(4L)
  included_n_reps <- matrix(c(0L, 4L), nrow = 1)

  call <- function() {
    bootstrap_histogram(3L, counts, mean_pmf_included_n_reps, included_n_reps,
                        matrix(0.5, nrow = 2, ncol = 2), random_seed = 7L)
  }
  assert_error(call(), "expected ERR_INVALID_INPUT for a study without a residual at a point", ERR_INVALID_INPUT)
  err <- tryCatch({ call(); NULL }, error = function(e) e)
  assert_true(identical(err$argument, "included_n_reps"), "expected the error on included_n_reps")
}

# Test helper: the studies' own counts (points x bins x studies) of the bootstrap fixture: two
# points, three bins, studies of 4 and 6 residuals (pool totals 10 and 10).
bootstrap_test_counts <- function() {
  counts <- array(0L, dim = c(2, 3, 2))
  counts[1, , 1] <- c(2L, 1L, 1L)
  counts[1, , 2] <- c(3L, 2L, 1L)
  counts[2, , 1] <- c(1L, 1L, 2L)
  counts[2, , 2] <- c(1L, 1L, 4L)
  counts
}

test_bootstrap_histogram_divergence_mode <- function() {
  # Call + shape for both divergence modes (the numbers are checked in the Fortran suite): each
  # returns a 2 x n_studies numeric interval with lower <= upper; omitting the mode equals the
  # default "divergence_pairwise_k2" exactly, which with two studies gives both the same interval.
  # An unknown mode is rejected before Fortran; per-study counts that do not add up, and a pool
  # total that is not the studies' sum, are invalid input.
  counts <- bootstrap_test_counts()
  mean_pmf_included_n_reps <- c(10L, 10L)
  included_n_reps <- matrix(c(4L, 4L, 6L, 6L), nrow = 2)
  seed_ci <- matrix(c(1.0, 0.0, 1.0, 0.0), nrow = 2)

  results <- list(
    default = bootstrap_histogram(20L, counts, mean_pmf_included_n_reps, included_n_reps, seed_ci,
                                  random_seed = 3L),
    symmetric = bootstrap_histogram(20L, counts, mean_pmf_included_n_reps, included_n_reps, seed_ci,
                                    random_seed = 3L, divergence_mode = "divergence_symmetric"),
    pairwise = bootstrap_histogram(20L, counts, mean_pmf_included_n_reps, included_n_reps, seed_ci,
                                   random_seed = 3L, divergence_mode = "divergence_pairwise_k2"))
  for (name in names(results)) {
    ci <- results[[name]]
    assert_true(is.numeric(ci) && identical(dim(ci), c(2L, 2L)), paste(name, ": numeric 2 x 2"))
    assert_true(all(is.finite(ci)) && all(ci[1, ] <= ci[2, ]), paste(name, ": lower <= upper"))
  }
  assert_true(identical(results$default, results$pairwise), "the default is divergence_pairwise_k2")
  assert_true(identical(results$pairwise[, 1], results$pairwise[, 2]), "pairwise: both studies alike")

  assert_error(bootstrap_histogram(3L, counts, mean_pmf_included_n_reps, included_n_reps,
                                   matrix(0.5, nrow = 2, ncol = 2), divergence_mode = "no_such_mode"),
               "expected an unknown divergence_mode to be rejected before Fortran")
  bad_counts <- counts
  bad_counts[1, 1, 1] <- 0L
  assert_error(bootstrap_histogram(3L, bad_counts, mean_pmf_included_n_reps, included_n_reps,
                                   matrix(0.5, nrow = 2, ncol = 2)),
               "expected ERR_INVALID_INPUT for per-study counts not adding up", ERR_INVALID_INPUT)
  assert_error(bootstrap_histogram(3L, counts, c(10L, 9L), included_n_reps, matrix(0.5, nrow = 2, ncol = 2)),
               "expected ERR_INVALID_INPUT for a pool total that is not the studies' sum", ERR_INVALID_INPUT)
}

test_calc_js_comp_test_candidate_bounds <- function() {
  # The bounds of the grid of test_generate_js_comp_test_candidates (8742 genes): the first
  # n_points 374 and floor(8742/(0.25*200)) = 174 neighbors at the floor; with
  # grid_min_points = 300 the neighbor bound is 116, with grid_points_factor = 8 the first
  # n_points is ceiling(8*sqrt(8742)) = 748.
  bounds <- calc_js_comp_test_candidate_bounds(8742L)
  assert_equal_int(as.integer(bounds$max_n_points_candidate), 374L, "expected max_n_points_candidate=374")
  assert_equal_int(as.integer(bounds$max_n_neighbors_candidate), 174L, "expected max_n_neighbors_candidate=174")

  candidates <- generate_js_comp_test_candidates(8742L)
  assert_true(bounds$max_n_neighbors_candidate >= max(candidates[2, ]),
              "max_n_neighbors_candidate must be a safe upper bound on the grid's n_neighbors")

  floor_300 <- calc_js_comp_test_candidate_bounds(8742L, grid_min_points = 300L)
  assert_equal_int(as.integer(c(floor_300$max_n_points_candidate, floor_300$max_n_neighbors_candidate)),
                   c(374L, 116L), "grid_min_points = 300")
  factor_8 <- calc_js_comp_test_candidate_bounds(8742L, grid_points_factor = 8.0)
  assert_equal_int(as.integer(factor_8$max_n_points_candidate), 748L, "grid_points_factor = 8")

  assert_error(calc_js_comp_test_candidate_bounds(0L),
               "expected ERR_INVALID_INPUT for max_n_genes_all_studies=0", ERR_INVALID_INPUT)
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

  # An out-of-range gene index must raise, confirming the wrapper's gene-index bounds check is wired through.
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

  # The routine validates its own arguments (test_calc_js_comp_test_n_top_k_jsds_validation).
  assert_error(calc_js_comp_test_n_top_k_jsds(0L),
               "expected ERR_INVALID_INPUT for n_bootstraps=0", ERR_INVALID_INPUT)
  assert_error(calc_js_comp_test_n_top_k_jsds(100L, two_sided_bootstrapping_significance_level = 100.5),
               "expected ERR_INVALID_INPUT for a significance level above 100", ERR_INVALID_INPUT)
}

test_run_js_comp_test_two_studies_hand_traceable <- function() {
  # Call-ability/shape check only (Fortran_Coding_Guides.pdf Sec 17.1) -- the same fixture as the
  # Fortran suite's own test_run_js_comp_test_two_studies_hand_traceable, which is where the exact
  # hand-derived closed-form JSD assertions live. Issue #187 removed the mandatory scalar `n_bins`
  # input (the bin count is now decided per reference point by an internal occupancy search) and
  # added `n_bins_per_point`/`max_n_bins_per_point` plus several per-point diagnostics;
  # `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts` are now always sized to the fixed 256-bin ceiling
  # (MAX_N_BINS) in their second dimension, of which only the leading `max_n_bins_per_point` bins
  # are meaningful -- a caller must slice `[, seq_len(max_n_bins_per_point), ]` themselves.
  gene_means <- matrix(c(1.0, 5.0, 1.0, 5.0), nrow = 2)

  residuals <- array(NaN, dim = c(2, 2, 2))
  residuals[, 1, 1] <- c(-1.0, 1.0)
  residuals[, 1, 2] <- c(-3.0, 3.0)

  x_star <- c(1.0)

  # min_residuals_per_bin=1 (its default is 10, which this tiny 4-residual fixture could never
  # satisfy) so the occupancy search actually succeeds, matching the Fortran fixture's own
  # override -- exercised here for call-ability, not to reproduce that test's exact numbers.
  result <- run_js_comp_test(n_neighbors = 1L, gene_means = gene_means,
                              residuals = residuals, x_star = x_star, min_residuals_per_bin = 1L)

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
  assert_equal_int(dim(result$pmfs), c(n_points, 256L, n_studies), "pmfs shape")
  assert_equal_int(dim(result$counts), c(n_points, 256L, n_studies), "counts shape")
  assert_equal_int(dim(result$included_n_reps), c(n_points, n_studies), "included_n_reps shape")
  assert_equal_int(dim(result$mean_pmf), c(n_points, 256L), "mean_pmf shape")
  assert_equal_int(dim(result$mean_pmf_counts), c(n_points, 256L), "mean_pmf_counts shape")
  assert_equal_int(length(result$mean_pmf_included_n_reps), n_points, "mean_pmf_included_n_reps shape")
  assert_equal_int(dim(result$js_divergences), c(n_points, n_studies), "js_divergences shape")
  assert_equal_int(dim(result$weights), c(n_points, n_studies), "weights shape")
  assert_equal_int(length(result$global_js_divergence), n_studies, "global_js_divergence shape")
  assert_true(is.null(result$p_values), "the final run has no permutation test (see run_js_comp_test_permutation)")

  # A caller must be able to slice down to the meaningful leading bins.
  max_n_bins <- as.integer(result$max_n_bins_per_point)
  assert_equal_int(dim(result$pmfs[, seq_len(max_n_bins), , drop = FALSE]), c(n_points, max_n_bins, n_studies),
                    "pmfs slices down to max_n_bins_per_point")

  assert_equal_int(as.integer(result$neighborhood_indices[, 1, 1]), 1L, "study 1 neighbor should be gene 1")
  assert_equal_int(as.integer(result$neighborhood_indices[, 1, 2]), 1L, "study 2 neighbor should be gene 1")

  # Two studies with the default divergence_mode (pairwise): both studies get the same value.
  assert_true(result$global_js_divergence[1] == result$global_js_divergence[2],
              "pairwise default: both studies share one global JSD")

  # bin_count_mode: every mode string is accepted (type/shape only; the bin counts themselves are
  # the Fortran suite's job), an unknown one is rejected before Fortran.
  for (mode in BIN_COUNT_MODES) {
    result_mode <- run_js_comp_test(n_neighbors = 1L, gene_means = gene_means,
                                    residuals = residuals, x_star = x_star,
                                    min_residuals_per_bin = 1L, bin_count_mode = mode)
    assert_equal_int(length(result_mode$n_bins_per_point), n_points, paste(mode, "n_bins_per_point length"))
    assert_true(is.logical(result_mode$occupancy_failed), paste(mode, "occupancy_failed must be logical"))
    assert_equal_int(dim(result_mode$pmfs), c(n_points, 256L, n_studies), paste(mode, "pmfs shape"))
    assert_equal_int(length(result_mode$global_js_divergence), n_studies, paste(mode, "global_js_divergence length"))
  }
  assert_error(run_js_comp_test(n_neighbors = 1L, gene_means = gene_means,
                                residuals = residuals, x_star = x_star,
                                bin_count_mode = "nonsense"),
               "expected an unknown bin_count_mode string to be rejected before Fortran")

  # divergence_mode: both strings accepted (type/shape), omitting it equals the default
  # "divergence_pairwise_k2" exactly, an unknown string is rejected before Fortran.
  for (mode in c("divergence_symmetric", "divergence_pairwise_k2")) {
    result_mode <- run_js_comp_test(n_neighbors = 1L, gene_means = gene_means, residuals = residuals,
                                    x_star = x_star, min_residuals_per_bin = 1L, divergence_mode = mode)
    assert_equal_int(length(result_mode$global_js_divergence), n_studies, paste(mode, "global_js_divergence length"))
    assert_equal_int(dim(result_mode$js_divergences), c(n_points, n_studies), paste(mode, "js_divergences shape"))
    if (mode == "divergence_pairwise_k2")
      assert_true(identical(result_mode$global_js_divergence, result$global_js_divergence), "default is pairwise")
  }
  assert_error(run_js_comp_test(n_neighbors = 1L, gene_means = gene_means, residuals = residuals, x_star = x_star,
                                min_residuals_per_bin = 1L, divergence_mode = "nonsense"),
               "expected an unknown divergence_mode string to be rejected before Fortran")
}

test_run_js_comp_test_expert_takes_the_callers_permutations <- function() {
  # run_js_comp_test sorts each study's gene means itself; run_js_comp_test_expert takes that order
  # from the caller (one column per study). Fixture 1, no ties: the expert tier with the ascending
  # order (order(), the only ascending order without ties) gives exactly what the plain call gives,
  # every output. Fixture 2: study 1's genes 1 and 2 tie at 1.0 and x_star = 0.5 lies below every
  # mean, so with one neighbor study 1's neighborhood is the first gene of the caller's order:
  # c(1, 2, 3) takes gene 1, c(2, 1, 3) gene 2
  # (test_run_js_comp_test_expert_uses_caller_tie_order in the Fortran suite).
  gene_means <- matrix(c(1.0, 5.0, 1.0, 5.0), nrow = 2)
  residuals <- array(NaN, dim = c(2, 2, 2))
  residuals[, 1, 1] <- c(-1.0, 1.0)
  residuals[, 1, 2] <- c(-3.0, 3.0)
  x_star <- c(1.0)

  plain <- run_js_comp_test(n_neighbors = 1L, gene_means = gene_means, residuals = residuals, x_star = x_star,
                            min_residuals_per_bin = 1L)
  ascending <- apply(gene_means, 2, order)
  storage.mode(ascending) <- "integer"
  expert <- run_js_comp_test_expert(n_neighbors = 1L, gene_means = gene_means, gene_means_perms = ascending,
                                    residuals = residuals, x_star = x_star, min_residuals_per_bin = 1L)
  assert_true(identical(sort(names(expert)), sort(names(plain))), "both tiers return the same outputs")
  for (key in names(plain)) {
    assert_true(identical(expert[[key]], plain[[key]]),
                paste(key, "expert tier with the ascending order differs from the plain call"))
  }

  tied_means <- matrix(c(1.0, 1.0, 5.0, 1.0, 5.0, 9.0), nrow = 3)
  tied_residuals <- array(0.0, dim = c(2, 3, 2))
  tied_residuals[, 1, 1] <- c(-1.0, 1.0)
  tied_residuals[, 2, 1] <- c(-3.0, 3.0)
  tied_residuals[, 3, 1] <- c(0.5, -0.5)
  tied_residuals[, 1, 2] <- c(-2.0, 2.0)
  tied_residuals[, 2, 2] <- c(0.25, -0.25)
  tied_residuals[, 3, 2] <- c(0.75, -0.75)
  for (case in list(list(order = c(1L, 2L, 3L), gene = 1L), list(order = c(2L, 1L, 3L), gene = 2L))) {
    perms <- cbind(case$order, c(1L, 2L, 3L))
    result <- run_js_comp_test_expert(n_neighbors = 1L, gene_means = tied_means, gene_means_perms = perms,
                                      residuals = tied_residuals, x_star = c(0.5), min_residuals_per_bin = 1L)
    label <- paste0("order c(", paste(case$order, collapse = ", "), "): ")
    assert_equal_int(as.integer(result$neighborhood_indices[1, 1, 1]), case$gene,
                     paste0(label, "study 1 should take gene ", case$gene))
    assert_equal_int(as.integer(result$neighborhood_indices[1, 1, 2]), 1L, paste0(label, "study 2 should take gene 1"))
  }
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
  # the FIRST (finest-resolution) candidate, (n_points, n_neighbors)=(200, 40): the default grid of
  # 2000 genes is that single candidate (ceiling(4*sqrt(2000)) = 179 clamped up to the floor 200,
  # floor(2000/(0.25*200)) = 40), and resets
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

  assert_equal_int(as.integer(result$n_points), 200L, "expected n_points=200")
  assert_equal_int(as.integer(result$n_neighbors), 40L, "expected n_neighbors=40")
  assert_false(result$plateau_established, "no candidate ever plateaus here")
  assert_equal_numeric(result$best_candidate_pair_confidence_interval[, 1], c(-1.0, -1.0), TOL, "study 1 CI reset")
  assert_equal_numeric(result$best_candidate_pair_confidence_interval[, 2], c(-1.0, -1.0), TOL, "study 2 CI reset")

  # ============================================================
  # Test 2 -- the GAMMA-decay grid collapses to exactly ONE candidate at
  # max_n_genes_all_studies=100 ((n_points, n_neighbors)=(200, 2): the floor 200 and
  # floor(100/(0.25*200)) = 2); with both admissibility
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

  assert_equal_int(as.integer(result2$n_points), 200L, "expected n_points=200")
  assert_equal_int(as.integer(result2$n_neighbors), 2L, "expected n_neighbors=2")
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
                "trace_delta_median", "trace_delta_max", "plateau_established", "trace_ci_overlap")) {
    assert_true(!is.null(result3[[key]]), paste0("missing expected output '", key, "'"))
  }
  assert_true(is.logical(result3$plateau_established), "plateau_established should be logical")
  n_admissible <- length(result3$trace_n_points)
  assert_true(n_admissible >= 1L, "expected at least one admissible candidate")
  assert_true(ncol(result3$trace_ci_width) == n_admissible, "trace_ci_width has n_admissible columns")
  assert_true(all(dim(result3$trace_ci_width_relative) == dim(result3$trace_ci_width)),
              "trace_ci_width_relative shape matches trace_ci_width")
  assert_true(all(dim(result3$trace_ci_overlap) == dim(result3$trace_ci_width)),
              "trace_ci_overlap shape matches trace_ci_width")
  # The first (and here, only) admissible candidate has no predecessor to diff against.
  assert_equal_numeric(result3$trace_delta[, 1], c(-1.0, -1.0), TOL, "first candidate's delta is the -1.0 sentinel")
  # Likewise no previous best to overlap with: the same -1.0 sentinel.
  assert_equal_numeric(result3$trace_ci_overlap[, 1], c(-1.0, -1.0), TOL, "first candidate's overlap is the -1.0 sentinel")

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
  # The n_points leading entries of n_bins_per_point must match the sole admissible candidate's
  # own trace_selected_n_bins column (test_param_search_final_n_bins_matches_selected_trace_column
  # is the Fortran suite's rigorous version of this same wiring check).
  n_points_2 <- result2$n_points
  assert_true(all(result2$n_bins_per_point[1:n_points_2] == result2$trace_selected_n_bins[1:n_points_2, 1]),
              "n_bins_per_point matches the sole admissible candidate's trace column")

  # ============================================================
  # Test 5 -- bin_count_mode: every mode string is accepted (type/shape only; the Stone bin counts
  # themselves are the Fortran suite's job, test_param_search_stone_matches_standalone), an
  # unknown one rejected. Reuses Test 2's single-candidate setup.
  # ============================================================
  for (mode in BIN_COUNT_MODES) {
    result_mode <- run_js_comp_test_parameter_search(gene_means_2, residuals_2,
                                                     n_bootstraps = 5L, join_method = "join_min",
                                                     max_n_points_candidate = bounds_2$max_n_points_candidate,
                                                     max_n_neighbors_candidate = bounds_2$max_n_neighbors_candidate,
                                                     min_residuals_per_bin = 0L, min_neighbor_overlap = 0.0,
                                                     random_seed = 1L, bin_count_mode = mode)
    n_points_mode <- result_mode$n_points
    assert_true(is.integer(result_mode$n_bins_per_point), paste(mode, "n_bins_per_point must be integer"))
    assert_equal_int(length(result_mode$n_bins_per_point), bounds_2$max_n_points_candidate,
                     paste(mode, "n_bins_per_point length"))
    assert_true(all(result_mode$n_bins_per_point[1:n_points_mode] >= 3 &
                    result_mode$n_bins_per_point[1:n_points_mode] <= 120),
                paste(mode, "bin counts within the default [m_min, m_max]"))
    assert_true(is.logical(result_mode$trace_occupancy_failed), paste(mode, "trace_occupancy_failed must be logical"))
    assert_true(ncol(result_mode$trace_selected_n_bins) >= 1L, paste(mode, "at least one traced candidate"))
    assert_true(all(dim(result_mode$trace_occupancy_failed) == dim(result_mode$trace_selected_n_bins)),
                paste(mode, "trace shapes"))
  }
  assert_error(run_js_comp_test_parameter_search(gene_means_2, residuals_2,
                                                 n_bootstraps = 5L, join_method = "join_min",
                                                 max_n_points_candidate = bounds_2$max_n_points_candidate,
                                                 max_n_neighbors_candidate = bounds_2$max_n_neighbors_candidate,
                                                 bin_count_mode = "nonsense"),
               "expected an unknown bin_count_mode string to be rejected before Fortran")

  # ============================================================
  # Test 6 -- neighborhood_size_mode: every mode string is accepted (type/shape only; the Fortran
  # suite's test_param_search_neighborhood_size_mode checks the grid it walks).
  # ============================================================
  for (mode in c("neighborhood_size_large", "neighborhood_size_small")) {
    result_size <- run_js_comp_test_parameter_search(gene_means, residuals, n_bootstraps = 10L, join_method = "join_min",
                                                     max_n_points_candidate = bounds$max_n_points_candidate,
                                                     max_n_neighbors_candidate = bounds$max_n_neighbors_candidate,
                                                     min_residuals_per_bin = 1000000L, random_seed = 1L,
                                                     neighborhood_size_mode = mode)
    assert_true(result_size$n_points >= 1L, paste(mode, "n_points"))
    assert_true(result_size$n_neighbors >= 1L, paste(mode, "n_neighbors"))
    assert_equal_int(length(result_size$n_bins_per_point), bounds$max_n_points_candidate, paste(mode, "length"))
  }
  assert_error(run_js_comp_test_parameter_search(gene_means, residuals, n_bootstraps = 10L, join_method = "join_min",
                                                 max_n_points_candidate = bounds$max_n_points_candidate,
                                                 max_n_neighbors_candidate = bounds$max_n_neighbors_candidate,
                                                 neighborhood_size_mode = "nonsense"),
               "expected an unknown neighborhood_size_mode string to be rejected before Fortran")

  # ============================================================
  # Test 7 -- delta_abs_tolerance is accepted (call + shape; test_param_search_delta_abs_tolerance
  # is the Fortran check), a negative one rejected. Reuses Test 2's single-candidate setup.
  # ============================================================
  result_tol <- run_js_comp_test_parameter_search(gene_means_2, residuals_2, n_bootstraps = 5L, join_method = "join_min",
                                                  max_n_points_candidate = bounds_2$max_n_points_candidate,
                                                  max_n_neighbors_candidate = bounds_2$max_n_neighbors_candidate,
                                                  min_residuals_per_bin = 0L, min_neighbor_overlap = 0.0,
                                                  random_seed = 1L, delta_abs_tolerance = 0.001)
  assert_true(result_tol$n_points >= 1L && result_tol$n_neighbors >= 1L, "tolerance: a candidate")
  assert_equal_int(dim(result_tol$trace_delta), c(2L, 1L), "tolerance: trace_delta shape")
  assert_error(run_js_comp_test_parameter_search(gene_means_2, residuals_2, n_bootstraps = 5L, join_method = "join_min",
                                                 max_n_points_candidate = bounds_2$max_n_points_candidate,
                                                 max_n_neighbors_candidate = bounds_2$max_n_neighbors_candidate,
                                                 delta_abs_tolerance = -0.1),
               "expected a negative delta_abs_tolerance to be rejected", ERR_INVALID_INPUT)

  # ============================================================
  # Test 8 -- divergence_mode: both strings accepted (call + shape; the values are the Fortran
  # suite's job), an unknown one rejected before Fortran. Reuses Test 2's setup.
  # ============================================================
  for (mode in c("divergence_symmetric", "divergence_pairwise_k2")) {
    result_div <- run_js_comp_test_parameter_search(gene_means_2, residuals_2, n_bootstraps = 5L,
                                                    join_method = "join_min",
                                                    max_n_points_candidate = bounds_2$max_n_points_candidate,
                                                    max_n_neighbors_candidate = bounds_2$max_n_neighbors_candidate,
                                                    min_residuals_per_bin = 0L, min_neighbor_overlap = 0.0,
                                                    random_seed = 1L, divergence_mode = mode)
    assert_equal_int(dim(result_div$best_candidate_pair_confidence_interval), c(2L, 2L), paste(mode, "CI shape"))
    assert_equal_int(dim(result_div$trace_global_js_divergence), c(2L, 1L), paste(mode, "trace shape"))
  }
  assert_error(run_js_comp_test_parameter_search(gene_means_2, residuals_2, n_bootstraps = 5L, join_method = "join_min",
                                                 max_n_points_candidate = bounds_2$max_n_points_candidate,
                                                 max_n_neighbors_candidate = bounds_2$max_n_neighbors_candidate,
                                                 divergence_mode = "no_such_mode"),
               "expected an unknown divergence_mode to be rejected before Fortran")

  # ============================================================
  # Test 9 -- grid_points_factor and grid_min_points: factor 6.7 and floor 240 give the 2000-gene
  # grid (300, 26), (240, 33); no candidate is admissible, so the first is returned. Capacities
  # sized for the default grid are too small for it: invalid input naming no argument.
  # ============================================================
  bounds_grid <- calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, grid_points_factor = 6.7,
                                                    grid_min_points = 240L)
  assert_equal_int(as.integer(c(bounds_grid$max_n_points_candidate, bounds_grid$max_n_neighbors_candidate)),
                   c(300L, 33L), "grid options: bounds")
  result_grid <- run_js_comp_test_parameter_search(gene_means, residuals, n_bootstraps = 10L, join_method = "join_min",
                                                   max_n_points_candidate = bounds_grid$max_n_points_candidate,
                                                   max_n_neighbors_candidate = bounds_grid$max_n_neighbors_candidate,
                                                   min_residuals_per_bin = 1000000L, random_seed = 1L,
                                                   grid_points_factor = 6.7, grid_min_points = 240L)
  assert_equal_int(as.integer(c(result_grid$n_points, result_grid$n_neighbors)), c(300L, 26L),
                   "grid options: the first candidate")
  undersized <- function() {
    run_js_comp_test_parameter_search(gene_means, residuals, n_bootstraps = 10L, join_method = "join_min",
                                      max_n_points_candidate = bounds$max_n_points_candidate,
                                      max_n_neighbors_candidate = bounds$max_n_neighbors_candidate,
                                      min_residuals_per_bin = 1000000L, random_seed = 1L,
                                      grid_points_factor = 6.7, grid_min_points = 240L)
  }
  assert_error(undersized(), "expected ERR_INVALID_INPUT for capacities below the optioned grid", ERR_INVALID_INPUT)
  err <- tryCatch({ undersized(); NULL }, error = function(e) e)
  assert_true(is.na(err$argument), "expected no argument named")
  assert_error(run_js_comp_test_parameter_search(gene_means, residuals, n_bootstraps = 10L, join_method = "join_min",
                                                 max_n_points_candidate = bounds$max_n_points_candidate,
                                                 max_n_neighbors_candidate = bounds$max_n_neighbors_candidate,
                                                 grid_points_factor = -1.0),
               "expected ERR_INVALID_INPUT for a negative grid_points_factor", ERR_INVALID_INPUT)
}

test_run_js_comp_test_parameter_search_rejects_undersized_bounds <- function() {
  # A max_n_points_candidate or max_n_neighbors_candidate one below what
  # calc_js_comp_test_candidate_bounds returns for the gene count (199 and 5 for 300 genes: the
  # default grid is the single candidate (200, 6)) is rejected with invalid input, checked by the
  # search against the grid options, instead of being written past
  # (test_param_search_rejects_undersized_max_n_points_candidate and
  # test_param_search_rejects_undersized_max_n_neighbors_candidate).
  n_genes <- 300L
  n_reps <- 3L
  n_studies <- 2L
  gene_means <- matrix(rep(as.numeric(seq_len(n_genes)), n_studies), nrow = n_genes, ncol = n_studies)
  residuals <- array(rep(c(-0.5, 0.0, 0.5), n_genes * n_studies), dim = c(n_reps, n_genes, n_studies))

  bounds <- calc_js_comp_test_candidate_bounds(n_genes)
  assert_error(run_js_comp_test_parameter_search(gene_means, residuals, n_bootstraps = 5L, join_method = "join_min",
                                                 max_n_points_candidate = bounds$max_n_points_candidate - 1L,
                                                 max_n_neighbors_candidate = bounds$max_n_neighbors_candidate,
                                                 min_neighbor_overlap = 0.0, random_seed = 1L),
               "expected ERR_INVALID_INPUT for max_n_points_candidate below its bound", ERR_INVALID_INPUT)
  assert_error(run_js_comp_test_parameter_search(gene_means, residuals, n_bootstraps = 5L, join_method = "join_min",
                                                 max_n_points_candidate = bounds$max_n_points_candidate,
                                                 max_n_neighbors_candidate = bounds$max_n_neighbors_candidate - 1L,
                                                 min_neighbor_overlap = 0.0, random_seed = 1L),
               "expected ERR_INVALID_INPUT for max_n_neighbors_candidate below its bound", ERR_INVALID_INPUT)
}

test_run_js_comp_test_rejects_reversed_residual_range_quantiles <- function() {
  # An upper residual-range quantile below the lower one (0.05 < 0.95) is rejected with invalid
  # input instead of collapsing every histogram into one bin
  # (test_run_js_comp_test_rejects_reversed_residual_range_quantiles in the Fortran suite).
  gene_means <- matrix(c(1.0, 5.0, 1.0, 5.0), nrow = 2)
  residuals <- array(NaN, dim = c(2, 2, 2))
  residuals[, 1, 1] <- c(-1.0, 1.0)
  residuals[, 1, 2] <- c(-3.0, 3.0)
  x_star <- c(1.0)

  call <- function() {
    run_js_comp_test(n_neighbors = 1L, gene_means = gene_means,
                     residuals = residuals, x_star = x_star,
                     min_residuals_per_bin = 1L, lower_residual_range_quantile = 0.95,
                     upper_residual_range_quantile = 0.05)
  }
  assert_error(call(), "expected ERR_INVALID_INPUT for upper_residual_range_quantile below the lower one",
               ERR_INVALID_INPUT)
  err <- tryCatch({ call(); NULL }, error = function(e) e)
  assert_true(identical(err$argument, "upper_residual_range_quantile"),
              "expected the error on upper_residual_range_quantile")
}

test_run_js_comp_test_rejects_nan_x_star <- function() {
  # A NaN reference point is rejected with ERR_NAN_INF on x_star, instead of silently taking the
  # lowest-mean genes as its neighborhood (test_run_js_comp_test_rejects_nan_x_star in the
  # Fortran suite).
  gene_means <- matrix(c(1.0, 5.0, 1.0, 5.0), nrow = 2)
  residuals <- array(NaN, dim = c(2, 2, 2))
  residuals[, 1, 1] <- c(-1.0, 1.0)
  residuals[, 1, 2] <- c(-3.0, 3.0)
  x_star <- c(1.0, NaN)

  call <- function() {
    run_js_comp_test(n_neighbors = 1L, gene_means = gene_means,
                     residuals = residuals, x_star = x_star,
                     min_residuals_per_bin = 1L)
  }
  assert_error(call(), "expected ERR_NAN_INF for a NaN reference point", ERR_NAN_INF)
  err <- tryCatch({ call(); NULL }, error = function(e) e)
  assert_true(identical(err$argument, "x_star"), "expected the error on x_star")
}

test_run_js_comp_test_rejects_empty_study <- function() {
  # A study whose gene means and residuals are all missing includes no residual at the reference
  # point, so the equal-weight consensus is undefined and the call fails with invalid input, not
  # tied to any one argument (test_run_js_comp_test_empty_study_sets_ierr in the Fortran suite).
  gene_means <- matrix(c(1.0, 5.0, NaN, NaN), nrow = 2)
  residuals <- array(NaN, dim = c(2, 2, 2))
  residuals[, 1, 1] <- c(-1.0, 1.0)
  residuals[, 2, 1] <- c(-3.0, 3.0)
  x_star <- c(1.0)

  call <- function() {
    run_js_comp_test(n_neighbors = 1L, gene_means = gene_means,
                     residuals = residuals, x_star = x_star,
                     min_residuals_per_bin = 1L)
  }
  assert_error(call(), "expected ERR_INVALID_INPUT for a study without any non-NaN residual", ERR_INVALID_INPUT)
  err <- tryCatch({ call(); NULL }, error = function(e) e)
  assert_true(is.na(err$argument), "expected no argument named")
}

test_generate_adaptive_js_comp_test_candidates <- function() {
  # Call-ability, return type/shape and non-error only -- the hand-derived sequences are the
  # Fortran suite's job (test_adaptive_candidates_*). G=2000, S=2 must not collapse to a single
  # candidate the way the fixed-k grid does at this size.
  candidates <- generate_adaptive_js_comp_test_candidates(2000L, 2L)
  assert_true(is.matrix(candidates), "expected a matrix")
  assert_true(is.integer(candidates), "expected an integer matrix")
  assert_equal_int(nrow(candidates), 3L, "expected 3 rows (k_start, k_step, k_max)")
  assert_true(ncol(candidates) >= 8 && ncol(candidates) <= 16, "expected 8..16 candidates")

  assert_error(generate_adaptive_js_comp_test_candidates(0L, 2L),
               "expected ERR_INVALID_INPUT for max_n_genes_all_studies=0", ERR_INVALID_INPUT)
  assert_error(generate_adaptive_js_comp_test_candidates(2000L, 0L),
               "expected ERR_INVALID_INPUT for n_studies=0", ERR_INVALID_INPUT)
}

test_calc_adaptive_js_comp_test_bounds <- function() {
  # Call-ability, return type and non-error only (numerics: test_adaptive_bounds_* in Fortran).
  max_n_points_candidate <- calc_adaptive_js_comp_test_bounds(2000L, 2L)
  assert_true(is.numeric(max_n_points_candidate) && length(max_n_points_candidate) == 1,
              "expected a single number")
  assert_true(max_n_points_candidate >= 1 && max_n_points_candidate <= 2000 * 2,
              "capacity must lie in [1, G*S]")

  assert_error(calc_adaptive_js_comp_test_bounds(0L, 2L),
               "expected ERR_INVALID_INPUT for max_n_genes_all_studies=0", ERR_INVALID_INPUT)
  assert_error(calc_adaptive_js_comp_test_bounds(2000L, 0L),
               "expected ERR_INVALID_INPUT for n_studies=0", ERR_INVALID_INPUT)
  # 50000 * 50000 = 2.5e9 does not fit a 32-bit integer (the pool size is a 32-bit extent downstream).
  assert_error(calc_adaptive_js_comp_test_bounds(50000L, 50000L),
               "expected ERR_INVALID_INPUT for a pool size beyond int32", ERR_INVALID_INPUT)
}

test_construct_adaptive_neighborhoods <- function() {
  # Call-ability, return types/shapes and non-error only -- the hand-derived constructions are
  # the Fortran suite's job (test_adaptive_*). Two studies of 6 genes, the second NaN-padded by
  # one gene, 10 replicates each; every output is trimmed to the same number of points.
  set.seed(217)
  gene_means <- matrix(runif(12, 0, 10), nrow = 6, ncol = 2)
  gene_means[6, 2] <- NaN
  residuals <- array(rnorm(120), dim = c(10, 6, 2))
  residuals[, 6, 2] <- NaN
  result <- construct_adaptive_neighborhoods(gene_means, residuals, 2L, 1L, 4L)
  assert_true(is.list(result), "expected a list")
  n_points <- length(result$x_star)
  assert_true(n_points >= 1 && n_points <= 12, "expected 1..G*S points")
  assert_true(is.matrix(result$pooled_neighborhood_range) && all(dim(result$pooled_neighborhood_range) == c(2, n_points)),
              "pooled_neighborhood_range must be 2 x n_points")
  assert_true(is.matrix(result$n_neighbors_per_point) && all(dim(result$n_neighbors_per_point) == c(2, n_points)),
              "n_neighbors_per_point must be n_studies x n_points")
  assert_equal_int(length(result$stop_reason), n_points, "stop_reason length")
  assert_equal_int(length(result$neighborhood_dispersion), n_points, "neighborhood_dispersion length")
  assert_equal_int(length(result$neighborhood_mad), n_points, "neighborhood_mad length")
  assert_true(all(result$stop_reason %in% c(STAT_ADAPTIVE_STOP_TAU, STAT_ADAPTIVE_STOP_K_MAX,
                                            STAT_ADAPTIVE_STOP_RESIDUAL_CAP, STAT_ADAPTIVE_STOP_EXHAUSTED,
                                            STAT_ADAPTIVE_STOP_ZERO_DISPERSION, STAT_ADAPTIVE_STOP_TOO_FEW_RESIDUALS)),
              "stop reasons are STAT_ADAPTIVE_STOP_* codes")
  assert_true(result$construction_status %in% c(ERR_OK, STAT_ADAPTIVE_CONSTRUCTION_TOO_FEW_MEANS,
                                                STAT_ADAPTIVE_CONSTRUCTION_EMPTY_STUDY_NEIGHBORHOOD,
                                                STAT_ADAPTIVE_CONSTRUCTION_TOO_FEW_RESIDUALS),
              "construction status is ERR_OK or a STAT_ADAPTIVE_CONSTRUCTION_* code")
  assert_true(result$max_n_neighbors >= 1, "max_n_neighbors >= 1")
  # structural: the per-study counts of a neighborhood partition its pooled range
  lengths <- result$pooled_neighborhood_range[2, ] - result$pooled_neighborhood_range[1, ] + 1L
  assert_true(all(colSums(result$n_neighbors_per_point) == lengths), "counts partition each range")

  optional <- construct_adaptive_neighborhoods(gene_means, residuals, 2L, 1L, 4L, tau = 0.5, mad_distance_factor = 0.5,
                                               max_pooled_residuals = 60L, min_study_neighbors = 1L)
  assert_true(length(optional$x_star) >= 1, "optionals accepted")

  assert_error(construct_adaptive_neighborhoods(gene_means, residuals, 3L, 1L, 2L),
               "expected ERR_INVALID_INPUT for k_max < k_start", ERR_INVALID_INPUT)
}

test_run_js_comp_test_adaptive <- function() {
  # Call-ability, return types/shapes and non-error only -- the hand-derived values and the
  # fixed-k equivalence are the Fortran suite's job (test_run_adaptive_*). construct first, then
  # run on its trimmed pooled ranges, as a caller does. Three interleaved studies of 8 genes,
  # the third NaN-padded by one gene, 6 replicates each.
  set.seed(2170)
  n_genes <- 8L
  n_studies <- 3L
  gene_means <- outer(seq_len(n_genes) - 1, 0.3 * (seq_len(n_studies) - 1), "+")
  gene_means[8, 3] <- NaN
  residuals <- array(rnorm(6 * n_genes * n_studies), dim = c(6, n_genes, n_studies))
  residuals[, 8, 3] <- NaN
  construction <- construct_adaptive_neighborhoods(gene_means, residuals, 6L, 3L, 12L)
  assert_true(construction$construction_status == ERR_OK, "expected a construction without empty studies")
  ranges <- construction$pooled_neighborhood_range
  n_points <- ncol(ranges)

  result <- run_js_comp_test_adaptive(gene_means, residuals, ranges, min_residuals_per_bin = 2L, m_min = 2L, m_max = 8L)
  assert_true(is.list(result), "expected a list")
  assert_true(is.matrix(result$n_neighbors_per_point) &&
                all(dim(result$n_neighbors_per_point) == c(n_studies, n_points)),
              "n_neighbors_per_point must be n_studies x n_points")
  assert_true(all(result$n_neighbors_per_point == construction$n_neighbors_per_point),
              "per-study counts equal the construction's")
  for (key in c("n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "occupancy_failed",
                "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy",
                "sturges_bins", "fd_bins", "mean_pmf_included_n_reps")) {
    assert_equal_int(length(result[[key]]), n_points, paste(key, "length"))
  }
  assert_true(is.numeric(result$max_n_bins_per_point), "max_n_bins_per_point must be numeric")
  assert_equal_int(dim(result$pmfs), c(n_points, 256L, n_studies), "pmfs shape")
  assert_equal_int(dim(result$counts), c(n_points, 256L, n_studies), "counts shape")
  assert_equal_int(dim(result$mean_pmf), c(n_points, 256L), "mean_pmf shape")
  assert_equal_int(dim(result$mean_pmf_counts), c(n_points, 256L), "mean_pmf_counts shape")
  for (key in c("included_n_reps", "js_divergences", "weights")) {
    assert_equal_int(dim(result[[key]]), c(n_points, n_studies), paste(key, "shape"))
  }
  assert_equal_int(length(result$global_js_divergence), n_studies, "global_js_divergence length")
  assert_true(is.null(result$p_values), "the final run has no permutation test")
  assert_true(all(is.finite(result$global_js_divergence)), "global JSD finite")

  # divergence_mode reaches Fortran: on the first two studies (with three or more both modes
  # coincide) the symmetric mode gives each study its own JSD against the consensus, the pairwise
  # mode one value shared by both, and omitting the mode equals the pairwise mode.
  gene_means_two <- gene_means[, 1:2]
  residuals_two <- residuals[, , 1:2]
  construction_two <- construct_adaptive_neighborhoods(gene_means_two, residuals_two, 6L, 3L, 12L)
  assert_true(construction_two$construction_status == ERR_OK, "two studies: expected a construction")
  ranges_two <- construction_two$pooled_neighborhood_range
  two_study_run <- function(...) {
    run_js_comp_test_adaptive(gene_means_two, residuals_two, ranges_two, min_residuals_per_bin = 2L, m_min = 2L,
                              m_max = 8L, ...)$global_js_divergence
  }
  pairwise <- two_study_run(divergence_mode = "divergence_pairwise_k2")
  assert_two_study_divergence_modes(two_study_run(divergence_mode = "divergence_symmetric"), pairwise,
                                    "run_js_comp_test_adaptive")
  assert_true(identical(two_study_run(), pairwise), "default is pairwise")
  assert_error(run_js_comp_test_adaptive(gene_means, residuals, ranges, divergence_mode = "nonsense"),
               "expected an unknown divergence_mode string to be rejected before Fortran")

  # bin_count_mode: every mode string is accepted (type/shape only), an unknown one rejected.
  for (mode in BIN_COUNT_MODES) {
    result_mode <- run_js_comp_test_adaptive(gene_means, residuals, ranges,
                                             min_residuals_per_bin = 2L, m_min = 2L, m_max = 8L, bin_count_mode = mode)
    assert_equal_int(length(result_mode$n_bins_per_point), n_points, paste(mode, "n_bins_per_point length"))
    assert_true(is.logical(result_mode$occupancy_failed), paste(mode, "occupancy_failed must be logical"))
    assert_true(all(result_mode$n_bins_per_point >= 2 & result_mode$n_bins_per_point <= 8),
                paste(mode, "bin counts within [m_min, m_max]"))
    assert_equal_int(length(result_mode$global_js_divergence), n_studies, paste(mode, "global_js_divergence length"))
  }
  assert_error(run_js_comp_test_adaptive(gene_means, residuals, ranges,
                                         bin_count_mode = "nonsense"),
               "expected an unknown bin_count_mode string to be rejected before Fortran")

  reversed_ranges <- ranges
  reversed_ranges[, 1] <- c(2L, 1L)
  assert_error(run_js_comp_test_adaptive(gene_means, residuals, reversed_ranges),
               "expected ERR_INVALID_INPUT for a reversed pooled range", ERR_INVALID_INPUT)
}

test_run_js_comp_test_adaptive_rejects_empty_study <- function() {
  # Study 2 has valid gene means, so both pooled ranges hold one of its genes, but all its
  # residuals are missing: it includes no residual at either reference point, so the
  # equal-weight consensus is undefined and the call fails with invalid input, not tied to any
  # one argument (test_run_adaptive_nan_residuals_only_study_sets_ierr in the Fortran suite).
  gene_means <- matrix(c(0.0, 2.0, 4.0, 1.0, 3.0, NaN), nrow = 3)
  residuals <- array(NaN, dim = c(2, 3, 2))
  residuals[, 1, 1] <- c(-1.0, 1.0)
  residuals[, 2, 1] <- c(0.5, NaN)
  residuals[, 3, 1] <- c(-0.5, 0.5)
  ranges <- matrix(c(1L, 3L, 4L, 5L), nrow = 2)

  call <- function() {
    run_js_comp_test_adaptive(gene_means, residuals, ranges, min_residuals_per_bin = 1L)
  }
  assert_error(call(), "expected ERR_INVALID_INPUT for a study without any non-NaN residual", ERR_INVALID_INPUT)
  err <- tryCatch({ call(); NULL }, error = function(e) e)
  assert_true(is.na(err$argument), "expected no argument named")
}

test_run_js_comp_test_adaptive_parameter_search <- function() {
  # Call-ability, return types/shapes and non-error only -- plateau selection, fallbacks,
  # candidate statuses and the construct + run round trip are the Fortran suite's job
  # (test_adaptive_search_*). Two interleaved studies of 500 genes, 3 replicates, so the adaptive
  # sequence has several candidates; gates relaxed so candidates get evaluated.
  set.seed(2171)
  n_genes <- 500L
  n_studies <- 2L
  gene_means <- outer(seq_len(n_genes), 0.5 * (seq_len(n_studies) - 1), "+")
  residuals <- array(rnorm(3 * n_genes * n_studies), dim = c(3, n_genes, n_studies))
  capacity <- calc_adaptive_js_comp_test_bounds(n_genes, n_studies)

  result <- run_js_comp_test_adaptive_parameter_search(gene_means, residuals, n_bootstraps = 10L,
                                                       join_method = "join_min", max_n_points_candidate = capacity,
                                                       min_residuals_per_bin = 0L, min_neighbor_overlap = 0.0,
                                                       plateau_mode = "plateau_both", tau = 0.1,
                                                       mad_distance_factor = 1.0, max_pooled_residuals = 0L,
                                                       min_study_neighbors = 1L, random_seed = 1L,
                                                       delta_abs_tolerance = 0.001, divergence_mode = "divergence_symmetric")
  assert_true(is.list(result), "expected a list")
  for (key in c("k_start", "k_step", "k_max", "n_points")) {
    assert_true(is.numeric(result[[key]]) && length(result[[key]]) == 1, paste(key, "must be a single number"))
  }
  assert_true(result$n_points >= 1 && result$n_points <= capacity, "n_points within the capacity")
  assert_true(result$k_start <= result$k_max, "k_start <= k_max")
  assert_true(is.logical(result$plateau_established), "plateau_established must be logical")
  assert_equal_int(length(result$n_bins_per_point), capacity, "n_bins_per_point length")
  assert_equal_int(dim(result$best_candidate_confidence_interval), c(2L, n_studies), "CI shape")

  # Two independent DM_RESULT_SIZE_IS counters: traces trimmed to the admissible candidates, the
  # candidate log to the candidates tried; both counters are dropped.
  assert_true(is.null(result$n_admissible_evaluated) && is.null(result$n_candidates_tried), "counters dropped")
  n_admissible <- length(result$trace_k_start)
  n_tried <- length(result$candidate_status)
  assert_true(n_admissible >= 1 && n_admissible <= n_tried && n_tried <= 16, "1 <= admissible <= tried <= 16")
  for (key in c("trace_k_step", "trace_k_max", "trace_n_points", "trace_delta_median", "trace_delta_max")) {
    assert_equal_int(length(result[[key]]), n_admissible, paste(key, "length"))
  }
  for (key in c("trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width",
                "trace_ci_width_relative", "trace_ci_overlap", "trace_delta")) {
    assert_equal_int(dim(result[[key]]), c(n_studies, n_admissible), paste(key, "shape"))
  }
  # The first admissible candidate has no previous best to overlap with: the -1.0 sentinel.
  assert_equal_numeric(as.numeric(result$trace_ci_overlap[, 1]), rep(-1.0, n_studies), TOL,
                        "first candidate's overlap is the -1.0 sentinel")
  for (key in c("trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals",
                "trace_shared_residual_range_low", "trace_shared_residual_range_high")) {
    assert_equal_int(dim(result[[key]]), c(capacity, n_admissible), paste(key, "shape"))
  }
  assert_equal_int(length(result$candidate_k_start), n_tried, "candidate_k_start length")
  assert_equal_int(length(result$candidate_n_points), n_tried, "candidate_n_points length")
  assert_true(all(result$candidate_status %in% c(ERR_OK, STAT_ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED,
                                                 STAT_ADAPTIVE_CANDIDATE_EMPTY_STUDY,
                                                 STAT_ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED,
                                                 STAT_ADAPTIVE_CANDIDATE_OVERLAP_FAILED,
                                                 STAT_ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED)),
              "status codes are ERR_OK or STAT_ADAPTIVE_CANDIDATE_* codes")

  assert_error(run_js_comp_test_adaptive_parameter_search(gene_means, residuals, n_bootstraps = 10L,
                                                          join_method = "join_min", max_n_points_candidate = 0L),
               "expected ERR_INVALID_INPUT for max_n_points_candidate=0", ERR_INVALID_INPUT)

  # bin_count_mode: every mode string is accepted (type/shape only), an unknown one rejected.
  for (mode in BIN_COUNT_MODES) {
    result_mode <- run_js_comp_test_adaptive_parameter_search(gene_means, residuals, n_bootstraps = 10L,
                                                              join_method = "join_min",
                                                              max_n_points_candidate = capacity,
                                                              min_residuals_per_bin = 0L, min_neighbor_overlap = 0.0,
                                                              random_seed = 1L, bin_count_mode = mode)
    assert_true(is.numeric(result_mode$n_points) && length(result_mode$n_points) == 1,
                paste(mode, "n_points must be a single number"))
    assert_equal_int(length(result_mode$n_bins_per_point), capacity, paste(mode, "n_bins_per_point length"))
    n_admissible_mode <- length(result_mode$trace_k_start)
    assert_equal_int(dim(result_mode$trace_selected_n_bins), c(capacity, n_admissible_mode),
                     paste(mode, "trace_selected_n_bins shape"))
    assert_true(is.logical(result_mode$trace_occupancy_failed), paste(mode, "trace_occupancy_failed must be logical"))
  }
  assert_error(run_js_comp_test_adaptive_parameter_search(gene_means, residuals, n_bootstraps = 10L,
                                                          join_method = "join_min", max_n_points_candidate = capacity,
                                                          bin_count_mode = "nonsense"),
               "expected an unknown bin_count_mode string to be rejected before Fortran")
  assert_error(run_js_comp_test_adaptive_parameter_search(gene_means, residuals, n_bootstraps = 10L,
                                                          join_method = "join_min", max_n_points_candidate = capacity,
                                                          divergence_mode = "nonsense"),
               "expected an unknown divergence_mode string to be rejected before Fortran")
}

# Three studies of 5, 4 and 4 samples (the last slot of studies 2 and 3 is NaN padding), 300
# genes spread over the expression range, as the Fortran suite's fixtures.
# Two studies: the symmetric divergence gives the two studies different values, the pairwise one a
# single shared value, and the two modes differ -- so a mode dropped by the binding fails.
assert_two_study_divergence_modes <- function(symmetric, pairwise, label) {
  assert_true(length(symmetric) == 2 && length(pairwise) == 2, paste(label, ": one value per study"))
  assert_true(symmetric[1] != symmetric[2], paste(label, ": symmetric, the two studies differ"))
  assert_true(pairwise[1] == pairwise[2], paste(label, ": pairwise, one shared value"))
  assert_true(symmetric[1] != pairwise[1], paste(label, ": the two modes give different values"))
}

label_test_expression <- function() {
  set.seed(9001)
  n_genes <- 300L
  n_reps <- c(5L, 4L, 4L)
  expression <- array(NaN, dim = c(5L, n_genes, 3L))
  centres <- 2.0 + 6.0 * (seq_len(n_genes) - 1) / (n_genes - 1)
  for (study in 1:3) {
    for (rep in seq_len(n_reps[study])) {
      expression[rep, , study] <- centres + rnorm(n_genes, 0.0, 0.5)
    }
  }
  expression
}

test_calc_js_comp_test_permutation_bounds <- function() {
  # Types first; then the exact values for the default and an optioned grid, so dropped grid options fail.
  bounds <- calc_js_comp_test_permutation_bounds(300L, 10L, 30L)
  for (key in c("max_n_points_candidate", "max_n_neighbors_candidate", "n_bootstrapping_top_k_jsds")) {
    assert_true(is.numeric(bounds[[key]]) && length(bounds[[key]]) == 1 && bounds[[key]] >= 1,
                paste(key, "must be a single number of at least 1"))
  }
  assert_error(calc_js_comp_test_permutation_bounds(300L, 10L, 301L),
               "expected ERR_INVALID_INPUT for n_neighbors above the gene count", ERR_INVALID_INPUT)
  # Exact values, by hand (G = 300, n_points = 10, n_neighbors = 30):
  # defaults: ceiling(4*sqrt(300)) = ceiling(69.28) = 70, clamped up to the floor 200; n_points_low = 200;
  #   points max(10, 200) = 200, neighbors max(30, min(300, floor(300/(0.25*200)) = 6)) = 30.
  # grid_points_factor = 20, grid_min_points = 300: ceiling(20*sqrt(300)) = ceiling(346.41) = 347,
  #   n_points_low = max(300, ceiling(0.2*347) = 70) = 300; points max(10, 347) = 347,
  #   neighbors max(30, floor(300/(0.25*300)) = 4) = 30. The heap size does not depend on the grid.
  assert_true(bounds$max_n_points_candidate == 200, "defaults: expected 200 points")
  assert_true(bounds$max_n_neighbors_candidate == 30, "defaults: expected 30 neighbors")
  optioned <- calc_js_comp_test_permutation_bounds(300L, 10L, 30L, grid_points_factor = 20.0, grid_min_points = 300L)
  assert_true(optioned$max_n_points_candidate == 347, "grid options: expected 347 points")
  assert_true(optioned$max_n_neighbors_candidate == 30, "grid options: expected 30 neighbors")
  assert_true(optioned$n_bootstrapping_top_k_jsds == bounds$n_bootstrapping_top_k_jsds,
              "grid options: the heap size does not depend on the grid")
  assert_error(calc_js_comp_test_permutation_bounds(300L, 10L, 30L, grid_min_points = 0L),
               "expected ERR_INVALID_INPUT for grid_min_points = 0", ERR_INVALID_INPUT)
}

test_calc_adaptive_js_comp_test_permutation_bounds <- function() {
  bounds <- calc_adaptive_js_comp_test_permutation_bounds(300L, 3L, 30L)
  for (key in c("max_n_points_candidate", "n_bootstrapping_top_k_jsds")) {
    assert_true(is.numeric(bounds[[key]]) && length(bounds[[key]]) == 1 && bounds[[key]] >= 1,
                paste(key, "must be a single number of at least 1"))
  }
  assert_error(calc_adaptive_js_comp_test_permutation_bounds(300L, 3L, 0L),
               "expected ERR_INVALID_INPUT for k_start = 0", ERR_INVALID_INPUT)
}

test_run_js_comp_test_on_sample_assignment <- function() {
  # Call-ability and return types/shapes only; the identity/swap equalities are the Fortran
  # suite's job.
  expression <- label_test_expression()
  assignment <- c(2L, 1L, 3:13)
  result <- run_js_comp_test_on_sample_assignment(expression, assignment, 10L, 30L)
  assert_true(is.list(result), "expected a list")
  assert_equal_int(length(result$global_js_divergence), 3L, "global_js_divergence length")
  assert_true(all(is.finite(result$global_js_divergence)), "global_js_divergence finite")
  assert_true(length(result$used_n_points) == 1 && length(result$used_n_neighbors) == 1, "used parameters scalars")

  repeated <- assignment
  repeated[2] <- repeated[1]
  assert_error(run_js_comp_test_on_sample_assignment(expression, repeated, 10L, 30L),
               "expected ERR_INVALID_INPUT for an assignment that is not a permutation", ERR_INVALID_INPUT)

  # divergence_mode reaches Fortran: the first two studies (9 samples), see
  # assert_two_study_divergence_modes.
  two_studies <- expression[, , 1:2]
  by_mode <- lapply(c(symmetric = "divergence_symmetric", pairwise = "divergence_pairwise_k2"), function(mode) {
    run_js_comp_test_on_sample_assignment(two_studies, 1:9, 10L, 30L, divergence_mode = mode)$global_js_divergence
  })
  assert_two_study_divergence_modes(by_mode$symmetric, by_mode$pairwise, "run_js_comp_test_on_sample_assignment")

  # The grid options reach the rerun search (they only matter with it). G = 300: the default grid
  # is the single candidate (200, 6) (ceiling(4*sqrt(300)) = 70 clamps up to the floor 200); with
  # grid_points_factor 17.32 and grid_min_points 240 it is (300, 4), (240, 5)
  # (ceiling(17.32*sqrt(300)) = 300, 0.8*300 = 240), so the rerun uses another n_points.
  default_rerun <- run_js_comp_test_on_sample_assignment(expression, 1:13, 10L, 30L, rerun_parameter_search = TRUE,
                                                         n_bootstraps = 20L)
  assert_true(default_rerun$used_n_points == 200, "default grid: expected 200 points")
  optioned_rerun <- run_js_comp_test_on_sample_assignment(expression, 1:13, 10L, 30L, rerun_parameter_search = TRUE,
                                                          n_bootstraps = 20L, grid_points_factor = 17.32,
                                                          grid_min_points = 240L)
  assert_true(optioned_rerun$used_n_points %in% c(300, 240), "grid options: expected 300 or 240 points")
  assert_error(run_js_comp_test_on_sample_assignment(expression, assignment, 10L, 30L, divergence_mode = "nonsense"),
               "expected an unknown divergence_mode string to be rejected before Fortran")
  assert_error(run_js_comp_test_on_sample_assignment(expression, assignment, 10L, 30L,
                                                     neighborhood_size_mode = "neighborhood_size_both"),
               "expected the removed neighborhood_size_mode 'both' to be rejected before Fortran")
}

test_run_js_comp_test_adaptive_on_sample_assignment <- function() {
  expression <- label_test_expression()
  result <- run_js_comp_test_adaptive_on_sample_assignment(expression, 1:13, 30L, 10L, 90L)
  assert_true(is.list(result), "expected a list")
  assert_equal_int(length(result$global_js_divergence), 3L, "global_js_divergence length")
  assert_true(all(is.finite(result$global_js_divergence)), "global_js_divergence finite")
  for (key in c("used_k_start", "used_k_step", "used_k_max", "n_points")) {
    assert_true(length(result[[key]]) == 1, paste(key, "must be a scalar"))
  }
  assert_error(run_js_comp_test_adaptive_on_sample_assignment(expression, 1:13, 2000L, 10L, 2000L),
               "expected ERR_INVALID_INPUT when no neighborhood is built", ERR_INVALID_INPUT)
  assert_error(run_js_comp_test_adaptive_on_sample_assignment(expression, 1:13, 30L, 10L, 90L,
                                                              divergence_mode = "nonsense"),
               "expected an unknown divergence_mode string to be rejected before Fortran")
  two_studies <- expression[, , 1:2]
  by_mode <- lapply(c(symmetric = "divergence_symmetric", pairwise = "divergence_pairwise_k2"), function(mode) {
    run_js_comp_test_adaptive_on_sample_assignment(two_studies, 1:9, 30L, 10L, 90L,
                                                   divergence_mode = mode)$global_js_divergence
  })
  assert_two_study_divergence_modes(by_mode$symmetric, by_mode$pairwise,
                                    "run_js_comp_test_adaptive_on_sample_assignment")
  # min_study_neighbors above every study's share of a neighborhood: construction status not ok.
  assert_error(run_js_comp_test_adaptive_on_sample_assignment(expression, 1:13, 30L, 10L, 90L,
                                                              min_study_neighbors = 1000L),
               "expected ERR_INVALID_INPUT when the construction status is not ok", ERR_INVALID_INPUT)
}

test_run_js_comp_test_permutation <- function() {
  expression <- label_test_expression()
  result <- run_js_comp_test_permutation(expression, 10L, 30L, 4L, random_seed = 3L)
  assert_true(is.list(result), "expected a list")
  assert_equal_int(length(result$global_js_divergence_observed), 3L, "observed length")
  assert_equal_int(length(result$p_values), 3L, "p_values length")
  assert_equal_int(dim(result$jsd_null), c(4L, 3L), "jsd_null shape")
  assert_true(length(result$n_failed_permutations) == 1, "n_failed_permutations scalar")
  assert_true(all(result$p_values > 0 & result$p_values <= 1), "p-values in (0, 1]")
  assert_error(run_js_comp_test_permutation(expression, 10L, 30L, 0L),
               "expected ERR_INVALID_INPUT for n_permutations = 0", ERR_INVALID_INPUT)
  # divergence_mode reaches Fortran: the observed values on the first two studies.
  two_studies <- expression[, , 1:2]
  by_mode <- lapply(c(symmetric = "divergence_symmetric", pairwise = "divergence_pairwise_k2"), function(mode) {
    run_js_comp_test_permutation(two_studies, 10L, 30L, 1L, random_seed = 3L,
                                 divergence_mode = mode)$global_js_divergence_observed
  })
  assert_two_study_divergence_modes(by_mode$symmetric, by_mode$pairwise, "run_js_comp_test_permutation")
  # The grid options reach the rerun search: the observed values equal the rerun on the identity
  # with the same options (whose n_points is checked in test_run_js_comp_test_on_sample_assignment)
  # and differ from the default grid's.
  optioned <- run_js_comp_test_permutation(expression, 10L, 30L, 1L, random_seed = 3L, rerun_parameter_search = TRUE,
                                           n_bootstraps = 20L, grid_points_factor = 17.32, grid_min_points = 240L)
  optioned_identity <- run_js_comp_test_on_sample_assignment(expression, 1:13, 10L, 30L,
                                                             rerun_parameter_search = TRUE, n_bootstraps = 20L,
                                                             grid_points_factor = 17.32, grid_min_points = 240L)
  default_identity <- run_js_comp_test_on_sample_assignment(expression, 1:13, 10L, 30L,
                                                            rerun_parameter_search = TRUE, n_bootstraps = 20L)
  assert_true(identical(optioned$global_js_divergence_observed, optioned_identity$global_js_divergence),
              "grid options: the observed values are the rerun's with the same options")
  assert_true(!identical(optioned$global_js_divergence_observed, default_identity$global_js_divergence),
              "grid options: another observed value than the default grid's")
  assert_error(run_js_comp_test_permutation(expression, 10L, 30L, 4L, divergence_mode = "nonsense"),
               "expected an unknown divergence_mode string to be rejected before Fortran")
}

test_run_js_comp_test_adaptive_permutation <- function() {
  expression <- label_test_expression()
  result <- run_js_comp_test_adaptive_permutation(expression, 30L, 10L, 90L, 3L, random_seed = 3L)
  assert_true(is.list(result), "expected a list")
  assert_equal_int(length(result$global_js_divergence_observed), 3L, "observed length")
  assert_equal_int(length(result$p_values), 3L, "p_values length")
  assert_equal_int(dim(result$jsd_null), c(3L, 3L), "jsd_null shape")
  assert_true(length(result$n_failed_permutations) == 1, "n_failed_permutations scalar")
  assert_error(run_js_comp_test_adaptive_permutation(expression, 30L, 10L, 90L, 0L),
               "expected ERR_INVALID_INPUT for n_permutations = 0", ERR_INVALID_INPUT)
  two_studies <- expression[, , 1:2]
  by_mode <- lapply(c(symmetric = "divergence_symmetric", pairwise = "divergence_pairwise_k2"), function(mode) {
    run_js_comp_test_adaptive_permutation(two_studies, 30L, 10L, 90L, 1L, random_seed = 3L,
                                          divergence_mode = mode)$global_js_divergence_observed
  })
  assert_two_study_divergence_modes(by_mode$symmetric, by_mode$pairwise, "run_js_comp_test_adaptive_permutation")
  assert_error(run_js_comp_test_adaptive_permutation(expression, 30L, 10L, 90L, 3L, divergence_mode = "nonsense"),
               "expected an unknown divergence_mode string to be rejected before Fortran")
}

run_all_tests()
