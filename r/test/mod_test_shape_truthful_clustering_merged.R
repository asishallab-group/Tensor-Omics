source("r/load_tensor_omics.R")
source("r/test_helpers.R")

STOP_REASON_FIXED_POINT <- 4L
STOP_REASON_ERROR <- 0L
MEMBER_ADDED_AT_STEP_NON_MEMBER <- -1L
MEMBER_ADDED_AT_STEP_SEED <- 0L

# D=2, N=7. Matches fixture_a in mod_test_shape_truthful_clustering.R exactly.
fixture_a <- function() {
  vectors <- cbind(rbind(seq(0, 4), rep(0, 5)), c(0, 1.5), c(0, 3.0))
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  list(vectors = vectors, kd_indices = kd_indices, dimension_order = dimension_order)
}

test_single_seed_matches_per_seed_kernel <- function() {
  fx <- fixture_a()
  seed_selection_mask <- rep(FALSE, 7)
  seed_selection_mask[1] <- TRUE

  res <- ensemble_identification_merged(fx$vectors, fx$kd_indices, fx$dimension_order, seed_selection_mask,
                                        k_min = 1, chordal_dist_max_as_prcnt_of_range = 0.1, d_max = 0, G_max = 1e10, RMSE_change_max = 1e10, o = 4)

  assert_true(res$ensemble_stop_reason[1] == STOP_REASON_FIXED_POINT)
  assert_true(abs(res$ensemble_growth_radii[1] - 1.0) < 1e-9)

  expected_mask <- rep(FALSE, 7)
  expected_mask[1:5] <- TRUE
  assert_true(all(res$ensemble_masks[, 1] == expected_mask))

  assert_true(all(res$ensemble_k_history[, 1] == c(2, 3, 4, 5)))
  assert_true(all(res$ensemble_accepted_history[, 1]))

  expected_step <- c(MEMBER_ADDED_AT_STEP_SEED, 1, 2, 3, 4, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                     MEMBER_ADDED_AT_STEP_NON_MEMBER)
  assert_true(all(res$ensemble_member_added_at_step[, 1] == expected_step))

  # Iteration 1's own bootstrap mask -- {seed=1, its one growth-radius neighbor=2}.
  expected_low_confidence <- rep(FALSE, 7)
  expected_low_confidence[1:2] <- TRUE
  assert_true(all(res$ensemble_low_confidence_masks[, 1] == expected_low_confidence))
}

test_two_independent_seeds <- function() {
  vectors <- cbind(rbind(seq(0, 4), rep(0, 5)), c(0, 1.5), c(0, 3.0))
  vectors <- cbind(vectors, vectors[, 1:7, drop = FALSE] + 100)
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)

  seed_selection_mask <- rep(FALSE, 14)
  seed_selection_mask[1] <- TRUE
  seed_selection_mask[8] <- TRUE

  res <- ensemble_identification_merged(vectors, kd_indices, dimension_order, seed_selection_mask,
                                        k_min = 1, chordal_dist_max_as_prcnt_of_range = 0.1, d_max = 0, G_max = 1e10, RMSE_change_max = 1e10, o = 4)

  assert_true(res$ensemble_stop_reason[1] == STOP_REASON_FIXED_POINT)
  assert_true(res$ensemble_stop_reason[2] == STOP_REASON_FIXED_POINT)

  expected_mask_1 <- rep(FALSE, 14)
  expected_mask_1[1:5] <- TRUE
  assert_true(all(res$ensemble_masks[, 1] == expected_mask_1))

  expected_mask_2 <- rep(FALSE, 14)
  expected_mask_2[8:12] <- TRUE
  assert_true(all(res$ensemble_masks[, 2] == expected_mask_2))

  assert_true(all(res$ensemble_k_history[, 1] == c(2, 3, 4, 5)))
  assert_true(all(res$ensemble_k_history[, 2] == c(2, 3, 4, 5)))

  expected_step_1 <- rep(MEMBER_ADDED_AT_STEP_NON_MEMBER, 14)
  expected_step_1[1] <- MEMBER_ADDED_AT_STEP_SEED
  expected_step_1[2:5] <- c(1, 2, 3, 4)
  assert_true(all(res$ensemble_member_added_at_step[, 1] == expected_step_1))

  expected_step_2 <- rep(MEMBER_ADDED_AT_STEP_NON_MEMBER, 14)
  expected_step_2[8] <- MEMBER_ADDED_AT_STEP_SEED
  expected_step_2[9:12] <- c(1, 2, 3, 4)
  assert_true(all(res$ensemble_member_added_at_step[, 2] == expected_step_2))

  # ensemble_U_first/ensemble_d_first: each column is its own seed's bootstrap basis, collinear
  # along the x-axis in both copies -- must not leak across columns.
  assert_true(res$ensemble_d_first[1] == 1)
  assert_true(res$ensemble_d_first[2] == 1)
  assert_true(abs(abs(res$ensemble_U_first[1, 1, 1]) - 1.0) < 1e-9)
  assert_true(abs(abs(res$ensemble_U_first[2, 1, 1]) - 0.0) < 1e-9)
  assert_true(abs(abs(res$ensemble_U_first[1, 1, 2]) - 1.0) < 1e-9)
  assert_true(abs(abs(res$ensemble_U_first[2, 1, 2]) - 0.0) < 1e-9)
}

test_zero_seeds <- function() {
  fx <- fixture_a()
  seed_selection_mask <- rep(FALSE, 7)

  res <- ensemble_identification_merged(fx$vectors, fx$kd_indices, fx$dimension_order, seed_selection_mask,
                                        k_min = 1, chordal_dist_max_as_prcnt_of_range = 0.1, d_max = 0, G_max = 1e10, RMSE_change_max = 1e10, o = 4)
  assert_true(all(dim(res$ensemble_masks) == c(7, 0)))
}

test_n_dimensions_too_small <- function() {
  vectors <- rbind(seq(0, 6))
  dimension_order <- c(1L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  seed_selection_mask <- rep(FALSE, 7)
  seed_selection_mask[1] <- TRUE

  assert_error(ensemble_identification_merged(vectors, kd_indices, dimension_order, seed_selection_mask,
                                              k_min = 1, chordal_dist_max_as_prcnt_of_range = 0.1, d_max = 0, G_max = 1e10, RMSE_change_max = 1e10, o = 4),
               "Expected error for n_dimensions=1", ERR_INVALID_INPUT)
}

run_all_tests()
