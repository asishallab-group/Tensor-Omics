source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# tox_shape_truthful_clustering_kernel's STOP_REASON_*/MEMBER_ADDED_AT_STEP_* parameters are
# not exposed as importable constants (no existing codegen mechanism does that for a plain
# result-code output, unlike a "mode" input) -- these mirror their literal Fortran values.
STOP_REASON_MAX_SIZE <- 1L
STOP_REASON_REJECTED_AFTER_STABLE <- 2L
STOP_REASON_REJECTED_IMMEDIATELY <- 3L
STOP_REASON_FIXED_POINT <- 4L
STOP_REASON_ERROR <- 0L
MEMBER_ADDED_AT_STEP_NON_MEMBER <- -1L
MEMBER_ADDED_AT_STEP_SEED <- 0L

# D=2, N=7. A 5-point line (0,0)..(4,0), plus two far-away points a growth radius of 1.0
# (k_min=1 from the seed's nearest neighbor at distance 1.0) never reaches.
fixture_a <- function() {
  vectors <- cbind(rbind(seq(0, 4), rep(0, 5)), c(0, 1.5), c(0, 3.0))
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  list(vectors = vectors, kd_indices = kd_indices, dimension_order = dimension_order)
}

# D=3, N=7. A 5-point x-axis line plus a branch at (1,1,0),(1,2,0) next to the 2nd x-axis
# point -- the branch is swept in already at growth iteration t=2.
fixture_b <- function() {
  vectors <- cbind(rbind(seq(0, 4), rep(0, 5), rep(0, 5)), c(1, 1, 0), c(1, 2, 0))
  dimension_order <- c(1L, 2L, 3L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  list(vectors = vectors, kd_indices = kd_indices, dimension_order = dimension_order)
}

# D=3, N=7. Same idea as fixture_b, but the branch sits next to the 3rd x-axis point, so it
# is only swept in at growth iteration t=3, after 2 accepted iterations.
fixture_c <- function() {
  vectors <- cbind(rbind(seq(0, 4), rep(0, 5), rep(0, 5)), c(2, 1, 0), c(2, 2, 0))
  dimension_order <- c(1L, 2L, 3L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  list(vectors = vectors, kd_indices = kd_indices, dimension_order = dimension_order)
}

test_natural_fixed_point <- function() {
  fx <- fixture_a()
  res <- ensemble_identification(fx$vectors, fx$kd_indices, fx$dimension_order, 1, k_min = 1,
                                 alpha_max = 0.1, d_max = 0, G_max = 1e10, o = 4)

  assert_true(res$stop_reason == STOP_REASON_FIXED_POINT)
  assert_true(abs(res$growth_radius - 1.0) < 1e-9)

  expected_mask <- rep(FALSE, 7)
  expected_mask[1:5] <- TRUE
  assert_true(all(res$final_ensemble_mask == expected_mask))

  assert_true(all(res$k_history == c(2, 3, 4, 5)))
  assert_true(all(res$accepted_history))
  assert_true(all(res$d_history == 1))

  expected_step <- c(MEMBER_ADDED_AT_STEP_SEED, 1, 2, 3, 4, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                     MEMBER_ADDED_AT_STEP_NON_MEMBER)
  assert_true(all(res$member_added_at_step == expected_step))

  # Iteration 1's own bootstrap mask -- {seed=1, its one growth-radius neighbor=2}.
  expected_low_confidence <- rep(FALSE, 7)
  expected_low_confidence[1:2] <- TRUE
  assert_true(all(res$low_confidence_mask == expected_low_confidence))
}

test_history_window_shifts <- function() {
  fx <- fixture_a()
  res <- ensemble_identification(fx$vectors, fx$kd_indices, fx$dimension_order, 1, k_min = 1,
                                 alpha_max = 0.1, d_max = 0, G_max = 1e10, o = 2)

  assert_true(res$stop_reason == STOP_REASON_FIXED_POINT)
  assert_true(all(res$k_history == c(4, 5)))
  assert_true(all(res$accepted_history))
}

test_max_size_at_bootstrap <- function() {
  fx <- fixture_a()
  res <- ensemble_identification(fx$vectors, fx$kd_indices, fx$dimension_order, 1, k_min = 1,
                                 alpha_max = 0.1, d_max = 0, G_max = 1e10, f_max = 0.2, o = 3)

  assert_true(res$stop_reason == STOP_REASON_MAX_SIZE)
  assert_true(!any(res$final_ensemble_mask))
  assert_true(all(res$k_history == 0))
  assert_true(!any(res$accepted_history))
  assert_true(all(res$member_added_at_step == MEMBER_ADDED_AT_STEP_NON_MEMBER))
  # Stop Condition 1 fires before observable is ever called at all: no genuine SVD ever
  # happened, so there is no real iteration-1 data to report.
  assert_true(!any(res$low_confidence_mask))
}

test_max_size_poisons_prior_accepts <- function() {
  fx <- fixture_a()
  res <- ensemble_identification(fx$vectors, fx$kd_indices, fx$dimension_order, 1, k_min = 1,
                                 alpha_max = 0.1, d_max = 0, G_max = 1e10, f_max = 0.35, o = 3)

  assert_true(res$stop_reason == STOP_REASON_MAX_SIZE)
  assert_true(!any(res$final_ensemble_mask))
  assert_true(all(res$k_history == 0))
  assert_true(all(res$member_added_at_step == MEMBER_ADDED_AT_STEP_NON_MEMBER))
  # The key behavior this output exists for: Stop Condition 1 wipes final_ensemble_mask and
  # every history array, but iteration 1's own genuinely-SVD-backed bootstrap mask --
  # {seed=1, its one growth-radius neighbor=2} -- survives that reset.
  expected_low_confidence <- rep(FALSE, 7)
  expected_low_confidence[1:2] <- TRUE
  assert_true(all(res$low_confidence_mask == expected_low_confidence))
}

test_rejected_immediately <- function() {
  fx <- fixture_b()
  res <- ensemble_identification(fx$vectors, fx$kd_indices, fx$dimension_order, 1, k_min = 1,
                                 alpha_max = 0.1, d_max = 0, G_max = 1e10, o = 2)

  assert_true(res$stop_reason == STOP_REASON_REJECTED_IMMEDIATELY)

  expected_mask <- rep(FALSE, 7)
  expected_mask[1:2] <- TRUE
  assert_true(all(res$final_ensemble_mask == expected_mask))

  assert_true(all(res$k_history == c(2, 4)))
  assert_true(all(res$accepted_history == c(TRUE, FALSE)))

  expected_step <- c(MEMBER_ADDED_AT_STEP_SEED, 1, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                     MEMBER_ADDED_AT_STEP_NON_MEMBER, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                     MEMBER_ADDED_AT_STEP_NON_MEMBER, MEMBER_ADDED_AT_STEP_NON_MEMBER)
  assert_true(all(res$member_added_at_step == expected_step))
}

test_rejected_after_stable <- function() {
  fx <- fixture_c()
  res <- ensemble_identification(fx$vectors, fx$kd_indices, fx$dimension_order, 1, k_min = 1,
                                 alpha_max = 0.1, d_max = 0, G_max = 1e10, o = 3)

  assert_true(res$stop_reason == STOP_REASON_REJECTED_AFTER_STABLE)

  expected_mask <- rep(FALSE, 7)
  expected_mask[1:3] <- TRUE
  assert_true(all(res$final_ensemble_mask == expected_mask))

  assert_true(all(res$k_history == c(2, 3, 5)))
  assert_true(all(res$accepted_history == c(TRUE, TRUE, FALSE)))

  expected_step <- c(MEMBER_ADDED_AT_STEP_SEED, 1, 2, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                     MEMBER_ADDED_AT_STEP_NON_MEMBER, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                     MEMBER_ADDED_AT_STEP_NON_MEMBER)
  assert_true(all(res$member_added_at_step == expected_step))
}

test_seed_index_out_of_range <- function() {
  fx <- fixture_a()
  assert_error(ensemble_identification(fx$vectors, fx$kd_indices, fx$dimension_order, 8, k_min = 1,
                                       alpha_max = 0.1, d_max = 0, G_max = 1e10, o = 3),
               "Expected error for seed_index > n_vectors", ERR_INVALID_INPUT)
}

test_o_zero <- function() {
  fx <- fixture_a()
  assert_error(ensemble_identification(fx$vectors, fx$kd_indices, fx$dimension_order, 1, k_min = 1,
                                       alpha_max = 0.1, d_max = 0, G_max = 1e10, o = 0),
               "Expected error for o=0", ERR_INVALID_INPUT)
}

test_n_dimensions_too_small <- function() {
  vectors <- rbind(seq(0, 6))
  dimension_order <- c(1L)
  kd_indices <- build_kd_index(vectors, dimension_order)

  assert_error(ensemble_identification(vectors, kd_indices, dimension_order, 1, k_min = 1,
                                       alpha_max = 0.1, d_max = 0, G_max = 1e10, o = 3),
               "Expected error for n_dimensions=1", ERR_INVALID_INPUT)
}

run_all_tests()
