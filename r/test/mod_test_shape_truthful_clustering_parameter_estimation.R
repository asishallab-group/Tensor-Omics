source("r/load_tensor_omics.R")
source("r/test_helpers.R")

TOL <- 1e-9

# =====================
# sample_estimator_anchors
# =====================
# 11 points, density labels equal to point index (1..11, already ascending -- the sort
# permutation is the identity). n_anchors=5 gives percentiles 20/40/60/80/100, ranks
# 3/5/7/9/11 exactly (no interpolation rounding needed).
test_sample_anchors_hand_computed <- function() {
  density_labels <- as.double(1:11)
  anchors <- sample_estimator_anchors(density_labels, n_anchors = 5)
  expected <- c(3, 5, 7, 9, 11)
  assert_true(all(anchors == expected))
}

# 5 points, n_anchors=5 (every point its own percentile mark): ranks 2/3/3/4/5 -- rank 3 is
# hit twice, so anchor_indices necessarily repeats an index. Documented, not a bug.
test_sample_anchors_duplicates_possible <- function() {
  density_labels <- as.double(1:5)
  anchors <- sample_estimator_anchors(density_labels, n_anchors = 5)
  expected <- c(2, 3, 3, 4, 5)
  assert_true(all(anchors == expected))
}

test_sample_anchors_invalid_n_anchors_zero <- function() {
  density_labels <- as.double(1:11)
  assert_error(sample_estimator_anchors(density_labels, n_anchors = 0),
               "Expected error for n_anchors < 1", ERR_INVALID_INPUT)
}

test_sample_anchors_invalid_n_anchors_too_large <- function() {
  density_labels <- as.double(1:11)
  assert_error(sample_estimator_anchors(density_labels, n_anchors = 12),
               "Expected error for n_anchors > n_vectors", ERR_INVALID_INPUT)
}

# =====================
# grow_estimator_anchor_clouds
# =====================
# D=2, N=7 points on a line, (0,0)..(6,0). Two anchors at the opposite ends, point 1 (x=0,
# 1-indexed) and point 7 (x=6).
line_fixture_7 <- function() {
  rbind(seq(0, 6), rep(0, 7))
}

# seed_max_set_size=100 (grow until every point is claimed). Every round is an exact distance
# tie between the two clouds' own nearest-unclaimed candidate (both always 1.0 apart on this
# evenly-spaced line) -- ties are broken by whichever EA is scanned first (anchor 1), so
# anchor 1 wins every single round. Anchor 2's cloud never grows past its own single point.
test_grow_clouds_symmetric_line_ties_favor_lower_index <- function() {
  vectors <- line_fixture_7()
  anchor_indices <- c(1L, 7L)
  result <- grow_estimator_anchor_clouds(vectors, anchor_indices, seed_max_set_size = 100.0)

  expected_sizes <- c(6, 1)
  assert_true(all(result$cloud_sizes == expected_sizes))

  expected_cloud_1 <- rep(FALSE, 7)
  expected_cloud_1[1:6] <- TRUE
  assert_true(all(result$cloud_masks[, 1] == expected_cloud_1))
  assert_true(result$cloud_masks[7, 2] && sum(result$cloud_masks[, 2]) == 1)
}

# Same fixture, seed_max_set_size=50 -> ceiling(0.5*7)=4 total claims: 2 anchors already
# present plus 2 more rounds, both won by EA 1.
test_grow_clouds_seed_max_set_size_stops_early <- function() {
  vectors <- line_fixture_7()
  anchor_indices <- c(1L, 7L)
  result <- grow_estimator_anchor_clouds(vectors, anchor_indices, seed_max_set_size = 50.0)

  expected_sizes <- c(3, 1)
  assert_true(all(result$cloud_sizes == expected_sizes))
  assert_true(sum(result$cloud_masks) == 4)
}

# Default seed_max_set_size (5.0): ceiling(0.05*7)=1, clamped up to n_anchors=2 itself --
# actual_max_claims never exceeds the anchor count, so no growth happens at all.
test_grow_clouds_default_seed_max_set_size_can_yield_no_growth <- function() {
  vectors <- line_fixture_7()
  anchor_indices <- c(1L, 7L)
  result <- grow_estimator_anchor_clouds(vectors, anchor_indices)

  expected_sizes <- c(1, 1)
  assert_true(all(result$cloud_sizes == expected_sizes))
}

test_grow_clouds_invalid_seed_max_set_size <- function() {
  vectors <- line_fixture_7()
  anchor_indices <- c(1L, 7L)
  assert_error(grow_estimator_anchor_clouds(vectors, anchor_indices, seed_max_set_size = 150.0),
               "Expected error for seed_max_set_size > 100", ERR_INVALID_INPUT)
}

# =====================
# estimate_stc_parameters
# =====================
# D=2, N=21, a perfectly collinear, evenly-spaced line (0,0)..(20,0). Every estimator
# anchor's grown cloud is itself a sub-interval of the same line, so every EA agrees exactly
# on d=1 and on tangent direction -- chordal_dist_max_as_prcnt_of_range and d_max must both
# come out at (or, for the former, an SVD-residual hair above) 0. This fixture is exactly
# symmetric, so sample_estimator_anchors_impl's own tie-break (ties resolved by ascending point
# index) is what pins the anchor set to [3,5,7,9,11] and, through it, k_min/density_quantile
# below to a single deterministic outcome. k_min/k_density/density_quantile/G_max are
# cross-checked against this exact, already-verified, fully deterministic kernel's own real
# output (no randomness anywhere here).
collinear_line_21 <- function() {
  vectors <- rbind(seq(0, 20), rep(0, 21))
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  list(vectors = vectors, kd_indices = kd_indices, dimension_order = dimension_order)
}

test_estimate_parameters_collinear_line <- function() {
  fx <- collinear_line_21()
  result <- estimate_stc_parameters(fx$vectors, fx$kd_indices, fx$dimension_order, seed_max_set_size = 50.0)

  assert_true(result$estimated_chordal_dist_max_as_prcnt_of_range < 1e-6)
  assert_true(abs(result$estimated_d_max - 0.0) < TOL)
  assert_true(abs(result$estimated_k_min - 2.0) < TOL)
  assert_true(abs(result$estimated_k_density - result$estimated_k_min) < TOL)
  assert_true(abs(result$estimated_density_quantile - 1.0) < TOL)
  assert_true(abs(result$estimated_G_max - 0.0) < TOL)
}

# seed_max_set_size=0: every EA's cloud stays size 1 -- zero clouds ever reach the size >= 2
# a genuine observable/SVD needs. Fewer than 2 usable EAs is a genuine, data-dependent
# runtime failure, not a validation error.
test_estimate_parameters_too_few_valid_eas <- function() {
  fx <- collinear_line_21()
  assert_error(estimate_stc_parameters(fx$vectors, fx$kd_indices, fx$dimension_order, seed_max_set_size = 0.0),
               "Expected error when fewer than 2 EAs ever grow past size 1", ERR_INTERNAL)
}

test_estimate_parameters_invalid_n_anchors <- function() {
  fx <- collinear_line_21()
  assert_error(estimate_stc_parameters(fx$vectors, fx$kd_indices, fx$dimension_order, n_anchors = 50),
               "Expected error for n_anchors > n_vectors", ERR_INVALID_INPUT)
}

test_estimate_parameters_invalid_seed_max_set_size <- function() {
  fx <- collinear_line_21()
  assert_error(estimate_stc_parameters(fx$vectors, fx$kd_indices, fx$dimension_order, seed_max_set_size = -1.0),
               "Expected error for seed_max_set_size < 0", ERR_INVALID_INPUT)
}

# The Fortran suite's test_estimate_parameters_omitted_n_anchors_is_clamped is a regression
# test for a crash that could only happen with n_anchors truly absent at the Fortran ABI
# boundary -- not reproducible here, since the R binding always resolves and passes
# n_anchors=5 explicitly (see misc/code_gen_footgun.md's third entry). What is worth covering
# from R: that this always-explicit default of 5 still gets validated normally on a dataset
# smaller than it.
test_estimate_parameters_default_n_anchors_too_large_for_dataset <- function() {
  vectors <- rbind(seq(0, 2), rep(0, 3))
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  assert_error(estimate_stc_parameters(vectors, kd_indices, dimension_order),
               "Expected error for the default n_anchors=5 on a 3-point dataset", ERR_INVALID_INPUT)
}

run_all_tests()
