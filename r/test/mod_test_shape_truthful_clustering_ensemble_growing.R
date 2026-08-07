source("r/load_tensor_omics.R")
source("r/test_helpers.R")

TOL <- 1e-9

line_fixture <- function(n = 11) {
  vectors <- rbind(seq(0, n - 1), rep(0, n))
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  list(vectors = vectors, kd_indices = kd_indices, dimension_order = dimension_order)
}

# =====================
# calc_ensemble_growth_radius
# =====================
# Seed at x=5 (1-based index 6). Its 4 nearest neighbors are x=4,6 (distance 1) and x=3,7
# (distance 2), sorted [1,1,2,2] -- even k_min=4, median = avg(2nd,3rd) = 1.5. Its 3 nearest
# are x=4,6 (distance 1) and one of x=3/x=7 (distance 2, a tie), sorted [1,1,2] -- odd k_min=3,
# median = the middle element = 1.
test_growth_radius_even_k <- function() {
  fx <- line_fixture(11)
  radius <- calc_ensemble_growth_radius(fx$vectors, fx$kd_indices, fx$dimension_order, 6, k_min = 4)
  assert_true(abs(radius - 1.5) < TOL)
}

test_growth_radius_odd_k <- function() {
  fx <- line_fixture(11)
  radius <- calc_ensemble_growth_radius(fx$vectors, fx$kd_indices, fx$dimension_order, 6, k_min = 3)
  assert_true(abs(radius - 1.0) < TOL)
}

test_growth_radius_seed_index_out_of_range <- function() {
  fx <- line_fixture(11)
  assert_error(calc_ensemble_growth_radius(fx$vectors, fx$kd_indices, fx$dimension_order, 12),
               "Expected error for seed_index > n_vectors", ERR_INVALID_INPUT)
}

test_growth_radius_k_min_too_large <- function() {
  fx <- line_fixture(11)
  assert_error(calc_ensemble_growth_radius(fx$vectors, fx$kd_indices, fx$dimension_order, 6, k_min = 11),
               "Expected error for k_min > n_vectors-1", ERR_INVALID_INPUT)
}

# =====================
# grow_ensemble
# =====================
# Same 11-point line fixture. With growth_radius=1.5: from a single member x=5 (index 6):
# covers x=4,5,6 (indices 5,6,7). From members x=4,5,6 (indices 5,6,7): the union covers
# x=3..7 (indices 4..8).
test_grow_ensemble_single_member <- function() {
  fx <- line_fixture(11)
  is_member_mask <- rep(FALSE, 11)
  is_member_mask[6] <- TRUE

  result <- grow_ensemble(fx$vectors, fx$kd_indices, fx$dimension_order, is_member_mask, 1.5)
  expected <- rep(FALSE, 11)
  expected[5:7] <- TRUE
  assert_true(all(result == expected))
}

test_grow_ensemble_multi_member_union <- function() {
  fx <- line_fixture(11)
  is_member_mask <- rep(FALSE, 11)
  is_member_mask[5:7] <- TRUE

  result <- grow_ensemble(fx$vectors, fx$kd_indices, fx$dimension_order, is_member_mask, 1.5)
  expected <- rep(FALSE, 11)
  expected[4:8] <- TRUE
  assert_true(all(result == expected))
}

# An all-FALSE is_member_mask is a well-defined degenerate case (nothing to grow from), not a
# validation error -- see the Fortran test's comment for why.
test_grow_ensemble_empty_ensemble_is_degenerate <- function() {
  fx <- line_fixture(11)
  is_member_mask <- rep(FALSE, 11)

  result <- grow_ensemble(fx$vectors, fx$kd_indices, fx$dimension_order, is_member_mask, 1.5)
  assert_true(!any(result))
}

test_grow_ensemble_negative_radius <- function() {
  fx <- line_fixture(11)
  is_member_mask <- rep(FALSE, 11)
  is_member_mask[6] <- TRUE

  assert_error(grow_ensemble(fx$vectors, fx$kd_indices, fx$dimension_order, is_member_mask, -1.5),
               "Expected error for a negative growth radius", ERR_INVALID_INPUT)
}

run_all_tests()
