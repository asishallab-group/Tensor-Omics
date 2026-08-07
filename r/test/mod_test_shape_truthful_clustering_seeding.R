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
# calculate_density_radius
# =====================
# Shared fixture: D=2, N=5 points on a line, (0,0),(1,0),(2,0),(3,0),(4,0). Mean = (2,0),
# so the mean-to-vector distances are [2,1,0,1,2], sorted ascending [0,1,1,2,2].
# Default 15th percentile: rank = 0.15*4+1 = 1.6 -> interpolate 0 and 1 at fraction 0.6 -> 0.6.
# 50th percentile: rank = 3.0 exactly -> 1.0.
radius_fixture <- function() {
  rbind(c(0, 1, 2, 3, 4), c(0, 0, 0, 0, 0))
}

test_density_radius_default_percentile <- function() {
  radius <- calculate_density_radius(radius_fixture())
  assert_true(abs(radius - 0.6) < TOL)
}

test_density_radius_custom_percentile <- function() {
  radius <- calculate_density_radius(radius_fixture(), mean_to_other_vecs_dist_quant = 0.5)
  assert_true(abs(radius - 1.0) < TOL)
}

test_density_radius_invalid_percentile <- function() {
  assert_error(calculate_density_radius(radius_fixture(), mean_to_other_vecs_dist_quant = 1.5),
               "Expected error for quantile > 1.0", ERR_INVALID_INPUT)
}

test_density_radius_single_vector <- function() {
  vectors <- rbind(5, 5)
  radius <- calculate_density_radius(vectors)
  assert_true(abs(radius - 0.0) < 1e-12)
}

# =====================
# density_labels
# =====================
test_density_labels_basic <- function() {
  fx <- line_fixture(11)
  labels <- density_labels(fx$vectors, fx$kd_indices, fx$dimension_order, 1.5)
  expected <- rep(3.0, 11)
  expected[1] <- 2.0
  expected[11] <- 2.0
  assert_true(all(abs(labels - expected) < 1e-12))
}

test_density_labels_zero_radius <- function() {
  fx <- line_fixture(11)
  labels <- density_labels(fx$vectors, fx$kd_indices, fx$dimension_order, 0.0)
  assert_true(all(abs(labels - 1.0) < 1e-12))
}

test_density_labels_invalid_kd_indices <- function() {
  fx <- line_fixture(11)
  bad_kd_indices <- fx$kd_indices
  bad_kd_indices[1] <- 12L
  assert_error(density_labels(fx$vectors, bad_kd_indices, fx$dimension_order, 1.5),
               "Expected error for kd_indices entry > n_vectors", ERR_INVALID_INPUT)
}

# =====================
# seeds
# =====================
two_clusters_fixture <- function() {
  vectors <- rbind(
    c(0.0, 0.1, 0.0, -0.1, 0.0, 10.0, 10.1, 10.0, 9.9, 10.0),
    c(0.0, 0.0, 0.1, 0.0, -0.1, 0.0, 0.0, 0.1, 0.0, -0.1)
  )
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  list(vectors = vectors, kd_indices = kd_indices, dimension_order = dimension_order)
}

test_seeds_two_separated_clusters <- function() {
  fx <- two_clusters_fixture()
  is_seed_mask <- seeds(fx$vectors, fx$kd_indices, fx$dimension_order)
  assert_true(sum(is_seed_mask) == 2)
  assert_true(any(is_seed_mask[1:5]))
  assert_true(any(is_seed_mask[6:10]))
}

single_cluster_fixture <- function() {
  vectors <- rbind(
    c(0.0, 0.1, 0.0, -0.1, 0.0),
    c(0.0, 0.0, 0.1, 0.0, -0.1)
  )
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  list(vectors = vectors, kd_indices = kd_indices, dimension_order = dimension_order)
}

# At the 100th percentile the density radius equals the cluster's own diameter (the
# farthest mean-to-vector distance), large enough to cover the whole cluster from a single
# pick. At the default (15th) percentile it would not -- see the Fortran test's comment for
# why that is a real property of the algorithm, not a bug.
test_seeds_single_cluster_one_seed <- function() {
  fx <- single_cluster_fixture()
  is_seed_mask <- seeds(fx$vectors, fx$kd_indices, fx$dimension_order, mean_to_other_vecs_dist_quant = 1.0)
  assert_true(sum(is_seed_mask) == 1)
}

test_seeds_invalid_percentile <- function() {
  fx <- single_cluster_fixture()
  assert_error(seeds(fx$vectors, fx$kd_indices, fx$dimension_order, mean_to_other_vecs_dist_quant = -0.1),
               "Expected error for negative mean_to_other_vecs_dist_quant", ERR_INVALID_INPUT)
}

run_all_tests()
