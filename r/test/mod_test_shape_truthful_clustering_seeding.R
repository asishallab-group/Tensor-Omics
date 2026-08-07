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
# density_labels
# =====================
# 3 points on a line, (0,0),(1,0),(3,0), k_density=2 (every other point). Bandwidth is the
# bandwidth_percentile-th percentile (default 68.27, the heuristic "1 SD" anchor) of each
# point's own k-NN distances, via calc_percentile_helper. Hand-computed (and cross-checked
# against an independent Python re-implementation of the same formula) expected densities --
# see the Fortran test's own comment for the full derivation.
test_density_labels_hand_computed <- function() {
  vectors <- rbind(c(0, 1, 3), c(0, 0, 0))
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)

  labels <- density_labels(vectors, kd_indices, dimension_order, k_density = 2)
  expected <- c(0.2434133437, 0.4702740525, 0.1795903795)
  assert_true(all(abs(labels - expected) < 1e-6))
}

# Same fixture, but bandwidth_percentile=50.0 (the median) instead of the default 68.27 --
# confirms the parameter actually changes the bandwidth and thus the resulting labels.
test_density_labels_bandwidth_percentile_median <- function() {
  vectors <- rbind(c(0, 1, 3), c(0, 0, 0))
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)

  labels <- density_labels(vectors, kd_indices, dimension_order, k_density = 2, bandwidth_percentile = 50.0)
  expected <- c(0.3017873425, 0.5385998637, 0.1940642069)
  assert_true(all(abs(labels - expected) < 1e-6))
}

test_density_labels_invalid_bandwidth_percentile <- function() {
  fx <- line_fixture(11)
  assert_error(density_labels(fx$vectors, fx$kd_indices, fx$dimension_order, k_density = 4, bandwidth_percentile = 101.0),
               "Expected error for bandwidth_percentile > 100", ERR_INVALID_INPUT)
}

# The center of an evenly-spaced plus shape has all 4 of its k_density=4 neighbors at the
# identical distance 0.1 -- the percentile-based bandwidth is exactly 0.1 (every distance
# equal), never zero, so this just confirms the label is a genuine, representable, strictly
# positive number.
test_density_labels_symmetric_neighborhood_does_not_underflow <- function() {
  vectors <- rbind(
    c(0.0, 0.1, 0.0, -0.1, 0.0),
    c(0.0, 0.0, 0.1, 0.0, -0.1)
  )
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)

  labels <- density_labels(vectors, kd_indices, dimension_order, k_density = 4)
  assert_true(labels[1] > 0.0)
}

# On an evenly-spaced 11-point line, every interior point's k_density=4 nearest neighbors
# form the identical distance pattern [1,1,2,2] by translation symmetry, so all interior
# points (3..9) must get exactly the same density label.
test_density_labels_uniform_interior_points_agree <- function() {
  fx <- line_fixture(11)
  labels <- density_labels(fx$vectors, fx$kd_indices, fx$dimension_order, k_density = 4)
  for (i in 4:9) {
    assert_true(abs(labels[i] - labels[3]) < 1e-9)
  }
}

# A dense cluster (spacing 0.1) and a sparse cluster (spacing 2.0), far enough apart that
# k_density=2 never crosses between them: the dense cluster's adaptive bandwidth is far
# smaller, so its members must get a strictly higher density label.
test_density_labels_dense_vs_sparse <- function() {
  vectors <- rbind(
    c(0.0, 0.1, 0.2, 100.0, 102.0, 104.0),
    c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  )
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)

  labels <- density_labels(vectors, kd_indices, dimension_order, k_density = 2)
  assert_true(labels[2] > labels[5])
}

test_density_labels_invalid_kd_indices <- function() {
  fx <- line_fixture(11)
  bad_kd_indices <- fx$kd_indices
  bad_kd_indices[1] <- 12L
  assert_error(density_labels(fx$vectors, bad_kd_indices, fx$dimension_order, k_density = 4),
               "Expected error for kd_indices entry > n_vectors", ERR_INVALID_INPUT)
}

test_density_labels_k_density_too_large <- function() {
  fx <- line_fixture(11)
  assert_error(density_labels(fx$vectors, fx$kd_indices, fx$dimension_order, k_density = 11),
               "Expected error for k_density > n_vectors - 1", ERR_INVALID_INPUT)
}

# =====================
# seeds
# =====================
# Two separated 2-point clusters, k_density=1: each point's single nearest neighbor is always
# its own cluster-mate (0.1 apart), never the other cluster (10 apart) -- see the Fortran
# test's own comment for why this stays a 2-point, k_density=1 fixture rather than a larger
# symmetric one.
two_clusters_fixture <- function() {
  vectors <- rbind(
    c(0.0, 0.1, 10.0, 10.1),
    c(0.0, 0.0, 0.0, 0.0)
  )
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  list(vectors = vectors, kd_indices = kd_indices, dimension_order = dimension_order)
}

test_seeds_two_separated_clusters <- function() {
  fx <- two_clusters_fixture()
  is_seed_mask <- seeds(fx$vectors, fx$kd_indices, fx$dimension_order, k_density = 1)
  assert_true(sum(is_seed_mask) == 2)
  assert_true(any(is_seed_mask[1:2]))
  assert_true(any(is_seed_mask[3:4]))
}

test_seeds_single_cluster_one_seed <- function() {
  vectors <- rbind(c(0.0, 0.1), c(0.0, 0.0))
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)

  is_seed_mask <- seeds(vectors, kd_indices, dimension_order, k_density = 1)
  assert_true(sum(is_seed_mask) == 1)
}

test_seeds_invalid_k_density <- function() {
  vectors <- rbind(
    c(0.0, 0.1, 0.0, -0.1, 0.0),
    c(0.0, 0.0, 0.1, 0.0, -0.1)
  )
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  assert_error(seeds(vectors, kd_indices, dimension_order, k_density = 0),
               "Expected error for k_density < 1", ERR_INVALID_INPUT)
}

# The Fortran suite's test_density_labels_omitted_k_density_is_clamped/
# test_seeds_omitted_k_density_is_clamped are regression tests for a crash that could only
# happen with k_density truly absent at the Fortran ABI boundary -- not reproducible here,
# since the R binding always resolves and passes k_density=30 explicitly, never a
# genuinely-absent optional (see misc/code_gen_footgun.md's third entry). What is worth
# covering from R: that this always-explicit default of 30 still gets validated normally (a
# clean, typed error, not a crash) on a dataset smaller than it.
test_seeds_default_k_density_too_large_for_dataset <- function() {
  vectors <- rbind(
    c(0.0, 0.1, 0.0, -0.1, 0.0),
    c(0.0, 0.0, 0.1, 0.0, -0.1)
  )
  dimension_order <- c(1L, 2L)
  kd_indices <- build_kd_index(vectors, dimension_order)
  assert_error(seeds(vectors, kd_indices, dimension_order),
               "Expected error for the default k_density=30 on a 5-point dataset", ERR_INVALID_INPUT)
}

run_all_tests()
