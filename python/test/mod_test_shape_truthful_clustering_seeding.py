"""
Python test suite for density_labels / seeds (tox_shape_truthful_clustering_seeding),
mirroring test/mod_test_shape_truthful_clustering_seeding.F90
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import density_labels, seeds, build_kd_index
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_INVALID_INPUT


def _line_fixture(n=11):
    """D=2, N points on a line: (0,0),(1,0),...,(n-1,0)."""
    vectors = np.zeros((2, n), dtype=np.float64, order='F')
    vectors[0, :] = np.arange(n, dtype=np.float64)
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    return vectors, kd_indices, dimension_order


# =====================
# density_labels
# =====================
# 3 points on a line, (0,0),(1,0),(3,0), k_density=2 (every other point). Bandwidth is the
# bandwidth_percentile-th percentile (default 68.27, the heuristic "1 SD" anchor) of each
# point's own k-NN distances, via calc_percentile_helper. Hand-computed (and cross-checked
# against an independent Python re-implementation of the same formula) expected densities --
# see the Fortran test's own comment for the full derivation.
def test_density_labels_hand_computed():
    vectors = np.array([[0.0, 1.0, 3.0], [0.0, 0.0, 0.0]], order='F')
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)

    labels = density_labels(vectors, kd_indices, dimension_order, k_density=2)
    expected = np.array([0.2434133437, 0.4702740525, 0.1795903795])
    assert np.allclose(labels, expected, atol=1e-6), f"expected {expected}, got {labels}"


# Same fixture, but bandwidth_percentile=50.0 (the median) instead of the default 68.27 --
# confirms the parameter actually changes the bandwidth and thus the resulting labels.
def test_density_labels_bandwidth_percentile_median():
    vectors = np.array([[0.0, 1.0, 3.0], [0.0, 0.0, 0.0]], order='F')
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)

    labels = density_labels(vectors, kd_indices, dimension_order, k_density=2, bandwidth_percentile=50.0)
    expected = np.array([0.3017873425, 0.5385998637, 0.1940642069])
    assert np.allclose(labels, expected, atol=1e-6), f"expected {expected}, got {labels}"


def test_density_labels_invalid_bandwidth_percentile():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    assert_error(lambda: density_labels(vectors, kd_indices, dimension_order, k_density=4, bandwidth_percentile=101.0),
                 "Expected error for bandwidth_percentile > 100", ERR_INVALID_INPUT)


def test_density_labels_symmetric_neighborhood_does_not_underflow():
    # The center of an evenly-spaced plus shape has all 4 of its k_density=4 neighbors at the
    # identical distance 0.1 -- the percentile-based bandwidth is exactly 0.1 (every distance
    # equal), never zero, so this just confirms the label is a genuine, representable,
    # strictly positive number.
    vectors = np.array([
        [0.0, 0.1, 0.0, -0.1, 0.0],
        [0.0, 0.0, 0.1, 0.0, -0.1],
    ], order='F')
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)

    labels = density_labels(vectors, kd_indices, dimension_order, k_density=4)
    assert labels[0] > 0.0, f"expected a strictly positive label, got {labels[0]}"


def test_density_labels_uniform_interior_points_agree():
    # On an evenly-spaced 11-point line, every interior point's k_density=4 nearest neighbors
    # form the identical distance pattern [1,1,2,2] by translation symmetry, so all interior
    # points (index 2..8, 0-based) must get exactly the same density label.
    vectors, kd_indices, dimension_order = _line_fixture(11)
    labels = density_labels(vectors, kd_indices, dimension_order, k_density=4)
    for i in range(3, 9):
        assert abs(labels[i] - labels[2]) < 1e-9, f"expected interior points to agree, got {labels}"


def test_density_labels_dense_vs_sparse():
    # A dense cluster (spacing 0.1) and a sparse cluster (spacing 2.0), far enough apart that
    # k_density=2 never crosses between them: the dense cluster's adaptive bandwidth is far
    # smaller, so its members must get a strictly higher density label.
    vectors = np.array([
        [0.0, 0.1, 0.2, 100.0, 102.0, 104.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ], order='F')
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)

    labels = density_labels(vectors, kd_indices, dimension_order, k_density=2)
    assert labels[1] > labels[4], f"expected the dense cluster's label to exceed the sparse cluster's, got {labels}"


def test_density_labels_invalid_kd_indices():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    bad_kd_indices = kd_indices.copy()
    bad_kd_indices[0] = 12
    assert_error(lambda: density_labels(vectors, bad_kd_indices, dimension_order, k_density=4),
                 "Expected error for kd_indices entry > n_vectors", ERR_INVALID_INPUT)


def test_density_labels_k_density_too_large():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    assert_error(lambda: density_labels(vectors, kd_indices, dimension_order, k_density=11),
                 "Expected error for k_density > n_vectors - 1", ERR_INVALID_INPUT)


# =====================
# seeds
# =====================
# Two separated 2-point clusters, k_density=1: each point's single nearest neighbor is always
# its own cluster-mate (0.1 apart), never the other cluster (10 apart) -- see the Fortran
# test's own comment for why this stays a 2-point, k_density=1 fixture rather than a larger
# symmetric one.
def _two_clusters_fixture():
    vectors = np.array([
        [0.0, 0.1, 10.0, 10.1],
        [0.0, 0.0, 0.0, 0.0],
    ], order='F')
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    return vectors, kd_indices, dimension_order


def test_seeds_two_separated_clusters():
    vectors, kd_indices, dimension_order = _two_clusters_fixture()
    is_seed_mask = seeds(vectors, kd_indices, dimension_order, k_density=1)
    assert int(np.sum(is_seed_mask)) == 2, f"expected 2 seeds, got {int(np.sum(is_seed_mask))}"
    assert np.any(is_seed_mask[0:2]), "cluster A (indices 0-1) should have a seed"
    assert np.any(is_seed_mask[2:4]), "cluster B (indices 2-3) should have a seed"


def test_seeds_single_cluster_one_seed():
    vectors = np.array([[0.0, 0.1], [0.0, 0.0]], order='F')
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)

    is_seed_mask = seeds(vectors, kd_indices, dimension_order, k_density=1)
    assert int(np.sum(is_seed_mask)) == 1, f"expected 1 seed, got {int(np.sum(is_seed_mask))}"


# The shared 11-point line fixture, k_density=4. A wider exclusion radius (100th percentile
# of the k_density distances, i.e. the farthest neighbor -- 2.0, vs. the default
# 50th-percentile median of 1.5) suppresses more of the line per seed, so fewer seeds are
# needed to cover it.
def test_seeds_exclusion_radius_percentile_widens_coverage():
    vectors, kd_indices, dimension_order = _line_fixture(11)

    mask_default = seeds(vectors, kd_indices, dimension_order, k_density=4)
    mask_wide = seeds(vectors, kd_indices, dimension_order, k_density=4, exclusion_radius_percentile=100.0)

    expected_default = np.zeros(11, dtype=np.bool_)
    expected_default[[1, 3, 5, 7, 9]] = True
    expected_wide = np.zeros(11, dtype=np.bool_)
    expected_wide[[0, 3, 7, 10]] = True

    assert np.array_equal(mask_default, expected_default), f"expected {expected_default}, got {mask_default}"
    assert np.array_equal(mask_wide, expected_wide), f"expected {expected_wide}, got {mask_wide}"
    assert int(np.sum(mask_wide)) < int(np.sum(mask_default)), \
        "a wider exclusion radius should need fewer seeds to cover the same line"


def test_seeds_invalid_exclusion_radius_percentile():
    vectors = np.array([
        [0.0, 0.1, 0.0, -0.1, 0.0],
        [0.0, 0.0, 0.1, 0.0, -0.1],
    ], order='F')
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    assert_error(lambda: seeds(vectors, kd_indices, dimension_order, k_density=4, exclusion_radius_percentile=101.0),
                 "Expected error for exclusion_radius_percentile > 100", ERR_INVALID_INPUT)


def test_seeds_invalid_k_density():
    vectors = np.array([
        [0.0, 0.1, 0.0, -0.1, 0.0],
        [0.0, 0.0, 0.1, 0.0, -0.1],
    ], order='F')
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    assert_error(lambda: seeds(vectors, kd_indices, dimension_order, k_density=0),
                 "Expected error for k_density < 1", ERR_INVALID_INPUT)


# The Fortran suite's test_density_labels_omitted_k_density_is_clamped/
# test_seeds_omitted_k_density_is_clamped are regression tests for a crash that could only
# happen with k_density truly *absent* at the Fortran ABI boundary -- not reproducible here,
# since the Python binding always resolves and passes k_density=30 explicitly, never a
# genuinely-absent optional (see misc/code_gen_footgun.md's third entry). What *is* worth
# covering from Python: that this always-explicit default of 30 still gets validated
# normally (a clean, typed error, not a crash) on a dataset smaller than it.
def test_seeds_default_k_density_too_large_for_dataset():
    vectors = np.array([
        [0.0, 0.1, 0.0, -0.1, 0.0],
        [0.0, 0.0, 0.1, 0.0, -0.1],
    ], order='F')
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    assert_error(lambda: seeds(vectors, kd_indices, dimension_order),
                 "Expected error for the default k_density=30 on a 5-point dataset", ERR_INVALID_INPUT)


if __name__ == "__main__":
    run_all_tests(globals().values())
