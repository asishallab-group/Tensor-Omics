"""
Python test suite for calculate_density_radius / density_labels / seeds
(tox_shape_truthful_clustering_seeding), mirroring test/mod_test_shape_truthful_clustering_seeding.F90
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import calculate_density_radius, density_labels, seeds, build_kd_index
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
# calculate_density_radius
# =====================
# Shared fixture: D=2, N=5 points on a line, (0,0),(1,0),(2,0),(3,0),(4,0). Mean = (2,0),
# so the mean-to-vector distances are [2,1,0,1,2], sorted ascending [0,1,1,2,2].
# Default 15th percentile: rank = 0.15*4+1 = 1.6 -> interpolate 0 and 1 at fraction 0.6 -> 0.6.
# 50th percentile: rank = 3.0 exactly -> 1.0.
def _radius_fixture():
    vectors = np.array([[0.0, 1.0, 2.0, 3.0, 4.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
    return vectors


def test_density_radius_default_percentile():
    radius = calculate_density_radius(_radius_fixture())
    assert abs(radius - 0.6) < 1e-9, f"expected 0.6, got {radius}"


def test_density_radius_custom_percentile():
    radius = calculate_density_radius(_radius_fixture(), mean_to_other_vecs_dist_quant=0.5)
    assert abs(radius - 1.0) < 1e-9, f"expected 1.0, got {radius}"


def test_density_radius_invalid_percentile():
    assert_error(lambda: calculate_density_radius(_radius_fixture(), mean_to_other_vecs_dist_quant=1.5),
                 "Expected error for quantile > 1.0", ERR_INVALID_INPUT)


def test_density_radius_single_vector():
    vectors = np.array([[5.0], [5.0]], order='F')
    radius = calculate_density_radius(vectors)
    assert abs(radius - 0.0) < 1e-12, f"expected 0.0, got {radius}"


# =====================
# density_labels
# =====================
def test_density_labels_basic():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    labels = density_labels(vectors, kd_indices, dimension_order, 1.5)
    expected = np.full(11, 3.0)
    expected[0] = 2.0
    expected[10] = 2.0
    assert np.allclose(labels, expected, atol=1e-12), f"expected {expected}, got {labels}"


def test_density_labels_zero_radius():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    labels = density_labels(vectors, kd_indices, dimension_order, 0.0)
    assert np.allclose(labels, 1.0, atol=1e-12), f"expected all 1.0, got {labels}"


def test_density_labels_invalid_kd_indices():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    bad_kd_indices = kd_indices.copy()
    bad_kd_indices[0] = 12
    assert_error(lambda: density_labels(vectors, bad_kd_indices, dimension_order, 1.5),
                 "Expected error for kd_indices entry > n_vectors", ERR_INVALID_INPUT)


# =====================
# seeds
# =====================
def _two_clusters_fixture():
    vectors = np.array([
        [0.0, 0.1, 0.0, -0.1, 0.0, 10.0, 10.1, 10.0, 9.9, 10.0],
        [0.0, 0.0, 0.1, 0.0, -0.1, 0.0, 0.0, 0.1, 0.0, -0.1],
    ], order='F')
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    return vectors, kd_indices, dimension_order


def test_seeds_two_separated_clusters():
    vectors, kd_indices, dimension_order = _two_clusters_fixture()
    is_seed_mask = seeds(vectors, kd_indices, dimension_order)
    assert int(np.sum(is_seed_mask)) == 2, f"expected 2 seeds, got {int(np.sum(is_seed_mask))}"
    assert np.any(is_seed_mask[0:5]), "cluster A (indices 0-4) should have a seed"
    assert np.any(is_seed_mask[5:10]), "cluster B (indices 5-9) should have a seed"


def _single_cluster_fixture():
    vectors = np.array([
        [0.0, 0.1, 0.0, -0.1, 0.0],
        [0.0, 0.0, 0.1, 0.0, -0.1],
    ], order='F')
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    return vectors, kd_indices, dimension_order


def test_seeds_single_cluster_one_seed():
    # At the 100th percentile the density radius equals the cluster's own diameter (the
    # farthest mean-to-vector distance), which is large enough to cover the whole cluster
    # from a single pick. At the default (15th) percentile it would not (see the Fortran
    # test's comment for why that is a real property of the algorithm, not a bug).
    vectors, kd_indices, dimension_order = _single_cluster_fixture()
    is_seed_mask = seeds(vectors, kd_indices, dimension_order, mean_to_other_vecs_dist_quant=1.0)
    assert int(np.sum(is_seed_mask)) == 1, f"expected 1 seed, got {int(np.sum(is_seed_mask))}"


def test_seeds_invalid_percentile():
    vectors, kd_indices, dimension_order = _single_cluster_fixture()
    assert_error(lambda: seeds(vectors, kd_indices, dimension_order, mean_to_other_vecs_dist_quant=-0.1),
                 "Expected error for negative mean_to_other_vecs_dist_quant", ERR_INVALID_INPUT)


if __name__ == "__main__":
    run_all_tests(globals().values())
