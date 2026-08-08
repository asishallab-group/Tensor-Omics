"""
Python test suite for calc_ensemble_growth_radius / grow_ensemble
(tox_shape_truthful_clustering_ensemble_growing), mirroring
test/mod_test_shape_truthful_clustering_ensemble_growing.F90
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import calc_ensemble_growth_radius, grow_ensemble, build_kd_index
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
# calc_ensemble_growth_radius
# =====================
# Seed at x=5 (index 6, 0-based index 5). Its 4 nearest neighbors are x=4,6 (distance 1) and
# x=3,7 (distance 2), sorted [1,1,2,2] -- even k_min=4, median = avg(2nd,3rd) = 1.5. Its 3
# nearest are x=4,6 (distance 1) and one of x=3/x=7 (distance 2, a tie), sorted [1,1,2] --
# odd k_min=3, median = the middle element = 1.
def test_growth_radius_even_k():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    radius = calc_ensemble_growth_radius(vectors, kd_indices, dimension_order, 6, k_min=4)
    assert abs(radius - 1.5) < 1e-9, f"expected 1.5, got {radius}"


def test_growth_radius_odd_k():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    radius = calc_ensemble_growth_radius(vectors, kd_indices, dimension_order, 6, k_min=3)
    assert abs(radius - 1.0) < 1e-9, f"expected 1.0, got {radius}"


def test_growth_radius_seed_index_out_of_range():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    assert_error(lambda: calc_ensemble_growth_radius(vectors, kd_indices, dimension_order, 12),
                 "Expected error for seed_index > n_vectors", ERR_INVALID_INPUT)


def test_growth_radius_k_min_too_large():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    assert_error(lambda: calc_ensemble_growth_radius(vectors, kd_indices, dimension_order, 6, k_min=11),
                 "Expected error for k_min > n_vectors-1", ERR_INVALID_INPUT)


# The Fortran suite's test_growth_radius_omitted_k_min_is_clamped is a regression test for a
# crash that could only happen with k_min truly *absent* at the Fortran ABI boundary -- not
# reproducible here, since the Python binding always resolves and passes k_min=30 explicitly
# (see misc/code_gen_footgun.md's third entry). What is worth covering from Python: that this
# always-explicit default still gets validated normally on a dataset smaller than it.
def test_growth_radius_default_k_min_too_large_for_dataset():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    assert_error(lambda: calc_ensemble_growth_radius(vectors, kd_indices, dimension_order, 6),
                 "Expected error for the default k_min=30 on an 11-point dataset", ERR_INVALID_INPUT)


# Same k_min=4 fixture, distances sorted [1,1,2,2]. radius_percentile=0.0 is the smallest
# value in sorted order (the nearest-neighbor distance, 1.0); radius_percentile=100.0 is the
# largest (the farthest of the k_min neighbors, 2.0) -- irrespective of k_min's parity.
def test_growth_radius_percentile_min():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    radius = calc_ensemble_growth_radius(vectors, kd_indices, dimension_order, 6, k_min=4, radius_percentile=0.0)
    assert abs(radius - 1.0) < 1e-9, f"expected 1.0, got {radius}"


def test_growth_radius_percentile_max():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    radius = calc_ensemble_growth_radius(vectors, kd_indices, dimension_order, 6, k_min=4, radius_percentile=100.0)
    assert abs(radius - 2.0) < 1e-9, f"expected 2.0, got {radius}"


def test_growth_radius_invalid_percentile():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    assert_error(lambda: calc_ensemble_growth_radius(vectors, kd_indices, dimension_order, 6,
                                                      k_min=4, radius_percentile=101.0),
                 "Expected error for radius_percentile > 100", ERR_INVALID_INPUT)


# =====================
# grow_ensemble
# =====================
# Same 11-point line fixture. With growth_radius=1.5: from a single member x=5 (index 6, 0-based
# 5): covers x=4,5,6 (0-based indices 4,5,6). From members x=4,5,6: the union covers x=3..7.
def test_grow_ensemble_single_member():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    is_member_mask = np.zeros(11, dtype=np.bool_)
    is_member_mask[5] = True

    result = grow_ensemble(vectors, kd_indices, dimension_order, is_member_mask, 1.5)
    expected = np.zeros(11, dtype=np.bool_)
    expected[4:7] = True
    assert np.array_equal(result, expected), f"expected {expected}, got {result}"


def test_grow_ensemble_multi_member_union():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    is_member_mask = np.zeros(11, dtype=np.bool_)
    is_member_mask[4:7] = True

    result = grow_ensemble(vectors, kd_indices, dimension_order, is_member_mask, 1.5)
    expected = np.zeros(11, dtype=np.bool_)
    expected[3:8] = True
    assert np.array_equal(result, expected), f"expected {expected}, got {result}"


def test_grow_ensemble_empty_ensemble_is_degenerate():
    # An all-False is_member_mask is a well-defined degenerate case (nothing to grow from),
    # not a validation error -- see the Fortran test's comment for why.
    vectors, kd_indices, dimension_order = _line_fixture(11)
    is_member_mask = np.zeros(11, dtype=np.bool_)

    result = grow_ensemble(vectors, kd_indices, dimension_order, is_member_mask, 1.5)
    assert not np.any(result), f"expected an all-False result, got {result}"


def test_grow_ensemble_negative_radius():
    vectors, kd_indices, dimension_order = _line_fixture(11)
    is_member_mask = np.zeros(11, dtype=np.bool_)
    is_member_mask[5] = True

    assert_error(lambda: grow_ensemble(vectors, kd_indices, dimension_order, is_member_mask, -1.5),
                 "Expected error for a negative growth radius", ERR_INVALID_INPUT)


if __name__ == "__main__":
    run_all_tests(globals().values())
