"""
Python test suite for ensemble_identification (tox_shape_truthful_clustering), mirroring
test/mod_test_shape_truthful_clustering.F90
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import ensemble_identification, build_kd_index
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_INVALID_INPUT

# tox_shape_truthful_clustering_kernel's STOP_REASON_*/MEMBER_ADDED_AT_STEP_* parameters are
# not exposed as importable constants (no existing codegen mechanism does that for a plain
# result-code output, unlike a "mode" input) -- these mirror their literal Fortran values.
STOP_REASON_MAX_SIZE = 1
STOP_REASON_REJECTED_AFTER_STABLE = 2
STOP_REASON_REJECTED_IMMEDIATELY = 3
STOP_REASON_FIXED_POINT = 4
STOP_REASON_ERROR = 0
MEMBER_ADDED_AT_STEP_NON_MEMBER = -1
MEMBER_ADDED_AT_STEP_SEED = 0


def _fixture_a():
    """D=2, N=7. A 5-point line (0,0)..(4,0), plus two far-away points a growth radius of
    1.0 (k_min=1 from the seed's nearest neighbor at distance 1.0) never reaches."""
    vectors = np.zeros((2, 7), dtype=np.float64, order='F')
    vectors[0, 0:5] = np.arange(5, dtype=np.float64)
    vectors[:, 5] = [0.0, 1.5]
    vectors[:, 6] = [0.0, 3.0]
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    return vectors, kd_indices, dimension_order


def _fixture_b():
    """D=3, N=7. A 5-point x-axis line plus a branch at (1,1,0),(1,2,0) next to the 2nd
    x-axis point -- the branch is swept in already at growth iteration t=2."""
    vectors = np.zeros((3, 7), dtype=np.float64, order='F')
    vectors[0, 0:5] = np.arange(5, dtype=np.float64)
    vectors[:, 5] = [1.0, 1.0, 0.0]
    vectors[:, 6] = [1.0, 2.0, 0.0]
    dimension_order = np.array([1, 2, 3], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    return vectors, kd_indices, dimension_order


def _fixture_c():
    """D=3, N=7. Same idea as fixture b, but the branch sits next to the 3rd x-axis point,
    so it is only swept in at growth iteration t=3, after 2 accepted iterations."""
    vectors = np.zeros((3, 7), dtype=np.float64, order='F')
    vectors[0, 0:5] = np.arange(5, dtype=np.float64)
    vectors[:, 5] = [2.0, 1.0, 0.0]
    vectors[:, 6] = [2.0, 2.0, 0.0]
    dimension_order = np.array([1, 2, 3], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    return vectors, kd_indices, dimension_order


def test_natural_fixed_point():
    vectors, kd_indices, dimension_order = _fixture_a()
    result = ensemble_identification(vectors, kd_indices, dimension_order, 1, 0.1, 0, 1.0e10, 4,
                                     k_min=1)

    assert result['stop_reason'] == STOP_REASON_FIXED_POINT
    assert abs(result['growth_radius'] - 1.0) < 1e-9

    expected_mask = np.zeros(7, dtype=np.bool_)
    expected_mask[0:5] = True
    assert np.array_equal(result['final_ensemble_mask'], expected_mask)

    assert np.array_equal(result['k_history'], np.array([2, 3, 4, 5], dtype=np.int32))
    assert np.all(result['accepted_history'])
    assert np.all(result['d_history'] == 1)

    expected_step = np.array([MEMBER_ADDED_AT_STEP_SEED, 1, 2, 3, 4, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                              MEMBER_ADDED_AT_STEP_NON_MEMBER], dtype=np.int32)
    assert np.array_equal(result['member_added_at_step'], expected_step)


def test_history_window_shifts():
    vectors, kd_indices, dimension_order = _fixture_a()
    result = ensemble_identification(vectors, kd_indices, dimension_order, 1, 0.1, 0, 1.0e10, 2,
                                     k_min=1)

    assert result['stop_reason'] == STOP_REASON_FIXED_POINT
    assert np.array_equal(result['k_history'], np.array([4, 5], dtype=np.int32))
    assert np.all(result['accepted_history'])


def test_max_size_at_bootstrap():
    vectors, kd_indices, dimension_order = _fixture_a()
    result = ensemble_identification(vectors, kd_indices, dimension_order, 1, 0.1, 0, 1.0e10, 3,
                                     k_min=1, f_max=0.2)

    assert result['stop_reason'] == STOP_REASON_MAX_SIZE
    assert not np.any(result['final_ensemble_mask'])
    assert np.array_equal(result['k_history'], np.zeros(3, dtype=np.int32))
    assert not np.any(result['accepted_history'])
    assert np.all(result['member_added_at_step'] == MEMBER_ADDED_AT_STEP_NON_MEMBER)


def test_max_size_poisons_prior_accepts():
    vectors, kd_indices, dimension_order = _fixture_a()
    result = ensemble_identification(vectors, kd_indices, dimension_order, 1, 0.1, 0, 1.0e10, 3,
                                     k_min=1, f_max=0.35)

    assert result['stop_reason'] == STOP_REASON_MAX_SIZE
    assert not np.any(result['final_ensemble_mask'])
    assert np.array_equal(result['k_history'], np.zeros(3, dtype=np.int32))
    assert np.all(result['member_added_at_step'] == MEMBER_ADDED_AT_STEP_NON_MEMBER)


def test_rejected_immediately():
    vectors, kd_indices, dimension_order = _fixture_b()
    result = ensemble_identification(vectors, kd_indices, dimension_order, 1, 0.1, 0, 1.0e10, 2,
                                     k_min=1)

    assert result['stop_reason'] == STOP_REASON_REJECTED_IMMEDIATELY

    expected_mask = np.zeros(7, dtype=np.bool_)
    expected_mask[0:2] = True
    assert np.array_equal(result['final_ensemble_mask'], expected_mask)

    assert np.array_equal(result['k_history'], np.array([2, 4], dtype=np.int32))
    assert np.array_equal(result['accepted_history'], np.array([True, False]))

    expected_step = np.array([MEMBER_ADDED_AT_STEP_SEED, 1, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                              MEMBER_ADDED_AT_STEP_NON_MEMBER, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                              MEMBER_ADDED_AT_STEP_NON_MEMBER, MEMBER_ADDED_AT_STEP_NON_MEMBER], dtype=np.int32)
    assert np.array_equal(result['member_added_at_step'], expected_step)


def test_rejected_after_stable():
    vectors, kd_indices, dimension_order = _fixture_c()
    result = ensemble_identification(vectors, kd_indices, dimension_order, 1, 0.1, 0, 1.0e10, 3,
                                     k_min=1)

    assert result['stop_reason'] == STOP_REASON_REJECTED_AFTER_STABLE

    expected_mask = np.zeros(7, dtype=np.bool_)
    expected_mask[0:3] = True
    assert np.array_equal(result['final_ensemble_mask'], expected_mask)

    assert np.array_equal(result['k_history'], np.array([2, 3, 5], dtype=np.int32))
    assert np.array_equal(result['accepted_history'], np.array([True, True, False]))

    expected_step = np.array([MEMBER_ADDED_AT_STEP_SEED, 1, 2, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                              MEMBER_ADDED_AT_STEP_NON_MEMBER, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                              MEMBER_ADDED_AT_STEP_NON_MEMBER], dtype=np.int32)
    assert np.array_equal(result['member_added_at_step'], expected_step)


def test_seed_index_out_of_range():
    vectors, kd_indices, dimension_order = _fixture_a()
    assert_error(lambda: ensemble_identification(vectors, kd_indices, dimension_order, 8, 0.1, 0, 1.0e10, 3,
                                                 k_min=1),
                 "Expected error for seed_index > n_vectors", ERR_INVALID_INPUT)


def test_o_zero():
    vectors, kd_indices, dimension_order = _fixture_a()
    assert_error(lambda: ensemble_identification(vectors, kd_indices, dimension_order, 1, 0.1, 0, 1.0e10, 0,
                                                 k_min=1),
                 "Expected error for o=0", ERR_INVALID_INPUT)


def test_n_dimensions_too_small():
    vectors = np.zeros((1, 7), dtype=np.float64, order='F')
    vectors[0, :] = np.arange(7, dtype=np.float64)
    dimension_order = np.array([1], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)

    assert_error(lambda: ensemble_identification(vectors, kd_indices, dimension_order, 1, 0.1, 0, 1.0e10, 3,
                                                 k_min=1),
                 "Expected error for n_dimensions=1", ERR_INVALID_INPUT)


if __name__ == "__main__":
    run_all_tests(globals().values())
