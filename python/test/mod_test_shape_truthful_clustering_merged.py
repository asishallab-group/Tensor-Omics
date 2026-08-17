"""
Python test suite for ensemble_identification_merged (tox_shape_truthful_clustering),
mirroring test/mod_test_shape_truthful_clustering_merged.F90
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import ensemble_identification_merged, build_kd_index
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_INVALID_INPUT

STOP_REASON_FIXED_POINT = 4
STOP_REASON_ERROR = 0
MEMBER_ADDED_AT_STEP_NON_MEMBER = -1
MEMBER_ADDED_AT_STEP_SEED = 0


def _fixture_a():
    """D=2, N=7. Matches fixture A in mod_test_shape_truthful_clustering.py exactly."""
    vectors = np.zeros((2, 7), dtype=np.float64, order='F')
    vectors[0, 0:5] = np.arange(5, dtype=np.float64)
    vectors[:, 5] = [0.0, 1.5]
    vectors[:, 6] = [0.0, 3.0]
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    return vectors, kd_indices, dimension_order


def test_single_seed_matches_per_seed_impl():
    vectors, kd_indices, dimension_order = _fixture_a()
    seed_selection_mask = np.zeros(7, dtype=np.bool_)
    seed_selection_mask[0] = True

    result = ensemble_identification_merged(vectors, kd_indices, dimension_order, seed_selection_mask,
                                            0.1, 0, 1.0e10, 1.0e10, 4, k_min=1)

    assert result['ensemble_stop_reason'][0] == STOP_REASON_FIXED_POINT
    assert abs(result['ensemble_growth_radii'][0] - 1.0) < 1e-9

    expected_mask = np.zeros(7, dtype=np.bool_)
    expected_mask[0:5] = True
    assert np.array_equal(result['ensemble_masks'][:, 0], expected_mask)

    assert np.array_equal(result['ensemble_k_history'][:, 0], np.array([2, 3, 4, 5], dtype=np.int32))
    assert np.all(result['ensemble_accepted_history'][:, 0])

    expected_step = np.array([MEMBER_ADDED_AT_STEP_SEED, 1, 2, 3, 4, MEMBER_ADDED_AT_STEP_NON_MEMBER,
                              MEMBER_ADDED_AT_STEP_NON_MEMBER], dtype=np.int32)
    assert np.array_equal(result['ensemble_member_added_at_step'][:, 0], expected_step)

    # Iteration 1's own bootstrap mask -- {seed=1, its one growth-radius neighbor=2}.
    expected_low_confidence = np.zeros(7, dtype=np.bool_)
    expected_low_confidence[0:2] = True
    assert np.array_equal(result['ensemble_low_confidence_masks'][:, 0], expected_low_confidence)


def test_two_independent_seeds():
    vectors = np.zeros((2, 14), dtype=np.float64, order='F')
    vectors[0, 0:5] = np.arange(5, dtype=np.float64)
    vectors[:, 5] = [0.0, 1.5]
    vectors[:, 6] = [0.0, 3.0]
    vectors[:, 7:14] = vectors[:, 0:7]
    vectors[0, 7:14] += 100.0
    vectors[1, 7:14] += 100.0
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)

    seed_selection_mask = np.zeros(14, dtype=np.bool_)
    seed_selection_mask[0] = True
    seed_selection_mask[7] = True

    result = ensemble_identification_merged(vectors, kd_indices, dimension_order, seed_selection_mask,
                                            0.1, 0, 1.0e10, 1.0e10, 4, k_min=1)

    assert result['ensemble_stop_reason'][0] == STOP_REASON_FIXED_POINT
    assert result['ensemble_stop_reason'][1] == STOP_REASON_FIXED_POINT

    expected_mask_1 = np.zeros(14, dtype=np.bool_)
    expected_mask_1[0:5] = True
    assert np.array_equal(result['ensemble_masks'][:, 0], expected_mask_1)

    expected_mask_2 = np.zeros(14, dtype=np.bool_)
    expected_mask_2[7:12] = True
    assert np.array_equal(result['ensemble_masks'][:, 1], expected_mask_2)

    assert np.array_equal(result['ensemble_k_history'][:, 0], np.array([2, 3, 4, 5], dtype=np.int32))
    assert np.array_equal(result['ensemble_k_history'][:, 1], np.array([2, 3, 4, 5], dtype=np.int32))

    expected_step_1 = np.full(14, MEMBER_ADDED_AT_STEP_NON_MEMBER, dtype=np.int32)
    expected_step_1[0] = MEMBER_ADDED_AT_STEP_SEED
    expected_step_1[1:5] = [1, 2, 3, 4]
    assert np.array_equal(result['ensemble_member_added_at_step'][:, 0], expected_step_1)

    expected_step_2 = np.full(14, MEMBER_ADDED_AT_STEP_NON_MEMBER, dtype=np.int32)
    expected_step_2[7] = MEMBER_ADDED_AT_STEP_SEED
    expected_step_2[8:12] = [1, 2, 3, 4]
    assert np.array_equal(result['ensemble_member_added_at_step'][:, 1], expected_step_2)

    # ensemble_U_first/ensemble_d_first: each column is its own seed's bootstrap basis,
    # collinear along the x-axis in both copies -- must not leak across columns.
    assert result['ensemble_d_first'][0] == 1
    assert result['ensemble_d_first'][1] == 1
    assert abs(abs(result['ensemble_U_first'][0, 0, 0]) - 1.0) < 1e-9
    assert abs(abs(result['ensemble_U_first'][1, 0, 0]) - 0.0) < 1e-9
    assert abs(abs(result['ensemble_U_first'][0, 0, 1]) - 1.0) < 1e-9
    assert abs(abs(result['ensemble_U_first'][1, 0, 1]) - 0.0) < 1e-9


def test_zero_seeds():
    vectors, kd_indices, dimension_order = _fixture_a()
    seed_selection_mask = np.zeros(7, dtype=np.bool_)

    result = ensemble_identification_merged(vectors, kd_indices, dimension_order, seed_selection_mask,
                                            0.1, 0, 1.0e10, 1.0e10, 4, k_min=1)
    assert result['ensemble_masks'].shape == (7, 0)


def test_n_dimensions_too_small():
    vectors = np.zeros((1, 7), dtype=np.float64, order='F')
    vectors[0, :] = np.arange(7, dtype=np.float64)
    dimension_order = np.array([1], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    seed_selection_mask = np.zeros(7, dtype=np.bool_)
    seed_selection_mask[0] = True

    assert_error(lambda: ensemble_identification_merged(vectors, kd_indices, dimension_order, seed_selection_mask,
                                                        0.1, 0, 1.0e10, 1.0e10, 4, k_min=1),
                 "Expected error for n_dimensions=1", ERR_INVALID_INPUT)


if __name__ == "__main__":
    run_all_tests(globals().values())
