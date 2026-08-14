"""
Python test suite for filter_ensembles_by_stop_condition / filter_ensembles_by_dimension /
filter_ensembles_by_var_explained / filter_ensembles (tox_shape_truthful_clustering_filter),
mirroring test/mod_test_shape_truthful_clustering_filter.F90
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (filter_ensembles_by_stop_condition, filter_ensembles_by_dimension,
                          filter_ensembles_by_var_explained, filter_ensembles)
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_INVALID_INPUT

STOP_REASON_MAX_SIZE = 1
STOP_REASON_REJECTED_AFTER_STABLE = 2
STOP_REASON_REJECTED_IMMEDIATELY = 3
STOP_REASON_FIXED_POINT = 4


# =====================
# filter_ensembles_by_stop_condition
# =====================
def test_filter_stop_condition_exact_exclusion():
    stop_reason = np.array([STOP_REASON_MAX_SIZE, STOP_REASON_REJECTED_AFTER_STABLE,
                            STOP_REASON_REJECTED_IMMEDIATELY, STOP_REASON_FIXED_POINT], dtype=np.int32)
    allowed = np.array([True, True, False, True], dtype=np.bool_)

    eligible = filter_ensembles_by_stop_condition(stop_reason, allowed)

    assert eligible[0]
    assert eligible[1]
    assert not eligible[2]
    assert eligible[3]


def test_filter_stop_condition_absent_is_noop():
    stop_reason = np.array([STOP_REASON_MAX_SIZE, STOP_REASON_REJECTED_AFTER_STABLE,
                            STOP_REASON_REJECTED_IMMEDIATELY, STOP_REASON_FIXED_POINT], dtype=np.int32)

    eligible = filter_ensembles_by_stop_condition(stop_reason)

    assert np.all(eligible)


def test_filter_stop_condition_all_four_values():
    for r in range(1, 5):
        stop_reason = np.array([r], dtype=np.int32)
        allowed = np.ones(4, dtype=np.bool_)
        allowed[r - 1] = False

        eligible = filter_ensembles_by_stop_condition(stop_reason, allowed)

        assert not eligible[0], f"disallowing value {r} should exclude it"


def test_filter_stop_condition_invalid_value():
    stop_reason = np.array([5], dtype=np.int32)
    assert_error(lambda: filter_ensembles_by_stop_condition(stop_reason),
                 "an out-of-range ensemble_stop_reason value must be rejected", ERR_INVALID_INPUT)


def test_filter_stop_condition_zero_ensembles():
    stop_reason = np.array([], dtype=np.int32)
    eligible = filter_ensembles_by_stop_condition(stop_reason)
    assert eligible.shape[0] == 0


# =====================
# filter_ensembles_by_dimension
# =====================
def test_filter_dimension_d_min_only():
    d_final = np.array([0, 1, 2, 3], dtype=np.int32)
    has_final = np.ones(4, dtype=np.bool_)

    eligible = filter_ensembles_by_dimension(3, d_final, has_final, d_min=2)

    assert not eligible[0]
    assert not eligible[1]
    assert eligible[2]
    assert eligible[3]


def test_filter_dimension_d_max_only():
    d_final = np.array([0, 1, 2, 3], dtype=np.int32)
    has_final = np.ones(4, dtype=np.bool_)

    eligible = filter_ensembles_by_dimension(3, d_final, has_final, d_max=1)

    assert eligible[0]
    assert eligible[1]
    assert not eligible[2]
    assert not eligible[3]


def test_filter_dimension_both_bounds():
    d_final = np.array([0, 1, 2, 3], dtype=np.int32)
    has_final = np.ones(4, dtype=np.bool_)

    eligible = filter_ensembles_by_dimension(3, d_final, has_final, d_min=1, d_max=2)

    assert not eligible[0]
    assert eligible[1]
    assert eligible[2]
    assert not eligible[3]


def test_filter_dimension_both_absent_is_noop():
    d_final = np.array([0, 3], dtype=np.int32)
    has_final = np.array([True, False], dtype=np.bool_)

    eligible = filter_ensembles_by_dimension(3, d_final, has_final)

    assert np.all(eligible)


def test_filter_dimension_no_final_excluded_once_bound_present():
    d_final = np.array([1], dtype=np.int32)
    has_final = np.array([False], dtype=np.bool_)

    eligible = filter_ensembles_by_dimension(3, d_final, has_final, d_min=0, d_max=3)

    assert not eligible[0]


def test_filter_dimension_invalid_n_dimensions():
    d_final = np.array([0], dtype=np.int32)
    has_final = np.array([True], dtype=np.bool_)
    assert_error(lambda: filter_ensembles_by_dimension(1, d_final, has_final, d_min=0),
                 "n_dimensions=1 must be rejected (minimum is 2)", ERR_INVALID_INPUT)


def test_filter_dimension_d_min_exceeds_d_max_still_computes():
    # d_min > d_max is not itself validated -- simply unsatisfiable, so every ensemble ends up
    # ineligible.
    d_final = np.array([0, 1, 2], dtype=np.int32)
    has_final = np.ones(3, dtype=np.bool_)

    eligible = filter_ensembles_by_dimension(2, d_final, has_final, d_min=2, d_max=1)

    assert not np.any(eligible)


# =====================
# filter_ensembles_by_var_explained
# =====================
def test_filter_var_explained_clean_fixture():
    # S=[10,1], k=2, d=1 -> eigenvalues [100,1] -> ve = 100/101 ~ 0.9901.
    S_final = np.array([[10.0], [1.0]], order='F')
    d_final = np.array([1], dtype=np.int32)
    k_final = np.array([2], dtype=np.int32)
    has_final = np.array([True], dtype=np.bool_)

    eligible = filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final, var_explained_min=0.9)

    assert eligible[0]


def test_filter_var_explained_threshold_at_boundary():
    S_final = np.array([[10.0], [1.0]], order='F')
    d_final = np.array([1], dtype=np.int32)
    k_final = np.array([2], dtype=np.int32)
    has_final = np.array([True], dtype=np.bool_)
    ve = 100.0 / 101.0

    assert filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final, var_explained_min=ve - 1e-9)[0]
    assert filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final, var_explained_min=ve)[0]
    assert not filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final, var_explained_min=ve + 1e-9)[0]


def test_filter_var_explained_k_le_one_guard():
    S_final = np.array([[10.0, 10.0], [1.0, 1.0]], order='F')
    d_final = np.array([1, 1], dtype=np.int32)
    k_final = np.array([1, 0], dtype=np.int32)
    has_final = np.array([True, True], dtype=np.bool_)

    eligible = filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final, var_explained_min=0.0)

    assert not eligible[0]
    assert not eligible[1]


def test_filter_var_explained_absent_is_noop():
    S_final = np.array([[10.0], [1.0]], order='F')
    d_final = np.array([1], dtype=np.int32)
    k_final = np.array([1], dtype=np.int32)  # degenerate k, irrelevant when absent
    has_final = np.array([True], dtype=np.bool_)

    eligible = filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final)

    assert eligible[0]


def test_filter_var_explained_invalid_threshold():
    S_final = np.array([[10.0], [1.0]], order='F')
    d_final = np.array([1], dtype=np.int32)
    k_final = np.array([2], dtype=np.int32)
    has_final = np.array([True], dtype=np.bool_)

    assert_error(lambda: filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final,
                                                            var_explained_min=1.5),
                 "var_explained_min > 1.0 must be rejected", ERR_INVALID_INPUT)


# =====================
# filter_ensembles (combined orchestrator)
# =====================
def test_filter_ensembles_combined_different_criteria():
    # D=2, o=1, 4 ensembles, each failing a different single criterion (or none).
    U = np.zeros((2, 2, 1, 4), order='F')
    mu = np.zeros((2, 1, 4), order='F')
    G = np.zeros((1, 4), order='F')
    k = np.full((1, 4), 2, dtype=np.int32, order='F')
    accepted = np.ones((1, 4), dtype=np.bool_, order='F')
    d = np.ones((1, 4), dtype=np.int32, order='F')
    d[0, 1] = 2  # ensemble 2 (0-indexed 1) fails d_max=1
    S = np.zeros((2, 1, 4), order='F')
    S[:, 0, 0] = [10.0, 1.0]
    S[:, 0, 1] = [10.0, 1.0]
    S[:, 0, 2] = [1.0, 10.0]  # ensemble 3 (0-indexed 2) fails var_explained_min=0.5
    S[:, 0, 3] = [10.0, 1.0]

    stop_reason = np.full(4, STOP_REASON_FIXED_POINT, dtype=np.int32)
    stop_reason[0] = STOP_REASON_REJECTED_IMMEDIATELY  # ensemble 1 fails stop condition
    allowed = np.array([True, True, False, True], dtype=np.bool_)

    result = filter_ensembles(U, d, S, mu, G, k, accepted, stop_reason,
                              allowed_stop_reasons=allowed, d_max=1, var_explained_min=0.5)

    assert not result['eligible_by_stop_condition'][0]
    assert result['eligible_by_dimension'][0]
    assert result['eligible_by_var_explained'][0]
    assert not result['eligible'][0]

    assert result['eligible_by_stop_condition'][1]
    assert not result['eligible_by_dimension'][1]
    assert result['eligible_by_var_explained'][1]
    assert not result['eligible'][1]

    assert result['eligible_by_stop_condition'][2]
    assert result['eligible_by_dimension'][2]
    assert not result['eligible_by_var_explained'][2]
    assert not result['eligible'][2]

    assert result['eligible_by_stop_condition'][3]
    assert result['eligible_by_dimension'][3]
    assert result['eligible_by_var_explained'][3]
    assert result['eligible'][3]


def test_filter_ensembles_all_omitted_is_noop():
    U = np.zeros((2, 2, 1, 2), order='F')
    d = np.zeros((1, 2), dtype=np.int32, order='F')
    S = np.zeros((2, 1, 2), order='F')
    mu = np.zeros((2, 1, 2), order='F')
    G = np.zeros((1, 2), order='F')
    k = np.zeros((1, 2), dtype=np.int32, order='F')
    accepted = np.zeros((1, 2), dtype=np.bool_, order='F')
    stop_reason = np.full(2, STOP_REASON_FIXED_POINT, dtype=np.int32)

    result = filter_ensembles(U, d, S, mu, G, k, accepted, stop_reason)

    assert np.all(result['eligible'])
    assert np.all(result['eligible_by_stop_condition'])
    assert np.all(result['eligible_by_dimension'])
    assert np.all(result['eligible_by_var_explained'])


if __name__ == "__main__":
    run_all_tests(globals().values())
