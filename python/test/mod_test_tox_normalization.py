"""
The `tox_normalization` suite: each published procedure can be called from Python, returns the
documented type and shape, and raises the documented error. Whether the numbers are right is the
Fortran suite's job (test/mod_test_tox_normalization.F90), as the coding guide's "Where to Test
What" asks.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    calc_fchange,
    calc_tiss_avg,
    log2_transformation,
    normalization_pipeline,
    normalize_by_std_dev,
    normalize_unit_length,
    quantile_normalization,
    root_mean_sq_normalization,
)
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import (
    ERR_DIVISION_BY_ZERO,
    ERR_EMPTY_INPUT,
    ERR_INVALID_INPUT,
    ERR_SIZE_MISMATCH,
)

N_REPLICATES, N_GENES = 6, 10
# two tissues of three replicates each
REPS_PER_TISSUE = np.array([3, 3], dtype=np.int32)


def _expr():
    """(replicates, genes), every gene's spread proportional to its level: the mean-sd trend the
    LOESS step fits, positive throughout, so `log2(x + 1)` accepts every normalized value."""
    return np.fromfunction(lambda i_replicate, i_gene: (i_gene + 1) * (10.0 + 0.5 * (i_replicate + 1)),
                           (N_REPLICATES, N_GENES))


def _assert_matrix(result, shape, name):
    assert isinstance(result, np.ndarray), f"{name}: expected a numpy array, got {type(result).__name__}"
    assert result.dtype == np.float64, f"{name}: expected float64, got {result.dtype}"
    assert result.shape == shape, f"{name}: expected shape {shape}, got {result.shape}"
    assert np.all(np.isfinite(result)), f"{name}: non-finite values in the result"


def test_normalize_unit_length():
    vector = np.array([3.0, 4.0], dtype=np.float64)
    assert normalize_unit_length(vector) is None, "the vector is modified in place, nothing is returned"
    # the only way to see the in-place write-back reach the caller's array
    assert np.isclose(np.linalg.norm(vector), 1.0), "the caller's array was not normalized in place"


def test_normalize_unit_length_rejects_zero_vector():
    assert_error(lambda: normalize_unit_length(np.zeros(3, dtype=np.float64)),
                 "a zero vector has no direction", ERR_DIVISION_BY_ZERO)


def test_normalize_unit_length_rejects_what_it_cannot_modify_in_place():
    # a list would be converted to a temporary array, and the result would be lost with it
    assert_error(lambda: normalize_unit_length([3.0, 4.0]), "a list cannot be modified in place")


def test_normalization_pipeline():
    for use_quantile in (False, True):
        result = normalization_pipeline(_expr(), REPS_PER_TISSUE, span=0.75, degree=2, use_quantile=use_quantile)
        _assert_matrix(result, (len(REPS_PER_TISSUE), N_GENES), f"use_quantile={use_quantile}")


def test_normalization_pipeline_rejects_reps_not_summing_to_the_replicates():
    assert_error(lambda: normalization_pipeline(_expr(), np.array([3, 2], dtype=np.int32)),
                 "reps_per_tissue sums to 5, but expr has 6 replicates", ERR_SIZE_MISMATCH)


def test_normalize_by_std_dev():
    _assert_matrix(normalize_by_std_dev(_expr()), (N_REPLICATES, N_GENES), "defaults")
    _assert_matrix(normalize_by_std_dev(_expr(), span=0.75, degree=1), (N_REPLICATES, N_GENES), "span and degree")


def test_normalize_by_std_dev_rejects_too_few_varying_genes():
    # the LOESS fit needs five genes that vary across replicates; four is one short
    assert_error(lambda: normalize_by_std_dev(_expr()[:, :4]),
                 "four genes are too few for the fit", ERR_INVALID_INPUT)


def test_root_mean_sq_normalization():
    _assert_matrix(root_mean_sq_normalization(_expr()), (N_REPLICATES, N_GENES), "root_mean_sq_normalization")


def test_quantile_normalization():
    result = quantile_normalization(_expr())
    assert isinstance(result, dict), f"expected a dict, got {type(result).__name__}"
    assert set(result) == {"normalized_expr", "rank_means"}, f"unexpected keys {sorted(result)}"
    _assert_matrix(result["normalized_expr"], (N_REPLICATES, N_GENES), "normalized_expr")
    _assert_matrix(result["rank_means"], (N_GENES,), "rank_means")


def test_log2_transformation():
    tissue_averages = calc_tiss_avg(REPS_PER_TISSUE, _expr())
    _assert_matrix(log2_transformation(tissue_averages), (len(REPS_PER_TISSUE), N_GENES), "log2_transformation")


def test_log2_transformation_rejects_values_at_minus_one():
    # log2(x + 1) is only defined above -1
    assert_error(lambda: log2_transformation(np.array([[1.0, -1.0]])),
                 "log2(0) is undefined", ERR_INVALID_INPUT)


def test_every_matrix_procedure_rejects_an_empty_expr():
    empty = np.empty((N_REPLICATES, 0), dtype=np.float64)
    for procedure in (normalize_by_std_dev, root_mean_sq_normalization, quantile_normalization,
                      log2_transformation):
        assert_error(lambda: procedure(empty), f"{procedure.__name__}: no genes", ERR_EMPTY_INPUT)


def test_calc_tiss_avg():
    _assert_matrix(calc_tiss_avg(REPS_PER_TISSUE, _expr()), (len(REPS_PER_TISSUE), N_GENES), "calc_tiss_avg")


def test_calc_tiss_avg_rejects_a_tissue_without_replicates():
    assert_error(lambda: calc_tiss_avg(np.array([3, 0, 3], dtype=np.int32), _expr()),
                 "a tissue with no replicates has no average", ERR_INVALID_INPUT)


def test_calc_fchange():
    tissue_averages = calc_tiss_avg(np.array([2, 2, 2], dtype=np.int32), _expr())
    # one control against two conditions
    result = calc_fchange(np.array([1, 1], dtype=np.int32), np.array([2, 3], dtype=np.int32), tissue_averages)
    _assert_matrix(result, (2, N_GENES), "calc_fchange")


def test_calc_fchange_rejects_a_tissue_past_the_last():
    tissue_averages = calc_tiss_avg(np.array([2, 2, 2], dtype=np.int32), _expr())
    assert_error(lambda: calc_fchange(np.array([1], dtype=np.int32), np.array([4], dtype=np.int32), tissue_averages),
                 "there are only three tissues", ERR_INVALID_INPUT)


def test_calc_fchange_rejects_unpaired_tissues():
    tissue_averages = calc_tiss_avg(np.array([2, 2, 2], dtype=np.int32), _expr())
    # the binding's own shape check, before the library is called
    assert_error(lambda: calc_fchange(np.array([1, 1], dtype=np.int32), np.array([2], dtype=np.int32), tissue_averages),
                 "two controls, one condition")


if __name__ == '__main__':
    run_all_tests(globals().values())
