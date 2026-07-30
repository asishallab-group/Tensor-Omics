"""
Test suite for TensorOmics normalization functions
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensor_omics
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    normalize_by_std_dev,
    quantile_normalization,
    log2_transformation,
    calc_tiss_avg,
    calc_fchange,
    normalize_unit_length
)
from test_helpers import run_all_tests, assert_error


def test_normalize_unit_length():
    TOL = 1e-12

    # -------------------------------
    # Case 1: Normal vector
    # -------------------------------
    vec = np.array([3.0, 4.0, -123.0], dtype=np.float64)
    expected = vec / np.linalg.norm(vec)
    normalize_unit_length(vec)
    assert np.allclose(vec, expected, atol=TOL), "normal vector: result mismatch"

    # -------------------------------
    # Case 2: Zero vector
    # -------------------------------
    vec = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    assert_error(lambda: normalize_unit_length(vec), "zero vector: should raise RuntimeError")

    # -------------------------------
    # Case 3: Already normalized
    # -------------------------------
    vec = np.array([0.6164770879765119, -0.42, 0.666], dtype=np.float64)
    expected = vec.copy()
    normalize_unit_length(vec)
    assert np.allclose(vec, expected, atol=0.0), "already normalized vector: result mismatch"

    # -------------------------------
    # Case 4: NaN
    # -------------------------------
    vec = np.array([0.6164770879765119, -0.42, np.nan], dtype=np.float64)
    assert_error(lambda: normalize_unit_length(vec), "NaN vector: should raise RuntimeError")

    # -------------------------------
    # Case 5: Infinity
    # -------------------------------
    vec = np.array([0.6164770879765119, -0.42, np.inf], dtype=np.float64)
    assert_error(lambda: normalize_unit_length(vec), "Infinity vector: should raise RuntimeError")


def test_stddev_example_1():
    # Input data - need at least 10 genes for LOESS
    mat = np.array([[j + i for j in range(1, 100)] for i in range(100)], dtype=np.float64, order="F")

    # Call tox function with LOESS parameters
    span = 0.75
    degree = 2
    result = normalize_by_std_dev(mat, span=span, degree=degree)

    for i in range(mat.shape[1]):
        col = mat[:, i]
        normalized_col = col / np.std(col, dtype=np.float64)
        assert np.allclose(result[:, i], normalized_col), "result doesn't match expected"


def test_stddev_example_2():
    # Need at least 10 genes for LOESS
    mat = np.array([[j * 1e6 + i*1e5 for j in range(10)] for i in range(10)], dtype=np.float64)

    span = 0.75
    degree = 2
    result = normalize_by_std_dev(mat, span=span, degree=degree)

    for i in range(mat.shape[1]):
        col = mat[:, i]
        normalized_col = col / np.std(col, dtype=np.float64)
        assert np.allclose(result[:, i], normalized_col), "result doesn't match expected"


def test_quantile_example_1():
    mat = np.array([[3.0, 1.0, 2.0], [6.0, 4.0, 5.0]], dtype=np.float64, order="F")

    result = quantile_normalization(mat)["normalized_expr"]

    # Check column distributions
    expected_col = np.array([2.5, 3.5, 4.5], dtype=np.float64)
    for j in range(mat.shape[0]):
        sorted_col = np.sort(result[j, :])
        assert np.allclose(sorted_col, expected_col), "result doesn't match expected"


def test_log2_example_1():
    """Example 1: Simple log2 transformation"""

    # Input data: [0, 3, 7, 15] → log2(x+1) = [0, 2, 3, 4]
    mat = np.array([[0.0, 3.0],
                    [7.0, 15.0]], dtype=np.float64)

    result = log2_transformation(mat)

    # Manual verification
    expected = np.log2(mat + 1)

    match = np.allclose(result, expected)
    assert match, "Results should match expected log2(x+1)"


def test_log2_example_2():
    """Example 2: Edge cases for log2"""

    # Edge cases: zeros, ones, large values
    mat = np.array([[0.0, 1.0, 1023.0]], dtype=np.float64)

    result = log2_transformation(mat)

    # Verify expected values
    expected_vals = [0.0, 1.0, 10.0]  # log2(1), log2(2), log2(1024)
    for i, expected in enumerate(expected_vals):
        assert abs(result[0,i] - expected) < 1e-10, f"Value {i} should be {expected}"


def test_calc_tiss_avg_example_1():
    """Example 1: Average tissue replicates"""

    # 2 genes × 6 samples (3 tissues, 2 replicates each)
    mat = np.array([[1.0, 7.0, 3.0, 9.0, 5.0, 11.0],   # Gene 1: samples 1-6
                    [2.0, 8.0, 4.0, 10.0, 6.0, 12.0]], dtype=np.float64).transpose().copy(order="F")  # Gene 2: samples 1-6

    group_starts = np.array([1, 3, 5], dtype=np.int32)  # 1-based indexing for Fortran
    group_counts = np.array([2, 2, 2], dtype=np.int32)

    result = calc_tiss_avg(group_counts, mat)

    # Manual verification

    # Verify the results match
    expected = np.array([[4.0, 6.0, 8.0], [5.0, 7.0, 9.0]]).transpose()
    match = np.allclose(result, expected)
    assert match, "Results should match expected averages"


def test_calc_tiss_avg_example_2():
    """Example 2: Unequal replicates"""

    # 2 genes × 7 samples (tissue1: 2 reps, tissue2: 3 reps, tissue3: 2 reps)
    mat = np.array([[1.0, 8.0, 2.0, 9.0, 3.0, 10.0, 4.0],   # Gene 1
                    [5.0, 12.0, 6.0, 13.0, 7.0, 14.0, 8.0]], dtype=np.float64).transpose().copy(order="F")  # Gene 2

    group_starts = np.array([1, 3, 6], dtype=np.int32)
    group_counts = np.array([2, 3, 2], dtype=np.int32)

    result = calc_tiss_avg(group_counts, mat)

    assert np.all(np.isfinite(result)), "All results should be finite"


def test_calc_fchange_example_1():
    """Example 1: Simple fold change"""

    # 2 genes × 2 samples (control, condition)
    mat = np.array([[1.0, 2.0],   # Gene 1: control=1, condition=2
                    [4.0, 8.0]], dtype=np.float64).transpose().copy(order="F")  # Gene 2: control=4, condition=8

    control_cols = np.array([1], dtype=np.int32)  # 1-based indexing
    condition_cols = np.array([2], dtype=np.int32)

    result = calc_fchange(control_cols, condition_cols, mat)

    # Verify expected values
    expected = np.array([[1.0], [4.0]]).transpose()  # Differences: 2-1=1, 8-4=4
    match = np.allclose(result, expected)
    assert match, "Results should match expected fold changes"


def test_calc_fchange_example_2():
    """Example 2: Multiple conditions vs same control"""

    # 2 genes × 3 samples (1 control, 2 conditions)
    mat = np.array([[2.0, 6.0, 10.0],   # Gene 1: control=2, condition1=6, condition2=10
                    [4.0, 16.0, 24.0]], dtype=np.float64).transpose().copy(order="F")  # Gene 2: control=4, condition1=16, condition2=24

    # Both pairs use column 1 as control, but different condition columns
    control_cols = np.array([1, 1], dtype=np.int32)      # Same control for both
    condition_cols = np.array([2, 3], dtype=np.int32)    # Different conditions

    result = calc_fchange(control_cols, condition_cols, mat)

    # Expected results
    expected = np.array([[4.0, 8.0], [12.0, 20.0]]).transpose()
    match = np.allclose(result, expected)
    assert match, "Results should match expected fold changes"


def test_error_handling():
    """Test error handling"""

    # Test NaN input
    mat_nan = np.array([[1.0, np.nan], [3.0, 4.0]], dtype=np.float64)

    assert_error(lambda: normalize_by_std_dev(mat_nan, span=0.75, degree=2), "Should have raised ValueError for NaN input")

    # Test infinite input
    mat_inf = np.array([[1.0, np.inf], [3.0, 4.0]], dtype=np.float64)

    assert_error(lambda: normalize_by_std_dev(mat_inf, span=0.75, degree=2), "Should have raised ValueError for infinite input")


if __name__ == '__main__':
    run_all_tests(globals().values())