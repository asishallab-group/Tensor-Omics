"""
Test suite for TensorOmics normalization functions
Uses wrapper functions from tensoromics_functions.py with tox_ prefix
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (
    tox_normalize_by_std_dev,
    tox_quantile_normalization,
    tox_log2_transformation,
    tox_calculate_tissue_averages,
    tox_calculate_fold_changes,
    tox_normalize_unit_length
)
from test_helpers import run_all_tests


def test_normalize_unit_length():
    TOL = 1e-12

    # -------------------------------
    # Case 1: Normal vector
    # -------------------------------
    vec = np.array([3.0, 4.0, -123.0], dtype=np.float64)
    expected = vec / np.linalg.norm(vec)
    result = tox_normalize_unit_length(vec.copy())
    assert np.allclose(result, expected, atol=TOL), "normal vector: result mismatch"

    # -------------------------------
    # Case 2: Zero vector
    # -------------------------------
    vec = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    try:
        tox_normalize_unit_length(vec.copy())
        raise AssertionError("zero vector: should raise RuntimeError")
    except RuntimeError:
        pass  # expected

    # -------------------------------
    # Case 3: Already normalized
    # -------------------------------
    vec = np.array([0.6164770879765119, -0.42, 0.666], dtype=np.float64)
    result = tox_normalize_unit_length(vec.copy())
    assert np.allclose(result, vec, atol=0.0), "already normalized vector: result mismatch"

    # -------------------------------
    # Case 4: NaN
    # -------------------------------
    vec = np.array([0.6164770879765119, -0.42, np.nan], dtype=np.float64)
    try:
        tox_normalize_unit_length(vec.copy())
        raise AssertionError("NaN vector: should raise RuntimeError")
    except RuntimeError:
        pass  # expected

    # -------------------------------
    # Case 5: Infinity
    # -------------------------------
    vec = np.array([0.6164770879765119, -0.42, np.inf], dtype=np.float64)
    try:
        tox_normalize_unit_length(vec.copy())
        raise AssertionError("Infinity vector: should raise RuntimeError")
    except RuntimeError:
        pass  # expected


def test_tox_stddev_example_1():
    # Input data - need at least 10 genes for LOESS
    mat = np.array([[j + i for j in range(1, 100)] for i in range(100)], dtype=np.float64, order="F")

    # Call tox function with LOESS parameters
    span = 0.75
    degree = 2
    result = tox_normalize_by_std_dev(mat, span=span, degree=degree)

    for i in range(mat.shape[1]):
        col = mat[:, i]
        normalized_col = col / np.std(col, dtype=np.float64)
        assert np.allclose(result[:, i], normalized_col), "result doesn't match expected"


def test_tox_stddev_example_2():
    # Need at least 10 genes for LOESS
    mat = np.array([[j * 1e6 + i*1e5 for j in range(10)] for i in range(10)], dtype=np.float64)

    span = 0.75
    degree = 2
    result = tox_normalize_by_std_dev(mat, span=span, degree=degree)

    for i in range(mat.shape[1]):
        col = mat[:, i]
        normalized_col = col / np.std(col, dtype=np.float64)
        assert np.allclose(result[:, i], normalized_col), "result doesn't match expected"


def test_tox_quantile_example_1():
    mat = np.array([[3.0, 1.0, 2.0], [6.0, 4.0, 5.0]], dtype=np.float64, order="F")

    result = tox_quantile_normalization(mat)

    # Check column distributions
    expected_col = np.array([2.5, 3.5, 4.5], dtype=np.float64)
    for j in range(mat.shape[0]):
        sorted_col = np.sort(result[j, :])
        assert np.allclose(sorted_col, expected_col), "result doesn't match expected"


def test_tox_log2_example_1():
    """Example 1: Simple log2 transformation"""
    print("="*50)
    print("TOX_LOG2 TRANSFORMATION EXAMPLE 1: Simple values")
    print("="*50)

    # Input data: [0, 3, 7, 15] → log2(x+1) = [0, 2, 3, 4]
    mat = np.array([[0.0, 3.0],
                    [7.0, 15.0]], dtype=np.float64)


    result = tox_log2_transformation(mat)

    # Manual verification
    print("\nManual verification (log2(x+1)):")
    expected = np.log2(mat + 1)

    match = np.allclose(result, expected)
    print(f"Match? {match}")
    assert match, "Results should match expected log2(x+1)"
    print("✓ Test passed!")


def test_tox_log2_example_2():
    """Example 2: Edge cases for log2"""
    print("="*50)
    print("TOX_LOG2 TRANSFORMATION EXAMPLE 2: Edge cases")
    print("="*50)

    # Edge cases: zeros, ones, large values
    mat = np.array([[0.0, 1.0, 1023.0]], dtype=np.float64)


    result = tox_log2_transformation(mat)

    print(f"log2(0+1)={result[0,0]:.3f}, log2(1+1)={result[0,1]:.3f}, log2(1023+1)={result[0,2]:.3f}")

    # Verify expected values
    expected_vals = [0.0, 1.0, 10.0]  # log2(1), log2(2), log2(1024)
    for i, expected in enumerate(expected_vals):
        assert abs(result[0,i] - expected) < 1e-10, f"Value {i} should be {expected}"

    print("✓ Test passed!")


def test_tox_calc_tiss_avg_example_1():
    """Example 1: Average tissue replicates"""
    print("="*50)
    print("TOX_TISSUE AVERAGING EXAMPLE 1: 3 tissues, 2 reps each")
    print("="*50)

    # 2 genes × 6 samples (3 tissues, 2 replicates each)
    mat = np.array([[1.0, 7.0, 3.0, 9.0, 5.0, 11.0],   # Gene 1: samples 1-6
                    [2.0, 8.0, 4.0, 10.0, 6.0, 12.0]], dtype=np.float64)  # Gene 2: samples 1-6

    print("Tissue groups: [1,3,5] start columns, [2,2,2] replicates each")

    group_starts = np.array([1, 3, 5], dtype=np.int32)  # 1-based indexing for Fortran
    group_counts = np.array([2, 2, 2], dtype=np.int32)

    result = tox_calculate_tissue_averages(mat, group_counts)

    # Manual verification
    print("\nManual verification:")
    print(f"Tissue 1 (samples 1-2): Gene1=mean([1,7])={np.mean([1,7]):.1f}, Gene2=mean([2,8])={np.mean([2,8]):.1f}")
    print(f"Tissue 2 (samples 3-4): Gene1=mean([3,9])={np.mean([3,9]):.1f}, Gene2=mean([4,10])={np.mean([4,10]):.1f}")
    print(f"Tissue 3 (samples 5-6): Gene1=mean([5,11])={np.mean([5,11]):.1f}, Gene2=mean([6,12])={np.mean([6,12]):.1f}")

    # Verify the results match
    expected = np.array([[4.0, 6.0, 8.0], [5.0, 7.0, 9.0]])
    match = np.allclose(result, expected)
    print(f"\nResults match expected? {match}")
    assert match, "Results should match expected averages"
    print("✓ Test passed!")


def test_tox_calc_tiss_avg_example_2():
    """Example 2: Unequal replicates"""
    print("="*50)
    print("TOX_TISSUE AVERAGING EXAMPLE 2: Unequal replicates")
    print("="*50)

    # 2 genes × 7 samples (tissue1: 2 reps, tissue2: 3 reps, tissue3: 2 reps)
    mat = np.array([[1.0, 8.0, 2.0, 9.0, 3.0, 10.0, 4.0],   # Gene 1
                    [5.0, 12.0, 6.0, 13.0, 7.0, 14.0, 8.0]], dtype=np.float64)  # Gene 2

    print("Tissue groups: [1,3,6] start columns, [2,3,2] replicates each")

    group_starts = np.array([1, 3, 6], dtype=np.int32)
    group_counts = np.array([2, 3, 2], dtype=np.int32)

    result = tox_calculate_tissue_averages(mat, group_counts)

    assert np.all(np.isfinite(result)), "All results should be finite"
    print("✓ Test passed!")


def test_tox_calc_fchange_example_1():
    """Example 1: Simple fold change"""
    print("="*50)
    print("TOX_FOLD CHANGE EXAMPLE 1: Simple condition vs control")
    print("="*50)

    # 2 genes × 2 samples (control, condition)
    mat = np.array([[1.0, 2.0],   # Gene 1: control=1, condition=2
                    [4.0, 8.0]], dtype=np.float64)  # Gene 2: control=4, condition=8

    print("Control column: 1, Condition column: 2")

    control_cols = np.array([1], dtype=np.int32)  # 1-based indexing
    condition_cols = np.array([2], dtype=np.int32)

    result = tox_calculate_fold_changes(mat, control_cols, condition_cols)

    print(f"\nOutput (fold changes): {result.flatten()}")

    # Manual verification
    print("\nManual verification (condition - control):")
    print(f"Gene 1: {mat[0,1]} - {mat[0,0]} = {mat[0,1] - mat[0,0]}")
    print(f"Gene 2: {mat[1,1]} - {mat[1,0]} = {mat[1,1] - mat[1,0]}")

    # Verify expected values
    expected = np.array([[1.0], [4.0]])  # Differences: 2-1=1, 8-4=4
    match = np.allclose(result, expected)
    print(f"Results match expected? {match}")
    assert match, "Results should match expected fold changes"
    print("✓ Test passed!")


def test_tox_calc_fchange_example_2():
    """Example 2: Multiple conditions vs same control"""
    print("="*50)
    print("TOX_FOLD CHANGE EXAMPLE 2: Multiple conditions vs same control")
    print("="*50)

    # 2 genes × 3 samples (1 control, 2 conditions)
    mat = np.array([[2.0, 6.0, 10.0],   # Gene 1: control=2, condition1=6, condition2=10
                    [4.0, 16.0, 24.0]], dtype=np.float64)  # Gene 2: control=4, condition1=16, condition2=24

    print("Control column: 1, Condition columns: 2 and 3")
    print("Pairs: (control=1,condition=2) and (control=1,condition=3)")

    # Both pairs use column 1 as control, but different condition columns
    control_cols = np.array([1, 1], dtype=np.int32)      # Same control for both
    condition_cols = np.array([2, 3], dtype=np.int32)    # Different conditions

    result = tox_calculate_fold_changes(mat, control_cols, condition_cols)

    print("\nManual verification:")
    print(f"Pair 1 (condition2-control1): Gene1={mat[0,1]}-{mat[0,0]}={mat[0,1]-mat[0,0]}, Gene2={mat[1,1]}-{mat[1,0]}={mat[1,1]-mat[1,0]}")
    print(f"Pair 2 (condition3-control1): Gene1={mat[0,2]}-{mat[0,0]}={mat[0,2]-mat[0,0]}, Gene2={mat[1,2]}-{mat[1,0]}={mat[1,2]-mat[1,0]}")

    # Expected results
    expected = np.array([[4.0, 8.0], [12.0, 20.0]])
    match = np.allclose(result, expected)
    print(f"\nResults match expected? {match}")
    assert match, "Results should match expected fold changes"
    print("✓ Test passed!")


def test_error_handling():
    """Test error handling"""
    print("="*50)
    print("TOX_ERROR HANDLING TESTS")
    print("="*50)

    # Test NaN input
    print("\n--- Test NaN input ---")
    mat_nan = np.array([[1.0, np.nan], [3.0, 4.0]], dtype=np.float64)

    try:
        tox_normalize_by_std_dev(mat_nan, span=0.75, degree=2)
        assert False, "Should have raised ValueError for NaN input"
    except ValueError as e:
        print(f"✓ Correctly caught NaN error: {e}")

    # Test infinite input
    print("\n--- Test infinite input ---")
    mat_inf = np.array([[1.0, np.inf], [3.0, 4.0]], dtype=np.float64)

    try:
        tox_normalize_by_std_dev(mat_inf, span=0.75, degree=2)
        assert False, "Should have raised ValueError for infinite input"
    except ValueError as e:
        print(f"✓ Correctly caught infinite error: {e}")

    # Test mismatched arrays
    print("\n--- Test mismatched arrays ---")
    mat = np.array([[1.0, 2.0]], dtype=np.float64)
    group_starts = np.array([1], dtype=np.int32)
    group_counts = np.array([1, 2], dtype=np.int32)  # Wrong length

    try:
        tox_calculate_tissue_averages(mat, group_counts)
        assert False, "Should have raised ValueError for mismatched arrays"
    except ValueError as e:
        print(f"✓ Correctly caught mismatch error: {e}")

    print("✓ All error handling tests passed!")


if __name__ == '__main__':
    run_all_tests([test_tox_stddev_example_1, test_tox_stddev_example_2, test_tox_quantile_example_1])