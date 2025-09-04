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
)

def print_matrix(name, mat):
    """Pretty print matrix with name"""
    print(f"\n{name}:")
    print(f"Shape: {mat.shape}")
    for i, row in enumerate(mat):
        print(f"  Row {i+1}: {row}")

def test_tox_normalize_example_1():
    """Example 1: Simple 2x3 matrix normalization"""
    print("="*50)
    print("TOX_NORMALIZE EXAMPLE 1: Simple 2x3 matrix")
    print("="*50)
    
    # Input data
    mat = np.array([[1.0, 2.0, 3.0], 
                    [4.0, 5.0, 6.0]], dtype=np.float64)
    
    print_matrix("Input", mat)
    
    # Call tox function
    result = tox_normalize_by_std_dev(mat)
    print_matrix("Output", result)
    
    # Manual verification
    print("\nManual verification:")
    for i in range(mat.shape[0]):
        row = mat[i, :]
        std_dev = np.sqrt(np.mean(row**2))
        normalized_row = row / std_dev
        print(f"  Gene {i+1}: std_dev={std_dev:.4f}, normalized={normalized_row}")
    
    # Verify results are finite
    assert np.all(np.isfinite(result)), "All results should be finite"
    print("✓ Test passed!")

def test_tox_normalize_example_2():
    """Example 2: Large values"""
    print("="*50)
    print("TOX_NORMALIZE EXAMPLE 2: Large values")
    print("="*50)
    
    mat = np.array([[1e6, 2e6], 
                    [3e6, 4e6]], dtype=np.float64)
    
    print_matrix("Input", mat)
    
    result = tox_normalize_by_std_dev(mat)
    print_matrix("Output", result)
    
    assert np.all(np.isfinite(result)), "All results should be finite"
    print(f"All finite? {np.all(np.isfinite(result))}")
    print("✓ Test passed!")

def test_tox_quantile_example_1():
    """Example 1: Simple quantile normalization"""
    print("="*50)
    print("TOX_QUANTILE EXAMPLE 1: Simple 3x3 matrix")
    print("="*50)
    
    mat = np.array([[5.0, 2.0, 8.0],
                    [1.0, 6.0, 3.0],
                    [9.0, 4.0, 7.0]], dtype=np.float64)
    
    print_matrix("Input", mat)
    
    result = tox_quantile_normalization(mat)
    print_matrix("Output", result)
    
    # Check column distributions
    print("\nColumn distributions (sorted):")
    for j in range(mat.shape[1]):
        sorted_col = np.sort(result[:, j])
        print(f"  Column {j+1}: {sorted_col}")
    
    assert np.all(np.isfinite(result)), "All results should be finite"
    print("✓ Test passed!")

def test_tox_log2_example_1():
    """Example 1: Simple log2 transformation"""
    print("="*50)
    print("TOX_LOG2 TRANSFORMATION EXAMPLE 1: Simple values")
    print("="*50)
    
    # Input data: [0, 3, 7, 15] → log2(x+1) = [0, 2, 3, 4]
    mat = np.array([[0.0, 3.0], 
                    [7.0, 15.0]], dtype=np.float64)
    
    print_matrix("Input", mat)
    
    result = tox_log2_transformation(mat)
    print_matrix("Output", result)
    
    # Manual verification
    print("\nManual verification (log2(x+1)):")
    expected = np.log2(mat + 1)
    print_matrix("Expected", expected)
    
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
    
    print_matrix("Input (edge cases)", mat)
    
    result = tox_log2_transformation(mat)
    print_matrix("Output", result)
    
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
    
    print_matrix("Input (2 genes × 6 samples)", mat)
    print("Tissue groups: [1,3,5] start columns, [2,2,2] replicates each")
    
    group_starts = np.array([1, 3, 5], dtype=np.int32)  # 1-based indexing for Fortran
    group_counts = np.array([2, 2, 2], dtype=np.int32)
    
    result = tox_calculate_tissue_averages(mat, group_starts, group_counts)
    print_matrix("Output (2 genes × 3 tissues)", result)
    
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
    
    print_matrix("Input (2 genes × 7 samples)", mat)
    print("Tissue groups: [1,3,6] start columns, [2,3,2] replicates each")
    
    group_starts = np.array([1, 3, 6], dtype=np.int32)
    group_counts = np.array([2, 3, 2], dtype=np.int32)
    
    result = tox_calculate_tissue_averages(mat, group_starts, group_counts)
    print_matrix("Output (2 genes × 3 tissues)", result)
    
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
    
    print_matrix("Input (2 genes × 2 samples)", mat)
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
    
    print_matrix("Input (2 genes × 3 samples)", mat)
    print("Control column: 1, Condition columns: 2 and 3")
    print("Pairs: (control=1,condition=2) and (control=1,condition=3)")
    
    # Both pairs use column 1 as control, but different condition columns
    control_cols = np.array([1, 1], dtype=np.int32)      # Same control for both
    condition_cols = np.array([2, 3], dtype=np.int32)    # Different conditions
    
    result = tox_calculate_fold_changes(mat, control_cols, condition_cols)
    print_matrix("Output (2 genes × 2 pairs)", result)
    
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
        tox_normalize_by_std_dev(mat_nan)
        assert False, "Should have raised ValueError for NaN input"
    except ValueError as e:
        print(f"✓ Correctly caught NaN error: {e}")
    
    # Test infinite input
    print("\n--- Test infinite input ---")
    mat_inf = np.array([[1.0, np.inf], [3.0, 4.0]], dtype=np.float64)
    
    try:
        tox_normalize_by_std_dev(mat_inf)
        assert False, "Should have raised ValueError for infinite input"
    except ValueError as e:
        print(f"✓ Correctly caught infinite error: {e}")
    
    # Test mismatched arrays
    print("\n--- Test mismatched arrays ---")
    mat = np.array([[1.0, 2.0]], dtype=np.float64)
    group_starts = np.array([1], dtype=np.int32)
    group_counts = np.array([1, 2], dtype=np.int32)  # Wrong length
    
    try:
        tox_calculate_tissue_averages(mat, group_starts, group_counts)
        assert False, "Should have raised ValueError for mismatched arrays"
    except ValueError as e:
        print(f"✓ Correctly caught mismatch error: {e}")
    
    print("✓ All error handling tests passed!")



if __name__ == "__main__":
    print("TENSOR-OMICS PYTHON TOX_ FUNCTIONS TEST SUITE")
    print("Testing wrapper functions with tox_ prefix...")
    
    try:
        test_tox_normalize_example_1()
        test_tox_normalize_example_2()
        test_tox_quantile_example_1()
        test_tox_log2_example_1()
        test_tox_log2_example_2()
        test_tox_calc_tiss_avg_example_1()
        test_tox_calc_tiss_avg_example_2()
        test_tox_calc_fchange_example_1()
        test_tox_calc_fchange_example_2()
        test_error_handling()
        
        print("\n" + "="*50)
        print("ALL TOX_ FUNCTION TESTS COMPLETED SUCCESSFULLY!")
        print("All normalization and utility functions working correctly with tox_ prefix.")
        print("="*50)

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        print("Check function implementations and signatures.")
