
"""
Comprehensive Python test suite for compute_edf function
Uses tensoromics_functions.py wrapper function (mirrors Fortran test suite)
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import compute_edf, compute_edf_expert

print("=== Testing compute_edf Python wrapper function ===")
print("Based on Fortran test suite with comprehensive test coverage\n")

# =====================
# Test 1: Simple EDF Test
# =====================
def test_edf_simple():
    """Test basic EDF computation with simple dataset"""
    print("[test_edf_simple] Basic EDF computation")
    
    values = np.array([1.0, 2.0, 2.0, 3.0, 3.0, 3.0])
    
    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']
    
    print(f"  Input values: {values}")
    print(f"  Number of unique values: {n_unique}")
    print(f"  Unique values: {unique_vals[:n_unique]}")
    print(f"  CDF values: {cdf_vals[:n_unique]}")
    
    # Expected: 1.0 appears once (1/6), 2.0 appears twice (3/6), 3.0 appears thrice (6/6)
    expected_unique = [1.0, 2.0, 3.0]
    expected_cdf = [1.0/6.0, 3.0/6.0, 6.0/6.0]
    
    assert n_unique == 3, f"Expected 3 unique values, got {n_unique}"
    
    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}: expected {expected_unique[i]}, got {unique_vals[i]}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}: expected {expected_cdf[i]}, got {cdf_vals[i]}"
    
    print("  ✓ Simple EDF test passed\n")


# =====================
# Test 2: All Unique Values
# =====================
def test_edf_all_unique():
    """Test EDF when all values are unique"""
    print("[test_edf_all_unique] All values are unique")
    
    values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    
    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']
    
    print(f"  Input values: {values}")
    print(f"  Number of unique values: {n_unique}")
    print(f"  Unique values: {unique_vals[:n_unique]}")
    print(f"  CDF values: {cdf_vals[:n_unique]}")
    
    # All values are unique, so CDF should be 0.2, 0.4, 0.6, 0.8, 1.0
    expected_unique = [1.0, 2.0, 3.0, 4.0, 5.0]
    expected_cdf = [0.2, 0.4, 0.6, 0.8, 1.0]
    
    assert n_unique == 5, f"Expected 5 unique values, got {n_unique}"
    
    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}"
    
    print("  ✓ All unique values test passed\n")


# =====================
# Test 3: All Same Values
# =====================
def test_edf_all_same():
    """Test EDF when all values are identical"""
    print("[test_edf_all_same] All values are identical")
    
    values = np.array([5.0, 5.0, 5.0, 5.0])
    
    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']
    
    print(f"  Input values: {values}")
    print(f"  Number of unique values: {n_unique}")
    print(f"  Unique values: {unique_vals[:n_unique]}")
    print(f"  CDF values: {cdf_vals[:n_unique]}")
    
    # Only one unique value with CDF = 1.0
    assert n_unique == 1, f"Expected 1 unique value, got {n_unique}"
    assert abs(unique_vals[0] - 5.0) < 1e-12, "Unique value should be 5.0"
    assert abs(cdf_vals[0] - 1.0) < 1e-12, "CDF should be 1.0"
    
    print("  ✓ All same values test passed\n")


# =====================
# Test 4: Duplicates
# =====================
def test_edf_duplicates():
    """Test EDF with various duplicate values"""
    print("[test_edf_duplicates] Various duplicate values")
    
    values = np.array([1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 4.0])
    
    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']
    
    print(f"  Input values: {values}")
    print(f"  Number of unique values: {n_unique}")
    print(f"  Unique values: {unique_vals[:n_unique]}")
    print(f"  CDF values: {cdf_vals[:n_unique]}")
    
    # Expected: 1.0 (2/7), 2.0 (3/7), 3.0 (6/7), 4.0 (7/7)
    expected_unique = [1.0, 2.0, 3.0, 4.0]
    expected_cdf = [2.0/7.0, 3.0/7.0, 6.0/7.0, 7.0/7.0]
    
    assert n_unique == 4, f"Expected 4 unique values, got {n_unique}"
    
    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}"
    
    print("  ✓ Duplicates test passed\n")


# =====================
# Test 5: Single Value
# =====================
def test_edf_single_value():
    """Test EDF with a single value"""
    print("[test_edf_single_value] Single value input")
    
    values = np.array([42.0])
    
    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']
    
    print(f"  Input values: {values}")
    print(f"  Number of unique values: {n_unique}")
    print(f"  Unique values: {unique_vals[:n_unique]}")
    print(f"  CDF values: {cdf_vals[:n_unique]}")
    
    assert n_unique == 1, f"Expected 1 unique value, got {n_unique}"
    assert abs(unique_vals[0] - 42.0) < 1e-12, "Unique value should be 42.0"
    assert abs(cdf_vals[0] - 1.0) < 1e-12, "CDF should be 1.0"
    
    print("  ✓ Single value test passed\n")


# =====================
# Test 6: Empty Input
# =====================
def test_edf_empty_input():
    """Test EDF with empty array (should fail with error 202)"""
    print("[test_edf_empty_input] Empty array input (should fail)")
    
    values = np.array([])
    
    try:
        result = compute_edf(values)
        print("  ✗ ERROR: Should have raised exception for empty input")
        assert False, "Expected RuntimeError for empty input"
    except RuntimeError as e:
        print(f"  ✓ Correctly raised exception: {e}\n")


# =====================
# Test 7: Large Dataset
# =====================
def test_edf_large_dataset():
    """Test EDF with a larger dataset"""
    print("[test_edf_large_dataset] Large dataset (1000 values)")
    
    # Create dataset with known distribution
    values = np.concatenate([
        np.full(250, 1.0),
        np.full(250, 2.0),
        np.full(250, 3.0),
        np.full(250, 4.0)
    ])
    
    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']
    
    print(f"  Input size: {len(values)}")
    print(f"  Number of unique values: {n_unique}")
    print(f"  Unique values: {unique_vals[:n_unique]}")
    print(f"  CDF values: {cdf_vals[:n_unique]}")
    
    # Expected: 1.0 (0.25), 2.0 (0.5), 3.0 (0.75), 4.0 (1.0)
    expected_unique = [1.0, 2.0, 3.0, 4.0]
    expected_cdf = [0.25, 0.5, 0.75, 1.0]
    
    assert n_unique == 4, f"Expected 4 unique values, got {n_unique}"
    
    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}"
    
    print("  ✓ Large dataset test passed\n")


# =====================
# Test 8: Negative Values
# =====================
def test_edf_negative_values():
    """Test EDF with negative values"""
    print("[test_edf_negative_values] Negative values")
    
    values = np.array([-3.0, -1.0, 0.0, 1.0, 3.0])
    
    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']
    
    print(f"  Input values: {values}")
    print(f"  Number of unique values: {n_unique}")
    print(f"  Unique values: {unique_vals[:n_unique]}")
    print(f"  CDF values: {cdf_vals[:n_unique]}")
    
    # All values are unique
    expected_unique = [-3.0, -1.0, 0.0, 1.0, 3.0]
    expected_cdf = [0.2, 0.4, 0.6, 0.8, 1.0]
    
    assert n_unique == 5, f"Expected 5 unique values, got {n_unique}"
    
    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}"
    
    print("  ✓ Negative values test passed\n")


# =====================
# Test 9: Unsorted Input
# =====================
def test_edf_unsorted_input():
    """Test EDF with unsorted input (function should handle internally)"""
    print("[test_edf_unsorted_input] Unsorted input values")
    
    values = np.array([3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0])
    
    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']
    
    print(f"  Input values: {values}")
    print(f"  Number of unique values: {n_unique}")
    print(f"  Unique values: {unique_vals[:n_unique]}")
    print(f"  CDF values: {cdf_vals[:n_unique]}")
    
    # Expected sorted unique values: 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 9.0
    # 1.0 appears twice (2/8), then each once
    expected_unique = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 9.0]
    expected_cdf = [2.0/8.0, 3.0/8.0, 4.0/8.0, 5.0/8.0, 6.0/8.0, 7.0/8.0, 8.0/8.0]
    
    assert n_unique == 7, f"Expected 7 unique values, got {n_unique}"
    
    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}: expected {expected_unique[i]}, got {unique_vals[i]}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}: expected {expected_cdf[i]}, got {cdf_vals[i]}"
    
    print("  ✓ Unsorted input test passed\n")


# =====================
# Test 10: Python List Input
# =====================
def test_edf_list_input():
    """Test EDF with Python list instead of numpy array"""
    print("[test_edf_list_input] Python list as input")
    
    values = [1.0, 2.0, 2.0, 3.0]
    
    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']
    
    print(f"  Input values: {values}")
    print(f"  Number of unique values: {n_unique}")
    print(f"  Unique values: {unique_vals[:n_unique]}")
    print(f"  CDF values: {cdf_vals[:n_unique]}")
    
    expected_unique = [1.0, 2.0, 3.0]
    expected_cdf = [0.25, 0.75, 1.0]
    
    assert n_unique == 3, f"Expected 3 unique values, got {n_unique}"
    
    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}"
    
    print("  ✓ Python list input test passed\n")


# =====================
# Run all tests
# =====================
if __name__ == "__main__":
    test_count = 0
    failed_count = 0
    
    tests = [
        test_edf_simple,
        test_edf_all_unique,
        test_edf_all_same,
        test_edf_duplicates,
        test_edf_single_value,
        test_edf_empty_input,
        test_edf_large_dataset,
        test_edf_negative_values,
        test_edf_unsorted_input,
        test_edf_list_input
    ]
    
    for test in tests:
        test_count += 1
        try:
            test()
        except AssertionError as e:
            failed_count += 1
            print(f"  ✗ Test failed: {e}\n")
        except Exception as e:
            failed_count += 1
            print(f"  ✗ Unexpected error: {e}\n")
    
    print("=" * 60)
    print(f"Test Summary: {test_count - failed_count}/{test_count} tests passed")
    if failed_count == 0:
        print("All tests passed! ✓")
    else:
        print(f"{failed_count} test(s) failed")
    print("=" * 60)
