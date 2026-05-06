
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
from test_helpers import run_all_tests, assert_error


# =====================
# Test 1: Simple EDF Test
# =====================
def test_edf_simple():
    """Test basic EDF computation with simple dataset"""
    values = np.array([1.0, 2.0, 2.0, 3.0, 3.0, 3.0])

    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']

    # Expected: 1.0 appears once (1/6), 2.0 appears twice (3/6), 3.0 appears thrice (6/6)
    expected_unique = [1.0, 2.0, 3.0]
    expected_cdf = [1.0/6.0, 3.0/6.0, 6.0/6.0]

    assert n_unique == 3, f"Expected 3 unique values, got {n_unique}"

    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}: expected {expected_unique[i]}, got {unique_vals[i]}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}: expected {expected_cdf[i]}, got {cdf_vals[i]}"


# =====================
# Test 2: All Unique Values
# =====================
def test_edf_all_unique():
    """Test EDF when all values are unique"""
    values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])

    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']

    # All values are unique, so CDF should be 0.2, 0.4, 0.6, 0.8, 1.0
    expected_unique = [1.0, 2.0, 3.0, 4.0, 5.0]
    expected_cdf = [0.2, 0.4, 0.6, 0.8, 1.0]

    assert n_unique == 5, f"Expected 5 unique values, got {n_unique}"

    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}"


# =====================
# Test 3: All Same Values
# =====================
def test_edf_all_same():
    """Test EDF when all values are identical"""
    values = np.array([5.0, 5.0, 5.0, 5.0])

    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']

    # Only one unique value with CDF = 1.0
    assert n_unique == 1, f"Expected 1 unique value, got {n_unique}"
    assert abs(unique_vals[0] - 5.0) < 1e-12, "Unique value should be 5.0"
    assert abs(cdf_vals[0] - 1.0) < 1e-12, "CDF should be 1.0"


# =====================
# Test 4: Duplicates
# =====================
def test_edf_duplicates():
    """Test EDF with various duplicate values"""
    values = np.array([1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 4.0])

    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']

    # Expected: 1.0 (2/7), 2.0 (3/7), 3.0 (6/7), 4.0 (7/7)
    expected_unique = [1.0, 2.0, 3.0, 4.0]
    expected_cdf = [2.0/7.0, 3.0/7.0, 6.0/7.0, 7.0/7.0]

    assert n_unique == 4, f"Expected 4 unique values, got {n_unique}"

    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}"


# =====================
# Test 5: Single Value
# =====================
def test_edf_single_value():
    """Test EDF with a single value"""
    values = np.array([42.0])

    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']

    assert n_unique == 1, f"Expected 1 unique value, got {n_unique}"
    assert abs(unique_vals[0] - 42.0) < 1e-12, "Unique value should be 42.0"
    assert abs(cdf_vals[0] - 1.0) < 1e-12, "CDF should be 1.0"


# =====================
# Test 6: Empty Input
# =====================
def test_edf_empty_input():
    """Test EDF with empty array (should fail with error 202)"""
    values = np.array([])

    assert_error(lambda: compute_edf(values), "Expected error for empty input")

# =====================
# Test 7: Large Dataset
# =====================
def test_edf_large_dataset():
    """Test EDF with a larger dataset"""
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

    # Expected: 1.0 (0.25), 2.0 (0.5), 3.0 (0.75), 4.0 (1.0)
    expected_unique = [1.0, 2.0, 3.0, 4.0]
    expected_cdf = [0.25, 0.5, 0.75, 1.0]

    assert n_unique == 4, f"Expected 4 unique values, got {n_unique}"

    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}"


# =====================
# Test 8: Negative Values
# =====================
def test_edf_negative_values():
    """Test EDF with negative values"""
    values = np.array([-3.0, -1.0, 0.0, 1.0, 3.0])

    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']

    # All values are unique
    expected_unique = [-3.0, -1.0, 0.0, 1.0, 3.0]
    expected_cdf = [0.2, 0.4, 0.6, 0.8, 1.0]

    assert n_unique == 5, f"Expected 5 unique values, got {n_unique}"

    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}"


# =====================
# Test 9: Unsorted Input
# =====================
def test_edf_unsorted_input():
    """Test EDF with unsorted input (function should handle internally)"""
    values = np.array([3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0])

    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']

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


# =====================
# Test 10: Python List Input
# =====================
def test_edf_list_input():
    """Test EDF with Python list instead of numpy array"""
    values = [1.0, 2.0, 2.0, 3.0]

    result = compute_edf(values)
    unique_vals = result['unique_values']
    cdf_vals = result['cdf_values']
    n_unique = result['n_unique']

    expected_unique = [1.0, 2.0, 3.0]
    expected_cdf = [0.25, 0.75, 1.0]

    assert n_unique == 3, f"Expected 3 unique values, got {n_unique}"

    for i in range(n_unique):
        assert abs(unique_vals[i] - expected_unique[i]) < 1e-12, \
            f"Unique value mismatch at index {i}"
        assert abs(cdf_vals[i] - expected_cdf[i]) < 1e-12, \
            f"CDF value mismatch at index {i}"


if __name__ == "__main__":
    run_all_tests(globals().values())
