"""
Test script for LOESS functions
Validation of the Python wrapper for plain and robust LOESS.
"""

import numpy as np
import sys
import os
import ctypes

# Path configuration to import your functions
# Adjust the path if your module is in a different directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (
    compute_empirical_p_values
)




def _assert_allclose(a, b, tol=1e-12, msg=""):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if not np.allclose(a, b, rtol=0.0, atol=tol, equal_nan=True):
        raise AssertionError(msg or f"Arrays not close.\nExpected: {b}\nGot: {a}")

def _naive_empirical_p_values(distribution, c_const):
    """
    Naive reference implementation that matches your documented behavior:
    - distribution values < 0 are "invalid": p=1 for those genes
    - distribution is clamped to 0 for building the empirical distribution D
    - Upper-tail one-sided p-value with >=
    """
    dist = np.asarray(distribution, dtype=np.float64)
    n = dist.size
    if n == 0:
        return np.asarray([], dtype=np.float64)

    # Build D (sorted, clamped negatives to 0)
    D = dist.copy()
    D[D < 0.0] = 0.0
    D_sorted = np.sort(D)

    denom = n + float(c_const)
    out = np.empty(n, dtype=np.float64)

    for i, d in enumerate(dist):
        if d < 0.0:
            out[i] = 1.0
        else:
            count_ge = int(np.sum(D_sorted >= d))
            out[i] = (count_ge + float(c_const)) / denom
    return out

print("=== Testing empirical p-values Python wrapper ===")
print("Based on Fortran test suite with comprehensive test coverage")

# =====================
# Tests for compute_empirical_p_values
# =====================

def test_empirical_p_values_basic():
    print("\n[test_empirical_p_values_basic] Basic empirical p-values calculation test")
    distribution = np.array([0.5, 1.2, 0.8, 0.3], dtype=np.float64)
    c_const = 1.0

    p_values = compute_empirical_p_values(distribution, c_const)

    # Verify p-values are within [0, 1]
    assert np.all(p_values >= 0.0) and np.all(p_values <= 1.0), "p-values must be in [0,1]"

    # Verify against naive implementation
    expected = _naive_empirical_p_values(distribution, c_const)
    _assert_allclose(p_values, expected, tol=1e-12, msg="Basic test: p-values mismatch vs naive reference")

    print("Basic empirical p-values calculation test passed ✓")

def test_empirical_p_values_all_zeros():
    print("\n[test_empirical_p_values_all_zeros] All zeros distribution test")
    distribution = np.array([0, 0, 0, 0, 0], dtype=np.float64)
    c_const = 1.0

    p_values = compute_empirical_p_values(distribution, c_const)

    # All elements >= 0, and distribution D is all zeros:
    # count_ge = n for d=0 => (n+c)/(n+c) = 1
    expected = np.ones_like(distribution, dtype=np.float64)
    _assert_allclose(p_values, expected, tol=0.0, msg="All zeros: expected all ones")

    print("All zeros distribution test passed ✓")

def test_empirical_p_values_negative_values():
    print("\n[test_empirical_p_values_negative_values] Negative values in distribution test")
    distribution = np.array([-0.5, 1.2, -0.8, 0.3], dtype=np.float64)
    c_const = 1.0

    p_values = compute_empirical_p_values(distribution, c_const)

    # Verify p-values for negative values are 1
    assert np.all(p_values[distribution < 0.0] == 1.0), "Negative inputs must return p=1"

    # Verify against naive implementation
    expected = _naive_empirical_p_values(distribution, c_const)
    _assert_allclose(p_values, expected, tol=1e-12, msg="Negative values: p-values mismatch vs naive reference")

    print("Negative values in distribution test passed ✓")

def test_empirical_p_values_large_distribution():
    print("\n[test_empirical_p_values_large_distribution] Large distribution test")
    rng = np.random.default_rng(42)  # reproducible
    distribution = rng.uniform(0.0, 10.0, size=1000).astype(np.float64)
    c_const = 1.0

    p_values = compute_empirical_p_values(distribution, c_const)

    # Verify p-values are within [0, 1]
    assert np.all(p_values >= 0.0) and np.all(p_values <= 1.0), "p-values must be in [0,1]"

    # Spot-check a handful against naive (avoid O(n^2) full check for 1000 if you want faster tests)
    expected = _naive_empirical_p_values(distribution, c_const)
    idx = np.array([0, 1, 2, 10, 123, 999], dtype=int)
    _assert_allclose(p_values[idx], expected[idx], tol=1e-12, msg="Large dist: spot-check mismatch vs naive")

    print("Large distribution test passed ✓")

# Optional extra tests (recommended)

def test_empirical_p_values_monotonicity_on_sorted_inputs():
    print("\n[test_empirical_p_values_monotonicity_on_sorted_inputs] Monotonicity sanity check")
    distribution = np.array([0.0, 0.5, 1.0, 2.0, 4.0], dtype=np.float64)
    c_const = 1.0

    p_values = compute_empirical_p_values(distribution, c_const)

    # For increasing d, p(d) should be non-increasing (upper tail)
    assert np.all(p_values[:-1] >= p_values[1:]), "p-values must be non-increasing as d increases"

    print("Monotonicity test passed ✓")
def test_empirical_p_values_extremes():
    print("\n[test_empirical_p_values_extremes] Extremes check")
    rdi = np.array([-1.0, 0.0, 10.0, 3.0], dtype=np.float64)
    c_const = 1.0

    p_values = compute_empirical_p_values(rdi, c_const)

    n = rdi.size
    denom = n + c_const

    # negative -> 1
    assert p_values[0] == 1.0, "negative -> p=1"

    # d=0 -> all clamped values >=0 -> count=n -> 1
    assert abs(p_values[1] - 1.0) <= 0.0, "d=0 -> p=1"

    # d=10 is IN the distribution and is the max -> count=1 -> (1+c)/(n+c)
    expected = (1.0 + c_const) / denom
    assert abs(p_values[2] - expected) < 1e-12, "d==max (10) -> p=(1+c)/(n+c)"

    # d=3 -> values >=3 are [3,10] -> count=2 -> (2+c)/(n+c)
    expected = (2.0 + c_const) / denom
    assert abs(p_values[3] - expected) < 1e-12, "d=3 -> p=(2+c)/(n+c)"

    print("Extremes test passed ✓")


# =====================
# Run all tests
# =====================

print("\n=================================================")
print("    EMPIRICAL P VALUE TESTS (PYTHON)")
print("=================================================\n")

test_empirical_p_values_basic()
test_empirical_p_values_all_zeros()
test_empirical_p_values_negative_values()
test_empirical_p_values_large_distribution()

# optional extras
test_empirical_p_values_monotonicity_on_sorted_inputs()
test_empirical_p_values_extremes()

print("\n=================================================")
print("             ALL TESTS COMPLETED")
print("=================================================")
print("If you see this message, all empirical p-value Python interface tests passed! ✓")