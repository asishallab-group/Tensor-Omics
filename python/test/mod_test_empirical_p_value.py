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
from tensor_omics import compute_scaled_distance_quantile as _compute_scaled_distance_quantile
from test_helpers import run_all_tests


def compute_scaled_distance_quantile(distribution, c_const):
    """Sort-prep in front of the raw Fortran routine, which takes rdi/sorted_rdi/perm.

    Negatives are invalid and clamped to zero; perm is the 1-based ascending
    permutation of that clamped array.
    """
    dist = np.ascontiguousarray(distribution, dtype=np.float64)
    if dist.size == 0:
        return np.ascontiguousarray([], dtype=np.float64)

    sorted_rdi = dist.copy()
    sorted_rdi[sorted_rdi < 0.0] = 0.0
    perm = (np.argsort(sorted_rdi, kind="mergesort").astype(np.int32) + 1)

    return _compute_scaled_distance_quantile(dist, sorted_rdi, perm, float(c_const))


def _assert_allclose(a, b, tol=1e-12, msg=""):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if not np.allclose(a, b, rtol=0.0, atol=tol, equal_nan=True):
        raise AssertionError(msg or f"Arrays not close.\nExpected: {b}\nGot: {a}")


def _naive_scaled_distance_quantile(distribution, c_const):
    """
    Naive reference implementation that matches your documented behavior:
    - distribution values < 0 are "invalid": quantile=1 for those genes
    - distribution is clamped to 0 for building the empirical distribution D
    - Upper-tail one-sided quantile with >=
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


# =====================
# Tests for compute_scaled_distance_quantile
# =====================

def test_empirical_p_values_basic():
    distribution = np.array([0.5, 1.2, 0.8, 0.3], dtype=np.float64)
    c_const = 1.0

    quantile = compute_scaled_distance_quantile(distribution, c_const)

    # Verify quantile is within [0, 1]
    assert np.all(quantile >= 0.0) and np.all(quantile <= 1.0), "quantile must be in [0,1]"

    # Verify against naive implementation
    expected = _naive_scaled_distance_quantile(distribution, c_const)
    _assert_allclose(quantile, expected, tol=1e-12, msg="Basic test: quantile mismatch vs naive reference")


def test_empirical_p_values_all_zeros():
    distribution = np.array([0, 0, 0, 0, 0], dtype=np.float64)
    c_const = 1.0

    quantile = compute_scaled_distance_quantile(distribution, c_const)

    # All elements >= 0, and distribution D is all zeros:
    # count_ge = n for d=0 => (n+c)/(n+c) = 1
    expected = np.ones_like(distribution, dtype=np.float64)
    _assert_allclose(quantile, expected, tol=0.0, msg="All zeros: expected all ones")


def test_empirical_p_values_negative_values():
    distribution = np.array([-0.5, 1.2, -0.8, 0.3], dtype=np.float64)
    c_const = 1.0

    quantile = compute_scaled_distance_quantile(distribution, c_const)

    # Verify quantile for negative values are 1
    assert np.all(quantile[distribution < 0.0] == 1.0), "Negative inputs must return quantile=1"

    # Verify against naive implementation
    expected = _naive_scaled_distance_quantile(distribution, c_const)
    _assert_allclose(quantile, expected, tol=1e-12, msg="Negative values: quantile mismatch vs naive reference")


def test_empirical_p_values_large_distribution():
    rng = np.random.default_rng(42)  # reproducible
    distribution = rng.uniform(0.0, 10.0, size=1000).astype(np.float64)
    c_const = 1.0

    quantile = compute_scaled_distance_quantile(distribution, c_const)

    # Verify quantile is within [0, 1]
    assert np.all(quantile >= 0.0) and np.all(quantile <= 1.0), "quantile must be in [0,1]"

    # Spot-check a handful against naive (avoid O(n^2) full check for 1000 if you want faster tests)
    expected = _naive_scaled_distance_quantile(distribution, c_const)
    idx = np.array([0, 1, 2, 10, 123, 999], dtype=int)
    _assert_allclose(quantile[idx], expected[idx], tol=1e-12, msg="Large dist: spot-check mismatch vs naive")


# Optional extra tests (recommended)

def test_empirical_p_values_monotonicity_on_sorted_inputs():
    distribution = np.array([0.0, 0.5, 1.0, 2.0, 4.0], dtype=np.float64)
    c_const = 1.0

    quantile = compute_scaled_distance_quantile(distribution, c_const)

    # For increasing d, quantile(d) should be non-increasing (upper tail)
    assert np.all(quantile[:-1] >= quantile[1:]), "quantile must be non-increasing as d increases"


def test_empirical_p_values_extremes():
    rdi = np.array([-1.0, 0.0, 10.0, 3.0], dtype=np.float64)
    c_const = 1.0

    quantile = compute_scaled_distance_quantile(rdi, c_const)

    n = rdi.size
    denom = n + c_const

    # negative -> 1
    assert quantile[0] == 1.0, "negative -> quantile=1"

    # d=0 -> all clamped values >=0 -> count=n -> 1
    assert abs(quantile[1] - 1.0) <= 0.0, "d=0 -> quantile=1"

    # d=10 is IN the distribution and is the max -> count=1 -> (1+c)/(n+c)
    expected = (1.0 + c_const) / denom
    assert abs(quantile[2] - expected) < 1e-12, "d==max (10) -> quantile=(1+c)/(n+c)"

    # d=3 -> values >=3 are [3,10] -> count=2 -> (2+c)/(n+c)
    expected = (2.0 + c_const) / denom
    assert abs(quantile[3] - expected) < 1e-12, "d=3 -> quantile=(2+c)/(n+c)"


if __name__ == "__main__":
    run_all_tests(globals().values())
