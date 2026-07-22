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
from tensor_omics import (
    tox_loess_required_workspace,
    loess_fit_plain,
    loess_fit_robust,
    loess
)
from test_helpers import run_all_tests


def test_workspace_calculation():
    """Validates that the workspace recommendation function returns consistent values."""

    # Typical parameters: univariate, 100 neighbours, factorization saved
    ws = tox_loess_required_workspace(n_dim=1, max_neighborhood_size=100, save_factorization=True)

    assert isinstance(ws, dict), "Should return a dictionary"
    assert ws["int_workspace_size"] > 0, "Integer workspace size should be positive"
    assert ws["real_workspace_size"] > 0, "Real workspace size should be positive"


def test_loess_plain_functionality():
    """Tests the low-level loess_fit_plain subroutine."""
    n = 20
    x = np.linspace(1, 10, n)
    # Generate a linear trend with minor Gaussian noise
    y = 2.0 * x + np.random.normal(0, 0.1, n) 
    w = np.ones(n)
    z = x.copy()

    yhat = loess_fit_plain(
        x=x, y=y, weights=w, eval_points=z.reshape(n, 1),
        span=0.5, degree=1, max_neighborhood_size=n,
        compute_influence=False, save_factorization=False
    )

    assert yhat.shape == (n,), "Output shape mismatch"
    assert not np.any(np.isnan(yhat)), "Output contains NaNs"


def test_loess_robust_functionality():
    """Tests the low-level loess_fit_robust subroutine with outlier suppression."""
    n = 20
    x = np.linspace(1, 10, n)
    y = 3.0 * x
    y[5] = 100.0  # Introduce an aggressive outlier

    w = np.ones(n)
    z = x.copy()
    # Additional arrays required specifically for the robust version
    rw = np.zeros(n, dtype=np.float64)
    ww = np.zeros(n, dtype=np.float64)
    res = np.zeros(n, dtype=np.float64)
    pi = np.zeros(n, dtype=np.int32)

    yhat = loess_fit_robust(
        x=x, y=y, weights=w, eval_points=z.reshape(n, 1),
        span=0.5, degree=1, max_neighborhood_size=n,
        compute_influence=False, save_factorization=False, n_iters=4,
        robust_weights=rw, combined_weights=ww,
        residuals=res, permutation_indices=pi
    )

    assert yhat.shape == (n,), "Output shape mismatch"
    # If robustness (bisquare reweighting) works, the outlier at index 5 
    # should be largely ignored, resulting in a value much lower than 100.
    assert yhat[5] < 50.0, f"Robust LOESS failed to suppress outlier: got {yhat[5]}"


def test_tox_loess_wrapper():
    """Tests the high-level wrapper that selects between plain and robust modes."""
    n = 30
    x = np.arange(n, dtype=np.float64)
    y = np.sin(x / 5.0)

    # Test Plain mode (mode="plain", n_iters=0)
    yhat_plain = loess(x, y, span=0.4, degree=1, mode="plain", n_iters=0)
    assert yhat_plain.shape == (n,), "Plain mode shape mismatch"

    # Test Robust mode (mode=1)
    yhat_robust = loess(x, y, span=0.4, degree=1, mode="robust", n_iters=2)
    assert yhat_robust.shape == (n,), "Robust mode shape mismatch"

    # Verify that results differ due to robust iterations
    assert not np.array_equal(yhat_plain, yhat_robust), "Plain and Robust results should not be identical"


if __name__ == "__main__":
    run_all_tests(globals().values())
