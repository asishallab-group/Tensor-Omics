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

    # The robust scratch arrays (robust/combined weights, residuals, permutation) are now
    # allocated by the wrapper, so the friendly binding no longer takes them.
    yhat = loess_fit_robust(
        x=x, y=y, weights=w, eval_points=z.reshape(n, 1),
        span=0.5, degree=1, max_neighborhood_size=n,
        compute_influence=False, save_factorization=False, n_iters=4
    )

    assert yhat.shape == (n,), "Output shape mismatch"
    # If robustness (bisquare reweighting) works, the outlier at index 5 
    # should be largely ignored, resulting in a value much lower than 100.
    assert yhat[5] < 50.0, f"Robust LOESS failed to suppress outlier: got {yhat[5]}"


def test_plain_and_robust_are_separate_entry_points():
    """The plain and robust fits are separate entry points, and disagree on noisy data."""
    n = 30
    x = np.arange(n, dtype=np.float64)
    y = np.sin(x / 5.0)
    y[7] += 5.0  # an outlier for the robust iterations to down-weight

    # the self-allocating entry points still take what the fit is *of*: the weights and the
    # points to evaluate at. Uniform weights, evaluated at the training x, is the common case.
    weights = np.ones(n)
    eval_points = x.reshape(n, 1)

    yhat_plain = loess_fit_plain(
        x=x, y=y, weights=weights, eval_points=eval_points,
        span=0.4, degree=1, max_neighborhood_size=n,
    )
    assert yhat_plain.shape == (n,), "Plain shape mismatch"

    yhat_robust = loess_fit_robust(
        x=x, y=y, weights=weights, eval_points=eval_points,
        span=0.4, degree=1, max_neighborhood_size=n, n_iters=2,
    )
    assert yhat_robust.shape == (n,), "Robust shape mismatch"

    # the robust fit suppresses the outlier, so the two cannot agree
    assert not np.array_equal(yhat_plain, yhat_robust), "Plain and Robust results should not be identical"


if __name__ == "__main__":
    run_all_tests(globals().values())
