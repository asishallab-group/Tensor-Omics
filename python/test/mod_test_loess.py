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
    tox_loess_required_workspace,
    loess_fit_plain,
    loess_fit_robust,
    tox_loess
)

def test_workspace_calculation():
    """Validates that the workspace recommendation function returns consistent values."""
    print("Testing tox_loess_required_workspace...")

    # Typical parameters: d=1 (univariate), nvmax=100, setlf=True
    ws = tox_loess_required_workspace(d=1, nvmax=100, setlf=True)

    assert isinstance(ws, dict), "Should return a dictionary"
    assert ws["liv"] > 0, "LIV (Integer Workspace) should be positive"
    assert ws["lv"] > 0, "LV (Real Workspace) should be positive"
    print(f" Workspace calculation passed (liv={ws['liv']}, lv={ws['lv']})")


def test_loess_plain_functionality():
    """Tests the low-level loess_fit_plain subroutine."""
    print("Testing loess_fit_plain...")
    n = 20
    x = np.linspace(1, 10, n)
    # Generate a linear trend with minor Gaussian noise
    y = 2.0 * x + np.random.normal(0, 0.1, n) 
    w = np.ones(n)
    z = x.copy()

    # Get required workspace sizes
    ws = tox_loess_required_workspace(d=1, nvmax=n, setlf=False)
    iv = np.zeros(ws["liv"], dtype=np.int32)
    wv = np.zeros(ws["lv"], dtype=np.float64)
    diagl = np.zeros(n, dtype=np.float64)

    print(" Workspace sizes - liv:", ws["liv"], "lv:", ws["lv"])

    yhat = loess_fit_plain(
        n=n, x=x, y=y, w=w, z=z, 
        span=0.5, degree=1, nvmax=n, 
        infl=False, setlf=False, 
        iv=iv, liv=ws["liv"], wv=wv, lv=ws["lv"], 
        diagl=diagl
    )

    assert yhat.shape == (n,), "Output shape mismatch"
    assert not np.any(np.isnan(yhat)), "Output contains NaNs"
    print(" loess_fit_plain passed.")


def test_loess_robust_functionality():
    """Tests the low-level loess_fit_robust subroutine with outlier suppression."""
    print("Testing loess_fit_robust...")
    n = 20
    x = np.linspace(1, 10, n)
    y = 3.0 * x
    y[5] = 100.0  # Introduce an aggressive outlier

    w = np.ones(n)
    z = x.copy()
    ws = tox_loess_required_workspace(d=1, nvmax=n, setlf=False)

    # Additional arrays required specifically for the robust version
    iv = np.zeros(ws["liv"], dtype=np.int32)
    wv = np.zeros(ws["lv"], dtype=np.float64)
    diagl = np.zeros(n, dtype=np.float64)
    rw = np.zeros(n, dtype=np.float64)
    ww = np.zeros(n, dtype=np.float64)
    res = np.zeros(n, dtype=np.float64)
    pi = np.zeros(n, dtype=np.int32)

    yhat = loess_fit_robust(
        n=n, x=x, y=y, w=w, z=z, 
        span=0.5, degree=1, nvmax=n, 
        infl=False, setlf=False, n_iters=4,
        iv=iv, liv=ws["liv"], wv=wv, lv=ws["lv"],
        diagl=diagl, rw=rw, ww=ww, res=res, pi=pi
    )

    assert yhat.shape == (n,), "Output shape mismatch"
    # If robustness (bisquare reweighting) works, the outlier at index 5 
    # should be largely ignored, resulting in a value much lower than 100.
    assert yhat[5] < 50.0, f"Robust LOESS failed to suppress outlier: got {yhat[5]}"
    print(" loess_fit_robust passed.")


def test_tox_loess_wrapper():
    """Tests the high-level wrapper that selects between plain and robust modes."""
    print("Testing tox_loess (High-level wrapper)...")
    n = 30
    x = np.arange(n, dtype=np.float64)
    y = np.sin(x / 5.0)

    # Test Plain mode (mode=0)
    yhat_plain = tox_loess(x, y, span=0.4, degree=1, mode=0)
    assert yhat_plain.shape == (n,), "Plain mode shape mismatch"

    # Test Robust mode (mode=1)
    yhat_robust = tox_loess(x, y, span=0.4, degree=1, mode=1, n_iters=2)
    assert yhat_robust.shape == (n,), "Robust mode shape mismatch"

    # Verify that results differ due to robust iterations
    assert not np.array_equal(yhat_plain, yhat_robust), "Plain and Robust results should not be identical"

    print(" tox_loess wrapper passed.")


def main():
    print("=================================================")
    print("      LOESS INTERFACE FUNCTIONAL TESTS")
    print("=================================================")
    print()

    try:
        test_workspace_calculation()
        test_loess_plain_functionality()
        test_loess_robust_functionality()
        test_tox_loess_wrapper()
        print("\n All LOESS integration tests passed successfully!")
    except Exception as e:
        print(f"\n Test failed with error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()