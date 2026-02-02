"""
Test script for clock hand angle functions
Python equivalent of the R and Fortran clock hand angle tests
"""

import numpy as np
import sys
import os
from math import pi as PI

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (
    tox_compute_contributions,
    tox_compute_all_contributions,
    tox_compute_baselines_factor_dependent,
    tox_perform_permutation_test,
    tox_compute_p_values
)

# Constants
TOL = 1e-12


def test_tox_compute_baselines_factor_dependent():
    """Test baseline computation wrapper across all supported modes and error cases."""

    factor = np.array([1.0, 3.0, 2.0, 4.0], dtype=np.float64)
    dependent = np.array([5.0, 7.0, 6.0, 8.0], dtype=np.float64)

    # RAW mode => zero baselines
    res_raw = tox_compute_baselines_factor_dependent(factor, dependent, mode="raw")
    assert np.isclose(res_raw['baseline_factor'], 0.0, atol=TOL)
    assert np.isclose(res_raw['baseline_dependent'], 0.0, atol=TOL)

    # MIN mode => min values
    res_min = tox_compute_baselines_factor_dependent(factor, dependent, mode="min")
    assert np.isclose(res_min['baseline_factor'], np.min(factor), atol=TOL)
    assert np.isclose(res_min['baseline_dependent'], np.min(dependent), atol=TOL)

    # MEAN mode => arithmetic mean
    res_mean = tox_compute_baselines_factor_dependent(factor, dependent, mode="mean")
    assert np.isclose(res_mean['baseline_factor'], np.mean(factor), atol=TOL)
    assert np.isclose(res_mean['baseline_dependent'], np.mean(dependent), atol=TOL)

    # Mismatched lengths should raise ValueError
    try:
        tox_compute_baselines_factor_dependent(factor, dependent[:-1], mode=1)
        raise AssertionError("Expected ValueError for mismatched lengths")
    except ValueError:
        pass

    # Invalid mode should bubble up as RuntimeError from Fortran layer
    try:
        tox_compute_baselines_factor_dependent(factor, dependent, mode="unknown_mode")
        raise AssertionError("Expected RuntimeError for invalid mode")
    except RuntimeError:
        pass

    print("✅ tox_compute_baselines_factor_dependent passed all tests.")


def test_compute_contributions():
    # -------------------------------
    # Case 1: RAW baseline
    # -------------------------------
    factor = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    dependent = np.array([2.0, 1.0, 0.0, -1.0], dtype=np.float64)
    local, total = tox_compute_contributions(factor, dependent, mode="raw").values()

    expected_local = factor * dependent
    expected_total = sum(expected_local)
    assert np.allclose(local, expected_local, atol=TOL), "Case 1 local contributions mismatch"
    assert abs(total - expected_total) < TOL, "Case 1 total contribution mismatch"

    # -------------------------------
    # Case 2: MIN baseline
    # -------------------------------
    factor = np.array([3.0, 5.0, 2.0, 4.0], dtype=np.float64)
    dependent = np.array([1.0, 2.0, 0.0, -1.0], dtype=np.float64)
    local, total = tox_compute_contributions(factor, dependent, mode="min").values()

    expected_local = (factor - factor.min()) * (dependent - dependent.min())
    expected_total = sum(expected_local)
    assert np.allclose(local, expected_local, atol=TOL), "Case 2 local contributions mismatch"
    assert abs(total - expected_total) < TOL, "Case 2 total contribution mismatch"

    # -------------------------------
    # Case 3: MEAN baseline
    # -------------------------------
    factor = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    dependent = np.array([4.0, 3.0, 2.0, 1.0], dtype=np.float64)
    local, total = tox_compute_contributions(factor, dependent, mode="mean").values()

    expected_local = (factor - factor.mean()) * (dependent - dependent.mean())
    expected_total = sum(expected_local)
    assert np.allclose(local, expected_local, atol=TOL), "Case 3 local contributions mismatch"
    assert abs(total - expected_total) < TOL, "Case 3 total contribution mismatch"

    print("✅ Compute contributions test passed.")


def test_compute_all_contributions():
    # -------------------------------
    # Case 1: MEAN baseline
    # -------------------------------
    # Factor trajectory: [1,2,3]
    # Dependent trajectory: [4,5,6]
    trajectories = np.empty((2, 1, 3), dtype=np.float64, order="F")
    trajectories[0, 0, :] = [1.0, 2.0, 3.0]
    trajectories[1, 0, :] = [4.0, 5.0, 6.0]

    factor_indices = np.array([1], dtype=np.int32, order="F")
    dependent_indices = np.array([2], dtype=np.int32, order="F")

    local, total = tox_compute_all_contributions(trajectories, factor_indices, dependent_indices, mode="mean").values()

    # Baselines: mean(factor)=2.0, mean(dependent)=5.0
    expected_local = np.array([1.0, 0.0, 1.0], dtype=np.float64, order="F")
    expected_total = 2.0

    assert np.allclose(local[:, 0, 0, 0], expected_local, atol=TOL), "Case 1 local contributions mismatch"
    assert abs(total[0, 0, 0] - expected_total) < TOL, "Case 1 total contribution mismatch"

    # -------------------------------
    # Case 2: MIN baseline
    # -------------------------------
    # Factor trajectory: [2,4,6]
    # Dependent trajectory: [1,3,5]
    trajectories[0, 0, :] = [2.0, 4.0, 6.0]
    trajectories[1, 0, :] = [5.0, 3.0, 5.0]

    factor_indices = np.array([1], dtype=np.int32, order="F")
    dependent_indices = np.array([2], dtype=np.int32, order="F")

    local, total = tox_compute_all_contributions(trajectories, factor_indices, dependent_indices, mode="min").values()

    # Baselines: min(factor)=2.0, min(dependent)=1.0
    expected_local = np.array([0.0, 0.0, 8.0], dtype=np.float64, order="F")
    expected_total = 8.0

    assert np.allclose(local[:, 0, 0, 0], expected_local, atol=TOL), "Case 2 local contributions mismatch"
    assert abs(total[0, 0, 0] - expected_total) < TOL, "Case 2 total contribution mismatch"

    print("✅ Compute all contributions test passed.")


def test_perform_permutation_test():
    # -------------------------------
    # Setup synthetic trajectories
    # -------------------------------
    n_factors, n_samples, n_timepoints = 2, 3, 3
    n_permutations = 3

    trajectories = np.empty((n_factors, n_samples, n_timepoints), dtype=np.float64, order="F")

    # Factor 1 values across samples/timepoints
    trajectories[0,0,:] = [1.0, 2.0, 3.0]
    trajectories[0,1,:] = [2.0, 4.0, 6.0]
    trajectories[0,2,:] = [3.0, 6.0, 9.0]

    # Dependent 2 values across samples/timepoints
    trajectories[1,0,:] = [4.0, 5.0, 6.0]
    trajectories[1,1,:] = [1.0, 3.0, 5.0]
    trajectories[1,2,:] = [2.0, 4.0, 6.0]

    factor_idx    = 1   # Fortran-style 1-based
    dependent_idx = 2
    sample_idx    = 1
    mode          = "mean"   # MEAN baseline

    # -------------------------------
    # Call wrapper with fixed seed
    # -------------------------------
    result = tox_perform_permutation_test(
        trajectories,
        factor_idx=factor_idx,
        dependent_idx=dependent_idx,
        sample_idx=sample_idx,
        mode=mode,
        n_permutations=n_permutations,
        random_seed=12345
    )

    local, total = result.values()

    # -------------------------------
    # Expected values for one known permutation
    # -------------------------------
    # Factor trajectory (sample 1): [1,2,3], mean=2.0
    # Suppose RNG picks dependent sample 2: [1,3,5], mean=3.0
    expected_local = np.array([
        (1.0-2.0)*(1.0-3.0),   # 2.0
        (2.0-2.0)*(3.0-3.0),   # 0.0
        (3.0-2.0)*(5.0-3.0)    # 2.0
    ], dtype=np.float64, order="F")
    expected_total = expected_local.sum()  # 4.0

    # -------------------------------
    # Assertions
    # -------------------------------
    # With fixed seed, first permutation is reproducible → check against sample 2
    assert np.allclose(local[:,0], expected_local, atol=TOL), "Permutation 1 local contributions mismatch"
    assert abs(total[0] - expected_total) < TOL, "Permutation 1 total contribution mismatch"

    # For permutations 2 and 3, RNG may pick sample 2 or 3.
    # We assert that totals are finite and contributions are not all zero.
    assert np.all(np.isfinite(total)), "Permutation totals must be finite"
    assert not np.allclose(local, 0.0, atol=TOL), "Local contributions should not all be zero"

    print("✅ Permutation test test passed.")


def test_compute_p_values():
    # -------------------------------
    # Case 1: Valid inputs
    # -------------------------------
    n_timepoints, n_permutations = 3, 4
    local_obs = np.array([2.0, 0.0, 2.0], dtype=np.float64, order="F")
    total_obs = 4.0

    local_perm = np.empty((n_timepoints, n_permutations), dtype=np.float64, order="F")
    total_perm = np.empty(n_permutations, dtype=np.float64, order="F")

    local_perm[:,0] = [1.0, 0.0, 1.0]; total_perm[0] = 2.0
    local_perm[:,1] = [2.0, 0.0, 2.0]; total_perm[1] = 4.0
    local_perm[:,2] = [3.0, 1.0, 3.0]; total_perm[2] = 7.0
    local_perm[:,3] = [0.0, 0.0, 0.0]; total_perm[3] = 0.0

    result = tox_compute_p_values(local_obs, total_obs, local_perm, total_perm)
    local_p = result["local_p_values"]
    total_p = result["total_p_value"]

    expected_local = np.array([0.5, 1.0, 0.5], dtype=np.float64)
    expected_total = 0.5

    assert np.allclose(local_p, expected_local, atol=TOL), "Valid local p-values mismatch"
    assert abs(total_p - expected_total) < TOL, "Valid total p-value mismatch"

    # -------------------------------
    # Case 2: NaN in observed contributions
    # -------------------------------
    local_obs_nan = np.array([2.0, 0.0, np.nan], dtype=np.float64, order="F")
    try:
        tox_compute_p_values(local_obs_nan, total_obs, local_perm, total_perm)
        raise AssertionError("Expected RuntimeError for NaN input")
    except RuntimeError:
        pass  # Expected error

    # -------------------------------
    # Case 3: Inf in permutation contributions
    # -------------------------------
    local_perm_inf = local_perm.copy(order="F")
    local_perm_inf[2,3] = np.inf
    try:
        tox_compute_p_values(local_obs, total_obs, local_perm_inf, total_perm)
        raise AssertionError("Expected RuntimeError for Infinity input")
    except RuntimeError:
        pass  # Expected error

    print("✅ Compute p values test passed.")


def main():
    print("=================================================")
    print("    TRAJECTORY CONTRIBUTION ANALYSIS PYTHON INTERFACE TESTS")
    print("=================================================")
    print()

    test_tox_compute_baselines_factor_dependent()
    test_compute_contributions()
    test_compute_all_contributions()
    test_perform_permutation_test()
    test_compute_p_values()


if __name__ == "__main__":
    main()
