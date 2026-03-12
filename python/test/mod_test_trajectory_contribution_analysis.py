"""
Test script for trajectory contribution analysis functions.
Python equivalent of the R and Fortran trajectory contribution tests.
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
    tox_compute_p_values,
    tox_compute_velocity_trajectory,
    tox_compute_acceleration_from_velocity_trajectory,
    tox_compute_velocity_trajectories,
    tox_compute_acceleration_from_velocity,
    tox_compute_velocity_acceleration_contributions,
    tox_compute_velocity_acceleration_contributions_expert
)


# Constants
TOL = 1e-12


def _expected_velocity(trajectories: np.ndarray) -> np.ndarray:
    trajectories = np.asarray(trajectories, dtype=np.float64)
    if trajectories.ndim != 3:
        raise ValueError("trajectories must be 3D (n_factors, n_samples, n_timepoints)")
    velocity = np.zeros_like(trajectories)
    if trajectories.shape[2] > 1:
        velocity[:, :, 1:] = trajectories[:, :, 1:] - trajectories[:, :, :-1]
    return velocity


def _expected_acceleration(velocity: np.ndarray) -> np.ndarray:
    velocity = np.asarray(velocity, dtype=np.float64)
    if velocity.ndim != 3:
        raise ValueError("velocity must be 3D (n_factors, n_samples, n_timepoints)")
    acceleration = np.zeros_like(velocity)
    if velocity.shape[2] > 2:
        acceleration[:, :, 2:] = velocity[:, :, 2:] - velocity[:, :, 1:-1]
    return acceleration


def test_tox_compute_velocity_trajectories():
    """Test velocity computation wrapper."""
    # trajectories shape: (n_factors=1, n_samples=2, n_timepoints=4)
    trajectories = np.array(
        [[[1.0, 2.0, 4.0, 7.0],
          [0.0, -1.0, -1.0,  0.0]]],
        dtype=np.float64,
    )

    print("Trajectories shape:", trajectories.shape)
    print("Trajectories:\n", trajectories)

    velocity = tox_compute_velocity_trajectories(trajectories)
    expected = _expected_velocity(trajectories)

    print("Computed velocity:\n", velocity)
    print("Expected velocity:\n", expected)
    print("Difference:\n", velocity - expected)

    assert velocity.shape == trajectories.shape
    assert np.allclose(velocity, expected, atol=TOL), f"Velocity mismatch! Max diff: {np.max(np.abs(velocity - expected))}"

    print("✅ tox_compute_velocity_trajectories passed.")


def test_tox_compute_acceleration_from_velocity():
    """Test acceleration computation wrapper."""
    # velocity shape: (n_factors=1, n_samples=2, n_timepoints=4)
    velocity = np.array(
        [[[0.0, 1.0, 2.0, 3.0],
          [0.0, -1.0, 0.0, 1.0]]],
        dtype=np.float64,
    )

    acceleration = tox_compute_acceleration_from_velocity(velocity)
    expected = _expected_acceleration(velocity)

    assert acceleration.shape == velocity.shape
    assert np.allclose(acceleration, expected, atol=TOL)

    print("✅ tox_compute_acceleration_from_velocity passed.")


def test_tox_compute_velocity_acceleration_contributions():
    """Test velocity and acceleration contribution computation wrapper."""
    # trajectories shape: (n_factors=2, n_samples=1, n_timepoints=4)
    trajectories = np.array(
        [[[1.0, 3.0, 6.0, 10.0]],   # Factor 1, Sample 1
         [[1.0, 2.0, 2.0, 1.0]]],   # Factor 2, Sample 1
        dtype=np.float64,
    )
    mode = "raw"

    print(f"Input trajectories shape: {trajectories.shape}")
    assert trajectories.shape == (2, 1, 4), f"Expected (2, 1, 4), got {trajectories.shape}"

    result = tox_compute_velocity_acceleration_contributions(trajectories, mode)

    C_vel = result["C_velocity"]
    series_vel = result["velocity_contribution_series"]
    C_acc = result["C_acceleration"]
    series_acc = result["acceleration_contribution_series"]

    # Output shapes: (n_samples, n_factors, n_factors, ...)
    assert C_vel.shape == (1, 2, 2), f"Expected (1, 2, 2), got {C_vel.shape}"
    assert series_vel.shape == (1, 2, 2, 4), f"Expected (1, 2, 2, 4), got {series_vel.shape}"
    assert C_acc.shape == (1, 2, 2), f"Expected (1, 2, 2), got {C_acc.shape}"
    assert series_acc.shape == (1, 2, 2, 4), f"Expected (1, 2, 2, 4), got {series_acc.shape}"

    print("\n=== DEBUG OUTPUT ===")
    print(f"C_vel:\n{C_vel}")
    print(f"series_vel shape: {series_vel.shape}")
    print(f"series_vel[0,0,1,:]:\n{series_vel[0,0,1,:]}")
    print(f"C_acc:\n{C_acc}")
    print(f"series_acc[0,0,1,:]:\n{series_acc[0,0,1,:]}")

    expected_velocity = _expected_velocity(trajectories)
    expected_acceleration = _expected_acceleration(expected_velocity)

    print(f"\nExpected velocity:\n{expected_velocity}")
    print(f"Expected acceleration:\n{expected_acceleration}")

    # Extract velocities: factor 0 and factor 1 for sample 0
    factor_velocity = expected_velocity[0, 0, 1:]      # Factor 0, timepoints 2-4
    dependent_velocity = expected_velocity[1, 0, 1:]   # Factor 1, timepoints 2-4
    factor_acceleration = expected_acceleration[0, 0, 2:]   # Factor 0, timepoints 3-4
    dependent_acceleration = expected_acceleration[1, 0, 2:]  # Factor 1, timepoints 3-4

    print(f"\nFactor velocity (t=2..4): {factor_velocity}")
    print(f"Dependent velocity (t=2..4): {dependent_velocity}")
    print(f"Factor acceleration (t=3..4): {factor_acceleration}")
    print(f"Dependent acceleration (t=3..4): {dependent_acceleration}")

    expected_series_vel = np.zeros(4)
    if factor_velocity.size > 0:
        vel_contribs = factor_velocity * dependent_velocity
        expected_total_vel = vel_contribs.sum()
        expected_series_vel[1:] = vel_contribs
        print(f"\nVelocity contributions: {vel_contribs}")
        print(f"Expected total velocity: {expected_total_vel}")
    else:
        expected_total_vel = 0.0

    expected_series_acc = np.zeros(4)
    if factor_acceleration.size > 0:
        acc_contribs = factor_acceleration * dependent_acceleration
        expected_total_acc = acc_contribs.sum()
        expected_series_acc[2:] = acc_contribs
        print(f"\nAcceleration contributions: {acc_contribs}")
        print(f"Expected total acceleration: {expected_total_acc}")
    else:
        expected_total_acc = 0.0

    print(f"\nAssertion check: C_vel[0, 0, 1] = {C_vel[0, 0, 1]} vs expected {expected_total_vel}")

    # Check sample 0, factor 0 → factor 1 contributions
    assert np.isclose(C_vel[0, 0, 1], expected_total_vel, atol=TOL), \
        f"C_vel[0,0,1] = {C_vel[0, 0, 1]}, expected {expected_total_vel}"
    assert np.allclose(series_vel[0, 0, 1, :], expected_series_vel, atol=TOL), \
        f"series_vel mismatch"
    assert np.isclose(C_acc[0, 0, 1], expected_total_acc, atol=TOL), \
        f"C_acc[0,0,1] = {C_acc[0, 0, 1]}, expected {expected_total_acc}"
    assert np.allclose(series_acc[0, 0, 1, :], expected_series_acc, atol=TOL), \
        f"series_acc mismatch"

    print("✅ tox_compute_velocity_acceleration_contributions passed.")


def test_tox_compute_velocity_acceleration_contributions_expert():
    trajectories = np.array([
        [[1.0, 2.0, 3.0, 4.0],
         [2.0, 4.0, 6.0, 8.0]]
    ], dtype=np.float64)

    result_base = tox_compute_velocity_acceleration_contributions(trajectories, "raw")
    result_expert = tox_compute_velocity_acceleration_contributions_expert(trajectories, "raw")

    for key in result_base:
        assert key in result_expert, f"Missing key {key} in expert result"
        assert result_base[key].shape == result_expert[key].shape
        assert np.allclose(result_base[key], result_expert[key], atol=TOL)

    print("✅ tox_compute_velocity_acceleration_contributions_expert passed.")


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


def test_tox_compute_velocity_trajectory():
    """Test single-trajectory velocity computation wrapper."""
    trajectory = np.array([1.0, 2.0, 4.0, 7.0], dtype=np.float64)
    velocity = tox_compute_velocity_trajectory(trajectory)

    expected = np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float64)

    assert velocity.shape == trajectory.shape
    assert np.allclose(velocity, expected, atol=TOL)

    # Dimensionality check
    try:
        tox_compute_velocity_trajectory(trajectory.reshape(1, -1))
        assert False, "Expected ValueError for 2D input"
    except ValueError:
        pass  # Expected

    print("✅ tox_compute_velocity_trajectory passed.")


def test_tox_compute_acceleration_from_velocity_trajectory():
    """Test single-trajectory acceleration computation wrapper."""
    velocity = np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float64)
    acceleration = tox_compute_acceleration_from_velocity_trajectory(velocity)

    expected = np.array([0.0, 0.0, 1.0, 1.0], dtype=np.float64)

    assert acceleration.shape == velocity.shape
    assert np.allclose(acceleration, expected, atol=TOL)

    # Dimensionality check
    try:
        tox_compute_acceleration_from_velocity_trajectory(velocity.reshape(1, -1))
        assert False, "Expected ValueError for 2D input"
    except ValueError:
        pass  # Expected

    print("✅ tox_compute_acceleration_from_velocity_trajectory passed.")


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
    test_tox_compute_velocity_trajectories()
    test_tox_compute_acceleration_from_velocity()
    test_tox_compute_velocity_acceleration_contributions()
    test_tox_compute_velocity_acceleration_contributions_expert()
    test_tox_compute_velocity_trajectory()
    test_tox_compute_acceleration_from_velocity_trajectory()

    


if __name__ == "__main__":
    main()
