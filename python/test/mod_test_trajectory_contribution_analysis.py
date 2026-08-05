"""
Test script for trajectory contribution analysis functions.
Python equivalent of the R and Fortran trajectory contribution tests.
"""

import numpy as np
import sys
import os
from math import pi as PI

# Add parent directory to path to import tensor_omics
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    compute_contributions,
    compute_all_contributions,
    compute_baselines_factor_dependent,
    perform_permutation_test,
    compute_p_values,
    compute_velocity_trajectory,
    compute_acceleration_from_velocity_trajectory,
    compute_velocity_trajectories,
    compute_acceleration_from_velocity,
    compute_velocity_acceleration_contributions,
    compute_velocity_acceleration_contributions_expert
)
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_INVALID_INPUT, ERR_NAN_INF


# Constants
TOL = 1e-12


def _expected_velocity(trajectories: np.ndarray) -> np.ndarray:
    """Successive differences along time, returned time-major and one step shorter.

    Positions come in as (n_factors, n_samples, n_timepoints); velocity comes back as
    (n_timepoints - 1, n_factors, n_samples).
    """
    trajectories = np.asarray(trajectories, dtype=np.float64)
    if trajectories.ndim != 3:
        raise ValueError("trajectories must be 3D (n_factors, n_samples, n_timepoints)")
    return np.moveaxis(np.diff(trajectories, axis=2), 2, 0)


def _expected_acceleration(velocity: np.ndarray) -> np.ndarray:
    """Successive differences of an already time-major velocity."""
    velocity = np.asarray(velocity, dtype=np.float64)
    if velocity.ndim != 3:
        raise ValueError("velocity must be 3D (n_timepoints - 1, n_factors, n_samples)")
    return np.diff(velocity, axis=0)


def test_tox_compute_velocity_trajectories():
    """Test velocity computation wrapper."""
    # trajectories shape: (n_factors=1, n_samples=2, n_timepoints=4)
    trajectories = np.array(
        [[[1.0, 2.0, 4.0, 7.0],
          [0.0, -1.0, -1.0,  0.0]]],
        dtype=np.float64,
    )

    velocity = compute_velocity_trajectories(trajectories)
    expected = _expected_velocity(trajectories)

    # one timepoint shorter, and time-major
    assert velocity.shape == (trajectories.shape[2] - 1,) + trajectories.shape[:2]
    assert np.allclose(velocity, expected, atol=TOL), f"Velocity mismatch! Max diff: {np.max(np.abs(velocity - expected))}"


def test_tox_compute_acceleration_from_velocity():
    """Test acceleration computation wrapper."""
    # velocity as the routine produces it: (n_timepoints - 1 = 3, n_factors = 1, n_samples = 2)
    n_timepoints = 4
    velocity = np.asfortranarray(np.array(
        [[[1.0, -1.0]],
         [[2.0, 0.0]],
         [[3.0, 1.0]]],
        dtype=np.float64,
    ))

    acceleration = compute_acceleration_from_velocity(velocity, n_timepoints)
    expected = _expected_acceleration(velocity)

    assert acceleration.shape == (n_timepoints - 2,) + velocity.shape[1:]
    assert np.allclose(acceleration, expected, atol=TOL)


def test_tox_compute_velocity_acceleration_contributions():
    """Test velocity and acceleration contribution computation wrapper."""
    # trajectories shape: (n_factors=2, n_samples=1, n_timepoints=4)
    trajectories = np.array(
        [  # factors
            [  # samples
                [1.0, 3.0, 6.0, 10.0]  # timepoints
            ],
            [
                [1.0, 2.0, 2.0, 1.0]  # timepoints
            ]
        ],
        dtype=np.float64,
    )
    n_factors, n_samples, n_timepoints = trajectories.shape
    mode = "raw"

    result = compute_velocity_acceleration_contributions(trajectories, mode)

    C_vel = result["contrib_velocity"]
    series_vel = result["velocity_contribution_series"]
    C_acc = result["contrib_acceleration"]
    series_acc = result["acceleration_contribution_series"]

    # Output shapes: (n_samples, n_factors, n_factors, ...)
    assert C_vel.shape == (n_factors, n_factors, n_samples), f"Expected ({n_samples}, {n_factors}, {n_factors}), got {C_vel.shape}"
    assert series_vel.shape == (n_timepoints, n_factors, n_factors, n_samples), f"Expected ({n_samples}, {n_factors}, {n_factors}, {n_timepoints}), got {series_vel.shape}"
    assert C_acc.shape == (n_factors, n_factors, n_samples), f"Expected ({n_samples}, {n_factors}, {n_factors}), got {C_acc.shape}"
    assert series_acc.shape == (n_timepoints, n_factors, n_factors, n_samples), f"Expected ({n_samples}, {n_factors}, {n_factors}, {n_timepoints}), got {series_acc.shape}"

    expected_velocity = _expected_velocity(trajectories)
    expected_acceleration = _expected_acceleration(expected_velocity)

    # velocity/acceleration are time-major: (n_timepoints - 1, n_factors, n_samples)
    factor_velocity = expected_velocity[:, 0, 0]           # Factor 0, timepoints 2-4
    dependent_velocity = expected_velocity[:, 1, 0]        # Factor 1, timepoints 2-4
    factor_acceleration = expected_acceleration[:, 0, 0]   # Factor 0, timepoints 3-4
    dependent_acceleration = expected_acceleration[:, 1, 0]  # Factor 1, timepoints 3-4

    expected_series_vel = np.zeros(4)
    vel_contribs = factor_velocity * dependent_velocity
    expected_vel_contribution = vel_contribs.sum()
    expected_series_vel[1:] = vel_contribs

    expected_series_acc = np.zeros(4)
    acc_contribs = factor_acceleration * dependent_acceleration
    expected_acc_contribution = acc_contribs.sum()
    expected_series_acc[2:] = acc_contribs

    # Check sample 0, factor 0 → factor 1 contributions
    assert np.isclose(C_vel[1, 0, 0], expected_vel_contribution, atol=TOL), \
        f"C_vel[1,0,0] = {C_vel[1, 0, 0]}, expected {expected_vel_contribution}"
    assert np.allclose(series_vel[:, 1, 0, 0], expected_series_vel, atol=TOL), \
        f"series_vel mismatch"
    assert np.isclose(C_acc[1, 0, 0], expected_acc_contribution, atol=TOL), \
        f"C_acc[1,0,0] = {C_acc[1, 0, 0]}, expected {expected_acc_contribution}"
    assert np.allclose(series_acc[:, 1, 0, 0], expected_series_acc, atol=TOL), \
        f"series_acc mismatch"


def test_tox_compute_velocity_acceleration_contributions_expert():
    trajectories = np.array([
        [[1.0, 2.0, 3.0, 4.0],
         [2.0, 4.0, 6.0, 8.0]]
    ], dtype=np.float64)

    result_base = compute_velocity_acceleration_contributions(trajectories, "raw")
    result_expert = compute_velocity_acceleration_contributions_expert(trajectories, "raw")

    for key in result_base:
        assert key in result_expert, f"Missing key {key} in expert result"
        assert result_base[key].shape == result_expert[key].shape, "Expert shape differs"
        assert np.allclose(result_base[key], result_expert[key], atol=TOL), "Expert result differs"


def test_tox_compute_baselines_factor_dependent():
    """Test baseline computation wrapper across all supported modes and error cases."""

    factor = np.array([1.0, 3.0, 2.0, 4.0], dtype=np.float64)
    dependent = np.array([5.0, 7.0, 6.0, 8.0], dtype=np.float64)

    # RAW mode => zero baselines
    res_raw = compute_baselines_factor_dependent(factor, dependent, baseline_mode="raw")
    assert np.isclose(res_raw['factor_baseline'], 0.0, atol=TOL)
    assert np.isclose(res_raw['dependent_baseline'], 0.0, atol=TOL)

    # MIN mode => min values
    res_min = compute_baselines_factor_dependent(factor, dependent, baseline_mode="min")
    assert np.isclose(res_min['factor_baseline'], np.min(factor), atol=TOL)
    assert np.isclose(res_min['dependent_baseline'], np.min(dependent), atol=TOL)

    # MEAN mode => arithmetic mean
    res_mean = compute_baselines_factor_dependent(factor, dependent, baseline_mode="mean")
    assert np.isclose(res_mean['factor_baseline'], np.mean(factor), atol=TOL)
    assert np.isclose(res_mean['dependent_baseline'], np.mean(dependent), atol=TOL)

    # Mismatched lengths should raise ValueError
    assert_error(lambda: compute_baselines_factor_dependent(factor, dependent[:-1], baseline_mode=1), "Expected ValueError for mismatched lengths")

    # Invalid mode should bubble up as RuntimeError from Fortran layer
    assert_error(lambda: compute_baselines_factor_dependent(factor, dependent, baseline_mode="unknown_mode"), "Expected RuntimeError for invalid mode", ERR_INVALID_INPUT)


def test_compute_contributions():
    # -------------------------------
    # Case 1: RAW baseline
    # -------------------------------
    factor = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    dependent = np.array([2.0, 1.0, 0.0, -1.0], dtype=np.float64)
    local, total = compute_contributions(factor, dependent, baseline_mode="raw").values()

    expected_local = factor * dependent
    expected_total = sum(expected_local)
    assert np.allclose(local, expected_local, atol=TOL), "Case 1 local contributions mismatch"
    assert abs(total - expected_total) < TOL, "Case 1 total contribution mismatch"

    # -------------------------------
    # Case 2: MIN baseline
    # -------------------------------
    factor = np.array([3.0, 5.0, 2.0, 4.0], dtype=np.float64)
    dependent = np.array([1.0, 2.0, 0.0, -1.0], dtype=np.float64)
    local, total = compute_contributions(factor, dependent, baseline_mode="min").values()

    expected_local = (factor - factor.min()) * (dependent - dependent.min())
    expected_total = sum(expected_local)
    assert np.allclose(local, expected_local, atol=TOL), "Case 2 local contributions mismatch"
    assert abs(total - expected_total) < TOL, "Case 2 total contribution mismatch"

    # -------------------------------
    # Case 3: MEAN baseline
    # -------------------------------
    factor = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    dependent = np.array([4.0, 3.0, 2.0, 1.0], dtype=np.float64)
    local, total = compute_contributions(factor, dependent, baseline_mode="mean").values()

    expected_local = (factor - factor.mean()) * (dependent - dependent.mean())
    expected_total = sum(expected_local)
    assert np.allclose(local, expected_local, atol=TOL), "Case 3 local contributions mismatch"
    assert abs(total - expected_total) < TOL, "Case 3 total contribution mismatch"


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

    local, total = compute_all_contributions(trajectories, factor_indices, dependent_indices, baseline_mode="mean").values()

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

    local, total = compute_all_contributions(trajectories, factor_indices, dependent_indices, baseline_mode="min").values()

    # Baselines: min(factor)=2.0, min(dependent)=1.0
    expected_local = np.array([0.0, 0.0, 8.0], dtype=np.float64, order="F")
    expected_total = 8.0

    assert np.allclose(local[:, 0, 0, 0], expected_local, atol=TOL), "Case 2 local contributions mismatch"
    assert abs(total[0, 0, 0] - expected_total) < TOL, "Case 2 total contribution mismatch"


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
    result = perform_permutation_test(
        trajectories,
        factor_idx=factor_idx,
        dependent_idx=dependent_idx,
        sample_idx=sample_idx,
        baseline_mode=mode,
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

    result = compute_p_values(local_obs, total_obs, local_perm, total_perm)
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
    assert_error(lambda: compute_p_values(local_obs_nan, total_obs, local_perm, total_perm), "Expected RuntimeError for NaN input", ERR_NAN_INF)

    # -------------------------------
    # Case 3: Inf in permutation contributions
    # -------------------------------
    local_perm_inf = local_perm.copy(order="F")
    local_perm_inf[2,3] = np.inf
    assert_error(lambda: compute_p_values(local_obs, total_obs, local_perm_inf, total_perm), "Expected RuntimeError for Infinity input", ERR_NAN_INF)


def test_tox_compute_velocity_trajectory():
    """Test single-trajectory velocity computation wrapper."""
    trajectory = np.array([1.0, 2.0, 4.0, 7.0], dtype=np.float64)
    velocity = compute_velocity_trajectory(trajectory)

    expected = np.array([1.0, 2.0, 3.0], dtype=np.float64)

    assert velocity.shape == (trajectory.shape[0] - 1,)
    assert np.allclose(velocity, expected, atol=TOL)

    # Dimensionality check
    assert_error(lambda: compute_velocity_trajectory(trajectory.reshape(1, -1)), "Expected ValueError for 2D input")


def test_tox_compute_acceleration_from_velocity_trajectory():
    """Test single-trajectory acceleration computation wrapper."""
    n_timepoints = 4
    velocity = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    acceleration = compute_acceleration_from_velocity_trajectory(velocity, n_timepoints)

    expected = np.array([1.0, 1.0], dtype=np.float64)

    assert acceleration.shape == (n_timepoints - 2,)
    assert np.allclose(acceleration, expected, atol=TOL)

    # Dimensionality check
    assert_error(lambda: compute_acceleration_from_velocity_trajectory(velocity.reshape(1, -1), n_timepoints), "Expected ValueError for 2D input")


if __name__ == "__main__":
    run_all_tests(globals().values())
