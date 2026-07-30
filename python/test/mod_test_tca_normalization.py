"""
Test suite for TensorOmics trajectory normalization functions
Tests the three normalization functions with various edge cases
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensor_omics
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from tensor_omics import (
    normalize_variable_timeseries,
    normalize_single_trajectory,
    normalize_all_trajectories
)
from test_helpers import run_all_tests, assert_error

TOL = 1e-12


def test_normalize_variable_timeseries():
    """Test normalize_variable_timeseries function"""

    # -------------------------------
    # Case 1: Normal time series
    # -------------------------------
    v = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)

    result = normalize_variable_timeseries(v)["v_norm"]

    expected = np.array([0.0, 0.25, 0.5, 0.75, 1.0], dtype=np.float64)
    assert np.allclose(result, expected, atol=TOL), "Normal case: result mismatch"
    assert np.min(result) >= 0.0 - TOL, "Values should be >= 0"
    assert np.max(result) <= 1.0 + TOL, "Values should be <= 1"

    # -------------------------------
    # Case 2: Constant vector
    # -------------------------------
    v = np.array([3.14, 3.14, 3.14, 3.14], dtype=np.float64)
    result = normalize_variable_timeseries(v)["v_norm"]
    res_stat = normalize_variable_timeseries(v)["status"]

    expected = np.zeros_like(v)
    assert np.equal(res_stat, 210), "Constant case: status code mismatch"
    assert np.allclose(result, expected, atol=TOL), "Constant case: result mismatch"

    # -------------------------------
    # Case 3: Negative values
    # -------------------------------
    v = np.array([-10.0, -5.0, 0.0, 5.0, 10.0], dtype=np.float64)
    result = normalize_variable_timeseries(v)["v_norm"]

    expected = np.array([0.0, 0.25, 0.5, 0.75, 1.0], dtype=np.float64)
    assert np.allclose(result, expected, atol=1e-8), "Negative values: result mismatch"

    # -------------------------------
    # Case 4: Single point
    # -------------------------------
    v = np.array([42.0], dtype=np.float64)
    result = normalize_variable_timeseries(v)["v_norm"]

    assert np.abs(result[0]) < TOL, "Single point should normalize to 0"

    # -------------------------------
    # Case 5: Zero vector
    # -------------------------------
    v = np.zeros(5, dtype=np.float64)
    v, status = normalize_variable_timeseries(v).values()
    assert np.all(v == 0.0), "Zero vector should remain zero"
    assert status != 0, "Empty vectord should indicate error"


def test_normalize_single_trajectory():
    """Test normalize_single_trajectory function"""

    n_factors = 2
    n_timepoints = 4

    # -------------------------------
    # Case 1: Simple 2x4 trajectory
    # -------------------------------
    trajectory = np.array([
        [11.0, 21.0],
        [12.0, 22.0],
        [13.0, 23.0],
        [14.0, 24.0]
    ], dtype=np.float64, order='F')

    result = normalize_single_trajectory(trajectory)["trajectory_norm"]

    # Check each factor is normalized across time
    for i_factor in range(n_factors):
        factor_values = result[:, i_factor]
        assert np.abs(np.min(factor_values)) < TOL, f"Factor {i_factor+1}: min should be 0"
        assert np.abs(np.max(factor_values) - 1.0) < TOL, f"Factor {i_factor+1}: max should be 1"

    # -------------------------------
    # Case 2: Larger 3x5 trajectory
    # -------------------------------
    n_factors = 3
    n_timepoints = 5
    trajectory = np.array([
        [1.0, 2.0, 5.0],
        [10.0, 4.0, 5.0],
        [100.0, 8.0, 5.0],
        [50.0, 16.0, 5.0],
        [25.0, 32.0, 5.0]
    ], dtype=np.float64, order='F')
    result = normalize_single_trajectory(trajectory)["trajectory_norm"]

    # Factor 0: Should have min=0 at time 0, max=1 at time 2
    assert np.abs(result[0, 0]) < TOL, "Factor 1: time 0 should be min"
    assert np.abs(result[2, 0] - 1.0) < TOL, "Factor 1: time 2 should be max"

    # Factor 1: Exponential - all values should be increasing
    for i_timepoint in range(1, n_timepoints):
        assert result[i_timepoint, 1] > result[i_timepoint - 1, 1] - TOL, f"Factor 2 should be increasing at time {i_timepoint}"

    # Factor 2: Constant - should normalize to all zeros
    assert np.all(np.abs(result[:, 2]) < TOL), "Constant factor should be all zeros"


def test_normalize_all_trajectories():
    """Test normalize_all_trajectories function"""

    # -------------------------------
    # Case 1: Small 2x2x3 dataset
    # -------------------------------
    trajectories = np.array([[[101.0, 102.0, 103.0],   # Factor 0, Sample 0
                              [111.0, 112.0, 113.0]],  # Factor 0, Sample 1

                             [[201.0, 202.0, 203.0],   # Factor 1, Sample 0
                              [211.0, 212.0, 213.0]]], # Factor 1, Sample 1
                            dtype=np.float64, order='F')

    result = normalize_all_trajectories(trajectories)["trajectories_norm"]

    # Check normalization for each (factor, sample) pair
    for i in range(trajectories.shape[0]):  # factors
        for j in range(trajectories.shape[1]):  # samples
            time_series = result[i, j, :]
            assert np.abs(np.min(time_series)) < TOL, f"Factor {i}, Sample {j}: min should be 0"
            assert np.abs(np.max(time_series) - 1.0) < TOL, f"Factor {i}, Sample {j}: max should be 1"


    # -------------------------------
    # Case 2: Random larger dataset
    # -------------------------------
    np.random.seed(42)
    trajectories = np.random.randn(5, 10, 20).astype(np.float64) * 100 + 500
    trajectories = np.ascontiguousarray(trajectories)

    result = normalize_all_trajectories(trajectories)["trajectories_norm"]

    # Verify all values in [0,1]
    assert np.all(result >= 0.0 - TOL), "All values should be >= 0"
    assert np.all(result <= 1.0 + TOL), "All values should be <= 1"

    # Verify each time series is properly normalized
    for i in range(result.shape[0]):
        for j in range(result.shape[1]):
            ts = result[i, j, :]
            # Check min and max are approximately 0 and 1 (allowing for epsilon)
            if not (np.abs(np.min(ts)) < 1e-6 or np.abs(np.max(ts) - 1.0) < 1e-6):
                # This can happen if all values are equal
                if np.std(trajectories[i, j, :]) > 1e-10:  # Not constant
                    assert False, f"Factor {i}, Sample {j}: not properly normalized"

    # -------------------------------
    # Case 3: Edge case with constant factors
    # -------------------------------
    trajectories = np.ones((3, 4, 5), dtype=np.float64) * 7.0
    trajectories[1, :, :] = np.arange(5).reshape(1, 1, 5)  # Make factor 1 vary

    result = normalize_all_trajectories(trajectories)["trajectories_norm"]

    # Factor 0 and 2 should be all zeros (constant)
    assert np.all(np.abs(result[0, :, :]) < TOL), "Constant factor 0 should be all zeros"
    assert np.all(np.abs(result[2, :, :]) < TOL), "Constant factor 2 should be all zeros"

    # Factor 1 should be properly normalized
    for j in range(trajectories.shape[1]):
        assert np.abs(np.min(result[1, j, :])) < TOL, f"Factor 1, Sample {j}: min should be 0"
        assert np.abs(np.max(result[1, j, :]) - 1.0) < TOL, f"Factor 1, Sample {j}: max should be 1"


def test_error_handling():
    """Test error conditions"""

    # -------------------------------
    # Case 1: NaN input
    # -------------------------------
    v = np.array([1.0, 2.0, np.nan, 4.0], dtype=np.float64)
    assert_error(lambda: normalize_variable_timeseries(v), "Should raise RuntimeError for NaN")

    # -------------------------------
    # Case 2: Infinity input
    # -------------------------------
    v = np.array([1.0, 2.0, np.inf, 4.0], dtype=np.float64)
    assert_error(lambda: normalize_variable_timeseries(v), "Should raise RuntimeError for infinity")

    # -------------------------------
    # Case 3: Empty array
    # -------------------------------
    v = np.array([], dtype=np.float64)
    assert_error(lambda: normalize_variable_timeseries(v), "Should raise RuntimeError for empty array")


if __name__ == "__main__":
    run_all_tests(globals().values())
