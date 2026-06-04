"""
Test suite for TensorOmics trajectory normalization functions
Tests the three normalization functions with various edge cases
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from tensoromics_functions import (
    tox_normalize_variable_timeseries,
    tox_normalize_single_trajectory,
    tox_normalize_all_trajectories
)

TOL = 1e-12


def print_array_info(name, arr):
    """Print array information"""
    print(f"\n{name}:")
    print(f"  Shape: {arr.shape}")
    print(f"  Dtype: {arr.dtype}")
    print(f"  Min: {np.min(arr):.6f}, Max: {np.max(arr):.6f}")
    if arr.ndim <= 2:
        print(f"  Values:\n{arr}")


def test_tox_normalize_variable_timeseries():
    """Test normalize_variable_timeseries function"""
    print("="*50)
    print("TEST: tox_normalize_variable_timeseries")
    print("="*50)

    # -------------------------------
    # Case 1: Normal time series
    # -------------------------------
    print("\nCase 1: Normal time series")
    v = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)
    print_array_info("Input", v)

    result = tox_normalize_variable_timeseries(v)["v_norm"]
    print_array_info("Output", result)

    expected = np.array([0.0, 0.25, 0.5, 0.75, 1.0], dtype=np.float64)
    assert np.allclose(result, expected, atol=TOL), "Normal case: result mismatch"
    assert np.min(result) >= 0.0 - TOL, "Values should be >= 0"
    assert np.max(result) <= 1.0 + TOL, "Values should be <= 1"
    print("✓ Normal case passed")

    # -------------------------------
    # Case 2: Constant vector
    # -------------------------------
    print("\nCase 2: Constant vector")
    v = np.array([3.14, 3.14, 3.14, 3.14], dtype=np.float64)
    result = tox_normalize_variable_timeseries(v)["v_norm"]
    res_stat = tox_normalize_variable_timeseries(v)["status"]
    print(f"Status code: {res_stat}")
    print_array_info("Output", result)

    expected = np.zeros_like(v)
    assert np.equal(res_stat, 210), "Constant case: status code mismatch"
    assert np.allclose(result, expected, atol=TOL), "Constant case: result mismatch"
    print("✓ Constant case passed")

    # -------------------------------
    # Case 3: Negative values
    # -------------------------------
    print("\nCase 3: Negative values")
    v = np.array([-10.0, -5.0, 0.0, 5.0, 10.0], dtype=np.float64)
    result = tox_normalize_variable_timeseries(v)["v_norm"]
    print_array_info("Output", result)

    expected = np.array([0.0, 0.25, 0.5, 0.75, 1.0], dtype=np.float64)
    assert np.allclose(result, expected, atol=1e-8), "Negative values: result mismatch"
    print("✓ Negative values case passed")

    # -------------------------------
    # Case 4: Single point
    # -------------------------------
    print("\nCase 4: Single point")
    v = np.array([42.0], dtype=np.float64)
    result = tox_normalize_variable_timeseries(v)["v_norm"]
    print_array_info("Output", result)

    assert np.abs(result[0]) < TOL, "Single point should normalize to 0"
    print("✓ Single point case passed")

    # -------------------------------
    # Case 5: Zero vector
    # -------------------------------
    print("\nCase 5: Zero vector")
    v = np.zeros(5, dtype=np.float64)
    try:
        result = tox_normalize_variable_timeseries(v)["v_norm"]
        assert np.all(v == 0.0), "Zero vector should remain zero"
        print("Zero vector case passed")
    except RuntimeError:
        print("Zero vector raised error")

    print("\n" + "="*50)
    print("✓ ALL tox_normalize_variable_timeseries TESTS PASSED")
    print("="*50)


def test_tox_normalize_single_trajectory():
    """Test normalize_single_trajectory function"""
    print("\n" + "="*50)
    print("TEST: tox_normalize_single_trajectory")
    print("="*50)

    n_factors = 2
    n_timepoints = 4

    # -------------------------------
    # Case 1: Simple 2x4 trajectory
    # -------------------------------
    print("\nCase 1: Simple 2x4 trajectory (2 factors, 4 timepoints)")
    trajectory = np.array([
        [11.0, 21.0],
        [12.0, 22.0],
        [13.0, 23.0],
        [14.0, 24.0]
    ], dtype=np.float64, order='F')
    print_array_info("Input", trajectory)

    result = tox_normalize_single_trajectory(trajectory)["traj_norm"]
    print_array_info("Output", result)

    # Check each factor is normalized across time
    for i_factor in range(n_factors):
        factor_values = result[:, i_factor]
        assert np.abs(np.min(factor_values)) < TOL, f"Factor {i_factor+1}: min should be 0"
        assert np.abs(np.max(factor_values) - 1.0) < TOL, f"Factor {i_factor+1}: max should be 1"
        print(f"✓ Factor {i_factor+1}: correctly normalized to [0,1]")

    # -------------------------------
    # Case 2: Larger 3x5 trajectory
    # -------------------------------
    n_factors = 3
    n_timepoints = 5
    print("\nCase 2: 3x5 trajectory with varying patterns")
    trajectory = np.array([
        [1.0, 2.0, 5.0],
        [10.0, 4.0, 5.0],
        [100.0, 8.0, 5.0],
        [50.0, 16.0, 5.0],
        [25.0, 32.0, 5.0]
    ], dtype=np.float64, order='F')
    result = tox_normalize_single_trajectory(trajectory)["traj_norm"]

    # Factor 0: Should have min=0 at time 0, max=1 at time 2
    assert np.abs(result[0, 0]) < TOL, "Factor 1: time 0 should be min"
    assert np.abs(result[2, 0] - 1.0) < TOL, "Factor 1: time 2 should be max"

    # Factor 1: Exponential - all values should be increasing
    for i_timepoint in range(1, n_timepoints):
        assert result[i_timepoint, 1] > result[i_timepoint - 1, 1] - TOL, f"Factor 2 should be increasing at time {i_timepoint}"

    # Factor 2: Constant - should normalize to all zeros
    assert np.all(np.abs(result[:, 2]) < TOL), "Constant factor should be all zeros"

    print("✓ All factors correctly normalized")

    print("\n" + "="*50)
    print("✓ ALL tox_normalize_single_trajectory TESTS PASSED")
    print("="*50)


def test_tox_normalize_all_trajectories():
    """Test normalize_all_trajectories function"""
    print("\n" + "="*50)
    print("TEST: tox_normalize_all_trajectories")
    print("="*50)

    # -------------------------------
    # Case 1: Small 2x2x3 dataset
    # -------------------------------
    print("\nCase 1: 2x2x3 dataset (2 factors, 2 samples, 3 timepoints)")
    trajectories = np.array([[[101.0, 102.0, 103.0],   # Factor 0, Sample 0
                              [111.0, 112.0, 113.0]],  # Factor 0, Sample 1

                             [[201.0, 202.0, 203.0],   # Factor 1, Sample 0
                              [211.0, 212.0, 213.0]]], # Factor 1, Sample 1
                            dtype=np.float64, order='F')
    print_array_info("Input", trajectories)

    result = tox_normalize_all_trajectories(trajectories)["traj_norm"]
    print_array_info("Output", result)

    # Check normalization for each (factor, sample) pair
    for i in range(trajectories.shape[0]):  # factors
        for j in range(trajectories.shape[1]):  # samples
            time_series = result[i, j, :]
            assert np.abs(np.min(time_series)) < TOL, f"Factor {i}, Sample {j}: min should be 0"
            assert np.abs(np.max(time_series) - 1.0) < TOL, f"Factor {i}, Sample {j}: max should be 1"

    print("✓ All (factor, sample) pairs correctly normalized across time")

    # -------------------------------
    # Case 2: Random larger dataset
    # -------------------------------
    print("\nCase 2: Random 5x10x20 dataset")
    np.random.seed(42)
    trajectories = np.random.randn(5, 10, 20).astype(np.float64) * 100 + 500
    trajectories = np.ascontiguousarray(trajectories)

    result = tox_normalize_all_trajectories(trajectories)["traj_norm"]

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

    print(f"✓ Random dataset normalized successfully")
    print(f"  Min of all results: {np.min(result):.6f}")
    print(f"  Max of all results: {np.max(result):.6f}")

    # -------------------------------
    # Case 3: Edge case with constant factors
    # -------------------------------
    print("\nCase 3: Dataset with constant factors")
    trajectories = np.ones((3, 4, 5), dtype=np.float64) * 7.0
    trajectories[1, :, :] = np.arange(5).reshape(1, 1, 5)  # Make factor 1 vary

    result = tox_normalize_all_trajectories(trajectories)["traj_norm"]

    # Factor 0 and 2 should be all zeros (constant)
    assert np.all(np.abs(result[0, :, :]) < TOL), "Constant factor 0 should be all zeros"
    assert np.all(np.abs(result[2, :, :]) < TOL), "Constant factor 2 should be all zeros"

    # Factor 1 should be properly normalized
    for j in range(trajectories.shape[1]):
        assert np.abs(np.min(result[1, j, :])) < TOL, f"Factor 1, Sample {j}: min should be 0"
        assert np.abs(np.max(result[1, j, :]) - 1.0) < TOL, f"Factor 1, Sample {j}: max should be 1"

    print("✓ Constant factors handled correctly")

    print("\n" + "="*50)
    print("✓ ALL tox_normalize_all_trajectories TESTS PASSED")
    print("="*50)


def test_error_handling():
    """Test error conditions"""
    print("\n" + "="*50)
    print("TEST: Error Handling")
    print("="*50)

    # -------------------------------
    # Case 1: NaN input
    # -------------------------------
    print("\nCase 1: NaN in input")
    v = np.array([1.0, 2.0, np.nan, 4.0], dtype=np.float64)
    try:
        tox_normalize_variable_timeseries(v)
        raise AssertionError("Should raise RuntimeError for NaN")
    except RuntimeError:
        print("✓ NaN input correctly raised error")

    # -------------------------------
    # Case 2: Infinity input
    # -------------------------------
    print("\nCase 2: Infinity in input")
    v = np.array([1.0, 2.0, np.inf, 4.0], dtype=np.float64)
    try:
        tox_normalize_variable_timeseries(v)
        raise AssertionError("Should raise RuntimeError for infinity")
    except RuntimeError:
        print("✓ Infinity input correctly raised error")

    # -------------------------------
    # Case 3: Empty array
    # -------------------------------
    print("\nCase 3: Empty array")
    v = np.array([], dtype=np.float64)
    try:
        tox_normalize_variable_timeseries(v)
        raise AssertionError("Should raise RuntimeError for empty array")
    except RuntimeError:
        print("✓ Empty array correctly raised error")

    print("\n" + "="*50)
    print("✓ ALL ERROR HANDLING TESTS PASSED")
    print("="*50)


def main():
    """Run all normalization tests"""
    print("TensorOmics Trajectory Normalization Test Suite")
    print("="*60)

    try:
        test_tox_normalize_variable_timeseries()
        test_tox_normalize_single_trajectory()
        test_tox_normalize_all_trajectories()
        test_error_handling()

        print("\n" + "="*60)
        print("SUCCESS: All normalization tests passed! 🎉")
        print("="*60)

    except Exception as e:
        print(f"\nFAILURE: Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())