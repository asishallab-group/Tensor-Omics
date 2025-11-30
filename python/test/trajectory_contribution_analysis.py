"""
Test script for clock hand angle functions
Python equivalent of the R and Fortran clock hand angle tests
"""

import numpy as np
import sys
import os
from math import pi as PI

os.environ['GFORTRAN_UNBUFFERED_ALL'] = '1'

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (
    tox_trajectory_contribution,
    tox_spike_contribution,
    tox_calc_contributions,
    tox_calc_contributions_expert,
    tox_process_trajectories,
    tox_process_trajectories_flat
)

# Constants
TOL = 1e-12
MODE_NORMAL = 1
MODE_RAP = 2


def test_tox_trajectory_contribution():
    """
    Test tox_trajectory_contribution for both modes using known vectors.
    """

    # Aligned vectors
    factor = np.array([1.0, 2.0, 3.0])
    dependent = np.array([2.0, 4.0, 6.0])
    dot_prod = np.dot(factor, dependent)
    magnitude = np.sqrt(np.sum(factor**2)) * np.sqrt(np.sum(dependent**2))

    # MODE_NORMAL
    expected_normal = dot_prod / magnitude
    result_normal = tox_trajectory_contribution(factor, dependent, mode=MODE_NORMAL)
    assert np.isclose(result_normal, expected_normal, atol=1e-12), f"MODE_NORMAL failed: {result_normal} vs {expected_normal}"

    # MODE_RAP
    expected_rap = np.arccos(expected_normal)
    result_rap = tox_trajectory_contribution(factor, dependent, mode=MODE_RAP)
    assert np.isclose(result_rap, expected_rap, atol=1e-12), f"MODE_RAP failed: {result_rap} vs {expected_rap}"

    print("✅ tox_trajectory_contribution passed all tests.")


def test_tox_spike_contribution():
    """
    Test tox_spike_contribution for both modes using known vectors.
    """

    # Aligned vectors
    factor = np.array([1.0, 2.0, 3.0])
    dependent = np.array([2.0, 4.0, 6.0])
    magnitude = np.sqrt(np.sum(factor**2)) * np.sqrt(np.sum(dependent**2))

    # MODE_NORMAL
    expected_normal = (factor * dependent) / magnitude
    result_normal = tox_spike_contribution(factor, dependent, mode=MODE_NORMAL)
    assert np.allclose(result_normal, expected_normal, atol=1e-12), f"MODE_NORMAL failed: {result_normal} vs {expected_normal}"

    # MODE_RAP
    expected_rap = np.arccos(expected_normal)
    result_rap = tox_spike_contribution(factor, dependent, mode=MODE_RAP)
    assert np.allclose(result_rap, expected_rap, atol=1e-12), f"MODE_RAP failed: {result_rap} vs {expected_rap}"

    print("✅ tox_spike_contribution passed all tests.")


def test_tox_calc_contributions():
    """
    Test both tox_calc_contributions and tox_calc_contributions_expert
    using a reproducible trajectory tensor with known alignment.
    """

    # Setup: T[i,j,k] = 100*i + 10*j + k
    n_factors, n_samples, n_timepoints = 2, 3, 4
    T = np.empty((n_factors, n_samples, n_timepoints), dtype=np.float64)
    for i in range(n_factors):
        for j in range(n_samples):
            for k in range(n_timepoints):
                T[i, j, k] = 100.0*(i+1) + 10.0*(j+1) + (k+1)

    # Aligned vectors: factor=1, dependent=0 (0-based)
    i_factor = 1
    dependent_idx = 0
    mode = MODE_NORMAL

    # Reference computation
    expected_spike = np.empty((n_timepoints, n_samples), dtype=np.float64)
    expected_integrated = np.empty(n_samples, dtype=np.float64)
    for j in range(n_samples):
        f = T[i_factor, j, :]
        d = T[dependent_idx, j, :]
        mag = np.sqrt(np.sum(f**2)) * np.sqrt(np.sum(d**2))
        expected_spike[:, j] = (f * d) / mag
        expected_integrated[j] = np.sum(f * d) / mag

    # Test calc_contributions_c
    res = tox_calc_contributions(T, i_factor, dependent_idx, mode)
    assert np.allclose(res["spikes"], expected_spike, atol=1e-12), "tox_calc_contributions: spike mismatch"
    assert np.allclose(res["trajectory"], expected_integrated, atol=1e-12), "tox_calc_contributions: integrated mismatch"

    # Test calc_contributions_expert_c
    temp_factor_vector = np.empty(n_timepoints, dtype=np.float64)
    temp_dependent_vector = np.empty(n_timepoints, dtype=np.float64)
    res_expert = tox_calc_contributions_expert(
        T, i_factor, dependent_idx, mode, temp_factor_vector, temp_dependent_vector
    )

    assert np.allclose(res_expert["spikes"], expected_spike, atol=1e-12), "tox_calc_contributions_expert: spike mismatch"
    assert np.allclose(res_expert["trajectory"], expected_integrated, atol=1e-12), "tox_calc_contributions_expert: integrated mismatch"

    print("✅ Both tox_calc_contributions wrappers passed all tests.")

def test_tox_process_trajectories():
    """
    Test tox_process_trajectories with multiple factors and known patterns.
    """
    # Create test data with predictable patterns
    n_factors, n_samples, n_timepoints = 4, 5, 6
    trajectories = np.random.RandomState(42).randn(n_factors, n_samples, n_timepoints)
    
    # Set factor mask: process factors 2, 3, 4 (exclude dependent_idx=1)
    factor_mask = np.array([False, True, True, True])
    dependent_idx = 1
    mode = MODE_NORMAL
    percentile = 90.0
    
    # Call the function
    result = tox_process_trajectories(
        trajectories, factor_mask, dependent_idx, mode, percentile
    )
    
    # Verify output dimensions
    n_processed = np.sum(factor_mask)
    if factor_mask[dependent_idx - 1]:
        n_processed -= 1
    
    assert result["integrated_contribs"].shape == (n_samples, n_processed), "Integrated contributions wrong shape"
    assert result["spike_contribs"].shape == (n_timepoints, n_samples, n_processed), "Spike contributions wrong shape"
    assert result["thresholds_integrated"].shape == (n_processed,), "Integrated thresholds wrong shape"
    assert result["thresholds_spike"].shape == (n_timepoints, n_processed), "Spike thresholds wrong shape"
    assert result["outliers_integrated"].shape == (n_samples, n_processed), "Integrated outliers wrong shape"
    assert result["outliers_spike"].shape == (n_timepoints, n_samples, n_processed), "Spike outliers wrong shape"
    
    # Verify data types
    assert result["integrated_contribs"].dtype == np.float64, "Integrated contribs should be float64"
    assert result["spike_contribs"].dtype == np.float64, "Spike contribs should be float64"
    assert result["thresholds_integrated"].dtype == np.float64, "Integrated thresholds should be float64"
    assert result["thresholds_spike"].dtype == np.float64, "Spike thresholds should be float64"
    assert result["outliers_integrated"].dtype == bool, "Integrated outliers should be boolean"
    assert result["outliers_spike"].dtype == bool, "Spike outliers should be boolean"
    
    print("✅ tox_process_trajectories passed all tests.")


def test_tox_process_trajectories_flat():
    """
    Test tox_process_trajectories_flat_alloc with global percentile for spike contributions.
    """
    # Create test data
    n_factors, n_samples, n_timepoints = 3, 4, 5
    trajectories = np.random.RandomState(123).randn(n_factors, n_samples, n_timepoints)
    
    # Set factor mask: process factors 2, 3 (exclude dependent_idx=1, fortran is 1 based)
    factor_mask = np.array([False, True, True])
    dependent_idx = 1
    mode = MODE_RAP
    percentile = 85.0
    
    # Call the function
    result = tox_process_trajectories_flat(
        trajectories, factor_mask, dependent_idx, mode, percentile
    )
    
    # Verify output dimensions
    n_processed = np.sum(factor_mask)
    if factor_mask[dependent_idx - 1]:
        n_processed -= 1

    assert result["integrated_contribs"].shape == (n_samples, n_processed), "Integrated contributions wrong shape"
    assert result["spike_contribs"].shape == (n_timepoints, n_samples, n_processed), "Spike contributions wrong shape"
    assert result["thresholds_integrated"].shape == (n_processed,), "Integrated thresholds wrong shape"
    assert result["thresholds_spike"].shape == (n_processed,), "Spike thresholds wrong shape (should be 1D)"
    assert result["outliers_integrated"].shape == (n_samples, n_processed), "Integrated outliers wrong shape"
    assert result["outliers_spike"].shape == (n_timepoints, n_samples, n_processed), "Spike outliers wrong shape"
    
    # Verify data types
    assert result["integrated_contribs"].dtype == np.float64, "Integrated contribs should be float64"
    assert result["spike_contribs"].dtype == np.float64, "Spike contribs should be float64"
    assert result["thresholds_integrated"].dtype == np.float64, "Integrated thresholds should be float64"
    assert result["thresholds_spike"].dtype == np.float64, "Spike thresholds should be float64"
    assert result["outliers_integrated"].dtype == bool, "Integrated outliers should be boolean"
    assert result["outliers_spike"].dtype == bool, "Spike outliers should be boolean"
    
    # For RAP mode, contributions should be angles in [0, pi]
    if mode == MODE_RAP:
        assert np.all(result["integrated_contribs"] >= 0.0) and np.all(result["integrated_contribs"] <= PI), \
            "RAP mode integrated contributions should be angles in [0, pi]"
        assert np.all(result["spike_contribs"] >= 0.0) and np.all(result["spike_contribs"] <= PI), \
            "RAP mode spike contributions should be angles in [0, pi]"
    
    print("✅ tox_process_trajectories_flat passed all tests.")


def test_process_trajectories_edge_cases():
    """
    Test edge cases for process trajectories functions.
    """
    # Test with minimum valid dimensions
    n_factors, n_samples, n_timepoints = 2, 2, 2
    trajectories = np.random.RandomState(456).randn(n_factors, n_samples, n_timepoints)
    
    # Only one factor to process (dependent_idx=0, process factor 1)
    factor_mask = np.array([False, True])
    dependent_idx = 1
    mode = MODE_NORMAL
    percentile = 95.0
    
    # Both functions should handle this case
    result_alloc = tox_process_trajectories(
        trajectories, factor_mask, dependent_idx, mode, percentile
    )
    result_flat = tox_process_trajectories_flat(
        trajectories, factor_mask, dependent_idx, mode, percentile
    )
    
    # Should process exactly one factor
    assert result_alloc["integrated_contribs"].shape[1] == 1, "Should process exactly one factor"
    assert result_flat["integrated_contribs"].shape[1] == 1, "Should process exactly one factor"
    
    print("✅ Process trajectories edge cases passed all tests.")


def test_process_trajectories_consistency():
    """
    Test that both process trajectories functions produce consistent results for the same input.
    """
    # Create identical test data
    n_factors, n_samples, n_timepoints = 3, 4, 5
    trajectories = np.random.RandomState(789).randn(n_factors, n_samples, n_timepoints)
    
    factor_mask = np.array([True, False, True])
    dependent_idx = 2
    mode = MODE_NORMAL
    percentile = 80.0
    
    # Call both functions
    result_alloc = tox_process_trajectories(
        trajectories, factor_mask, dependent_idx, mode, percentile
    )
    result_flat = tox_process_trajectories_flat(
        trajectories, factor_mask, dependent_idx, mode, percentile
    )
    
    # Contributions should be identical (same calculation method)
    assert np.allclose(result_alloc["integrated_contribs"], result_flat["integrated_contribs"]), \
        "Integrated contributions should be identical"
    assert np.allclose(result_alloc["spike_contribs"], result_flat["spike_contribs"]), \
        "Spike contributions should be identical"
    
    # Outliers should be the same (same data, different threshold methods)
    # Note: outliers might differ due to different threshold calculation methods
    
    print("✅ Process trajectories consistency test passed.")

def main():
    print("=================================================")
    print("    TRAJECTORY CONTRIBUTION ANALYSIS PYTHON INTERFACE TESTS")
    print("=================================================")
    print()

    test_tox_spike_contribution()
    test_tox_trajectory_contribution()
    test_tox_calc_contributions()
    test_tox_process_trajectories()
    test_tox_process_trajectories_flat()
    test_process_trajectories_edge_cases()
    test_process_trajectories_consistency()

if __name__ == "__main__":
    main()
