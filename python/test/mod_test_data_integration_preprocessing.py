"""
Comprehensive Python test suite for data integration preprocessing functions in tensoromics.
Uses tensor_omics.py wrapper function (mirrors Fortran test suite)
"""
import numpy as np
import sys
import os
from pathlib import Path

# Add parent directory to path to import tensor_omics
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from test_helpers import run_all_tests
from tensor_omics import euclidean_distance, distance_to_centroid

# Import your functions
from tensor_omics import (
    compute_gene_means,
    compute_residuals,
    pool_study_means,
    construct_neighborhoods,
    pool_means_expert,
    calc_neighborhood_size,
    construct_neighborhoods
)


def test_compute_gene_means():
    """Test compute_gene_means function"""

    # Test 1: Basic 3x3 matrix
    expr = np.array([[10, 20, 30],
                     [12, 22, 32],
                     [14, 24, 34]], dtype=np.float64, order='F')

    means = compute_gene_means(expr)

    # Verify calculations
    expected = np.array([12.0, 22.0, 32.0])  # (10+12+14)/3, (20+22+24)/3, (30+32+34)/3
    np.testing.assert_array_almost_equal(means, expected, decimal=10)
    assert means.shape == (3,)

    # Test 2: With NaN values
    expr_with_nan = np.array([[10, 20, np.nan],
                              [12, np.nan, 32],
                              [14, 24, 34]], dtype=np.float64, order='F')

    means_with_nan = compute_gene_means(expr_with_nan)

    # Gene 1: (10+12+14)/3 = 12
    # Gene 2: (20+24)/2 = 22 (NaN excluded)
    # Gene 3: (32+34)/2 = 33 (NaN excluded)
    expected_nan = np.array([12.0, 22.0, 33.0])
    np.testing.assert_array_almost_equal(means_with_nan, expected_nan, decimal=10)



def test_compute_residuals():
    """Test compute_residuals function"""

    expr = np.array([[10, 20, 30],
                     [12, 22, 32],
                     [14, 24, 34]], dtype=np.float64, order='F')

    means = np.array([12.0, 22.0, 32.0], dtype=np.float64, order='F')

    residuals = compute_residuals(expr, means)

    # Verify calculations: residuals = expr - means
    expected = np.array([[10-12, 20-22, 30-32],
                         [12-12, 22-22, 32-32],
                         [14-12, 24-22, 34-32]], dtype=np.float64)

    np.testing.assert_array_almost_equal(residuals, expected, decimal=10)
    assert residuals.shape == (3, 3)

    # Test with NaN values
    expr_nan = np.array([[10, np.nan, 30],
                         [12, 22, np.nan],
                         [np.nan, 24, 34]], dtype=np.float64, order='F')

    means_nan = np.array([11.0, 23.0, 32.0], dtype=np.float64, order='F')

    residuals_nan = compute_residuals(expr_nan, means_nan)

    # Check that residuals at NaN positions are also NaN
    assert np.isnan(residuals_nan[0, 1])
    assert np.isnan(residuals_nan[1, 2])
    assert np.isnan(residuals_nan[2, 0])



def test_pool_means():
    """Test pool_study_means function"""

    # Test 1: Basic pooling
    mean_S1 = np.array([1.0, 3.0, 5.0, 7.0, 9.0], dtype=np.float64, order='F')
    mean_S2 = np.array([2.0, 4.0, 6.0, 8.0, 10.0], dtype=np.float64, order='F')
    n_points = 3

    def pooled(mean_S1, mean_S2, n_points):
        """The expert entry point, given the pooling and sorting the other one does."""
        pool = np.concatenate([mean_S1, mean_S2])
        perm = (np.argsort(pool, kind="mergesort").astype(np.int32) + 1)
        return pool_means_expert(pool, perm, n_points)

    for pool in (pool_study_means, pooled):
        result = pool(mean_S1, mean_S2, n_points)

        # Verify n_pool (no NaN values, so should be 10)
        assert result['n_pool'] == 10

        # Verify x_star has correct length
        assert len(result['x_star']) == n_points

        # Test 2: With NaN values
        mean_S1_nan = np.array([1.0, np.nan, 5.0, 7.0, 9.0], dtype=np.float64, order='F')
        mean_S2_nan = np.array([2.0, 4.0, np.nan, 8.0, 10.0], dtype=np.float64, order='F')

        result_nan = pool(mean_S1_nan, mean_S2_nan, n_points)

        # Should exclude NaN values: 4 from S1 + 4 from S2 = 8
        assert result_nan['n_pool'] == 8



def test_construct_neighborhoods():
    """Test construct_neighborhoods function"""

    # Create test data
    n_points = 3
    x_star = np.array([10.0, 20.0, 30.0], dtype=np.float64, order='F')

    # 5 genes, 2 replicates
    mean_S = np.array([8.0, 12.0, 18.0, 22.0, 28.0], dtype=np.float64, order='F')
    resid_S = np.array([[1.0, -1.0, 2.0, -2.0, 3.0],
                        [-1.0, 1.0, -2.0, 2.0, -3.0]],
                       dtype=np.float64, order='F')

    n_pool = 10
    neighborhood_size = 50  # Explicit size

    expected_n_neighbors = calc_neighborhood_size(n_pool, n_points, mean_S, neighborhood_size)
    result = construct_neighborhoods(x_star, mean_S, resid_S, expected_n_neighbors)

    # Basic assertions
    assert result['neighborhood_residuals'].shape[2] == n_points
    assert result['neighborhood_residuals'].shape[1] == expected_n_neighbors
    assert result['neighborhood_residuals'].shape[0] == 2
    assert result['neighborhood_indices'].shape[1] == n_points
    assert result['neighborhood_indices'].shape[0] == expected_n_neighbors

    # Check that indices are 1-based (from Fortran)
    # Should be >= 1 (valid) or -1 (invalid)
    indices = result['neighborhood_indices']
    assert np.all(indices[indices != -1] >= 1)

    # Test with automatic neighborhood size
    expected_n_neighbors = calc_neighborhood_size(n_pool, n_points, mean_S)
    result_auto = construct_neighborhoods(x_star, mean_S, resid_S, expected_n_neighbors)
    assert result_auto['neighborhood_residuals'].shape[1] == expected_n_neighbors



def test_integration():
    """Test integration of multiple functions"""

    # Create test expression data
    expr_S1 = np.array([[10.0, 20.0, 30.0, 40.0],
                        [12.0, 22.0, 32.0, 42.0],
                        [14.0, 24.0, 34.0, 44.0]],
                       dtype=np.float64, order='F')

    expr_S2 = np.array([[15.0, 25.0, 35.0, 45.0],
                        [17.0, 27.0, 37.0, 47.0],
                        [19.0, 29.0, 39.0, 49.0],
                        [21.0, 31.0, 41.0, 51.0]],
                       dtype=np.float64, order='F')

    # Step 1: Compute means for both studies
    means_S1 = compute_gene_means(expr_S1)
    means_S2 = compute_gene_means(expr_S2)

    assert means_S1.shape == (4,)
    assert means_S2.shape == (4,)

    # Step 2: Compute residuals
    resid_S1 = compute_residuals(expr_S1, means_S1)
    resid_S2 = compute_residuals(expr_S2, means_S2)

    assert resid_S1.shape == expr_S1.shape
    assert resid_S2.shape == expr_S2.shape

    # Step 3: Pool means
    n_points = 5
    pool_result = pool_study_means(means_S1, means_S2, n_points)

    assert pool_result['n_pool'] == 8  # 4 + 4, no NaN values
    assert len(pool_result['x_star']) == n_points

    # Step 4: Construct neighborhoods
    n_neighbors = calc_neighborhood_size(
        pool_result['n_pool'], n_points, means_S1, desired_size=50)
    neighborhoods = construct_neighborhoods(
        pool_result['x_star'], means_S1, resid_S1, n_neighbors
    )


if __name__ == "__main__":
    run_all_tests(globals().values())
