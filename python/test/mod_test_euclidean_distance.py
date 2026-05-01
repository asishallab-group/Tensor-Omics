#!/usr/bin/env python3
"""
Comprehensive Python test suite for Euclidean distance functions
Uses tensoromics_functions.py wrapper functions (mirrors R euclidean_distance.R tests)
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import tox_euclidean_distance, tox_distance_to_centroid
from test_helpers import run_all_tests, assert_error


# =====================
# Tests for euclidean_distance
# =====================

def test_euclidean_distance_3d():
    """Test simple 3D vectors"""
    vec1 = np.array([1.0, 2.0, 3.0])
    vec2 = np.array([4.0, 5.0, 6.0])

    result = tox_euclidean_distance(vec1, vec2)
    expected = np.linalg.norm(vec1 - vec2)

    assert np.isclose(result, expected)


def test_euclidean_distance_to_origin():
    """Test 2D vector to origin (3-4-5 triangle)"""
    vec1 = np.array([3.0, 4.0])
    vec2 = np.array([0.0, 0.0])

    result = tox_euclidean_distance(vec1, vec2)
    expected = 5.0

    assert np.isclose(result, expected)


def test_euclidean_distance_identical():
    """Test identical vectors"""
    vec1 = np.array([1.5, 2.7, 3.9])
    vec2 = vec1.copy()  # identical

    result = tox_euclidean_distance(vec1, vec2)
    expected = 0.0

    assert np.isclose(result, expected)


def test_euclidean_distance_high_dimensional():
    """Test high-dimensional vectors"""
    # 100-dimensional vectors
    d = 100
    vec1 = np.arange(1, d+1, dtype=np.float64)  # [1, 2, 3, ..., 100]
    vec2 = np.arange(2, d+2, dtype=np.float64)  # [2, 3, 4, ..., 101] (shift by 1)

    result = tox_euclidean_distance(vec1, vec2)
    expected = np.sqrt(d)  # sqrt(100 * 1^2) = 10
    assert np.isclose(result, expected)


def test_euclidean_distance_invalid_inputs():
    """Test invalid inputs (should throw errors)"""
    # Test different lengths
    assert_error(lambda: tox_euclidean_distance(np.array([1, 2]), np.array([1, 2, 3])), "Expected error")

    # Test empty vectors
    assert_error(lambda: tox_euclidean_distance(np.array([]), np.array([])), "Expected error")

    # Test non-numeric input
    assert_error(lambda: tox_euclidean_distance(np.array(["a", "b"]), np.array([1, 2])), "Expected error")


def test_euclidean_distance_1d():
    """Test single-dimensional vectors"""
    vec1 = np.array([5.0])
    vec2 = np.array([2.0])

    result = tox_euclidean_distance(vec1, vec2)
    expected = 3.0

    assert np.isclose(result, expected)


# =====================
# Tests for distance_to_centroid
# =====================
def test_distance_to_centroid_basic():
    """Test basic distance to centroid functionality"""

    # Gene expression data (genes as columns, dimensions as rows)
    # Gene 1: [1, 0, 0] - Family 1
    # Gene 2: [0, 1, 0] - Family 1  
    # Gene 3: [3, 0, 0] - Family 2
    # Gene 4: [0, 3, 0] - Family 2
    genes = np.array([
        [1.0, 0.0, 0.0],  # Gene 1
        [0.0, 1.0, 0.0],  # Gene 2
        [3.0, 0.0, 0.0],  # Gene 3
        [0.0, 3.0, 0.0]   # Gene 4
    ], dtype=np.float64).T.copy(order="F")

    # Family centroids
    # Family 1 centroid: [0.5, 0.5, 0.0]
    # Family 2 centroid: [1.5, 1.5, 0.0]
    centroids = np.array([
        [0.5, 0.5, 0.0],  # Family 1
        [1.5, 1.5, 0.0]   # Family 2
    ], dtype=np.float64).T.copy(order="F")

    # Gene-to-family mapping (1-based)
    gene_to_fam = np.array([1, 1, 2, 2], dtype=np.int32)

    result = tox_distance_to_centroid(genes, centroids, gene_to_fam)

    # Expected distances
    # Gene 1: [1,0,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
    # Gene 2: [0,1,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
    # Gene 3: [3,0,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
    # Gene 4: [0,3,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
    expected = np.array([np.sqrt(0.5**2 + 0.5**2), np.sqrt(0.5**2 + 0.5**2),
                         np.sqrt(1.5**2 + 1.5**2), np.sqrt(1.5**2 + 1.5**2)])

    assert all(np.isclose(result, expected)), "Expected distances mismatch"


def test_distance_to_centroid_invalid_families():
    """Test handling invalid family indices (should return -1 for invalid genes)"""
    # Gene data
    genes = np.array([
        [1.0, 2.0],  # Gene 1
        [3.0, 4.0],  # Gene 2
        [5.0, 6.0]   # Gene 3
    ], dtype=np.float64).T.copy(order="F")

    # Centroids
    centroids = np.array([
        [0.0, 0.0],  # Family 1
        [1.0, 1.0]   # Family 2
    ], dtype=np.float64).T.copy(order="F")

    # Mixed family assignments: valid (1), invalid (3), no family (0)
    gene_to_fam = np.array([1, 0, 0], dtype=np.int32)  # family 0 = no assignment

    result = tox_distance_to_centroid(genes, centroids, gene_to_fam)

    expected = [np.sqrt(sum((genes[:, 0] - centroids[:, 0]) ** 2)), -1, -1]
    assert all(np.isclose(result, expected)), "Expected distances mismatch"


def test_distance_to_centroid_performance():
    """Test performance with realistic genomic data size"""

    n_genes = 1000
    n_families = 50
    n_tissues = 20

    # Generate random-like data
    np.random.seed(12345)
    genes = np.random.randn(n_tissues * n_genes).reshape((n_tissues, n_genes), order="F")
    centroids = np.random.randn(n_tissues * n_families).reshape((n_tissues, n_families), order="F")
    gene_to_fam = np.random.randint(1, n_families + 1, n_genes, dtype=np.int32)

    # Time the operation
    import time
    start_time = time.time()

    result = tox_distance_to_centroid(genes, centroids, gene_to_fam)

    end_time = time.time()
    elapsed = end_time - start_time

    assert elapsed < 1, "Should be faster than a second"

    # Verify all distances are positive
    assert np.all(result >= 0)
    assert len(result) == n_genes


def test_distance_to_centroid_input_validation():
    """Test input validation for distance_to_centroid"""

    genes = np.array([
        [1.0, 2.0],  # Gene 1
        [3.0, 4.0],  # Gene 2
        [5.0, 6.0]   # Gene 3
    ], dtype=np.float64).T.copy(order="F")

    # Centroids
    centroids = np.array([
        [0.0, 0.0],  # Family 1
        [1.0, 1.0]   # Family 2
    ], dtype=np.float64).T.copy(order="F")

    gene_to_fam = np.array([1, 0, 0], dtype=np.int32)  # family 0 = no assignment

    assert_error(lambda: tox_distance_to_centroid(genes, centroids, -gene_to_fam), "Expected error")


if __name__ == "__main__":
    run_all_tests(globals().values())
