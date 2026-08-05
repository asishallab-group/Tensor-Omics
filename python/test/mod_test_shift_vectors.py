#!/usr/bin/env python3
"""
Comprehensive Python test suite for shift vector field (mirrors Fortran unit tests)
Uses the modular tensor_omics module
"""

import numpy as np
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))
from tensor_omics import compute_shift_vector_field
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_EMPTY_INPUT, ERR_INVALID_INPUT


# 1. Test correct mapping between families and genes
def test_correct_family_mapping():
    expression_vectors = np.array([
        [1.0, 4.0, 7.0],
        [2.0, 5.0, 8.0],
        [3.0, 6.0, 9.0]
    ], dtype=np.float64)
    family_centroids = np.array([
        [5.0, 4.0, 3.0],
        [4.0, 3.0, 2.0],
        [3.0, 2.0, 1.0]
    ], dtype=np.float64)
    gene_to_centroid = np.array([1, 2, 3], dtype=np.int32)
    shift_vectors = compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
    shift_vectors = np.vstack([shift_vectors[:, 0, :], shift_vectors[:, 1, :]])
    # Expected: rows 0..2 = centroid, rows 3..5 = shift
    expected_centroids = np.array([
        [5.0, 4.0, 3.0],
        [4.0, 3.0, 2.0],
        [3.0, 2.0, 1.0],
        [-4.0, 0.0, 4.0],
        [-2.0, 2.0, 6.0],
        [0.0, 4.0, 8.0]
    ])
    assert shift_vectors.shape == (6, 3)
    np.testing.assert_allclose(shift_vectors, expected_centroids, atol=1e-12)


# 2. Test for invalid family id mapping raising error
def test_invalid_family_mapping():
    expression_vectors = np.array([
        [1.0, 4.0],
        [2.0, 5.0],
        [3.0, 6.0]
    ], dtype=np.float64)
    family_centroids = np.array([
        [5.0, 4.0],
        [4.0, 3.0],
        [3.0, 2.0]
    ], dtype=np.float64)
    # gene_to_centroid contains invalid mapping (3)
    gene_to_centroid = np.array([1, 3], dtype=np.int32)
    # Check that error is raised
    assert_error(lambda: compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid), "Expected error", ERR_INVALID_INPUT)


# 3. Test for zero distance between paralog and centroid
def test_zero_distance():
    expression_vectors = np.array([
        [1.0, 4.0],
        [2.0, 5.0],
        [3.0, 6.0]
    ], dtype=np.float64)
    family_centroids = np.array([
        [1.0, 4.0],
        [2.0, 5.0],
        [3.0, 6.0]
    ], dtype=np.float64)
    gene_to_centroid = np.array([1, 2], dtype=np.int32)
    shift_vectors = compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
    shift_vectors = np.vstack([shift_vectors[:, 0, :], shift_vectors[:, 1, :]])
    expected_shift_vectors = np.array([
        [1.0, 4.0],
        [2.0, 5.0],
        [3.0, 6.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0]
    ], dtype=np.float64)
    np.testing.assert_allclose(shift_vectors, expected_shift_vectors, atol=1e-12)


# 4. Test for multiple genes per family centroid
def test_multiple_genes_per_family():
    expression_vectors = np.array([
        [1.0, 3.0, 5.0, 7.0],
        [2.0, 4.0, 6.0, 8.0]
    ], dtype=np.float64)
    family_centroids = np.array([
        [10.0, 30.0],
        [20.0, 40.0]
    ], dtype=np.float64)
    gene_to_centroid = np.array([1, 2, 1, 2], dtype=np.int32)
    shift_vectors = compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
    shift_vectors = np.vstack([shift_vectors[:, 0, :], shift_vectors[:, 1, :]])
    expected_shift_vectors = np.array([
        [10.0, 30.0, 10.0, 30.0],
        [20.0, 40.0, 20.0, 40.0],
        [-9.0, -27.0, -5.0, -23.0],
        [-18.0, -36.0, -14.0, -32.0]
    ], dtype=np.float64)
    np.testing.assert_allclose(shift_vectors, expected_shift_vectors, atol=1e-12)


# 5. Test for single gene per family centroid
def test_single_gene_per_family():
    expression_vectors = np.array([
        [1.0, 3.0, 5.0, 7.0],
        [2.0, 4.0, 6.0, 8.0]
    ], dtype=np.float64)
    family_centroids = np.array([
        [10.0, 30.0, 50.0, 70.0],
        [20.0, 40.0, 60.0, 80.0]
    ], dtype=np.float64)
    gene_to_centroid = np.array([1, 2, 3, 4], dtype=np.int32)
    shift_vectors = compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
    shift_vectors = np.vstack([shift_vectors[:, 0, :], shift_vectors[:, 1, :]])
    expected_shift_vectors = np.array([
        [10.0, 30.0, 50.0, 70.0],
        [20.0, 40.0, 60.0, 80.0],
        [-9.0, -27.0, -45.0, -63.0],
        [-18.0, -36.0, -54.0, -72.0]
    ], dtype=np.float64)
    np.testing.assert_allclose(shift_vectors, expected_shift_vectors, atol=1e-12)


# 6. Test for dimension edge cases (0 genes with dimension 1 and 1 family)
def test_dimension_edge_cases():
    expression_vectors = np.empty((1, 0), dtype=np.float64)
    family_centroids = np.empty((1, 1), dtype=np.float64)
    gene_to_centroid = np.empty((0,), dtype=np.int32)
    # Check that error is raised
    assert_error(lambda: compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid), "Expected error", ERR_EMPTY_INPUT)


if __name__ == "__main__":
    run_all_tests(globals().values())
