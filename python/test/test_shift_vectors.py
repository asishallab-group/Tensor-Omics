#!/usr/bin/env python3
"""
Comprehensive Python test suite for shift vector field (mirrors Fortran unit tests)
Uses the modular tensoromics_functions module
"""

import numpy as np
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))
from tensoromics_functions import tox_compute_shift_vector_field

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
    shift_vectors = tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
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
    print("test_correct_family_mapping passed")

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
    error_raised = False
    # Check that error is raised
    try:
        tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
    except ValueError as e:
        error_raised = True
        assert "gene_to_centroid contains invalid family IDs" in str(e)
    assert error_raised, "Expected RuntimeError was not raised"
    print("test_invalid_family_mapping passed")

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
    shift_vectors = tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
    expected_shift_vectors = np.array([
        [1.0, 4.0],
        [2.0, 5.0],
        [3.0, 6.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0]
    ], dtype=np.float64)
    np.testing.assert_allclose(shift_vectors, expected_shift_vectors, atol=1e-12)
    print("test_zero_distance passed")

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
    shift_vectors = tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
    expected_shift_vectors = np.array([
        [10.0, 30.0, 10.0, 30.0],
        [20.0, 40.0, 20.0, 40.0],
        [-9.0, -27.0, -5.0, -23.0],
        [-18.0, -36.0, -14.0, -32.0]
    ], dtype=np.float64)
    np.testing.assert_allclose(shift_vectors, expected_shift_vectors, atol=1e-12)
    print("test_multiple_genes_per_family passed")

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
    shift_vectors = tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
    expected_shift_vectors = np.array([
        [10.0, 30.0, 50.0, 70.0],
        [20.0, 40.0, 60.0, 80.0],
        [-9.0, -27.0, -45.0, -63.0],
        [-18.0, -36.0, -54.0, -72.0]
    ], dtype=np.float64)
    np.testing.assert_allclose(shift_vectors, expected_shift_vectors, atol=1e-12)
    print("test_single_gene_per_family passed")

# 6. Test for dimension edge cases (0 genes with dimension 1 and 1 family)
def test_dimension_edge_cases():
    expression_vectors = np.empty((1, 0), dtype=np.float64)
    family_centroids = np.empty((1, 1), dtype=np.float64)
    gene_to_centroid = np.empty((0,), dtype=np.int32)
    # Check that error is raised
    try:
        tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
    except RuntimeError as e:
        error_raised = True
        assert "Empty input arrays provided." in str(e)
    assert error_raised, "Expected RuntimeError was not raised"
    print("test_dimension_edge_cases passed")

def main():
    print("=================================================")
    print("    SHIFT VECTOR FIELD FULL PYTHON INTERFACE TESTS")
    print("=================================================\n")
    test_correct_family_mapping()
    test_invalid_family_mapping()
    test_zero_distance()
    test_multiple_genes_per_family()
    test_single_gene_per_family()
    test_dimension_edge_cases()
    print("=================================================")
    print("             ALL TESTS COMPLETED")
    print("=================================================")
    print("If you see this message, all shift vector Python interface tests passed! ✓")

if __name__ == "__main__":
    main()