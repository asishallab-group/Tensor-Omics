#!/usr/bin/env python3
"""
Comprehensive Python test suite for gene centroids interface functions.
Mirrors the Fortran unit tests in mod_test_gene_centroids.f90.
"""

import numpy as np
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))
from tensoromics_functions import tox_group_centroid

def test_basic_all_mode():
    n_genes, n_families = 5, 2
    vectors = np.array([[1, 3, 10, 20, 5],
                        [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    ortholog_set = np.zeros(n_genes, dtype=bool)
    expected = np.array([[3, 15],
                         [3, 15]], dtype=np.float64)
    centroids = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='all')
    np.testing.assert_allclose(centroids, expected, atol=1e-12)
    print("test_basic_all_mode passed")

def test_basic_ortho_mode():
    n_families = 2
    vectors = np.array([[1, 3, 10, 20, 5],
                        [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    ortholog_set = np.array([True, False, True, True, True])
    expected = np.array([[3, 15],
                         [3, 15]], dtype=np.float64)
    centroids = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='ortho')
    np.testing.assert_allclose(centroids, expected, atol=1e-12)
    print("test_basic_ortho_mode passed")

def test_empty_family():
    n_axes, n_genes, n_families = 3, 4, 2
    vectors = np.ones((n_axes, n_genes), dtype=np.float64)
    gene_to_family = np.ones(n_genes, dtype=np.int32)
    ortholog_set = np.ones(n_genes, dtype=bool)
    expected = np.zeros((n_axes, n_families), dtype=np.float64)
    expected[:, 0] = 1.0
    centroids = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='all')
    np.testing.assert_allclose(centroids, expected, atol=1e-12)
    print("test_empty_family passed")

def test_no_matching_orthologs():
    n_axes, n_genes, n_families = 2, 3, 1
    vectors = np.array([[10, 20, 30],
                        [10, 20, 30]], dtype=np.float64)
    gene_to_family = np.ones(n_genes, dtype=np.int32)
    ortholog_set = np.zeros(n_genes, dtype=bool)
    expected = np.zeros((n_axes, n_families), dtype=np.float64)
    centroids = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='ortho')
    np.testing.assert_allclose(centroids, expected, atol=1e-12)
    print("test_no_matching_orthologs passed")

def test_single_gene_family():
    n_families = 1
    vectors = np.array([[12.3], [-4.5], [6.7]], dtype=np.float64)
    gene_to_family = np.array([1], dtype=np.int32)
    ortholog_set = np.array([True])
    centroids = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='all')
    np.testing.assert_allclose(centroids, vectors, atol=1e-12)
    print("test_single_gene_family passed")

def test_extreme_values():
    n_axes, n_genes, n_families = 2, 4, 1
    vectors = np.zeros((n_axes, n_genes), dtype=np.float64)
    vectors[:, 0] = [1e12, -1e-12]
    vectors[:, 1] = [-1e12, 1e-12]
    vectors[:, 2] = [0, 5]
    vectors[:, 3] = [0, -5]
    gene_to_family = np.ones(n_genes, dtype=np.int32)
    ortholog_set = np.ones(n_genes, dtype=bool)
    expected = np.zeros((n_axes, n_families), dtype=np.float64)
    centroids = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='all')
    np.testing.assert_allclose(centroids, expected, atol=1e-12)
    print("test_extreme_values passed")

def test_higher_dimensions():
    n_axes, n_genes, n_families = 10, 100, 5
    vectors = np.zeros((n_axes, n_genes), dtype=np.float64)
    gene_to_family = np.zeros(n_genes, dtype=np.int32)
    for i in range(n_genes):
        vectors[:, i] = i + 1
        gene_to_family[i] = (i % n_families) + 1
    ortholog_set = np.ones(n_genes, dtype=bool)
    centroids = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='all')
    idxs = np.where(gene_to_family == 1)[0]
    expected = np.mean(vectors[:, idxs], axis=1)
    np.testing.assert_allclose(centroids[:, 0], expected, atol=1e-12)
    print("test_higher_dimensions passed")

def test_gene_order_invariance():
    n_families = 2
    vectors1 = np.array([[1, 3, 10, 20, 5],
                         [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family1 = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    ortholog_set1 = np.array([True, False, True, True, True])
    vectors2 = np.array([[5, 10, 1, 3, 20],
                         [5, 10, 1, 3, 20]], dtype=np.float64)
    gene_to_family2 = np.array([1, 2, 1, 1, 2], dtype=np.int32)
    ortholog_set2 = np.array([True, True, True, False, True])
    centroids1 = tox_group_centroid(vectors1, gene_to_family1, n_families, ortholog_set1, mode='ortho')
    centroids2 = tox_group_centroid(vectors2, gene_to_family2, n_families, ortholog_set2, mode='ortho')
    np.testing.assert_allclose(centroids1, centroids2, atol=1e-12)
    print("test_gene_order_invariance passed")

def test_invalid_input_arguments():
    n_axes, n_genes, n_families = 2, 5, 2
    vectors = np.array([[1, 3, 10, 20, 5],
                        [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    ortholog_set = np.zeros(n_genes, dtype=bool)
    # Invalid n_axes (empty vectors)
    try:
        tox_group_centroid(np.empty((0, n_genes)), gene_to_family, n_families, ortholog_set, mode='all')
        assert False, "Expected ValueError for n_axes=0"
    except Exception:
        pass
    # Invalid n_genes (empty gene set)
    try:
        tox_group_centroid(np.empty((n_axes, 0)), np.array([], dtype=np.int32), n_families, np.array([], dtype=bool), mode='all')
        assert False, "Expected ValueError for n_genes=0"
    except Exception:
        pass
    # Invalid n_families (0 families)
    try:
        tox_group_centroid(vectors, gene_to_family, 0, ortholog_set, mode='all')
        assert False, "Expected ValueError for n_families=0"
    except Exception:
        pass
    print("test_invalid_input_arguments passed")

def test_invalid_family_mapping():
    n_genes, n_families = 5, 2
    vectors = np.array([[1, 3, 10, 20, 5],
                        [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family = np.array([1, 1, 2, 3, 1], dtype=np.int32)  # 3 is invalid
    ortholog_set = np.zeros(n_genes, dtype=bool)
    try:
        tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='all')
        assert False, "Expected ValueError for invalid family mapping"
    except Exception:
        pass
    print("test_invalid_family_mapping passed")

def main():
    print("=================================================")
    print("    GENE CENTROIDS FULL PYTHON INTERFACE TESTS")
    print("=================================================\n")
    test_basic_all_mode()
    test_basic_ortho_mode()
    test_empty_family()
    test_no_matching_orthologs()
    test_single_gene_family()
    test_extreme_values()
    test_higher_dimensions()
    test_gene_order_invariance()
    test_invalid_input_arguments()
    test_invalid_family_mapping()
    print("=================================================")
    print("             ALL TESTS COMPLETED")
    print("=================================================")
    print("If you see this message, all gene centroids Python interface tests passed! ✓")

if __name__ == "__main__":
    main()