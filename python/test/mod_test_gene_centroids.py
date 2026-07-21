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
from tensor_omics import group_centroid
from test_helpers import run_all_tests, assert_error


def test_basic_all_mode():
    n_families = 2
    vectors = np.array([[1, 3, 10, 20, 5],
                        [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    mode = 'group_all'
    expected = np.array([[3, 15],
                         [3, 15]], dtype=np.float64)
    centroids = group_centroid(vectors, gene_to_family, n_families, mode)
    np.testing.assert_allclose(centroids, expected, atol=1e-12)


def test_basic_orthologs_mode():
    n_families = 2
    vectors = np.array([[1, 3, 10, 20, 5],
                        [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    mode = 'group_orthologs'
    ortholog_set = np.array([True, False, True, True, True])
    expected = np.array([[3, 15],
                         [3, 15]], dtype=np.float64)
    centroids = group_centroid(vectors, gene_to_family, n_families, mode, ortholog_set)
    np.testing.assert_allclose(centroids, expected, atol=1e-12)


def test_empty_family():
    n_axes, n_genes, n_families = 3, 4, 2
    vectors = np.ones((n_axes, n_genes), dtype=np.float64)
    gene_to_family = np.ones(n_genes, dtype=np.int32)
    mode = 'group_all'
    expected = np.zeros((n_axes, n_families), dtype=np.float64)
    expected[:, 0] = 1.0
    centroids = group_centroid(vectors, gene_to_family, n_families, mode)
    np.testing.assert_allclose(centroids, expected, atol=1e-12)


def test_no_matching_orthologs():
    n_axes, n_genes, n_families = 2, 3, 1
    vectors = np.array([[10, 20, 30],
                        [10, 20, 30]], dtype=np.float64)
    gene_to_family = np.ones(n_genes, dtype=np.int32)
    mode = 'group_orthologs'
    ortholog_set = np.zeros(n_genes, dtype=bool)
    expected = np.zeros((n_axes, n_families), dtype=np.float64)
    centroids = group_centroid(vectors, gene_to_family, n_families,  mode, ortholog_set)
    np.testing.assert_allclose(centroids, expected, atol=1e-12)


def test_single_gene_family():
    n_families = 1
    vectors = np.array([[12.3], [-4.5], [6.7]], dtype=np.float64)
    gene_to_family = np.array([1], dtype=np.int32)
    mode = 'group_all'
    centroids = group_centroid(vectors, gene_to_family, n_families, mode)
    np.testing.assert_allclose(centroids, vectors, atol=1e-12)


def test_extreme_values():
    n_axes, n_genes, n_families = 2, 4, 1
    vectors = np.zeros((n_axes, n_genes), dtype=np.float64)
    vectors[:, 0] = [1e12, -1e-12]
    vectors[:, 1] = [-1e12, 1e-12]
    vectors[:, 2] = [0, 5]
    vectors[:, 3] = [0, -5]
    gene_to_family = np.ones(n_genes, dtype=np.int32)
    mode = 'group_all'
    expected = np.zeros((n_axes, n_families), dtype=np.float64)
    centroids = group_centroid(vectors, gene_to_family, n_families, mode)
    np.testing.assert_allclose(centroids, expected, atol=1e-12)


def test_higher_dimensions():
    n_axes, n_genes, n_families = 10, 100, 5
    vectors = np.zeros((n_axes, n_genes), dtype=np.float64)
    gene_to_family = np.zeros(n_genes, dtype=np.int32)
    for i in range(n_genes):
        vectors[:, i] = i + 1
        gene_to_family[i] = (i % n_families) + 1
    mode = 'group_all'
    centroids = group_centroid(vectors, gene_to_family, n_families, mode)
    idxs = np.where(gene_to_family == 1)[0]
    expected = np.mean(vectors[:, idxs], axis=1)
    np.testing.assert_allclose(centroids[:, 0], expected, atol=1e-12)


def test_gene_order_invariance():
    n_families = 2
    vectors1 = np.array([[1, 3, 10, 20, 5],
                         [1, 3, 10, 20, 5]], dtype=np.float64)
    mode = 'group_orthologs'
    gene_to_family1 = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    ortholog_set1 = np.array([True, False, True, True, True])
    vectors2 = np.array([[5, 10, 1, 3, 20],
                         [5, 10, 1, 3, 20]], dtype=np.float64)
    gene_to_family2 = np.array([1, 2, 1, 1, 2], dtype=np.int32)
    ortholog_set2 = np.array([True, True, True, False, True])
    centroids1 = group_centroid(vectors1, gene_to_family1, n_families, mode, ortholog_set1)
    centroids2 = group_centroid(vectors2, gene_to_family2, n_families, mode, ortholog_set2)
    np.testing.assert_allclose(centroids1, centroids2, atol=1e-12)


def test_invalid_input_arguments():
    n_axes, n_genes, n_families = 2, 5, 2
    vectors = np.array([[1, 3, 10, 20, 5],
                        [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    mode = 'group_all'
    # Invalid n_axes (empty vectors)
    assert_error(lambda: group_centroid(np.empty((0, n_genes)), gene_to_family, n_families, mode), "Expected ValueError for n_axes=0")
    # Invalid n_genes (empty gene set)
    assert_error(lambda: group_centroid(np.empty((n_axes, 0)), np.array([], dtype=np.int32), n_families, mode), "Expected ValueError for n_genes=0")
    # Invalid n_families (0 families)
    assert_error(lambda: group_centroid(vectors, gene_to_family, 0, mode), "Expected ValueError for n_families=0")


def test_invalid_family_mapping():
    n_families = 2
    vectors = np.array([[1, 3, 10, 20, 5],
                        [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family = np.array([1, 1, 2, 3, 1], dtype=np.int32)  # 3 is invalid
    mode = 'group_all'
    assert_error(lambda: group_centroid(vectors, gene_to_family, n_families, mode), "Expected ValueError for invalid family mapping")


def test_invalid_mode_string():
    n_families = 2
    vectors = np.array([[1, 3, 10, 20, 5],
                        [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    mode = ''
    assert_error(lambda: group_centroid(vectors, gene_to_family, n_families, mode), "Expected ValueError for invalid family mapping")

    mode = 'invalid_mode'
    assert_error(lambda: group_centroid(vectors, gene_to_family, n_families, mode), "Expected ValueError for mode string")


def test_missing_ortholog_set():
    n_families = 2
    vectors = np.array([[1, 3, 10, 20, 5],
                        [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    mode = 'group_orthologs'
    assert_error(lambda: group_centroid(vectors, gene_to_family, n_families, mode), "Expected ValueError for missing ortholog set")


def test_present_ortholog_set_in_all_mode():
    n_families = 2
    vectors = np.array([[1, 3, 10, 20, 5],
                        [1, 3, 10, 20, 5]], dtype=np.float64)
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    mode = 'group_all'
    ortholog_set = np.array([False, False, True, True, True])
    expected = np.array([[3, 15],
                         [3, 15]], dtype=np.float64)
    centroids = group_centroid(vectors, gene_to_family, n_families, mode, ortholog_set)
    np.testing.assert_allclose(centroids, expected, atol=1e-12)


if __name__ == "__main__":
    run_all_tests(globals().values())
