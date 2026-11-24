"""
The `tox_euclidean_distance` suite: `euclidean_distance` and `distance_to_centroid` can be called
from Python, return the documented types and shapes, and raise the documented errors, naming the
argument. Whether the numbers are right is the Fortran suite's job
(test/mod_test_tox_euclidean_distance/), as the coding guide's "Where to Test What" asks; the one
value each procedure is checked for here is there to show that the binding hands over the whole
vectors, and the family ids 1-based.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import euclidean_distance, distance_to_centroid
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_EMPTY_INPUT, ERR_INVALID_INPUT, ERR_NAN_INF

N_TISSUES, N_GENES, N_FAMILIES = 2, 4, 2
#: gene_to_fam of a gene in no family
UNASSIGNED = 0
#: the distance such a gene gets
NO_DISTANCE = -1.0


def _genes():
    """(tissues, genes): family 1's centroid plus [3, 4], family 2's plus [6, 8], unassigned,
    family 2's centroid itself."""
    return np.array([[3.0, 4.0], [16.0, 28.0], [5.0, 5.0], [10.0, 20.0]]).T.copy(order='F')


def _centroids():
    """(tissues, families): [0, 0] and [10, 20]."""
    return np.array([[0.0, 0.0], [10.0, 20.0]]).T.copy(order='F')


def _gene_to_fam():
    return np.array([1, 2, UNASSIGNED, 2], dtype=np.int32)


# =====================
# euclidean_distance
# =====================

def test_euclidean_distance():
    result = euclidean_distance(np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0]))
    assert isinstance(result, float), f"expected a float, got {type(result).__name__}"


def test_euclidean_distance_takes_the_whole_vectors():
    # the one value: [3, 4] to the origin is 5 only if both elements arrive
    result = euclidean_distance(np.array([3.0, 4.0]), np.array([0.0, 0.0]))
    assert result == 5.0, f"expected 5.0, got {result}"


def test_euclidean_distance_rejects_vectors_of_different_lengths():
    # the binding's own shape check, before the library is called
    assert_error(lambda: euclidean_distance(np.array([1.0, 2.0]), np.array([1.0, 2.0, 3.0])),
                 "vectors of different lengths", argument="vec2")


def test_euclidean_distance_rejects_empty_vectors():
    # n_elements is read off vec1, so that is the argument named
    assert_error(lambda: euclidean_distance(np.array([]), np.array([])),
                 "empty vectors", ERR_EMPTY_INPUT, argument="vec1")


def test_euclidean_distance_rejects_non_numeric_input():
    assert_error(lambda: euclidean_distance(np.array(["a", "b"]), np.array([1.0, 2.0])),
                 "strings", argument="vec1")


def test_euclidean_distance_rejects_nan():
    assert_error(lambda: euclidean_distance(np.array([1.0, np.nan]), np.array([1.0, 2.0])),
                 "a NaN in vec1", ERR_NAN_INF, argument="vec1")


def test_euclidean_distance_rejects_inf():
    assert_error(lambda: euclidean_distance(np.array([1.0, 2.0]), np.array([1.0, -np.inf])),
                 "an infinite value in vec2", ERR_NAN_INF, argument="vec2")


# =====================
# distance_to_centroid
# =====================

def test_distance_to_centroid():
    result = distance_to_centroid(_genes(), _centroids(), _gene_to_fam())
    assert isinstance(result, np.ndarray), f"expected a numpy array, got {type(result).__name__}"
    assert result.dtype == np.float64, f"expected float64, got {result.dtype}"
    assert result.shape == (N_GENES,), f"expected shape {(N_GENES,)}, got {result.shape}"
    assert not result.flags.writeable, "a result is a value, it should be read-only"


def test_distance_to_centroid_takes_family_ids_1_based():
    # the one value: id 2 is the second column of centroids, and 0 has no centroid
    result = distance_to_centroid(_genes(), _centroids(), _gene_to_fam())
    assert np.array_equal(result, [5.0, 10.0, NO_DISTANCE, 0.0]), f"got {result}"


def test_distance_to_centroid_rejects_a_family_id_past_n_families():
    gene_to_fam = _gene_to_fam()
    gene_to_fam[-1] = N_FAMILIES + 1
    assert_error(lambda: distance_to_centroid(_genes(), _centroids(), gene_to_fam),
                 "a family id past n_families", ERR_INVALID_INPUT, argument="gene_to_fam")


def test_distance_to_centroid_rejects_a_negative_family_id():
    gene_to_fam = _gene_to_fam()
    gene_to_fam[0] = -1
    assert_error(lambda: distance_to_centroid(_genes(), _centroids(), gene_to_fam),
                 "a negative family id", ERR_INVALID_INPUT, argument="gene_to_fam")


def test_distance_to_centroid_rejects_no_genes():
    # n_genes is read off genes, so that is the argument named
    assert_error(lambda: distance_to_centroid(np.empty((N_TISSUES, 0), order='F'), _centroids(),
                                              np.empty(0, dtype=np.int32)),
                 "no genes", ERR_EMPTY_INPUT, argument="genes")


def test_distance_to_centroid_rejects_no_tissues():
    assert_error(lambda: distance_to_centroid(np.empty((0, N_GENES), order='F'), np.empty((0, N_FAMILIES), order='F'),
                                              _gene_to_fam()),
                 "no tissues", ERR_EMPTY_INPUT, argument="genes")


def test_distance_to_centroid_rejects_no_families():
    # n_families is read off centroids
    assert_error(lambda: distance_to_centroid(_genes(), np.empty((N_TISSUES, 0), order='F'),
                                              np.zeros(N_GENES, dtype=np.int32)),
                 "no families", ERR_EMPTY_INPUT, argument="centroids")


def test_distance_to_centroid_rejects_nan_expression():
    genes = _genes()
    genes[-1, -1] = np.nan
    assert_error(lambda: distance_to_centroid(genes, _centroids(), _gene_to_fam()),
                 "a NaN expression", ERR_NAN_INF, argument="genes")


def test_distance_to_centroid_rejects_infinite_centroid():
    centroids = _centroids()
    centroids[-1, -1] = np.inf
    assert_error(lambda: distance_to_centroid(_genes(), centroids, _gene_to_fam()),
                 "an infinite centroid", ERR_NAN_INF, argument="centroids")


def test_distance_to_centroid_rejects_centroids_of_other_tissues():
    # the binding's own shape check, before the library is called
    assert_error(lambda: distance_to_centroid(_genes(), np.ones((N_TISSUES + 1, N_FAMILIES), order='F'),
                                              _gene_to_fam()),
                 "centroids over another number of tissues", argument="centroids")


def test_distance_to_centroid_rejects_a_mapping_of_other_genes():
    # the binding's own shape check, before the library is called
    assert_error(lambda: distance_to_centroid(_genes(), _centroids(), np.ones(N_GENES + 1, dtype=np.int32)),
                 "a mapping of another number of genes", argument="gene_to_fam")


if __name__ == "__main__":
    run_all_tests(globals().values())
