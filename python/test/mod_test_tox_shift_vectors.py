"""
The `tox_shift_vectors` suite: `compute_shift_vector_field` can be called from Python, returns
the documented type and shape, and raises the documented error, naming the argument. Whether
the numbers are right is the Fortran suite's job (test/mod_test_tox_shift_vectors/), as the
coding guide's "Where to Test What" asks; the one value a test checks here is there to show
that the binding hands the arrays over in the right order and the family ids 1-based.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import compute_shift_vector_field
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_EMPTY_INPUT, ERR_INVALID_INPUT, ERR_NAN_INF

N_TISSUES, N_GENES, N_FAMILIES = 3, 4, 2
#: gene_to_fam of a gene in no family
UNASSIGNED = 0


def _expression_vectors():
    """(tissues, genes), every gene different."""
    return np.arange(1.0, N_TISSUES * N_GENES + 1.0).reshape((N_TISSUES, N_GENES), order='F')


def _family_centroids():
    """(tissues, families): family f is 10*f in every tissue."""
    return np.tile(np.array([10.0, 20.0]), (N_TISSUES, 1))


def _gene_to_fam():
    return np.array([1, 2, UNASSIGNED, 2], dtype=np.int32)


def test_compute_shift_vector_field():
    result = compute_shift_vector_field(_expression_vectors(), _family_centroids(), _gene_to_fam())
    assert isinstance(result, np.ndarray), f"expected a numpy array, got {type(result).__name__}"
    assert result.dtype == np.float64, f"expected float64, got {result.dtype}"
    assert result.shape == (N_TISSUES, 2, N_GENES), f"expected shape {(N_TISSUES, 2, N_GENES)}, got {result.shape}"
    assert result.flags.f_contiguous, "expected column-major (order='F')"
    assert not result.flags.writeable, "a result is a value, it should be read-only"


def test_compute_shift_vector_field_takes_family_ids_1_based():
    # the one value: id 2 is the second column of family_centroids, and [:, 0, g] the centroid
    result = compute_shift_vector_field(_expression_vectors(), _family_centroids(), _gene_to_fam())
    assert np.array_equal(result[:, 0, 1], [20.0] * N_TISSUES), f"centroid of gene 2: got {result[:, 0, 1]}"


def test_compute_shift_vector_field_rejects_a_family_id_past_n_families():
    gene_to_fam = _gene_to_fam()
    gene_to_fam[-1] = N_FAMILIES + 1
    assert_error(lambda: compute_shift_vector_field(_expression_vectors(), _family_centroids(), gene_to_fam),
                 "a family id past n_families", ERR_INVALID_INPUT, argument="gene_to_fam")


def test_compute_shift_vector_field_rejects_a_negative_family_id():
    gene_to_fam = _gene_to_fam()
    gene_to_fam[0] = -1
    assert_error(lambda: compute_shift_vector_field(_expression_vectors(), _family_centroids(), gene_to_fam),
                 "a negative family id", ERR_INVALID_INPUT, argument="gene_to_fam")


def test_compute_shift_vector_field_rejects_no_genes():
    # n_genes is read off expression_vectors, so that is the argument named
    assert_error(lambda: compute_shift_vector_field(np.empty((N_TISSUES, 0), order='F'), _family_centroids(),
                                                    np.empty(0, dtype=np.int32)),
                 "no genes", ERR_EMPTY_INPUT, argument="expression_vectors")


def test_compute_shift_vector_field_rejects_no_tissues():
    assert_error(lambda: compute_shift_vector_field(np.empty((0, N_GENES), order='F'), np.empty((0, N_FAMILIES), order='F'),
                                                    _gene_to_fam()),
                 "no tissues", ERR_EMPTY_INPUT, argument="expression_vectors")


def test_compute_shift_vector_field_rejects_no_families():
    # n_families is read off family_centroids
    assert_error(lambda: compute_shift_vector_field(_expression_vectors(), np.empty((N_TISSUES, 0), order='F'),
                                                    np.zeros(N_GENES, dtype=np.int32)),
                 "no families", ERR_EMPTY_INPUT, argument="family_centroids")


def test_compute_shift_vector_field_rejects_nan_expression():
    expression_vectors = _expression_vectors()
    expression_vectors[-1, -1] = np.nan
    assert_error(lambda: compute_shift_vector_field(expression_vectors, _family_centroids(), _gene_to_fam()),
                 "a NaN expression", ERR_NAN_INF, argument="expression_vectors")


def test_compute_shift_vector_field_rejects_infinite_centroid():
    family_centroids = _family_centroids()
    family_centroids[-1, -1] = np.inf
    assert_error(lambda: compute_shift_vector_field(_expression_vectors(), family_centroids, _gene_to_fam()),
                 "an infinite centroid", ERR_NAN_INF, argument="family_centroids")


def test_compute_shift_vector_field_rejects_centroids_of_other_tissues():
    # the binding's own shape check, before the library is called
    assert_error(lambda: compute_shift_vector_field(_expression_vectors(), np.ones((N_TISSUES + 1, N_FAMILIES)),
                                                    _gene_to_fam()),
                 "centroids over another number of tissues", argument="family_centroids")


def test_compute_shift_vector_field_rejects_a_mapping_of_other_genes():
    # the binding's own shape check, before the library is called
    assert_error(lambda: compute_shift_vector_field(_expression_vectors(), _family_centroids(),
                                                    np.ones(N_GENES + 1, dtype=np.int32)),
                 "a mapping of another number of genes", argument="gene_to_fam")


if __name__ == "__main__":
    run_all_tests(globals().values())
