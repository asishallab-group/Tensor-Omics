"""
The `tox_gene_centroids` suite: each published procedure can be called from Python, returns the
documented type and shape, and raises the documented error. Whether the numbers are right is the
Fortran suite's job (test/mod_test_tox_gene_centroids.F90), as the coding guide's "Where to Test
What" asks.

The one value per procedure picks single columns, so no arithmetic is involved: it only shows that
gene indices and family ids cross the binding 1-based and that the axes stay along the rows.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    group_centroid_all,
    group_centroid_orthologs,
    mean_vector,
)
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import (
    ERR_EMPTY_INPUT,
    ERR_INVALID_INPUT,
    ERR_NAN_INF,
)

N_AXES, N_GENES = 2, 4


def _assert_result_array(result, shape, dtype, name, fortran_order):
    assert isinstance(result, np.ndarray), f"{name}: expected a numpy array, got {type(result).__name__}"
    assert result.dtype == dtype, f"{name}: expected {np.dtype(dtype)}, got {result.dtype}"
    assert result.shape == shape, f"{name}: expected shape {shape}, got {result.shape}"
    assert not result.flags.writeable, f"{name}: a result is a value, it should be read-only"
    if fortran_order:
        assert result.flags.f_contiguous, f"{name}: expected column-major (order='F')"
    else:
        assert result.flags.c_contiguous, f"{name}: expected a contiguous array"


def _expression_vectors():
    """(axes, genes) in C order, so the binding has to convert it; every column is distinct."""
    return np.array([[1.0, 2.0, 3.0, 4.0],
                     [10.0, 20.0, 30.0, 40.0]])


# -----------------------------------------------------------------------------------------------
# mean_vector
# -----------------------------------------------------------------------------------------------

def test_mean_vector():
    centroid = mean_vector(_expression_vectors(), [3])
    _assert_result_array(centroid, (N_AXES,), np.float64, "centroid", fortran_order=False)
    # the one value: a single 1-based index selects that column, the axes down its rows
    assert np.array_equal(centroid, [3.0, 30.0]), f"centroid: got {centroid}"


def test_mean_vector_accepts_an_empty_selection():
    # no genes selected is allowed (n_selected_genes >= 0), not an error
    centroid = mean_vector(_expression_vectors(), np.empty(0, dtype=np.int32))
    _assert_result_array(centroid, (N_AXES,), np.float64, "centroid", fortran_order=False)


def test_mean_vector_rejects_an_index_out_of_range():
    # indices are 1-based: 0 and n_genes + 1 both lie outside
    assert_error(lambda: mean_vector(_expression_vectors(), [0]), "gene index 0", ERR_INVALID_INPUT)
    assert_error(lambda: mean_vector(_expression_vectors(), [N_GENES + 1]), "gene index n_genes + 1",
                 ERR_INVALID_INPUT)


def test_mean_vector_rejects_more_indices_than_genes():
    # n_selected_genes, read off gene_indices, may not exceed n_genes
    assert_error(lambda: mean_vector(_expression_vectors(), np.ones(N_GENES + 1, dtype=np.int32)),
                 "more indices than genes", ERR_INVALID_INPUT)


def test_mean_vector_rejects_no_axes():
    assert_error(lambda: mean_vector(np.empty((0, N_GENES)), [1]), "no axes", ERR_EMPTY_INPUT)


def test_mean_vector_rejects_nan():
    expression_vectors = _expression_vectors()
    expression_vectors[1, 2] = np.nan
    assert_error(lambda: mean_vector(expression_vectors, [1]), "a NaN expression", ERR_NAN_INF)


def test_mean_vector_rejects_a_vector():
    # the binding's own rank check, before the library is called
    assert_error(lambda: mean_vector(np.ones(3), [1]), "expression_vectors must be 2-D")


# -----------------------------------------------------------------------------------------------
# group_centroid_orthologs
# -----------------------------------------------------------------------------------------------

def test_group_centroid_orthologs():
    # gene 2 is unassigned (0); gene 4 is in family 2 but not an ortholog
    gene_to_family = [2, 0, 1, 2]
    ortholog_set = [True, True, True, False]
    centroid_matrix = group_centroid_orthologs(_expression_vectors(), gene_to_family, 2, ortholog_set)
    _assert_result_array(centroid_matrix, (N_AXES, 2), np.float64, "centroid_matrix", fortran_order=True)
    # the one value: each family keeps exactly one ortholog, so its centroid is that gene's column --
    # family 1 is gene 3, family 2 is gene 1
    assert np.array_equal(centroid_matrix, [[3.0, 1.0], [30.0, 10.0]]), f"centroid_matrix: got {centroid_matrix}"


def test_group_centroid_orthologs_rejects_a_family_out_of_range():
    ortholog_set = np.ones(N_GENES, dtype=bool)
    assert_error(lambda: group_centroid_orthologs(_expression_vectors(), [1, 3, 1, 1], 2, ortholog_set),
                 "family id n_families + 1", ERR_INVALID_INPUT)
    assert_error(lambda: group_centroid_orthologs(_expression_vectors(), [1, -1, 1, 1], 2, ortholog_set),
                 "a negative family id", ERR_INVALID_INPUT)


def test_group_centroid_orthologs_rejects_no_families():
    assert_error(lambda: group_centroid_orthologs(_expression_vectors(), np.zeros(N_GENES, dtype=np.int32), 0,
                                                  np.ones(N_GENES, dtype=bool)),
                 "no families", ERR_EMPTY_INPUT)


def test_group_centroid_orthologs_rejects_no_genes():
    assert_error(lambda: group_centroid_orthologs(np.empty((N_AXES, 0)), np.empty(0, dtype=np.int32), 1,
                                                  np.empty(0, dtype=bool)),
                 "no genes", ERR_EMPTY_INPUT)


def test_group_centroid_orthologs_rejects_infinity():
    expression_vectors = _expression_vectors()
    expression_vectors[0, 3] = np.inf
    assert_error(lambda: group_centroid_orthologs(expression_vectors, np.ones(N_GENES, dtype=np.int32), 1,
                                                  np.ones(N_GENES, dtype=bool)),
                 "an infinite expression", ERR_NAN_INF)


def test_group_centroid_orthologs_rejects_a_short_ortholog_set():
    # the binding's own extent check, before the library is called
    assert_error(lambda: group_centroid_orthologs(_expression_vectors(), np.ones(N_GENES, dtype=np.int32), 1,
                                                  np.ones(N_GENES - 1, dtype=bool)),
                 "ortholog_set shorter than n_genes")


# -----------------------------------------------------------------------------------------------
# group_centroid_all
# -----------------------------------------------------------------------------------------------

def test_group_centroid_all():
    # genes 2 and 4 are unassigned (0)
    centroid_matrix = group_centroid_all(_expression_vectors(), [2, 0, 1, 0], 2)
    _assert_result_array(centroid_matrix, (N_AXES, 2), np.float64, "centroid_matrix", fortran_order=True)
    # the one value: each family holds one gene, so its centroid is that gene's column --
    # family 1 is gene 3, family 2 is gene 1
    assert np.array_equal(centroid_matrix, [[3.0, 1.0], [30.0, 10.0]]), f"centroid_matrix: got {centroid_matrix}"


def test_group_centroid_all_rejects_a_family_out_of_range():
    assert_error(lambda: group_centroid_all(_expression_vectors(), [1, 1, 3, 1], 2),
                 "family id n_families + 1", ERR_INVALID_INPUT)


def test_group_centroid_all_rejects_no_families():
    assert_error(lambda: group_centroid_all(_expression_vectors(), np.zeros(N_GENES, dtype=np.int32), 0),
                 "no families", ERR_EMPTY_INPUT)


def test_group_centroid_all_rejects_no_axes():
    assert_error(lambda: group_centroid_all(np.empty((0, N_GENES)), np.ones(N_GENES, dtype=np.int32), 1),
                 "no axes", ERR_EMPTY_INPUT)


def test_group_centroid_all_rejects_nan():
    expression_vectors = _expression_vectors()
    expression_vectors[1, 0] = np.nan
    assert_error(lambda: group_centroid_all(expression_vectors, np.ones(N_GENES, dtype=np.int32), 1),
                 "a NaN expression", ERR_NAN_INF)


def test_group_centroid_all_rejects_a_long_gene_to_family():
    # the binding's own extent check, before the library is called
    assert_error(lambda: group_centroid_all(_expression_vectors(), np.ones(N_GENES + 1, dtype=np.int32), 1),
                 "gene_to_family longer than n_genes")


if __name__ == '__main__':
    run_all_tests(globals().values())
