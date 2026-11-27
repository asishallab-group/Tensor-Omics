"""
The `tox_data_tools` suite: the four readers and filters this module exports
(`read_gene_ids_from_tsv_file`, `read_expression_vectors_tsv`, `read_orthofinder_file`,
`get_unassigned_mask`) can be called from Python, return the documented types and shapes, and
raise the documented errors; the binding's own refusals name the argument. Whether the values
are right is the Fortran suite's job (test/mod_test_tox_data_tools/), as the coding guide's
"Where to Test What" asks; the one value each reader is checked for here is there to show that
the binding hands over the column, header and value-column arguments as they were meant.

Every file is written into a temporary directory that the case removes again.

Not called here: `read_gene_ids_from_tsv_file` and `read_expression_vectors_tsv` on a missing
file. Both close a Fortran unit they never opened, which can close the process's stderr or abort
it (see the batch's report).
"""

import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (read_gene_ids_from_tsv_file, read_expression_vectors_tsv, read_orthofinder_file,
                          get_unassigned_mask)
from test_helpers import run_all_tests, assert_error, assert_refused_by_name
from tensor_omics.error_handling import ERR_FILE_OPEN, ERR_INVALID_INPUT, ERR_SIZE_MISMATCH

#: the gene_to_fam of a gene in no family
UNASSIGNED = 0


def _write(directory, name, content):
    """Writes `content` to `name` in `directory`, byte for byte, and returns the path."""
    path = os.path.join(directory, name)
    with open(path, "w", newline="") as file:
        file.write(content)
    return path


def _expression_file(directory):
    """Three genes, a header, and the values of g2 and g1 swapped against `gene_ids` order."""
    return _write(directory, "expression.tsv",
                  "id\ta\tb\n"
                  "g2\t1.5\t2.5\n"
                  "g1\t3.5\t4.5\n"
                  "g3\t5.5\t6.5\n")


def _orthogroups_file(directory):
    """Two families: OG1 lists g3 and g1, OG2 lists g2; g4 is in none."""
    return _write(directory, "Orthogroups.tsv",
                  "Orthogroup\thuman\n"
                  "OG1\tg3, g1\n"
                  "OG2\tg2\n")


# =====================
# read_gene_ids_from_tsv_file
# =====================

def test_read_gene_ids_from_tsv_file():
    with tempfile.TemporaryDirectory() as directory:
        tsv = _write(directory, "genes.tsv", "h\tid\nx\tg1\ny\tg22\n")
        result = read_gene_ids_from_tsv_file(tsv, 8, 2, 1, 2)
    assert isinstance(result, list), f"expected a list, got {type(result).__name__}"
    assert all(isinstance(gene, str) for gene in result), f"expected strings, got {result!r}"
    # the one value: column 2 after one header line, so gene_col and n_header_rows arrive
    assert result == ["g1", "g22"], f"got {result!r}"


def test_read_gene_ids_from_tsv_file_reads_a_file_of_no_gene():
    # moved from mod_test_tox_data.py's test_an_array_of_no_strings (#221 B3)
    with tempfile.TemporaryDirectory() as directory:
        tsv = _write(directory, "no_genes.tsv", "id\tv1\n")
        assert read_gene_ids_from_tsv_file(tsv, 8, 0, 1, 1) == []


def test_read_gene_ids_from_tsv_file_refuses_more_rows_than_ids():
    with tempfile.TemporaryDirectory() as directory:
        tsv = _write(directory, "genes.tsv", "id\ng1\ng2\ng3\n")
        assert_error(lambda: read_gene_ids_from_tsv_file(tsv, 8, 2, 1, 1),
                     "three rows for two ids", ERR_SIZE_MISMATCH)


def test_read_gene_ids_from_tsv_file_refuses_a_negative_length_or_count_by_name():
    # moved from mod_test_tox_data.py: both size the strings the binding allocates before the
    # file is even opened (#215)
    assert_refused_by_name(lambda: read_gene_ids_from_tsv_file("absent.test.tsv", -1, 3, 1, 1),
                           "gene_ids_strlen = -1", "gene_ids_strlen")
    assert_refused_by_name(lambda: read_gene_ids_from_tsv_file("absent.test.tsv", 8, -1, 1, 1),
                           "n_gene_ids_elements = -1", "n_gene_ids_elements")


def test_read_gene_ids_from_tsv_file_refuses_wrong_types_by_name():
    assert_refused_by_name(lambda: read_gene_ids_from_tsv_file(42, 8, 3, 1, 1), "a number for a filename",
                           "filename")
    assert_refused_by_name(lambda: read_gene_ids_from_tsv_file("absent.test.tsv", 8, 3, "one", 1),
                           "a str for n_header_rows", "n_header_rows")
    assert_refused_by_name(lambda: read_gene_ids_from_tsv_file("absent.test.tsv", 8, 3, 1, 1.5),
                           "a fraction for gene_col", "gene_col")


# =====================
# read_expression_vectors_tsv
# =====================

def test_read_expression_vectors_tsv():
    expression = np.zeros((2, 3), order='F')
    with tempfile.TemporaryDirectory() as directory:
        result = read_expression_vectors_tsv([_expression_file(directory)], ["g1", "g2", "g3"], expression,
                                             1, 1, np.array([3, 2], dtype=np.int32), 1)
    assert result is None, f"expected None (expression_vectors is filled in place), got {result!r}"
    assert expression.shape == (2, 3) and expression.dtype == np.float64, f"got {expression.shape}, {expression.dtype}"
    # the one value: g1 is the file's second row, and value_cols [3, 2] puts column 3 first
    assert np.array_equal(expression[:, 0], [4.5, 3.5]), f"got {expression[:, 0]}"


def test_read_expression_vectors_tsv_reads_no_file():
    # moved from mod_test_tox_data.py's test_an_array_of_no_strings (#221 B3): no file and no gene
    assert read_expression_vectors_tsv([], [], np.ones((1, 0), order='F'), 1, 1, np.array([2], np.int32), 1) is None


def test_read_expression_vectors_tsv_takes_a_delimiter():
    expression = np.zeros((1, 1), order='F')
    with tempfile.TemporaryDirectory() as directory:
        csv = _write(directory, "expression.csv", "id,a\ng1,2.5\n")
        read_expression_vectors_tsv([csv], ["g1"], expression, 1, 1, np.array([2], np.int32), 1, delimiter=",")
    assert expression[0, 0] == 2.5, f"got {expression[0, 0]}"


def test_read_expression_vectors_tsv_refuses_an_array_it_cannot_fill_in_place_by_name():
    cols = np.array([2], np.int32)
    # 2 x 2: an array of no element is column-major and row-major at once
    assert_refused_by_name(lambda: read_expression_vectors_tsv([], ["g1", "g2"], np.ones((2, 2), order='C'),
                                                               1, 1, cols, 1),
                           "a row-major array", "expression_vectors")
    assert_refused_by_name(lambda: read_expression_vectors_tsv([], [], np.ones((1, 0), np.float32, order='F'),
                                                               1, 1, cols, 1),
                           "a float32 array", "expression_vectors")
    assert_refused_by_name(lambda: read_expression_vectors_tsv([], [], np.ones(0), 1, 1, cols, 1),
                           "one dimension", "expression_vectors")
    assert_refused_by_name(lambda: read_expression_vectors_tsv([], [], [[1.0]], 1, 1, cols, 1),
                           "a list", "expression_vectors")


def test_read_expression_vectors_tsv_refuses_wrong_types_by_name():
    expression = np.ones((1, 0), order='F')
    assert_refused_by_name(lambda: read_expression_vectors_tsv("one.tsv", [], expression, 1, 1, [2], 1),
                           "a single str for file_list", "file_list")
    assert_refused_by_name(lambda: read_expression_vectors_tsv([], [], expression, 1, 1, [2.5], 1),
                           "a fraction in value_cols", "value_cols")
    assert_refused_by_name(lambda: read_expression_vectors_tsv([], [], expression, 1, 1, [2], 1, delimiter=", "),
                           "a delimiter of two characters", "delimiter")


def test_read_expression_vectors_tsv_refuses_no_value_column():
    assert_error(lambda: read_expression_vectors_tsv([], [], np.ones((1, 0), order='F'), 1, 1,
                                                     np.array([], np.int32), 1),
                 "no value column", ERR_INVALID_INPUT)


def test_read_expression_vectors_tsv_refuses_a_column_per_gene_too_few():
    assert_error(lambda: read_expression_vectors_tsv([], ["g1", "g2"], np.ones((1, 1), order='F'), 1, 1,
                                                     np.array([2], np.int32), 1),
                 "one column for two genes", ERR_SIZE_MISMATCH)


def test_read_expression_vectors_tsv_refuses_nan():
    expression = np.zeros((1, 1), order='F')
    with tempfile.TemporaryDirectory() as directory:
        tsv = _write(directory, "expression.tsv", "id\ta\ng1\tNaN\n")
        assert_error(lambda: read_expression_vectors_tsv([tsv], ["g1"], expression, 1, 1, np.array([2], np.int32), 1),
                     "a NaN in the file", ERR_INVALID_INPUT)


# =====================
# read_orthofinder_file
# =====================

def test_read_orthofinder_file():
    with tempfile.TemporaryDirectory() as directory:
        result = read_orthofinder_file(_orthogroups_file(directory), ["g1", "g2", "g3", "g4"], 8, 2, 4)
    assert isinstance(result, dict) and set(result) == {"family_ids", "gene_to_fam"}, f"got {result!r}"
    family_ids, gene_to_fam = result["family_ids"], result["gene_to_fam"]
    assert isinstance(family_ids, list) and all(isinstance(family, str) for family in family_ids), \
        f"expected a list of strings, got {family_ids!r}"
    assert isinstance(gene_to_fam, np.ndarray) and gene_to_fam.dtype == np.int32, \
        f"expected an int32 array, got {type(gene_to_fam).__name__} {getattr(gene_to_fam, 'dtype', '')}"
    assert gene_to_fam.shape == (4,), f"got shape {gene_to_fam.shape}"
    assert not gene_to_fam.flags.writeable, "a result is a value, it should be read-only"
    # the one value: family indices are 1-based, and a gene in no family gets the sentinel
    assert family_ids == ["OG1", "OG2"], f"got {family_ids!r}"
    assert np.array_equal(gene_to_fam, [1, 2, 1, UNASSIGNED]), f"got {gene_to_fam}"


def test_read_orthofinder_file_reads_a_file_of_no_family():
    # moved from mod_test_tox_data.py's test_an_array_of_no_strings (#221 B3)
    with tempfile.TemporaryDirectory() as directory:
        tsv = _write(directory, "no_families.tsv", "Orthogroup\tsp1\n")
        families = read_orthofinder_file(tsv, [], 8, 0, 0)
    assert families["family_ids"] == [] and families["gene_to_fam"].size == 0, families


def test_read_orthofinder_file_refuses_a_missing_file():
    with tempfile.TemporaryDirectory() as directory:
        missing = os.path.join(directory, "missing.tsv")
        assert_error(lambda: read_orthofinder_file(missing, ["g1"], 8, 1, 1), "a missing file", ERR_FILE_OPEN)


def test_read_orthofinder_file_refuses_a_negative_length_or_count_by_name():
    for argument, sizes in (("family_ids_strlen", (-1, 1, 1)), ("n_family_ids_elements", (8, -1, 1)),
                            ("n_gene_to_fam_elements", (8, 1, -1))):
        assert_refused_by_name(lambda: read_orthofinder_file("absent.test.tsv", ["g1"], *sizes),
                               f"{argument} = -1", argument)


def test_read_orthofinder_file_refuses_wrong_types_by_name():
    assert_refused_by_name(lambda: read_orthofinder_file("absent.test.tsv", "g1", 8, 1, 1),
                           "a single str for gene_ids", "gene_ids")
    assert_refused_by_name(lambda: read_orthofinder_file(None, ["g1"], 8, 1, 1), "None for a filename", "filename")


# =====================
# get_unassigned_mask
# =====================

def test_get_unassigned_mask():
    result = get_unassigned_mask(np.array([2, UNASSIGNED, 1], dtype=np.int32))
    assert isinstance(result, dict) and set(result) == {"mask", "n_genes_kept"}, f"got {result!r}"
    mask, n_genes_kept = result["mask"], result["n_genes_kept"]
    assert isinstance(mask, np.ndarray) and mask.dtype == np.bool_, f"expected a bool array, got {mask!r}"
    assert mask.shape == (3,), f"got shape {mask.shape}"
    assert not mask.flags.writeable, "a result is a value, it should be read-only"
    assert isinstance(n_genes_kept, int), f"expected an int, got {type(n_genes_kept).__name__}"
    # the one value: the sentinel's gene is the one dropped
    assert mask.tolist() == [True, False, True] and n_genes_kept == 2, f"got {mask}, {n_genes_kept}"


def test_get_unassigned_mask_of_no_gene():
    result = get_unassigned_mask(np.array([], dtype=np.int32))
    assert result["mask"].shape == (0,) and result["n_genes_kept"] == 0, f"got {result!r}"


def test_get_unassigned_mask_refuses_wrong_types_by_name():
    assert_refused_by_name(lambda: get_unassigned_mask([1.5, 2.0]), "a fraction", "gene_to_fam")
    assert_refused_by_name(lambda: get_unassigned_mask(np.array([[1, 2]], dtype=np.int32)), "two dimensions",
                           "gene_to_fam")
    assert_refused_by_name(lambda: get_unassigned_mask(["OG1"]), "strings", "gene_to_fam")


if __name__ == "__main__":
    run_all_tests(globals().values())
