"""
Python test suite for tox_stc_csv (serialize_stc_points_as_csv,
serialize_stc_ensemble_overlap_as_csv, serialize_stc_super_ensembles_as_tsv), mirroring
test/mod_test_shape_truthful_clustering_csv.F90
"""

import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    serialize_stc_points_as_csv,
    serialize_stc_ensemble_overlap_as_csv,
    serialize_stc_super_ensembles_as_tsv,
)
from test_helpers import run_all_tests


def _fixture():
    """N=4, 2 ensembles: {1,2,3} (seed=1) and {2,3,4} (seed=4), overlapping on {2,3} --
    Overlap Coefficient 2/3 -- merged into one super-ensemble. Point 1 is also flagged
    low-confidence. Matches mod_test_shape_truthful_clustering_json.py's own fixture."""
    seed_selection_mask = np.array([True, False, False, True], dtype=np.bool_)
    ensemble_masks = np.zeros((4, 2), dtype=np.bool_, order='F')
    ensemble_masks[0:3, 0] = True
    ensemble_masks[1:4, 1] = True
    ensemble_low_confidence_masks = np.zeros((4, 2), dtype=np.bool_, order='F')
    ensemble_low_confidence_masks[0, 0] = True
    super_ensembles = np.asfortranarray([[1, 0], [2, 0]], dtype=np.int32)
    return seed_selection_mask, ensemble_masks, ensemble_low_confidence_masks, super_ensembles


def test_points_csv_membership():
    seed_selection_mask, ensemble_masks, ensemble_low_confidence_masks, super_ensembles = _fixture()

    filename = "test_stc_points_py.csv"
    serialize_stc_points_as_csv(filename, 1, seed_selection_mask, ensemble_masks,
                                ensemble_low_confidence_masks, super_ensembles)

    with open(filename, "r") as f:
        lines = [line.rstrip("\n") for line in f.readlines()]
    os.remove(filename)

    assert lines[0] == "row,ensembles,super_ensembles,low_confidence_ensembles,seed_of"
    assert lines[1] == '1,"1","1","1","1"'
    assert lines[2] == '2,"1,2","1","",""'
    assert lines[3] == '3,"1,2","1","",""'
    assert lines[4] == '4,"2","1","","2"'


def test_points_csv_zero_ensembles():
    seed_selection_mask = np.array([False, False], dtype=np.bool_)
    ensemble_masks = np.zeros((2, 0), dtype=np.bool_, order='F')
    ensemble_low_confidence_masks = np.zeros((2, 0), dtype=np.bool_, order='F')
    super_ensembles = np.zeros((2, 0), dtype=np.int32, order='F')

    filename = "test_stc_points_zero_py.csv"
    serialize_stc_points_as_csv(filename, 0, seed_selection_mask, ensemble_masks,
                                ensemble_low_confidence_masks, super_ensembles)

    with open(filename, "r") as f:
        lines = [line.rstrip("\n") for line in f.readlines()]
    os.remove(filename)

    assert lines[1] == '1,"","","",""'


def test_overlap_csv_pairwise_coefficient():
    _, ensemble_masks, _, _ = _fixture()

    filename = "test_stc_overlap_py.csv"
    serialize_stc_ensemble_overlap_as_csv(filename, ensemble_masks)

    with open(filename, "r") as f:
        lines = [line.rstrip("\n") for line in f.readlines()]
    os.remove(filename)

    assert lines[0] == "ensemble_a,ensemble_b,overlap_coefficient"
    assert lines[1] == "1,2,6.6666666666666663E-001"


def test_overlap_csv_no_intersecting_pairs():
    ensemble_masks = np.zeros((4, 2), dtype=np.bool_, order='F')
    ensemble_masks[0:2, 0] = True
    ensemble_masks[2:4, 1] = True

    filename = "test_stc_overlap_disjoint_py.csv"
    serialize_stc_ensemble_overlap_as_csv(filename, ensemble_masks)

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert content.strip() == "ensemble_a,ensemble_b,overlap_coefficient"


def test_super_ensembles_tsv_gene_family_format():
    super_ensembles = np.asfortranarray([[1], [2], [0]], dtype=np.int32)

    filename = "test_stc_super_ensembles_py.tsv"
    serialize_stc_super_ensembles_as_tsv(filename, super_ensembles)

    with open(filename, "r") as f:
        lines = [line.rstrip("\n") for line in f.readlines()]
    os.remove(filename)

    assert lines[0] == "1\t1,2"


def test_super_ensembles_tsv_empty_when_zero_groups():
    super_ensembles = np.zeros((2, 0), dtype=np.int32, order='F')

    filename = "test_stc_super_ensembles_empty_py.tsv"
    serialize_stc_super_ensembles_as_tsv(filename, super_ensembles)

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert content.strip() == ""


if __name__ == "__main__":
    run_all_tests(globals().values())
