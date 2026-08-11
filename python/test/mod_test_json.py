"""
Python test suite for serialize_tox_data_as_flyer_json (tox_flyer_json), mirroring
test/mod_test_json.f90
"""

import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import serialize_tox_data_as_flyer_json
from test_helpers import run_all_tests


def test_flyer_serialization():
    filename = "test_flyer_serialization.json"

    tissues = ["Adipose", "Thyroid"]
    family_ids = ["EMPTYFAM1", "OG0000000", "OG0000001", "EMPTYFAM2"]
    gene_ids = [
        "NP_001243379.1",
        "UNASSIGNED.1",
        "XP_038381480.1",
        "NP_001243386.1",
        "XP_038421312.1",
    ]
    gene_species = [
        "Canis_lupus_protein.1",
        "Whatever_protein",
        "Canis_lupus_protein.2",
        "Canis_lupus_protein.3",
        "Canis_lupus_protein.4",
    ]
    gene_types = ["ortholog", "ortholog", "paralog", "ortholog", "paralog"]

    gene_to_fam = np.array([2, 0, 3, 2, 3], dtype=np.int32)
    sorted_gene_to_fam_perm = np.array([2, 1, 4, 3, 5], dtype=np.int32)
    gene_outliers = np.array([True, False, False, False, True], dtype=np.bool_)

    genes = np.asfortranarray([
        [0.0541530138142222, 3.1415926535897932, 1.115620102135, 0.0060563875402208, 0.0060563875402208],
        [0.0041979991981664, 2.7182818284590452, 0.308001439923219, 0.147377684857407, 0.0041979991981664],
    ], dtype=np.float64)

    centroids = np.asfortranarray([
        [0.0, 0.03010470067722, 0.56083824483761, 0.0],
        [0.0, 0.07578784202779, 0.15609971956069, 0.0],
    ], dtype=np.float64)

    serialize_tox_data_as_flyer_json(
        filename, tissues, family_ids, centroids, gene_ids, genes,
        gene_to_fam, sorted_gene_to_fam_perm, gene_outliers, gene_species, gene_types,
    )

    assert os.path.exists(filename), "Output JSON file does not exist"

    with open(filename, "r") as f:
        text = f.read().strip()

    expected_fragment = (
        '{'
            '"tissues":["Adipose","Thyroid"],'
            '"families":['
                '{"family":"EMPTYFAM1","gene_indices":[],"centroid":[0.0000000000000000E+000,0.0000000000000000E+000]},'
                '{"family":"OG0000000","gene_indices":[1,4],"centroid":[3.0104700677220000E-002,7.5787842027790001E-002]},'
                '{"family":"OG0000001","gene_indices":[3,5],"centroid":[5.6083824483761002E-001,1.5609971956068999E-001]},'
                '{"family":"EMPTYFAM2","gene_indices":[],"centroid":[0.0000000000000000E+000,0.0000000000000000E+000]}'
            '],'
            '"genes":['
                '{"coordinates":[5.4153013814222203E-002,4.1979991981664000E-003],"id":"NP_001243379.1","family":"OG0000000","species":"Canis_lupus_protein.1","is_outlier":true,"type":"ortholog"},'
                '{"coordinates":[3.1415926535897931E+000,2.7182818284590451E+000],"id":"UNASSIGNED.1","family":null,"species":"Whatever_protein","is_outlier":false,"type":"ortholog"},'
                '{"coordinates":[1.1156201021350001E+000,3.0800143992321899E-001],"id":"XP_038381480.1","family":"OG0000001","species":"Canis_lupus_protein.2","is_outlier":false,"type":"paralog"},'
                '{"coordinates":[6.0563875402208003E-003,1.4737768485740699E-001],"id":"NP_001243386.1","family":"OG0000000","species":"Canis_lupus_protein.3","is_outlier":false,"type":"ortholog"},'
                '{"coordinates":[6.0563875402208003E-003,4.1979991981664000E-003],"id":"XP_038421312.1","family":"OG0000001","species":"Canis_lupus_protein.4","is_outlier":true,"type":"paralog"}'
            ']'
        '}'
    )

    assert expected_fragment == text, "JSON output does not match expected fragment"

    os.remove(filename)


if __name__ == "__main__":
    run_all_tests(globals().values())
