"""Contract tests for save_flyer_json: the call, the file it writes, and its documented errors.

The values themselves are tested in Fortran (test/mod_test_flyer_json.F90); here the file is
compared with the same hand-written expected file and read back by an independent JSON parser.
"""
import sys
import os
import json
import math
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
# the generated package loads the shared library itself

from tensor_omics import save_flyer_json
from tensor_omics.error_handling import ToxError, ERR_NAN_INF, ERR_INVALID_INPUT, ERR_FILE_OPEN
from test_helpers import run_all_tests, assert_error

EXPECTED_FILE = os.path.join(os.path.dirname(__file__), '..', '..', 'test', 'test_files', 'flyer_expected.json')


def reference_data():
    """The data set of the Fortran suite and of the expected file, one gene per row before transposing."""
    return dict(
        expression_vectors=np.array([[1.0, 0.5, -2.0],
                                     [0.0, -0.0, 0.25],
                                     [1024.0, -0.125, 3.0],
                                     [1.5, 2.5, -3.5],
                                     [0.75, 8.0, -16.0]]).T,
        family_centroids=np.array([[0.375, 4.25, -8.125],
                                   [512.5, 0.1875, 0.5],
                                   [-0.0, 0.0, 0.0625]]).T,
        gene_to_fam=[2, 1, 2, 0, 1],
        is_outlier=[False, True, False, False, True],
        axis_labels=["liver", "brain", "heart"],
        family_ids=["F1", "F2", "F3"],
        gene_ids=["g1", "g2", "g\u00e8ne", "g4", "g5"],
        gene_species=["human", "mouse", "human", "fly", 'say "hi"'],
        gene_types=["ortholog", "paralog", "ortholog", "ortholog", "a\\b"],
    )


def assert_tox_error(call, code, argument, msg):
    """Fail unless `call` raises a ToxError with this code, naming this argument."""
    try:
        call()
    except ToxError as error:
        assert error.code == code, f"{msg}: expected code {code}, got {error.code}: {error}"
        assert error.argument == argument, f"{msg}: expected argument '{argument}', got '{error.argument}': {error}"
        assert f"'{argument}'" in str(error), f"{msg}: the message does not name '{argument}': {error}"
    else:
        raise AssertionError(f"{msg}: nothing was raised")


def test_reference_file_matches_expected():
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "flyer.json")
        save_flyer_json(filename=path, **reference_data())
        with open(path, "rb") as written, open(EXPECTED_FILE, "rb") as expected:
            assert written.read() == expected.read(), "the file differs from the expected file"


def test_file_parses_with_every_field():
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "flyer.json")
        save_flyer_json(filename=path, **reference_data())
        with open(path, encoding="utf-8") as stream:
            document = json.load(stream)

    assert document["kind"] == "tox-flyer" and document["version"] == 1, "kind and version"
    assert document["tissues"] == ["liver", "brain", "heart"], "tissues"

    families = document["families"]
    assert [family["family"] for family in families] == ["F1", "F2", "F3"], "family order"
    assert [family["gene_indices"] for family in families] == [[2, 5], [1, 3], []], "1-based gene indices"
    assert families[0]["centroid"] == [0.375, 4.25, -8.125], "centroid values"
    assert math.copysign(1.0, families[2]["centroid"][0]) == -1.0, "a negative zero keeps its sign"

    genes = document["genes"]
    assert len(genes) == 5, "one object per gene"
    for gene in genes:
        assert set(gene) == {"coordinates", "id", "family", "species", "is_outlier", "type"}, "all six keys"
    assert genes[2]["coordinates"] == [1024.0, -0.125, 3.0], "coordinates"
    assert [gene["id"] for gene in genes] == ["g1", "g2", "g\u00e8ne", "g4", "g5"], "ids, non-ASCII included"
    assert [gene["family"] for gene in genes] == ["F2", "F1", "F2", None, "F1"], "family ids, null for none"
    assert [gene["is_outlier"] for gene in genes] == [False, True, False, False, True], "outlier flags"
    assert all(type(gene["is_outlier"]) is bool for gene in genes), "outlier flags are JSON booleans"
    assert genes[4]["species"] == 'say "hi"' and genes[4]["type"] == "a\\b", "escaped strings"


def test_non_ascii_ids_round_trip():
    data = reference_data()
    data["gene_ids"] = ["基因", "\U0001F600", "été", "∅", "g5"]
    data["family_ids"] = ["族1", "F2", "F3"]
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "flyer.json")
        save_flyer_json(filename=path, **data)
        with open(path, encoding="utf-8") as stream:
            document = json.load(stream)
    assert [gene["id"] for gene in document["genes"]] == data["gene_ids"], "non-ASCII ids read back"
    assert document["genes"][1]["family"] == "族1", "non-ASCII family id read back"


def test_documented_errors_name_the_argument():
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "flyer.json")

        data = reference_data()
        data["expression_vectors"][1, 2] = np.nan
        assert_tox_error(lambda: save_flyer_json(filename=path, **data), ERR_NAN_INF, "expression_vectors", "NaN")

        data = reference_data()
        data["family_centroids"][0, 1] = np.inf
        assert_tox_error(lambda: save_flyer_json(filename=path, **data), ERR_NAN_INF, "family_centroids", "Inf")

        data = reference_data()
        data["gene_to_fam"] = [2, 1, 4, 0, 1]
        assert_tox_error(lambda: save_flyer_json(filename=path, **data), ERR_INVALID_INPUT, "gene_to_fam",
                         "family index above the number of families")

        data = reference_data()
        data["axis_labels"] = ["liver", "brain", "liver"]
        assert_tox_error(lambda: save_flyer_json(filename=path, **data), ERR_INVALID_INPUT, "axis_labels",
                         "repeated axis label")

        data = reference_data()
        data["gene_ids"] = ["g1", "g2", "g1", "g4", "g5"]
        assert_tox_error(lambda: save_flyer_json(filename=path, **data), ERR_INVALID_INPUT, "gene_ids",
                         "repeated gene id")

        data = reference_data()
        data["family_ids"] = ["F1", "", "F3"]
        assert_tox_error(lambda: save_flyer_json(filename=path, **data), ERR_INVALID_INPUT, "family_ids",
                         "empty family id")

        assert not os.path.exists(path), "a refused call creates no file"


def test_existing_file_is_refused():
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "flyer.json")
        with open(path, "w") as stream:
            stream.write("keep me")
        assert_tox_error(lambda: save_flyer_json(filename=path, **reference_data()), ERR_FILE_OPEN, "filename",
                         "existing file")
        with open(path) as stream:
            assert stream.read() == "keep me", "the existing file is untouched"


def test_shape_mismatch_is_refused_before_the_call():
    data = reference_data()
    data["gene_ids"] = ["g1", "g2", "g3", "g4"]
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "flyer.json")
        assert_error(lambda: save_flyer_json(filename=path, **data), "gene_ids shorter than the genes")
        assert not os.path.exists(path), "a refused call creates no file"


if __name__ == "__main__":
    run_all_tests(list(globals().values()))
