#!/usr/bin/env python3
"""The input contract of the Python binding: every case of the table, accepted and refused.

The table is the one both bindings follow (helper/codegen/design/language-layers.md, "what an
argument accepts"); r/test/mod_test_input_contract.R is the same list in R. A refusal is the
binding's own: a ToxArgumentError naming the argument, with no tox code (#221).
"""

import os
import pathlib
import re
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    calc_quantile,
    compute_family_scaling,
    deserialize_char_helper,
    deserialize_int_helper,
    euclidean_distance,
    group_centroid_orthologs,
    linkage_clustering,
    read_expression_vectors_tsv,
    serialize_char_helper,
    serialize_int_helper,
    validate_expression_data,
    validate_string_array_uniqueness,
)
from tensor_omics.error_handling import ERR_INVALID_INPUT, ToxArgumentError
from test_helpers import run_all_tests, assert_error

EXPR = np.array([[1.0, 2.0, 3.0, 4.0]])          # one axis, four genes
FAMILY = [1, 1, 1, 1]
MASK = [True, False, True, False]
SCRATCH = tempfile.mkdtemp(prefix="tox_contract_")


def refused(func, argument, pattern):
    """`func` raises the binding's own refusal of `argument`, its message matching `pattern`."""
    try:
        func()
    except ToxArgumentError as error:
        assert error.code is None, f"a refusal carries no code, but got {error.code}"
        assert error.argument == argument, f"expected '{argument}' named, got {error.argument!r}: {error}"
        assert re.search(pattern, str(error)), f"{str(error)!r} does not match {pattern!r}"
        return
    raise AssertionError(f"'{argument}' was not refused")


def centroid(gene_to_family=FAMILY, n_families=1, mask=MASK, expr=EXPR):
    return group_centroid_orthologs(expr, gene_to_family, n_families, mask)


def path(name):
    return os.path.join(SCRATCH, name)


# =====================================================================================
# integer(int32)
# =====================================================================================

def test_integer_accepts_python_and_numpy_integers_and_0d():
    for n_families in (1, np.int64(1), np.int32(1), np.uint8(1), np.array(1)):
        assert np.allclose(centroid(n_families=n_families), [[2.0]])
    for gene_to_family in ([1, 1, 1, 1], np.int64(FAMILY), np.uint16(FAMILY), np.int32(FAMILY)):
        assert np.allclose(centroid(gene_to_family), [[2.0]])


def test_integer_accepts_whole_floats_and_negative_zero():
    assert np.allclose(centroid([1.0, 1.0, 1.0, 1.0], n_families=1.0), [[2.0]])
    assert np.allclose(centroid(n_families=np.float64(1.0)), [[2.0]])
    # -0.0 is 0: a gene in no family, which the mask leaves out anyway
    assert np.allclose(centroid([1.0, -0.0, 1.0, 1.0]), [[2.0]])


def test_integer_refuses_bools():
    refused(lambda: centroid(n_families=True), "n_families", "not True \\(bool\\)")
    refused(lambda: centroid([True, True, True, True]), "gene_to_family", "dtype bool")


def test_integer_refuses_fractions_nan_and_infinities():
    # these were truncated or wrapped before (#221 N1, G-S4)
    refused(lambda: centroid([1.7, 1.2, 1.9, 1.0]), "gene_to_family", "element 0 is 1.7")
    refused(lambda: centroid([1.0, np.nan, 1.0, 1.0]), "gene_to_family", "element 1 is nan")
    refused(lambda: centroid([1.0, 1.0, np.inf, 1.0]), "gene_to_family", "element 2 is inf")
    refused(lambda: centroid(n_families=1.5), "n_families", "but it is 1.5")
    refused(lambda: centroid(n_families=float("nan")), "n_families", "but it is nan")


def test_integer_refuses_out_of_range_before_narrowing():
    # wrapped modulo 2**32 before, which also slipped past DM_MAX (#221 N2, G-S4)
    refused(lambda: centroid(np.int64([2**32 + 1] * 4)), "gene_to_family", "4294967297")
    refused(lambda: centroid(n_families=2**32 + 1), "n_families", "4294967297")
    refused(lambda: centroid(n_families=2**31), "n_families", "2147483648")
    refused(lambda: centroid(n_families=2**70), "n_families", "not 1180591620717411303424")


def test_integer_refuses_strings_objects_datetimes_and_complex():
    refused(lambda: centroid(["1", "1", "1", "1"]), "gene_to_family", "dtype <U1")
    refused(lambda: centroid([1, None, 1, 1]), "gene_to_family", "dtype object")
    refused(lambda: centroid(np.array(["2020-01-01"] * 4, dtype="datetime64[D]")),
            "gene_to_family", "datetime64")
    refused(lambda: centroid(np.ones(4, dtype=complex)), "gene_to_family", "complex128")
    refused(lambda: centroid(n_families="1"), "n_families", "not '1' \\(str\\)")


def test_integer_refuses_a_masked_array_with_masked_elements():
    masked = np.ma.masked_array(FAMILY, mask=[0, 0, 1, 0])
    refused(lambda: centroid(masked), "gene_to_family", "masked elements")
    # one with nothing masked is just its data
    assert np.allclose(centroid(np.ma.masked_array(FAMILY, mask=False)), [[2.0]])


def test_a_uint32_bit_mask_word_is_refused_with_the_view_hint():
    refused(lambda: centroid(np.uint32([2**31, 1, 1, 1])), "gene_to_family", "view\\(np.int32\\)")


# =====================================================================================
# real(real64)
# =====================================================================================

def test_real_accepts_ints_floats_and_their_numpy_kinds():
    for vec1 in ([3, 4], [3.0, 4.0], np.int32([3, 4]), np.uint8([3, 4]), np.float32([3, 4])):
        assert np.isclose(euclidean_distance(vec1, [0, 0]), 5.0)
    for level in (0.5, np.float64(0.5), np.float32(0.5), np.array(0.5)):
        assert np.isclose(calc_quantile([1.0, 2.0, 3.0], level), 2.0)
    assert np.isclose(calc_quantile([1.0, 2.0, 3.0], 1), 3.0)


def test_real_accepts_an_integer_only_where_float64_holds_it_exactly():
    # 2**53 + 1 reached Fortran as 2**53, and a long double was rounded (#221 WP8)
    big = 2**53
    assert euclidean_distance([big, 0], [0, 0]) == big
    assert euclidean_distance(np.int64([big, 0]), [0, 0]) == big
    assert euclidean_distance(np.uint64([2**63, 0]), [0, 0]) == 2**63
    assert calc_quantile([1.0, 2.0, 3.0], np.int64(1)) == 3.0
    for inexact in ([big + 1, 0], np.int64([big + 1, 0]), np.uint64([big + 1, 0])):
        refused(lambda: euclidean_distance(inexact, [0, 0]), "vec1",
                "that a float64 holds exactly, but element 0 is 9007199254740993, "
                r"which it would round to 9007199254740992\.0")
    refused(lambda: euclidean_distance(np.uint64([2**64 - 1, 0]), [0, 0]), "vec1",
            "element 0 is 18446744073709551615")
    for level in (big + 1, np.int64(big + 1), np.uint64(big + 1)):
        refused(lambda: calc_quantile([1.0, 2.0], level), "level", "it is 9007199254740993")
    if np.finfo(np.longdouble).nmant > 52:  # an 80- or 128-bit long double
        inexact = np.longdouble(1) + np.finfo(np.longdouble).eps
        refused(lambda: euclidean_distance(np.array([inexact, 0]), [0, 0]), "vec1", "element 0 is 1.0000")
        assert euclidean_distance(np.longdouble([3, 4]), [0, 0]) == 5.0


def test_real_leaves_nan_to_the_procedure():
    # DM_ALLOW_NAN decides, so the binding passes NaN on and the Fortran refuses it with a code
    error = None
    try:
        euclidean_distance([np.nan, 4.0], [0.0, 0.0])
    except Exception as caught:
        error = caught
    assert error is not None and not isinstance(error, ToxArgumentError) and error.code is not None


def test_real_refuses_bools_strings_complex_datetimes_and_objects():
    refused(lambda: euclidean_distance([True, False], [0, 0]), "vec1", "dtype bool")
    # numpy would parse these as numbers
    refused(lambda: euclidean_distance(["1.5", "2"], [0, 0]), "vec1", "dtype <U3")
    refused(lambda: euclidean_distance([1 + 0j, 2], [0, 0]), "vec1", "complex")
    # numpy would count the days since 1970
    refused(lambda: euclidean_distance(np.array(["2020-01-01", "2020-01-02"], dtype="datetime64[D]"),
                                       [0, 0]), "vec1", "datetime64")
    refused(lambda: euclidean_distance([1.0, None], [0, 0]), "vec1", "dtype object")
    refused(lambda: calc_quantile([1.0, 2.0], True), "level", "bool")
    refused(lambda: calc_quantile([1.0, 2.0], "0.5"), "level", "str")


def test_real_refuses_a_masked_array_with_masked_elements():
    # the masked 1e6 was read before (#221 F7)
    masked = np.ma.masked_array([1.0, 2.0, 1e6], mask=[0, 0, 1])
    refused(lambda: euclidean_distance(masked, np.zeros(3)), "vec1", "masked elements")


# =====================================================================================
# logical
# =====================================================================================

def test_logical_accepts_bools_and_numbers_that_are_0_or_1():
    for mask in (MASK, np.array(MASK), [1, 0, 1, 0], np.int32([1, 0, 1, 0]), [1.0, 0.0, 1.0, -0.0]):
        assert np.allclose(centroid(mask=mask), [[2.0]])
    # a negative value is reported exactly when the check is asked for
    negative = np.array([[-1.0, 1.0]])
    for flag in (True, np.True_, 1, np.array(True), 1.0):
        assert_error(lambda: validate_expression_data(negative, flag), f"checked for {flag!r}",
                     ERR_INVALID_INPUT)
    for flag in (False, np.False_, 0, 0.0):
        assert validate_expression_data(negative, flag) is None


def test_logical_refuses_other_numbers_nan_none_and_strings():
    # these were all True before, and None False (#221 item 5, N5)
    refused(lambda: centroid(mask=[1, 2, 3, 4]), "ortholog_set", "element 1 is 2")
    refused(lambda: centroid(mask=[0.5, 0, 0, 0]), "ortholog_set", "element 0 is 0.5")
    refused(lambda: centroid(mask=[1.0, np.nan, 0.0, 0.0]), "ortholog_set", "element 1 is nan")
    refused(lambda: centroid(mask=[True, None, True, False]), "ortholog_set", "dtype object")
    refused(lambda: centroid(mask=["no", "", "", ""]), "ortholog_set", "dtype <U2")
    refused(lambda: validate_expression_data(EXPR, 2), "check_non_negative", "but it is 2")
    refused(lambda: validate_expression_data(EXPR, "False"), "check_non_negative", "str")
    refused(lambda: validate_expression_data(EXPR, None), "check_non_negative", "NoneType")


# =====================================================================================
# a string
# =====================================================================================

def test_a_string_accepts_str_and_a_path():
    serialize_char_helper(["x"], path("str.bin"))
    serialize_char_helper(["y"], pathlib.Path(SCRATCH) / "path.bin")
    assert deserialize_char_helper(pathlib.Path(SCRATCH) / "path.bin").tolist() == ["y"]


def test_a_string_refuses_bytes_none_sequences_and_numbers():
    # each wrote a file of its str() before (#221 N3, G-S5)
    refused(lambda: serialize_char_helper(["x"], path("b.bin").encode()), "filename", "bytes")
    refused(lambda: serialize_char_helper(["x"], None), "filename", "NoneType")
    refused(lambda: serialize_char_helper(["x"], [path("f1.bin"), path("f2.bin")]), "filename", "list")
    refused(lambda: serialize_char_helper(["x"], 42), "filename", "int")
    assert sorted(os.listdir(SCRATCH)) == ["path.bin", "str.bin"]


def test_a_string_refuses_a_nul_and_a_lone_surrogate():
    # a NUL wrote a file named up to it before; a surrogate raised a bare UnicodeEncodeError
    refused(lambda: serialize_char_helper(["x"], path("nul\x00evil.bin")), "filename", "NUL")
    refused(lambda: serialize_char_helper(["x"], path("\udcff.bin")), "filename", "surrogates")
    assert not os.path.exists(path("nul"))


def test_a_string_longer_than_its_declared_length_is_refused_not_cut():
    # `delimiter='::'` reached Fortran as ':' before
    refused(lambda: read_expression_vectors_tsv([path("none.tsv")], ["g"], np.zeros((1, 1), order="F"),
                                                0, 1, [2], 1, delimiter="::"),
            "delimiter", "2 bytes long in UTF-8, but it is declared 1 long")
    refused(lambda: read_expression_vectors_tsv([path("none.tsv")], ["g"], np.zeros((1, 1), order="F"),
                                                0, 1, [2], 1, delimiter="\u00e9"),
            "delimiter", "2 bytes long")


# =====================================================================================
# an array of strings
# =====================================================================================

def test_strings_accept_a_list_tuple_ndarray_series_and_categorical():
    for strings in (["a", "b"], ("a", "b"), np.array(["a", "b"]), np.array(["a", "b"], dtype=object)):
        assert validate_string_array_uniqueness(strings) is None
    assert validate_string_array_uniqueness(pd.Series(["a", "b"])) is None
    assert validate_string_array_uniqueness(pd.Categorical(["a", "b"])) is None


def test_strings_refuse_a_bare_str_and_any_element_that_is_not_a_str():
    # 'abc' was taken as ['a', 'b', 'c'], and None and 1 as "None" and "1" (#221 N4)
    refused(lambda: validate_string_array_uniqueness("abc"), "str_arr", "single str; wrap it in a list")
    refused(lambda: validate_string_array_uniqueness([None, "a"]), "str_arr", "element 0 is NoneType")
    refused(lambda: validate_string_array_uniqueness([1, "1"]), "str_arr", "element 0 is int")
    refused(lambda: validate_string_array_uniqueness(42), "str_arr", "int")
    refused(lambda: validate_string_array_uniqueness({"a", "b"}), "str_arr", "set, which has no order")
    refused(lambda: validate_string_array_uniqueness(np.array([["a", "b"]])), "str_arr", "1 dimension")


def test_strings_refuse_bytes_with_the_decode_hint():
    refused(lambda: validate_string_array_uniqueness(np.array([b"a", b"b"])), "str_arr", "astype\\(str\\)")
    refused(lambda: validate_string_array_uniqueness([b"a"]), "str_arr", "bytes")


def test_strings_refuse_a_nul_and_a_lone_surrogate_naming_the_element():
    refused(lambda: validate_string_array_uniqueness(["a", "b\x00c"]), "str_arr", "NUL character \\(element 1\\)")
    refused(lambda: validate_string_array_uniqueness(["a", "\udcff"]), "str_arr", "Unicode \\(element 1\\)")


def test_a_serialized_array_of_strings_takes_any_rank_and_refuses_non_strings():
    serialize_char_helper([["a", "b"], ["c", "d"]], path("grid.bin"))
    assert deserialize_char_helper(path("grid.bin")).tolist() == [["a", "b"], ["c", "d"]]
    refused(lambda: serialize_char_helper([1, 2], path("ints.bin")), "arr", "element 0 is int")


# =====================================================================================
# trailing blanks: Fortran semantics, accepted and dropped
# =====================================================================================

def test_trailing_blanks_are_accepted_and_not_significant():
    # to Fortran 'a ' is 'a', so they are duplicates
    assert_error(lambda: validate_string_array_uniqueness(["a", "a "]), "duplicates", ERR_INVALID_INPUT)
    # and never returned; leading blanks survive
    serialize_char_helper(["ab  ", " c"], path("blanks.bin"))
    assert deserialize_char_helper(path("blanks.bin")).tolist() == ["ab", " c"]


def test_an_empty_array_of_strings_comes_back_as_str():
    serialize_char_helper(np.array([], dtype=str), path("empty.bin"))
    back = deserialize_char_helper(path("empty.bin"))
    assert back.dtype.kind == "U" and back.size == 0


# =====================================================================================
# modes
# =====================================================================================

DISTANCES = np.asfortranarray([[0.0, 1.0, 2.0], [1.0, 0.0, 1.0], [2.0, 1.0, 0.0]])


def test_a_mode_matches_whole_with_the_case_rule_kept():
    assert linkage_clustering(DISTANCES.copy(order="F"), "ward")["heights"].size == 2
    assert linkage_clustering(DISTANCES.copy(order="F"), "WARD")["heights"].size == 2
    # trailing blanks are Fortran semantics, as for every string; a leading one is part of it
    assert linkage_clustering(DISTANCES.copy(order="F"), "ward  ")["heights"].size == 2
    refused(lambda: linkage_clustering(DISTANCES.copy(order="F"), " ward"), "method", "not ' ward'")


def test_a_mode_is_never_truncated_and_the_message_lists_the_modes():
    # 'weightedXYZ...' ran the weighted method before (#221 G-I1); the refusal named nothing (G-I2)
    refused(lambda: linkage_clustering(DISTANCES.copy(order="F"), "weightedXYZ_not_a_method"),
            "method", "must be one of 'average', .*'weighted'.*not 'weightedXYZ_not_a_method'")
    refused(lambda: compute_family_scaling(1, [1.0, 2.0], [1, 1], mode="not_a_mode"),
            "mode", "one of 'plain', 'robust'")
    refused(lambda: linkage_clustering(DISTANCES.copy(order="F"), 5), "method", "not int")


# =====================================================================================
# rank
# =====================================================================================

def test_rank_must_be_the_declared_one_and_nothing_is_flattened():
    refused(lambda: euclidean_distance([[3.0, 4.0]], [0.0, 0.0]), "vec1", "must have 1 dimension, but has 2")
    refused(lambda: centroid(expr=[1.0, 2.0, 3.0, 4.0]), "expression_vectors", "2 dimensions, but has 1")
    refused(lambda: centroid(n_families=[1]), "n_families", "single value, not an array of shape \\(1,\\)")


def test_a_dataframe_is_a_matrix():
    assert np.allclose(centroid(expr=pd.DataFrame(EXPR), gene_to_family=pd.Series(FAMILY),
                                mask=pd.Series(MASK)), [[2.0]])


def test_a_serialized_array_takes_any_rank_and_checks_its_values():
    serialize_int_helper(np.array([[1.0, 2.0], [3.0, 4.0]]), path("ints2d.bin"))
    back = deserialize_int_helper(path("ints2d.bin"))
    assert back.dtype == np.int32 and back.tolist() == [[1, 2], [3, 4]]
    refused(lambda: serialize_int_helper([1.5, 2.0], path("frac.bin")), "arr", "element 0 is 1.5")


if __name__ == "__main__":
    run_all_tests(globals().values())
