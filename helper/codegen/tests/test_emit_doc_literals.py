"""Rendering the Fortran literals in a doc for the language that will read it.

The author writes on a Fortran kernel and the `DM_` range macros quote that Fortran back
verbatim, so a bound arrives as `0_int32` and a default as `.false.`. Carried into a Python
docstring unchanged, the second one contradicts the signature three lines above it.
"""

import pytest

from codegen.emit.doc_literals import render


class TestKindSuffixes:
    @pytest.mark.parametrize("fortran, expected", [
        ("The minimum valid value is `0_int32`.", "The minimum valid value is `0`."),
        ("The maximum valid value is `1.0_real64`.", "The maximum valid value is `1.0`."),
        ("of shape (max(0_int32, n - 1),)", "of shape (max(0, n - 1),)"),
        ("a tolerance of `1.0e-9_real64`", "a tolerance of `1.0e-9`"),
        ("`2_int64` and `3.5_real32`", "`2` and `3.5`"),
    ])
    def test_the_kind_goes_and_the_number_stays(self, fortran, expected):
        assert render(fortran, "python") == expected

    def test_an_identifier_that_merely_ends_that_way_is_untouched(self):
        # `n_int32` is a name, not a kinded literal -- the digit boundary is what tells them apart
        assert render("the extent `n_int32`", "python") == "the extent `n_int32`"


class TestLogicals:
    def test_python_gets_python_booleans(self):
        assert render("The default value is `.false.`.", "python") == "The default value is `False`."

    def test_r_gets_r_booleans(self):
        assert render("The default value is `.false.`.", "r") == "The default value is `FALSE`."

    def test_either_case(self):
        assert render(".TRUE. if 1 else .False.", "r") == "TRUE if 1 else FALSE"


class TestProseIsNotRewritten:
    def test_ordinary_text_survives_unchanged(self):
        text = "Number of genes in the study, which must match `gene_to_fam`."
        assert render(text, "python") == text

    def test_a_bare_number_is_left_alone(self):
        assert render("at most 32 chunks", "python") == "at most 32 chunks"


class TestFordBlockTags:
    @pytest.mark.parametrize("tag", ["@note", "@endnote", "@warning", "@endwarning", "  @note  "])
    def test_a_delimiter_is_recognised(self, tag):
        from codegen.emit.doc_literals import is_ford_block_tag

        assert is_ford_block_tag(tag)

    @pytest.mark.parametrize("text", ["@note the caller owns this", "note that", "", "a @note here"])
    def test_prose_is_not(self, text):
        from codegen.emit.doc_literals import is_ford_block_tag

        # only a line that is the delimiter and nothing else -- prose keeps its words
        assert not is_ford_block_tag(text)
