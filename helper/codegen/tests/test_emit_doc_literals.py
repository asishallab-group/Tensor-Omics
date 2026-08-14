"""Rendering the Fortran literals in a doc for the language that will read it.

The author writes on a Fortran impl and the `DM_` range macros quote that Fortran back
verbatim, so a bound arrives as `0_int32` and a default as `.false.`. Carried into a Python
docstring unchanged, the second one contradicts the signature three lines above it.
"""

import pytest

from codegen.emit.doc_literals import render, render_mode_default


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


class TestInlineMaths:
    def test_r_gets_rd_maths(self):
        # `\(` is not an Rd macro; `\eqn` is, and it passes the LaTeX through to the manual
        assert render(r"the ratio \( \frac{a}{b} \)", "r") == r"the ratio \eqn{\frac{a}{b}}"

    def test_python_is_left_alone(self):
        # numpydoc has no equivalent worth inventing here, and the text is harmless
        text = r"the ratio \( \frac{a}{b} \)"
        assert render(text, "python") == text

    def test_several_on_one_line_each_convert(self):
        assert render(r"\( a \) and \( b \)", "r") == r"\eqn{a} and \eqn{b}"


class TestModeDefaults:
    """`DM_DEFAULT` on a mode argument quotes back the integer the author wrote, which is the
    one value a caller of the binding must not pass -- the binding takes the mode string."""

    SENTENCE = "The default value is `1`."

    def test_python_gets_its_own_quotes(self):
        assert (render_mode_default(self.SENTENCE, "robust", "python")
                == "The default value is `'robust'`.")

    def test_r_gets_its_own_quotes(self):
        # matching the type line above it, which lists the modes as "plain", "robust"
        assert (render_mode_default(self.SENTENCE, "robust", "r")
                == 'The default value is `"robust"`.')

    def test_the_c_wrapper_gets_a_fortran_literal(self):
        assert (render_mode_default(self.SENTENCE, "robust", "fortran")
                == "The default value is `'robust'`.")

    def test_the_kinded_form_is_rewritten_too(self):
        # the C wrapper keeps Fortran literals, so the sentence still carries its kind there
        assert (render_mode_default("The default value is `1_int32`.", "plain", "fortran")
                == "The default value is `'plain'`.")

    @pytest.mark.parametrize("text", [
        "The minimum valid value is `1`.",
        "The maximum valid value is `1`.",
        "Defaults to the value 1 when omitted.",
        "The default value of the neighbour count is `1`.",
    ])
    def test_only_the_macro_s_own_sentence_is_touched(self, text):
        # the wording is DM_DEFAULT's, which is what makes rewriting it the generator
        # editing its own output rather than the author's prose
        assert render_mode_default(text, "robust", "python") == text

    def test_surrounding_prose_survives(self):
        text = "Mode for LOESS fitting. The default value is `1`. See the table below."
        assert render_mode_default(text, "robust", "python") == (
            "Mode for LOESS fitting. The default value is `'robust'`. See the table below."
        )
