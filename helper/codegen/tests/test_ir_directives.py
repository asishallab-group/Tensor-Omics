"""Directive parsing.

Patterns are built from the real src/macros.h, so these tests fail if a macro is
reworded into something the parser no longer recognises -- which is the whole point of
deriving the patterns from the definitions.
"""

import pytest

from codegen.frontend.macros import (
    MacroTable,
    MissingMacroError,
    build_directive_patterns,
)
from codegen.ir.directives import (
    AllowInfinite,
    AllowNan,
    Default,
    DirectiveError,
    DirectiveParser,
    Directives,
    Maximum,
    Minimum,
    OptionalOutput,
    OutputFrom,
    OutputFromMode,
    RequiredIfMode,
    ResultSizeIs,
    Sentinel,
)
from codegen.ir.doc import Doc

from conftest import REPO_ROOT


@pytest.fixture(scope="module")
def macros():
    return MacroTable(REPO_ROOT / "src/macros.h", include_paths=(REPO_ROOT,))


@pytest.fixture(scope="module")
def parser(macros):
    return DirectiveParser(build_directive_patterns(macros))


def documented(macros, *invocations):
    """A Doc whose lines are the given macro invocations, expanded as Ford would see them."""
    return Doc.parse([macros.expand(line) for line in invocations])


class TestBuildDirectivePatterns:
    def test_patterns_are_built_from_the_real_header(self, macros):
        patterns = build_directive_patterns(macros)

        assert patterns.default.search(macros.expand("DM_DEFAULT(1)")) is not None

    def test_a_header_without_the_doc_macros_is_rejected(self, tmp_path):
        header = tmp_path / "macros.h"
        header.write_text("#define M_NAN nan\n")

        with pytest.raises(MissingMacroError, match=r"does not define:.*\bDM_DEFAULT\b"):
            build_directive_patterns(MacroTable(header))


class TestDefault:
    @pytest.mark.parametrize(
        "expression",
        ["0.1_real64", "300_int32", "PI", "achar(9)", ".false.", "'x'", "-1_int32"],
    )
    def test_the_expression_is_captured_verbatim(self, macros, parser, expression):
        directives = parser.parse(documented(macros, f"DM_DEFAULT({expression})"))

        assert directives.default == Default(expression, line_number=None)
        assert directives.has_default

    def test_a_default_among_other_prose(self, macros, parser):
        doc = Doc.parse(
            ["LOESS span parameter.", macros.expand("DM_DEFAULT(0.7_real64)"), "See notes."]
        )

        assert parser.parse(doc).default.expression == "0.7_real64"

    def test_the_group_does_not_swallow_a_later_quoted_word(self, macros, parser):
        # The macro body ends in '`.', so a greedy group runs to the *last* '`.' on the
        # line and captures "0.7_real64`. Same as `PI" instead.
        line = macros.expand("DM_DEFAULT(0.7_real64)") + " Same as `PI`."

        assert parser.parse(Doc.parse([line])).default.expression == "0.7_real64"

    def test_no_default_when_absent(self, parser):
        assert parser.parse(Doc.parse(["just prose"])).default is None

    def test_similar_prose_is_not_a_default(self, parser):
        assert parser.parse(Doc.parse(["The default value is unspecified"])).default is None

    def test_the_line_number_is_recorded(self, macros, parser):
        doc = Doc.parse(["prose", macros.expand("DM_DEFAULT(1)")], first_line_number=10)

        assert parser.parse(doc).default.line_number == 11

    def test_duplicates_are_rejected(self, macros, parser):
        doc = documented(macros, "DM_DEFAULT(1)", "DM_DEFAULT(2)")

        with pytest.raises(DirectiveError, match="duplicate DM_DEFAULT directive"):
            parser.parse(doc)


class TestRequiredIfMode:
    def test_fields_are_captured(self, macros, parser):
        directives = parser.parse(
            documented(
                macros,
                "DM_REQUIRED_IF_MODE(pattern_mode, tox_paralog_analysis, MODE_SUBFUNC_PATTERN)",
            )
        )

        assert directives.required_if_mode == RequiredIfMode(
            mode_arg="pattern_mode",
            module="tox_paralog_analysis",
            mode_param="MODE_SUBFUNC_PATTERN",
            line_number=None,
        )

    def test_absent_when_not_documented(self, parser):
        assert parser.parse(Doc.parse(["prose"])).required_if_mode is None


class TestOptionalOutput:
    def test_recognised(self, macros, parser):
        directives = parser.parse(documented(macros, "DM_OPTIONAL_OUTPUT"))

        assert directives.optional_output == OptionalOutput(line_number=None)
        assert directives.is_optional_output

    def test_absent_when_not_documented(self, parser):
        assert not parser.parse(Doc.parse(["prose"])).is_optional_output


class TestResultSizeIs:
    @pytest.mark.parametrize("argument", ["n_results", "num_matches"])
    def test_the_argument_is_captured(self, macros, parser, argument):
        directives = parser.parse(documented(macros, f"DM_RESULT_SIZE_IS({argument})"))

        assert directives.result_size_is == ResultSizeIs(argument, line_number=None)


class TestOutputFrom:
    def test_auto(self, macros, parser):
        directives = parser.parse(
            documented(macros, "DM_OUTPUT_FROM(count, mask_chunk_count, tox_paralog_analysis, AUTO)")
        )

        assert directives.output_from == OutputFrom(
            argument="count",
            procedure="mask_chunk_count",
            module="tox_paralog_analysis",
            mode=OutputFromMode.AUTO,
            line_number=None,
        )
        assert directives.output_from.is_automatic

    def test_just_info(self, macros, parser):
        directives = parser.parse(
            documented(
                macros, "DM_OUTPUT_FROM(count, mask_chunk_count, tox_paralog_analysis, JUST_INFO)"
            )
        )

        assert directives.output_from.mode is OutputFromMode.JUST_INFO
        assert not directives.output_from.is_automatic

    def test_the_two_modes_are_told_apart(self, macros, parser):
        # both expand to the same sentence bar the leading clause
        auto = parser.parse(documented(macros, "DM_OUTPUT_FROM(a, p, m, AUTO)"))
        info = parser.parse(documented(macros, "DM_OUTPUT_FROM(a, p, m, JUST_INFO)"))

        assert auto.output_from.mode is OutputFromMode.AUTO
        assert info.output_from.mode is OutputFromMode.JUST_INFO


class TestRange:
    @pytest.mark.parametrize(
        "expression", ["0.0_real64", "n_values", "above(0.0_real64)", "PI", "1_int32"]
    )
    def test_minimum_is_captured(self, macros, parser, expression):
        directives = parser.parse(documented(macros, f"DM_MIN({expression})"))

        assert directives.minimum == Minimum(expression, line_number=None)
        assert directives.has_range

    @pytest.mark.parametrize(
        "expression", ["100.0_real64", "n_genes", "below(1.0_real64)", "PI"]
    )
    def test_maximum_is_captured(self, macros, parser, expression):
        directives = parser.parse(documented(macros, f"DM_MAX({expression})"))

        assert directives.maximum == Maximum(expression, line_number=None)
        assert directives.has_range

    def test_minimum_and_maximum_together(self, macros, parser):
        directives = parser.parse(
            documented(macros, "DM_MIN(0.0_real64)", "DM_MAX(100.0_real64)")
        )

        assert directives.minimum.expression == "0.0_real64"
        assert directives.maximum.expression == "100.0_real64"

    def test_min_and_max_are_told_apart(self, macros, parser):
        # both are "The <bound> valid value is `EXPR`." -- a lazy group must not let one
        # sentence be read as the other
        assert parser.parse(documented(macros, "DM_MIN(1_int32)")).maximum is None
        assert parser.parse(documented(macros, "DM_MAX(1_int32)")).minimum is None

    @pytest.mark.parametrize("expression", ["0_int32", "-1_int32"])
    def test_sentinel_is_captured(self, macros, parser, expression):
        directives = parser.parse(documented(macros, f"DM_SENTINEL({expression})"))

        assert directives.sentinel == Sentinel(expression, line_number=None)
        assert directives.has_range

    def test_no_range_when_absent(self, parser):
        directives = parser.parse(Doc.parse(["just prose"]))

        assert not directives.has_range
        assert directives.minimum is None
        assert directives.maximum is None
        assert directives.sentinel is None

    def test_duplicate_minimum_is_rejected(self, macros, parser):
        with pytest.raises(DirectiveError, match="duplicate DM_MIN directive"):
            parser.parse(documented(macros, "DM_MIN(1)", "DM_MIN(2)"))


class TestFiniteness:
    def test_allow_nan_is_recognised(self, macros, parser):
        directives = parser.parse(documented(macros, "DM_ALLOW_NAN"))

        assert directives.allow_nan == AllowNan(line_number=None)
        assert directives.allows_nan

    def test_allow_infinite_is_recognised(self, macros, parser):
        directives = parser.parse(documented(macros, "DM_ALLOW_INFINITE"))

        assert directives.allow_infinite == AllowInfinite(line_number=None)
        assert directives.allows_infinite

    def test_the_two_are_told_apart(self, macros, parser):
        assert parser.parse(documented(macros, "DM_ALLOW_NAN")).allow_infinite is None
        assert parser.parse(documented(macros, "DM_ALLOW_INFINITE")).allow_nan is None

    def test_absent_by_default(self, parser):
        directives = parser.parse(Doc.parse(["prose"]))

        assert not directives.allows_nan
        assert not directives.allows_infinite


class TestDirectivesTogether:
    def test_several_directives_on_one_entity(self, macros, parser):
        doc = documented(macros, "DM_DEFAULT(1_int32)", "DM_RESULT_SIZE_IS(n_results)")

        directives = parser.parse(doc)

        assert directives.default.expression == "1_int32"
        assert directives.result_size_is.argument == "n_results"
        assert len(directives.all) == 2

    def test_directives_are_found_in_table_cells(self, macros, parser):
        doc = Doc.parse(["| Note |", "|------|", f"| {macros.expand('DM_DEFAULT(1)')} |"])

        assert parser.parse(doc).default.expression == "1"

    def test_a_default_and_a_required_if_mode_contradict(self, macros, parser):
        doc = documented(
            macros,
            "DM_DEFAULT(1_int32)",
            "DM_REQUIRED_IF_MODE(mode, tox_x, MODE_Y)",
        )

        with pytest.raises(DirectiveError, match="cannot have both a default and a mode"):
            parser.parse(doc)

    def test_empty_directives_are_falsy(self, parser):
        assert not parser.parse(Doc.parse(["prose"]))
        assert not Directives()

    def test_populated_directives_are_truthy(self, macros, parser):
        assert parser.parse(documented(macros, "DM_DEFAULT(1)"))

    def test_all_lists_only_what_was_found(self, macros, parser):
        directives = parser.parse(documented(macros, "DM_OPTIONAL_OUTPUT"))

        assert directives.all == (OptionalOutput(line_number=None),)

    def test_an_undocumented_entity_yields_empty_directives(self, parser):
        assert parser.parse(Doc.parse([])) == Directives()
