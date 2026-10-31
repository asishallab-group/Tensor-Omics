import io
import re
from contextlib import redirect_stdout

import pytest

from codegen.config import CONVENTIONS, Paths
from codegen.frontend.ford_frontend import FordFrontend
from codegen.frontend.macros import (
    DOC_MACROS,
    EXPORT_MACRO,
    OPENMP_SWITCH,
    MacroTable,
    MissingMacroError,
    export_category,
    require_omp_macros,
)

from conftest import REPO_ROOT


@pytest.fixture(scope="module")
def macros():
    return MacroTable(REPO_ROOT / "src/macros.h", include_paths=(REPO_ROOT,))


@pytest.fixture
def header(tmp_path):
    def write(text):
        path = tmp_path / "macros.h"
        path.write_text(text)
        return MacroTable(path)

    return write


def test_defines_are_collected(macros):
    assert "M_ALLOCATE" in macros
    assert "DM_DEFAULT" in macros
    assert "NOT_A_MACRO" not in macros


def test_included_headers_are_followed(macros):
    # src/macros.h opens with #include <authors.h> from the repo root
    assert "AUTHOR_FRANZ_ERIC_SILL" in macros.names


def test_loading_writes_nothing_to_stdout(header):
    table = header("#define M_NAN quiet_nan\n")

    buffer = io.StringIO()
    with redirect_stdout(buffer):
        table.expand("M_NAN")

    assert buffer.getvalue() == ""


def test_expand_resolves_an_object_macro(macros):
    assert macros.expand("M_NAN") == "ieee_value(1.0_real64, ieee_quiet_nan)"


def test_expand_resolves_a_function_macro(header):
    table = header("#define DM_DEFAULT(VAL) The default value is `VAL`.\n")

    assert table.expand("DM_DEFAULT(0.7_real64)") == "The default value is `0.7_real64`."


def test_expand_leaves_unknown_text_alone(header):
    table = header("#define M_X y\n")

    assert table.expand("plain prose") == "plain prose"


def test_pattern_escapes_the_macro_body_but_not_the_arguments(header):
    table = header("#define DM_DEFAULT(VAL) The default value is `VAL`.\n")

    pattern = table.pattern("DM_DEFAULT((?P<value>.*))")

    # the argument stays a live group
    assert "(?P<value>.*)" in pattern
    assert re.compile(pattern).match("The default value is `0.7_real64`.") is not None
    # while the body's full stop is literal and no longer matches any character
    assert re.compile(pattern).match("The default value is `0.7_real64`X") is None


def test_pattern_does_not_match_similar_prose(header):
    table = header("#define DM_DEFAULT(VAL) The default value is `VAL`.\n")
    compiled = table.compiled("DM_DEFAULT((?P<value>.*))")

    assert compiled.match("The default value is unspecified") is None


def test_compiled_captures_the_argument(header):
    table = header("#define DM_DEFAULT(VAL) The default value is `VAL`.\n")

    match = table.compiled("DM_DEFAULT((?P<value>.*))").match(
        "The default value is `0.7_real64`."
    )

    assert match.group("value") == "0.7_real64"


def test_pattern_survives_regex_metacharacters_in_the_body(header):
    # '.', '*', '(', ')' and '[[...]]' are all regex metacharacters
    table = header("#define DM_NOTE(X) See [[mod(module):X]] (*important*).\n")
    compiled = table.compiled("DM_NOTE((?P<item>.*))")

    assert compiled.match("See [[mod(module):foo]] (*important*).").group("item") == "foo"
    # the metacharacters must be literal, not a wildcard match
    assert compiled.match("See XXmodXmoduleXXfooXX XXimportantXXX") is None


def test_escaped_and_plain_tables_do_not_interfere(header):
    table = header("#define DM_DEFAULT(VAL) The default value is `VAL`.\n")

    table.pattern("DM_DEFAULT((?P<value>.*))")

    # escaping mutates macro bodies in place, so a shared preprocessor would corrupt this
    assert table.expand("DM_DEFAULT(1)") == "The default value is `1`."


def test_token_pasting_selects_the_variant(macros):
    auto = macros.expand("DM_OUTPUT_FROM(n_work, sizes_alloc, tox_x, AUTO)")
    just_info = macros.expand("DM_OUTPUT_FROM(n_work, sizes_alloc, tox_x, JUST_INFO)")

    assert auto.startswith("It is *VERY IMPORTANT*")
    assert just_info.startswith("It is recommended")
    assert "[[tox_x(module):sizes_alloc]]" in auto


class TestRealDocMacros:
    """The DM_ macros in src/macros.h must expand and round-trip into patterns."""

    def test_every_known_doc_macro_is_defined(self, macros):
        missing = [name for name in DOC_MACROS.all() if name not in macros]

        assert missing == []

    def test_default(self, macros):
        match = macros.compiled("DM_DEFAULT((?P<value>.*))").match(
            macros.expand("DM_DEFAULT(0.7_real64)")
        )

        assert match.group("value") == "0.7_real64"

    def test_required_if_mode(self, macros):
        template = "DM_REQUIRED_IF_MODE((?P<mode_arg>.*), (?P<module>.*), (?P<mode_param>.*))"
        expanded = macros.expand("DM_REQUIRED_IF_MODE(link_method, tox_clustering, METHOD_WARD)")

        match = macros.compiled(template).match(expanded)

        assert match.group("mode_arg") == "link_method"
        assert match.group("module") == "tox_clustering"
        assert match.group("mode_param") == "METHOD_WARD"

    def test_result_size_is(self, macros):
        match = macros.compiled("DM_RESULT_SIZE_IS((?P<argument>.*))").match(
            macros.expand("DM_RESULT_SIZE_IS(n_results)")
        )

        assert match.group("argument") == "n_results"

    def test_output_from(self, macros):
        template = (
            "DM_OUTPUT_FROM((?P<argument>.*), (?P<procedure>.*), (?P<module>.*), AUTO)"
        )
        expanded = macros.expand("DM_OUTPUT_FROM(n_work, sizes_alloc, tox_x, AUTO)")

        match = macros.compiled(template).match(expanded)

        assert match.group("argument") == "n_work"
        assert match.group("procedure") == "sizes_alloc"
        assert match.group("module") == "tox_x"

    def test_doc_macros_expand_to_prose_not_code(self, macros):
        # a doc macro that still looked like a macro call would leak into the docs
        expanded = macros.expand("DM_DEFAULT(0.7_real64)")

        assert "DM_" not in expanded

    def test_doc_macros_avoid_apostrophes(self, macros):
        # pcpp tokenises a lone apostrophe as an unterminated character literal
        for invocation in (
            "DM_DEFAULT(1)",
            "DM_MIN(1_int32)",
            "DM_RESULT_SIZE_IS(n)",
            "DM_OUTPUT_FROM(a, p, m, AUTO)",
            "DM_OUTPUT_FROM(a, p, m, JUST_INFO)",
            "DM_REQUIRED_IF_MODE(a, m, MODE_X)",
        ):
            assert "'" not in macros.expand(invocation), invocation


class TestExportCategory:
    """M_EXPORT_C is the single source of the export marker."""

    def test_the_macro_is_defined(self, macros):
        assert EXPORT_MACRO in macros

    def test_the_category_is_read_from_the_macro(self, macros):
        assert export_category(macros) == "C-binding"

    def test_it_expands_to_a_ford_category_tag(self, macros):
        # so Ford still parses it into meta.category
        assert macros.expand(EXPORT_MACRO).strip().lower().startswith("category:")

    def test_a_header_without_the_macro_is_reported(self, tmp_path):
        header = tmp_path / "macros.h"
        header.write_text("#define M_NAN nan\n")

        with pytest.raises(MissingMacroError, match=EXPORT_MACRO):
            export_category(MacroTable(header))

    def test_a_macro_that_is_not_a_category_tag_is_reported(self, tmp_path):
        header = tmp_path / "macros.h"
        header.write_text(f"#define {EXPORT_MACRO} not a category\n")

        with pytest.raises(MissingMacroError, match="category: <name>"):
            export_category(MacroTable(header))

    def test_the_value_follows_a_changed_macro(self, tmp_path):
        header = tmp_path / "macros.h"
        header.write_text(f"#define {EXPORT_MACRO} category: C-export\n")

        assert export_category(MacroTable(header)) == "C-export"



@pytest.fixture(scope="module")
def with_openmp():
    """The real header as a compiler with OpenMP reads it."""
    return MacroTable(
        REPO_ROOT / "src/macros.h", include_paths=(REPO_ROOT,), defines=(OPENMP_SWITCH,)
    )


def header_without(tmp_path, name, branches=None):
    """The real header without the definitions of `name` in `branches` (0: with OpenMP, 1: without),
    or without every definition of it."""
    lines = (REPO_ROOT / "src/macros.h").read_text().splitlines(keepends=True)
    definitions = [i for i, line in enumerate(lines) if line.startswith(f"#define {name}")]
    if branches is None:
        dropped = set(definitions)
    else:
        assert len(definitions) == 2, "one definition per branch of #ifdef _OPENMP"
        dropped = {definitions[branch] for branch in branches}
    assert dropped, f"the header defines {name}"
    header = tmp_path / "macros.h"
    header.write_text("".join(line for i, line in enumerate(lines) if i not in dropped))
    return MacroTable(header, include_paths=(REPO_ROOT,))


#: Defined in each branch of the header's `#ifdef _OPENMP`, differently
OMP_BRANCH_MACRO_NAMES = [CONVENTIONS.omp_pure_macro, CONVENTIONS.omp_n_threads_macro]
#: Every OpenMP macro the generator requires: those, and the parallel loop's, defined once
OMP_MACRO_NAMES = OMP_BRANCH_MACRO_NAMES + [CONVENTIONS.omp_parallel_macro]


class TestOmpMacros:
    """The OpenMP macros: the generator relies on their names, the compilers on their bodies."""

    def test_the_real_header_defines_them(self, macros):
        require_omp_macros(macros)

    def test_the_purity_follows_the_switch(self, macros, with_openmp):
        assert macros.expand(CONVENTIONS.omp_pure_macro).strip() == "pure"
        assert with_openmp.expand(CONVENTIONS.omp_pure_macro).strip() == "impure"

    def test_the_thread_count_follows_the_switch(self, macros, with_openmp):
        invocation = f"{CONVENTIONS.omp_n_threads_macro}(n_omp_threads)"
        assert macros.expand(invocation).strip() == "1_int32"
        assert with_openmp.expand(invocation).strip() == "n_omp_threads"

    def test_the_parallel_loop_is_one_directive_in_both_builds(self, macros, with_openmp):
        # the `!` stays in the source line, so the expansion is a comment without OpenMP
        line = "!M_OMP_PARALLEL(n_items > 1) shared(values) private(i_thread)"
        directive = (
            "!$omp parallel do default(none) num_threads(n_omp_threads)"
            " if(n_omp_threads > 1 .and. (n_items > 1)) schedule(dynamic)"
            " shared(values) private(i_thread)"
        )
        assert macros.expand(line) == directive
        assert with_openmp.expand(line) == directive

    @pytest.mark.parametrize("name", OMP_MACRO_NAMES)
    def test_a_header_without_one_is_reported(self, tmp_path, name):
        with pytest.raises(MissingMacroError, match=rf"{name} without {OPENMP_SWITCH}"):
            require_omp_macros(header_without(tmp_path, name))

    @pytest.mark.parametrize("name", OMP_BRANCH_MACRO_NAMES)
    def test_one_missing_only_with_openmp_is_reported(self, tmp_path, name):
        with pytest.raises(MissingMacroError, match=rf"{name} with {OPENMP_SWITCH}"):
            require_omp_macros(header_without(tmp_path, name, branches=(0,)))

    def test_the_parallel_loop_defined_only_without_openmp_is_reported(self, tmp_path):
        # defined once in the real header, but a header that moves it into one branch is
        # checked in both, like the other two
        name = CONVENTIONS.omp_parallel_macro
        text = (REPO_ROOT / "src/macros.h").read_text()
        definition = next(line for line in text.splitlines(keepends=True)
                          if line.startswith(f"#define {name}"))
        text = text.replace(definition, "").replace("#else\n", f"#else\n{definition}", 1)
        header = tmp_path / "macros.h"
        header.write_text(text)

        with pytest.raises(MissingMacroError, match=rf"{name} with {OPENMP_SWITCH}"):
            require_omp_macros(MacroTable(header, include_paths=(REPO_ROOT,)))

    @pytest.mark.parametrize("name", OMP_MACRO_NAMES)
    def test_the_generator_refuses_such_a_header(self, tmp_path, name):
        header = header_without(tmp_path, name)

        with pytest.raises(MissingMacroError, match=name):
            FordFrontend(Paths(root=REPO_ROOT, macros_header=header.header))
