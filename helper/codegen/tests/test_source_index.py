import re
from pathlib import Path

import pytest

from codegen.frontend.source_index import SourceFile, SourceIndex

from conftest import REPO_ROOT

SOURCE = """\
#include <src/macros.h>

!> Module with normalization routines.
module tox_normalization
    use safeguard
    implicit none

contains

    !> category: C-interface
    !| Normalizes an input vector to unit length in-place
    pure subroutine normalize_unit_length(vector, n_dims, ierr)
        integer(int32), intent(in) :: n_dims
            !! number of elements in `vector`
        real(real64), dimension(n_dims), intent(inout) :: vector
            !! Vector that will be normalized to unit length
        integer(int32), intent(out) :: ierr
            !! Error code

        real(real64) :: vector_norm

        ! n_dims is used here, but this is not a declaration
        call validate_dimension_size(n_dims, ierr)
        vector_norm = norm(vector)
    end subroutine normalize_unit_length

    !> A function
    pure integer(int32) function get_gene_index(gene_id, genes) result(idx)
        character(len=*), intent(in) :: gene_id
            !! the gene to look for
        integer(int32) :: idx
    end function get_gene_index

    elemental impure subroutine odd_one(a)
        real(real64), intent(in) :: a
    end subroutine odd_one
end module tox_normalization
"""


@pytest.fixture
def source():
    return SourceFile(Path("src/tox/tox_normalization.F90"), SOURCE)


def line_of(text):
    """The 1-based line of the first line containing `text`, for readable expectations."""
    for number, line in enumerate(SOURCE.split("\n"), start=1):
        if text in line:
            return number
    raise AssertionError(f"{text!r} not in the fixture")


class TestModules:
    def test_a_module_is_found(self, source):
        assert source.module_line("tox_normalization") == line_of("module tox_normalization")

    def test_lookup_is_case_insensitive(self, source):
        assert source.module_line("TOX_NORMALIZATION") is not None

    def test_an_absent_module_is_none(self, source):
        assert source.module_line("tox_absent") is None

    def test_end_module_is_not_the_declaration(self, source):
        assert source.module_line("tox_normalization") < line_of("end module")


class TestProcedures:
    def test_a_subroutine_is_found(self, source):
        assert source.procedure_line("normalize_unit_length") == line_of(
            "pure subroutine normalize_unit_length"
        )

    def test_a_function_with_a_type_prefix_and_result_is_found(self, source):
        assert source.procedure_line("get_gene_index") == line_of("function get_gene_index")

    def test_several_prefixes_are_handled(self, source):
        assert source.procedure_line("odd_one") == line_of("elemental impure subroutine odd_one")

    def test_a_call_is_not_mistaken_for_a_declaration(self, source):
        # 'call validate_dimension_size(...)' must not look like a declaration
        assert source.procedure_line("validate_dimension_size") is None

    def test_an_absent_procedure_is_none(self, source):
        assert source.procedure_line("nope") is None

    def test_a_similarly_named_procedure_is_not_matched(self, source):
        assert source.procedure_line("normalize") is None
        assert source.procedure_line("normalize_unit") is None

    def test_the_end_statement_is_not_the_declaration(self, source):
        assert source.procedure_line("normalize_unit_length") < line_of(
            "end subroutine normalize_unit_length"
        )


class TestProcedureSpans:
    def test_a_span_runs_from_the_statement_to_its_end(self, source):
        start, end = source.procedure_span("normalize_unit_length")

        assert start == line_of("pure subroutine normalize_unit_length")
        assert end == line_of("end subroutine normalize_unit_length")

    def test_an_absent_procedure_has_no_span(self, source):
        assert source.procedure_span("nope") is None


class TestArguments:
    @pytest.mark.parametrize(
        "argument, declaration",
        [
            ("n_dims", "integer(int32), intent(in) :: n_dims"),
            ("vector", "real(real64), dimension(n_dims), intent(inout) :: vector"),
            ("ierr", "integer(int32), intent(out) :: ierr"),
        ],
    )
    def test_arguments_are_found(self, source, argument, declaration):
        assert source.argument_line("normalize_unit_length", argument) == line_of(declaration)

    def test_a_use_in_the_body_is_not_a_declaration(self, source):
        # n_dims appears in a call below its declaration; the declaration must win
        assert source.argument_line("normalize_unit_length", "n_dims") == line_of(
            "integer(int32), intent(in) :: n_dims"
        )

    def test_a_name_only_inside_a_dimension_attribute_is_not_a_declaration(self, source):
        # 'dimension(n_dims)' is left of the '::', so it must not match
        assert source.argument_line("normalize_unit_length", "n_dims") != line_of(
            "real(real64), dimension(n_dims), intent(inout) :: vector"
        )

    def test_a_name_only_in_a_doc_comment_is_not_a_declaration(self, source):
        # the doc for n_dims mentions `vector`
        assert source.argument_line("normalize_unit_length", "vector") == line_of(
            "real(real64), dimension(n_dims), intent(inout) :: vector"
        )

    def test_the_search_is_confined_to_the_procedure(self, source):
        # gene_id is declared in get_gene_index, not in normalize_unit_length
        assert source.argument_line("normalize_unit_length", "gene_id") is None
        assert source.argument_line("get_gene_index", "gene_id") is not None

    def test_an_absent_argument_is_none(self, source):
        assert source.argument_line("normalize_unit_length", "nope") is None

    def test_an_argument_of_an_absent_procedure_is_none(self, source):
        assert source.argument_line("nope", "n_dims") is None

    def test_a_local_variable_is_found_too(self, source):
        # the index does not know what is a dummy and what is a local; that is fine
        assert source.argument_line("normalize_unit_length", "vector_norm") is not None


class TestDocLines:
    def test_argument_documentation_starts_after_the_declaration(self, source):
        assert source.argument_doc_line("normalize_unit_length", "n_dims") == line_of(
            "number of elements in `vector`"
        )

    def test_documentation_of_an_absent_argument_is_none(self, source):
        assert source.argument_doc_line("normalize_unit_length", "nope") is None


class TestVariables:
    def test_a_module_parameter_is_found(self):
        source = SourceFile(
            Path("x.F90"),
            "module tox_errors\n"
            "    integer(int32), parameter :: ERR_OK = 0\n"
            "        !! no error\n"
            "end module tox_errors\n",
        )

        assert source.variable_line("ERR_OK") == 2
        assert source.variable_doc_line("ERR_OK") == 3


class TestSourceIndex:
    def test_locations_carry_the_path_and_line(self, tmp_path):
        path = tmp_path / "m.F90"
        path.write_text(SOURCE)
        index = SourceIndex(tmp_path)

        location = index.procedure(Path("m.F90"), "normalize_unit_length")

        assert location.file == Path("m.F90")
        assert location.line == line_of("pure subroutine normalize_unit_length")

    def test_a_missing_file_yields_a_location_without_a_line(self, tmp_path):
        index = SourceIndex(tmp_path)

        location = index.procedure(Path("absent.F90"), "p")

        assert location.file == Path("absent.F90")
        assert location.line is None

    def test_an_unfindable_declaration_yields_no_line_rather_than_a_guess(self, tmp_path):
        path = tmp_path / "m.F90"
        path.write_text(SOURCE)
        index = SourceIndex(tmp_path)

        assert index.procedure(Path("m.F90"), "absent").line is None

    def test_files_are_read_once(self, tmp_path):
        path = tmp_path / "m.F90"
        path.write_text(SOURCE)
        index = SourceIndex(tmp_path)

        first = index.file(Path("m.F90"))
        second = index.file(Path("m.F90"))

        assert first is second


#: Every Fortran source in the repo, so odd real signatures are exercised
REAL_SOURCES = sorted((REPO_ROOT / "src").rglob("*.[fF]90"))

#: Matches a procedure statement the way the index does, to find what *should* be found
DECLARATION_RE = re.compile(
    r"^\s*(?:[\w()=*:,\s]*?\s)??\b(?:subroutine|function)\s+(\w+)\s*\(", re.MULTILINE
)


@pytest.fixture(scope="module")
def errors_file():
    return SourceFile(REPO_ROOT / "src/tox/tox_errors.F90")


class TestAgainstRealSources:
    """The index must work on the sources as written, not only on a tailored fixture."""

    def test_the_tox_errors_module_is_found(self, errors_file):
        line = errors_file.module_line("tox_errors")

        assert line is not None
        assert errors_file.lines[line - 1].strip() == "module tox_errors"

    @pytest.mark.parametrize(
        "name",
        ["create_err_code", "get_err_code", "get_err_arg_pos", "set_err", "is_err", "set_ok"],
    )
    def test_real_procedures_are_found_on_their_own_statement(self, errors_file, name):
        line = errors_file.procedure_line(name)

        assert line is not None, f"{name} not found"
        text = errors_file.lines[line - 1]
        assert re.search(rf"\b(?:subroutine|function)\s+{name}\b", text), text

    def test_a_real_argument_is_found(self, errors_file):
        line = errors_file.argument_line("create_err_code", "arg_pos")

        assert line is not None
        assert "::" in errors_file.lines[line - 1]
        assert "arg_pos" in errors_file.lines[line - 1]

    def test_a_real_parameter_is_found_with_its_documentation(self, errors_file):
        line = errors_file.variable_line("ERR_INVALID_INPUT")

        assert line is not None
        assert "ERR_INVALID_INPUT = 201" in errors_file.lines[line - 1]
        assert "invalid input" in errors_file.lines[errors_file.variable_doc_line("ERR_INVALID_INPUT") - 1]

    def test_every_procedure_in_tox_errors_is_locatable(self, errors_file):
        names = self._declared_names(errors_file.path.read_text())

        assert names, "no procedures found in tox_errors"
        assert [name for name in names if errors_file.procedure_line(name) is None] == []

    @pytest.mark.parametrize("path", REAL_SOURCES, ids=lambda p: p.name)
    def test_every_procedure_in_every_source_is_locatable(self, path):
        # the prefix pattern has to cope with every signature actually written, not
        # only the ones the fixture happens to contain
        source = SourceFile(path)
        names = self._declared_names(path.read_text(errors="replace"))

        unlocatable = [name for name in names if source.procedure_line(name) is None]

        assert unlocatable == [], f"{path.name}: cannot locate {unlocatable}"

    @pytest.mark.parametrize("path", REAL_SOURCES, ids=lambda p: p.name)
    def test_a_located_line_really_declares_the_procedure(self, path):
        # a wrong line is worse than no line, so check every hit lands on a statement
        source = SourceFile(path)

        for name in self._declared_names(path.read_text(errors="replace")):
            line = source.procedure_line(name)
            text = source.lines[line - 1]
            assert re.search(rf"\b(?:subroutine|function)\s+{re.escape(name)}\b", text), (
                f"{path.name}:{line} does not declare {name}: {text!r}"
            )

    @staticmethod
    def _declared_names(text: str) -> set[str]:
        return {
            name
            for name in DECLARATION_RE.findall(text)
            if not name.lower().startswith("end")
        }
