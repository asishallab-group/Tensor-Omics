"""The Ford frontend, run against the fixture modules.

These are slow: each parse runs Ford over the fixture sources. The project is parsed once
per session and shared, because what is being checked is the translation, not the parse.
"""

import os
from pathlib import Path

import pytest

from codegen.config import Paths
from codegen.diagnostics import DiagnosticBag
from codegen.frontend.ford_frontend import FordFrontend
from codegen.ir.directives import OutputFromMode
from codegen.ir.roles import analyse_project
from codegen.ir.types import BaseType, Intent
from codegen.ir.validate import validate_project

from conftest import REPO_ROOT

FIXTURE_SRC = Path("helper/codegen/tests/fixtures/src")


def parse(src_dir=FIXTURE_SRC, bag=None):
    bag = bag if bag is not None else DiagnosticBag()
    paths = Paths(root=REPO_ROOT, src_dir=src_dir)
    return FordFrontend(paths, bag).parse(), bag


@pytest.fixture(scope="session")
def parsed():
    result, bag = parse()
    analyse_project(result.project, bag)
    validate_project(result.project, bag)
    return result, bag


@pytest.fixture(scope="session")
def project(parsed):
    return parsed[0].project


@pytest.fixture(scope="session")
def bag(parsed):
    return parsed[1]


class TestParsingProducesCleanIR:
    def test_the_fixtures_produce_no_diagnostics(self, bag):
        # the fixtures are the specification: they must be exemplary
        assert bag.render() == ""

    def test_both_modules_are_found(self, project):
        assert [m.name for m in project] == ["fx_basics", "fx_edges"]

    def test_the_arg_pos_factor_comes_from_the_macro(self, parsed):
        assert parsed[0].arg_pos_factor == 10000

    def test_a_module_records_what_it_uses(self, project):
        # the generated re-export modules are built from these
        assert "tox_errors" in project.module("fx_edges").uses


class TestExportSelection:
    def test_only_tagged_procedures_are_exported(self, project):
        module = project.module("fx_basics")

        assert "fx_internal" in [p.name for p in module]
        assert "fx_internal" not in [p.name for p in module.exported_procedures]

    def test_an_internal_procedure_is_not_held_to_the_contract(self, bag):
        # fx_internal has a deferred-length string and an untabled mode argument
        assert bag.render() == ""


class TestTypes:
    def test_a_kinded_real(self, project):
        argument = project.procedure("fx_basics", "fx_normalize").argument("vector")

        assert argument.type.base is BaseType.REAL
        assert argument.type.kind == "real64"

    def test_a_default_logical_needs_conversion(self, project):
        argument = project.procedure("fx_basics", "fx_optionals").argument("use_quantile")

        assert argument.type.is_logical
        assert not argument.type.is_c_bool
        assert argument.type.needs_conversion

    def test_a_c_bool_logical_needs_no_conversion(self, project):
        argument = project.procedure("fx_edges", "fx_c_bool_flag").argument("flag")

        assert argument.type.is_c_bool
        assert not argument.type.needs_conversion

    @pytest.mark.parametrize(
        "name, length, is_constant",
        [("assumed", "*", False), ("fixed", "8", True), ("sized", "n_chars", False)],
    )
    def test_character_lengths(self, project, name, length, is_constant):
        argument = project.procedure("fx_edges", "fx_strings").argument(name)

        assert argument.type.is_character
        assert str(argument.type.length) == length
        assert argument.type.length.is_constant is is_constant

    def test_a_character_vector_keeps_length_and_rank_apart(self, project):
        argument = project.procedure("fx_edges", "fx_strings").argument("names")

        assert str(argument.type.length) == "16"
        assert argument.rank == 1


class TestDimensions:
    def test_a_dimension_attribute_is_read(self, project):
        # Ford leaves .dimension empty and puts dimension(n_dims) in attribs
        argument = project.procedure("fx_basics", "fx_normalize").argument("vector")

        assert argument.dimension.extents == ("n_dims",)

    def test_a_matrix(self, project):
        argument = project.procedure("fx_basics", "fx_sum_matrix").argument("matrix")

        assert argument.dimension.extents == ("n_rows", "n_cols")
        assert argument.rank == 2

    def test_a_scalar_has_no_extents(self, project):
        assert project.procedure("fx_basics", "fx_normalize").argument("n_dims").is_scalar

    def test_assumed_shape(self, project):
        argument = project.procedure("fx_edges", "fx_serialized").argument("data")

        assert argument.dimension.is_assumed_shape


class TestDimensionAttributes:
    """Reading dimension(...) out of Ford's attribute list."""

    @pytest.mark.parametrize(
        "attribute, expected",
        [
            ("dimension(n)", ("n",)),
            ("dimension(n, m)", ("n", "m")),
            ("dimension(:)", (":",)),
            ("dimension(*)", ("*",)),
            # nested calls: a \(([^)]*)\) pattern stops at the first bracket and
            # returns something unbalanced. All of these appear in the real sources.
            ("dimension(max(0, n_timepoints - 1))", ("max(0, n_timepoints - 1)",)),
            (
                "dimension(max(0_int32, n_timepoints - 1), n_factors, n_samples)",
                ("max(0_int32, n_timepoints - 1)", "n_factors", "n_samples"),
            ),
            ("dimension(clen, product(orig_shape))", ("clen", "product(orig_shape)")),
            ("dimension((n_reps_S1 + n_reps_S2)*n_neighbors)", ("(n_reps_S1 + n_reps_S2)*n_neighbors",)),
            ("DIMENSION(N)", ("N",)),
            ("dimension  ( n )", ("n",)),
            # not a dimension attribute at all
            ("allocatable", None),
            ("target", None),
        ],
    )
    def test_extents_are_extracted_with_bracket_matching(self, attribute, expected):
        from codegen.frontend.ford_frontend import _dimension_attribute
        from codegen.ir.types import Dimension

        extracted = _dimension_attribute(attribute)

        if expected is None:
            assert extracted is None
        else:
            assert Dimension.parse(extracted).extents == expected

    def test_an_unclosed_attribute_yields_nothing_rather_than_a_broken_string(self):
        from codegen.frontend.ford_frontend import _dimension_attribute

        assert _dimension_attribute("dimension(max(0, n") is None


class TestTheWholeRealSource:
    """The frontend must cope with every construct actually written, not just fixtures.

    This is the broadest check there is: parse src/ and require it to come out clean.
    """

    @pytest.fixture(scope="session")
    def real_project(self):
        result, bag = parse(src_dir=Path("src"))
        return result.project, bag

    def test_the_real_sources_parse_without_diagnostics(self, real_project):
        _, bag = real_project

        assert bag.render() == ""

    def test_every_module_is_found(self, real_project):
        project, _ = real_project

        assert len(project) > 30
        assert project.module("tox_errors") is not None
        assert project.module("tox_conversions") is not None

    def test_every_argument_has_a_type(self, real_project):
        # a type the frontend cannot read comes back as None and would break the ABI
        project, _ = real_project

        untyped = [
            f"{m.name}.{p.name}.{a.name}"
            for m in project
            for p in m
            for a in p.arguments
            if a.type is None
        ]

        assert untyped == []

    def test_derived_types_keep_their_name(self, real_project):
        # Ford renders it as a link to its documentation page
        project, _ = real_project

        derived = [
            a.type
            for m in project
            for p in m
            for a in p.arguments
            if a.type is not None and a.type.derived_name
        ]

        assert derived, "no derived types found; this test needs a new example"
        assert all("<" not in t.derived_name for t in derived)

    def test_generic_interfaces_and_derived_types_are_carried(self, real_project):
        """Nothing is generated from either, so nothing more of them is modelled -- but a link
        naming one has to resolve (`validate.check_doc_links`)."""
        project, _ = real_project

        sort = project.module("f42_sort_impl")
        assert {d.name for d in sort.declarations} >= {"sort_array", "sort_array_heapsort"}
        assert all(d.kind == "interface" for d in sort.declarations)

        hashmap = project.module("f42_xxh3_hashmap")
        assert {(d.name, d.kind) for d in hashmap.declarations} == {
            ("hashmap_type", "type"),
            ("hashset_type", "type"),
        }

    def test_an_interface_s_own_documentation_is_read(self, real_project):
        """It belongs to no procedure and no module, so if it is not read here the links in
        it are read nowhere."""
        project, _ = real_project

        interface = next(
            d for d in project.module("f42_sort_impl").declarations
            if d.name == "sort_array_heapsort"
        )

        assert "worst-case" in interface.doc.text
        assert [str(link) for link in interface.doc.links] == [
            "[[f42_sort_impl(module):sort_array(interface)]]"
        ]

    def test_a_private_routine_is_not_carried(self, real_project):
        """Ford documents only what a module makes public, and the link check follows it: a
        link to a private name dangles in the rendered docs exactly as this says it does."""
        project, _ = real_project
        hashmap = project.module("f42_xxh3_hashmap")

        assert "resize_hashmap" not in hashmap.public_names
        assert "hashmap_get" in hashmap.public_names


class TestIntentsAndAttributes:
    @pytest.mark.parametrize(
        "name, intent",
        [("n_dims", Intent.IN), ("vector", Intent.INOUT), ("ierr", Intent.OUT)],
    )
    def test_intents(self, project, name, intent):
        assert project.procedure("fx_basics", "fx_normalize").argument(name).intent is intent

    def test_optional_is_read(self, project):
        procedure = project.procedure("fx_basics", "fx_optionals")

        assert procedure.argument("span").optional
        assert not procedure.argument("values").optional


class TestFunctions:
    def test_a_function_has_a_result(self, project):
        procedure = project.procedure("fx_basics", "fx_count_positive")

        assert procedure.is_function
        assert procedure.result.name == "n_positive"
        assert procedure.result.is_result

    def test_a_result_is_an_output(self, project):
        # Fortran does not say so, but a result is one by definition
        assert project.procedure("fx_basics", "fx_count_positive").result.intent is Intent.OUT

    def test_a_subroutine_has_no_result(self, project):
        assert project.procedure("fx_basics", "fx_normalize").result is None


class TestDocumentation:
    def test_argument_documentation_is_read(self, project):
        argument = project.procedure("fx_basics", "fx_normalize").argument("n_dims")

        assert argument.doc.summary == "number of elements in `vector`"

    def test_the_summary_meta_tag_is_read(self, project):
        procedure = project.procedure("fx_basics", "fx_normalize")

        assert procedure.meta.summary == "Normalizes a vector to unit length in-place"

    def test_the_author_meta_tag_is_stripped_of_markup(self, project):
        # Ford renders author as an anchor element for its HTML pages
        author = project.procedure("fx_basics", "fx_normalize").meta.author

        assert author == "A Developer"
        assert "<" not in author

    def test_a_mode_table_survives_into_the_doc(self, project):
        argument = project.procedure("fx_basics", "fx_modes").argument("mode")

        assert len(argument.doc.tables) == 1


class TestDirectives:
    @pytest.mark.parametrize(
        "name, expression",
        [("span", "0.1_real64"), ("max_iter", "300_int32"), ("use_quantile", ".false.")],
    )
    def test_defaults_are_read(self, project, name, expression):
        argument = project.procedure("fx_basics", "fx_optionals").argument(name)

        assert argument.directives.default.expression == expression

    def test_output_from_is_read(self, project):
        argument = project.procedure("fx_edges", "fx_cluster_expert").argument("n_work")

        directive = argument.directives.output_from
        assert directive.procedure == "fx_work_size"
        assert directive.module == "fx_edges"
        assert directive.mode is OutputFromMode.AUTO

    def test_result_size_resolves_to_the_counting_argument(self, project):
        procedure = project.procedure("fx_edges", "fx_masked")

        assert procedure.argument("results").roles.result_size_arg is procedure.argument("n_results")


class TestRolesOnRealParses:
    def test_extents_link_to_their_arrays(self, project):
        procedure = project.procedure("fx_basics", "fx_sum_matrix")

        assert [a.name for a in procedure.argument("n_cols").roles.extent_of] == [
            "matrix",
            "weights",
        ]

    def test_a_temporary_is_recognised(self, project):
        assert project.procedure("fx_basics", "fx_optionals").argument("tmp_work").roles.is_temporary

    def test_mode_values_are_read_from_the_table(self, project):
        mode = project.procedure("fx_basics", "fx_modes").argument("mode").roles.mode

        assert [v.string for v in mode.values] == ["mean", "median"]
        assert mode.values[0].parameter == "MODE_MEAN"
        assert mode.values[0].module == "fx_basics"

    def test_a_method_argument_uses_the_method_prefix(self, project):
        mode = project.procedure("fx_basics", "fx_modes").argument("link_method").roles.mode

        assert mode.alias == "method"
        assert [v.string for v in mode.values] == ["ward", "single"]

    def test_a_shape_argument_links_to_its_array(self, project):
        procedure = project.procedure("fx_edges", "fx_serialized")

        assert procedure.argument("data_shape").roles.shape_of is procedure.argument("data")

    def test_a_mask_count_links_to_its_mask(self, project):
        procedure = project.procedure("fx_edges", "fx_masked")

        assert procedure.argument("n_selected_genes").roles.mask_count_of is procedure.argument(
            "genes_selection_mask"
        )


class TestParameters:
    def test_module_parameters_are_read_with_their_expressions(self, project):
        module = project.module("fx_basics")

        assert [(p.name, p.expression) for p in module.parameters] == [
            ("MODE_MEAN", "1"),
            ("MODE_MEDIAN", "2"),
            ("METHOD_WARD", "1"),
            ("METHOD_SINGLE", "2"),
        ]

    def test_parameters_evaluate_to_constants(self, project):
        values = project.constant_values()

        assert values["MODE_MEAN"] == 1
        assert values["MODE_MEDIAN"] == 2
        assert values["METHOD_WARD"] == 1
        assert values["METHOD_SINGLE"] == 2
        assert values["MODE_GROUP_ORTHOLOGS"] == 1


class TestSourceLocations:
    def test_a_procedure_points_at_its_statement(self, project):
        location = project.procedure("fx_basics", "fx_normalize").location

        assert location.file.name == "fx_basics.F90"
        source = (REPO_ROOT / FIXTURE_SRC / "fx_basics.F90").read_text().split("\n")
        assert "subroutine fx_normalize" in source[location.line - 1]

    def test_an_argument_points_at_its_declaration(self, project):
        location = project.procedure("fx_basics", "fx_normalize").argument("vector").location

        source = (REPO_ROOT / FIXTURE_SRC / "fx_basics.F90").read_text().split("\n")
        assert ":: vector" in source[location.line - 1]

    def test_locations_are_original_lines_not_preprocessed_ones(self, project):
        # the include of macros.h injects some forty lines before Ford sees line 1
        location = project.module("fx_basics").location
        source = (REPO_ROOT / FIXTURE_SRC / "fx_basics.F90").read_text().split("\n")

        assert source[location.line - 1].strip() == "module fx_basics"


class TestWorkingDirectoryIndependence:
    def test_parsing_does_not_depend_on_the_working_directory(self, tmp_path):
        # fpm.toml configures the preprocessor with '-I.', so a wrong cwd silently
        # drops every DM_ directive rather than failing
        before = os.getcwd()
        try:
            os.chdir(tmp_path)
            result, bag = parse()
        finally:
            os.chdir(before)

        argument = result.project.procedure("fx_basics", "fx_optionals").argument("span")
        assert argument.directives.default.expression == "0.1_real64"
        assert bag.errors == ()

    def test_the_working_directory_is_restored(self):
        before = os.getcwd()

        parse()

        assert os.getcwd() == before


class TestUnexpandedMacros:
    def test_a_file_that_forgets_the_include_is_reported(self, tmp_path):
        # without the include the directive is silently ignored and the wrapper comes
        # out quietly wrong, so it must be an error
        src = tmp_path / "src"
        src.mkdir()
        (src / "fx_forgot.F90").write_text(
            "!> summary: forgot the include\n"
            "module fx_forgot\n"
            "    use, intrinsic :: iso_fortran_env, only: real64, int32\n"
            "    implicit none\n"
            "contains\n"
            "    !> category: C-binding\n"
            "    !| summary: p\n"
            "    !| author: A\n"
            "    subroutine fx_p(span, ierr)\n"
            "        real(real64), intent(in), optional :: span\n"
            "            !! a span.\n"
            "            !! DM_DEFAULT(0.1_real64)\n"
            "        integer(int32), intent(out) :: ierr\n"
            "            !! Error code\n"
            "    end subroutine fx_p\n"
            "end module fx_forgot\n"
        )

        _, bag = parse(src_dir=src)

        assert len(bag.errors) == 1
        error = bag.errors[0]
        assert "'DM_DEFAULT' looks like a documentation macro but never expanded" in error.message
        assert "#include" in error.note
        assert error.location.line is not None

    def test_a_misspelt_doc_macro_is_an_error(self, tmp_path):
        # a DM_ typo is silent: the default is simply never found
        _, bag = parse(src_dir=_module_with(tmp_path, "!! DM_DEFAULR(0.1_real64)"))

        assert len(bag.errors) == 1
        assert "'DM_DEFAULR'" in bag.errors[0].message
        assert "documentation macro" in bag.errors[0].message

    def test_a_misspelt_code_macro_is_a_warning(self, tmp_path):
        # an M_ typo leaves the name standing in the docs: wrong, but visible, and
        # nothing is mis-generated from it
        _, bag = parse(src_dir=_module_with(tmp_path, "!! see `M_GENE_TO_FAM_SENTINELL`"))

        assert bag.errors == ()
        assert len(bag.warnings) == 1
        assert "'M_GENE_TO_FAM_SENTINELL'" in bag.warnings[0].message
        assert "a macro from a macro header" in bag.warnings[0].message

    def test_a_misspelt_custom_macro_is_a_warning(self, tmp_path):
        _, bag = parse(src_dir=_module_with(tmp_path, "!! default: `CM_SPAN_DEFAUL`"))

        assert bag.errors == ()
        assert "custom macro defined in the source file itself" in bag.warnings[0].message

    def test_a_correctly_spelt_macro_is_not_reported(self, tmp_path):
        # M_GENE_TO_FAM_SENTINEL is defined, so it expands and leaves nothing behind
        _, bag = parse(src_dir=_module_with(tmp_path, "!! see `M_GENE_TO_FAM_SENTINEL`"))

        assert bag.render() == ""

    def test_several_typos_on_one_line_are_all_reported(self, tmp_path):
        _, bag = parse(src_dir=_module_with(tmp_path, "!! `M_NOPE` and `CM_ALSO_NOPE`"))

        assert len(bag.warnings) == 2


def _module_with(tmp_path, doc_line):
    """A minimal exported procedure carrying `doc_line` on its argument."""
    src = tmp_path / "src"
    src.mkdir()
    (src / "fx_typo.F90").write_text(
        "#include <src/macros.h>\n"
        "\n"
        "!> summary: typo fixture\n"
        "module fx_typo\n"
        "    use, intrinsic :: iso_fortran_env, only: real64, int32\n"
        "    implicit none\n"
        "contains\n"
        "    !> category: C-binding\n"
        "    !| summary: p\n"
        "    !| author: A\n"
        "    subroutine fx_p(span, ierr)\n"
        "        real(real64), intent(in) :: span\n"
        "            !! a span.\n"
        f"            {doc_line}\n"
        "        integer(int32), intent(out) :: ierr\n"
        "            !! Error code\n"
        "    end subroutine fx_p\n"
        "end module fx_typo\n"
    )
    return src


class TestAllocatableLocals:
    """The generator never reads a body, so an allocation is seen through its declaration."""

    def parse_impl(self, tmp_path, declarations):
        src = tmp_path / "src"
        src.mkdir()
        (src / "fx_alloc_impl.F90").write_text(
            "!> summary: an implementation module\n"
            "module fx_alloc_impl\n"
            "    use, intrinsic :: iso_fortran_env, only: real64, int32\n"
            "    implicit none\n"
            "contains\n"
            "    !> summary: p\n"
            "    !| author: A\n"
            "    subroutine fx_p_impl(n)\n"
            "        integer(int32), intent(in) :: n\n"
            "            !! a count\n"
            f"{declarations}"
            "    end subroutine fx_p_impl\n"
            "end module fx_alloc_impl\n"
        )
        parsed, bag = parse(src_dir=src)
        return parsed.project.procedure("fx_alloc_impl", "fx_p_impl"), bag

    def test_an_allocatable_local_is_recorded(self, tmp_path):
        procedure, _ = self.parse_impl(
            tmp_path,
            "        real(real64), dimension(:), allocatable :: workspace, scratch\n",
        )

        assert procedure.allocatable_locals == ("workspace", "scratch")

    def test_the_rule_fires_on_it(self, tmp_path):
        procedure, bag = self.parse_impl(
            tmp_path,
            "        real(real64), dimension(:), allocatable :: workspace\n",
        )
        validate_project(procedure.module.project, bag)

        assert any("allocates: 'workspace'" in e.message for e in bag.errors), bag.render()

    def test_a_pointer_local_is_not_an_allocation(self, tmp_path):
        procedure, _ = self.parse_impl(
            tmp_path, "        real(real64), dimension(:), pointer :: view\n"
        )

        assert procedure.allocatable_locals == ()

    def test_a_plain_local_is_not_an_allocation(self, tmp_path):
        procedure, _ = self.parse_impl(
            tmp_path, "        integer(int32) :: i\n"
        )

        assert procedure.allocatable_locals == ()

    def test_an_argument_is_never_a_local(self, tmp_path):
        # `n` is a dummy; only what is left after Ford splits them off counts
        procedure, _ = self.parse_impl(tmp_path, "")

        assert procedure.allocatable_locals == ()


class TestTypesFordCannotBeReadFrom:
    """The two refusals in `_type`, reached directly rather than through a parse.

    `_type` sees only Ford's rendered `full_type` string, so a stand-in carrying that one
    attribute exercises it exactly as a parse would. Note the deliberately fresh bag: the
    module-wide `bag` fixture is the *parse's* bag, shared for the session.
    """

    def _type_of(self, full_type, report=True):
        from types import SimpleNamespace

        from codegen.diagnostics import SourceLocation

        bag = DiagnosticBag()
        variable = SimpleNamespace(name="x", full_type=full_type, kind=None)
        result = FordFrontend(Paths(), diagnostics=bag)._type(
            variable, SourceLocation(file=Path("src/m.F90")), report=report
        )
        return result, [d.message for d in bag.errors]

    def test_a_type_spec_that_does_not_parse_is_reported(self):
        result, errors = self._type_of("???")

        assert result is None
        assert "cannot read the type of 'x'" in errors[0]

    def test_a_type_the_generator_does_not_support_is_reported(self):
        result, errors = self._type_of("procedure")

        assert result is None
        assert len(errors) == 1

    def test_a_type_that_parses_but_cannot_be_built_is_reported(self):
        # `real` with no kind: the base type is fine, but a default kind is processor
        # dependent and so has no defensible C mapping
        result, errors = self._type_of("real")

        assert result is None
        assert "has an unusable type 'real'" in errors[0]
        assert "a kind is required" in errors[0]

    def test_report_false_stays_quiet(self):
        # `_parameter` reads types this way for entities that need not be interoperable
        result, errors = self._type_of("real", report=False)

        assert result is None
        assert errors == []


class TestFordGivingUpOnAFile:
    """Ford warns and carries on; the frontend has to turn that into a diagnostic.

    Reached through `_report_unparsed` with the text Ford actually writes, rather than by
    corrupting a fixture source: the point is the reporting, and a real parse of a broken
    file is slow and leaves the failure mode to Ford's version.
    """

    def _report(self, captured):
        bag = DiagnosticBag()
        FordFrontend(Paths(), diagnostics=bag)._report_unparsed(captured)
        return [d.message for d in bag.errors]

    def test_a_complaint_becomes_an_error(self):
        errors = self._report(
            "ERROR in file 'tox_normalization_impl.F90': END statement outside of any nesting"
        )

        assert len(errors) == 1
        assert "Ford could not read a source file" in errors[0]
        assert "tox_normalization_impl.F90" in errors[0]

    def test_fords_own_give_up_warning_counts_too(self):
        # the per-file `except` in fortran_project warns with this exact wording; a generic
        # "Warning:" does not count, which is the whole point of matching the two phrases
        errors = self._report("Warning: Error parsing src/m.F90.")

        assert len(errors) == 1
        assert "m.F90" in errors[0]

    def test_ordinary_narration_is_not_an_error(self):
        # whatever Ford prints about its own progress is not a complaint about a source
        errors = self._report("Reading file src/m.F90\nCorrelating information\n\n")

        assert errors == []

    def test_an_unknown_metadata_key_is_not_a_dropped_file(self):
        # Ford warns about these on every run of this repository and keeps the file. Taking
        # them for parse failures is what turned six benign warnings into a red CI.
        errors = self._report(
            "Warning: In 'f42_sort_impl.F90:unknown': Ignoring unknown Ford metadata key "
            "'helper'"
        )

        assert errors == []

    def test_colour_does_not_decide_what_counts(self):
        # rich colours its output at a terminal and not in CI, so a matcher that sees the
        # escape codes behaves differently in the two -- which is exactly how this broke
        coloured = (
            "\x1b[1;31mWarning:\x1b[0m In \x1b[32m'f42_sort_impl.F90:unknown'\x1b[0m: "
            "Ignoring unknown Ford metadata key \x1b[32m'helper'\x1b[0m"
        )

        assert self._report(coloured) == []
        assert len(self._report("\x1b[1;31mERROR in file 'a.F90'\x1b[0m: broken")) == 1

    def test_a_message_rich_wrapped_across_lines_is_still_one_complaint(self):
        # rich hard-wraps to the console width, so a per-line matcher sees only a fragment
        errors = self._report("ERROR in file\n'tox_normalization_impl.F90': END statement\noutside")

        assert len(errors) == 1
        assert "tox_normalization_impl.F90" in errors[0]

    def test_each_complaint_is_reported_separately(self):
        errors = self._report("ERROR in file 'a.F90': bad\nERROR in file 'b.F90': worse")

        assert len(errors) == 2

    def test_one_file_is_reported_once_however_often_ford_says_it(self):
        # its reader raises `ERROR in file`, then the per-file except above warns
        # `Error parsing` about the same file
        errors = self._report(
            "ERROR in file 'm.F90': bad Warning: Error parsing src/m.F90."
        )

        assert len(errors) == 1

    def test_the_detail_stops_at_ford_going_back_to_narrating(self):
        errors = self._report(
            "ERROR in file 'm.F90': bad Preprocessing /some/other/file.F90"
        )

        assert "Preprocessing" not in errors[0]
