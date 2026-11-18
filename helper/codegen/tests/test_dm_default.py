"""`DM_DEFAULT` is allowed only on an `optional` argument.

The rule: on a required scalar it is an error, read off the declaration as written (the
implementation's, or a hand-written export's), reported once at the directive's line. The note
says what to do for the case the argument is in. An array keeps its own refusal
(`test_ir_validate`), whether or not it is optional.
"""

import textwrap

import pytest

from codegen.config import CONVENTIONS, Paths
from codegen.diagnostics import DiagnosticBag
from codegen.frontend.ford_frontend import FordFrontend
from codegen.generate import generate
from codegen.ir.directives import Default, Directives, OutputFrom, OutputFromMode
from codegen.ir.entities import Meta
from codegen.ir.roles import analyse_project
from codegen.ir.types import Intent
from codegen.ir.validate import check_defaults_are_on_optionals, validate_project
from codegen.synthesize import synthesize_wrappers

import builders as b
from conftest import REPO_ROOT

IMPL_META = Meta(summary="Scale the values", author="AUTHOR")
REFUSED = "argument '{}' has a DM_DEFAULT but is not optional"

DERIVED_NOTE = "Python and R derive '{}' and never take it, so the default is never used; drop it"
OUTPUT_NOTE = "an output takes no value from its caller, so there is nothing to default; drop it"
EXTENT_NOTE = "an extent cannot be optional; drop the DM_DEFAULT"
INOUT_NOTE = ("an in-out argument cannot be optional, so it cannot have a default either; drop "
              "the DM_DEFAULT")
OPTIONAL_NOTE = ("declare it `optional` and resolve it with M_DEFAULT_VAL (the C, Python and R "
                 "bindings come out the same, and Fortran callers get the default too), or drop "
                 "the DM_DEFAULT and state a recommended value in the prose")


def default(expression="2.0_real64", **others):
    return Directives(default=Default(expression), **others)


def scaled(*extra):
    """`scale_impl(n, values, ..., scaled)`: `n` is read off `values`."""
    return (
        b.integer("n", doc="count"),
        b.real("values", Intent.IN, "(n)", doc="values"),
        *extra,
        b.real("scaled", Intent.OUT, "(n)", doc="scaled values"),
    )


def impl_project(*arguments):
    impl = b.procedure("scale_impl", *arguments, meta=IMPL_META)
    return b.project(b.module("tox_demo_impl", impl, doc="demo", path="src/tox/tox_demo_impl.F90"))


def refusals(project):
    """The DM_DEFAULT refusals of the project as the pipeline validates it: (message, note)."""
    bag = DiagnosticBag()
    synthesis = synthesize_wrappers(project, CONVENTIONS, bag)
    analyse_project(synthesis.project, bag, CONVENTIONS)
    validate_project(synthesis.project, bag, CONVENTIONS)
    return [(d.message, d.note) for d in bag.errors if "has a DM_DEFAULT" in d.message]


# -- where it may go ---------------------------------------------------------------------------


class TestWhereItGoes:
    def test_on_an_optional_argument_it_is_accepted(self):
        factor = b.real("factor", doc="factor", optional=True, directives=default())
        assert refusals(impl_project(*scaled(factor))) == []

    def test_on_a_required_argument_it_is_refused(self):
        factor = b.real("factor", doc="factor", directives=default())
        assert [m for m, _ in refusals(impl_project(*scaled(factor)))] == [REFUSED.format("factor")]

    def test_an_array_keeps_its_own_refusal_only(self):
        weights = b.real("weights", Intent.IN, "(n)", doc="weights", directives=default())
        bag = DiagnosticBag()
        synthesis = synthesize_wrappers(impl_project(*scaled(weights)), CONVENTIONS, bag)
        analyse_project(synthesis.project, bag, CONVENTIONS)
        validate_project(synthesis.project, bag, CONVENTIONS)

        messages = [d.message for d in bag.errors]
        assert REFUSED.format("weights") not in messages
        assert "array argument 'weights' has a DM_DEFAULT" in messages

    def test_a_hand_written_export_is_checked(self):
        # M_EXPORT_C: no implementation, no synthesized wrapper, the declaration is the export's
        exported = b.procedure(
            "gd_scaled",
            *scaled(b.real("factor", doc="factor", directives=default())),
            b.ierr(),
        )
        project = b.project(b.module("gd_export", exported, doc="demo"))
        assert exported.is_exported  # the premise

        assert refusals(project) == [(REFUSED.format("factor"), OPTIONAL_NOTE)]

    def test_an_internal_helper_is_not_checked(self):
        # nothing is generated from it, so its documentation reaches no caller
        helper = b.procedure("scale_helper", *scaled(b.real("factor", doc="f", directives=default())),
                             meta=Meta())
        impl = b.procedure("scale_impl", *scaled(), meta=IMPL_META)
        project = b.project(b.module("tox_demo_impl", impl, helper, doc="demo",
                                     path="src/tox/tox_demo_impl.F90"))
        assert not helper.is_exported  # the premise

        assert refusals(project) == []


# -- the note: first match ---------------------------------------------------------------------


class TestTheNote:
    def test_an_extent_the_bindings_read_off_an_array(self):
        n = b.integer("n", doc="count", directives=default("3_int32"))
        arguments = (n,) + scaled()[1:]
        assert refusals(impl_project(*arguments)) == [(REFUSED.format("n"), DERIVED_NOTE.format("n"))]

    def test_a_mask_count(self):
        arguments = scaled(
            b.logical("values_mask", Intent.IN, "(n)", kind="c_bool", doc="which values"),
            b.integer("n_selected_values", doc="how many", directives=default("1_int32")),
        )
        assert refusals(impl_project(*arguments)) == [
            (REFUSED.format("n_selected_values"), DERIVED_NOTE.format("n_selected_values")),
        ]

    def test_a_value_a_producer_computes(self):
        work_size = b.integer("work_size", doc="work size", directives=default(
            "10_int32",
            output_from=OutputFrom("work_size", "work_size_of", "tox_demo", OutputFromMode.AUTO),
        ))
        arguments = scaled(work_size, b.real("tmp_work", Intent.OUT, "(work_size)", doc="work"))
        assert refusals(impl_project(*arguments)) == [
            (REFUSED.format("work_size"), DERIVED_NOTE.format("work_size")),
        ]

    def test_an_extent_only_an_output_carries(self):
        # the caller passes it, but an optional extent is refused too: drop the default
        n = b.integer("n", doc="count", directives=default("3_int32"))
        filled = b.real("filled", Intent.OUT, "(n)", doc="filled")
        arguments = (n, b.real("value", doc="value"), filled)
        assert refusals(impl_project(*arguments)) == [(REFUSED.format("n"), EXTENT_NOTE)]

    def test_a_just_info_size_of_a_work_array(self):
        # JUST_INFO: the caller looks it up and passes it, and it is an extent all the same
        work_size = b.integer("work_size", doc="work size", directives=default(
            "10_int32",
            output_from=OutputFrom("work_size", "work_size_of", "tox_demo", OutputFromMode.JUST_INFO),
        ))
        arguments = scaled(work_size, b.real("tmp_work", Intent.OUT, "(work_size)", doc="work"))
        assert refusals(impl_project(*arguments)) == [(REFUSED.format("work_size"), EXTENT_NOTE)]

    def test_a_just_info_value_that_sizes_nothing(self):
        # the caller passes it and it may be optional, so the default would be used
        threshold = b.integer("threshold", doc="threshold", directives=default(
            "10_int32",
            output_from=OutputFrom("threshold", "threshold_of", "tox_demo", OutputFromMode.JUST_INFO),
        ))
        assert refusals(impl_project(*scaled(threshold))) == [(REFUSED.format("threshold"),
                                                              OPTIONAL_NOTE)]

    def test_an_output(self):
        total = b.real("total", Intent.OUT, doc="total", directives=default("0.0_real64"))
        assert refusals(impl_project(*scaled(total))) == [(REFUSED.format("total"), OUTPUT_NOTE)]

    def test_an_in_out_argument(self):
        running = b.real("running", Intent.INOUT, doc="running", directives=default("1.0_real64"))
        assert refusals(impl_project(*scaled(running))) == [(REFUSED.format("running"), INOUT_NOTE)]

    @pytest.mark.parametrize("argument", [
        b.real("factor", doc="factor", directives=default()),
        b.logical("negate", doc="negate", directives=default(".false.")),
        b.character("separator", length="1", doc="separator", directives=default("','")),
    ], ids=["real", "logical", "character"])
    def test_any_other_input(self, argument):
        assert refusals(impl_project(*scaled(argument))) == [(REFUSED.format(argument.name),
                                                             OPTIONAL_NOTE)]


# -- from real source --------------------------------------------------------------------------

SOURCE = textwrap.dedent("""\
    #include <src/macros.h>

    !> Defaults on required arguments, in the shapes that reach the generated wrappers.
    module gd_default_impl
        use, intrinsic :: iso_fortran_env, only: real64, int32
        M_IMPLICIT_NONE
        private
        public :: gd_pair_impl, gd_split_impl, gd_fine_impl
        public :: MODE_A, MODE_B
        integer(int32), parameter :: MODE_A = 1_int32
            !! a
        integer(int32), parameter :: MODE_B = 2_int32
            !! b
    contains
        !> summary: scale, with a work array, so there are two tiers
        !| author: A Developer
        pure subroutine gd_pair_impl(n, values, factor, tmp_work, scaled)
            integer(int32), intent(in) :: n
                !! number of values
            real(real64), intent(in) :: values(n)
                !! the values
            real(real64), intent(in) :: factor
                !! the factor
                !! DM_DEFAULT(2.0_real64)
            real(real64), intent(out) :: tmp_work(n)
                !! scratch
            real(real64), intent(out) :: scaled(n)
                !! the scaled values
            tmp_work = values
            scaled = tmp_work * factor
        end subroutine gd_pair_impl

        !> summary: split mode, a required argument with a mode it is needed in and a default
        !| author: A Developer
        pure subroutine gd_split_impl(n, values, pattern_mode, gain, total)
            integer(int32), intent(in) :: n
                !! number of values
            real(real64), intent(in) :: values(n)
                !! the values
            integer(int32), intent(in) :: pattern_mode
                !! the mode
                !!
                !! | Mode | Value | Procedure |
                !! |------|-------|-----------|
                !! | a | [[gd_default_impl(module):MODE_A(variable)]] | gd_split_a |
                !! | b | [[gd_default_impl(module):MODE_B(variable)]] | gd_split_b |
            real(real64), intent(in) :: gain
                !! a gain only mode a needs
                !! DM_REQUIRED_IF_MODE(pattern_mode, gd_default_impl, MODE_A)
                !! DM_DEFAULT(0.5_real64)
            real(real64), intent(out) :: total
                !! their sum
            total = sum(values) * gain
        end subroutine gd_split_impl

        !> summary: an optional factor with a default, as it should be
        !| author: A Developer
        pure subroutine gd_fine_impl(n, values, factor, scaled)
            integer(int32), intent(in) :: n
                !! number of values
            real(real64), intent(in) :: values(n)
                !! the values
            real(real64), intent(in), optional :: factor
                !! the factor
                !! DM_DEFAULT(2.0_real64)
            real(real64), intent(out) :: scaled(n)
                !! the scaled values
            real(real64) :: factor_value
            M_DEFAULT_VAL(factor, factor_value, 2.0_real64)
            scaled = values * factor_value
        end subroutine gd_fine_impl
    end module gd_default_impl
    """)


def line_of(text, source=SOURCE, after=""):
    """The 1-based line of the first `text` below the first `after`."""
    lines = source.splitlines()
    start = next(i for i, line in enumerate(lines) if after in line)
    return next(i for i, line in enumerate(lines) if i >= start and text in line) + 1


@pytest.fixture(scope="module")
def real_source(tmp_path_factory):
    """`SOURCE` parsed by Ford, synthesized, and generated: (synthesis, diagnostics)."""
    src = tmp_path_factory.mktemp("dm_default") / "src" / "tox"
    src.mkdir(parents=True)
    (src / "gd_default_impl.F90").write_text(SOURCE)
    paths = Paths(root=REPO_ROOT, src_dir=src.parent)
    parsed = FordFrontend(paths, DiagnosticBag(), CONVENTIONS).parse()
    synthesis = synthesize_wrappers(parsed.project, CONVENTIONS, DiagnosticBag())
    result = generate(paths=paths, targets=("fortran",), parsed=parsed)
    return synthesis, result.diagnostics


def refused(diagnostics, name):
    return [d for d in diagnostics.errors if d.message == REFUSED.format(name)]


class TestRealSource:
    def test_reported_once_though_both_tiers_carry_it(self, real_source):
        synthesis, diagnostics = real_source
        tiers = [p for m in synthesis.project if m.name == "gd_default"
                 for p in m.procedures if p.name in ("gd_pair", "gd_pair_expert")]
        # the premise: both tiers were synthesized, each with the defaulted argument
        assert sorted(p.name for p in tiers) == ["gd_pair", "gd_pair_expert"]
        assert all(p.argument("factor").directives.has_default for p in tiers)

        [error] = refused(diagnostics, "factor")
        assert error.note == OPTIONAL_NOTE

    def test_at_the_directive_s_line(self, real_source):
        _, diagnostics = real_source
        [error] = refused(diagnostics, "factor")
        assert error.location.file.name == "gd_default_impl.F90"
        assert error.location.line == line_of("DM_DEFAULT(2.0_real64)", after="gd_pair_impl(")

    def test_a_split_mode_argument_is_refused(self, real_source):
        # the split wrapper takes it as optional because of its default and hands it to the
        # implementation's mandatory dummy: an omitted `gain` would crash a Fortran caller
        synthesis, diagnostics = real_source
        split = next(p for m in synthesis.project if m.name == "gd_default"
                     for p in m.procedures if p.name == "gd_split_a")
        assert split.argument("gain").optional  # the premise

        [error] = refused(diagnostics, "gain")
        assert error.location.line == line_of("DM_DEFAULT(0.5_real64)")
        assert error.note == OPTIONAL_NOTE

    def test_an_optional_argument_is_accepted(self, real_source):
        _, diagnostics = real_source
        assert all("gd_fine" not in str(d.context) for d in diagnostics.errors), diagnostics.render()
        assert len([d for d in diagnostics.errors if "has a DM_DEFAULT" in d.message]) == 2


def test_the_check_reads_the_declaration_not_the_wrappers():
    # a synthesized module is cloned from its implementation: checked there, not again
    factor = b.real("factor", doc="factor", directives=default())
    bag = DiagnosticBag()
    synthesis = synthesize_wrappers(impl_project(*scaled(factor)), CONVENTIONS, bag)
    check_defaults_are_on_optionals(synthesis.project, bag, CONVENTIONS)

    [error] = bag.errors
    assert "scale_impl" in error.render()
