"""Emitting the generated validating wrapper module.

Text-level checks: that the wrapper sets ierr, validates each input the right way (an
extent by its size, a bounded real by its range, an unbounded real by finiteness alone),
bails on error, and calls the implementation. Compilation is exercised end to end elsewhere.
"""

from codegen.config import CONVENTIONS
from codegen.diagnostics import DiagnosticBag
from codegen.emit.fortran_wrapper import FortranWrapperEmitter
from codegen.ir.roles import analyse_project
from codegen.ir.directives import Directives, OutputFrom, OutputFromMode
from codegen.ir.types import Intent
from codegen.synthesize import synthesize_wrappers

from builders import C_BINDING

from test_synthesize import alloc_impl_module, ierr_impl_module, impl_module


def emitted(module_source):
    project = synthesize_wrappers(module_source).project
    analyse_project(project, DiagnosticBag(), CONVENTIONS)
    generated = project.module("tox_demo")
    return FortranWrapperEmitter(project=project).module(generated)


def emitted_with_info(module_source):
    """Like `emitted`, but threading the per-wrapper info the mode-split path needs."""
    from codegen.emit.fortran_wrapper import WrapperInfo

    synthesis = synthesize_wrappers(module_source)
    analyse_project(synthesis.project, DiagnosticBag(), CONVENTIONS)
    info = {}
    for spec in synthesis.specs:
        for wrapper in (spec.validating, spec.allocating):
            if wrapper is not None:
                info[wrapper.name.lower()] = WrapperInfo(spec.impl.name, spec.mode_fix)
    return FortranWrapperEmitter(project=synthesis.project).module(
        synthesis.project.module("tox_demo"), info
    )


def demo():
    from builders import project

    return emitted(project(impl_module()))


def paired_demo():
    from builders import project

    return emitted(project(alloc_impl_module()))


def alloc_body():
    """Just the crunch subroutine text."""
    text = paired_demo()
    start = text.index("subroutine crunch(")
    return text[start : text.index("end subroutine crunch")]


class TestPurity:
    """A wrapper is `pure` exactly when everything it calls is.

    So a caller who wants one inside `do concurrent` does not have to reach past the
    wrapper to the implementation, giving up validation to get it.
    """

    def module_of(self, impl, **kwargs):
        from builders import Meta, integer, module, procedure, project, real

        return project(
            module(
                "tox_demo_impl",
                procedure(
                    impl,
                    real("values", Intent.INOUT, "(n)", doc="the data"),
                    integer("n", Intent.IN, doc="length"),
                    meta=Meta(summary="Scale", author="AUTHOR"),
                    **kwargs,
                ),
            )
        )

    def test_a_pure_implementation_gives_a_pure_wrapper(self):
        text = emitted(self.module_of("scale_impl", is_pure=True))

        assert "pure subroutine scale(" in text

    def test_an_impure_one_does_not(self):
        text = emitted(self.module_of("scale_impl"))

        assert "pure subroutine scale(" not in text
        assert "subroutine scale(" in text

    def test_an_impure_producer_makes_the_wrapper_impure(self):
        # the wrapper calls the recommend routine itself, so its purity counts too --
        # the implementation being pure is necessary and not sufficient
        from builders import Meta, integer, module, procedure, project, real

        source = project(
            module(
                "tox_demo_impl",
                procedure(
                    "work_size",
                    integer("n", Intent.IN, doc="length"),
                    integer("wsize", Intent.OUT, doc="recommended size"),
                    meta=C_BINDING,
                ),
                procedure(
                    "crunch_impl",
                    real("values", Intent.IN, "(n)", doc="the data"),
                    integer("n", Intent.IN, doc="length"),
                    real("tmp_work", Intent.OUT, "(wsize)", doc="work buffer"),
                    integer(
                        "wsize",
                        Intent.IN,
                        directives=Directives(
                            output_from=OutputFrom(
                                "wsize", "work_size", "tox_demo_impl", OutputFromMode.AUTO
                            )
                        ),
                        doc="work size",
                    ),
                    meta=Meta(summary="Crunch", author="AUTHOR"),
                    is_pure=True,
                ),
            )
        )

        text = emitted(source)

        assert "subroutine crunch(" in text
        assert "pure subroutine crunch(" not in text


class TestModuleShell:
    def test_module_header_and_footer(self):
        text = demo()
        assert "module tox_demo" in text
        assert "end module tox_demo" in text
        assert "#include <src/macros.h>" in text
        assert "M_IMPLICIT_NONE" in text

    def test_uses_the_impl_and_errors(self):
        text = demo()
        assert "use tox_demo_impl, only: scale_vector_impl" in text
        assert "use tox_errors, only:" in text
        assert "use, intrinsic :: iso_fortran_env, only:" in text

    def test_the_wrapper_is_public(self):
        assert "public :: scale_vector" in demo()

    def test_no_position_is_cleared_when_the_impl_reports_nothing(self):
        # an implementation with no ierr cannot have packed a position, so neither the call nor the
        # import belongs here
        assert "clear_err_arg_pos" not in demo()

    def test_the_module_carries_a_generated_doc(self):
        # a non-empty module doc avoids a Ford "module has no documentation" warning
        assert "do not edit" in demo()


class TestModuleDocumentation:
    """A generated module is the published API, so it carries the author's own module doc."""

    def emitted(self):
        from builders import module, project

        impl = impl_module()
        documented = module(
            "tox_demo_impl",
            *impl.procedures,
            doc=["Scaling vectors in place.", "", "A second paragraph."],
            path="src/tox/tox_demo_impl.F90",
        )
        return emitted(project(documented))

    def test_the_author_s_documentation_is_carried_verbatim(self):
        text = self.emitted()
        assert "!> Scaling vectors in place." in text
        assert "!| A second paragraph." in text

    def test_the_author_s_first_line_is_the_summary(self):
        """No `summary:` tag of the emitter's own -- it would displace the author's opening
        line everywhere Ford shows a one-liner, to say what the note below already says."""
        header = self.emitted().split("module tox_demo")[0].splitlines()
        assert not any(line.startswith("!> summary:") for line in header)
        assert header[2] == "!> Scaling vectors in place."

    def test_the_generated_note_comes_last(self):
        header = self.emitted().split("module tox_demo")[0].splitlines()
        assert header[-1] == (
            "!| Generated from [[tox_demo_impl(module)]]; do not edit -- regenerate instead."
        )

    def test_an_undocumented_module_still_gets_the_note(self):
        # `impl_module` has no module doc, so the note is the whole of it -- and has to open
        # with `!>` rather than the `!|` continuation, or Ford reads no documentation at all
        assert demo().splitlines()[2].startswith("!> Generated from ")

    def test_the_synthesised_ierr_is_documented(self):
        # ... and a doc on the synthesised ierr avoids the argument-level warning
        assert "Error code; zero on success" in demo()


class TestValidation:
    def test_sets_ierr_ok_first(self):
        assert "call set_ok(ierr)" in demo()

    def test_an_extent_is_validated_by_its_size(self):
        # n is the extent of vector(n); it is the 2nd dummy
        assert "call validate_dimension_size(n, ierr, arg_pos=2_int32)" in demo()

    def test_a_bounded_real_is_range_checked(self):
        # factor carries DM_MIN(0.0_real64); it is the 3rd dummy
        assert (
            "call validate_in_range_real(factor, ierr, arg_pos=3_int32, min=0.0_real64)"
            in demo()
        )

    def test_an_unbounded_real_array_is_finiteness_checked(self):
        # vector has no bound, but finiteness is the default contract; count is its extent
        assert (
            "call validate_all_in_range_real(vector, n, ierr, arg_pos=1_int32)" in demo()
        )

    def test_an_extent_with_a_range_uses_the_range_not_the_dimension_check(self):
        # n_selected is the extent of indices, but its DM_MIN(0)/DM_MAX permits an empty
        # selection, which validate_dimension_size would reject
        from builders import project
        from test_synthesize import count_extent_impl_module

        text = emitted(project(count_extent_impl_module()))
        assert (
            "call validate_in_range_int(n_selected, ierr, arg_pos=3_int32, min=0_int32, max=n_values)"
            in text
        )
        assert "validate_dimension_size(n_selected" not in text
        # a plain extent still uses the dimension check
        assert "call validate_dimension_size(n_values, ierr, arg_pos=2_int32)" in text

    def test_bails_before_calling_the_impl(self):
        text = demo()
        assert "if (is_err(ierr)) return" in text
        assert text.index("if (is_err(ierr)) return") < text.index("call scale_vector_impl(")


class TestImplCall:
    def test_calls_the_impl_with_the_arguments_unchanged(self):
        text = demo()
        assert "call scale_vector_impl(&" in text
        assert "vector = vector" in text
        assert "factor = factor" in text

    def test_does_not_pass_ierr_to_the_impl(self):
        # the implementation has no ierr; the call must not invent one
        call = demo().split("call scale_vector_impl(&", 1)[1]
        assert "ierr = ierr" not in call


class TestTmpPermutation:
    def bodies(self):
        from builders import project
        from test_synthesize import tmp_perm_impl_module

        text = emitted(project(tmp_perm_impl_module()))
        start = text.index("subroutine rank(")
        return text[start : text.index("end subroutine rank")], text

    def test_the_tmp_permutation_is_allocated(self):
        body, _ = self.bodies()
        assert "M_ALLOCATE(tmp_column_perm(n))" in body
        assert "M_ALLOCATE(tmp_column(n))" in body

    def test_it_is_not_seeded_or_sorted_by_the_wrapper(self):
        # a tmp_ permutation is the implementation's own scratch: the implementation seeds and sorts it, so
        # the allocating wrapper must only allocate it -- no init_perm, no sort, no f42_sort
        body, text = self.bodies()
        assert "init_perm" not in body
        assert "sort_array_heapsort" not in body
        assert "use f42_sort_impl" not in text


class TestSuffixCollisions:
    def bodies(self):
        from builders import project
        from test_synthesize import tmp_suffix_collision_impl_module

        text = emitted(project(tmp_suffix_collision_impl_module()))
        start = text.index("subroutine work(")
        return text[start : text.index("end subroutine work")], text

    def test_a_tmp_shape_is_allocated_not_derived(self):
        body, _ = self.bodies()
        assert "M_ALLOCATE(tmp_shape(2))" in body
        # dropped from the wrapper's signature (a work array), not carried as a dummy
        signature = body.split(")", 1)[0]
        assert "tmp_shape" not in signature

    def test_a_tmp_mask_is_allocated_not_derived(self):
        body, _ = self.bodies()
        assert "M_ALLOCATE(tmp_selection_mask(n))" in body

    def test_no_seeding_or_perm_helpers(self):
        body, text = self.bodies()
        assert "init_perm" not in body
        assert "use f42_sort_impl" not in text


class TestImplThatDeclaresIerr:
    def demo_ierr(self):
        from builders import project

        return emitted(project(ierr_impl_module()))

    def test_the_wrapper_passes_ierr_to_a_impl_that_declares_one(self):
        # an implementation that propagates a sub-helper's error takes ierr; the wrapper must pass it
        call = self.demo_ierr().split("call risky_impl(&", 1)[1]
        assert "ierr = ierr" in call

    def test_the_wrapper_still_validates_and_bails(self):
        text = self.demo_ierr()
        assert "call set_ok(ierr)" in text
        assert "if (is_err(ierr)) return" in text

    def test_the_argument_position_the_impl_packed_is_cleared(self):
        # whatever position came back numbers the implementation's dummy list -- or a private
        # helper's, further down -- so it names nothing the wrapper's caller passed
        text = self.demo_ierr()
        after_call = text.split("call risky_impl(&", 1)[1]
        assert "call clear_err_arg_pos(ierr)" in after_call

    def test_and_the_clear_is_imported(self):
        assert "clear_err_arg_pos" in self.demo_ierr().split("contains", 1)[0]


class TestModeSplitEmit:
    def emitted(self):
        from builders import project
        from test_synthesize import mode_split_impl_module

        return emitted_with_info(project(mode_split_impl_module()))

    def dosage(self):
        text = self.emitted()
        start = text.index("subroutine detect_dosage_effect(")
        return text[start : text.index("end subroutine detect_dosage_effect")]

    def subfunc(self):
        text = self.emitted()
        start = text.index("subroutine detect_subfunctionalization(")
        return text[start : text.index("end subroutine detect_subfunctionalization")]

    def test_each_mode_becomes_a_subroutine(self):
        text = self.emitted()
        assert "subroutine detect_dosage_effect(" in text
        assert "subroutine detect_subfunctionalization(" in text
        # no single runtime-mode procedure
        assert "subroutine detect_patterns(" not in text

    def test_the_impl_call_fixes_the_mode(self):
        dosage = self.dosage()
        assert "call detect_patterns_impl(&" in dosage
        assert "pattern_mode = MODE_DOSAGE" in dosage

    def test_the_required_argument_is_passed_in_its_mode(self):
        dosage = self.dosage()
        assert "threshold = threshold" in dosage

    def test_the_argument_of_another_mode_is_omitted(self):
        # detect_subfunctionalization takes no threshold, so its call must not name it
        subfunc = self.subfunc()
        assert "pattern_mode = MODE_SUBFUNC" in subfunc
        assert "threshold" not in subfunc

    def test_the_fixed_parameter_is_imported(self):
        header = self.emitted().split("contains", 1)[0]
        assert "MODE_DOSAGE" in header
        assert "MODE_SUBFUNC" in header


class TestMaskCountConvention:
    def emitted(self):
        from builders import project
        from test_synthesize import mask_count_impl_module

        return emitted(project(mask_count_impl_module()))

    def test_the_count_is_checked_against_the_mask(self):
        # the wrapper, not the implementation, guards count(mask) == n_selected -- so the implementation stays
        # pure and needs no ierr for it
        text = self.emitted()
        assert (
            "if (count(vecs_selection_mask, kind=int32) /= n_selected_vecs) "
            "call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)" in text
        )

    def test_it_imports_the_error_symbols(self):
        header = self.emitted().split("contains", 1)[0]
        assert "set_err_once" in header
        assert "ERR_INVALID_INPUT" in header


class TestDistanceMatrixConvention:
    def emitted(self):
        from builders import project
        from test_synthesize import distance_matrix_impl_module

        return emitted(project(distance_matrix_impl_module()))

    def test_the_matrix_is_structurally_validated(self):
        assert (
            "call validate_distance_matrix(distances, n_points, ierr, arg_pos=1_int32)"
            in self.emitted()
        )

    def test_the_matrix_opts_out_of_the_finiteness_check(self):
        # validate_distance_matrix subsumes the real finiteness check, which would otherwise
        # run on the same argument
        assert "validate_all_in_range_real(distances" not in self.emitted()


class TestAllocatingWrapper:
    def test_the_module_imports_what_the_alloc_needs(self):
        text = paired_demo()
        assert "use f42_sort_impl, only: init_perm, sort_array_heapsort" in text
        # the recommend routine lives in the implementation module, imported with the implementation
        assert "work_size" in text.split("contains")[0]
        assert "ERR_ALLOC_FAIL" in text

    def test_the_taken_over_arguments_become_allocatable_locals(self):
        body = alloc_body()
        assert "integer(int32), dimension(:), allocatable :: values_perm" in body
        assert "real(real64), dimension(:), allocatable :: tmp_scratch" in body
        assert "real(real64), dimension(:), allocatable :: tmp_work" in body
        assert "integer(int32) :: wsize" in body  # the computed size is a scalar local

    def test_the_recommend_routine_is_called_into_the_size(self):
        body = alloc_body()
        assert "call work_size(&" in body
        assert "n = n" in body
        assert "wsize = wsize" in body

    def test_the_work_arrays_are_allocated(self):
        body = alloc_body()
        assert "M_ALLOCATE(tmp_scratch(n))" in body
        assert "M_ALLOCATE(tmp_work(wsize))" in body
        assert "M_ALLOCATE(values_perm(n))" in body

    def test_the_permutation_is_seeded_and_sorted(self):
        body = alloc_body()
        assert "call init_perm(values_perm)" in body
        assert "call sort_array_heapsort(values, values_perm)" in body

    def test_it_calls_the_impl_with_every_argument(self):
        body = alloc_body()
        call = body.split("call crunch_impl(&", 1)[1]
        for name in ("values", "n", "values_perm", "tmp_scratch", "tmp_work", "wsize"):
            assert f"{name} = {name}" in call
        assert "ierr = ierr" not in call

    def test_recommend_before_allocate_before_sort(self):
        body = alloc_body()
        assert body.index("call work_size(") < body.index("M_ALLOCATE(tmp_work(wsize))")
        assert body.index("M_ALLOCATE(values_perm(n))") < body.index(
            "call sort_array_heapsort(values, values_perm)"
        )


class TestPrologue:
    """A routine the wrapper runs before the implementation, which may handle the call itself."""

    def emitted_with(self, scratch=False):
        from builders import project
        from test_synthesize import prologue_impl_module

        return emitted(project(prologue_impl_module(self._guard(scratch))))

    @staticmethod
    def _guard(scratch):
        """The prologue's dummies; with `scratch` it also takes the implementation's work array."""
        if not scratch:
            return None
        from codegen.ir.types import Intent

        from builders import ierr, integer, logical, real

        return (
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length"),
            real("tmp_scratch", Intent.OUT, "(n)", doc="scratch it fills itself"),
            real("result", Intent.OUT, "(n)", doc="the answer"),
            logical("handled", Intent.OUT, doc="whether the call was dealt with"),
            ierr(),
        )

    def body(self, text, name):
        start = text.index(f"subroutine {name}(")
        return text[start : text.index(f"end subroutine {name}")]

    def test_it_runs_only_in_the_allocating_wrapper(self):
        # a prologue is the sugar that tier adds; the validating wrapper is the tier that
        # lets a caller do the same preparation themselves, so it must not run there
        text = self.emitted_with()

        assert "call guard(&" in self.body(text, "crunch")
        assert "call guard(&" not in self.body(text, "crunch_expert")
        assert "logical :: handled" not in self.body(text, "crunch_expert")

    def test_it_returns_early_when_the_prologue_handled_the_call(self):
        body = self.body(self.emitted_with(), "crunch")

        assert "logical :: handled" in body
        assert "if (handled) return" in body
        # and the implementation is only reached afterwards
        assert body.index("if (handled) return") < body.index("call crunch_impl(")

    def test_it_runs_below_the_work_arrays_it_prepares(self):
        # preparing them is what it is for, so there is nothing to hand it any earlier
        body = self.body(self.emitted_with(scratch=True), "crunch")

        assert body.index("M_ALLOCATE(tmp_scratch(n))") < body.index("call guard(&")
        assert body.index("call guard(&") < body.index("call crunch_impl(")

    def test_it_runs_below_them_even_when_it_takes_none(self):
        # the placement does not depend on what this prologue happens to take: it is where
        # the tier's preparation belongs, so it is the same for every prologue
        body = self.body(self.emitted_with(), "crunch")

        assert body.index("M_ALLOCATE(tmp_scratch(n))") < body.index("call guard(&")

    def call_to(self, body, name):
        """Just the `call <name>(...)` statement -- the whole body proves nothing.

        Every actual is spelled `x = x`, so an assertion over the whole body is satisfied by
        the *impl* call two lines below and holds even when the prologue call is empty.
        """
        start = body.index(f"call {name}(")
        return body[start : body.index(")", body.index("&\n", start))]

    def test_it_is_handed_the_work_array_it_asked_for(self):

        body = self.body(self.emitted_with(scratch=True), "crunch")

        assert "tmp_scratch = tmp_scratch" in self.call_to(body, "guard")

    def test_an_ordinary_prologue_is_handed_only_what_it_asked_for(self):
        # the negative control for the assertion above: the plain prologue takes no work
        # array, so its call must not mention one even though the body around it does

        body = self.body(self.emitted_with(), "crunch")

        assert "tmp_scratch" not in self.call_to(body, "guard")
        assert "tmp_scratch" in body

    def test_the_prologue_is_imported(self):

        header = self.emitted_with().split("contains", 1)[0]
        assert "guard" in header


class TestReexportModule:
    def text(self):
        from builders import project

        from test_synthesize import split_family_modules

        synthesis = synthesize_wrappers(project(*split_family_modules()))
        analyse_project(synthesis.project, DiagnosticBag(), CONVENTIONS)
        return FortranWrapperEmitter(project=synthesis.project).module(
            synthesis.project.module("tox_demo")
        )

    def test_it_is_a_module_of_nothing_but_uses(self):
        text = self.text()
        assert "module tox_demo" in text
        assert "use tox_demo_left" in text
        assert "use tox_demo_right" in text
        assert "end module tox_demo" in text
        assert "contains" not in text

    def test_nothing_is_hidden(self):
        # re-exporting is the point: no `only` list to filter what a child made public,
        # and no `private` to hide it again
        text = self.text()
        assert "only:" not in text
        assert "private" not in text

    def test_it_declares_nothing_so_needs_no_macros(self):
        text = self.text()
        assert "#include" not in text
        assert "M_IMPLICIT_NONE" not in text

    def test_it_carries_a_generated_doc(self):
        assert "do not edit" in self.text()


class TestMaterialisedProducerInput:
    """A recommend routine takes by value what the wrapper takes optionally.

    Forwarding the optional straight into a mandatory dummy is not a Fortran program: when
    the caller omits it, the routine reads an absent argument. The wrapper has to resolve
    the documented default into a local first.
    """

    def text(self):
        from builders import project

        from test_synthesize import optional_producer_input_impl_module

        return emitted(project(optional_producer_input_impl_module()))

    def test_the_default_is_resolved_into_a_local(self):
        text = self.text()
        assert "logical :: exact_value" in text
        assert "M_DEFAULT_VAL(exact, exact_value, .false.)" in text

    def test_the_recommend_call_takes_the_local(self):
        assert "exact = exact_value" in self.text()

    def test_the_impl_still_takes_the_argument_itself(self):
        # the implementation declares it optional, so it is forwarded as it came
        assert "exact = exact,&" in self.text()


class TestModeMembership:
    """A runtime-mode argument is checked against the values its own table names."""

    def text(self, module_source):
        from builders import project

        return emitted(project(module_source))

    def test_the_wrapper_rejects_a_value_the_table_does_not_name(self):
        from test_synthesize import mode_impl_module

        text = self.text(mode_impl_module())
        assert ("if (mode /= MODE_FAST .and. mode /= MODE_EXACT) "
                "call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=2_int32)") in text

    def test_the_parameters_it_compares_against_are_imported(self):
        from test_synthesize import mode_impl_module

        text = self.text(mode_impl_module()).split("contains", 1)[0]
        assert "MODE_FAST" in text and "MODE_EXACT" in text

    def test_a_mode_split_wrapper_has_no_mode_to_check(self):
        from builders import project

        from test_synthesize import mode_split_impl_module

        # it *is* its mode; there is no argument left to compare
        assert "/= MODE_" not in emitted_with_info(project(mode_split_impl_module()))


class TestAPrologueThatBuildsThePermutation:
    """`<base>_perm` is seeded and heapsorted by the allocating wrapper -- unless the prologue
    says it builds that permutation itself, which is how a family gets a non-default ordering
    without giving up the expert tier's ability to take one from the caller."""

    def emitted_with(self, perm_intent):
        from codegen.ir.types import Intent

        from builders import ierr, integer, logical, project, real
        from test_synthesize import prologue_impl_module

        guard = (
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length"),
            logical("handled", Intent.OUT, doc="dealt with"),
            ierr(),
        )
        if perm_intent is not None:
            guard = guard[:2] + (
                integer("values_perm", perm_intent, "(n)", doc="ordered how this family wants"),
            ) + guard[2:]
        impl = (
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length"),
            integer("values_perm", Intent.OUT, "(n)", doc="values, ordered"),
            real("result", Intent.OUT, "(n)", doc="the answer"),
        )
        return emitted(project(prologue_impl_module(guard, impl)))

    def body(self, text):
        start = text.index("subroutine crunch(")
        return text[start : text.index("end subroutine crunch")]

    def test_by_default_the_wrapper_seeds_and_sorts_it(self):
        text = self.emitted_with(perm_intent=None)

        body = self.body(text)
        assert "M_ALLOCATE(values_perm(n))" in body
        assert "call init_perm(values_perm)" in body
        assert "call sort_array_heapsort(values, values_perm)" in body
        assert "use f42_sort_impl" in text

    def test_a_prologue_that_declares_it_intent_out_takes_it_over(self):
        from codegen.ir.types import Intent

        text = self.emitted_with(perm_intent=Intent.OUT)

        body = self.body(text)
        # still allocated for it, but neither seeded nor sorted -- and it is handed over
        assert "M_ALLOCATE(values_perm(n))" in body
        assert "init_perm" not in body
        assert "sort_array_heapsort" not in body
        assert "values_perm = values_perm" in body
        # ... so the helpers are not imported either
        assert "use f42_sort_impl" not in text

    def test_intent_inout_refines_the_default_instead_of_replacing_it(self):
        # an in-out prologue dummy reads what it is given, so the wrapper still builds it
        from codegen.ir.types import Intent

        text = self.emitted_with(perm_intent=Intent.INOUT)

        body = self.body(text)
        assert "call sort_array_heapsort(values, values_perm)" in body
        assert body.index("call sort_array_heapsort(values, values_perm)") < body.index("call guard(&")


class TestAPrologueThatProducesAImplInput:
    """Prologue `intent(out)` over impl `intent(in)`: the value never leaves the wrapper,
    so the allocating wrapper drops it from its signature and declares it as a local."""

    def emitted_with(self, extra_impl, extra_guard):
        from codegen.ir.types import Intent

        from builders import ierr, integer, logical, project, real
        from test_synthesize import prologue_impl_module

        guard = (
            integer("n", Intent.IN, doc="length"),
            extra_guard,
            real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
            logical("handled", Intent.OUT, doc="dealt with"),
            ierr(),
        )
        impl = (
            integer("n", Intent.IN, doc="length"),
            extra_impl,
            real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
            real("result", Intent.OUT, "(n)", doc="the answer"),
        )
        return emitted(project(prologue_impl_module(guard, impl)))

    def bodies(self, text):
        alloc = text[text.index("subroutine crunch(") : text.index("end subroutine crunch")]
        foo = text[text.index("subroutine crunch_expert(") : text.index("end subroutine crunch_expert")]
        return alloc, foo

    def test_a_scalar_becomes_a_plain_local(self):
        from codegen.ir.types import Intent

        from builders import integer

        alloc, foo = self.bodies(self.emitted_with(
            integer("room", Intent.IN, doc="read by the implementation"),
            integer("room", Intent.OUT, doc="written by the prologue"),
        ))

        assert "integer(int32) :: room" in alloc     # a local, not allocated
        assert "M_ALLOCATE(room" not in alloc
        assert "room" not in alloc[: alloc.index(")")]  # gone from the signature
        # ... but the expert tier still takes it, because it has no prologue
        assert "intent(in) :: room" in foo

    def test_an_array_becomes_an_allocated_local(self):
        from codegen.ir.types import Intent

        from builders import integer

        alloc, foo = self.bodies(self.emitted_with(
            integer("ranking", Intent.IN, "(n)", doc="read by the implementation"),
            integer("ranking", Intent.OUT, "(n)", doc="built by the prologue"),
        ))

        assert "integer(int32), dimension(:), allocatable :: ranking" in alloc
        assert "M_ALLOCATE(ranking(n))" in alloc
        assert "ranking" not in alloc[: alloc.index(")")]
        assert "intent(in) :: ranking" in foo

    def test_it_is_allocated_before_the_prologue_that_fills_it(self):
        from codegen.ir.types import Intent

        from builders import integer

        alloc, _ = self.bodies(self.emitted_with(
            integer("ranking", Intent.IN, "(n)", doc="read by the implementation"),
            integer("ranking", Intent.OUT, "(n)", doc="built by the prologue"),
        ))

        assert alloc.index("M_ALLOCATE(ranking(n))") < alloc.index("call guard(&")
        assert alloc.index("call guard(&") < alloc.index("call crunch_impl(")


class TestValidationCanBeCompiledOut:
    """`NO_INPUT_VALIDATION` removes the checks for a caller who has already established
    that the inputs are good. What must survive it is anything that is not a check."""

    def body(self):
        return alloc_body()

    def test_the_checks_sit_behind_the_guard(self):
        body = self.body()

        assert "#ifndef NO_INPUT_VALIDATION" in body
        assert "#endif" in body
        guarded = body[body.index("#ifndef NO_INPUT_VALIDATION") : body.index("#endif")]
        assert "validate_" in guarded
        assert "if (is_err(ierr)) return" in guarded

    def test_set_ok_stays_outside_it(self):
        # not a check: it is what leaves ierr defined when nothing goes wrong, and a build
        # without validation still reports the runtime errors an implementation raises
        body = self.body()

        assert body.index("call set_ok(ierr)") < body.index("#ifndef NO_INPUT_VALIDATION")

    def test_the_impl_call_stays_outside_it(self):
        body = self.body()

        assert body.index("#endif") < body.index("call crunch_impl(")

    def test_the_directives_are_at_column_zero(self):
        # gfortran's preprocessor rejects an indented directive outright, and a body is
        # rendered on its own writer before being nested into its module's
        text = paired_demo()

        for line in text.splitlines():
            if line.lstrip().startswith("#"):
                assert line == line.lstrip(), repr(line)

    def test_a_wrapper_with_nothing_to_check_emits_no_guard(self):
        from builders import project
        from test_synthesize import impl_module

        # an empty #ifndef/#endif pair would be noise, and there is no is_err to test either
        text = emitted(project(impl_module()))
        for procedure in ("scale_vector",):
            assert "call set_ok(ierr)" in text


class TestAPrologueWithArgumentsOfItsOwn:
    """What a prologue derives *from* is the allocating tier's vocabulary: a threshold comes
    from a percentile, and the implementation takes the threshold."""

    def built(self):
        from codegen.ir.types import Intent

        from builders import ierr, integer, logical, project, real
        from test_synthesize import prologue_impl_module

        guard = (
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length"),
            real("percentile", Intent.IN, doc="where to put the threshold"),
            real("threshold", Intent.OUT, doc="derived here"),
            real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
            logical("handled", Intent.OUT, doc="dealt with"),
            ierr(),
        )
        impl = (
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length"),
            real("threshold", Intent.IN, doc="what counts as an outlier"),
            real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
            real("result", Intent.OUT, "(n)", doc="the answer"),
        )
        return emitted(project(prologue_impl_module(guard, impl)))

    def body(self, text, name):
        return text[text.index(f"subroutine {name}(") : text.index(f"end subroutine {name}")]

    def test_the_allocating_wrapper_takes_it(self):
        alloc = self.body(self.built(), "crunch")

        assert "real(real64), intent(in) :: percentile" in alloc
        assert "percentile = percentile" in alloc

    def test_the_expert_wrapper_does_not(self):
        # it takes the derived value directly, so what it was derived from means nothing here
        expert = self.body(self.built(), "crunch_expert")

        assert "percentile" not in expert
        assert "intent(in) :: threshold" in expert

    def test_the_derived_value_is_a_local_of_the_allocating_wrapper(self):
        alloc = self.body(self.built(), "crunch")

        assert "real(real64) :: threshold" in alloc      # a local, not a dummy
        assert "threshold" not in alloc[: alloc.index(")")]
        assert alloc.index("call guard(&") < alloc.index("call crunch_impl(")

    def test_the_prologue_argument_comes_after_the_implementations_and_before_ierr(self):
        alloc = self.body(self.built(), "crunch")
        signature = alloc[: alloc.index(")")]

        assert signature.index("result") < signature.index("percentile")
        assert signature.index("percentile") < signature.index("ierr")

    def test_a_prologue_argument_is_validated_like_any_other(self):
        # it reaches the wrapper's own dummy list, so the range checks apply to it too
        alloc = self.body(self.built(), "crunch")

        assert "call validate_in_range_real(percentile, ierr" in alloc
