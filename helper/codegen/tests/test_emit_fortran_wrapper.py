"""Emitting the generated validating wrapper module.

Text-level checks: that the wrapper sets ierr, validates each input the right way (an
extent by its size, a bounded real by its range, an unbounded real by finiteness alone),
bails on error, and calls the kernel. Compilation is exercised end to end elsewhere.
"""

from codegen.config import CONVENTIONS
from codegen.diagnostics import DiagnosticBag
from codegen.emit.fortran_wrapper import FortranWrapperEmitter
from codegen.ir.roles import analyse_project
from codegen.synthesize import synthesize_wrappers

from test_synthesize import alloc_kernel_module, ierr_kernel_module, kernel_module


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
                info[wrapper.name.lower()] = WrapperInfo(spec.kernel.name, spec.mode_fix)
    return FortranWrapperEmitter(project=synthesis.project).module(
        synthesis.project.module("tox_demo"), info
    )


def demo():
    from builders import project

    return emitted(project(kernel_module()))


def demo_alloc():
    from builders import project

    return emitted(project(alloc_kernel_module()))


def alloc_body():
    """Just the crunch_alloc subroutine text."""
    text = demo_alloc()
    start = text.index("subroutine crunch_alloc(")
    return text[start : text.index("end subroutine crunch_alloc")]


class TestModuleShell:
    def test_module_header_and_footer(self):
        text = demo()
        assert "module tox_demo" in text
        assert "end module tox_demo" in text
        assert "#include <src/macros.h>" in text
        assert "M_IMPLICIT_NONE" in text

    def test_uses_the_kernel_and_errors(self):
        text = demo()
        assert "use tox_demo_kernel, only: scale_vector_kernel" in text
        assert "use tox_errors, only:" in text
        assert "use, intrinsic :: iso_fortran_env, only:" in text

    def test_the_wrapper_is_public(self):
        assert "public :: scale_vector" in demo()

    def test_no_position_is_cleared_when_the_kernel_reports_nothing(self):
        # a kernel with no ierr cannot have packed a position, so neither the call nor the
        # import belongs here
        assert "clear_err_arg_pos" not in demo()

    def test_the_module_carries_a_generated_doc(self):
        # a non-empty module doc avoids a Ford "module has no documentation" warning
        assert "do not edit" in demo()

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
        from test_synthesize import count_extent_kernel_module

        text = emitted(project(count_extent_kernel_module()))
        assert (
            "call validate_in_range_int(n_selected, ierr, arg_pos=3_int32, min=0_int32, max=n_values)"
            in text
        )
        assert "validate_dimension_size(n_selected" not in text
        # a plain extent still uses the dimension check
        assert "call validate_dimension_size(n_values, ierr, arg_pos=2_int32)" in text

    def test_bails_before_calling_the_kernel(self):
        text = demo()
        assert "if (is_err(ierr)) return" in text
        assert text.index("if (is_err(ierr)) return") < text.index("call scale_vector_kernel(")


class TestKernelCall:
    def test_calls_the_kernel_with_the_arguments_unchanged(self):
        text = demo()
        assert "call scale_vector_kernel(&" in text
        assert "vector = vector" in text
        assert "factor = factor" in text

    def test_does_not_pass_ierr_to_the_kernel(self):
        # the kernel has no ierr; the call must not invent one
        call = demo().split("call scale_vector_kernel(&", 1)[1]
        assert "ierr = ierr" not in call


class TestTmpPermutation:
    def bodies(self):
        from builders import project
        from test_synthesize import tmp_perm_kernel_module

        text = emitted(project(tmp_perm_kernel_module()))
        start = text.index("subroutine rank_alloc(")
        return text[start : text.index("end subroutine rank_alloc")], text

    def test_the_tmp_permutation_is_allocated(self):
        body, _ = self.bodies()
        assert "M_ALLOCATE(tmp_column_perm(n))" in body
        assert "M_ALLOCATE(tmp_column(n))" in body

    def test_it_is_not_seeded_or_sorted_by_the_wrapper(self):
        # a tmp_ permutation is the kernel's own scratch: the kernel seeds and sorts it, so
        # the allocating wrapper must only allocate it -- no init_perm, no sort, no f42_sort
        body, text = self.bodies()
        assert "init_perm" not in body
        assert "sort_array_heapsort" not in body
        assert "use f42_sort" not in text


class TestSuffixCollisions:
    def bodies(self):
        from builders import project
        from test_synthesize import tmp_suffix_collision_kernel_module

        text = emitted(project(tmp_suffix_collision_kernel_module()))
        start = text.index("subroutine work_alloc(")
        return text[start : text.index("end subroutine work_alloc")], text

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
        assert "use f42_sort" not in text


class TestKernelThatDeclaresIerr:
    def demo_ierr(self):
        from builders import project

        return emitted(project(ierr_kernel_module()))

    def test_the_wrapper_passes_ierr_to_a_kernel_that_declares_one(self):
        # a kernel that propagates a sub-helper's error takes ierr; the wrapper must pass it
        call = self.demo_ierr().split("call risky_kernel(&", 1)[1]
        assert "ierr = ierr" in call

    def test_the_wrapper_still_validates_and_bails(self):
        text = self.demo_ierr()
        assert "call set_ok(ierr)" in text
        assert "if (is_err(ierr)) return" in text

    def test_the_argument_position_the_kernel_packed_is_cleared(self):
        # whatever position came back numbers the kernel's dummy list -- or a private
        # helper's, further down -- so it names nothing the wrapper's caller passed
        text = self.demo_ierr()
        after_call = text.split("call risky_kernel(&", 1)[1]
        assert "call clear_err_arg_pos(ierr)" in after_call

    def test_and_the_clear_is_imported(self):
        assert "clear_err_arg_pos" in self.demo_ierr().split("contains", 1)[0]


class TestModeSplitEmit:
    def emitted(self):
        from builders import project
        from test_synthesize import mode_split_kernel_module

        return emitted_with_info(project(mode_split_kernel_module()))

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

    def test_the_kernel_call_fixes_the_mode(self):
        dosage = self.dosage()
        assert "call detect_patterns_kernel(&" in dosage
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
        from test_synthesize import mask_count_kernel_module

        return emitted(project(mask_count_kernel_module()))

    def test_the_count_is_checked_against_the_mask(self):
        # the wrapper, not the kernel, guards count(mask) == n_selected -- so the kernel stays
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
        from test_synthesize import distance_matrix_kernel_module

        return emitted(project(distance_matrix_kernel_module()))

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
        text = demo_alloc()
        assert "use f42_sort, only: init_perm, sort_array_heapsort" in text
        # the recommend routine lives in the kernel module, imported with the kernel
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

    def test_it_calls_the_kernel_with_every_argument(self):
        body = alloc_body()
        call = body.split("call crunch_kernel(&", 1)[1]
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
    """A routine the wrapper runs before the kernel, which may handle the call itself."""

    def emitted_with(self, scope):
        from builders import project
        from test_synthesize import prologue_kernel_module

        return emitted(project(prologue_kernel_module(scope)))

    def body(self, text, name):
        start = text.index(f"subroutine {name}(")
        return text[start : text.index(f"end subroutine {name}")]

    def test_both_runs_it_in_each_wrapper(self):
        from codegen.ir.directives import PrologueScope

        text = self.emitted_with(PrologueScope.BOTH)

        for name in ("crunch", "crunch_alloc"):
            body = self.body(text, name)
            assert "call guard(&" in body, name
            assert "handled = handled" in body, name

    def test_it_returns_early_when_the_prologue_handled_the_call(self):
        from codegen.ir.directives import PrologueScope

        body = self.body(self.emitted_with(PrologueScope.BOTH), "crunch")

        assert "logical :: handled" in body
        assert "if (handled) return" in body
        # and the kernel is only reached afterwards
        assert body.index("if (handled) return") < body.index("call crunch_kernel(")

    def test_alloc_scope_runs_it_only_in_the_allocating_wrapper(self):
        from codegen.ir.directives import PrologueScope

        text = self.emitted_with(PrologueScope.ALLOC)

        assert "call guard(&" in self.body(text, "crunch_alloc")
        assert "call guard(&" not in self.body(text, "crunch")

    def test_expert_scope_runs_it_only_in_the_validating_wrapper(self):
        from codegen.ir.directives import PrologueScope

        text = self.emitted_with(PrologueScope.EXPERT)

        assert "call guard(&" in self.body(text, "crunch")
        assert "call guard(&" not in self.body(text, "crunch_alloc")

    def test_it_runs_before_the_allocation_it_may_make_unnecessary(self):
        from codegen.ir.directives import PrologueScope

        body = self.body(self.emitted_with(PrologueScope.ALLOC), "crunch_alloc")

        assert body.index("call guard(&") < body.index("M_ALLOCATE(tmp_scratch(n))")

    def test_the_prologue_is_imported(self):
        from codegen.ir.directives import PrologueScope

        header = self.emitted_with(PrologueScope.BOTH).split("contains", 1)[0]
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

        from test_synthesize import optional_producer_input_kernel_module

        return emitted(project(optional_producer_input_kernel_module()))

    def test_the_default_is_resolved_into_a_local(self):
        text = self.text()
        assert "logical :: exact_value" in text
        assert "M_DEFAULT_VAL(exact, exact_value, .false.)" in text

    def test_the_recommend_call_takes_the_local(self):
        assert "exact = exact_value" in self.text()

    def test_the_kernel_still_takes_the_argument_itself(self):
        # the kernel declares it optional, so it is forwarded as it came
        assert "exact = exact,&" in self.text()
