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

from test_synthesize import alloc_kernel_module, kernel_module


def emitted(module_source):
    project = synthesize_wrappers(module_source).project
    analyse_project(project, DiagnosticBag(), CONVENTIONS)
    generated = project.module("tox_demo")
    return FortranWrapperEmitter().module(generated)


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


class TestAllocatingWrapper:
    def test_the_module_imports_what_the_alloc_needs(self):
        text = demo_alloc()
        assert "use f42_utils, only: init_perm, sort_array_heapsort" in text
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
