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

from test_synthesize import kernel_module


def emitted(module_source):
    project = synthesize_wrappers(module_source).project
    analyse_project(project, DiagnosticBag(), CONVENTIONS)
    generated = project.module("tox_demo")
    return FortranWrapperEmitter().module(generated)


def demo():
    from builders import project

    return emitted(project(kernel_module()))


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
