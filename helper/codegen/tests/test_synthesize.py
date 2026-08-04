"""Synthesising the validating wrapper from a kernel.

A kernel carries no `ierr` and no validation; the wrapper the generator builds from it
adds both and is the exported entry point. These tests hand-build a kernel and check the
injected wrapper, without going near Ford.
"""

from codegen.ir.directives import Directives, Minimum
from codegen.ir.entities import Meta
from codegen.ir.types import Intent
from codegen.synthesize import synthesize_wrappers

from builders import integer, module, procedure, project, real


def kernel_module():
    """A one-kernel module, as an author would write it: no ierr, no validation."""
    return module(
        "tox_demo_kernel",
        procedure(
            "scale_vector_kernel",
            real("vector", Intent.INOUT, "(n)", doc="the data to scale in place"),
            integer("n", Intent.IN, doc="length of `vector`"),
            real(
                "factor",
                Intent.IN,
                directives=Directives(minimum=Minimum("0.0_real64")),
                doc="scale factor",
            ),
            meta=Meta(summary="Scale a vector", author="AUTHOR"),  # a kernel: not exported
        ),
    )


class TestSynthesis:
    def test_a_generated_module_is_injected(self):
        result = synthesize_wrappers(project(kernel_module()))

        assert result.project.module("tox_demo") is not None
        # the kernel module stays, inert
        assert result.project.module("tox_demo_kernel") is not None

    def test_the_wrapper_drops_the_kernel_suffix(self):
        result = synthesize_wrappers(project(kernel_module()))

        wrapper = result.project.procedure("tox_demo", "scale_vector")
        assert wrapper is not None
        assert result.project.procedure("tox_demo", "scale_vector_kernel") is None

    def test_the_wrapper_is_exported(self):
        result = synthesize_wrappers(project(kernel_module()))

        wrapper = result.project.procedure("tox_demo", "scale_vector")
        assert wrapper.is_exported

    def test_the_wrapper_takes_the_kernel_arguments_plus_ierr(self):
        result = synthesize_wrappers(project(kernel_module()))

        wrapper = result.project.procedure("tox_demo", "scale_vector")
        names = [a.name for a in wrapper.arguments]
        assert names == ["vector", "n", "factor", "ierr"]
        assert wrapper.argument("ierr").intent is Intent.OUT

    def test_the_kernel_keeps_no_ierr(self):
        # a kernel that already had an ierr would not get a second one; this one had none
        result = synthesize_wrappers(project(kernel_module()))
        kernel = result.project.procedure("tox_demo_kernel", "scale_vector_kernel")

        assert kernel.argument("ierr") is None

    def test_the_spec_links_wrapper_to_kernel(self):
        result = synthesize_wrappers(project(kernel_module()))

        (spec,) = result.specs
        assert spec.kernel.name == "scale_vector_kernel"
        assert spec.validating.name == "scale_vector"
        assert spec.module_name == "tox_demo"

    def test_a_non_kernel_procedure_is_not_wrapped(self):
        # a recommend routine sitting in the kernel module has no _kernel suffix
        mod = module(
            "tox_demo_kernel",
            procedure(
                "scale_vector_kernel",
                real("vector", Intent.INOUT, "(n)"),
                integer("n", Intent.IN),
                meta=Meta(),
            ),
            procedure(
                "recommended_factor",
                integer("n", Intent.IN),
                real("factor", Intent.OUT),
                meta=Meta(),
            ),
        )
        result = synthesize_wrappers(project(mod))

        assert result.project.procedure("tox_demo", "scale_vector") is not None
        assert result.project.procedure("tox_demo", "recommended_factor") is None
        assert len(result.specs) == 1
