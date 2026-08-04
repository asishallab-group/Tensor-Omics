"""Synthesising the validating wrapper from a kernel.

A kernel carries no `ierr` and no validation; the wrapper the generator builds from it
adds both and is the exported entry point. These tests hand-build a kernel and check the
injected wrapper, without going near Ford.
"""

from codegen.ir.directives import (
    Directives,
    Minimum,
    OutputFrom,
    OutputFromMode,
    RequiredIfMode,
)
from codegen.ir.entities import Meta
from codegen.ir.types import Intent
from codegen.synthesize import synthesize_wrappers

from builders import ierr, integer, module, procedure, project, real


def mode_split_kernel_module():
    """A kernel with a mode argument whose table names a procedure per mode."""
    mode_table = [
        "Which pattern to detect.",
        "| Mode | Value | Procedure |",
        "|------|-------|-----------|",
        "| dosage effect | [[tox_demo_kernel(module):MODE_DOSAGE(variable)]] | detect_dosage_effect |",
        "| subfunctionalisation | [[tox_demo_kernel(module):MODE_SUBFUNC(variable)]] | detect_subfunctionalization |",
    ]
    return module(
        "tox_demo_kernel",
        procedure(
            "detect_patterns_kernel",
            real("genes", Intent.IN, "(n)", doc="gene values"),
            integer("n", Intent.IN, doc="number of genes"),
            integer("pattern_mode", Intent.IN, doc=mode_table),
            real(
                "threshold",
                Intent.IN,
                optional=True,
                directives=Directives(
                    required_if_mode=RequiredIfMode(
                        "pattern_mode", "tox_demo_kernel", "MODE_DOSAGE"
                    )
                ),
                doc="dosage threshold",
            ),
            real("result", Intent.OUT, "(n)", doc="pattern score per gene"),
            meta=Meta(summary="Detect paralog patterns", author="AUTHOR"),
        ),
    )

C_BINDING = Meta(summary="a summary", author="AUTHOR", category="C-binding")


def alloc_kernel_module():
    """A kernel that needs work arrays: a scratch buffer, a permutation, and a
    recommend-sized buffer -- plus the recommend routine it is sized by."""
    return module(
        "tox_demo_kernel",
        procedure(  # a recommend routine: exported, but not a kernel
            "work_size",
            integer("n", Intent.IN, doc="length"),
            integer("wsize", Intent.OUT, doc="recommended work size"),
            meta=C_BINDING,
        ),
        procedure(
            "crunch_kernel",
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length of `values`"),
            integer("values_perm", Intent.IN, "(n)", doc="ascending permutation of `values`"),
            real("tmp_scratch", Intent.OUT, "(n)", doc="scratch space"),
            real("tmp_work", Intent.OUT, "(wsize)", doc="work buffer"),
            integer(
                "wsize",
                Intent.IN,
                directives=Directives(
                    output_from=OutputFrom(
                        "wsize", "work_size", "tox_demo_kernel", OutputFromMode.AUTO
                    )
                ),
                doc="work-buffer size",
            ),
            meta=Meta(summary="Crunch the data", author="AUTHOR"),  # a kernel: not exported
        ),
    )


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


def ierr_kernel_module():
    """A kernel that itself declares `ierr` -- it propagates a sub-helper's failure."""
    return module(
        "tox_demo_kernel",
        procedure(
            "risky_kernel",
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length of `values`"),
            real("result", Intent.OUT, "(n)", doc="the result"),
            ierr(),
            meta=Meta(summary="A fallible kernel", author="AUTHOR"),
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

    def test_a_kernel_that_declares_ierr_gets_no_second_one(self):
        result = synthesize_wrappers(project(ierr_kernel_module()))

        foo = result.project.procedure("tox_demo", "risky")
        names = [a.name for a in foo.arguments]
        assert names == ["values", "n", "result", "ierr"]

    def test_the_spec_links_wrapper_to_kernel(self):
        result = synthesize_wrappers(project(kernel_module()))

        (spec,) = result.specs
        assert spec.kernel.name == "scale_vector_kernel"
        assert spec.validating.name == "scale_vector"
        assert spec.module_name == "tox_demo"

    def test_no_allocating_wrapper_without_work_arrays(self):
        # scale_vector_kernel takes no tmp_/perm/computed argument
        result = synthesize_wrappers(project(kernel_module()))

        assert result.project.procedure("tox_demo", "scale_vector_alloc") is None
        assert result.specs[0].allocating is None

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


class TestAllocatingSynthesis:
    def test_an_allocating_wrapper_is_generated(self):
        result = synthesize_wrappers(project(alloc_kernel_module()))

        assert result.project.procedure("tox_demo", "crunch_alloc") is not None
        assert result.specs[0].allocating is not None

    def test_it_drops_the_work_arrays_permutation_and_computed_size(self):
        result = synthesize_wrappers(project(alloc_kernel_module()))

        alloc = result.project.procedure("tox_demo", "crunch_alloc")
        names = [a.name for a in alloc.arguments]
        # values + n survive; values_perm, tmp_scratch, tmp_work, wsize are taken over
        assert names == ["values", "n", "ierr"]

    def test_the_validating_wrapper_keeps_everything(self):
        result = synthesize_wrappers(project(alloc_kernel_module()))

        foo = result.project.procedure("tox_demo", "crunch")
        names = [a.name for a in foo.arguments]
        assert names == [
            "values", "n", "values_perm", "tmp_scratch", "tmp_work", "wsize", "ierr"
        ]

    def test_both_variants_are_exported(self):
        result = synthesize_wrappers(project(alloc_kernel_module()))

        assert result.project.procedure("tox_demo", "crunch").is_exported
        assert result.project.procedure("tox_demo", "crunch_alloc").is_exported


class TestModeSplitSynthesis:
    def result(self):
        return synthesize_wrappers(project(mode_split_kernel_module()))

    def test_one_procedure_per_mode_named_from_the_table(self):
        project_ = self.result().project

        assert project_.procedure("tox_demo", "detect_dosage_effect") is not None
        assert project_.procedure("tox_demo", "detect_subfunctionalization") is not None
        # the mode argument is gone, so there is no single runtime-mode procedure
        assert project_.procedure("tox_demo", "detect_patterns") is None

    def test_the_mode_argument_is_dropped(self):
        dosage = self.result().project.procedure("tox_demo", "detect_dosage_effect")

        assert "pattern_mode" not in [a.name for a in dosage.arguments]

    def test_a_required_in_this_mode_argument_is_mandatory(self):
        dosage = self.result().project.procedure("tox_demo", "detect_dosage_effect")

        threshold = dosage.argument("threshold")
        assert threshold is not None
        assert not threshold.optional
        # the now-unconditional requirement drops the directive
        assert threshold.directives.required_if_mode is None

    def test_a_required_in_another_mode_argument_is_absent(self):
        subfunc = self.result().project.procedure(
            "tox_demo", "detect_subfunctionalization"
        )

        assert subfunc.argument("threshold") is None
        assert [a.name for a in subfunc.arguments] == ["genes", "n", "result", "ierr"]

    def test_the_mode_fix_is_recorded_on_the_spec(self):
        specs = self.result().specs
        dosage = next(s for s in specs if s.validating.name == "detect_dosage_effect")

        assert dosage.mode_fix.argument == "pattern_mode"
        assert dosage.mode_fix.parameter == "MODE_DOSAGE"
        assert dosage.mode_fix.module == "tox_demo_kernel"
        assert dosage.kernel.name == "detect_patterns_kernel"

    def test_an_ordinary_kernel_has_no_mode_fix(self):
        (spec,) = synthesize_wrappers(project(kernel_module())).specs

        assert spec.mode_fix is None
