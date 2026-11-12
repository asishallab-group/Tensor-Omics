"""Which argument an error blames when the Fortran one is derived (`emit.error_arguments`).

An error names a position in the Fortran dummy list; the binding names the argument the
caller passed in its place. That is only true of an argument the caller did pass.
"""

from codegen.config import CONVENTIONS
from codegen.diagnostics import DiagnosticBag
from codegen.emit.error_arguments import source_of, sources_of
from codegen.ir.roles import analyse
from codegen.ir.types import Intent

import builders as b


def analysed(*arguments):
    procedure = b.procedure("demo", *arguments, b.ierr())
    analyse(procedure, DiagnosticBag(), CONVENTIONS)
    return procedure


class TestAnExtent:
    def test_it_names_the_input_it_is_read_off(self):
        procedure = analysed(
            b.integer("n"),
            b.real("values", Intent.IN, "(n)", doc="the data"),
        )
        assert source_of(procedure.argument("n")) == "values"

    def test_it_names_an_input_even_when_an_output_is_declared_first(self):
        # `n_tissues` sized `log_transformed_expr` before `reps_per_tissue` in
        # normalization_pipeline, and an empty `reps_per_tissue` blamed the output
        procedure = analysed(
            b.integer("n"),
            b.real("result", Intent.OUT, "(n)", doc="the result"),
            b.integer("counts", Intent.IN, "(n)", doc="the counts"),
        )
        extent = procedure.argument("n")
        assert extent.roles.extent_of[0].name == "result"  # the premise
        assert source_of(extent) == "counts"

    def test_one_sizing_outputs_alone_is_the_callers_own(self):
        # the caller passes it, so the error names it, not an output they never wrote
        procedure = analysed(
            b.integer("n_vectors"),
            b.real("status", Intent.OUT, "(n_vectors)", doc="the status"),
        )
        extent = procedure.argument("n_vectors")
        assert extent.roles.is_extent and not extent.roles.is_derived  # the premise
        assert source_of(extent) is None
        assert sources_of(procedure) == [None, None, None]


class TestAMaskCount:
    def test_one_sizing_outputs_alone_names_its_mask(self):
        # `n_selected_vectors` sizes only outputs, but the caller passed the mask it counts
        procedure = analysed(
            b.integer("n_vectors"),
            b.logical("vectors_selection_mask", Intent.IN, "(n_vectors)", kind="c_bool",
                      doc="the mask"),
            b.integer("n_selected_vectors"),
            b.real("result", Intent.OUT, "(n_selected_vectors)", doc="the result"),
        )
        count = procedure.argument("n_selected_vectors")
        assert count.roles.is_extent and not count.roles.readable_owners  # the premise
        assert source_of(count) == "vectors_selection_mask"

    def test_one_sizing_nothing_names_its_mask(self):
        procedure = analysed(
            b.integer("n_vectors"),
            b.logical("vectors_selection_mask", Intent.IN, "(n_vectors)", kind="c_bool",
                      doc="the mask"),
            b.integer("n_selected_vectors"),
        )
        assert source_of(procedure.argument("n_selected_vectors")) == "vectors_selection_mask"
