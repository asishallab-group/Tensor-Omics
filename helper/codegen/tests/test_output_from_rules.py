"""The `DM_OUTPUT_FROM` rules that settle what a producer call is and what it is given.

Rule K -- a producer call is taken over whole or kept whole -- with its refusal N1 and its
assertion; the producer-input table's grammar and rule O (`OMITTED`); rule T (a producer input is
given a value of its own kind); and the optional pass-through of P1. The cases that run are in
`test_end_to_end_impl.py` (`fixtures/impl/fx_sized_impl.F90`).
"""

import pytest

from codegen.config import CONVENTIONS
from codegen.diagnostics import DiagnosticBag
from codegen.ir.directives import Directives
from codegen.ir.doc import Doc
from codegen.ir.types import Intent
from codegen.synthesize import binding_computes, synthesize_wrappers, taken_over_arguments

from builders import integer, module, parameter, procedure, project, real
from test_output_from_producers import (
    IMPL, analysed, auto, body, consumer, fortran, producers, table,
)


def norms_producer():
    """`large_norms_impl`: one call gives an array and a size that depends on its values."""
    return procedure(
        "large_norms_impl",
        integer("n_dims", doc="rows"),
        integer("n_items", doc="columns"),
        real("x", Intent.IN, "(n_dims, n_items)", doc="the matrix"),
        real("norms", Intent.OUT, "(n_items)", doc="the norm of each column"),
        integer("n_large", Intent.OUT, doc="how many exceed one"),
        meta=IMPL, is_pure=True,
    )


def norms_consumer(name, last, *, norms_doc="the norms", n_large_doc="how many", leading=(),
                   parameters=()):
    """A consumer of both outputs of one `large_norms` call, `last` sized by `n_large`."""
    return module("tox_cons_impl", procedure(
        name,
        *leading,
        integer("n_dims", doc="rows"),
        integer("n_items", doc="columns"),
        real("x", Intent.IN, "(n_dims, n_items)", doc="the matrix"),
        real("norms", Intent.IN, "(n_items)", doc=norms_doc,
             directives=auto("norms", "large_norms", "tox_prod")),
        integer("n_large", doc=n_large_doc, directives=auto("n_large", "large_norms", "tox_prod")),
        last,
        real("total", Intent.OUT, doc="the result"),
        meta=IMPL, is_pure=True,
    ), parameters=parameters)


def errors(bag):
    return [d for d in bag.errors]


def located(argument, lines, first_line=100):
    """`argument` documented by `lines` numbered from `first_line`, as a parsed source's are:
    a `table(...)`'s header is then `first_line + 1` and its first row `first_line + 3`."""
    argument.doc = Doc.parse(lines, first_line_number=first_line)
    return argument


class TestACallIsKeptWhole:
    """Rule K: where `foo` must keep one value a call fills, it keeps every value that call fills."""

    def test_data_and_a_size_that_sizes_a_result_are_both_kept(self):
        synthesis, bag = analysed(producers(norms_producer()), norms_consumer(
            "picked_impl", real("picked", Intent.OUT, "(n_large)", doc="the picked")))

        assert errors(bag) == []
        names = [a.name for a in synthesis.project.procedure("tox_cons", "picked").arguments]
        assert {"norms", "n_large"} <= set(names)
        # nothing is left to take over: one wrapper, which calls no producer in Fortran
        assert synthesis.project.procedure("tox_cons", "picked_expert") is None
        assert "call large_norms(" not in fortran(synthesis)

    def test_the_bindings_compute_both_in_one_call(self):
        synthesis, _ = analysed(producers(norms_producer()), norms_consumer(
            "picked_impl", real("picked", Intent.OUT, "(n_large)", doc="the picked")))
        wrapper = synthesis.project.procedure("tox_cons", "picked")
        norms, n_large = wrapper.argument("norms"), wrapper.argument("n_large")

        assert binding_computes(norms, wrapper) and binding_computes(n_large, wrapper)
        assert norms.roles.computed_from.call == n_large.roles.computed_from.call

    def test_a_call_nothing_keeps_is_taken_over_whole(self):
        synthesis, bag = analysed(producers(norms_producer()), norms_consumer(
            "large_total_impl", real("tmp_large", Intent.OUT, "(n_large)", doc="scratch")))

        assert errors(bag) == []
        plain = synthesis.project.procedure("tox_cons", "large_total")
        assert plain.argument("norms") is None and plain.argument("n_large") is None
        assert body(fortran(synthesis), "large_total").count("call large_norms(") == 1

    def test_the_collapse_a_size_of_a_kept_input(self):
        # `n_large` sizes the weights the caller passes: K keeps the whole call, nothing is left
        # to take over, and the pair becomes one wrapper whose binding computes both (Q10)
        synthesis, bag = analysed(producers(norms_producer()), norms_consumer(
            "weighted_impl", real("w", Intent.IN, "(n_large)", doc="weights")))

        assert errors(bag) == []
        assert synthesis.project.procedure("tox_cons", "weighted_expert") is None
        wrapper = synthesis.project.procedure("tox_cons", "weighted")
        assert binding_computes(wrapper.argument("norms"), wrapper)
        assert binding_computes(wrapper.argument("n_large"), wrapper)

    def test_one_call_spelt_two_ways_is_one_call(self):
        # K2: `4_int32` and the parameter FOUR are one constant, so one call, kept whole
        cons = module(
            "tox_cons_impl",
            procedure(
                "spell_impl",
                integer("a_size", doc=table(("n_items", "4_int32")),
                        directives=auto("first", "two_sizes", "tox_prod")),
                real("y", Intent.OUT, "(a_size)", doc="the result"),
                integer("b_size", doc=table(("n_items", "FOUR")),
                        directives=auto("second", "two_sizes", "tox_prod")),
                real("tmp_b", Intent.OUT, "(b_size)", doc="scratch"),
                meta=IMPL, is_pure=True,
            ),
            parameters=(parameter("FOUR", "4_int32"),),
        )
        synthesis, bag = analysed(producers(), cons)

        assert errors(bag) == []
        assert synthesis.project.procedure("tox_cons", "spell").argument("b_size") is not None

    def test_two_tables_are_two_calls(self):
        # K6: `n_a` and `n_b` are two calls: the one sizing the result is kept, the other taken
        synthesis, bag = analysed(producers(), consumer(
            "tables_impl",
            integer("n_a", doc="items of the result"),
            integer("n_b", doc="items of the scratch"),
            integer("a_size", doc=table(("n_items", "n_a")),
                    directives=auto("first", "two_sizes", "tox_prod")),
            real("y", Intent.OUT, "(a_size)", doc="the result"),
            integer("b_size", doc=table(("n_items", "n_b")),
                    directives=auto("second", "two_sizes", "tox_prod")),
            real("tmp_b", Intent.OUT, "(b_size)", doc="scratch"),
        ))

        assert errors(bag) == []
        plain = synthesis.project.procedure("tox_cons", "tables")
        assert plain.argument("a_size") is not None and plain.argument("b_size") is None

    def test_an_identity_row_is_no_row(self):
        # `| n_items | N_ITEMS |` names the argument name-matching finds: one call with the other
        synthesis, bag = analysed(producers(), consumer(
            "same_impl",
            integer("n_items", doc="items"),
            integer("a_size", doc=table(("n_items", "N_ITEMS")),
                    directives=auto("first", "two_sizes", "tox_prod")),
            real("y", Intent.OUT, "(a_size)", doc="the result"),
            integer("b_size", doc="the scratch size",
                    directives=auto("second", "two_sizes", "tox_prod")),
            real("tmp_b", Intent.OUT, "(b_size)", doc="scratch"),
        ))

        assert errors(bag) == []
        assert synthesis.project.procedure("tox_cons", "same").argument("b_size") is not None

    def test_one_evaluator_keys_the_calls_everywhere(self):
        source = project(producers(), consumer("x_impl", integer("n", doc="n")))
        synthesis = synthesize_wrappers(source, CONVENTIONS, DiagnosticBag())

        assert synthesis.project.evaluator() is source.evaluator()


class TestTheAssertion:
    """After analysis, no `plan.call` group of a pair is kept in part: an internal error."""

    def test_a_call_keyed_apart_from_its_plans_is_reported(self, monkeypatch):
        # break the key so each value is its own call: `n_large` is kept for the result it
        # sizes, `norms` taken over, while their plans are one call
        import codegen.synthesize as synthesize

        monkeypatch.setattr(synthesize, "producer_call_key",
                            lambda argument, *rest: argument.name.lower())
        _, bag = analysed(producers(norms_producer()), norms_consumer(
            "picked_impl", real("picked", Intent.OUT, "(n_large)", doc="the picked")))

        [error] = [d for d in bag.errors if "internal error" in d.message]
        assert "'picked' keeps 'n_large' but takes over 'norms'" in error.message
        assert "'large_norms'" in error.message
        assert "rule K" in error.note


class TestARowNamesAnInput:
    """N1: a producer-input table row names an input of the producer, not an output or `ierr`."""

    def outrow(self, row=("n_large", "1_int32"), **kwargs):
        cons = norms_consumer(
            "outrow_impl", real("picked", Intent.OUT, "(n_large)", doc="the picked"),
            norms_doc=table(row), **kwargs)
        located(cons.procedures[0].argument("norms"), table(row))
        return cons

    def test_a_row_naming_an_output_is_refused_at_its_line(self):
        _, bag = analysed(producers(norms_producer()), self.outrow())

        [error] = bag.errors
        assert error.message == ("the producer-input table of 'norms' maps 'n_large', which is "
                                 "an output of 'large_norms', not an input")
        assert error.note == ("a row supplies an input the producer reads; DM_OUTPUT_FROM names "
                              "the output it takes")
        assert error.location.line == 103  # the row, not the directive

    def test_a_row_naming_ierr_is_refused(self):
        _, bag = analysed(producers(norms_producer()), self.outrow(("ierr", "1_int32")))

        assert [e.message for e in bag.errors] == [
            "the producer-input table of 'norms' maps 'ierr', which is an output of "
            "'large_norms', not an input"]

    def test_a_mode_split_reports_it_once(self):
        # each per-mode wrapper finds the row again; the cause is one row of the implementation
        modes = [
            "Which way.",
            "| Mode | Value | Procedure |",
            "|------|-------|-----------|",
            "| one | [[tox_cons_impl(module):MODE_ONE(variable)]] | outrow_one |",
            "| two | [[tox_cons_impl(module):MODE_TWO(variable)]] | outrow_two |",
        ]
        cons = self.outrow(leading=(integer("mode", doc=modes),), parameters=(
            parameter("MODE_ONE", "1_int32"), parameter("MODE_TWO", "2_int32")))
        synthesis, bag = analysed(producers(norms_producer()), cons)

        # the premise: two wrappers, each resolving the table
        assert synthesis.project.procedure("tox_cons", "outrow_one") is not None
        assert synthesis.project.procedure("tox_cons", "outrow_two") is not None

        assert [e.message for e in bag.errors] == [
            "the producer-input table of 'norms' maps 'n_large', which is an output of "
            "'large_norms', not an input"]
