"""The `DM_OUTPUT_FROM` rules that settle what a producer call is and what it is given.

Rule K -- a producer call is taken over whole or kept whole -- with its refusal N1 and its
assertion; the producer-input table's grammar and rule O (`OMITTED`); rule T (a producer input is
given a value of its own kind); and the optional pass-through of P1. The cases that run are in
`test_end_to_end_impl.py` (`fixtures/impl/fx_sized_impl.F90`).
"""

import pytest

from codegen.config import CONVENTIONS
from codegen.diagnostics import DiagnosticBag
from codegen.ir.directives import (
    AllowInfinite, AllowNan, Default, Directives, Equals, Maximum, Minimum, OutputFrom,
    OutputFromMode, Prologue, RequiredIfMode, Sentinel,
)
from codegen.ir.validate import validate_project
from codegen.ir.doc import Doc
from codegen.ir.types import Intent
from codegen.synthesize import binding_computes, synthesize_wrappers, taken_over_arguments

from builders import character, integer, logical, module, parameter, procedure, project, real
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


def optional_producers():
    """`tox_prod_impl` with producers whose inputs are optional: one without a default, one with,
    and one required in a mode chosen at run time."""
    modes = [
        "Which way.",
        "| Mode | Value |",
        "|------|-------|",
        "| plain | [[tox_prod_impl(module):MODE_PLAIN(variable)]] |",
        "| gained | [[tox_prod_impl(module):MODE_GAINED(variable)]] |",
    ]
    return module(
        "tox_prod_impl",
        procedure(
            "padded_size_impl",
            integer("n_items", doc="items"),
            integer("padding", doc="extra slots", optional=True),
            integer("padded", Intent.OUT, doc="the size"),
            meta=IMPL, is_pure=True,
        ),
        procedure(
            "strided_size_impl",
            integer("n_items", doc="items"),
            integer("stride", doc="slots per item", optional=True,
                    directives=Directives(default=Default("3_int32"))),
            integer("strided", Intent.OUT, doc="the size"),
            meta=IMPL, is_pure=True,
        ),
        procedure(
            "gained_size_impl",
            integer("n_items", doc="items"),
            integer("mode", doc=modes),
            integer("gain", doc="the gain", optional=True, directives=Directives(
                required_if_mode=RequiredIfMode("mode", "tox_prod_impl", "MODE_GAINED"))),
            integer("gained", Intent.OUT, doc="the size"),
            meta=IMPL, is_pure=True,
        ),
        parameters=(parameter("MODE_PLAIN", "1_int32"), parameter("MODE_GAINED", "2_int32")),
    )


def omitting(rows, producer="padded_size", output="padded", *extra, name="omit_impl",
             sized="tmp_buf", leading=(), parameters=(), mode_dm="AUTO"):
    """A consumer whose `work_size` is computed from `producer` with the table `rows`, located."""
    mode = OutputFromMode.AUTO if mode_dm == "AUTO" else OutputFromMode.JUST_INFO
    work = integer("work_size", doc="the size", directives=Directives(
        output_from=OutputFrom(output, producer, "tox_prod", mode)))
    located(work, table(*rows))
    last = (real(sized, Intent.OUT, "(work_size)", doc="scratch") if sized else
            real("total", Intent.OUT, doc="the total"))
    return module("tox_cons_impl", procedure(
        name, *leading, integer("n_items", doc="items"), *extra, work, last,
        meta=IMPL, is_pure=True,
    ), parameters=parameters)


class TestOmitted:
    """Rule O: a cell `OMITTED` leaves the producer input out of the call, in every tier."""

    def test_no_call_passes_it(self):
        synthesis, bag = analysed(optional_producers(), omitting([("padding", "OMITTED")]))

        assert errors(bag) == []
        plan = synthesis.project.procedure("tox_cons", "omit_expert").argument(
            "work_size").roles.computed_from
        assert [supply.name for supply in plan.inputs] == ["n_items"]
        call = body(fortran(synthesis), "omit").split("call padded_size(", 1)[1].split(")", 1)[0]
        assert "padding" not in call and "n_items = n_items" in call

    def test_a_row_overrides_name_matching(self):
        # the consumer has a `padding` of its own; the row still leaves the producer's out
        synthesis, bag = analysed(optional_producers(), omitting(
            [("padding", "OMITTED")], "padded_size", "padded", integer("padding", doc="its own")))

        assert errors(bag) == []
        plan = synthesis.project.procedure("tox_cons", "omit_expert").argument(
            "work_size").roles.computed_from
        assert [supply.name for supply in plan.inputs] == ["n_items"]

    def test_a_defaulted_optional_may_be_omitted(self):
        _, bag = analysed(optional_producers(), omitting(
            [("stride", "OMITTED")], "strided_size", "strided"))

        assert errors(bag) == []

    def test_omitting_and_supplying_are_two_calls(self):
        cons = module("tox_cons_impl", procedure(
            "twice_impl",
            integer("n_items", doc="items"),
            located(integer("a_size", doc="a", directives=auto("strided", "strided_size",
                                                                 "tox_prod")),
                    table(("stride", "OMITTED"))),
            real("tmp_a", Intent.OUT, "(a_size)", doc="scratch"),
            located(integer("b_size", doc="b", directives=auto("strided", "strided_size",
                                                                 "tox_prod")),
                    table(("stride", "2_int32"))),
            real("tmp_b", Intent.OUT, "(b_size)", doc="scratch"),
            meta=IMPL, is_pure=True,
        ))
        synthesis, bag = analysed(optional_producers(), cons)

        assert errors(bag) == []
        assert body(fortran(synthesis), "twice").count("call strided_size(") == 2

    def test_a_required_input_is_refused_at_its_row(self):
        _, bag = analysed(optional_producers(), omitting([("n_items", "OMITTED")]))

        [error] = bag.errors
        assert error.message == ("the producer-input table of 'work_size' omits 'n_items', "
                                 "which 'padded_size' requires")
        assert error.note == ("only an optional input may be left out; supply 'n_items' by an "
                              "argument the caller passes, or a constant")
        assert error.location.line == 103

    def test_an_input_a_runtime_mode_requires_is_refused(self):
        _, bag = analysed(optional_producers(), omitting(
            [("gain", "OMITTED")], "gained_size", "gained", integer("mode", doc=[
                "Which way.", "| Mode | Value |", "|------|-------|",
                "| plain | [[tox_cons_impl(module):MODE_PLAIN(variable)]] |",
                "| gained | [[tox_cons_impl(module):MODE_GAINED(variable)]] |"]),
            parameters=(parameter("MODE_PLAIN", "1_int32"), parameter("MODE_GAINED", "2_int32"))))

        assert [(e.message, e.note) for e in bag.errors] == [(
            "the producer-input table of 'work_size' omits 'gain', which 'gained_size' requires "
            "in mode 'gained'",
            "the mode is chosen at run time, so the input may be needed: supply 'gain'")]

    @pytest.mark.parametrize("cell", ["Omitted", "omitted", "'OMITTED'", '"OMITTED"', "`OMITTED`"])
    def test_any_other_spelling_is_refused(self, cell):
        _, bag = analysed(optional_producers(), omitting([("padding", cell)]))

        [error] = bag.errors
        assert error.message == (f"the producer-input table of 'work_size' supplies 'padding' by "
                                 f"'{cell}', which is not the keyword OMITTED")
        assert error.note.startswith("the keyword is written bare and in capitals, as AUTO is")
        assert error.location.line == 103

    def test_an_argument_of_that_name_is_refused(self):
        _, bag = analysed(optional_producers(), omitting(
            [("padding", "OMITTED")], "padded_size", "padded", integer("omitted", doc="a count")))

        [error] = bag.errors
        assert error.message == (
            "the producer-input table of 'work_size' supplies 'padding' by 'OMITTED', which is the "
            "keyword, but 'omitted' is also an argument of 'omit_impl'")
        assert error.note.startswith("rename")

    def test_a_constant_of_that_name_is_refused(self):
        _, bag = analysed(optional_producers(), omitting(
            [("padding", "OMITTED")], parameters=(parameter("OMITTED", "0_int32"),)))

        [error] = bag.errors
        assert error.message == (
            "the producer-input table of 'work_size' supplies 'padding' by 'OMITTED', which is the "
            "keyword, but 'omitted' is also a constant")

    def test_an_optional_nothing_supplies_is_still_an_error_whose_note_offers_the_row(self):
        cons = omitting([], name="omit_impl")
        synthesis, bag = analysed(optional_producers(), cons)

        [error] = bag.errors
        assert error.message == ("'work_size' is computed from 'padded_size', but its input(s) "
                                 "padding have no argument of the same name in 'omit_impl'")
        assert error.note.startswith(
            "'padding' is optional, so the row `| padding | OMITTED |` leaves it out of the call")
        assert "| padding        | OMITTED     |" in error.note


class TestTheTableGrammar:
    """One row per input, an input in each, the same rules on a JUST_INFO table."""

    def test_two_rows_for_one_input_are_refused_at_the_second(self):
        _, bag = analysed(optional_producers(), omitting(
            [("padding", "n_items"), ("padding", "1_int32")]))

        [error] = bag.errors
        assert error.message == ("the producer-input table of 'work_size' has two rows for "
                                 "'padding' (rows 1 and 2)")
        assert error.location.line == 104
        assert error.note

    def test_a_just_info_table_is_checked_too(self):
        _, bag = analysed(optional_producers(), omitting(
            [("n_items", "OMITTED")], sized="", mode_dm="JUST_INFO"))

        assert [e.message for e in bag.errors] == [
            "the producer-input table of 'work_size' omits 'n_items', which 'padded_size' requires"]

    def test_one_error_per_cell(self):
        _, bag = analysed(optional_producers(), omitting(
            [("padding", "Omitted"), ("padded", "1_int32")]))

        assert len(bag.errors) == 2


class TestWhatNothingSupplies:
    """The unmatched input at the directive's line; a cell naming nothing at its row (4a); a
    split mode as a supplier (4b, T); each naming the implementation, not the wrapper (4c)."""

    def split(self, row=None):
        modes = [
            "Which way.",
            "| Mode | Value | Procedure |",
            "|------|-------|-----------|",
            "| one | [[tox_cons_impl(module):MODE_ONE(variable)]] | scaled_one |",
            "| two | [[tox_cons_impl(module):MODE_TWO(variable)]] | scaled_two |",
        ]
        work = integer("work_size", doc="the size", directives=Directives(
            output_from=OutputFrom("scaled", "scaled_size", "tox_prod", OutputFromMode.AUTO, 7)))
        located(work, table(row or ("factor", "mode")))
        return module("tox_cons_impl", procedure(
            "scaled_impl", integer("mode", doc=modes), integer("n_items", doc="items"), work,
            real("tmp_buf", Intent.OUT, "(work_size)", doc="scratch"),
            meta=IMPL, is_pure=True,
        ), parameters=(parameter("MODE_ONE", "1_int32"), parameter("MODE_TWO", "2_int32")))

    def test_a_split_mode_is_refused_once_per_implementation(self):
        _, bag = analysed(producers(), self.split())

        [error] = bag.errors
        assert error.message == ("'work_size' is computed from 'scaled_size', whose input 'factor' "
                                 "would be 'mode', the mode 'scaled_impl' splits into one wrapper "
                                 "per mode")
        assert error.location.line == 7  # the directive's
        assert "each per-mode wrapper fixes the mode" in error.note

    def test_a_cell_naming_nothing_is_at_its_row(self):
        # 4a: `2*scale` is neither an argument nor a constant; the row, not the directive
        _, bag = analysed(producers(), self.split(("factor", "2*scale")))

        [error] = bag.errors
        assert error.message == ("the producer-input table of 'work_size' supplies 'factor' by "
                                 "'2*scale', which is neither an argument of 'scaled_impl' nor a "
                                 "constant")
        assert error.location.line == 103
        assert "an expression is not computed" in error.note

    def test_an_unmatched_input_is_at_the_directive_once_naming_the_implementation(self):
        _, bag = analysed(producers(), self.split(("n_items", "n_items")))

        [error] = bag.errors
        assert error.message == ("'work_size' is computed from 'scaled_size', but its input(s) "
                                 "factor have no argument of the same name in 'scaled_impl'")
        assert error.location.line == 7


def method_doc(module_name, *names):
    rows = [f"| {n.lower()} | [[{module_name}(module):METHOD_{n}(variable)]] |" for n in names]
    return ["How.", "| Method | Value |", "|--------|-------|", *rows]


def kinded_producers(two_is=2):
    """`tox_prod_impl`: producers whose inputs are of several kinds, a mode among them."""
    return module(
        "tox_prod_impl",
        procedure("method_size_impl",
                  integer("n_items", doc="items"),
                  integer("method", doc=method_doc("tox_prod_impl", "ONE", "TWO", "THREE")),
                  integer("sized", Intent.OUT, doc="the size"), meta=IMPL, is_pure=True),
        procedure("factor_size_impl",
                  integer("n_items", doc="items"),
                  integer("factor", doc="a factor"),
                  integer("sized", Intent.OUT, doc="the size"), meta=IMPL, is_pure=True),
        procedure("ratio_size_impl",
                  integer("n_items", doc="items"),
                  real("ratio", doc="a ratio"),
                  integer("sized", Intent.OUT, doc="the size"), meta=IMPL, is_pure=True),
        procedure("flag_size_impl",
                  integer("n_items", doc="items"),
                  logical("flag", kind="c_bool", doc="a flag"),
                  integer("sized", Intent.OUT, doc="the size"), meta=IMPL, is_pure=True),
        procedure("label_size_impl",
                  integer("n_items", doc="items"),
                  character("label", length="8", doc="a label"),
                  integer("sized", Intent.OUT, doc="the size"), meta=IMPL, is_pure=True),
        procedure("extent_size_impl",
                  integer("n", doc="values"),
                  real("x", Intent.IN, "(n)", doc="the values"),
                  integer("sized", Intent.OUT, doc="the size"), meta=IMPL, is_pure=True),
        parameters=(parameter("METHOD_ONE", "1_int32"), parameter("METHOD_TWO", f"{two_is}_int32"),
                    parameter("METHOD_THREE", "3_int32")),
    )


def kind_consumer(producer, *own, rows=(), fortran=True, parameters=()):
    """A consumer of `producer`'s `sized`: through `foo` (a work array) or the bindings (a result)."""
    work = integer("work_size", doc="the size", directives=Directives(output_from=OutputFrom(
        "sized", producer, "tox_prod", OutputFromMode.AUTO, 9)))
    if rows:
        work.doc = Doc.parse(table(*rows))
    sized = (real("tmp_buf", Intent.OUT, "(work_size)", doc="scratch") if fortran else
             real("filled", Intent.OUT, "(work_size)", doc="the result"))
    return module("tox_cons_impl", procedure(
        "use_impl", integer("n_items", doc="items"), *own, work, sized, meta=IMPL, is_pure=True,
    ), parameters=parameters)


def refusals(*modules):
    _, bag = analysed(*modules)
    return [(e.message, e.note, e.location.line) for e in bag.errors]


class TestAModeMeetsOnlyAMode:
    """Rule T 1: a mode and a non-mode never meet; tables agree; a constant names a parameter."""

    def test_a_mode_into_a_non_mode_input(self):
        [(message, note, line)] = refusals(kinded_producers(), kind_consumer(
            "factor_size", integer("method", doc=method_doc("tox_cons_impl", "ONE")),
            rows=[("factor", "method")],
            parameters=(parameter("METHOD_ONE", "1_int32"),)))

        assert message == ("'work_size' is computed from 'factor_size', whose input 'factor' would "
                           "be the mode 'method', which Python and R pass as its string and Fortran "
                           "as its code")
        assert note.startswith("a mode is a choice, not a number (Section 5.10)")
        assert line == 9

    def test_a_non_mode_into_a_mode_input(self):
        [(message, note, _)] = refusals(kinded_producers(), kind_consumer(
            "method_size", integer("k", doc="a number"), rows=[("method", "k")]))

        assert message == ("'work_size' is computed from 'method_size', whose input 'method' is a "
                           "mode, but would be 'k', which is not")
        assert "(METHOD_ONE, METHOD_TWO, METHOD_THREE)" in note

    def test_a_constant_value_into_a_mode_input(self):
        # Q8: by parameter, not by value -- `2_int32` is METHOD_TWO's code, and still refused
        [(message, _, _)] = refusals(kinded_producers(), kind_consumer(
            "method_size", rows=[("method", "2_int32")]))

        assert message == ("'work_size' is computed from 'method_size', whose input 'method' is a "
                           "mode, but would be '2_int32', which is not")

    @pytest.mark.parametrize("fortran", [True, False])
    def test_a_mode_parameter_into_a_mode_input(self, fortran):
        synthesis, bag = analysed(kinded_producers(), kind_consumer(
            "method_size", rows=[("method", "METHOD_TWO")], fortran=fortran))

        assert errors(bag) == []
        tier = "use_expert" if fortran else "use"
        [supply] = [s for s in synthesis.project.procedure("tox_cons", tier).argument(
            "work_size").roles.computed_from.inputs if s.name == "method"]
        assert (supply.constant, supply.mode.string) == (2, "two")

    @pytest.mark.parametrize("fortran", [True, False])
    def test_the_producer_may_accept_more_modes(self, fortran):
        assert refusals(kinded_producers(), kind_consumer(
            "method_size", integer("method", doc=method_doc("tox_cons_impl", "ONE", "TWO")),
            fortran=fortran, parameters=(parameter("METHOD_ONE", "1_int32"),
                                         parameter("METHOD_TWO", "2_int32")))) == []

    def test_a_mode_the_producer_does_not_accept(self):
        [(message, note, _)] = refusals(kinded_producers(), kind_consumer(
            "method_size", integer("method", doc=method_doc("tox_cons_impl", "ONE", "FOUR")),
            parameters=(parameter("METHOD_ONE", "1_int32"), parameter("METHOD_FOUR", "4_int32"))))

        assert message == ("'work_size' is computed from 'method_size', whose mode input 'method' "
                           "accepts 'one', 'two', 'three', but would be the mode 'method', which "
                           "also accepts 'four'")
        assert note.startswith("a mode passes through as its string in Python and R")

    @pytest.mark.parametrize("fortran", [True, False])
    def test_a_mode_with_another_value_is_refused_by_its_own_module(self, fortran):
        # q_collide: METHOD_TWO is 5 in the producer's module, 2 in the consumer's. The project's
        # flat table of constants has one METHOD_TWO, so only evaluating each in its own module
        # sees that the two disagree
        [(message, _, _)] = refusals(kinded_producers(two_is=5), kind_consumer(
            "method_size", integer("method", doc=method_doc("tox_cons_impl", "ONE", "TWO")),
            fortran=fortran, parameters=(parameter("METHOD_ONE", "1_int32"),
                                         parameter("METHOD_TWO", "2_int32"))))

        assert message.endswith("but would be the mode 'method', which gives 'two' another value")

    def test_a_mode_is_not_computed(self):
        cons = module("tox_cons_impl", procedure(
            "use_impl", integer("n_items", doc="items"),
            integer("method", doc=method_doc("tox_cons_impl", "ONE"), directives=Directives(
                output_from=OutputFrom("sized", "factor_size", "tox_prod", OutputFromMode.AUTO,
                                       9))),
            integer("factor", doc="a factor"),
            real("total", Intent.OUT, doc="the total"), meta=IMPL, is_pure=True,
        ), parameters=(parameter("METHOD_ONE", "1_int32"),))
        [(message, note, line)] = refusals(kinded_producers(), cons)

        assert message == "'method' is a mode, so it cannot be computed with 'factor_size'"
        assert note.startswith("a mode is a choice the caller makes")
        assert line == 9

    def test_a_mode_table_on_an_output_is_refused(self):
        from codegen.ir.validate import validate_project

        cons = module("tox_cons_impl", procedure(
            "use_impl", integer("n_items", doc="items"),
            integer("method", Intent.OUT, doc=method_doc("tox_cons_impl", "ONE")),
            meta=IMPL, is_pure=True,
        ), parameters=(parameter("METHOD_ONE", "1_int32"),))
        synthesis, bag = analysed(cons)
        validate_project(synthesis.project, bag)

        assert [(e.message, e.note) for e in bag.errors
                if "mode table" in e.message] == [(
            "'method' is an output, but has a mode table, whose values the wrapper checks before "
            "the implementation sets it", "drop the table, and document the codes in prose")]


class TestAnArgumentOfItsKind:
    """Rule T 2: base type and rank on every route; kind, length and extents where Fortran calls."""

    def test_rank_and_type(self):
        [(message, note, _)] = refusals(kinded_producers(), kind_consumer(
            "factor_size", real("x", Intent.IN, "(n_items)", doc="values"),
            rows=[("factor", "x")], fortran=False))

        assert message == ("'work_size' is computed from 'factor_size', whose input 'factor' is "
                           "integer(int32), a scalar, but would be 'x', real(real64), rank 1")
        assert note == ("supply it by an argument of its type, kind and rank, or a constant in "
                        "the table")

    def test_a_logical_kind_on_the_fortran_route_only(self):
        flag = logical("flag", doc="a default logical")
        [(message, note, _)] = refusals(kinded_producers(), kind_consumer("flag_size", flag))

        assert message == ("'work_size' is computed from 'flag_size', whose input 'flag' is "
                           "logical(c_bool), a scalar, but would be 'flag', logical, a scalar")
        assert note.endswith("so the caller passes it")  # the JUST_INFO way out
        assert refusals(kinded_producers(), kind_consumer(
            "flag_size", logical("flag", doc="a default logical"), fortran=False)) == []

    def test_a_length(self):
        [(message, note, _)] = refusals(kinded_producers(), kind_consumer(
            "label_size", character("label", length="4", doc="a label")))

        assert message == ("'work_size' is computed from 'label_size', whose input 'label' is "
                           "character(len=8), but would be 'label', character(len=4)")
        assert note.startswith("a fixed-length input takes a value of that length")
        assert refusals(kinded_producers(), kind_consumer(
            "label_size", character("label", length="8", doc="a label"))) == []
        # a binding reads the length off the value
        assert refusals(kinded_producers(), kind_consumer(
            "label_size", character("label", length="4", doc="a label"), fortran=False)) == []

    def test_extents_on_the_fortran_route_only(self):
        # q_extent: `n = k` beside `x = y(m)` reads past `y` for k > m
        own = (integer("m", doc="values"), integer("k", doc="a count"),
               real("y", Intent.IN, "(m)", doc="the values"))
        [(message, note, _)] = refusals(kinded_producers(), kind_consumer(
            "extent_size", *own, rows=[("n", "k"), ("x", "y")]))

        assert message == ("'work_size' is computed from 'extent_size', whose input 'x' has shape "
                           "(n), with 'n' given 'k', but would be 'y', of shape (m)")
        assert note.startswith("the generated wrapper passes 'n = k' beside 'x = y', so a "
                               "different 'm' reads past 'y': supply 'n' by what sizes 'y'")
        assert refusals(kinded_producers(), kind_consumer(
            "extent_size", *own, rows=[("n", "k"), ("x", "y")], fortran=False)) == []
        assert refusals(kinded_producers(), kind_consumer(
            "extent_size", *own, rows=[("n", "m"), ("x", "y")])) == []

    def test_one_error_per_failing_input(self):
        two = module("tox_prod_impl", procedure(
            "two_size_impl", integer("a", doc="a"), integer("b", doc="b"),
            integer("sized", Intent.OUT, doc="the size"), meta=IMPL, is_pure=True))
        cons = kind_consumer("two_size", real("p", doc="p"), real("q", doc="q"),
                             rows=[("a", "p"), ("b", "q")])

        assert len(refusals(two, cons)) == 2


class TestAConstantOfItsType:
    """Rule T 3: a constant supplies a scalar input of its own type, and fits its kind."""

    @pytest.mark.parametrize("producer, row, why", [
        ("ratio_size", ("ratio", "3"), "but the table gives it the constant 3, an integer"),
        ("factor_size", ("factor", ".true."), "but the table gives it the constant .true., a logical"),
        ("factor_size", ("factor", "3000000000_int64"),
         "but the table gives it the constant 3000000000_int64, which does not fit"),
        ("label_size", ("label", "'abcdefghij'"),
         "but the table gives it the constant 'abcdefghij', which is longer than 8"),
        ("extent_size", ("x", "1.0_real64"),
         "rank 1, but the table gives it the constant 1.0_real64, a scalar"),
    ])
    @pytest.mark.parametrize("fortran", [True, False])
    def test_refused_on_either_route(self, producer, row, why, fortran):
        rows = [row] + ([("n", "n_items")] if producer == "extent_size" else [])
        [(message, note, _)] = refusals(kinded_producers(), kind_consumer(
            producer, rows=rows, fortran=fortran))

        assert message.endswith(why)
        assert note == ("a constant supplies a scalar input of its own type, written as Fortran "
                        "writes it (3.0_real64, .false., 'abc')")

    def test_a_string_is_quoted_in_fortran(self):
        synthesis, bag = analysed(kinded_producers(), kind_consumer(
            "label_size", rows=[("label", "'abc'")]))

        assert errors(bag) == []
        assert "label = 'abc     '" in fortran(synthesis)  # padded to its length of 8

    def test_a_real_and_a_logical_of_their_own_type(self):
        assert refusals(kinded_producers(), kind_consumer(
            "ratio_size", rows=[("ratio", "3.0_real64")])) == []
        assert refusals(kinded_producers(), kind_consumer(
            "flag_size", rows=[("flag", ".false.")])) == []


def validated(*modules):
    synthesis, bag = analysed(*modules)
    validate_project(synthesis.project, bag)
    return synthesis, bag


def staged(*arguments, prologue=None, guard_arguments=()):
    """`staged_impl`: a consumer whose `limit` `fx`-style `scaled_size` computes, with
    `arguments` after it, and optionally a prologue `guard` taking `guard_arguments`."""
    directives = Directives(prologue=Prologue("guard", "tox_cons_impl")) if prologue else Directives()
    procedures = [procedure(
        "staged_impl",
        integer("n_items", doc="items"),
        integer("factor", doc="a factor"),
        integer("limit", doc="a limit", directives=auto("scaled", "scaled_size", "tox_prod")),
        *arguments,
        real("total", Intent.OUT, doc="the total"),
        meta=IMPL, is_pure=True, directives=directives,
    )]
    if prologue:
        procedures.append(procedure(
            "guard", *guard_arguments,
            logical("handled", Intent.OUT, doc="never"),
            integer("ierr", Intent.INOUT, doc="left as it is"),
            meta=IMPL, is_pure=True,
        ))
    return module("tox_cons_impl", *procedures)


class TestTheChecksAreStaged:
    """Rule V's emission: the late bound alone runs below the producers, in its own block."""

    def test_a_bound_naming_a_produced_value_waits_for_the_producers(self):
        synthesis, bag = analysed(producers(), staged(
            real("x", Intent.IN, "(n_items)", doc="values",
                 directives=Directives(maximum=Maximum("real(limit, real64)")))))

        assert errors(bag) == []
        text = body(fortran(synthesis), "staged")
        top, after = text.split("call scaled_size(", 1)
        assert "call validate_all_in_range_real(x, n_items, ierr, arg_pos=3_int32)" in top
        assert "max=real(limit" not in top
        assert ("call validate_all_in_range_real(x, n_items, ierr, arg_pos=3_int32, "
                "max=real(limit, real64), allow_nan=.true._c_bool, allow_infinite=.true._c_bool)"
                in after.split("call staged_impl(")[0])
        # its own block, after the producer's error return
        assert "end if\n#ifndef NO_INPUT_VALIDATION" in after

    def test_the_expert_checks_everything_at_the_top(self):
        synthesis, _ = analysed(producers(), staged(
            real("x", Intent.IN, "(n_items)", doc="values",
                 directives=Directives(maximum=Maximum("real(limit, real64)")))))
        text = body(fortran(synthesis), "staged_expert")

        assert "max=real(limit, real64))" in text
        assert "allow_nan" not in text

    def test_an_integer_whose_only_bound_waits_is_not_checked_at_the_top(self):
        synthesis, _ = analysed(producers(), staged(
            integer("k", doc="a count", directives=Directives(maximum=Maximum("limit")))))
        top = body(fortran(synthesis), "staged").split("call scaled_size(", 1)[0]

        assert "validate_in_range_int(k" not in top


class TestWhatABoundMayName:
    """Rule V's refusals: a name set at none of the three points, at the bound's line."""

    def refused(self, *arguments, **kwargs):
        _, bag = validated(producers(), staged(*arguments, **kwargs))
        return [(e.message, e.note, e.location.line) for e in bag.errors]

    def test_a_permutation_foo_builds(self):
        [(message, note, line)] = self.refused(
            real("x", Intent.IN, "(n_items)", doc="values"),
            integer("x_perm", Intent.IN, "(n_items)", doc="an order"),
            integer("k", doc="an index", directives=Directives(minimum=Minimum("x_perm(1)", 12))))

        assert message == "'k' is checked against 'x_perm', which 'staged' builds itself"
        assert note == ("a bound names an argument the caller passes, a value a DM_OUTPUT_FROM "
                        "producer computes, a value the prologue writes, or a constant")
        assert line == 12

    def test_a_work_array_and_its_size(self):
        [(message, note, _)] = self.refused(
            real("tmp_w", Intent.OUT, "(n_items)", doc="scratch"),
            integer("k", doc="an index", directives=Directives(maximum=Maximum("size(tmp_w)"))))

        assert message == "'k' is checked against 'tmp_w', which 'staged' builds itself"
        assert note.endswith("; for a work array's size, name its extent")

    def test_an_output_as_a_value(self):
        [(message, _, _)] = self.refused(
            integer("k", doc="a value", directives=Directives(maximum=Maximum("total"))))

        assert message == ("'k' is checked against 'total', an output, which has no value until "
                           "the implementation sets it")

    def test_the_size_of_an_output_foo_keeps_is_fine(self):
        assert self.refused(
            real("out", Intent.OUT, "(n_items)", doc="a result"),
            integer("k", doc="a value", directives=Directives(maximum=Maximum("size(out)")))) == []

    def test_an_argument_only_the_prologue_has(self):
        [(message, _, _)] = self.refused(
            integer("k", doc="a value", directives=Directives(maximum=Maximum("fraction"))),
            prologue=True, guard_arguments=(real("fraction", doc="a fraction"),))

        assert message == ("'k' is checked against 'fraction', which only the prologue takes, so "
                           "'staged_expert' has no 'fraction'")

    def test_a_sentinel_naming_a_value_set_later(self):
        [(message, note, _)] = self.refused(
            integer("k", doc="a value", directives=Directives(
                minimum=Minimum("0_int32"), sentinel=Sentinel("limit"))))

        assert message == "'k' has a DM_SENTINEL naming 'limit', which 'staged' computes"
        assert note.endswith("; a sentinel is a fixed marker")

    def test_a_bound_after_the_prologue_on_what_the_prologue_writes(self):
        # fv_self: the caller's y is within the bound, the prologue's is not, and the error
        # would blame the caller after their y had been overwritten
        [(message, _, _)] = self.refused(
            real("level", doc="a level"),
            real("y", Intent.INOUT, doc="refined", directives=Directives(maximum=Maximum("level"))),
            prologue=True, guard_arguments=(real("level", Intent.OUT, doc="the level"),
                                            real("y", Intent.INOUT, doc="refined")))

        assert message == ("'y' is checked against 'level' once the prologue has run, but the "
                           "prologue writes 'y' too, so its error would blame the caller for the "
                           "prologue's value")

    def test_what_each_point_sets_is_accepted(self):
        assert self.refused(
            real("level", doc="a level"),
            integer("y", doc="between", directives=Directives(
                minimum=Minimum("limit"), maximum=Maximum("int(level, int32)"))),
            integer("m", doc="the limit", directives=Directives(equals=Equals("limit"))),
            prologue=True, guard_arguments=(real("level", Intent.OUT, doc="the level"),)) == []


class TestAComputedLengthOfAPassedString:
    """3b: a value a producer computes is not the length of a string the caller passes."""

    def test_refused_at_the_directive(self):
        cons = module("tox_cons_impl", procedure(
            "len_impl",
            integer("n_items", doc="items"),
            integer("str_len", doc="the length", directives=Directives(output_from=OutputFrom(
                "scaled", "scaled_size", "tox_prod", OutputFromMode.AUTO, 8))),
            character("names", Intent.IN, "(n_items)", length="str_len", doc="the names"),
            integer("factor", doc="a factor"),
            real("total", Intent.OUT, doc="the total"),
            meta=IMPL, is_pure=True,
        ))
        _, bag = validated(producers(), cons)

        [error] = [e for e in bag.errors if "length" in e.message]
        assert error.message == ("'str_len' is computed with 'scaled_size', but it is the length "
                                 "of 'names', which the caller passes: the strings carry their own "
                                 "length")
        assert error.note.startswith("a computed length is for strings the wrapper allocates")
        assert error.location.line == 8

    def test_an_output_s_computed_length_stays(self):
        cons = module("tox_cons_impl", procedure(
            "len_impl",
            integer("n_items", doc="items"),
            integer("str_len", doc="the length",
                    directives=auto("scaled", "scaled_size", "tox_prod")),
            integer("factor", doc="a factor"),
            character("names", Intent.OUT, "(n_items)", length="str_len", doc="the names"),
            meta=IMPL, is_pure=True,
        ))
        _, bag = validated(producers(), cons)

        assert [e.message for e in bag.errors if "length" in e.message] == []


class TestAPrologueKeepsAProducedValue:
    """Item 5: a prologue's intent(out) does not write over a value a producer computes; an
    intent(inout) one refines it (Q4)."""

    def test_an_intent_out_is_refused_once_at_the_directive(self):
        # the AUTO size of a work array: one error, not this and the read-first rule's too
        cons = staged(real("tmp_w", Intent.OUT, "(limit)", doc="scratch"), prologue=True,
                      guard_arguments=(integer("limit", Intent.OUT, doc="overwritten"),))
        limit = cons.procedures[0].argument("limit")
        limit.directives = Directives(output_from=OutputFrom(
            "scaled", "scaled_size", "tox_prod", OutputFromMode.AUTO, 21))
        _, bag = validated(producers(), cons)

        [error] = bag.errors
        assert error.location.line == 21
        assert error.message == ("the prologue 'guard' writes 'limit' (intent(out)), which "
                                 "'staged' computes with 'scaled_size' first, so the producer's "
                                 "value is discarded")
        assert error.note.startswith("a value has one source")

    def test_an_intent_inout_refines_it(self):
        _, bag = validated(producers(), staged(
            integer("y", doc="a value"),
            prologue=True, guard_arguments=(integer("limit", Intent.INOUT, doc="refined"),)))

        assert bag.errors == ()


class TestAKeptValueOfAPairIsReported:
    """Once per implementation, also where `foo` keeps the value: the plain half resolves into a
    throwaway bag, and must not mark the report as made, or the expert's is suppressed too."""

    def pair(self, rows, producer="padded_size", output="padded"):
        work = integer("out_size", doc="the size", directives=Directives(
            output_from=OutputFrom(output, producer, "tox_prod", OutputFromMode.AUTO, 7)))
        located(work, table(*rows))
        return module("tox_cons_impl", procedure(
            "pair_impl", integer("n_items", doc="items"),
            real("tmp_buf", Intent.OUT, "(n_items)", doc="scratch, which makes this a pair"),
            work, real("filled", Intent.OUT, "(out_size)", doc="the result"),
            meta=IMPL, is_pure=True,
        ))

    @pytest.mark.parametrize("rows, message", [
        ([("n_items", "OMITTED")],
         "the producer-input table of 'out_size' omits 'n_items', which 'padded_size' requires"),
        ([("padded", "n_items")],
         "the producer-input table of 'out_size' maps 'padded', which is an output of "
         "'padded_size', not an input"),
        ([("padding", "nothing_here")],
         "the producer-input table of 'out_size' supplies 'padding' by 'nothing_here', which is "
         "neither an argument of 'pair_impl' nor a constant"),
        ([], "'out_size' is computed from 'padded_size', but its input(s) padding have no "
             "argument of the same name in 'pair_impl'"),
    ])
    def test_each_refusal_is_reported_once(self, rows, message):
        synthesis, bag = analysed(optional_producers(), self.pair(rows))

        # the premise: a pair, whose plain tier keeps `out_size`
        assert synthesis.project.procedure("tox_cons", "pair_expert") is not None
        assert synthesis.project.procedure("tox_cons", "pair").argument("out_size") is not None
        assert [e.message for e in bag.errors] == [message]


class TestAConstantRowIsReportedAtItsRow:
    """M2: a refusal caused by one table row is reported at that row, once per implementation,
    T's constant rows included -- a split consumer finds the row in every per-mode wrapper."""

    def split(self, row):
        modes = [
            "Which way.",
            "| Mode | Value | Procedure |",
            "|------|-------|-----------|",
            "| one | [[tox_cons_impl(module):MODE_ONE(variable)]] | use_one |",
            "| two | [[tox_cons_impl(module):MODE_TWO(variable)]] | use_two |",
        ]
        work = integer("work_size", doc="the size", directives=Directives(output_from=OutputFrom(
            "sized", "factor_size", "tox_prod", OutputFromMode.AUTO, 9)))
        located(work, table(row))
        return module("tox_cons_impl", procedure(
            "use_impl", integer("mode", doc=modes), integer("n_items", doc="items"), work,
            real("tmp_buf", Intent.OUT, "(work_size)", doc="scratch"), meta=IMPL, is_pure=True,
        ), parameters=(parameter("MODE_ONE", "1_int32"), parameter("MODE_TWO", "2_int32")))

    def test_once_at_the_row(self):
        synthesis, bag = analysed(kinded_producers(), self.split(("factor", ".true.")))

        assert synthesis.project.procedure("tox_cons", "use_one") is not None  # the premise
        [error] = bag.errors
        assert error.message.endswith("but the table gives it the constant .true., a logical")
        assert error.location.line == 103

    def test_a_mode_input_given_a_number_too(self):
        cons = self.split(("method", "2_int32"))
        work = cons.procedures[0].argument("work_size")
        work.directives = Directives(output_from=OutputFrom(
            "sized", "method_size", "tox_prod", OutputFromMode.AUTO, 9))
        _, bag = analysed(kinded_producers(), cons)

        [error] = bag.errors
        assert error.location.line == 103


class TestAStringConstantOfAFixedLength:
    """M1 and N4: a string constant shorter than a fixed length is padded to it, as Fortran
    assigns and the bindings pass it; a longer one is refused, for a length named by a
    parameter too."""

    def producers(self, length):
        return module(
            "tox_prod_impl",
            procedure("label_size_impl", integer("n_items", doc="items"),
                      character("label", length=length, doc="a label"),
                      integer("sized", Intent.OUT, doc="the size"), meta=IMPL, is_pure=True),
            parameters=(parameter("LABEL_LEN", "6_int32"),),
        )

    @pytest.mark.parametrize("length", ["6", "LABEL_LEN"])
    def test_a_shorter_one_is_padded(self, length):
        synthesis, bag = analysed(self.producers(length), kind_consumer(
            "label_size", rows=[("label", "'abc'")]))

        assert errors(bag) == []
        assert "label = 'abc   '" in fortran(synthesis)

    @pytest.mark.parametrize("length", ["6", "LABEL_LEN"])
    def test_a_longer_one_is_refused(self, length):
        [(message, _, _)] = refusals(self.producers(length), kind_consumer(
            "label_size", rows=[("label", "'abcdefg'")]))

        assert message.endswith(f"which is longer than 6")


class TestNothingLeftToCheckAtTheTop:
    """N3: a real allowing NaN and infinity whose bounds all wait has no call at the top."""

    def test_no_empty_call(self):
        synthesis, bag = analysed(producers(), staged(
            real("y", doc="a value", directives=Directives(
                maximum=Maximum("real(limit, real64)"), allow_nan=AllowNan(),
                allow_infinite=AllowInfinite()))))

        assert errors(bag) == []
        top, after = body(fortran(synthesis), "staged").split("call scaled_size(", 1)
        assert "validate_in_range_real(y" not in top
        assert "validate_in_range_real(y, ierr, arg_pos=3_int32, max=real(limit, real64)" in after


class TestABoundOfALoneWrapper:
    """S1: a lone wrapper's bound naming an output reads what nothing set, in every tier."""

    def test_an_output_named_as_a_value_is_refused(self):
        cons = module("tox_cons_impl", procedure(
            "lone_impl",
            integer("k", doc="a value", directives=Directives(maximum=Maximum("n_out", 5))),
            integer("n_out", Intent.OUT, doc="a count"),
            meta=IMPL, is_pure=True,
        ))
        synthesis, bag = validated(cons)

        assert synthesis.project.procedure("tox_cons", "lone_expert") is None  # the premise
        [error] = bag.errors
        assert error.message == ("'k' is checked against 'n_out', an output, which has no value "
                                 "until the implementation sets it")
        assert error.location.line == 5


class TestANamedLengthIsTheProducersOwn:
    """A producer input's `len=LABEL_LEN` is evaluated in the producer's own module: the
    consumer's module may declare a `LABEL_LEN` of its own, which the project's flat table of
    constants would give instead (`fx/m1clash` of the second review)."""

    def pair(self, producer_len, consumer_len):
        producer = module(
            "tox_prod_impl",
            procedure("label_size_impl", integer("n_items", doc="items"),
                      character("label", length="LABEL_LEN", doc="a label"),
                      integer("sized", Intent.OUT, doc="the size"), meta=IMPL, is_pure=True),
            parameters=(parameter("LABEL_LEN", f"{producer_len}_int32"),),
        )
        work = integer("work_size", doc="the size", directives=Directives(output_from=OutputFrom(
            "sized", "label_size", "tox_prod", OutputFromMode.AUTO, 9)))
        work.doc = Doc.parse(table(("label", "'abc'")))
        # named to sort after the producer's module, so its LABEL_LEN is the flat table's
        consumer = module("tox_zz_impl", procedure(
            "use_impl", integer("n_items", doc="items"), work,
            real("tmp_buf", Intent.OUT, "(work_size)", doc="scratch"), meta=IMPL, is_pure=True,
        ), parameters=(parameter("LABEL_LEN", f"{consumer_len}_int32"),))
        return producer, consumer

    def test_a_longer_one_of_the_consumer_s_does_not_refuse_it(self):
        synthesis, bag = analysed(*self.pair(6, 2))

        assert errors(bag) == []
        assert "label = 'abc   '" in fortran(synthesis, "tox_zz")

    def test_a_shorter_one_of_the_producer_s_refuses_it(self):
        [(message, _, _)] = refusals(*self.pair(2, 6))

        assert message.endswith("which is longer than 2")
