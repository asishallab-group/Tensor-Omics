"""`DM_OUTPUT_FROM` naming a producer that is itself a generated wrapper.

The producer is a published entry point: a hand-written `M_EXPORT_C` procedure, or the plain-named
wrapper an implementation generates. A generated one validates its inputs, so it can fail, and
reports the failure by its own argument's position; these tests pin what the consumer makes of
that, and the rest of what calling one directly takes (a same-module producer, its outputs nobody
takes, the extents its binding derives, its purity). The same cases are compiled and run through
Python and R by `test_end_to_end_impl.py` (`fixtures/impl/fx_producers_impl.F90`).
"""

import pytest

from codegen.config import CONVENTIONS
from codegen.diagnostics import DiagnosticBag
from codegen.emit.fortran_wrapper import FortranWrapperEmitter, wrapper_info
from codegen.ir.directives import Directives, Minimum, OutputFrom, OutputFromMode
from codegen.ir.entities import Meta
from codegen.ir.roles import analyse_project
from codegen.ir.types import Intent
from codegen.synthesize import synthesize_wrappers

from builders import C_BINDING, integer, module, procedure, project, real

#: What an implementation carries: no export, it is not published itself
IMPL = Meta(summary="an implementation", author="AUTHOR")


def auto(output, producer, producer_module):
    return Directives(output_from=OutputFrom(output, producer, producer_module, OutputFromMode.AUTO))


def table(*rows):
    lines = ["Where it comes from.", "| Producer input | Supplied by |", "|---|---|"]
    return lines + [f"| {left} | {right} |" for left, right in rows]


def producers(*extra, pure=True):
    """`tox_prod_impl`: implementations whose generated wrappers in `tox_prod` are producers."""
    return module(
        "tox_prod_impl",
        procedure(
            "scaled_size_impl",
            integer("n_items", doc="items", directives=Directives(minimum=Minimum("1_int32"))),
            integer("factor", doc="factor", directives=Directives(minimum=Minimum("1_int32"))),
            integer("scaled", Intent.OUT, doc="the size"),
            meta=IMPL, is_pure=pure,
        ),
        procedure(
            "extent_size_impl",
            integer("n", doc="values"),
            real("x", Intent.IN, "(n)", doc="the values"),
            integer("extent", Intent.OUT, doc="the size"),
            meta=IMPL, is_pure=True,
        ),
        procedure(
            "two_sizes_impl",
            integer("n_items", doc="items"),
            integer("first", Intent.OUT, doc="one size"),
            integer("second", Intent.OUT, doc="another"),
            meta=IMPL, is_pure=True,
        ),
        procedure(
            "sized_and_rows_impl",
            integer("n_items", doc="items"),
            integer("sized", Intent.OUT, doc="the size"),
            real("rows", Intent.OUT, "(n_items)", doc="an array output"),
            meta=IMPL, is_pure=True,
        ),
        *extra,
    )


def split_producer(pure=True):
    """`tox_split_impl`: an implementation split into one wrapper per mode, in `tox_split`."""
    modes = [
        "Which multiple.",
        "| Mode | Value | Procedure |",
        "|------|-------|-----------|",
        "| double | [[tox_split_impl(module):MODE_DOUBLE(variable)]] | split_size_double |",
        "| triple | [[tox_split_impl(module):MODE_TRIPLE(variable)]] | split_size_triple |",
    ]
    return module(
        "tox_split_impl",
        procedure(
            "split_size_impl",
            integer("n_items", doc="items"),
            integer("size_mode", doc=modes),
            integer("sized", Intent.OUT, doc="the size"),
            meta=IMPL, is_pure=pure,
        ),
    )


def consumer(name, *arguments, module_name="tox_cons_impl"):
    """An implementation module holding one consumer implementation."""
    return module(module_name, procedure(name, *arguments, meta=IMPL, is_pure=True))


def work_consumer(producer="scaled_size", producer_module="tox_prod", output="scaled",
                  doc=table(("factor", "scale")), name="work_total_impl"):
    """A consumer whose AUTO value sizes a work array: the generated wrapper calls the producer."""
    return consumer(
        name,
        integer("n_items", doc="items"),
        real("x", Intent.IN, "(n_items)", doc="values"),
        integer("scale", doc="scale"),
        integer("work_size", doc=doc, directives=auto(output, producer, producer_module)),
        real("tmp_buf", Intent.OUT, "(work_size)", doc="scratch"),
        real("total", Intent.OUT, doc="the result"),
    )


def result_consumer(producer="scaled_size", output="scaled", doc=table(("factor", "scale"))):
    """A consumer whose AUTO value sizes a result: Python and R call the producer's binding."""
    return consumer(
        "fill_impl",
        integer("n_items", doc="items"),
        real("x", Intent.IN, "(n_items)", doc="values"),
        integer("scale", doc="scale"),
        integer("out_size", doc=doc, directives=auto(output, producer, "tox_prod")),
        real("filled", Intent.OUT, "(out_size)", doc="the result"),
    )


def analysed(*modules):
    bag = DiagnosticBag()
    synthesis = synthesize_wrappers(project(*modules), CONVENTIONS, bag)
    analyse_project(synthesis.project, bag, CONVENTIONS)
    return synthesis, bag


def fortran(synthesis, module_name="tox_cons"):
    emitter = FortranWrapperEmitter(project=synthesis.project, wrapper_info=wrapper_info(synthesis))
    info = wrapper_info(synthesis)
    return emitter.module(synthesis.project.module(module_name), info.get(module_name))


def body(text, name):
    start = text.index(f"subroutine {name}(")
    return text[start:text.index(f"end subroutine {name}\n", start)]


def call_of(text, callee):
    start = text.index(f"call {callee}(&")
    return text[start:text.index(")\n", start)]


class TestTheProducerIsAPublishedEntryPoint:
    """Rule 1: `foo`, the plain-named wrapper in the generated module -- never `foo_impl`."""

    def test_the_generated_wrapper_is_the_producer(self):
        synthesis, bag = analysed(producers(), work_consumer())

        plan = synthesis.project.procedure("tox_cons", "work_total_expert") \
            .argument("work_size").roles.computed_from
        assert bag.errors == ()
        assert (plan.producer.module.name, plan.producer.name) == ("tox_prod", "scaled_size")

    def test_naming_the_implementation_is_refused_and_names_its_wrapper(self):
        _, bag = analysed(producers(), work_consumer("scaled_size_impl", "tox_prod_impl"))

        error = next(e for e in bag.errors if "not exported" in e.message)
        assert "'scaled_size_impl'" in error.message
        assert error.note == ("an implementation is published as its generated wrapper: name "
                              "'scaled_size' in 'tox_prod'")

    def test_the_implementation_with_its_generated_module_is_refused_and_names_it(self):
        # `scaled_size_impl` in `tox_prod`: there is no such procedure there
        _, bag = analysed(producers(), work_consumer("scaled_size_impl", "tox_prod"))

        error = next(e for e in bag.errors if "does not exist" in e.message)
        assert error.note == ("an implementation is published as its generated wrapper: name "
                              "'scaled_size' in 'tox_prod'")

    def test_the_wrapper_with_the_implementation_module_is_refused_and_named(self):
        # `scaled_size` in `tox_prod_impl`: the wrapper's name with the implementation's module
        _, bag = analysed(producers(), work_consumer("scaled_size", "tox_prod_impl"))

        error = next(e for e in bag.errors if "does not exist" in e.message)
        assert error.note == ("'scaled_size' is generated into 'tox_prod', not into "
                              "'tox_prod_impl': name it there")

    def test_a_dangling_link_beside_a_missing_producer_is_still_reported(self):
        # the refusal of the producer covers the directive; the link the macro renders into the
        # argument's documentation dangles too, and is reported as itself
        from codegen.ir.validate import validate_project

        link = ["Size, produced by [[tox_prod_impl(module):scaled_size]]."]
        cons = work_consumer("scaled_size", "tox_prod_impl", doc=link + table(("factor", "scale"))[1:])
        synthesis, bag = analysed(producers(), cons)
        validate_project(synthesis.project, bag)

        assert any("does not exist" in e.message for e in bag.errors)
        assert any("documentation links to 'scaled_size', which 'tox_prod_impl' does not publish"
                   in e.message for e in bag.errors)

    def test_the_expert_with_the_implementation_module_is_pointed_at_the_plain_one(self):
        middle = procedure(
            "middle_impl",
            integer("n_items", doc="items"),
            integer("scale", doc="scale"),
            integer("work_size", doc=table(("factor", "scale")),
                    directives=auto("scaled", "scaled_size", "tox_prod")),
            real("tmp_buf", Intent.OUT, "(work_size)", doc="scratch"),
            integer("middle_size", Intent.OUT, doc="the size"),
            meta=IMPL, is_pure=True,
        )
        cons = work_consumer("middle_expert", "tox_prod_impl", "middle_size", doc=table())
        _, bag = analysed(producers(middle), cons)

        error = next(e for e in bag.errors if "does not exist" in e.message)
        assert error.note == ("'middle' is generated into 'tox_prod', not into "
                              "'tox_prod_impl': name it there")

    def test_the_expert_tier_is_refused_and_names_the_plain_one(self):
        # `middle` takes a work array over, so it has a `middle_expert` -- published, but not
        # what a producer call may name
        middle = procedure(
            "middle_impl",
            integer("n_items", doc="items"),
            integer("scale", doc="scale"),
            integer("work_size", doc=table(("factor", "scale")),
                    directives=auto("scaled", "scaled_size", "tox_prod")),
            real("tmp_buf", Intent.OUT, "(work_size)", doc="scratch"),
            integer("middle_size", Intent.OUT, doc="the size"),
            meta=IMPL, is_pure=True,
        )
        cons = work_consumer("middle_expert", "tox_prod", "middle_size", doc=table())
        _, bag = analysed(producers(middle), cons)

        error = next(e for e in bag.errors if "middle_expert" in e.message)
        assert "the expert tier" in error.message
        assert error.note == ("name the plain-named wrapper, which prepares it itself: "
                              "'middle' in 'tox_prod'")

    def test_a_hand_written_expert_name_is_the_author_s_own(self):
        hand = module("tox_hand", procedure(
            "sizer_expert", integer("n_items", doc="n"), integer("scaled", Intent.OUT, doc="s"),
            integer("ierr", Intent.OUT, doc="e")))
        _, bag = analysed(hand, work_consumer("sizer_expert", "tox_hand", doc=table()))

        assert bag.errors == ()

    def test_a_mode_split_implementation_names_its_per_mode_wrappers(self):
        cons = work_consumer("split_size_impl", "tox_split_impl", "sized", doc=table())
        _, bag = analysed(split_producer(), cons)

        error = next(e for e in bag.errors if "not exported" in e.message)
        assert error.note == ("'split_size_impl' splits by mode, into one wrapper per mode: name "
                              "the one you mean in 'tox_split' ('split_size_double', "
                              "'split_size_triple')")

    def test_a_per_mode_wrapper_is_a_producer(self):
        cons = work_consumer("split_size_double", "tox_split", "sized", doc=table())
        synthesis, bag = analysed(split_producer(), cons)

        assert bag.errors == ()
        assert "call split_size_double(&" in fortran(synthesis)

    def test_a_hand_written_procedure_is_told_to_export_itself(self):
        hand = module(
            "tox_hand",
            procedure("sizer", integer("n_items", doc="n"), integer("scaled", Intent.OUT, doc="s"),
                      integer("ierr", Intent.OUT, doc="e"), meta=Meta()),
        )
        _, bag = analysed(hand, work_consumer("sizer", "tox_hand", doc=table(("n_items", "n_items"))))

        error = next(e for e in bag.errors if "not exported" in e.message)
        assert error.note == "mark 'sizer' M_EXPORT_C to export it"

    def test_just_info_names_the_wrapper_too(self):
        # documentation only, but the link it renders is what a reader follows
        cons = consumer(
            "info_impl",
            integer("n_items", doc="items"),
            integer("out_size", doc="size", directives=Directives(output_from=OutputFrom(
                "scaled", "scaled_size_impl", "tox_prod_impl", OutputFromMode.JUST_INFO))),
            real("ones", Intent.OUT, "(out_size)", doc="ones"),
        )
        _, bag = analysed(producers(), cons)

        error = next(e for e in bag.errors if "not exported" in e.message)
        assert "name 'scaled_size' in 'tox_prod'" in error.note


class TestEachProducerCallIsChecked:
    """Rule 2: `is_err` after each call, its position mapped to the dummy that supplied it."""

    def two_producers(self):
        return consumer(
            "two_impl",
            integer("n_items", doc="items"),
            real("x", Intent.IN, "(n_items)", doc="values"),
            integer("scale", doc="scale"),
            integer("first_size", doc="size", directives=auto("first", "two_sizes", "tox_prod")),
            real("tmp_first", Intent.OUT, "(first_size)", doc="scratch"),
            integer("second_size", doc=table(("factor", "scale")),
                    directives=auto("scaled", "scaled_size", "tox_prod")),
            real("tmp_second", Intent.OUT, "(second_size)", doc="scratch"),
            real("total", Intent.OUT, doc="the result"),
        )

    def test_the_first_producer_is_checked_before_the_second_is_called(self):
        synthesis, bag = analysed(producers(), self.two_producers())
        text = body(fortran(synthesis), "two")

        assert bag.errors == ()
        first, second = text.index("call two_sizes(&"), text.index("call scaled_size(&")
        assert first < text.index("if (is_err(ierr)) then", first) < second
        assert text.count("if (is_err(ierr)) then") == 2

    def test_a_position_is_mapped_to_the_dummy_that_supplied_it(self):
        # scaled_size(n_items, factor, ...) is passed two(n_items, x, scale, total, ierr)'s
        # n_items (1 -> 1) and scale (2 -> 3)
        synthesis, _ = analysed(producers(), self.two_producers())
        text = body(fortran(synthesis), "two")
        check = text[text.index("call scaled_size(&"):]
        check = check[check.index("if (is_err(ierr)) then"):check.index("end if")]

        assert "select case (get_err_arg_pos(ierr))" in check
        # n_items keeps its position: kept, not rewritten -- the default would clear it
        assert "case (1_int32)\n" in check
        assert "map_err_arg_pos(ierr, 1_int32" not in check
        assert "call map_err_arg_pos(ierr, 2_int32, 3_int32)" in check
        assert "case default\n" in check and "call clear_err_arg_pos(ierr)" in check
        assert check.rstrip().endswith("return")

    def test_a_constant_supplies_nothing_to_map_so_the_position_is_cleared(self):
        cons = work_consumer(doc=table(("factor", "2_int32")))
        synthesis, bag = analysed(producers(), cons)
        text = body(fortran(synthesis), "work_total")
        check = text[text.index("call scaled_size(&"):]
        check = check[check.index("if (is_err(ierr)) then"):check.index("end if")]

        assert bag.errors == ()
        assert "case (1_int32)" in check  # n_items is still a dummy
        assert "case (2_int32)" not in check  # factor, the constant, falls to the default
        assert "case default\n" in check

    def test_the_helpers_are_imported(self):
        synthesis, _ = analysed(producers(), work_consumer())
        header = fortran(synthesis).split("contains", 1)[0]

        for name in ("get_err_arg_pos", "map_err_arg_pos", "clear_err_arg_pos"):
            assert name in header

    def test_one_producer_given_two_sets_of_inputs_is_called_once_for_each(self):
        cons = consumer(
            "twice_impl",
            integer("n_items", doc="items"),
            integer("scale_a", doc="scale"),
            integer("scale_b", doc="scale"),
            integer("size_a", doc=table(("factor", "scale_a")),
                    directives=auto("scaled", "scaled_size", "tox_prod")),
            real("tmp_a", Intent.OUT, "(size_a)", doc="scratch"),
            integer("size_b", doc=table(("factor", "scale_b")),
                    directives=auto("scaled", "scaled_size", "tox_prod")),
            real("tmp_b", Intent.OUT, "(size_b)", doc="scratch"),
            real("total", Intent.OUT, doc="the result"),
        )
        synthesis, bag = analysed(producers(), cons)
        text = body(fortran(synthesis), "twice")

        assert bag.errors == ()
        assert text.count("call scaled_size(&") == 2
        first, second = text.split("call scaled_size(&")[1:]
        assert "factor = scale_a" in first and "scaled = size_a" in first
        assert "factor = scale_b" in second and "scaled = size_b" in second

    def test_two_arguments_taking_one_call_s_same_output_are_refused_on_either_route(self):
        # one call, `factor | scale` for both: one value, which the Fortran call wrote once
        for sized in ("tmp", "res"):
            cons = consumer(
                "same_impl",
                integer("n_items", doc="items"),
                integer("scale", doc="scale"),
                integer("size_a", doc=table(("factor", "scale")),
                        directives=auto("scaled", "scaled_size", "tox_prod")),
                real(f"{sized}_a", Intent.OUT, "(size_a)", doc="sized"),
                integer("size_b", doc=table(("factor", "scale")),
                        directives=auto("scaled", "scaled_size", "tox_prod")),
                real(f"{sized}_b", Intent.OUT, "(size_b)", doc="sized"),
            )
            _, bag = analysed(producers(), cons)

            [error] = [e for e in bag.errors if "takes the same output" in e.message]
            assert error.message == ("'size_b' takes the same output 'scaled' of 'scaled_size', "
                                     "given the same inputs, as 'size_a'")
            assert error.note == ("one call fills both with one value: use 'size_a' alone, "
                                  "since 'size_b' would hold the same value")

    def test_a_producer_without_ierr_is_not_checked(self):
        hand = module(
            "tox_hand",
            procedure("sizer", integer("n_items", doc="n"), integer("scaled", Intent.OUT, doc="s")),
        )
        synthesis, bag = analysed(hand, work_consumer("sizer", "tox_hand", doc=table()))
        text = body(fortran(synthesis), "work_total")

        assert bag.errors == ()
        after = text[text.index("call sizer(&"):text.index("M_ALLOCATE(tmp_buf")]
        assert "is_err" not in after


class TestTheBindingsNameTheirOwnArgument:
    """Rule 3: Python and R raise a producer's error again, in the consumer's words."""

    def binding(self):
        from codegen.abi.c_abi import build_project

        synthesis, bag = analysed(producers(), result_consumer())
        binding = build_project(synthesis.project, bag)
        assert bag.errors == (), bag.render()
        return next(w for m in binding for w in m if w.stripped_name == "fill")

    def test_python_maps_the_producer_s_arguments_to_its_own_positions(self):
        from codegen.emit.python_ctypes import PythonEmitter

        text = PythonEmitter().function(self.binding())
        # fill(n_items, x, scale, out_size, filled, ierr): factor is supplied by scale (3rd)
        assert "    out_size = scaled_size(n_items=n_items, factor=scale)\n" in text
        assert ('raise _consumer_error(_error, {"n_items": 1, "factor": 3}, _FILL_ARGUMENTS, '
                "_FILL_ARGUMENT_SOURCES) from None") in text
        assert text.index("try:") < text.index("scaled_size(") < text.index("except _ToxError")

    def test_r_maps_them_too(self):
        from codegen.emit.r_wrapper import RWrapperEmitter

        text = RWrapperEmitter().function(self.binding())
        assert ("out_size <- tryCatch(scaled_size(n_items = length(x), factor = scale), "
                "tox_error = function(e) .tox_consumer_error(e, c(n_items = 1L, factor = 3L), "
                ".arguments, .sources))") in text
        # what the producer's error is decoded against is set before the producer is called
        assert text.index(".arguments <- c(") < text.index("tryCatch(")
        assert text.count(".arguments <- c(") == 1


class TestPurity:
    """Rule 4: a generated producer is as pure as `_is_pure` emits it, asked recursively."""

    def test_a_pure_generated_producer_leaves_its_consumer_pure(self):
        synthesis, _ = analysed(producers(), work_consumer())

        assert "pure subroutine work_total(" in fortran(synthesis)

    def test_an_impure_generated_producer_makes_it_impure(self):
        synthesis, _ = analysed(producers(pure=False), work_consumer())
        text = fortran(synthesis)

        assert "subroutine work_total(" in text
        assert "pure subroutine work_total(" not in text

    def test_a_producer_that_is_a_consumer_is_asked_in_turn(self):
        # tox_prod's own work_total calls scaled_size; tox_cons's calls that one's wrapper
        middle = procedure(
            "middle_impl",
            integer("n_items", doc="items"),
            integer("scale", doc="scale"),
            integer("work_size", doc=table(("factor", "scale")),
                    directives=auto("scaled", "scaled_size", "tox_prod")),
            real("tmp_buf", Intent.OUT, "(work_size)", doc="scratch"),
            integer("middle_size", Intent.OUT, doc="the size"),
            meta=IMPL, is_pure=True,
        )
        cons = work_consumer("middle", "tox_prod", "middle_size", doc=table())
        for pure, expected in ((True, True), (False, False)):
            synthesis, bag = analysed(producers(middle, pure=pure), cons)
            assert bag.errors == (), bag.render()
            assert ("pure subroutine work_total(" in fortran(synthesis)) is expected

    def test_a_producer_that_reaches_back_to_its_consumer_ends(self):
        # work_total sizes itself by its own wrapper: refused (`TestACycle`), but nothing
        # proves that pure either, and asking must not recurse until Python gives up
        cons = work_consumer("work_total", "tox_cons", "total", doc=table())
        synthesis, _ = analysed(producers(), cons)

        text = fortran(synthesis)  # no RecursionError
        assert "pure subroutine work_total(" not in text

    def test_a_mode_split_producer_in_another_module_is_seen_through(self):
        # its implementation is named by the mode table, not by the wrapper's name, so it is
        # found only with the other module's wrapper info at hand
        cons = work_consumer("split_size_double", "tox_split", "sized", doc=table())
        for pure in (True, False):
            synthesis, bag = analysed(split_producer(pure=pure), cons)
            assert bag.errors == (), bag.render()
            assert ("pure subroutine work_total(" in fortran(synthesis)) is pure


class TestAnInputTheProducerRefines:
    """An `intent(inout)` producer input -- a cap the producer lowers in place.

    Where the generated Fortran calls the producer, it is passed through a required
    `intent(inout)` argument of the consumer: that is what the call can write and what carries
    the refined value on. Python and R adopt the value themselves.
    """

    MESSAGE = ("'work_size' is computed from 'capped_size', which refines its input 'cap' in "
               "place, so it must be passed through a required intent(inout) argument; it is "
               "given {}")

    def capping(self):
        return module("tox_hand", procedure(
            "capped_size",
            integer("n_items", doc="items"),
            integer("cap", Intent.INOUT, doc="a cap it lowers"),
            integer("sized", Intent.OUT, doc="the size"),
            integer("ierr", Intent.OUT, doc="error"),
        ))

    def capped(self, cap=None, sized_result=False, doc=None):
        size = integer("work_size", doc=doc or "size",
                       directives=auto("sized", "capped_size", "tox_hand"))
        sized = (real("res", Intent.OUT, "(work_size)", doc="the result") if sized_result
                 else real("tmp_buf", Intent.OUT, "(work_size)", doc="scratch"))
        caps = [cap] if cap is not None else []
        return consumer("capped_impl", integer("n_items", doc="items"), *caps, size, sized)

    def refusals(self, bag):
        return [e for e in bag.errors if "refines its input" in e.message]

    def refused_as(self, cons, given):
        _, bag = analysed(self.capping(), cons)

        [error] = self.refusals(bag)
        assert error.message == self.MESSAGE.format(given)
        assert error.note == (
            "the generated wrapper calls 'capped_size' itself, and writes 'cap' back through "
            "what supplies it: supply it by a required intent(inout) argument of the "
            "implementation, which the refined value then reaches, and the caller with it")

    def test_an_intent_in_argument_is_refused(self):
        self.refused_as(self.capped(integer("cap", Intent.IN, doc="a cap")), "intent(in) 'cap'")

    def test_a_constant_is_refused(self):
        # as written in the table, not as evaluated
        self.refused_as(self.capped(doc=table(("cap", "5_int32"))), "the constant 5_int32")

    def test_an_optional_with_a_default_is_refused(self):
        # its value reaches the call as a local, which would swallow the refined one
        from codegen.ir.directives import Default

        cap = integer("cap", Intent.IN, doc="a cap", optional=True,
                      directives=Directives(default=Default("20_int32")))
        self.refused_as(self.capped(cap), "optional intent(in) 'cap'")

    def test_an_optional_without_one_is_refused(self):
        cap = integer("cap", Intent.IN, doc="a cap", optional=True)
        self.refused_as(self.capped(cap), "optional intent(in) 'cap'")

    def test_an_optional_intent_inout_is_refused_too(self):
        cap = integer("cap", Intent.INOUT, doc="a cap", optional=True)
        self.refused_as(self.capped(cap), "optional intent(inout) 'cap'")

    def test_a_required_intent_inout_argument_is_passed_to_it(self):
        synthesis, bag = analysed(self.capping(),
                                  self.capped(integer("cap", Intent.INOUT, doc="a cap")))

        assert bag.errors == ()
        assert "cap = cap" in call_of(body(fortran(synthesis), "capped"), "capped_size")

    def test_the_bindings_adopt_the_refined_value_instead(self):
        # nothing in Fortran calls it: Python and R take the producer's returned cap
        _, bag = analysed(self.capping(),
                          self.capped(integer("cap", Intent.IN, doc="a cap"), sized_result=True))

        assert bag.errors == ()


class TestARefinedValueFeedsOneCall:
    """A value a producer refines in place feeds that one producer call, and no other.

    Another call reading it, or refining it again, is refused on either route: what it would
    see depended on the call order. One call filling several outputs from it is one call.
    """

    def producers(self):
        return module("tox_hand",
            procedure("capped_size", integer("n_items", doc="items"),
                      integer("cap", Intent.INOUT, doc="a cap it lowers"),
                      integer("sized", Intent.OUT, doc="the size"),
                      integer("ierr", Intent.OUT, doc="error")),
            procedure("capped_pair", integer("cap", Intent.INOUT, doc="a cap it lowers"),
                      integer("first", Intent.OUT, doc="a size"),
                      integer("second", Intent.OUT, doc="another"),
                      integer("ierr", Intent.OUT, doc="error")),
            procedure("plain_size", integer("base", doc="a base"),
                      integer("sized", Intent.OUT, doc="the size"),
                      integer("ierr", Intent.OUT, doc="error")),
        )

    def consumer_of(self, *sized, binding=False):
        """`cap` a required intent(inout) argument; each of `sized` is (name, producer, output,
        table rows), sizing a work array (the Fortran route) or a result (the binding route)."""
        arguments = [integer("n_items", doc="items"), integer("cap", Intent.INOUT, doc="a cap")]
        for name, producer, output, rows in sized:
            arguments.append(integer(name, doc=table(*rows),
                                     directives=auto(output, producer, "tox_hand")))
            arguments.append(real(("res_" if binding else "tmp_") + name, Intent.OUT, f"({name})",
                                  doc="sized"))
        return consumer("refining_impl", *arguments)

    def refusals(self, bag):
        return [e.message for e in bag.errors if "is refined in place by" in e.message]

    REFINER = ("s1", "capped_size", "sized", [])

    def test_a_second_refiner_is_refused_on_the_fortran_route(self):
        cons = self.consumer_of(self.REFINER, ("s2", "capped_pair", "first", []))
        _, bag = analysed(self.producers(), cons)

        assert self.refusals(bag) == [
            "'cap' is refined in place by 'capped_size' (computing 's1'), so it may feed no "
            "other producer call, but 'capped_pair' (computing 's2') refines it too"]

    def test_another_reader_is_refused_whatever_the_order(self):
        reader = ("s0", "plain_size", "sized", [("base", "cap")])
        for sized in ((reader, self.REFINER), (self.REFINER, reader)):
            _, bag = analysed(self.producers(), self.consumer_of(*sized))

            assert self.refusals(bag) == [
                "'cap' is refined in place by 'capped_size' (computing 's1'), so it may feed no "
                "other producer call, but 'plain_size' (computing 's0') reads it too"]

    def test_both_are_refused_on_the_binding_route_too(self):
        reader = ("s0", "plain_size", "sized", [("base", "cap")])
        for other, verb in ((reader, "reads"), (("s2", "capped_pair", "first", []), "refines")):
            _, bag = analysed(self.producers(),
                              self.consumer_of(self.REFINER, other, binding=True))

            assert len(self.refusals(bag)) == 1 and self.refusals(bag)[0].endswith(f"{verb} it too")

    def test_a_hand_written_consumer_beside_its_own_expert_named_procedure_is_checked(self):
        # `refining` and an unrelated hand-written `refining_expert`: not a generated pair, so
        # nothing else checks `refining`'s arguments
        def hand(binding_name):
            return procedure(binding_name, integer("n_items", doc="items"),
                             integer("total", Intent.OUT, doc="a total"),
                             integer("ierr", Intent.OUT, doc="error"))
        reader = ("s0", "plain_size", "sized", [("base", "cap")])
        arguments = [integer("n_items", doc="items"), integer("cap", Intent.INOUT, doc="a cap")]
        for name, producer, output, rows in (self.REFINER, reader):
            arguments.append(integer(name, doc=table(*rows),
                                     directives=auto(output, producer, "tox_hand")))
            arguments.append(real("res_" + name, Intent.OUT, f"({name})", doc="sized"))
        arguments.append(integer("ierr", Intent.OUT, doc="error"))
        cons = module("tox_cons", procedure("refining", *arguments), hand("refining_expert"))
        _, bag = analysed(self.producers(), cons)

        assert len(self.refusals(bag)) == 1

    def test_one_call_filling_two_outputs_from_it_is_one_call(self):
        for binding in (False, True):
            _, bag = analysed(self.producers(), self.consumer_of(
                ("s1", "capped_pair", "first", []), ("s2", "capped_pair", "second", []),
                binding=binding))

            assert bag.errors == ()


class TestChainingIsRefused:
    """A producer input supplied by another AUTO value of the consumer: one producer feeding
    another. Rejected (`design/impl-layer.md`, "Rejected -- chaining producers"), on either
    route, a value given back to its own producer included."""

    def sizes(self):
        return module("tox_base",
            procedure("base_size", integer("n_items", doc="items"),
                      integer("base", Intent.OUT, doc="a size"),
                      integer("ierr", Intent.OUT, doc="error")),
            procedure("plus_size", integer("other", doc="a size"),
                      integer("plus", Intent.OUT, doc="a size"),
                      integer("ierr", Intent.OUT, doc="error")),
        )

    def consumer_of(self, *sized, binding=False):
        arguments = [integer("n_items", doc="items")]
        for name, producer, output, rows in sized:
            arguments.append(integer(name, doc=table(*rows),
                                     directives=auto(output, producer, "tox_base")))
            arguments.append(real(("res_" if binding else "tmp_") + name, Intent.OUT, f"({name})",
                                  doc="sized"))
        return consumer("chained_impl", *arguments)

    def chained(self, bag):
        return [(e.message, e.note) for e in bag.errors if "a producer computes here" in e.message]

    def test_one_fed_by_another_is_refused_on_either_route_in_either_order(self):
        base = ("s0", "base_size", "base", [])
        plus = ("s1", "plus_size", "plus", [("other", "s0")])
        for binding in (False, True):
            for sized in ((base, plus), (plus, base)):
                _, bag = analysed(self.sizes(), self.consumer_of(*sized, binding=binding))

                [(message, note)] = self.chained(bag)
                assert message == ("'s1' is computed from 'plus_size', whose input 'other' would "
                                   "be 's0', which a producer computes here too")
                assert note == ("producers do not feed each other: size for an upper bound "
                                "instead (the sizing formulas are monotone, so a buffer sized for "
                                "the bound fits), or take the value as an argument of the caller's")

    def test_a_refined_input_is_told_to_be_a_required_intent_inout_argument(self):
        sizes = module("tox_base",
            procedure("base_size", integer("n_items", doc="items"),
                      integer("base", Intent.OUT, doc="a size"),
                      integer("ierr", Intent.OUT, doc="error")),
            procedure("capped_size", integer("cap", Intent.INOUT, doc="a cap it lowers"),
                      integer("sized", Intent.OUT, doc="a size"),
                      integer("ierr", Intent.OUT, doc="error")))
        cons = self.consumer_of(("s0", "base_size", "base", []),
                                ("s1", "capped_size", "sized", [("cap", "s0")]))
        _, bag = analysed(sizes, cons)

        [(_, note)] = self.chained(bag)
        assert note.endswith(" -- a required intent(inout) one, since the producer refines it "
                             "in place")

    def test_one_given_what_it_fills_is_refused(self):
        _, bag = analysed(self.sizes(),
                          self.consumer_of(("s1", "plus_size", "plus", [("other", "s1")])))

        [(message, _)] = self.chained(bag)
        assert message == ("'s1' is computed from 'plus_size', whose input 'other' would be "
                           "'s1', which a producer computes here too")


class TestACycle:
    """A producer that, through its producers, comes back to its consumer is refused.

    Each call would make the next until the stack ran out -- in the generated Fortran, and in
    Python and R through the producers' bindings.
    """

    def sized_by(self, name, producer, output, size_name="k"):
        return procedure(
            name,
            integer("n_items", doc="items"),
            integer(size_name, doc="size", directives=auto(output, producer, "tox_cyc")),
            real("tmp_w", Intent.OUT, f"({size_name})", doc="scratch"),
            integer(f"{name[:-5]}_size", Intent.OUT, doc="a size"),
            meta=IMPL, is_pure=True,
        )

    def cycles(self, bag):
        return [e for e in bag.errors if "comes back to it" in e.message]

    def test_one_sizing_itself(self):
        _, bag = analysed(module("tox_cyc_impl", self.sized_by("self_impl", "self", "self_size")))

        [error] = self.cycles(bag)
        assert error.message == ("'k' of 'self' is computed from a producer that comes back to "
                                 "it: self -> self")

    def test_two_sizing_each_other(self):
        source = module("tox_cyc_impl",
                        self.sized_by("a_impl", "b", "b_size"),
                        self.sized_by("b_impl", "a", "a_size", size_name="j"))
        _, bag = analysed(source)

        [error] = self.cycles(bag)  # one cycle, reported once
        assert error.message == ("'k' of 'a' is computed from a producer that comes back to "
                                 "it: a -> b -> a")

    def test_two_the_bindings_call_are_refused_too(self):
        # hand-written: nothing in Fortran calls either, but Python and R would, forever
        def sizes_result(name, producer, output):
            return procedure(
                name,
                integer("n_items", doc="items"),
                integer("k", doc="size", directives=auto(output, producer, "tox_hand")),
                real("res", Intent.OUT, "(k)", doc="the result"),
                integer(f"{name}_size", Intent.OUT, doc="a size"),
                integer("ierr", Intent.OUT, doc="error"),
            )
        _, bag = analysed(module("tox_hand", sizes_result("p", "q", "q_size"),
                                 sizes_result("q", "p", "p_size")))

        assert [e.message.split(": ", 1)[1] for e in self.cycles(bag)] == ["p -> q -> p"]

    def test_a_chain_that_ends_is_not_a_cycle(self):
        middle = self.sized_by("middle_impl", "end", "end_size")
        end = procedure("end_impl", integer("n_items", doc="items"),
                        integer("end_size", Intent.OUT, doc="a size"), meta=IMPL, is_pure=True)
        cons = module("tox_top_impl", procedure(
            "top_impl", integer("n_items", doc="items"),
            integer("k", doc="size", directives=auto("middle_size", "middle", "tox_cyc")),
            real("tmp_w", Intent.OUT, "(k)", doc="scratch"), meta=IMPL, is_pure=True))
        _, bag = analysed(module("tox_cyc_impl", middle, end), cons)

        assert bag.errors == ()


class TestASameModuleProducer:
    """Rule 5: a producer generated into the consumer's own module is not `use`d."""

    def test_the_module_does_not_use_itself(self):
        same = producers(procedure(
            "local_total_impl",
            integer("n_items", doc="items"),
            real("x", Intent.IN, "(n_items)", doc="values"),
            integer("scale", doc="scale"),
            integer("work_size", doc=table(("factor", "scale")),
                    directives=auto("scaled", "scaled_size", "tox_prod")),
            real("tmp_buf", Intent.OUT, "(work_size)", doc="scratch"),
            real("total", Intent.OUT, doc="the result"),
            meta=IMPL, is_pure=True,
        ))
        synthesis, bag = analysed(same)
        text = fortran(synthesis, "tox_prod")

        assert bag.errors == ()
        assert "use tox_prod," not in text and "use tox_prod\n" not in text
        assert "call scaled_size(&" in body(text, "local_total")


class TestAnOutputNothingTakes:
    """Rule 6: a scalar goes into a local; an array is refused."""

    def test_a_scalar_is_passed_a_local(self):
        synthesis, bag = analysed(producers(), work_consumer("two_sizes", output="first", doc=table()))
        text = body(fortran(synthesis), "work_total")

        assert bag.errors == ()
        assert "integer(int32) :: second_unused" in text
        assert "second = second_unused" in call_of(text, "two_sizes")
        assert "first = work_size" in call_of(text, "two_sizes")

    def test_its_kind_is_imported_where_no_argument_has_it(self):
        # an integer-only consumer drops a logical(c_bool) and a real64 output
        from builders import logical

        odd = producers(procedure(
            "odd_outputs_impl",
            integer("n_items", doc="items"),
            integer("sized", Intent.OUT, doc="the size"),
            logical("flag", Intent.OUT, kind="c_bool", doc="a flag"),
            real("weight", Intent.OUT, doc="a weight"),
            meta=IMPL, is_pure=True,
        ))
        cons = consumer(
            "count_impl",
            integer("n_items", doc="items"),
            integer("work_size", doc="size", directives=auto("sized", "odd_outputs", "tox_prod")),
            integer("tmp_buf", Intent.OUT, "(work_size)", doc="scratch"),
            integer("total", Intent.OUT, doc="the result"),
        )
        synthesis, bag = analysed(odd, cons)
        text = fortran(synthesis)
        header = text.split("contains", 1)[0]

        assert bag.errors == ()
        assert "logical(c_bool) :: flag_unused" in text and "real(real64) :: weight_unused" in text
        assert "use, intrinsic :: iso_c_binding, only: c_bool" in header
        assert "use, intrinsic :: iso_fortran_env, only: int32, real64" in header
        assert "use f42_safeguard" in header  # a C kind: the safeguard comes first

    def test_two_of_one_name_get_a_local_each(self):
        # two_sizes' `second` is an integer, other_two's a real
        other = producers(procedure(
            "other_two_impl",
            integer("n_items", doc="items"),
            integer("first", Intent.OUT, doc="one size"),
            real("second", Intent.OUT, doc="a real"),
            meta=IMPL, is_pure=True,
        ))
        cons = consumer(
            "both_impl",
            integer("n_items", doc="items"),
            integer("a_size", doc="size", directives=auto("first", "two_sizes", "tox_prod")),
            real("tmp_a", Intent.OUT, "(a_size)", doc="scratch"),
            integer("b_size", doc="size", directives=auto("first", "other_two", "tox_prod")),
            real("tmp_b", Intent.OUT, "(b_size)", doc="scratch"),
            real("total", Intent.OUT, doc="the result"),
        )
        synthesis, bag = analysed(other, cons)
        text = body(fortran(synthesis), "both")

        assert bag.errors == ()
        assert "integer(int32) :: second_unused" in text
        assert "real(real64) :: other_two_second_unused" in text
        assert "second = second_unused" in call_of(text, "two_sizes")
        assert "second = other_two_second_unused" in call_of(text, "other_two")

    def test_an_array_is_refused(self):
        _, bag = analysed(producers(), work_consumer("sized_and_rows", output="sized", doc=table()))

        assert any("'rows'" in e.message and "nothing here takes" in e.message for e in bag.errors)

    def test_an_array_the_bindings_receive_anyway_is_not(self):
        # Python and R take every output the producer's binding returns; nothing is dropped
        _, bag = analysed(producers(), result_consumer("sized_and_rows", "sized", doc=table()))

        assert bag.errors == ()


class TestADerivedProducerInput:
    """Rule 7: an input the producer's binding derives resolves like any other.

    Marked derived, so Python and R still leave it out; an error only where a generated
    Fortran wrapper calls the producer, which has to pass it.
    """

    def extent_consumer(self, doc, sized_result=False):
        size = integer("extent_size", doc=doc, directives=auto("extent", "extent_size", "tox_prod"))
        sized = (real("res", Intent.OUT, "(extent_size)", doc="the result") if sized_result
                 else real("tmp_buf", Intent.OUT, "(extent_size)", doc="scratch"))
        return consumer(
            "extent_total_impl",
            integer("m", doc="values"),
            real("x", Intent.IN, "(m)", doc="the values"),
            size, sized,
            real("total", Intent.OUT, doc="the result"),
        )

    def test_a_table_row_supplies_it(self):
        synthesis, bag = analysed(producers(), self.extent_consumer(table(("n", "m"))))
        plan = synthesis.project.procedure("tox_cons", "extent_total_expert") \
            .argument("extent_size").roles.computed_from

        assert bag.errors == ()
        assert [(s.name, s.argument, s.derived) for s in plan.inputs] == [
            ("n", "m", True), ("x", "x", False)]
        assert [s.name for s in plan.binding_inputs] == ["x"]
        assert "n = m" in call_of(body(fortran(synthesis), "extent_total"), "extent_size")

    def test_without_a_row_it_is_an_error_where_fortran_calls_the_producer(self):
        _, bag = analysed(producers(), self.extent_consumer(table()))

        assert any("input(s) n have no argument of the same name" in e.message
                   for e in bag.errors)

    def test_without_a_row_the_bindings_need_nothing(self):
        # the producer's binding reads its `n` off `x`, and nothing in Fortran calls it
        synthesis, bag = analysed(producers(), self.extent_consumer(table(), sized_result=True))
        plan = synthesis.project.procedure("tox_cons", "extent_total") \
            .argument("extent_size").roles.computed_from

        assert bag.errors == ()
        assert [s.name for s in plan.binding_inputs] == ["x"]

    def test_a_hand_written_consumer_needs_nothing_either(self):
        hand = module("tox_hand", procedure(
            "extent_total",
            integer("m", doc="values"),
            real("x", Intent.IN, "(m)", doc="the values"),
            integer("extent_size", doc="size", directives=auto("extent", "extent_size", "tox_prod")),
            real("tmp_buf", Intent.OUT, "(extent_size)", doc="scratch"),
            integer("ierr", Intent.OUT, doc="error"),
            meta=C_BINDING,
        ))
        _, bag = analysed(producers(), hand)

        assert bag.errors == ()


class TestAShapeTheBindingsReadOffAnArray:
    """A producer input that is a shape the consumer's bindings read off its array.

    R spelled it by name, `arr_shape`, which is not a variable there -- the shim reads it off
    `arr` -- so the call ended in "object 'arr_shape' not found". Python has the variable.
    """

    def test_r_reads_it_off_the_array_as_the_shim_does(self):
        from codegen.abi.c_abi import build_project
        from codegen.emit.r_wrapper import RWrapperEmitter

        counter = procedure(
            "count_vectors",
            integer("n_dims", doc="rank"),
            integer("arr_shape", Intent.IN, "(n_dims)", doc="a shape"),
            integer("n_vectors", Intent.OUT, doc="how many vectors"),
            integer("ierr", Intent.OUT, doc="error"),
        )
        norm = procedure(
            "norm",
            integer("n_elements", doc="elements"),
            real("arr", Intent.IN, "(n_elements)", doc="the array"),
            integer("n_dims", doc="rank"),
            integer("arr_shape", Intent.IN, "(n_dims)", doc="its shape"),
            integer("n_vectors", doc="vectors",
                    directives=auto("n_vectors", "count_vectors", "tox_hand")),
            integer("status", Intent.OUT, "(n_vectors)", doc="per vector"),
            integer("ierr", Intent.OUT, doc="error"),
        )
        bag = DiagnosticBag()
        source = project(module("tox_hand", counter, norm))
        analyse_project(source, bag)
        binding = build_project(source, bag)
        assert bag.errors == (), bag.render()
        wrapper = next(w for m in binding for w in m if w.stripped_name == "norm")

        text = RWrapperEmitter().function(wrapper)

        assert "count_vectors(arr_shape = if (is.null(dim(arr))) length(arr) else dim(arr))" in text
