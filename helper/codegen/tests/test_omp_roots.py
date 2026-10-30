"""The `n_omp_threads` convention (#218): an OpenMP root, its wrappers, and what is refused.

A root is an implementation that takes `n_omp_threads` and holds the parallel loop. Its wrappers
carry `M_OMP_PURE` where they would carry `pure`, check that the count is at least 1, document it
with the generator's one text, and -- the allocating one -- size the per-thread work arrays and the
count it hands the implementation with `M_OMP_N_THREADS`, one without OpenMP. The argument and
the two macros are refused everywhere the wrappers could not do that.

Text-level checks on hand-built IR; the fixture `fx_omp_impl` is compiled and run, with and
without OpenMP, by the end-to-end tests.
"""

import re
from pathlib import Path

import pytest

import builders as b
from codegen.abi.c_abi import build_project
from codegen.config import CONVENTIONS, Paths
from codegen.diagnostics import DiagnosticBag, SourceLocation
from codegen.emit.fortran_c import FortranCEmitter
from codegen.emit.fortran_wrapper import FortranWrapperEmitter, WrapperInfo
from codegen.emit.python_ctypes import PythonEmitter
from codegen.emit.r_wrapper import RWrapperEmitter
from codegen.frontend.ford_frontend import FordFrontend
from codegen.frontend.macros import MacroTable, build_directive_patterns
from codegen.ir.directives import (
    DirectiveParser,
    Directives,
    Minimum,
    Prologue,
    RequiredIfMode,
)
from codegen.ir.doc import Doc
from codegen.ir.entities import MacroUse, Meta
from codegen.ir.roles import analyse_project
from codegen.ir.types import Intent
from codegen.ir.validate import validate_project
from codegen.synthesize import _OMP_THREADS_DOC, synthesize_wrappers
from conftest import REPO_ROOT

IMPL_META = Meta(summary="Sum the values", author="AUTHOR")

#: The per-thread work array's declaration, and the root's own lines in the "source"
_STATEMENT_LINE = 10


def omp_threads(**kwargs):
    """The argument as a root declares it: undocumented, the generator documents it."""
    return b.integer("n_omp_threads", kwargs.pop("intent", Intent.IN), **kwargs)


def root(*arguments, name="sums_impl", prefix="M_OMP_PURE", is_pure=True, **kwargs):
    """A root implementation. `arguments` replaces the default signature."""
    if not arguments:
        arguments = (
            b.integer("n", Intent.IN, doc="rows"),
            b.real("values", Intent.IN, "(n)", doc="the data"),
            b.real("tmp_w", Intent.OUT, "(n, n_omp_threads)", doc="scratch per thread"),
            b.real("total", Intent.OUT, doc="their sum"),
            omp_threads(),
        )
    return b.procedure(
        name, *arguments, meta=IMPL_META, is_pure=is_pure, statement_prefix=prefix,
        location=SourceLocation(Path("src/tox/tox_demo_impl.F90"), _STATEMENT_LINE), **kwargs,
    )


def lone_root(threads=None, **kwargs):
    """A root without work arrays: one wrapper, which takes nothing over. `threads` replaces
    its `n_omp_threads`."""
    return root(
        b.integer("n", Intent.IN, doc="rows"),
        b.real("values", Intent.IN, "(n)", doc="the data"),
        b.real("total", Intent.OUT, doc="their sum"),
        threads if threads is not None else omp_threads(),
        **kwargs,
    )


def impl_module(*procedures, **kwargs):
    return b.module("tox_demo_impl", *(procedures or (root(),)), doc="Demo roots.", **kwargs)


def pipeline(*modules):
    """Synthesise, analyse and validate, as a real run does; the synthesis and the diagnostics."""
    bag = DiagnosticBag()
    synthesis = synthesize_wrappers(b.project(*modules), CONVENTIONS, bag)
    analyse_project(synthesis.project, bag)
    validate_project(synthesis.project, bag)
    return synthesis, bag


def errors(*modules):
    _, bag = pipeline(*modules)
    return [d.message for d in bag.errors]


def emitted(*modules, generated="tox_demo"):
    synthesis, bag = pipeline(*modules)
    assert bag.errors == (), bag.render()
    info = {}
    for spec in synthesis.specs:
        for wrapper in (spec.validating, spec.allocating):
            if wrapper is not None:
                info[wrapper.name.lower()] = WrapperInfo(spec.impl.name, spec.mode_fix)
    return FortranWrapperEmitter(project=synthesis.project).module(
        synthesis.project.module(generated), info
    )


def subroutine(text, name):
    """One subroutine of an emitted module, from its statement to its end."""
    start = text.index(f"subroutine {name}(")
    return text[text.rindex("\n", 0, start) + 1: text.index(f"end subroutine {name}\n", start)]


def bound(*modules):
    synthesis, bag = pipeline(*modules)
    binding = build_project(synthesis.project, bag)
    assert bag.errors == (), bag.render()
    return next(module for module in binding if module.stripped_name == "tox_demo")


# -- what the wrappers make of it --------------------------------------------------------


class TestThePrefix:
    def test_both_tiers_of_a_pure_root_carry_the_macro(self):
        text = emitted(impl_module())

        assert "M_OMP_PURE subroutine sums(" in text
        assert "M_OMP_PURE subroutine sums_expert(" in text
        assert "pure subroutine" not in text.replace("M_OMP_PURE subroutine", "")

    def test_a_pure_implementation_without_the_argument_keeps_plain_pure(self):
        text = emitted(impl_module(b.procedure(
            "sums_impl", b.integer("n", Intent.IN, doc="rows"),
            b.real("values", Intent.IN, "(n)", doc="the data"),
            b.real("total", Intent.OUT, doc="their sum"), meta=IMPL_META, is_pure=True,
        )))

        assert "pure subroutine sums(" in text
        assert "M_OMP_PURE" not in text

    def test_an_impure_root_gets_no_prefix(self):
        text = emitted(impl_module(root(prefix="", is_pure=False)))

        assert re.search(r"^    subroutine sums\(", text, re.MULTILINE)
        assert re.search(r"^    subroutine sums_expert\(", text, re.MULTILINE)
        assert "M_OMP_PURE" not in text

    def test_the_c_wrapper_stays_impure(self):
        text = FortranCEmitter().module(bound(impl_module()))

        assert "M_OMP_PURE" not in text
        assert "pure subroutine" not in text


class TestTheImpliedMinimum:
    """At least 1 in both tiers, as an injected `DM_MIN(1_int32)`, inside the validation guard."""

    def test_both_tiers_check_it_at_their_own_position(self):
        text = emitted(impl_module())

        # foo drops tmp_w, which comes before it: 4th there, 5th in the expert tier
        assert ("call validate_in_range_int(n_omp_threads, ierr, arg_pos=4_int32, min=1_int32)"
                in subroutine(text, "sums"))
        assert ("call validate_in_range_int(n_omp_threads, ierr, arg_pos=5_int32, min=1_int32)"
                in subroutine(text, "sums_expert"))

    def test_the_range_replaces_the_extent_check_in_the_expert_tier(self):
        # it is an extent there (of tmp_w), and a documented range wins over the dimension
        # check, which would report 0 as ERR_EMPTY_INPUT rather than ERR_INVALID_INPUT
        assert "validate_dimension_size(n_omp_threads" not in emitted(impl_module())

    def test_it_sits_inside_the_validation_guard(self):
        # FES 2026-10-30 (Q8): like every check of a value the caller knows
        for name in ("sums", "sums_expert"):
            body = subroutine(emitted(impl_module()), name)
            guard = body.index("#ifndef NO_INPUT_VALIDATION")
            check = body.index("validate_in_range_int(n_omp_threads")
            assert guard < check < body.index("#endif")

    def test_the_lone_wrapper_checks_it_too(self):
        assert ("call validate_in_range_int(n_omp_threads, ierr, arg_pos=4_int32, min=1_int32)"
                in emitted(impl_module(lone_root())))


class TestTheThreadCount:
    def test_the_allocating_wrapper_sizes_the_thread_axis_with_the_macro(self):
        body = subroutine(emitted(impl_module()), "sums")

        assert "M_ALLOCATE(tmp_w(n, M_OMP_N_THREADS(n_omp_threads)))" in body

    def test_and_hands_the_implementation_the_same_count(self):
        body = subroutine(emitted(impl_module()), "sums")

        assert "n_omp_threads = M_OMP_N_THREADS(n_omp_threads)" in body

    def test_the_expert_tier_passes_the_caller_s_count_on(self):
        body = subroutine(emitted(impl_module()), "sums_expert")

        assert "n_omp_threads = n_omp_threads" in body
        assert "M_OMP_N_THREADS" not in body

    def test_a_lone_wrapper_passes_it_on(self):
        text = emitted(impl_module(lone_root()))

        assert "subroutine sums_expert(" not in text
        assert "n_omp_threads = n_omp_threads&" in subroutine(text, "sums")
        assert "M_OMP_N_THREADS" not in text
        assert "M_OMP_PURE subroutine sums(" in text

    def test_an_array_whose_last_extent_is_another_keeps_it(self):
        text = emitted(impl_module(root(
            b.integer("n", Intent.IN, doc="rows"),
            b.real("values", Intent.IN, "(n)", doc="the data"),
            b.real("tmp_w", Intent.OUT, "(n, n_omp_threads)", doc="scratch per thread"),
            b.real("tmp_shared", Intent.OUT, "(n)", doc="shared scratch"),
            b.real("total", Intent.OUT, doc="their sum"),
            omp_threads(),
        )))

        assert "M_ALLOCATE(tmp_shared(n))" in text


class TestAModeSplitRoot:
    """`group_centroid` is one: two modes, each with both tiers, four wrappers."""

    def module(self):
        mode_table = [
            "Which statistic.",
            "| Mode | Value | Procedure |",
            "|------|-------|-----------|",
            "| the mean | [[tox_demo_impl(module):MODE_MEAN(variable)]] | sums_mean |",
            "| the median | [[tox_demo_impl(module):MODE_MEDIAN(variable)]] | sums_median |",
        ]
        return impl_module(root(
            b.integer("n", Intent.IN, doc="rows"),
            b.real("values", Intent.IN, "(n)", doc="the data"),
            b.integer("stat_mode", Intent.IN, doc=mode_table),
            b.real("weight", Intent.IN, optional=True, doc="median weight", directives=Directives(
                required_if_mode=RequiredIfMode("stat_mode", "tox_demo_impl", "MODE_MEDIAN"))),
            b.real("tmp_w", Intent.OUT, "(n, n_omp_threads)", doc="scratch per thread"),
            b.real("total", Intent.OUT, doc="their sum"),
            omp_threads(),
        ), parameters=(b.parameter("MODE_MEAN", "1"), b.parameter("MODE_MEDIAN", "2")))

    def test_four_wrappers_each_with_the_prefix(self):
        text = emitted(self.module())

        for name in ("sums_mean", "sums_mean_expert", "sums_median", "sums_median_expert"):
            assert f"M_OMP_PURE subroutine {name}(" in text

    def test_each_allocating_one_sizes_and_passes_with_the_macro(self):
        text = emitted(self.module())

        for name in ("sums_mean", "sums_median"):
            body = subroutine(text, name)
            assert "M_ALLOCATE(tmp_w(n, M_OMP_N_THREADS(n_omp_threads)))" in body
            assert "n_omp_threads = M_OMP_N_THREADS(n_omp_threads)" in body
            assert "validate_in_range_int(n_omp_threads" in body
        for name in ("sums_mean_expert", "sums_median_expert"):
            assert "n_omp_threads = n_omp_threads" in subroutine(text, name)


class TestTheBindings:
    """Required in every language, in the implementation's order, with no default (DECIDED 3)."""

    def test_python_takes_it_after_the_data(self):
        text = PythonEmitter().module(bound(impl_module()))

        assert re.search(r"def sums\(\s*values,\s*n_omp_threads,?\s*\)", text)

    def test_r_takes_it_after_the_data(self):
        text = RWrapperEmitter().module(bound(impl_module()))

        assert re.search(r"sums <- function\(values, n_omp_threads\)", text)

    def test_neither_sees_the_work_array_or_an_expert_tier(self):
        synthesis, _ = pipeline(impl_module())

        assert "sums_expert" in synthesis.expert_only_in_fortran


class TestTheDocumentation:
    """One text, the generator's, reaching every language with the minimum it checks."""

    def test_it_reaches_the_fortran_wrappers(self):
        text = emitted(impl_module())

        for line in _OMP_THREADS_DOC:
            assert text.count(f"!! {line}") == 2  # foo and foo_expert

    def test_it_reaches_c_python_and_r(self):
        module = bound(impl_module())
        texts = {
            "C": FortranCEmitter().module(module),
            "Python": PythonEmitter().module(module),
            "R": RWrapperEmitter().module(module),
        }
        prose, minimum = _OMP_THREADS_DOC[:-1], _OMP_THREADS_DOC[-1]
        for language, text in texts.items():
            for line in prose:
                assert line in text, (language, line)
        # the minimum as every DM_MIN reaches each language: Python and R show the value
        # without the Fortran kind (`emit.doc_literals`)
        assert minimum in texts["C"]
        for language in ("Python", "R"):
            assert "The minimum valid value is `1`." in texts[language], language

    def test_its_minimum_is_what_dm_min_writes(self):
        # the last line stands for `DM_MIN(1_int32)`: the header's own expansion, which the
        # directive parser reads back as the minimum the wrappers check
        macros = MacroTable(REPO_ROOT / "src/macros.h", include_paths=(REPO_ROOT,))
        assert _OMP_THREADS_DOC[-1] == macros.expand("DM_MIN(1_int32)").strip()

        parser = DirectiveParser(build_directive_patterns(macros))
        directives = parser.parse(Doc.parse(list(_OMP_THREADS_DOC)))
        assert directives == Directives(minimum=Minimum("1_int32"))

    def test_the_wrapper_carries_that_minimum(self):
        synthesis, _ = pipeline(impl_module())
        argument = synthesis.project.procedure("tox_demo", "sums").argument("n_omp_threads")

        assert argument.directives == Directives(minimum=Minimum("1_int32"))
        assert argument.doc.text.split("\n")[-1] == _OMP_THREADS_DOC[-1]


class TestTheRole:
    def test_it_is_marked_on_the_wrapper_argument(self):
        synthesis, _ = pipeline(impl_module())
        procedure = synthesis.project.procedure("tox_demo", "sums")

        assert procedure.argument("n_omp_threads").roles.is_omp_threads
        assert not procedure.argument("n").roles.is_omp_threads


# -- what is refused ------------------------------------------------------------------


class TestAcceptedRoots:
    @pytest.mark.parametrize("prefix", ["M_OMP_PURE", "", "impure", "M_OMP_PURE elemental",
                                        "impure elemental", "recursive"])
    def test_these_statements_are_fine(self, prefix):
        assert errors(impl_module(root(prefix=prefix))) == []

    def test_a_root_using_both_macros_where_they_belong(self):
        uses = (
            MacroUse("M_OMP_PURE", SourceLocation(None, _STATEMENT_LINE), "sums_impl"),
            MacroUse("M_OMP_PARALLEL", SourceLocation(None, _STATEMENT_LINE + 9), "sums_impl"),
        )
        assert errors(impl_module(macro_uses=uses)) == []


class TestWhereItMayAppear:
    """R1: only on an implementation and the wrappers generated from it."""

    def test_an_export_may_not_take_it(self):
        export = b.module("tox_plain", b.procedure(
            "hand_written", b.integer("n", Intent.IN, doc="rows"),
            b.integer("n_omp_threads", Intent.IN, doc="threads"), b.ierr(),
        ), doc="Hand-written.")

        assert errors(export) == [
            "argument 'n_omp_threads' of 'hand_written', which is not an implementation"
        ]

    def test_nor_return_it(self):
        export = b.module("tox_plain", b.procedure(
            "threads_query", b.ierr(),
            result=b.integer("n_omp_threads", Intent.OUT, doc="threads", is_result=True),
        ), doc="Hand-written.")

        assert errors(export) == [
            "result 'n_omp_threads' of 'threads_query', which is not an implementation"
        ]

    def test_a_recommend_routine_beside_a_root_may_not_take_it(self):
        recommend = b.procedure(
            "work_size", b.integer("n_omp_threads", Intent.IN, doc="threads"),
            b.integer("size", Intent.OUT, doc="size"),
        )

        assert errors(impl_module(root(), recommend)) == [
            "argument 'n_omp_threads' of 'work_size', which is not an implementation"
        ]

    def test_a_prologue_may_not_take_it(self):
        prologue = b.procedure(
            "guard",
            b.integer("n", Intent.IN, doc="rows"),
            omp_threads(doc="threads"),
            b.logical("handled", Intent.OUT, doc="whether the call was dealt with"),
            meta=Meta(summary="Guard", author="AUTHOR"),
        )
        guarded = root(directives=Directives(prologue=Prologue("guard", "tox_demo_impl")))

        assert errors(impl_module(prologue, guarded)) == [
            "prologue 'guard' takes 'n_omp_threads'"
        ]

    def test_an_implementation_function_may_not_return_it(self):
        function = b.procedure(
            "count_impl", b.integer("n", Intent.IN, doc="rows"),
            result=b.integer("n_omp_threads", Intent.OUT, doc="threads", is_result=True),
            meta=IMPL_META,
        )

        assert "result of 'count_impl' is named 'n_omp_threads'" in errors(impl_module(function))


class TestItsDeclaration:
    """R2: a required, scalar `integer(int32)`, `intent(in)`."""

    @pytest.mark.parametrize("argument, problem", [
        (omp_threads(optional=True), "optional"),
        (omp_threads(dimension="(2)"), "an array of rank 1"),
        (omp_threads(kind="int64"), "integer(int64)"),
        (omp_threads(intent=Intent.INOUT), "intent(inout)"),
    ])
    def test_anything_else_is_refused(self, argument, problem):
        assert f"'n_omp_threads' is {problem}" in errors(impl_module(lone_root(argument)))


class TestWhatItMaySize:
    """R3: only the last extent of a `tmp_` array, as the bare name."""

    def refused(self, *extra):
        return errors(impl_module(root(
            b.integer("n", Intent.IN, doc="rows"),
            b.real("values", Intent.IN, "(n)", doc="the data"),
            b.real("total", Intent.OUT, doc="their sum"),
            omp_threads(),
            *extra,
        )))

    def test_a_non_last_extent(self):
        assert self.refused(b.real("tmp_w", Intent.OUT, "(n_omp_threads, n)", doc="scratch")) == [
            "'n_omp_threads' is used in the extent 'n_omp_threads' of 'tmp_w'"
        ]

    def test_an_array_that_is_not_a_work_array(self):
        assert self.refused(b.real("per_thread", Intent.OUT, "(n, n_omp_threads)", doc="out")) == [
            "'n_omp_threads' is used in the extent 'n_omp_threads' of 'per_thread'"
        ]

    @pytest.mark.parametrize("name, extents", [
        ("tmp_w", "(n, 2*n_omp_threads)"),
        ("out", "(n_omp_threads + 1)"),
    ])
    def test_inside_an_extent_expression(self, name, extents):
        found = self.refused(b.real(name, Intent.OUT, extents, doc="sized by it"))

        assert len(found) == 1 and found[0].startswith(
            f"'n_omp_threads' is used in the extent "
        ) and found[0].endswith(f"of '{name}'")

    def test_a_character_length(self):
        found = self.refused(b.character("label", Intent.OUT, length="n_omp_threads", doc="a"))

        assert found == ["'n_omp_threads' is used in the length 'n_omp_threads' of 'label'"]


class TestItsDocumentation:
    """R4: the generator writes it; an author's text or directive is refused, not dropped."""

    def test_its_own_text(self):
        source = impl_module(lone_root(omp_threads(doc="how many threads")))

        assert errors(source) == ["'n_omp_threads' has documentation of its own"]

    def test_a_directive(self):
        threads = omp_threads(doc="The minimum valid value is `2_int32`.",
                              directives=Directives(minimum=Minimum("2_int32")))

        assert errors(impl_module(lone_root(threads))) == ["'n_omp_threads' carries DM_MIN"]


class TestItsPurity:
    """R5: a root is `M_OMP_PURE` or impure, read from the statement as written."""

    def test_plain_pure(self):
        assert errors(impl_module(root(prefix="pure"))) == [
            "OpenMP root 'sums_impl' is declared pure"
        ]

    def test_elemental_without_impure(self):
        assert errors(impl_module(root(prefix="elemental"))) == [
            "OpenMP root 'sums_impl' is elemental, which is pure"
        ]

    def test_a_statement_that_could_not_be_read(self):
        assert errors(impl_module(root(prefix=None))) == [
            "the statement of 'sums_impl' could not be read"
        ]


class TestWhereTheMacrosMayBeWritten:
    """R5: `M_OMP_PURE` on a root's statement, `!M_OMP_PARALLEL(` inside a root -- nowhere else."""

    def plain_impl(self):
        return b.procedure(
            "scale_impl", b.integer("n", Intent.IN, doc="rows"),
            b.real("values", Intent.INOUT, "(n)", doc="the data"),
            meta=IMPL_META, is_pure=True, statement_prefix="M_OMP_PURE",
            location=SourceLocation(None, 30),
        )

    def test_m_omp_pure_on_an_implementation_without_the_argument(self):
        uses = (MacroUse("M_OMP_PURE", SourceLocation(None, 30), "scale_impl"),)

        assert errors(impl_module(root(), self.plain_impl(), macro_uses=uses)) == [
            "M_OMP_PURE is written in 'scale_impl', and only an OpenMP root may carry it"
        ]

    def test_a_parallel_loop_in_an_implementation_without_the_argument(self):
        uses = (MacroUse("M_OMP_PARALLEL", SourceLocation(None, 35), "scale_impl"),)

        assert errors(impl_module(root(), self.plain_impl(), macro_uses=uses)) == [
            "!M_OMP_PARALLEL(...) is written in 'scale_impl', and only an OpenMP root may carry it"
        ]

    def test_m_omp_pure_in_the_body_of_a_root(self):
        uses = (MacroUse("M_OMP_PURE", SourceLocation(None, _STATEMENT_LINE + 3), "sums_impl"),)

        assert errors(impl_module(macro_uses=uses)) == [
            "M_OMP_PURE is written in the body of 'sums_impl', and only an OpenMP root may carry it"
        ]

    def test_outside_every_procedure(self):
        uses = (MacroUse("M_OMP_PARALLEL", SourceLocation(None, 3), None),)

        assert errors(impl_module(macro_uses=uses)) == [
            "!M_OMP_PARALLEL(...) is written outside every procedure, and only an OpenMP root "
            "may carry it"
        ]

    @pytest.mark.parametrize("macro, written", [
        ("M_OMP_PURE", "M_OMP_PURE"), ("M_OMP_PARALLEL", "!M_OMP_PARALLEL(...)"),
    ])
    def test_in_a_module_that_is_not_an_implementation_module(self, macro, written):
        export = b.module("tox_plain", b.procedure(
            "hand_written", b.integer("n", Intent.IN, doc="rows"), b.ierr(),
            location=SourceLocation(None, 5),
        ), doc="Hand-written.", macro_uses=(MacroUse(macro, SourceLocation(None, 5),
                                                     "hand_written"),))

        assert errors(export) == [
            f"{written} is written in 'hand_written', and only an OpenMP root may carry it"
        ]


# -- the frontend reads both from the source the author edits ---------------------------


@pytest.fixture(scope="module")
def fixture_impls():
    bag = DiagnosticBag()
    parsed = FordFrontend(
        Paths(root=REPO_ROOT, src_dir=Path("helper/codegen/tests/fixtures/impl")), bag
    ).parse()
    assert bag.render() == ""
    return parsed.project


class TestTheFrontend:
    def test_it_reads_the_macro_before_the_preprocessor_replaces_it(self, fixture_impls):
        root_ = fixture_impls.procedure("fx_omp_impl", "fx_column_sums_impl")

        assert root_.statement_prefix == "M_OMP_PURE"
        assert root_.is_pure  # what Ford reads, after the preprocessor

    def test_it_finds_each_macro_in_the_procedure_that_holds_it(self, fixture_impls):
        uses = fixture_impls.module("fx_omp_impl").macro_uses

        assert [(use.macro, use.procedure) for use in uses] == [
            ("M_OMP_PURE", "fx_column_sums_impl"), ("M_OMP_PARALLEL", "fx_column_sums_impl"),
            ("M_OMP_PURE", "fx_column_maxima_impl"), ("M_OMP_PARALLEL", "fx_column_maxima_impl"),
        ]
        statement = fixture_impls.procedure("fx_omp_impl", "fx_column_sums_impl").location
        assert uses[0].location == statement
