"""The Python and R bindings must ask their callers for the same things.

Both emitters decide, independently, which arguments a caller supplies and which the
wrapper works out. They share the `roles` model but not the code, so a change to a role
lands in whichever emitter it was written for and leaves the other on the old rule. That
has happened repeatedly, and the symptoms are bad: an argument both asked for *and*
derived, silently overwriting what the caller passed; or a whole target that no longer
compiles.

Nothing here checks that a particular decision is right -- the other suites do that. It
checks only that the two targets agree, which is the property no single-target test can
see.
"""

from pathlib import Path

import pytest

from codegen.abi.c_abi import build_project
from codegen.config import CONVENTIONS, Paths
from codegen.diagnostics import DiagnosticBag
from codegen.emit.doc_literals import is_ford_block_tag
from codegen.emit.python_ctypes import PythonEmitter
from codegen.emit.c_call import CCallEmitter
from codegen.emit.r_wrapper import RWrapperEmitter
from codegen.frontend.ford_frontend import FordFrontend
from codegen.ir.directives import Directives, ResultSizeIs
from codegen.ir.doc import DocTable
from codegen.ir.entities import Project
from codegen.ir.roles import analyse_project
from codegen.ir.types import Intent
from codegen.ir.validate import validate_project
from codegen.synthesize import synthesize_wrappers
from conftest import REPO_ROOT

import builders as b

FIXTURE_SRC = Path("helper/codegen/tests/fixtures/src")


def _build(src: Path):
    """Parse, synthesise, then run the pipeline -- the order `generate.py` uses.

    Synthesis is not optional here, and skipping it was quietly halving what these tests
    covered: without it the project holds only the *hand-written* exports, so every generated
    wrapper -- which is most of the published API -- went unchecked for the very disagreement
    this file exists to catch. It also made `validate` see a project with no generated
    modules in it, so a documentation link into one resolved to nothing.
    """
    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=src), bag).parse()
    project = synthesize_wrappers(parsed.project, CONVENTIONS, bag).project
    analyse_project(project, bag)
    validate_project(project, bag)
    built = build_project(project, bag)
    assert bag.errors == (), bag.render()
    return built


@pytest.fixture(scope="module")
def binding():
    return _build(FIXTURE_SRC)


@pytest.fixture(scope="module")
def project_binding():
    """The real sources. The fixtures are deliberately small, and every divergence found
    so far was in `src/` -- the serde helpers, the outlier work arrays -- so checking the
    fixtures alone would have caught none of them."""
    return _build(Path("src"))


def _wrappers(binding):
    for module in binding:
        for wrapper in module:
            yield wrapper


def _names(arguments):
    return [argument.name for argument in arguments]


class TestTheTargetsAgree:
    def test_every_wrapper_asks_for_the_same_arguments(self, binding):
        python, r = PythonEmitter(), RWrapperEmitter()

        divergent = {
            wrapper.stripped_name: (_names(python._inputs(wrapper)), _names(r._inputs(wrapper)))
            for wrapper in _wrappers(binding)
            if _names(python._inputs(wrapper)) != _names(r._inputs(wrapper))
        }

        assert not divergent, (
            "Python and R disagree about what the caller supplies:\n"
            + "\n".join(f"  {name}:\n    python: {py}\n    r:      {r_}"
                        for name, (py, r_) in sorted(divergent.items()))
        )


    def test_r_calls_c_with_the_arguments_c_declares(self, binding):
        """The R wrapper and the C++ function it calls must agree, argument for argument.

        They are decided in two places -- `r_wrapper._call_inputs` and `c_call._inputs` --
        so a role change that reaches one and not the other produces an R function that
        calls `.Call("name_call", ...)` with the wrong arguments. C is compiled separately, so
        nothing catches it until the call happens.
        """
        r, c = RWrapperEmitter(), CCallEmitter()

        divergent = {
            wrapper.stripped_name: (_names(r._call_inputs(wrapper)), _names(c._inputs(wrapper)))
            for wrapper in _wrappers(binding)
            if _names(r._call_inputs(wrapper)) != _names(c._inputs(wrapper))
        }

        assert not divergent, (
            "the R wrapper and the C function disagree:\n"
            + "\n".join(f"  {name}:\n    r passes:    {passed}\n    c expects:   {declared}"
                        for name, (passed, declared) in sorted(divergent.items()))
        )

    def test_no_argument_is_both_asked_for_and_derived(self, binding):
        """Deriving an argument the signature also asks for overwrites what was passed.

        It is the failure mode that reads as a mystery: the caller supplies a value, the
        wrapper quietly replaces it with one computed from something that may not even
        exist yet.
        """
        emitter = PythonEmitter()

        overwritten = {}
        for wrapper in _wrappers(binding):
            asked = set(_names(emitter._inputs(wrapper)))
            derived = {
                argument.name for argument in wrapper
                if argument.name in asked
                and emitter._extent_expression(argument, wrapper)
            }
            if derived:
                overwritten[wrapper.stripped_name] = sorted(derived)

        assert not overwritten, f"asked for and then overwritten: {overwritten}"


class TestTheTargetsAgreeOnTheRealProject:
    def test_the_generated_wrappers_are_in_scope(self, project_binding):
        """Assert the premise, or the rest of this class silently stops meaning much.

        Skipping synthesis left 31 of the 132 wrappers here -- the hand-written exports only --
        and every test below still passed, over a quarter of the API.
        """
        published = {w.stripped_name for w in _wrappers(project_binding)}

        assert len(published) > 100
        assert "loess_fit_plain" in published, "a generated wrapper"
        assert "loess_fit_plain_expert" in published, "a generated expert tier"
        assert "get_array_metadata" in published, "a hand-written export"

    def test_every_wrapper_asks_for_the_same_arguments(self, project_binding):
        python, r = PythonEmitter(), RWrapperEmitter()

        divergent = {
            wrapper.stripped_name: (_names(python._inputs(wrapper)), _names(r._inputs(wrapper)))
            for wrapper in _wrappers(project_binding)
            if _names(python._inputs(wrapper)) != _names(r._inputs(wrapper))
        }

        assert not divergent, (
            "Python and R disagree about what the caller supplies:\n"
            + "\n".join(f"  {name}:\n    python: {py}\n    r:      {r_}"
                        for name, (py, r_) in sorted(divergent.items()))
        )


    def test_r_calls_c_with_the_arguments_c_declares(self, project_binding):
        """The R wrapper and the C++ function it calls must agree, argument for argument.

        They are decided in two places -- `r_wrapper._call_inputs` and `c_call._inputs` --
        so a role change that reaches one and not the other produces an R function that
        calls `.Call("name_call", ...)` with the wrong arguments. C is compiled separately, so
        nothing catches it until the call happens.
        """
        r, c = RWrapperEmitter(), CCallEmitter()

        divergent = {
            wrapper.stripped_name: (_names(r._call_inputs(wrapper)), _names(c._inputs(wrapper)))
            for wrapper in _wrappers(project_binding)
            if _names(r._call_inputs(wrapper)) != _names(c._inputs(wrapper))
        }

        assert not divergent, (
            "the R wrapper and the C function disagree:\n"
            + "\n".join(f"  {name}:\n    r passes:    {passed}\n    c expects:   {declared}"
                        for name, (passed, declared) in sorted(divergent.items()))
        )

    def test_no_argument_is_both_asked_for_and_derived(self, project_binding):
        emitter = PythonEmitter()

        overwritten = {}
        for wrapper in _wrappers(project_binding):
            asked = set(_names(emitter._inputs(wrapper)))
            derived = {
                argument.name for argument in wrapper
                if argument.name in asked
                and emitter._extent_expression(argument, wrapper)
            }
            if derived:
                overwritten[wrapper.stripped_name] = sorted(derived)

        assert not overwritten, f"asked for and then overwritten: {overwritten}"


class TestTierNotes:
    """What each half of a published pair says about the other."""

    def note(self, is_expert, prepares=("`values_perm`",), how=("seeds `values_perm`",)):
        from codegen.emit.doc_tiers import TierNote

        return TierNote("crunch_expert" if not is_expert else "crunch", prepares, how,
                        is_expert=is_expert)

    def test_the_plain_half_says_what_it_does_and_offers_the_other(self):
        lines = self.note(is_expert=False).lines("`{}`")

        assert lines[0] == "This entry point seeds `values_perm`."
        assert lines[1] == "Call `crunch_expert` to do that yourself."

    def test_the_expert_half_says_what_it_leaves_to_you(self):
        lines = self.note(is_expert=True).lines("`{}`")

        assert lines[0] == "The expert entry point: you supply `values_perm` yourself."
        assert lines[1] == "`crunch` seeds `values_perm`."

    def test_a_prologue_alone_changes_no_argument_so_it_reads_differently(self):
        lines = self.note(is_expert=True, prepares=(),
                          how=("runs `guard` first, which may answer the call outright",)
                          ).lines("`{}`")

        assert lines == [
            "The expert entry point: `crunch` runs `guard` first, which may answer the "
            "call outright; this one does not."
        ]

    def test_r_marks_identifiers_its_own_way(self):
        # the parts are stored with backticks; each language renders them as it marks code,
        # and the rest of the generated roxygen uses \code{} for text the emitter wrote
        lines = self.note(is_expert=False).lines("\\code{{{}}}")

        assert lines[0] == "This entry point seeds \\code{values_perm}."
        assert lines[1] == "Call \\code{crunch_expert} to do that yourself."

    def test_only_pairs_that_both_reach_the_language_get_a_note(self, real_binding=None):
        # over the real tree: a note exists for exactly the published _expert functions
        from codegen.config import CONVENTIONS
        from codegen.diagnostics import DiagnosticBag
        from codegen.emit.doc_tiers import build

        binding, synthesis = _real_published()
        notes = build(binding, synthesis.specs, CONVENTIONS)
        published = {w.stripped_name for m in binding for w in m}
        experts = {name for name in published if name.endswith(CONVENTIONS.expert_suffix)}
        assert experts, "no expert tier survived, so this asserts nothing"
        for expert in experts:
            assert expert in notes, expert
            assert expert[: -len(CONVENTIONS.expert_suffix)] in notes, expert
        assert set(notes) == experts | {e[: -len(CONVENTIONS.expert_suffix)] for e in experts}


def _real_published():
    """The real tree's binding as Python and R publish it, plus its synthesis."""
    from pathlib import Path

    from codegen.abi.c_abi import build_project
    from codegen.config import CONVENTIONS, Paths
    from codegen.diagnostics import DiagnosticBag
    from codegen.frontend.ford_frontend import FordFrontend
    from codegen.generate import _published_to_the_languages
    from codegen.ir.roles import analyse_project
    from codegen.synthesize import synthesize_wrappers

    from conftest import REPO_ROOT

    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=Path("src")), bag).parse()
    synthesis = synthesize_wrappers(parsed.project, CONVENTIONS, bag)
    analyse_project(synthesis.project, bag)
    binding = build_project(synthesis.project, bag, CONVENTIONS)
    return _published_to_the_languages(binding, synthesis), synthesis


class TestTheTargetsDocumentTheSame:
    """Both targets render the same Fortran doc, and for a long time only one of them
    rendered all of it: roxygen wrote the summary *and* the body, numpydoc wrote the summary
    and went straight to Parameters. So every generated Python function was missing whatever
    the author wrote below `summary:` -- which is usually where the contract is, and there is
    no other place a Python caller could have read it.

    Checked on the real sources rather than the fixtures, for the reason `project_binding`
    gives: the fixtures are small and this kind of loss only shows on real prose.
    """

    def _body_lines(self, wrapper):
        """The prose blocks of a procedure's doc, tables and Ford tags dropped -- what both
        renderers are supposed to carry over."""
        return [
            block.text.strip()
            for block in wrapper.doc
            if not isinstance(block, DocTable)
            and not is_ford_block_tag(block.text)
            and block.text.strip()
        ]

    def test_the_doc_body_reaches_python_wherever_it_reaches_r(self, project_binding):
        python, r = PythonEmitter(), RWrapperEmitter()

        missing = []
        for wrapper in _wrappers(project_binding):
            body = self._body_lines(wrapper)
            if not body:
                continue
            py = python.function(wrapper)
            roxygen = r.function(wrapper)
            for line in body:
                # a sentence long enough to be prose, not a fragment reformatted per target
                probe = line.split(".")[0][:60]
                if len(probe) < 25 or probe not in roxygen:
                    continue
                if probe not in py:
                    missing.append((wrapper.stripped_name, probe))

        assert not missing, "documented in R and not in Python:\n" + "\n".join(
            f"  {name}: {probe!r}" for name, probe in missing[:10]
        )

    def test_the_scale_direction_is_not_lost(self, project_binding):
        """The worked example from the report. `compute_tissue_versatility` normalises to
        [0, 1] and only the doc body says which end is which, so without it the returned
        number cannot be read at all."""
        wrapper = next(
            w for w in _wrappers(project_binding)
            if w.stripped_name == "compute_tissue_versatility"
        )

        text = PythonEmitter().function(wrapper)

        assert "0 means uniform expression" in text
        # and it lands where numpydoc puts an extended summary: after the summary line,
        # before Parameters
        assert text.index("0 means uniform expression") < text.index("Parameters")


class TestConstantSizesInBothBindings:
    """Neither binding can name a Fortran module constant, so both use its value (#221 M2)."""

    def module(self):
        import builders as b
        from codegen.abi.c_abi import build_module
        from codegen.ir.roles import analyse
        from codegen.ir.types import Intent

        module = b.module(
            "m",
            b.procedure(
                "p",
                b.integer("n", Intent.IN, doc="points"),
                b.real("points", Intent.IN, "(N_AXES, n)", doc="the points"),
                b.character("name", Intent.IN, length="MAX_LEN", doc="a name"),
                b.real("sums", Intent.OUT, "(N_AXES)", doc="per axis"),
                b.ierr(),
            ),
            parameters=(b.parameter("N_AXES", "2"), b.parameter("MAX_LEN", "8")),
        )
        b.project(module)
        bag = DiagnosticBag()
        analyse(module.procedure("p"), bag)
        return build_module(module, bag)

    def test_python_allocates_and_pads_with_the_values(self):
        import ast

        text = PythonEmitter().module(self.module())

        assert "np.empty((2,)" in text
        assert ".ljust(8)" in text
        # the code reads no such name; only what a reader is shown keeps it
        read = {node.id for node in ast.walk(ast.parse(text)) if isinstance(node, ast.Name)}
        assert not read & {"N_AXES", "MAX_LEN"}

    def test_python_checks_the_constant_axis_against_its_value(self):
        text = PythonEmitter().module(self.module())

        assert "if points.shape[0] != 2:" in text

    def test_the_documentation_shows_the_constant_by_name(self):
        # what the author wrote, not the value it stands for (#221 WP1 review N5)
        text = PythonEmitter().module(self.module())

        assert "points : np.ndarray[np.float64] of shape (N_AXES, n,)" in text

    def test_the_r_wrapper_names_no_constant(self):
        import re

        text = RWrapperEmitter().module(self.module())
        code = re.sub(r'"[^"]*"', '""', "\n".join(
            line for line in text.splitlines() if not line.lstrip().startswith("#")))

        assert "N_AXES" not in code and "MAX_LEN" not in code
        assert "if (dim(points)[1] != 2L)" in code


class TestTheTargetsTrimTheSameAxis:
    """#221 B8: `DM_RESULT_SIZE_IS` counts along the LAST axis, at any rank, in both languages.

    R took `utils::head(x, n)`, which on a matrix keeps the first n rows -- the padding
    columns included -- while Python's `[..., :n]` kept the first n columns.
    """

    @pytest.mark.parametrize("rank", [1, 2, 3])
    def test_the_count_trims_the_last_axis(self, rank):
        extents = ", ".join(["2"] * (rank - 1) + ["n"])
        found = b.real("found", Intent.OUT, f"({extents})",
                       directives=Directives(result_size_is=ResultSizeIs("n_found")))
        procedure = b.procedure("p", b.integer("n"), found, b.integer("n_found", Intent.OUT),
                                b.ierr())
        project = Project([b.module("m", procedure)])
        bag = DiagnosticBag()
        analyse_project(project, bag)
        (wrapper,) = _wrappers(build_project(project, bag))
        result = next(a for a in wrapper if a.name == "found")

        python = PythonEmitter()._result_expression(result, wrapper)
        r = RWrapperEmitter()._result_expression(result, wrapper)

        assert python == "found[..., :n_found.value]"
        if rank == 1:
            assert r == "utils::head(.result$found, .result$n_found)"
        else:
            leading = ", " * (rank - 1)
            assert r == f".result$found[{leading}seq_len(.result$n_found), drop = FALSE]"


FIXTURE_IMPL = Path("helper/codegen/tests/fixtures/impl")


@pytest.fixture(scope="module")
def impl_published():
    """The implementation fixtures' binding as Python and R publish it, the empty experts gone.

    These are where the producers that fill an array live, and with them every AUTO value a
    generated expert computes in one tier and takes as a parameter in another (E2).
    """
    from codegen.generate import _published_to_the_languages

    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=FIXTURE_IMPL), bag).parse()
    synthesis = synthesize_wrappers(parsed.project, CONVENTIONS, bag)
    analyse_project(synthesis.project, bag)
    validate_project(synthesis.project, bag)
    built = build_project(synthesis.project, bag)
    assert bag.errors == (), bag.render()
    return _published_to_the_languages(built, synthesis)


def _call_arity(r_function: str) -> tuple[str, int]:
    """The `.Call` an R function makes: the entry point's name, and how many arguments it passes."""
    import re

    match = re.search(r'\.Call\("(\w+)"((?:, [\w.]+)*)\)', r_function)
    assert match, r_function
    return match.group(1), len([a for a in match.group(2).split(", ") if a])


class TestEveryPublishedBindingAgrees:
    """For every published binding of every fixture: what R passes is what the shim registers,
    nothing a caller passes is computed or derived over, and nothing is both asked and derived.

    The E2 attack's prototype broke every one of these with the sites it did not route: R took
    `norms` and then overwrote it with the producer's, and called a `.Call` registered with two
    arguments with three. `test_parameters_reach_the_call.py` runs them, a sentinel for every
    parameter.
    """

    @pytest.fixture(params=["src", "impl"])
    def published(self, request, binding, impl_published):
        return binding if request.param == "src" else impl_published

    def test_r_passes_what_the_shim_registers(self, published):
        import re

        r, c = RWrapperEmitter(), CCallEmitter()
        registered = dict(re.findall(r'\{"(\w+)", \(DL_FUNC\) &\w+, (\d+)\}',
                                     c.registration(list(published))))
        divergent = {}
        for wrapper in _wrappers(published):
            name, passed = _call_arity(r.function(wrapper.with_constant_sizes_evaluated()))
            if int(registered[name]) != passed:
                divergent[name] = (passed, registered[name])

        assert not divergent, f"R passes, the shim registers: {divergent}"

    def test_no_parameter_is_computed_or_derived_over(self, published):
        from codegen.emit.producer_calls import auto_plans

        python, r = PythonEmitter(), RWrapperEmitter()
        overwritten = {}
        for wrapper in _wrappers(published):
            computed = {argument.name for argument, _ in auto_plans(wrapper)}
            derived = {a.name for a in wrapper if python._extent_expression(a, wrapper)}
            for language, emitter in (("python", python), ("r", r)):
                asked = set(_names(emitter._inputs(wrapper)))
                if asked & (computed | derived):
                    overwritten[(language, wrapper.stripped_name)] = sorted(
                        asked & (computed | derived))

        assert not overwritten, f"asked for and then overwritten: {overwritten}"

    def test_python_and_r_ask_for_the_same_arguments(self, published):
        python, r = PythonEmitter(), RWrapperEmitter()

        divergent = {
            wrapper.stripped_name: (_names(python._inputs(wrapper)), _names(r._inputs(wrapper)))
            for wrapper in _wrappers(published)
            if _names(python._inputs(wrapper)) != _names(r._inputs(wrapper))
        }

        assert not divergent, divergent

    def test_the_premise_a_published_expert_takes_an_array_a_producer_fills(self, impl_published):
        python = PythonEmitter()
        expert = next(w for w in _wrappers(impl_published)
                      if w.stripped_name == "fx_norm_rank_expert")

        assert "norms" in _names(python._inputs(expert))


class TestTierNotesListWhatTheLanguagePublishes:
    """R5: the Python and R note names the parameters that language's expert takes and its plain
    tier does not -- not the C arguments, whose work arrays and sizes both tiers of a binding
    prepare themselves -- and says a value a producer computes before the permutation sorted by
    it."""

    def notes(self, impl_published):
        from codegen.emit.doc_tiers import build
        from codegen.frontend.ford_frontend import FordFrontend

        bag = DiagnosticBag()
        parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=FIXTURE_IMPL), bag).parse()
        specs = synthesize_wrappers(parsed.project, CONVENTIONS, bag).specs
        return build(impl_published, specs, CONVENTIONS)

    def test_a_computed_value_comes_before_the_permutation(self, impl_published):
        lines = self.notes(impl_published)["fx_norm_rank_expert"].lines("`{}`")

        assert lines == [
            "The expert entry point: you supply `norms` and `norms_perm` yourself.",
            "`fx_norm_rank` computes `norms` with `fx_column_norms` and seeds `norms_perm` and "
            "sorts it by `norms`.",
        ]

    def test_work_arrays_and_sizes_are_not_listed(self, impl_published):
        notes = self.notes(impl_published)

        assert notes["fx_pos_sum_expert"].lines("`{}`")[0] == \
            "The expert entry point: you supply `values_perm` yourself."
        assert notes["fx_clip_expert"].lines("`{}`") == [
            "The expert entry point: you supply `level` yourself.",
            "`fx_clip` computes `level` for you and runs `fx_clip_guard` first, which may answer "
            "the call outright.",
        ]
