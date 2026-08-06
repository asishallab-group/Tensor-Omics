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
from codegen.config import Paths
from codegen.diagnostics import DiagnosticBag
from codegen.emit.python_ctypes import PythonEmitter
from codegen.emit.c_call import CCallEmitter
from codegen.emit.r_wrapper import RWrapperEmitter
from codegen.frontend.ford_frontend import FordFrontend
from codegen.ir.roles import analyse_project
from codegen.ir.validate import validate_project
from conftest import REPO_ROOT

FIXTURE_SRC = Path("helper/codegen/tests/fixtures/src")


def _build(src: Path):
    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=src), bag).parse()
    analyse_project(parsed.project, bag)
    validate_project(parsed.project, bag)
    built = build_project(parsed.project, bag)
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
        experts = {name for name in published if name.endswith(CONVENTIONS.expert_infix)}
        assert experts, "no expert tier survived, so this asserts nothing"
        for expert in experts:
            assert expert in notes, expert
            assert expert[: -len(CONVENTIONS.expert_infix)] in notes, expert
        assert set(notes) == experts | {e[: -len(CONVENTIONS.expert_infix)] for e in experts}


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
