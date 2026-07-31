"""The Python and R interfaces must ask their callers for the same things.

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
def interface():
    return _build(FIXTURE_SRC)


@pytest.fixture(scope="module")
def project_interface():
    """The real sources. The fixtures are deliberately small, and every divergence found
    so far was in `src/` -- the serde helpers, the outlier work arrays -- so checking the
    fixtures alone would have caught none of them."""
    return _build(Path("src"))


def _wrappers(interface):
    for module in interface:
        for wrapper in module:
            yield wrapper


def _names(arguments):
    return [argument.name for argument in arguments]


class TestTheTargetsAgree:
    def test_every_wrapper_asks_for_the_same_arguments(self, interface):
        python, r = PythonEmitter(), RWrapperEmitter()

        divergent = {
            wrapper.stripped_name: (_names(python._inputs(wrapper)), _names(r._inputs(wrapper)))
            for wrapper in _wrappers(interface)
            if _names(python._inputs(wrapper)) != _names(r._inputs(wrapper))
        }

        assert not divergent, (
            "Python and R disagree about what the caller supplies:\n"
            + "\n".join(f"  {name}:\n    python: {py}\n    r:      {r_}"
                        for name, (py, r_) in sorted(divergent.items()))
        )


    def test_r_calls_c_with_the_arguments_c_declares(self, interface):
        """The R wrapper and the C++ function it calls must agree, argument for argument.

        They are decided in two places -- `r_wrapper._call_inputs` and `c_call._inputs` --
        so a role change that reaches one and not the other produces an R function that
        calls `.Call("name_call", ...)` with the wrong arguments. C is compiled separately, so
        nothing catches it until the call happens.
        """
        r, c = RWrapperEmitter(), CCallEmitter()

        divergent = {
            wrapper.stripped_name: (_names(r._call_inputs(wrapper)), _names(c._inputs(wrapper)))
            for wrapper in _wrappers(interface)
            if _names(r._call_inputs(wrapper)) != _names(c._inputs(wrapper))
        }

        assert not divergent, (
            "the R wrapper and the C function disagree:\n"
            + "\n".join(f"  {name}:\n    r passes:    {passed}\n    c expects:   {declared}"
                        for name, (passed, declared) in sorted(divergent.items()))
        )

    def test_no_argument_is_both_asked_for_and_derived(self, interface):
        """Deriving an argument the signature also asks for overwrites what was passed.

        It is the failure mode that reads as a mystery: the caller supplies a value, the
        wrapper quietly replaces it with one computed from something that may not even
        exist yet.
        """
        emitter = PythonEmitter()

        overwritten = {}
        for wrapper in _wrappers(interface):
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
    def test_every_wrapper_asks_for_the_same_arguments(self, project_interface):
        python, r = PythonEmitter(), RWrapperEmitter()

        divergent = {
            wrapper.stripped_name: (_names(python._inputs(wrapper)), _names(r._inputs(wrapper)))
            for wrapper in _wrappers(project_interface)
            if _names(python._inputs(wrapper)) != _names(r._inputs(wrapper))
        }

        assert not divergent, (
            "Python and R disagree about what the caller supplies:\n"
            + "\n".join(f"  {name}:\n    python: {py}\n    r:      {r_}"
                        for name, (py, r_) in sorted(divergent.items()))
        )


    def test_r_calls_c_with_the_arguments_c_declares(self, project_interface):
        """The R wrapper and the C++ function it calls must agree, argument for argument.

        They are decided in two places -- `r_wrapper._call_inputs` and `c_call._inputs` --
        so a role change that reaches one and not the other produces an R function that
        calls `.Call("name_call", ...)` with the wrong arguments. C is compiled separately, so
        nothing catches it until the call happens.
        """
        r, c = RWrapperEmitter(), CCallEmitter()

        divergent = {
            wrapper.stripped_name: (_names(r._call_inputs(wrapper)), _names(c._inputs(wrapper)))
            for wrapper in _wrappers(project_interface)
            if _names(r._call_inputs(wrapper)) != _names(c._inputs(wrapper))
        }

        assert not divergent, (
            "the R wrapper and the C function disagree:\n"
            + "\n".join(f"  {name}:\n    r passes:    {passed}\n    c expects:   {declared}"
                        for name, (passed, declared) in sorted(divergent.items()))
        )

    def test_no_argument_is_both_asked_for_and_derived(self, project_interface):
        emitter = PythonEmitter()

        overwritten = {}
        for wrapper in _wrappers(project_interface):
            asked = set(_names(emitter._inputs(wrapper)))
            derived = {
                argument.name for argument in wrapper
                if argument.name in asked
                and emitter._extent_expression(argument, wrapper)
            }
            if derived:
                overwritten[wrapper.stripped_name] = sorted(derived)

        assert not overwritten, f"asked for and then overwritten: {overwritten}"
