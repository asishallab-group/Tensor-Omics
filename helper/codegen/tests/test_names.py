"""The name table holds every name the generated code introduces.

`ir.validate` refuses a published argument named like a name the generated code introduces
itself, and it reads those names from `names`. That only protects anything while the table is
complete -- an emitter that starts declaring a new local without adding it would bring the
silent clash back. So this reads the generated code back, for the fixtures and for the real
sources, and checks that every name it binds beside the arguments is one the table knows.
"""

import ast
import re
from pathlib import Path

import pytest

from codegen import names
from codegen.abi.c_abi import build_project
from codegen.config import CONVENTIONS, Paths
from codegen.diagnostics import DiagnosticBag
from codegen.emit.c_call import CCallEmitter
from codegen.emit.fortran_c import FortranCEmitter
from codegen.emit.fortran_wrapper import FortranWrapperEmitter, WrapperInfo
from codegen.emit.python_ctypes import PythonEmitter
from codegen.frontend.ford_frontend import FordFrontend
from codegen.ir.roles import analyse_project
from codegen.ir.validate import validate_project
from codegen.synthesize import synthesize_wrappers
from conftest import REPO_ROOT

SOURCES = {
    "fixtures": Path("helper/codegen/tests/fixtures/src"),
    "fixture impls": Path("helper/codegen/tests/fixtures/impl"),
    "src": Path("src"),
}


@pytest.fixture(scope="module", params=sorted(SOURCES))
def pipeline(request):
    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=SOURCES[request.param]), bag).parse()
    synthesis = synthesize_wrappers(parsed.project, CONVENTIONS, bag)
    analyse_project(synthesis.project, bag)
    validate_project(synthesis.project, bag)
    binding = build_project(synthesis.project, bag)
    assert bag.errors == (), bag.render()
    return synthesis, binding


def _wrappers(binding):
    return [wrapper for module in binding for wrapper in module]


def _built_from(name: str, arguments) -> bool:
    """Whether `name` is `<argument><suffix>` for one of the table's suffixes."""
    return names.clashing_local(name, arguments) is not None


class TestPython:
    def test_a_module_imports_only_what_the_table_names(self, pipeline):
        _, binding = pipeline
        emitter = PythonEmitter()
        for module in binding:
            tree = ast.parse(emitter.module(module))
            for node in tree.body:
                if isinstance(node, (ast.Import, ast.ImportFrom)):
                    for alias in node.names:
                        assert (alias.asname or alias.name) in names.PYTHON_BOUND, alias.name

    def test_a_function_binds_nothing_but_its_arguments_and_the_table(self, pipeline):
        _, binding = pipeline
        emitter = PythonEmitter()
        exported = {wrapper.stripped_name for wrapper in _wrappers(binding)}
        for wrapper in _wrappers(binding):
            function = ast.parse(emitter.function(wrapper.with_constant_sizes_evaluated()))
            arguments = {argument.name for argument in wrapper}
            for node in ast.walk(function.body[0]):
                bound = []
                if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
                    bound.append(node.id)
                elif isinstance(node, ast.ExceptHandler) and node.name:
                    bound.append(node.name)
                elif isinstance(node, ast.ImportFrom):
                    # a producer in another module, imported where it is called
                    bound += [a.asname or a.name for a in node.names if a.name in exported]
                    assert all(a.name in exported for a in node.names)
                for name in bound:
                    if name.startswith("_") or name in arguments or name in exported:
                        continue  # a Fortran name cannot start with `_`
                    assert name in names.PYTHON_BOUND, (wrapper.name, name)


#: A C declaration in a shim: `SEXP x`, `int n = ...`, `unsigned char* x_c = ...`, `int x_dim[2]`
_C_DECLARATION = re.compile(
    r"\b(?:SEXP|int|unsigned char|char|double|float|Rcomplex|size_t|int(?:8|16|32|64)_t)\b"
    r"\s*\*?\s*(\w+)\s*(?=[=\[;,)])"
)


class TestTheRShim:
    def test_every_local_is_an_argument_or_built_from_the_table(self, pipeline):
        _, binding = pipeline
        emitter = CCallEmitter()
        for wrapper in _wrappers(binding):
            wrapper = wrapper.with_constant_sizes_evaluated()
            text = emitter.function(wrapper)
            arguments = [argument.name for argument in wrapper]
            for name in _C_DECLARATION.findall(text):
                if name.startswith("_") or name == "i" or name in arguments:
                    continue  # `_out`, `_nms`; a loop index in its own block scope
                if name in names.GENERATED_LOCALS[names.R_SHIM]:
                    continue
                assert _built_from(name, arguments), (wrapper.name, name)


#: A Fortran declaration's entity: `real(real64), dimension(:) :: x_f`
_FORTRAN_DECLARATION = re.compile(r"^\s*(?!use\b|public\b|private\b)[^!\n]*::\s*(\w+)", re.M)


class TestTheFortranCWrapper:
    def test_every_local_is_an_argument_or_built_from_the_table(self, pipeline):
        _, binding = pipeline
        emitter = FortranCEmitter()
        for wrapper in _wrappers(binding):
            text = emitter.wrapper(wrapper)
            arguments = [argument.name for argument in wrapper]
            for name in _FORTRAN_DECLARATION.findall(text):
                if name.lower() in {a.lower() for a in arguments}:
                    continue
                assert _built_from(name, arguments), (wrapper.name, name)


class TestTheFortranWrappers:
    def test_every_local_is_an_argument_or_built_from_the_table(self, pipeline):
        synthesis, _ = pipeline
        emitter = FortranWrapperEmitter(project=synthesis.project)
        impls = {}
        info = {}
        for spec in synthesis.specs:
            for wrapper in (spec.validating, spec.allocating):
                if wrapper is not None:
                    impls[wrapper.name.lower()] = spec
                    info.setdefault(spec.module_name, {})[wrapper.name.lower()] = (
                        WrapperInfo(spec.impl.name, spec.mode_fix)
                    )
        for module in synthesis.project:
            if module.name not in synthesis.generated_names or not module.procedures:
                continue
            text = emitter.module(module, info.get(module.name))
            for procedure, body in _subroutines(text):
                spec = impls[procedure.lower()]
                # what the allocating tier takes over becomes a local of the same name
                arguments = [a.name for a in spec.impl.arguments]
                wrapper = spec.allocating if spec.allocating is not None and (
                    spec.allocating.name.lower() == procedure.lower()) else spec.validating
                arguments += [a.name for a in wrapper.arguments]
                if spec.prologue is not None:
                    arguments += [a.name for a in spec.prologue.arguments]
                for name in _FORTRAN_DECLARATION.findall(body):
                    if name.lower() in {a.lower() for a in arguments}:
                        continue
                    if name.lower() in {n.lower() for n in names.GENERATED_LOCALS[names.FORTRAN]}:
                        continue
                    assert _built_from(name, arguments), (procedure, name)


def _subroutines(text: str):
    for match in re.finditer(
        r"subroutine (\w+)\(.*?end subroutine \1", text, re.S | re.I
    ):
        yield match.group(1), match.group(0)


class TestTheTableItself:
    def test_no_suffix_is_a_name_on_its_own(self):
        # a suffix must start with `_`, or `<x><suffix>` would not read as built from `x`
        assert all(suffix.startswith("_") for suffix in names.GENERATED_SUFFIXES)

    def test_the_excluded_names_are_not_in_it(self):
        # main publishes these, and none of them breaks a binding (see `names`)
        for name in ("value", "mode", "vector", "array", "expression", "factor", "points",
                     "weights", "match", "length"):
            assert names.language_word(name) is None, name

    def test_the_handled_flag_is_the_prologue_convention(self):
        assert names.FORTRAN_HANDLED == CONVENTIONS.prologue_handled_arg
