"""Every name the generated code binds beside the arguments is one `names` accounts for.

No argument may meet a name the generated code introduces itself. Python and the R binding
make that impossible rather than checking it: everything they bind starts with `_` (Python, the
R shim's C) or `.` (the R wrapper), which no Fortran name can. A Fortran wrapper cannot, so its
locals are `<argument><suffix>` from `names`, and `emit.name_clashes` refuses an argument where
the wrapper emits the same name. Both only hold while the emitters keep to it -- an emitter
that starts binding a plain `status` would bring the silent clash back -- so this reads the
generated code back, for the fixtures (a prologue among them) and for the real sources.
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
from codegen.emit.r_wrapper import RWrapperEmitter
from codegen.frontend.ford_frontend import FordFrontend
from codegen.ir.roles import analyse_project
from codegen.ir.validate import validate_project
from codegen.synthesize import synthesize_wrappers
from conftest import REPO_ROOT

SOURCES = {
    "fixtures": Path("helper/codegen/tests/fixtures/src"),
    # fx_prologue_impl is the prologue: `handled` and the prologue-produced local
    "fixture impls": Path("helper/codegen/tests/fixtures/impl"),
    "src": Path("src"),
}

#: The suffixes of a Fortran wrapper's locals
FORTRAN_SUFFIXES = (names.CONVERTED, names.MODE_VALUE, names.MODE_POINTER, names.VIEW_ADDRESS,
                    names.DEFAULT_VALUE)
#: What a Fortran C wrapper module declares beside its wrappers
FORTRAN_C_MODULE_VARIABLES = (names.EMPTY_VIEW_ANCHOR,)


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


def test_the_fixture_impls_have_a_prologue():
    # the premise of the prologue path being read back at all
    text = (REPO_ROOT / SOURCES["fixture impls"] / "fx_prologue_impl.F90").read_text()
    assert "DM_PROLOGUE(" in text


def _wrappers(binding):
    return [wrapper for module in binding for wrapper in module]


def _fortran_local(name: str, arguments) -> bool:
    """Whether `name` is `<argument><suffix>` for a Fortran wrapper's suffix."""
    lowered = name.lower()
    return any(
        lowered == names.local_name(argument, suffix).lower()
        for argument in arguments for suffix in FORTRAN_SUFFIXES
    )


class TestPython:
    def test_a_module_binds_its_imports_under_underscored_names(self, pipeline):
        _, binding = pipeline
        emitter = PythonEmitter()
        for module in binding:
            tree = ast.parse(emitter.module(module))
            for node in tree.body:
                if isinstance(node, (ast.Import, ast.ImportFrom)):
                    for alias in node.names:
                        bound = alias.asname or alias.name
                        assert bound in names.PYTHON_BOUND and bound.startswith("_"), bound

    def test_a_function_binds_nothing_but_its_arguments_and_underscored_names(self, pipeline):
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
                    # a producer in another module, imported where it is called; `validate`
                    # refuses an argument named like one
                    assert all(a.name in exported for a in node.names)
                for name in bound:
                    assert name.startswith("_") or name in arguments, (wrapper.name, name)


#: A C declaration in a shim: `SEXP x`, `int n = ...`, `unsigned char* _c_x = ...`, `int _dim_x[2]`
_C_DECLARATION = re.compile(
    r"\b(?:SEXP|int|unsigned char|char|double|float|Rcomplex|size_t|int(?:8|16|32|64)_t)\b"
    r"\s*\*?\s*(\w+)\s*(?=[=\[;,)])"
)


class TestTheRShim:
    def test_every_local_is_an_argument_or_underscored(self, pipeline):
        _, binding = pipeline
        emitter = CCallEmitter()
        for wrapper in _wrappers(binding):
            wrapper = wrapper.with_constant_sizes_evaluated()
            arguments = [argument.name for argument in wrapper]
            for name in _C_DECLARATION.findall(emitter.function(wrapper)):
                assert name.startswith("_") or name in arguments, (wrapper.name, name)


#: An assignment in an R wrapper: `x <- ...`
_R_ASSIGNMENT = re.compile(r"^\s*([A-Za-z.][\w.]*)\s*<-")


class TestTheRWrapper:
    def test_every_local_is_an_argument_or_dotted(self, pipeline):
        _, binding = pipeline
        emitter = RWrapperEmitter()
        for module in binding:
            wrappers = {wrapper.stripped_name: wrapper for wrapper in module}
            current = None
            for line in emitter.module(module).splitlines():
                match = _R_ASSIGNMENT.match(line)
                if match is None:
                    continue
                name = match.group(1)
                if not line.startswith((" ", "\t")):
                    # a top-level assignment is a wrapper, and what follows is its body
                    assert name in wrappers, name
                    current = wrappers[name]
                    continue
                arguments = {argument.name for argument in current}
                assert name.startswith(".") or name in arguments, (current.name, name)


#: A Fortran declaration's entity: `real(real64), dimension(:) :: x_f`
_FORTRAN_DECLARATION = re.compile(r"^\s*(?!use\b|public\b|private\b)[^!\n]*::\s*(\w+)", re.M)


#: A declaration in a module's specification part: what is not a `use` or an access list
_MODULE_VARIABLE = re.compile(r"^(?![ \t]*(?:use|public|private)\b)[^!\n]*::\s*(\w+)", re.M)


class TestTheFortranCWrapper:
    def test_every_local_is_an_argument_or_built_from_the_table(self, pipeline):
        _, binding = pipeline
        emitter = FortranCEmitter()
        for wrapper in _wrappers(binding):
            arguments = [argument.name for argument in wrapper]
            for name in _FORTRAN_DECLARATION.findall(emitter.wrapper(wrapper)):
                if name.lower() in {a.lower() for a in arguments}:
                    continue
                assert _fortran_local(name, arguments), (wrapper.name, name)


    def test_every_module_variable_is_from_the_table(self, pipeline):
        _, binding = pipeline
        emitter = FortranCEmitter()
        for module in binding:
            header = emitter.module(module).partition("contains")[0]
            for name in _MODULE_VARIABLE.findall(header):
                assert name in FORTRAN_C_MODULE_VARIABLES, (module.name, name)


class TestTheFortranWrappers:
    def test_every_local_is_an_argument_or_built_from_the_table(self, pipeline):
        synthesis, _ = pipeline
        emitter = FortranWrapperEmitter(project=synthesis.project)
        specs, info = {}, {}
        for spec in synthesis.specs:
            for wrapper in (spec.validating, spec.allocating):
                if wrapper is not None:
                    specs[wrapper.name.lower()] = spec
                    info.setdefault(spec.module_name, {})[wrapper.name.lower()] = (
                        WrapperInfo(spec.impl.name, spec.mode_fix)
                    )
        for module in synthesis.project:
            if module.name not in synthesis.generated_names or not module.procedures:
                continue
            text = emitter.module(module, info.get(module.name))
            for procedure, body in _subroutines(text):
                spec = specs[procedure.lower()]
                # what the allocating tier takes over becomes a local of the same name
                arguments = [a.name for a in spec.impl.arguments]
                arguments += [a.name for a in module.procedure(procedure).arguments]
                for name in _FORTRAN_DECLARATION.findall(body):
                    if name.lower() in {a.lower() for a in arguments}:
                        continue
                    if spec.prologue is not None and name.lower() == names.FORTRAN_HANDLED:
                        continue
                    assert _fortran_local(name, arguments), (procedure, name)


def _subroutines(text: str):
    for match in re.finditer(r"subroutine (\w+)\(.*?end subroutine \1", text, re.S | re.I):
        yield match.group(1), match.group(0)


class TestTheTableItself:
    def test_the_excluded_names_are_not_in_it(self):
        # main publishes these, and none of them breaks a binding (see `names`)
        for name in ("value", "mode", "vector", "array", "expression", "factor", "points",
                     "weights", "match", "length", "count"):
            assert names.language_word(name) is None, name

    def test_the_handled_flag_is_the_prologue_convention(self):
        assert names.FORTRAN_HANDLED == CONVENTIONS.prologue_handled_arg

    def test_every_shim_local_is_underscored(self):
        for tag in (names.SCALAR_VALUE, names.OPTIONAL_POINTER, names.OPTIONAL_SIZE, names.DIM,
                    names.OPTIONAL_DIM_ATTRIBUTE, names.BUFFER, names.INOUT_COPY):
            # underscore, then a lower-case tag: legal C at block scope whatever the case
            assert names.shim_local("K", tag).startswith("_") and tag[0].islower()
        assert names.C_PROTECT_COUNT.startswith("_") and names.C_LOOP_INDEX.startswith("_")
