"""What `tensor_omics` binds at its top level (#221 C).

The package re-exports every procedure and everything the error module defines -- the codes,
`STAT_*` included, the exception classes, `ARG_POS_FACTOR` and `check_err_code` -- from the one
list the error module is emitted from, so a code added to `tox_errors` is exported with no
other edit. A procedure or module that would take one of those names, or a submodule's, is a
generator error. Drives `generate` on hand-built projects, as `test_name_clashes` does.
"""

import importlib.util
from pathlib import Path

import pytest

from codegen.diagnostics import DiagnosticBag, SourceLocation
from codegen.emit.errors_python import PythonErrorEmitter, public_names
from codegen.frontend.ford_frontend import ParsedProject
from codegen.frontend.macros import MacroTable
from codegen.generate import generate
from codegen.ir.entities import Project
from codegen.ir.errors import ErrorCatalogue
from codegen.ir.types import Intent

import builders as b
from conftest import REPO_ROOT

MACROS = MacroTable(REPO_ROOT / "src/macros.h", include_paths=(REPO_ROOT,))

CODES = (
    b.parameter("ERR_OK", "0", doc="no error, operation successful"),
    b.parameter("ERR_INVALID_INPUT", "201", doc="invalid input arguments"),
    b.parameter("ERR_DIVISION_BY_ZERO", "210", doc="division by zero"),
    b.parameter("ERR_ALLOC_FAIL", "301", doc="memory allocation failed"),
)


def errors_module(*codes):
    return b.module("tox_errors", parameters=CODES + codes)


def exported(name, module="tox_demo"):
    """A module holding one hand-written export called `name`."""
    return b.module(module, b.procedure(
        name, b.integer("n"), b.real("x", Intent.IN, "(n)"), b.ierr(),
    ))


def run(*modules, extra_codes=()):
    parsed = ParsedProject(project=Project([errors_module(*extra_codes), *modules]),
                           macros=MACROS, arg_pos_factor=10000)
    return generate(targets=("python",), parsed=parsed)


def written(result, name):
    return next(f.content for f in result.files if f.path.name == name)


def errors(*modules):
    return [d.message for d in run(*modules).diagnostics.errors]


class TestTheExports:
    def test_init_imports_and_lists_every_public_name_of_the_error_module(self):
        result = run(exported("demo_sum"))
        assert result.ok, result.diagnostics.render()
        init = written(result, "__init__.py")

        imports = init.split("from .error_handling import (")[1].split(")")[0]
        listed = init.split("__all__ = [")[1]
        for name in ("ARG_POS_FACTOR", "ToxError", "ToxInputError", "ToxMemoryError",
                     "ToxArgumentError", "ERR_OK", "ERR_DIVISION_BY_ZERO", "check_err_code"):
            assert f"    {name},\n" in imports
            assert f'    "{name}",\n' in listed
        assert '    "demo_sum",\n' in listed

    def test_a_code_added_to_the_table_is_exported_with_no_other_edit(self):
        # #220 brings the first STAT_* codes; they must follow without touching the generator
        result = run(exported("demo_sum"), extra_codes=(
            b.parameter("STAT_NO_STABLE_DIRECTION", "1", doc="no stable direction"),
            b.parameter("ERR_SOMETHING_NEW", "211", doc="something new went wrong"),
        ))
        assert result.ok, result.diagnostics.render()
        init = written(result, "__init__.py")
        for name in ("STAT_NO_STABLE_DIRECTION", "ERR_SOMETHING_NEW"):
            assert f"    {name},\n" in init
            assert f'    "{name}",\n' in init

    def test_the_list_is_exactly_what_the_error_module_defines(self, tmp_path):
        """The helper is the module's public surface: nothing more, nothing it lacks."""
        bag = DiagnosticBag()
        catalogue = ErrorCatalogue.from_module(errors_module(
            b.parameter("STAT_CONVERGED", "1", doc="converged")), bag)
        assert bag.errors == (), bag.render()
        path = tmp_path / "error_handling.py"
        path.write_text(PythonErrorEmitter(catalogue).module())
        spec = importlib.util.spec_from_file_location("surface_error_handling", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        defined = {name for name in vars(module) if not name.startswith("_")}
        assert defined == set(public_names(catalogue))
        assert len(public_names(catalogue)) == len(defined)


class TestAClashIsAnError:
    @pytest.mark.parametrize("name", ["check_err_code", "err_ok", "Err_Division_By_Zero",
                                      "toxinputerror", "arg_pos_factor"])
    def test_a_procedure_named_like_an_export_of_the_error_module(self, name):
        found = errors(exported(name))
        assert len(found) == 1, found
        assert "of the error module and the procedure" in found[0]
        assert f"the procedure '{name}' of 'tox_demo'" in found[0]
        assert "top level of the Python package" in found[0]

    @pytest.mark.parametrize("name", ["library", "error_handling", "Library"])
    def test_a_procedure_named_like_a_submodule_of_the_package(self, name):
        found = errors(exported(name))
        assert found == [
            f"the package's own submodule '{name.lower()}' and the procedure '{name}' of "
            "'tox_demo' have the same name at the top level of the Python package"
        ]

    def test_a_procedure_named_like_its_own_module(self):
        # `from .tox_demo import tox_demo` would replace the submodule by the function
        assert errors(exported("tox_demo")) == [
            "the module 'tox_demo' and the procedure 'tox_demo' of 'tox_demo' have the same "
            "name at the top level of the Python package"
        ]

    def test_a_module_whose_python_file_is_the_error_module_s(self):
        # it would be written to error_handling.py, over the error module
        found = errors(exported("demo_sum", module="error_handling"))
        assert found == [
            "the package's own submodule 'error_handling' and the module 'error_handling' have "
            "the same name at the top level of the Python package"
        ]

    def test_names_that_only_look_alike_are_fine(self):
        assert errors(exported("check_err_codes"), exported("err_ok_count", module="tox_other"),
                      exported("libraries", module="tox_third")) == []

    def test_the_error_points_at_the_procedure(self):
        where = SourceLocation(file=Path("src/tox/tox_demo.F90"), line=12)
        demo = b.module("tox_demo", b.procedure(
            "check_err_code", b.integer("n"), b.real("x", Intent.IN, "(n)"), b.ierr(),
            location=where,
        ))
        diagnostic = run(demo).diagnostics.errors[0]
        assert diagnostic.location == where
        assert diagnostic.note and "ignoring case" in diagnostic.note
