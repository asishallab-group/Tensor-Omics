"""The Python error module.

Generated Python is tested by running it: each test execs the emitted module and calls
it. Asserting on the text would check that it looks right, which is not the question.
"""

import importlib.util

import pytest

from codegen_reworked.diagnostics import DiagnosticBag
from codegen_reworked.emit.errors_python import PythonErrorEmitter
from codegen_reworked.ir.errors import ErrorCatalogue

import builders as b


def code(name, value, doc):
    return b.parameter(name, str(value), doc=doc)


@pytest.fixture(scope="module")
def module(tmp_path_factory):
    """The emitted module, imported."""
    bag = DiagnosticBag()
    source = b.module(
        "tox_errors",
        parameters=(
            code("ERR_OK", 0, "no error, operation successful"),
            code("ERR_FILE_OPEN", 101, "could not open file"),
            code("ERR_INVALID_INPUT", 201, "invalid input arguments"),
            code("ERR_ALLOC_FAIL", 301, "memory allocation failed"),
            code("ERR_UNIT_NOT_CONNECTED", 5002, "unit not connected"),
            code("ERR_INTERNAL", 9001, "unexpected internal state"),
            code("STAT_CONVERGED", 1, "the algorithm converged"),
        ),
    )
    catalogue = ErrorCatalogue.from_module(source, bag)
    assert bag.errors == (), bag.render()

    path = tmp_path_factory.mktemp("errors") / "error_handling.py"
    path.write_text(PythonErrorEmitter(catalogue).module())

    spec = importlib.util.spec_from_file_location("error_handling", path)
    emitted = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(emitted)
    return emitted


class TestItRuns:
    def test_the_emitted_module_imports(self, module):
        assert module.check_err_code is not None

    def test_the_factor_comes_from_the_catalogue(self, module):
        assert module.ARG_POS_FACTOR == 10000

    def test_the_codes_are_exported_as_constants(self, module):
        assert module.ERR_INVALID_INPUT == 201
        assert module.ERR_OK == 0


class TestSuccess:
    def test_ok_does_not_raise_and_returns_nothing(self, module):
        # ierr never reaches a caller: a successful call returns its results only
        assert module.check_err_code(0) is None

    def test_an_argument_position_on_ok_is_still_ok(self, module):
        # is_err compares only the code, so must this
        assert module.check_err_code(10000 * 2) is None


class TestStatuses:
    def test_a_status_does_not_raise(self, module):
        # a status is an outcome, not a failure
        assert module.check_err_code(1) == "STAT_CONVERGED"

    def test_a_status_is_not_in_the_error_tables(self, module):
        assert 1 not in module._MESSAGES


class TestErrors:
    def test_an_error_raises_with_the_fortran_documentation_as_its_message(self, module):
        with pytest.raises(module.ToxError, match="invalid input arguments"):
            module.check_err_code(201)

    def test_the_exception_carries_the_code_and_its_name(self, module):
        with pytest.raises(module.ToxError) as excinfo:
            module.check_err_code(201)

        assert excinfo.value.code == 201
        assert excinfo.value.name == "ERR_INVALID_INPUT"
        assert excinfo.value.argument is None

    def test_an_unmapped_code_still_raises(self, module):
        with pytest.raises(module.ToxError, match="unmapped error code 777"):
            module.check_err_code(777)


class TestArgumentPositions:
    def test_the_offending_argument_is_named(self, module):
        # the position is packed into the code, and the wrapper knows its own arguments
        with pytest.raises(module.ToxError, match=r"argument 'n_dims'"):
            module.check_err_code(10000 * 2 + 201, arguments=("vector", "n_dims", "ierr"))

    def test_the_name_is_on_the_exception(self, module):
        with pytest.raises(module.ToxError) as excinfo:
            module.check_err_code(10000 * 1 + 201, arguments=("vector", "n_dims"))

        assert excinfo.value.argument == "vector"

    def test_without_names_the_position_is_reported(self, module):
        with pytest.raises(module.ToxError, match=r"\(argument 2\)"):
            module.check_err_code(10000 * 2 + 201)

    def test_a_position_beyond_the_names_falls_back_to_the_number(self, module):
        # rather than an IndexError from inside the error handler
        with pytest.raises(module.ToxError, match=r"\(argument 9\)"):
            module.check_err_code(10000 * 9 + 201, arguments=("vector",))

    def test_position_zero_names_no_argument(self, module):
        with pytest.raises(module.ToxError) as excinfo:
            module.check_err_code(201, arguments=("vector",))

        assert excinfo.value.argument is None
        # the message is 'invalid input arguments', so look for the parenthesised suffix
        assert "(argument" not in str(excinfo.value)


class TestExceptionHierarchy:
    @pytest.mark.parametrize(
        "code, expected",
        [
            (101, "ToxIOError"),
            (201, "ToxInputError"),
            (301, "ToxMemoryError"),
            (5002, "ToxRuntimeError"),
            (9001, "ToxInternalError"),
        ],
    )
    def test_the_range_decides_the_type(self, module, code, expected):
        with pytest.raises(module.ToxError) as excinfo:
            module.check_err_code(code)

        assert type(excinfo.value).__name__ == expected

    def test_every_type_is_catchable_as_one(self, module):
        # so a caller who does not care which kind can still catch them all
        for code in (101, 201, 301, 5002, 9001, 777):
            with pytest.raises(module.ToxError):
                module.check_err_code(code)

    def test_they_are_runtime_errors(self, module):
        assert issubclass(module.ToxError, RuntimeError)

    def test_only_the_groups_present_get_a_class(self, tmp_path):
        bag = DiagnosticBag()
        source = b.module("tox_errors", parameters=(code("ERR_INVALID_INPUT", 201, "bad"),))
        text = PythonErrorEmitter(ErrorCatalogue.from_module(source, bag)).module()

        assert "class ToxInputError" in text
        assert "class ToxIOError" not in text


@pytest.fixture(scope="session")
def real(tmp_path_factory):
    """The error module emitted from the real tox_errors, imported."""
    from pathlib import Path

    from codegen_reworked.config import Paths
    from codegen_reworked.frontend.ford_frontend import FordFrontend
    from conftest import REPO_ROOT

    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=Path("src/tox")), bag).parse()
    catalogue = ErrorCatalogue.from_module(
        parsed.project.module("tox_errors"),
        bag,
        parsed.project.constant_values(),
        arg_pos_factor=parsed.arg_pos_factor,
    )
    path = tmp_path_factory.mktemp("real_errors") / "error_handling.py"
    path.write_text(PythonErrorEmitter(catalogue).module())

    spec = importlib.util.spec_from_file_location("real_error_handling", path)
    emitted = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(emitted)
    return emitted, catalogue


class TestAgainstTheRealCatalogue:
    """Emitted from the real tox_errors, then run."""

    def test_it_runs(self, real):
        module, _ = real

        assert module.check_err_code(0) is None

    def test_every_real_code_is_present_and_raises(self, real):
        module, catalogue = real

        for error in catalogue.errors:
            assert getattr(module, error.name) == error.value
            with pytest.raises(module.ToxError) as excinfo:
                module.check_err_code(error.value)
            assert excinfo.value.name == error.name

    def test_the_messages_are_the_fortran_documentation(self, real):
        module, _ = real

        with pytest.raises(module.ToxError, match="null pointer dereference"):
            module.check_err_code(module.ERR_POINTER_NULL)
