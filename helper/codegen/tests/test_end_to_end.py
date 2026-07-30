"""The whole pipeline: generate, build a shared library, and call it from Python.

Everything else checks a layer. This checks that the layers agree -- that the Python the
generator writes calls the C wrapper the generator writes, and gets the answer the
Fortran actually computed. Nothing short of running it can show that.

Skipped without gfortran.
"""

import inspect
import os
import shutil
import subprocess
import sys

import pytest

from codegen.abi.c_abi import build_project
from codegen.config import Paths
from codegen.diagnostics import DiagnosticBag
from codegen.emit.errors_python import PythonErrorEmitter
from codegen.emit.fortran_c import FortranCEmitter
from codegen.emit.python_ctypes import PythonEmitter
from codegen.frontend.ford_frontend import FordFrontend
from codegen.ir.errors import ErrorCatalogue
from codegen.ir.roles import analyse_project
from codegen.ir.validate import validate_project

from conftest import REPO_ROOT

from pathlib import Path

GFORTRAN = shutil.which("gfortran")
FIXTURE_SRC = Path("helper/codegen/tests/fixtures/src")

pytestmark = pytest.mark.skipif(GFORTRAN is None, reason="gfortran is not installed")


def _generate(out: Path) -> None:
    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=FIXTURE_SRC), bag).parse()
    analyse_project(parsed.project, bag)
    validate_project(parsed.project, bag)
    interface = build_project(parsed.project, bag)
    assert bag.errors == (), bag.render()

    for module in interface:
        (out / f"{module.name}.F90").write_text(FortranCEmitter().module(module))

    emitter = PythonEmitter(library="build/libfixtures.so")
    package = out / "tensor_omics"
    package.mkdir(exist_ok=True)
    (package / "library.py").write_text(emitter.library_module())
    (package / "__init__.py").write_text(emitter.package_init(list(interface)))
    for module in interface:
        (package / f"{module.stripped_name}.py").write_text(emitter.module(module))

    errors = DiagnosticBag()
    real = FordFrontend(Paths(root=REPO_ROOT, src_dir=Path("src/tox")), errors).parse()
    catalogue = ErrorCatalogue.from_module(
        real.project.module("tox_errors"),
        errors,
        real.project.constant_values(),
        arg_pos_factor=real.arg_pos_factor,
    )
    (package / "error_handling.py").write_text(PythonErrorEmitter(catalogue).module())


def _build(out: Path) -> Path:
    for name in ("fx_basics.F90", "fx_edges.F90"):
        shutil.copy(REPO_ROOT / FIXTURE_SRC / name, out / name)

    flags = ["-cpp", "-I.", "-std=f2018", "-ffree-line-length-none", "-fPIC", f"-J{out}"]
    sources = [
        REPO_ROOT / "src/tox/tox_errors.F90",
        REPO_ROOT / "src/tox/tox_conversions.F90",
        REPO_ROOT / "src/safeguard.F90",
        out / "fx_basics.F90",
        out / "fx_edges.F90",
        out / "fx_basics_c.F90",
        out / "fx_edges_c.F90",
    ]
    objects = []
    for source in sources:
        obj = out / f"{Path(source).stem}.o"
        result = subprocess.run(
            [GFORTRAN, *flags, "-c", str(source), "-o", str(obj)],
            cwd=REPO_ROOT, capture_output=True, text=True,
        )
        assert result.returncode == 0, f"{source}:\n{result.stderr}"
        objects.append(str(obj))

    library = out / "build" / "libfixtures.so"
    library.parent.mkdir(exist_ok=True)
    result = subprocess.run(
        [GFORTRAN, "-shared", "-o", str(library), *objects],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    return library


@pytest.fixture(scope="session")
def tox(tmp_path_factory):
    """The generated package, built and imported."""
    out = tmp_path_factory.mktemp("e2e")
    _generate(out)
    library = _build(out)

    os.environ["TENSOR_OMICS_LIBRARY"] = str(library)
    sys.path.insert(0, str(out))
    try:
        import tensor_omics
    finally:
        sys.path.remove(str(out))
    return tensor_omics


@pytest.fixture
def np():
    return pytest.importorskip("numpy")


class TestItRuns:
    def test_the_package_imports_and_exports_the_wrappers(self, tox):
        assert "fx_normalize" in tox.__all__
        assert "fx_sum_matrix" in tox.__all__

    def test_the_alloc_variant_takes_the_plain_name(self, tox):
        assert "fx_cluster" in tox.__all__

    def test_a_wrapper_using_auto_output_from_is_generated(self, tox):
        # fx_cluster (the expert twin) gets n_work via DM_OUTPUT_FROM(AUTO): it calls
        # fx_work_size itself, so the caller never supplies the work size
        assert "fx_cluster_expert" in tox.__all__


class TestAutoOutputFrom:
    def test_the_work_size_is_never_asked_for(self, tox):
        import inspect

        # n_work is computed by calling fx_work_size, not a parameter
        assert list(inspect.signature(tox.fx_cluster_expert).parameters) == ["values"]

    def test_it_computes_the_work_size_and_runs(self, tox, np):
        # fx_work_size returns 2*n_values; fx_cluster returns count(values > 0)
        assert tox.fx_cluster_expert(np.array([-1.0, 2, 3, -4, 5])) == 3

    def test_the_producer_is_callable_directly_too(self, tox):
        assert tox.fx_work_size(4) == 8

    def test_a_producer_in_another_module_is_asked_for_too(self, tox):
        import inspect

        # fx_cross_module lives in fx_basics, fx_work_size in fx_edges
        assert list(inspect.signature(tox.fx_cross_module).parameters) == ["values"]

    def test_a_cross_module_producer_runs(self, tox, np):
        # the work array is sized by fx_edges::fx_work_size, and the sum comes back
        assert tox.fx_cross_module(np.array([1.0, 2.0, 3.0])) == pytest.approx(6.0)

    def test_a_renamed_producer_input_is_asked_for_by_the_producer_name(self, tox):
        # the consumer calls it n_samples, the producer n_values; the keyword has to be
        # the producer's, or the call is a TypeError
        source = (Path(tox.__file__).parent / "fx_edges.py").read_text()
        body = source.split("def fx_renamed_input")[1]
        assert "fx_work_size(n_values=n_samples)" in body

    def test_a_renamed_producer_input_runs(self, tox, np):
        # fx_work_size returns 2*n, and fx_renamed_input fills the work array with it
        assert list(inspect.signature(tox.fx_renamed_input).parameters) == ["samples"]
        assert tox.fx_renamed_input(np.array([1.0, 2.0, 3.0])) is None

    def test_a_constant_producer_input_is_passed_outright(self, tox, np):
        # the consumer has no argument for `factor`, so the table gives it a value
        assert list(inspect.signature(tox.fx_constant_input).parameters) == ["samples"]
        source = (Path(tox.__file__).parent / "fx_edges.py").read_text()
        body = source.split("def fx_constant_input")[1]
        assert "fx_scaled_size(n_values=n_samples, factor=3)" in body
        assert tox.fx_constant_input(np.array([1.0, 2.0])) is None

    def test_the_cross_module_import_is_deferred(self, tox):
        # a module-level import would be circular the moment two modules size each
        # other's outputs, so it sits inside the function
        source = (Path(tox.__file__).parent / "fx_basics.py").read_text()
        body = source.split("def fx_cross_module")[1]
        assert "from .fx_edges import fx_work_size" in body


class TestResultsAreRight:
    def test_an_inout_array_is_modified_in_place_without_a_copy(self, tox, np):
        # the whole point of not converting an intent(inout) array
        vector = np.array([3.0, 4.0], dtype=np.float64)

        assert tox.fx_normalize(vector) is None
        assert np.allclose(vector, [0.6, 0.8])

    def test_the_extent_is_never_asked_for(self, tox, np):
        # n_dims comes from the array; asking for it would be asking twice
        import inspect

        assert list(inspect.signature(tox.fx_normalize).parameters) == ["vector"]

    def test_a_column_major_matrix_arrives_intact(self, tox, np):
        matrix = np.asfortranarray([[1.0, 2.0], [3.0, 4.0]])
        weights = np.array([10.0, 100.0])

        # 10*(1+3) + 100*(2+4)
        assert tox.fx_sum_matrix(matrix, weights) == pytest.approx(640.0)

    def test_a_c_ordered_matrix_is_converted_rather_than_rejected(self, tox, np):
        matrix = np.array([[1.0, 2.0], [3.0, 4.0]])  # C order

        assert tox.fx_sum_matrix(matrix, np.array([10.0, 100.0])) == pytest.approx(640.0)

    def test_an_output_extent_is_sourced_from_the_input_not_the_output(self, tox):
        # n_out is an extent of both the intent(out) `averaged` and the intent(in)
        # `group_sizes`, with the output declared first. It must be read off the input:
        # the output is allocated *from* n_out and does not exist yet.
        # Regression -- the emitter used to read `averaged.shape[0]`, an UnboundLocalError.
        source = (Path(tox.__file__).parent / "fx_basics.py").read_text()
        body = source.split("def fx_grouped_output")[1]
        assert "n_out = group_sizes.shape[0]" in body

    def test_a_shared_output_extent_wrapper_runs(self, tox, np):
        # averaged[i_out, i_col] = values[i_col] * group_sizes[i_out]
        averaged = tox.fx_grouped_output(
            np.array([1.0, 2.0]), np.array([3, 5], dtype=np.int32))
        assert averaged.shape == (2, 2)
        assert np.allclose(averaged, [[3.0, 6.0], [5.0, 10.0]])

    def test_a_list_is_accepted_as_an_input_array(self, tox, np):
        assert tox.fx_count_positive([-1.0, 2.0, 3.0]) == 2

    def test_a_function_result_comes_back_as_the_return_value(self, tox, np):
        assert tox.fx_count_positive(np.array([-1.0, 2.0, 3.0])) == 2


class TestModes:
    def test_a_mode_travels_as_a_string_and_selects_the_branch(self, tox, np):
        values = np.array([1.0, 2.0, 30.0])

        assert tox.fx_modes(values, "mean", "ward") == pytest.approx(11.0)
        assert tox.fx_modes(values, "median", "ward") == pytest.approx(2.0)

    def test_an_unknown_mode_is_rejected_by_the_wrapper(self, tox, np):
        with pytest.raises(tox.ToxError):
            tox.fx_modes(np.array([1.0]), "nonsense", "ward")


class TestCharacters:
    def test_a_scalar_string_comes_back_decoded(self, tox, np):
        result = tox.fx_labels(np.array([1.0, -2.0, 3.0]))

        assert result["label"] == "summary"
        assert isinstance(result["label"], str)

    def test_a_string_vector_comes_back_as_a_list_of_str(self, tox, np):
        result = tox.fx_labels(np.array([1.0, -2.0, 3.0]))

        assert result["labels"] == ["pos", "nonpos", "pos"]
        assert all(isinstance(label, str) for label in result["labels"])

    def test_a_short_string_in_a_long_buffer_is_not_read_past(self, tox, np):
        # 'pos' is three characters in an eight-character buffer; the null padding is
        # what stops the trailing bytes being read back as part of the string
        labels = tox.fx_labels(np.array([1.0]))["labels"]

        assert labels == ["pos"]

    def test_a_string_vector_goes_in(self, tox):
        names = ["alpha", "beta", "alpha", "gamma"]

        assert tox.fx_count_matching(names, "alpha") == 2
        assert tox.fx_count_matching(names, "beta") == 1
        assert tox.fx_count_matching(names, "zzz") == 0

    def test_an_assumed_length_scalar_carries_its_own_length(self, tox):
        # len=* : numpy sizes the S dtype from the string, and the wrapper reads it back
        assert tox.fx_count_matching(["ab", "abc", "ab"], "ab") == 2
        assert tox.fx_count_matching(["ab", "abc", "ab"], "abc") == 1

    def test_a_non_iterable_for_a_string_vector_names_the_argument(self, tox):
        # str() of anything succeeds, so the failure is a non-iterable, not a bad element
        with pytest.raises(TypeError, match="'names' must be a sequence of strings"):
            tox.fx_count_matching(42, "x")


class TestShapeChecking:
    def test_a_shared_extent_that_disagrees_is_caught(self, tox, np):
        # Fortran cannot make this check: it would be a wrong answer or a segfault
        matrix = np.asfortranarray([[1.0, 2.0], [3.0, 4.0]])

        with pytest.raises(ValueError, match="implies n_cols == 2"):
            tox.fx_sum_matrix(matrix, np.array([1.0, 2.0, 3.0]))

    def test_the_message_names_both_arguments(self, tox, np):
        matrix = np.asfortranarray([[1.0, 2.0], [3.0, 4.0]])

        with pytest.raises(ValueError, match="'weights'.*'matrix'"):
            tox.fx_sum_matrix(matrix, np.array([1.0, 2.0, 3.0]))


class TestBadInputsAreNamed:
    def test_a_value_that_cannot_be_an_array_names_the_argument(self, tox):
        # numpy raises "could not convert string to float", which does not say which
        # argument was wrong -- it cannot know
        with pytest.raises(TypeError, match="'values' must be an array of np.float64"):
            tox.fx_count_positive("hello")

    def test_the_underlying_reason_is_kept(self, tox):
        with pytest.raises(TypeError, match="could not convert"):
            tox.fx_count_positive("hello")

    def test_a_wrong_type_in_a_second_argument_names_that_one(self, tox, np):
        matrix = np.asfortranarray([[1.0, 2.0], [3.0, 4.0]])

        with pytest.raises(TypeError, match="'weights'"):
            tox.fx_sum_matrix(matrix, "oops")

    def test_an_inout_argument_is_rejected_rather_than_converted(self, tox, np):
        # converting would copy, and the caller would never see the modification
        with pytest.raises(TypeError, match="'vector' is modified in place"):
            tox.fx_normalize([3.0, 4.0])

    def test_an_inout_argument_of_the_wrong_dtype_is_rejected(self, tox, np):
        with pytest.raises(TypeError, match="'vector' is modified in place"):
            tox.fx_normalize(np.array([3, 4], dtype=np.int32))


class TestErrors:
    def test_a_fortran_error_arrives_as_an_exception(self, tox, np):
        # ierr never reaches the caller
        with pytest.raises(tox.ToxError, match="Division by zero"):
            tox.fx_normalize(np.zeros(3))

    def test_the_message_is_the_fortran_documentation(self, tox, np):
        with pytest.raises(tox.ToxError) as excinfo:
            tox.fx_normalize(np.zeros(3))

        assert excinfo.value.name == "ERR_DIVISION_BY_ZERO"

    def test_a_successful_call_returns_no_error_code(self, tox, np):
        assert tox.fx_normalize(np.array([3.0, 4.0])) is None
