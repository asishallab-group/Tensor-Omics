"""The R pipeline: generate, build, and drive it through R.

The R analogue of `test_end_to_end.py`. It generates the C++ and R from the fixtures,
compiles the C++ against the fixture Fortran library with `R CMD SHLIB`, and runs a
driver script that sources everything and checks the answers and the classed conditions.

Skipped without gfortran, R and Rcpp.
"""

import os
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

from codegen.abi.c_abi import build_project
from codegen.config import Paths
from codegen.diagnostics import DiagnosticBag
from codegen.emit.errors_r import RErrorEmitter
from codegen.emit.fortran_c import FortranCEmitter
from codegen.emit.r_wrapper import RWrapperEmitter
from codegen.emit.rcpp import RcppEmitter
from codegen.frontend.ford_frontend import FordFrontend
from codegen.ir.errors import ErrorCatalogue
from codegen.ir.roles import analyse_project
from codegen.ir.validate import validate_project

from conftest import REPO_ROOT

GFORTRAN = shutil.which("gfortran")
RSCRIPT = shutil.which("Rscript")
FIXTURE_SRC = Path("helper/codegen/tests/fixtures/src")


def _has_rcpp() -> bool:
    if RSCRIPT is None:
        return False
    result = subprocess.run(
        [RSCRIPT, "-e", 'quit(status = !requireNamespace("Rcpp", quietly = TRUE))'],
        capture_output=True,
    )
    return result.returncode == 0


pytestmark = pytest.mark.skipif(
    GFORTRAN is None or RSCRIPT is None or not _has_rcpp(),
    reason="needs gfortran, R and Rcpp",
)


def _generate(out: Path) -> None:
    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=FIXTURE_SRC), bag).parse()
    analyse_project(parsed.project, bag)
    validate_project(parsed.project, bag)
    interface = build_project(parsed.project, bag)
    assert bag.errors == (), bag.render()

    fortran, cpp, r = FortranCEmitter(), RcppEmitter(), RWrapperEmitter()
    (out / "tox_marshal.h").write_text(cpp.marshal_header_content())
    (out / "tox_validate.R").write_text(r.validators())
    for module in interface:
        (out / f"{module.name}.F90").write_text(fortran.module(module))
        (out / f"{module.stripped_name}.cpp").write_text(cpp.module(module))
        (out / f"{module.stripped_name}.R").write_text(r.module(module))

    errors = DiagnosticBag()
    real = FordFrontend(Paths(root=REPO_ROOT, src_dir=Path("src/tox")), errors).parse()
    catalogue = ErrorCatalogue.from_module(
        real.project.module("tox_errors"), errors,
        real.project.constant_values(), arg_pos_factor=real.arg_pos_factor,
    )
    (out / "error_handling.R").write_text(RErrorEmitter(catalogue).module())


def _build_fortran(out: Path) -> Path:
    for name in ("fx_basics.F90", "fx_edges.F90"):
        shutil.copy(REPO_ROOT / FIXTURE_SRC / name, out / name)
    flags = ["-cpp", "-I.", "-std=f2018", "-ffree-line-length-none", "-fPIC", f"-J{out}"]
    objects = []
    for source in (
        REPO_ROOT / "src/tox/tox_errors.F90",
        REPO_ROOT / "src/tox/tox_conversions.F90",
        REPO_ROOT / "src/safeguard.F90",
        out / "fx_basics.F90", out / "fx_edges.F90",
        out / "fx_basics_c.F90", out / "fx_edges_c.F90",
    ):
        obj = out / f"{Path(source).stem}.o"
        result = subprocess.run(
            [GFORTRAN, *flags, "-c", str(source), "-o", str(obj)],
            cwd=REPO_ROOT, capture_output=True, text=True,
        )
        assert result.returncode == 0, f"{source}:\n{result.stderr}"
        objects.append(str(obj))
    library = out / "libfixtures.so"
    result = subprocess.run(
        [GFORTRAN, "-shared", "-o", str(library), *objects],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    return library


@pytest.fixture(scope="session")
def built(tmp_path_factory):
    out = tmp_path_factory.mktemp("e2e_r")
    _generate(out)
    _build_fortran(out)
    return out


def run_r(built: Path, body: str) -> str:
    """Source the generated R, run `body`, return its stdout."""
    script = textwrap.dedent(f"""
        suppressMessages({{
          Rcpp::sourceCpp("fx_basics.cpp")
          Rcpp::sourceCpp("fx_edges.cpp")
        }})
        source("error_handling.R"); source("tox_validate.R")
        source("fx_basics.R"); source("fx_edges.R")
        {body}
    """)
    environ = dict(os.environ)
    environ["PKG_CPPFLAGS"] = f"-I{built}"
    environ["PKG_LIBS"] = f"-L{built} -lfixtures -Wl,-rpath,{built}"
    result = subprocess.run(
        [RSCRIPT, "-e", script], cwd=built, capture_output=True, text=True, env=environ,
    )
    assert result.returncode == 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    return result.stdout


class TestResultsAreRight:
    def test_an_inout_array_comes_back_and_the_caller_is_untouched(self, built):
        # R is copy-on-modify: the modification is returned, the caller's vector stays
        out = run_r(built, """
            v <- c(3, 4)
            r <- fx_normalize(v)
            cat(all.equal(r, c(0.6, 0.8)), all(v == c(3, 4)))
        """)
        assert out.split() == ["TRUE", "TRUE"]

    def test_a_column_major_matrix_passes_without_a_copy(self, built):
        out = run_r(built, 'cat(fx_sum_matrix(matrix(c(1,2,3,4), 2), c(10, 100)))')
        assert float(out) == pytest.approx(730.0)

    def test_a_function_result_comes_back(self, built):
        out = run_r(built, 'cat(fx_count_positive(c(-1, 2, 3)))')
        assert int(out) == 2

    def test_a_default_is_applied_when_the_argument_is_omitted(self, built):
        # fx_optionals has span/max_iter/use_quantile with DM_DEFAULT
        out = run_r(built, 'cat(is.null(fx_optionals(c(1, 2, 3))))')
        assert out.strip() == "TRUE"


class TestModes:
    def test_a_mode_string_selects_the_branch(self, built):
        out = run_r(built, """
            cat(fx_modes(c(1,2,30), "mean", "ward"),
                fx_modes(c(1,2,30), "median", "ward"))
        """)
        assert [float(x) for x in out.split()] == [pytest.approx(11.0), pytest.approx(2.0)]

    def test_an_unknown_mode_is_a_classed_type_error(self, built):
        out = run_r(built, """
            cat(tryCatch(fx_modes(c(1,2,3), "nonsense", "ward"),
                         tox_type_error = function(e) "caught"))
        """)
        assert out.strip() == "caught"


class TestCharacters:
    def test_a_string_vector_comes_back(self, built):
        out = run_r(built, 'cat(fx_labels(c(1, -2, 3))$labels)')
        assert out.split() == ["pos", "nonpos", "pos"]

    def test_a_scalar_string_comes_back_unwrapped(self, built):
        out = run_r(built, 'cat(fx_labels(c(1, -2, 3))$label)')
        assert out.strip() == "summary"

    def test_a_string_vector_goes_in(self, built):
        out = run_r(built, 'cat(fx_count_matching(c("a", "b", "a"), "a"))')
        assert int(out) == 2


class TestOptionals:
    def test_a_nullable_optional_may_be_given_or_omitted(self, built):
        out = run_r(built, """
            given <- fx_nullable(c(1,2,3), "ungrouped", c(1L,1L,2L))
            omitted <- fx_nullable(c(1,2,3), "ungrouped")
            cat(is.null(given), is.null(omitted))
        """)
        assert out.split() == ["TRUE", "TRUE"]


class TestClassedConditions:
    def test_a_shape_mismatch_is_a_tox_shape_error_naming_both(self, built):
        out = run_r(built, """
            cat(tryCatch(fx_sum_matrix(matrix(c(1,2,3,4), 2), c(1,2,3)),
                         tox_shape_error = function(e) conditionMessage(e)))
        """)
        assert "'weights'" in out and "'matrix'" in out

    def test_a_wrong_type_is_a_tox_type_error_naming_the_argument(self, built):
        out = run_r(built, """
            cat(tryCatch(fx_normalize("hello"),
                         tox_type_error = function(e) conditionMessage(e)))
        """)
        assert "'vector'" in out

    def test_a_fortran_error_is_a_classed_condition(self, built):
        out = run_r(built, """
            cat(tryCatch(fx_normalize(numeric(3)),
                         tox_input_error = function(e) e$name))
        """)
        assert out.strip() == "ERR_DIVISION_BY_ZERO"

    def test_every_tox_error_is_catchable_as_the_parent(self, built):
        out = run_r(built, """
            cat(tryCatch(fx_normalize(numeric(3)), tox_error = function(e) "caught"))
        """)
        assert out.strip() == "caught"

    def test_an_na_in_an_integer_input_is_rejected(self, built):
        # fx_nullable's ortholog_set is integer; NA_integer_ is INT_MIN to Fortran
        out = run_r(built, """
            cat(tryCatch(fx_nullable(c(1,2,3), "ungrouped", c(1L, NA, 2L)),
                         tox_na_error = function(e) "caught"))
        """)
        assert out.strip() == "caught"
