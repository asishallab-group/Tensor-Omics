"""The R pipeline: generate, build, and drive it through R.

The R analogue of `test_end_to_end.py`. It generates the C and R from the fixtures and
compiles the Fortran fixtures *and* the C `.Call` shims into one `libfixtures.so` -- the
same bundling fpm does for libtensor-omics.so in production -- then runs a driver script
that `dyn.load`s it, sources everything, and checks the answers and the classed conditions.

Skipped without gfortran, gcc and R. No Rcpp -- the shims are pure C.
"""

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
from codegen.emit.c_call import CCallEmitter
from codegen.frontend.ford_frontend import FordFrontend
from codegen.ir.errors import ErrorCatalogue
from codegen.ir.roles import analyse_project
from codegen.ir.validate import validate_project

from conftest import REPO_ROOT

GFORTRAN = shutil.which("gfortran")
GCC = shutil.which("gcc") or shutil.which("cc")
RSCRIPT = shutil.which("Rscript")
RBIN = shutil.which("R")
FIXTURE_SRC = Path("helper/codegen/tests/fixtures/src")


pytestmark = pytest.mark.skipif(
    GFORTRAN is None or GCC is None or RSCRIPT is None or RBIN is None,
    reason="needs gfortran, gcc and R",
)


def _generate(out: Path) -> None:
    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=FIXTURE_SRC), bag).parse()
    analyse_project(parsed.project, bag)
    validate_project(parsed.project, bag)
    binding = build_project(parsed.project, bag)
    assert bag.errors == (), bag.render()

    fortran, c, r = FortranCEmitter(), CCallEmitter(), RWrapperEmitter()
    modules = list(binding)
    (out / "tox_marshal.h").write_text(c.marshal_header_content())
    (out / "init.c").write_text(c.registration(modules))
    (out / "tox_validate.R").write_text(r.validators())
    for module in modules:
        (out / f"{module.name}.F90").write_text(fortran.module(module))
        (out / f"{module.stripped_name}.c").write_text(c.module(module))
        (out / f"{module.stripped_name}.R").write_text(r.module(module))

    errors = DiagnosticBag()
    real = FordFrontend(Paths(root=REPO_ROOT, src_dir=Path("src/f42")), errors).parse()
    catalogue = ErrorCatalogue.from_module(
        real.project.module("tox_errors"), errors,
        real.project.constant_values(), arg_pos_factor=real.arg_pos_factor,
    )
    (out / "error_handling.R").write_text(RErrorEmitter(catalogue).module())


def _r_cppflags() -> list[str]:
    out = subprocess.run([RBIN, "CMD", "config", "--cppflags"],
                         capture_output=True, text=True, check=True).stdout
    return out.split()


def _build_lib(out: Path, with_r: bool = True) -> Path:
    """Compile the Fortran fixtures and the generated C `.Call` shims into one
    `libfixtures.so` -- mirroring the bundled production build, where fpm links the R shims
    into libtensor-omics.so. With `with_r=False` the shims are compiled with
    `-DNO_R_BINDING`, so they are empty objects and the library has no R entry points."""
    lib_dir = out / ("lib_r" if with_r else "lib_nor")
    lib_dir.mkdir(exist_ok=True)
    for name in ("fx_basics.F90", "fx_edges.F90"):
        shutil.copy(REPO_ROOT / FIXTURE_SRC / name, out / name)

    objects = []
    fflags = ["-cpp", "-I.", "-std=f2018", "-ffree-line-length-none", "-fPIC", f"-J{lib_dir}"]
    for source in (
        REPO_ROOT / "src/f42/tox_errors.F90",
        REPO_ROOT / "src/f42/tox_conversions.F90",
        REPO_ROOT / "src/f42/f42_safeguard.F90",
        out / "fx_basics.F90", out / "fx_edges.F90",
        out / "fx_basics_c.F90", out / "fx_edges_c.F90",
    ):
        obj = lib_dir / f"{Path(source).stem}.o"
        result = subprocess.run([GFORTRAN, *fflags, "-c", str(source), "-o", str(obj)],
                                cwd=REPO_ROOT, capture_output=True, text=True)
        assert result.returncode == 0, f"{source}:\n{result.stderr}"
        objects.append(str(obj))

    # the C `.Call` shims -> distinctly-named objects (fx_basics.c would otherwise collide
    # with the Fortran fx_basics.o). R symbols stay undefined, resolved at dyn.load.
    cflags = ["-fPIC", f"-I{out}", *_r_cppflags()] + ([] if with_r else ["-DNO_R_BINDING"])
    for source in sorted(out.glob("*.c")):
        obj = lib_dir / f"{source.stem}_shim.o"
        result = subprocess.run([GCC, *cflags, "-c", str(source), "-o", str(obj)],
                                cwd=REPO_ROOT, capture_output=True, text=True)
        assert result.returncode == 0, f"{source}:\n{result.stderr}"
        objects.append(str(obj))

    library = lib_dir / "libfixtures.so"
    result = subprocess.run([GFORTRAN, "-shared", "-o", str(library), *objects],
                            cwd=REPO_ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return library


@pytest.fixture(scope="session")
def built(tmp_path_factory):
    out = tmp_path_factory.mktemp("e2e_r")
    _generate(out)
    library = _build_lib(out, with_r=True)
    shutil.copy(library, out / "libfixtures.so")   # a stable path for run_r to dyn.load
    return out


def run_r(built: Path, body: str) -> str:
    """Load the bundled library, source the generated R, run `body`, return its stdout."""
    script = textwrap.dedent(f"""
        dyn.load("libfixtures.so")
        source("error_handling.R"); source("tox_validate.R")
        source("fx_basics.R"); source("fx_edges.R")
        {body}
    """)
    result = subprocess.run(
        [RSCRIPT, "-e", script], cwd=built, capture_output=True, text=True,
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

    def test_an_inout_scalar_is_capped_and_returned(self, built):
        # the modified scalar comes back from its <name>_v local, not the SEXP argument
        out = run_r(built, 'cat(fx_cap_value(15L), fx_cap_value(3L))')
        assert [int(x) for x in out.split()] == [10, 3]

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

    def test_the_padding_is_trimmed_rather_than_returned(self, built):
        # 'pos' arrives in an eight-wide blank-padded buffer and 'summary' in a sixteen-wide
        # one. `identical` does not trim, so this is the assertion `cat` cannot make.
        out = run_r(built, """
            r <- fx_labels(c(1, -2, 3))
            cat(identical(r$labels, c("pos", "nonpos", "pos")),
                identical(r$label, "summary"))
        """)
        assert out.split() == ["TRUE", "TRUE"]

    def test_a_string_vector_goes_in(self, built):
        out = run_r(built, 'cat(fx_count_matching(c("a", "b", "a"), "a"))')
        assert int(out) == 2

    def test_an_omitted_string_reaches_the_callee_as_absent(self, built):
        # C passes null, the wrapper leaves the view disassociated, and F2018 15.5.2.12
        # makes a disassociated pointer actual an absent optional dummy
        out = run_r(built, """
            cat(fx_optional_strings(),
                fx_optional_strings(tag = "x"),
                fx_optional_strings(extras = c("a", "b", "c")),
                fx_optional_strings(tag = "   "))
        """)
        assert out.split() == ["0", "1", "3", "0"]


class TestOptionals:
    def test_a_nullable_optional_may_be_given_or_omitted(self, built):
        out = run_r(built, """
            given <- fx_nullable(c(1,2,3), "ungrouped", c(1L,1L,2L))
            omitted <- fx_nullable(c(1,2,3), "ungrouped")
            cat(is.null(given), is.null(omitted))
        """)
        assert out.split() == ["TRUE", "TRUE"]


class TestAutoOutputFrom:
    def test_the_work_size_is_computed_by_calling_the_producer(self, built):
        # fx_cluster_expert calls fx_work_size itself; the caller passes only values
        out = run_r(built, 'cat(fx_cluster_expert(c(-1, 2, 3, -4, 5)))')
        assert int(out) == 3

    def test_the_producer_is_a_wrapper_in_its_own_right(self, built):
        out = run_r(built, 'cat(fx_work_size(4))')
        assert int(out) == 8


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


class TestNoRBinding:
    def test_no_r_binding_build_drops_the_call_shims(self, tmp_path_factory):
        """Built with NO_R_BINDING, the library keeps the Fortran C ABI but has no R
        `.Call` entry points -- the shims compiled to empty objects, needing no R headers."""
        out = tmp_path_factory.mktemp("e2e_r_nor")
        _generate(out)
        library = _build_lib(out, with_r=False)
        syms = subprocess.run(["nm", "-D", str(library)],
                              capture_output=True, text=True).stdout
        assert "fx_normalize_call" not in syms   # the R shim is gone
        assert "fx_normalize_c" in syms          # the Fortran C ABI stays


class TestPythonCanLoadTheBundledLibrary:
    def test_ctypes_loads_the_r_bundled_library(self, built):
        """The one library carries the R shims (undefined R symbols), but must still load
        into a non-R host under the *default* eager binding: every R symbol is marked weak,
        so it resolves to null (Python never calls the R code). This guards the fix that lets
        Python's ctypes load libtensor-omics.so once it also contains the R binding -- and,
        because the load is eager, that the weak set is complete (a strong R symbol would
        make this fail)."""
        import ctypes
        lib = ctypes.CDLL(str(built / "libfixtures.so"))   # default mode == eager
        assert hasattr(lib, "fx_normalize_c")   # a Fortran C-ABI symbol is callable
