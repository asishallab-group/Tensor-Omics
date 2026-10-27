"""The implementation layer through Python and R: generate everything, build one library, call it.

`test_end_to_end.py` and `test_end_to_end_r.py` drive the hand-written fixture exports, and
`test_end_to_end_fortran.py` stops at the Fortran tiers. This builds the implementation fixtures
the way fpm builds the real library -- the generated Fortran wrappers, the C wrappers and the
R `.Call` shims all in one `libtensor_omics.so` -- and calls them from the generated Python and
R. It is what shows that a declaration the wrappers can compile (an extent or a character
length named by a module constant) also works in the languages that cannot name the constant.

Python runs in a subprocess, so this library and the one `test_end_to_end.py` loads into the
test process never meet. Skipped without gfortran, gcc and R.
"""

import json
import os
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from codegen.config import Paths
from codegen.diagnostics import DiagnosticBag
from codegen.emit.errors_python import PythonErrorEmitter
from codegen.emit.errors_r import RErrorEmitter
from codegen.frontend.ford_frontend import FordFrontend
from codegen.generate import generate
from codegen.ir.errors import ErrorCatalogue

from conftest import REPO_ROOT
from test_end_to_end_fortran import DEPENDENCIES, FIXTURE_IMPL
from test_end_to_end_r import GCC, GFORTRAN, RBIN, RSCRIPT, _r_cppflags

pytestmark = pytest.mark.skipif(
    GFORTRAN is None or GCC is None or RSCRIPT is None or RBIN is None,
    reason="needs gfortran, gcc and R",
)

#: The fixture project has no `tox_errors`, so the error handling is generated from the real one
_NO_ERROR_MODULE = "no 'tox_errors' module"

#: The generated code is held to -std=f2018; the f42 sources are compiled the way fpm compiles
#: them, without it (see `test_end_to_end_fortran.DEPENDENCY_FLAGS`)
FORTRAN_FLAGS = ["-cpp", "-I.", "-std=f2018", "-ffree-line-length-none", "-fPIC", "-g",
                 "-fcheck=all"]
DEPENDENCY_FLAGS = [flag for flag in FORTRAN_FLAGS if flag != "-std=f2018"]


def _generate(out: Path) -> None:
    paths = Paths(
        root=REPO_ROOT, src_dir=REPO_ROOT / FIXTURE_IMPL, generated_dir=out / "fortran",
        c_binding_dir=out / "c", r_binding_dir=out / "shim",
        python_out_dir=out / "py" / "fx_impl", r_out_dir=out / "r",
    )
    result = generate(paths=paths, targets=("fortran", "c", "python", "r"),
                      library="build/libtensor_omics.so")
    unexpected = [d for d in result.diagnostics if _NO_ERROR_MODULE not in d.message]
    assert unexpected == [], "\n".join(d.render() for d in unexpected)
    for file in result.files:
        file.write()

    errors = DiagnosticBag()
    real = FordFrontend(Paths(root=REPO_ROOT, src_dir=Path("src/f42")), errors).parse()
    catalogue = ErrorCatalogue.from_module(
        real.project.module("tox_errors"), errors,
        real.project.constant_values(), arg_pos_factor=real.arg_pos_factor,
    )
    (out / "py" / "fx_impl" / "error_handling.py").write_text(PythonErrorEmitter(catalogue).module())
    (out / "r" / "error_handling.R").write_text(RErrorEmitter(catalogue).module())


def _compile_fortran(source: Path, out: Path, flags=FORTRAN_FLAGS) -> str:
    obj = out / "obj" / f"{source.stem}.o"
    result = subprocess.run(
        [GFORTRAN, *flags, f"-J{out / 'obj'}", "-c", str(source), "-o", str(obj)],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, f"{source}:\n{result.stderr}"
    return str(obj)


def _build(out: Path) -> Path:
    (out / "obj").mkdir()
    dependencies = [REPO_ROOT / d for d in DEPENDENCIES]
    objects = [_compile_fortran(source, out, DEPENDENCY_FLAGS) for source in dependencies]
    sources = sorted((REPO_ROOT / FIXTURE_IMPL).glob("*_impl.F90"))
    sources += sorted((out / "fortran").rglob("*.F90"))   # the wrappers before their C wrappers
    sources += sorted((out / "c").glob("*.F90"))
    objects += [_compile_fortran(source, out) for source in sources]

    # the R shims go into the same library, as fpm bundles them into libtensor_omics.so
    cflags = ["-fPIC", f"-I{out / 'shim'}", *_r_cppflags()]
    for source in sorted((out / "shim").glob("*.c")):
        obj = out / "obj" / f"{source.stem}_shim.o"
        result = subprocess.run([GCC, *cflags, "-c", str(source), "-o", str(obj)],
                                cwd=REPO_ROOT, capture_output=True, text=True)
        assert result.returncode == 0, f"{source}:\n{result.stderr}"
        objects.append(str(obj))

    library = out / "libtensor_omics.so"
    result = subprocess.run([GFORTRAN, "-shared", "-o", str(library), *objects],
                            cwd=REPO_ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return library


@pytest.fixture(scope="session")
def built(tmp_path_factory):
    out = tmp_path_factory.mktemp("e2e_impl")
    _generate(out)
    _build(out)
    return out


def run_python(built: Path, body: str):
    """Import the generated package in a fresh interpreter, run `body`, return what it prints as JSON."""
    script = "import json\nimport numpy as np\nimport fx_impl as tox\n" + textwrap.dedent(body)
    env = dict(os.environ, TENSOR_OMICS_LIBRARY=str(built / "libtensor_omics.so"),
               PYTHONPATH=str(built / "py"))
    result = subprocess.run([sys.executable, "-c", script], cwd=built, env=env,
                            capture_output=True, text=True)
    assert result.returncode == 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    return json.loads(result.stdout)


def run_r(built: Path, body: str) -> list[str]:
    modules = " ".join(f'source("r/{p.name}");' for p in sorted((built / "r").glob("fx_*.R")))
    script = textwrap.dedent(f"""
        dyn.load("libtensor_omics.so")
        source("r/error_handling.R"); source("r/tox_validate.R")
        {modules}
    """) + textwrap.dedent(body)
    result = subprocess.run([RSCRIPT, "-e", script], cwd=built, capture_output=True, text=True)
    assert result.returncode == 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    return result.stdout.split()


class TestConstantSizes:
    """A module constant in an extent or a length works in the languages that cannot name it.

    The shim allocated `Rf_allocVector(REALSXP, FX_N_AXES)`, which does not compile -- and
    with the shims bundled into the library, nothing built -- and Python raised NameError.
    """

    def test_an_input_and_an_output_extent_in_python(self, built):
        out = run_python(built, """
            points = np.asfortranarray([[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]])
            print(json.dumps(list(tox.axis_sums(points))))
        """)

        assert out == [6.0, 60.0]

    def test_an_input_and_an_output_extent_in_r(self, built):
        out = run_r(built, "cat(axis_sums(matrix(c(1, 10, 2, 20, 3, 30), 2)))")

        assert out == ["6", "60"]

    def test_a_constant_character_length_in_python(self, built):
        out = run_python(built, """
            result = tox.label_length("abc")
            print(json.dumps([int(result["n_chars"]), result["tag"]]))
        """)

        assert out == [3, "fx"]

    def test_a_constant_character_length_in_r(self, built):
        out = run_r(built, 'r <- label_length("abc"); cat(r$n_chars, r$tag)')

        assert out == ["3", "fx"]


class TestTheRestOfTheFixture:
    def test_an_assumed_shape_input_and_an_impure_elemental_in_python(self, built):
        out = run_python(built, """
            print(json.dumps([tox.sum_values([1.0, 2.0, 4.0]), tox.note_value(2.5)]))
        """)

        assert out == [7.0, 2.5]

    def test_an_assumed_shape_input_and_an_impure_elemental_in_r(self, built):
        out = run_r(built, "cat(sum_values(c(1, 2, 4)), note_value(2.5))")

        assert out == ["7", "2.5"]


class TestAPrologueBesideAnOwnIerr:
    """#221 B2: the prologue's own `fraction` comes before the implementation's own `ierr`.

    It used to be appended after it, so `fx_clip_c` took `..., ierr, fraction` and the
    languages, which pass by position into C, would have handed one for the other.
    """

    def test_python_clips_at_the_level_the_prologue_derives(self, built):
        out = run_python(built, "print(json.dumps(list(tox.fx_clip([1.0, 2.0, 4.0], 0.5))))")

        assert out == [1.0, 2.0, 2.0]

    def test_r_clips_at_the_level_the_prologue_derives(self, built):
        assert run_r(built, "cat(fx_clip(c(1, 2, 4), 0.5))") == ["1", "2", "2"]


class TestAnOptionalModeWithoutADefault:
    """#221 B10: an optional mode with no default, omitted, is absent to the implementation.

    The C wrapper took `c_char_as_view` of the null pointer it was handed, which killed Python
    with SIGSEGV and R with "caught segfault" -- and had it survived, it would have passed its
    local integer on, so the implementation would have seen the mode as present.
    """

    def test_python_omits_it(self, built):
        out = run_python(built, "print(json.dumps(list(tox.fx_scale_by_method([1.0, 4.0]))))")

        assert out == [1.0, 4.0]

    def test_python_passes_it(self, built):
        out = run_python(built, """
            print(json.dumps([list(tox.fx_scale_by_method([1.0, 4.0], method)) for method in
                              ("double", "halve")]))
        """)

        assert out == [[2.0, 8.0], [0.5, 2.0]]

    def test_r_omits_it(self, built):
        assert run_r(built, "cat(fx_scale_by_method(c(1, 4)))") == ["1", "4"]

    def test_r_passes_it(self, built):
        out = run_r(built, 'cat(fx_scale_by_method(c(1, 4), "double"), '
                           'fx_scale_by_method(c(1, 4), "halve"))')

        assert out == ["2", "8", "0.5", "2"]

