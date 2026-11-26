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
from codegen.generate import _clean, generate
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


def _generate(out: Path, impl_dir: Path = REPO_ROOT / FIXTURE_IMPL) -> None:
    """Generate every target from `impl_dir` into `out`, cleaning first as a real run does."""
    paths = Paths(
        root=REPO_ROOT, src_dir=impl_dir, generated_dir=out / "fortran",
        c_binding_dir=out / "c", r_binding_dir=out / "shim",
        python_out_dir=out / "py" / "fx_impl", r_out_dir=out / "r",
    )
    targets = ("fortran", "c", "python", "r")
    result = generate(paths=paths, targets=targets, library="build/libtensor_omics.so")
    unexpected = [d for d in result.diagnostics if _NO_ERROR_MODULE not in d.message]
    assert unexpected == [], "\n".join(d.render() for d in unexpected)
    # `generate_and_write` would clean and write, but writes nothing here: the missing
    # tox_errors is an error, and its stand-in is written below
    _clean(targets, paths)
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


def _build(out: Path, impl_dir: Path = REPO_ROOT / FIXTURE_IMPL) -> Path:
    (out / "obj").mkdir()
    dependencies = [REPO_ROOT / d for d in DEPENDENCIES]
    objects = [_compile_fortran(source, out, DEPENDENCY_FLAGS) for source in dependencies]
    sources = sorted(impl_dir.glob("*_impl.F90"))
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


class TestAnExtentTheBindingsCompute:
    """#221 B11: `(n - 10)/3 + 5` is 3 for n = 3 in Fortran, which truncates -7/3 to -2.

    Python's `//` floors it to -3, and an output allocated with that is one element short: the
    Fortran then writes past it. The shim's C truncates as Fortran does.
    """

    def test_python_truncates_toward_zero(self, built):
        out = run_python(built, """
            print(json.dumps([list(map(int, tox.truncated_count(n))) for n in (3, 13)]))
        """)

        assert out == [[1, 2, 3], [1, 2, 3, 4, 5, 6]]

    def test_r_truncates_toward_zero(self, built):
        out = run_r(built, "cat(length(truncated_count(3L)), length(truncated_count(13L)))")

        assert out == ["3", "6"]


#: A refusal by the binding, as Python reports it: its class, whether it carries no code, the
#: argument, and the message
_PY_REFUSED = """
            def refused(call):
                try:
                    call()
                except tox.error_handling.ToxInputError as error:
                    return [type(error).__name__, error.code, error.argument, str(error)]
                return None
"""  # indented as the bodies below are, which `run_python` dedents together with it

#: The same for R: the first class, the code (NA for the binding's own) and the argument
_R_REFUSED = """
            refused <- function(expr) {
              e <- tryCatch({ expr; NULL }, tox_input_error = function(e) e)
              if (is.null(e)) return(cat("none\\n"))
              cat(class(e)[1], if (is.na(e$code)) "NA" else e$code, e$argument, "\\n")
            }
"""


class TestEveryInputAxisIsCompared:
    """#221 A-shape: an input axis that is wrong is refused before Fortran reads past it.

    On main each of these was garbage: Fortran's explicit-shape dummies take an address, and an
    axis that is not a bare extent shared with another input was compared with nothing.
    """

    def test_a_constant_axis_in_python(self, built):
        out = run_python(built, _PY_REFUSED + """
            print(json.dumps([refused(lambda: tox.pair_sums(np.ones((3, 1)))),
                              list(tox.pair_sums(np.ones((3, 2))))]))
        """)

        name, code, argument, message = out[0]
        assert (name, code, argument) == ("ToxArgumentError", None, "pairs")
        assert "along axis 1, but its declared shape fixes it at 2" in message
        assert out[1] == [2.0, 2.0, 2.0]

    def test_a_module_constant_axis_in_python(self, built):
        # WP1 review round 2: `axis_sums(points(FX_N_AXES, n))` with a 3-row points over-read
        out = run_python(built, _PY_REFUSED + """
            print(json.dumps(refused(lambda: tox.axis_sums(np.ones((3, 2))))))
        """)

        assert out[:3] == ["ToxArgumentError", None, "points"]
        assert "is declared with shape (FX_N_AXES, n_points)" in out[3]

    def test_an_optional_array_is_compared_when_given(self, built):
        out = run_python(built, _PY_REFUSED + """
            print(json.dumps([refused(lambda: tox.weighted_total([1.0, 2.0, 3.0], weights=[1.0])),
                              tox.weighted_total([1.0, 2.0, 3.0]),
                              tox.weighted_total([1.0, 2.0, 3.0], weights=[1.0, 0.0, 2.0])]))
        """)

        assert out[0][:3] == ["ToxArgumentError", None, "weights"]
        assert "'values' implies n == 3" in out[0][3]
        assert out[1:] == [6.0, 7.0]

    def test_an_input_sized_by_a_mask_count(self, built):
        out = run_python(built, _PY_REFUSED + """
            mask = [True, False, True]
            print(json.dumps([refused(lambda: tox.selected_total([1.0, 2.0, 3.0], mask, [1.0, 1.0, 1.0])),
                              tox.selected_total([1.0, 2.0, 3.0], mask, [1.0, 2.0])]))
        """)

        assert out[0][:3] == ["ToxArgumentError", None, "weights"]
        assert "'values_mask' implies n_selected_values == 2" in out[0][3]
        assert out[1] == 7.0

    def test_an_axis_the_declaration_computes_is_compared_in_c(self, built):
        # C is given the actual extent and compares it with `n + 1`: ERR_SIZE_MISMATCH, naming x
        out = run_python(built, _PY_REFUSED + """
            x = np.array([1.0, 2.0, 3.0, 99.0])
            print(json.dumps([refused(lambda: tox.shifted_sum([1.0, 1.0, 1.0], x[:3])),
                              tox.shifted_sum([1.0, 1.0, 1.0], x),
                              tox.error_handling.ERR_SIZE_MISMATCH]))
        """)

        name, code, argument, message = out[0]
        assert (name, code, argument) == ("ToxInputError", out[2], "x")
        # a code carries no text: the binding adds what the declaration says x must be
        assert message.endswith("(argument 'x'); 'x' is declared with shape (n + 1)")
        assert out[1] == 108.0

    def test_an_axis_reading_another_array_s_contents(self, built):
        out = run_python(built, _PY_REFUSED + """
            print(json.dumps([refused(lambda: tox.grouped_sums([2, 2], [1.0, 2.0])),
                              list(tox.grouped_sums([2, 1], [1.0, 2.0, 4.0]))]))
        """)

        assert out[0][:3] == ["ToxInputError", 206, "values"]
        assert out[0][3].endswith("'values' is declared with shape (sum(group_sizes))")
        assert out[1] == [3.0, 4.0]

    def test_an_optional_axis_the_declaration_computes(self, built):
        out = run_python(built, _PY_REFUSED + """
            print(json.dumps([tox.optional_shift([1.0, 1.0]),
                              tox.optional_shift([1.0, 1.0], x=[1.0, 1.0, 1.0]),
                              refused(lambda: tox.optional_shift([1.0, 1.0], x=[1.0, 1.0]))]))
        """)

        assert out[:2] == [2.0, 5.0]
        assert out[2][:3] == ["ToxInputError", 206, "x"]

    def test_r_refuses_each_the_same_way(self, built):
        out = run_r(built, _R_REFUSED + """
            refused(pair_sums(matrix(1, 3, 1)))
            refused(axis_sums(matrix(1, 3, 2)))
            refused(weighted_total(c(1, 2, 3), weights = 1))
            refused(selected_total(c(1, 2, 3), c(TRUE, FALSE, TRUE), c(1, 1, 1)))
            refused(shifted_sum(c(1, 1, 1), c(1, 2, 3)))
            refused(grouped_sums(c(2L, 2L), c(1, 2)))
            refused(optional_shift(c(1, 1), x = c(1, 1)))
            cat(pair_sums(matrix(1, 3, 2)), weighted_total(c(1, 2, 3)),
                selected_total(c(1, 2, 3), c(TRUE, FALSE, TRUE), c(1, 2)),
                shifted_sum(c(1, 1, 1), c(1, 2, 3, 99)), grouped_sums(c(2L, 1L), c(1, 2, 4)),
                optional_shift(c(1, 1)), optional_shift(c(1, 1), x = c(1, 1, 1)), "\\n")
        """)

        assert out == [
            "tox_shape_error", "NA", "pairs",
            "tox_shape_error", "NA", "points",
            "tox_shape_error", "NA", "weights",
            "tox_shape_error", "NA", "weights",
            "tox_input_error", "206", "x",
            "tox_input_error", "206", "values",
            "tox_input_error", "206", "x",
            "2", "2", "2", "6", "7", "108", "3", "4", "2", "5",
        ]

    def test_r_states_the_declared_shape_of_a_size_mismatch(self, built):
        out = run_r(built, """
            cat(tryCatch(grouped_sums(c(2L, 2L), c(1, 2)),
                         tox_input_error = function(e) conditionMessage(e)), "\\n")
        """)

        assert " ".join(out).endswith("'values' is declared with shape (sum(group_sizes))")

    def test_a_null_array_an_extent_reads_is_reported_not_dereferenced(self, built):
        # the C wrapper used to declare `values(sum(group_sizes))`, evaluated on entry -- before
        # the null check of `group_sizes` (attack-shape F3): now ERR_POINTER_NULL comes back
        out = run_python(built, """
            import ctypes
            lib = ctypes.CDLL(str(tox.library.load_library()._name))
            n, n_values, ierr = ctypes.c_int(2), ctypes.c_int(2), ctypes.c_int(0)
            values = (ctypes.c_double * 2)(1.0, 2.0)
            sums = (ctypes.c_double * 2)()
            lib.grouped_sums_c(ctypes.byref(n), None, values, ctypes.byref(n_values), sums,
                               ctypes.byref(ierr))
            print(json.dumps([ierr.value, tox.error_handling.ERR_POINTER_NULL]))
        """)

        assert out[0] == out[1]


class TestAnExtentThatMustEqual:
    """#203: `n_rows` must equal `sum(group_sizes)`; the caller passes neither, the binding reads it.

    A mismatch is ERR_SIZE_MISMATCH at the extent, blamed on `values` and saying what it had to
    equal; an empty input stays ERR_EMPTY_INPUT; a bad count is the count's own error.
    """

    def test_python(self, built):
        out = run_python(built, _PY_REFUSED + """
            print(json.dumps([
                refused(lambda: tox.grouped_rows([2, 2], [1.0, 2.0, 3.0])),
                refused(lambda: tox.grouped_rows([5, -3], [1.0, 2.0])),
                refused(lambda: tox.grouped_rows([], [1.0, 2.0])),
                list(tox.grouped_rows([2, 1], [1.0, 2.0, 4.0])),
            ]))
        """)

        mismatch, negative, empty, sums = out
        assert mismatch == ["ToxInputError", 206, "values",
                            "Array size mismatch (argument 'values', via 'n_rows'); "
                            "'n_rows' must equal sum(group_sizes)"]
        assert negative[:3] == ["ToxInputError", 201, "group_sizes"]
        assert empty[:3] == ["ToxInputError", 202, "group_sizes"]
        assert sums == [3.0, 4.0]

    def test_r(self, built):
        out = run_r(built, _R_REFUSED + """
            refused(grouped_rows(c(2L, 2L), c(1, 2, 3)))
            refused(grouped_rows(c(5L, -3L), c(1, 2)))
            refused(grouped_rows(integer(), c(1, 2)))
            cat(grouped_rows(c(2L, 1L), c(1, 2, 4)), "\\n")
            cat(tryCatch(grouped_rows(c(2L, 2L), c(1, 2, 3)),
                         tox_input_error = function(e) conditionMessage(e)), "\\n")
        """)

        assert out[:9] == [
            "tox_input_error", "206", "values",
            "tox_input_error", "201", "group_sizes",
            "tox_input_error", "202", "group_sizes",
        ]
        assert out[9:11] == ["3", "4"]
        assert " ".join(out[11:]) == ("Array size mismatch (argument 'values', via 'n_rows'); "
                                      "'n_rows' must equal sum(group_sizes)")


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
        out = run_python(built, """
            result = tox.fx_scale_by_method([1.0, 4.0])
            print(json.dumps([list(result["scaled"]), int(result["method_given"])]))
        """)

        assert out == [[1.0, 4.0], 0]

    def test_python_passes_it(self, built):
        out = run_python(built, """
            results = [tox.fx_scale_by_method([1.0, 4.0], method) for method in ("double", "halve")]
            print(json.dumps([[list(r["scaled"]), int(r["method_given"])] for r in results]))
        """)

        assert out == [[[2.0, 8.0], 1], [[0.5, 2.0], 1]]

    def test_r_omits_it(self, built):
        out = run_r(built, "r <- fx_scale_by_method(c(1, 4)); cat(r$scaled, r$method_given)")

        assert out == ["1", "4", "0"]

    def test_r_passes_it(self, built):
        out = run_r(built, """
            for (method in c("double", "halve")) {
                r <- fx_scale_by_method(c(1, 4), method); cat(r$scaled, r$method_given, "")
            }
        """)

        assert out == ["2", "8", "1", "0.5", "2", "1"]


class TestFunctions:
    """#221 B1: a function implementation, its result published right before `ierr`.

    The wrapper called it as a subroutine, which does not compile, and nothing said so.
    """

    def test_python_returns_each_result(self, built):
        out = run_python(built, """
            print(json.dumps([
                tox.fx_total([1.0, 2.0, 5.0]), list(tox.fx_doubled([1.0, 2.0])),
                bool(tox.fx_all_positive([1.0, -2.0])), tox.fx_halve(5.0),
                tox.fx_share([1.0, 2.0, 5.0]), tox.fx_smallest([3.0, 1.0, 2.0]),
            ]))
        """)

        assert out == [8.0, [2.0, 4.0], False, 2.5, 0.125, 1.0]

    def test_r_returns_each_result(self, built):
        out = run_r(built, """
            cat(fx_total(c(1, 2, 5)), fx_doubled(c(1, 2)), fx_all_positive(c(1, -2)),
                fx_halve(5), fx_share(c(1, 2, 5)), fx_smallest(c(3, 1, 2)))
        """)

        assert out == ["8", "2", "4", "FALSE", "2.5", "0.125", "1"]

    def test_the_implementation_s_own_error_reaches_python_and_r(self, built):
        # ERR_DIVISION_BY_ZERO (210), raised by the implementation itself
        out = run_python(built, """
            try:
                tox.fx_share([0.0, 0.0])
            except tox.error_handling.ToxError as error:
                print(json.dumps(error.code))
        """)
        assert out == 210
        assert run_r(built, """
            r <- tryCatch(fx_share(c(0, 0)), tox_error = function(e) e$code); cat(r)
        """) == ["210"]

    def test_c_takes_the_result_right_before_ierr(self, built):
        # by position, as a C caller passes it: fx_share_c(n_values, values, share, ierr)
        out = run_python(built, """
            import ctypes, os
            library = ctypes.CDLL(os.environ["TENSOR_OMICS_LIBRARY"])
            n_values = ctypes.c_int(3)
            values = (ctypes.c_double * 3)(1.0, 2.0, 5.0)
            share, ierr = ctypes.c_double(-1.0), ctypes.c_int(-1)
            library.fx_share_c(ctypes.byref(n_values), values, ctypes.byref(share),
                               ctypes.byref(ierr))
            print(json.dumps([share.value, ierr.value]))
        """)

        assert out == [0.125, 0]

    def test_a_complex_output_comes_back_whole_in_python_c_and_r(self, built):
        # #221 WP8: Python boxed it in one c_double, so the imaginary part was lost and Fortran
        # wrote 8 bytes past the box; C and R had it right already
        out = run_python(built, """
            z, w = tox.fx_complex(1.0, -2.0), tox.fx_make_complex(1.0, -2.0)
            print(json.dumps([type(z).__name__, z.real, z.imag, type(w).__name__, w.real, w.imag]))
        """)
        assert out == ["complex", 1.0, -2.0, "complex", 1.0, -2.0]

        # by position, as a C caller passes it: two doubles, real part first
        out = run_python(built, """
            import ctypes, os
            library = ctypes.CDLL(os.environ["TENSOR_OMICS_LIBRARY"])
            re, im = ctypes.c_double(1.0), ctypes.c_double(-2.0)
            z, w, ierr = (ctypes.c_double * 2)(), (ctypes.c_double * 2)(), ctypes.c_int(-1)
            library.fx_complex_c(ctypes.byref(re), ctypes.byref(im), z, ctypes.byref(ierr))
            library.fx_make_complex_c(ctypes.byref(re), ctypes.byref(im), w, ctypes.byref(ierr))
            print(json.dumps([list(z), list(w), ierr.value]))
        """)
        assert out == [[1.0, -2.0], [1.0, -2.0], 0]

        assert run_r(built, """
            z <- fx_complex(1, -2); w <- fx_make_complex(1, -2)
            cat(class(z), format(z), class(w), format(w))
        """) == ["complex", "1-2i", "complex", "1-2i"]

    def test_a_complex_input_crosses_whole_as_in_c_and_r(self, built):
        # #221 WP8 fix 1: Python boxed a complex input in one c_double, so a complex raised a bare
        # TypeError, and a float's imaginary part was read from whatever lay past the box; an
        # intent(inout) one was passed as the caller's own object
        out = run_python(built, """
            print(json.dumps([
                tox.fx_parts(3 - 4j), tox.fx_parts(3.0), tox.fx_parts(np.complex64(1 + 2j)),
                repr(tox.fx_offset(1.0, -2.0)), repr(tox.fx_offset(1.0, -2.0, 1j)),
                repr(tox.fx_offset(1.0, -2.0, offset=2)), repr(tox.fx_twice(3 - 4j)),
                repr(tox.fx_twice(5)),
            ]))
        """)
        assert out == [{"re": 3.0, "im": -4.0}, {"re": 3.0, "im": 0.0}, {"re": 1.0, "im": 2.0},
                       "(1-2j)", "(1-1j)", "(3-2j)", "(6-8j)", "(10+0j)"]

        # refused by name, as every contract refusal is
        out = run_python(built, """
            from fx_impl.error_handling import ToxArgumentError
            caught = []
            for call in (lambda: tox.fx_parts("3"), lambda: tox.fx_parts(True),
                         lambda: tox.fx_parts(None), lambda: tox.fx_twice([1j]),
                         lambda: tox.fx_offset(1.0, 2.0, "1j"), lambda: tox.fx_parts(2**53 + 1)):
                try:
                    call()
                    caught.append(None)
                except ToxArgumentError as error:
                    caught.append(error.argument)
            print(json.dumps(caught))
        """)
        assert out == ["z", "z", "z", "z", "offset", "z"]

        # C by position, and R: the same numbers
        out = run_python(built, """
            import ctypes, os
            library = ctypes.CDLL(os.environ["TENSOR_OMICS_LIBRARY"])
            z, re, im, ierr = (ctypes.c_double * 2)(3.0, -4.0), ctypes.c_double(), ctypes.c_double(), ctypes.c_int(-1)
            library.fx_parts_c(z, ctypes.byref(re), ctypes.byref(im), ctypes.byref(ierr))
            library.fx_twice_c(z, ctypes.byref(ierr))
            print(json.dumps([re.value, im.value, list(z), ierr.value]))
        """)
        assert out == [3.0, -4.0, [6.0, -8.0], 0]
        assert run_r(built, """
            p <- fx_parts(3-4i); cat(p$re, p$im, format(fx_offset(1, -2)), format(fx_offset(1, -2, 1i)),
                                    format(fx_twice(3-4i)))
        """) == ["3", "-4", "1-2i", "1-1i", "6-8i"]


#: Where `_generate` writes the Fortran wrappers, the C wrappers, the R shims, Python and R
_TARGET_DIRS = ("fortran", "c", "shim", "py", "r")


@pytest.fixture(scope="module")
def regenerated(tmp_path_factory):
    """The fixture implementations generated, one of them deleted, and generated again."""
    impl_dir = tmp_path_factory.mktemp("impl")
    for source in (REPO_ROOT / FIXTURE_IMPL).glob("*_impl.F90"):
        shutil.copy(source, impl_dir)
    out = tmp_path_factory.mktemp("e2e_deleted")
    (out / "r").mkdir()
    (out / "r" / "DESCRIPTION").write_text("Package: hand-written\n")
    _generate(out, impl_dir)
    assert (out / "fortran" / "fx_modes.F90").is_file()   # the premise: it was generated

    (impl_dir / "fx_modes_impl.F90").unlink()
    _generate(out, impl_dir)
    return out, impl_dir


class TestADeletedImplementationLeavesNothingBehind:
    """#221 B21: regenerating after deleting an implementation removes everything it made.

    The fortran target cleaned only the wrappers of the implementations that still existed, so
    a deleted one's wrapper stayed in `src/generated` -- where fpm compiles everything -- and
    broke every build after with "Unable to find source for module dependency".
    """

    def test_no_file_of_it_remains_in_any_target(self, regenerated):
        out, _ = regenerated

        generated = [p for target in _TARGET_DIRS for p in (out / target).rglob("*") if p.is_file()]
        left = [p.relative_to(out).as_posix() for p in generated
                if "fx_modes" in p.name or "fx_scale_by_method" in p.read_text()]
        assert generated and left == []

    def test_a_hand_written_file_beside_the_r_wrappers_stays(self, regenerated):
        out, _ = regenerated

        assert (out / "r" / "DESCRIPTION").read_text() == "Package: hand-written\n"

    def test_the_library_still_builds(self, regenerated):
        out, impl_dir = regenerated

        assert _build(out, impl_dir).is_file()


#: A producer's error as Python reports it: its class, its code's name, the argument and the
#: message. None where the call returned, and then what it returned.
_PY_PRODUCER_ERROR = """
            def raised(call):
                try:
                    return ["returned", call()]
                except tox.error_handling.ToxError as error:
                    return [type(error).__name__, error.name, error.argument, str(error)]
"""  # indented as the bodies below are, which `run_python` dedents together with it

def r_raised(built: Path, call: str) -> list[str]:
    """What `call` raises in R: the first class, the code's name, the argument, the message.

    `["returned", <the value, space-separated>]` where it returns instead.
    """
    out = run_r(built, f"""
        e <- tryCatch({{ v <- {call}; NULL }}, tox_error = function(e) e)
        if (is.null(e)) cat("returned", paste(v, collapse = " "), sep = "|") else
            cat(class(e)[1], e$name, e$argument, conditionMessage(e), sep = "|")
    """)
    return " ".join(out).split("|")


class TestAGeneratedWrapperAsProducer:
    """`DM_OUTPUT_FROM` naming an implementation's generated wrapper (`fx_producers_impl`).

    The producer validates its inputs, so its error carries the position of its own argument. A
    consumer that passed it through unchanged named the consumer's argument at that position --
    `x` for a bad `scale` -- and a second producer's `set_ok` swallowed the first one's error.
    """

    def test_python_names_the_argument_a_fortran_called_producer_was_given(self, built):
        # fx_work_total's wrapper calls fx_scaled_size itself, `scale` as its `factor`
        out = run_python(built, _PY_PRODUCER_ERROR + """
            print(json.dumps([raised(lambda: tox.fx_work_total([1.0, 2.0], 0)),
                              raised(lambda: tox.fx_work_total([1.0, 2.0], 3))]))
        """)

        assert out[0][:3] == ["ToxInputError", "ERR_INVALID_INPUT", "scale"]
        assert out[0][3].endswith("(argument 'scale')")
        assert out[1] == ["returned", 9.0]

    def test_r_names_the_argument_a_fortran_called_producer_was_given(self, built):
        assert r_raised(built, "fx_work_total(c(1, 2), 0L)") == [
            "tox_input_error", "ERR_INVALID_INPUT", "scale",
            "invalid input arguments (argument 'scale')"]
        assert r_raised(built, "fx_work_total(c(1, 2), 3L)") == ["returned", "9"]

    def test_the_first_of_two_producers_to_fail_is_the_error(self, built):
        # fx_capped_size refuses 7, then fx_scaled_size would succeed and reset ierr
        out = run_python(built, _PY_PRODUCER_ERROR + """
            print(json.dumps([raised(lambda: tox.fx_two_producers([1.0], 7, 2)),
                              raised(lambda: tox.fx_two_producers([1.0], 4, 2))]))
        """)

        assert out[0][:3] == ["ToxInputError", "ERR_INVALID_INPUT", "n_count"]
        assert out[1] == ["returned", 1.0 + 5 + 2]

    def test_a_producer_in_the_consumer_s_own_module(self, built):
        out = run_python(built, _PY_PRODUCER_ERROR + """
            print(json.dumps([raised(lambda: tox.fx_local_total([1.0, 2.0], 2)),
                              raised(lambda: tox.fx_local_total([1.0, 2.0], 0))]))
        """)

        assert out[0] == ["returned", 3.0 + 4]
        assert out[1][2] == "scale"

    def test_python_names_the_argument_a_binding_called_producer_was_given(self, built):
        # fx_fill's result is sized by fx_scaled_size, which Python calls before fx_fill
        out = run_python(built, _PY_PRODUCER_ERROR + """
            print(json.dumps([raised(lambda: tox.fx_fill([1.0, 2.0], 0)),
                              raised(lambda: tox.fx_fill([], 2)),
                              raised(lambda: list(tox.fx_fill([1.0, 2.0], 2)))]))
        """)

        assert out[0][:3] == ["ToxInputError", "ERR_INVALID_INPUT", "scale"]
        assert out[0][3].endswith("(argument 'scale')")
        # the producer's `n_items` is the consumer's, which the caller wrote as `x`: the
        # producer refuses it (DM_MIN) before fx_fill itself is called
        assert out[1][:3] == ["ToxInputError", "ERR_INVALID_INPUT", "x"]
        assert out[1][3].endswith("(argument 'x', via 'n_items')")
        assert out[2] == ["returned", [2.0, 2.0, 2.0, 2.0]]

    def test_r_names_the_argument_a_binding_called_producer_was_given(self, built):
        assert r_raised(built, "fx_fill(c(1, 2), 0L)") == [
            "tox_input_error", "ERR_INVALID_INPUT", "scale",
            "invalid input arguments (argument 'scale')"]
        assert r_raised(built, "fx_fill(numeric(0), 2L)")[2:] == [
            "x", "invalid input arguments (argument 'x', via 'n_items')"]
        assert r_raised(built, "fx_fill(c(1, 2), 2L)") == ["returned", "2 2 2 2"]

    def test_r_reports_the_first_of_two_producers_to_fail(self, built):
        assert r_raised(built, "fx_two_producers(1, 7L, 2L)")[:3] == [
            "tox_input_error", "ERR_INVALID_INPUT", "n_count"]

    def test_an_extent_the_producer_calls_by_another_name(self, built):
        # fx_extent_size(n, x(n)) for fx_extent_total(m, x(m)): the table row `n | m`
        assert run_python(built, """
            print(json.dumps(tox.fx_extent_total([1.0, 2.0, 3.0])))
        """) == 6.0 + 4
        assert run_r(built, "cat(fx_extent_total(c(1, 2, 3)))") == ["10"]

    def test_an_output_nothing_takes(self, built):
        # fx_two_sizes also outputs `second_size`, which fx_first_of_two has no use for
        assert run_python(built, """
            print(json.dumps(tox.fx_first_of_two([1.0, 2.0])))
        """) == 3.0 + 2

    def test_outputs_nothing_takes_of_other_kinds_and_of_one_name(self, built):
        # fx_two_sizes' integer `second_size` and fx_flagged_size's real one, and its c_bool flag:
        # each dropped into a local of its own type, whose kind the module imports
        assert run_python(built, """
            print(json.dumps(tox.fx_both_firsts([1.0, 2.0])))
        """) == 3.0 + 2 + 3

    def test_one_producer_given_two_inputs_is_called_twice(self, built):
        # fx_scaled_size sizes tmp_a with scale_a and tmp_b with scale_b: two calls, two sizes
        assert run_python(built, """
            print(json.dumps(tox.fx_two_scales([1.0, 2.0], 1, 5)))
        """) == 3.0 + 2 + 1000 * 10
        assert run_r(built, "cat(fx_two_scales(c(1, 2), 1L, 5L))") == ["10005"]

    def test_a_shape_the_bindings_read_off_an_array(self, built):
        # fx_vector_count takes `arr_shape`, which neither language asks of the caller
        assert run_python(built, """
            print(json.dumps([list(map(int, tox.fx_status_per_vector(np.ones((2, 3, 4))))),
                              list(map(int, tox.fx_status_per_vector(np.ones(3))))]))
        """) == [[2] * 12, [3]]
        assert run_r(built, """
            cat(length(fx_status_per_vector(array(1, c(2, 3, 4)))),
                fx_status_per_vector(c(1, 2, 3)))
        """) == ["12", "3"]


class TestJustInfoLinksThePublishedWrapper:
    """`DM_OUTPUT_FROM(..., JUST_INFO)` names the producer as `AUTO` does, by its wrapper.

    So the link every layer renders is to the entry point a reader can call -- the generated
    `fx_scaled_size` in `fx_producers` -- not to the unvalidated implementation behind it. The
    link is checked like any (`ir.validate.check_doc_links`), which `_generate` asserts passed.
    """

    LINK = "produced by [[fx_producers(module):fx_scaled_size]]."

    def info_lines(self, path):
        return [line.strip() for line in path.read_text().splitlines() if "It is recommended" in line]

    def test_the_fortran_and_c_wrappers_link_it_for_ford(self, built):
        for path in (built / "fortran" / "fx_sized.F90", built / "c" / "fx_sized_c.F90"):
            lines = self.info_lines(path)
            assert lines and all(line.endswith(self.LINK) for line in lines), (path, lines)

    def test_python_and_r_link_their_own_binding(self, built):
        python = self.info_lines(built / "py" / "fx_impl" / "fx_sized.py")
        r = self.info_lines(built / "r" / "fx_sized.R")
        assert python == ["It is recommended to compute this argument from the `scaled_size` "
                          "output produced by :func:`tensor_omics.fx_scaled_size`."]
        assert r == ["#'   It is recommended to compute this argument from the `scaled_size` "
                     "output produced by \\code{\\link{fx_scaled_size}}."]

    def test_the_caller_passes_the_size(self, built):
        assert run_python(built, """
            print(json.dumps(list(tox.fx_info([1.0], 3))))
        """) == [1.0, 1.0, 1.0]


#: A 2 x 3 matrix whose column norms are 5, 1 and 2.5: unsorted, and one of them not above one
_NORMS_X_PY = """
            x = np.asfortranarray([[3.0, 0.6, 1.5], [4.0, 0.8, 2.0]])
"""  # indented as the bodies below are, which `run_python` dedents together with it
_NORMS_X = "matrix(c(3, 4, 0.6, 0.8, 1.5, 2), 2)"
_NORMS_X_R = f"""
            x <- {_NORMS_X}
"""


class TestAnArrayAProducerFills:
    """An array a producer fills (`fx_norm_rank`, `fx_large_first`), through Python and R.

    The plain tier computes it in Fortran. The expert tier takes it from the caller, as Fortran's
    and C's do (E2) -- unless the same producer call also gives a size the binding allocates by,
    which it then computes whole.
    """

    def test_the_plain_call_in_python(self, built):
        out = run_python(built, _NORMS_X_PY + """
            result = tox.fx_norm_rank(x, 2)
            print(json.dumps([list(result["ranked"]), result["total"]]))
        """)

        assert out == [[1.0, 2.5, 5.0], 6.0]

    def test_the_plain_call_in_r(self, built):
        out = run_r(built, _NORMS_X_R + """
            r <- fx_norm_rank(x, 2L); cat(r$ranked, r$total)
        """)

        assert out == ["1", "2.5", "5", "6"]

    def test_the_expert_takes_the_caller_s_norms_in_python(self, built):
        # sentinels no producer would compute, read-only as a result handed back would be:
        # they reach the implementation unchanged
        out = run_python(built, _NORMS_X_PY + """
            norms = np.array([9.0, 7.0, 8.0]); norms.flags.writeable = False
            result = tox.fx_norm_rank_expert(x, 2, norms, [2, 3, 1])
            print(json.dumps(list(result["ranked"])))
        """)

        assert out == [7.0, 8.0, 9.0]

    def test_the_expert_takes_the_caller_s_norms_in_r(self, built):
        out = run_r(built, _NORMS_X_R + """
            cat(fx_norm_rank_expert(x, 2L, c(9, 7, 8), c(2L, 3L, 1L))$ranked)
        """)

        assert out == ["7", "8", "9"]

    def test_the_expert_signature(self, built):
        out = run_python(built, """
            import inspect
            print(json.dumps([list(inspect.signature(f).parameters) for f in (
                tox.fx_norm_rank_expert, tox.fx_large_first_expert, tox.fx_rounded_total_expert,
                tox.fx_bounded_total_expert)]))
        """)

        assert out == [["x", "scale", "norms", "norms_perm"], ["scale", "x", "norms_perm"],
                       ["x", "rounded"], ["x", "scale", "limit"]]

    def test_the_docstring_names_what_the_expert_takes(self, built):
        # R5: the Python parameters, not the C arguments (`work_size`, `tmp_buf`)
        out = run_python(built, """
            print(json.dumps(tox.fx_norm_rank_expert.__doc__))
        """)

        assert "you supply `norms` and `norms_perm` yourself." in out
        assert ("`fx_norm_rank` computes `norms` with `fx_column_norms` and seeds `norms_perm` "
                "and sorts it by `norms`.") in out

    def test_nan_in_the_caller_s_norms_is_refused_naming_them_in_python(self, built):
        out = run_python(built, _PY_PRODUCER_ERROR + _NORMS_X_PY + """
            print(json.dumps(raised(lambda: tox.fx_norm_rank_expert(
                x, 2, [9.0, float("nan"), 8.0], [2, 3, 1]))))
        """)

        assert out[2] == "norms", out

    def test_na_in_the_caller_s_norms_is_refused_naming_them_in_r(self, built):
        out = r_raised(built, f"fx_norm_rank_expert({_NORMS_X}, 2L, c(9, NA, 8), c(2L, 3L, 1L))")

        assert out[2] == "norms", out

    def test_a_wrong_length_is_refused(self, built):
        out = run_python(built, _PY_REFUSED + _NORMS_X_PY + """
            print(json.dumps(refused(lambda: tox.fx_norm_rank_expert(x, 2, [9.0, 7.0], [2, 1]))))
        """)

        assert out[:3] == ["ToxArgumentError", None, "norms"]

    def test_a_producer_error_names_the_caller_s_argument(self, built):
        # fx_column_norms takes at most three rows; four reach it from the plain tier's Fortran
        out = run_python(built, _PY_PRODUCER_ERROR + """
            print(json.dumps(raised(lambda: tox.fx_norm_rank(np.ones((4, 3), order="F"), 2))))
        """)

        assert out[:3] == ["ToxInputError", "ERR_INVALID_INPUT", "x"], out
        assert out[3].endswith("(argument 'x', via 'n_dims')"), out
        assert r_raised(built, "fx_norm_rank(matrix(1, 4, 3), 2L)")[1:] == [
            "ERR_INVALID_INPUT", "x", "invalid input arguments (argument 'x', via 'n_dims')"]


class TestDataAndASizeFromOneCall:
    """`fx_large_first`: one producer call gives `norms` and the size `n_large` of a work array.

    The expert computes both, whole, from that one call: a size can depend on the values the call
    computes, so it is not paired with data the caller passed. `norms` is declared before `x`,
    and a second producer, called first, is given `n_items`: it is read off `x`, which the caller
    passed, before any producer runs (R4). It was read off `norms`, which a producer fills --
    after the producers, so Python raised `UnboundLocalError`, and R read a global of that name.
    """

    #: 2.5 + 5, the norms above one, plus a thousand for each of the 3 x 2 slots of scratch
    TOTAL = 7.5 + 6000.0

    def test_the_expert_in_python(self, built):
        out = run_python(built, _NORMS_X_PY + """
            print(json.dumps([tox.fx_large_first_expert(2, x, [2, 3, 1]), tox.fx_large_first(2, x)]))
        """)

        assert out == [self.TOTAL, self.TOTAL]

    def test_the_expert_in_r(self, built):
        out = run_r(built, _NORMS_X_R + """
            cat(fx_large_first_expert(2L, x, c(2L, 3L, 1L)), fx_large_first(2L, x))
        """)

        assert out == [str(self.TOTAL)] * 2

    def test_r_reads_no_global_of_the_array_s_or_the_extent_s_name(self, built):
        out = run_r(built, _NORMS_X_R + """
            norms <- c(100, 200, 300, 400, 500); n_items <- 99L; n_large <- 7L
            cat(fx_large_first_expert(2L, x, c(2L, 3L, 1L)))
        """)

        assert out == [str(self.TOTAL)]


class TestAnIntegerArrayTheCallerPasses:
    """`fx_rounded_total_expert` takes the integer array a producer fills in `foo`.

    The R shim checks it as it converts it, like any integer array a caller passes: the routing
    the E2 attack's prototype missed let a double through, and NA reached Fortran as INT_MIN.
    """

    def test_python(self, built):
        out = run_python(built, """
            print(json.dumps([tox.fx_rounded_total([1.2, 2.7, -0.4]),
                              tox.fx_rounded_total_expert([1.2, 2.7, -0.4], [1, 3, 0])]))
        """)

        assert out == [4.0, 4.0]

    def test_r_takes_whole_doubles_and_refuses_na(self, built):
        assert r_raised(built, "fx_rounded_total_expert(c(1.2, 2.7, -0.4), c(1, 3, 0))") == [
            "returned", "4"]
        assert r_raised(built, "fx_rounded_total_expert(c(1.2, 2.7, -0.4), c(1L, NA, 0L))")[2] \
            == "rounded"
        assert r_raised(built, "fx_rounded_total_expert(c(1.2, 2.7, -0.4), c(1, 2.5, 0))")[2] \
            == "rounded"


class TestAValueOnlyABoundNames:
    """`fx_bounded_total`: `limit` sizes nothing, so the Python and R expert takes it as a
    parameter, and is published for it (E2). Only the expert is called: `foo` checks `x`
    against `limit` before its producer has computed it, a defect of its own."""

    def test_the_expert_validates_by_the_caller_s_limit(self, built):
        out = run_python(built, _PY_PRODUCER_ERROR + """
            print(json.dumps([raised(lambda: tox.fx_bounded_total_expert([1.0, 2.0], 9, 3)),
                              raised(lambda: tox.fx_bounded_total_expert([1.0, 5.0], 9, 3))]))
        """)

        assert out[0] == ["returned", 3.0]
        assert out[1][:3] == ["ToxInputError", "ERR_INVALID_INPUT", "x"]


def function_body(text: str, start: str, end: str) -> str:
    """The text of one generated function, from `start` to the next `end` after it."""
    head = text.index(start)
    tail = text.find(end, head + len(start))
    return text[head:] if tail < 0 else text[head:tail]


class TestACallKeptWhole:
    """Rule K: one producer call gives `norms` and `n_large`, which sizes something the caller
    passes or receives, so the call is kept whole -- the wrapper computes neither, and Python and
    R make the call once. Before, `foo` took `norms` over and called the producer again for it.
    """

    def test_python(self, built):
        out = run_python(built, _PY_PRODUCER_ERROR + """
            x = np.array([[2.0, 0.5, 3.0]])
            print(json.dumps([tox.fx_picked_norms(x).tolist(),
                              raised(lambda: tox.fx_weighted_norms(x, np.array([10.0, 20.0]))),
                              raised(lambda: tox.fx_weighted_norms(x, np.array([10.0])))]))
        """)

        assert out[0] == [2.0, 3.0]
        assert out[1] == ["returned", 2.0 + 0.5 + 3.0 + 30.0]
        assert out[2][:3] == ["ToxArgumentError", None, "w"], out[2]

    def test_r(self, built):
        x = "matrix(c(2, 0.5, 3), nrow = 1)"
        assert r_raised(built, f"fx_picked_norms({x})") == ["returned", "2 3"]
        assert r_raised(built, f"fx_weighted_norms({x}, c(10, 20))") == ["returned", "35.5"]
        assert r_raised(built, f"fx_weighted_norms({x}, c(10))")[:3] == [
            "tox_shape_error", "NA", "w"]

    def test_each_binding_makes_the_call_once_and_the_wrapper_none(self, built):
        python = (built / "py" / "fx_impl" / "fx_sized.py").read_text()
        r = (built / "r" / "fx_sized.R").read_text()
        fortran = (built / "fortran" / "fx_sized.F90").read_text()
        for name in ("fx_picked_norms", "fx_weighted_norms"):
            assert function_body(python, f"def {name}(", "\ndef ").count("fx_large_norms(") == 1
            assert function_body(r, f"{name} <- function(", "\n}").count("fx_large_norms(") == 1
            assert "call fx_large_norms(" not in function_body(
                fortran, f"subroutine {name}(", f"end subroutine {name}")
        # the negative control: a call nothing keeps is taken over, so `foo` makes it
        assert "call fx_large_norms(" in function_body(
            fortran, "subroutine fx_large_total(", "end subroutine fx_large_total")

    def test_nothing_is_left_to_take_over_so_there_is_one_wrapper(self, built):
        out = run_python(built, """
            import inspect
            print(json.dumps([list(inspect.signature(tox.fx_weighted_norms).parameters),
                              hasattr(tox, "fx_weighted_norms_expert"),
                              hasattr(tox, "fx_picked_norms_expert")]))
        """)

        assert out == [["x", "w"], False, False]


class TestAnOmittedProducerInput:
    """Rule O: `| padding | OMITTED |` leaves the input out of every call, so the producer does
    what it does for a caller who leaves it out -- `fx_padded_size` adds nothing, and
    `fx_strided_size` takes its stride of three. On the Fortran route (a work array, through `foo`)
    and the binding route (a result) alike; omitting in one table and supplying in another are
    two calls."""

    def test_python(self, built):
        out = run_python(built, """
            x = np.array([1.0, 2.0, 3.0])
            print(json.dumps([tox.fx_unpadded_total(x), len(tox.fx_unpadded(x)),
                              tox.fx_unstrided_total(x), len(tox.fx_unstrided(x)),
                              tox.fx_omitted_once(x)]))
        """)

        assert out == [6 + 3, 3, 6 + 9, 9, 6 + 9 + 6]

    def test_r(self, built):
        out = run_r(built, """
            x <- c(1, 2, 3)
            cat(fx_unpadded_total(x), length(fx_unpadded(x)), fx_unstrided_total(x),
                length(fx_unstrided(x)), fx_omitted_once(x))
        """)

        assert out == ["9", "3", "15", "9", "21"]

    def test_no_call_passes_it(self, built):
        python = (built / "py" / "fx_impl" / "fx_sized.py").read_text()
        r = (built / "r" / "fx_sized.R").read_text()
        fortran = (built / "fortran" / "fx_sized.F90").read_text()
        assert "padding" not in function_body(python, "def fx_unpadded(", "\ndef ").split(
            "fx_padded_size(", 1)[1].split(")", 1)[0]
        assert "stride" not in function_body(r, "fx_unstrided <- function(", "\n}").split(
            "fx_strided_size(", 1)[1].split(")", 1)[0]
        assert "stride =" not in function_body(
            fortran, "subroutine fx_unstrided_total(", "end subroutine fx_unstrided_total")
        # the negative control: a consumer that supplies it passes it
        assert "padding = padding" in function_body(
            fortran, "subroutine fx_padded_total(", "end subroutine fx_padded_total")


class TestAModeReachesAProducerAsAMode:
    """Rule T: a mode argument supplies a producer's mode input whose table names its modes, and
    one more (`fx_method_size` accepts `three`); a constant naming one of the input's mode
    parameters is its code in Fortran and its string in Python and R. Before, Python and R passed
    `METHOD_TWO` as the integer 2, which the producer's binding refused."""

    def test_python(self, built):
        out = run_python(built, """
            x = np.array([1.0, 2.0, 3.0])
            print(json.dumps([tox.fx_method_total(x, "one"), tox.fx_method_total(x, "two"),
                              len(tox.fx_method_filled(x, "one")), len(tox.fx_method_filled(x, "two")),
                              tox.fx_second_total(x), len(tox.fx_second(x))]))
        """)

        assert out == [6 + 3, 6 + 6, 3, 6, 6 + 6, 6]

    def test_r(self, built):
        out = run_r(built, """
            x <- c(1, 2, 3)
            cat(fx_method_total(x, "one"), fx_method_total(x, "two"), length(fx_method_filled(x, "one")),
                length(fx_method_filled(x, "two")), fx_second_total(x), length(fx_second(x)))
        """)

        assert out == ["9", "12", "3", "6", "12", "6"]

    def test_the_constant_crosses_as_its_code_and_its_string(self, built):
        python = (built / "py" / "fx_impl" / "fx_sized.py").read_text()
        r = (built / "r" / "fx_sized.R").read_text()
        fortran = (built / "fortran" / "fx_sized.F90").read_text()
        assert "fx_method_size(n_items=n_items, method='two')" in function_body(
            python, "def fx_second(", "\ndef ")
        assert 'method = "two"' in function_body(r, "fx_second <- function(", "\n}")
        assert "method = 2_int32" in function_body(
            fortran, "subroutine fx_second_total(", "end subroutine fx_second_total")


FIXTURE_CARRIER = Path("helper/codegen/tests/fixtures/carrier")


@pytest.fixture(scope="module")
def carrier_built(tmp_path_factory):
    """`fixtures/carrier` built into a library of its own, R shims included."""
    out = tmp_path_factory.mktemp("e2e_carrier")
    _generate(out, REPO_ROOT / FIXTURE_CARRIER)
    _build(out, REPO_ROOT / FIXTURE_CARRIER)
    return out


class TestAValueThatSizesOnlyAnInput:
    """The carrier case of E2: in the expert, `n_cap` is read off `y`, which carries it.

    Before, the Python and R expert computed it through the producer and refused a `y` of
    another length -- and of any length where the producer refuses its input. Now it reads it
    off `y` and accepts any length, as Fortran's and C's expert always have. The plain tier
    still computes it, and checks `y` against it.
    """

    def test_python_reads_it_off_the_array(self, carrier_built):
        out = run_python(carrier_built, _PY_PRODUCER_ERROR + """
            print(json.dumps([
                raised(lambda: tox.fx_capacity_total_expert([3.0, 1.0], [2, 1], [1.0, 1.0])),
                # six items: the producer refuses more than five, and is not called
                raised(lambda: tox.fx_capacity_total_expert([6.0, 5, 4, 3, 2, 1], [6, 5, 4, 3, 2, 1], [10.0])),
                raised(lambda: tox.fx_capacity_total([3.0, 1.0], [1.0, 1.0, 1.0])),
                raised(lambda: tox.fx_capacity_total([3.0, 1.0], [1.0, 1.0])),
            ]))
        """)

        assert out[0] == ["returned", 1.0 + 2.0]
        assert out[1] == ["returned", 1.0 + 10.0]
        assert out[2] == ["returned", 1.0 + 3.0]
        assert out[3][:3] == ["ToxArgumentError", None, "y"], out[3]

    def test_the_expert_asks_for_no_capacity(self, carrier_built):
        out = run_python(carrier_built, """
            import inspect
            print(json.dumps(list(inspect.signature(tox.fx_capacity_total_expert).parameters)))
        """)

        assert out == ["x", "x_perm", "y"]

    def test_r_reads_it_off_the_array(self, carrier_built):
        assert r_raised(carrier_built,
                        "fx_capacity_total_expert(c(3, 1), c(2L, 1L), c(1, 1))") == ["returned", "3"]
        # six items: the producer refuses more than five, and is not called
        assert r_raised(carrier_built, "fx_capacity_total_expert(c(6, 5, 4, 3, 2, 1), "
                                       "c(6L, 5L, 4L, 3L, 2L, 1L), c(10))") == ["returned", "11"]

    def test_r_computes_it_in_the_plain_tier_and_checks_the_array(self, carrier_built):
        # the shim takes the capacity R computed, and reads no second one off `y`: it compiled
        assert r_raised(carrier_built, "fx_capacity_total(c(3, 1), c(1, 1, 1))") == ["returned", "4"]
        wrong = r_raised(carrier_built, "fx_capacity_total(c(3, 1), c(1, 1))")
        assert wrong[:3] == ["tox_shape_error", "NA", "y"], wrong

    def test_the_r_expert_asks_for_no_capacity(self, carrier_built):
        assert run_r(carrier_built, 'cat(names(formals(fx_capacity_total_expert)))') == [
            "x", "x_perm", "y"]

    def test_python_computes_it_in_the_plain_tier_and_reads_nothing_off_the_array(self):
        # computed through the producer, then `y` is checked against it: no `n_cap = y.shape[0]`
        # first, which the producer's value overwrote (3a). The expert, which reads it off `y`,
        # is the negative control
        from codegen.generate import generate

        paths = Paths(root=REPO_ROOT, src_dir=REPO_ROOT / FIXTURE_CARRIER)
        result = generate(paths=paths, targets=("python",), library="build/libtensor_omics.so")
        [module] = [f.content for f in result.files if f.path.name == "fx_carrier.py"]
        bodies = {name: module.split(f"def {name}(", 1)[1].split("\ndef ", 1)[0]
                  for name in ("fx_capacity_total", "fx_capacity_total_expert",
                               "fx_inout_capacity")}

        assert "n_cap = y.shape[0]" not in bodies["fx_capacity_total"]
        assert "n_cap = fx_room_for(" in bodies["fx_capacity_total"]
        assert "n_cap = y.shape[0]" not in bodies["fx_inout_capacity"]
        assert "n_cap = y.shape[0]" in bodies["fx_capacity_total_expert"]


class TestAComputedCapacityOfAnInOutArray:
    """A lone wrapper whose computed capacity sizes an `intent(inout)` array (3a).

    Nothing is taken over, so there is one wrapper, whose binding computes `n_cap` through the
    producer and checks `y` against it. Its R shim declared `n_cap` twice -- the `.Call`
    parameter R computed, and a second one read off `y` -- and did not compile.
    """

    def test_python(self, carrier_built):
        out = run_python(carrier_built, _PY_PRODUCER_ERROR + """
            y = np.zeros(3)
            first = raised(lambda: tox.fx_inout_capacity([1.0, 2.0], y))
            print(json.dumps([first, y.tolist(),
                              raised(lambda: tox.fx_inout_capacity([1.0, 2.0], np.zeros(2)))]))
        """)

        assert out[0] == ["returned", None]
        assert out[1] == [3.0, 3.0, 3.0]
        assert out[2][:3] == ["ToxArgumentError", None, "y"], out[2]

    def test_r(self, carrier_built):
        assert r_raised(carrier_built, "fx_inout_capacity(c(1, 2), c(0, 0, 0))") == [
            "returned", "3 3 3"]
        wrong = r_raised(carrier_built, "fx_inout_capacity(c(1, 2), c(0, 0))")
        assert wrong[:3] == ["tox_shape_error", "NA", "y"], wrong
