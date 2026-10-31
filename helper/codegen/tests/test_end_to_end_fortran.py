"""The implementation layer, end to end: generate a wrapper, compile it, and run it.

The Fortran sibling of `test_end_to_end.py` (Python) and `test_end_to_end_r.py` (R). Every
other generator test stops at the emitted text. This one compiles the emitted `foo` and
`foo_expert` against the real `tox_errors` and `f42_sort` and runs them -- the only thing
that can show that validation rejects what it should, with the argument position the caller
actually sees, and that `foo` allocates, sorts and computes what the implementation computes.
The C wrappers are compiled too, because a declaration they cannot compile (an extent naming
a constant they do not import) is as much a broken wrapper as a wrong answer.

Skipped without gfortran.
"""

import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

from codegen.config import Paths
from codegen.generate import generate

from conftest import REPO_ROOT

GFORTRAN = shutil.which("gfortran")
FIXTURE_IMPL = Path("helper/codegen/tests/fixtures/impl")

pytestmark = pytest.mark.skipif(GFORTRAN is None, reason="gfortran is not installed")

#: What the generated Fortran is held to: the project's own flags plus -std=f2018 for
#: conformance (as in test_emit_fortran_c.py) and -fcheck=all, so a work array allocated to
#: the wrong extent is a bounds error rather than a plausible answer.
WRAPPER_FLAGS = ["-cpp", "-I.", "-std=f2018", "-ffree-line-length-none", "-g", "-fcheck=all"]

#: The same without -std=f2018. The f42 sources are compiled the way fpm compiles them (it
#: passes no -std at all). They are not generated code, so the generator's conformance bar
#: is not theirs to meet.
DEPENDENCY_FLAGS = [flag for flag in WRAPPER_FLAGS if flag != "-std=f2018"]

#: Every module the wrapper `use`s, in dependency order. The real sources, not stubs: a stub
#: could not show that `sort_array_heapsort` resolves for the array and the *allocatable*
#: permutation the wrapper hands it.
DEPENDENCIES = (
    "src/f42/f42_safeguard.F90",         # first: tox_errors names c_bool and depends on it
    "src/f42/tox_errors.F90",
    "src/f42/tox_conversions.F90",       # c_char_as_view, for a C wrapper's mode string
    "src/f42/utils/f42_math_impl.F90",   # f42_sort uses clamp
    "src/f42/utils/f42_sort_impl.F90",   # init_perm, sort_array_heapsort
)


def _compile(source, out: Path, flags) -> Path:
    obj = out / f"{Path(source).stem}.o"
    result = subprocess.run(
        [GFORTRAN, *flags, f"-J{out}", "-c", str(source), "-o", str(obj)],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, f"{source}:\n{result.stderr}"
    return obj


@pytest.fixture(scope="session")
def built(tmp_path_factory):
    """Generate the wrappers for the implementation fixture and compile everything they need."""
    out = tmp_path_factory.mktemp("e2e_fortran")
    impl_dir = REPO_ROOT / FIXTURE_IMPL

    result = generate(
        paths=Paths(root=REPO_ROOT, src_dir=impl_dir, generated_dir=out,
                    c_binding_dir=out / "bindings" / "c"),
        targets=("fortran", "c"),
    )
    # the fixtures are the specification, so they must generate cleanly
    assert result.diagnostics.render() == ""
    for file in result.files:
        file.write()

    for dependency in DEPENDENCIES:
        _compile(REPO_ROOT / dependency, out, DEPENDENCY_FLAGS)
    for impl in sorted(impl_dir.glob("*_impl.F90")):
        _compile(impl, out, WRAPPER_FLAGS)
    # the Fortran wrappers first: each C wrapper uses the one it wraps
    for file in sorted(result.files, key=lambda file: file.path.name.endswith("_c.F90")):
        _compile(file.path, out, WRAPPER_FLAGS)
    return out


#: The declarations every probe shares, so a test body is the two lines that matter
_PROBE = """\
program probe
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fx_ranks, only: rank_scores, rank_scores_expert
    use fx_ranks_impl, only: rank_scores_impl
    implicit none
    integer(int32), parameter :: N = 4
    real(real64) :: scores(N) = [3.0_real64, 1.0_real64, 4.0_real64, 2.0_real64]
    real(real64) :: ranked(N), shifted(N), by_hand(N)
    ! the permutation that sorts `scores` ascending, for calling the implementation by hand
    integer(int32) :: perm(N) = [2, 4, 1, 3]
    integer(int32) :: ierr

{body}
end program probe
"""


#: The probe for `fx_shapes_impl`: the extent forms and purities a wrapper must reproduce
_SHAPES_PROBE = """\
program probe
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fx_shapes, only: sum_values, axis_sums, note_value
    use fx_shapes_impl, only: FX_N_AXES, n_noted
    implicit none
    real(real64) :: values(3) = [1.0_real64, 2.0_real64, 4.0_real64]
    real(real64) :: points(FX_N_AXES, 3) = reshape([1, 10, 2, 20, 3, 30] * 1.0_real64, [FX_N_AXES, 3])
    real(real64) :: scalar, sums(FX_N_AXES)
    integer(int32) :: ierr

{body}
end program probe
"""


#: The probe for `fx_functions_impl`: functions, published as subroutines with the result
#: right before `ierr`, the call passing every argument by position to pin that order
_FUNCTIONS_PROBE = """\
program probe
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use fx_functions, only: fx_total, fx_doubled, fx_all_positive, fx_halve, fx_share, &
        fx_smallest, fx_smallest_expert
    implicit none
    real(real64) :: values(3) = [1.0_real64, 2.0_real64, 5.0_real64]
    real(real64) :: scalar, doubled(3), copy(3)
    logical(c_bool) :: flag
    integer(int32) :: ierr

{body}
end program probe
"""


def run_fortran(built: Path, tmp_path: Path, body: str, probe: str = _PROBE) -> list[str]:
    """Compile a probe around `body`, link it against `built`, return its output words."""
    source = tmp_path / "probe.F90"
    source.write_text(probe.format(body=textwrap.indent(textwrap.dedent(body), "    ")))
    exe = tmp_path / "probe"
    compiled = subprocess.run(
        [GFORTRAN, *WRAPPER_FLAGS, f"-I{built}", f"-J{tmp_path}",
         "-o", str(exe), str(source), *sorted(str(o) for o in built.glob("*.o"))],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert compiled.returncode == 0, compiled.stderr
    run = subprocess.run([str(exe)], cwd=tmp_path, capture_output=True, text=True)
    assert run.returncode == 0, f"stdout:\n{run.stdout}\nstderr:\n{run.stderr}"
    return run.stdout.split()


class TestItCompiles:
    def test_the_emitted_wrapper_is_valid_fortran_2018(self, built):
        # the fixture compiled it with -std=f2018; this is the assertion that names the fact
        assert (built / "fx_ranks.o").exists()


class TestTheAllocatingWrapper:
    def test_it_allocates_sorts_and_computes(self, built, tmp_path):
        # scores arrive unsorted; without init_perm + heapsort this would be 2.00 0.00 3.00 1.00
        out = run_fortran(built, tmp_path, """
            call rank_scores(N, scores, 1.0_real64, ranked, ierr)
            print '(i0, 4(1x, f6.2))', ierr, ranked
        """)

        assert out == ["0", "0.00", "1.00", "2.00", "3.00"]

    def test_it_matches_the_impl_called_by_hand(self, built, tmp_path):
        # the same computation with the permutation and work array supplied by the caller
        out = run_fortran(built, tmp_path, """
            call rank_scores(N, scores, 1.0_real64, ranked, ierr)
            call rank_scores_impl(N, scores, perm, shifted, 1.0_real64, by_hand)
            print '(l1)', all(ranked == by_hand)
        """)

        assert out == ["T"]


class TestTheValidatingWrapper:
    def test_it_passes_a_valid_call_through_to_the_impl(self, built, tmp_path):
        out = run_fortran(built, tmp_path, """
            call rank_scores_expert(N, scores, perm, shifted, 1.0_real64, ranked, ierr)
            print '(i0, 4(1x, f6.2))', ierr, ranked
        """)

        assert out == ["0", "0.00", "1.00", "2.00", "3.00"]


class TestValidationRejects:
    def test_a_value_below_the_documented_minimum(self, built, tmp_path):
        # ERR_INVALID_INPUT (201) packed with arg_pos 3: min_score is rank_scores's
        # third dummy, and 3*M_ERR_ARG_POS_FACTOR + 201 == 30201
        out = run_fortran(built, tmp_path, """
            call rank_scores(N, scores, -1.0_real64, ranked, ierr)
            print '(i0)', ierr
        """)

        assert out == ["30201"]

    def test_and_the_two_wrappers_number_their_arguments_independently(self, built, tmp_path):
        # the same min_score is the fifth dummy of rank_scores_expert and the third of rank_scores,
        # which drops the work array and the permutation ahead of it
        out = run_fortran(built, tmp_path, """
            call rank_scores_expert(N, scores, perm, shifted, -1.0_real64, ranked, ierr)
            print '(i0)', ierr
            call rank_scores(N, scores, -1.0_real64, ranked, ierr)
            print '(i0)', ierr
        """)

        assert out == ["50201", "30201"]

    def test_a_non_finite_value_though_nothing_asked_for_it(self, built, tmp_path):
        # finiteness is the framework's default contract, not an opt-in: ERR_NAN_INF (204)
        # at arg_pos 2, the scores array
        out = run_fortran(built, tmp_path, """
            scores(3) = ieee_value(0.0_real64, ieee_quiet_nan)
            call rank_scores(N, scores, 1.0_real64, ranked, ierr)
            print '(i0)', ierr
        """)

        assert out == ["20204"]

    def test_an_empty_extent(self, built, tmp_path):
        # ERR_EMPTY_INPUT (202) at arg_pos 1, from the automatic validate_dimension_size
        out = run_fortran(built, tmp_path, """
            call rank_scores(0_int32, scores(1:0), 1.0_real64, ranked(1:0), ierr)
            print '(i0)', ierr
        """)

        assert out == ["10202"]


class TestExtentForms:
    """What the wrappers must reproduce of a declaration they redeclare (issue #221, B13)."""

    def test_an_assumed_shape_input_is_validated_with_its_own_size(self, built, tmp_path):
        # the validator was handed the extent `(:)` verbatim, which is not Fortran
        out = run_fortran(built, tmp_path, """
            call sum_values(values, scalar, ierr)
            print '(i0, 1x, f6.2)', ierr, scalar
            values(2) = ieee_value(0.0_real64, ieee_quiet_nan)
            call sum_values(values, scalar, ierr)
            print '(i0)', ierr
        """, probe=_SHAPES_PROBE)

        # ERR_NAN_INF (204) at arg_pos 1: the whole array was checked, not none of it
        assert out == ["0", "7.00", "10204"]

    def test_an_extent_naming_a_module_constant_compiles_and_runs(self, built, tmp_path):
        # the wrapper and the C wrapper redeclare `points(FX_N_AXES, n_points)`, so both
        # import the constant; neither compiled before
        out = run_fortran(built, tmp_path, """
            call axis_sums(3_int32, points, sums, ierr)
            print '(i0, 2(1x, f6.2))', ierr, sums
        """, probe=_SHAPES_PROBE)

        assert out == ["0", "6.00", "60.00"]
        assert (built / "fx_shapes_c.o").exists()


class TestPurity:
    def test_an_impure_elemental_implementation_gets_an_impure_wrapper(self, built, tmp_path):
        # a `pure` wrapper calling it was a compile error; the call must also really happen
        out = run_fortran(built, tmp_path, """
            call note_value(2.5_real64, scalar, ierr)
            call note_value(3.5_real64, scalar, ierr)
            print '(i0, 1x, i0, 1x, f6.2)', ierr, n_noted, scalar
        """, probe=_SHAPES_PROBE)

        assert out == ["0", "2", "3.50"]


class TestFunctions:
    """A function implementation is published as subroutines: its result an `intent(out)`
    argument right before `ierr` (#221 B1). The wrapper used to `call` it, which does not
    compile, and nothing said so."""

    def test_a_scalar_an_array_and_a_logical_result(self, built, tmp_path):
        out = run_fortran(built, tmp_path, """
            call fx_total(3_int32, values, scalar, ierr)
            print '(i0, 1x, f6.2)', ierr, scalar
            call fx_doubled(3_int32, values, doubled, ierr)
            print '(i0, 3(1x, f6.2))', ierr, doubled
            call fx_all_positive(3_int32, values, flag, ierr)
            print '(i0, 1x, l1)', ierr, flag
        """, probe=_FUNCTIONS_PROBE)

        assert out == ["0", "8.00", "0", "2.00", "4.00", "10.00", "0", "T"]

    def test_an_elemental_function_is_wrapped_like_an_elemental_subroutine(self, built, tmp_path):
        out = run_fortran(built, tmp_path, """
            call fx_halve(5.0_real64, scalar, ierr)
            print '(i0, 1x, f6.2)', ierr, scalar
        """, probe=_FUNCTIONS_PROBE)

        assert out == ["0", "2.50"]

    def test_the_result_goes_before_the_implementation_s_own_ierr(self, built, tmp_path):
        # ERR_DIVISION_BY_ZERO (210), the implementation's own, its position cleared
        out = run_fortran(built, tmp_path, """
            call fx_share(3_int32, values, scalar, ierr)
            print '(i0, 1x, f6.3)', ierr, scalar
            values = 0.0_real64
            call fx_share(3_int32, values, scalar, ierr)
            print '(i0)', ierr
        """, probe=_FUNCTIONS_PROBE)

        assert out == ["0", "0.125", "210"]

    def test_both_tiers_publish_the_result(self, built, tmp_path):
        # `fx_smallest` allocates the work array; the expert tier is handed one
        out = run_fortran(built, tmp_path, """
            call fx_smallest(3_int32, values, scalar, ierr)
            print '(i0, 1x, f6.2)', ierr, scalar
            call fx_smallest_expert(3_int32, values, copy, scalar, ierr)
            print '(i0, 1x, f6.2, 1x, f6.2)', ierr, scalar, copy(3)
        """, probe=_FUNCTIONS_PROBE)

        assert out == ["0", "1.00", "0", "1.00", "5.00"]

    def test_the_wrappers_are_pure_where_the_function_is(self, built, tmp_path):
        text = (built / "fx_functions.F90").read_text()

        assert "pure subroutine fx_total(" in text
        assert "\n    subroutine fx_share(" in text  # it has an intent(out) `ierr`


#: The probe for `fx_omp_impl`: OpenMP roots, compiled here without OpenMP. `in_pure` calls the
#: plain wrapper from a pure procedure, which compiles only while the wrapper is pure
_OMP_PROBE = """\
program probe
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use fx_omp, only: fx_column_sums, fx_column_sums_expert, fx_column_maxima
    implicit none
    real(real64) :: values(2, 3) = reshape([1, 2, 3, 4, 5, 6] * 1.0_real64, [2, 3])
    real(real64) :: sums(3), maxima(3), work(2, 4)
    integer(int32) :: columns, ierr

{body}
contains
    pure subroutine in_pure(values, sums, columns, ierr)
        real(real64), intent(in) :: values(2, 3)
        real(real64), intent(out) :: sums(3)
        integer(int32), intent(out) :: columns, ierr
        call fx_column_sums(2_int32, 3_int32, values, sums, 4_int32, columns, ierr)
    end subroutine in_pure
end program probe
"""


class TestAnOpenMPRootWithoutOpenMP:
    """#218: a root's wrappers, built without OpenMP. The plain one allocates one thread column
    whatever the caller asks for, and is pure; both tiers refuse a count below 1."""

    def test_the_plain_wrapper_allocates_one_column_and_is_pure(self, built, tmp_path):
        out = run_fortran(built, tmp_path, """
            call in_pure(values, sums, columns, ierr)
            print '(i0, 1x, i0, 3(1x, f6.2))', ierr, columns, sums
        """, probe=_OMP_PROBE)

        assert out == ["0", "1", "3.00", "7.00", "11.00"]

    def test_the_expert_tier_takes_the_caller_s_columns(self, built, tmp_path):
        out = run_fortran(built, tmp_path, """
            call fx_column_sums_expert(2_int32, 3_int32, values, sums, 4_int32, columns, work, ierr)
            print '(i0, 1x, i0, 3(1x, f6.2))', ierr, columns, sums
        """, probe=_OMP_PROBE)

        assert out == ["0", "4", "3.00", "7.00", "11.00"]

    def test_zero_threads_is_refused_by_every_tier(self, built, tmp_path):
        # ERR_INVALID_INPUT (201) at n_omp_threads' position, 5, in each wrapper
        out = run_fortran(built, tmp_path, """
            call fx_column_sums(2_int32, 3_int32, values, sums, 0_int32, columns, ierr)
            print '(i0)', ierr
            call fx_column_sums_expert(2_int32, 3_int32, values, sums, 0_int32, columns, work, ierr)
            print '(i0)', ierr
            call fx_column_maxima(2_int32, 3_int32, values, maxima, 0_int32, ierr)
            print '(i0)', ierr
        """, probe=_OMP_PROBE)

        assert out == ["50201", "50201", "50201"]
