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


#: The probe for `fx_extents_impl`'s extent that must equal a sum (#203)
_EQUALS_PROBE = """\
program probe
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use fx_extents, only: grouped_rows
    use tox_errors, only: get_err_code, get_err_arg_pos
    implicit none
    real(real64) :: values(4) = [1.0_real64, 2.0_real64, 4.0_real64, 8.0_real64]
    real(real64) :: sums(3)
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


class TestAnExtentThatMustEqual:
    """#203: `n_rows` must equal `sum(group_sizes)`, and keeps its dimension check."""

    def test_each_input_is_blamed_for_its_own_fault(self, built, tmp_path):
        out = run_fortran(built, tmp_path, """
            call grouped_rows(3_int32, 2_int32, [2_int32, 1_int32], values, sums, ierr)
            print '(i0, 2(1x, f6.2))', ierr, sums(1:2)
            call grouped_rows(4_int32, 2_int32, [2_int32, 1_int32], values, sums, ierr)
            print '(i0, 1x, i0)', get_err_code(ierr), get_err_arg_pos(ierr)
            call grouped_rows(0_int32, 2_int32, [2_int32, 1_int32], values, sums, ierr)
            print '(i0, 1x, i0)', get_err_code(ierr), get_err_arg_pos(ierr)
            call grouped_rows(2_int32, 0_int32, [2_int32, 1_int32], values, sums, ierr)
            print '(i0, 1x, i0)', get_err_code(ierr), get_err_arg_pos(ierr)
            call grouped_rows(2_int32, 2_int32, [5_int32, -3_int32], values, sums, ierr)
            print '(i0, 1x, i0)', get_err_code(ierr), get_err_arg_pos(ierr)
        """, probe=_EQUALS_PROBE)

        assert out == [
            "0", "3.00", "4.00",
            "206", "1",     # a plain mismatch, at n_rows
            "202", "1",     # n_rows = 0 is empty, not a mismatch
            "202", "2",     # n_groups = 0: the comparison runs after it, so not a mismatch either
            "201", "3",     # [5, -3] sums to n_rows, and the count is refused
        ]


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


#: The probe for `fx_positions_impl`: each implementation (and the prologue) reports the
#: position the probe sets, and `report` prints what the wrapper returned for every one in turn
_POSITIONS_PROBE = """\
program probe
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: get_err_code, get_err_arg_pos
    use fx_positions_impl, only: fx_report_pos, fx_report_prologue_pos
    use fx_positions, only: fx_pos_sum, fx_pos_sum_expert, fx_pos_group_a, fx_pos_group_a_expert, &
        fx_pos_group_b, fx_pos_group_b_expert, fx_pos_clip, fx_pos_clip_expert, fx_pos_share, &
        fx_pos_share_expert
    implicit none
    integer(int32), parameter :: N = 3
    real(real64) :: values(N) = [3.0_real64, 1.0_real64, 2.0_real64]
    real(real64) :: total, scratch(N), work(2*N), clipped(N)
    integer(int32) :: perm(N) = [2, 3, 1]
    logical(c_bool) :: mask(N) = .true._c_bool
    integer(int32) :: ierr, k

{body}
contains
    !> one word per position: the wrapper's position, or `ok` where nothing was reported
    subroutine report(k)
        integer(int32), intent(in) :: k
        if (ierr == 0) then
            write(*, '(a)', advance='no') ' ok'
        else if (get_err_code(ierr) /= 201) then
            write(*, '(a, i0)', advance='no') ' code', get_err_code(ierr)
        else
            write(*, '(1x, i0)', advance='no') get_err_arg_pos(ierr)
        end if
        if (k < 0) print *
    end subroutine report
end program probe
"""


class TestAnImplementationPositionIsTranslated:
    """The position an implementation or a prologue reports numbers its OWN dummy list, and
    every wrapper returns the one that numbers ITS list: kept where it passed one of its own
    dummies, 0 for everything it prepared itself (codegen_guide.md 5.14).

    Each line is one wrapper, the implementation reporting 0 (nothing), then every position
    from 1 to one past its last dummy. The error code itself, ERR_INVALID_INPUT, always survives.
    """

    def run(self, built, tmp_path, call, last):
        out = run_fortran(built, tmp_path, f"""
            do k = 0, {last}
                fx_report_pos = k
                {call}
                call report(merge(-1_int32, k, k == {last}))
            end do
        """, probe=_POSITIONS_PROBE)
        return out

    def test_the_allocating_tier_numbers_its_own_list(self, built, tmp_path):
        # fx_pos_sum(n, values, values_mask, n_selected_values, scale, total, ierr): the
        # permutation, the recommend-sized work_size and the work array are its locals
        call = "call fx_pos_sum(N, values, mask, N, 2.0_real64, total, ierr)"
        assert self.run(built, tmp_path, call, 11) == \
            ["ok", "1", "2", "0", "0", "0", "3", "4", "5", "6", "0", "0"]

    def test_an_absent_optional_still_maps(self, built, tmp_path):
        call = "call fx_pos_sum(N, values, mask, N, total=total, ierr=ierr)"
        assert self.run(built, tmp_path, call, 11) == \
            ["ok", "1", "2", "0", "0", "0", "3", "4", "5", "6", "0", "0"]

    def test_the_expert_tier_keeps_every_position_but_ierr(self, built, tmp_path):
        call = "call fx_pos_sum_expert(N, values, perm, 2*N, work, mask, N, 2.0_real64, total, ierr)"
        assert self.run(built, tmp_path, call, 11) == \
            ["ok", "1", "2", "3", "4", "5", "6", "7", "8", "9", "0", "0"]

    def test_a_mode_split_numbers_each_mode_and_tier_its_own_way(self, built, tmp_path):
        # the mode is fixed in all four, and mode b does not pass `weight` at all
        calls = {
            "call fx_pos_group_a(N, values, 2.0_real64, total, ierr)":
                ["ok", "1", "2", "0", "0", "3", "4", "0", "0"],
            "call fx_pos_group_a_expert(N, values, scratch, 2.0_real64, total, ierr)":
                ["ok", "1", "2", "0", "3", "4", "5", "0", "0"],
            "call fx_pos_group_b(N, values, total, ierr)":
                ["ok", "1", "2", "0", "0", "0", "3", "0", "0"],
            "call fx_pos_group_b_expert(N, values, scratch, total, ierr)":
                ["ok", "1", "2", "0", "3", "0", "4", "0", "0"],
        }
        for call, expected in calls.items():
            assert self.run(built, tmp_path, call, 8) == expected, call

    def test_behind_a_prologue(self, built, tmp_path):
        # fx_pos_clip(n, values, clipped, fraction, ierr): the level the prologue produced and
        # the work array are its locals
        call = "call fx_pos_clip(N, values, clipped, 0.5_real64, ierr)"
        assert self.run(built, tmp_path, call, 7) == ["ok", "1", "2", "0", "0", "3", "0", "0"]
        call = "call fx_pos_clip_expert(N, values, 1.0_real64, scratch, clipped, ierr)"
        assert self.run(built, tmp_path, call, 7) == ["ok", "1", "2", "3", "4", "5", "0", "0"]

    def test_the_prologue_s_own_position(self, built, tmp_path):
        # fx_pos_clip_guard(n, values, fraction, level, handled, ierr): its `fraction` is the
        # wrapper's fourth argument, `level` a local, `handled` never an argument
        out = run_fortran(built, tmp_path, """
            fx_report_pos = 0
            do k = 0, 7
                fx_report_prologue_pos = k
                call fx_pos_clip(N, values, clipped, 0.5_real64, ierr)
                call report(merge(-1_int32, k, k == 7))
            end do
        """, probe=_POSITIONS_PROBE)

        assert out == ["ok", "1", "2", "4", "0", "0", "0", "0"]

    def test_a_function_s_result_moves_nothing(self, built, tmp_path):
        call = "call fx_pos_share(N, values, total, ierr)"
        assert self.run(built, tmp_path, call, 5) == ["ok", "1", "2", "0", "0", "0"]
        call = "call fx_pos_share_expert(N, values, scratch, total, ierr)"
        assert self.run(built, tmp_path, call, 5) == ["ok", "1", "2", "3", "0", "0"]


#: The probe for the arrays a producer fills (`fx_sized_impl`): a 2 x 3 matrix whose column
#: norms are 5, 1 and 2.5, so one call gives them unsorted and one of them is not above one
_SIZED_PROBE = """\
program probe
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: get_err_code, get_err_arg_pos
    use fx_sized, only: fx_norm_rank, fx_norm_rank_expert, fx_large_total, fx_default_scale, &
        fx_padded_total, fx_wide_stride_total
    implicit none
    real(real64) :: x(2, 3) = reshape([3.0_real64, 4.0_real64, 0.6_real64, 0.8_real64, &
                                       1.5_real64, 2.0_real64], [2, 3])
    real(real64) :: ranked(3), total
    real(real64) :: values(3) = [1.0_real64, 2.0_real64, 4.0_real64]
    integer(int32) :: ierr

{body}
end program probe
"""


class TestAnArrayAProducerFills:
    """`foo` allocates an array right above the producer call that fills it (R1).

    It used to allocate it after every producer call, so the producer wrote into an unallocated
    array: under -fcheck=all, which these probes are compiled with, that is a runtime error.
    """

    def test_foo_fills_sorts_and_passes_it(self, built, tmp_path):
        # norms 5, 1, 2.5 sorted ascending through norms_perm; the scratch is 3 columns x 2
        out = run_fortran(built, tmp_path, """
            call fx_norm_rank(2, 3, x, 2, ranked, total, ierr)
            print '(i0, 3(1x, f6.2), 1x, f6.2)', ierr, ranked, total
        """, probe=_SIZED_PROBE)

        assert out == ["0", "1.00", "2.50", "5.00", "6.00"]

    def test_the_expert_takes_it_from_the_caller(self, built, tmp_path):
        out = run_fortran(built, tmp_path, """
            block
                real(real64) :: scratch(6)
                call fx_norm_rank_expert(2, 3, x, 2, [9.0_real64, 7.0_real64, 8.0_real64], &
                    [2, 3, 1], 6, scratch, ranked, total, ierr)
            end block
            print '(i0, 3(1x, f6.2))', ierr, ranked
        """, probe=_SIZED_PROBE)

        assert out == ["0", "7.00", "8.00", "9.00"]

    def test_a_producer_error_names_the_caller_s_argument(self, built, tmp_path):
        # fx_column_norms takes at most three rows, a bound fx_norm_rank leaves out: its
        # ERR_INVALID_INPUT at its own n_dims (1) is fx_norm_rank's n_dims, also its first
        out = run_fortran(built, tmp_path, """
            block
                real(real64) :: tall(4, 3) = 1.0_real64
                call fx_norm_rank(4, 3, tall, 2, ranked, total, ierr)
            end block
            print '(i0, 1x, i0)', get_err_code(ierr), get_err_arg_pos(ierr)
        """, probe=_SIZED_PROBE)

        assert out == ["201", "1"]

    def test_one_call_gives_the_array_and_a_size_its_values_decide(self, built, tmp_path):
        # 5 + 1 + 2.5, plus a thousand for each of the two norms above one
        out = run_fortran(built, tmp_path, """
            call fx_large_total(2, 3, x, total, ierr)
            print '(i0, 1x, f8.2)', ierr, total
        """, probe=_SIZED_PROBE)

        assert out == ["0", "2008.50"]

    def test_a_zero_extent(self, built, tmp_path):
        # no columns at all: the array, the size and the scratch are all empty
        out = run_fortran(built, tmp_path, """
            call fx_large_total(2, 0, x(:, 1:0), total, ierr)
            print '(i0, 1x, f8.2)', ierr, total
        """, probe=_SIZED_PROBE)

        assert out == ["0", "0.00"]


class TestAnOptionalSupplier:
    """P1: an optional producer input's supplier, resolved or passed through as it stands."""

    def test_its_default_is_resolved_under_the_name_the_table_gives_it(self):
        # `factor | scale`: left out, `scale` is its default 2, so the scratch is 3 x 2
        out = run_fortran(self.built, self.tmp_path, """
            call fx_default_scale(3, values, total=total, ierr=ierr)
            print '(i0, 1x, f6.2)', ierr, total
            call fx_default_scale(3, values, 3, total, ierr)
            print '(i0, 1x, f6.2)', ierr, total
        """, probe=_SIZED_PROBE)

        assert out == ["0", "13.00", "0", "16.00"]

    def test_an_absent_one_stays_absent_through_to_the_producer(self):
        out = run_fortran(self.built, self.tmp_path, """
            call fx_padded_total(3, values, total=total, ierr=ierr)
            print '(i0, 1x, f6.2)', ierr, total
            call fx_padded_total(3, values, 2, total, ierr)
            print '(i0, 1x, f6.2)', ierr, total
        """, probe=_SIZED_PROBE)

        assert out == ["0", "10.00", "0", "12.00"]

    def test_its_default_is_resolved_where_the_producer_s_input_is_optional_too(self):
        # `stride` defaults to 4 here and to 3 in the producer: left out, `foo` passes the 4 it
        # resolved, so the work array has the 12 slots the implementation expects. Passed
        # through absent, the producer sized it for 3 (9 slots), and 3 more were expected than
        # there were: 9 + 3000
        out = run_fortran(self.built, self.tmp_path, """
            call fx_wide_stride_total(3, values, total=total, ierr=ierr)
            print '(i0, 1x, f8.2)', ierr, total
            call fx_wide_stride_total(3, values, 2, total, ierr)
            print '(i0, 1x, f8.2)', ierr, total
        """, probe=_SIZED_PROBE)

        assert out == ["0", "12.00", "0", "6.00"]

    @pytest.fixture(autouse=True)
    def _paths(self, built, tmp_path):
        self.built, self.tmp_path = built, tmp_path


#: The probe for `fx_staged_impl` and `fx_bounded_total`: bounds naming what `foo` sets later
_STAGED_PROBE = """\
program probe
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tox_errors, only: get_err_code, get_err_arg_pos
    use fx_sized, only: fx_bounded_total
    use fx_staged, only: fx_handled, fx_index, fx_two_points, fx_feeds, fx_refined
    implicit none
    real(real64) :: total, nan
    integer(int32) :: ierr, itotal
    nan = ieee_value(1.0_real64, ieee_quiet_nan)

{body}
end program probe
"""

#: prints an error's code and position, or 0 0
_REPORT = "print '(i0, 1x, i0)', get_err_code(ierr), get_err_arg_pos(ierr)"


class TestABoundIsCheckedOnceItsValuesAreSet:
    """Rule V: in `foo`, a check runs once every value it names is set -- at the top (P0), after
    the producers (P1), or after the prologue (P2) -- and what the caller passes is checked at the
    top all the same. Only the late bound moves, in a second call: the rest of the argument's
    checks, finiteness included, stay at the top, so the prologue never reads it unchecked."""

    def run(self, body):
        return run_fortran(self.built, self.tmp_path, body, probe=_STAGED_PROBE)

    @pytest.mark.parametrize("init", [[], ["-finit-integer=-1"], ["-finit-integer=1000"]])
    def test_after_the_producers_whatever_the_local_held_before(self, init):
        # `x`'s DM_MAX(real(limit, real64)) was checked before `fx_scaled_size` set `limit`:
        # with -finit-integer=-1 valid input was refused, with 1000 invalid input accepted
        built = self.built
        if init:
            built = self.tmp_path / "init"
            built.mkdir()
            for obj in self.built.glob("*.o"):
                if obj.name != "fx_sized.o":
                    (built / obj.name).symlink_to(obj)
            for mod in self.built.glob("*.mod"):
                (built / mod.name).symlink_to(mod)
            source = next(self.built.rglob("fx_sized.F90"))
            _compile(source, built, [*WRAPPER_FLAGS, *init])
        out = run_fortran(built, self.tmp_path, f"""
            call fx_bounded_total(2, [1.0_real64, 2.0_real64], 3, total, ierr)
            {_REPORT}
            call fx_bounded_total(2, [1.0_real64, 7.0_real64], 3, total, ierr)
            {_REPORT}
        """, probe=_STAGED_PROBE)

        assert out == ["0", "0", "201", "2"]

    def test_the_latest_point_a_check_names_wins(self):
        # limit 2 (a producer), level 100 (the prologue): y between them is checked after the
        # prologue, m's equality with the limit after the producers
        out = self.run(f"""
            call fx_two_points(n_items=1, factor=2, y=50, w=3, m=2, total=itotal, ierr=ierr)
            {_REPORT}
            call fx_two_points(n_items=1, factor=2, y=1, w=3, m=2, total=itotal, ierr=ierr)
            {_REPORT}
            call fx_two_points(n_items=1, factor=2, y=50, w=3, m=3, total=itotal, ierr=ierr)
            {_REPORT}
        """)

        assert out == ["0", "0", "201", "3", "206", "5"]

    def test_finiteness_is_checked_before_a_prologue_that_answers_the_call(self):
        # the prologue answers a single item itself; y's bound names the level it writes, but
        # y's finiteness does not wait for it: NaN is refused on the answered path too
        out = self.run(f"""
            call fx_handled(n_items=1, x=[1.0_real64], y=nan, total=total, ierr=ierr)
            {_REPORT}
            call fx_handled(n_items=2, x=[1.0_real64, 2.0_real64], y=5.0_real64, total=total, ierr=ierr)
            {_REPORT}
        """)

        assert out == ["204", "3", "201", "3"]

    def test_an_index_is_checked_before_the_prologue_reads_it(self):
        # k's DM_MIN(1_int32) stays at the top; its DM_MAX(level) waits for the prologue, which
        # reads x(k) -- out of bounds for k = 0, were the whole check moved below it
        out = self.run(f"""
            call fx_index(n_items=3, x=[1.0_real64, 2.0_real64, 3.0_real64], k=2, total=total, ierr=ierr)
            {_REPORT}
            print '(f6.2)', total
            call fx_index(n_items=3, x=[1.0_real64, 2.0_real64, 3.0_real64], k=0, total=total, ierr=ierr)
            {_REPORT}
        """)

        assert out == ["0", "0", "2.00", "201", "3"]

    def test_a_bound_on_an_input_the_producer_reads(self):
        # factor's DM_MIN(2_int32) is the caller's to meet, at the top; the producer refuses 0
        # itself, its position mapped to factor's
        out = self.run(f"""
            call fx_feeds(n_items=3, factor=1, total=itotal, ierr=ierr)
            {_REPORT}
            call fx_feeds(n_items=3, factor=0, total=itotal, ierr=ierr)
            {_REPORT}
            call fx_feeds(n_items=3, factor=2, total=itotal, ierr=ierr)
            {_REPORT}
            print '(i0)', itotal
        """)

        assert out == ["201", "2", "201", "2", "0", "0", "8"]

    def test_a_value_the_prologue_refines_is_checked_as_refined(self):
        # item 5 allows a prologue to refine (intent(inout)) a produced value: limit 3 * 4,
        # halved to 6. y's DM_MAX(limit) runs after the prologue, against the 6
        out = self.run(f"""
            call fx_refined(n_items=3, factor=4, y=5, total=itotal, ierr=ierr)
            {_REPORT}
            print '(i0)', itotal
            call fx_refined(n_items=3, factor=4, y=7, total=itotal, ierr=ierr)
            {_REPORT}
        """)

        assert out == ["0", "0", "11", "201", "3"]

    @pytest.fixture(autouse=True)
    def _paths(self, built, tmp_path):
        self.built, self.tmp_path = built, tmp_path
