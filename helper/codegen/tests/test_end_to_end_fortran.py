"""The implementation layer, end to end: generate a wrapper, compile it, and run it.

The Fortran sibling of `test_end_to_end.py` (Python) and `test_end_to_end_r.py` (R). Every
other generator test stops at the emitted text. This one compiles the emitted `foo` and
`foo_expert` against the real `tox_errors` and `f42_sort` and runs them -- the only thing
that can show that validation rejects what it should, with the argument position the caller
actually sees, and that `foo` allocates, sorts and computes what the implementation computes.

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
#: passes no -std at all): f42_math uses the F2023 `reduce()` locality spec and f42_sort
#: declares `array(n)` before `n` is typed, a GNU extension. Neither is generated code, so
#: holding them to the generator's conformance bar would only break this test.
DEPENDENCY_FLAGS = [flag for flag in WRAPPER_FLAGS if flag != "-std=f2018"]

#: Every module the wrapper `use`s, in dependency order. The real sources, not stubs: a stub
#: could not show that `sort_array_heapsort` resolves for the array and the *allocatable*
#: permutation the wrapper hands it. Three files, 0.35s.
DEPENDENCIES = (
    "src/f42/tox_errors.F90",
    "src/f42/utils/f42_math.F90",   # f42_sort uses clamp
    "src/f42/utils/f42_sort.F90",   # init_perm, sort_array_heapsort
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
        paths=Paths(root=REPO_ROOT, src_dir=impl_dir, generated_dir=out),
        targets=("fortran",),
    )
    # the fixtures are the specification, so they must generate cleanly
    assert result.diagnostics.render() == ""
    for file in result.files:
        file.write()

    for dependency in DEPENDENCIES:
        _compile(REPO_ROOT / dependency, out, DEPENDENCY_FLAGS)
    _compile(impl_dir / "fx_ranks_impl.F90", out, WRAPPER_FLAGS)
    for file in result.files:
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


def run_fortran(built: Path, tmp_path: Path, body: str) -> list[str]:
    """Compile a probe around `body`, link it against `built`, return its output words."""
    source = tmp_path / "probe.F90"
    source.write_text(_PROBE.format(body=textwrap.indent(textwrap.dedent(body), "    ")))
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
