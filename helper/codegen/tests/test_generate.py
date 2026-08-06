"""The pipeline and CLI: parse to written files.

Checks the orchestration -- that every target is produced, that a clean run writes and a
dirty one does not, that `--check` reports without writing. The emitters themselves are
tested elsewhere; here the question is only that they are wired together correctly.
"""

from pathlib import Path

import pytest

from codegen import cli
from codegen.config import CONVENTIONS, Paths
from codegen.diagnostics import DiagnosticBag
from codegen.frontend.ford_frontend import FordFrontend
from codegen.generate import generate, generate_and_write

from conftest import REPO_ROOT

FIXTURE_SRC = Path("helper/codegen/tests/fixtures/src")
#: the real tox_errors lives here, so a full generate can build the error module
# The whole source tree: a DM_DEFAULT may reference a parameter from another package
# (max_angle's default is f42_utils's PI), and constant resolution only sees the modules
# that are parsed. The generator is always run on the full tree, so these end-to-end tests
# are too. Parsing it per test is slow, so the `real_project` fixture parses it once and the
# tests reuse it through generate(parsed=...); only the CLI tests, which go through argv,
# parse for themselves.
REAL_SRC = Path("src")


def paths(root: Path, src: Path = FIXTURE_SRC) -> Paths:
    return Paths(root=root, src_dir=src)


@pytest.fixture(scope="module")
def real_project():
    """Parse the whole real source once, for the tests that generate from it.

    Ford is the one slow stage, and the rest of the pipeline is idempotent on the parsed
    project (roles are re-derived, not accumulated), so every generate() test can pass this
    in as `parsed=` instead of re-running Ford. Only the CLI tests, which drive argv, parse
    for themselves.
    """
    bag = DiagnosticBag()
    parsed = FordFrontend(paths(REPO_ROOT, REAL_SRC), bag, CONVENTIONS).parse()
    assert bag.errors == (), bag.render()
    return parsed


class TestGenerate:
    def test_the_fixture_run_lacks_tox_errors_and_says_so(self):
        # the fixtures have no tox_errors module, so the error module cannot be built
        result = generate(paths(REPO_ROOT))

        assert not result.ok
        assert any("tox_errors" in e.message for e in result.diagnostics.errors)

    def test_a_source_tree_with_tox_errors_generates_cleanly(self, real_project):
        result = generate(paths(REPO_ROOT, REAL_SRC), parsed=real_project)

        assert result.ok, result.diagnostics.render()

    def test_every_target_is_produced(self, real_project):
        result = generate(paths(REPO_ROOT, REAL_SRC), parsed=real_project)

        suffixes = {file.path.suffix for file in result.files}
        assert {".py", ".h", ".R"} <= suffixes

    def test_the_error_module_is_generated_for_each_language(self, real_project):
        result = generate(paths(REPO_ROOT, REAL_SRC), parsed=real_project)

        names = {file.path.name for file in result.files}
        assert "error_handling.py" in names
        assert "error_handling.R" in names

    def test_a_single_target_produces_only_that_target(self, real_project):
        result = generate(paths(REPO_ROOT, REAL_SRC), targets=("python",), parsed=real_project)

        suffixes = {file.path.suffix for file in result.files}
        assert suffixes == {".py"}

    def test_auto_output_from_is_generated_not_skipped(self, real_project):
        # every DM_OUTPUT_FROM(AUTO) in the tree resolves, so none is skipped with a warning
        result = generate(paths(REPO_ROOT, REAL_SRC), targets=("python",), parsed=real_project)

        assert not any("DM_OUTPUT_FROM" in w.message for w in result.diagnostics.warnings)

    def test_files_are_not_written_by_generate(self, isolated_repo, real_project):
        result = generate(paths(isolated_repo, REAL_SRC), targets=("python",), parsed=real_project)

        # generate() builds content only; nothing reaches disk
        assert not (isolated_repo / "python").exists()
        assert result.files

    def test_nothing_is_written_when_there_are_errors(self, isolated_repo):
        # a fixture run errors on the missing tox_errors, so no file should appear
        result = generate_and_write(paths(isolated_repo, FIXTURE_SRC), targets=("python",))

        assert not result.ok
        assert not (isolated_repo / "python").exists()


class TestGenerateAndWrite:
    def test_it_writes_the_files(self, isolated_repo, real_project):
        result = generate_and_write(
            paths(isolated_repo, REAL_SRC), targets=("python",), parsed=real_project
        )

        assert result.ok
        assert (isolated_repo / "python" / "tensor_omics" / "error_handling.py").is_file()

    def test_clean_removes_a_stale_file(self, isolated_repo, real_project):
        out = isolated_repo / "python" / "tensor_omics"
        out.mkdir(parents=True)
        stale = out / "was_exported_once.py"
        stale.write_text("# a wrapper for a procedure no longer exported")

        generate_and_write(
            paths(isolated_repo, REAL_SRC), targets=("python",), clean=True, parsed=real_project
        )

        assert not stale.exists()

    def test_no_clean_keeps_a_stale_file(self, isolated_repo, real_project):
        out = isolated_repo / "python" / "tensor_omics"
        out.mkdir(parents=True)
        stale = out / "leftover.py"
        stale.write_text("# left alone")

        generate_and_write(
            paths(isolated_repo, REAL_SRC), targets=("python",), clean=False, parsed=real_project
        )

        assert stale.exists()


class TestCli:
    def test_check_writes_nothing_and_reports_changes(self, isolated_repo, capsys):
        code = cli.main(["--root", str(isolated_repo), "--src", str(REAL_SRC),
                         "--target", "python", "--check", "--no-color"])

        assert code == 1  # nothing generated yet, so files would change
        assert not (isolated_repo / "python").exists()
        assert "would change" in capsys.readouterr().err

    def test_check_is_up_to_date_after_a_write(self, isolated_repo):
        common = ["--root", str(isolated_repo), "--src", str(REAL_SRC), "--target", "python"]
        assert cli.main(common) == 0

        assert cli.main(common + ["--check"]) == 0

    def test_a_run_that_errors_exits_non_zero(self, isolated_repo, capsys):
        # the fixtures lack tox_errors
        code = cli.main(["--root", str(isolated_repo), "--src", str(FIXTURE_SRC),
                         "--no-color"])

        assert code == 1
        assert "nothing was written" in capsys.readouterr().err

    def test_a_clean_run_reports_what_it_wrote(self, isolated_repo, capsys):
        code = cli.main(["--root", str(isolated_repo), "--src", str(REAL_SRC),
                         "--target", "python", "--no-color"])

        assert code == 0
        assert "generated" in capsys.readouterr().err


@pytest.fixture
def isolated_repo(tmp_path):
    """A tmp dir that looks enough like the repo for a generate run.

    The frontend reads fpm.toml, src and macros.h from the root, so those are linked in;
    output goes to the tmp dir, leaving the real tree untouched.
    """
    for name in ("fpm.toml", "src", "helper"):
        (tmp_path / name).symlink_to(REPO_ROOT / name)
    return tmp_path


class TestAParseThatFailed:
    """A source the frontend could not read stops the run, rather than half-building an IR.

    An argument whose type has no mapping is left with no type at all. Carrying that into the
    semantic pass raised an AttributeError over the top of the diagnostic the author was
    about to read -- a traceback instead of a report, for an ordinary authoring mistake.
    """

    def source(self, tmp_path):
        src = tmp_path / "src"
        src.mkdir()
        (src / "fx_untyped.F90").write_text(
            "#include <src/macros.h>\n"
            "!> summary: a module\n"
            "module fx_untyped\n"
            "    use, intrinsic :: iso_fortran_env, only: int32\n"
            "    M_IMPLICIT_NONE\n"
            "contains\n"
            "    !> M_EXPORT_C\n"
            "    !| summary: p\n"
            "    !| author: A\n"
            "    subroutine fx_p(factor, ierr)\n"
            "        double precision, intent(in) :: factor\n"
            "            !! a type the generator has no C mapping for\n"
            "        integer(int32), intent(out) :: ierr\n"
            "            !! Error code\n"
            "    end subroutine fx_p\n"
            "end module fx_untyped\n"
        )
        return src

    def test_it_reports_instead_of_raising(self, tmp_path):
        result = generate(Paths(root=REPO_ROOT, src_dir=self.source(tmp_path)))

        assert not result.ok
        assert any("unsupported type" in d.message for d in result.diagnostics.errors)

    def test_nothing_is_generated_from_it(self, tmp_path):
        result = generate(Paths(root=REPO_ROOT, src_dir=self.source(tmp_path)))

        assert result.files == []
