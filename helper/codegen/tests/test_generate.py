"""The pipeline and CLI: parse to written files.

Checks the orchestration -- that every target is produced, that a clean run writes and a
dirty one does not, that `--check` reports without writing. The emitters themselves are
tested elsewhere; here the question is only that they are wired together correctly.
"""

from pathlib import Path

import pytest

from codegen import cli
from codegen.config import Paths
from codegen.generate import generate, generate_and_write

from conftest import REPO_ROOT

FIXTURE_SRC = Path("helper/codegen/tests/fixtures/src")
#: the real tox_errors lives here, so a full generate can build the error module
REAL_SRC = Path("src/tox")


def paths(root: Path, src: Path = FIXTURE_SRC) -> Paths:
    return Paths(root=root, src_dir=src)


class TestGenerate:
    def test_the_fixture_run_lacks_tox_errors_and_says_so(self):
        # the fixtures have no tox_errors module, so the error module cannot be built
        result = generate(paths(REPO_ROOT))

        assert not result.ok
        assert any("tox_errors" in e.message for e in result.diagnostics.errors)

    def test_a_source_tree_with_tox_errors_generates_cleanly(self):
        # src/tox has tox_errors but no C-interface tags: scaffolding only, no errors
        result = generate(paths(REPO_ROOT, REAL_SRC))

        assert result.ok, result.diagnostics.render()

    def test_every_target_is_produced(self):
        result = generate(paths(REPO_ROOT, REAL_SRC))

        suffixes = {file.path.suffix for file in result.files}
        assert {".py", ".h", ".R"} <= suffixes

    def test_the_error_module_is_generated_for_each_language(self):
        result = generate(paths(REPO_ROOT, REAL_SRC))

        names = {file.path.name for file in result.files}
        assert "error_handling.py" in names
        assert "error_handling.R" in names

    def test_a_single_target_produces_only_that_target(self):
        result = generate(paths(REPO_ROOT, REAL_SRC), targets=("python",))

        suffixes = {file.path.suffix for file in result.files}
        assert suffixes == {".py"}

    def test_the_auto_output_from_warning_is_surfaced(self):
        # fx_cluster's expert twin needs DM_OUTPUT_FROM(AUTO); it is skipped, not emitted
        result = generate(paths(REPO_ROOT), targets=("c",))

        assert any("DM_OUTPUT_FROM" in w.message for w in result.diagnostics.warnings)

    def test_files_are_not_written_by_generate(self, isolated_repo):
        result = generate(paths(isolated_repo, REAL_SRC), targets=("python",))

        # generate() builds content only; nothing reaches disk
        assert not (isolated_repo / "python").exists()
        assert result.files

    def test_nothing_is_written_when_there_are_errors(self, isolated_repo):
        # a fixture run errors on the missing tox_errors, so no file should appear
        result = generate_and_write(paths(isolated_repo, FIXTURE_SRC), targets=("python",))

        assert not result.ok
        assert not (isolated_repo / "python").exists()


class TestGenerateAndWrite:
    def test_it_writes_the_files(self, isolated_repo):
        result = generate_and_write(paths(isolated_repo, REAL_SRC), targets=("python",))

        assert result.ok
        assert (isolated_repo / "python" / "tensor_omics" / "error_handling.py").is_file()

    def test_clean_removes_a_stale_file(self, isolated_repo):
        out = isolated_repo / "python" / "tensor_omics"
        out.mkdir(parents=True)
        stale = out / "was_exported_once.py"
        stale.write_text("# a wrapper for a procedure no longer exported")

        generate_and_write(paths(isolated_repo, REAL_SRC), targets=("python",), clean=True)

        assert not stale.exists()

    def test_no_clean_keeps_a_stale_file(self, isolated_repo):
        out = isolated_repo / "python" / "tensor_omics"
        out.mkdir(parents=True)
        stale = out / "leftover.py"
        stale.write_text("# left alone")

        generate_and_write(paths(isolated_repo, REAL_SRC), targets=("python",), clean=False)

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
