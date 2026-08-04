"""The fortran target: synthesis wired into `generate`, and the migration-safe cleaning.

Uses the `parsed=` seam to drive `generate` from a hand-built kernel project without Ford,
and a temp source tree to check that only the generator's own files are cleaned.
"""

from codegen.config import Paths
from codegen.frontend.ford_frontend import ParsedProject
from codegen.generate import _clean, generate
from codegen.ir.entities import Meta, Project
from codegen.synthesize import generated_wrapper_paths

from builders import integer, module, procedure
from test_synthesize import alloc_kernel_module, kernel_module


def parsed_of(*modules):
    # macros is unused by the fortran target (only the error catalogue reads it)
    return ParsedProject(project=Project(modules), macros=None, arg_pos_factor=10000)


class TestFortranTarget:
    def test_the_wrapper_file_is_generated(self):
        result = generate(targets=("fortran",), parsed=parsed_of(kernel_module()))

        assert result.ok, result.diagnostics.render()
        paths = [f.path.as_posix() for f in result.files]
        assert any(p.endswith("src/tox/tox_demo.F90") for p in paths)

    def test_the_content_is_the_wrapper_module(self):
        result = generate(targets=("fortran",), parsed=parsed_of(alloc_kernel_module()))

        (file,) = [f for f in result.files if f.path.name == "tox_demo.F90"]
        assert "module tox_demo" in file.content
        assert "subroutine crunch(" in file.content
        assert "subroutine crunch_alloc(" in file.content

    def test_nothing_is_generated_without_kernels(self):
        plain = module("tox_plain", procedure("thing", integer("n"), meta=Meta()))
        result = generate(targets=("fortran",), parsed=parsed_of(plain))

        assert result.files == []


class TestGeneratedWrapperPaths:
    def _kernel_tree(self, tmp_path):
        kernel = tmp_path / "src/kernel"
        kernel.mkdir(parents=True)
        (kernel / "tox_loess_kernel.F90").write_text("module tox_loess_kernel\nend module\n")
        # a support module in the kernel tree that is not itself a kernel module
        (kernel / "shared.F90").write_text("module shared\nend module\n")
        return tmp_path

    def test_only_kernel_modules_map_to_tox_files(self, tmp_path):
        root = self._kernel_tree(tmp_path)

        result = generated_wrapper_paths(Paths(root=root))

        assert [p.name for p in result] == ["tox_loess.F90"]
        assert result[0] == root / "src/tox/tox_loess.F90"

    def test_none_when_there_is_no_kernel_tree(self, tmp_path):
        assert generated_wrapper_paths(Paths(root=tmp_path)) == []


class TestCleanIsMigrationSafe:
    def test_clean_removes_only_the_generated_wrapper(self, tmp_path):
        (tmp_path / "src/kernel").mkdir(parents=True)
        (tmp_path / "src/kernel/tox_loess_kernel.F90").write_text(
            "module tox_loess_kernel\nend module\n"
        )
        tox = tmp_path / "src/tox"
        tox.mkdir(parents=True)
        (tox / "tox_loess.F90").write_text("generated")     # the generator owns this
        (tox / "tox_errors.F90").write_text("hand-written")  # must survive the migration

        _clean(("fortran",), Paths(root=tmp_path))

        assert not (tox / "tox_loess.F90").exists()
        assert (tox / "tox_errors.F90").exists()
