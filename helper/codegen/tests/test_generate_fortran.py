"""The fortran target: synthesis wired into `generate`, and the migration-safe cleaning.

Uses the `parsed=` seam to drive `generate` from a hand-built impl project without Ford,
and a temp source tree to check that only the generator's own files are cleaned.
"""

from codegen.config import Paths
from codegen.frontend.ford_frontend import ParsedProject
from codegen.generate import _clean, generate
from codegen.ir.entities import Meta, Project
from codegen.synthesize import generated_wrapper_paths

from builders import integer, module, procedure
from test_synthesize import alloc_impl_module, impl_module


def parsed_of(*modules):
    # macros is unused by the fortran target (only the error catalogue reads it)
    return ParsedProject(project=Project(modules), macros=None, arg_pos_factor=10000)


class TestFortranTarget:
    def test_the_wrapper_file_is_generated(self):
        result = generate(targets=("fortran",), parsed=parsed_of(impl_module()))

        assert result.ok, result.diagnostics.render()
        paths = [f.path.as_posix() for f in result.files]
        assert any(p.endswith("src/generated/tox/tox_demo.F90") for p in paths)

    def test_the_content_is_the_wrapper_module(self):
        result = generate(targets=("fortran",), parsed=parsed_of(alloc_impl_module()))

        (file,) = [f for f in result.files if f.path.name == "tox_demo.F90"]
        assert "module tox_demo" in file.content
        assert "subroutine crunch_expert(" in file.content
        assert "subroutine crunch(" in file.content

    def test_nothing_is_generated_without_implementations(self):
        plain = module("tox_plain", procedure("thing", integer("n"), meta=Meta()))
        result = generate(targets=("fortran",), parsed=parsed_of(plain))

        assert result.files == []


class TestGeneratedWrapperPaths:
    def _source_tree(self, tmp_path):
        tox = tmp_path / "src/tox"
        tox.mkdir(parents=True)
        (tox / "tox_loess_impl.F90").write_text("module tox_loess_impl\nend module\n")
        # a hand-written module beside it that is not an implementation
        (tox / "shared.F90").write_text("module shared\nend module\n")
        return tmp_path

    def test_only_impl_modules_are_owned(self, tmp_path):
        root = self._source_tree(tmp_path)

        result = generated_wrapper_paths(Paths(root=root))

        assert [p.name for p in result] == ["tox_loess.F90"]
        assert result[0] == root / "src/generated/tox/tox_loess.F90"

    def test_the_layer_is_mirrored_rather_than_named(self, tmp_path):
        # the same file under a different layer generates under that layer, with nothing in
        # the rule knowing what a layer is -- which is what lets f42 write implementations
        f42 = tmp_path / "src/f42/utils"
        f42.mkdir(parents=True)
        (f42 / "f42_stats_impl.F90").write_text("module f42_stats_impl\nend module\n")

        (result,) = generated_wrapper_paths(Paths(root=tmp_path))

        assert result == tmp_path / "src/generated/f42/utils/f42_stats.F90"

    def test_the_generated_tree_is_never_claimed_as_a_source(self, tmp_path):
        # a wrapper the generator wrote must not itself look like an implementation
        out = tmp_path / "src/generated/tox"
        out.mkdir(parents=True)
        (out / "tox_loess_impl.F90").write_text("module tox_loess_impl\nend module\n")

        assert generated_wrapper_paths(Paths(root=tmp_path)) == []

    def test_none_when_there_is_no_source_tree(self, tmp_path):
        assert generated_wrapper_paths(Paths(root=tmp_path)) == []


class TestCleanRemovesWhatItOwns:
    def test_clean_removes_only_the_files_the_generator_owns(self, tmp_path):
        (tmp_path / "src/tox").mkdir(parents=True)
        (tmp_path / "src/tox/tox_loess_impl.F90").write_text(
            "module tox_loess_impl\nend module\n"
        )
        tox = tmp_path / "src/generated/tox"
        tox.mkdir(parents=True)
        (tox / "tox_loess.F90").write_text("generated")   # one implementation, so one owned file
        # anything else in the tree is not the generator's to delete, even though nothing
        # hand-written lives there any more
        (tox / "stray.F90").write_text("not ours")

        _clean(("fortran",), Paths(root=tmp_path))

        assert not (tox / "tox_loess.F90").exists()
        assert (tox / "stray.F90").exists()
