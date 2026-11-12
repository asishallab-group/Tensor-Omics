"""The fortran target: synthesis wired into `generate`, and the migration-safe cleaning.

Uses the `parsed=` seam to drive `generate` from a hand-built impl project without Ford,
and a temp source tree to check that only the generator's own files are cleaned.
"""

from pathlib import Path

from codegen.config import Paths
from codegen.frontend.ford_frontend import ParsedProject
from codegen.generate import _clean, generate
from codegen.ir.entities import Meta, Project
from codegen.synthesize import generated_path_for, generated_wrapper_paths

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


class TestTheMirror:
    """Where an implementation's wrappers are written (`generated_path_for`)."""

    def test_only_an_impl_module_has_wrappers(self, tmp_path):
        paths = Paths(root=tmp_path)

        assert generated_path_for(Path("src/tox/tox_loess_impl.F90"), paths) == (
            tmp_path / "src/generated/tox/tox_loess.F90"
        )
        # a hand-written module beside it that is not an implementation
        assert generated_path_for(Path("src/tox/shared.F90"), paths) is None

    def test_the_layer_is_mirrored_rather_than_named(self, tmp_path):
        # the same file under a different layer generates under that layer, with nothing in
        # the rule knowing what a layer is -- which is what lets f42 write implementations
        result = generated_path_for(Path("src/f42/utils/f42_stats_impl.F90"), Paths(root=tmp_path))

        assert result == tmp_path / "src/generated/f42/utils/f42_stats.F90"

    def test_the_generated_tree_is_never_claimed_as_a_source(self, tmp_path):
        # a wrapper the generator wrote must not itself look like an implementation
        source = Path("src/generated/tox/tox_loess_impl.F90")

        assert generated_path_for(source, Paths(root=tmp_path)) is None


def _tree(root: Path, *files: str) -> None:
    for file in files:
        path = root / file
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("! test file\n")


class TestGeneratedWrapperPaths:
    """What the fortran target owns: every wrapper module on disk, whatever made it (#221 B21)."""

    def test_a_wrapper_whose_implementation_is_gone_is_owned(self, tmp_path):
        # the orphan is the case that matters: it still `use`s a module nothing builds any
        # more, and fpm compiles everything under src/generated
        _tree(tmp_path, "src/tox/tox_loess_impl.F90",
              "src/generated/tox/tox_loess.F90",
              "src/generated/tox/tox_deleted.F90",
              "src/generated/f42/utils/f42_stats.F90")

        result = generated_wrapper_paths(Paths(root=tmp_path))

        assert [p.relative_to(tmp_path).as_posix() for p in result] == [
            "src/generated/f42/utils/f42_stats.F90",
            "src/generated/tox/tox_deleted.F90",
            "src/generated/tox/tox_loess.F90",
        ]

    def test_the_binding_trees_are_not_the_fortran_targets(self, tmp_path):
        # they belong to the c and r targets, so `--target fortran` must leave them alone
        _tree(tmp_path, "src/tox/tox_loess_impl.F90",
              "src/generated/tox/tox_loess.F90",
              "src/generated/bindings/c/tox_loess_c.F90",
              "src/generated/bindings/r/tox_loess.c")

        (result,) = generated_wrapper_paths(Paths(root=tmp_path))

        assert result == tmp_path / "src/generated/tox/tox_loess.F90"

    def test_none_when_there_is_no_source_tree(self, tmp_path):
        # a mistyped --src must not take every wrapper with it
        _tree(tmp_path, "src/generated/tox/tox_loess.F90")

        assert generated_wrapper_paths(Paths(root=tmp_path, src_dir=Path("scr"))) == []

    def test_none_when_nothing_was_generated_yet(self, tmp_path):
        _tree(tmp_path, "src/tox/tox_loess_impl.F90")

        assert generated_wrapper_paths(Paths(root=tmp_path)) == []


class TestCleanRemovesWhatItOwns:
    def test_an_orphaned_wrapper_is_removed(self, tmp_path):
        _tree(tmp_path, "src/tox/tox_loess_impl.F90",
              "src/generated/tox/tox_loess.F90",
              "src/generated/tox/tox_deleted.F90")

        _clean(("fortran",), Paths(root=tmp_path))

        assert not (tmp_path / "src/generated/tox/tox_loess.F90").exists()
        assert not (tmp_path / "src/generated/tox/tox_deleted.F90").exists()

    def test_the_sources_and_the_binding_trees_stay(self, tmp_path):
        kept = ("src/tox/tox_loess_impl.F90", "src/tox/shared.F90",
                "src/generated/bindings/c/tox_loess_c.F90", "src/generated/bindings/r/init.c")
        _tree(tmp_path, *kept, "src/generated/tox/tox_loess.F90")

        _clean(("fortran",), Paths(root=tmp_path))

        assert all((tmp_path / file).exists() for file in kept)
