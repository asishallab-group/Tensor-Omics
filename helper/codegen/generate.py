"""The pipeline, end to end: sources in, generated bindings out.

Wires the stages together in the one order they can run: parse, work out roles, validate,
build the C ABI, emit. Everything a caller might want to drive separately -- a different
source tree, a subset of targets, a dry run -- is a parameter, so the CLI is a thin shell
over this and the tests can call it directly.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path

from .abi.c_abi import build_project
from .abi.model import CBinding
from .config import CONVENTIONS, Conventions, Paths
from .diagnostics import DiagnosticBag
from .emit.errors_python import PythonErrorEmitter
from .emit.errors_r import RErrorEmitter
from .emit.fortran_c import FortranCEmitter
from .emit.fortran_wrapper import FortranWrapperEmitter
from .emit.python_ctypes import PythonEmitter
from .emit.r_wrapper import RWrapperEmitter
from .emit.c_call import CCallEmitter
from .emit.vscode_snippets import SnippetEmitter
from .frontend.ford_frontend import FordFrontend, ParsedProject
from .ir.errors import ErrorCatalogue
from .ir.roles import analyse_project
from .ir.validate import validate_project
from .synthesize import (
    SynthesisResult,
    generated_wrapper_paths,
    synthesize_wrappers,
)

#: The module the error catalogue is read from
ERROR_MODULE = "tox_errors"


@dataclass
class GeneratedFile:
    path: Path
    content: str

    def write(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(self.content)


@dataclass
class Result:
    diagnostics: DiagnosticBag
    files: list[GeneratedFile] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.diagnostics.errors


def generate(
    paths: Paths = Paths(),
    targets: tuple[str, ...] = ("fortran", "c", "python", "r", "snippets"),
    conventions: Conventions = CONVENTIONS,
    library: str = "build/libtensor-omics.so",
    parsed: ParsedProject | None = None,
) -> Result:
    """Build the generated files, without writing them.

    Returns them in `Result.files` so a caller can inspect or diff before committing to
    disk. Diagnostics are collected throughout; `Result.ok` is false if any are errors,
    and no files should be written in that case.

    `parsed` lets a caller supply an already-parsed project and skip the Ford run -- the
    one slow, side-effect-free stage. The rest of the pipeline is idempotent on it, so a
    test that exercises many generate options over the same real source can parse it once
    and reuse it. When omitted, the source at `paths` is parsed as usual.
    """
    diagnostics = DiagnosticBag()
    if parsed is None:
        parsed = FordFrontend(paths, diagnostics, conventions).parse()

    # Synthesise the wrappers each kernel implies and inject them before the semantic pass,
    # so the C / Python / R targets wrap them like any procedure read from source and the
    # single kernel parse is the one source of truth. Runs unconditionally -- a project with
    # no kernels comes back unchanged, so this is a no-op until the first kernel exists.
    synthesis = synthesize_wrappers(parsed.project, conventions)
    project = synthesis.project

    analyse_project(project, diagnostics, conventions)
    validate_project(project, diagnostics, conventions)

    binding = build_project(project, diagnostics, conventions)

    files: list[GeneratedFile] = []
    if "fortran" in targets:
        files += _fortran_files(project, synthesis, paths, conventions)
    if "c" in targets:
        files += _c_files(binding, paths)
    if any(target in targets for target in ("python", "r", "snippets")):
        catalogue = _catalogue(project, parsed.arg_pos_factor, diagnostics, conventions)
        if "python" in targets:
            files += _python_files(binding, catalogue, paths, library)
        if "r" in targets:
            files += _r_files(binding, catalogue, paths)
        if "snippets" in targets:
            files += _snippets_files(binding, catalogue, paths)

    return Result(diagnostics=diagnostics, files=files)


def generate_and_write(
    paths: Paths = Paths(),
    targets: tuple[str, ...] = ("fortran", "c", "python", "r", "snippets"),
    conventions: Conventions = CONVENTIONS,
    library: str = "build/libtensor-omics.so",
    clean: bool = True,
    parsed: ParsedProject | None = None,
) -> Result:
    """Generate and, if there are no errors, write the files to disk.

    `clean` removes each target's output first, so a routine that stops being exported does
    not leave a stale wrapper behind. `parsed` is passed through to `generate` to reuse an
    already-parsed project (see there).
    """
    result = generate(paths, targets, conventions, library, parsed=parsed)
    if not result.ok:
        return result

    if clean:
        _clean(targets, paths, conventions)
    for file in result.files:
        file.write()
    return result


# -- per target -----------------------------------------------------------------


def _fortran_files(
    project, synthesis: SynthesisResult, paths: Paths, conventions: Conventions
) -> list[GeneratedFile]:
    emitter = FortranWrapperEmitter(
        conventions, macros_header=str(paths.macros_header), project=project
    )
    out = paths.resolve(paths.tox_out_dir)
    generated = {spec.module_name for spec in synthesis.specs}
    return [
        GeneratedFile(out / f"{module.name}.F90", emitter.module(module))
        for module in project
        if module.name in generated
    ]


def _c_files(binding: CBinding, paths: Paths) -> list[GeneratedFile]:
    emitter = FortranCEmitter()
    out = paths.resolve(paths.c_binding_dir)
    return [
        GeneratedFile(out / f"{module.name}.F90", emitter.module(module))
        for module in binding
    ]


def _python_files(binding: CBinding, catalogue, paths: Paths,
                  library: str) -> list[GeneratedFile]:
    emitter = PythonEmitter(library=library)
    out = paths.resolve(paths.python_out_dir)
    files = [
        GeneratedFile(out / "library.py", emitter.library_module()),
        GeneratedFile(out / "__init__.py", emitter.package_init(list(binding))),
        GeneratedFile(out / "error_handling.py", PythonErrorEmitter(catalogue).module()),
    ]
    files += [
        GeneratedFile(out / f"{module.stripped_name}.py", emitter.module(module))
        for module in binding
    ]
    return files


def _r_files(binding: CBinding, catalogue, paths: Paths) -> list[GeneratedFile]:
    emitter, wrapper = CCallEmitter(), RWrapperEmitter()
    # the C `.Call` shims live under src/ so fpm compiles them into the one
    # libtensor-omics.so (mirroring the generated src/c_binding/*.F90); the R-language
    # wrappers live in the R package tree.
    csrc = paths.resolve(paths.r_binding_dir)
    out = paths.resolve(paths.r_out_dir)
    modules = list(binding)
    files = [
        GeneratedFile(csrc / "tox_marshal.h", emitter.marshal_header_content()),
        GeneratedFile(csrc / "init.c", emitter.registration(modules)),
        GeneratedFile(out / "tox_validate.R", wrapper.validators()),
        GeneratedFile(out / "error_handling.R", RErrorEmitter(catalogue).module()),
    ]
    for module in modules:
        files.append(
            GeneratedFile(csrc / f"{module.stripped_name}.c", emitter.module(module))
        )
        files.append(
            GeneratedFile(out / f"{module.stripped_name}.R", wrapper.module(module))
        )
    return files


def _snippets_files(binding: CBinding, catalogue, paths: Paths) -> list[GeneratedFile]:
    out = paths.resolve(paths.snippets_dir)
    files = SnippetEmitter().snippets_files(binding, catalogue)
    return [GeneratedFile(out / name, content) for name, content in files.items()]


def _catalogue(project, arg_pos_factor: int, diagnostics: DiagnosticBag,
               conventions: Conventions) -> ErrorCatalogue:
    module = project.module(ERROR_MODULE)
    if module is None:
        diagnostics.error(
            f"no '{ERROR_MODULE}' module, so the error handling cannot be generated"
        )
        return ErrorCatalogue((), conventions, arg_pos_factor)
    return ErrorCatalogue.from_module(
        module,
        diagnostics,
        project.constant_values(),
        conventions,
        arg_pos_factor=arg_pos_factor,
    )


def _clean(targets: tuple[str, ...], paths: Paths,
           conventions: Conventions = CONVENTIONS) -> None:
    directories = []
    globs = []
    files = []
    if "c" in targets:
        directories.append(paths.resolve(paths.c_binding_dir))
    if "python" in targets:
        directories.append(paths.resolve(paths.python_out_dir))
    if "r" in targets:
        directories.append(paths.resolve(paths.r_binding_dir))   # the C .Call shims
        # the R wrappers sit directly in r_out_dir; remove only the generated `.R` so a
        # hand-written DESCRIPTION or NAMESPACE alongside them is left alone
        globs.append((paths.resolve(paths.r_out_dir), "*.R"))
    if "fortran" in targets:
        # src/tox holds hand-written modules during the migration, so remove only the
        # files the generator owns (one per kernel), never the whole directory
        files += generated_wrapper_paths(paths, conventions)
    for directory in directories:
        if directory.is_dir():
            shutil.rmtree(directory)
    for directory, pattern in globs:
        for path in directory.glob(pattern):
            path.unlink()
    for path in files:
        if path.exists():
            path.unlink()
