"""The pipeline, end to end: sources in, generated interfaces out.

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
from .abi.model import CInterface
from .config import CONVENTIONS, Conventions, Paths
from .diagnostics import DiagnosticBag
from .emit.errors_python import PythonErrorEmitter
from .emit.errors_r import RErrorEmitter
from .emit.fortran_c import FortranCEmitter
from .emit.python_ctypes import PythonEmitter
from .emit.r_wrapper import RWrapperEmitter
from .emit.rcpp import RcppEmitter
from .emit.vscode_snippets import SnippetEmitter
from .frontend.ford_frontend import FordFrontend, ParsedProject
from .ir.errors import ErrorCatalogue
from .ir.roles import analyse_project
from .ir.validate import validate_project

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
    targets: tuple[str, ...] = ("c", "python", "r", "snippets"),
    conventions: Conventions = CONVENTIONS,
    library: str = "build/libtensor-omics.so",
    parsed: ParsedProject | None = None,
) -> Result:
    """Build the interface files, without writing them.

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

    analyse_project(parsed.project, diagnostics, conventions)
    validate_project(parsed.project, diagnostics, conventions)

    interface = build_project(parsed.project, diagnostics, conventions)
    catalogue = _catalogue(parsed, diagnostics, conventions)

    files: list[GeneratedFile] = []
    if "c" in targets:
        files += _c_files(interface, paths)
    if "python" in targets:
        files += _python_files(interface, catalogue, paths, library)
    if "r" in targets:
        files += _r_files(interface, catalogue, paths)
    if "snippets" in targets:
        files += _snippets_files(interface, catalogue, paths)

    return Result(diagnostics=diagnostics, files=files)


def generate_and_write(
    paths: Paths = Paths(),
    targets: tuple[str, ...] = ("c", "python", "r", "snippets"),
    conventions: Conventions = CONVENTIONS,
    library: str = "build/libtensor-omics.so",
    clean: bool = True,
    parsed: ParsedProject | None = None,
) -> Result:
    """Generate and, if there are no errors, write the files to disk.

    `clean` removes each target's output directory first, so a routine that stops being
    exported does not leave a stale wrapper behind. `parsed` is passed through to
    `generate` to reuse an already-parsed project (see there).
    """
    result = generate(paths, targets, conventions, library, parsed=parsed)
    if not result.ok:
        return result

    if clean:
        _clean(targets, paths)
    for file in result.files:
        file.write()
    return result


# -- per target -----------------------------------------------------------------


def _c_files(interface: CInterface, paths: Paths) -> list[GeneratedFile]:
    emitter = FortranCEmitter()
    out = paths.resolve(paths.c_interface_dir)
    return [
        GeneratedFile(out / f"{module.name}.F90", emitter.module(module))
        for module in interface
    ]


def _python_files(interface: CInterface, catalogue, paths: Paths,
                  library: str) -> list[GeneratedFile]:
    emitter = PythonEmitter(library=library)
    out = paths.resolve(paths.python_out_dir)
    files = [
        GeneratedFile(out / "library.py", emitter.library_module()),
        GeneratedFile(out / "__init__.py", emitter.package_init(list(interface))),
        GeneratedFile(out / "error_handling.py", PythonErrorEmitter(catalogue).module()),
    ]
    files += [
        GeneratedFile(out / f"{module.stripped_name}.py", emitter.module(module))
        for module in interface
    ]
    return files


def _r_files(interface: CInterface, catalogue, paths: Paths) -> list[GeneratedFile]:
    rcpp, wrapper = RcppEmitter(), RWrapperEmitter()
    out = paths.resolve(paths.rcpp_out_dir)
    files = [
        GeneratedFile(out / "src" / "tox_marshal.h", rcpp.marshal_header_content()),
        GeneratedFile(out / "R" / "tox_validate.R", wrapper.validators()),
        GeneratedFile(out / "R" / "error_handling.R", RErrorEmitter(catalogue).module()),
    ]
    for module in interface:
        files.append(
            GeneratedFile(out / "src" / f"{module.stripped_name}.cpp", rcpp.module(module))
        )
        files.append(
            GeneratedFile(out / "R" / f"{module.stripped_name}.R", wrapper.module(module))
        )
    return files


def _snippets_files(interface: CInterface, catalogue, paths: Paths) -> list[GeneratedFile]:
    out = paths.resolve(paths.snippets_dir)
    files = SnippetEmitter().snippets_files(interface, catalogue)
    return [GeneratedFile(out / name, content) for name, content in files.items()]


def _catalogue(parsed: ParsedProject, diagnostics: DiagnosticBag,
               conventions: Conventions) -> ErrorCatalogue:
    module = parsed.project.module(ERROR_MODULE)
    if module is None:
        diagnostics.error(
            f"no '{ERROR_MODULE}' module, so the error handling cannot be generated"
        )
        return ErrorCatalogue((), conventions, parsed.arg_pos_factor)
    return ErrorCatalogue.from_module(
        module,
        diagnostics,
        parsed.project.constant_values(),
        conventions,
        arg_pos_factor=parsed.arg_pos_factor,
    )


def _clean(targets: tuple[str, ...], paths: Paths) -> None:
    directories = []
    if "c" in targets:
        directories.append(paths.resolve(paths.c_interface_dir))
    if "python" in targets:
        directories.append(paths.resolve(paths.python_out_dir))
    if "r" in targets:
        # only the generated subdirectories, so a hand-written DESCRIPTION or NAMESPACE
        # alongside them is left alone
        base = paths.resolve(paths.rcpp_out_dir)
        directories += [base / "src", base / "R"]
    for directory in directories:
        if directory.is_dir():
            shutil.rmtree(directory)
