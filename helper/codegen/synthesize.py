"""Synthesise the wrappers a kernel implies.

An author writes only the kernel -- `loess_fit_kernel` in module `tox_loess_kernel`, with
no `ierr` and no input validation. From it the generator builds the validating wrapper
`loess_fit` (validate the inputs, then call the kernel) and, when the kernel needs work
arrays, the allocating wrapper `loess_fit_alloc` (allocate the work arrays, build and sort
the permutations, call the recommend routines, then call the kernel). Both land in a
generated module `tox_loess`, carry the export category, and the rest of the pipeline
wraps them to C/Python/R exactly as it does any procedure read from source.

The wrappers are built here as IR and injected into the project *before* the semantic
pass, so a single parse of the kernel is the one source of truth: the Fortran wrappers and
their bindings both come from it, with no round-trip through Ford. The kernel modules stay
in the project but are inert -- nothing generates from an unexported procedure.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

from .config import CONVENTIONS, Conventions, Paths
from .diagnostics import DiagnosticBag
from .ir.doc import Doc
from .ir.entities import Argument, Meta, Module, Procedure, Project
from .ir.roles import analyse
from .ir.types import BaseType, FortranType, Intent

#: Documentation for the `ierr` a wrapper carries when the kernel declares none
_ERROR_DOC = "Error code; zero on success, non-zero on failure."
#: Documentation body of a generated wrapper module (the emitter adds the summary above it)
_MODULE_DOC = "Generated from the kernel; do not edit -- regenerate instead."


def generated_wrapper_paths(
    paths: Paths, conventions: Conventions = CONVENTIONS
) -> list[Path]:
    """The `src/tox` files the generator owns, one per kernel module found on disk.

    Derived from the kernel tree so the Ford frontend (which excludes them from the parse)
    and the cleaner (which removes them before rewriting) agree on the set without parsing
    -- and so `src/tox` can hold hand-written modules during the migration without either
    touching them.
    """
    kernel_dir = paths.resolve(paths.kernel_src_dir)
    out_dir = paths.resolve(paths.tox_out_dir)
    suffix = conventions.kernel_suffix
    if not kernel_dir.is_dir():
        return []
    generated = []
    for path in sorted(kernel_dir.rglob("*.[Ff]90")):
        stem = path.stem
        if stem.lower().endswith(suffix):
            generated.append(out_dir / f"{stem[: -len(suffix)]}.F90")
    return generated


# -- what the allocating wrapper takes over from the caller ----------------------
#
# These are the arguments a caller should not have to supply: they are dropped from the
# allocating wrapper's signature and become locals it prepares itself. The predicates are
# name/directive based, so they hold both before the semantic pass (synthesis, to shape the
# signature) and after it (emission, to fill the body).


def is_temporary(argument: Argument, conventions: Conventions = CONVENTIONS) -> bool:
    """A `tmp_` work array: allocated silently, never returned."""
    return argument.name.lower().startswith(conventions.temporary_prefix)


def is_permutation(argument: Argument, conventions: Conventions = CONVENTIONS) -> bool:
    """A `<base>_perm` vector: seeded and sorted in the allocating wrapper."""
    return argument.name.lower().endswith(conventions.perm_suffix)


def is_computed(argument: Argument) -> bool:
    """Sized by a recommend routine (`DM_OUTPUT_FROM(..., AUTO)`)."""
    directive = argument.directives.output_from
    return directive is not None and directive.is_automatic


def is_taken_over(argument: Argument, conventions: Conventions = CONVENTIONS) -> bool:
    return (
        is_temporary(argument, conventions)
        or is_permutation(argument, conventions)
        or is_computed(argument)
    )


@dataclass(frozen=True)
class ModeFix:
    """The mode a per-mode wrapper fixes: the argument dropped, and what it is set to."""

    #: The kernel's mode argument, dropped from the wrapper's signature
    argument: str
    #: The parameter the mode is fixed to in the kernel call, e.g. `MODE_DOSAGE_PATTERN`
    parameter: str
    #: The module the parameter lives in, so the wrapper can import it
    module: str


@dataclass(frozen=True)
class WrapperSpec:
    """A synthesised wrapper set and the kernel it came from."""

    kernel: Procedure
    #: The generated module the wrappers live in
    module_name: str
    validating: Procedure
    allocating: Procedure | None = None
    #: Set when this is one of a mode-split kernel's per-mode wrappers
    mode_fix: ModeFix | None = None


@dataclass(frozen=True)
class SynthesisResult:
    project: Project
    specs: tuple[WrapperSpec, ...]


def synthesize_wrappers(
    project: Project,
    conventions: Conventions = CONVENTIONS,
    diagnostics: DiagnosticBag | None = None,
) -> SynthesisResult:
    """Build the wrappers every kernel implies and inject them into `project`.

    A kernel is a `<name>_kernel` procedure in a `<module>_kernel` module. Its wrappers go
    into a generated module named by dropping the module's `_kernel` suffix, so
    `tox_loess_kernel` yields `tox_loess`.

    A kernel is analysed here, before its wrappers are built, so its mode table is known: a
    mode argument whose table names a procedure per mode makes the kernel expand into one
    wrapper (pair) per mode instead of a single procedure taking the mode. The kernel is not
    exported, so the later project-wide pass leaves it alone; the wrappers, cloned from it,
    are analysed there like any procedure.
    """
    diagnostics = diagnostics if diagnostics is not None else DiagnosticBag()
    suffix = conventions.kernel_suffix
    by_module: dict[str, list[Procedure]] = {}
    sources: dict[str, Module] = {}
    specs: list[WrapperSpec] = []

    for module in project:
        if not module.name.lower().endswith(suffix):
            continue
        generated_name = module.name[: -len(suffix)]
        for kernel in module.procedures:
            if not kernel.name.lower().endswith(suffix):
                continue  # a recommend routine or private helper, not a kernel

            analyse(kernel, diagnostics, conventions)
            mode_argument = _split_mode_argument(kernel)

            for base, arguments, mode_fix in _variants(kernel, mode_argument, conventions):
                validating, allocating = _wrappers_for(base, arguments, kernel, conventions)
                wrappers = [validating] + ([allocating] if allocating is not None else [])
                by_module.setdefault(generated_name, []).extend(wrappers)
                sources[generated_name] = module
                specs.append(
                    WrapperSpec(
                        kernel=kernel,
                        module_name=generated_name,
                        validating=validating,
                        allocating=allocating,
                        mode_fix=mode_fix,
                    )
                )

    generated_modules = [
        Module(
            name,
            procedures=procedures,
            # a fresh doc rather than the kernel module's own, which describes the kernel and
            # not this API; the emitter adds a "Wrappers for <kernel>" summary above it
            doc=Doc.parse([_MODULE_DOC]),
            meta=Meta(),
            location=sources[name].location,
        )
        for name, procedures in by_module.items()
    ]

    augmented = Project(list(project.modules) + generated_modules, name=project.name)
    return SynthesisResult(project=augmented, specs=tuple(specs))


def _base_name(kernel: Procedure, conventions: Conventions) -> str:
    return kernel.name[: -len(conventions.kernel_suffix)]


def _split_mode_argument(kernel: Procedure) -> Argument | None:
    """The kernel's mode argument whose table opts into per-mode splitting, or None.

    Requires the kernel to have been analysed, so the mode table is on its roles.
    """
    for argument in kernel.arguments:
        roles = argument.roles
        if roles is not None and roles.mode is not None and roles.mode.is_split:
            return argument
    return None


def _variants(
    kernel: Procedure, mode_argument: Argument | None, conventions: Conventions
) -> list[tuple[str, list[Argument], ModeFix | None]]:
    """The (base name, exposed kernel arguments, mode fix) of each wrapper to generate.

    One entry for an ordinary kernel; one per mode value for a mode-split kernel.
    """
    if mode_argument is None:
        return [(_base_name(kernel, conventions), list(kernel.arguments), None)]
    return [
        (
            value.procedure_name,
            _mode_kept_arguments(kernel, mode_argument, value, conventions),
            ModeFix(mode_argument.name, value.parameter, value.module),
        )
        for value in mode_argument.roles.mode.values
    ]


def _mode_kept_arguments(
    kernel: Procedure, mode_argument: Argument, value, conventions: Conventions
) -> list[Argument]:
    """The kernel arguments a per-mode wrapper exposes.

    The mode argument is dropped (fixed to `value.parameter` in the call). An argument that
    `DM_REQUIRED_IF_MODE` ties to this mode becomes a mandatory dummy; one tied to another
    mode is dropped; everything else is carried through.
    """
    kept: list[Argument] = []
    for argument in kernel.arguments:
        if argument is mode_argument:
            continue
        required = argument.directives.required_if_mode
        if required is not None and required.mode_arg.lower() == mode_argument.name.lower():
            if required.mode_param.lower() == value.parameter.lower():
                kept.append(_as_required(argument))
            # otherwise it belongs to another mode and is absent here
        else:
            kept.append(argument.with_name(argument.name))
    return kept


def _as_required(argument: Argument) -> Argument:
    """A copy of a mode-scoped optional, as this mode's wrapper takes it.

    Its `DM_REQUIRED_IF_MODE` directive is dropped: this wrapper is the mode, so the
    conditional-requirement relation no longer applies. What is left is a plain argument,
    mandatory -- unless it carries a `DM_DEFAULT`, which the binding supplies when the
    caller omits it. Such an argument is never *required*; there the directive only scopes
    it to its mode (absent from the others), and it stays optional.
    """
    return Argument(
        argument.name,
        argument.type,
        dimension=argument.dimension,
        intent=argument.intent,
        optional=argument.directives.has_default,
        doc=argument.doc,
        directives=replace(argument.directives, required_if_mode=None),
        attributes=argument.attributes,
        location=argument.location,
    )


def _wrappers_for(
    base: str, arguments: list[Argument], kernel: Procedure, conventions: Conventions
) -> tuple[Procedure, Procedure | None]:
    """The validating wrapper, and the allocating one when the exposed arguments need it."""
    foo_arguments = [argument.with_name(argument.name) for argument in arguments]
    _append_error(foo_arguments, kernel, conventions)
    validating = _wrapper(
        base + conventions.validating_suffix, foo_arguments, kernel, conventions
    )

    kept = [a for a in arguments if not is_taken_over(a, conventions)]
    allocating = None
    if len(kept) != len(arguments):
        alloc_arguments = [argument.with_name(argument.name) for argument in kept]
        _append_error(alloc_arguments, kernel, conventions)
        allocating = _wrapper(
            base + conventions.alloc_suffix, alloc_arguments, kernel, conventions
        )
    return validating, allocating


def _wrapper(
    name: str, arguments: list[Argument], kernel: Procedure, conventions: Conventions
) -> Procedure:
    meta = Meta(
        summary=kernel.meta.summary,
        author=kernel.meta.author,
        category=conventions.c_binding_category,
    )
    return Procedure(
        name,
        arguments=arguments,
        doc=kernel.doc,
        meta=meta,
        location=kernel.location,
        conventions=conventions,
    )


def _append_error(
    arguments: list[Argument], kernel: Procedure, conventions: Conventions
) -> None:
    if kernel.argument(conventions.error_arg) is None:
        arguments.append(
            Argument(
                conventions.error_arg,
                FortranType(BaseType.INTEGER, "int32"),
                intent=Intent.OUT,
                doc=Doc.parse([_ERROR_DOC]),
            )
        )
