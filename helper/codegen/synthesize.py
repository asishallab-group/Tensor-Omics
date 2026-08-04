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

from dataclasses import dataclass

from .config import CONVENTIONS, Conventions
from .ir.entities import Argument, Meta, Module, Procedure, Project
from .ir.types import BaseType, FortranType, Intent


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
class WrapperSpec:
    """A synthesised wrapper set and the kernel it came from."""

    kernel: Procedure
    #: The generated module the wrappers live in
    module_name: str
    validating: Procedure
    allocating: Procedure | None = None


@dataclass(frozen=True)
class SynthesisResult:
    project: Project
    specs: tuple[WrapperSpec, ...]


def synthesize_wrappers(
    project: Project, conventions: Conventions = CONVENTIONS
) -> SynthesisResult:
    """Build the wrappers every kernel implies and inject them into `project`.

    A kernel is a `<name>_kernel` procedure in a `<module>_kernel` module. Its wrappers go
    into a generated module named by dropping the module's `_kernel` suffix, so
    `tox_loess_kernel` yields `tox_loess`.
    """
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

            validating = _validating_wrapper(kernel, conventions)
            wrappers = [validating]
            allocating = _allocating_wrapper(kernel, conventions)
            if allocating is not None:
                wrappers.append(allocating)

            by_module.setdefault(generated_name, []).extend(wrappers)
            sources[generated_name] = module
            specs.append(
                WrapperSpec(
                    kernel=kernel,
                    module_name=generated_name,
                    validating=validating,
                    allocating=allocating,
                )
            )

    generated_modules = [
        Module(
            name,
            procedures=procedures,
            doc=sources[name].doc,
            meta=Meta(),
            location=sources[name].location,
        )
        for name, procedures in by_module.items()
    ]

    augmented = Project(list(project.modules) + generated_modules, name=project.name)
    return SynthesisResult(project=augmented, specs=tuple(specs))


def _base_name(kernel: Procedure, conventions: Conventions) -> str:
    return kernel.name[: -len(conventions.kernel_suffix)]


def _validating_wrapper(kernel: Procedure, conventions: Conventions) -> Procedure:
    """`foo`: the kernel's arguments plus `ierr`, validating then calling the kernel."""
    name = _base_name(kernel, conventions) + conventions.validating_suffix
    arguments = [argument.with_name(argument.name) for argument in kernel.arguments]
    _append_error(arguments, kernel, conventions)
    return _wrapper(name, arguments, kernel, conventions)


def _allocating_wrapper(kernel: Procedure, conventions: Conventions) -> Procedure | None:
    """`foo_alloc`: the kernel's caller-facing arguments plus `ierr`.

    The work arrays, permutations, and recommend-sized values are dropped -- the wrapper
    prepares them itself. None when the kernel takes none of those, so there is nothing to
    allocate and the validating wrapper is the only entry point.
    """
    kept = [a for a in kernel.arguments if not is_taken_over(a, conventions)]
    if len(kept) == len(kernel.arguments):
        return None

    name = _base_name(kernel, conventions) + conventions.alloc_suffix
    arguments = [argument.with_name(argument.name) for argument in kept]
    _append_error(arguments, kernel, conventions)
    return _wrapper(name, arguments, kernel, conventions)


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
            )
        )
