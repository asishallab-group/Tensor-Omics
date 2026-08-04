"""Synthesise the wrappers a kernel implies.

An author writes only the kernel -- `loess_fit_kernel` in module `tox_loess_kernel`, with
no `ierr` and no input validation. From it the generator builds the validating wrapper
`loess_fit` (validate the inputs, then call the kernel) and, when the kernel needs work
arrays, the allocating wrapper `loess_fit_alloc`. Both land in a generated module
`tox_loess`, carry the export category, and the rest of the pipeline wraps them to
C/Python/R exactly as it does any procedure read from source.

The wrappers are built here as IR and injected into the project *before* the semantic
pass, so a single parse of the kernel is the one source of truth: the Fortran wrappers and
their bindings both come from it, with no round-trip through Ford. The kernel modules stay
in the project but are inert -- nothing generates from an unexported procedure.

Only the validating wrapper is synthesised for now; the allocating wrapper follows.
"""

from __future__ import annotations

from dataclasses import dataclass

from .config import CONVENTIONS, Conventions
from .ir.entities import Argument, Meta, Module, Procedure, Project
from .ir.types import BaseType, FortranType, Intent


@dataclass(frozen=True)
class WrapperSpec:
    """A synthesised wrapper set and the kernel it came from.

    Carried alongside the augmented project so the Fortran emitter has the kernel to hand
    without re-deriving it, and so a later stage can reason about what was generated.
    """

    kernel: Procedure
    #: The generated module the wrappers live in
    module_name: str
    validating: Procedure


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
            by_module.setdefault(generated_name, []).append(validating)
            sources[generated_name] = module
            specs.append(
                WrapperSpec(
                    kernel=kernel, module_name=generated_name, validating=validating
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


def _validating_wrapper(kernel: Procedure, conventions: Conventions) -> Procedure:
    """`foo`: the kernel's arguments plus `ierr`, validating then calling the kernel."""
    base = kernel.name[: -len(conventions.kernel_suffix)]
    arguments = [argument.with_name(argument.name) for argument in kernel.arguments]
    if kernel.argument(conventions.error_arg) is None:
        arguments.append(_error_argument(conventions))

    meta = Meta(
        summary=kernel.meta.summary,
        author=kernel.meta.author,
        category=conventions.c_binding_category,
    )
    return Procedure(
        base,
        arguments=arguments,
        doc=kernel.doc,
        meta=meta,
        location=kernel.location,
        conventions=conventions,
    )


def _error_argument(conventions: Conventions) -> Argument:
    return Argument(
        conventions.error_arg,
        FortranType(BaseType.INTEGER, "int32"),
        intent=Intent.OUT,
    )
