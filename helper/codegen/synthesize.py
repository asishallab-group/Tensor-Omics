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

import re
from collections.abc import Sequence
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
#: Documentation body of a generated re-export module
_REEXPORT_DOC = (
    "Generated from the kernel tree; do not edit -- regenerate instead. "
    "Use this module to reach the whole family; the split into the modules below is an "
    "implementation detail."
)


def generated_wrapper_paths(
    paths: Paths, conventions: Conventions = CONVENTIONS
) -> list[Path]:
    """The `src/generated/tox` files the generator owns, one per kernel module on disk.

    Derived from the kernel tree so the Ford frontend (which excludes them from the parse)
    and the cleaner (which removes them before rewriting) agree on the set without parsing
    -- without either of them having to parse anything.
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


#: Identifiers named in an extent expression, to see which arguments an extent depends on
_EXTENT_IDENTIFIER_RE = re.compile(r"[A-Za-z_]\w*")


def _sizes_a_kept_argument(
    argument: Argument, arguments: Sequence[Argument], conventions: Conventions
) -> bool:
    """Whether an argument the caller still passes or receives is sized by `argument`."""
    name = argument.name.lower()
    for other in arguments:
        if other is argument or is_taken_over(other, conventions):
            continue
        for extent in other.dimension.extents:
            if name in {i.lower() for i in _EXTENT_IDENTIFIER_RE.findall(extent)}:
                return True
    return False


def taken_over_arguments(
    arguments: Sequence[Argument], conventions: Conventions = CONVENTIONS
) -> list[Argument]:
    """The arguments the allocating wrapper prepares itself, rather than taking from the caller.

    Work arrays always qualify -- they are scratch. A `<base>_perm` qualifies only while the
    wrapper can actually build it, which means the `<base>` it orders is an argument too:
    `values_perm` is seeded and sorted against `values`, but a name that merely ends in `_perm`
    and orders something the kernel never receives is the caller's own data. A recommend-sized
    value qualifies as well, though only while nothing the caller still sees depends on it: an
    argument that also gives an extent of a returned array has to stay in the signature, or
    neither the caller nor the binding could size what comes back.
    """
    names = {argument.name.lower() for argument in arguments}
    taken = []
    for argument in arguments:
        if is_temporary(argument, conventions):
            taken.append(argument)
        elif is_permutation(argument, conventions):
            base = argument.name.lower()[: -len(conventions.perm_suffix)]
            if base in names:
                taken.append(argument)
        elif is_computed(argument) and not _sizes_a_kept_argument(
            argument, arguments, conventions
        ):
            taken.append(argument)
    return taken


def sorted_permutations(
    taken: Sequence[Argument], prologue=None, conventions: Conventions = CONVENTIONS
) -> list[Argument]:
    """The permutations the allocating wrapper seeds and heapsorts itself.

    Every `<base>_perm` it took over, except the ones something else builds:

    - a `tmp_` permutation is the kernel's own scratch, seeded and ordered inside it;
    - a permutation the prologue declares `intent(out)` is the prologue's. That is what makes
      a non-default ordering expressible without giving up the expert tier: name the argument
      `<base>_perm` so an expert caller can still supply their own order, and let the prologue
      build it however this family needs. `intent(inout)` means the opposite -- the prologue
      refines an order, so the default sort still runs first and hands it one.

    Shared by the emitter (which emits the calls) and `validate` (which checks nothing the
    sort reads is written afterwards), so the two cannot disagree about what gets sorted.
    """
    built = {
        dummy.name.lower()
        for dummy in (prologue.arguments if prologue is not None else ())
        if dummy.intent is Intent.OUT
    }
    return [
        argument
        for argument in taken
        if is_permutation(argument, conventions)
        and not is_temporary(argument, conventions)
        and argument.name.lower() not in built
    ]


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
    #: names of the generated modules that only re-export their children (see
    #: `_reexport_modules`); they hold no procedures, so no `WrapperSpec` names them
    reexports: tuple[str, ...] = ()

    @property
    def generated_names(self) -> set[str]:
        """Every module this synthesis put into the project."""
        return {spec.module_name for spec in self.specs} | set(self.reexports)


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
    kernel_modules = [m for m in project if m.name.lower().endswith(suffix)]

    for module in kernel_modules:
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
    reexports = _reexport_modules(kernel_modules, set(by_module), conventions)
    generated_modules += reexports

    augmented = Project(list(project.modules) + generated_modules, name=project.name)
    return SynthesisResult(
        project=augmented,
        specs=tuple(specs),
        reexports=tuple(module.name for module in reexports),
    )


def _reexport_modules(
    kernel_modules: Sequence[Module], wrapper_modules: set[str], conventions: Conventions
) -> list[Module]:
    """The generated counterparts of the kernel tree's re-export modules.

    A family too big for one file is split into several kernel modules and gathered by a
    parent that holds no procedures of its own and only `use`s its children -- the shape f42
    uses for its serde tree. That parent generates the matching parent over the *wrappers*,
    so `use tox_data_integration` still reaches the whole family and the split stays an
    implementation detail of the kernel tree.

    Only children that generate something are re-exported: a kernel module used for its
    constants or its recommend routines has no generated counterpart to `use`. Parents are
    resolved to a fixed point, so a parent may gather other parents.
    """
    suffix = conventions.kernel_suffix
    candidates = {
        module.name[: -len(suffix)]: sorted(
            {
                used[: -len(suffix)]
                for used in module.uses
                if used.lower().endswith(suffix)
            }
        )
        for module in kernel_modules
        if not module.procedures and module.uses
    }
    if not candidates:
        return []

    locations = {
        module.name[: -len(suffix)]: module.location for module in kernel_modules
    }
    generated = set(wrapper_modules)
    while True:
        resolved = {
            name
            for name, children in candidates.items()
            if any(child in generated for child in children)
        }
        if resolved <= generated:
            break
        generated |= resolved

    return [
        Module(
            name,
            doc=Doc.parse([_REEXPORT_DOC]),
            meta=Meta(),
            location=locations[name],
            uses=[child for child in children if child in generated],
        )
        for name, children in sorted(candidates.items())
        if name in generated
    ]


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


def _mode_prose_line(argument: Argument) -> int | None:
    """The line `DM_REQUIRED_IF_MODE` wrote on this argument, if it can be identified.

    The directive's line number comes from a declaration-line + 1 heuristic in the source
    index, so it is confirmed against the line's own text before anything is deleted: a line
    that does not name the mode parameter is not the macro's, and survives. If the numbering
    assumption ever breaks, the failure is prose that stays rather than prose that vanishes.
    """
    directive = argument.directives.required_if_mode
    if directive is None or directive.line_number is None:
        return None
    for block in argument.doc.blocks:
        if block.line_number == directive.line_number and directive.mode_param in block.text:
            return directive.line_number
    return None


def _as_required(argument: Argument) -> Argument:
    """A copy of a mode-scoped optional, as this mode's wrapper takes it.

    Its `DM_REQUIRED_IF_MODE` directive is dropped: this wrapper is the mode, so the
    conditional-requirement relation no longer applies. What is left is a plain argument,
    mandatory -- unless it carries a `DM_DEFAULT`, which the binding supplies when the
    caller omits it. Such an argument is never *required*; there the directive only scopes
    it to its mode (absent from the others), and it stays optional.

    The prose the macro wrote goes with the directive. It conditions on a mode argument this
    wrapper does not have, and calls an argument the split has just made mandatory
    "optional" -- a Python reader sees `rdi_threshold : float` with "This optional argument
    needs to be passed if..." directly beneath it.
    """
    return Argument(
        argument.name,
        argument.type,
        dimension=argument.dimension,
        intent=argument.intent,
        optional=argument.directives.has_default,
        doc=argument.doc.without_line(_mode_prose_line(argument)),
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

    taken = taken_over_arguments(arguments, conventions)
    kept = [a for a in arguments if a not in taken]
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
