"""Emitting the generated wrapper modules.

One `tox_<name>.F90` per kernel module, holding the validating wrapper `foo` and, when the
kernel needs work arrays, the allocating wrapper `foo_alloc`, for every kernel. Unlike the
C wrappers these are ordinary library sources -- no `bind(C)`, no `NO_C_BINDING` guard, no
pointer marshalling.

A validating wrapper: set `ierr` ok (the kernel it calls has none and reports nothing),
validate every input against its documented range -- and, for reals, the framework's
default finiteness contract -- bail on error, then call the kernel unchanged.

An allocating wrapper takes over the work the caller should not do: it drops the kernel's
work arrays, permutations and recommend-sized values from its signature, validates what is
left, calls the recommend routines to size the work arrays, allocates them, seeds and sorts
the permutations, then calls the kernel directly (it has just built the permutation, so
there is nothing left for the validating wrapper to re-check).

The allocating wrapper reads the kernel's full picture -- every argument, its extents, and
its `DM_OUTPUT_FROM` plan -- off its sibling validating wrapper, which carries them all.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..config import CONVENTIONS, Conventions
from ..ir.entities import Argument, Module, Procedure
from ..ir.types import BaseType
from ..render import Writer
from ..synthesize import (
    ModeFix,
    is_computed,
    is_permutation,
    is_temporary,
    taken_over_arguments,
)
from .doc_ford import render_doc
from .fortran_c import _EXTENT_IDENTIFIER_RE, _chunks, _product


#: The exclusive-bound helpers a range expression may wrap a bound in, and where they live.
#: Named by their own module rather than by the `f42_utils` parent that re-exports them, so
#: the import says where the thing is defined -- as a resolved module constant already does.
_BOUND_HELPERS = {"above": "f42_math", "below": "f42_math"}

#: What seeding and sorting a permutation needs, and where those live
_PERMUTATION_HELPERS = {"init_perm": "f42_sort", "sort_array_heapsort": "f42_sort"}

#: The intrinsic-module functions a range expression may call, and where they come from.
#: A bound is often "how many of these values are usable", which is a NaN test away.
_BOUND_INTRINSICS = {
    "ieee_is_nan": "ieee_arithmetic",
    "ieee_is_finite": "ieee_arithmetic",
}

#: Modules that must be imported `use, intrinsic ::`, not plain `use`
_INTRINSIC_MODULES = frozenset({"ieee_arithmetic"})


@dataclass(frozen=True)
class WrapperInfo:
    """What the emitter needs about a generated wrapper beyond its own signature.

    The kernel it calls (its name is not derivable for a mode-split wrapper, which is named
    from the mode table) and, for a per-mode wrapper, the mode it fixes.
    """

    kernel_name: str
    mode_fix: ModeFix | None = None


class FortranWrapperEmitter:
    def __init__(
        self,
        conventions: Conventions = CONVENTIONS,
        macros_header: str = "src/macros.h",
        project=None,
    ):
        self.conventions = conventions
        self.macros_header = macros_header
        #: used to resolve module constants named in a `DM_MIN`/`DM_MAX` expression, so
        #: they can be imported; None in a unit test whose bounds use no constants
        self.project = project
        #: per-wrapper info (kernel name, mode fix), keyed by lower-case procedure name;
        #: set for the duration of a `module()` call. Empty in a unit test that passes none,
        #: which drives the ordinary (non-split) path via the fallbacks below.
        self._wrapper_info: dict[str, WrapperInfo] = {}

    # -- module -----------------------------------------------------------------

    def module(self, module: Module, wrapper_info: dict[str, WrapperInfo] | None = None) -> str:
        if not module.procedures and module.uses:
            return self.reexport_module(module)
        self._wrapper_info = wrapper_info or {}
        writer = Writer()
        writer.line(f"#include <{self.macros_header}>")
        writer.blank()

        link = f"[[{self._kernel_module(module)}(module)]]"
        writer.block(
            render_doc(module.doc, kind="module", summary=f"Wrappers for {link}")
        )
        writer.line(f"module {module.name}")
        with writer.indent():
            self._module_uses(writer, module)
            # Not a bare `implicit none`: a mistyped kernel call would otherwise compile as
            # an implicit external and fail only at link. Generated code is exactly where
            # that must surface at the generator's output, not in a build log.
            writer.line("M_IMPLICIT_NONE")
            writer.line("private")
            writer.blank()
            for procedure in module:
                writer.line(f"public :: {procedure.name}")
        writer.blank()
        writer.line("contains")
        writer.blank()

        for procedure in module:
            with writer.indent():
                writer.block(self.subroutine(procedure, module))
            writer.blank()

        writer.line(f"end module {module.name}")
        return writer.render(trailing_newline=True)

    def reexport_module(self, module: Module) -> str:
        """A parent that gathers a split family, `use`ing its children and nothing else.

        No `only` list and no `private`: re-exporting is the whole point, so whatever a child
        makes public this module makes public too. It declares nothing itself, so it needs
        neither the macros header nor `implicit none`.
        """
        writer = Writer()
        summary = f"Gathers the {module.name} family"
        writer.block(render_doc(module.doc, kind="module", summary=summary))
        writer.line(f"module {module.name}")
        with writer.indent():
            for used in module.uses:
                writer.line(f"use {used}")
        writer.line(f"end module {module.name}")
        return writer.render(trailing_newline=True)

    def _kernel_module(self, module: Module) -> str:
        return f"{module.name}{self.conventions.kernel_suffix}"

    def _module_uses(self, writer: Writer, module: Module) -> None:
        kernel_module = self._kernel_module(module)

        # everything imported besides the intrinsic kinds and tox_errors, keyed by module:
        # the recommend routines, the f42_utils permutation/bound helpers, and any module
        # constant a range expression names
        extra: dict[str, set[str]] = {}
        for producer_module, names in self._producers(module).items():
            extra.setdefault(producer_module, set()).update(names)
        for prologue_module, names in self._prologues(module).items():
            extra.setdefault(prologue_module, set()).update(names)
        if self._has_permutations(module):
            for name, source in _PERMUTATION_HELPERS.items():
                extra.setdefault(source, set()).add(name)
        for helper_module, names in self._bound_imports(module).items():
            extra.setdefault(helper_module, set()).update(names)
        for procedure in module:
            fix = self._mode_fix(procedure)
            if fix is not None:  # the parameter a per-mode wrapper fixes the mode to
                extra.setdefault(fix.module, set()).add(fix.parameter)
            # and the parameters a runtime-mode wrapper checks its argument against
            for argument in procedure.arguments:
                roles = argument.roles
                if roles is None or roles.mode is None or roles.mode.is_split:
                    continue
                for value in roles.mode.values:
                    extra.setdefault(value.module, set()).add(value.parameter)

        # the kernels, plus any recommend routine / constant that lives in the kernel module
        kernel_names = {self._kernel_name(p, module) for p in module}
        kernel_names |= extra.pop(kernel_module, set())
        for chunk in _chunks(sorted(kernel_names), 4):
            writer.line(f"use {kernel_module}, only: {', '.join(chunk)}")

        kinds = sorted({a.type.kind for p in module for a in p.arguments if a.type.kind})
        for chunk in _chunks(kinds, 6):
            writer.line(f"use, intrinsic :: iso_fortran_env, only: {', '.join(chunk)}")

        for other_module, names in sorted(extra.items()):
            intrinsic = ", intrinsic :: " if other_module in _INTRINSIC_MODULES else " "
            for chunk in _chunks(sorted(names), 4):
                writer.line(f"use{intrinsic}{other_module}, only: {', '.join(chunk)}")

        for chunk in _chunks(self._error_imports(module), 4):
            writer.line(f"use tox_errors, only: {', '.join(chunk)}")

    def _bound_imports(self, module: Module) -> dict[str, set[str]]:
        """Modules and names a range expression refers to and so must import.

        A bound may be `above(0.0_real64)` (helpers from f42_utils), call an IEEE inquiry
        (`count(.not. ieee_is_nan(x))`, for a bound that is "how many of these are usable"),
        or name a module constant like `PI`. Identifiers that are the wrapper's own
        arguments, Fortran intrinsics, or literals need nothing.
        """
        identifiers: set[str] = set()
        for procedure in module:
            for argument in procedure.arguments:
                directives = argument.directives
                for directive in (directives.minimum, directives.maximum, directives.sentinel):
                    if directive is not None:
                        identifiers.update(
                            i.lower() for i in _EXTENT_IDENTIFIER_RE.findall(directive.expression)
                        )
        if not identifiers:
            return {}

        result: dict[str, set[str]] = {}
        for helper, source in _BOUND_HELPERS.items():
            if helper in identifiers:
                result.setdefault(source, set()).add(helper)
        for name, source in _BOUND_INTRINSICS.items():
            if name in identifiers:
                result.setdefault(source, set()).add(name)

        dummies = {a.name.lower() for p in module for a in p.arguments}
        if self.project is not None:
            parameters = {p.name.lower(): p for p in self.project.parameters}
            for identifier in identifiers - dummies:
                parameter = parameters.get(identifier)
                if parameter is not None and parameter.parent is not None:
                    result.setdefault(parameter.parent.name, set()).add(parameter.name)
        return result

    def _error_imports(self, module: Module) -> list[str]:
        names = {"set_ok", "is_err"}
        if self._emits_arg_pos_clear(module):
            names.add("clear_err_arg_pos")
        for procedure in module:
            for argument in procedure.arguments:
                validator = self._validator_for(argument)
                if validator is not None:
                    names.add(validator)
        if any(self._is_allocating(p, module) for p in module):
            # M_ALLOCATE reports through set_err / ERR_ALLOC_FAIL
            names |= {"set_err", "ERR_ALLOC_FAIL"}
        if any(self._is_distance_matrix(a) for p in module for a in p.arguments):
            names.add("validate_distance_matrix")
        # both the mask-count convention and the mode-membership check report this way
        if any(
            (a.roles is not None and a.roles.mask_count_of is not None)
            or self._mode_membership_check(a, 1) is not None
            for p in module
            for a in p.arguments
        ):
            names |= {"set_err_once", "ERR_INVALID_INPUT"}
        ordered = ["set_ok", "is_err"]
        return ordered + sorted(names - set(ordered))

    def _prologues(self, module: Module) -> dict[str, set[str]]:
        """The prologue procedures this module's wrappers call, by the module they live in."""
        found: dict[str, set[str]] = {}
        for procedure in module:
            prologue = self._prologue_for(
                self._sibling(module, procedure)
                if self._is_allocating(procedure, module)
                else procedure,
                module,
                is_allocating=self._is_allocating(procedure, module),
            )
            if prologue is not None and prologue.module is not None:
                found.setdefault(prologue.module.name, set()).add(prologue.name)
        return found

    def _producers(self, module: Module) -> dict[str, set[str]]:
        producers: dict[str, set[str]] = {}
        for alloc in module:
            if not self._is_allocating(alloc, module):
                continue
            foo = self._sibling(module, alloc)
            for argument in taken_over_arguments(
                self._kernel_arguments(foo), self.conventions
            ):
                if not is_computed(argument):
                    continue
                plan = foo.argument(argument.name).roles.computed_from
                producers.setdefault(plan.producer.module.name, set()).add(
                    plan.producer.name
                )
        return producers

    def _has_permutations(self, module: Module) -> bool:
        return any(
            is_permutation(a, self.conventions) and not is_temporary(a, self.conventions)
            for alloc in module
            if self._is_allocating(alloc, module)
            for a in self._kernel_arguments(self._sibling(module, alloc))
        )

    # -- subroutine -------------------------------------------------------------

    def _is_allocating(self, procedure: Procedure, module: Module) -> bool:
        """Whether `procedure` is emitted with the allocating body.

        True only for a real allocating wrapper: an `_alloc` variant that has a validating
        sibling to take work arrays over from. A self-allocating *pipeline* kernel whose
        public name merely ends in `_alloc` (so `is_alloc_variant` is true, e.g. `loess_alloc`)
        has no such sibling -- it prepares its own work arrays internally -- and is emitted
        with the ordinary validating body, its foo being a thin validate-then-call shim.
        """
        return procedure.is_alloc_variant and self._sibling(module, procedure) is not None

    def subroutine(self, procedure: Procedure, module: Module) -> str:
        writer = Writer()
        kernel = self._kernel_name(procedure, module)
        link = f"[[{self._kernel_module(module)}(module):{kernel}]]"
        summary = (
            f"Allocates its work arrays, then calls {link}."
            if self._is_allocating(procedure, module)
            else f"Validates its inputs, then calls {link}."
        )
        writer.block(render_doc(procedure.doc, kind="procedure", summary=summary))
        self._signature(writer, procedure)
        with writer.indent():
            if self._is_allocating(procedure, module):
                self._allocating_body(writer, procedure, module)
            else:
                self._validating_body(writer, procedure, module)
        writer.line(f"end subroutine {procedure.name}")
        return writer.render()

    def _validating_body(
        self, writer: Writer, foo: Procedure, module: Module
    ) -> None:
        prologue = self._prologue_for(foo, module, is_allocating=False)
        self._declarations(writer, foo.arguments)
        self._prologue_local(writer, prologue)
        writer.blank()
        self._validation(writer, foo)
        writer.blank()
        self._prologue_call(writer, prologue, foo, module)
        self._kernel_call(writer, foo, module)

    def _allocating_body(
        self, writer: Writer, alloc: Procedure, module: Module
    ) -> None:
        foo = self._sibling(module, alloc)
        taken = taken_over_arguments(self._kernel_arguments(foo), self.conventions)

        prologue = self._prologue_for(foo, module, is_allocating=True)
        materialized = self._materialized_producer_inputs(foo, alloc, taken)
        self._declarations(writer, alloc.arguments)
        self._locals(writer, taken)
        self._materialized_locals(writer, materialized)
        self._prologue_local(writer, prologue)
        writer.blank()
        self._validation(writer, alloc)
        writer.blank()
        self._materialized_defaults(writer, materialized)
        self._prologue_call(writer, prologue, alloc, module)
        self._recommend_calls(writer, foo, alloc, taken, materialized)
        self._allocations(writer, taken)
        self._permutations(writer, taken)
        writer.blank()
        self._kernel_call(writer, foo, module)

    # -- prologue ---------------------------------------------------------------

    def _prologue_for(self, foo: Procedure, module: Module, is_allocating: bool):
        """The prologue this wrapper runs, resolved to the procedure it names.

        Declared on the kernel, since it is the kernel's own contract; which wrappers run it
        is the directive's scope. Returns None when there is none, or when the project is
        unavailable (a unit test), so the ordinary path is unaffected.
        """
        if self.project is None:
            return None
        kernel = self.project.procedure(
            self._kernel_module(module), self._kernel_name(foo, module)
        )
        if kernel is None:
            return None
        directive = kernel.directives.prologue
        if directive is None or not directive.runs_in(is_allocating):
            return None
        return self.project.procedure(directive.module, directive.procedure)

    def _prologue_local(self, writer: Writer, prologue) -> None:
        if prologue is not None:
            writer.line(f"logical :: {self.conventions.prologue_handled_arg}")

    def _prologue_call(self, writer: Writer, prologue, wrapper: Procedure, module: Module) -> None:
        """Run the prologue, and return early if it handled the call itself.

        Its dummies are supplied by name from the wrapper's own arguments, as a recommend
        routine's are; `handled` and `ierr` come from the wrapper. An argument the wrapper
        does not have is left out, so a prologue may take fewer arguments than the kernel.
        """
        if prologue is None:
            return
        error = self.conventions.error_arg
        handled = self.conventions.prologue_handled_arg
        actuals = []
        for dummy in prologue.arguments:
            lowered = dummy.name.lower()
            if lowered == error.lower():
                actuals.append((dummy.name, error))
            elif lowered == handled.lower():
                actuals.append((dummy.name, handled))
            elif wrapper.argument(dummy.name) is not None:
                actuals.append((dummy.name, dummy.name))
        self._call(writer, prologue.name, actuals)
        self._clear_arg_pos(writer)
        writer.line(f"if (is_err({error})) return")
        writer.line(f"if ({handled}) return")
        writer.blank()

    # -- names ------------------------------------------------------------------

    def _base(self, procedure: Procedure) -> str:
        name = procedure.name
        for suffix in (self.conventions.alloc_suffix, self.conventions.validating_suffix):
            if suffix and name.lower().endswith(suffix):
                return name[: -len(suffix)]
        return name

    def _kernel_name(self, procedure: Procedure, module: Module) -> str:
        # a mode-split wrapper is named from the mode table, so its kernel is not derivable
        # from its name; take it from the spec when available
        info = self._wrapper_info.get(procedure.name.lower())
        if info is not None:
            return info.kernel_name
        return f"{self._base(procedure)}{self.conventions.kernel_suffix}"

    def _mode_fix(self, procedure: Procedure) -> ModeFix | None:
        info = self._wrapper_info.get(procedure.name.lower())
        return info.mode_fix if info is not None else None

    def _sibling(self, module: Module, alloc: Procedure) -> Procedure:
        return module.procedure(self._base(alloc) + self.conventions.validating_suffix)

    def _kernel_arguments(self, foo: Procedure) -> list[Argument]:
        """The kernel's own arguments minus ierr -- the ones that may be taken over.

        Used to work out the allocations and permutations; ierr is never among them, so
        dropping it here is harmless and independent of whether the kernel declares one.
        """
        error = self.conventions.error_arg.lower()
        return [a for a in foo.arguments if a.name.lower() != error]

    def _kernel_declares_error(self, foo: Procedure, module: Module) -> bool:
        """Whether the kernel takes `ierr`, so the wrapper has a position to clear."""
        error = self.conventions.error_arg.lower()
        return any(
            a.name.lower() == error for a in self._kernel_call_arguments(foo, module)
        )

    def _emits_arg_pos_clear(self, module: Module) -> bool:
        """Whether any body here calls out to hand-written code that returns `ierr`.

        The one source of truth for both the emitted call and the `use` list -- deciding the
        two separately is how a missing import gets shipped.
        """
        return any(
            self._kernel_declares_error(procedure, module)
            for procedure in module
            if not self._is_allocating(procedure, module)
        ) or bool(self._prologues(module))

    def _kernel_call_arguments(self, foo: Procedure, module: Module) -> list[Argument]:
        """The kernel's own signature, for the call.

        Includes `ierr` exactly when the kernel declares one: a kernel that propagates a
        sub-helper's failure takes `ierr`, a pure one does not, and the wrapper must pass on
        whichever the kernel actually has. Resolved from the kernel itself when the project
        is available; without it (a unit test) the kernel has no ierr, since synthesis only
        appends one to the wrapper.
        """
        if self.project is not None:
            kernel = self.project.procedure(
                self._kernel_module(module), self._kernel_name(foo, module)
            )
            if kernel is not None:
                return list(kernel.arguments)
        return self._kernel_arguments(foo)

    # -- declarations -----------------------------------------------------------

    def _signature(self, writer: Writer, procedure: Procedure) -> None:
        writer.line(f"subroutine {procedure.name}(&")
        names = [a.name for a in procedure.arguments]
        with writer.indent(2):
            for name in names[:-1]:
                writer.line(f"{name},&")
            writer.line(f"{names[-1]}&")
        with writer.indent():
            writer.line(")")

    def _declarations(self, writer: Writer, arguments) -> None:
        for argument in self._declaration_order(arguments):
            writer.line(self._declaration(argument))
            if argument.doc:
                with writer.indent():
                    writer.block(render_doc(argument.doc, kind="argument"))

    @staticmethod
    def _declaration_order(arguments) -> list[Argument]:
        """An extent is declared before the array whose extent it is (gfortran -std=f2018
        rejects the reverse). Everything else keeps the author's order."""
        referenced = {
            identifier.lower()
            for argument in arguments
            for extent in argument.dimension.extents
            for identifier in _EXTENT_IDENTIFIER_RE.findall(extent)
        }
        extents = [a for a in arguments if a.name.lower() in referenced]
        rest = [a for a in arguments if a.name.lower() not in referenced]
        return [*extents, *rest]

    def _declaration(self, argument: Argument) -> str:
        attributes = []
        if argument.dimension:
            attributes.append(f"dimension({', '.join(argument.dimension.extents)})")
        attributes.append(f"intent({argument.intent.value})")
        if argument.optional:
            attributes.append("optional")
        return f"{argument.type}, {', '.join(attributes)} :: {argument.name}"

    def _locals(self, writer: Writer, taken) -> None:
        """The taken-over arguments, as locals the allocating wrapper prepares itself."""
        for argument in taken:
            if argument.is_array:
                deferred = ", ".join(":" for _ in argument.dimension.extents)
                writer.line(
                    f"{argument.type}, dimension({deferred}), allocatable :: {argument.name}"
                )
            else:
                writer.line(f"{argument.type} :: {argument.name}")

    # -- validation -------------------------------------------------------------

    def _validation(self, writer: Writer, procedure: Procedure) -> None:
        error = self.conventions.error_arg
        writer.line(f"call set_ok({error})")
        # scalars (extents, bounds) before arrays, as the hand-written wrappers do; arg_pos
        # stays the argument's true position, so the ordering is presentational only
        calls = [
            (argument.is_array, self._validation_call(argument, position))
            for position, argument in enumerate(procedure.arguments, start=1)
            if argument.name.lower() != error.lower()
        ]
        for is_array in (False, True):
            for array, call in calls:
                if array is is_array and call:
                    writer.line(call)
        for line in self._convention_checks(procedure):
            writer.line(line)
        writer.line(f"if (is_err({error})) return")

    def _convention_checks(self, procedure: Procedure) -> list[str]:
        """Checks driven by a naming convention rather than a per-argument range.

        A distance-matrix argument is validated for structure; a `n_selected_<x>` count is
        validated against `count(<x>_mask)`. Both come from information the wrapper already has
        (the argument's roles / shape), so the kernel need not carry them or an `ierr`.
        """
        error = self.conventions.error_arg
        lines: list[str] = []
        for position, argument in enumerate(procedure.arguments, start=1):
            if self._is_distance_matrix(argument):
                extent = argument.dimension.extents[0]
                lines.append(
                    f"call validate_distance_matrix({argument.name}, {extent}, {error}, "
                    f"arg_pos={position}_int32)"
                )
        for position, argument in enumerate(procedure.arguments, start=1):
            roles = argument.roles
            if roles is None or roles.mask_count_of is None:
                continue
            mask = roles.mask_count_of
            lines.append(
                f"if (count({mask.name}, kind=int32) /= {argument.name}) "
                f"call set_err_once({error}, ERR_INVALID_INPUT, arg_pos={position}_int32)"
            )
        for position, argument in enumerate(procedure.arguments, start=1):
            check = self._mode_membership_check(argument, position)
            if check is not None:
                lines.append(check)
        return lines

    def _mode_membership_check(self, argument: Argument, position: int) -> str | None:
        """Reject a mode value that is not one of the ones its table names.

        The generator parsed that table, so it knows the accepted parameters exactly -- which
        makes this a better answer than a hand-written `DM_MIN`/`DM_MAX` pair: the accepted
        set need not be contiguous, and a range annotation goes stale the moment a mode is
        added, while this is regenerated from the table.

        A chain of `/=` rather than `all(x /= [...])`: the array constructor would build a
        temporary on every call, for a check that is a handful of integer comparisons.

        An optional mode is guarded by `present()`. Absent means "use the documented default",
        which is by construction one of the accepted values -- and reading an absent optional
        to find that out is not a Fortran program.

        Only for a mode the wrapper still takes. A mode-split wrapper *is* its mode and has no
        such argument, and the C layer rejects an unknown mode string before Fortran is
        entered -- but a direct Fortran caller reaches this, and nothing else stops it.
        """
        roles = argument.roles
        if roles is None or roles.mode is None or roles.mode.is_split:
            return None
        values = roles.mode.values
        if not values:
            return None
        error = self.conventions.error_arg
        comparisons = " .and. ".join(
            f"{argument.name} /= {value.parameter}" for value in values
        )
        check = (
            f"if ({comparisons}) "
            f"call set_err_once({error}, ERR_INVALID_INPUT, arg_pos={position}_int32)"
        )
        if argument.optional:
            return f"if (present({argument.name})) then; {check}; end if"
        return check

    def _is_distance_matrix(self, argument: Argument) -> bool:
        """A square real matrix named by the distance-matrix convention.

        It is validated for distance-matrix structure (symmetry, non-negativity, zero diagonal)
        by `validate_distance_matrix`, which no per-argument range validator expresses, so it also
        opts out of the finiteness contract.
        """
        if argument.type.base is not BaseType.REAL or not argument.intent.is_input:
            return False
        extents = argument.dimension.extents
        if len(extents) != 2 or extents[0] != extents[1]:
            return False
        name = argument.name.lower()
        return any(
            name == suffix or name.endswith(f"_{suffix}")
            for suffix in self.conventions.distance_matrix_suffixes
        )

    def _validator_for(self, argument: Argument) -> str | None:
        """The `tox_errors` validator this argument needs, or None."""
        if not argument.intent.is_input:
            return None  # an output carries no value to check
        if self._is_distance_matrix(argument):
            return None  # validated by validate_distance_matrix, emitted as a convention check
        if (
            argument.roles is not None
            and argument.roles.is_extent
            and not argument.directives.has_range
        ):
            # an extent is validated as a dimension by default, but a documented range wins:
            # it lets an extent that is really a count permit zero (an empty selection), which
            # validate_dimension_size rejects
            return "validate_dimension_size"
        base = argument.type.base
        if base is BaseType.REAL:
            # finiteness is the default contract, so a real input is always checked --
            # unless it opts out of both failure modes and states no bound
            directives = argument.directives
            if (
                not directives.has_range
                and directives.allows_nan
                and directives.allows_infinite
            ):
                return None
            return "validate_all_in_range_real" if argument.is_array else "validate_in_range_real"
        if base is BaseType.INTEGER:
            if not argument.directives.has_range:
                return None  # an integer with no bound has nothing to check
            return "validate_all_in_range_int" if argument.is_array else "validate_in_range_int"
        return None

    def _validation_call(self, argument: Argument, position: int) -> str:
        validator = self._validator_for(argument)
        if validator is None:
            return ""
        error = self.conventions.error_arg
        arg_pos = f"arg_pos={position}_int32"

        if validator == "validate_dimension_size":
            return f"call validate_dimension_size({argument.name}, {error}, {arg_pos})"

        if argument.is_array:
            count = _product(argument.dimension.extents)
            head = f"{argument.name}, {count}, {error}"
        else:
            head = f"{argument.name}, {error}"

        keywords = [arg_pos]
        directives = argument.directives
        if directives.minimum is not None:
            keywords.append(f"min={directives.minimum.expression}")
        if directives.maximum is not None:
            keywords.append(f"max={directives.maximum.expression}")
        if directives.sentinel is not None:
            keywords.append(f"sentinel={directives.sentinel.expression}")
        if argument.type.base is BaseType.REAL:
            if directives.allows_nan:
                keywords.append("allow_nan=.true.")
            if directives.allows_infinite:
                keywords.append("allow_infinite=.true.")

        return f"call {validator}({head}, {', '.join(keywords)})"

    # -- allocation, recommend routines, permutations ---------------------------

    def _materialized_producer_inputs(
        self, foo: Procedure, alloc: Procedure, taken
    ) -> dict[str, Argument]:
        """Wrapper arguments a recommend routine needs by value, keyed by name.

        A recommend routine is called directly, so an input the wrapper takes optionally
        cannot simply be forwarded when the routine's own dummy is mandatory: an absent
        optional passed there is not a Fortran program. Such an input has a `DM_DEFAULT`
        (that is what makes it optional in the first place), so the wrapper resolves it into
        a local first -- the argument when the caller gave one, the documented default when
        not -- and hands the routine that.
        """
        materialized: dict[str, Argument] = {}
        for argument in taken:
            if not is_computed(argument):
                continue
            roles = foo.argument(argument.name).roles
            plan = roles.computed_from if roles is not None else None
            if plan is None:
                continue
            for dummy in plan.producer.arguments:
                if dummy.optional:
                    continue
                supplied = alloc.argument(dummy.name)
                if supplied is None or not supplied.optional:
                    continue
                if supplied.directives.default is None:
                    continue
                materialized[supplied.name] = supplied
        return materialized

    @staticmethod
    def _materialized_name(name: str) -> str:
        return f"{name}_value"

    def _materialized_locals(self, writer: Writer, materialized: dict[str, Argument]) -> None:
        for name, argument in materialized.items():
            writer.line(f"{argument.type} :: {self._materialized_name(name)}")

    def _materialized_defaults(
        self, writer: Writer, materialized: dict[str, Argument]
    ) -> None:
        if not materialized:
            return
        for name, argument in materialized.items():
            default = argument.directives.default.expression
            writer.line(
                f"M_DEFAULT_VAL({name}, {self._materialized_name(name)}, {default})"
            )
        writer.blank()

    def _recommend_calls(
        self, writer: Writer, foo: Procedure, alloc: Procedure, taken,
        materialized: dict[str, Argument] | None = None,
    ) -> None:
        """Call each recommend routine once, into the sizes it produces.

        The routine is called directly (its Fortran, not a binding), so every one of its
        dummies is supplied: its inputs from the `DM_OUTPUT_FROM` plan or by name from the
        wrapper's own arguments, its outputs into the local sizes, its ierr as ierr.
        """
        computed = [a for a in taken if is_computed(a)]
        if not computed:
            return

        groups: dict[str, list[Argument]] = {}
        for argument in computed:
            plan = foo.argument(argument.name).roles.computed_from
            groups.setdefault(plan.producer.name, []).append(argument)

        reports_error = False
        for consumer_arguments in groups.values():
            plan = foo.argument(consumer_arguments[0].name).roles.computed_from
            producer = plan.producer
            produced = {
                foo.argument(a.name).roles.computed_from.output.name.lower(): a.name
                for a in consumer_arguments
            }
            supplied = {pi.name.lower(): pi for pi in plan.inputs}

            actuals = []
            for dummy in producer.arguments:
                lowered = dummy.name.lower()
                if lowered == self.conventions.error_arg.lower():
                    actuals.append((dummy.name, self.conventions.error_arg))
                    reports_error = True
                elif lowered in produced:
                    actuals.append((dummy.name, produced[lowered]))
                elif lowered in supplied:
                    supply = supplied[lowered]
                    value = (
                        supply.argument
                        if supply.argument is not None
                        else _fortran_literal(supply.constant, dummy)
                    )
                    actuals.append((dummy.name, self._actual(value, materialized)))
                elif alloc.argument(dummy.name) is not None:
                    # a producer input the wrapper carries under the same name -- e.g. an
                    # extent the producer's own binding would derive, which a direct call
                    # must still supply
                    actuals.append((dummy.name, self._actual(dummy.name, materialized)))
            self._call(writer, producer.name, actuals)

        if reports_error:
            writer.line(f"if (is_err({self.conventions.error_arg})) return")

    def _actual(self, name: str, materialized) -> str:
        """The name a producer call passes: the materialised local where there is one."""
        if materialized and name in materialized:
            return self._materialized_name(name)
        return name

    def _allocations(self, writer: Writer, taken) -> None:
        for argument in taken:
            if argument.is_array:
                extents = ", ".join(argument.dimension.extents)
                writer.line(f"M_ALLOCATE({argument.name}({extents}))")

    def _permutations(self, writer: Writer, taken) -> None:
        for argument in taken:
            if is_permutation(argument, self.conventions) and not is_temporary(
                argument, self.conventions
            ):
                base = argument.name[: -len(self.conventions.perm_suffix)]
                writer.line(f"call init_perm({argument.name})")
                writer.line(f"call sort_array_heapsort({base}, {argument.name})")

    # -- the kernel call --------------------------------------------------------

    def _kernel_call(self, writer: Writer, foo: Procedure, module: Module) -> None:
        kernel = self._kernel_name(foo, module)
        arguments = self._kernel_call_arguments(foo, module)
        mode_fix = self._mode_fix(foo)
        # the validating wrapper carries exactly the kernel arguments this variant supplies
        # (as dummies, or -- in the allocating body -- as the locals it prepared); anything
        # else is a mode fixed to its parameter, or an argument another mode takes and this
        # one omits (an absent optional)
        present = {a.name.lower() for a in foo.arguments}
        actuals = []
        for argument in arguments:
            lowered = argument.name.lower()
            if mode_fix is not None and lowered == mode_fix.argument.lower():
                actuals.append((argument.name, mode_fix.parameter))
            elif lowered in present:
                actuals.append((argument.name, argument.name))
        self._call(writer, kernel, actuals)
        if self._kernel_declares_error(foo, module):
            self._clear_arg_pos(writer)

    def _clear_arg_pos(self, writer: Writer) -> None:
        """Drop the position an `ierr` carried in from a hand-written procedure.

        It is that procedure's own numbering -- or a private helper's, further down, since a
        position propagates unchanged through every call that does not rewrite it -- so it
        names nothing in this wrapper's dummy list. The code itself is still true.

        Safe because `ierr` is provably OK when the call is made: validation ends with
        `if (is_err(ierr)) return`, and so do the prologue and recommend calls. Nothing this
        clears can be a position the wrapper itself set.
        """
        writer.line(f"call clear_err_arg_pos({self.conventions.error_arg})")

    def _call(self, writer: Writer, name: str, actuals) -> None:
        writer.line(f"call {name}(&")
        with writer.indent():
            for index, (dummy, value) in enumerate(actuals):
                separator = "&" if index == len(actuals) - 1 else ",&"
                writer.line(f"{dummy} = {value}{separator}")
        writer.line(")")


def _fortran_literal(value: object, argument: Argument) -> str:
    """A constant producer input, rendered back as a kinded Fortran literal.

    The kind is the producer dummy's own, so `1` handed to an `integer(int32)` input comes
    out `1_int32` and `.false.` for a logical.
    """
    if isinstance(value, bool):
        return ".true." if value else ".false."
    kind = argument.type.kind
    if isinstance(value, int):
        return f"{value}_{kind}" if kind else str(value)
    if isinstance(value, float):
        return f"{value!r}_{kind}" if kind else repr(value)
    return str(value)
