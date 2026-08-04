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

from ..config import CONVENTIONS, Conventions
from ..ir.entities import Argument, Module, Procedure
from ..ir.types import BaseType
from ..render import Writer
from ..synthesize import is_computed, is_permutation, is_taken_over, is_temporary
from .doc_ford import render_doc
from .fortran_c import _EXTENT_IDENTIFIER_RE, _chunks, _product


class FortranWrapperEmitter:
    def __init__(
        self, conventions: Conventions = CONVENTIONS, macros_header: str = "src/macros.h"
    ):
        self.conventions = conventions
        self.macros_header = macros_header

    # -- module -----------------------------------------------------------------

    def module(self, module: Module) -> str:
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

    def _kernel_module(self, module: Module) -> str:
        return f"{module.name}{self.conventions.kernel_suffix}"

    def _module_uses(self, writer: Writer, module: Module) -> None:
        kernel_module = self._kernel_module(module)
        producers = self._producers(module)

        # the kernels, plus any recommend routine that lives in the kernel module
        kernel_names = {self._kernel_name(p, module) for p in module}
        kernel_names |= producers.get(kernel_module, set())
        for chunk in _chunks(sorted(kernel_names), 4):
            writer.line(f"use {kernel_module}, only: {', '.join(chunk)}")

        # recommend routines an allocating wrapper calls that live elsewhere
        for producer_module, names in sorted(producers.items()):
            if producer_module == kernel_module:
                continue  # merged into the kernel import above
            for chunk in _chunks(sorted(names), 4):
                writer.line(f"use {producer_module}, only: {', '.join(chunk)}")

        kinds = sorted({a.type.kind for p in module for a in p.arguments if a.type.kind})
        for chunk in _chunks(kinds, 6):
            writer.line(f"use, intrinsic :: iso_fortran_env, only: {', '.join(chunk)}")

        if self._has_permutations(module):
            writer.line("use f42_utils, only: init_perm, sort_array_heapsort")

        for chunk in _chunks(self._error_imports(module), 4):
            writer.line(f"use tox_errors, only: {', '.join(chunk)}")

    def _error_imports(self, module: Module) -> list[str]:
        names = {"set_ok", "is_err"}
        for procedure in module:
            for argument in procedure.arguments:
                validator = self._validator_for(argument)
                if validator is not None:
                    names.add(validator)
        if any(p.is_alloc_variant for p in module):
            # M_ALLOCATE reports through set_err / ERR_ALLOC_FAIL
            names |= {"set_err", "ERR_ALLOC_FAIL"}
        ordered = ["set_ok", "is_err"]
        return ordered + sorted(names - set(ordered))

    def _producers(self, module: Module) -> dict[str, set[str]]:
        producers: dict[str, set[str]] = {}
        for alloc in module:
            if not alloc.is_alloc_variant:
                continue
            foo = self._sibling(module, alloc)
            for argument in self._kernel_arguments(foo):
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
            if alloc.is_alloc_variant
            for a in self._kernel_arguments(self._sibling(module, alloc))
        )

    # -- subroutine -------------------------------------------------------------

    def subroutine(self, procedure: Procedure, module: Module) -> str:
        writer = Writer()
        kernel = self._kernel_name(procedure, module)
        link = f"[[{self._kernel_module(module)}(module):{kernel}]]"
        summary = (
            f"Allocates its work arrays, then calls {link}."
            if procedure.is_alloc_variant
            else f"Validates its inputs, then calls {link}."
        )
        writer.block(render_doc(procedure.doc, kind="procedure", summary=summary))
        self._signature(writer, procedure)
        with writer.indent():
            if procedure.is_alloc_variant:
                self._allocating_body(writer, procedure, module)
            else:
                self._validating_body(writer, procedure, module)
        writer.line(f"end subroutine {procedure.name}")
        return writer.render()

    def _validating_body(
        self, writer: Writer, foo: Procedure, module: Module
    ) -> None:
        self._declarations(writer, foo.arguments)
        writer.blank()
        self._validation(writer, foo)
        writer.blank()
        self._kernel_call(writer, self._kernel_name(foo, module), self._kernel_arguments(foo))

    def _allocating_body(
        self, writer: Writer, alloc: Procedure, module: Module
    ) -> None:
        foo = self._sibling(module, alloc)
        kernel_arguments = self._kernel_arguments(foo)
        taken = [a for a in kernel_arguments if is_taken_over(a, self.conventions)]

        self._declarations(writer, alloc.arguments)
        self._locals(writer, taken)
        writer.blank()
        self._validation(writer, alloc)
        writer.blank()
        self._recommend_calls(writer, foo, alloc, taken)
        self._allocations(writer, taken)
        self._permutations(writer, taken)
        writer.blank()
        self._kernel_call(writer, self._kernel_name(foo, module), kernel_arguments)

    # -- names ------------------------------------------------------------------

    def _base(self, procedure: Procedure) -> str:
        name = procedure.name
        for suffix in (self.conventions.alloc_suffix, self.conventions.validating_suffix):
            if suffix and name.lower().endswith(suffix):
                return name[: -len(suffix)]
        return name

    def _kernel_name(self, procedure: Procedure, module: Module) -> str:
        return f"{self._base(procedure)}{self.conventions.kernel_suffix}"

    def _sibling(self, module: Module, alloc: Procedure) -> Procedure:
        return module.procedure(self._base(alloc) + self.conventions.validating_suffix)

    def _kernel_arguments(self, foo: Procedure) -> list[Argument]:
        """The kernel's own arguments: everything the validating wrapper carries but ierr."""
        error = self.conventions.error_arg.lower()
        return [a for a in foo.arguments if a.name.lower() != error]

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
        writer.line(f"if (is_err({error})) return")

    def _validator_for(self, argument: Argument) -> str | None:
        """The `tox_errors` validator this argument needs, or None."""
        if not argument.intent.is_input:
            return None  # an output carries no value to check
        if argument.roles is not None and argument.roles.is_extent:
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

    def _recommend_calls(
        self, writer: Writer, foo: Procedure, alloc: Procedure, taken
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
                    actuals.append((dummy.name, value))
                elif alloc.argument(dummy.name) is not None:
                    # a producer input the wrapper carries under the same name -- e.g. an
                    # extent the producer's own binding would derive, which a direct call
                    # must still supply
                    actuals.append((dummy.name, dummy.name))
            self._call(writer, producer.name, actuals)

        if reports_error:
            writer.line(f"if (is_err({self.conventions.error_arg})) return")

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

    def _kernel_call(self, writer: Writer, kernel: str, arguments) -> None:
        # the kernel has no ierr (it does no reporting), so it is not passed on
        self._call(writer, kernel, [(a.name, a.name) for a in arguments])

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
