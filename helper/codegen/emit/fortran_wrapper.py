"""Emitting the generated wrapper modules.

One `tox_<name>.F90` per kernel module, holding the validating wrapper `foo` (and, later,
the allocating wrapper `foo_alloc`) for every kernel. Unlike the C wrappers these are
ordinary library sources -- no `bind(C)`, no `NO_C_BINDING` guard, no pointer marshalling.

A validating wrapper's job, in order:

1. set `ierr` ok, because the kernel it calls has none and reports nothing
2. validate every input against the range its documentation states -- and, for reals,
   against the framework's default finiteness contract
3. bail if anything failed
4. call the kernel with the arguments unchanged

Finiteness is the default: a real input is checked for NaN/Inf whether or not it carries a
`DM_MIN`/`DM_MAX`, unless it opts out with `DM_ALLOW_NAN` / `DM_ALLOW_INFINITE`.
"""

from __future__ import annotations

from ..config import CONVENTIONS, Conventions
from ..ir.entities import Argument, Module, Procedure
from ..ir.types import BaseType
from ..render import Writer
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
        """Emit `module tox_<name>` with a wrapper per procedure it holds."""
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

    def _kernel_name(self, procedure: Procedure) -> str:
        return f"{procedure.name}{self.conventions.kernel_suffix}"

    def _module_uses(self, writer: Writer, module: Module) -> None:
        kernel_module = self._kernel_module(module)
        kernel_names = sorted({self._kernel_name(p) for p in module})
        for chunk in _chunks(kernel_names, 4):
            writer.line(f"use {kernel_module}, only: {', '.join(chunk)}")

        kinds = sorted({a.type.kind for p in module for a in p.arguments if a.type.kind})
        for chunk in _chunks(kinds, 6):
            writer.line(f"use, intrinsic :: iso_fortran_env, only: {', '.join(chunk)}")

        errors = self._error_imports(module)
        for chunk in _chunks(errors, 4):
            writer.line(f"use tox_errors, only: {', '.join(chunk)}")

    def _error_imports(self, module: Module) -> list[str]:
        names = {"set_ok", "is_err"}
        for procedure in module:
            for argument in procedure.arguments:
                validator = self._validator_for(argument)
                if validator is not None:
                    names.add(validator)
        # a stable order, with the always-present two first
        ordered = ["set_ok", "is_err"]
        ordered += sorted(names - set(ordered))
        return ordered

    # -- subroutine -------------------------------------------------------------

    def subroutine(self, procedure: Procedure, module: Module) -> str:
        writer = Writer()
        kernel = self._kernel_name(procedure)
        link = f"[[{self._kernel_module(module)}(module):{kernel}]]"
        writer.block(
            render_doc(
                procedure.doc,
                kind="procedure",
                summary=f"Validates its inputs, then calls {link}.",
            )
        )
        self._signature(writer, procedure)
        with writer.indent():
            self._declarations(writer, procedure)
            writer.blank()
            self._validation(writer, procedure)
            writer.blank()
            self._kernel_call(writer, procedure, kernel)
        writer.line(f"end subroutine {procedure.name}")
        return writer.render()

    def _signature(self, writer: Writer, procedure: Procedure) -> None:
        writer.line(f"subroutine {procedure.name}(&")
        names = [a.name for a in procedure.arguments]
        with writer.indent(2):
            for name in names[:-1]:
                writer.line(f"{name},&")
            writer.line(f"{names[-1]}&")
        with writer.indent():
            writer.line(")")

    def _declarations(self, writer: Writer, procedure: Procedure) -> None:
        for argument in self._declaration_order(procedure.arguments):
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

    # -- body -------------------------------------------------------------------

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

    def _kernel_call(self, writer: Writer, procedure: Procedure, kernel: str) -> None:
        # The kernel has no ierr (it does no reporting), so it is not passed on.
        error = self.conventions.error_arg
        actuals = [a for a in procedure.arguments if a.name.lower() != error.lower()]
        writer.line(f"call {kernel}(&")
        with writer.indent():
            for index, argument in enumerate(actuals):
                separator = "&" if index == len(actuals) - 1 else ",&"
                writer.line(f"{argument.name} = {argument.name}{separator}")
        writer.line(")")
