"""Emitting the Fortran C wrapper modules.

One `<module>_c.F90` per module that exports anything, wrapped in `#ifndef
NO_C_INTERFACE` so the whole C interface -- and with it the `safeguard` module -- can be
compiled out with a single directive.

The wrapper's job, in order:

1. check `ierr` itself, since it is the only channel for reporting anything
2. set it ok, because a procedure that declares no `ierr` will never set it
3. null-check the arguments, in an order `c_loc` tolerates (see `CWrapper.validation_order`)
4. convert what C cannot express: strings, logicals, mode strings
5. call the procedure
6. convert the outputs back
"""

from __future__ import annotations

import re

from ..abi.model import CArgument, Conversion, CWrapper, CWrapperModule, Origin
from ..config import CONVENTIONS, Conventions
from ..ir.types import BaseType, Intent
from ..render import Writer
from .doc_ford import render_doc

#: A Fortran identifier inside an extent expression. A leading digit is required to be
#: absent so an integer literal's kind suffix (`0_int32` -> `_int32`) is the only spurious
#: match, and that never collides with a real argument name.
_EXTENT_IDENTIFIER_RE = re.compile(r"[A-Za-z_]\w*")

#: Suffix of the local a converted argument is passed to the callee as
CONVERTED_SUFFIX = "_f"
#: Suffix of the local holding a mode string's integer parameter
MODE_SUFFIX = "_mode_f"
#: Suffix of the pointer standing in for a nullable optional
POINTER_SUFFIX = "_p"


def _product(extents) -> str:
    """The product of some extents, each parenthesised.

    An extent may be an expression (`n_timepoints - 1`), and `*` binds tighter than `-`
    in both C++ and Fortran: joining raw gives `n_timepoints - 1 * n_factors`, which is
    a different number and, for a buffer size, a heap overrun waiting to happen.
    """
    parts = [e if e.isidentifier() or e.isdigit() else f"({e})" for e in extents]
    return " * ".join(parts) if parts else "1"

class FortranCEmitter:
    def __init__(self, conventions: Conventions = CONVENTIONS,
                 macros_header: str = "src/macros.h"):
        self.conventions = conventions
        self.macros_header = macros_header

    # -- module -----------------------------------------------------------------

    def module(self, module: CWrapperModule) -> str:
        writer = Writer()
        writer.line("#ifndef NO_C_INTERFACE")
        writer.line(f"#include <{self.macros_header}>")
        writer.blank()

        link = f"[[{module.stripped_name}(module)]]"
        writer.block(render_doc(module.doc, kind="module",
                                summary=f"C-wrappers for {link}"))
        writer.line(f"module {module.name}")
        with writer.indent():
            self._module_uses(writer, module)
            # Not a bare `implicit none`: that constrains only variables, so a call to
            # a procedure that does not exist compiles as an implicit external and fails
            # at link time. Generated code is exactly where that must not be possible --
            # a typo here is the generator's fault, and should surface at its output
            # rather than in someone's build log.
            writer.line("M_IMPLICIT_NONE")
            writer.line("private")
            writer.blank()
            for wrapper in module:
                writer.line(f"public :: {wrapper.name}")
        writer.blank()
        writer.line("contains")
        writer.blank()

        for wrapper in module:
            with writer.indent():
                writer.block(self.wrapper(wrapper))
            writer.blank()

        writer.line(f"end module {module.name}")
        writer.line("#endif")
        return writer.render(trailing_newline=True)

    def _module_uses(self, writer: Writer, module: CWrapperModule) -> None:
        # safeguard lives only here: nothing outside the generated C interface needs it,
        # which is the point of giving the wrappers their own modules
        writer.line("use safeguard")

        kinds = sorted({a.type.kind for w in module for a in w if a.type.kind})
        kinds = sorted(set(kinds) | {"c_loc", "c_associated"})
        for chunk in _chunks(kinds, 6):
            writer.line(f"use, intrinsic :: iso_c_binding, only: {', '.join(chunk)}")

        conversions = sorted({
            name for wrapper in module for name in _conversion_helpers(wrapper)
        })
        if conversions:
            for chunk in _chunks(conversions, 3):
                writer.line(f"use tox_conversions, only: {', '.join(chunk)}")

        errors = ["set_ok", "set_err", "is_err", "ERR_POINTER_NULL"]
        if any(a.conversion is Conversion.MODE for w in module for a in w):
            errors.append("ERR_INVALID_INPUT")
        if _needs_allocation(module):
            errors.append("ERR_ALLOC_FAIL")
        writer.line(f"use tox_errors, only: {', '.join(errors)}")

    # -- wrapper ----------------------------------------------------------------

    def wrapper(self, wrapper: CWrapper) -> str:
        writer = Writer()
        procedure = wrapper.procedure
        link = f"[[{wrapper.module_name}(module):{procedure.name}({procedure.kind})]]"
        writer.block(
            render_doc(wrapper.doc, kind="procedure", summary=f"C-wrapper for {link}")
        )

        writer.line(f"subroutine {wrapper.name}(&")
        with writer.indent(2):
            names = [a.name for a in wrapper.arguments]
            for name in names[:-1]:
                writer.line(f"{name},&")
            writer.line(f"{names[-1]}&")
        with writer.indent():
            writer.line(f') bind(C, name="{wrapper.name}")')

        with writer.indent():
            writer.line(f"use {wrapper.module_name}, only: {procedure.name}")
            self._mode_uses(writer, wrapper)
            writer.blank()
            self._declarations(writer, wrapper)
            self._locals(writer, wrapper)
            writer.blank()
            self._validation(writer, wrapper)
            writer.blank()
            self._input_conversions(writer, wrapper)
            self._call(writer, wrapper)
            self._output_conversions(writer, wrapper)

        writer.line(f"end subroutine {wrapper.name}")
        return writer.render()

    def _mode_uses(self, writer: Writer, wrapper: CWrapper) -> None:
        by_module: dict[str, set[str]] = {}
        for argument in wrapper:
            if argument.mode is None:
                continue
            for value in argument.mode.values:
                by_module.setdefault(value.module, set()).add(value.parameter)
        for module_name in sorted(by_module):
            names = ", ".join(sorted(by_module[module_name]))
            writer.line(f"use {module_name}, only: {names}")

    def _declarations(self, writer: Writer, wrapper: CWrapper) -> None:
        for argument in self._declaration_order(wrapper):
            writer.line(self._declaration(argument))
            if argument.doc:
                with writer.indent():
                    writer.block(render_doc(argument.doc, kind="argument"))

    @staticmethod
    def _declaration_order(wrapper: CWrapper) -> list[CArgument]:
        """Declarations ordered so an extent is typed before the array that uses it.

        Referring to a symbol that is typed further down the specification part is a GNU
        extension, not standard Fortran: gfortran -std=f2018 rejects `dimension(n)` above
        `integer :: n` outright. So anything named in someone's extents is hoisted, and
        everything else keeps the order the author wrote.

        An extent may be an expression -- `max(0_int32, n_timepoints - 1)`,
        `sum(reps_per_tissue)` -- so identifiers are pulled out of it rather than comparing
        the whole string to an argument name; otherwise the argument the expression depends
        on is not recognised as referenced and is left declared below the array using it.
        """
        referenced = {
            identifier.lower()
            for argument in wrapper
            for extent in argument.dimension.extents
            for identifier in _EXTENT_IDENTIFIER_RE.findall(extent)
        }
        extents = [a for a in wrapper if a.name.lower() in referenced]
        rest = [a for a in wrapper if a.name.lower() not in referenced]
        return [*extents, *rest]

    def _declaration(self, argument: CArgument) -> str:
        attributes = [f"intent({argument.intent.value})"]
        if argument.dimension:
            attributes.insert(0, f"dimension({', '.join(argument.dimension.extents)})")
        if argument.optional:
            # An OPTIONAL dummy of a bind(C) procedure is absent exactly when C passes a
            # null pointer (TS 29113, in F2018). The C prototype stays a plain pointer,
            # so this costs the ABI nothing, and the wrapper can hand the argument
            # straight to the callee: an optional associated with an absent optional is
            # itself absent. No branch, no descriptor, no pointer games.
            attributes.append("optional")
        else:
            # c_loc needs a target. An optional gets none: it is never null-checked,
            # and c_loc of an absent argument is not allowed anyway.
            attributes.append("target")
        return f"{argument.type}, {', '.join(attributes)} :: {argument.name}"

    def _locals(self, writer: Writer, wrapper: CWrapper) -> None:
        for argument in wrapper:
            for declaration in self._local_declarations(argument):
                writer.line(declaration)

    def _local_declarations(self, argument: CArgument) -> list[str]:
        source = argument.source
        if argument.conversion is Conversion.MODE:
            return [f"integer(int32) :: {argument.name}{MODE_SUFFIX}"]

        # A nullable optional that needs a converted local cannot use an automatic one:
        # when C passes null the argument is absent, and the local must then be *absent*
        # to the callee too. An unallocated allocatable actual is an absent optional dummy
        # (F2018), so the converted local is allocatable and simply left unallocated. Its
        # shape is deferred -- one `:` per rank of the source, so a rank-2 mask does not
        # come out declared rank-1.
        nullable = argument.optional
        deferred = (
            ", ".join(":" for _ in source.dimension.extents)
            if source is not None and source.dimension else ""
        )

        if argument.conversion is Conversion.LOGICAL:
            base = "logical"
            if nullable:
                base += f", dimension({deferred}), allocatable" if deferred else ", allocatable"
            elif source.dimension:
                base += f", dimension({', '.join(source.dimension.extents)})"
            return [f"{base} :: {argument.name}{CONVERTED_SUFFIX}"]

        if argument.conversion is Conversion.CHARACTER:
            # The local takes the callee's declared length, so it is a valid actual
            # argument whatever C sent. A deferred length (`len=*` callee) forces an
            # allocatable local -- and allocatable requires a deferred shape too, so the
            # explicit extents give way to `:` per rank (an explicit-shape allocatable is
            # rejected). An assumed-length *output* array is then allocated before the call.
            length = source.type.length
            allocatable = length.is_assumed or nullable
            declared = "len=:" if length.is_assumed else f"len={length}"
            base = f"character({declared})"
            if allocatable:
                base += ", allocatable"
            if source.dimension:
                base += f", dimension({deferred})" if allocatable else \
                    f", dimension({', '.join(source.dimension.extents)})"
            return [f"{base} :: {argument.name}{CONVERTED_SUFFIX}"]

        return []

    # -- body -------------------------------------------------------------------

    def _validation(self, writer: Writer, wrapper: CWrapper) -> None:
        writer.line("M_CHECK_IERR_NON_NULL")
        # A procedure that declares no ierr never sets one, so the wrapper must
        writer.line(f"call set_ok({wrapper.error_argument.name})")

        for argument in wrapper.validation_order:
            if argument.is_error:
                continue
            if argument.is_scalar:
                writer.line(f"M_CHECK_NON_NULL({argument.name})")
            else:
                writer.line(
                    f"M_CHECK_ARRAY_NON_NULL({argument.name}, {self._size_of(argument)})"
                )

    def _size_of(self, argument: CArgument) -> str:
        """The element count that guards `c_loc` of an array."""
        if argument.shape_arg is not None:
            # Its shape travels separately, so the count is that array's product. The
            # ordering has already made the shape argument safe to read.
            return f"product({argument.shape_arg})"
        return _product(argument.size_extents)

    def _input_conversions(self, writer: Writer, wrapper: CWrapper) -> None:
        wrote = False
        for argument in wrapper:
            if not argument.needs_conversion:
                continue
            if argument.intent is Intent.OUT and argument.conversion is not Conversion.MODE:
                # No value flows in, but an allocatable local still has to be allocated
                # before the call so the intent(out) callee has something to write into:
                # a nullable optional (allocated only when present), or an assumed-length
                # character array (allocated at the C-supplied string length).
                block = self._allocate_output_local(argument, wrapper)
                if block:
                    writer.block(block)
                    wrote = True
                continue
            block = self._input_conversion(argument, wrapper)
            if block:
                writer.block(block)
                wrote = True
        if wrote:
            writer.blank()

    def _allocate_output_local(self, argument: CArgument, wrapper: CWrapper) -> str:
        source = argument.source
        assumed_char = (
            argument.conversion is Conversion.CHARACTER and source.type.length.is_assumed
        )
        # An automatic local (fixed-length character, or a non-optional logical) is sized by
        # its declaration; only an allocatable one needs allocating here.
        if not argument.optional and not assumed_char:
            return ""
        # The C-visible extents (real names, synthesised for an assumed-shape source rather
        # than the source's own `:`). For a character the leading extent is the string
        # length, which `character(len=...)` carries, so it is dropped from the shape.
        extents = list(argument.dimension.extents)
        if argument.conversion is Conversion.CHARACTER and extents:
            extents = extents[1:]
        target = f"{argument.name}{CONVERTED_SUFFIX}"
        if extents:
            target += f"({', '.join(extents)})"
        if argument.conversion is Conversion.CHARACTER:
            allocate = f"allocate(character(len={self._strlen_name(argument, wrapper)}) :: {target})"
        else:
            allocate = f"allocate({target})"
        if argument.optional:
            return f"if (present({argument.name})) {allocate}"
        return allocate

    @staticmethod
    def _strlen_name(argument: CArgument, wrapper: CWrapper) -> str:
        """The C argument carrying this character's length (see abi `_character_length`)."""
        for other in wrapper:
            if other.origin is Origin.STRLEN and other.sizes == argument.name:
                return other.name
        return f"{argument.name}_strlen"

    def _conversion_source(self, argument: CArgument) -> str:
        """The wrapper's own array, sliced when its length travels separately.

        An argument with a shape argument is assumed-size in the wrapper, and an
        assumed-size array can be neither assigned whole nor passed to an assumed-shape
        dummy. Its element count is the product of the shape argument, exactly as when it
        is passed to the callee unconverted.
        """
        if argument.shape_arg is None:
            return argument.name
        if argument.type.is_character:
            # a character buffer carries its length as the leading extent, so the items
            # are the second dimension
            return f"{argument.name}(:, 1:product({argument.shape_arg}))"
        return f"{argument.name}(1:product({argument.shape_arg}))"

    def _input_conversion(self, argument: CArgument, wrapper: CWrapper) -> str:
        writer = Writer()
        name = argument.name
        error = wrapper.error_argument.name

        match argument.conversion:
            case Conversion.LOGICAL:
                # logical(c_bool) and default logical differ in kind, so intrinsic
                # assignment does the conversion. A nullable optional is only converted
                # when present -- absent leaves the allocatable local unallocated, which
                # the callee reads as absent.
                source_expression = self._conversion_source(argument)
                if argument.optional:
                    writer.line(
                        f"if (present({name})) {name}{CONVERTED_SUFFIX} = {source_expression}"
                    )
                else:
                    writer.line(f"{name}{CONVERTED_SUFFIX} = {source_expression}")
            case Conversion.CHARACTER:
                writer.block(self._guard_optional(argument, self._character_in(argument, error)))
            case Conversion.MODE:
                writer.block(self._mode_in(argument, error))
        return writer.render()

    def _guard_optional(self, argument: CArgument, block: str) -> str:
        """Wrap a multi-line conversion in `if (present(...)) then` for a nullable optional."""
        if not argument.optional:
            return block
        writer = Writer()
        writer.line(f"if (present({argument.name})) then")
        with writer.indent():
            writer.block(block)
        writer.line("end if")
        return writer.render()

    def _character_in(self, argument: CArgument, error: str) -> str:
        writer = Writer()
        name = argument.name
        source = argument.source
        helper = "c_char_1d_as_string" if source.rank == 0 else "c_char_2d_as_string"
        from_buffer = self._conversion_source(argument)

        if source.type.length.is_assumed:
            writer.line(f"call {helper}({from_buffer}, {name}{CONVERTED_SUFFIX}, {error})")
        else:
            # The callee's length is fixed, so convert into a deferred-length local and
            # let assignment pad or truncate to exactly that length
            writer.line(f"block")
            with writer.indent():
                declared = "character(len=:), allocatable"
                if source.rank:
                    declared += ", dimension(:)"
                writer.line(f"{declared} :: converted")
                writer.line(f"call {helper}({name}, converted, {error})")
                writer.line(f"if (is_err({error})) return")
                writer.line(f"{name}{CONVERTED_SUFFIX} = converted")
            writer.line("end block")
            return writer.render()

        writer.line(f"if (is_err({error})) return")
        return writer.render()

    def _mode_in(self, argument: CArgument, error: str) -> str:
        """Recover the mode string from the c_char buffer and map it to its parameter.

        The local lives in a block construct: it exists only for this conversion, and a
        block keeps it next to the code that uses it instead of in a declaration section
        away from any context.
        """
        name = argument.name
        local = f"{name}{CONVERTED_SUFFIX}"

        writer = Writer()
        writer.line("block")
        with writer.indent():
            writer.line(f"character(len=:), allocatable :: {local}")
            writer.line(f"call c_char_1d_as_string({name}, {local}, {error})")
            writer.line(f"if (is_err({error})) return")
            writer.blank()
            writer.line(f"select case ({local})")
            with writer.indent():
                for value in argument.mode.values:
                    writer.line(f'case ("{value.string}")')
                    with writer.indent():
                        writer.line(f"{name}{MODE_SUFFIX} = {value.parameter}")
                writer.line("case default")
                with writer.indent():
                    writer.line(f"call set_err({error}, ERR_INVALID_INPUT)")
                    writer.line("return")
            writer.line("end select")
        writer.line("end block")
        return writer.render()

    def _call(self, writer: Writer, wrapper: CWrapper) -> None:
        procedure = wrapper.procedure
        actuals = []
        for argument in wrapper:
            if argument.origin in (Origin.EXTENT, Origin.STRLEN):
                # Invented for C; the callee learns the size from its own declaration
                if argument.sizes and _is_declared_extent(argument, wrapper):
                    actuals.append((argument.name, argument.name))
                continue
            if argument.origin is Origin.ERROR and argument.source is None:
                continue
            if argument.origin is Origin.RESULT:
                continue
            actuals.append((argument.name, self._actual(argument)))

        result = next((a for a in wrapper if a.origin is Origin.RESULT), None)
        prefix = "call" if result is None else f"{self._actual(result)} ="

        writer.line(f"{prefix} {procedure.name}(&")
        with writer.indent():
            for index, (dummy, actual) in enumerate(actuals):
                separator = "&" if index == len(actuals) - 1 else ",&"
                writer.line(f"{dummy} = {actual}{separator}")
        writer.line(")")

    def _actual(self, argument: CArgument) -> str:
        if argument.conversion is Conversion.MODE:
            return f"{argument.name}{MODE_SUFFIX}"
        if argument.needs_conversion:
            return f"{argument.name}{CONVERTED_SUFFIX}"
        if argument.shape_arg is not None:
            # The wrapper has it as assumed size, which cannot be passed to the callee's
            # assumed-shape dummy. Slicing to its actual element count gives the shape
            # back. The count is the product of the shape argument, not its size -- its
            # size is the rank.
            return f"{argument.name}(1:product({argument.shape_arg}))"
        return argument.name

    def _output_conversions(self, writer: Writer, wrapper: CWrapper) -> None:
        blocks = []
        for argument in wrapper:
            if not argument.needs_conversion or argument.intent is Intent.IN:
                continue
            if argument.conversion is Conversion.MODE:
                continue
            blocks.append(self._output_conversion(argument))
        if blocks:
            writer.blank()
            for block in blocks:
                writer.block(block)

    def _output_conversion(self, argument: CArgument) -> str:
        name = argument.name
        source = argument.source
        # the same slice the input direction needs: an assumed-size array can be neither
        # assigned whole nor passed to an assumed-shape dummy
        target = self._conversion_source(argument)
        if argument.conversion is Conversion.LOGICAL:
            assignment = f"{target} = {name}{CONVERTED_SUFFIX}"
        else:
            helper = "string_as_c_char_1d" if source.rank == 0 else "string_as_c_char_2d"
            assignment = f"call {helper}({name}{CONVERTED_SUFFIX}, {target})"
        # A nullable optional is only written back when present -- absent means C passed
        # null, so there is nothing to convert into.
        if argument.optional:
            return f"if (present({name})) {assignment}"
        return assignment


def _is_declared_extent(argument: CArgument, wrapper: CWrapper) -> bool:
    """Whether the callee has an argument of this name to receive."""
    procedure = wrapper.procedure
    return procedure.argument(argument.name) is not None


def _conversion_helpers(wrapper: CWrapper) -> set[str]:
    helpers: set[str] = set()
    for argument in wrapper:
        if argument.conversion is Conversion.CHARACTER:
            rank = argument.source.rank
            helpers.add("c_char_1d_as_string" if rank == 0 else "c_char_2d_as_string")
            if argument.intent is not Intent.IN:
                helpers.add("string_as_c_char_1d" if rank == 0 else "string_as_c_char_2d")
        elif argument.conversion is Conversion.MODE:
            helpers.add("c_char_1d_as_string")
    return helpers


def _needs_allocation(module: CWrapperModule) -> bool:
    return any(
        a.conversion is Conversion.CHARACTER for wrapper in module for a in wrapper
    )


def _chunks(items, size):
    items = list(items)
    for start in range(0, len(items), size):
        yield items[start : start + size]
