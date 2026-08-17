"""Emitting the Fortran C wrapper modules.

One `<module>_c.F90` per module that exports anything, wrapped in `#ifndef
NO_C_BINDING` so the whole C binding -- and with it the `f42_safeguard` module -- can be
compiled out with a single directive.

The wrapper's job, in order:

1. check `ierr` itself, since it is the only channel for reporting anything
2. set it ok, because a procedure that declares no `ierr` will never set it
3. null-check the arguments, in an order `c_loc` tolerates (see `CWrapper.validation_order`)
4. prepare what C cannot express: a pointer view of a string buffer, a logical work
   array, the integer a mode string names
5. call the procedure
6. convert the outputs back -- only a logical, since a string view *is* the caller's buffer
"""

from __future__ import annotations

import re

from ..abi.model import CArgument, Conversion, CWrapper, CWrapperModule, Origin
from ..config import CONVENTIONS, Conventions
from ..ir.types import Intent
from ..render import Writer
from .doc_ford import render_doc
from .doc_literals import render_mode_default

#: A Fortran identifier inside an extent expression. A leading digit is required to be
#: absent so an integer literal's kind suffix (`0_int32` -> `_int32`) is the only spurious
#: match, and that never collides with a real argument name.
_EXTENT_IDENTIFIER_RE = re.compile(r"[A-Za-z_]\w*")

#: Suffix of the local a converted argument is passed to the callee as
CONVERTED_SUFFIX = "_f"
#: Suffix of the local holding a mode string's integer parameter
MODE_SUFFIX = "_mode_f"


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
        writer.line("#ifndef NO_C_BINDING")
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
        # f42_safeguard lives only here: nothing outside the generated C binding needs it,
        # which is the point of giving the wrappers their own modules
        writer.line("use f42_safeguard")

        kinds = sorted({a.type.kind for w in module for a in w if a.type.kind})
        kinds = set(kinds) | {"c_loc", "c_associated"}
        if any(a.conversion is Conversion.CHARACTER for w in module for a in w):
            # a character crosses as a pointer view of the caller's buffer, remapped here.
            # `M_IMPLICIT_NONE` is `implicit none (type, external)`, so an unimported
            # intrinsic is a compile error rather than a silent implicit external.
            kinds.add("c_f_pointer")
        kinds = sorted(kinds)
        for chunk in _chunks(kinds, 6):
            writer.line(f"use, intrinsic :: iso_c_binding, only: {', '.join(chunk)}")

        conversions = sorted({
            name for wrapper in module for name in _conversion_helpers(wrapper)
        })
        if conversions:
            for chunk in _chunks(conversions, 3):
                writer.line(f"use tox_conversions, only: {', '.join(chunk)}")

        # what the body will actually name -- an unused `only:` name is not a compile error
        # in any Fortran, so a fixed list goes stale silently
        errors = ["set_ok", "set_err"]
        errors.append("ERR_POINTER_NULL")
        if any(a.conversion is Conversion.MODE for w in module for a in w):
            errors.append("ERR_INVALID_INPUT")
        # Only a converted argument ever reaches `_allocate_local` in the body, so asking it
        # about the others answered a question the wrapper never puts -- and left
        # `ERR_ALLOC_FAIL` imported by modules that allocate nothing.
        if any(a.needs_conversion and self._allocate_local(a, w) for w in module for a in w):
            errors.append("ERR_ALLOC_FAIL")  # M_ALLOCATE names it on the failure path
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
                # C takes the mode as a string, so the documented integer default names a
                # value this layer's caller cannot pass. The plain Fortran wrapper keeps it.
                text = render_doc(argument.doc, kind="argument")
                if argument.mode_default is not None:
                    text = render_mode_default(text, argument.mode_default, "fortran")
                with writer.indent():
                    writer.block(text)

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
            if argument.conversion is Conversion.CHARACTER:
                # A character is remapped, not copied, so the wrapper takes `c_loc` of it
                # under `present(...)` -- and `c_loc` needs a TARGET. `optional, target`
                # together are legal on a bind(C) dummy and cost the ABI nothing; without
                # the target both gfortran and ifx reject the c_loc outright.
                attributes.append("target")
        else:
            # c_loc needs a target, for the null check and for the pointer view.
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

        # A logical work array whose size the *caller* chooses is allocatable, and is
        # allocated through `M_ALLOCATE`. An automatic array is the obvious spelling and was
        # the original one, but it is stack storage sized at run time: a 10-million-element
        # mask is 40 MB, which overflows a default stack -- ifx segfaults on it, and gfortran
        # only escapes by quietly moving large automatics to the heap. `M_ALLOCATE` puts it
        # on the heap on every compiler and turns a failure into `ERR_ALLOC_FAIL` returned to
        # the caller, which is the whole point of a wrapper that owns memory.
        #
        # A nullable optional needs it for a second reason: when C passes null the argument
        # is absent, and the local must be *absent* to the callee too. An unallocated
        # allocatable actual is an absent optional dummy (F2018 15.5.2.12), so it is simply
        # left unallocated. Its shape is deferred -- one `:` per rank of the source, so a
        # rank-2 mask does not come out declared rank-1. A character reaches the same end
        # through the other half of that rule: a *disassociated pointer* is absent too.
        nullable = argument.optional
        deferred = (
            ", ".join(":" for _ in source.dimension.extents)
            if source is not None and source.dimension else ""
        )

        if argument.conversion is Conversion.LOGICAL:
            # A scalar stays automatic: four bytes, and its size is not the caller's to pick.
            base = "logical"
            if deferred:
                base += f", dimension({deferred}), allocatable"
            elif nullable:
                base += ", allocatable"
            return [f"{base} :: {argument.name}{CONVERTED_SUFFIX}"]

        if argument.conversion is Conversion.CHARACTER:
            # A string does not get converted at all: the local is a *view* of the caller's
            # buffer, remapped with `c_f_pointer`. `character(len=n), dimension(m)` and
            # `character(len=1), dimension(n, m)` have identical layout -- contiguous,
            # column-major, string i at bytes [(i-1)n+1 .. i*n] -- so the remap is pure
            # reinterpretation, in either direction. See design/c-layer.md.
            #
            # The length is the buffer's leading C extent: the synthesised `<name>_strlen`
            # for a `len=*` callee, the callee's own `len=` when it is fixed. It must NEVER
            # be spelled `len=:`: a deferred-length remap compiles clean on both gfortran
            # and ifx and then yields an empty string on one and a segfault on the other.
            base = f"character(len={argument.dimension.extents[0]}), pointer"
            if source.dimension:
                base += f", dimension({deferred})"
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
            # The local comes first in either direction: an intent(out) callee needs
            # something to write into, and an input conversion needs somewhere to convert
            # into. Nothing below may run before the local exists.
            block = self._prepare_local(argument, wrapper)
            if block:
                writer.block(block)
                wrote = True
            if argument.intent is Intent.OUT and argument.conversion is not Conversion.MODE:
                continue
            block = self._input_conversion(argument, wrapper)
            if block:
                writer.block(block)
                wrote = True
        if wrote:
            writer.blank()

    def _prepare_local(self, argument: CArgument, wrapper: CWrapper) -> str:
        """Give the converted local something to be, before the callee is handed it.

        A character is *remapped* onto the caller's buffer; anything else that needs one is
        allocated. Both run whatever the intent -- an `intent(out)` view has to exist before
        the call, or the callee writes nowhere.
        """
        if argument.conversion is Conversion.CHARACTER:
            return self._remap_local(argument)
        return self._allocate_local(argument, wrapper)

    def _remap_local(self, argument: CArgument) -> str:
        """`c_f_pointer` making the local a view of the caller's own string buffer.

        No copy, no allocation, and nothing that can fail -- so no `is_err` check follows it,
        and the one statement serves every intent. The shape drops the leading C extent,
        which is the string length the pointer's own `character(len=...)` carries; a scalar
        is left with no shape at all.
        """
        name = argument.name
        local = f"{name}{CONVERTED_SUFFIX}"
        if argument.shape_arg is not None:
            # the C extent is `*`; how many strings there really are is the product of the
            # shape that travels beside them
            extents = [f"product({argument.shape_arg})"]
        else:
            extents = list(argument.dimension.extents)[1:]
        shape = f", [{', '.join(extents)}]" if extents else ""
        remap = f"call c_f_pointer(c_loc({name}), {local}{shape})"
        if not argument.optional:
            return remap
        # C passed null -> the argument is absent -> the local must be absent to the callee
        # too, which for a pointer means *disassociated* (F2018 15.5.2.12). An undefined
        # pointer is not a disassociated one, so it is nullified explicitly -- as a
        # statement, never as an initialiser, which would give the local an implicit SAVE.
        return f"nullify({local})\nif (present({name})) {remap}"

    def _allocate_local(self, argument: CArgument, wrapper: CWrapper) -> str:
        """`M_ALLOCATE` for a converted local, or empty where the local is not allocatable.

        Every logical work array sized by something C supplied is allocated here rather than
        declared as an automatic array, so the storage is on the heap on every compiler and a
        failure comes back as `ERR_ALLOC_FAIL` instead of a stack overflow. Only a scalar
        logical stays automatic -- four bytes, and not the caller's size to choose. A
        character allocates nothing at all: it is a view, see `_remap_local`.
        """
        source = argument.source
        if argument.conversion in (Conversion.MODE, Conversion.CHARACTER):
            return ""
        allocatable = (source is not None and bool(source.dimension)) or argument.optional
        if not allocatable:
            return ""
        # The C-visible extents (real names, synthesised for an assumed-shape source rather
        # than the source's own `:`).
        extents = list(argument.dimension.extents)
        if argument.shape_arg is not None:
            # the C extent is `*`; how many items there really are is the product of the
            # shape that travels beside them
            extents = [f"product({argument.shape_arg})"]
        target = f"{argument.name}{CONVERTED_SUFFIX}"
        if extents:
            target += f"({', '.join(extents)})"
        allocate = f"M_ALLOCATE({target})"
        if not argument.optional:
            return allocate
        # M_ALLOCATE expands to several statements, so the guard has to be a block rather
        # than the one-line `if (present(x)) allocate(...)` this replaced.
        return f"if (present({argument.name})) then\n    {allocate}\nend if"

    def _conversion_source(self, argument: CArgument) -> str:
        """The wrapper's own array, sliced when its length travels separately.

        An argument with a shape argument is assumed-size in the wrapper, and an
        assumed-size array can be neither assigned whole nor passed to an assumed-shape
        dummy. Its element count is the product of the shape argument, exactly as when it
        is passed to the callee unconverted.
        """
        if argument.shape_arg is None:
            return argument.name
        return f"{argument.name}(1:product({argument.shape_arg}))"

    def _input_conversion(self, argument: CArgument, wrapper: CWrapper) -> str:
        """Read a value in, for the conversions that still copy one.

        A character is not among them: `_remap_local` has already pointed the local at the
        caller's bytes, and there is nothing left to read.
        """
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
            case Conversion.MODE:
                writer.block(self._mode_in(argument, error))
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
            # A view, not a copy: nothing to allocate and nothing that can fail, so no ierr
            # check follows. It ends at the first NUL where the caller wrote one and spans
            # the buffer where it did not, so a C caller's "ward\0\0\0\0" and a generated
            # binding's "ward    " both reach the same case -- the second because Fortran
            # blank-pads the shorter operand of the comparison.
            writer.line(f"character(len=:), pointer :: {local}")
            writer.line(f"{local} => c_char_as_view({name})")
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
        """Copy back what was converted on the way in.

        Only a logical: a character local is a view *of* the caller's buffer, so the callee
        has already written there and copying it back would be copying it onto itself.
        """
        blocks = []
        for argument in wrapper:
            if argument.conversion is not Conversion.LOGICAL:
                continue
            if argument.intent is Intent.IN:
                continue
            blocks.append(self._output_conversion(argument))
        if blocks:
            writer.blank()
            for block in blocks:
                writer.block(block)

    def _output_conversion(self, argument: CArgument) -> str:
        name = argument.name
        # the same slice the input direction needs: an assumed-size array can be neither
        # assigned whole nor passed to an assumed-shape dummy
        target = self._conversion_source(argument)
        assignment = f"{target} = {name}{CONVERTED_SUFFIX}"
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
    """The tox_conversions helpers this wrapper calls.

    Only a mode string is left. A plain character is remapped, not converted, so a module
    full of strings imports nothing from `tox_conversions` in either direction.
    """
    helpers: set[str] = set()
    for argument in wrapper:
        if argument.conversion is Conversion.MODE:
            helpers.add("c_char_as_view")
    return helpers


# `is_err` used to be imported wherever a wrapper read a value in through a conversion that
# could fail -- which meant a mode string, the only one that allocated. Nothing in this layer
# converts an input any more: a character is remapped, a logical is assigned, and a mode takes
# a view, none of which can fail. An unused `only:` name is not a compile error in any Fortran,
# so the import would have gone stale silently; `TestImportsMatchTheBody` is what caught it.


def _chunks(items, size):
    items = list(items)
    for start in range(0, len(items), size):
        yield items[start : start + size]
