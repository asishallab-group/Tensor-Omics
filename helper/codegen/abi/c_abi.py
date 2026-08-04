"""Turning an IR procedure into the C ABI's view of it.

Every decision about what C sees is made here and nowhere else:

- what the symbol is called
- which arguments C is given that Fortran does not have, and what they are called
- what each Fortran type becomes on the C side
- what has to be converted, and what can be passed straight through

The emitters then render this. That is the whole point of the split: three targets that
must agree cannot each be left to work it out.
"""

from __future__ import annotations

from ..config import CONVENTIONS, Conventions
from ..diagnostics import DiagnosticBag
from ..ir.constants import ConstantError, ConstantEvaluator
from ..ir.doc import Doc
from ..ir.entities import Argument, Module, Procedure, Project
from ..ir.types import BaseType, CharacterLength, Dimension, FortranType, Intent
from .model import CArgument, CBinding, Conversion, CWrapper, CWrapperModule, Origin

#: Fortran kind -> the iso_c_binding kind C interoperates through.
#:
#: Explicit rather than inferred: a kind whose C counterpart is not written down here has
#: no defensible mapping, and guessing one produces a wrapper that compiles and lies.
KIND_MAP = {
    BaseType.INTEGER: {
        "int8": "c_int8_t",
        "int16": "c_int16_t",
        "int32": "c_int",
        "int64": "c_int64_t",
        "c_int": "c_int",
        "c_int8_t": "c_int8_t",
        "c_int16_t": "c_int16_t",
        "c_int32_t": "c_int32_t",
        "c_int64_t": "c_int64_t",
        "c_size_t": "c_size_t",
    },
    BaseType.REAL: {
        "real32": "c_float",
        "real64": "c_double",
        "c_float": "c_float",
        "c_double": "c_double",
    },
    BaseType.COMPLEX: {
        "real32": "c_float_complex",
        "real64": "c_double_complex",
        "c_float_complex": "c_float_complex",
        "c_double_complex": "c_double_complex",
    },
}

#: The C kind of a logical. Agreed ABI: C passes a real bool rather than an int, so the
#: binding languages hand over a boolean and not a 0/1 integer. The copy into a
#: default logical stays -- c_bool is one byte and a default logical is four -- but it is
#: one copy instead of a conversion at both ends. A logical already declared c_bool needs
#: no copy at all.
LOGICAL_C_KIND = "c_bool"

#: The C kind of a character. Its length becomes the leading extent.
CHARACTER_C_KIND = "c_char"


def build_project(project: Project, diagnostics: DiagnosticBag,
                  conventions: Conventions = CONVENTIONS,
                  evaluator: ConstantEvaluator | None = None) -> CBinding:
    """Build the C binding of every module that exports anything."""
    if evaluator is None:
        # The project's own parameters, so DM_DEFAULT(PI) resolves
        evaluator = ConstantEvaluator(project.constant_values())
    modules = []
    for module in project:
        if not module.has_exports:
            # A module with nothing to export produces no file at all
            continue
        modules.append(build_module(module, diagnostics, conventions, evaluator))
    return CBinding(tuple(modules))


def build_module(module: Module, diagnostics: DiagnosticBag,
                 conventions: Conventions = CONVENTIONS,
                 evaluator: ConstantEvaluator | None = None) -> CWrapperModule:
    wrappers = tuple(
        build_wrapper(procedure, diagnostics, conventions, evaluator)
        for procedure in module.exported_procedures
    )
    return CWrapperModule(
        name=f"{module.name}{conventions.c_suffix}", module=module, wrappers=wrappers
    )


def build_wrapper(procedure: Procedure, diagnostics: DiagnosticBag,
                  conventions: Conventions = CONVENTIONS,
                  evaluator: ConstantEvaluator | None = None) -> CWrapper:
    return _Builder(procedure, diagnostics, conventions, evaluator).build()


def c_symbol_name(procedure: Procedure, conventions: Conventions = CONVENTIONS) -> str:
    """The exported symbol for `procedure`.

    `<p>_alloc` is the one callers want, so it takes the plain name `<p>_c`. Its
    non-allocating twin `<p>` is the expert entry point and becomes `<p>_expert_c` -- but
    only where such a twin exists, otherwise a lone `<p>` would be needlessly renamed.
    """
    return f"{stripped_name(procedure, conventions)}{conventions.c_suffix}"


def stripped_name(procedure: Procedure, conventions: Conventions = CONVENTIONS) -> str:
    """The wrapper name without the `_c` suffix, as the binding languages call it."""
    name = procedure.name
    if procedure.is_alloc_variant:
        return name[: -len(conventions.alloc_suffix)]
    if procedure.has_alloc_sibling:
        return f"{name}{conventions.expert_infix}"
    return name


class _Builder:
    def __init__(self, procedure: Procedure, diagnostics: DiagnosticBag,
                 conventions: Conventions, evaluator: ConstantEvaluator | None = None):
        self.procedure = procedure
        self.diagnostics = diagnostics
        self.conventions = conventions
        self.evaluator = evaluator or ConstantEvaluator()
        self.taken = {a.name.lower() for a in procedure.arguments}
        if procedure.result is not None:
            self.taken.add(procedure.result.name.lower())

    def build(self) -> CWrapper:
        arguments: list[CArgument] = []
        for argument in self.procedure.arguments:
            arguments.extend(self._argument(argument))

        # A function is exposed as a subroutine: C receives the result as an argument.
        # It goes last so that adding a result to a subroutine does not renumber
        # everything a caller already passes.
        if self.procedure.result is not None:
            arguments.extend(self._argument(self.procedure.result, is_result=True))

        # Only when the procedure has none of its own: one it already declares was
        # emitted above, in the place the author put it.
        if not self.procedure.has_error_argument:
            arguments.append(self._error_argument())

        return CWrapper(
            name=c_symbol_name(self.procedure, self.conventions),
            stripped_name=stripped_name(self.procedure, self.conventions),
            procedure=self.procedure,
            arguments=tuple(arguments),
        )

    # -- arguments --------------------------------------------------------------

    def _argument(self, argument: Argument, is_result: bool = False) -> list[CArgument]:
        """One Fortran argument, plus whatever C needs alongside it."""
        if argument.type is None:
            return []

        # The error argument of the procedure itself is kept where it is; the wrapper
        # does not add a second one.
        extents: list[str] = list(argument.dimension.extents)
        extras: list[CArgument] = []

        for axis, extent in enumerate(extents):
            if extent != ":":
                continue
            # Assumed shape is not interoperable, so C is given the extent explicitly
            extra = self._extent_argument(argument, axis)
            extents[axis] = extra.name
            extras.append(extra)

        shape_arg = argument.roles.shape_arg if argument.roles else None
        if shape_arg is not None:
            # The shape travels separately, so the array itself is assumed size and the
            # synthesised extents are unnecessary
            extents = ["*"]
            extras = []

        conversion = Conversion.NONE
        mode = argument.roles.mode if argument.roles else None
        c_type = self._c_type(argument)
        if c_type is None:
            return []

        if mode is not None:
            # The integer the callee compares against parameters is a string to C
            length = str(max(mode.max_string_length, 1))
            c_type = FortranType(
                BaseType.CHARACTER,
                kind=CHARACTER_C_KIND,
                length=CharacterLength.parse("1"),
            )
            extents.insert(0, length)
            conversion = Conversion.MODE
        elif argument.type.is_character:
            length, strlen = self._character_length(argument)
            if strlen is not None:
                extras.insert(0, strlen)
            extents.insert(0, length)
            conversion = Conversion.CHARACTER
        elif argument.type.is_logical and not argument.type.is_c_bool:
            conversion = Conversion.LOGICAL

        size_extents = tuple(e for e in extents if e != "*")

        c_argument = CArgument(
            name=argument.name,
            type=c_type,
            dimension=Dimension(tuple(extents)),
            intent=argument.intent or Intent.INOUT,
            origin=self._origin_of(argument, is_result),
            conversion=conversion,
            optional=self._is_nullable(argument),
            doc=argument.doc,
            source=argument,
            mode=mode,
            shape_arg=shape_arg.name if shape_arg is not None else None,
            size_extents=size_extents,
            default=self._default_of(argument),
            has_default=argument.directives.has_default,
        )

        # The synthesised arguments follow the one they belong to, so the C signature
        # reads in the same order as the Fortran declaration.
        return [c_argument, *extras]

    def _extent_argument(self, argument: Argument, axis: int) -> CArgument:
        suffix = f"_dim_{axis + 1}" if argument.rank > 1 else ""
        name = self._unique(f"n_{argument.name}_elements{suffix}", argument)
        ordinal = f"{axis + 1}. dimension of" if argument.rank > 1 else "number of elements in"
        return CArgument(
            name=name,
            type=FortranType(BaseType.INTEGER, kind="c_int"),
            dimension=Dimension(),
            intent=Intent.IN,
            origin=Origin.EXTENT,
            doc=Doc.parse([f"{ordinal} `{argument.name}`"]),
            sizes=argument.name,
            axis=axis,
        )

    def _character_length(self, argument: Argument):
        """The leading extent carrying a character's length, and its argument if invented.

        A character's length is not part of its shape in Fortran, but in C it is just the
        fastest-moving extent of the buffer. `len=*` means the callee is told at run time,
        so C has to be told too.
        """
        length = argument.type.length
        if length.is_assumed:
            name = self._unique(f"{argument.name}_strlen", argument)
            extra = CArgument(
                name=name,
                type=FortranType(BaseType.INTEGER, kind="c_int"),
                dimension=Dimension(),
                intent=Intent.IN,
                origin=Origin.STRLEN,
                doc=Doc.parse([f"length of the strings in `{argument.name}`"]),
                sizes=argument.name,
            )
            return name, extra
        return str(length), None

    def _default_of(self, argument: Argument):
        """The value an omitted optional takes, evaluated now.

        The binding languages have to pass it, so it must be a value and not an
        expression by the time they are emitted. A default that will not evaluate is
        reported here, once, rather than by each emitter in turn.
        """
        directive = argument.directives.default
        if directive is None:
            return None
        try:
            value = self.evaluator.evaluate(directive.expression)
        except ConstantError as error:
            self.diagnostics.error(
                f"the default of '{argument.name}' is not a constant: {error}",
                entity=argument,
                note=(
                    "a default must be evaluable at generation time, because Python and "
                    "R pass it when the argument is omitted"
                ),
            )
            return None

        mode = argument.roles.mode if argument.roles else None
        if mode is not None:
            # A mode is a string to the bindings, so its default is the mode string whose
            # parameter has this integer value -- not the integer, which would be passed
            # through verbatim and rejected as an unknown mode.
            return self._mode_string_for(value, mode, argument)
        return value

    def _mode_string_for(self, value, mode, argument: Argument):
        """The mode string whose parameter evaluates to `value`, for a defaulted mode."""
        for entry in mode.values:
            try:
                if self.evaluator.evaluate(entry.parameter) == value:
                    return entry.string
            except ConstantError:
                continue
        self.diagnostics.error(
            f"the default of '{argument.name}' ({value!r}) matches no value in its mode table",
            entity=argument,
            note="a defaulted mode must default to one of its tabulated values",
        )
        return value

    @staticmethod
    def _is_nullable(argument: Argument) -> bool:
        """Whether C may pass a null pointer to mean "not given".

        An optional with a `DM_DEFAULT` is *not* nullable: the binding languages know
        the default and pass it, so C always receives a value. That is deliberate, from
        issue #131 -- every nullable optional is another branch in the wrapper, and the
        default has to be applied somewhere regardless. Applying it in Python and R keeps
        the wrapper flat and the default in one place.

        What stays nullable is an optional with no default: one required only in a
        certain mode, or an output the caller may decline.
        """
        return argument.optional and not argument.directives.has_default

    def _origin_of(self, argument: Argument, is_result: bool) -> Origin:
        if is_result:
            return Origin.RESULT
        if argument.name.lower() == self.conventions.error_arg:
            # The procedure's own error argument. Marked as such so nothing has to
            # recognise it by name again, and so a second one is not invented for it.
            return Origin.ERROR
        return Origin.ARGUMENT

    def _error_argument(self) -> CArgument:
        """The error code for a procedure that reports none of its own.

        The null checks report through it, so a wrapper without one could not tell the
        caller that a pointer was null. It goes last, where an added argument disturbs
        the fewest callers.
        """
        name = self._unique(self.conventions.error_arg, self.procedure)
        return CArgument(
            name=name,
            type=FortranType(BaseType.INTEGER, kind="c_int"),
            dimension=Dimension(),
            intent=Intent.OUT,
            origin=Origin.ERROR,
            doc=Doc.parse(["Error code"]),
        )

    # -- types ------------------------------------------------------------------

    def _c_type(self, argument: Argument) -> FortranType | None:
        fortran = argument.type

        if fortran.is_character:
            return FortranType(
                BaseType.CHARACTER,
                kind=CHARACTER_C_KIND,
                length=CharacterLength.parse("1"),
            )

        if fortran.is_logical:
            return FortranType(BaseType.LOGICAL, kind=LOGICAL_C_KIND)

        kinds = KIND_MAP.get(fortran.base)
        if kinds is None:
            self.diagnostics.error(
                f"'{argument.name}' is {fortran.base.value}, which has no C mapping",
                entity=argument,
            )
            return None

        c_kind = kinds.get((fortran.kind or "").lower())
        if c_kind is None:
            self.diagnostics.error(
                f"'{argument.name}' has kind '{fortran.kind}', which has no known C "
                f"counterpart",
                entity=argument,
                note=(
                    f"known {fortran.base.value} kinds: {', '.join(sorted(kinds))}"
                ),
            )
            return None

        return FortranType(fortran.base, kind=c_kind)

    # -- names ------------------------------------------------------------------

    def _unique(self, name: str, entity) -> str:
        """A synthesised name, reported rather than silently mangled if it collides.

        A collision means the wrapper would have two arguments of one name, which no
        amount of renaming makes less confusing than saying so.
        """
        if name.lower() in self.taken:
            self.diagnostics.error(
                f"the wrapper needs an argument called '{name}', but "
                f"'{self.procedure.name}' already has one",
                entity=entity,
                note=(
                    "this name is reserved for the argument the C binding adds; "
                    "rename the Fortran one"
                ),
            )
        self.taken.add(name.lower())
        return name
