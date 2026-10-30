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

import re

from ..config import CONVENTIONS, Conventions
from ..diagnostics import DiagnosticBag
from ..ir.constants import ConstantError, ConstantEvaluator, fits_integer_kind, integer_kind_value
from ..ir.doc import Doc
from ..ir.entities import Argument, Module, Procedure, Project
from ..ir.types import BaseType, CharacterLength, Dimension, FortranType, Intent
from ..ir.validate import bindings_cannot_carry
from ..names import extent_argument_name, strlen_argument_name
from .model import (
    CArgument,
    CBinding,
    Conversion,
    CWrapper,
    CWrapperModule,
    Origin,
    SizeConstant,
)

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

#: An identifier in an extent. A kind suffix (`0_int32` -> `_int32`) is the only spurious
#: match, and it names no module constant.
_IDENTIFIER_RE = re.compile(r"[A-Za-z_]\w*")

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


def reads_a_dummy(extent: str, procedure: Procedure) -> bool:
    """Whether an extent names a dummy of `procedure` anywhere in it, a bare name included."""
    dummies = {dummy.name.lower() for dummy in procedure.arguments}
    return any(identifier.lower() in dummies for identifier in _IDENTIFIER_RE.findall(extent))


def c_symbol_name(procedure: Procedure, conventions: Conventions = CONVENTIONS) -> str:
    """The exported symbol for `procedure`: what the bindings call it, plus `_c`."""
    return f"{stripped_name(procedure, conventions)}{conventions.c_suffix}"


def stripped_name(procedure: Procedure, conventions: Conventions = CONVENTIONS) -> str:
    """The wrapper name without the `_c` suffix, as the binding languages call it.

    Which is the Fortran name, unchanged. Synthesis names the tiers `foo` and `foo_expert`
    outright, and a hand-written export is published as written.

    It was not always the identity. Until f42 converted (2026-08-11) this translated a
    hand-written `foo` / `foo_alloc` pair into the published `foo` / `foo_expert`, which is
    how the framework spelled the two tiers before any of them were generated. No source
    writes that pair now, so the translation went with it: a rule that renames procedures
    behind the author's back has to earn its keep, and one with no remaining customer cannot.
    An author who writes `foo_alloc` today gets a procedure published as `foo_alloc`.

    `Conventions.alloc_suffix` outlives this. `validate` still refuses an implementation whose
    base name ends in it, which is what stops the shape being reintroduced where it would do
    real damage -- beside a generated `foo` that already means the allocating tier.
    """
    return procedure.name


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
        # A function is exposed as a subroutine: C receives the result as an argument, right
        # before the error code -- before the procedure's own `ierr` where it declares one,
        # else last, ahead of the `ierr` added below -- so that adding a result to a
        # subroutine does not move `ierr` for existing callers, and the signature ends in it
        # (#221 B1; it was appended after an own `ierr`, `hw_ownierr_c(n, x, ierr, total)`).
        arguments: list[CArgument] = []
        result = self.procedure.result
        for argument in self.procedure.arguments:
            if result is not None and argument.name.lower() == self.conventions.error_arg:
                arguments.extend(self._argument(result, is_result=True))
                result = None
            arguments.extend(self._argument(argument))
        if result is not None:
            arguments.extend(self._argument(result, is_result=True))

        # Only when the procedure has none of its own: one it already declares was
        # emitted above, in the place the author put it.
        if not self.procedure.has_error_argument:
            arguments.append(self._error_argument())
        self._check_error_is_last(arguments)

        return CWrapper(
            name=c_symbol_name(self.procedure, self.conventions),
            stripped_name=stripped_name(self.procedure, self.conventions),
            procedure=self.procedure,
            arguments=tuple(arguments),
            size_constants=self._size_constants(arguments),
        )

    def _check_error_is_last(self, arguments: list[CArgument]) -> None:
        """The C signature ends in the error code, as every published signature does.

        `validate` checks the Fortran order, and a procedure that passes it gets a C
        signature that ends in `ierr` too: a function's result is placed before it (`build`).
        Checked here all the same, because this is where C's order is decided, and an added
        argument is exactly what would move it.
        """
        error = next((a for a in arguments if a.origin is Origin.ERROR), None)
        if error is None or arguments[-1] is error:
            return
        own = self.procedure.argument(self.conventions.error_arg)
        if own is not None and self.procedure.arguments[-1] is not own:
            return  # `validate` reports the Fortran order, which this one only follows
        self.diagnostics.error(
            f"the C wrapper of '{self.procedure.name}' does not end in '{error.name}'",
            entity=self.procedure,
            note="every published signature ends in the error code",
        )

    def _size_constants(self, arguments: list[CArgument]) -> tuple[SizeConstant, ...]:
        """The module constants named in the extents C sees, where each is defined, and its value.

        An extent that is not an argument of this wrapper may be a module constant
        (`x(N_FIXED, n)`, a character's `len=MAX_LEN`). The C wrapper redeclares every
        argument and imports only the procedure it calls, so without these imports such a
        declaration names something undefined. Python and R cannot name it at all, so they
        are given its value, which therefore has to evaluate to an integer.
        """
        module = self.procedure.module
        project = module.project if module is not None else None
        if project is None:
            return ()
        dummies = {argument.name.lower() for argument in arguments}
        parameters = {parameter.name.lower(): parameter for parameter in project.parameters}
        found: dict[tuple[str, str], SizeConstant] = {}
        for argument in arguments:
            for extent in argument.dimension.extents:
                for identifier in _IDENTIFIER_RE.findall(extent):
                    parameter = parameters.get(identifier.lower())
                    if (identifier.lower() in dummies or parameter is None
                            or parameter.parent is None):
                        continue
                    key = (parameter.parent.name, parameter.name)
                    if key in found:
                        continue
                    value = self._constant_value(parameter)
                    if value is None:
                        self.diagnostics.error(
                            f"'{argument.name}' is sized by the module constant "
                            f"'{parameter.name}', which the generator cannot evaluate to an "
                            f"integer",
                            entity=argument.source or self.procedure,
                            note=(
                                "Python and R allocate and check the argument with the "
                                "constant's value, so it must be an integer constant "
                                "expression the generator can evaluate"
                            ),
                        )
                        continue
                    found[key] = SizeConstant(*key, value)
        return tuple(found[key] for key in sorted(found))

    def _constant_value(self, parameter) -> int | None:
        """A module parameter's integer value: from the project's table, else its own expression."""
        for expression in (parameter.name, parameter.expression):
            try:
                value = self.evaluator.evaluate(expression)
            except ConstantError:
                continue
            if isinstance(value, int) and not isinstance(value, bool):
                return value
        return None

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
        # what a reader is shown: an assumed-shape axis by the extent C is given, any other as
        # the author wrote it, which the axes passed below do not change
        declared = list(extents)

        shape_arg = argument.roles.shape_arg if argument.roles else None
        if shape_arg is not None:
            # The shape travels separately, so the array itself is assumed size and the
            # synthesised extents are unnecessary
            extents = ["*"]
            declared = ["*"]
            extras = []
        elif self._reads_its_own_extents(argument, is_result):
            for axis, extent in enumerate(extents):
                if not self._is_computed_extent(extent):
                    continue
                if not self._can_check(argument, extent):
                    continue
                # C is given the actual extent of an axis the declaration computes, declares
                # the array with it and compares the two (`CArgument.checked_against`)
                extra = self._extent_argument(argument, axis, checked_against=extent)
                extents[axis] = extra.name
                extras.append(extra)

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
            declared.insert(0, length)
            conversion = Conversion.MODE
        elif argument.type.is_character:
            length, strlen = self._character_length(argument)
            if strlen is not None:
                extras.insert(0, strlen)
            extents.insert(0, length)
            declared.insert(0, length)
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
            declared_extents=tuple(declared),
        )

        # The synthesised arguments follow the one they belong to, so the C signature
        # reads in the same order as the Fortran declaration.
        return [c_argument, *extras]

    def _extent_argument(self, argument: Argument, axis: int,
                         checked_against: str | None = None) -> CArgument:
        name = self._unique(extent_argument_name(argument.name, axis, argument.rank), argument)
        ordinal = f"{axis + 1}. dimension of" if argument.rank > 1 else "number of elements in"
        doc = f"{ordinal} `{argument.name}`"
        if checked_against is not None:
            doc += f", which must be `{checked_against}`"
        return CArgument(
            name=name,
            type=FortranType(BaseType.INTEGER, kind="c_int"),
            dimension=Dimension(),
            intent=Intent.IN,
            origin=Origin.EXTENT,
            doc=Doc.parse([doc]),
            sizes=argument.name,
            axis=axis,
            checked_against=checked_against,
        )

    # -- the extents C is given and the C wrapper checks (#221 A-shape) --------------------------

    def _reads_its_own_extents(self, argument: Argument, is_result: bool) -> bool:
        """Whether an argument's extents are read off what the caller passes, not computed.

        An input array, then: an output or a work array is allocated *from* its extents, by the
        caller or the binding, so there is nothing to compare them with. Only Python and R read
        the shape of what they are handed; C hands a pointer, and the C wrapper compares what it
        is told about an axis the declaration computes.
        """
        roles = argument.roles
        return (not is_result and argument.is_array and argument.intent is not None
                and argument.intent.is_input and not (roles and roles.is_temporary))

    def _dummies(self) -> dict[str, Argument]:
        return {dummy.name.lower(): dummy for dummy in self.procedure.arguments}

    def _is_computed_extent(self, extent: str) -> bool:
        """Whether an extent is computed from other dummies, rather than named or constant.

        A bare dummy name is read off the arrays by the binding and compared there, and an
        extent naming no dummy -- a literal, a module constant -- is a value the binding checks
        exactly. Everything else needs Fortran to evaluate it: `sum(reps_per_tissue)`, `n + 1`.
        """
        if extent.lower() in self._dummies():
            return False
        return reads_a_dummy(extent, self.procedure)

    def _can_check(self, argument: Argument, extent: str) -> bool:
        """Whether the C wrapper can evaluate `extent` once its pointers are checked; else refused.

        The check runs after every null check, so what it reads is safe to read -- unless it reads
        what has no defined value at all: an optional, which may be absent, or an `intent(out)`
        dummy, which is the procedure's to set. A dummy the C wrapper converts (a string, a mode,
        a default-kind logical) is not, under its name there, the value the declaration means.
        And an array whose own extent is computed would have to be checked first, an order this
        does not attempt: none of these occur on main, so they are refused instead.
        """
        dummies = self._dummies()
        for identifier in sorted({i.lower() for i in _IDENTIFIER_RE.findall(extent)}):
            dummy = dummies.get(identifier)
            if dummy is None or dummy.type is None:
                continue
            problem = None
            if dummy.optional:
                problem = f"it reads the optional '{dummy.name}', which may be absent"
            elif dummy.intent is Intent.OUT:
                problem = f"it reads '{dummy.name}', which is an output"
            elif (dummy.type.is_character or (dummy.roles and dummy.roles.is_mode)
                  or (dummy.type.is_logical and not dummy.type.is_c_bool)):
                problem = (f"it reads '{dummy.name}', which the C wrapper converts and so does "
                           f"not hold under that name")
            elif dummy.is_array and any(self._is_computed_extent(e) for e in dummy.dimension):
                problem = (f"it reads '{dummy.name}', whose own extent is computed too, so "
                           f"one check would have to wait for the other")
            if problem is not None:
                self.diagnostics.error(
                    f"'{argument.name}' of '{self.procedure.name}' is sized by '{extent}', "
                    f"which its C wrapper cannot check: {problem}",
                    entity=argument,
                    note=(
                        "C is given the actual extent of an input axis the declaration "
                        "computes, and the C wrapper compares it with the declaration before "
                        "the call; size the axis by an argument instead"
                    ),
                )
                return False
        return True

    def _character_length(self, argument: Argument):
        """The leading extent carrying a character's length, and its argument if invented.

        A character's length is not part of its shape in Fortran, but in C it is just the
        fastest-moving extent of the buffer. `len=*` means the callee is told at run time,
        so C has to be told too.
        """
        length = argument.type.length
        if length.is_assumed:
            name = self._unique(strlen_argument_name(argument.name), argument)
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

        kind = integer_kind_value(argument.type.kind)
        if (argument.type.base is BaseType.INTEGER and kind is not None
                and isinstance(value, int) and not isinstance(value, bool)
                and not fits_integer_kind(value, kind)):
            self.diagnostics.error(
                f"the default of '{argument.name}' ({value}) does not fit its kind, "
                f"{argument.type}",
                entity=argument,
                note=(
                    "Python and R would pass a value no Fortran caller can, wrapped to "
                    "the kind's range"
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
            if bindings_cannot_carry(fortran):
                return None  # `validate` has refused it already (`integer(8)`, `real(4)`)
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
