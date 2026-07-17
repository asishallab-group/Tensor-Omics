"""Rules an exported procedure must satisfy.

Only exported procedures are checked. An internal routine is held to none of this,
because nothing is generated from it -- it may use deferred-length strings, derived types
and untyped intents as freely as Fortran allows.

Each rule exists because generating correct code is otherwise impossible, and each says
what to do instead. A rule that merely encodes taste belongs in a review, not here.
"""

from __future__ import annotations

from ..config import CONVENTIONS, Conventions
from ..diagnostics import DiagnosticBag
from .entities import Argument, Module, Procedure, Project
from .types import BaseType, Intent

#: Attributes on a dummy argument that have no C interoperable form
NON_INTEROPERABLE_ATTRIBUTES = ("allocatable", "pointer")

#: The greatest character rank tox_conversions can convert. A scalar character maps to a
#: 1-D c_char buffer and a character vector to a 2-D one, which is where
#: c_char_1d_as_string / c_char_2d_as_string stop.
MAX_CHARACTER_RANK = 1


def validate_project(project: Project, diagnostics: DiagnosticBag,
                     conventions: Conventions = CONVENTIONS) -> None:
    for module in project:
        validate_module(module, diagnostics, conventions)


def validate_module(module: Module, diagnostics: DiagnosticBag,
                    conventions: Conventions = CONVENTIONS) -> None:
    if not module.doc:
        diagnostics.warn("module has no documentation", entity=module)
    for procedure in module.exported_procedures:
        validate_procedure(procedure, diagnostics, conventions)


def validate_procedure(procedure: Procedure, diagnostics: DiagnosticBag,
                       conventions: Conventions = CONVENTIONS) -> None:
    _Validator(procedure, diagnostics, conventions).run()


class _Validator:
    def __init__(self, procedure: Procedure, diagnostics: DiagnosticBag,
                 conventions: Conventions):
        self.procedure = procedure
        self.diagnostics = diagnostics
        self.conventions = conventions
        self.arguments = tuple(procedure.arguments)
        if procedure.result is not None:
            self.arguments += (procedure.result,)

    def error(self, message, argument, note=None):
        self.diagnostics.error(message, entity=argument, note=note)

    def run(self) -> None:
        self._check_meta()
        for argument in self.arguments:
            self._check_argument(argument)

    def _check_meta(self) -> None:
        if not self.procedure.meta.summary:
            self.diagnostics.warn(
                "procedure has no summary meta tag",
                entity=self.procedure,
                note="add '!> summary: ...' -- it becomes the docstring in Python and R",
            )
        if not self.procedure.meta.author:
            self.diagnostics.warn("procedure has no author meta tag", entity=self.procedure)

    def _check_argument(self, argument: Argument) -> None:
        self._check_documented(argument)
        self._check_intent(argument)
        self._check_attributes(argument)
        self._check_type(argument)
        self._check_temporary(argument)
        self._check_mode(argument)
        self._check_shape_argument(argument)
        self._check_extent(argument)
        self._check_error_argument(argument)

    def _check_documented(self, argument: Argument) -> None:
        if not argument.doc and not argument.is_result:
            self.diagnostics.warn(
                f"argument '{argument.name}' has no documentation",
                entity=argument,
                note="it is inherited by the C wrapper and by the Python and R docstrings",
            )

    def _check_intent(self, argument: Argument) -> None:
        if argument.is_result:
            return
        if argument.intent is None:
            self.error(
                f"argument '{argument.name}' has no intent",
                argument,
                note=(
                    "an exported argument needs an explicit intent: it decides constness "
                    "in C and R, and whether the argument is an input, an output or both"
                ),
            )

    def _check_attributes(self, argument: Argument) -> None:
        for attribute in NON_INTEROPERABLE_ATTRIBUTES:
            if argument.has_attribute(attribute):
                self.error(
                    f"argument '{argument.name}' is '{attribute}', which has no C "
                    f"interoperable form",
                    argument,
                    note=(
                        "pass the data and its extents as separate arguments, and let the "
                        "_alloc variant own the allocation"
                    ),
                )

    def _check_type(self, argument: Argument) -> None:
        if argument.type.base is BaseType.DERIVED:
            self.error(
                f"argument '{argument.name}' is a derived type, which the generator "
                f"cannot map to C yet",
                argument,
            )
            return

        if argument.type.is_character:
            self._check_character(argument)

    def _check_character(self, argument: Argument) -> None:
        length = argument.type.length
        if length.is_deferred:
            self.error(
                f"argument '{argument.name}' has a deferred length (len=:)",
                argument,
                note=(
                    "a deferred length implies allocatable or pointer, neither of which "
                    "is interoperable; use len=* to have the length passed in, or a fixed "
                    "len=n"
                ),
            )

        if argument.rank > MAX_CHARACTER_RANK:
            self.error(
                f"argument '{argument.name}' is a rank-{argument.rank} character array, "
                f"and tox_conversions can convert up to rank {MAX_CHARACTER_RANK}",
                argument,
                note=(
                    "a character carries its length as a leading C extent, so a rank-2 "
                    "character array would need a rank-3 c_char conversion"
                ),
            )

    def _check_temporary(self, argument: Argument) -> None:
        if argument.roles is None or not argument.roles.is_temporary:
            return
        if argument.intent is Intent.IN:
            self.error(
                f"temporary '{argument.name}' is intent(in)",
                argument,
                note=(
                    "a work array is allocated by the caller in Python and R, so it is "
                    "intent(out) or intent(inout); drop the "
                    f"'{self.conventions.temporary_prefix}' prefix if it is really an input"
                ),
            )

    def _check_mode(self, argument: Argument) -> None:
        if self.conventions.mode_alias_of(argument.name) is None:
            return
        # A missing or malformed table is already reported by roles.analyse
        if argument.type.base is not BaseType.INTEGER or argument.is_array:
            self.error(
                f"mode argument '{argument.name}' is not a scalar integer",
                argument,
                note=(
                    "a mode is a scalar integer compared against the MODE_/METHOD_ "
                    "parameters; the C wrapper is what receives it as a string"
                ),
            )

    def _check_shape_argument(self, argument: Argument) -> None:
        roles = argument.roles
        if roles is None:
            return

        if roles.is_shape_arg:
            self._check_is_a_usable_shape_argument(argument)
            owner = roles.shape_of
            if owner.rank != 1:
                self.error(
                    f"'{owner.name}' has a shape argument but is rank {owner.rank}",
                    owner,
                    note=(
                        "an argument whose shape is passed separately is flat: declare it "
                        "rank-1 and let "
                        f"'{argument.name}' carry the shape"
                    ),
                )

    def _check_is_a_usable_shape_argument(self, argument: Argument) -> None:
        if argument.intent is not Intent.IN:
            self.error(
                f"shape argument '{argument.name}' is not intent(in)",
                argument,
                note="a shape describes an argument, so it is only ever read",
            )
        if argument.type.base is not BaseType.INTEGER or argument.rank != 1:
            self.error(
                f"shape argument '{argument.name}' is not a rank-1 integer array",
                argument,
                note="a shape is the list of extents, one integer per dimension",
            )
        if argument.optional:
            self.error(
                f"shape argument '{argument.name}' is optional",
                argument,
                note=_NON_OPTIONAL_NOTE,
            )

    def _check_extent(self, argument: Argument) -> None:
        roles = argument.roles
        if roles is None or not roles.is_extent:
            return
        if argument.optional:
            owners = ", ".join(f"'{owner.name}'" for owner in roles.extent_of)
            self.error(
                f"extent argument '{argument.name}' is optional, but {owners} needs it",
                argument,
                note=_NON_OPTIONAL_NOTE,
            )

    def _check_error_argument(self, argument: Argument) -> None:
        if argument.name.lower() != self.conventions.error_arg:
            return
        if argument.type.base is not BaseType.INTEGER or argument.is_array:
            self.error(
                f"'{argument.name}' is not a scalar integer", argument,
                note="the error argument carries an encoded tox_errors code",
            )
        elif argument.intent is Intent.IN:
            self.error(
                f"'{argument.name}' is intent(in), so no error can be reported through it",
                argument,
            )


_NON_OPTIONAL_NOTE = (
    "the C wrapper reads it before it may take c_loc of the arrays it describes, so it "
    "has to be there; see the null validation order in the generator README"
)
