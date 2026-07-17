"""The semantic pass: what each argument *means*.

A signature says an argument is a scalar integer called `n_genes`. It does not say that
it is the first extent of `expr`, and so must not be asked of the caller in Python. That
meaning comes from the naming conventions and the `DM_` directives, and working it out is
this module's only job.

It is a procedure-level pass because most of the relations are between arguments: an
extent belongs to the arrays that declare it, a `<x>_shape` to the `<x>` it describes, an
`n_selected_<x>` to the mask it counts. `analyse` resolves them all at once and attaches
an `ArgumentRoles` to each argument.

Diagnostics are collected rather than raised, so one run reports every broken argument.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from ..config import CONVENTIONS, Conventions
from ..diagnostics import DiagnosticBag, SourceLocation
from .doc import DocTable, FordLink
from .entities import Argument, Procedure
from .types import BaseType, Intent


@dataclass(frozen=True)
class ModeValue:
    """One row of a mode table: a parameter, and the string C passes for it."""

    #: The Fortran parameter, e.g. `METHOD_WARD`
    parameter: str
    #: The module the parameter lives in
    module: str
    #: What C passes instead of the integer, e.g. `ward`
    string: str
    #: The prose from the first column
    description: str = ""


@dataclass(frozen=True)
class ModeTable:
    """The values a mode argument accepts, taken from its documentation table."""

    #: `mode` or `method`, as spelled in the table header
    alias: str
    values: tuple[ModeValue, ...]
    #: The table it came from, so an emitter can rewrite it in place
    source: DocTable | None = None

    @property
    def max_string_length(self) -> int:
        """The longest mode string, which sets the character length in the C wrapper."""
        return max((len(value.string) for value in self.values), default=0)

    def value_for(self, string: str) -> ModeValue | None:
        for value in self.values:
            if value.string == string:
                return value
        return None


@dataclass
class ArgumentRoles:
    """Everything the conventions and directives say about one argument."""

    #: `tmp_` prefixed: a work array, allocated silently and never returned
    is_temporary: bool = False

    #: The modes this argument accepts, when it is a mode argument
    mode: ModeTable | None = None

    #: Arrays that name this scalar as one of their extents
    extent_of: tuple[Argument, ...] = ()
    #: The array this `<x>_shape` argument describes
    shape_of: Argument | None = None
    #: The `<x>_shape` argument describing this array
    shape_arg: Argument | None = None
    #: The mask this `n_selected_<x>` counts
    mask_count_of: Argument | None = None
    #: The `n_selected_<x>` counting this mask
    count_arg: Argument | None = None
    #: The scalar naming how many leading elements of this array carry results
    result_size_arg: Argument | None = None

    @property
    def is_mode(self) -> bool:
        return self.mode is not None

    @property
    def is_extent(self) -> bool:
        return bool(self.extent_of)

    @property
    def is_shape_arg(self) -> bool:
        return self.shape_of is not None

    @property
    def has_shape_arg(self) -> bool:
        return self.shape_arg is not None

    @property
    def is_mask_count(self) -> bool:
        return self.mask_count_of is not None

    @property
    def is_derived(self) -> bool:
        """Whether the interfacing languages can work this argument out themselves.

        Such an argument is not asked of the caller: it comes from another argument.
        """
        return self.is_extent or self.is_shape_arg or self.is_mask_count


def analyse(procedure: Procedure, diagnostics: DiagnosticBag,
            conventions: Conventions = CONVENTIONS) -> None:
    """Work out and attach the roles of every argument of `procedure`."""
    _Analyser(procedure, diagnostics, conventions).run()


def analyse_project(project, diagnostics: DiagnosticBag,
                    conventions: Conventions = CONVENTIONS) -> None:
    """Analyse every procedure marked for export.

    Only those. Analysis reports what it cannot interpret -- a mode argument with no
    table of values, say -- and an internal routine is held to no such contract, because
    nothing is generated from it.
    """
    for module in project:
        for procedure in module.exported_procedures:
            analyse(procedure, diagnostics, conventions)


class _Analyser:
    def __init__(self, procedure: Procedure, diagnostics: DiagnosticBag,
                 conventions: Conventions):
        self.procedure = procedure
        self.diagnostics = diagnostics
        self.conventions = conventions
        self.arguments = tuple(procedure.arguments)
        if procedure.result is not None:
            self.arguments += (procedure.result,)

    def run(self) -> None:
        for argument in self.arguments:
            argument.roles = ArgumentRoles(
                is_temporary=self._is_temporary(argument),
                mode=self._mode_table(argument),
            )

        # Relations need every argument to exist first
        self._link_extents()
        self._link_shape_arguments()
        self._link_mask_counts()
        self._link_result_sizes()

    def _is_temporary(self, argument: Argument) -> bool:
        return argument.name.lower().startswith(self.conventions.temporary_prefix)

    # -- modes ------------------------------------------------------------------

    def _mode_table(self, argument: Argument) -> ModeTable | None:
        alias = self.conventions.mode_alias_of(argument.name)
        if alias is None:
            return None

        tables = [t for t in argument.doc.tables if self._looks_like_a_mode_table(t)]
        if not tables:
            self.diagnostics.error(
                f"mode argument '{argument.name}' has no table of accepted values",
                entity=argument,
                note=_mode_table_example(alias),
            )
            return None
        if len(tables) > 1:
            self.diagnostics.error(
                f"mode argument '{argument.name}' has {len(tables)} tables of accepted values",
                entity=argument,
                note="a mode argument documents its values in exactly one table",
            )
            return None

        table = tables[0]
        header_alias = table.header_text[0].strip().lower()
        if header_alias != alias:
            # The argument's own name is authoritative. Letting the header decide means
            # a 'link_method' argument documented with a 'Mode' header and MODE_ prefixed
            # parameters is accepted, and nothing ever points out the inconsistency.
            self.diagnostics.error(
                f"the table of '{argument.name}' is headed '{table.header_text[0].strip()}' "
                f"but '{argument.name}' is a {alias} argument",
                entity=argument,
                location=_location_of(argument, table.line_number),
                note=(
                    f"head the table '{alias.capitalize()}' and name its parameters "
                    f"{alias.upper()}_<NAME>, or rename the argument"
                ),
            )
            return None

        values = self._mode_values(argument, table, alias)
        if values is None:
            return None
        if not values:
            self.diagnostics.error(
                f"mode argument '{argument.name}' has a table with no values",
                entity=argument,
                note=_mode_table_example(alias),
            )
            return None

        return ModeTable(alias=alias, values=tuple(values), source=table)

    def _looks_like_a_mode_table(self, table: DocTable) -> bool:
        """A two-column table headed by a mode alias and `Value`.

        Recognition is deliberately shallow. Whether the *rows* are well formed is
        checked afterwards and reported, rather than quietly disqualifying the table --
        the old parser bailed out of mode detection on a malformed row, and the author
        then got 'no mode table' pointing away from the actual mistake.
        """
        if table.n_columns != 2:
            return False
        header = [text.strip().lower() for text in table.header_text]
        return header[0] in self.conventions.mode_aliases and header[1] == "value"

    def _mode_values(self, argument, table: DocTable, header_alias: str):
        expected_prefix = f"{header_alias.upper()}_"
        pattern = re.compile(rf"{expected_prefix}(?P<name>[A-Z_0-9]+)\Z")
        values: list[ModeValue] = []

        for row_index, (description, value) in enumerate(table.rows, start=1):
            # An empty Value cell documents a mode that has no parameter yet
            if not value.text.strip():
                continue

            links = value.links
            if len(links) != 1 or value.text != str(links[0]):
                self.diagnostics.error(
                    f"row {row_index} of the {header_alias} table of '{argument.name}' "
                    f"does not name a parameter",
                    entity=argument,
                    location=_location_of(argument, value.line_number),
                    note=_mode_table_example(header_alias),
                )
                return None

            link = links[0]
            if link.item is None:
                self.diagnostics.error(
                    f"row {row_index} of the {header_alias} table of '{argument.name}' "
                    f"links to module '{link.component}' rather than to a parameter in it",
                    entity=argument,
                    location=_location_of(argument, value.line_number),
                    note=_mode_table_example(header_alias),
                )
                return None

            match = pattern.match(link.item)
            if match is None:
                self.diagnostics.error(
                    f"'{link.item}' is not a valid {header_alias} parameter",
                    entity=argument,
                    location=_location_of(argument, value.line_number),
                    note=(
                        f"a {header_alias} parameter is named {expected_prefix}<NAME> in "
                        f"upper case, and C is passed '<name>' in lower case"
                    ),
                )
                return None

            values.append(
                ModeValue(
                    parameter=link.item,
                    module=link.component,
                    string=match.group("name").lower(),
                    description=description.text.strip(),
                )
            )

        duplicates = _duplicates(value.string for value in values)
        if duplicates:
            self.diagnostics.error(
                f"the {header_alias} table of '{argument.name}' maps "
                f"{', '.join(sorted(duplicates))} more than once",
                entity=argument,
                note="each mode string must identify exactly one parameter",
            )
            return None

        return values

    # -- relations between arguments --------------------------------------------

    def _link_extents(self) -> None:
        for argument in self.arguments:
            if not argument.is_scalar or argument.type.base is not BaseType.INTEGER:
                continue
            name = argument.name.lower()
            owners = tuple(
                other
                for other in self.arguments
                if other is not argument
                and any(extent.lower() == name for extent in other.dimension)
            )
            if owners:
                argument.roles.extent_of = owners

    def _link_shape_arguments(self) -> None:
        suffix = self.conventions.shape_suffix
        for argument in self.arguments:
            if not argument.name.lower().endswith(suffix):
                continue
            owner = self.procedure.argument(argument.name[: -len(suffix)])
            if owner is None:
                continue
            argument.roles.shape_of = owner
            owner.roles.shape_arg = argument

    def _link_mask_counts(self) -> None:
        prefix = self.conventions.mask_count_prefix
        for argument in self.arguments:
            name = argument.name.lower()
            if not name.startswith(prefix):
                continue
            owner_name = name[len(prefix):]
            if not owner_name:
                continue
            mask = self._find_mask_for(owner_name)
            if mask is None:
                continue
            argument.roles.mask_count_of = mask
            mask.roles.count_arg = argument

    def _find_mask_for(self, owner_name: str) -> Argument | None:
        for candidate in self.arguments:
            if self.conventions.mask_arg_name_of(candidate.name) == owner_name:
                return candidate
        return None

    def _link_result_sizes(self) -> None:
        for argument in self.arguments:
            directive = argument.directives.result_size_is
            if directive is None:
                continue
            counter = self.procedure.argument(directive.argument)
            if counter is None:
                self.diagnostics.error(
                    f"'{argument.name}' says its result size is '{directive.argument}', "
                    f"which is not an argument of '{self.procedure.name}'",
                    entity=argument,
                    location=_location_of(argument, directive.line_number),
                )
                continue
            if counter.is_array or counter.type.base is not BaseType.INTEGER:
                self.diagnostics.error(
                    f"'{directive.argument}' cannot be the result size of "
                    f"'{argument.name}': a result size is a scalar integer",
                    entity=argument,
                    location=_location_of(argument, directive.line_number),
                )
                continue
            argument.roles.result_size_arg = counter


def _duplicates(values) -> set:
    seen, repeated = set(), set()
    for value in values:
        if value in seen:
            repeated.add(value)
        seen.add(value)
    return repeated


def _location_of(argument: Argument, line_number: int | None) -> SourceLocation:
    if line_number is None:
        return argument.location
    return SourceLocation(argument.location.file, line_number)


def _mode_table_example(alias: str) -> str:
    title = alias.capitalize()
    return (
        "document the accepted values in a table, like:\n"
        f"!! | {title} | Value |\n"
        "!! |------|-------|\n"
        f"!! | what it does | [[a_module(module):{alias.upper()}_BLA(variable)]] |"
    )
