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
from .constants import ConstantError, ConstantEvaluator
from .doc import DocTable, FordLink
from .entities import Argument, Procedure
from .types import BaseType, Intent

#: A generated procedure name in a mode table's third column: a Fortran identifier.
_PROCEDURE_NAME_RE = re.compile(r"[A-Za-z]\w*\Z")


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
    #: The per-mode procedure name from the optional third column; set only when the mode
    #: table opts into per-mode splitting
    procedure_name: str = ""


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

    @property
    def is_split(self) -> bool:
        """Whether the table names a procedure per mode, opting into per-mode wrappers.

        The third column is all-or-nothing (the reader rejects a partial one), so any value
        carrying a procedure name means every one does.
        """
        return any(value.procedure_name for value in self.values)

    def value_for(self, string: str) -> ModeValue | None:
        for value in self.values:
            if value.string == string:
                return value
        return None


@dataclass(frozen=True)
class ProducerInput:
    """One input of a `DM_OUTPUT_FROM` producer, and where its value comes from."""

    #: The producer's own parameter name, which is what the call must use as a keyword
    name: str
    #: The consumer argument supplying it, or None when it is a constant
    argument: str | None = None
    #: The value, when the producer input is a constant the consumer has no argument for
    constant: object = None

    @property
    def is_constant(self) -> bool:
        return self.argument is None


@dataclass(frozen=True)
class OutputFromPlan:
    """How an argument is obtained by calling another procedure.

    For `DM_OUTPUT_FROM(count, mask_chunk_count, ..., AUTO)`: call `mask_chunk_count`,
    supply its inputs from the consumer's own arguments, and take its `count` output. The
    binding languages call the producer's *generated wrapper*, so its own error
    checking and result handling come for free.
    """

    #: The procedure to call
    producer: Procedure
    #: The producer's output argument whose value fills the consumer argument
    output: Argument
    #: Where each of the producer's inputs comes from, in the producer's own order
    inputs: tuple[ProducerInput, ...]
    #: AUTO -> the languages call it; JUST_INFO -> the doc only tells the caller to
    is_automatic: bool

    @property
    def inout_feedback(self) -> tuple[ProducerInput, ...]:
        """The producer inputs it refines in place, that a consumer argument supplied.

        An `intent(inout)` producer input is returned alongside the annotated output with a
        possibly different value -- `calc_work_arr_paralog_subsets_size` caps the
        `max_subset_size` it is handed to what actually fits. After the AUTO call the
        consumer must adopt that returned value, not the raw one it passed in, or it works
        with a stale figure. Constants are never fed back (nothing supplied them).
        """
        fed_back = []
        for supply in self.inputs:
            if supply.argument is None:
                continue
            producer_input = self.producer.argument(supply.name)
            intent = producer_input.intent if producer_input else None
            if intent and intent.is_input and intent.is_output:
                fed_back.append(supply)
        return tuple(fed_back)


@dataclass
class ArgumentRoles:
    """Everything the conventions and directives say about one argument."""

    #: `tmp_` prefixed: a work array, allocated silently and never returned
    is_temporary: bool = False

    #: How this argument is computed from another procedure, if it is
    computed_from: OutputFromPlan | None = None

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
    def is_inferable_extent(self) -> bool:
        """Whether an extent can actually be read off an argument the caller supplies.

        An extent named only by `intent(out)` arrays -- the length of a result whose size
        the caller alone knows -- cannot be. Neither can one named only by work arrays:
        those are allocated *from* the extent, so reading it back off them is circular.
        Treating either as derived would drop it from the signature and then leave nothing
        to compute it from, so it stays a parameter.
        """
        def readable(owner) -> bool:
            if owner.roles and owner.roles.is_temporary:
                return False
            if owner.roles and owner.roles.shape_arg is not None:
                # its extents travel in a shape argument, so the count is their product
                # whichever way the array itself flows
                return True
            return owner.intent.is_input

        return any(readable(owner) for owner in self.extent_of)

    @property
    def is_shape_arg(self) -> bool:
        return self.shape_of is not None

    @property
    def is_inferable_shape_arg(self) -> bool:
        """Whether the shape can be read off the array it describes.

        Only when that array is an input. Describing an `intent(out)` array, the shape is
        what the caller states the result should be, so it has to be asked for.
        """
        return self.shape_of is not None and self.shape_of.intent.is_input

    @property
    def has_shape_arg(self) -> bool:
        return self.shape_arg is not None

    @property
    def is_mask_count(self) -> bool:
        return self.mask_count_of is not None

    @property
    def is_computed(self) -> bool:
        """Computed by calling another procedure (AUTO). JUST_INFO does not count -- the
        caller still supplies it, the documentation only says where to get it."""
        return self.computed_from is not None and self.computed_from.is_automatic

    @property
    def is_derived(self) -> bool:
        """Whether the binding languages can work this argument out themselves.

        Such an argument is not asked of the caller: it comes from another argument, or
        from calling another procedure.
        """
        return (self.is_inferable_extent or self.is_inferable_shape_arg
                or self.is_mask_count or self.is_computed)


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

    # output_from is a cross-procedure relation, so it is resolved once every procedure
    # has its own roles and the producer can be looked up. The evaluator is built once:
    # a producer-input table may supply a constant, which has to be a value by the time
    # an emitter passes it on.
    evaluator = ConstantEvaluator(project.constant_values())
    for module in project:
        for procedure in module.exported_procedures:
            _resolve_output_from(procedure, project, diagnostics, conventions, evaluator)


def _resolve_output_from(consumer: Procedure, project, diagnostics: DiagnosticBag,
                         conventions: Conventions = CONVENTIONS,
                         evaluator: ConstantEvaluator | None = None) -> None:
    """Attach an `OutputFromPlan` to every argument documented with `DM_OUTPUT_FROM`.

    Each of the producer's inputs is supplied by the consumer argument of the same name.
    Where the two disagree, the argument documents the mapping in a table:

        !! | Producer input | Supplied by |
        !! |----------------|-------------|
        !! | n_elements     | n_genes     |

    A cell that is not an argument of the consumer is evaluated as a Fortran constant, so
    a producer input the consumer simply has no argument for -- `n_dim` on a routine that
    is always one-dimensional -- can be given a value outright.

    A producer input that is neither name-matched, nor mapped to an argument, nor a
    constant is an error naming it. The producer may live in another module; the Python
    emitter imports it where it is called, and R has every wrapper in one environment.
    """
    for argument in consumer.arguments:
        directive = argument.directives.output_from
        if directive is None:
            continue

        producer = project.procedure(directive.module, directive.procedure)
        if producer is None:
            diagnostics.error(
                f"'{argument.name}' is computed from '{directive.procedure}' in "
                f"'{directive.module}', which does not exist",
                entity=argument,
            )
            continue

        if not producer.is_exported:
            diagnostics.error(
                f"'{argument.name}' is computed from '{producer.name}', which is not "
                f"exported, so there is no wrapper to call",
                entity=argument,
                note=(
                    f"add 'category: {CONVENTIONS.c_binding_category}' to "
                    f"'{producer.name}'"
                ),
            )
            continue

        output = producer.argument(directive.argument)
        if output is None or not (output.intent and output.intent.is_output):
            diagnostics.error(
                f"'{directive.argument}' is not an output of '{producer.name}'",
                entity=argument,
            )
            continue

        if not directive.is_automatic:
            # JUST_INFO is documentation only: the caller still supplies the argument, and
            # nothing calls the producer. So the consumer need not be able to feed the
            # producer's inputs -- a pointer to a producer whose inputs live elsewhere (a
            # filter that needs data the consumer never sees) is exactly what JUST_INFO is
            # for. Record the relation, resolving no inputs; the producer existing, being
            # exported, and owning the named output (checked above) is all that must hold.
            argument.roles.computed_from = OutputFromPlan(
                producer=producer, output=output, inputs=(), is_automatic=False,
            )
            continue

        renames = _producer_input_table(argument, consumer, producer,
                                        diagnostics, conventions)
        if renames is None:
            continue

        inputs = []
        unmatched = []
        for producer_input in producer.arguments:
            if producer_input is output or not producer_input.intent.is_input:
                continue
            if producer_input.name.lower() == "ierr":
                continue
            if producer_input.roles and producer_input.roles.is_derived:
                # the producer's own wrapper derives this input (an extent read off one of
                # its array arguments, a shape, a mask count, or a nested producer) and
                # drops it from its signature, so the call must not -- and need not --
                # supply it: the producer recomputes it from the array inputs it is passed
                continue
            supplier = renames.get(producer_input.name.lower(), producer_input.name)
            match = consumer.argument(supplier)
            if match is not None:
                inputs.append(ProducerInput(producer_input.name, argument=match.name))
                continue
            # not an argument of the consumer: a table may supply it as a constant, which
            # is how a producer input the consumer simply does not have gets a value
            if supplier != producer_input.name and evaluator is not None:
                try:
                    value = evaluator.evaluate(supplier)
                except ConstantError:
                    pass
                else:
                    inputs.append(ProducerInput(producer_input.name, constant=value))
                    continue
            unmatched.append(producer_input.name)

        if unmatched:
            diagnostics.error(
                f"'{argument.name}' is computed from '{producer.name}', but its "
                f"input(s) {', '.join(sorted(unmatched))} have no argument of the same "
                f"name in '{consumer.name}'",
                entity=argument,
                note=_producer_input_table_example(conventions),
            )
            continue

        argument.roles.computed_from = OutputFromPlan(
            producer=producer,
            output=output,
            inputs=tuple(inputs),
            is_automatic=directive.is_automatic,
        )



def _looks_like_a_producer_input_table(table: DocTable, conventions: Conventions) -> bool:
    """A two-column table headed `Producer input` and `Supplied by`.

    Shallow, like the mode-table check: a table that looks like this one is treated as
    one, and anything wrong with its rows is reported rather than quietly disqualifying
    it and leaving the author with a misleading 'no matching argument' further down.
    """
    if table.n_columns != 2:
        return False
    header = [text.strip().lower() for text in table.header_text]
    return (header[0] == conventions.producer_input_header
            and header[1] == conventions.producer_supplied_by_header)


def _producer_input_table(argument: Argument, consumer: Procedure, producer: Procedure,
                          diagnostics: DiagnosticBag,
                          conventions: Conventions) -> dict[str, str] | None:
    """The producer-input -> consumer-argument overrides this argument documents.

    An empty mapping when there is no table. `None` when the table is there but wrong, so
    the caller stops rather than reporting a pile of consequential name mismatches.
    """
    tables = [t for t in argument.doc.tables
              if _looks_like_a_producer_input_table(t, conventions)]
    if not tables:
        return {}
    if len(tables) > 1:
        diagnostics.error(
            f"'{argument.name}' has {len(tables)} tables mapping the inputs of "
            f"'{producer.name}'",
            entity=argument,
            note="the mapping goes in exactly one table",
        )
        return None

    mapping: dict[str, str] = {}
    for row_index, (producer_input, supplied_by) in enumerate(tables[0].rows, start=1):
        wanted = producer_input.text.strip()
        given = supplied_by.text.strip()
        if not wanted or not given:
            diagnostics.error(
                f"row {row_index} of the producer-input table of '{argument.name}' "
                f"has an empty cell",
                entity=argument,
                note=_producer_input_table_example(conventions),
            )
            return None
        if producer.argument(wanted) is None:
            diagnostics.error(
                f"the producer-input table of '{argument.name}' maps '{wanted}', which "
                f"is not an argument of '{producer.name}'",
                entity=argument,
            )
            return None
        mapping[wanted.lower()] = given
    return mapping


def _producer_input_table_example(conventions: Conventions) -> str:
    left = conventions.producer_input_header.title()
    right = conventions.producer_supplied_by_header.capitalize()
    return (
        "name-matching is tried first; where the names differ, or where the consumer has\n"
        "no argument for an input at all, map it in a table -- the right-hand cell is an\n"
        "argument of the consumer, or a Fortran constant:\n"
        f"!! | {left} | {right} |\n"
        f"!! |----------------|-------------|\n"
        f"!! | n_elements     | n_genes     |\n"
        f"!! | n_dim          | 1_int32     |"
    )

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
        """A mode table: headed by a mode alias and `Value`, and optionally `Procedure`.

        Recognition is deliberately shallow. Whether the *rows* are well formed is
        checked afterwards and reported, rather than quietly disqualifying the table --
        the old parser bailed out of mode detection on a malformed row, and the author
        then got 'no mode table' pointing away from the actual mistake. The optional third
        column opts the argument into per-mode splitting.
        """
        if table.n_columns not in (2, 3):
            return False
        header = [text.strip().lower() for text in table.header_text]
        if header[0] not in self.conventions.mode_aliases or header[1] != "value":
            return False
        if table.n_columns == 3 and header[2] != self.conventions.mode_procedure_header:
            return False
        return True

    def _mode_values(self, argument, table: DocTable, header_alias: str):
        expected_prefix = f"{header_alias.upper()}_"
        pattern = re.compile(rf"{expected_prefix}(?P<name>[A-Z_0-9]+)\Z")
        has_procedure = table.n_columns == 3
        values: list[ModeValue] = []

        for row_index, row in enumerate(table.rows, start=1):
            description, value = row[0], row[1]
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

            procedure_name = ""
            if has_procedure:
                procedure_name = row[2].text.strip()
                if not _PROCEDURE_NAME_RE.match(procedure_name):
                    self.diagnostics.error(
                        f"row {row_index} of the {header_alias} table of '{argument.name}' "
                        f"has no valid procedure name in its third column",
                        entity=argument,
                        location=_location_of(argument, value.line_number),
                        note=(
                            "the third column names the generated procedure for that mode, "
                            "an identifier such as `detect_dosage_effect`"
                        ),
                    )
                    return None

            values.append(
                ModeValue(
                    parameter=link.item,
                    module=link.component,
                    string=match.group("name").lower(),
                    description=description.text.strip(),
                    procedure_name=procedure_name,
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

        if has_procedure:
            repeated = _duplicates(v.procedure_name for v in values)
            if repeated:
                self.diagnostics.error(
                    f"the {header_alias} table of '{argument.name}' names procedure(s) "
                    f"{', '.join(sorted(repeated))} more than once",
                    entity=argument,
                    note="each mode's procedure name must be unique",
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
