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
from .constants import ConstantError, ConstantEvaluator, fits_integer_kind
from .doc import DocTable, FordLink
from .entities import Argument, Procedure
from .types import BaseType, FortranType, Intent
from .validate import _kind_value, bindings_cannot_carry

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
    #: Whether the producer's own wrapper derives it (an extent, a shape, a mask count, a
    #: nested producer's value) and so leaves it out of what Python and R ask for. Only a
    #: generated Fortran wrapper, which calls the producer directly, passes it.
    derived: bool = False
    #: The mode, where the input is a mode and a cell names one of its mode parameters (rule
    #: T): Fortran passes `constant`, its code evaluated in its own module, Python and R its string
    mode: ModeValue | None = None

    @property
    def is_constant(self) -> bool:
        return self.argument is None


@dataclass(frozen=True)
class Supplier:
    """What one `Supplied by` cell resolves to: the key a producer call is made of (rule K).

    `kind` is `argument` (`value` the consumer's argument, lower-cased), `constant` (`value` the
    evaluated constant with its Python type, so `.true.` is not `1`), `omitted` (rule O: no call
    passes the input), or `unresolved` (`value` the text, which `_resolve_output_from` reports).
    """

    kind: str
    value: object


#: The `Supplied by` cell that leaves a producer input out of the call (rule O), written bare and
#: in capitals as `AUTO` is
OMITTED = "OMITTED"


def resolve_cell(text: str, has_argument, evaluator: ConstantEvaluator | None) -> Supplier:
    """Resolve one `Supplied by` cell: `OMITTED`, the consumer's argument of that name, else a
    constant. `OMITTED` is its own kind of supplier (`omitted`), so a call omitting an input is
    another call than one supplying it.

    The one resolution of a cell. `roles` builds the plans with it, and synthesis keys a producer
    call with it (`producer_call_key`) before any plan exists, with the project's one evaluator
    (`Project.evaluator`), so the two group calls alike. `has_argument(name)` says whether the
    consumer has an argument of that name.
    """
    cell = text.strip()
    if cell == OMITTED:
        return Supplier("omitted", None)
    if has_argument(cell):
        return Supplier("argument", cell.lower())
    if evaluator is not None:
        try:
            value = evaluator.evaluate(cell)
        except ConstantError:
            pass
        else:
            return Supplier("constant", (type(value).__name__, value))
    return Supplier("unresolved", cell.lower())


def producer_input_rows(argument: Argument, conventions: Conventions = CONVENTIONS):
    """The `(producer input, cell)` rows of `argument`'s producer-input tables, as written."""
    return [
        (wanted.text.strip(), given.text.strip())
        for table in argument.doc.tables
        if _looks_like_a_producer_input_table(table, conventions)
        for wanted, given in table.rows
    ]


def producer_call_key(argument: Argument, arguments, conventions: Conventions = CONVENTIONS,
                      evaluator: ConstantEvaluator | None = None) -> tuple:
    """The producer call a `DM_OUTPUT_FROM` names, as synthesis keys it (rule K).

    The producer, and each table row resolved as `_resolve_output_from` resolves it
    (`resolve_cell`) against `arguments`, the consumer's; a row whose cell resolves to the
    name-matched argument is no row (`| n_items | N_ITEMS |`). Name-matched inputs need no part
    in it: one consumer matches them alike in every plan. Two arguments share a key exactly
    where their plans share `OutputFromPlan.call`, which `_check_calls_kept_whole` asserts.
    """
    directive = argument.directives.output_from
    names = {other.name.lower() for other in arguments}
    rows = set()
    for wanted, given in producer_input_rows(argument, conventions):
        supplier = resolve_cell(given, lambda name: name.lower() in names, evaluator)
        if supplier == Supplier("argument", wanted.lower()):
            continue
        rows.add((wanted.lower(), supplier))
    return directive.module.lower(), directive.procedure.lower(), frozenset(rows)


@dataclass(frozen=True)
class OutputFromPlan:
    """How an argument is obtained by calling another procedure.

    For `DM_OUTPUT_FROM(work_array_size, calc_paralog_subsets_bitmasks_size, ..., AUTO)`: call
    `calc_paralog_subsets_bitmasks_size`, supply its inputs from the consumer's own arguments,
    and take its `work_array_size` output. The
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
    def call(self) -> tuple:
        """What makes two plans one producer call: the producer, and every input it is given.

        One call fills every argument it outputs into -- but only arguments whose plans pass it
        the same inputs. Two arguments a table supplies differently (`a | a1`, `a | a2`) are two
        calls, each into its own argument. Synthesis keys a call before any plan exists, through
        the resolver this is built with (`producer_call_key`), so `foo` keeps or takes over a
        call whole (rule K). Python and R group their calls by the call they render, which can
        only merge plans whose binding inputs agree, never split one.
        """
        return id(self.producer), self.inputs

    @property
    def binding_inputs(self) -> tuple[ProducerInput, ...]:
        """What Python and R pass when they call the producer: every input but the derived.

        The producer's own binding works a derived input out from the arrays it is handed,
        so it is not one of its parameters there. A Fortran call passes them all.
        """
        return tuple(supply for supply in self.inputs if not supply.derived)

    @property
    def inout_feedback(self) -> tuple[ProducerInput, ...]:
        """The producer inputs it refines in place, that a consumer argument supplied.

        An `intent(inout)` producer input is returned alongside the annotated output with a
        possibly different value -- `calc_paralog_subsets_bitmasks_size` caps the
        `max_subset_size` it is handed to what actually fits. After the AUTO call the
        consumer must adopt that returned value, not the raw one it passed in, or it works
        with a stale figure. Constants are never fed back (nothing supplied them).
        """
        fed_back = []
        for supply in self.binding_inputs:
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

    #: `n_omp_threads`: the most threads an OpenMP root's loop may use, and the last extent of
    #: its per-thread work arrays. The allocating wrapper hands the implementation
    #: `M_OMP_N_THREADS(n_omp_threads)` for it, one without OpenMP; every other tier passes it on
    is_omp_threads: bool = False

    #: Named `bitmask(s)` or `*_bitmask(s)`: it holds bit-mask words, every int32 bit pattern
    #: of which is a word (`Conventions.is_bitmask_name`; `validate` refuses it on anything but
    #: an `integer(int32)` array)
    is_bitmask: bool = False

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
        return bool(self.readable_owners)

    @property
    def readable_owners(self) -> tuple[Argument, ...]:
        """The arrays naming this extent that it can be read off (`is_inferable_extent`)."""
        def readable(owner) -> bool:
            if owner.roles and owner.roles.is_temporary:
                return False
            if owner.roles and owner.roles.shape_arg is not None:
                # its extents travel in a shape argument, so the count is their product
                # whichever way the array itself flows
                return True
            # a missing intent is `validate`'s to report, and `validate` runs after this;
            # dereferencing it here ends the run in a traceback instead
            return owner.intent is not None and owner.intent.is_input

        return tuple(owner for owner in self.extent_of if readable(owner))

    @property
    def is_shape_arg(self) -> bool:
        return self.shape_of is not None

    @property
    def is_inferable_shape_arg(self) -> bool:
        """Whether the shape can be read off the array it describes.

        Only when that array is an input. Describing an `intent(out)` array, the shape is
        what the caller states the result should be, so it has to be asked for.
        """
        return (
            self.shape_of is not None
            and self.shape_of.intent is not None
            and self.shape_of.intent.is_input
        )

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
    evaluator = project.evaluator()
    # what the rules whose cause sits in an implementation's source have said, so each says it
    # once per implementation, not once per wrapper of a mode split (`_OncePerImplementation`)
    reported: set = set()
    for module in project:
        for procedure in module.exported_procedures:
            # once per generated pair: its expert half carries every argument the plain half has,
            # and resolves each against the implementation's own, so it reports for both. The
            # plain half's plans are still attached, for its binding, but what is wrong with them
            # is said once, where the whole picture is.
            expert = module.procedure(procedure.name + conventions.expert_suffix)
            plain_half = (expert is not None
                          and plain_of_an_expert(expert, project, conventions) is procedure)
            bag = DiagnosticBag() if plain_half else diagnostics
            implementation = _implementation_of(procedure, project, conventions)
            # the plain half's reports are thrown away, so it marks nothing as reported: else
            # where `foo` keeps a value, its expert's report would be suppressed with them
            once = _OncePerImplementation(bag, set() if plain_half else reported,
                                          id(implementation or procedure))
            _resolve_output_from(procedure, project, bag, conventions, evaluator, once)
            if not plain_half:
                _check_refined_values(procedure, project, diagnostics, conventions)
                _check_repeated_outputs(procedure, diagnostics)
                plain = plain_of_an_expert(procedure, project, conventions)
                if plain is not None:
                    _check_calls_kept_whole(procedure, plain, diagnostics)
    _check_producer_cycles(project, diagnostics, conventions)


class _OncePerImplementation:
    """Reports a diagnostic once per implementation: once, not once per wrapper of a mode split.

    For the rules whose cause sits in the implementation's own source -- a table row, the
    supplier a directive names -- and which every per-mode wrapper would otherwise find again,
    the same message at the same line. Keyed by the implementation, the location and the
    message, which name the implementation rather than the wrapper (4c), so they agree.
    """

    def __init__(self, diagnostics: DiagnosticBag, reported: set, implementation: int):
        self.diagnostics = diagnostics
        self.reported = reported
        self.implementation = implementation

    def error(self, message: str, entity=None, location: SourceLocation | None = None,
              note: str | None = None) -> None:
        where = location or (entity.location if entity is not None else None)
        key = (self.implementation, str(where), message)
        if key in self.reported:
            return
        self.reported.add(key)
        self.diagnostics.error(message, entity=entity, location=location, note=note)


def _check_calls_kept_whole(expert: Procedure, plain: Procedure,
                            diagnostics: DiagnosticBag) -> None:
    """Assert rule K: no producer call of a generated pair is kept in part and taken over in part.

    Synthesis decided `plain`'s signature by `producer_call_key` before any plan existed; this
    checks, on the plans, that it grouped the calls as `OutputFromPlan.call` does. With the
    table rows refused that the plans would ignore (N1), no source reaches it: a failure is the
    generator's, reported once per pair at the first kept value's directive.
    """
    for arguments in _producer_calls(expert):
        kept = [a for a in arguments if plain.argument(a.name) is not None]
        taken = [a for a in arguments if plain.argument(a.name) is None]
        if not kept or not taken:
            continue
        plan = kept[0].roles.computed_from
        directive = kept[0].directives.output_from
        diagnostics.error(
            f"internal error: '{plain.name}' keeps {', '.join(repr(a.name) for a in kept)} but "
            f"takes over {', '.join(repr(a.name) for a in taken)}, which one call of "
            f"'{plan.producer.name}' fills (pair '{plain.name}' / '{expert.name}')",
            entity=kept[0],
            location=_location_of(kept[0], directive.line_number if directive else None),
            note=("a producer call is taken over whole or kept whole (rule K); synthesis keyed "
                  "this call differently from its plans, which is a defect of the generator"),
        )


def _resolve_output_from(consumer: Procedure, project, diagnostics: DiagnosticBag,
                         conventions: Conventions = CONVENTIONS,
                         evaluator: ConstantEvaluator | None = None,
                         once: _OncePerImplementation | None = None) -> None:
    """Attach an `OutputFromPlan` to every argument documented with `DM_OUTPUT_FROM`.

    The producer is a published entry point: a hand-written `M_EXPORT_C` procedure, or the
    generated wrapper of an implementation -- the plain-named one, in the generated module.
    An implementation itself is not one; the error names the wrapper to spell instead.

    Each of the producer's inputs is supplied by the consumer argument of the same name.
    Where the two disagree, the argument documents the mapping in a table:

        !! | Producer input | Supplied by |
        !! |----------------|-------------|
        !! | n_elements     | n_genes     |

    A cell that is not an argument of the consumer is evaluated as a Fortran constant, so
    a producer input the consumer simply has no argument for -- `n_dim` on a routine that
    is always one-dimensional -- can be given a value outright.

    A producer input that is neither name-matched, nor mapped to an argument, nor a
    constant is an error naming it. An input the producer's binding derives (an extent read
    off one of its arrays) is resolved the same way and marked derived: Python and R leave
    it out of their call, and only a generated Fortran wrapper calling the producer passes
    it -- so only there is one that resolves to nothing an error. The producer may live in
    another module; the Python emitter imports it where it is called, and R has every
    wrapper in one environment.
    """
    if once is None:
        once = _OncePerImplementation(diagnostics, set(), id(consumer))
    called_here: dict[tuple, list[Argument]] = {}
    # where `consumer` is a generated pair's expert: the plain half, which `foo` is
    plain = plain_of_an_expert(consumer, project, conventions)
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
                note=_missing_producer_note(directive, project, conventions),
            )
            continue

        if not producer.is_exported:
            diagnostics.error(
                f"'{argument.name}' is computed from '{producer.name}', which is not "
                f"exported, so there is no wrapper to call",
                entity=argument,
                note=_not_exported_note(producer, project, conventions),
            )
            continue

        producer_plain = plain_of_an_expert(producer, project, conventions)
        if producer_plain is not None:
            diagnostics.error(
                f"'{argument.name}' is computed from '{producer.name}', the expert tier, "
                f"which leaves preparing its work to the caller",
                entity=argument,
                note=(f"name the plain-named wrapper, which prepares it itself: "
                      f"'{producer_plain.name}' in '{producer_plain.module.name}'"),
            )
            continue

        output = producer.argument(directive.argument)
        if output is None or not (output.intent and output.intent.is_output):
            diagnostics.error(
                f"'{directive.argument}' is not an output of '{producer.name}'",
                entity=argument,
            )
            continue

        fortran_calls = _is_called_from_fortran(consumer, argument, project, conventions)
        if _refuses_the_output(argument, directive, producer, output, fortran_calls,
                               diagnostics):
            continue

        renames = _producer_input_table(argument, consumer, producer,
                                        diagnostics, conventions, once, evaluator, project)
        if renames is None:
            continue

        if not directive.is_automatic:
            # JUST_INFO is documentation only: the caller still supplies the argument, and
            # nothing calls the producer. So the consumer need not be able to feed the
            # producer's inputs -- a pointer to a producer whose inputs live elsewhere (a
            # filter that needs data the consumer never sees) is exactly what JUST_INFO is
            # for. Record the relation, resolving no inputs; the producer existing, being
            # exported, and owning the named output (checked above) is all that must hold.
            # Its table is still checked by the row rules above: Ford publishes it.
            argument.roles.computed_from = OutputFromPlan(
                producer=producer, output=output, inputs=(), is_automatic=False,
            )
            continue

        inputs = []
        unmatched = []
        # (producer input, cell): a row naming neither an argument nor a constant (4a)
        unresolved = []
        # producer inputs a mode the implementation splits would supply (T, 4b)
        split_supplied = []
        split = _split_mode_of(consumer, project, conventions)
        # (producer input, supplier) where the supplier is an argument of `foo` alone: the
        # prologue's own, which the implementation never has
        prologue_only = []
        for producer_input in producer.arguments:
            if producer_input is output or producer_input.intent is None:
                continue
            if not producer_input.intent.is_input:
                continue
            if producer_input.name.lower() == "ierr":
                continue
            # the producer's own wrapper derives this input (an extent read off one of its
            # array arguments, a shape, a mask count, or a nested producer) and drops it from
            # its signature, so Python and R must not supply it -- but a Fortran call of the
            # producer has to, so it is resolved like any other input
            # a nested producer's value counts before its own plan is resolved: the directive
            # says so already, and the consumers are resolved in no particular order
            nested = producer_input.directives.output_from
            derived = bool(producer_input.roles and producer_input.roles.is_derived) or bool(
                nested is not None and nested.is_automatic)
            cell = renames.get(producer_input.name.lower())
            supplier = cell if cell is not None else producer_input.name
            # a row is resolved as synthesis keys the call (`resolve_cell`): the consumer's
            # argument, else a constant, which is how a producer input the consumer simply
            # does not have gets a value; without a row, the argument of the input's name
            if cell is not None:
                resolved = resolve_cell(cell, lambda name: consumer.argument(name) is not None,
                                        evaluator)
            elif consumer.argument(supplier) is not None:
                resolved = Supplier("argument", supplier.lower())
            else:
                resolved = Supplier("unresolved", supplier.lower())
            if resolved.kind == "argument":
                inputs.append(ProducerInput(producer_input.name,
                                            argument=consumer.argument(resolved.value).name,
                                            derived=derived))
                continue
            if resolved.kind == "constant":
                # a cell naming one of a mode input's own parameters is that mode (T, Q8): its
                # code from its own module, its string for Python and R
                mode = _mode_named(producer_input, cell)
                code = _mode_code(mode, project, evaluator) if mode is not None else None
                inputs.append(ProducerInput(
                    producer_input.name, constant=resolved.value[1] if code is None else code,
                    derived=derived, mode=mode))
                continue
            if resolved.kind == "omitted":
                continue  # rule O: no call passes it, in any tier
            if derived and not fortran_calls:
                continue  # nothing that calls the producer here passes it
            if plain is not None and plain.argument(supplier) is not None:
                prologue_only.append((producer_input.name, plain.argument(supplier).name))
                continue
            if split is not None and supplier.lower() == split.lower():
                split_supplied.append(producer_input.name)
            elif cell is not None:
                unresolved.append((producer_input.name, cell))
            else:
                unmatched.append(producer_input.name)

        if _refuses_what_nothing_supplies(argument, directive, consumer, producer, split,
                                          split_supplied, unresolved, unmatched, renames,
                                          fortran_calls, project, conventions, once):
            continue

        if _refuses_an_axis(argument, directive, consumer, producer, output, inputs,
                            fortran_calls, evaluator, diagnostics):
            continue
        if _refuses_chaining(argument, consumer, producer, inputs, fortran_calls, diagnostics):
            continue
        # before P1, which would refuse some of the same suppliers less precisely: an input
        # the producer refines has the narrower rule of its own
        if fortran_calls and _refuses_a_refined_input(argument, consumer, producer, inputs,
                                                      renames, diagnostics, conventions):
            continue
        if _refuses_a_supplier(argument, directive, consumer, plain, producer, inputs,
                               prologue_only, fortran_calls, project, conventions, diagnostics):
            continue
        if _refuses_a_kind(argument, directive, consumer, producer, inputs, renames,
                           fortran_calls, project, evaluator, diagnostics, once):
            continue

        argument.roles.computed_from = OutputFromPlan(
            producer=producer,
            output=output,
            inputs=tuple(inputs),
            is_automatic=directive.is_automatic,
        )
        if fortran_calls:
            called_here.setdefault(argument.roles.computed_from.call, []).append(argument)

    for arguments in called_here.values():
        _check_unused_outputs(arguments, diagnostics, conventions)


def _check_producer_cycles(project, diagnostics: DiagnosticBag,
                           conventions: Conventions) -> None:
    """Refuse a producer that, through its own producers, comes back to its consumer.

    Calling it recurses until the stack runs out: the generated wrapper calls its producer, which
    calls its own, and so on, and the bindings do the same through the producers' bindings. One
    procedure sizing an argument by itself is the shortest such chain.

    A generated pair is one procedure here: a producer is always the plain-named wrapper, and what
    it calls is what either tier of it names, the work it takes over sitting on its expert half.
    """
    edges: dict[int, list[tuple[Procedure, Argument]]] = {}
    nodes: dict[int, Procedure] = {}
    for module in project:
        for consumer in module.exported_procedures:
            entry = plain_of_an_expert(consumer, project, conventions) or consumer
            for argument in consumer.arguments:
                plan = argument.roles.computed_from if argument.roles else None
                if plan is None or not plan.is_automatic:
                    continue
                nodes[id(entry)] = entry
                nodes.setdefault(id(plan.producer), plan.producer)
                callees = edges.setdefault(id(entry), [])
                if all(callee is not plan.producer for callee, _ in callees):
                    callees.append((plan.producer, argument))

    reported: set[frozenset] = set()
    finished: set[int] = set()

    def visit(node: Procedure, chain: list[tuple[Procedure, Argument]]) -> None:
        # `chain`: each procedure on the way here, with its argument that calls the next one
        on_chain = [id(procedure) for procedure, _ in chain]
        for callee, argument in edges.get(id(node), []):
            here = chain + [(node, argument)]
            if id(callee) == id(node) or id(callee) in on_chain:
                cycle = here[(on_chain + [id(node)]).index(id(callee)):]
                key = frozenset(id(procedure) for procedure, _ in cycle)
                if key not in reported:
                    reported.add(key)
                    _report_cycle(cycle, diagnostics)
            elif id(callee) not in finished:
                visit(callee, here)
        finished.add(id(node))

    for node_id, node in nodes.items():
        if node_id not in finished:
            visit(node, [])


def _report_cycle(cycle: list[tuple[Procedure, Argument]], diagnostics: DiagnosticBag) -> None:
    """One diagnostic per cycle, on the argument where its chain starts, naming the whole chain.

    `cycle` holds each procedure on it with its argument that calls the next one; the last one
    calls the first again.
    """
    first, argument = cycle[0]
    chain = " -> ".join([procedure.name for procedure, _ in cycle] + [first.name])
    diagnostics.error(
        f"'{argument.name}' of '{first.name}' is computed from a producer that comes back to "
        f"it: {chain}",
        entity=argument,
        note=(
            "each call would make the next until the stack runs out; size one of them by "
            "something that does not depend on the others"
        ),
    )


#: What every refusal of a producer input's supplier says to do instead (rule P1)
PRODUCER_INPUT_NOTE = "supply it by an argument the caller passes, or a constant in the table"

#: Appended to a note where the generated `foo` calls the producer and the caller could supply
#: the value instead: JUST_INFO keeps the argument, and calls nothing
_JUST_INFO_WAY_OUT = (
    "; or, where the caller can supply the value, name the producer with JUST_INFO instead "
    "of AUTO, so the caller passes it"
)


def _way_out(fortran_calls: bool) -> str:
    return _JUST_INFO_WAY_OUT if fortran_calls else ""


def _rank_text(argument: Argument) -> str:
    return "a scalar" if argument.is_scalar else f"rank {argument.rank}"


def _logical_kind(fortran: FortranType) -> str | None:
    return (fortran.kind or "").lower() or None


def _types_agree(consumer: FortranType, produced: FortranType, fortran_calls: bool) -> bool:
    """Whether a consumer argument of type `consumer` can take a producer output of `produced`.

    Where the generated Fortran calls the producer it passes the consumer's local by keyword,
    so the two must be one Fortran type: the base type, a numeric kind by value (`c_int` is
    `int32`), a logical's kind by name (`c_bool` is not the default kind, which the bindings
    convert). Where a binding calls it, the value crosses as the C type the binding passes, so
    the two logical kinds are one.
    """
    if consumer.base is not produced.base:
        return False
    if consumer.is_character:
        return True
    if consumer.is_logical:
        return not fortran_calls or _logical_kind(consumer) == _logical_kind(produced)
    ours, theirs = _kind_value(consumer.kind), _kind_value(produced.kind)
    if ours is None or theirs is None:
        return (consumer.kind or "").lower() == (produced.kind or "").lower()
    return ours == theirs


def _refuses_the_output(argument: Argument, directive, producer: Procedure, output: Argument,
                        fortran_calls: bool, diagnostics: DiagnosticBag) -> bool:
    """Refuse a consumer argument that cannot hold the producer output it is computed from.

    The rows of the design's R2 that need only the two declarations, each with its own message,
    at the directive's line, and only the first that applies (row c, which needs the producer's
    inputs resolved, is `_refuses_an_axis`):

    - (e) a character array: its length would be one more thing to agree;
    - (i) an output that is not `intent(out)`: the producer would read what nothing set;
    - (g) a rank that differs, and (f) a type that differs -- on a JUST_INFO pointer too, since
      the caller is told to pass that output here; silent where `validate` refuses the type
      itself, so a type is reported once;
    - (a) an array `foo` allocates that is assumed-shape or assumed-size, with no extent to
      allocate it by;
    - (b) a producer output that is assumed-shape or assumed-size, or that it fills only partly
      (`DM_RESULT_SIZE_IS`), which an array declared with its own extents cannot hold. The
      partial fill is refused where something allocates the array with its declared extents:
      `foo`, and a binding of an explicit-shape consumer; an assumed-shape one takes what comes.

    True when something was refused.
    """
    location = _location_of(argument, directive.line_number)
    automatic = directive.is_automatic
    computed = f"'{argument.name}' is computed from '{producer.name}'"
    way_out = _way_out(fortran_calls)

    def refuse(message: str, note: str) -> bool:
        diagnostics.error(message, entity=argument, location=location, note=note)
        return True

    if automatic and argument.roles is not None and argument.roles.is_mode:
        return refuse(
            f"'{argument.name}' is a mode, so it cannot be computed with '{producer.name}'",
            "a mode is a choice the caller makes (Section 5.10); it comes from the caller, never "
            "from a producer" + way_out,
        )
    if automatic and argument.is_array and argument.type.is_character:
        return refuse(
            f"'{argument.name}' is a character array computed from '{producer.name}'",
            "a producer does not fill a character array: its length would be one more thing "
            "to agree with the producer's" + way_out,
        )
    if automatic and output.intent is not Intent.OUT:
        return refuse(
            f"{computed}, whose output '{output.name}' is intent({output.intent.value}), not "
            "intent(out)",
            "a producer's output is intent(out): "
            + ("the generated wrapper would pass it a local nothing has set, which the "
               "producer then reads" if fortran_calls else
               "nothing here has a value to hand it, which the producer would read first")
            + way_out,
        )
    if argument.rank != output.rank:
        return refuse(
            f"'{argument.name}' is {_rank_text(argument)}, but '{producer.name}' gives "
            f"'{output.name}' as {_rank_text(output)}",
            f"declare it with the rank '{output.name}' has",
        )
    silent = any(t.base is BaseType.DERIVED or bindings_cannot_carry(t)
                 for t in (argument.type, output.type))
    if not silent and not _types_agree(argument.type, output.type, fortran_calls):
        return refuse(
            f"'{argument.name}' is {argument.type}, but '{producer.name}' gives "
            f"'{output.name}' as {output.type}",
            f"declare it as '{output.name}' is declared"
            + (", kind and all: the generated wrapper passes it to the producer by keyword"
               if fortran_calls else ""),
        )
    if not automatic or not argument.is_array:
        return False
    if fortran_calls and not argument.dimension.is_explicit_shape:
        shape = "assumed-shape" if argument.dimension.is_assumed_shape else "assumed-size"
        return refuse(
            f"{computed}, and the generated wrapper allocates it before the call, but it is "
            f"{shape}, so there is no extent to allocate it with",
            f"declare its extents as '{output.name}' has them, in this procedure's names"
            + way_out,
        )
    if not output.dimension.is_explicit_shape:
        shape = "assumed-shape" if output.dimension.is_assumed_shape else "assumed-size"
        return refuse(
            f"{computed}, whose output '{output.name}' is {shape}",
            "an output a consumer takes is declared with its extents: whatever calls the "
            "producer allocates the array with the consumer's, which an output that takes its "
            "shape from the actual argument never names" + way_out,
        )
    allocated = fortran_calls or argument.dimension.is_explicit_shape
    counter = output.directives.result_size_is
    if allocated and counter is not None:
        return refuse(
            f"{computed}, which fills only the first '{counter.argument}' elements of "
            f"'{output.name}'",
            "the implementation would be handed elements the producer never wrote; have the "
            "producer fill the whole array" + way_out,
        )
    return False


def _extent_tokens(text: str, rename) -> list[tuple[str, object]] | None:
    """`text` as comparable tokens: names lower-cased and renamed, numbers as values.

    `rename(name)` gives a name's replacement -- a name, a value, or None where the name has
    nothing to compare by -- and the whole extent is then None.
    """
    from ..synthesize import extent_tokens  # local: synthesize imports this module

    tokens = []
    for kind, text_ in extent_tokens(text):
        if kind == "number":
            literal = text_.split("_")[0]
            tokens.append(("number", int(literal) if literal.isdigit() else literal.lower()))
        elif kind == "name":
            replacement = rename(text_)
            if replacement is None:
                return None
            tokens.append(replacement)
        else:
            tokens.append((kind, text_))
    return tokens


def _same_extent(left: list, right: list) -> bool:
    """Whether two token lists are one extent: as trees where both parse, else token by token."""
    from ..abi.extents import ExtentError, parse  # local: the ABI layer sits above the IR

    def text(tokens) -> str:
        return " ".join(str(value) for _, value in tokens)

    try:
        return parse(text(left)) == parse(text(right))
    except ExtentError:
        return left == right


def _renamers(consumer: Procedure, producer: Procedure, inputs, evaluator):
    """`(theirs, mine)`: how `_extent_tokens` reads a name of the producer's declarations, and one
    of the consumer's, so the two compare.

    A producer input is replaced by what supplies it here; a name that is no input is a constant
    of the producer's own scope, and one that is no argument of the consumer a constant of its --
    each evaluated where it is declared, so `5_int32` is `5`. `theirs` gives None for an input
    nothing here supplies (one only the producer's binding derives): nothing spells it.
    """
    supplied = {supply.name.lower(): supply for supply in inputs}
    producer_inputs = {dummy.name.lower() for dummy in producer.arguments}
    ours = {dummy.name.lower() for dummy in consumer.arguments}

    def constant(name: str):
        if evaluator is None:
            return None
        try:
            value = evaluator.evaluate(name)
        except ConstantError:
            return None
        return ("number", value) if isinstance(value, int) else None

    def theirs(name: str):
        supply = supplied.get(name.lower())
        if supply is not None:
            if supply.argument is not None:
                return ("name", supply.argument.lower())
            value = supply.constant
            return ("number", int(value)) if isinstance(value, int) else ("name", str(value))
        if name.lower() in producer_inputs:
            return None  # an input nothing here supplies
        return constant(name) or ("name", name.lower())

    def mine(name: str):
        if name.lower() in ours:
            return ("name", name.lower())
        return constant(name) or ("name", name.lower())

    return theirs, mine


def _refuses_an_axis(argument: Argument, directive, consumer: Procedure, producer: Procedure,
                     output: Argument, inputs: list[ProducerInput], fortran_calls: bool,
                     evaluator: ConstantEvaluator | None, diagnostics: DiagnosticBag) -> bool:
    """Refuse an AUTO array whose extent differs from the producer output's (R2 row c).

    The producer writes the shape it declares into an array allocated with the consumer's, so
    each axis must be the producer's extent with every producer input replaced by what supplies
    it here -- all at once, so a table swapping two (`n | m`, `m | n`) swaps them, inside
    `size(a, k)` too. A name that is no input is a constant of the producer's own scope, and one
    that is no argument of the consumer a constant of its: each evaluated where it is declared,
    so `5_int32` is `5` and a named parameter compares by its value. Names compare in any case.
    Where `foo` allocates it, and a binding of an explicit-shape consumer; an assumed-shape one
    takes the shape it is given. An axis naming an input nothing here supplies (one only the
    producer's binding derives, on the binding route) is not compared: there is nothing to
    spell it by. Rank is row g's, and a shape the producer does not declare row b's. True when
    something was refused.
    """
    if not (directive.is_automatic and argument.is_array):
        return False
    if not (fortran_calls or argument.dimension.is_explicit_shape):
        return False
    theirs, mine = _renamers(consumer, producer, inputs, evaluator)
    for axis, (extent, produced) in enumerate(
            zip(argument.dimension.extents, output.dimension.extents), start=1):
        expected = _extent_tokens(produced, theirs)
        if expected is None:
            continue
        given = _extent_tokens(extent, mine)
        if _same_extent(given, expected):
            continue
        spelled = " ".join(str(value) for _, value in expected)
        diagnostics.error(
            f"'{argument.name}' is computed from '{producer.name}', but its axis {axis} is "
            f"'{extent}', where '{output.name}' is declared with '{produced}', which is "
            f"'{spelled}' here",
            entity=argument,
            location=_location_of(argument, directive.line_number),
            note=(
                "the producer writes the shape it declares into the array, so spell each "
                "extent as it does, in this procedure's names: its inputs replaced by what "
                "supplies them" + _way_out(fortran_calls)
            ),
        )
        return True
    return False


def _mode_named(producer_input: Argument, cell: str | None) -> ModeValue | None:
    """The mode of `producer_input`'s own table whose parameter the cell names, or None."""
    table = producer_input.roles.mode if producer_input.roles is not None else None
    if table is None or cell is None:
        return None
    return next((value for value in table.values
                 if value.parameter.lower() == cell.strip().lower()), None)


def _mode_code(value: ModeValue, project, evaluator: ConstantEvaluator | None):
    """A mode's code: its parameter evaluated in the module its table row links (T).

    Not looked up by name in the project's flat table of constants, where two modules' parameters
    of one name are one entry. None where it cannot be evaluated.
    """
    module = project.module(value.module) if project is not None else None
    if module is None or evaluator is None:
        return None
    for parameter in module.parameters:
        if parameter.name.lower() == value.parameter.lower() and parameter.expression:
            try:
                return evaluator.evaluate(parameter.expression)
            except ConstantError:
                return None
    return None


def _same_mode(left: ModeValue, right: ModeValue, project, evaluator) -> bool:
    """Whether two rows of mode tables are one mode with one value: one parameter, or two that
    evaluate, each in its own module, to one code."""
    if (left.module.lower(), left.parameter.lower()) == (right.module.lower(),
                                                         right.parameter.lower()):
        return True
    code = _mode_code(left, project, evaluator)
    return code is not None and code == _mode_code(right, project, evaluator)


def _quoted(strings) -> str:
    return ", ".join(f"'{string}'" for string in strings)


#: The note every type, kind or rank refusal of an argument supplier carries (T)
_SUPPLY_ITS_KIND = "supply it by an argument of its type, kind and rank, or a constant in the table"

#: The note every refusal of a constant supplier carries (T)
_CONSTANT_OF_ITS_TYPE = ("a constant supplies a scalar input of its own type, written as Fortran "
                         "writes it (3.0_real64, .false., 'abc')")


def _refuses_a_kind(argument: Argument, directive, consumer: Procedure, producer: Procedure,
                    inputs: list[ProducerInput], renames, fortran_calls: bool, project,
                    evaluator: ConstantEvaluator | None, diagnostics: DiagnosticBag,
                    once: _OncePerImplementation | None = None) -> bool:
    """Refuse a producer input given a value of another kind than its own (rule T).

    Judged where the generated Fortran calls the producer for every input it passes, whether or
    not a binding calls it too; where only a binding does, for `binding_inputs`. Silent where
    `validate` refuses the type itself (`bindings_cannot_carry`, a derived type), so a type is
    reported once. One error per failing input, at most one cause each, in this order:

    1. modes: a mode and a non-mode never meet; a mode argument supplies only a mode input whose
       table names each of its modes with the same value, each evaluated in its own module (the
       producer may accept more); a constant given to a mode input names one of its parameters;
    2. an argument's base type and rank -- and where the Fortran calls the producer, its kind
       (numeric by value, logical by name), a fixed length's value, and an explicit shape's
       extents with every producer input replaced by what supplies it (`_renamers`). The
       binding route carries one kind per base type, reads a string's length and an array's
       extents off the value, so only base type and rank are compared there;
    3. a constant: a scalar input, of its own type -- a logical for a logical, an integer that
       fits the kind, a real for a real (an integer is refused, not converted), a string no
       longer than a fixed length.

    An argument supplier at the directive's line, once per pair; a constant, which a table row
    gives, at that row, once per implementation, as every refusal one row causes. It runs after P1,
    so a supplier P1 refused is not judged here. True when something was refused.
    """
    once = once or _OncePerImplementation(diagnostics, set(), id(consumer))
    location = _location_of(argument, directive.line_number)
    computed = f"'{argument.name}' is computed from '{producer.name}'"
    theirs, mine = _renamers(consumer, producer, inputs, evaluator)
    refused = False
    for supply in inputs:
        if supply.derived and not fortran_calls:
            continue  # not passed where only a binding calls the producer
        wanted = producer.argument(supply.name)
        if wanted is None:
            continue
        if supply.is_constant:
            cause = _a_constant_of_another_kind(supply, wanted, renames, evaluator)
        else:
            cause = _an_argument_of_another_kind(supply, wanted, consumer.argument(supply.argument),
                                                 inputs, fortran_calls, theirs, mine, project,
                                                 evaluator)
        if cause is None:
            continue
        message, note = cause
        if supply.is_constant:
            once.error(f"{computed}, whose {message}", entity=argument, note=note,
                       location=_location_of(argument, renames.lines.get(supply.name.lower())))
        else:
            diagnostics.error(f"{computed}, whose {message}", entity=argument,
                              location=location, note=note + _way_out(fortran_calls))
        refused = True
    return refused


def _carried(*types: FortranType) -> bool:
    """Whether every type is one `validate` lets through, so T may judge them."""
    return not any(t.base is BaseType.DERIVED or bindings_cannot_carry(t) for t in types)


def _mode_parameters(table: ModeTable) -> str:
    return ", ".join(value.parameter for value in table.values)


def _an_argument_of_another_kind(supply, wanted: Argument, given: Argument | None, inputs,
                                 fortran_calls: bool, theirs, mine, project, evaluator):
    """Why the argument `given` cannot be passed for the producer input `wanted` (T), as
    (message, note), or None."""
    if given is None or not _carried(wanted.type, given.type):
        return None
    wanted_mode = wanted.roles.mode if wanted.roles is not None else None
    given_mode = given.roles.mode if given.roles is not None else None
    if given_mode is not None and wanted_mode is None:
        return (f"input '{wanted.name}' would be the mode '{given.name}', which Python and R pass "
                "as its string and Fortran as its code",
                "a mode is a choice, not a number (Section 5.10): it reaches a producer only as a "
                f"mode, whose table names the same modes. Supply '{wanted.name}' by an argument "
                "or a constant in the table")
    if wanted_mode is not None and given_mode is None:
        return (f"input '{wanted.name}' is a mode, but would be '{given.name}', which is not",
                "supply a mode input by a mode argument whose table names the same modes, or by "
                f"one of its mode parameters in the table ({_mode_parameters(wanted_mode)})")
    if wanted_mode is not None:
        accepts = _quoted(value.string for value in wanted_mode.values)
        missing = [value.string for value in given_mode.values
                   if wanted_mode.value_for(value.string) is None]
        other = [value.string for value in given_mode.values
                 if wanted_mode.value_for(value.string) is not None
                 and not _same_mode(value, wanted_mode.value_for(value.string), project,
                                    evaluator)]
        why = (f"which also accepts {_quoted(missing)}" if missing else
               f"which gives {_quoted(other)} another value" if other else None)
        if why is None:
            return None
        return (f"mode input '{wanted.name}' accepts {accepts}, but would be the mode "
                f"'{given.name}', {why}",
                "a mode passes through as its string in Python and R and as its code in Fortran, "
                "so the producer's table must name every mode of the supplier with the same value")
    if (given.rank != wanted.rank
            or not (_types_agree(wanted.type, given.type, True) if fortran_calls
                    else wanted.type.base is given.type.base)):
        return (f"input '{wanted.name}' is {wanted.type}, {_rank_text(wanted)}, but would be "
                f"'{given.name}', {given.type}, {_rank_text(given)}", _SUPPLY_ITS_KIND)
    if not fortran_calls:
        return None
    length = wanted.type.length
    if length is not None and not length.is_assumed:
        expected = _extent_tokens(length.expr, theirs)
        given_length = given.type.length
        same = (given_length is not None and not given_length.is_assumed and expected is not None
                and _same_extent(_extent_tokens(given_length.expr, mine), expected))
        if expected is not None and not same:
            return (f"input '{wanted.name}' is {wanted.type}, but would be '{given.name}', "
                    f"{given.type}",
                    "a fixed-length input takes a value of that length; declare the producer's "
                    f"input len=*, or give '{given.name}' its length")
    if wanted.is_array and wanted.dimension.is_explicit_shape:
        for extent, own in zip(wanted.dimension.extents, given.dimension.extents):
            expected = _extent_tokens(extent, theirs)
            if expected is None:
                continue  # it names an input nothing here supplies
            if _same_extent(_extent_tokens(own, mine), expected):
                continue
            from ..synthesize import size_names_in  # local: synthesize imports this module

            supplied = {s.name.lower(): s for s in inputs}
            named = [n.name for n in size_names_in(extent) if n.name.lower() in supplied]
            given_as = {n: supplied[n.lower()].argument or _as_written(supplied[n.lower()])
                        for n in named}
            with_ = (", with " + ", ".join(f"'{n}' given '{s}'" for n, s in given_as.items())
                     if given_as else "")
            passes = ", ".join(f"'{n} = {s}'" for n, s in given_as.items())
            return (f"input '{wanted.name}' has shape ({', '.join(wanted.dimension.extents)})"
                    f"{with_}, but would be '{given.name}', of shape "
                    f"({', '.join(given.dimension.extents)})",
                    f"the generated wrapper passes {passes or 'its extents'} beside "
                    f"'{wanted.name} = {given.name}', so a different '{own}' reads past "
                    f"'{given.name}': supply {_quoted(named) or 'its extents'} by what sizes "
                    f"'{given.name}'")
    return None


def _as_written(supply: ProducerInput) -> str:
    """How a constant supply reads in a message: as Fortran would write it."""
    value = supply.constant
    if isinstance(value, bool):
        return ".true." if value else ".false."
    return str(value)


def fixed_length(argument: Argument, evaluator: ConstantEvaluator | None) -> int | None:
    """A character argument's fixed length -- a literal, or a constant (`len=LABEL_LEN`) evaluated
    where the argument is declared -- or None for `len=*` or a length it cannot evaluate.

    Not by name in the project's flat table of constants, where two modules' parameters of one
    name are one entry: a consumer's `LABEL_LEN` would size the producer's input. A generated
    wrapper is declared with its implementation's names, so its implementation module counts as
    its own (`_in_its_module`).
    """
    length = argument.type.length if argument.type is not None else None
    if length is None or length.expr is None:
        return None
    if length.is_constant:
        return int(length.expr)
    value = _in_its_module(length.expr, argument, evaluator)
    return value if isinstance(value, int) and not isinstance(value, bool) else None


def _in_its_module(expression: str, argument: Argument, evaluator: ConstantEvaluator | None):
    """`expression` evaluated with the parameters of the module `argument` is declared in -- and,
    for a generated wrapper, of its implementation module -- before the project's others; None
    where it is not a constant."""
    if evaluator is None:
        return None
    procedure = argument.parent
    module = procedure.module if procedure is not None else None
    project = module.project if module is not None else None
    constants = dict(evaluator.constants)
    if project is not None:
        own = [module, project.module(module.name + CONVENTIONS.impl_suffix)]
        local = ConstantEvaluator(constants, evaluator.intrinsics)
        parameters = [p for m in own if m is not None for p in m.parameters if p.expression]
        for _ in parameters:  # as many rounds as there are parameters: one may name another
            for parameter in parameters:
                try:
                    local.constants[parameter.name.lower()] = local.evaluate_with_kind(
                        parameter.expression)
                except ConstantError:
                    continue
        constants = local.constants
    try:
        return ConstantEvaluator(constants, evaluator.intrinsics).evaluate(expression)
    except ConstantError:
        return None


def _a_constant_of_another_kind(supply, wanted: Argument, renames,
                                evaluator: ConstantEvaluator | None = None):
    """Why the table's constant cannot be passed for the producer input `wanted` (T), as
    (message, note), or None."""
    if not _carried(wanted.type):
        return None
    cell = renames.get(wanted.name.lower(), _as_written(supply))
    wanted_mode = wanted.roles.mode if wanted.roles is not None else None
    if wanted_mode is not None:
        if supply.mode is not None:
            return None
        return (f"input '{wanted.name}' is a mode, but would be '{cell}', which is not",
                "supply a mode input by a mode argument whose table names the same modes, or by "
                f"one of its mode parameters in the table ({_mode_parameters(wanted_mode)})")
    given = f"input '{wanted.name}' is {wanted.type}"
    if wanted.is_array:
        return (f"{given}, {_rank_text(wanted)}, but the table gives it the constant {cell}, a "
                "scalar", _CONSTANT_OF_ITS_TYPE)
    value = supply.constant
    base = wanted.type.base
    kinds = {bool: "a logical", int: "an integer", float: "a real", str: "a string",
             complex: "a complex"}
    what = kinds.get(type(value), "of another type")
    expected = {BaseType.LOGICAL: bool, BaseType.INTEGER: int, BaseType.REAL: float,
                BaseType.CHARACTER: str, BaseType.COMPLEX: complex}.get(base)
    if expected is None or type(value) is not expected:
        return (f"{given}, but the table gives it the constant {cell}, {what}",
                _CONSTANT_OF_ITS_TYPE)
    if base is BaseType.INTEGER:
        if not fits_integer_kind(value, _kind_value(wanted.type.kind) or 4):
            return (f"{given}, but the table gives it the constant {cell}, which does not fit",
                    _CONSTANT_OF_ITS_TYPE)
    length = fixed_length(wanted, evaluator) if base is BaseType.CHARACTER else None
    if length is not None and len(value) > length:
        return (f"{given}, but the table gives it the constant {cell}, which is longer than "
                f"{length}", _CONSTANT_OF_ITS_TYPE)
    return None


def _implementation_of(wrapper: Procedure, project, conventions: Conventions):
    """The implementation a generated wrapper calls, or None for any other procedure.

    `<base>_impl` beside it in the `_impl` module, or -- for a per-mode wrapper, named from its
    implementation's mode table -- the implementation whose table names it.
    """
    module = wrapper.module
    impl_module = project.module(module.name + conventions.impl_suffix) if module else None
    if impl_module is None:
        return None
    base = wrapper.name
    if base.lower().endswith(conventions.expert_suffix):
        base = base[: -len(conventions.expert_suffix)]
    impl = impl_module.procedure(base + conventions.impl_suffix)
    if impl is not None:
        return impl
    for candidate in impl_module.procedures:
        for dummy in candidate.arguments:
            mode = dummy.roles.mode if dummy.roles is not None else None
            if mode is not None and any(v.procedure_name.lower() == base.lower()
                                        for v in mode.values):
                return candidate
    return None


def _prologue_written(wrapper: Procedure, project, conventions: Conventions) -> set[str]:
    """What the prologue of a generated wrapper's implementation writes, lower-cased."""
    impl = _implementation_of(wrapper, project, conventions)
    directive = impl.directives.prologue if impl is not None else None
    prologue = (project.procedure(directive.module, directive.procedure)
                if directive is not None else None)
    if prologue is None:
        return set()
    return {dummy.name.lower() for dummy in prologue.arguments
            if dummy.intent is not None and dummy.intent.is_output}


def _refuses_a_supplier(argument: Argument, directive, consumer: Procedure, plain,
                        producer: Procedure, inputs: list[ProducerInput], prologue_only: list,
                        fortran_calls: bool, project, conventions: Conventions,
                        diagnostics: DiagnosticBag) -> bool:
    """Refuse a producer input that nothing has supplied by the time the producer is called (P1).

    A producer input that a call passes is supplied by a scalar constant, or by an input of the
    implementation that the caller passes, in every tier that calls the producer: in Fortran an
    `intent(in)` or `intent(inout)` dummy `foo` keeps -- required; optional with a `DM_DEFAULT`,
    which `foo` resolves before the call; or optional where the producer's input is optional
    too, passed through so an absent value stays absent -- and in a binding a parameter, or what
    the binding reads off one before the producers run. Refused, each with its own message: a
    work array, a permutation `foo` seeds and sorts, an output, an optional without a default
    feeding a required input, and an argument only the prologue has.

    Not judged here, as their own rules: a value another producer computes (chaining,
    `_refuses_chaining`), and one the prologue writes (`validate`'s rule that nothing above the
    prologue reads what it writes, which then still finds this plan). The supplier is resolved
    among the implementation's arguments, which the expert half of a pair has, so a pair is
    checked once, there. An input only the producer's binding derives is passed by a Fortran
    call alone. True when something was refused.
    """
    location = _location_of(argument, directive.line_number)
    note = PRODUCER_INPUT_NOTE + _way_out(fortran_calls)

    def refuse(name: str, supplier: str, why: str) -> bool:
        diagnostics.error(
            f"'{argument.name}' is computed from '{producer.name}', whose input '{name}' would "
            f"be '{supplier}', {why}",
            entity=argument, location=location, note=note,
        )
        return True

    for name, supplier in prologue_only:
        return refuse(name, supplier, "which only the prologue takes: it is the prologue's "
                                      "argument, not the implementation's")
    written = _prologue_written(consumer, project, conventions) if plain is not None else set()
    for supply in inputs:
        if supply.argument is None or (supply.derived and not fortran_calls):
            continue
        dummy = consumer.argument(supply.argument)
        if dummy is None:
            continue
        own = dummy.directives.output_from
        if own is not None and own.is_automatic:
            continue  # chaining: `_refuses_chaining`
        if dummy.name.lower() in written:
            continue  # the prologue's rule: `validate._check_prologue_outputs_are_not_read_first`
        if dummy.name.lower().startswith(conventions.temporary_prefix):
            return refuse(supply.name, dummy.name,
                          "a work array, which holds nothing before the implementation fills it")
        if (plain is not None and plain.argument(dummy.name) is None
                and dummy.name.lower().endswith(conventions.perm_suffix)):
            return refuse(supply.name, dummy.name,
                          f"a permutation '{plain.name}' seeds and sorts only after its "
                          "producers have run")
        if dummy.intent is Intent.OUT:
            return refuse(supply.name, dummy.name,
                          "an output, which has no value before the implementation sets it")
        wanted = producer.argument(supply.name)
        if (dummy.optional and dummy.directives.default is None
                and wanted is not None and not wanted.optional):
            return refuse(supply.name, dummy.name,
                          "an optional without a DM_DEFAULT, which reaches the producer's "
                          "required input absent when the caller leaves it out")
    return False


def _refuses_chaining(argument: Argument, consumer: Procedure, producer: Procedure,
                      inputs: list[ProducerInput], fortran_calls: bool,
                      diagnostics: DiagnosticBag) -> bool:
    """Refuse a producer input supplied by another AUTO `DM_OUTPUT_FROM` value of the consumer.

    Chaining producers, so one feeds another, is rejected (`design/impl-layer.md`, "Rejected --
    chaining producers"): every producer input comes from the caller, as an argument or a
    constant. On either route, and for a value given back to its own producer too -- but only
    for an input that is passed: one the producer's binding derives is passed by a Fortran call
    alone. True when something was refused.
    """
    refused = False
    for supply in inputs:
        if supply.derived and not fortran_calls:
            continue
        supplier = consumer.argument(supply.argument) if supply.argument is not None else None
        directive = supplier.directives.output_from if supplier is not None else None
        if directive is None or not directive.is_automatic:
            continue
        diagnostics.error(
            f"'{argument.name}' is computed from '{producer.name}', whose input "
            f"'{supply.name}' would be '{supplier.name}', which a producer computes here too",
            entity=argument,
            # at the directive, like every rule about what a DM_OUTPUT_FROM may name
            location=_location_of(argument, argument.directives.output_from.line_number),
            note=(
                "producers do not feed each other: size for an upper bound instead (the sizing "
                "formulas are monotone, so a buffer sized for the bound fits), or take the value "
                "as an argument of the caller's"
                + (" -- a required intent(inout) one, since the producer refines it in place"
                   if fortran_calls and producer.argument(supply.name).intent is Intent.INOUT
                   else "")
            ),
        )
        refused = True
    return refused


def _refuses_a_refined_input(argument: Argument, consumer: Procedure, producer: Procedure,
                             inputs: list[ProducerInput], renames: dict[str, str],
                             diagnostics: DiagnosticBag, conventions: Conventions) -> bool:
    """Refuse an `intent(inout)` producer input not supplied by a required `intent(inout)` dummy.

    The rule, where the generated Fortran wrapper calls the producer itself: an input the
    producer refines in place is passed through a required `intent(inout)` argument of the
    consumer, and nothing else. That is what the call can write, and what carries the refined
    value on: to the caller, and to the allocations and the implementation call after it (no
    other producer call may read it, `_check_refined_values`). An `intent(in)` argument or a
    constant does not compile there, and the local an optional is materialised into would
    swallow the refined value. Python and R adopt the refined value themselves
    (`inout_feedback`), so where only they call the producer nothing is refused. True when
    something was.
    """
    plain = consumer.module.procedure(consumer.name[: -len(conventions.expert_suffix)])
    refused = False
    for supply in inputs:
        if producer.argument(supply.name).intent is not Intent.INOUT:
            continue
        dummy = plain.argument(supply.argument) if supply.argument is not None else None
        if dummy is not None and dummy.intent is Intent.INOUT and not dummy.optional:
            continue
        if supply.argument is None:
            given = f"the constant {renames.get(supply.name.lower(), supply.constant)}"
        elif dummy is None:
            given = f"'{supply.argument}', a work array or other value the wrapper prepares"
        else:
            intent = f"intent({dummy.intent.value})" if dummy.intent else "no intent"
            given = f"{'optional ' if dummy.optional else ''}{intent} '{supply.argument}'"
        diagnostics.error(
            f"'{argument.name}' is computed from '{producer.name}', which refines its input "
            f"'{supply.name}' in place, so it must be passed through a required intent(inout) "
            f"argument; it is given {given}",
            entity=argument,
            note=(
                f"the generated wrapper calls '{producer.name}' itself, and writes "
                f"'{supply.name}' back through what supplies it: supply it by a required "
                "intent(inout) argument of the implementation, which the refined value then "
                "reaches, and the caller with it"
            ),
        )
        refused = True
    return refused


def _producer_calls(consumer: Procedure) -> list[list[Argument]]:
    """The consumer's AUTO arguments, grouped by the producer call that fills them."""
    calls: dict[tuple, list[Argument]] = {}
    for argument in consumer.arguments:
        plan = argument.roles.computed_from if argument.roles else None
        if plan is not None and plan.is_automatic:
            calls.setdefault(plan.call, []).append(argument)
    return list(calls.values())


def _check_repeated_outputs(consumer: Procedure, diagnostics: DiagnosticBag) -> None:
    """Refuse two arguments taking the same output of one producer call.

    One call -- the same producer, given the same inputs -- fills each output once, so the two
    would hold one value; the generated Fortran wrapper wrote only one of them. Either route.
    """
    for arguments in _producer_calls(consumer):
        taken: dict[str, Argument] = {}
        for argument in arguments:
            plan = argument.roles.computed_from
            first = taken.setdefault(plan.output.name.lower(), argument)
            if first is argument:
                continue
            diagnostics.error(
                f"'{argument.name}' takes the same output '{plan.output.name}' of "
                f"'{plan.producer.name}', given the same inputs, as '{first.name}'",
                entity=argument,
                note=(f"one call fills both with one value: use '{first.name}' alone, since "
                      f"'{argument.name}' would hold the same value"),
            )


def _check_refined_values(consumer: Procedure, project, diagnostics: DiagnosticBag,
                          conventions: Conventions) -> None:
    """Refuse a value a producer refines in place that feeds any other producer call too.

    The rule: a value a producer refines -- an `intent(inout)` input -- feeds
    that one producer call and no other, on either route. Whether another call reading it saw
    it before or after the refinement would depend on the call order, and a second refiner would
    refine what the first already had; so neither is allowed. One call filling several outputs
    from it is still one call.
    """
    calls = _producer_calls(consumer)

    def refined_by(arguments) -> set[str]:
        plan = arguments[0].roles.computed_from
        return {s.argument.lower() for s in plan.inputs if s.argument is not None
                and plan.producer.argument(s.name).intent is Intent.INOUT}

    def reads(arguments) -> set[str]:
        return {s.argument.lower() for s in arguments[0].roles.computed_from.inputs
                if s.argument is not None}

    reported: set[tuple] = set()
    for refiner in calls:
        for value in sorted(refined_by(refiner)):
            for other in calls:
                if other is refiner or value not in reads(other):
                    continue
                key = (value, frozenset({id(refiner[0]), id(other[0])}))
                if key in reported:
                    continue
                reported.add(key)
                name = consumer.argument(value).name
                refining = refiner[0].roles.computed_from.producer.name
                reader = other[0].roles.computed_from.producer.name
                refines = value in refined_by(other)
                verb = "refines it too" if refines else "reads it too"
                # where the Fortran wrapper calls a refiner, its input must be a required
                # intent(inout) argument (`_refuses_a_refined_input`); Python and R take any
                fortran = refines and _is_called_from_fortran(consumer, other[0], project,
                                                              conventions)
                diagnostics.error(
                    f"'{name}' is refined in place by '{refining}' (computing "
                    f"'{refiner[0].name}'), so it may feed no other producer call, but "
                    f"'{reader}' (computing '{other[0].name}') {verb}",
                    entity=other[0],
                    note=(
                        f"a value a producer refines feeds that producer alone: give "
                        f"'{reader}' its own argument of the caller's for that input"
                        + (" -- a required intent(inout) one, since it refines it as well"
                           if fortran else "")
                    ),
                )


def _is_called_from_fortran(consumer: Procedure, argument: Argument, project,
                            conventions: Conventions) -> bool:
    """Whether a generated Fortran wrapper calls the producer of `argument` itself.

    It does where the allocating tier took the argument over: `consumer` is the expert half
    of a generated pair, and the plain-named half beside it no longer has the argument.
    Everywhere else the caller supplies the value, or Python and R compute it through the
    producer's binding.
    """
    module = consumer.module
    if module is None or project.module(module.name + conventions.impl_suffix) is None:
        return False  # a hand-written module: nothing generated calls anything here
    suffix = conventions.expert_suffix
    if not consumer.name.lower().endswith(suffix):
        return False
    plain = module.procedure(consumer.name[: -len(suffix)])
    return plain is not None and plain.argument(argument.name) is None


def unused_outputs(arguments: list[Argument],
                   conventions: Conventions = CONVENTIONS) -> list[Argument]:
    """The outputs of the producer `arguments` are computed from that none of them takes.

    `arguments` are a consumer's AUTO values from one producer, which one call fills. A
    direct Fortran call must pass the rest of the producer's outputs all the same.
    """
    plan = arguments[0].roles.computed_from
    used = {a.roles.computed_from.output.name.lower() for a in arguments}
    error = conventions.error_arg.lower()
    return [
        dummy for dummy in plan.producer.arguments
        if dummy.intent is Intent.OUT and dummy.name.lower() not in used | {error}
    ]


def fits_a_local(argument: Argument) -> bool:
    """Whether a wrapper can declare a local of this output's type as it is spelled.

    A scalar can, as long as nothing in its declaration names the producer's other
    arguments: an array's extents do, and so does a string's length unless it is a literal.
    """
    length = argument.type.length
    return not argument.is_array and (length is None or length.is_constant)


def _check_unused_outputs(arguments: list[Argument], diagnostics: DiagnosticBag,
                          conventions: Conventions) -> None:
    """Refuse an output of a Fortran-called producer that nothing takes and no local can hold.

    A direct call has to pass every output. One the consumer takes goes into the value it
    computes; a scalar it does not take goes into a local and is dropped (`fits_a_local`).
    An array cannot: its local would need the producer's extents, spelled in the consumer's
    arguments.
    """
    producer = arguments[0].roles.computed_from.producer
    for dummy in unused_outputs(arguments, conventions):
        if fits_a_local(dummy):
            continue
        diagnostics.error(
            f"'{arguments[0].name}' is computed from '{producer.name}', whose output "
            f"'{dummy.name}' nothing here takes, and it is not a scalar the generated "
            f"wrapper can drop",
            entity=arguments[0],
            note=(
                "the generated wrapper calls the producer itself and must pass every output; "
                "a scalar it does not need goes into a local, but an array (or a string whose "
                "length is not a literal) would need the producer's own sizes, which nothing "
                "here can spell"
            ),
        )


def _missing_producer_note(directive, project, conventions: Conventions) -> str | None:
    """Where the producer a directive names is, when it names the right thing in the wrong place.

    The two mix-ups of an implementation and its wrapper: `foo_impl` with the generated module,
    and `foo` with the implementation module.
    """
    suffix = conventions.impl_suffix
    implementation = project.procedure(directive.module + suffix, directive.procedure)
    if implementation is not None:
        return _not_exported_note(implementation, project, conventions)
    if directive.module.lower().endswith(suffix):
        generated = directive.module[: -len(suffix)]
        wrapper = project.procedure(generated, directive.procedure)
        if wrapper is not None:
            # the plain-named one, even for `foo_expert`: that is what a producer names
            wrapper = plain_of_an_expert(wrapper, project, conventions) or wrapper
            return (f"'{wrapper.name}' is generated into '{generated}', not into "
                    f"'{directive.module}': name it there")
    return None


def plain_of_an_expert(producer: Procedure, project, conventions: Conventions):
    """The plain-named wrapper beside `producer`, if `producer` is a generated `foo_expert`."""
    module = producer.module
    suffix = conventions.expert_suffix
    if module is None or not producer.name.lower().endswith(suffix):
        return None
    if project.module(module.name + conventions.impl_suffix) is None:
        return None  # a hand-written pair: its `_expert` is the author's own name
    return module.procedure(producer.name[: -len(suffix)])


def _not_exported_note(producer: Procedure, project, conventions: Conventions) -> str:
    """What to name instead of a producer that is not exported.

    An implementation never is, and marking it `M_EXPORT_C` is refused as well: what is
    published is its generated wrapper, so that is what the directive names. A hand-written
    procedure is exported by the marker.
    """
    suffix = conventions.impl_suffix
    module = producer.module.name if producer.module is not None else ""
    if not (module.lower().endswith(suffix) and producer.name.lower().endswith(suffix)):
        return f"mark '{producer.name}' {conventions.export_marker} to export it"
    generated, base = module[: -len(suffix)], producer.name[: -len(suffix)]
    if project.procedure(generated, base) is not None:
        return (f"an implementation is published as its generated wrapper: name '{base}' "
                f"in '{generated}'")
    per_mode = [
        value.procedure_name
        for argument in producer.arguments
        if argument.roles is not None and argument.roles.mode is not None
        and argument.roles.mode.is_split
        for value in argument.roles.mode.values
    ]
    if per_mode:
        names = ", ".join(f"'{name}'" for name in per_mode)
        return (f"'{producer.name}' splits by mode, into one wrapper per mode: name the one "
                f"you mean in '{generated}' ({names})")
    return f"an implementation is published as its generated wrapper, in '{generated}'"


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


class _ProducerInputTable(dict):
    """A producer-input table: producer input (lower-cased) -> its `Supplied by` cell, as written.

    `lines` gives each row's line, where a diagnostic about the row points.
    """

    def __init__(self):
        super().__init__()
        self.lines: dict[str, int | None] = {}


def _producer_input_table(argument: Argument, consumer: Procedure, producer: Procedure,
                          diagnostics: DiagnosticBag, conventions: Conventions,
                          once: _OncePerImplementation | None = None,
                          evaluator: ConstantEvaluator | None = None,
                          project=None) -> _ProducerInputTable | None:
    """The producer-input -> consumer-argument overrides this argument documents.

    An empty mapping when there is no table. `None` when the table is there but wrong, so
    the caller stops rather than reporting a pile of consequential name mismatches. The rows
    are checked by one grammar, on a `JUST_INFO` table too, which Ford publishes. Each refusal
    below is at the row's line, once per implementation, one error per cell:

    - a row names an input of the producer (N1): one naming an output or `ierr` would be ignored
      by the plans while synthesis keys the call by it (rule K);
    - one row per input: of two, the last won silently;
    - `OMITTED` (rule O) is the keyword only bare and in capitals, as `AUTO` is; it leaves out
      only an input the producer may be called without -- an optional one, not one it requires
      in a mode chosen at run time -- and it is refused where the consumer has an argument, or
      the evaluator a constant, of that name, which the cell would otherwise mean.
    """
    once = once or _OncePerImplementation(diagnostics, set(), id(consumer))
    tables = [t for t in argument.doc.tables
              if _looks_like_a_producer_input_table(t, conventions)]
    mapping = _ProducerInputTable()
    if not tables:
        return mapping
    if len(tables) > 1:
        diagnostics.error(
            f"'{argument.name}' has {len(tables)} tables mapping the inputs of "
            f"'{producer.name}'",
            entity=argument,
            note="the mapping goes in exactly one table",
        )
        return None

    implementation = _implementation_name(consumer, project, conventions)
    table_of = f"the producer-input table of '{argument.name}'"
    rows: dict[str, int] = {}
    refused = False
    for row_index, (producer_input, supplied_by) in enumerate(tables[0].rows, start=1):
        wanted = producer_input.text.strip()
        given = supplied_by.text.strip()
        location = _location_of(argument, producer_input.line_number)

        def refuse(message: str, note: str) -> None:
            once.error(message, entity=argument, location=location, note=note)

        if not wanted or not given:
            diagnostics.error(
                f"row {row_index} of the producer-input table of '{argument.name}' "
                f"has an empty cell",
                entity=argument,
                note=_producer_input_table_example(conventions),
            )
            return None
        dummy = producer.argument(wanted)
        if dummy is None:
            diagnostics.error(
                f"the producer-input table of '{argument.name}' maps '{wanted}', which "
                f"is not an argument of '{producer.name}'",
                entity=argument,
            )
            return None
        if (dummy.intent is None or not dummy.intent.is_input
                or dummy.name.lower() == conventions.error_arg.lower()):
            refuse(f"{table_of} maps '{wanted}', which is an output of '{producer.name}', not "
                   f"an input",
                   "a row supplies an input the producer reads; DM_OUTPUT_FROM names the output "
                   "it takes")
            refused = True
            continue
        first = rows.setdefault(wanted.lower(), row_index)
        if first != row_index:
            refuse(f"{table_of} has two rows for '{wanted}' (rows {first} and {row_index})",
                   "a producer input is supplied once: keep the row you mean")
            refused = True
            continue
        if given != OMITTED and given.strip("'\"`").strip().lower() == OMITTED.lower():
            refuse(f"{table_of} supplies '{wanted}' by '{given}', which is not the keyword "
                   f"{OMITTED}",
                   f"the keyword is written bare and in capitals, as AUTO is: {OMITTED}")
            refused = True
            continue
        if given == OMITTED:
            why = _refuses_to_omit(dummy, producer, consumer, evaluator, implementation)
            if why is not None:
                refuse(f"{table_of} {why[0]}", why[1])
                refused = True
                continue
        mapping[wanted.lower()] = given
        mapping.lines[wanted.lower()] = producer_input.line_number
    return None if refused else mapping


def _refuses_to_omit(dummy: Argument, producer: Procedure, consumer: Procedure,
                     evaluator: ConstantEvaluator | None,
                     implementation: str) -> tuple[str, str] | None:
    """Why an `OMITTED` row may not leave the producer input `dummy` out, as (message, note)."""
    name = dummy.name
    if consumer.argument(OMITTED) is not None:
        return (f"supplies '{name}' by '{OMITTED}', which is the keyword, but "
                f"'{consumer.argument(OMITTED).name}' is also an argument of '{implementation}'",
                f"rename the argument, so the cell means one thing")
    if evaluator is not None:
        try:
            evaluator.evaluate(OMITTED)
        except ConstantError:
            pass
        else:
            return (f"supplies '{name}' by '{OMITTED}', which is the keyword, but "
                    f"'{OMITTED.lower()}' is also a constant",
                    "rename the constant, so the cell means one thing")
    if not dummy.optional:
        return (f"omits '{name}', which '{producer.name}' requires",
                f"only an optional input may be left out; supply '{name}' by an argument the "
                "caller passes, or a constant")
    required = dummy.directives.required_if_mode
    # a mode chosen at run time: a split producer's per-mode wrapper carries no such directive --
    # its input is required there (refused above) or not an argument at all
    if required is not None:
        selector = producer.argument(required.mode_arg)
        table = selector.roles.mode if selector is not None and selector.roles else None
        mode = required.mode_param
        value = next((v for v in table.values if v.parameter.lower() == mode.lower()),
                     None) if table is not None else None
        return (f"omits '{name}', which '{producer.name}' requires in mode "
                f"'{value.string if value is not None else mode}'",
                f"the mode is chosen at run time, so the input may be needed: supply '{name}'")
    return None


def _implementation_name(consumer: Procedure, project, conventions: Conventions) -> str:
    """What a diagnostic calls the consumer: the implementation the author wrote, where the
    consumer is a generated wrapper (4c), rather than a wrapper the author never sees."""
    implementation = (_implementation_of(consumer, project, conventions)
                      if project is not None else None)
    return implementation.name if implementation is not None else consumer.name


def _split_mode_of(consumer: Procedure, project, conventions: Conventions) -> str | None:
    """The mode argument the implementation of `consumer` splits into per-mode wrappers, where
    `consumer` is one of them -- which has no such argument, its mode fixed there."""
    implementation = _implementation_of(consumer, project, conventions)
    if implementation is None:
        return None
    for dummy in implementation.arguments:
        mode = dummy.roles.mode if dummy.roles is not None else None
        if mode is not None and mode.is_split and consumer.argument(dummy.name) is None:
            return dummy.name
    return None


def _refuses_what_nothing_supplies(argument: Argument, directive, consumer: Procedure,
                                   producer: Procedure, split: str | None, split_supplied,
                                   unresolved, unmatched, renames, fortran_calls: bool,
                                   project, conventions: Conventions,
                                   once: _OncePerImplementation) -> bool:
    """Refuse a producer input that resolves to nothing, each cause with its own message.

    Once per implementation, one error per input or cell:

    - a mode the implementation splits (T, Q1 -- 4b): each per-mode wrapper fixes it, and
      it is no argument there;
    - a `Supplied by` cell that is neither an argument nor a constant (4a), at its row;
    - an input with no row and no argument of its name, at the directive's line; where it is
      optional, the note offers the `OMITTED` row.

    The messages name the implementation, not the generated wrapper (4c). True when something
    was refused.
    """
    from .validate import _nearly  # local: validate is imported here at the top, by name

    location = _location_of(argument, directive.line_number)
    implementation = _implementation_name(consumer, project, conventions)
    computed = f"'{argument.name}' is computed from '{producer.name}'"
    for name in split_supplied:
        once.error(
            f"{computed}, whose input '{name}' would be '{split}', the mode "
            f"'{implementation}' splits into one wrapper per mode",
            entity=argument, location=location,
            note=(f"each per-mode wrapper fixes the mode, and passes it to no producer: supply "
                  f"'{name}' by an argument the caller passes, or a constant in the table"
                  + _way_out(fortran_calls)),
        )
    arguments = {dummy.name.lower() for dummy in consumer.arguments}
    for name, cell in unresolved:
        near = _nearly(cell.strip().lower(), arguments)
        once.error(
            f"the producer-input table of '{argument.name}' supplies '{name}' by '{cell}', "
            f"which is neither an argument of '{implementation}' nor a constant",
            entity=argument,
            location=_location_of(argument, renames.lines.get(name.lower())),
            note=("a cell is an argument the caller passes, or a constant the generator can "
                  "evaluate (1_int32, .false., a named parameter); an expression is not "
                  "computed: give the implementation an argument the caller passes for the value"
                  + (f"; did you mean '{near}'?" if near is not None else "")),
        )
    if unmatched:
        optional = [name for name in sorted(unmatched) if producer.argument(name).optional]
        offer = "".join(
            f"'{name}' is optional, so the row `| {name} | {OMITTED} |` leaves it out of the "
            "call; " for name in optional)
        once.error(
            f"{computed}, but its input(s) {', '.join(sorted(unmatched))} have no argument of "
            f"the same name in '{implementation}'",
            entity=argument, location=location,
            note=offer + _producer_input_table_example(conventions),
        )
    return bool(split_supplied or unresolved or unmatched)


def _producer_input_table_example(conventions: Conventions) -> str:
    left = conventions.producer_input_header.title()
    right = conventions.producer_supplied_by_header.capitalize()
    return (
        "name-matching is tried first; where the names differ, or where the consumer has\n"
        "no argument for an input at all, map it in a table -- the right-hand cell is an\n"
        "argument of the implementation that the caller passes, a Fortran constant, or\n"
        f"{OMITTED}, which leaves an optional input out of the call:\n"
        f"!! | {left} | {right} |\n"
        f"!! |----------------|-------------|\n"
        f"!! | n_elements     | n_genes     |\n"
        f"!! | n_dim          | 1_int32     |\n"
        f"!! | padding        | {OMITTED}     |"
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
                is_omp_threads=self.conventions.is_omp_threads_arg(argument.name),
                is_bitmask=self.conventions.is_bitmask_name(argument.name),
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
            self._warn_unread_mode_table(argument)
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

    def _warn_unread_mode_table(self, argument: Argument) -> None:
        """Warn where a table of accepted values sits on an argument that is not a mode.

        Mode detection goes by the argument's *name*, so a table like this is never read: the
        binding passes the raw integer instead of the mode string, and a `Procedure` column
        splits nothing. Both are silent -- the same silence that makes an unexpanded `DM_`
        macro an error -- but a table with a `Value` column need not be a mode table at all,
        so this stays a warning rather than an error.
        """
        for table in argument.doc.tables:
            if table.n_columns not in (2, 3):
                continue
            header = [text.strip().lower() for text in table.header_text]
            if header[1] != "value" and header[0] not in self.conventions.mode_aliases:
                continue
            aliases = " / ".join(f"<name>_{a}" for a in self.conventions.mode_aliases)
            self.diagnostics.warn(
                f"'{argument.name}' documents a table of values but is not a mode argument",
                entity=argument,
                location=_location_of(argument, table.line_number),
                note=(
                    f"a mode argument is named for what it selects ({aliases}, or the alias "
                    "alone); named otherwise the table is documentation only, and the value "
                    "crosses the binding as an integer rather than its name"
                ),
            )
            return

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
