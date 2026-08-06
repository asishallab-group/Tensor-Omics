"""Rules an exported procedure must satisfy.

Only exported procedures are checked. An internal routine is held to none of this,
because nothing is generated from it -- it may use deferred-length strings, derived types
and untyped intents as freely as Fortran allows.

Each rule exists because generating correct code is otherwise impossible, and each says
what to do instead. A rule that merely encodes taste belongs in a review, not here.
"""

from __future__ import annotations

import re

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
    if module.name.lower().endswith(conventions.kernel_suffix):
        validate_kernel_module(module, diagnostics, conventions)
    for procedure in module.exported_procedures:
        validate_procedure(procedure, diagnostics, conventions)


def validate_kernel_module(module: Module, diagnostics: DiagnosticBag,
                           conventions: Conventions = CONVENTIONS) -> None:
    """Rules a hand-written kernel module must satisfy.

    These hold for the *unexported* procedures the rules above skip, because a kernel is
    read for what it makes the generator write rather than wrapped itself. They are what
    keeps the kernel tree a tree of kernels: every entry point published through a
    generated wrapper, every allocation owned by the generated `_alloc`.
    """
    for procedure in module.procedures:
        _check_kernel_allocates(procedure, diagnostics)
        if not procedure.name.lower().endswith(conventions.kernel_suffix):
            continue  # a recommend routine or a private helper, not a kernel
        _check_kernel_is_not_exported(procedure, diagnostics, conventions)
        _check_kernel_is_not_named_alloc(procedure, diagnostics, conventions)
        _check_prologue(procedure, diagnostics, conventions)


def _check_kernel_allocates(procedure: Procedure, diagnostics: DiagnosticBag) -> None:
    """A kernel module allocates nothing: the generated `_alloc` wrapper owns the memory.

    A kernel that allocates for itself hides the allocation from the caller who is meant to
    be able to avoid it -- the expert tier exists precisely so a caller can hand in reused
    buffers -- and it hides `ERR_ALLOC_FAIL` behind a wrapper whose name promises no
    allocation. The rule covers every procedure in the module, not only the kernels: a
    kernel that allocates nothing itself but calls a helper that does is no better off.

    What to do instead: declare the buffer as a `tmp_` dummy. It disappears from the
    allocating wrapper's signature and is allocated there, which is the whole convention.
    Where its size is not an expression over the other arguments, `DM_OUTPUT_FROM(..., AUTO)`
    names the routine that computes it.
    """
    if not procedure.allocatable_locals:
        return
    names = ", ".join(f"'{name}'" for name in procedure.allocatable_locals)
    diagnostics.error(
        f"'{procedure.name}' allocates: {names}",
        entity=procedure,
        note=(
            "a kernel module owns no memory -- pass the buffer as a 'tmp_' dummy and the "
            "generated '_alloc' wrapper allocates it; size it with "
            "'DM_OUTPUT_FROM(..., AUTO)' where no expression over the arguments will do"
        ),
    )


def _check_kernel_is_not_exported(procedure: Procedure, diagnostics: DiagnosticBag,
                                  conventions: Conventions) -> None:
    """A kernel is reached through its wrapper, never directly.

    Exporting it publishes a second entry point that skips every check the wrapper exists
    to make, under a name (`foo_kernel`) that a binding caller cannot tell apart from the
    validated `foo` beside it. Support procedures in the same module -- the recommend
    routines a `DM_OUTPUT_FROM` producer names -- are exported as usual; they are not
    kernels and have no wrapper.
    """
    if not procedure.is_exported:
        return
    diagnostics.error(
        f"kernel '{procedure.name}' is exported",
        entity=procedure,
        note=(
            f"drop the '{conventions.export_marker}': the generated wrapper is what the "
            "bindings call, and exporting the kernel beside it publishes an unvalidated "
            "twin of the same procedure"
        ),
    )


def _check_kernel_is_not_named_alloc(procedure: Procedure, diagnostics: DiagnosticBag,
                                     conventions: Conventions) -> None:
    """`_alloc` is the generator's suffix, so a kernel may not claim it.

    A `foo_alloc_kernel` generates `foo_alloc` as its *validating* wrapper -- a procedure
    that allocates nothing wearing the name of the one that does, with no expert tier
    behind it. Name the kernel for what it computes and let its `tmp_` arguments decide
    whether an allocating wrapper is generated at all.
    """
    base = procedure.name.lower()[: -len(conventions.kernel_suffix)]
    if not base.endswith(conventions.alloc_suffix):
        return
    diagnostics.error(
        f"kernel '{procedure.name}' is named for the allocating wrapper",
        entity=procedure,
        note=(
            f"'{conventions.alloc_suffix}' is the generator's suffix: name the kernel "
            "for what it computes, and the allocating wrapper is generated from its "
            f"'{conventions.temporary_prefix}' arguments"
        ),
    )


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
        self._check_optional_output(argument)
        self._check_attributes(argument)
        self._check_type(argument)
        self._check_temporary(argument)
        self._check_mode(argument)
        self._check_required_if_mode(argument)
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

    def _check_optional_output(self, argument: Argument) -> None:
        """An optional output generates a binding that contradicts its own signature.

        The Fortran says the argument may be absent; every generated binding allocates it,
        passes it and returns it regardless. Nothing says so, so the declaration and the
        wrapper disagree and the reader has no way to tell which is true.

        Honouring it instead is what cannot be generated well: Python would return a dict
        whose keys vary per call and R a list of varying length, so a caller could not
        write against the signature at all. The shape that expresses the same intent is an
        optional *input* flag plus a `tmp_` work array -- what `loess_fit_plain` does with
        `compute_influence`: the work stays skippable and the return type stays fixed.

        Work arrays are exempt: they are dropped from the allocating wrapper and never
        reach a caller.
        """
        if argument.intent is not Intent.OUT or not argument.optional:
            return
        if argument.roles is not None and argument.roles.is_temporary:
            return
        self.error(
            f"argument '{argument.name}' is an optional output",
            argument,
            note=(
                "the binding languages would have to vary their return shape per call. "
                "Express it as an optional input flag plus a 'tmp_' work array, the way "
                "'compute_influence' does"
            ),
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

    def _check_required_if_mode(self, argument: Argument) -> None:
        """A default and a required-in mode contradict, unless the kernel splits per mode.

        Where the mode is resolved at runtime the argument is always passed on -- the binding
        supplies the default -- so "required in that mode" says nothing. Where the mode table
        names a procedure per value, the same pairing is meaningful: the directive scopes the
        argument to its mode (absent from the others) and the default applies within it.
        """
        required = argument.directives.required_if_mode
        if required is None or not argument.directives.has_default:
            return
        mode = self.procedure.argument(required.mode_arg)
        roles = mode.roles if mode is not None else None
        if roles is not None and roles.mode is not None and roles.mode.is_split:
            return
        self.error(
            f"argument '{argument.name}' has both a default and a mode it is required in",
            argument,
            note=(
                "DM_REQUIRED_IF_MODE is for optionals that have no default; an argument "
                "with a default is always passed on -- unless the mode table names a "
                "procedure per value, where the directive scopes the argument to its mode"
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


def _check_prologue(kernel: Procedure, diagnostics: DiagnosticBag,
                    conventions: Conventions) -> None:
    """A `DM_PROLOGUE` must name something the wrapper can actually call.

    Nothing checked any of this before, and every way of getting it wrong failed silently:
    the emitter is the only consumer of the directive, and an emitter has no line to point
    at. A prologue that named nothing produced a wrapper with no prologue; a dummy that
    matched nothing was dropped from the keyword call; a prologue with no `handled` left the
    wrapper branching on an undefined logical, which compiles.
    """
    directive = kernel.directives.prologue
    if directive is None:
        return
    project = kernel.module.project if kernel.module is not None else None
    if project is None:
        return  # a hand-built procedure in a unit test; the real pipeline always has one

    prologue = project.procedure(directive.module, directive.procedure)
    if prologue is None:
        diagnostics.error(
            f"prologue '{directive.module}:{directive.procedure}' does not exist",
            entity=kernel,
            note=(
                "the wrapper would be generated without a prologue at all, which is why "
                "this is an error rather than a warning"
            ),
        )
        return

    _check_prologue_reports_handled(kernel, prologue, diagnostics, conventions)
    _check_prologue_runs_somewhere(kernel, prologue, diagnostics, conventions)
    _check_prologue_arguments_resolve(kernel, prologue, diagnostics, conventions)
    _check_prologue_outputs_are_not_read_first(kernel, prologue, diagnostics, conventions)


def _check_prologue_runs_somewhere(kernel: Procedure, prologue: Procedure,
                                   diagnostics: DiagnosticBag,
                                   conventions: Conventions) -> None:
    """A prologue scoped to a wrapper the kernel does not generate never runs.

    An allocating wrapper exists only where the kernel has something to take over -- a work
    array, a permutation, a recommend-sized value. Scope a prologue to `ALLOC` on a kernel
    with none of those and it is attached to a procedure that is never generated: no call is
    emitted anywhere, and nothing said so.
    """
    from ..synthesize import taken_over_arguments

    directive = kernel.directives.prologue
    if directive.runs_in(False):
        return  # it runs in the validating wrapper, which every kernel has
    if taken_over_arguments(kernel.arguments, conventions):
        return
    diagnostics.error(
        f"prologue '{prologue.name}' is scoped to the allocating wrapper, which "
        f"'{kernel.name}' does not generate",
        entity=kernel,
        note=(
            f"an allocating wrapper appears only where the kernel has work arrays to take "
            f"over, and this one has none -- so the prologue would never be called"
        ),
    )


def _check_prologue_reports_handled(kernel: Procedure, prologue: Procedure,
                                    diagnostics: DiagnosticBag,
                                    conventions: Conventions) -> None:
    """The wrapper returns early on `handled`, so the prologue has to set it.

    The generated wrapper declares `handled` and branches on it unconditionally. A prologue
    that does not take it never assigns it, and the branch reads an undefined logical --
    which no compiler is obliged to reject and gfortran does not.
    """
    name = conventions.prologue_handled_arg
    handled = prologue.argument(name)
    if handled is None:
        diagnostics.error(
            f"prologue '{prologue.name}' has no '{name}' argument",
            entity=kernel,
            note=(
                f"every wrapper returns early on `if ({name})`, so a prologue that never "
                f"sets it leaves that branch reading an undefined value; declare "
                f"`logical, intent(out) :: {name}` and set it on every path"
            ),
        )
        return
    if handled.type.base is not BaseType.LOGICAL or handled.is_array:
        diagnostics.error(
            f"prologue argument '{name}' is not a scalar logical", handled,
            note="it is the flag the wrapper returns early on",
        )
    elif handled.intent is not Intent.OUT:
        diagnostics.error(
            f"prologue argument '{name}' is not intent(out)", handled,
            note="the prologue reports through it; the wrapper never supplies a value",
        )


def _check_prologue_arguments_resolve(kernel: Procedure, prologue: Procedure,
                                      diagnostics: DiagnosticBag,
                                      conventions: Conventions) -> None:
    """Every prologue dummy is supplied by name, so every name has to be there.

    The wrapper supplies `handled` and `ierr` itself and everything else from the kernel's
    own arguments -- either a dummy it passes on, or a work array it prepared. A name that
    is neither used to vanish from the generated call, which Fortran then rejected in code
    the author did not write.

    There is no rename table, deliberately: a producer needs one because it is a published
    routine whose own parameter names cannot move, whereas a prologue is internal to the
    kernel module and its dummy can simply be renamed to match.
    """
    supplied = {argument.name.lower() for argument in kernel.arguments}
    supplied |= {conventions.error_arg.lower(), conventions.prologue_handled_arg.lower()}
    supplied -= _dropped_by_the_mode_split(kernel)
    for dummy in prologue.arguments:
        if dummy.name.lower() in supplied:
            continue
        diagnostics.error(
            f"prologue argument '{dummy.name}' names nothing every wrapper of "
            f"'{kernel.name}' has",
            dummy,
            note=(
                f"a prologue's dummies are supplied by name from the wrapper's arguments; "
                f"rename this one to whatever '{kernel.name}' calls the same thing, or -- "
                f"if it is scoped to one mode -- take something every mode has"
            ),
        )


def _dropped_by_the_mode_split(kernel: Procedure) -> set[str]:
    """The kernel's arguments that some generated wrapper will not have.

    A prologue is declared once and runs in every wrapper, so it may only name what every
    wrapper has. That is all of them unless the kernel splits per mode, where the mode
    argument is fixed in the call and dropped, and a `DM_REQUIRED_IF_MODE` argument belongs
    to its own mode and is absent from the rest (`synthesize._mode_kept_arguments`).

    Checking against the kernel alone was wrong in both directions: it accepted the mode
    argument, which no wrapper has, and rejected every `DM_REQUIRED_IF_MODE` argument, which
    on an unsplit mode every wrapper does have -- there the directive only says the argument
    is required at run time, and the wrapper still passes it on.
    """
    mode = next(
        (
            argument
            for argument in kernel.arguments
            if argument.roles is not None
            and argument.roles.mode is not None
            and argument.roles.mode.is_split
        ),
        None,
    )
    if mode is None:
        return set()
    dropped = {mode.name.lower()}
    for argument in kernel.arguments:
        required = argument.directives.required_if_mode
        if required is not None and required.mode_arg.lower() == mode.name.lower():
            dropped.add(argument.name.lower())
    return dropped


def _check_prologue_outputs_are_not_read_first(kernel: Procedure, prologue: Procedure,
                                               diagnostics: DiagnosticBag,
                                               conventions: Conventions) -> None:
    """A late prologue may not produce something the setup above it already read.

    A prologue that takes one of the wrapper's work arrays runs below the allocations, since
    there would be nothing to hand it above them. Anything it writes is therefore undefined
    while those allocations and the recommend calls run -- and a name resolves the same
    either way, so reading one early compiles and computes rubbish.
    """
    # local imports: synthesize imports this module
    from ..synthesize import is_permutation, is_temporary, taken_over_arguments

    if not prologue_runs_late(kernel, prologue, conventions):
        return

    written = {
        dummy.name.lower()
        for dummy in prologue.arguments
        if dummy.intent is not None and dummy.intent.is_output
    }
    written -= {conventions.error_arg.lower(), conventions.prologue_handled_arg.lower()}
    if not written:
        return

    taken = taken_over_arguments(kernel.arguments, conventions)

    # `M_ALLOCATE(tmp_x(<extent>))`
    for argument in taken:
        for extent in argument.dimension.extents:
            for identifier in _IDENTIFIER_RE.findall(extent):
                if identifier.lower() not in written:
                    continue
                diagnostics.error(
                    f"'{argument.name}' is sized by '{identifier}', which the prologue "
                    f"only fills afterwards",
                    argument,
                    note=_READ_FIRST_NOTE,
                )

    # `call sort_array_heapsort(<base>, <base>_perm)` -- the permutation is built from the
    # data as it stands there, so a prologue that then rewrites the data leaves an order
    # that describes something else. Nothing crashes; the kernel is handed a wrong answer.
    for argument in taken:
        if not is_permutation(argument, conventions) or is_temporary(argument, conventions):
            continue
        base = argument.name[: -len(conventions.perm_suffix)]
        if base.lower() not in written:
            continue
        diagnostics.error(
            f"'{argument.name}' is sorted against '{base}' before the prologue rewrites it",
            argument,
            note=_READ_FIRST_NOTE,
        )

    # `call <recommend>(...)` -- its inputs are read where it is called, above the prologue.
    # The plan is resolved on the generated wrapper rather than here: `analyse` runs on a
    # kernel without a project, so it can never look a producer up, and a kernel argument's
    # `computed_from` is always None. Reading it here would have made this half dead code.
    for argument in _wrapper_arguments_with_plans(kernel, conventions):
        plan = argument.roles.computed_from
        for supply in plan.inputs:
            if supply.argument is None or supply.argument.lower() not in written:
                continue
            diagnostics.error(
                f"'{plan.producer.name}' is passed '{supply.argument}', which the prologue "
                f"only fills afterwards",
                argument,
                note=_READ_FIRST_NOTE,
            )


def _wrapper_arguments_with_plans(kernel: Procedure, conventions: Conventions):
    """The generated validating wrapper's arguments that a producer fills, if it is built.

    `DM_OUTPUT_FROM` is only resolved into an `OutputFromPlan` by the project-wide pass, and
    that runs on the wrappers, not on the kernels. So the plans live on `foo`, which is
    reachable from here: the kernel `x_kernel` in module `tox_y_kernel` becomes `x` in
    `tox_y`. Absent in a unit test that validates a kernel module on its own.
    """
    module = kernel.module
    project = module.project if module is not None else None
    if project is None:
        return []
    suffix = conventions.kernel_suffix
    foo = project.procedure(module.name[: -len(suffix)], kernel.name[: -len(suffix)])
    if foo is None:
        return []
    return [
        argument
        for argument in foo.arguments
        if argument.roles is not None and argument.roles.computed_from is not None
    ]


def prologue_runs_late(kernel: Procedure, prologue: Procedure,
                       conventions: Conventions = CONVENTIONS) -> bool:
    """Whether the allocating wrapper emits its prologue below the allocations.

    Mirrors `emit.fortran_wrapper._prologue_needs_locals`, on the kernel side: a prologue
    that takes a work array cannot precede the allocation of that work array.
    """
    from ..synthesize import taken_over_arguments

    directive = kernel.directives.prologue
    if directive is None or not directive.runs_in(True):
        return False
    taken = {a.name.lower() for a in taken_over_arguments(kernel.arguments, conventions)}
    return any(dummy.name.lower() in taken for dummy in prologue.arguments)


#: Identifiers named in an extent expression
_IDENTIFIER_RE = re.compile(r"[A-Za-z_]\w*")

_READ_FIRST_NOTE = (
    "this prologue takes a work array, so it runs below the allocations; give it the value "
    "as a 'tmp_' argument instead, or compute the size from something the caller passes"
)

_NON_OPTIONAL_NOTE = (
    "the C wrapper reads it before it may take c_loc of the arrays it describes, so it "
    "has to be there; see the null validation order in the generator README"
)
