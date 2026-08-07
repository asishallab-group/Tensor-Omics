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
    if module.name.lower().endswith(conventions.impl_suffix):
        validate_impl_module(module, diagnostics, conventions)
    for procedure in module.exported_procedures:
        validate_procedure(procedure, diagnostics, conventions)


def validate_impl_module(module: Module, diagnostics: DiagnosticBag,
                         conventions: Conventions = CONVENTIONS) -> None:
    """Rules a hand-written implementation module must satisfy.

    These hold for the *unexported* procedures the rules above skip, because an
    implementation is read for what it makes the generator write rather than wrapped itself.
    They are what keeps an implementation module a module of implementations: every entry
    point published through a generated wrapper, every allocation owned by that wrapper.

    The module name is the whole trigger, so these rules follow the `_impl` suffix wherever
    it is written -- under `src/tox`, under `src/f42`, anywhere. That is the price of
    dropping the fixed directory the rules used to be scoped to, and it makes `_impl` a
    reserved suffix for the entire source tree.
    """
    _check_module_is_named_for_its_file(module, diagnostics)
    _check_impl_imports(module, diagnostics, conventions)
    for procedure in module.procedures:
        _check_impl_allocates(procedure, diagnostics)
        if not procedure.name.lower().endswith(conventions.impl_suffix):
            continue  # a recommend routine or a private helper, not an implementation
        _check_impl_is_not_exported(procedure, diagnostics, conventions)
        _check_impl_is_not_named_for_a_wrapper(procedure, diagnostics, conventions)
        _check_prologue(procedure, diagnostics, conventions)


def _check_module_is_named_for_its_file(module: Module, diagnostics: DiagnosticBag) -> None:
    """An implementation module lives in a file named after it.

    Two independent things derive from that name and they have to agree. Synthesis triggers
    on the *module* name and writes `tox_x.F90`; the cleaner and the Ford exclusion scan
    *file* names, without parsing, to know what the generator owns. Let the two diverge and
    the generator writes a wrapper that nothing cleans and nothing excludes -- so the next
    run parses its own output and defines the module twice.

    Skipped where there is no file to check: a hand-built project in a unit test.
    """
    path = module.location.file
    if path is None or path.stem.lower() == module.name.lower():
        return
    diagnostics.error(
        f"implementation module '{module.name}' is in '{path.name}'",
        entity=module,
        note=(
            "the generated wrapper's path is derived from the file name and its contents "
            f"from the module name -- name the file '{module.name}.F90' so the two agree"
        ),
    )


def _check_impl_allocates(procedure: Procedure, diagnostics: DiagnosticBag) -> None:
    """An implementation module allocates nothing: its generated wrapper owns the memory.

    An implementation that allocates for itself hides the allocation from the caller who is
    meant to be able to avoid it -- the expert tier exists precisely so a caller can hand in
    reused buffers -- and it hides `ERR_ALLOC_FAIL` behind a name that promises no
    allocation. The rule covers every procedure in the module, not only the implementations:
    one that allocates nothing itself but calls a helper that does is no better off. It holds
    across modules too, which is `_check_impl_imports`' doing: nothing else could see it.

    A local and an `allocatable` dummy count alike. The dummy is the subtler of the two --
    it looks like the caller's memory, and the caller does receive it, but whoever fills it
    is the one who allocated it, and a `tmp_` argument already expresses that without an
    allocatable in the signature (which the C layer could not carry anyway).

    What to do instead: declare the buffer as a `tmp_` dummy. It disappears from the
    allocating wrapper's signature and is allocated there, which is the whole convention.
    Where its size is not an expression over the other arguments, `DM_OUTPUT_FROM(..., AUTO)`
    names the routine that computes it.
    """
    dummies = [a.name for a in procedure.arguments if "allocatable" in a.attributes]
    offenders = list(procedure.allocatable_locals) + dummies
    if not offenders:
        return
    names = ", ".join(f"'{name}'" for name in offenders)
    diagnostics.error(
        f"'{procedure.name}' allocates: {names}",
        entity=procedure,
        note=(
            "an implementation module owns no memory -- pass the buffer as a 'tmp_' dummy "
            "and the generated wrapper allocates it; size it with "
            "'DM_OUTPUT_FROM(..., AUTO)' where no expression over the arguments will do"
        ),
    )


def _check_impl_imports(module: Module, diagnostics: DiagnosticBag,
                        conventions: Conventions) -> None:
    """An implementation reaches only other implementations and the listed infrastructure.

    This is what makes `_check_impl_allocates` mean anything beyond one file. That check
    reads declarations, so it sees a procedure allocating for itself and nothing else; an
    implementation calling a helper in some other module that allocates would pass it. Bound
    the imports and the property becomes checkable: another `_impl` module is held to the
    same rule, and `Conventions.impl_import_whitelist` is the curated rest.

    It fixes the direction too. Nothing else stops an implementation from `use`ing a
    *generated* wrapper -- a layering inversion, and within one family a module cycle that
    surfaces as a build error with no hint of what caused it.

    Procedure-level imports count. Fortran lets a `use` sit inside a procedure, and a rule
    about what a module may reach that only reads the module header is one indentation level
    away from being bypassed.
    """
    allowed = {name.lower() for name in conventions.impl_import_whitelist}
    suffix = conventions.impl_suffix
    for entity, used_by in [(module, module.uses)] + [
        (procedure, procedure.uses) for procedure in module.procedures
    ]:
        for used in used_by:
            lowered = used.lower()
            if lowered in allowed or lowered.endswith(suffix):
                continue
            diagnostics.error(
                f"implementation module uses '{used}'",
                entity=entity,
                note=(
                    "an implementation may use another implementation module or one of "
                    f"{', '.join(repr(n) for n in conventions.impl_import_whitelist)} -- "
                    "that bound is what keeps it allocation-free and keeps it from reaching "
                    "a generated wrapper. Add it to 'impl_import_whitelist' if it really is "
                    "infrastructure, or take what you need as an argument"
                ),
            )


def _check_impl_is_not_exported(procedure: Procedure, diagnostics: DiagnosticBag,
                                conventions: Conventions) -> None:
    """An implementation is reached through its wrapper, never directly.

    Exporting it publishes a second entry point that skips every check the wrapper exists
    to make, under a name (`foo_impl`) that a binding caller cannot tell apart from the
    validated `foo` beside it. Support procedures in the same module -- the recommend
    routines a `DM_OUTPUT_FROM` producer names -- are exported as usual; they are not
    implementations and have no wrapper.
    """
    if not procedure.is_exported:
        return
    diagnostics.error(
        f"implementation '{procedure.name}' is exported",
        entity=procedure,
        note=(
            f"drop the '{conventions.export_marker}': the generated wrapper is what the "
            "bindings call, and exporting the implementation beside it publishes an "
            "unvalidated twin of the same procedure"
        ),
    )


#: The suffixes an implementation's base name may not carry, and what each would collide
#: with. `_expert` is the generator's own: `foo_expert_impl` beside `foo_impl` makes two
#: procedures called `foo_expert`, and the emitter would strip the suffix and call
#: `foo_impl` from the wrong one -- wrong code that compiles, because `foo_impl` exists.
#: `_alloc` is the hand-written pair convention `abi.c_abi.stripped_name` still reads, so
#: `foo_alloc_impl` would generate `foo_alloc` and publish it to Python and R as `foo`,
#: colliding with a real `foo` in the same family.
_RESERVED_BASE_SUFFIXES = ("expert_suffix", "alloc_suffix")


def _check_impl_is_not_named_for_a_wrapper(procedure: Procedure, diagnostics: DiagnosticBag,
                                           conventions: Conventions) -> None:
    """An implementation may not be named for one of the wrappers it generates.

    Name it for what it computes, and let its `tmp_` arguments decide whether a second tier
    is generated at all -- the suffix is never the author's to choose.
    """
    base = procedure.name.lower()[: -len(conventions.impl_suffix)]
    for field in _RESERVED_BASE_SUFFIXES:
        suffix = getattr(conventions, field)
        if not base.endswith(suffix):
            continue
        diagnostics.error(
            f"implementation '{procedure.name}' is named for a wrapper",
            entity=procedure,
            note=(
                f"'{suffix}' is the generator's to add: name the implementation for what it "
                "computes, and the wrappers are generated from its "
                f"'{conventions.temporary_prefix}' arguments"
            ),
        )
        return


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
                        "generated wrapper own the allocation"
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
        """A default and a required-in mode contradict, unless the implementation splits per mode.

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


def _check_prologue(impl: Procedure, diagnostics: DiagnosticBag,
                    conventions: Conventions) -> None:
    """A `DM_PROLOGUE` must name something the wrapper can actually call.

    Nothing checked any of this before, and every way of getting it wrong failed silently:
    the emitter is the only consumer of the directive, and an emitter has no line to point
    at. A prologue that named nothing produced a wrapper with no prologue; a dummy that
    matched nothing was dropped from the keyword call; a prologue with no `handled` left the
    wrapper branching on an undefined logical, which compiles.
    """
    directive = impl.directives.prologue
    if directive is None:
        return
    project = impl.module.project if impl.module is not None else None
    if project is None:
        return  # a hand-built procedure in a unit test; the real pipeline always has one

    prologue = project.procedure(directive.module, directive.procedure)
    if prologue is None:
        diagnostics.error(
            f"prologue '{directive.module}:{directive.procedure}' does not exist",
            entity=impl,
            note=(
                "the wrapper would be generated without a prologue at all, which is why "
                "this is an error rather than a warning"
            ),
        )
        return

    _check_prologue_reports_handled(impl, prologue, diagnostics, conventions)
    _check_prologue_runs_somewhere(impl, prologue, diagnostics, conventions)
    _check_prologue_arguments_resolve(impl, prologue, diagnostics, conventions)
    _check_prologue_writes_what_it_may(impl, prologue, diagnostics, conventions)
    _check_prologue_outputs_are_not_read_first(impl, prologue, diagnostics, conventions)


def _check_prologue_runs_somewhere(impl: Procedure, prologue: Procedure,
                                   diagnostics: DiagnosticBag,
                                   conventions: Conventions) -> None:
    """A prologue on an implementation with no allocating wrapper never runs.

    An allocating wrapper exists where the implementation has something to take over -- a work array,
    a permutation, a recommend-sized value -- or where the prologue asks for an argument of
    its own. A prologue with neither is attached to a procedure that is never generated: no
    call is emitted anywhere, and nothing said so.
    """
    from ..synthesize import prologue_only_arguments, taken_over_arguments

    if taken_over_arguments(impl.arguments, conventions, prologue):
        return
    if prologue_only_arguments(prologue, impl.arguments, conventions):
        return  # the wrapper exists to carry those, whether or not anything is allocated
    diagnostics.error(
        f"prologue '{prologue.name}' is on an implementation that generates no allocating wrapper",
        entity=impl,
        note=(
            "a prologue is what the allocating wrapper prepares beyond allocating, and that "
            "wrapper appears only where the implementation has work arrays to take over. This one "
            "has none, so the prologue would never be called -- if the work belongs to every "
            "caller of the implementation, put it at the top of the implementation instead"
        ),
    )


def _check_prologue_reports_handled(impl: Procedure, prologue: Procedure,
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
            entity=impl,
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


def _check_prologue_arguments_resolve(impl: Procedure, prologue: Procedure,
                                      diagnostics: DiagnosticBag,
                                      conventions: Conventions) -> None:
    """Every prologue dummy is supplied by name, so every name has to mean something.

    A name the implementation has is passed straight on. A name it does not becomes an argument of
    the allocating wrapper -- what the prologue derives from, which is the allocating tier's
    own vocabulary and no business of the implementation's: a threshold's `percentile`.

    Which leaves a misspelling nowhere to be caught, since it now reads as a new argument.
    So a name that is *nearly* one the implementation has is refused: `n_gene` beside `n_genes` is a
    typo, and silently turning it into something the caller must pass would mean the prologue
    and the implementation working from different numbers.

    There is no rename table, deliberately: a producer needs one because it is a published
    routine whose own parameter names cannot move, whereas a prologue is internal to the
    impl module and its dummy can simply be renamed to match.
    """
    supplied = {argument.name.lower() for argument in impl.arguments}
    reserved = {conventions.error_arg.lower(), conventions.prologue_handled_arg.lower()}
    dropped = _dropped_by_the_mode_split(impl)

    for dummy in prologue.arguments:
        name = dummy.name.lower()
        if name in reserved:
            continue
        if name in dropped:
            diagnostics.error(
                f"prologue argument '{dummy.name}' is not on every wrapper of "
                f"'{impl.name}'",
                dummy,
                note=(
                    "the implementation splits per mode, and this argument is the mode or belongs to "
                    "one of them, so the wrappers for the other modes do not have it -- take "
                    "something every mode has"
                ),
            )
            continue
        if name in supplied:
            continue
        near = _nearly(name, supplied - dropped)
        if near is not None:
            diagnostics.error(
                f"prologue argument '{dummy.name}' looks like a misspelling of '{near}'",
                dummy,
                note=(
                    f"a name '{impl.name}' does not have becomes an argument of the "
                    f"allocating wrapper, which is right for something the prologue derives "
                    f"from -- but '{dummy.name}' is one edit from '{near}', and the two would "
                    f"then be different values. Rename it either way"
                ),
            )


def _nearly(name: str, candidates) -> str | None:
    """A candidate one edit away from `name`, or None.

    One edit only: two names that differ by more than that are two names. Cheap enough to
    run over a handful of arguments, and it is the last thing standing between a typo and a
    silently invented argument.
    """
    for candidate in sorted(candidates):
        if abs(len(candidate) - len(name)) > 1:
            continue
        if _within_one_edit(name, candidate):
            return candidate
    return None


def _within_one_edit(a: str, b: str) -> bool:
    if a == b:
        return True
    if len(a) > len(b):
        a, b = b, a
    if len(b) - len(a) > 1:
        return False
    for index, (x, y) in enumerate(zip(a, b)):
        if x != y:
            if len(a) == len(b):  # a substitution
                return a[index + 1:] == b[index + 1:]
            return a[index:] == b[index + 1:]  # an insertion in the longer one
    return True


def _dropped_by_the_mode_split(impl: Procedure) -> set[str]:
    """The implementation's arguments that some generated wrapper will not have.

    A prologue is declared once and runs in every wrapper, so it may only name what every
    wrapper has. That is all of them unless the implementation splits per mode, where the mode
    argument is fixed in the call and dropped, and a `DM_REQUIRED_IF_MODE` argument belongs
    to its own mode and is absent from the rest (`synthesize._mode_kept_arguments`).

    Checking against the implementation alone was wrong in both directions: it accepted the mode
    argument, which no wrapper has, and rejected every `DM_REQUIRED_IF_MODE` argument, which
    on an unsplit mode every wrapper does have -- there the directive only says the argument
    is required at run time, and the wrapper still passes it on.
    """
    mode = next(
        (
            argument
            for argument in impl.arguments
            if argument.roles is not None
            and argument.roles.mode is not None
            and argument.roles.mode.is_split
        ),
        None,
    )
    if mode is None:
        return set()
    dropped = {mode.name.lower()}
    for argument in impl.arguments:
        required = argument.directives.required_if_mode
        if required is not None and required.mode_arg.lower() == mode.name.lower():
            dropped.add(argument.name.lower())
    return dropped


def _check_prologue_writes_what_it_may(impl: Procedure, prologue: Procedure,
                                       diagnostics: DiagnosticBag,
                                       conventions: Conventions) -> None:
    """A prologue writing what the implementation reads makes it a local -- unless it cannot be one.

    Prologue `intent(out)` over impl `intent(in)` means the value never leaves the wrapper,
    so it is dropped from the allocating signature and declared as a local
    (`synthesize.is_prologue_produced`). That fails for the one argument that cannot be
    dropped: one that also gives an extent of something the caller still passes or receives,
    which has to stay in the signature or nothing could size what comes back. It is then an
    `intent(in)` dummy the wrapper would have to hand to something that writes it, which
    gfortran rejects -- in code the author never wrote.

    The other pairings are meaningful and left alone. Impl `intent(out)` means the two are
    alternative producers of one output: the prologue's value is what the caller gets on the
    `handled` path and the implementation overwrites it otherwise, exactly what `loess_degenerate_fit`
    did with `fitted_values`. Impl `intent(inout)` means the caller supplies a value the
    prologue then refines.
    """
    from ..synthesize import taken_over_arguments  # local: synthesize imports this module

    taken = {a.name.lower() for a in taken_over_arguments(impl.arguments, conventions,
                                                          prologue)}
    for dummy in prologue.arguments:
        if dummy.intent is None or not dummy.intent.is_output:
            continue
        argument = impl.argument(dummy.name)
        if argument is None or argument.intent is not Intent.IN:
            continue
        if argument.name.lower() in taken:
            continue  # dropped from the signature and made a local, which is the point
        diagnostics.error(
            f"prologue writes '{dummy.name}', which '{impl.name}' reads and the wrapper "
            f"cannot make a local",
            argument,
            note=(
                "it would become a local, but it also sizes something the caller still "
                "passes, so it has to stay a dummy -- and Fortran will not hand an "
                "intent(in) dummy to something that writes it. Size that argument from "
                "something the caller passes instead"
            ),
        )


def _check_prologue_outputs_are_not_read_first(impl: Procedure, prologue: Procedure,
                                               diagnostics: DiagnosticBag,
                                               conventions: Conventions) -> None:
    """A prologue may not produce something the setup above it already read.

    A prologue runs below the work arrays it exists to prepare, so anything it writes is
    undefined while the recommend calls, the allocations and the permutation sorts run --
    and a name resolves the same either way, so reading one early compiles and computes
    rubbish rather than failing.
    """
    # local imports: synthesize imports this module
    from ..synthesize import sorted_permutations, taken_over_arguments

    written = {
        dummy.name.lower()
        for dummy in prologue.arguments
        if dummy.intent is not None and dummy.intent.is_output
    }
    written -= {conventions.error_arg.lower(), conventions.prologue_handled_arg.lower()}
    if not written:
        return

    taken = taken_over_arguments(impl.arguments, conventions, prologue)

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
    # that describes something else. Nothing crashes; the implementation is handed a wrong answer.
    # A permutation the prologue builds itself is not sorted here, so it reads nothing.
    for argument in sorted_permutations(taken, prologue, conventions):
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
    # impl without a project, so it can never look a producer up, and an implementation argument's
    # `computed_from` is always None. Reading it here would have made this half dead code.
    for argument in _wrapper_arguments_with_plans(impl, conventions):
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


def _wrapper_arguments_with_plans(impl: Procedure, conventions: Conventions):
    """The expert wrapper's arguments that a producer fills, if one is built.

    `DM_OUTPUT_FROM` is only resolved into an `OutputFromPlan` by the project-wide pass, and
    that runs on the wrappers, not on the implementations. So the plans live on the wrapper
    that still *has* the recommend-sized arguments, which is the expert tier: the
    implementation `x_impl` in module `tox_y_impl` becomes `x_expert` in `tox_y`, and only
    `x` where there is no second tier.

    Asking for `x` first would find the allocating wrapper -- which has had exactly these
    arguments taken over, so the plan list comes back empty and this whole rule quietly
    stops firing. Absent in a unit test that validates an implementation module on its own.
    """
    module = impl.module
    project = module.project if module is not None else None
    if project is None:
        return []
    suffix = conventions.impl_suffix
    generated, base = module.name[: -len(suffix)], impl.name[: -len(suffix)]
    foo = project.procedure(generated, base + conventions.expert_suffix)
    if foo is None:
        foo = project.procedure(generated, base)
    if foo is None:
        return []
    return [
        argument
        for argument in foo.arguments
        if argument.roles is not None and argument.roles.computed_from is not None
    ]


#: Identifiers named in an extent expression
_IDENTIFIER_RE = re.compile(r"[A-Za-z_]\w*")

_READ_FIRST_NOTE = (
    "a prologue runs below the work arrays it prepares; take this value as a 'tmp_' argument "
    "instead, or derive it from something the caller passes"
)

_NON_OPTIONAL_NOTE = (
    "the C wrapper reads it before it may take c_loc of the arrays it describes, so it "
    "has to be there; see the null validation order in the generator README"
)
