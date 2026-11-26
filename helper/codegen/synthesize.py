"""Synthesise the wrappers an implementation implies.

An author writes only the implementation -- `loess_fit_impl` in module `tox_loess_impl`,
with no `ierr` and no input validation. From it the generator builds the entry point
`loess_fit`: validate the inputs, then call the implementation. Where the implementation
has work arrays to take over, that entry point does more -- allocate them, call the
recommend routines, build and sort the permutations, run the prologue -- and the plain
validate-and-call wrapper becomes `loess_fit_expert`, for a caller who wants to supply all
of that themselves. Both land in a generated module `tox_loess`, carry the export category,
and the rest of the pipeline wraps them to C/Python/R exactly as it does any procedure read
from source.

The wrappers are built here as IR and injected into the project *before* the semantic pass,
so a single parse of the implementation is the one source of truth: the Fortran wrappers and
their bindings both come from it, with no round-trip through Ford. The implementation
modules stay in the project but are inert -- nothing generates from an unexported procedure.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, replace
from pathlib import Path

from .config import CONVENTIONS, Conventions, Paths
from .diagnostics import DiagnosticBag
from .ir.doc import Doc
from .ir.entities import Argument, Meta, Module, Procedure, Project
from .ir.roles import analyse, plain_of_an_expert, producer_call_key
from .ir.types import BaseType, FortranType, Intent
from .ir.validate import function_result_problems

#: Documentation for the `ierr` a wrapper carries when the implementation declares none
_ERROR_DOC = "Error code; zero on success, non-zero on failure."


def generated_path_for(
    source: Path, paths: Paths, conventions: Conventions = CONVENTIONS
) -> Path | None:
    """Where the wrappers for the implementation file `source` are written, or None.

    The rule is a mirror: `src/<rest>/<module>_impl.F90` generates
    `src/generated/<rest>/<module>.F90`. Nothing in it knows about `tox` or `f42`, so an
    implementation anywhere under `src/` generates beside its own layer rather than into one
    directory the rule would have to name -- which is what lets f42 write implementations
    without the generator learning a second case.

    None where `source` is not an implementation at all, or lies outside the source tree
    (a hand-built test project, whose modules have no file). A file already under
    `generated_dir` is never claimed: that is the generator's own output.

    `source` may be either spelling. The IR carries paths relative to the root, so that a
    diagnostic can quote one; a caller holding an absolute path should not have to know
    that, so both are accepted and joined the one way `Paths` joins anything.
    """
    src = paths.resolve(paths.src_dir)
    out = paths.resolve(paths.generated_dir)
    suffix = conventions.impl_suffix
    if not source.stem.lower().endswith(suffix):
        return None
    source = paths.resolve(source)
    try:
        # both sides resolved: the frontend hands back paths from Ford's own resolution,
        # while `paths.resolve` only joins -- so on a checkout reached through a symlink (or
        # this repository's `/mnt/c` mount) the two spellings of the same directory would
        # not compare equal and every wrapper would look misplaced
        relative = source.resolve().relative_to(src.resolve())
    except (ValueError, OSError):
        return None
    if out.resolve() in source.resolve().parents:
        return None
    return out / relative.parent / f"{source.stem[: -len(suffix)]}.F90"


def generated_wrapper_paths(paths: Paths) -> list[Path]:
    """Every wrapper module on disk under `generated_dir`: what the fortran target cleans.

    Read off the disk, not off the implementations that exist now. Derived from those, the set
    missed the wrapper of an implementation that had been deleted, and fpm compiles all of
    `src/generated`: the orphan, still `use`-ing the module that went with it, broke every
    build after, regenerated or not (#221 B21). Reading the disk is safe by the rule `Paths`
    states -- nothing hand-written lives under `generated_dir` -- and the fortran target owns
    every `.F90` there outside the binding trees, which the c and r targets clean as theirs.
    So `--target fortran` leaves the C wrappers and the R shims alone.

    A missing source tree owns nothing, as it did when the set was derived from it: a mistyped
    `--src` must not take every wrapper with it.
    """
    src = paths.resolve(paths.src_dir)
    out = paths.resolve(paths.generated_dir)
    if not src.is_dir() or not out.is_dir():
        return []
    # resolved, as in `generated_path_for`: two spellings of one directory must compare equal
    bindings = {paths.resolve(d).resolve() for d in (paths.c_binding_dir, paths.r_binding_dir)}
    return [
        path
        for path in sorted(out.rglob("*.F90"))
        if bindings.isdisjoint(path.resolve().parents)
    ]


# -- what the wrappers publish ---------------------------------------------------


def published_result(impl: Procedure, conventions: Conventions = CONVENTIONS) -> Argument | None:
    """An implementation function's result, as its wrappers publish it, or None.

    Both tiers are subroutines, so the result becomes one of their `intent(out)` arguments,
    named and documented as the author declared it, and the wrapper assigns the call to it
    (`total = foo_impl(...)`). It is an ordinary argument from here on -- validated, bound and
    returned like any output -- which is why the copy no longer says it is a result: the
    wrapper has none. None for a subroutine, and for a result the wrappers cannot carry under
    its name (`validate.function_result_problems`, which says why).
    """
    result = impl.result
    if result is None or function_result_problems(impl, conventions):
        return None
    published = Argument(
        result.name,
        result.type,
        dimension=result.dimension,
        intent=Intent.OUT,
        doc=result.doc,
        directives=result.directives,
        attributes=result.attributes,
        location=result.location,
    )
    published.parent = impl
    return published


def published_arguments(impl: Procedure, conventions: Conventions = CONVENTIONS) -> list[Argument]:
    """The implementation's arguments as its wrappers start from: a function's result included.

    The result goes right before `ierr` (#221 B1, FES's decision): before the implementation's
    own `ierr` where it declares one, else last, which is where the wrapper's `ierr` is then
    appended -- so the published signature ends in `..., result, ierr` either way, as an
    exported function's C wrapper does (`abi.c_abi`).
    """
    arguments = list(impl.arguments)
    result = published_result(impl, conventions)
    if result is None:
        return arguments
    error = next(
        (index for index, argument in enumerate(arguments)
         if argument.name.lower() == conventions.error_arg.lower()),
        len(arguments),
    )
    return arguments[:error] + [result] + arguments[error:]


# -- what the allocating wrapper takes over from the caller ----------------------
#
# These are the arguments a caller should not have to supply: they are dropped from the
# allocating wrapper's signature and become locals it prepares itself. The predicates are
# name/directive based, so they hold both before the semantic pass (synthesis, to shape the
# signature) and after it (emission, to fill the body).


def is_temporary(argument: Argument, conventions: Conventions = CONVENTIONS) -> bool:
    """A `tmp_` work array: allocated silently, never returned."""
    return argument.name.lower().startswith(conventions.temporary_prefix)


def is_permutation(argument: Argument, conventions: Conventions = CONVENTIONS) -> bool:
    """A `<base>_perm` vector: seeded and sorted in the allocating wrapper."""
    return argument.name.lower().endswith(conventions.perm_suffix)


def is_computed(argument: Argument) -> bool:
    """Sized by a recommend routine (`DM_OUTPUT_FROM(..., AUTO)`)."""
    directive = argument.directives.output_from
    return directive is not None and directive.is_automatic


def is_taken_over(argument: Argument, conventions: Conventions = CONVENTIONS) -> bool:
    return (
        is_temporary(argument, conventions)
        or is_permutation(argument, conventions)
        or is_computed(argument)
    )


# -- what sizes what ------------------------------------------------------------------


@dataclass(frozen=True)
class SizeName:
    """An identifier an extent or a character length names, and how it names it."""

    name: str
    #: The array operand of `size()`: the declaration reads its shape there, not its value
    is_size_operand: bool = False


#: One token of an extent or a length: a number (its kind with it), a string, a name, a symbol
_SIZE_TOKEN_RE = re.compile(
    r"\s*(?:(?P<number>(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?(?:_\w+)?)"
    r"|(?P<string>'[^']*'|\"[^\"]*\")"
    r"|(?P<name>[A-Za-z]\w*)"
    r"|(?P<symbol>\S))"
)


def extent_tokens(text: str) -> list[tuple[str, str]]:
    """An extent or a length expression as `(kind, text)` tokens: `number` (its kind suffix
    with it), `string`, `name` or `symbol`. Blanks are dropped."""
    tokens = []
    position = 0
    while position < len(text):
        match = _SIZE_TOKEN_RE.match(text, position)
        if match is None:
            break  # only blanks are left
        position = match.end()
        tokens.append((match.lastgroup, match.group(match.lastgroup)))
    return tokens


def size_names_in(text: str) -> list[SizeName]:
    """The identifiers an extent or a length expression names, in order, each marked as a value
    or as the array operand of `size()`.

    The one answer to "what does this declaration size itself by", which every caller asks
    rather than scanning the text its own way: `_sizes_a_kept_argument` (where a `size()` operand
    counts, since the kept declaration names it), `binding_computes` (where it does not: the
    binding reads the shape of an array the caller passed) and `abi.c_abi.reads_a_dummy`. A
    number names nothing, its kind suffix included (`5_int32`), nor does a string, the `size`
    intrinsic itself, or a keyword (`dim=`).
    """
    tokens = extent_tokens(text)

    def token(index: int) -> tuple[str, str]:
        return tokens[index] if 0 <= index < len(tokens) else ("", "")

    found = []
    for index, (kind, value) in enumerate(tokens):
        if kind != "name":
            continue
        following = token(index + 1)[1]
        if value.lower() == "size" and following == "(":
            continue  # the intrinsic
        if following == "=" and token(index + 2)[1] != "=":
            continue  # a keyword, `size(a, dim=1)`
        operand = (token(index - 1)[1] == "(" and token(index - 2)[0] == "name"
                   and token(index - 2)[1].lower() == "size")
        found.append(SizeName(value, operand))
    return found


def size_names(argument: Argument) -> list[SizeName]:
    """`size_names_in` over every extent of `argument`, then its character length."""
    texts = list(argument.dimension.extents)
    length = argument.type.length if argument.type is not None else None
    if length is not None and length.expr is not None:
        texts.append(length.expr)
    return [named for text in texts for named in size_names_in(text)]


def _sizes_a_kept_argument(
    argument: Argument, arguments: Sequence[Argument], conventions: Conventions
) -> bool:
    """Whether an argument the caller still passes or receives is sized by `argument`.

    By any name in an extent or the length, a `size()` operand included: the kept argument's
    declaration names it, so the wrapper has to have it.
    """
    name = argument.name.lower()
    return any(
        named.name.lower() == name
        for other in arguments
        if other is not argument and not is_taken_over(other, conventions)
        for named in size_names(other)
    )


def is_prologue_produced(argument: Argument, prologue=None) -> bool:
    """Whether the prologue writes this argument and the implementation only reads it.

    Then nobody outside the wrapper is involved: the prologue produces the value, the
    implementation consumes it, and the caller neither supplies nor receives it. Both intents
    say so already, so this needs no naming convention on top -- and unlike a `tmp_` name, it
    lets the implementation keep the honest `intent(in)` for something it only reads.

    `foo` still takes it, because `foo` has no prologue: the expert tier is where a caller
    supplies the value themselves. Exactly the relation `<base>_perm` and the recommend-sized
    arguments already have with their two tiers.
    """
    if prologue is None or argument.intent is not Intent.IN:
        return False
    dummy = prologue.argument(argument.name)
    return dummy is not None and dummy.intent is Intent.OUT


def taken_over_arguments(
    arguments: Sequence[Argument],
    conventions: Conventions = CONVENTIONS,
    prologue=None,
    evaluator=None,
) -> list[Argument]:
    """The arguments the allocating wrapper prepares itself, rather than taking from the
    caller.

    Work arrays always qualify -- they are scratch. A `<base>_perm` qualifies only while the
    wrapper can actually build it, which means the `<base>` it orders is an argument too:
    `values_perm` is seeded and sorted against `values`, but a name that merely ends in `_perm`
    and orders something the implementation never receives is the caller's own data. A
    recommend-sized value qualifies as well, and so does one the prologue produces for the
    implementation to read -- both only while nothing the caller still sees depends on them: an
    argument that also gives an extent of a returned array has to stay in the signature, or
    neither the caller nor the binding could size what comes back.

    **A producer call is taken over whole or kept whole (rule K).** Where one value a call
    fills has to stay, every value that call fills stays, so the caller makes the one call and
    `foo` none. A call is keyed by `roles.producer_call_key`, the resolver the plans are built
    with, over `evaluator` -- the project's one (`Project.evaluator`), which every caller passes,
    so `foo`'s signature and its body agree. A JUST_INFO value belongs to no call. No closure
    is needed: a kept value names no other call's AUTO value in its declaration (an AUTO extent
    is spelled in suppliers, a supplier is never AUTO, and a computed length of a string the
    caller passes is refused).
    """
    names = {argument.name.lower() for argument in arguments}
    kept_calls = {
        producer_call_key(argument, arguments, conventions, evaluator)
        for argument in arguments
        if is_computed(argument) and _sizes_a_kept_argument(argument, arguments, conventions)
    }
    taken = []
    for argument in arguments:
        if is_temporary(argument, conventions):
            taken.append(argument)
        elif is_permutation(argument, conventions):
            base = argument.name.lower()[: -len(conventions.perm_suffix)]
            if base in names:
                taken.append(argument)
        elif is_computed(argument):
            if producer_call_key(argument, arguments, conventions, evaluator) not in kept_calls:
                taken.append(argument)
        elif is_prologue_produced(argument, prologue) and (
            not _sizes_a_kept_argument(argument, arguments, conventions)
        ):
            taken.append(argument)
    return taken


def sorted_permutations(
    taken: Sequence[Argument], prologue=None, conventions: Conventions = CONVENTIONS
) -> list[Argument]:
    """The permutations the allocating wrapper seeds and heapsorts itself.

    Every `<base>_perm` it took over, except the ones something else builds:

    - a `tmp_` permutation is the implementation's own scratch, seeded and ordered inside it;
    - a permutation the prologue declares `intent(out)` is the prologue's. That is what makes
      a non-default ordering expressible without giving up the expert tier: name the argument
      `<base>_perm` so an expert caller can still supply their own order, and let the prologue
      build it however this family needs. `intent(inout)` means the opposite -- the prologue
      refines an order, so the default sort still runs first and hands it one.

    Shared by the emitter (which emits the calls) and `validate` (which checks nothing the
    sort reads is written afterwards), so the two cannot disagree about what gets sorted.
    """
    built = {
        dummy.name.lower()
        for dummy in (prologue.arguments if prologue is not None else ())
        if dummy.intent is Intent.OUT
    }
    return [
        argument
        for argument in taken
        if is_permutation(argument, conventions)
        and not is_temporary(argument, conventions)
        and argument.name.lower() not in built
    ]


@dataclass(frozen=True)
class ModeFix:
    """The mode a per-mode wrapper fixes: the argument dropped, and what it is set to."""

    #: The implementation's mode argument, dropped from the wrapper's signature
    argument: str
    #: The parameter the mode is fixed to in the call, e.g. `MODE_DOSAGE_PATTERN`
    parameter: str
    #: The module the parameter lives in, so the wrapper can import it
    module: str


@dataclass(frozen=True)
class WrapperSpec:
    """A synthesised wrapper set and the implementation it came from."""

    impl: Procedure
    #: The generated module the wrappers live in
    module_name: str
    #: The wrapper that only validates. Named `<base>_expert` when `allocating` exists, and
    #: `<base>` when it is the only wrapper -- there being no second tier to distinguish it
    #: from. Read the name off the procedure; do not assume the suffix either way.
    validating: Procedure
    #: The wrapper that also allocates, when there is anything to take over. Always named
    #: `<base>`: the plain name belongs to the entry point a caller should reach for first.
    allocating: Procedure | None = None
    #: Set when this is one of a mode-split impl's per-mode wrappers
    mode_fix: ModeFix | None = None
    #: The procedure this implementation's `DM_PROLOGUE` names, resolved
    prologue: Procedure | None = None


@dataclass(frozen=True)
class SynthesisResult:
    project: Project
    specs: tuple[WrapperSpec, ...]
    #: names of the generated modules that only re-export their children (see
    #: `_reexport_modules`); they hold no procedures, so no `WrapperSpec` names them
    reexports: tuple[str, ...] = ()

    @property
    def expert_only_in_fortran(self) -> set[str]:
        """The names of the `foo_expert` wrappers that stay Fortran-and-C only.

        The name `foo_expert` promises control over what reaches the implementation. It can
        only keep that promise where `foo` does something an expert caller would want to do
        differently (`alloc_does_more`): sort a permutation, run a prologue, or compute a value
        the Python and R expert takes as a parameter. Where `foo` merely validates and
        allocates, the binding languages allocate for *both* tiers anyway, so `foo_expert`
        would be the same call under a name claiming otherwise.

        Fortran and C keep both: there the expert tier really does hand the buffers over, so
        a caller who wants to reuse them still can. Asked after the semantic pass, which is
        what resolves the producers `alloc_does_more` looks at.
        """
        return {
            spec.validating.name.lower()
            for spec in self.specs
            if spec.allocating is not None and not alloc_does_more(spec)
        }

    @property
    def generated_names(self) -> set[str]:
        """Every module this synthesis put into the project."""
        return {spec.module_name for spec in self.specs} | set(self.reexports)


def synthesize_wrappers(
    project: Project,
    conventions: Conventions = CONVENTIONS,
    diagnostics: DiagnosticBag | None = None,
) -> SynthesisResult:
    """Build the wrappers every implementation implies and inject them into `project`.

    An implementation is a `<name>_impl` procedure in a `<module>_impl` module. Its wrappers
    go into a generated module named by dropping the module's `_impl` suffix, so
    `tox_loess_impl` yields `tox_loess`. The module name is the whole trigger -- where the
    file sits decides only where its wrappers are *written* (`generated_path_for`), which is
    what lets an implementation live under any layer of `src/`.

    It is analysed here, before its wrappers are built, so its mode table is known: a mode
    argument whose table names a procedure per mode makes it expand into one wrapper (pair)
    per mode instead of a single procedure taking the mode. The implementation is not
    exported, so the later project-wide pass leaves it alone; the wrappers, cloned from it,
    are analysed there like any procedure.
    """
    diagnostics = diagnostics if diagnostics is not None else DiagnosticBag()
    evaluator = project.evaluator()
    suffix = conventions.impl_suffix
    by_module: dict[str, list[Procedure]] = {}
    sources: dict[str, Module] = {}
    specs: list[WrapperSpec] = []
    impl_modules = [m for m in project if m.name.lower().endswith(suffix)]

    for module in impl_modules:
        generated_name = module.name[: -len(suffix)]
        for impl in module.procedures:
            if not impl.name.lower().endswith(suffix):
                continue  # a recommend routine or private helper, not an implementation

            analyse(impl, diagnostics, conventions)
            prologue = _prologue_of(impl, project)
            mode_argument = _split_mode_argument(impl)

            for base, arguments, mode_fix in _variants(impl, mode_argument, conventions):
                validating, allocating = _wrappers_for(
                    base, arguments, impl, conventions, prologue, evaluator
                )
                # the plain entry point first: it leads the generated `public ::` block
                # and the subroutine order, which is the order a reader meets the tiers in
                wrappers = ([allocating] if allocating is not None else []) + [validating]
                by_module.setdefault(generated_name, []).extend(wrappers)
                sources[generated_name] = module
                specs.append(
                    WrapperSpec(
                        impl=impl,
                        module_name=generated_name,
                        validating=validating,
                        allocating=allocating,
                        mode_fix=mode_fix,
                        prologue=prologue,
                    )
                )

    generated_modules = [
        Module(
            name,
            procedures=procedures,
            # the implementation module's own documentation, verbatim. A generated module is
            # the published API -- what Python imports, what the R help pages are built from,
            # what a Fortran caller `use`s -- so what the author wrote about the family is what
            # a reader of any of those gets. It used to get a "do not edit" banner instead, on
            # the argument that an implementation module's doc describes the implementation;
            # the answer to that is that the doc has to be written as API prose, not that the
            # API should go undocumented. Each emitter adds its own generated-from note.
            doc=sources[name].doc,
            meta=sources[name].meta,
            location=sources[name].location,
        )
        for name, procedures in by_module.items()
    ]
    reexports = _reexport_modules(impl_modules, set(by_module), conventions)
    generated_modules += reexports

    # the same constants, so the same evaluator: what keyed the calls here keys them there
    augmented = Project(list(project.modules) + generated_modules, name=project.name,
                        evaluator=evaluator)
    return SynthesisResult(
        project=augmented,
        specs=tuple(specs),
        reexports=tuple(module.name for module in reexports),
    )


def _reexport_modules(
    impl_modules: Sequence[Module], wrapper_modules: set[str], conventions: Conventions
) -> list[Module]:
    """The generated counterparts of the implementation tree's re-export modules.

    A family too big for one file is split into several implementation modules and gathered by
    a parent that holds no procedures of its own and only `use`s its children -- the shape f42
    uses for its serde tree. That parent generates the matching parent over the *wrappers*, so
    `use tox_data_integration` still reaches the whole family and the split stays an
    implementation detail of the implementation tree.

    Only children that generate something are re-exported: an implementation module used for
    its constants or its recommend routines has no generated counterpart to `use`. Parents are
    resolved to a fixed point, so a parent may gather other parents.
    """
    suffix = conventions.impl_suffix
    candidates = {
        module.name[: -len(suffix)]: sorted(
            {
                used[: -len(suffix)]
                for used in module.uses
                if used.lower().endswith(suffix)
            }
        )
        for module in impl_modules
        if not module.procedures and module.uses
    }
    if not candidates:
        return []

    parents = {module.name[: -len(suffix)]: module for module in impl_modules}
    generated = set(wrapper_modules)
    while True:
        resolved = {
            name
            for name, children in candidates.items()
            if any(child in generated for child in children)
        }
        if resolved <= generated:
            break
        generated |= resolved

    return [
        Module(
            name,
            doc=parents[name].doc,
            meta=parents[name].meta,
            location=parents[name].location,
            uses=[child for child in children if child in generated],
        )
        for name, children in sorted(candidates.items())
        if name in generated
    ]


def alloc_does_more(spec: WrapperSpec, conventions: Conventions = CONVENTIONS) -> bool:
    """Whether the allocating wrapper does anything beyond validating and allocating.

    Three things qualify, and they are what an expert caller would be taking over: a
    permutation it seeds and heapsorts, which fixes one ordering; a prologue, which may derive
    an input or decide the call is degenerate and answer it directly; and a value the Python
    and R plain tier computes -- through `foo`, or through the binding -- that their expert
    takes as a parameter (E2): data reaching the implementation, such as an array a producer
    fills, is control over what it is handed. Sizing a buffer through `DM_OUTPUT_FROM` does not
    qualify -- that is part of allocating it, and the binding does it in every tier
    (`binding_computes`). Neither does validating, which both wrappers do.

    Evaluated after the semantic pass, where it is read (`expert_only_in_fortran`): whether a
    value is a parameter of the expert depends on which producer call fills it, and the
    producers are resolved there.
    """
    taken = taken_over_arguments(spec.validating.arguments, conventions, spec.prologue,
                                 project_evaluator(spec.impl))
    if spec.prologue is not None or sorted_permutations(taken, spec.prologue, conventions):
        return True
    expert = spec.validating
    # every AUTO value of a pair's expert is computed by its plain tier: by `foo` where it
    # took the value over, by the binding where `foo` kept it
    return any(
        automatic_plan(argument) is not None and not is_derived(argument, expert, conventions)
        for argument in expert.arguments
    )


# -- what a binding computes (E2) and where it reads an extent (R4) -----------------------
#
# The binding facts, which `roles` cannot state: `roles.is_computed` and `roles.is_derived`
# are what the Fortran `foo` does, and in a generated pair the Python and R expert tier does
# less. Every emitter and ABI rule that decides what a binding asks for, computes or reads
# asks these, never the IR facts, so that no two of them can disagree (`design/
# language-layers.md`, "What a binding computes"). Each takes the tier -- a `CWrapper` or the
# `Procedure` it wraps -- since one argument is computed in one tier and passed in another.
# They read resolved producer plans, so they hold after the semantic pass.


def _tier(wrapper) -> Procedure:
    """The IR procedure of a binding tier: a `CWrapper`'s, or the procedure itself."""
    return getattr(wrapper, "procedure", wrapper)


def _ir_argument(argument):
    """The IR argument behind a binding's (`CArgument.source`), or the argument itself."""
    return argument.source if hasattr(argument, "source") else argument


def automatic_plan(argument):
    """The `DM_OUTPUT_FROM(..., AUTO)` plan resolved for `argument`, or None.

    What a producer call is -- the producer, its inputs, the output taken. Whether a binding
    makes that call is `binding_computes`'s to say.
    """
    argument = _ir_argument(argument)
    roles = argument.roles if argument is not None else None
    plan = roles.computed_from if roles is not None else None
    return plan if plan is not None and plan.is_automatic else None


def is_generated_expert(procedure: Procedure, conventions: Conventions = CONVENTIONS) -> bool:
    """Whether `procedure` is the expert half of a generated pair."""
    module = procedure.module
    project = module.project if module is not None else None
    return project is not None and plain_of_an_expert(procedure, project, conventions) is not None


def _computed_by_the_expert(expert: Procedure, conventions: Conventions) -> set[str]:
    """The AUTO values a generated expert's binding computes, lower-cased.

    One that sizes a buffer the binding allocates -- its name is in an extent or the length
    of a work array or an output, other than as the array operand of `size()`, which reads the
    shape of an array the caller passed -- and every value the same producer call fills: a
    size can depend on the values that call computes (a count, a pack), so the two come from
    one call or the size need not fit the data.
    """
    allocated = [a for a in expert.arguments
                 if is_temporary(a, conventions) or a.intent is Intent.OUT]
    named = {named.name.lower() for a in allocated for named in size_names(a)
             if not named.is_size_operand}
    automatic = [(a, automatic_plan(a)) for a in expert.arguments]
    automatic = [(a, plan) for a, plan in automatic if plan is not None]
    calls = {plan.call for a, plan in automatic if a.name.lower() in named}
    return {a.name.lower() for a, plan in automatic if plan.call in calls}


def binding_computes(argument, wrapper, conventions: Conventions = CONVENTIONS) -> bool:
    """Whether the Python and R binding of `wrapper` computes `argument` through its producer (E2).

    In the expert tier of a generated pair, only a value that sizes a buffer the binding
    allocates, and whatever one call fills with it (`_computed_by_the_expert`); every other
    AUTO value is the caller's to pass there, as it is in Fortran and C. In every other tier --
    the plain one, a lone wrapper, a hand-written consumer -- every AUTO value it has.
    """
    argument, procedure = _ir_argument(argument), _tier(wrapper)
    if automatic_plan(argument) is None:
        return False
    if not is_generated_expert(procedure, conventions):
        return True
    return argument.name.lower() in _computed_by_the_expert(procedure, conventions)


def carries_computed(owner, wrapper, conventions: Conventions = CONVENTIONS) -> bool:
    """Whether the binding fills `owner` itself, so an extent read off it waits for the producers.

    An input the binding computes, or an output whose shape argument it computes (serde's
    `arr`, sized by the `arr_shape` its producer returns).
    """
    owner, procedure = _ir_argument(owner), _tier(wrapper)
    if owner.intent is not None and owner.intent.is_input:
        return binding_computes(owner, procedure, conventions)
    shape = owner.roles.shape_arg if owner.roles is not None else None
    return shape is not None and binding_computes(shape, procedure, conventions)


def _read_by_a_producer(argument: Argument, procedure: Procedure,
                        conventions: Conventions) -> bool:
    """Whether a producer the binding calls is passed `argument` (its `binding_inputs`)."""
    name = argument.name.lower()
    for other in procedure.arguments:
        plan = automatic_plan(other)
        if plan is None or not binding_computes(other, procedure, conventions):
            continue
        if any(supply.argument is not None and supply.argument.lower() == name
               for supply in plan.binding_inputs):
            return True
    return False


def readable_owners(argument, wrapper, conventions: Conventions = CONVENTIONS) -> tuple:
    """The arrays the binding of `wrapper` reads the extent `argument` off, the first one first (R4).

    1. A required array the caller passes, whenever one carries it: read before the producers,
       so any of them may be passed it. These come first; an optional one follows them, and is
       compared with it when given.
    2. Otherwise an array the binding fills itself (`carries_computed`), read after the
       producers -- but only where no producer is passed the extent, which would need it first.
    3. Otherwise none: the extent is asked of the caller, as one only outputs carry is, and no
       error names an array the caller never passed (`error_arguments.source_of`).

    An optional array never provides. Empty for anything that is no extent; `roles.readable_owners`
    is the IR's answer, which counts every input alike.
    """
    argument, procedure = _ir_argument(argument), _tier(wrapper)
    roles = argument.roles if argument is not None else None
    if roles is None:
        return ()
    passed, computed = [], []
    for owner in roles.readable_owners:
        (computed if carries_computed(owner, procedure, conventions) else passed).append(owner)
    required = [owner for owner in passed if not owner.optional]
    optional = [owner for owner in passed if owner.optional]
    if required:
        return tuple(required + optional)
    if computed:
        if _read_by_a_producer(argument, procedure, conventions):
            return ()
        return tuple(computed + optional)
    return tuple(optional)  # `validate` refuses an extent only optional arrays carry


def read_off_a_carrier(argument, wrapper, conventions: Conventions = CONVENTIONS) -> bool:
    """Whether the binding of `wrapper` reads the extent `argument` off an array that carries it.

    Only where it does not compute it: a tier whose binding computes a value through its
    producer passes C what the producer gave, and the carriers are checked against that
    (`shape_checks`), so nothing also reads it off them -- a value is computed or read, never
    both. The two sites that emit a derivation ask this (`c_call`, `python_ctypes`);
    `readable_owners` itself is unchanged, so an error at the extent still names the carrier
    (`error_arguments.source_of`).
    """
    return (not binding_computes(argument, wrapper, conventions)
            and bool(readable_owners(argument, wrapper, conventions)))


def is_derived(argument, wrapper, conventions: Conventions = CONVENTIONS) -> bool:
    """Whether the binding of `wrapper` works `argument` out itself, rather than asking for it.

    Derived in this tier: computed through its producer (`binding_computes`), read off a
    carrier (`readable_owners`), a shape read off its array, or a mask's count. The tier-aware
    counterpart of `roles.is_derived`, which says what the Fortran `foo` does.
    """
    argument, procedure = _ir_argument(argument), _tier(wrapper)
    roles = argument.roles if argument is not None else None
    if roles is None:
        return False
    return (binding_computes(argument, procedure, conventions)
            or bool(readable_owners(argument, procedure, conventions))
            or roles.is_inferable_shape_arg or roles.is_mask_count)


def project_evaluator(procedure: Procedure):
    """The one constant evaluator of the project `procedure` is in, or None outside one (a unit
    test's hand-built procedure). What every caller of `taken_over_arguments` passes."""
    module = procedure.module if procedure is not None else None
    project = module.project if module is not None else None
    return project.evaluator() if project is not None else None


def _prologue_of(impl: Procedure, project: Project):
    """The procedure an implementation's `DM_PROLOGUE` names, or None.

    Resolved here as well as in the emitter, because the drop set depends on it and the
    signature is shaped here while the body is written there -- both call
    `taken_over_arguments` with it, so the two cannot disagree about what the caller passes.
    A name that resolves to nothing is `validate`'s to report; here it simply means no
    prologue, which is what an unannotated impl gets anyway.
    """
    directive = impl.directives.prologue
    if directive is None:
        return None
    return project.procedure(directive.module, directive.procedure)


def _base_name(impl: Procedure, conventions: Conventions) -> str:
    return impl.name[: -len(conventions.impl_suffix)]


def _split_mode_argument(impl: Procedure) -> Argument | None:
    """The implementation's mode argument whose table opts into per-mode splitting, or None.

    Requires the implementation to have been analysed, so the mode table is on its roles.
    """
    for argument in impl.arguments:
        roles = argument.roles
        if roles is not None and roles.mode is not None and roles.mode.is_split:
            return argument
    return None


def _variants(
    impl: Procedure, mode_argument: Argument | None, conventions: Conventions
) -> list[tuple[str, list[Argument], ModeFix | None]]:
    """The (base name, exposed impl arguments, mode fix) of each wrapper to generate.

    One entry for an ordinary impl; one per mode value for a mode-split impl.
    """
    if mode_argument is None:
        return [(_base_name(impl, conventions), published_arguments(impl, conventions), None)]
    return [
        (
            value.procedure_name,
            _mode_kept_arguments(impl, mode_argument, value, conventions),
            ModeFix(mode_argument.name, value.parameter, value.module),
        )
        for value in mode_argument.roles.mode.values
    ]


def _mode_kept_arguments(
    impl: Procedure, mode_argument: Argument, value, conventions: Conventions
) -> list[Argument]:
    """The implementation arguments a per-mode wrapper exposes.

    The mode argument is dropped (fixed to `value.parameter` in the call). An argument that
    `DM_REQUIRED_IF_MODE` ties to this mode becomes a mandatory dummy; one tied to another
    mode is dropped; everything else is carried through.
    """
    kept: list[Argument] = []
    for argument in published_arguments(impl, conventions):
        if argument is mode_argument:
            continue
        required = argument.directives.required_if_mode
        if required is not None and required.mode_arg.lower() == mode_argument.name.lower():
            if required.mode_param.lower() == value.parameter.lower():
                kept.append(_as_required(argument))
            # otherwise it belongs to another mode and is absent here
        else:
            kept.append(argument.with_name(argument.name))
    return kept


def _mode_prose_line(argument: Argument) -> int | None:
    """The line `DM_REQUIRED_IF_MODE` wrote on this argument, if it can be identified.

    The directive's line number comes from a declaration-line + 1 heuristic in the source
    index, so it is confirmed against the line's own text before anything is deleted: a line
    that does not name the mode parameter is not the macro's, and survives. If the numbering
    assumption ever breaks, the failure is prose that stays rather than prose that vanishes.
    """
    directive = argument.directives.required_if_mode
    if directive is None or directive.line_number is None:
        return None
    for block in argument.doc.blocks:
        if block.line_number == directive.line_number and directive.mode_param in block.text:
            return directive.line_number
    return None


def _as_required(argument: Argument) -> Argument:
    """A copy of a mode-scoped optional, as this mode's wrapper takes it.

    Its `DM_REQUIRED_IF_MODE` directive is dropped: this wrapper is the mode, so the
    conditional-requirement relation no longer applies. What is left is a plain argument,
    mandatory -- unless it carries a `DM_DEFAULT`, which the binding supplies when the
    caller omits it. Such an argument is never *required*; there the directive only scopes
    it to its mode (absent from the others), and it stays optional.

    The prose the macro wrote goes with the directive. It conditions on a mode argument this
    wrapper does not have, and calls an argument the split has just made mandatory
    "optional" -- a Python reader sees `rdi_threshold : float` with "This optional argument
    needs to be passed if..." directly beneath it.
    """
    return Argument(
        argument.name,
        argument.type,
        dimension=argument.dimension,
        intent=argument.intent,
        optional=argument.directives.has_default,
        doc=argument.doc.without_line(_mode_prose_line(argument)),
        directives=replace(argument.directives, required_if_mode=None),
        attributes=argument.attributes,
        location=argument.location,
    )


def prologue_only_arguments(
    prologue, arguments: Sequence[Argument], conventions: Conventions = CONVENTIONS
) -> list[Argument]:
    """The prologue's own dummies -- the ones the implementation knows nothing about.

    A prologue derives something for the implementation, and what it derives it *from* is often
    not the implementation's business: a threshold comes from a quantile level, and the
    implementation takes the threshold. That level has to come from somewhere, so it
    becomes an argument of the allocating wrapper -- of that one alone, because only it runs
    the prologue. The expert tier takes the derived value directly and would have no use for
    it.

    `handled` and `ierr` are the wrapper's own and never count.
    """
    if prologue is None:
        return []
    known = {argument.name.lower() for argument in arguments}
    known |= {conventions.error_arg.lower(), conventions.prologue_handled_arg.lower()}
    return [
        dummy.with_name(dummy.name)
        for dummy in prologue.arguments
        if dummy.name.lower() not in known
    ]


def _wrappers_for(
    base: str, arguments: list[Argument], impl: Procedure, conventions: Conventions,
    prologue=None, evaluator=None,
) -> tuple[Procedure, Procedure | None]:
    """The validating wrapper, and the allocating one when the two would differ.

    They differ when the allocating one takes something over (a work array, a permutation, a
    recommend-sized value) or when the prologue asks for something of its own. Either way the
    caller sees a different signature, which is the whole reason for two entry points; where
    neither happens there is one wrapper and nothing to choose between.

    **The plain name always goes to the entry point a caller should reach for first.** With
    two of them that is the allocating one, and the validating one becomes `<base>_expert`;
    with one it is the validating one, because there is no second tier for the name to
    distinguish it from. Naming the lone wrapper `<base>_expert` unconditionally would
    generate, compile and bind perfectly well -- and quietly rename every entry point that
    has no work arrays to take over.
    """
    taken = taken_over_arguments(arguments, conventions, prologue, evaluator)
    extra = prologue_only_arguments(prologue, arguments, conventions)
    paired = bool(taken or extra)

    foo_arguments = [argument.with_name(argument.name) for argument in arguments]
    _append_error(foo_arguments, impl, conventions)
    validating = _wrapper(
        base + (conventions.expert_suffix if paired else ""),
        foo_arguments, impl, conventions,
    )
    if not paired:
        return validating, None

    # the prologue's own arguments go after the implementation's, before `ierr`: they are
    # the allocating tier's own vocabulary, not part of what the implementation asked for.
    # Before the implementation's own `ierr` too, where it has one, and before a function's
    # result, which is published right before `ierr`, so the signature still ends in both.
    kept = [a.with_name(a.name) for a in arguments if a not in taken]
    ending = {conventions.error_arg.lower()}
    if impl.result is not None:
        ending.add(impl.result.name.lower())
    error = next(
        (index for index, a in enumerate(kept) if a.name.lower() in ending),
        len(kept),
    )
    alloc_arguments = kept[:error] + extra + kept[error:]
    _append_error(alloc_arguments, impl, conventions)
    return validating, _wrapper(base, alloc_arguments, impl, conventions)


def _wrapper(
    name: str, arguments: list[Argument], impl: Procedure, conventions: Conventions
) -> Procedure:
    meta = Meta(
        summary=impl.meta.summary,
        author=impl.meta.author,
        category=conventions.c_binding_category,
    )
    return Procedure(
        name,
        arguments=arguments,
        doc=impl.doc,
        meta=meta,
        location=impl.location,
        conventions=conventions,
    )


def _append_error(
    arguments: list[Argument], impl: Procedure, conventions: Conventions
) -> None:
    if impl.argument(conventions.error_arg) is None:
        arguments.append(
            Argument(
                conventions.error_arg,
                FortranType(BaseType.INTEGER, "int32"),
                intent=Intent.OUT,
                doc=Doc.parse([_ERROR_DOC]),
            )
        )
