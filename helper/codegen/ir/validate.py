"""Rules an exported procedure must satisfy.

Only exported procedures are checked. An internal routine is held to none of this,
because nothing is generated from it -- it may use deferred-length strings, derived types
and untyped intents as freely as Fortran allows.

Each rule exists because generating correct code is otherwise impossible, and each says
what to do instead. A rule that merely encodes taste belongs in a review, not here.
"""

from __future__ import annotations

import re

from .. import names
from ..config import CONVENTIONS, Conventions
from ..diagnostics import DiagnosticBag
from .constants import KIND_VALUES
from .doc import EXTERNAL_PREFIX
from .entities import Argument, Module, Procedure, Project
from .types import BaseType, FortranType, Intent, split_top_level

#: Attributes on a dummy argument that have no C interoperable form
NON_INTEROPERABLE_ATTRIBUTES = ("allocatable", "pointer")

#: The one kind, by value, each numeric type may have in an exported signature. R holds an
#: integer as int32 and a number as real64 (a complex as two of them), and the generated
#: validators exist for those kinds only; anything else reaches R as the wrong storage -- an
#: int64 output is written eight bytes per element into a four-byte INTSXP -- and the
#: validator call does not compile. Compared by value, so `c_int` is `int32` and `c_double`
#: is `real64`, and `integer(8)` is int64; a kind whose value the generator does not know is
#: the C ABI's to refuse. Logicals are compared by name (`SUPPORTED_LOGICAL_KINDS`).
SUPPORTED_KIND_VALUES = {BaseType.INTEGER: 4, BaseType.REAL: 8, BaseType.COMPLEX: 8}

#: The kinds a logical may have, by *name*: `c_bool`, which C passes as it is, and the default
#: kind, which the C wrapper converts. Not by value: any other spelling -- `logical(1)`, though
#: it is c_bool's value -- would be converted into a default logical and handed to a dummy of
#: another kind, which does not compile.
SUPPORTED_LOGICAL_KINDS = (None, "c_bool")


def bindings_cannot_carry(fortran: FortranType) -> bool:
    """Whether the bindings cannot carry a value of this type (#221 B4).

    A numeric kind is compared by value with `SUPPORTED_KIND_VALUES`, a logical kind by name
    with `SUPPORTED_LOGICAL_KINDS`. A numeric kind whose value is unknown is not this rule's:
    the C ABI refuses it as having no C counterpart. Shared with `abi.c_abi`, which stays
    silent about a kind this already refuses, so one cause gets one error.
    """
    if fortran.is_logical:
        return ((fortran.kind or "").lower() or None) not in SUPPORTED_LOGICAL_KINDS
    supported = SUPPORTED_KIND_VALUES.get(fortran.base)
    if supported is None:
        return False
    value = _kind_value(fortran.kind)
    return value is not None and value != supported


def function_result_problems(procedure: Procedure,
                             conventions: Conventions = CONVENTIONS) -> list[tuple[str, str]]:
    """Why a function's result cannot be published as an output argument, as (message, note).

    Every published signature is a subroutine's: an implementation function gets Fortran tiers
    that are subroutines, and an exported function a C wrapper that is one. Its result becomes
    an `intent(out)` argument right before `ierr` (#221 B1), under the result's own name, so it
    is held to everything an `intent(out)` dummy is -- and these rules come on top, because a
    result can break them where a dummy cannot:

    - no `result(...)` clause: the result is then named after the function, and every wrapper
      that calls the function names it too, so the dummy and the procedure it calls would be
      one name (gfortran: "Symbol ... conflicts"; ifx: error #6401). The stripped name is no
      way out, being the wrapper's own;
    - named like `ierr`: the two would be one dummy (gfortran: "Duplicate symbol");
    - named with the `tmp_` prefix: a work array, which the bindings allocate and never return,
      and which the allocating tier takes over -- the function's only output, swallowed;
    - an assumed length (`len=*`): the result's length comes from the caller's declaration,
      which no binding writes. A deferred one is refused as on any dummy.

    Empty for a subroutine. `synthesize.published_result` leaves a result with any of these
    out of the wrappers, so nothing downstream builds the signature this refuses.
    """
    result = procedure.result
    if result is None:
        return []
    problems = []
    if result.name.lower() == procedure.name.lower():
        problems.append((
            f"function '{procedure.name}' has no result(...) clause",
            "its result is published as an output argument under its own name, and the "
            "wrapper that calls the function names the function too -- named after it, the "
            "two are one name and the wrapper does not compile. Name the result: "
            f"`function {procedure.name}(...) result(<name>)`",
        ))
    elif result.name.lower() == conventions.error_arg.lower():
        problems.append((
            f"result '{result.name}' is named like the error code",
            f"the result is published as an output argument right before "
            f"'{conventions.error_arg}', which every published signature ends in, so the two "
            "would be one name. Name the result for what it is",
        ))
    elif result.name.lower().startswith(conventions.temporary_prefix):
        problems.append((
            f"result '{result.name}' has the work-array prefix "
            f"'{conventions.temporary_prefix}'",
            "a work array is allocated by the bindings and never returned, so the function's "
            "result would be thrown away; drop the prefix",
        ))
    length = result.type.length if result.type is not None else None
    if length is not None and length.is_assumed:
        problems.append((
            f"result '{result.name}' has an assumed length (len=*)",
            "an assumed-length result takes its length from the caller's declaration, which "
            "no binding writes (and the form is obsolescent); declare a fixed length, len=n",
        ))
    return problems


def _check_function_result(procedure: Procedure, diagnostics: DiagnosticBag,
                           conventions: Conventions) -> None:
    for message, note in function_result_problems(procedure, conventions):
        diagnostics.error(message, entity=procedure.result, note=note)


#: The greatest character rank the C layer can carry. A character's length is a leading C
#: extent, so a scalar character is a 1-D c_char buffer and a character vector a 2-D one.
#: A rank-2 character array would need a rank-3 buffer, and the wrapper's pointer view of
#: it would have to be rank 2 -- which is a wire format, not a limit of one helper.
MAX_CHARACTER_RANK = 1


def validate_project(project: Project, diagnostics: DiagnosticBag,
                     conventions: Conventions = CONVENTIONS) -> None:
    for module in project:
        validate_module(module, diagnostics, conventions)
    check_doc_links(project, diagnostics)
    check_whitelist_is_allocation_free(project, diagnostics, conventions)
    check_c_kinds_depend_on_the_safeguard(project, diagnostics)
    check_extents_are_declared_first(project, diagnostics)
    check_exported_names_are_unique(project, diagnostics)


def check_exported_names_are_unique(project: Project, diagnostics: DiagnosticBag) -> None:
    """No two exported procedures share a name, ignoring case.

    The name is the C symbol (`foo_c`), the R shim's `.Call` entry (`foo_call`), the name the
    Python package re-exports and the name in R's one environment. Two modules exporting
    `foo` link as `multiple definition of 'foo_c'`; and were it not for the link, Python and R
    would quietly let the second shadow the first. Ignoring case because a Fortran caller
    cannot tell `foo` from `FOO` either. Within one module the same can come from a mode
    table, whose third column names a generated procedure that something else may already be
    called.
    """
    seen: dict[str, Procedure] = {}
    for module in project:
        for procedure in module.exported_procedures:
            first = seen.setdefault(procedure.name.lower(), procedure)
            if first is procedure:
                continue
            where = (
                f"both '{first.module.name}' and '{module.name}' export"
                if first.module is not module else f"'{module.name}' exports twice"
            )
            diagnostics.error(
                f"{where} a procedure called '{procedure.name}'",
                entity=procedure,
                note=(
                    "an exported name is the C symbol, the R `.Call` entry and the name the "
                    "Python package and the R environment publish, so it has to be unique "
                    "across the whole library, ignoring case; rename one of them"
                ),
            )


def allocations_in(procedure: Procedure) -> list[str]:
    """What a procedure allocates: its `allocatable` locals and its `allocatable` dummies.

    One definition, because two rules ask it -- `_check_impl_allocates` of an implementation,
    and `check_whitelist_is_allocation_free` of the infrastructure an implementation may
    reach. Asking it twice in two spellings is how the two would come to disagree about what
    counts, and the second rule exists precisely to close the gap the first leaves.

    A `pointer` local does not count; only `allocatable` does, which is a complete proxy here
    because `M_ALLOCATE` needs one.
    """
    dummies = [a.name for a in procedure.arguments if "allocatable" in a.attributes]
    return list(procedure.allocatable_locals) + dummies


#: The module whose compile-time guards say which C kind a platform is missing.
SAFEGUARD_MODULE = "f42_safeguard"


def check_c_kinds_depend_on_the_safeguard(project: Project, diagnostics: DiagnosticBag) -> None:
    """A module naming a C kind `use`s `f42_safeguard`.

    `iso_c_binding` gives a C kind the value -1 where the platform has no counterpart, and a
    negative kind is a compile error at every declaration naming it. `f42_safeguard` turns that
    into one error saying *which* kind is missing -- but only if it is compiled first, and the
    only way to order a Fortran build is a dependency. A module that names a C kind without
    depending on the safeguard compiles before it and fails with the compiler's message instead.

    Checked because it has already regressed once. Before the `c_bool` sweep the only modules
    naming a C kind were the generated C wrappers, every one of which had the `use`; the sweep
    put C kinds in 39 more compilation units and left the dependency behind, silently. The rule
    was in the safeguard's own comment the whole time.

    Reads argument kinds, which is what the IR carries. A module whose only C kind is on a
    module variable or a local is invisible here -- `f42_config` was exactly that case -- so
    this narrows the gap rather than closing it.
    """
    # Nothing to depend on, nothing to require: a project parsed without the safeguard -- a
    # unit-test fixture tree -- cannot satisfy this and is not wrong for that.
    if not any(m.name.lower() == SAFEGUARD_MODULE for m in project):
        return

    #: a synthesized wrapper is `<impl>` beside the `<impl>_impl` it came from. Its `use` list
    #: is the emitter's to write and is not in the IR at all, so checking it here would report
    #: what the emitter has already done. `emit/fortran_wrapper` handles those.
    synthesized = {
        m.name[: -len("_impl")].lower()
        for m in project
        if m.name.lower().endswith("_impl")
    }
    for module in project:
        if module.name.lower() == SAFEGUARD_MODULE or module.name.lower() in synthesized:
            continue
        offenders = sorted({
            argument.type.kind
            for procedure in module.procedures
            for argument in procedure.arguments
            if argument.type.kind and argument.type.kind.lower().startswith("c_")
        })
        if not offenders:
            continue
        reached = {u.lower() for u in module.uses}
        reached.update(u.lower() for p in module.procedures for u in p.uses)
        if SAFEGUARD_MODULE in reached:
            continue
        kinds = ", ".join(f"'{k}'" for k in offenders)
        diagnostics.error(
            f"module '{module.name}' names {kinds} without using '{SAFEGUARD_MODULE}'",
            entity=module,
            note=(
                f"add 'use {SAFEGUARD_MODULE}'. Its guards report which C kind a platform is "
                "missing, and they only get to speak if it is compiled first -- which a "
                "dependency is the only way to arrange. Without it this module compiles first "
                "and fails with the compiler's message about the kind instead"
            ),
        )


def check_extents_are_declared_first(project: Project, diagnostics: DiagnosticBag) -> None:
    """An extent -- or a character length -- is declared above the argument it sizes.

    `real(real64) :: array(n)` written above `integer(int32) :: n` is a GNU extension, not
    standard Fortran: gfortran under `-std=` reports `Symbol 'n' is used before it is typed`,
    then falls back to implicit typing for the array and buries the real message under a
    cascade. Writing both in one statement -- `:: a(n), b(n), n` -- is rejected the same way,
    so an extent on the *same* line as the array it sizes is a finding too.

    The generator already emits its own wrappers in this order (`emit/fortran_c` and
    `emit/fortran_wrapper` both sort extents ahead of their arrays). This is the same rule
    turned back on the hand-written sources, which had drifted from it in five files until
    the diagnostics profile gained `-std=f2023` and stopped compiling.

    Every procedure, not only the exported ones: this is about whether the file compiles, not
    about what crosses into C. Argument *declaration* order is what matters and it is
    independent of the dummy-argument list, so the line each declaration was found on is the
    only thing that can answer it -- an argument whose line is unknown is skipped rather than
    guessed at.
    """
    for module in project:
        for procedure in module.procedures:
            lines = {
                argument.name.lower(): argument.location.line
                for argument in procedure.arguments
                if argument.location is not None and argument.location.line is not None
            }
            for argument in procedure.arguments:
                if argument.location is None or argument.location.line is None:
                    continue
                declared = argument.location.line
                # A character length sizes the argument exactly as an extent does, and
                # `character(len=clen) :: a(n)` above `integer :: clen` is rejected for the
                # same reason -- `asserts.F90` had precisely that shape.
                sizes = list(argument.dimension.extents)
                length = argument.type.length
                if length is not None and not length.is_constant and length.expr:
                    sizes.append(length.expr)
                for extent in sizes:
                    for identifier in _IDENTIFIER_RE.findall(extent):
                        sizer = lines.get(identifier.lower())
                        if sizer is None or sizer < declared:
                            continue
                        # Not "below" vs "in the same statement": the extent's own
                        # location is looked up by searching for its name, and in exactly
                        # this broken ordering it matches the `array(n)` above its real
                        # declaration. One wording that holds either way.
                        diagnostics.error(
                            f"'{argument.name}' is sized by '{identifier}', which is not "
                            f"declared above it",
                            entity=argument,
                            note=(
                                "declare the extent first. Referring to a symbol that is "
                                "typed later in the specification part is a GNU extension "
                                "that gfortran rejects under -std=, and the array then takes "
                                "an implicit type, which buries the real error"
                            ),
                        )


def check_whitelist_is_allocation_free(project: Project, diagnostics: DiagnosticBag,
                                       conventions: Conventions = CONVENTIONS) -> None:
    """Every module an implementation may `use` allocates nothing.

    `_check_impl_allocates` reads declarations, so it sees a procedure allocating for itself
    and nothing else. What makes it hold *across* modules is `_check_impl_imports`, which
    bounds an implementation to other implementations and `impl_import_whitelist` -- the first
    are held to the same rule, and the second was a curated list nothing verified. This is the
    verification. Without it the bound is an assertion: whitelist a module, let it grow an
    allocation, and every implementation reaching it allocates with nothing to say so.

    It is checkable now and was not before. The list carried `f42_utils` while its `f42_stats`
    child still had hand-written `_alloc` halves, and `tox_conversions` until its string
    conversions became pointer views -- so this shipped as a standing failure twice over, and
    was twice deferred. Both are resolved, so the cheap module-level form is enough.

    Only what the project parsed is checked. The intrinsic modules are on the list too and
    have no source here; an entry that names nothing is not an error, because the import rule
    already refuses an implementation that reaches for a module that does not exist.
    """
    whitelisted = {name.lower() for name in conventions.impl_import_whitelist}
    for module in project:
        if module.name.lower() not in whitelisted:
            continue
        for procedure in module.procedures:
            offenders = allocations_in(procedure)
            if not offenders:
                continue
            names = ", ".join(f"'{name}'" for name in offenders)
            diagnostics.error(
                f"whitelisted module '{module.name}' allocates in '{procedure.name}': {names}",
                entity=procedure,
                note=(
                    "an implementation may use this module, and an implementation allocates "
                    "nothing -- a rule that only reaches one file unless everything it may "
                    "reach is allocation-free too. Either take the allocation out (a pointer "
                    "view into the caller's memory is usually what was wanted), or drop the "
                    f"module from 'impl_import_whitelist' -- a module that has to be "
                    "whitelisted may be one that should be an implementation"
                ),
            )


def check_doc_links(project: Project, diagnostics: DiagnosticBag) -> None:
    """Every Ford `[[...]]` cross-reference names something that exists.

    A dangling link is the failure with no symptom. Ford renders it as the literal text an
    author typed; `emit/doc_links.py` falls back to plain code, deliberately, because an R
    `\\link` to a missing topic fails `R CMD check`. So a link left pointing at a renamed
    module keeps generating, keeps compiling, keeps passing -- and quietly stops being a link
    in four languages at once. Converting f42 to `_impl` modules broke twenty-three of them in
    one pass and nothing said a word for a month.

    Errors, not warnings. The generator's documentation *is* the API's documentation in four
    languages, and a link is the one piece of it that can be checked mechanically -- so it is,
    at the same weight as a signature the C layer cannot express. It also has the property that
    decides these: nothing downstream can tell a broken link from a working one, so nothing
    downstream will ever report it.

    Resolved against the whole project, generated modules included -- a link into
    `[[f42_binary_search_tree(module):build_bst_index]]` is *correct*: it names the wrapper a
    reader can call, and `doc_links` turns it into that binding's name.

    What it cannot see: documentation attached to something the IR does not model -- a derived
    type's components, a private routine's own comment. Those are read by nobody downstream
    either, so the gap costs a check on prose that reaches no output.
    """
    modules = {module.name.lower(): module for module in project}
    # One authored line, one diagnostic. A generated wrapper carries the implementation's
    # documentation verbatim, so a single bad link in an argument comment reaches here once
    # from the implementation and once per wrapper generated from it -- four times, for a
    # procedure with both tiers. They all point at the same file and line, which is what makes
    # them the same finding.
    reported: set[tuple] = set()
    for module in project:
        for entity, doc in _documented(module):
            for link in doc.links:
                where = (entity.location.file, entity.location.line, str(link))
                if where in reported:
                    continue
                reported.add(where)
                _check_link(link, entity, modules, diagnostics)


def _documented(module: Module):
    """Every entity of `module` that carries documentation, with its doc."""
    yield module, module.doc
    for procedure in module.procedures:
        yield procedure, procedure.doc
        for argument in procedure.arguments:
            yield argument, argument.doc
        if procedure.result is not None:
            yield procedure.result, procedure.result.doc
    for parameter in module.parameters:
        yield parameter, parameter.doc
    for declaration in module.declarations:
        yield declaration, declaration.doc


def _check_link(link, entity, modules, diagnostics: DiagnosticBag) -> None:
    # `[[zlib(extmodule)]]` names an entity of another project entirely. Nothing here can
    # resolve one, and reporting it missing would be a lie.
    if (link.component_type or "").lower().startswith(EXTERNAL_PREFIX):
        return
    named = modules.get(link.component.lower())
    if named is None:
        # a link whose component is not a module names an entity directly -- Ford resolves
        # that against the whole project, so accept it wherever any module publishes the name
        if link.item is None and any(
            link.component.lower() in m.public_names for m in modules.values()
        ):
            return
        diagnostics.error(
            f"documentation links to '{link.component}', which is not a module here",
            entity=entity,
            note=(
                f"{link} resolves to nothing: Ford renders it as this literal text, and the "
                "Python and R bindings fall back to plain code. Name the module that defines "
                "it, or drop the link and write the name in backticks"
            ),
        )
        return
    if link.item is not None and link.item.lower() not in named.public_names:
        diagnostics.error(
            f"documentation links to '{link.item}', which '{named.name}' does not publish",
            entity=entity,
            note=(
                f"{link} resolves to nothing -- the name is private, or it moved. Ford "
                "documents only what a module makes public, so a link to a private routine "
                "dangles as surely as a link to one that never existed"
            ),
        )


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
        # the wrappers check the result they publish as any argument; what they cannot see is
        # a result they had to leave out (an exported one is `_Validator`'s, as any export)
        if not procedure.is_exported:
            _check_function_result(procedure, diagnostics, conventions)
        _check_taken_over_is_not_an_output(procedure, diagnostics, conventions)
        _check_prologue(procedure, diagnostics, conventions)


def _check_module_is_named_for_its_file(module: Module, diagnostics: DiagnosticBag) -> None:
    """An implementation module lives in a file named after it.

    Two independent things derive from that name and they have to agree. Synthesis triggers
    on the *module* name and names the generated module after it; the wrapper's *path* is
    mirrored from the *file* name. Let the two diverge and `module tox_x` is written into a
    file named for something else -- or, in a file that is not an `_impl` file, nowhere.

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
    offenders = allocations_in(procedure)
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
#: `_alloc` is the hand-written pair convention the framework used before any wrapper was
#: generated. Nothing reads it any more, which is exactly why it stays reserved: `foo_alloc`
#: beside a generated `foo` would read to every author as the allocating tier while being an
#: ordinary second procedure, and the shape is one nobody should reintroduce.
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


def _check_taken_over_is_not_an_output(impl: Procedure, diagnostics: DiagnosticBag,
                                       conventions: Conventions) -> None:
    """What the allocating wrapper prepares itself cannot also be what the implementation returns.

    `foo` takes over a `<base>_perm` whose `<base>` is an argument (it seeds and heapsorts
    it) and a `DM_OUTPUT_FROM(..., AUTO)` value (it calls the producer), declaring each as a
    local and handing it in (`synthesize.taken_over_arguments`). An `intent(out)` one is the
    implementation's *result* -- an argsort, a count -- and taking it over swallows it: `foo`
    computes it and throws it away, and Python and R return nothing. Only the expert tier would
    return it. A `tmp_` work array is scratch and stays exempt, and a prologue only takes over
    what the implementation reads.
    """
    from ..synthesize import (  # local: synthesize imports this module
        _prologue_of,
        is_computed,
        is_permutation,
        is_temporary,
        published_arguments,
        taken_over_arguments,
    )

    project = impl.module.project if impl.module is not None else None
    prologue = _prologue_of(impl, project) if project is not None else None
    arguments = published_arguments(impl, conventions)
    for argument in taken_over_arguments(arguments, conventions, prologue):
        if argument.intent is not Intent.OUT or is_temporary(argument, conventions):
            continue
        if is_permutation(argument, conventions):
            base = argument.name[: -len(conventions.perm_suffix)]
            how = f"builds it from '{base}' itself"
            hint = (
                f"a '<base>{conventions.perm_suffix}' whose base is an argument is seeded and "
                "sorted by the allocating wrapper"
            )
        elif is_computed(argument):
            producer = argument.directives.output_from.procedure
            how = f"computes it with '{producer}' itself"
            hint = "a DM_OUTPUT_FROM(..., AUTO) value is computed by the allocating wrapper"
        else:
            continue
        diagnostics.error(
            f"'{argument.name}' is an output, but the allocating wrapper {how}",
            entity=argument,
            note=(
                f"{hint}, which would then discard what the implementation computes. Name the "
                "result for what it is (`order`, `ranks`) or drop the directive, or make it "
                "intent(in) or intent(inout) if the implementation really is handed it"
            ),
        )


def _what(argument: Argument) -> str:
    """How a diagnostic names it: a function's result is published as an argument, but the
    author wrote a result."""
    return "result" if argument.is_result else "argument"


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
        self._check_procedure_name()
        self._check_error_is_last()
        _check_function_result(self.procedure, self.diagnostics, self.conventions)
        for argument in self.arguments:
            self._check_argument(argument)

    def _check_procedure_name(self) -> None:
        """A procedure is a Python and an R function of the same name, so it is held to their words.

        A keyword is a SyntaxError (and `import tensor_omics` fails with it) or an R parse
        error; a builtin's name hides the builtin for the whole generated module, where the
        other wrappers call it. Its C symbols carry a suffix (`foo_c`, `foo_call`), which no C
        keyword has, so C's list does not apply.
        """
        word = names.language_word(self.procedure.name, languages=("Python", "R"))
        if word is None:
            return
        self.diagnostics.error(
            f"procedure '{self.procedure.name}' is {word}",
            entity=self.procedure,
            note=(
                "a procedure is published as a Python and an R function of the same name; "
                "rename it (an implementation, its `_impl` base name)"
            ),
        )

    def _check_error_is_last(self) -> None:
        """The error code ends the published signature.

        Every published signature -- each Fortran tier, the C wrapper, and through them the
        Python and R call order -- ends in `ierr`, where it is added when a procedure has
        none of its own. One that declares its own somewhere else would publish it in the
        middle, and only for this procedure: a caller passing arguments by position gets them
        shifted by one. A synthesised wrapper takes the implementation's order, so this also
        reaches an implementation's `ierr`.
        """
        error = self.procedure.argument(self.conventions.error_arg)
        if error is None or self.procedure.arguments[-1] is error:
            return
        self.error(
            f"'{error.name}' is not the last argument of '{self.procedure.name}'",
            error,
            note=(
                "every published signature ends in the error code, which is where the "
                "generator adds it to a procedure that has none -- declare it last"
            ),
        )

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
        self._check_name(argument)
        self._check_optional_output(argument)
        self._check_optional_inout(argument)
        self._check_attributes(argument)
        self._check_type(argument)
        self._check_inout_conversion(argument)
        self._check_temporary(argument)
        self._check_mode(argument)
        self._check_required_if_mode(argument)
        self._check_default_on_an_array(argument)
        self._check_shape_argument(argument)
        self._check_extent(argument)
        self._check_mask_count(argument)
        self._check_result_size(argument)
        self._check_lower_bounds(argument)
        self._check_error_argument(argument)

    def _check_name(self, argument: Argument) -> None:
        """A published name every binding language takes, and no procedure the wrapper names.

        The name is a Python parameter, an R parameter, a parameter of the R `.Call` shim and
        a dummy of every generated wrapper, so the words those languages reserve are refused
        (`names.language_word`). So are the procedures a generated wrapper names -- itself, its
        C symbol, and a `DM_OUTPUT_FROM` producer it calls -- which a dummy would hide. The
        locals and imports of the generated code are another matter: the R shim's and
        Python's start with `_` and so meet no argument, and the Fortran wrappers' are checked
        where they are emitted (`emit.name_clashes`). A function's result is published as an
        output argument, so it is held to the same.
        """
        name = argument.name
        what = _what(argument)
        word = names.language_word(name)
        if word is not None:
            self.error(
                f"{what} '{name}' is {word}",
                argument,
                note=(
                    "a published argument is a parameter in Python, R and the C of the R "
                    "binding, so it cannot be a word one of them reserves; rename it"
                ),
            )
        called = self._called_procedure(name)
        if argument.is_result and name.lower() == self.procedure.name.lower():
            called = None  # no `result(...)` clause, which `function_result_problems` says
        if called is not None:
            self.error(
                f"{what} '{name}' has the name of {called}",
                argument,
                note=(
                    "the generated wrapper names that procedure, and inside it the argument "
                    "would hide it; rename the argument"
                ),
            )

    def _called_procedure(self, name: str) -> str | None:
        """What a generated wrapper of this procedure names that is called `name`, or None."""
        lowered = name.lower()
        procedure = self.procedure.name
        if lowered == procedure.lower():
            return f"the procedure '{procedure}' itself"
        if lowered == f"{procedure}{self.conventions.c_suffix}".lower():
            return f"the C wrapper of '{procedure}'"
        for argument in self.arguments:
            directive = argument.directives.output_from
            if directive is not None and directive.procedure.lower() == lowered:
                return f"'{directive.procedure}', which computes '{argument.name}'"
        return None

    def _check_optional_inout(self, argument: Argument) -> None:
        """An optional argument that is modified in place cannot be bound.

        Python passes a temporary for an optional scalar and keeps no result for it, so the
        new value is lost; R points straight into the caller's vector, which breaks R's
        copy-on-modify for every other object sharing it, and returns a stale copy. There is
        no shape that honours both "may be absent" and "is returned" in the bindings either
        -- the same reason an optional output is refused.
        """
        if argument.intent is not Intent.INOUT or not argument.optional:
            return
        self.error(
            f"argument '{argument.name}' is optional and intent(inout)",
            argument,
            note=(
                "the bindings cannot return a value the caller may omit. Make it a required "
                "intent(inout), or an optional intent(in) beside a required output"
            ),
        )

    def _check_inout_conversion(self, argument: Argument) -> None:
        """A string, or a logical array, modified in place cannot be bound.

        Both cross into Fortran through a buffer the binding builds -- R's logicals are four
        bytes wide and its strings are not blank-padded -- and nothing copies that buffer
        back: the R shim names its result after the input it already has (a compile error
        that fails the library), and Python discards the buffer, so the modification is lost.
        """
        if argument.intent is not Intent.INOUT:
            return
        if argument.type.is_character:
            what = "a character"
        elif argument.type.is_logical and argument.is_array:
            what = "a logical array"
        else:
            return
        self.error(
            f"argument '{argument.name}' is {what} that is intent(inout)",
            argument,
            note=(
                "the bindings convert it into a buffer of their own and cannot hand the "
                "modification back; take it intent(in) and return the result as a separate "
                "intent(out) argument"
            ),
        )

    def _check_default_on_an_array(self, argument: Argument) -> None:
        """A `DM_DEFAULT` is one value, and an array is not.

        The binding passes the default in the argument's place, so an omitted array arrives
        as a one-element array and fails the shape check against its extent.
        """
        if not argument.is_array or not argument.directives.has_default:
            return
        self.error(
            f"array argument '{argument.name}' has a DM_DEFAULT",
            argument,
            note=(
                "a default is a single value the bindings pass in the argument's place, so "
                "an omitted array would arrive with one element. Leave it optional without a "
                "default and let the implementation fill it in when it is absent"
            ),
        )

    def _check_documented(self, argument: Argument) -> None:
        # a result too: it is published as an output argument, documented like any (#221 B1)
        if not argument.doc:
            self.diagnostics.warn(
                f"{_what(argument)} '{argument.name}' has no documentation",
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
                    f"{_what(argument)} '{argument.name}' is '{attribute}', which has no C "
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
        self._check_kind(argument)

    def _check_kind(self, argument: Argument) -> None:
        """Only the kinds R and the validators can carry (`bindings_cannot_carry`)."""
        fortran = argument.type
        if not bindings_cannot_carry(fortran):
            return
        self.error(
            f"{_what(argument)} '{argument.name}' is {fortran}, which the bindings cannot carry",
            argument,
            note=(
                "R holds integers as int32 and numbers as real64, and the generated "
                "validators take only those kinds: declare it integer(int32), real(real64), "
                "complex(real64), logical(c_bool) or a default logical"
            ),
        )

    def _check_character(self, argument: Argument) -> None:
        length = argument.type.length
        if length.is_deferred:
            self.error(
                f"{_what(argument)} '{argument.name}' has a deferred length (len=:)",
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
                f"and the C layer carries up to rank {MAX_CHARACTER_RANK}",
                argument,
                note=(
                    "a character carries its length as a leading C extent, so a rank-2 "
                    "character array would need a rank-3 c_char buffer"
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
        owners = ", ".join(f"'{owner.name}'" for owner in roles.extent_of)
        if argument.optional:
            self.error(
                f"extent argument '{argument.name}' is optional, but {owners} needs it",
                argument,
                note=_NON_OPTIONAL_NOTE,
            )
        if argument.intent is not None and argument.intent is not Intent.IN:
            # Python and R read it off the arrays it sizes and pass that value; an extent the
            # procedure also writes is neither theirs to read nor theirs to return.
            self.error(
                f"extent argument '{argument.name}' is intent({argument.intent.value})",
                argument,
                note=(
                    f"an extent sizes {owners} on the way in, so it is intent(in); return "
                    "a count the procedure computes as a separate intent(out) argument"
                ),
            )
        readable = roles.readable_owners
        if readable and all(owner.optional for owner in readable):
            # the bindings drop it from the signature, since an array says what it is, and
            # then have nothing to read it off when that array is omitted
            names_ = ", ".join(f"'{owner.name}'" for owner in readable)
            self.error(
                f"extent argument '{argument.name}' can only be read off {names_}, which "
                f"may be omitted",
                argument,
                note=(
                    "Python and R read an extent off the arrays it sizes and do not ask for it; "
                    "make one of them required, or size the optional array by an extent that a "
                    "required argument also has"
                ),
            )

    def _check_mask_count(self, argument: Argument) -> None:
        """An `n_selected_<x>` is counted from its mask, so both are inputs, and both present.

        The bindings count the mask's `.true.` values and the wrapper checks the count against
        it before the call. With the mask an output there is nothing to count yet; with the
        count an output the check reads a value the caller never set (Python's
        UnboundLocalError, a compile error in the R shim, and a spurious
        `ERR_INVALID_INPUT` from Fortran). An omitted mask has nothing to count either.
        """
        roles = argument.roles
        if roles is None or roles.mask_count_of is None:
            return
        mask = roles.mask_count_of
        if mask.optional:
            self.error(
                f"mask '{mask.name}' is optional, but '{argument.name}' is counted from it",
                argument,
                note=(
                    "an omitted mask leaves nothing to count, in any layer. Make the mask "
                    f"required, or name the count something other than "
                    f"'{self.conventions.mask_count_prefix}<x>'"
                ),
            )
        if mask.intent is not Intent.IN or argument.intent is not Intent.IN:
            self.error(
                f"'{argument.name}' counts mask '{mask.name}', so both must be intent(in) "
                f"(they are intent({_intent(mask)}) and intent({_intent(argument)}))",
                argument,
                note=(
                    "the count is derived from the mask before the call, which needs both on "
                    "the way in. To return a selection and its size, name the count something "
                    f"other than '{self.conventions.mask_count_prefix}<x>'"
                ),
            )

    def _check_result_size(self, argument: Argument) -> None:
        """`DM_RESULT_SIZE_IS` names what the procedure reports back, not what it was given.

        The bindings trim the result to the count after the call, reading it as an output;
        an `intent(in)` count is the caller's own value and is never returned.
        """
        roles = argument.roles
        counter = roles.result_size_arg if roles is not None else None
        if counter is None or counter.intent is not Intent.IN:
            return
        self.error(
            f"the result size of '{argument.name}' is '{counter.name}', which is intent(in)",
            argument,
            note=(
                "DM_RESULT_SIZE_IS names the count the procedure fills in, as an output. A "
                f"size the caller chooses sizes the result directly: declare it "
                f"'{argument.name}({counter.name})'"
            ),
        )

    def _check_lower_bounds(self, argument: Argument) -> None:
        """An axis that does not start at 1 (`x(0:n)`, `x(0:)`).

        Every layer counts an axis from its extent alone: the validator's element count,
        the C wrapper's null check and the bindings' shapes would all read `0:n` as an
        expression and paste it where a number belongs. A lower bound of 1 is the default,
        and the frontend already drops it where it is written out.
        """
        bounded = [
            extent for extent in argument.dimension.extents
            if extent != ":" and len(split_top_level(extent, ":")) > 1
        ]
        if bounded:
            self.error(
                f"argument '{argument.name}' has a lower bound other than 1 "
                f"({', '.join(bounded)})",
                argument,
                note=(
                    "every layer sizes an axis from its extent alone, so an axis must start "
                    "at 1: declare the extent (`x(n + 1)`) and shift the index inside the "
                    "implementation"
                ),
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
    from ..synthesize import prologue_only_arguments, published_arguments, taken_over_arguments

    arguments = published_arguments(impl, conventions)
    if taken_over_arguments(arguments, conventions, prologue):
        return
    if prologue_only_arguments(prologue, arguments, conventions):
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

    And it declares it a default `logical` (`fortran_wrapper._prologue_local`), which it passes
    by keyword: argument association does not convert kinds, so a prologue's `logical(c_bool)`
    passed the generator and failed to compile (#221 WP8).
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
    elif handled.type.kind is not None:
        diagnostics.error(
            f"prologue '{prologue.name}' declares '{name}' as {handled.type}, not as a default "
            f"logical", handled,
            note=(
                f"the wrapper declares its flag `logical :: {name}` and passes it, and argument "
                f"association does not convert kinds, so this would not compile (\"Type mismatch "
                f"in argument '{name}'\"); declare `logical, intent(out) :: {name}`. It never "
                f"crosses to C, so it needs no C kind"
            ),
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
    own vocabulary and no business of the implementation's: a threshold's `quantile_level`.

    Which leaves a misspelling nowhere to be caught, since it now reads as a new argument.
    So a name that is *nearly* one the implementation has is refused: `n_gene` beside `n_genes` is a
    typo, and silently turning it into something the caller must pass would mean the prologue
    and the implementation working from different numbers.

    There is no rename table, deliberately: a producer needs one because it is a published
    routine whose own parameter names cannot move, whereas a prologue is internal to the
    impl module and its dummy can simply be renamed to match.
    """
    from ..synthesize import published_arguments  # local: synthesize imports this module

    # a function's result too: the wrappers publish it, so the prologue may answer it
    supplied = {argument.name.lower() for argument in published_arguments(impl, conventions)}
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
    # local: synthesize imports this module
    from ..synthesize import published_arguments, taken_over_arguments

    arguments = published_arguments(impl, conventions)
    taken = {a.name.lower() for a in taken_over_arguments(arguments, conventions, prologue)}
    by_name = {a.name.lower(): a for a in arguments}
    for dummy in prologue.arguments:
        if dummy.intent is None or not dummy.intent.is_output:
            continue
        argument = by_name.get(dummy.name.lower())
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
    from ..synthesize import published_arguments, sorted_permutations, taken_over_arguments

    written = {
        dummy.name.lower()
        for dummy in prologue.arguments
        if dummy.intent is not None and dummy.intent.is_output
    }
    written -= {conventions.error_arg.lower(), conventions.prologue_handled_arg.lower()}
    if not written:
        return

    taken = taken_over_arguments(published_arguments(impl, conventions), conventions,
                                 prologue)

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


def _kind_value(kind: str | None) -> int | None:
    """A kind's value, from its name (`int32`, `c_double`) or as written (`8`); None if unknown."""
    text = (kind or "").strip().lower()
    if text.isdigit():
        return int(text)
    return KIND_VALUES.get(text)


def _intent(argument: Argument) -> str:
    return argument.intent.value if argument.intent is not None else "unknown"

_READ_FIRST_NOTE = (
    "a prologue runs below the work arrays it prepares; take this value as a 'tmp_' argument "
    "instead, or derive it from something the caller passes"
)

_NON_OPTIONAL_NOTE = (
    "the C wrapper reads it before it may take c_loc of the arrays it describes, so it "
    "has to be there; see the null validation order in the generator README"
)
