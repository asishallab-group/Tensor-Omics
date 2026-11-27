"""Emitting the generated wrapper modules.

One `tox_<name>.F90` per implementation module, holding `foo` and -- where the implementation
has work arrays to take over -- `foo_expert` beside it. Unlike the C wrappers these are
ordinary library sources: no `bind(C)`, no `NO_C_BINDING` guard, no pointer marshalling.

A validating body: set `ierr` ok (the implementation it calls has none and reports nothing),
validate every input against its documented range -- and, for reals, the framework's default
finiteness contract -- bail on error, then call the implementation unchanged. That is the
whole of `foo_expert`, and the whole of a lone `foo` with nothing to take over.

An allocating body does the work the caller should not have to: it drops the implementation's
work arrays, permutations and recommend-sized values from its signature, validates what is
left, calls the recommend routines to size the work arrays, allocates them, seeds and sorts
the permutations, runs the prologue, then calls the implementation directly (it has just
built the permutation, so there is nothing left for the expert tier to re-check).

It reads the implementation's full picture -- every argument, its extents, and its
`DM_OUTPUT_FROM` plan -- off its `foo_expert` sibling, which carries them all.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

from ..config import CONVENTIONS, Conventions
from ..ir.doc import Doc, DocLine
from ..ir.entities import Argument, Module, Procedure
from ..ir.roles import fixed_length, unused_outputs
from ..ir.types import BaseType
from ..ir.validate import is_omp_root
from ..names import DEFAULT_VALUE, UNUSED_OUTPUT, local_name
from ..render import Writer
from ..synthesize import (
    ModeFix,
    is_computed,
    is_permutation,
    is_temporary,
    sorted_permutations,
    taken_over_arguments,
)
from .doc_ford import render_doc
from .fortran_c import _EXTENT_IDENTIFIER_RE, _chunks, _product


#: The exclusive-bound helpers a range expression may wrap a bound in, and where they live.
#: Named by their own implementation module rather than by any parent re-exporting them, so
#: the import says where the thing is defined -- as a resolved module constant already does.
_BOUND_HELPERS = {"above": "f42_math_impl", "below": "f42_math_impl"}

#: What seeding and sorting a permutation needs, and where those live
_PERMUTATION_HELPERS = {
    "init_perm": "f42_sort_impl",
    "sort_array_heapsort": "f42_sort_impl",
}

#: The intrinsic-module functions a range expression may call, and where they come from.
#: A bound is often "how many of these values are usable", which is a NaN test away.
_BOUND_INTRINSICS = {
    "ieee_is_nan": "ieee_arithmetic",
    "ieee_is_finite": "ieee_arithmetic",
}

#: Modules that must be imported `use, intrinsic ::`, not plain `use`
_INTRINSIC_MODULES = frozenset({"ieee_arithmetic"})


@dataclass(frozen=True)
class WrapperInfo:
    """What the emitter needs about a generated wrapper beyond its own signature.

    The implementation it calls (its name is not derivable for a mode-split wrapper, which is
    named from the mode table) and, for a per-mode wrapper, the mode it fixes.
    """

    impl_name: str
    mode_fix: ModeFix | None = None


def wrapper_info(synthesis) -> dict[str, dict[str, WrapperInfo]]:
    """The `WrapperInfo` of every generated wrapper, by generated module, then by wrapper name."""
    info: dict[str, dict[str, WrapperInfo]] = {}
    for spec in synthesis.specs:
        module_info = info.setdefault(spec.module_name, {})
        for wrapper in (spec.validating, spec.allocating):
            if wrapper is not None:
                module_info[wrapper.name.lower()] = WrapperInfo(spec.impl.name, spec.mode_fix)
    return info


class Purity(IntEnum):
    """How pure a generated wrapper can be declared, weakest first, so that `min` picks the weakest.

    `OMP_PURE` is an OpenMP root's `M_OMP_PURE`: pure without OpenMP, impure with it. A wrapper
    is as pure as the least pure thing it calls, so one that calls a root -- its implementation,
    or a `DM_OUTPUT_FROM` producer that is a root's wrapper or calls one -- is `M_OMP_PURE` too:
    a plain `pure` would not compile with OpenMP. Impure beats both.
    """

    IMPURE = 0
    OMP_PURE = 1
    PURE = 2


class FortranWrapperEmitter:
    def __init__(
        self,
        conventions: Conventions = CONVENTIONS,
        macros_header: str = "src/macros.h",
        project=None,
        wrapper_info: dict[str, dict[str, WrapperInfo]] | None = None,
    ):
        self.conventions = conventions
        self.macros_header = macros_header
        #: used to resolve module constants named in a `DM_MIN`/`DM_MAX` expression, so
        #: they can be imported; None in a unit test whose bounds use no constants
        self.project = project
        #: per-wrapper info (impl name, mode fix), keyed by lower-case procedure name;
        #: set for the duration of a `module()` call. Empty in a unit test that passes none,
        #: which drives the ordinary (non-split) path via the fallbacks below.
        self._wrapper_info: dict[str, WrapperInfo] = {}
        #: the same for every generated module (`wrapper_info`), for a wrapper this one calls
        #: in another module -- a producer, whose purity is asked of it in turn
        self._project_wrapper_info = wrapper_info or {}

    # -- module -----------------------------------------------------------------

    def module(self, module: Module, wrapper_info: dict[str, WrapperInfo] | None = None) -> str:
        if not module.procedures and module.uses:
            return self.reexport_module(module)
        self._wrapper_info = wrapper_info or {}
        writer = Writer()
        writer.line(f"#include <{self.macros_header}>")
        writer.blank()

        writer.block(self._module_doc(module))
        writer.line(f"module {module.name}")
        with writer.indent():
            self._module_uses(writer, module)
            # Not a bare `implicit none`: a mistyped impl call would otherwise compile as
            # an implicit external and fail only at link. Generated code is exactly where
            # that must surface at the generator's output, not in a build log.
            writer.line("M_IMPLICIT_NONE")
            writer.line("private")
            writer.blank()
            for procedure in module:
                writer.line(f"public :: {procedure.name}")
        writer.blank()
        writer.line("contains")
        writer.blank()

        for procedure in module:
            with writer.indent():
                writer.block(self.subroutine(procedure, module))
            writer.blank()

        writer.line(f"end module {module.name}")
        return writer.render(trailing_newline=True)

    def reexport_module(self, module: Module) -> str:
        """A parent that gathers a split family, `use`ing its children and nothing else.

        No `only` list and no `private`: re-exporting is the whole point, so whatever a child
        makes public this module makes public too. It declares nothing itself, so it needs
        neither the macros header nor `implicit none`.
        """
        writer = Writer()
        writer.block(self._module_doc(module))
        writer.line(f"module {module.name}")
        with writer.indent():
            for used in module.uses:
                writer.line(f"use {used}")
        writer.line(f"end module {module.name}")
        return writer.render(trailing_newline=True)

    def _module_doc(self, module: Module) -> str:
        """The module's documentation, with a note saying where the file came from.

        The documentation is the implementation module's own, unchanged: a generated module is
        the published API, so the author's prose about the family is what a reader of it wants.
        The note is added here rather than carried in the IR because it is about *this file* --
        Python and R say the same thing in their own words, and would otherwise print it twice.

        No `summary:` tag either. The doc's own first line is the summary; a "Wrappers for
        <impl>" tag above it would displace that everywhere Ford shows a one-liner, to say what
        the note at the bottom of the same comment already says.
        """
        impl = self._impl_module(module)
        note = f"Generated from [[{impl}(module)]]; do not edit -- regenerate instead."
        blocks = list(module.doc.blocks)
        if blocks:
            blocks.append(DocLine.parse(""))
        blocks.append(DocLine.parse(note))
        return render_doc(Doc(tuple(blocks)), kind="module")

    def _impl_module(self, module: Module) -> str:
        return f"{module.name}{self.conventions.impl_suffix}"

    def _module_uses(self, writer: Writer, module: Module) -> None:
        impl_module = self._impl_module(module)

        # everything imported besides the intrinsic kinds and tox_errors, keyed by module:
        # the recommend routines, the permutation/bound helpers, and any module
        # constant a range expression names
        extra: dict[str, set[str]] = {}
        for producer_module, names in self._producers(module).items():
            if producer_module.lower() == module.name.lower():
                continue  # generated into this very module: in scope, and no module uses itself
            extra.setdefault(producer_module, set()).update(names)
        for prologue_module, names in self._prologues(module).items():
            extra.setdefault(prologue_module, set()).update(names)
        if self._has_permutations(module):
            for name, source in _PERMUTATION_HELPERS.items():
                extra.setdefault(source, set()).add(name)
        for helper_module, names in self._bound_imports(module).items():
            extra.setdefault(helper_module, set()).update(names)
        for procedure in module:
            fix = self._mode_fix(procedure)
            if fix is not None:  # the parameter a per-mode wrapper fixes the mode to
                extra.setdefault(fix.module, set()).add(fix.parameter)
            # and the parameters a runtime-mode wrapper checks its argument against
            for argument in procedure.arguments:
                roles = argument.roles
                if roles is None or roles.mode is None or roles.mode.is_split:
                    continue
                for value in roles.mode.values:
                    extra.setdefault(value.module, set()).add(value.parameter)

        kinds = {a.type.kind for p in module for a in p.arguments if a.type.kind}
        # kinds the emitter itself writes, which no argument need carry: a module whose
        # arguments are all real still spells `allow_nan=.true._c_bool`
        for procedure in module:
            kinds |= self._emitted_literal_kinds(procedure, module)
        # and of the locals a producer's dropped outputs go into, which no argument declares
        kinds |= self._unused_output_kinds(module)

        # A module that names a C kind depends on f42_safeguard, whose compile-time guards say
        # which kind is missing where the platform has none. The dependency is what puts
        # safeguard at the front of the build: without it the module compiles first and fails
        # with the compiler's message about the kind instead of the guard's about which one.
        if any(_kind_module(kind) == "iso_c_binding" for kind in kinds):
            writer.line("use f42_safeguard")

        # the implementations, plus any recommend routine / constant that lives in the
        # implementation module
        impl_names = {self._impl_name(p, module) for p in module}
        impl_names |= extra.pop(impl_module, set())
        for chunk in _chunks(sorted(impl_names), 4):
            writer.line(f"use {impl_module}, only: {', '.join(chunk)}")

        by_owner: dict[str, set[str]] = {}
        for kind in kinds:
            by_owner.setdefault(_kind_module(kind), set()).add(kind)
        for owner in sorted(by_owner):
            for chunk in _chunks(sorted(by_owner[owner]), 6):
                writer.line(f"use, intrinsic :: {owner}, only: {', '.join(chunk)}")

        for other_module, names in sorted(extra.items()):
            intrinsic = ", intrinsic :: " if other_module in _INTRINSIC_MODULES else " "
            for chunk in _chunks(sorted(names), 4):
                writer.line(f"use{intrinsic}{other_module}, only: {', '.join(chunk)}")

        for chunk in _chunks(self._error_imports(module), 4):
            writer.line(f"use tox_errors, only: {', '.join(chunk)}")

    def _bound_imports(self, module: Module) -> dict[str, set[str]]:
        """Modules and names a range expression refers to and so must import.

        A bound may be `above(0.0_real64)` (helpers from `f42_math_impl`), call an IEEE inquiry
        (`count(.not. ieee_is_nan(x))`, for a bound that is "how many of these are usable"),
        or name a module constant like `PI`. Identifiers that are the wrapper's own
        arguments, Fortran intrinsics, or literals need nothing.

        A module constant may also size an argument -- `x(N_FIXED, n)`, `len=MAX_LEN` -- and
        the wrapper redeclares every argument, so those are imported too. Before they were,
        such a wrapper compiled only where a range directive happened to name the same
        constant.
        """
        identifiers: set[str] = set()
        sizes: set[str] = set()
        for procedure in module:
            for argument in procedure.arguments:
                directives = argument.directives
                for directive in (directives.minimum, directives.maximum, directives.sentinel,
                                  directives.equals):
                    if directive is not None:
                        identifiers.update(
                            i.lower() for i in _EXTENT_IDENTIFIER_RE.findall(directive.expression)
                        )
                sizes.update(i.lower() for i in _size_identifiers(argument))
        if not identifiers and not sizes:
            return {}

        result: dict[str, set[str]] = {}
        for helper, source in _BOUND_HELPERS.items():
            if helper in identifiers:
                result.setdefault(source, set()).add(helper)
        for name, source in _BOUND_INTRINSICS.items():
            if name in identifiers:
                result.setdefault(source, set()).add(name)

        dummies = {a.name.lower() for p in module for a in p.arguments}
        if self.project is not None:
            parameters = {p.name.lower(): p for p in self.project.parameters}
            for identifier in (identifiers | sizes) - dummies:
                parameter = parameters.get(identifier)
                if parameter is not None and parameter.parent is not None:
                    result.setdefault(parameter.parent.name, set()).add(parameter.name)
        return result

    def _error_imports(self, module: Module) -> list[str]:
        names = {"set_ok", "is_err"}
        names |= self._position_imports(module)
        for procedure in module:
            for argument in procedure.arguments:
                validator = self._validator_for(argument)
                if validator is not None:
                    names.add(validator)
        if any(self._is_allocating(p, module) for p in module):
            # M_ALLOCATE reports through set_err / ERR_ALLOC_FAIL
            names |= {"set_err", "ERR_ALLOC_FAIL"}
        if any(self._is_distance_matrix(a) for p in module for a in p.arguments):
            names.add("validate_distance_matrix")
        # both the mask-count convention and the mode-membership check report this way
        if any(
            (a.roles is not None and a.roles.mask_count_of is not None)
            or self._mode_membership_check(a, 1) is not None
            for p in module
            for a in p.arguments
        ):
            names |= {"set_err_once", "ERR_INVALID_INPUT"}
        # a `DM_EQUALS` comparison, which no validator call names
        if any(a.directives.equals is not None for p in module for a in p.arguments):
            names |= {"set_err_once", "ERR_SIZE_MISMATCH"}
        ordered = ["set_ok", "is_err"]
        return ordered + sorted(names - set(ordered))

    def _prologues(self, module: Module) -> dict[str, set[str]]:
        """The prologue procedures this module's wrappers call, by the module they live in.

        Only the allocating wrappers have one, so only they are asked.
        """
        found: dict[str, set[str]] = {}
        for procedure in module:
            if not self._is_allocating(procedure, module):
                continue
            prologue = self._prologue_for(self._sibling(module, procedure), module)
            if prologue is not None and prologue.module is not None:
                found.setdefault(prologue.module.name, set()).add(prologue.name)
        return found

    def _producers(self, module: Module) -> dict[str, set[str]]:
        producers: dict[str, set[str]] = {}
        for _, foo, group in self._producer_calls(module):
            plan = group[0].roles.computed_from
            producers.setdefault(plan.producer.module.name, set()).add(plan.producer.name)
        return producers

    def _producer_calls(self, module: Module):
        """Every producer call an allocating wrapper here makes: (alloc, foo, group).

        `group` is the expert sibling's arguments one call fills, in their order -- each
        producer is called once per distinct set of inputs, for all of the outputs a wrapper
        takes over from it with those inputs (`_producer_groups`).
        """
        for alloc in module:
            if not self._is_allocating(alloc, module):
                continue
            foo = self._sibling(module, alloc)
            taken = taken_over_arguments(
                self._impl_arguments(foo), self.conventions, self._prologue_for(foo, module),
                self._evaluator(),
            )
            for group in self._producer_groups(foo, taken):
                yield alloc, foo, group

    @staticmethod
    def _producer_groups(foo: Procedure, taken) -> list[list[Argument]]:
        """The taken-over AUTO arguments of `foo`, grouped by the producer call that fills them.

        One call per producer and set of inputs (`OutputFromPlan.call`): a producer whose inputs
        two arguments' tables supply differently is called once for each. In argument order: no
        producer is given another's output (`roles._refuses_chaining`), so the order is free.
        """
        groups: dict[tuple, list[Argument]] = {}
        for argument in taken:
            if not is_computed(argument):
                continue
            own = foo.argument(argument.name)
            groups.setdefault(own.roles.computed_from.call, []).append(own)
        return list(groups.values())

    def _error_returns(self, module: Module):
        """The position pairs of every `_error_return` a body here emits, one list per call.

        After each producer, prologue and implementation call that takes `ierr`, in the order
        the bodies make them. The one source of truth for both the emitted blocks and the `use`
        list -- deciding the two separately is how a missing import gets shipped.
        """
        for procedure in module:
            foo = procedure
            if self._is_allocating(procedure, module):
                foo = self._sibling(module, procedure)
                prologue = self._prologue_for(foo, module)
                taken = taken_over_arguments(
                    self._impl_arguments(foo), self.conventions, prologue, self._evaluator()
                )
                for group in self._producer_groups(foo, taken):
                    if group[0].roles.computed_from.producer.has_error_argument:
                        yield self._producer_positions(group, procedure)
                if prologue is not None:
                    actuals = self._prologue_actuals(prologue, procedure, taken)
                    yield self._call_positions(prologue.arguments, actuals, procedure)
            if self._impl_declares_error(foo, module):
                yield self._impl_positions(foo, module, procedure)

    def _position_imports(self, module: Module) -> set[str]:
        """What the blocks `_error_return` writes here name from `tox_errors`.

        `clear_err_arg_pos` for any block, `get_err_arg_pos` once a position maps at all, and
        `set_err_arg_pos` only where one moves.
        """
        names: set[str] = set()
        for positions in self._error_returns(module):
            names.add("clear_err_arg_pos")
            if positions:
                names.add("get_err_arg_pos")
            if any(old != new for old, new in positions):
                names.add("set_err_arg_pos")
        return names

    def _has_permutations(self, module: Module) -> bool:
        """Whether anything here is seeded and sorted, so `f42_sort_impl` is worth importing.

        Asks the same question the emission does, prologue included -- a family whose only
        permutation the prologue builds needs neither helper.
        """
        for alloc in module:
            if not self._is_allocating(alloc, module):
                continue
            foo = self._sibling(module, alloc)
            prologue = self._prologue_for(foo, module)
            taken = taken_over_arguments(
                self._impl_arguments(foo), self.conventions, prologue, self._evaluator()
            )
            if sorted_permutations(taken, prologue, self.conventions):
                return True
        return False

    # -- subroutine -------------------------------------------------------------

    def _is_allocating(self, procedure: Procedure, module: Module) -> bool:
        """Whether `procedure` is emitted with the allocating body.

        Exactly the wrappers that have an expert sibling to take the work arrays over from.
        A procedure with no `<name>_expert` beside it is the only entry point there is, so
        it takes nothing over and is emitted with the ordinary validating body -- which is
        also why a lone wrapper wears the plain name (`synthesize._wrappers_for`).
        """
        return self._sibling(module, procedure) is not None

    def subroutine(self, procedure: Procedure, module: Module) -> str:
        writer = Writer()
        impl = self._impl_name(procedure, module)
        link = f"[[{self._impl_module(module)}(module):{impl}]]"
        # Which tier this is, not merely what it does: both validate, and a reader of the
        # generated manual has to be able to tell `foo` from `foo_expert` by its summary.
        if self._is_allocating(procedure, module):
            summary = (
                f"Validates its inputs, prepares what {link} needs, then calls it. "
                f"The entry point to reach for first; see "
                f"[[{module.name}(module):{procedure.name}{self.conventions.expert_suffix}]] "
                "to prepare it yourself."
            )
        elif self._plain_of(module, procedure) is not None:
            summary = (
                f"Validates its inputs, then calls {link} with what you supply. "
                "The expert entry point: it allocates nothing and prepares nothing; "
                f"[[{module.name}(module):{self._base(procedure)}]] does both."
            )
        else:
            summary = f"Validates its inputs, then calls {link}."
        writer.block(render_doc(procedure.doc, kind="procedure", summary=summary))
        self._signature(writer, procedure, module)
        with writer.indent():
            if self._is_allocating(procedure, module):
                self._allocating_body(writer, procedure, module)
            else:
                self._validating_body(writer, procedure, module)
        writer.line(f"end subroutine {procedure.name}")
        return writer.render()

    def _validating_body(
        self, writer: Writer, foo: Procedure, module: Module
    ) -> None:
        self._declarations(writer, foo.arguments)
        writer.blank()
        self._validation(writer, foo)
        writer.blank()
        self._impl_call(writer, foo, module, foo)

    def _allocating_body(
        self, writer: Writer, alloc: Procedure, module: Module
    ) -> None:
        foo = self._sibling(module, alloc)
        prologue = self._prologue_for(foo, module)
        taken = taken_over_arguments(
            self._impl_arguments(foo), self.conventions, prologue, self._evaluator()
        )
        materialized = self._materialized_producer_inputs(foo, alloc, taken)
        self._declarations(writer, alloc.arguments)
        self._locals(writer, taken)
        self._unused_output_locals(writer, foo, taken)
        self._materialized_locals(writer, materialized)
        self._prologue_local(writer, prologue)
        writer.blank()
        # rule V: a check runs once every value it names is set -- at the top (P0), after the
        # producers (P1), or after the prologue (P2)
        late = self._validation(writer, alloc, self._late_values(taken, prologue))
        writer.blank()
        self._materialized_defaults(writer, materialized)
        self._recommend_calls(writer, foo, alloc, taken, materialized)
        self._late_checks(writer, late[1])
        self._allocations(writer, taken)
        self._permutations(writer, taken, prologue)
        writer.blank()
        # Below the work arrays, because preparing them is what a prologue is for -- it is the
        # sugar this tier adds over the expert one, which lets a caller prepare them instead.
        self._prologue_call(writer, prologue, alloc, module, taken)
        self._late_checks(writer, late[2])
        # mapped to `alloc`'s own dummy list: `foo` is the expert sibling, whose work arrays
        # and recommend-sized values are locals here
        self._impl_call(writer, foo, module, alloc)

    def _evaluator(self):
        """The project's one constant evaluator, which keys a producer call as synthesis did
        (rule K); None in a unit test without a project."""
        return self.project.evaluator() if self.project is not None else None

    # -- prologue ---------------------------------------------------------------

    def _prologue_for(self, foo: Procedure, module: Module):
        """The prologue the allocating wrapper runs, resolved to the procedure it names.

        Declared on the implementation, since it is the implementation's own contract. Returns
        None when there is none, or when the project is unavailable (a unit test), so the
        ordinary path is unaffected.
        """
        impl = self._impl(foo, module)
        if impl is None:
            return None
        directive = impl.directives.prologue
        if directive is None:
            return None
        return self.project.procedure(directive.module, directive.procedure)

    def _prologue_local(self, writer: Writer, prologue) -> None:
        if prologue is not None:
            writer.line(f"logical :: {self.conventions.prologue_handled_arg}")

    def _prologue_call(self, writer: Writer, prologue, wrapper: Procedure, module: Module,
                       taken=()) -> None:
        """Run the prologue, and return early if it handled the call itself.

        Its dummies are supplied by name -- from the wrapper's own arguments, or from the
        work arrays the allocating wrapper took over and prepared as locals; `handled` and
        `ierr` come from the wrapper. A dummy that matches none of these is refused by
        `validate._check_prologue`, so nothing is silently left out of the call here.
        """
        if prologue is None:
            return
        actuals = self._prologue_actuals(prologue, wrapper, taken)
        self._call(writer, prologue.name, actuals)
        self._error_return(writer, self._call_positions(prologue.arguments, actuals, wrapper))
        writer.line(f"if ({self.conventions.prologue_handled_arg}) return")
        writer.blank()

    def _prologue_actuals(self, prologue, wrapper: Procedure, taken) -> list[tuple[str, str]]:
        """The prologue call's `(dummy, actual)` pairs, as `_prologue_call` makes it."""
        error = self.conventions.error_arg
        handled = self.conventions.prologue_handled_arg
        local = {argument.name.lower() for argument in taken}
        actuals = []
        for dummy in prologue.arguments:
            lowered = dummy.name.lower()
            if lowered == error.lower():
                actuals.append((dummy.name, error))
            elif lowered == handled.lower():
                actuals.append((dummy.name, handled))
            elif wrapper.argument(dummy.name) is not None or lowered in local:
                actuals.append((dummy.name, dummy.name))
        return actuals

    def _purity(self, procedure: Procedure, module: Module,
                visiting: frozenset = frozenset()) -> Purity:
        """How pure this wrapper can be: as pure as the least pure thing it calls.

        A caller who wants a `pure` procedure -- to call it from their own pure code, which
        may call only pure procedures -- should not have to reach past the wrapper to the
        implementation to get one. The wrapper is pure exactly when the implementation is
        and so is everything the wrapper adds, so the property is derived rather than
        declared: `tox_errors`' validators and `set_ok`/`is_err` are pure, and so are
        `init_perm` and every `sort_array_heapsort` specific, which leaves the
        implementation itself, the `DM_OUTPUT_FROM` producers and the prologue.

        A producer that is itself a generated wrapper is pure by this same rule, so it is
        asked in turn (`_purity_of_call`); `visiting` holds the wrappers already being asked, so
        a chain that comes back to one of them ends instead of recursing forever.

        An OpenMP root is pure only without OpenMP (`Purity.OMP_PURE`), and so is everything that
        calls one: a root's own wrappers, and a consumer whose producer is a root's wrapper or
        reaches one, however far down. Such a consumer is `M_OMP_PURE` even where all else it
        calls is plain `pure`, since `pure` would not compile with OpenMP.

        (`allocate(..., stat=)` is permitted in a pure procedure, so the allocating tier is
        no less eligible than the expert one.)

        False without a project to resolve them -- a unit test. Guessing `pure` wrongly is a
        compile error in generated code; guessing it away only costs the caller.
        """
        if self.project is None:
            return Purity.IMPURE
        impl = self._impl(procedure, module)
        if impl is None:
            return Purity.IMPURE
        visiting = visiting | {self._key(procedure, module)}
        called = [impl]
        foo = self._sibling(module, procedure) or procedure
        called.append(self._prologue_for(foo, module))
        for argument in foo.arguments:
            roles = argument.roles
            plan = roles.computed_from if roles is not None else None
            if plan is not None:
                called.append(plan.producer)
        return min(self._purity_of_call(p, visiting) for p in called if p is not None)

    def _purity_of_call(self, called: Procedure, visiting: frozenset) -> Purity:
        """How pure a procedure the wrapper calls is: as declared, or as generated.

        Declared pure, it is plain pure unless it is an OpenMP root, which is declared
        `M_OMP_PURE` (`validate._check_omp_root_purity`). A generated wrapper is built impure
        (`synthesize`) and emitted pure only by `_purity`, so its declaration says nothing; it is
        asked the same question instead. One already being asked is a cycle, which nothing
        proves pure.
        """
        if called.is_pure:
            return Purity.OMP_PURE if is_omp_root(called, self.conventions) else Purity.PURE
        owner = called.module
        if owner is None or self.project.module(self._impl_module(owner)) is None:
            return Purity.IMPURE  # hand-written, and not declared pure
        if self._key(called, owner) in visiting:
            return Purity.IMPURE
        return self._purity(called, owner, visiting)

    @staticmethod
    def _key(procedure: Procedure, module: Module) -> tuple[str, str]:
        return module.name.lower(), procedure.name.lower()

    # -- names ------------------------------------------------------------------

    def _base(self, procedure: Procedure) -> str:
        """The name both tiers of a pair share: `foo` for `foo` and for `foo_expert`.

        Only the expert suffix is stripped, because it is the only one the generator adds.
        `validate` refuses an implementation whose base name already ends in it, so no
        wrapper reaches here whose own name would be shortened by mistake.
        """
        name = procedure.name
        suffix = self.conventions.expert_suffix
        return name[: -len(suffix)] if name.lower().endswith(suffix) else name

    def _impl_name(self, procedure: Procedure, module: Module) -> str:
        # a mode-split wrapper is named from the mode table, so its implementation is not
        # derivable from its name; take it from the spec when available
        info = self._info(procedure, module)
        if info is not None:
            return info.impl_name
        return f"{self._base(procedure)}{self.conventions.impl_suffix}"

    def _info(self, procedure: Procedure, module: Module) -> WrapperInfo | None:
        """The info of a wrapper in `module`: the one being emitted, or any other generated one."""
        info = self._project_wrapper_info.get(module.name, self._wrapper_info)
        return info.get(procedure.name.lower())

    def _mode_fix(self, procedure: Procedure) -> ModeFix | None:
        info = self._wrapper_info.get(procedure.name.lower())
        return info.mode_fix if info is not None else None

    def _sibling(self, module: Module, allocating: Procedure) -> Procedure | None:
        """The expert half of `allocating`, which carries the full signature to read.

        `foo_expert` is where the implementation's whole argument list survives -- the work
        arrays, the permutations and the recommend-sized values that `foo` drops -- so the
        allocating body reads its picture from there rather than from its own dummies.
        """
        return module.procedure(allocating.name + self.conventions.expert_suffix)

    def _plain_of(self, module: Module, expert: Procedure) -> Procedure | None:
        """The allocating half `expert` is the expert of, or None if it is a lone wrapper."""
        base = self._base(expert)
        return module.procedure(base) if base != expert.name else None

    def _impl_arguments(self, foo: Procedure) -> list[Argument]:
        """The implementation's own arguments minus ierr -- the ones that may be taken over.

        Used to work out the allocations and permutations; ierr is never among them, so
        dropping it here is harmless and independent of whether the implementation declares
        one.
        """
        error = self.conventions.error_arg.lower()
        return [a for a in foo.arguments if a.name.lower() != error]

    def _impl_declares_error(self, foo: Procedure, module: Module) -> bool:
        """Whether the implementation takes `ierr`, so the wrapper has a position to translate."""
        error = self.conventions.error_arg.lower()
        return any(
            a.name.lower() == error for a in self._impl_call_arguments(foo, module)
        )

    def _impl_call_arguments(self, foo: Procedure, module: Module) -> list[Argument]:
        """The implementation's own signature, for the call.

        Includes `ierr` exactly when the implementation declares one: an implementation that
        propagates a sub-helper's failure takes `ierr`, a pure one does not, and the wrapper
        must pass on whichever the implementation actually has. Resolved from the
        implementation itself when the project is available; without it (a unit test) the
        implementation has no ierr, since synthesis only appends one to the wrapper.
        """
        impl = self._impl(foo, module)
        if impl is not None:
            return list(impl.arguments)
        return self._impl_arguments(foo)

    def _impl(self, wrapper: Procedure, module: Module) -> Procedure | None:
        """The implementation `wrapper` calls, or None without a project (a unit test)."""
        if self.project is None:
            return None
        return self.project.procedure(self._impl_module(module), self._impl_name(wrapper, module))

    # -- declarations -----------------------------------------------------------

    def _signature(self, writer: Writer, procedure: Procedure, module: Module) -> None:
        # `M_OMP_PURE` says `pure` or `impure` per build, for a wrapper that calls an OpenMP
        # root, whose purity is the build's too: a plain `pure` would not compile with OpenMP
        prefix = {
            Purity.IMPURE: "",
            Purity.OMP_PURE: f"{self.conventions.omp_pure_macro} ",
            Purity.PURE: "pure ",
        }[self._purity(procedure, module)]
        writer.line(f"{prefix}subroutine {procedure.name}(&")
        names = [a.name for a in procedure.arguments]
        with writer.indent(2):
            for name in names[:-1]:
                writer.line(f"{name},&")
            writer.line(f"{names[-1]}&")
        with writer.indent():
            writer.line(")")

    def _declarations(self, writer: Writer, arguments) -> None:
        for argument in self._declaration_order(arguments):
            writer.line(self._declaration(argument))
            if argument.doc:
                with writer.indent():
                    writer.block(render_doc(argument.doc, kind="argument"))

    @staticmethod
    def _declaration_order(arguments) -> list[Argument]:
        """An extent is declared before the array whose extent it is (gfortran -std=f2018
        rejects the reverse). Everything else keeps the author's order."""
        referenced = {
            identifier.lower()
            for argument in arguments
            for extent in argument.dimension.extents
            for identifier in _EXTENT_IDENTIFIER_RE.findall(extent)
        }
        extents = [a for a in arguments if a.name.lower() in referenced]
        rest = [a for a in arguments if a.name.lower() not in referenced]
        return [*extents, *rest]

    def _declaration(self, argument: Argument) -> str:
        attributes = []
        if argument.dimension:
            attributes.append(f"dimension({', '.join(argument.dimension.extents)})")
        attributes.append(f"intent({argument.intent.value})")
        if argument.optional:
            attributes.append("optional")
        return f"{argument.type}, {', '.join(attributes)} :: {argument.name}"

    def _locals(self, writer: Writer, taken) -> None:
        """The taken-over arguments, as locals the allocating wrapper prepares itself."""
        for argument in taken:
            if argument.is_array:
                deferred = ", ".join(":" for _ in argument.dimension.extents)
                writer.line(
                    f"{argument.type}, dimension({deferred}), allocatable :: {argument.name}"
                )
            else:
                writer.line(f"{argument.type} :: {argument.name}")

    # -- validation -------------------------------------------------------------

    def _late_values(self, taken, prologue) -> tuple[set[str], set[str]]:
        """What the allocating `foo` sets after the top, lower-cased (rule V): the AUTO values it
        computes (P1, after the producers) and what the prologue writes (P2, after it)."""
        produced = {argument.name.lower() for argument in taken if is_computed(argument)}
        reserved = {self.conventions.error_arg.lower(),
                    self.conventions.prologue_handled_arg.lower()}
        written = ({dummy.name.lower() for dummy in prologue.arguments
                    if dummy.intent is not None and dummy.intent.is_output} - reserved
                   if prologue is not None else set())
        return produced, written

    @staticmethod
    def _point(expression: str, late: tuple[set[str], set[str]] | None) -> int:
        """Where a check naming `expression` runs: the latest point any of its names needs."""
        if late is None:
            return 0
        from ..synthesize import size_names_in  # local, as the other synthesize helpers here

        names = {named.name.lower() for named in size_names_in(expression)}
        produced, written = late
        return 2 if names & written else 1 if names & produced else 0

    def _deferred(self, argument: Argument, late) -> dict[str, int]:
        """The `min`/`max` keywords of `argument`'s validator call that wait for P1 or P2."""
        directives = argument.directives
        found = {}
        for keyword, directive in (("min", directives.minimum), ("max", directives.maximum)):
            if directive is not None:
                point = self._point(directive.expression, late)
                if point:
                    found[keyword] = point
        return found

    def _late_checks(self, writer: Writer, lines: list[str]) -> None:
        """A point's deferred checks, in a validation block of their own (rule V)."""
        if not lines:
            return
        writer.directive("#ifndef NO_INPUT_VALIDATION")
        for line in lines:
            writer.line(line)
        writer.line(f"if (is_err({self.conventions.error_arg})) return")
        writer.directive("#endif")
        writer.blank()

    def _validation(self, writer: Writer, procedure: Procedure,
                    late: tuple[set[str], set[str]] | None = None) -> dict[int, list[str]]:
        """Check every input, then bail -- all of it behind `NO_INPUT_VALIDATION`.

        `set_ok` stays outside the guard. It is not a check: it is what leaves `ierr` defined
        on the path where nothing goes wrong, and a build without validation still has to
        report the runtime errors an implementation raises.

        Compiling the checks out is for a caller who has already established that the inputs
        are good -- an inner loop calling the same wrapper a million times over data it
        produced itself. It buys the checks' cost and gives up every diagnostic the framework
        offers, which is why it is a whole-build switch rather than a per-call argument: it is
        not a decision to take one call site at a time.

        In the allocating `foo` (`late` given), a `DM_MIN` or `DM_MAX` naming a value set later
        is split off its argument's call (rule V): the call stays here with every other keyword
        -- finiteness, constant bounds, bounds naming what the caller passes -- so nothing reads
        the argument unchecked, and the deferred keywords run in a second call at their point
        (P1, P2), with the sentinel and, for a real, finiteness allowed, which was checked here
        already. A `DM_EQUALS` moves whole. Returned, by point, for the caller to emit there.
        """
        error = self.conventions.error_arg
        writer.line(f"call set_ok({error})")
        later: dict[int, dict[str, list[str]]] = {1: {"scalars": [], "arrays": []},
                                                  2: {"scalars": [], "arrays": []}}
        # scalars (extents, bounds) before arrays, as the hand-written wrappers do; arg_pos
        # stays the argument's true position, so the ordering is presentational only
        calls = []
        for position, argument in enumerate(procedure.arguments, start=1):
            if argument.name.lower() == error.lower():
                continue
            deferred = self._deferred(argument, late)
            calls.append((argument.is_array, self._validation_call(argument, position,
                                                                   frozenset(deferred))))
            for point in sorted(set(deferred.values())):
                keys = frozenset(k for k, at in deferred.items() if at == point)
                later[point]["arrays" if argument.is_array else "scalars"].append(
                    self._deferred_call(argument, position, keys))
        scalars = [c for array, c in calls if not array and c]
        # after every scalar check, so an empty dimension is still ERR_EMPTY_INPUT rather than
        # a mismatch of its empty sum; before the arrays, which may be sized by the argument
        for point, check in self._equality_checks(procedure, late):
            (later[point]["scalars"] if point else scalars).append(check)
        emitted = scalars + [c for array, c in calls if array and c]
        emitted += self._convention_checks(procedure)
        deferred_lines = {point: lines["scalars"] + lines["arrays"]
                          for point, lines in later.items()}
        if not emitted:
            return deferred_lines  # nothing to guard, and no `is_err` to test either
        writer.directive("#ifndef NO_INPUT_VALIDATION")
        for line in emitted:
            writer.line(line)
        writer.line(f"if (is_err({error})) return")
        writer.directive("#endif")
        return deferred_lines

    def _equality_checks(self, procedure: Procedure, late=None) -> list[tuple[int, str]]:
        """`DM_EQUALS(EXPR)`: the argument must equal `EXPR`, else `ERR_SIZE_MISMATCH` at it --
        each with the point it runs at (rule V), the top unless `EXPR` names a later value."""
        error = self.conventions.error_arg
        return [
            (self._point(argument.directives.equals.expression, late),
             f"if ({argument.name} /= ({argument.directives.equals.expression})) "
             f"call set_err_once({error}, ERR_SIZE_MISMATCH, arg_pos={position}_int32)")
            for position, argument in enumerate(procedure.arguments, start=1)
            if argument.directives.equals is not None
        ]

    def _convention_checks(self, procedure: Procedure) -> list[str]:
        """Checks driven by a naming convention rather than a per-argument range.

        A distance-matrix argument is validated for structure; a `n_selected_<x>` count is
        validated against `count(<x>_mask)`. Both come from information the wrapper already has
        (the argument's roles / shape), so the implementation need not carry them or an `ierr`.
        """
        error = self.conventions.error_arg
        lines: list[str] = []
        for position, argument in enumerate(procedure.arguments, start=1):
            if self._is_distance_matrix(argument):
                extent = argument.dimension.extents[0]
                lines.append(
                    f"call validate_distance_matrix({argument.name}, {extent}, {error}, "
                    f"arg_pos={position}_int32)"
                )
        for position, argument in enumerate(procedure.arguments, start=1):
            roles = argument.roles
            if roles is None or roles.mask_count_of is None:
                continue
            mask = roles.mask_count_of
            lines.append(
                f"if (count({mask.name}, kind=int32) /= {argument.name}) "
                f"call set_err_once({error}, ERR_INVALID_INPUT, arg_pos={position}_int32)"
            )
        for position, argument in enumerate(procedure.arguments, start=1):
            check = self._mode_membership_check(argument, position)
            if check is not None:
                lines.append(check)
        return lines

    def _mode_membership_check(self, argument: Argument, position: int) -> str | None:
        """Reject a mode value that is not one of the ones its table names.

        The generator parsed that table, so it knows the accepted parameters exactly -- which
        makes this a better answer than a hand-written `DM_MIN`/`DM_MAX` pair: the accepted
        set need not be contiguous, and a range annotation goes stale the moment a mode is
        added, while this is regenerated from the table.

        A chain of `/=` rather than `all(x /= [...])`: the array constructor would build a
        temporary on every call, for a check that is a handful of integer comparisons.

        An optional mode is guarded by `present()`. Absent means "use the documented default",
        which is by construction one of the accepted values -- and reading an absent optional
        to find that out is not a Fortran program.

        Only for a mode the wrapper still takes. A mode-split wrapper *is* its mode and has no
        such argument, and the C layer rejects an unknown mode string before Fortran is
        entered -- but a direct Fortran caller reaches this, and nothing else stops it.
        """
        roles = argument.roles
        if roles is None or roles.mode is None or roles.mode.is_split:
            return None
        values = roles.mode.values
        if not values:
            return None
        error = self.conventions.error_arg
        comparisons = " .and. ".join(
            f"{argument.name} /= {value.parameter}" for value in values
        )
        check = (
            f"if ({comparisons}) "
            f"call set_err_once({error}, ERR_INVALID_INPUT, arg_pos={position}_int32)"
        )
        if argument.optional:
            return f"if (present({argument.name})) then; {check}; end if"
        return check

    def _is_distance_matrix(self, argument: Argument) -> bool:
        """A square real matrix named by the distance-matrix convention.

        It is validated for distance-matrix structure (symmetry, non-negativity, zero diagonal)
        by `validate_distance_matrix`, which no per-argument range validator expresses, so it
        also opts out of the finiteness contract.
        """
        if argument.type.base is not BaseType.REAL or not argument.intent.is_input:
            return False
        extents = argument.dimension.extents
        if len(extents) != 2 or extents[0] != extents[1]:
            return False
        name = argument.name.lower()
        return any(
            name == suffix or name.endswith(f"_{suffix}")
            for suffix in self.conventions.distance_matrix_suffixes
        )

    def _validator_for(self, argument: Argument) -> str | None:
        """The `tox_errors` validator this argument needs, or None."""
        if not argument.intent.is_input:
            return None  # an output carries no value to check
        if self._is_distance_matrix(argument):
            return None  # validated by validate_distance_matrix, emitted as a convention check
        if (
            argument.roles is not None
            and argument.roles.is_extent
            and not argument.directives.has_range
        ):
            # an extent is validated as a dimension by default, but a documented range wins:
            # it lets an extent that is really a count permit zero (an empty selection), which
            # validate_dimension_size rejects
            return "validate_dimension_size"
        base = argument.type.base
        if base is BaseType.REAL:
            # finiteness is the default contract, so a real input is always checked --
            # unless it opts out of both failure modes and states no bound
            directives = argument.directives
            if (
                not directives.has_range
                and directives.allows_nan
                and directives.allows_infinite
            ):
                return None
            return "validate_all_in_range_real" if argument.is_array else "validate_in_range_real"
        if base is BaseType.INTEGER:
            if not argument.directives.has_range:
                return None  # an integer with no bound has nothing to check
            return "validate_all_in_range_int" if argument.is_array else "validate_in_range_int"
        return None

    def _emitted_literal_kinds(self, procedure: Procedure, module: Module | None = None) -> set[str]:
        """Kinds the emitter writes as literals, which no argument need carry.

        The import list is built from the arguments, so a kind that only ever appears in
        something the emitter invents is invisible to it -- and an unimported kind is a
        compile error in code nobody wrote. Three sources, and this has to stay in step with
        them: the `allow_nan=`/`allow_infinite=` keywords below, the ones a deferred check on a
        real carries (rule V), and a `DM_OUTPUT_FROM` producer input supplied as a constant.
        """
        kinds: set[str] = set()
        late = self._late_of(procedure, module) if module is not None else None
        for argument in procedure.arguments:
            if (argument.type.base is BaseType.REAL and self._validator_for(argument) is not None
                    and self._deferred(argument, late)):
                kinds.add(_VALIDATOR_LOGICAL_KIND)
            directives = argument.directives
            if (argument.type.base is BaseType.REAL
                    and self._validator_for(argument) is not None
                    and (directives.allows_nan or directives.allows_infinite)):
                kinds.add(_VALIDATOR_LOGICAL_KIND)
            plan = argument.roles.computed_from if argument.roles else None
            if plan is None:
                continue
            supplied = {pi.name.lower(): pi for pi in plan.inputs}
            for dummy in plan.producer.arguments:
                supply = supplied.get(dummy.name.lower())
                if supply is None or supply.argument is not None:
                    continue
                # rendered in the input's kind (`_fortran_literal`), whatever its type; a
                # string carries none
                if dummy.type.kind and not isinstance(supply.constant, str):
                    kinds.add(dummy.type.kind)
        return kinds

    def _late_of(self, procedure: Procedure, module: Module):
        """`_late_values` of `procedure` where it is an allocating `foo`, else None."""
        if not self._is_allocating(procedure, module):
            return None
        foo = self._sibling(module, procedure)
        prologue = self._prologue_for(foo, module)
        taken = taken_over_arguments(self._impl_arguments(foo), self.conventions, prologue,
                                     self._evaluator())
        return self._late_values(taken, prologue)

    def _deferred_call(self, argument: Argument, position: int, keys: frozenset) -> str:
        """The second validator call of a check deferred to P1 or P2 (rule V): the deferred
        bounds, the sentinel, and for a real finiteness allowed -- the top already checked it."""
        validator = self._validator_for(argument)
        directives = argument.directives
        keywords = [f"arg_pos={position}_int32"]
        if "min" in keys:
            keywords.append(f"min={directives.minimum.expression}")
        if "max" in keys:
            keywords.append(f"max={directives.maximum.expression}")
        if directives.sentinel is not None:
            keywords.append(f"sentinel={directives.sentinel.expression}")
        if argument.type.base is BaseType.REAL:
            keywords.append(f"allow_nan=.true._{_VALIDATOR_LOGICAL_KIND}")
            keywords.append(f"allow_infinite=.true._{_VALIDATOR_LOGICAL_KIND}")
        return f"call {validator}({self._validation_head(argument)}, {', '.join(keywords)})"

    def _validation_head(self, argument: Argument) -> str:
        """A validator call's leading arguments: the value, an array's element count, `ierr`."""
        error = self.conventions.error_arg
        if argument.is_array:
            # an assumed-shape axis has no extent to name, so its size is asked of the array
            count = _product(
                f"size({argument.name}, {axis})" if extent == ":" else extent
                for axis, extent in enumerate(argument.dimension.extents, start=1)
            )
            return f"{argument.name}, {count}, {error}"
        return f"{argument.name}, {error}"

    def _validation_call(self, argument: Argument, position: int,
                         deferred: frozenset = frozenset()) -> str:
        validator = self._validator_for(argument)
        if validator is None:
            return ""
        error = self.conventions.error_arg
        arg_pos = f"arg_pos={position}_int32"

        if validator == "validate_dimension_size":
            return f"call validate_dimension_size({argument.name}, {error}, {arg_pos})"

        head = self._validation_head(argument)
        keywords = [arg_pos]
        directives = argument.directives
        if directives.minimum is not None and "min" not in deferred:
            keywords.append(f"min={directives.minimum.expression}")
        if directives.maximum is not None and "max" not in deferred:
            keywords.append(f"max={directives.maximum.expression}")
        unbounded = not any(k.startswith(("min=", "max=")) for k in keywords)
        unchecked = argument.type.base is BaseType.INTEGER or (
            argument.type.base is BaseType.REAL and directives.allows_nan
            and directives.allows_infinite)
        if deferred and unbounded and unchecked:
            return ""  # every check it has waits for a later point: nothing to check here
        if directives.sentinel is not None:
            keywords.append(f"sentinel={directives.sentinel.expression}")
        if argument.type.base is BaseType.REAL:
            # `tox_errors` declares these two `logical(c_bool)`, and argument association
            # does not convert kinds -- so the emitter's own literal has to carry it.
            if directives.allows_nan:
                keywords.append(f"allow_nan=.true._{_VALIDATOR_LOGICAL_KIND}")
            if directives.allows_infinite:
                keywords.append(f"allow_infinite=.true._{_VALIDATOR_LOGICAL_KIND}")

        return f"call {validator}({head}, {', '.join(keywords)})"

    # -- allocation, recommend routines, permutations ---------------------------

    def _materialized_producer_inputs(
        self, foo: Procedure, alloc: Procedure, taken
    ) -> dict[str, Argument]:
        """Wrapper arguments a recommend routine needs by value, keyed by name.

        A recommend routine is called directly, so an input the wrapper takes optionally with a
        `DM_DEFAULT` is resolved into a local first -- the argument when the caller gave one,
        the documented default when not -- and the routine is handed that (P1). Whether or not
        the routine's own dummy is optional: an absent optional passed to a mandatory one is not
        a Fortran program, and passed to an optional one it applied the *producer's* default,
        while the implementation applies the consumer's -- a work array sized for a stride of
        three, used with one of four. Python and R pass the consumer's published default too. An
        optional without a default is passed through, so absent stays absent.
        """
        materialized: dict[str, Argument] = {}
        for argument in taken:
            if not is_computed(argument):
                continue
            roles = foo.argument(argument.name).roles
            plan = roles.computed_from if roles is not None else None
            if plan is None:
                continue
            for supply in plan.inputs:
                dummy = plan.producer.argument(supply.name)
                if supply.argument is None or dummy is None:
                    continue
                # what supplies it, which a table may name differently from the producer
                supplied = alloc.argument(supply.argument)
                if supplied is None or not supplied.optional:
                    continue
                if supplied.directives.default is None:
                    continue
                materialized[supplied.name] = supplied
        return materialized

    @staticmethod
    def _materialized_name(name: str) -> str:
        return local_name(name, DEFAULT_VALUE)

    def _materialized_locals(self, writer: Writer, materialized: dict[str, Argument]) -> None:
        for name, argument in materialized.items():
            writer.line(f"{argument.type} :: {self._materialized_name(name)}")

    def _materialized_defaults(
        self, writer: Writer, materialized: dict[str, Argument]
    ) -> None:
        if not materialized:
            return
        for name, argument in materialized.items():
            default = argument.directives.default.expression
            writer.line(
                f"M_DEFAULT_VAL({name}, {self._materialized_name(name)}, {default})"
            )
        writer.blank()

    def _recommend_calls(
        self, writer: Writer, foo: Procedure, alloc: Procedure, taken,
        materialized: dict[str, Argument] | None = None,
    ) -> None:
        """Call each recommend routine once per distinct set of inputs, into the sizes it produces.

        The routine is called directly (its Fortran, not a binding), so every one of its
        dummies is supplied: its inputs from the `DM_OUTPUT_FROM` plan -- the ones its own
        binding derives included -- its outputs into the local sizes, or into a local that is
        dropped where nothing here takes one (`_unused_output_locals`), its ierr as ierr.
        Each call is checked before the next is made, since each starts by setting ierr ok.
        """
        locals_ = self._unused_locals(foo, taken)
        for index, group in enumerate(self._producer_groups(foo, taken)):
            plan = group[0].roles.computed_from
            producer = plan.producer
            produced = {a.roles.computed_from.output.name.lower(): a.name for a in group}
            unused = {d.name.lower(): name for (i, d, name) in locals_ if i == index}
            supplied = {pi.name.lower(): pi for pi in plan.inputs}

            actuals = []
            for dummy in producer.arguments:
                lowered = dummy.name.lower()
                if lowered == self.conventions.error_arg.lower():
                    actuals.append((dummy.name, self.conventions.error_arg))
                elif lowered in produced:
                    actuals.append((dummy.name, produced[lowered]))
                elif lowered in unused:
                    actuals.append((dummy.name, unused[lowered]))
                elif lowered in supplied:
                    supply = supplied[lowered]
                    value = (
                        supply.argument
                        if supply.argument is not None
                        else _fortran_literal(supply.constant, dummy, self._evaluator())
                    )
                    actuals.append((dummy.name, self._actual(value, materialized)))
            # an array this call fills is allocated right above it: below the earlier calls'
            # error checks, so a failed call returns before it, and above the call, which writes
            # into it. Its extents name only what `foo` keeps and constants (`roles` row c and
            # rule P1), so they are known here.
            self._allocations(writer, [argument for argument in group if argument.is_array],
                              produced=True)
            self._call(writer, producer.name, actuals)
            if producer.has_error_argument:
                self._error_return(writer, self._producer_positions(group, alloc))

    def _producer_positions(self, group: list[Argument], alloc: Procedure) -> list[tuple[int, int]]:
        """(producer position, wrapper position) for each producer input a dummy supplied.

        A producer reports a bad input by its own position. Where that input is one of this
        wrapper's dummies, passed straight through, the position is translated to the dummy's,
        so the error names what the caller passed; everything else (a constant, a value the
        wrapper computed) names nothing the caller wrote.
        """
        plan = group[0].roles.computed_from
        producer_position = {
            dummy.name.lower(): index
            for index, dummy in enumerate(plan.producer.arguments, start=1)
        }
        own = {dummy.name.lower(): index for index, dummy in enumerate(alloc.arguments, start=1)}
        return [
            (producer_position[supply.name.lower()], own[supply.argument.lower()])
            for supply in plan.inputs
            if supply.argument is not None and supply.argument.lower() in own
        ]

    def _error_return(self, writer: Writer, positions: list[tuple[int, int]]) -> None:
        """Return on a callee's error, its position translated into this wrapper's dummy list.

        The one block written after every producer, prologue and implementation call that
        takes `ierr`. `positions` pairs a callee position with the position of the wrapper's
        own dummy that supplied it (`_call_positions`, `_producer_positions`): one it shares is
        kept, one that moved is rewritten, and everything else -- a work array or permutation
        the wrapper prepared, a recommend-sized or prologue-produced local, a mode fixed to its
        parameter, an argument another mode takes, a position outside the callee's list --
        names nothing the caller wrote and is cleared: 0, "not argument related".

        Guarded by `is_err`, so an OK `ierr` is left exactly as the callee returned it, and so
        the `select` only runs on an error. That the callee numbers its *own* dummy list is
        its author's job (codegen_guide.md 5.14); this only translates.
        """
        error = self.conventions.error_arg
        writer.line(f"if (is_err({error})) then")
        with writer.indent():
            if positions:
                writer.line(f"select case (get_err_arg_pos({error}))")
                kept = sorted(old for old, new in positions if old == new)
                if kept:
                    writer.line(f"case ({_position_ranges(kept)})")
                    with writer.indent():
                        writer.line("! the same position here")
                for old, new in sorted(positions):
                    if old == new:
                        continue
                    writer.line(f"case ({old}_int32)")
                    with writer.indent():
                        writer.line(f"call set_err_arg_pos({error}, {new}_int32)")
                writer.line("case default")
                with writer.indent():
                    writer.line(f"call clear_err_arg_pos({error})")
                writer.line("end select")
            else:
                writer.line(f"call clear_err_arg_pos({error})")
            writer.line("return")
        writer.line("end if")

    def _call_positions(self, callee_arguments, actuals, wrapper: Procedure) -> list[tuple[int, int]]:
        """(callee position, wrapper position) for each callee dummy passed one of `wrapper`'s own.

        Keyed on the actual as the call was built, not on the dummy's name: only an actual that
        IS a dummy of the wrapper being emitted maps -- for `foo`, `foo`'s own list, not its
        expert sibling's. A fixed mode is passed its parameter and a prepared work array a
        local, so neither maps. `ierr` and the prologue's `handled` never do: an error is not
        about either.
        """
        never = {self.conventions.error_arg.lower(), self.conventions.prologue_handled_arg.lower()}
        callee_position = {
            dummy.name.lower(): index for index, dummy in enumerate(callee_arguments, start=1)
        }
        own = {dummy.name.lower(): index for index, dummy in enumerate(wrapper.arguments, start=1)}
        return [
            (callee_position[dummy.lower()], own[value.lower()])
            for dummy, value in actuals
            if dummy.lower() not in never and value.lower() not in never and value.lower() in own
        ]

    def _unused_locals(self, foo: Procedure, taken) -> list[tuple[int, Argument, str]]:
        """The local each producer output the wrapper must pass and then drops goes into.

        As `(producer call, output, local)`: one local per output, named `<output>_unused`, or
        `<producer>_<output>_unused` where another producer's dropped output already has that
        name -- the two may be of different types -- and `<producer>_<output>_<n>_unused` for the
        same producer called again with other inputs. Only scalars reach here; `roles` refuses any
        other output nothing takes.
        """
        found = []
        taken_names: set[str] = set()
        for index, group in enumerate(self._producer_groups(foo, taken)):
            producer = group[0].roles.computed_from.producer
            for dummy in unused_outputs(group, self.conventions):
                name = local_name(dummy.name, UNUSED_OUTPUT)
                if name.lower() in taken_names:
                    name = local_name(f"{producer.name}_{dummy.name}", UNUSED_OUTPUT)
                count = 1
                while name.lower() in taken_names:  # the same producer, called again
                    count += 1
                    name = local_name(f"{producer.name}_{dummy.name}_{count}", UNUSED_OUTPUT)
                taken_names.add(name.lower())
                found.append((index, dummy, name))
        return found

    def _unused_output_locals(self, writer: Writer, foo: Procedure, taken) -> None:
        for _, dummy, name in self._unused_locals(foo, taken):
            writer.line(f"{dummy.type} :: {name}")

    def _unused_output_kinds(self, module: Module) -> set[str]:
        """The kinds of the locals `_unused_output_locals` declares, which `use` must import."""
        return {
            dummy.type.kind
            for _, _, group in self._producer_calls(module)
            for dummy in unused_outputs(group, self.conventions)
            if dummy.type.kind
        }

    def _actual(self, name: str, materialized) -> str:
        """The name a producer call passes: the materialised local where there is one."""
        if materialized and name in materialized:
            return self._materialized_name(name)
        return name

    def _allocations(self, writer: Writer, taken, produced: bool = False) -> None:
        """Allocate the taken-over arrays: the ones a producer fills (`produced`), or the rest.

        An array a producer fills is allocated right above the call that fills it
        (`_recommend_calls`); every other one -- a work array, a permutation -- after all the
        producer calls, since its extent may be a size one of them computes.
        """
        for argument in taken:
            if argument.is_array and is_computed(argument) is produced:
                extents = list(argument.dimension.extents)
                # a per-thread work array: as many columns as the count the implementation is
                # handed (`_impl_call`), which is one without OpenMP
                if self.conventions.is_omp_threads_arg(extents[-1]):
                    extents[-1] = self._omp_n_threads(extents[-1])
                writer.line(f"M_ALLOCATE({argument.name}({', '.join(extents)}))")

    def _omp_n_threads(self, name: str) -> str:
        """The thread count the allocating wrapper sizes and passes: `M_OMP_N_THREADS(name)`."""
        return f"{self.conventions.omp_n_threads_macro}({name})"

    def _is_omp_threads(self, argument: Argument) -> bool:
        roles = argument.roles
        if roles is not None:
            return roles.is_omp_threads
        return self.conventions.is_omp_threads_arg(argument.name)

    def _permutations(self, writer: Writer, taken, prologue) -> None:
        for argument in sorted_permutations(taken, prologue, self.conventions):
            base = argument.name[: -len(self.conventions.perm_suffix)]
            writer.line(f"call init_perm({argument.name})")
            writer.line(f"call sort_array_heapsort({base}, {argument.name})")

    # -- the implementation call --------------------------------------------------------

    def _impl_call(self, writer: Writer, foo: Procedure, module: Module,
                   wrapper: Procedure) -> None:
        """Call the implementation, translating the position it returns into `wrapper`'s list.

        `foo` carries the implementation's full picture (the expert sibling, in the allocating
        body); `wrapper` is the procedure being emitted, whose dummy list the caller numbers.
        """
        impl = self._impl_name(foo, module)
        actuals = self._impl_actuals(foo, module, wrapper)
        # A function is called for its value, which goes to the argument its result is
        # published as (`synthesize.published_result`): of the same name, on both tiers.
        implementation = self._impl(foo, module)
        result = implementation.result if implementation is not None else None
        published = result is not None and foo.argument(result.name) is not None
        self._call(writer, impl, actuals, result=result.name if published else None)
        if self._impl_declares_error(foo, module):
            self._error_return(writer, self._impl_positions(foo, module, wrapper))

    def _impl_actuals(self, foo: Procedure, module: Module,
                      wrapper: Procedure) -> list[tuple[str, str]]:
        """The implementation call's `(dummy, actual)` pairs, as `wrapper` makes the call.

        The allocating body (`wrapper` is not `foo`, its expert sibling) hands an OpenMP root
        the thread count it sized the work arrays with rather than the caller's
        (`_allocations`); that actual is no dummy of `wrapper`, so its position does not map.
        Every other tier passes the caller's on: its caller sized the arrays, or there are none.
        """
        allocating = wrapper is not foo
        arguments = self._impl_call_arguments(foo, module)
        mode_fix = self._mode_fix(foo)
        # the validating wrapper carries exactly the implementation arguments this variant
        # supplies (as dummies, or -- in the allocating body -- as the locals it prepared);
        # anything else is a mode fixed to its parameter, or an argument another mode takes and
        # this one omits (an absent optional)
        present = {a.name.lower() for a in foo.arguments}
        actuals = []
        for argument in arguments:
            lowered = argument.name.lower()
            if mode_fix is not None and lowered == mode_fix.argument.lower():
                actuals.append((argument.name, mode_fix.parameter))
            elif lowered in present:
                value = argument.name
                if allocating and self._is_omp_threads(argument):
                    value = self._omp_n_threads(argument.name)
                actuals.append((argument.name, value))
        return actuals

    def _impl_positions(self, foo: Procedure, module: Module,
                        wrapper: Procedure) -> list[tuple[int, int]]:
        """`_call_positions` for the implementation call `_impl_call` makes."""
        return self._call_positions(
            self._impl_call_arguments(foo, module),
            self._impl_actuals(foo, module, wrapper),
            wrapper,
        )

    def _call(self, writer: Writer, name: str, actuals, result: str | None = None) -> None:
        """`call name(dummy = value, ...)`, or `result = name(...)` for a function."""
        writer.line(f"call {name}(&" if result is None else f"{result} = {name}(&")
        with writer.indent():
            for index, (dummy, value) in enumerate(actuals):
                separator = "&" if index == len(actuals) - 1 else ",&"
                writer.line(f"{dummy} = {value}{separator}")
        writer.line(")")


def _position_ranges(positions: list[int]) -> str:
    """Sorted positions as `case` selectors, a contiguous run as one range: `1_int32:4_int32, 6_int32`."""
    runs: list[list[int]] = []
    for position in positions:
        if runs and position == runs[-1][-1] + 1:
            runs[-1].append(position)
        else:
            runs.append([position])
    return ", ".join(
        f"{run[0]}_int32" if len(run) == 1 else f"{run[0]}_int32:{run[-1]}_int32" for run in runs
    )


def _size_identifiers(argument: Argument) -> list[str]:
    """The identifiers an argument's declaration is sized by: its extents and its length."""
    texts = list(argument.dimension.extents)
    length = argument.type.length
    if length is not None and length.expr is not None:
        texts.append(length.expr)
    return [identifier for text in texts for identifier in _EXTENT_IDENTIFIER_RE.findall(text)]


#: The kind `tox_errors` declares its logical dummies with, which is what the emitter's own
#: `allow_nan=`/`allow_infinite=` literals must carry.
_VALIDATOR_LOGICAL_KIND = "c_bool"

#: Kind constants that live in `iso_c_binding` rather than `iso_fortran_env`. Splitting them
#: matters from the moment any argument is `logical(c_bool)`: emitting one `only:` list from
#: `iso_fortran_env` put `c_bool` in it, and gfortran rejects the module outright -- then
#: reports every dummy of that kind as having no implicit type, which reads like ten separate
#: faults rather than one bad `use`.
_C_BINDING_KINDS = frozenset({
    "c_bool", "c_char", "c_signed_char", "c_short", "c_int", "c_long", "c_long_long",
    "c_float", "c_double", "c_long_double", "c_float_complex", "c_double_complex",
    "c_size_t", "c_int8_t", "c_int16_t", "c_int32_t", "c_int64_t", "c_intptr_t",
    "c_ptr", "c_funptr",
})


def _kind_module(kind: str) -> str:
    """The intrinsic module a kind constant comes from."""
    return "iso_c_binding" if kind in _C_BINDING_KINDS else "iso_fortran_env"


def _fortran_literal(value: object, argument: Argument, evaluator=None) -> str:
    """A constant producer input, rendered back as a kinded Fortran literal.

    The kind is the producer dummy's own, so `1` handed to an `integer(int32)` input comes
    out `1_int32`, and `.false.` handed to a `logical(c_bool)` one comes out
    `.false._c_bool`. A logical was the exception until logicals stopped being default-kind:
    argument association does not convert kinds, so an unkinded literal against a
    `logical(c_bool)` dummy is a compile error rather than a silent widening.
    """
    kind = argument.type.kind
    if isinstance(value, bool):
        literal = ".true." if value else ".false."
        return f"{literal}_{kind}" if kind else literal
    if isinstance(value, int):
        return f"{value}_{kind}" if kind else str(value)
    if isinstance(value, float):
        return f"{value!r}_{kind}" if kind else repr(value)
    if isinstance(value, str):
        # a string constant, quoted as Fortran quotes it: unquoted, it was a name (rule T); and
        # padded to a fixed length, as Fortran assigns one and the bindings pass it -- a shorter
        # actual is no Fortran program (a longer one `roles` refuses)
        length = fixed_length(argument, evaluator)
        padded = value.ljust(length) if length is not None else value
        return "'" + padded.replace("'", "''") + "'"
    return str(value)
