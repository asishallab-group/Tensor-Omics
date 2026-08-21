"""What to say about a procedure that reaches a binding language as two functions.

`foo` and `foo_expert` are not fast and slow variants, and a caller reading either docstring
has no way to tell them apart otherwise: they compute the same thing, and where the plain one
is preferable is not a matter of taste. The difference is that `foo` *derives* what
`foo_expert` takes from you -- a permutation in a particular order, a value a prologue
computes, a degenerate-input policy -- so the choice is about control, and only a caller who
wants that control should reach past `foo`.

The generator knows exactly what the difference is, so it says so rather than leaving the two
docstrings identical but for a `Generated from ...` line naming a Fortran procedure the
reader has never seen.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from ..abi.model import CBinding, CWrapper


@dataclass(frozen=True)
class TierNote:
    """The sentences one half of a pair says about the other."""

    #: the sibling as the binding language spells it
    sibling: str
    #: the arguments the expert half takes and the plain half prepares
    prepares: tuple[str, ...]
    #: how the plain half prepares them, as clauses -- and what a prologue decides
    how: tuple[str, ...]
    #: True on the expert half, which takes those things instead of preparing them
    is_expert: bool

    def lines(self, code: str) -> list[str]:
        """The note, with `code` wrapping each identifier as this language marks code.

        The parts are built with identifiers in backticks -- one marker, whatever the target
        -- and rendered here: Python keeps the backticks, R turns them into `\\code{}`, which
        is what the rest of its generated roxygen uses for text the emitter wrote. Author
        prose keeps its own backticks; this is only the emitter's own sentence.

        `how` reads as a predicate of the *allocating* half, so the same clauses serve both
        sides: this entry point does X, or the other one does and this one leaves it to you.
        """
        mark = lambda text: re.sub(r"`([^`]+)`", lambda m: code.format(m.group(1)), text)
        sibling = code.format(self.sibling)
        how = mark(_join(self.how))
        if self.is_expert:
            if not self.prepares:  # a prologue and nothing else: no argument changes hands
                return [f"The expert entry point: {sibling} {how}; this one does not."]
            return [
                f"The expert entry point: you supply {mark(_join(self.prepares))} yourself.",
                f"{sibling} {how}.",
            ]
        return [
            f"This entry point {how}.",
            f"Call {sibling} to do that yourself.",
        ]


def _join(parts) -> str:
    parts = tuple(parts)
    if len(parts) == 1:
        return parts[0]
    return ", ".join(parts[:-1]) + f" and {parts[-1]}"


def build(binding: CBinding, specs=(), conventions=None) -> dict[str, TierNote]:
    """A note per published function that has a sibling tier, keyed by its published name.

    Only pairs that *both* reach the language get one: where the expert half was dropped
    (it prepared nothing the caller could have done differently) there is one function and
    no choice to explain.
    """
    from ..config import CONVENTIONS

    conventions = conventions or CONVENTIONS
    by_impl_name = {spec.validating.name.lower(): spec for spec in specs}
    notes: dict[str, TierNote] = {}
    for module in binding:
        published = {wrapper.stripped_name: wrapper for wrapper in module}
        for name, wrapper in published.items():
            expert = f"{name}{conventions.expert_suffix}"
            if expert not in published:
                continue
            prepares, how = _prepares(wrapper, published[expert], by_impl_name, conventions)
            if not how:
                continue
            notes[name] = TierNote(expert, prepares, how, is_expert=False)
            notes[expert] = TierNote(name, prepares, how, is_expert=True)
    return notes


def _prepares(plain: CWrapper, expert: CWrapper, by_impl_name, conventions):
    """What the plain half does that the expert half leaves to the caller.

    The arguments only the expert half asks for, named -- and, for a generated pair, what is
    actually done to them, which the signature alone cannot say. A hand-written pair gets the
    argument names only: the generator does not read bodies, so it will not claim to know
    how a hand-written entry point prepares them.
    """
    theirs = {argument.name.lower() for argument in expert.arguments}
    ours = {argument.name.lower() for argument in plain.arguments}
    extra = [a.name for a in expert.arguments if a.name.lower() in theirs - ours]
    names = tuple(f"`{name}`" for name in extra)

    spec = by_impl_name.get(expert.procedure.name.lower())
    if spec is None:
        # a hand-written pair: the generator does not read bodies, so it names the arguments
        # and will not claim to know how the author's own entry point prepares them
        return (names, (f"prepares {_join(names)}",)) if names else ((), ())

    how = []
    ordered = set()
    for argument in _sorted_permutations_of(spec, conventions):
        base = argument.name[: -len(conventions.perm_suffix)]
        how.append(f"seeds `{argument.name}` and sorts it by `{base}`")
        ordered.add(argument.name.lower())
    rest = [f"`{name}`" for name in extra if name.lower() not in ordered]
    if rest:
        how.append(f"computes {_join(rest)} for you")
    if spec.prologue is not None:
        how.append(f"runs `{spec.prologue.name}` first, which may answer the call outright")
    return names, tuple(how)


def _sorted_permutations_of(spec, conventions):
    """The permutations the allocating wrapper seeds and heapsorts for this pair.

    The prologue is passed on, so one the prologue builds instead is correctly left out --
    that one is described by the prologue sentence, not as an ordering the wrapper chose.
    """
    from ..synthesize import sorted_permutations, taken_over_arguments

    taken = taken_over_arguments(spec.validating.arguments, conventions, spec.prologue)
    return sorted_permutations(taken, spec.prologue, conventions)


def _names(names) -> str:
    return _join(tuple(f"`{name}`" for name in names))
