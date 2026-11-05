"""Which argument a caller actually passed, for an error that blames a derived one.

An `ierr` names a position in the Fortran procedure's dummy list, and that list contains
arguments no Python or R caller writes: an extent the binding reads off the array, a
`<x>_shape`, an `n_selected_<x>` counted from a mask. Reporting `n_elements` to someone who
called `euclidean_distance(vec1, vec2)` names nothing they can go and fix.

The generator knows where each of those came from -- that is how the binding derives them in
the first place -- so it can hand the error layer the caller's word alongside the Fortran one.
An argument with no caller-visible source (a work array, `ierr`) keeps only its Fortran name;
those are managed entirely by the generated layers and should never be blamed at all.
"""

from __future__ import annotations

from ..ir.entities import Argument, Procedure


def source_of(argument: Argument) -> str | None:
    """The caller-visible argument this one is derived from, or None if it is not derived."""
    roles = argument.roles
    if roles is None:
        return None
    if roles.extent_of:
        # an extent may size several arrays; they are checked against each other before the
        # call, so by the time Fortran can complain they agree and naming the first is true
        return roles.extent_of[0].name
    if roles.shape_of is not None:
        return roles.shape_of.name
    if roles.mask_count_of is not None:
        return roles.mask_count_of.name
    return None


def sources_of(procedure: Procedure) -> list[str | None]:
    """`source_of` for every argument, positionally, so an `arg_pos` indexes straight in."""
    return [source_of(argument) for argument in procedure.arguments]


def shapes_of(wrapper, language: str) -> list[str | None]:
    """The declared shape of each array the C wrapper compares, positionally, else None.

    The C wrapper reports a computed axis that does not match its declaration as
    `ERR_SIZE_MISMATCH` at the array's position (`CWrapper.size_checks`), and a code carries
    no text. The binding knows the declaration, so the error can say what the array must be:
    `expr` needs `sum(reps_per_tissue)` rows. Written the way `language` writes a literal.
    """
    from .doc_literals import render
    from .shape_checks import declared_shape

    checked = {check.array.name.lower(): check.array for check in wrapper.size_checks}
    return [
        render(declared_shape(checked[argument.name.lower()]), language)
        if argument.name.lower() in checked else None
        for argument in wrapper.procedure.arguments
    ]
