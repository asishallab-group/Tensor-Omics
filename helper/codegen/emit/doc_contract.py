"""The sentences the generator adds to a parameter's documentation.

Everything else a generated docstring says is the author's, rendered (`doc_literals`). These
are the generator's own, because what they describe is the binding's: a string crosses into
Fortran blank-padded, so a trailing blank is not part of the value -- Fortran semantics, like
1-based indexing, which the contract documents rather than refuses -- and a string declared
with a length cannot be longer, which is refused rather than cut short (#221). And a
`DM_EQUALS` on an extent the binding reads off an array, whose own documentation never reaches
the caller (#203). Plain text, which numpydoc and roxygen both render as it stands.
"""

from __future__ import annotations

from ..abi.model import CArgument, Conversion
from .doc_literals import render

#: Said of every string a caller passes
TRAILING_BLANKS = "Trailing blanks are not part of the value, as in Fortran: they are dropped."


def string_input_notes(argument: CArgument) -> list[str]:
    """The contract's sentences for a string the caller passes; none for anything else.

    A mode gets the blanks sentence too -- `'ward '` selects `ward` -- but no length: its values
    are listed, and one longer than them is not a mode.
    """
    if not argument.type.is_character or not argument.intent.is_input:
        return []
    notes = [TRAILING_BLANKS]
    if argument.conversion is Conversion.MODE:
        return notes
    length = argument.source.type.length if argument.source else None
    if length is not None and not length.is_assumed:
        unit = "byte" if str(length) == "1" else "bytes"
        notes.append(f"At most {length} {unit} long in UTF-8: a longer value is refused, "
                     "not cut short.")
    return notes


#: Where each language starts counting axes
_FIRST_AXIS = {"python": 0, "r": 1}


def equality_notes(argument: CArgument, language: str) -> list[str]:
    """What a `DM_EQUALS` on an extent read off this array says, where the extent is hidden.

    The binding derives such an extent from the array, so it has no parameter, and the
    "Must equal" its own documentation carries goes with it. The array it is read off is
    what the caller passes in its place, and what an error then blames (`source_of`), so
    the sentence goes there: "`n_replicates`, its axis 0, must equal `sum(reps_per_tissue)`."
    The axis is counted as `language` counts it: from 0 in Python, as numpy's `shape[0]`, and
    from 1 in R, as `dim(x)[1]`. An extent the caller passes keeps the sentence in its own
    documentation, and nothing is added.
    """
    array = argument.source
    procedure = array.parent if array is not None else None
    if procedure is None or not argument.is_array:
        return []
    notes = []
    for extent in procedure.arguments:
        equals = extent.directives.equals
        roles = extent.roles
        if equals is None or roles is None or not roles.is_derived:
            continue
        owners = roles.readable_owners
        if not owners or owners[0] is not array:
            continue
        axes = [
            str(axis)
            for axis, declared in enumerate(array.dimension.extents, start=_FIRST_AXIS[language])
            if declared.lower() == extent.name.lower()
        ]
        which = f"its axis {axes[0]}" if len(axes) == 1 else f"its axes {' and '.join(axes)}"
        notes.append(render(f"`{extent.name}`, {which}, must equal `{equals.expression}`.",
                            language))
    return notes
