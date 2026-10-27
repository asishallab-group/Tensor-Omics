"""The sentences the input contract adds to a string argument's documentation (#221).

Everything else a generated docstring says is the author's, rendered (`doc_literals`). These
are the generator's own, because what they describe is the binding's: a string crosses into
Fortran blank-padded, so a trailing blank is not part of the value -- Fortran semantics, like
1-based indexing, which the contract documents rather than refuses -- and a string declared
with a length cannot be longer, which is refused rather than cut short. Plain text, which
numpydoc and roxygen both render as it stands.
"""

from __future__ import annotations

from ..abi.model import CArgument, Conversion

#: Said of every string a caller passes
TRAILING_BLANKS = "Trailing blanks are not part of the value, as in Fortran: they are dropped."


def string_input_notes(argument: CArgument) -> list[str]:
    """The contract's sentences for a string the caller passes; none for anything else.

    A mode is left out: its values are listed, and a blank is in none of them.
    """
    if (not argument.type.is_character or argument.conversion is Conversion.MODE
            or not argument.intent.is_input):
        return []
    notes = [TRAILING_BLANKS]
    length = argument.source.type.length if argument.source else None
    if length is not None and not length.is_assumed:
        unit = "byte" if str(length) == "1" else "bytes"
        notes.append(f"At most {length} {unit} long in UTF-8: a longer value is refused, "
                     "not cut short.")
    return notes
