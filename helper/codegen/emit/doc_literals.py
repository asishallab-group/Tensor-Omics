"""Rendering a Fortran doc for the language that will read it.

The author writes documentation on a Fortran impl, and the `DM_` range macros quote the
author's own Fortran back verbatim -- so a bound arrives as `0_int32` and a default as
`.false.`. Carried into a Python docstring or an R help page unchanged, those are not merely
ugly: `.false.` sits two lines under a signature that says `compute_influence=False`, and
`1_int32` is offered as the default of an argument the binding types as `str`.

Only literals are translated. Prose is the author's and is never rewritten.
"""

from __future__ import annotations

import re

#: A numeric literal carrying a Fortran kind suffix: `0_int32`, `1.0_real64`, `1.0e-9_real64`
_KINDED_NUMBER = re.compile(
    r"\b(\d+(?:\.\d*)?(?:[eEdD][-+]?\d+)?)_(?:int8|int16|int32|int64|real32|real64|real128)\b"
)

#: Ford writes inline maths as `\\( ... \\)`. Rd writes it as `\\eqn{...}` and passes the
#: LaTeX inside straight to the PDF manual, so the author's `\\frac{}{}` survives the trip --
#: whereas `\\(` is not an Rd macro at all and would fail R CMD check.
_FORD_INLINE_MATH = re.compile(r"\\\(\s*(.+?)\s*\\\)")

_LOGICALS = {
    "python": {".true.": "True", ".false.": "False"},
    "r": {".true.": "TRUE", ".false.": "FALSE"},
}


def render(text: str, language: str) -> str:
    """`text` with its Fortran literals written the way `language` writes them."""
    text = _KINDED_NUMBER.sub(r"\1", text)
    for fortran, native in _LOGICALS[language].items():
        text = re.sub(re.escape(fortran), native, text, flags=re.IGNORECASE)
    if language == "r":
        text = _FORD_INLINE_MATH.sub(r"\\eqn{\1}", text)
    return text


#: Ford's block delimiters. They mark a note or a warning in the Fortran documentation and
#: mean nothing to numpydoc or roxygen -- `@endnote` is not even a roxygen tag. The prose
#: between them is the author's and is kept; only the delimiters go.
_FORD_BLOCK_TAGS = frozenset(
    {"@note", "@endnote", "@warning", "@endwarning", "@bug", "@endbug", "@todo", "@endtodo"}
)


def is_ford_block_tag(text: str) -> bool:
    """Whether a doc line is one of Ford's block delimiters and nothing else."""
    return text.strip().lower() in _FORD_BLOCK_TAGS
