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

#: The literal inside the sentence `DM_DEFAULT(EXPR)` expands to. Anchoring on the whole
#: sentence is what keeps this off the author's own prose: the wording is the macro's, so
#: rewriting what it quoted back is the generator editing its own output.
_DEFAULT_LITERAL = re.compile(r"(?<=The default value is )`[^`]*`")

#: How each language writes a string literal, so the default agrees with the type line
#: above it -- which already lists the accepted modes in the same style.
_STRING_QUOTE = {"python": "'", "r": '"', "fortran": "'"}


def render_mode_default(text: str, mode_string: str, language: str) -> str:
    """`text` with the `DM_DEFAULT` literal rewritten as the mode string.

    A mode is an integer in the Fortran the author wrote and a string everywhere a binding
    reads it, so the macro quotes back the one value such a caller must *not* pass:
    `mode : str, one of 'plain' | 'robust', optional, default 'robust'` used to be followed
    by ``The default value is `1`.`` Only the layers that type the mode as a string call
    this; the generated Fortran wrapper keeps the integer, because there it is one.
    """
    quote = _STRING_QUOTE[language]
    return _DEFAULT_LITERAL.sub(f"`{quote}{mode_string}{quote}`", text)


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
