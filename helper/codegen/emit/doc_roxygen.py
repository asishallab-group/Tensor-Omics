"""Rendering a wrapper's documentation as roxygen2.

The R analogue of `doc_numpydoc`: the Fortran documentation is the source, rendered as the
`#'` comment block roxygen2 turns into an `.Rd` help page. Nothing is written here that an
author did not write in the Fortran.
"""

from __future__ import annotations

from ..abi.model import CArgument, Conversion, CWrapper
from ..ir.doc import Doc, DocTable
from .doc_links import render_spans as _spans
from .doc_literals import is_ford_block_tag, render as _render, render_mode_default
from ..ir.types import BaseType
from ..render import Writer


def r_type_of(argument: CArgument) -> str:
    """How to describe an argument's type to an R caller."""
    if argument.conversion is Conversion.MODE:
        options = ", ".join(f'"{v.string}"' for v in argument.mode.values)
        return f"a string, one of {options}"

    if argument.type.is_character:
        return "a character vector" if argument.rank > 1 else "a string"

    noun = {
        BaseType.INTEGER: "integer",
        BaseType.REAL: "numeric",
        BaseType.COMPLEX: "complex",
        BaseType.LOGICAL: "logical",
    }[argument.type.base]

    if argument.is_scalar:
        return f"a {noun} scalar"
    if argument.rank == 1:
        return f"a {noun} vector"
    if argument.rank == 2:
        return f"a {noun} matrix"
    return f"a {noun} array of rank {argument.rank}"


def render_roxygen(wrapper: CWrapper, emitter) -> str:
    """The roxygen block for a generated R function."""
    writer = Writer()
    procedure = wrapper.procedure

    title = _render(procedure.meta.summary or wrapper.stripped_name, "r")
    writer.line(f"#' {title}")
    writer.line("#'")

    resolver = getattr(emitter, "links", None)
    notes = _notes(wrapper.doc, resolver)
    if notes:
        for line in notes:
            writer.line(f"#' {line}" if line else "#'")
        writer.line("#'")

    # the procedure, not just the module: an error message names an argument from
    # *its* dummy list -- including the extents, work arrays and ierr that no caller of this
    # binding passes, and, on the expert tier, several this one does not either. Naming it is
    # what lets a reader take "(argument 'n_dscale_elements')" back to a signature that has one.
    #
    # And it goes here, in the description, rather than after @return: roxygen gives untagged
    # text to whichever tag precedes it, so trailing it after @return filed the provenance
    # under Value, where it described the return value of nothing.
    origin = (
        f"{procedure.module.name}::{procedure.name}" if procedure.module else procedure.name
    )
    writer.line(
        f"#' Generated from the Fortran procedure \\code{{{origin}}}, whose argument names"
    )
    writer.line("#' are the ones an error message reports.")
    writer.line("#'")

    note = getattr(emitter, "tiers", {}).get(wrapper.stripped_name)
    if note is not None:
        for line in note.lines("\\code{{{}}}"):
            writer.line(f"#' {line}")
        writer.line("#'")

    for argument in emitter._inputs(wrapper):
        first, *rest = _param_lines(argument, resolver)
        writer.line(f"#' @param {argument.name} {first}")
        # roxygen continues a tag until the next one, so the rest of the description is
        # indented under it rather than dropped
        for line in rest:
            writer.line(f"#'   {line}")

    for line in _return_lines(wrapper, emitter):
        writer.line(line)
    writer.line("#' @export")
    return writer.render()


def _param_lines(argument: CArgument, resolver=None) -> list[str]:
    """The `@param` text, as the line the tag opens with plus any continuation lines.

    An argument's description is often several source lines; taking only the first silently
    dropped the rest, which is how "use `mask_chunk_count` for calculation" stopped reaching
    an R reader.
    """
    type_ = r_type_of(argument)
    lines = _prose_lines(argument.doc, resolver, argument.mode_default)
    if not lines:
        return [type_]
    return [f"{type_}. {lines[0]}", *lines[1:]]


def _prose_lines(doc: Doc, resolver=None, mode_default: str | None = None) -> list[str]:
    """Every non-table line of a doc, blanks trimmed from both ends."""
    lines = [
        _render(_spans(block, resolver, "r"), "r")
        for block in doc
        if not isinstance(block, DocTable) and not is_ford_block_tag(block.text)
    ]
    if mode_default is not None:
        lines = [render_mode_default(line, mode_default, "r") for line in lines]
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    return lines


def _return_lines(wrapper: CWrapper, emitter) -> list[str]:
    """The `@return` tag, as the lines it spans.

    Says what R hands back, in R's terms: the type first, as `@param` does, then the
    argument's own description. A list of several outputs documents each element with
    `\\item`, which is how `\\value` describes a list's components -- naming them and
    nothing else threw away every word the author wrote about them.
    """
    outputs = emitter._outputs(wrapper)
    if not outputs:
        return ["#' @return invisibly `NULL`; called for its effect."]

    if len(outputs) == 1:
        first, *rest = _param_lines(outputs[0], getattr(emitter, "links", None))
        return [f"#' @return {first}", *(f"#'   {line}" for line in rest)]

    lines = ["#' @return a named list with elements:"]
    for argument in outputs:
        first, *rest = _param_lines(argument, getattr(emitter, "links", None))
        lines.append(f"#'   \\item{{{argument.name}}}{{{first}")
        lines += [f"#'     {line}" for line in rest]
        lines[-1] += "}"
    return lines


def _notes(doc: Doc, resolver=None) -> list[str]:
    """Prose from the procedure doc, tables dropped (values are stated on the @param)."""
    lines = []
    for block in doc:
        if isinstance(block, DocTable) or is_ford_block_tag(block.text):
            continue
        lines.append(_render(_spans(block, resolver, "r"), "r"))
    while lines and not lines[-1]:
        lines.pop()
    return lines


def _first_line(doc: Doc, resolver=None) -> str:
    for block in doc:
        if not isinstance(block, DocTable) and block.text.strip():
            return _render(_spans(block, resolver, "r"), "r")
    return ""
