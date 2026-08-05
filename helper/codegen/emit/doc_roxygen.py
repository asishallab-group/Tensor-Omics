"""Rendering a wrapper's documentation as roxygen2.

The R analogue of `doc_numpydoc`: the Fortran documentation is the source, rendered as the
`#'` comment block roxygen2 turns into an `.Rd` help page. Nothing is written here that an
author did not write in the Fortran.
"""

from __future__ import annotations

from ..abi.model import CArgument, Conversion, CWrapper
from ..ir.doc import Doc, DocTable
from .doc_literals import render as _render
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

    notes = _notes(wrapper.doc)
    if notes:
        for line in notes:
            writer.line(f"#' {line}" if line else "#'")
        writer.line("#'")

    for argument in emitter._inputs(wrapper):
        writer.line(f"#' @param {argument.name} {_param_line(argument)}")

    writer.line(f"#' @return {_return_line(wrapper, emitter)}")
    link = f"{procedure.module.name}::{procedure.name}" if procedure.module else procedure.name
    writer.line("#'")
    writer.line(f"#' Generated from the Fortran procedure \\code{{{link}}}.")
    writer.line("#' @export")
    return writer.render()


def _param_line(argument: CArgument) -> str:
    description = _first_line(argument.doc)
    type_ = r_type_of(argument)
    if description:
        return f"{type_}. {description}"
    return type_


def _return_line(wrapper: CWrapper, emitter) -> str:
    outputs = emitter._outputs(wrapper)
    if not outputs:
        return "invisibly `NULL`; called for its effect."
    if len(outputs) == 1:
        description = _first_line(outputs[0].doc)
        return description or r_type_of(outputs[0])
    names = ", ".join(f"`{a.name}`" for a in outputs)
    return f"a named list with elements {names}."


def _notes(doc: Doc) -> list[str]:
    """Prose from the procedure doc, tables dropped (values are stated on the @param)."""
    lines = []
    for block in doc:
        if isinstance(block, DocTable):
            continue
        lines.append(_render(block.text, "r"))
    while lines and not lines[-1]:
        lines.pop()
    return lines


def _first_line(doc: Doc) -> str:
    for block in doc:
        if not isinstance(block, DocTable) and block.text.strip():
            return _render(block.text, "r")
    return ""
