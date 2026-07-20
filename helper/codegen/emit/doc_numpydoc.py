"""Rendering a wrapper's documentation as numpydoc.

The Fortran documentation is the source: a procedure's `summary:` becomes the opening
line, an argument's `!!` comment becomes its parameter description. Nothing is written
here that an author did not write there, which is what stops the two drifting.
"""

from __future__ import annotations

from ..abi.model import CArgument, Conversion, CWrapper
from ..ir.doc import Doc, DocTable
from ..ir.types import BaseType, Intent
from ..render import Writer


def python_type_of(argument: CArgument) -> str:
    """How to describe an argument's type to a Python caller."""
    if argument.conversion is Conversion.MODE:
        options = " | ".join(f"'{v.string}'" for v in argument.mode.values)
        return f"str, one of {options}"

    if argument.type.is_character:
        # The length is a C extent, not a numpy axis
        rank = argument.rank - 1
        return "str" if rank <= 0 else f"sequence of str, of length {argument.dimension.extents[1]}"

    scalar = {
        BaseType.INTEGER: "int",
        BaseType.REAL: "float",
        BaseType.COMPLEX: "complex",
        BaseType.LOGICAL: "bool",
    }[argument.type.base]

    if argument.is_scalar:
        return scalar

    from .python_ctypes import dtype_of

    shape = ", ".join(argument.dimension.extents)
    order = ", column-major (order='F')" if argument.rank > 1 else ""
    return f"np.ndarray[{dtype_of(argument)}] of shape ({shape},){order}"


def render_docstring(wrapper: CWrapper) -> str:
    """The whole numpydoc docstring for a generated function."""
    from .python_ctypes import PythonEmitter

    emitter = PythonEmitter()
    writer = Writer()
    # A raw docstring: the argument descriptions carry LaTeX (`\(`, `\frac`), which is an
    # invalid escape sequence in a plain string and warns (a hard error in a future Python).
    # The closing `"""` is always on its own line, so no trailing backslash can abut it.
    writer.line('r"""' + (wrapper.procedure.meta.summary or wrapper.stripped_name))

    inputs = emitter._inputs(wrapper)
    outputs = emitter._outputs(wrapper)

    if inputs:
        writer.blank()
        writer.line("Parameters")
        writer.line("----------")
        for argument in inputs:
            _parameter(writer, argument, emitter)

    writer.blank()
    writer.line("Returns")
    writer.line("-------")
    if not outputs:
        writer.line("None")
    elif len(outputs) == 1:
        _result(writer, outputs[0])
    else:
        # several outputs come back as a dict, so document it as one
        writer.line("dict")
        with writer.indent():
            writer.line("with keys:")
            writer.blank()
            for argument in outputs:
                _result(writer, argument)

    _raises(writer, wrapper)
    _notes(writer, wrapper)

    writer.line('"""')
    return writer.render()


def _parameter(writer: Writer, argument: CArgument, emitter) -> None:
    annotations = []
    if argument.intent is Intent.INOUT:
        annotations.append("modified in place")
    if argument.optional:
        annotations.append("optional")
    elif argument.has_default:
        annotations.append(f"optional, default {argument.default!r}")

    suffix = ", " + ", ".join(annotations) if annotations else ""
    writer.line(f"{argument.name} : {python_type_of(argument)}{suffix}")
    with writer.indent():
        _description(writer, argument.doc)


def _result(writer: Writer, argument: CArgument) -> None:
    writer.line(f"{argument.name} : {python_type_of(argument)}")
    with writer.indent():
        _description(writer, argument.doc)


def _description(writer: Writer, doc: Doc) -> None:
    if not doc:
        return
    for block in doc:
        if isinstance(block, DocTable):
            # A mode table is already stated as the accepted values in the type line
            continue
        writer.line(block.text)


def _raises(writer: Writer, wrapper: CWrapper) -> None:
    writer.blank()
    writer.line("Raises")
    writer.line("------")
    writer.line("ToxError")
    with writer.indent():
        writer.line("If the underlying Fortran reports an error.")


def _notes(writer: Writer, wrapper: CWrapper) -> None:
    procedure = wrapper.procedure
    link = f"{procedure.module.name}::{procedure.name}" if procedure.module else procedure.name
    writer.blank()
    writer.line("Notes")
    writer.line("-----")
    writer.line(f"Generated from the Fortran procedure `{link}`.")
