"""Rendering a wrapper's documentation as numpydoc.

The Fortran documentation is the source: a procedure's `summary:` becomes the opening
line, an argument's `!!` comment becomes its parameter description. Nothing is written
here that an author did not write there, which is what stops the two drifting.
"""

from __future__ import annotations

from ..abi.model import CArgument, Conversion, CWrapper
from ..ir.doc import Doc, DocTable
from .doc_links import render_spans as _spans
from .doc_literals import is_ford_block_tag, render as _render
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
    return _render(f"np.ndarray[{dtype_of(argument)}] of shape ({shape},){order}", "python")


def render_docstring(wrapper: CWrapper, emitter=None) -> str:
    """The whole numpydoc docstring for a generated function.

    `emitter` carries the link resolver; without one (a unit test) a Ford link renders as
    plain code rather than a cross-reference, which is the honest fallback anyway.
    """
    from .python_ctypes import PythonEmitter

    emitter = emitter if emitter is not None else PythonEmitter()
    writer = Writer()
    # A raw docstring: the argument descriptions carry LaTeX (`\(`, `\frac`), which is an
    # invalid escape sequence in a plain string and warns (a hard error in a future Python).
    # The closing `"""` is always on its own line, so no trailing backslash can abut it.
    writer.line('r"""' + _render(wrapper.procedure.meta.summary or wrapper.stripped_name, "python"))

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
        _result(writer, outputs[0], emitter)
    else:
        # several outputs come back as a dict, so document it as one
        writer.line("dict")
        with writer.indent():
            writer.line("with keys:")
            writer.blank()
            for argument in outputs:
                _result(writer, argument, emitter)

    _raises(writer, wrapper)
    _notes(writer, wrapper, emitter)

    writer.line('"""')
    return writer.render()


def _parameter(writer: Writer, argument: CArgument, emitter) -> None:
    annotations = []
    # an array really is mutated through the buffer the caller passed; a scalar cannot be --
    # Python has no mutable int, so the binding passes by value and returns the new one
    if argument.intent is Intent.INOUT and not argument.is_scalar:
        annotations.append("modified in place")
    if argument.optional:
        annotations.append("optional")
    elif argument.has_default:
        annotations.append(f"optional, default {argument.default!r}")

    suffix = ", " + ", ".join(annotations) if annotations else ""
    writer.line(f"{argument.name} : {python_type_of(argument)}{suffix}")
    with writer.indent():
        _description(writer, argument.doc, emitter)


def _result(writer: Writer, argument: CArgument, emitter=None) -> None:
    annotation = ""
    # said here rather than left for the caller to discover from a ValueError. Only the
    # arrays: a scalar and a decoded string are immutable in Python anyway
    if argument.is_array and not argument.type.is_character:
        annotation = ", read-only"
    writer.line(f"{argument.name} : {python_type_of(argument)}{annotation}")
    with writer.indent():
        _description(writer, argument.doc, emitter)
        if annotation:
            writer.line("A result is a value; call `.copy()` to obtain a modifiable array.")


def _resolver(emitter):
    """The link resolver an emitter carries, absent in a unit test that builds one by hand."""
    return getattr(emitter, "links", None)


def _description(writer: Writer, doc: Doc, emitter=None) -> None:
    if not doc:
        return
    for block in doc:
        if isinstance(block, DocTable):
            # A mode table is already stated as the accepted values in the type line
            continue
        if is_ford_block_tag(block.text):
            continue
        writer.line(_render(_spans(block, _resolver(emitter), "python"), "python"))


def _raises(writer: Writer, wrapper: CWrapper) -> None:
    writer.blank()
    writer.line("Raises")
    writer.line("------")
    writer.line("ToxError")
    with writer.indent():
        writer.line("If the underlying Fortran reports an error.")


def _notes(writer: Writer, wrapper: CWrapper, emitter=None) -> None:
    procedure = wrapper.procedure
    # the procedure, not just the module: an error message names an argument from
    # *its* dummy list -- including the extents, work arrays and ierr that no caller of this
    # binding passes, and, on the expert tier, several this one does not either. Naming it is
    # what lets a reader take "(argument 'n_dscale_elements')" back to a signature that has one.
    origin = (
        f"{procedure.module.name}::{procedure.name}" if procedure.module else procedure.name
    )
    writer.blank()
    writer.line("Notes")
    writer.line("-----")
    writer.line(f"Generated from the Fortran procedure `{origin}`, whose argument names are")
    writer.line("the ones an error message reports.")
    note = getattr(emitter, "tiers", {}).get(wrapper.stripped_name)
    if note is not None:
        writer.blank()
        for line in note.lines("`{}`"):
            writer.line(line)
