"""Concise constructors for IR fixtures.

Tests should read like the Fortran they stand for, not like constructor calls:

    integer("n_dims", Intent.IN)
    real("vector", "(n_dims)", Intent.INOUT)
    character("name", length="*")

These exist only because the IR is constructible without Ford. That is the point of the
entity layer, and this module is the smallest demonstration of it.
"""

from pathlib import Path

from codegen.diagnostics import SourceLocation
from codegen.ir.directives import Directives
from codegen.ir.doc import Doc
from codegen.ir.entities import (
    Argument,
    Declaration,
    Meta,
    Module,
    Parameter,
    Procedure,
    Project,
)
from codegen.ir.types import (
    BaseType,
    CharacterLength,
    Dimension,
    FortranType,
    Intent,
)

C_BINDING = Meta(summary="a summary", author="AUTHOR", category="C-binding")


def _argument(type_, name, dimension, intent, doc, **kwargs):
    return Argument(
        name=name,
        type=type_,
        dimension=Dimension.parse(dimension) if isinstance(dimension, str) else dimension,
        intent=intent,
        doc=Doc.parse(doc if isinstance(doc, list) else [doc]) if doc else Doc(),
        **kwargs,
    )


def integer(name, intent=Intent.IN, dimension="", kind="int32", doc="", **kwargs):
    return _argument(FortranType(BaseType.INTEGER, kind=kind), name, dimension, intent, doc, **kwargs)


def real(name, intent=Intent.IN, dimension="", kind="real64", doc="", **kwargs):
    return _argument(FortranType(BaseType.REAL, kind=kind), name, dimension, intent, doc, **kwargs)


def complex_(name, intent=Intent.IN, dimension="", kind="real64", doc="", **kwargs):
    return _argument(FortranType(BaseType.COMPLEX, kind=kind), name, dimension, intent, doc, **kwargs)


def logical(name, intent=Intent.IN, dimension="", kind=None, doc="", **kwargs):
    return _argument(FortranType(BaseType.LOGICAL, kind=kind), name, dimension, intent, doc, **kwargs)


def character(name, intent=Intent.IN, dimension="", length="*", kind=None, doc="", **kwargs):
    type_ = FortranType(BaseType.CHARACTER, kind=kind, length=CharacterLength.parse(length))
    return _argument(type_, name, dimension, intent, doc, **kwargs)


def ierr(name="ierr"):
    return integer(name, Intent.OUT, doc="Error code")


def procedure(name, *arguments, result=None, meta=C_BINDING, doc="", directives=None, **kwargs):
    return Procedure(
        name=name,
        arguments=arguments,
        result=result,
        meta=meta,
        doc=Doc.parse(doc if isinstance(doc, list) else [doc]) if doc else Doc(),
        directives=directives or Directives(),
        **kwargs,
    )


def module(name, *procedures, parameters=(), doc="", path=None, **kwargs):
    """A module. `path` gives it a source file, as a parsed one always has.

    Only the tests that go through `generate` need it -- the generated wrappers are written
    beside the implementation they came from, so that target has to know where it was.
    """
    if path is not None:
        kwargs["location"] = SourceLocation(file=Path(path))
    return Module(
        name=name,
        procedures=procedures,
        parameters=parameters,
        doc=Doc.parse(doc if isinstance(doc, list) else [doc]) if doc else Doc(),
        **kwargs,
    )


def parameter(name, expression, kind="int32", base=BaseType.INTEGER, doc=""):
    return Parameter(
        name=name,
        type=FortranType(base, kind=kind),
        expression=expression,
        doc=Doc.parse([doc]) if doc else Doc(),
    )


def declaration(name, kind="interface", doc=""):
    """A public generic interface or derived type: a name and its documentation, nothing more."""
    return Declaration(
        name=name,
        kind=kind,
        doc=Doc.parse(doc if isinstance(doc, list) else [doc]) if doc else Doc(),
    )


def project(*modules):
    return Project(modules)
