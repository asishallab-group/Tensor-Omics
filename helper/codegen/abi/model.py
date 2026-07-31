"""The C wrapper model: a procedure as the C ABI presents it.

This sits between the IR and the emitters. The IR says what the Fortran is; this says what
C sees; the emitters render it. Keeping it a layer of its own is what makes the three
targets agree by construction -- the Fortran wrapper, the Python `ctypes` binding and the
R C `.Call` binding all read the *same* answer to "what does this function look like from
C", rather than each working it out again and drifting.

Everything here is derived. `c_abi.build_wrapper` is the only thing that should construct
it, and nothing downstream may modify it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from ..ir.doc import Doc
from ..ir.entities import Argument, Module, Procedure
from ..ir.roles import ModeTable
from ..ir.types import Dimension, FortranType, Intent


class Origin(Enum):
    """Where a C argument came from."""

    #: A dummy argument of the wrapped procedure, one to one
    ARGUMENT = "argument"
    #: The result of a function, which C receives as an output argument
    RESULT = "result"
    #: An extent invented because the Fortran was assumed-shape, which C cannot express
    EXTENT = "extent"
    #: The length of a `character(len=*)`, which C has no way of knowing
    STRLEN = "strlen"
    #: The error code, invented when the procedure has none of its own
    ERROR = "error"


class Conversion(Enum):
    """What must happen to a value between C and the callee."""

    NONE = "none"
    #: logical(c_bool) <-> default logical
    LOGICAL = "logical"
    #: c_char buffer <-> character(len=...)
    CHARACTER = "character"
    #: a mode string <-> the integer parameter it names
    MODE = "mode"


@dataclass(frozen=True)
class CArgument:
    """One argument of a C wrapper."""

    name: str
    #: The type as C sees it, already mapped (c_int, c_double, c_bool, c_char, ...)
    type: FortranType
    #: The extents as C sees them: a character's length is the leading one
    dimension: Dimension
    intent: Intent
    origin: Origin
    conversion: Conversion = Conversion.NONE
    optional: bool = False
    doc: Doc = field(default_factory=Doc)
    #: The IR argument this came from, absent for a wholly invented one
    source: Argument | None = None
    #: For EXTENT and STRLEN: the argument whose size this carries
    sizes: str | None = None
    #: For EXTENT: which dimension of `sizes` it is, 0-based
    axis: int | None = None
    #: For a mode argument: the values it accepts
    mode: ModeTable | None = None
    #: For an array whose shape travels separately: the name of that shape argument
    shape_arg: str | None = None
    #: Extents to multiply for the element count, empty for a scalar
    size_extents: tuple[str, ...] = ()
    #: The value an omitted optional takes, evaluated from DM_DEFAULT. Evaluated here
    #: rather than in each emitter: Python and R must agree on it, and a constant
    #: expression means the same thing in both.
    default: object | None = None
    #: Whether a default was documented, which `default is None` cannot say -- the
    #: default may legitimately be None-like
    has_default: bool = False

    @property
    def rank(self) -> int:
        return self.dimension.rank

    @property
    def is_scalar(self) -> bool:
        return self.dimension.is_scalar

    @property
    def is_array(self) -> bool:
        return not self.dimension.is_scalar

    @property
    def is_synthesised(self) -> bool:
        """Whether C sees an argument the Fortran procedure does not have."""
        return self.origin in (Origin.EXTENT, Origin.STRLEN, Origin.ERROR)

    @property
    def is_error(self) -> bool:
        return self.origin is Origin.ERROR

    @property
    def needs_conversion(self) -> bool:
        return self.conversion is not Conversion.NONE

    @property
    def is_temporary(self) -> bool:
        return bool(self.source is not None and self.source.roles and self.source.roles.is_temporary)

    def __repr__(self) -> str:
        return f"CArgument({self.name!r}, {self.origin.value})"


@dataclass(frozen=True)
class CWrapper:
    """A procedure as C sees it."""

    #: The exported symbol, e.g. `fx_cluster_c`
    name: str
    #: The symbol without the `_c` suffix, which the binding languages use
    stripped_name: str
    #: The procedure being wrapped
    procedure: Procedure
    arguments: tuple[CArgument, ...]

    @property
    def module_name(self) -> str:
        return self.procedure.module.name if self.procedure.module else ""

    @property
    def doc(self) -> Doc:
        return self.procedure.doc

    def argument(self, name: str) -> CArgument | None:
        lowered = name.lower()
        for argument in self.arguments:
            if argument.name.lower() == lowered:
                return argument
        return None

    @property
    def error_argument(self) -> CArgument | None:
        for argument in self.arguments:
            if argument.origin is Origin.ERROR or argument.name.lower() == "ierr":
                return argument
        return None

    @property
    def scalars(self) -> tuple[CArgument, ...]:
        return tuple(a for a in self.arguments if a.is_scalar)

    @property
    def arrays(self) -> tuple[CArgument, ...]:
        return tuple(a for a in self.arguments if a.is_array)

    @property
    def validation_order(self) -> tuple[CArgument, ...]:
        """The arguments to null-check, in the order they may safely be checked.

        `c_loc` may not be taken of a zero-size target, so an array cannot be checked
        before the extents that size it are known to be readable. Hence: the error code
        first, since it is the only channel for reporting anything; then every scalar,
        which includes every extent and every shape length; then the arrays.

        Within the arrays, one whose shape travels in a separate argument has to wait for
        that argument -- its element count is the product of that array's contents.

        Optionals are excluded throughout: a null pointer is how C says "not present", so
        checking one would reject the very thing it is meant to allow. This is why an
        extent or shape argument may not be optional; `ir.validate` enforces that.
        """
        error = self.error_argument
        checkable = [
            a for a in self.arguments if not a.optional and a is not error
        ]

        scalars = [a for a in checkable if a.is_scalar]
        arrays = [a for a in checkable if a.is_array]

        shape_args = {a.shape_arg for a in arrays if a.shape_arg}
        described_last = sorted(arrays, key=lambda a: a.name not in shape_args)

        ordered = []
        if error is not None:
            ordered.append(error)
        ordered.extend(scalars)
        ordered.extend(described_last)
        return tuple(ordered)

    def __iter__(self):
        return iter(self.arguments)


@dataclass(frozen=True)
class CWrapperModule:
    """The wrappers of one module, in the module that will hold them."""

    #: e.g. `fx_edges_c`
    name: str
    #: The module being wrapped
    module: Module
    wrappers: tuple[CWrapper, ...]

    @property
    def stripped_name(self) -> str:
        return self.module.name

    @property
    def doc(self) -> Doc:
        return self.module.doc

    def __iter__(self):
        return iter(self.wrappers)

    def __len__(self) -> int:
        return len(self.wrappers)


@dataclass(frozen=True)
class CBinding:
    """Every wrapper the project exports."""

    modules: tuple[CWrapperModule, ...]

    def module(self, name: str) -> CWrapperModule | None:
        for module in self.modules:
            if module.name == name or module.stripped_name == name:
                return module
        return None

    def __iter__(self):
        return iter(self.modules)

    def __len__(self) -> int:
        return len(self.modules)
