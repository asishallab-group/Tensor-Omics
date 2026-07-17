"""The Fortran type system, as much of it as the generator needs.

Deliberately narrow: this models what appears in the signature of a procedure marked for
export, not Fortran at large.

Two departures from the previous generator's model:

- `intent` is not part of the type. It is a property of the argument that has the type,
  and conflating them meant the result variable's intent had to be patched in by mutating
  a type after construction.
- `dimension` is not part of the type either. What the C ABI does with rank -- synthesising
  extent arguments, folding a character length in as a leading extent -- is an ABI
  decision, and keeping it out of the type stops that decision leaking in here.

Everything in this module is immutable and constructible without Ford.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum


class Intent(Enum):
    """The intent of a dummy argument."""

    IN = "in"
    OUT = "out"
    INOUT = "inout"

    @property
    def is_input(self) -> bool:
        return self in (Intent.IN, Intent.INOUT)

    @property
    def is_output(self) -> bool:
        return self in (Intent.OUT, Intent.INOUT)


class BaseType(Enum):
    """The intrinsic type of an entity.

    An enum rather than a string, so a `match` over it can be checked for completeness
    and an unsupported type fails where it is introduced rather than deep in an emitter.
    """

    INTEGER = "integer"
    REAL = "real"
    COMPLEX = "complex"
    LOGICAL = "logical"
    CHARACTER = "character"
    DERIVED = "type"

    @classmethod
    def parse(cls, text: str) -> BaseType:
        try:
            return cls(text.strip().lower())
        except ValueError:
            raise UnsupportedTypeError(f"unsupported type '{text}'") from None


class UnsupportedTypeError(Exception):
    """A type the generator has no mapping for."""


class LengthKind(Enum):
    """How a `character` entity states its length."""

    #: `character(len=*)` -- taken from the actual argument, so C must supply it
    ASSUMED = "assumed"
    #: `character(len=:)` -- implies allocatable/pointer, never interoperable
    DEFERRED = "deferred"
    #: `character(len=n)` or `character(len=123)`
    EXPLICIT = "explicit"


_INTEGER_LITERAL_RE = re.compile(r"\d+\Z")


@dataclass(frozen=True)
class CharacterLength:
    """The `len` type parameter of a `character` entity."""

    kind: LengthKind
    #: The length expression, present only for `LengthKind.EXPLICIT`
    expr: str | None = None

    def __post_init__(self):
        if (self.kind is LengthKind.EXPLICIT) != (self.expr is not None):
            raise ValueError("an explicit length needs an expression, and only it does")

    @classmethod
    def parse(cls, text: str) -> CharacterLength:
        text = text.strip()
        if text == "*":
            return cls(LengthKind.ASSUMED)
        if text == ":":
            return cls(LengthKind.DEFERRED)
        return cls(LengthKind.EXPLICIT, text)

    @property
    def is_assumed(self) -> bool:
        return self.kind is LengthKind.ASSUMED

    @property
    def is_deferred(self) -> bool:
        return self.kind is LengthKind.DEFERRED

    @property
    def is_constant(self) -> bool:
        """True for `len=123`, false for `len=n`.

        A constant length needs no extra argument from C; a variable one refers to
        another dummy argument that carries it.
        """
        return self.expr is not None and _INTEGER_LITERAL_RE.match(self.expr) is not None

    def __str__(self) -> str:
        return {
            LengthKind.ASSUMED: "*",
            LengthKind.DEFERRED: ":",
        }.get(self.kind, self.expr)


@dataclass(frozen=True)
class FortranType:
    """An intrinsic Fortran type: what it is, and its kind or length parameter."""

    base: BaseType
    kind: str | None = None
    length: CharacterLength | None = None
    #: For `BaseType.DERIVED`, the name of the type
    derived_name: str | None = None

    def __post_init__(self):
        if self.base is BaseType.CHARACTER:
            if self.length is None:
                raise ValueError("a character type needs a length")
            if self.kind is not None and self.kind.lower() not in ("c_char",):
                raise ValueError(f"unsupported character kind '{self.kind}'")
        elif self.length is not None:
            raise ValueError(f"a length is meaningless for '{self.base.value}'")

        if self.base in (BaseType.INTEGER, BaseType.REAL, BaseType.COMPLEX) and not self.kind:
            # Default kinds are processor-dependent, so an unkinded numeric argument has
            # no defensible C mapping. The codebase always kinds them.
            raise ValueError(f"a kind is required for '{self.base.value}'")

        if (self.base is BaseType.DERIVED) != (self.derived_name is not None):
            raise ValueError("a derived type needs a name, and only it does")

    @property
    def is_character(self) -> bool:
        return self.base is BaseType.CHARACTER

    @property
    def is_logical(self) -> bool:
        return self.base is BaseType.LOGICAL

    @property
    def is_c_bool(self) -> bool:
        """True when the entity is already declared with the C interoperable kind.

        Such a logical needs no conversion in the wrapper -- the one case where the
        copy can be skipped.
        """
        return self.is_logical and (self.kind or "").lower() == "c_bool"

    @property
    def needs_conversion(self) -> bool:
        """Whether a C-facing value of this type must be converted for the callee.

        Characters always: C has no Fortran string. Logicals unless already `c_bool`.
        """
        if self.is_character:
            return True
        return self.is_logical and not self.is_c_bool

    def __str__(self) -> str:
        if self.base is BaseType.DERIVED:
            return f"type({self.derived_name})"
        if self.is_character:
            kind = f", kind={self.kind}" if self.kind else ""
            return f"character(len={self.length}{kind})"
        if self.kind:
            return f"{self.base.value}({self.kind})"
        return self.base.value


@dataclass(frozen=True)
class Dimension:
    """The extents of an entity, outermost-last as written in Fortran.

    A scalar has no extents. Each extent is kept as the source expression rather than a
    parsed value, because most are names of other dummy arguments and the emitters pass
    them through verbatim.
    """

    extents: tuple[str, ...] = ()

    def __post_init__(self):
        object.__setattr__(self, "extents", tuple(e.strip() for e in self.extents))

    @classmethod
    def parse(cls, text: str) -> Dimension:
        """Parse a dimension list, with or without its enclosing parentheses.

        Splitting respects nesting, so `dimension(int(n, int32), m)` is two extents and
        not three. The previous generator split on every comma and carried a TODO about
        exactly this producing wrong extents.
        """
        text = text.strip()
        if text.startswith("(") and text.endswith(")"):
            text = text[1:-1]
        if not text.strip():
            return cls()
        return cls(tuple(split_top_level(text)))

    @property
    def rank(self) -> int:
        return len(self.extents)

    @property
    def is_scalar(self) -> bool:
        return not self.extents

    @property
    def is_assumed_shape(self) -> bool:
        """`dimension(:)` -- not interoperable, so the ABI gives it explicit extents."""
        return any(e == ":" for e in self.extents)

    @property
    def is_assumed_size(self) -> bool:
        """`dimension(*)` -- interoperable, the callee is told the size some other way."""
        return bool(self.extents) and self.extents[-1] == "*"

    @property
    def is_explicit_shape(self) -> bool:
        return bool(self.extents) and not self.is_assumed_shape and not self.is_assumed_size

    def index_of(self, extent: str) -> int | None:
        """Position of `extent`, or None. Used to map an extent back to its argument."""
        try:
            return self.extents.index(extent)
        except ValueError:
            return None

    def with_extents(self, extents) -> Dimension:
        return Dimension(tuple(extents))

    def __iter__(self):
        return iter(self.extents)

    def __len__(self) -> int:
        return len(self.extents)

    def __getitem__(self, item):
        result = self.extents[item]
        return Dimension(result) if isinstance(item, slice) else result

    def __bool__(self) -> bool:
        return bool(self.extents)

    def __str__(self) -> str:
        return f"({', '.join(self.extents)})" if self.extents else ""


def split_top_level(text: str, separator: str = ",") -> list[str]:
    """Split on `separator`, ignoring separators nested in brackets or strings.

    `int(n, int32), m` splits into two parts, not three.
    """
    parts: list[str] = []
    depth = 0
    quote: str | None = None
    current: list[str] = []

    for char in text:
        if quote is not None:
            current.append(char)
            if char == quote:
                quote = None
            continue

        if char in "\"'":
            quote = char
        elif char in "([":
            depth += 1
        elif char in ")]":
            depth -= 1
            if depth < 0:
                raise ValueError(f"unbalanced brackets in '{text}'")
        elif char == separator and depth == 0:
            parts.append("".join(current).strip())
            current = []
            continue

        current.append(char)

    if depth != 0:
        raise ValueError(f"unbalanced brackets in '{text}'")
    if quote is not None:
        raise ValueError(f"unterminated string in '{text}'")

    parts.append("".join(current).strip())
    return parts
