"""Evaluation of Fortran constant expressions.

`DM_DEFAULT(...)` states the value an optional argument takes when omitted, and the
binding languages have to pass that value on, so the generator must know it at
generation time. The contract is that it is a *constant expression*: literals, module
parameters, and intrinsic calls over them.

This replaces the previous `eval_expr`, which rewrote the expression with a chain of
`str.replace` calls and handed the result to `eval`. That is unsound as substring
surgery, not merely unsafe:

    acos(1.0_real64)   ->  np.arcnp.cos(1.0)   ('acos'->'np.arccos', then 'cos'->'np.cos')
    .false.            ->  ''                  (eval('') raises SyntaxError)
    spike_threshold    ->  s3.14159...ke_threshold  ('pi' replaced inside a name)

Here the expression is tokenised and parsed, names resolve through an explicit symbol
table, and only whitelisted intrinsics are callable. Anything unsupported is reported
as an error naming the expression rather than raising NameError from inside `eval`.

The arithmetic is Fortran's, not Python's, because the value the bindings pass has to be
the one a Fortran caller gets: integers keep their kind (`Integer`) and overflowing it is
an error, integer division and `**` with a negative exponent truncate toward zero without
going through a float, `nint`/`anint` round half away from zero, and a negative real
raised to a real power is an error rather than a complex number, and a sign binds below `**`
(`-2**2` is -4). Real kinds are the one thing not tracked: every real is read as real64, so
an unkinded `0.95` -- a real32 0.949999988 to the compiler -- is 0.95 here.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Callable, Mapping


class ConstantError(Exception):
    """An expression that is not a constant the generator can evaluate."""


class Kind:
    """A kind name used as an argument, e.g. the `int32` in `int(x, int32)`.

    Kinds carry no value. They exist so intrinsics can accept and ignore them rather
    than the name failing to resolve as a constant.
    """

    __slots__ = ("name",)

    def __init__(self, name: str):
        self.name = name

    def __repr__(self) -> str:
        return f"Kind({self.name!r})"

    def __eq__(self, other) -> bool:
        return isinstance(other, Kind) and other.name == self.name


#: The intrinsics whose integer result has the kind their `kind` argument names; every other
#: integer-valued intrinsic returns the kind of its integer arguments (`abs`, `max`, `mod`) or
#: the default (`len`, `iachar`)
_KIND_ARGUMENT_INTRINSICS = frozenset({"int", "nint", "floor", "ceiling"})

#: Kind names that may appear as a literal suffix or an intrinsic argument, with their
#: values. The values are the ones every supported compiler (gfortran, ifx, nvfortran)
#: gives them, which is what an integer's range is read from: kind 4 is 32 bits.
KIND_VALUES: Mapping[str, int] = {
    "int8": 1, "int16": 2, "int32": 4, "int64": 8,
    "real32": 4, "real64": 8, "real128": 16,
    "c_bool": 1, "c_char": 1, "c_int": 4, "c_int8_t": 1, "c_int16_t": 2, "c_int32_t": 4,
    "c_int64_t": 8, "c_size_t": 8,
    "c_double": 8, "c_float": 4, "c_double_complex": 8, "c_float_complex": 4,
}
KIND_NAMES = frozenset(KIND_VALUES)

#: The kind of an integer literal written without one (`2`), and of what `len()` returns
DEFAULT_INTEGER_KIND = 4

#: The integer kinds there are; each holds `8 * kind` bits
_INTEGER_KINDS = frozenset({1, 2, 4, 8, 16})


def integer_kind_value(kind: str | int | None) -> int | None:
    """The value of an integer kind given by name (`int32`) or number (`4`), or None."""
    if isinstance(kind, int) and not isinstance(kind, bool):
        return kind if kind in _INTEGER_KINDS else None
    if isinstance(kind, str):
        text = kind.strip().lower()
        if text.isdigit():
            return integer_kind_value(int(text))
        value = KIND_VALUES.get(text)
        return value if value in _INTEGER_KINDS else None
    return None


def fits_integer_kind(value: int, kind: int) -> bool:
    """Whether `value` is representable in an integer of kind `kind`."""
    bits = 8 * kind
    return -(2 ** (bits - 1)) <= value < 2 ** (bits - 1)


class Integer(int):
    """An integer constant that knows its kind, as every Fortran integer does.

    Fortran evaluates `2147483647 + 1` in the kind of its operands, and overflowing it is
    an error the compiler reports; Python's unbounded `int` would carry on and hand the
    bindings a value no Fortran caller can get. So every integer the evaluator produces
    carries its kind, and constructing one out of range raises.
    """

    def __new__(cls, value, kind: int = DEFAULT_INTEGER_KIND):
        value = int(value)
        if kind not in _INTEGER_KINDS:
            raise ConstantError(f"there is no integer kind {kind}")
        if not fits_integer_kind(value, kind):
            raise ConstantError(
                f"{value} overflows an integer of kind {kind} ({8 * kind} bits)"
            )
        self = super().__new__(cls, value)
        self.kind = kind
        return self

    def __repr__(self) -> str:
        return int.__repr__(self)


def _is_integer(value) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _kind_of(value) -> int:
    return value.kind if isinstance(value, Integer) else DEFAULT_INTEGER_KIND


def _truncating_division(a: int, b: int) -> int:
    """`a / b` for integers as Fortran does it: exact, truncated toward zero.

    Never through a float: `int(a / b)` rounds `9007199254740993 / 1` to an even number.
    """
    if b == 0:
        raise ConstantError("division by zero")
    quotient = abs(a) // abs(b)
    return -quotient if (a < 0) != (b < 0) else quotient


def _round_half_away(x) -> int:
    """`nint`: to the nearest integer, a half away from zero.

    Not Python's `round`, which rounds a half to even (`nint(2.5)` would be 2), and not
    `floor(x + 0.5)`, whose addition rounds: `0.49999999999999994 + 0.5` is 1.0. The
    distance to the floor is exact, so comparing it with a half is too.
    """
    magnitude = abs(float(x))
    whole = math.floor(magnitude)
    if magnitude - whole >= 0.5:
        whole += 1
    return -whole if x < 0 else whole


def _fortran_mod(a, b):
    if b == 0:
        raise ConstantError("mod by zero")
    if _is_integer(a) and _is_integer(b):
        return a - b * _truncating_division(a, b)
    return math.fmod(a, b)


def _drop_kinds(args):
    return [a for a in args if not isinstance(a, Kind)]


def _kind_argument(args) -> int:
    """The integer kind an `int(x, kind)`-style call asks for, or the default."""
    for argument in args[1:]:
        if isinstance(argument, Kind):
            kind = integer_kind_value(argument.name)
        elif _is_integer(argument):
            kind = integer_kind_value(int(argument))
        else:
            continue
        if kind is None:
            raise ConstantError(f"'{argument}' is not an integer kind")
        return kind
    return DEFAULT_INTEGER_KIND


def _int(*args):
    return Integer(int(args[0]), _kind_argument(args))


def _nint(*args):
    return Integer(_round_half_away(args[0]), _kind_argument(args))


def _anint(x, *rest):
    return math.copysign(float(_round_half_away(x)), x)


def _floor(*args):
    return Integer(math.floor(args[0]), _kind_argument(args))


def _ceiling(*args):
    return Integer(math.ceil(args[0]), _kind_argument(args))


def _real(*args):
    return float(_drop_kinds(args)[0])


def _fortran_modulo(a, b):
    # modulo(a, b) takes the sign of the divisor (Python's % does too); mod() takes a's.
    if b == 0:
        raise ConstantError("modulo by zero")
    return a % b


def _sign(a, b):
    # sign(a, b) = |a| with the sign of b; preserves a's int/float type. b == 0 -> +|a|.
    return abs(a) if b >= 0 else -abs(a)


#: The intrinsics a default value may use. The set is generous on purpose (FES, 2026-07-29):
#: any elemental numeric/character intrinsic that has a well-defined constant value is fair
#: game, so a default like PI = 4*atan(1) evaluates rather than being rejected. Type-model
#: inquiries (huge/tiny/epsilon) are intentionally left out: their value depends on the kind
#: of the argument, and only an integer's kind is tracked (`Integer`). Every real is read
#: as real64, so `epsilon(1.0)` would resolve to real64's epsilon, not the real32 one
#: Fortran gives; `huge` of an integer could be evaluated, but no default needs it.
INTRINSICS: Mapping[str, Callable] = {
    "achar": lambda i, *rest: chr(int(i)),
    "char": lambda i, *rest: chr(int(i)),
    "iachar": lambda c, *rest: ord(c),
    "ichar": lambda c, *rest: ord(c),
    "int": _int,
    "real": _real,
    "dble": lambda x: float(x),
    "nint": _nint,
    "floor": _floor,
    "ceiling": _ceiling,
    "aint": lambda x, *rest: float(math.trunc(float(x))),
    "anint": _anint,
    "abs": abs,
    "max": max,
    "min": min,
    "mod": _fortran_mod,
    "modulo": _fortran_modulo,
    "sign": _sign,
    "len": len,
    "len_trim": lambda s: len(s.rstrip()),
    "trim": lambda s: s.rstrip(),
    "sqrt": lambda x: math.sqrt(float(x)),
    "exp": lambda x: math.exp(float(x)),
    "log": lambda x: math.log(float(x)),
    "log10": lambda x: math.log10(float(x)),
    "hypot": lambda x, y: math.hypot(float(x), float(y)),
    "sin": lambda x: math.sin(float(x)),
    "cos": lambda x: math.cos(float(x)),
    "tan": lambda x: math.tan(float(x)),
    "asin": lambda x: math.asin(float(x)),
    "acos": lambda x: math.acos(float(x)),
    "atan": lambda x: math.atan(float(x)),
    "atan2": lambda y, x: math.atan2(float(y), float(x)),
    "sinh": lambda x: math.sinh(float(x)),
    "cosh": lambda x: math.cosh(float(x)),
    "tanh": lambda x: math.tanh(float(x)),
    "asinh": lambda x: math.asinh(float(x)),
    "acosh": lambda x: math.acosh(float(x)),
    "atanh": lambda x: math.atanh(float(x)),
}


_TOKEN_RE = re.compile(
    r"""
      (?P<space>\s+)
    | (?P<logical>\.(?:true|false)\.(?:_(?P<logical_kind>[A-Za-z]\w*|\d+))?)
    | (?P<number>
          (?: \d+\.\d* (?:[edED][+-]?\d+)?
            | \.\d+    (?:[edED][+-]?\d+)?
            | \d+      (?:[edED][+-]?\d+)
            | \d+
          )
          (?:_(?P<kind_suffix>[A-Za-z]\w*|\d+))?
      )
    | (?P<string>'(?:[^']|'')*' | "(?:[^"]|"")*")
    | (?P<power>\*\*)
    | (?P<name>[A-Za-z]\w*)
    | (?P<op>[-+*/(),])
    """,
    # Fortran is case-insensitive: .TRUE., 1.0D0 and ACHAR are all legal spellings
    re.VERBOSE | re.IGNORECASE,
)


@dataclass(frozen=True)
class Token:
    kind: str
    text: str
    position: int


def tokenize(expression: str) -> list[Token]:
    tokens: list[Token] = []
    position = 0
    while position < len(expression):
        match = _TOKEN_RE.match(expression, position)
        if match is None:
            raise ConstantError(
                f"cannot parse '{expression}': unexpected '{expression[position]}' "
                f"at position {position}"
            )
        kind = match.lastgroup
        # lastgroup reports the innermost group, so recover the outer token kind
        for candidate in ("space", "logical", "number", "string", "power", "name", "op"):
            if match.group(candidate) is not None:
                kind = candidate
                break
        if kind != "space":
            tokens.append(Token(kind, match.group(), position))
        position = match.end()
    return tokens


class ConstantEvaluator:
    """Evaluates Fortran constant expressions against a table of module parameters."""

    def __init__(
        self,
        constants: Mapping[str, object] | None = None,
        intrinsics: Mapping[str, Callable] | None = None,
    ):
        # Fortran is case-insensitive, so the table is keyed in lower case
        self.constants = {
            name.lower(): value for name, value in (constants or {}).items()
        }
        self.intrinsics = dict(INTRINSICS if intrinsics is None else intrinsics)

    def evaluate(self, expression: str):
        """Evaluate `expression`, returning a Python int, float, bool or str."""
        value = self.evaluate_with_kind(expression)
        return int(value) if isinstance(value, Integer) else value

    def evaluate_with_kind(self, expression: str):
        """Evaluate `expression`, keeping an integer's kind (`Integer`).

        For a value that is evaluated again in later expressions -- a module parameter --
        so `MAX_ITER * 2` is computed in `MAX_ITER`'s kind.
        """
        if not expression or not expression.strip():
            raise ConstantError("empty expression")

        tokens = tokenize(expression)
        if not tokens:
            raise ConstantError(f"'{expression}' has no value")

        parser = _Parser(tokens, expression, self.constants, self.intrinsics)
        value = parser.parse_expression()
        parser.expect_end()
        if isinstance(value, Kind):
            raise ConstantError(f"'{expression}' is a kind, not a value")
        return value


class _Parser:
    """Recursive descent over the constant-expression grammar.

        expression  := add_operand (('+' | '-') add_operand)*
        add_operand := ('+' | '-') add_operand | term
        term        := factor (('*' | '/') signed)*
        factor      := primary ('**' signed)?        -- right associative
        signed      := ('+' | '-') signed | factor
        primary     := literal | name | name '(' args ')' | '(' expression ')'

    A sign binds as in Fortran (F2018 10.1.2.4): it applies to a whole add-operand, below
    `*` and `**`, so `-2**2` is -4 and `-a*b` is `-(a*b)`. A sign straight after `*`, `/` or
    `**` (`2*-3`, `2**-1`) is a gfortran and ifx extension, accepted the way both read it: on
    the operand that follows, `**` included, so `2*-3**2` is -18 and `2.0**-2**2` is 1/16.
    """

    def __init__(self, tokens, source, constants, intrinsics):
        self.tokens = tokens
        self.source = source
        self.constants = constants
        self.intrinsics = intrinsics
        self.index = 0

    @property
    def current(self) -> Token | None:
        return self.tokens[self.index] if self.index < len(self.tokens) else None

    def advance(self) -> Token:
        token = self.current
        self.index += 1
        return token

    def accept(self, kind: str, text: str | None = None) -> Token | None:
        token = self.current
        if token is not None and token.kind == kind and (text is None or token.text == text):
            return self.advance()
        return None

    def expect_end(self) -> None:
        if self.current is not None:
            raise ConstantError(
                f"cannot parse '{self.source}': unexpected '{self.current.text}'"
            )

    def parse_expression(self):
        value = self.parse_add_operand()
        while (token := self.accept("op", "+")) or (token := self.accept("op", "-")):
            right = self.parse_add_operand()
            value = self._arith(token.text, value, right)
        return value

    def parse_add_operand(self):
        if self.accept("op", "-"):
            return self._negate(self.parse_add_operand())
        if self.accept("op", "+"):
            return self.parse_add_operand()
        return self.parse_term()

    def parse_term(self):
        value = self.parse_factor()
        while (token := self.accept("op", "*")) or (token := self.accept("op", "/")):
            right = self.parse_signed()
            value = self._arith(token.text, value, right)
        return value

    def parse_factor(self):
        value = self.parse_primary()
        if self.accept("power"):
            # right associative: 2**3**2 is 2**(3**2)
            return self._arith("**", value, self.parse_signed())
        return value

    def parse_signed(self):
        if self.accept("op", "-"):
            return self._negate(self.parse_signed())
        if self.accept("op", "+"):
            return self.parse_signed()
        return self.parse_factor()

    def parse_primary(self):
        token = self.current
        if token is None:
            raise ConstantError(f"'{self.source}' ends unexpectedly")

        if self.accept("op", "("):
            value = self.parse_expression()
            if not self.accept("op", ")"):
                raise ConstantError(f"'{self.source}' is missing a closing bracket")
            return value

        if token.kind == "number":
            self.advance()
            return self._number(token)

        if token.kind == "logical":
            self.advance()
            # `.true._c_bool`: the kind changes the storage, not the value
            _, _, kind = token.text.partition("._")
            if kind and not (kind.isdigit() or kind.lower() in KIND_NAMES):
                raise ConstantError(f"'{token.text}' has an unknown kind '{kind}'")
            return token.text.lower().startswith(".true.")

        if token.kind == "string":
            self.advance()
            return self._string(token)

        if token.kind == "name":
            self.advance()
            if self.accept("op", "("):
                return self._call(token)
            return self._name(token)

        raise ConstantError(f"cannot parse '{self.source}': unexpected '{token.text}'")

    def _call(self, token: Token):
        name = token.text.lower()
        arguments = []
        if not self.accept("op", ")"):
            arguments.append(self.parse_expression())
            while self.accept("op", ","):
                arguments.append(self.parse_expression())
            if not self.accept("op", ")"):
                raise ConstantError(f"'{self.source}' is missing a closing bracket")

        if name not in self.intrinsics:
            raise ConstantError(
                f"'{token.text}' is not an intrinsic the generator can evaluate; "
                f"a default value must be a constant expression using one of: "
                f"{', '.join(sorted(self.intrinsics))}"
            )
        try:
            value = self.intrinsics[name](*arguments)
            if _is_integer(value) and name not in _KIND_ARGUMENT_INTRINSICS:
                kinds = [_kind_of(a) for a in arguments if _is_integer(a)]
                value = Integer(value, max(kinds, default=DEFAULT_INTEGER_KIND))
        except Exception as error:  # a ConstantError from the intrinsic included
            raise ConstantError(f"cannot evaluate '{self.source}': {error}") from None
        return self._checked(value)

    def _name(self, token: Token):
        name = token.text.lower()
        if name in self.constants:
            value = self.constants[name]
            if _is_integer(value) and not isinstance(value, Integer):
                # A table built by `Project.constant_values` carries each parameter's
                # declared kind. One built by hand does not: its value gets the default
                # kind, or int64 when it is too wide for that.
                value = Integer(value, DEFAULT_INTEGER_KIND
                                if fits_integer_kind(value, DEFAULT_INTEGER_KIND) else 8)
            return value
        if name in KIND_NAMES:
            return Kind(name)
        raise ConstantError(
            f"'{token.text}' is not a known constant; a default value must be a "
            f"literal or a module parameter"
        )

    def _number(self, token: Token):
        text = token.text
        kind = None
        if "_" in text:
            text, _, kind = text.rpartition("_")
            if not (kind.isdigit() or kind.lower() in KIND_NAMES):
                raise ConstantError(f"'{token.text}' has an unknown kind '{kind}'")

        normalised = re.sub(r"[dD]", "e", text)
        if re.search(r"[.eE]", normalised):
            return self._checked(float(normalised))
        integer_kind = DEFAULT_INTEGER_KIND if kind is None else integer_kind_value(kind)
        if integer_kind is None:
            raise ConstantError(f"'{token.text}' has a kind that is not an integer kind")
        try:
            return Integer(int(normalised), integer_kind)
        except ConstantError as error:
            raise ConstantError(f"cannot evaluate '{self.source}': {error}") from None

    @staticmethod
    def _string(token: Token):
        quote = token.text[0]
        return token.text[1:-1].replace(quote * 2, quote)

    def _arith(self, operator: str, left, right):
        if isinstance(left, Kind) or isinstance(right, Kind):
            raise ConstantError(f"cannot use a kind in arithmetic in '{self.source}'")
        integers = _is_integer(left) and _is_integer(right)
        try:
            match operator:
                case "+":
                    value = left + right
                case "-":
                    value = left - right
                case "*":
                    value = left * right
                case "/":
                    value = (_truncating_division(left, right) if integers
                             else left / right)
                case "**":
                    value = self._power(left, right, integers)
            if integers:
                # Fortran computes in the wider of the two kinds, and overflow is an error
                value = Integer(value, max(_kind_of(left), _kind_of(right)))
        except ConstantError as error:
            raise ConstantError(f"cannot evaluate '{self.source}': {error}") from None
        except ZeroDivisionError:
            raise ConstantError(f"division by zero in '{self.source}'") from None
        except (TypeError, OverflowError) as error:
            raise ConstantError(f"cannot evaluate '{self.source}': {error}") from None
        return self._checked(value)

    def _power(self, base, exponent, integers: bool):
        if integers:
            if exponent < 0:
                # `1 / base**|exponent|` in integer arithmetic: zero unless |base| is 1
                if base == 0:
                    raise ConstantError("zero raised to a negative power")
                return base ** -exponent if abs(base) == 1 else 0
            if abs(base) > 1 and exponent >= 8 * max(_INTEGER_KINDS):
                # overflows every kind; computing it first could take all the memory
                raise ConstantError(f"{base}**{exponent} overflows every integer kind")
            return base**exponent
        if isinstance(exponent, float) and base < 0:
            # Python's answer is a complex number; Fortran's is an error
            raise ConstantError("a negative real raised to a real power has no real value")
        return base**exponent

    def _checked(self, value):
        """`value`, if it is one a Fortran real or integer can hold."""
        if isinstance(value, complex):
            raise ConstantError(f"'{self.source}' has no real value")
        if isinstance(value, float) and not math.isfinite(value):
            raise ConstantError(f"'{self.source}' overflows a real")
        return value

    def _negate(self, value):
        if isinstance(value, Integer):
            try:
                return Integer(-int(value), value.kind)
            except ConstantError as error:
                raise ConstantError(f"cannot evaluate '{self.source}': {error}") from None
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return -value
        raise ConstantError(f"cannot negate a non-numeric value in '{self.source}'")
