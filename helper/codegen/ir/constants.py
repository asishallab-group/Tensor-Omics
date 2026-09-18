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


#: Kind names that may appear as a literal suffix or an intrinsic argument.
KIND_NAMES = frozenset(
    {
        "int8", "int16", "int32", "int64",
        "real32", "real64", "real128",
        "c_bool", "c_char", "c_int", "c_int32_t", "c_int64_t",
        "c_double", "c_float", "c_double_complex", "c_size_t",
    }
)


def _fortran_mod(a, b):
    if b == 0:
        raise ConstantError("mod by zero")
    return a - b * int(a / b)


def _drop_kinds(args):
    return [a for a in args if not isinstance(a, Kind)]


def _int(*args):
    return int(_drop_kinds(args)[0])


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
#: of the argument, which this evaluator does not track, so they would resolve unsoundly.
INTRINSICS: Mapping[str, Callable] = {
    "achar": lambda i, *rest: chr(int(i)),
    "char": lambda i, *rest: chr(int(i)),
    "iachar": lambda c, *rest: ord(c),
    "ichar": lambda c, *rest: ord(c),
    "int": _int,
    "real": _real,
    "dble": lambda x: float(x),
    "nint": lambda x, *rest: int(round(float(x))),
    "floor": lambda x, *rest: math.floor(float(x)),
    "ceiling": lambda x, *rest: math.ceil(float(x)),
    "aint": lambda x, *rest: float(math.trunc(float(x))),
    "anint": lambda x, *rest: float(round(float(x))),
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
    | (?P<logical>\.(?:true|false)\.)
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

        expression := term (('+' | '-') term)*
        term       := factor (('*' | '/') factor)*
        factor     := unary ('**' factor)?          -- right associative
        unary      := ('+' | '-') unary | primary
        primary    := literal | name | name '(' args ')' | '(' expression ')'
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
        value = self.parse_term()
        while (token := self.accept("op", "+")) or (token := self.accept("op", "-")):
            right = self.parse_term()
            value = self._arith(token.text, value, right)
        return value

    def parse_term(self):
        value = self.parse_factor()
        while (token := self.accept("op", "*")) or (token := self.accept("op", "/")):
            right = self.parse_factor()
            value = self._arith(token.text, value, right)
        return value

    def parse_factor(self):
        value = self.parse_unary()
        if self.accept("power"):
            # right associative: 2**3**2 is 2**(3**2)
            return self._arith("**", value, self.parse_factor())
        return value

    def parse_unary(self):
        if self.accept("op", "-"):
            return self._negate(self.parse_unary())
        if self.accept("op", "+"):
            return self.parse_unary()
        return self.parse_primary()

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
            return token.text.lower() == ".true."

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
            return self.intrinsics[name](*arguments)
        except ConstantError:
            raise
        except Exception as error:
            raise ConstantError(f"cannot evaluate '{self.source}': {error}") from None

    def _name(self, token: Token):
        name = token.text.lower()
        if name in self.constants:
            return self.constants[name]
        if name in KIND_NAMES:
            return Kind(name)
        raise ConstantError(
            f"'{token.text}' is not a known constant; a default value must be a "
            f"literal or a module parameter"
        )

    def _number(self, token: Token):
        text = token.text
        if "_" in text:
            text, _, kind = text.rpartition("_")
            if not (kind.isdigit() or kind.lower() in KIND_NAMES):
                raise ConstantError(f"'{token.text}' has an unknown kind '{kind}'")

        normalised = re.sub(r"[dD]", "e", text)
        if re.search(r"[.eE]", normalised):
            return float(normalised)
        return int(normalised)

    @staticmethod
    def _string(token: Token):
        quote = token.text[0]
        return token.text[1:-1].replace(quote * 2, quote)

    def _arith(self, operator: str, left, right):
        if isinstance(left, Kind) or isinstance(right, Kind):
            raise ConstantError(f"cannot use a kind in arithmetic in '{self.source}'")
        try:
            match operator:
                case "+":
                    return left + right
                case "-":
                    return left - right
                case "*":
                    return left * right
                case "/":
                    # Fortran integer division truncates
                    if isinstance(left, int) and isinstance(right, int):
                        if right == 0:
                            raise ConstantError(f"division by zero in '{self.source}'")
                        return int(left / right)
                    return left / right
                case "**":
                    return left**right
        except ConstantError:
            raise
        except ZeroDivisionError:
            raise ConstantError(f"division by zero in '{self.source}'") from None
        except TypeError as error:
            raise ConstantError(f"cannot evaluate '{self.source}': {error}") from None

    def _negate(self, value):
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return -value
        raise ConstantError(f"cannot negate a non-numeric value in '{self.source}'")
