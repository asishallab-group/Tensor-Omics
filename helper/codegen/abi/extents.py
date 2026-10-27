"""The extents a binding computes itself, parsed once and written out for Python and C.

Python and the R shim allocate every output and work array before the call, and pad every
fixed-length string, so they evaluate the Fortran extent that sizes it -- `n_points - 1`,
`max(0_int32, n_timepoints - 1)`, `size(gene_to_fam)`. They used to do it with two regex
translators that passed whatever they did not recognise through unchanged (#221 B11), so an
extent whose meaning differs between the languages compiled and sized the buffer wrongly: `/`
floors in Python and truncates in Fortran (`-7/2` is -4 against -3), `mod` and `nint` round
differently, and `**` did not even compile in the shim.

So the grammar is one whose meaning is identical in Fortran, Python and C, and nothing else is
accepted (`check_binding_extents` refuses it, pointing at `DM_OUTPUT_FROM`, which is how a size
that needs Fortran to compute it is obtained):

- argument names and module constants (by the time a binding reads an extent the constant is its
  value, `CWrapper.with_constant_sizes_evaluated`)
- integer literals, a kind suffix included (`0_int32`)
- `+`, `-`, `*`, a sign at the start of an expression as Fortran places it, and parentheses
- integer `/`, truncating toward zero as Fortran does -- C's `/` on `int` already does; Python's
  `//` floors, so it is written out exactly (`_python_division`)
- `max` and `min` of two or more operands, and `size(a)` / `size(a, k)`

`parse` builds the tree; `to_python` and `to_c` write it out; `evaluate` gives a constant extent
its value; `names` says what it reads.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from typing import Callable

from ..ir.types import Intent
from .model import Conversion

# -- the tree ---------------------------------------------------------------------------------


@dataclass(frozen=True)
class Name:
    name: str


@dataclass(frozen=True)
class Integer:
    value: int


@dataclass(frozen=True)
class Negation:
    operand: object


@dataclass(frozen=True)
class Binary:
    #: one of `+ - * /`
    operator: str
    left: object
    right: object


@dataclass(frozen=True)
class Extremum:
    #: `max` or `min`
    function: str
    operands: tuple


@dataclass(frozen=True)
class Size:
    array: str
    #: the 1-based dimension, or None for the element count
    dimension: int | None


class ExtentError(ValueError):
    """An extent outside the grammar, with what in it is not allowed."""


# -- parsing ----------------------------------------------------------------------------------

_TOKEN = re.compile(
    r"\s*(?:"
    r"(?P<real>\d+\.\d*|\.\d+|\d+[eEdD])"     # a real literal, refused
    r"|(?P<integer>\d+)(?:_(?P<kind>\w+))?"
    r"|(?P<name>[A-Za-z]\w*)"
    r"|(?P<power>\*\*)"
    r"|(?P<symbol>[-+*/(),])"
    r"|(?P<other>\S)"
    r")"
)


def _tokens(text: str) -> list[tuple[str, str]]:
    tokens = []
    position = 0
    text = text.rstrip()
    while position < len(text):
        match = _TOKEN.match(text, position)
        position = match.end()
        kind = match.lastgroup if match.lastgroup != "kind" else "integer"
        if kind == "real":
            raise ExtentError(f"the real literal '{match.group('real')}'")
        if kind == "power":
            raise ExtentError("the operator '**'")
        if kind == "other":
            raise ExtentError(f"'{match.group('other')}'")
        tokens.append((kind, match.group(kind)))
    return tokens


class _Parser:
    """Fortran's level-2 expressions (F2018 10.1.2.4), the subset `parse` documents."""

    def __init__(self, text: str):
        self.text = text
        self.tokens = _tokens(text)
        self.position = 0

    def parse(self):
        if not self.tokens:
            raise ExtentError("an empty extent")
        node = self._expression()
        if self.position != len(self.tokens):
            raise ExtentError(f"'{self.tokens[self.position][1]}' after a complete expression")
        return node

    def _peek(self) -> tuple[str, str] | None:
        return self.tokens[self.position] if self.position < len(self.tokens) else None

    def _take(self, symbol: str | None = None) -> tuple[str, str]:
        token = self._peek()
        if token is None:
            raise ExtentError("an expression that ends early")
        if symbol is not None and token != ("symbol", symbol):
            raise ExtentError(f"'{token[1]}' where '{symbol}' belongs")
        self.position += 1
        return token

    def _at(self, *symbols: str) -> bool:
        token = self._peek()
        return token is not None and token[0] == "symbol" and token[1] in symbols

    def _expression(self):
        # a sign applies to the first add-operand, as in Fortran: `-a*b` is -(a*b)
        negate = False
        if self._at("+", "-"):
            negate = self._take()[1] == "-"
        node = self._term()
        if negate:
            node = Negation(node)
        while self._at("+", "-"):
            operator = self._take()[1]
            node = Binary(operator, node, self._term())
        return node

    def _term(self):
        node = self._primary()
        while self._at("*", "/"):
            operator = self._take()[1]
            node = Binary(operator, node, self._primary())
        return node

    def _primary(self):
        kind, value = self._take()
        if kind == "integer":
            return Integer(int(value))
        if kind == "symbol" and value == "(":
            node = self._expression()
            self._take(")")
            return node
        if kind == "symbol":
            # `a * -b` is an extension, not Fortran, and the sign is where it would be
            raise ExtentError(f"'{value}' where an operand belongs")
        if not self._at("("):
            return Name(value)
        function = value.lower()
        self._take("(")
        if function in ("max", "min"):
            operands = [self._expression()]
            while self._at(","):
                self._take(",")
                operands.append(self._expression())
            self._take(")")
            if len(operands) < 2:
                raise ExtentError(f"'{value}' of a single operand")
            return Extremum(function, tuple(operands))
        if function == "size":
            kind, array = self._take()
            if kind != "name" or self._at("("):
                raise ExtentError(f"'size' of something other than an argument")
            dimension = None
            if self._at(","):
                self._take(",")
                kind, literal = self._take()
                if kind != "integer" or int(literal) < 1:
                    raise ExtentError("a 'size' dimension that is not a positive literal")
                dimension = int(literal)
            self._take(")")
            return Size(array, dimension)
        raise ExtentError(f"a call to '{value}'")


@lru_cache(maxsize=None)
def parse(text: str):
    """The tree of an extent, or ExtentError naming what the grammar does not allow."""
    return _Parser(text).parse()


# -- reading the tree -------------------------------------------------------------------------


def names(node) -> set[str]:
    """Every identifier the extent reads: the names in it and the arrays it asks the size of."""
    match node:
        case Name(name):
            return {name}
        case Size(array, _):
            return {array}
        case Integer():
            return set()
        case Negation(operand):
            return names(operand)
        case Binary(_, left, right):
            return names(left) | names(right)
        case Extremum(_, operands):
            return set().union(*(names(operand) for operand in operands))
    raise TypeError(node)


def evaluate(node, values: dict[str, int] | None = None) -> int:
    """The extent's value, with Fortran's integer arithmetic; `values` gives the names."""
    values = values or {}
    match node:
        case Integer(value):
            return value
        case Name(name):
            if name not in values:
                raise ExtentError(f"'{name}', which has no value here")
            return values[name]
        case Size():
            raise ExtentError("'size', which has no value before the call")
        case Negation(operand):
            return -evaluate(operand, values)
        case Binary(operator, left, right):
            a, b = evaluate(left, values), evaluate(right, values)
            if operator == "+":
                return a + b
            if operator == "-":
                return a - b
            if operator == "*":
                return a * b
            if b == 0:
                raise ExtentError("a division by zero")
            quotient = abs(a) // abs(b)
            return quotient if (a < 0) == (b < 0) else -quotient
        case Extremum(function, operands):
            evaluated = [evaluate(operand, values) for operand in operands]
            return max(evaluated) if function == "max" else min(evaluated)
    raise TypeError(node)


# -- writing it out ---------------------------------------------------------------------------

#: How tightly each node binds, for parenthesising only where the meaning needs it
_ATOM, _PRODUCT, _SUM = 3, 2, 1


def _precedence(node) -> int:
    if isinstance(node, Binary):
        return _PRODUCT if node.operator in "*/" else _SUM
    if isinstance(node, Negation):
        return _SUM
    return _ATOM


class _Writer:
    def __init__(self, name: Callable[[str], str], size: Callable[[str, int | None], str],
                 extremum: Callable[[str, list[str]], str],
                 division: Callable[[str, str], str], division_is_atomic: bool):
        self.name = name
        self.size = size
        self.extremum = extremum
        self.division = division
        #: whether `division` writes something parenthesised already, as Python's does
        self.division_is_atomic = division_is_atomic

    def write(self, node) -> str:
        match node:
            case Integer(value):
                return str(value)
            case Name(name):
                return self.name(name)
            case Size(array, dimension):
                return self.size(array, dimension)
            case Negation(operand):
                return f"-{self._operand(operand, _PRODUCT)}"
            case Extremum(function, operands):
                return self.extremum(function, [self.write(operand) for operand in operands])
            case Binary("/", left, right):
                return self.division(self._operand(left, _ATOM), self._operand(right, _ATOM))
            case Binary(operator, left, right):
                precedence = _precedence(node)
                # left-associative: a right operand binding as tightly as the operator keeps
                # its parentheses (`a - (b - c)`, `a * (b / c)`), a left one needs none
                return (f"{self._operand(left, precedence)} {operator} "
                        f"{self._operand(right, precedence + 1)}")
        raise TypeError(node)

    def _operand(self, node, minimum: int) -> str:
        text = self.write(node)
        precedence = _precedence(node)
        if self.division_is_atomic and isinstance(node, Binary) and node.operator == "/":
            precedence = _ATOM
        if precedence < minimum:
            return f"({text})"
        return text


def _python_division(left: str, right: str) -> str:
    """Fortran's integer `/`, which truncates toward zero, in Python, whose `//` floors.

    The two agree when the operands have the same sign; otherwise the quotient of the magnitudes
    is negated. Exact for any size, which `int(a / b)` through a float is not.
    """
    return f"({left} // {right} if ({left} >= 0) == ({right} > 0) else -(-{left} // {right}))"


def to_python(text: str, name: Callable[[str], str] = lambda n: n) -> str:
    """An extent as a Python expression over the generated function's variables."""
    def size(array: str, dimension: int | None) -> str:
        # no dimension: the element count, whatever the rank; numpy's axes are 0-based
        return f"{name(array)}.size" if dimension is None else f"{name(array)}.shape[{dimension - 1}]"

    def extremum(function: str, operands: list[str]) -> str:
        return f"{function}({', '.join(operands)})"

    return _Writer(name, size, extremum, _python_division, True).write(parse(text))


def to_c(text: str, name: Callable[[str], str], size: Callable[[str, int | None], str]) -> str:
    """An extent as a C `int` expression; `size` writes `size(a[, k])` for the caller's arrays.

    C's `/` on `int` truncates toward zero, as Fortran's does. `max` and `min` become the
    two-operand helpers of `tox_marshal.h`, nested for more operands.
    """
    def extremum(function: str, operands: list[str]) -> str:
        text = operands[-1]
        for operand in reversed(operands[:-1]):
            text = f"tox_i{function}({operand}, {text})"
        return text

    def division(left: str, right: str) -> str:
        return f"{left} / {right}"

    return _Writer(name, size, extremum, division, False).write(parse(text))


# -- what the bindings compute, refused where the grammar cannot say it ----------------------


def binding_computed_extents(wrapper):
    """Every extent Python and the R shim evaluate for `wrapper`, as (argument, extent) pairs.

    What they allocate before the call -- an output, a work array -- they size by all of its
    extents, and a string by its length, input strings included, which are padded to it. An
    array whose shape travels in a shape argument is sized by that argument's product instead.
    An input array's axes are read off it, and compared: with an argument, which needs nothing
    computed, or with a constant (`fields(n, 2, m)`, `x(N_AXES)`), which is evaluated here too.
    An axis the declaration computes from arguments is not among them: C is given its extent
    and the C wrapper compares it (`CArgument.checked_against`), so it names that.
    """
    for argument in wrapper:
        if argument.is_error:
            continue
        allocated = argument.intent is Intent.OUT or argument.is_temporary
        character = argument.conversion is Conversion.CHARACTER
        if allocated and argument.is_array and argument.shape_arg is None:
            for extent in argument.dimension.extents:
                yield argument, extent
            continue
        if character and (allocated or argument.intent.is_input):
            yield argument, argument.dimension.extents[0]
        if (argument.is_array and argument.intent.is_input and argument.shape_arg is None
                and argument.conversion is not Conversion.MODE):
            axes = argument.dimension.extents[1 if character else 0:]
            for extent in axes:
                if wrapper.argument(extent) is None:
                    yield argument, extent


def check_binding_extents(binding, diagnostics) -> None:
    """Refuse every extent a binding computes that it cannot compute as Fortran does.

    Run on the binding as Python and R publish it, since an expert tier that reaches neither
    computes nothing there. A name in such an extent must be an argument the binding has in
    hand -- an input, never optional -- or a module constant: an output or an omitted optional
    has no value yet, and anything else (a module variable) none the binding can name.
    """
    for module in binding:
        for wrapper in module:
            constants = {constant.name.lower() for constant in wrapper.size_constants}
            for argument, extent in binding_computed_extents(wrapper):
                problem = _binding_extent_problem(extent, wrapper, constants)
                if problem is None:
                    continue
                diagnostics.error(
                    f"'{argument.name}' of '{wrapper.procedure.name}' is sized by '{extent}', "
                    f"which Python and R cannot compute: {problem}",
                    entity=argument.source or wrapper.procedure,
                    note=(
                        "the bindings allocate and pad by such an extent before the call, so it "
                        "may use argument names, module constants, integer literals, + - * /, "
                        "parentheses, max, min and size; obtain any other size from a "
                        "DM_OUTPUT_FROM producer"
                    ),
                )


def _binding_extent_problem(extent: str, wrapper, constants: set[str]) -> str | None:
    try:
        node = parse(extent)
    except ExtentError as error:
        return f"it uses {error}"
    for identifier in sorted(names(node)):
        if identifier.lower() in constants:
            continue
        argument = wrapper.argument(identifier)
        if argument is None:
            return f"'{identifier}' is neither an argument nor a module constant"
        if argument.optional:
            return f"it reads the optional '{identifier}', which may be absent"
        if not argument.intent.is_input:
            return f"it reads '{identifier}', which is an output"
    return None
