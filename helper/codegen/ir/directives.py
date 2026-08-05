"""Documentation directives: the `DM_*` macros, recognised and typed.

A `DM_*` macro has already expanded to prose by the time the generator sees it, so
recognising one means matching that prose. The patterns come from expanding the macro
definitions themselves (see `frontend.macros`) and are injected here, which keeps `ir`
free of the preprocessor while the meaning of a directive stays in the IR.

The previous generator kept raw `re.Match` objects in a dict keyed by string, so every
use site read `meta["required_if_mode"].group("mode_var_name")` -- untyped, and silently
None when the key was absent. Here each directive is a value with named fields, and an
argument exposes them through `Directives`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum

from .doc import Doc, DocLine


class DirectiveError(Exception):
    """A directive that cannot be understood or contradicts another."""

    def __init__(self, message: str, line_number: int | None = None, note: str | None = None):
        super().__init__(message)
        self.line_number = line_number
        self.note = note


@dataclass(frozen=True)
class Default:
    """`DM_DEFAULT(VALUE)` -- the value an omitted optional argument takes.

    `expression` is the Fortran source text as written. It is evaluated later, where the
    module constants it may refer to are known.
    """

    expression: str
    line_number: int | None = None


@dataclass(frozen=True)
class RequiredIfMode:
    """`DM_REQUIRED_IF_MODE(MODE_ARG, MODULE, MODE_PARAM)`.

    An optional argument with no default that is required in exactly one mode.
    """

    mode_arg: str
    module: str
    mode_param: str
    line_number: int | None = None


@dataclass(frozen=True)
class OptionalOutput:
    """`DM_OPTIONAL_OUTPUT` -- an output the caller may decline to receive."""

    line_number: int | None = None


@dataclass(frozen=True)
class ResultSizeIs:
    """`DM_RESULT_SIZE_IS(ARGUMENT)` -- how many leading elements carry results."""

    argument: str
    line_number: int | None = None


class OutputFromMode(Enum):
    #: The binding languages call the procedure themselves
    AUTO = "auto"
    #: The caller has to call it; the generator only documents the fact
    JUST_INFO = "just_info"


@dataclass(frozen=True)
class OutputFrom:
    """`DM_OUTPUT_FROM(ARGUMENT, PROCEDURE, MODULE, MODE)`.

    This argument's value is produced by another procedure, typically a work array size
    that cannot be foreseen.
    """

    argument: str
    procedure: str
    module: str
    mode: OutputFromMode
    line_number: int | None = None

    @property
    def is_automatic(self) -> bool:
        return self.mode is OutputFromMode.AUTO


@dataclass(frozen=True)
class Minimum:
    """`DM_MIN(EXPR)` -- the inclusive lower bound of a numeric argument.

    `expression` is Fortran source, evaluated where the constants it refers to are known.
    Wrap it in `above(...)` for an exclusive bound.
    """

    expression: str
    line_number: int | None = None


@dataclass(frozen=True)
class Maximum:
    """`DM_MAX(EXPR)` -- the inclusive upper bound of a numeric argument.

    Wrap `expression` in `below(...)` for an exclusive bound.
    """

    expression: str
    line_number: int | None = None


@dataclass(frozen=True)
class Sentinel:
    """`DM_SENTINEL(EXPR)` -- a value accepted regardless of the range (e.g. a marker)."""

    expression: str
    line_number: int | None = None


@dataclass(frozen=True)
class AllowNan:
    """`DM_ALLOW_NAN` -- opt a real argument out of the default NaN rejection."""

    line_number: int | None = None


@dataclass(frozen=True)
class AllowInfinite:
    """`DM_ALLOW_INFINITE` -- opt a real argument out of the default infinity rejection."""

    line_number: int | None = None


Directive = (
    Default
    | RequiredIfMode
    | OptionalOutput
    | ResultSizeIs
    | OutputFrom
    | Minimum
    | Maximum
    | Sentinel
    | AllowNan
    | AllowInfinite
)


@dataclass(frozen=True)
class Directives:
    """The directives found on one documented entity."""

    default: Default | None = None
    required_if_mode: RequiredIfMode | None = None
    optional_output: OptionalOutput | None = None
    result_size_is: ResultSizeIs | None = None
    output_from: OutputFrom | None = None
    minimum: Minimum | None = None
    maximum: Maximum | None = None
    sentinel: Sentinel | None = None
    allow_nan: AllowNan | None = None
    allow_infinite: AllowInfinite | None = None

    @property
    def all(self) -> tuple[Directive, ...]:
        found = (
            self.default,
            self.required_if_mode,
            self.optional_output,
            self.result_size_is,
            self.output_from,
            self.minimum,
            self.maximum,
            self.sentinel,
            self.allow_nan,
            self.allow_infinite,
        )
        return tuple(directive for directive in found if directive is not None)

    @property
    def has_default(self) -> bool:
        return self.default is not None

    @property
    def is_optional_output(self) -> bool:
        return self.optional_output is not None

    @property
    def has_range(self) -> bool:
        """A `DM_MIN`/`DM_MAX`/`DM_SENTINEL` bound was documented."""
        return (
            self.minimum is not None
            or self.maximum is not None
            or self.sentinel is not None
        )

    @property
    def allows_nan(self) -> bool:
        return self.allow_nan is not None

    @property
    def allows_infinite(self) -> bool:
        return self.allow_infinite is not None

    def __bool__(self) -> bool:
        return bool(self.all)


#: Argument fragments used to build the recognition patterns. Lazy, so a line carrying a
#: directive plus further prose does not let a group swallow the rest of the line.
def _group(name: str) -> str:
    return f"(?P<{name}>.+?)"


@dataclass(frozen=True)
class DirectivePatterns:
    """Compiled patterns for each directive, built from the macro definitions."""

    default: re.Pattern
    required_if_mode: re.Pattern
    optional_output: re.Pattern
    result_size_is: re.Pattern
    output_from_auto: re.Pattern
    output_from_just_info: re.Pattern
    minimum: re.Pattern
    maximum: re.Pattern
    sentinel: re.Pattern
    allow_nan: re.Pattern
    allow_infinite: re.Pattern

    @staticmethod
    def templates() -> dict[str, str]:
        """Macro invocations whose arguments are regex groups, keyed by field name."""
        return {
            "default": f"DM_DEFAULT({_group('expression')})",
            "required_if_mode": (
                f"DM_REQUIRED_IF_MODE({_group('mode_arg')}, {_group('module')}, "
                f"{_group('mode_param')})"
            ),
            "optional_output": "DM_OPTIONAL_OUTPUT",
            "result_size_is": f"DM_RESULT_SIZE_IS({_group('argument')})",
            "output_from_auto": (
                f"DM_OUTPUT_FROM({_group('argument')}, {_group('procedure')}, "
                f"{_group('module')}, AUTO)"
            ),
            "output_from_just_info": (
                f"DM_OUTPUT_FROM({_group('argument')}, {_group('procedure')}, "
                f"{_group('module')}, JUST_INFO)"
            ),
            "minimum": f"DM_MIN({_group('expression')})",
            "maximum": f"DM_MAX({_group('expression')})",
            "sentinel": f"DM_SENTINEL({_group('expression')})",
            "allow_nan": "DM_ALLOW_NAN",
            "allow_infinite": "DM_ALLOW_INFINITE",
        }


class DirectiveParser:
    """Finds directives in a `Doc`."""

    def __init__(self, patterns: DirectivePatterns):
        self.patterns = patterns

    def parse(self, doc: Doc) -> Directives:
        found: dict[str, Directive] = {}

        for line in doc.all_lines:
            for field, directive in self._parse_line(line):
                if field in found:
                    raise DirectiveError(
                        f"duplicate {_macro_name(field)} directive",
                        line.line_number,
                        note="an entity can carry each directive only once",
                    )
                found[field] = directive

        self._reject_contradictions(found)
        return Directives(**found)

    def _parse_line(self, line: DocLine):
        text = line.text
        number = line.line_number

        if (match := self.patterns.default.search(text)) is not None:
            yield "default", Default(match.group("expression").strip(), number)

        if (match := self.patterns.required_if_mode.search(text)) is not None:
            yield "required_if_mode", RequiredIfMode(
                mode_arg=match.group("mode_arg").strip(),
                module=match.group("module").strip(),
                mode_param=match.group("mode_param").strip(),
                line_number=number,
            )

        if self.patterns.optional_output.search(text) is not None:
            yield "optional_output", OptionalOutput(number)

        if (match := self.patterns.result_size_is.search(text)) is not None:
            yield "result_size_is", ResultSizeIs(match.group("argument").strip(), number)

        for pattern, mode in (
            (self.patterns.output_from_auto, OutputFromMode.AUTO),
            (self.patterns.output_from_just_info, OutputFromMode.JUST_INFO),
        ):
            if (match := pattern.search(text)) is not None:
                yield "output_from", OutputFrom(
                    argument=match.group("argument").strip(),
                    procedure=match.group("procedure").strip(),
                    module=match.group("module").strip(),
                    mode=mode,
                    line_number=number,
                )
                break

        if (match := self.patterns.minimum.search(text)) is not None:
            yield "minimum", Minimum(match.group("expression").strip(), number)

        if (match := self.patterns.maximum.search(text)) is not None:
            yield "maximum", Maximum(match.group("expression").strip(), number)

        if (match := self.patterns.sentinel.search(text)) is not None:
            yield "sentinel", Sentinel(match.group("expression").strip(), number)

        if self.patterns.allow_nan.search(text) is not None:
            yield "allow_nan", AllowNan(number)

        if self.patterns.allow_infinite.search(text) is not None:
            yield "allow_infinite", AllowInfinite(number)

    @staticmethod
    def _reject_contradictions(found: dict[str, Directive]) -> None:
        # `DM_DEFAULT` together with `DM_REQUIRED_IF_MODE` is contradictory only where the
        # mode is resolved at runtime: there the argument is always passed on, so "required
        # in that mode" says nothing. In a mode-*split* kernel the same pairing is meaningful
        # -- the directive scopes the argument to its mode, and the default applies within it
        # -- and whether the kernel splits is not knowable from the documentation alone. The
        # rule therefore lives in the semantic pass, which has the mode table.
        return


_MACRO_NAMES = {
    "default": "DM_DEFAULT",
    "required_if_mode": "DM_REQUIRED_IF_MODE",
    "optional_output": "DM_OPTIONAL_OUTPUT",
    "result_size_is": "DM_RESULT_SIZE_IS",
    "output_from": "DM_OUTPUT_FROM",
    "minimum": "DM_MIN",
    "maximum": "DM_MAX",
    "sentinel": "DM_SENTINEL",
    "allow_nan": "DM_ALLOW_NAN",
    "allow_infinite": "DM_ALLOW_INFINITE",
}


def _macro_name(field: str) -> str:
    return _MACRO_NAMES.get(field, field)
