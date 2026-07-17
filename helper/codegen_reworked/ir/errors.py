"""The error code catalogue, read from `tox_errors`.

`tox_errors` is the single source of truth for what can go wrong and what it is called.
The generator reads its parameters and their documentation so the Python and R error
modules are derived from it rather than restated by hand -- the R mapping was maintained
separately, which is why `ERR_UNIT_NOT_CONNECTED` carries a comment about keeping its
value "for compatibility with existing R mapping".

Two things the catalogue must model:

- `ERR_*` codes are failures and raise in the interfacing languages. `STAT_*` codes are
  outcomes, not failures, and must never raise. None exist yet; the split is by prefix so
  the first one added needs no change here.
- An error code encodes the argument that caused it: `ierr = 10000*arg_pos + error`. So a
  raised error can name the offending argument, given the procedure's argument list.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from ..config import CONVENTIONS, Conventions
from ..diagnostics import DiagnosticBag
from .constants import ConstantError, ConstantEvaluator
from .doc import Doc
from .entities import Module

#: `create_err_code` in tox_errors packs the argument position in as
#: `10000*arg_pos + error`, and `get_err_code`/`get_err_arg_pos` unpack it again.
#: Mirrored here rather than parsed out of the Fortran body; a test guards the drift.
ARG_POS_FACTOR = 10000


class ErrorGroup(Enum):
    """The coarse grouping `tox_errors` documents through its numeric ranges.

    Used to give the interfacing languages a small hierarchy of exception types, rather
    than one type per code or a single opaque error.
    """

    IO = "io"
    INPUT = "input"
    MEMORY = "memory"
    RUNTIME = "runtime"
    INTERNAL = "internal"
    OTHER = "other"

    @classmethod
    def of(cls, value: int) -> ErrorGroup:
        for low, high, group in _GROUP_RANGES:
            if low <= value <= high:
                return group
        return cls.OTHER


#: Ranges as documented by the section headings in tox_errors
_GROUP_RANGES = (
    (100, 199, ErrorGroup.IO),
    (200, 299, ErrorGroup.INPUT),
    (300, 399, ErrorGroup.MEMORY),
    (5000, 5999, ErrorGroup.RUNTIME),
    (9000, 9999, ErrorGroup.INTERNAL),
)


@dataclass(frozen=True)
class ErrorCode:
    """One `ERR_*` or `STAT_*` parameter."""

    name: str
    value: int
    #: The parameter's documentation, which becomes the message text
    doc: Doc = Doc()
    is_status: bool = False

    @property
    def message(self) -> str:
        return self.doc.summary or self.name

    @property
    def group(self) -> ErrorGroup:
        return ErrorGroup.of(self.value)

    @property
    def is_ok(self) -> bool:
        return self.value == 0


@dataclass(frozen=True)
class DecodedError:
    """An `ierr` value taken apart."""

    code: int
    arg_pos: int
    error: ErrorCode | None = None

    @property
    def is_error(self) -> bool:
        """Mirrors `is_err`: the code, not the whole value, decides."""
        return self.code != 0

    @property
    def is_argument_related(self) -> bool:
        return self.arg_pos > 0


class ErrorCatalogue:
    """Every error and status code, keyed by name and by value."""

    def __init__(self, codes=(), conventions: Conventions = CONVENTIONS):
        self.conventions = conventions
        self.codes = tuple(sorted(codes, key=lambda c: c.value))
        self._by_name = {code.name.lower(): code for code in self.codes}
        self._by_value = {}
        for code in self.codes:
            self._by_value.setdefault(code.value, code)

    @classmethod
    def from_module(
        cls,
        module: Module,
        diagnostics: DiagnosticBag,
        constants: dict | None = None,
        conventions: Conventions = CONVENTIONS,
    ) -> ErrorCatalogue:
        """Read the catalogue out of the parsed `tox_errors` module."""
        evaluator = ConstantEvaluator(constants or {})
        codes = []

        for parameter in module.parameters:
            name = parameter.name
            is_error = name.startswith(conventions.error_code_prefix)
            is_status = name.startswith(conventions.status_code_prefix)
            if not (is_error or is_status):
                continue

            try:
                value = evaluator.evaluate(parameter.expression)
            except ConstantError as error:
                diagnostics.error(
                    f"error code '{name}' does not have a constant value: {error}",
                    entity=parameter,
                )
                continue

            if not isinstance(value, int) or isinstance(value, bool):
                diagnostics.error(
                    f"error code '{name}' is not an integer", entity=parameter
                )
                continue

            if not parameter.doc.summary:
                diagnostics.warn(
                    f"error code '{name}' has no documentation to use as its message",
                    entity=parameter,
                )

            codes.append(
                ErrorCode(name=name, value=value, doc=parameter.doc, is_status=is_status)
            )

        catalogue = cls(codes, conventions)
        catalogue._report_duplicates(diagnostics, module)
        return catalogue

    def _report_duplicates(self, diagnostics: DiagnosticBag, module: Module) -> None:
        seen: dict[int, ErrorCode] = {}
        for code in self.codes:
            if code.value in seen:
                diagnostics.error(
                    f"'{code.name}' and '{seen[code.value].name}' share the value "
                    f"{code.value}, so an error cannot be told from the other",
                    entity=module.parameter(code.name) or module,
                )
            else:
                seen[code.value] = code

    @property
    def errors(self) -> tuple[ErrorCode, ...]:
        """Failures. These raise in the interfacing languages."""
        return tuple(c for c in self.codes if not c.is_status and not c.is_ok)

    @property
    def statuses(self) -> tuple[ErrorCode, ...]:
        """Outcomes. These never raise."""
        return tuple(c for c in self.codes if c.is_status)

    @property
    def ok(self) -> ErrorCode | None:
        return self._by_name.get(self.conventions.ok_code.lower())

    def code(self, name: str) -> ErrorCode | None:
        return self._by_name.get(name.lower())

    def by_value(self, value: int) -> ErrorCode | None:
        return self._by_value.get(value)

    def decode(self, ierr: int) -> DecodedError:
        """Take an `ierr` apart the way `get_err_code`/`get_err_arg_pos` do."""
        code = ierr % ARG_POS_FACTOR
        arg_pos = ierr // ARG_POS_FACTOR
        return DecodedError(code=code, arg_pos=arg_pos, error=self.by_value(code))

    def encode(self, name: str, arg_pos: int = 0) -> int:
        """The inverse of `decode`, mirroring `create_err_code`. For tests."""
        code = self.code(name)
        if code is None:
            raise KeyError(f"no error code named '{name}'")
        return ARG_POS_FACTOR * arg_pos + code.value

    def groups(self) -> tuple[ErrorGroup, ...]:
        """The groups actually present, in a stable order."""
        present = {code.group for code in self.errors}
        return tuple(group for group in ErrorGroup if group in present)

    def __iter__(self):
        return iter(self.codes)

    def __len__(self) -> int:
        return len(self.codes)
