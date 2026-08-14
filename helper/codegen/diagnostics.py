"""Diagnostics: how the generator reports problems in Fortran sources.

The generator's users are Fortran authors, not its own developers, so a problem in a
source file must read as a message about that file -- not as a Python traceback from
somewhere inside an emitter. Every check therefore produces a `Diagnostic` carrying the
source location and the entity chain it was found in.

Collecting into a `DiagnosticBag` rather than raising on the first problem means one run
reports every broken procedure, instead of making the author re-run once per mistake.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


class Severity(Enum):
    WARNING = "warning"
    ERROR = "error"


class CodegenError(Exception):
    """Raised when generation cannot proceed because sources contain errors."""

    def __init__(self, message: str, diagnostics: tuple[Diagnostic, ...] = ()):
        super().__init__(message)
        self.diagnostics = diagnostics


@dataclass(frozen=True)
class SourceLocation:
    """Where a diagnostic points. Rendered as a clickable `file:line`."""

    file: Path | None = None
    line: int | None = None

    def __str__(self) -> str:
        if self.file is None:
            return "<unknown>"
        if self.line is None:
            return str(self.file)
        return f"{self.file}:{self.line}"


@dataclass(frozen=True)
class Context:
    """The entity chain a diagnostic was found in, innermost first.

    Each element is a (kind, name) pair, e.g. ("argument", "vector"). Built from an IR
    entity by `Context.of`, which follows `parent` links -- the same chain the old
    generator printed, but kept as data so it can be asserted on in tests.

    A link with no name is skipped rather than rendered as `in project ''`. The Ford
    frontend names the project, but IR assembled by hand -- every fixture, and every test
    that builds a `Project` to get the parent links -- has no name to give and should not
    have to invent one to keep diagnostics readable.
    """

    entities: tuple[tuple[str, str], ...] = ()

    @classmethod
    def of(cls, entity) -> Context:
        chain: list[tuple[str, str]] = []
        current = entity
        while current is not None:
            kind = getattr(current, "entity_kind", None)
            name = getattr(current, "name", None)
            if kind is not None and name:
                chain.append((kind, str(name)))
            current = getattr(current, "parent", None)
        return cls(tuple(chain))

    def __str__(self) -> str:
        return " in ".join(f"{kind} '{name}'" for kind, name in self.entities)

    def __bool__(self) -> bool:
        return bool(self.entities)


@dataclass(frozen=True)
class Diagnostic:
    severity: Severity
    message: str
    location: SourceLocation = SourceLocation()
    context: Context = Context()
    note: str | None = None

    def render(self, color: bool = False) -> str:
        paint = _painter(color)
        severity_color = _RED if self.severity is Severity.ERROR else _YELLOW

        lines = [f"{paint(self.severity.value + ':', severity_color)} {self.message}"]
        if self.location.file is not None:
            lines.append(f"  --> {self.location}")
        if self.context:
            lines.append(f"  {paint(str(self.context), _DIM)}")
        if self.note:
            note_label = paint("note:", _CYAN)
            indented = str(self.note).strip("\n").replace("\n", "\n  ")
            lines.append(f"  {note_label} {indented}")
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.render(color=False)


class DiagnosticBag:
    """Collects diagnostics so one run reports every problem it can find."""

    def __init__(self) -> None:
        self._diagnostics: list[Diagnostic] = []

    def add(self, diagnostic: Diagnostic) -> Diagnostic:
        self._diagnostics.append(diagnostic)
        return diagnostic

    def error(self, message: str, entity=None, location: SourceLocation | None = None,
              note: str | None = None) -> Diagnostic:
        return self.add(self._make(Severity.ERROR, message, entity, location, note))

    def warn(self, message: str, entity=None, location: SourceLocation | None = None,
             note: str | None = None) -> Diagnostic:
        return self.add(self._make(Severity.WARNING, message, entity, location, note))

    @staticmethod
    def _make(severity: Severity, message: str, entity, location: SourceLocation | None,
              note: str | None) -> Diagnostic:
        if location is None:
            location = getattr(entity, "location", None) or SourceLocation()
        return Diagnostic(
            severity=severity,
            message=message,
            location=location,
            context=Context.of(entity) if entity is not None else Context(),
            note=note,
        )

    @property
    def diagnostics(self) -> tuple[Diagnostic, ...]:
        return tuple(self._diagnostics)

    @property
    def errors(self) -> tuple[Diagnostic, ...]:
        return tuple(d for d in self._diagnostics if d.severity is Severity.ERROR)

    @property
    def warnings(self) -> tuple[Diagnostic, ...]:
        return tuple(d for d in self._diagnostics if d.severity is Severity.WARNING)

    def extend(self, other: DiagnosticBag) -> None:
        self._diagnostics.extend(other._diagnostics)

    def render(self, color: bool = False) -> str:
        return "\n\n".join(d.render(color) for d in self._diagnostics)

    def raise_if_errors(self) -> None:
        """Abort generation if any error was collected.

        Warnings never abort: an undocumented procedure should not stop a build.
        """
        errors = self.errors
        if not errors:
            return
        plural = "s" if len(errors) > 1 else ""
        raise CodegenError(
            f"{len(errors)} error{plural} in Fortran sources:\n\n"
            + "\n\n".join(e.render(color=_color_enabled()) for e in errors),
            diagnostics=errors,
        )

    def __len__(self) -> int:
        return len(self._diagnostics)

    def __iter__(self):
        return iter(self._diagnostics)


_RED = "\033[38;5;196m"
_YELLOW = "\033[38;5;226m"
_CYAN = "\033[38;5;44m"
_DIM = "\033[38;5;245m"
_RESET = "\033[0m"


def _color_enabled() -> bool:
    if os.environ.get("NO_COLOR"):
        return False
    return sys.stderr.isatty()


def _painter(color: bool):
    if not color:
        return lambda text, _: text
    return lambda text, code: f"{code}{text}{_RESET}"
