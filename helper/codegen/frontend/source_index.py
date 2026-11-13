"""Finding where things are written, so a diagnostic can point at a line.

Ford does not keep the line an entity was declared on -- `num_lines` is a count, not a
position -- so the location has to be recovered. Recovering it by searching the *original*
file is not a workaround for that so much as the only correct option: Ford parses
preprocessed text, in which `#include <src/macros.h>` alone has injected some forty lines
before the first declaration. A line number from Ford would refer to a file nobody has.

So this searches the file the author actually edits, anchored on the entity's own name.
Every lookup returns None rather than a guess when it cannot find the declaration: a
diagnostic pointing at the wrong line is worse than one pointing at the file alone.
"""

from __future__ import annotations

import re
from pathlib import Path

from ..diagnostics import SourceLocation

#: Prefixes a procedure statement may carry: `pure subroutine`, `elemental impure
#: subroutine`, `pure integer(int32) function`, `recursive function`, and so on.
_PROCEDURE_PREFIX = r"(?:[\w()=*:,\s]*?\s)??"

#: Any subroutine or function statement, naming what it declares
_ANY_PROCEDURE_RE = re.compile(
    rf"^\s*{_PROCEDURE_PREFIX}\b(?:subroutine|function)\s+(?P<name>\w+)\s*[(\s]",
    re.IGNORECASE,
)

_END_RE = re.compile(r"^\s*end\s+(?:subroutine|function|module)\b", re.IGNORECASE)
_COMMENT_RE = re.compile(r"^\s*!")

#: `_PROCEDURE_PREFIX` is a lazy optional around a lazy star, so *rejecting* a line costs the
#: engine a retry at every prefix length -- 22 us on a long line against 2 us for a plain
#: search. Since the pattern cannot match a line that does not contain the word at all, a
#: substring test first skips that work on the ~99% of lines that are not declarations.
_PROCEDURE_WORDS = ("subroutine", "function")


def _word(name: str) -> str:
    return rf"(?<![\w]){re.escape(name)}(?![\w])"


class SourceFile:
    """One Fortran file, searched by entity name. Line numbers are 1-based."""

    def __init__(self, path: Path, text: str | None = None):
        self.path = Path(path)
        if text is None:
            text = self.path.read_text(errors="replace")
        self.lines = text.split("\n")
        #: `lines` never changes after this, so every lookup below is a constant of
        #: (file, name) and is answered once. The frontend asks for the same procedure
        #: once per argument it declares, which is where the repetition came from.
        self._procedure_lines: dict[str, int | None] = {}
        self._procedure_spans: dict[str, tuple[int, int] | None] = {}
        self._argument_lines: dict[tuple[str, str], int | None] = {}
        self._variable_lines: dict[str, int | None] = {}
        #: the line ranges of every procedure, from statement to `end`, found once on demand
        self._containers: list[tuple[int, int]] | None = None

    def location(self, line: int | None = None) -> SourceLocation:
        return SourceLocation(self.path, line)

    def module_line(self, name: str) -> int | None:
        pattern = re.compile(rf"^\s*module\s+{_word(name)}\s*$", re.IGNORECASE)
        return self._first_match(pattern)

    def procedure_line(self, name: str) -> int | None:
        """The line a subroutine or function statement declares `name` on.

        Anchored on `subroutine`/`function` immediately followed by the name, so a call
        to the procedure elsewhere in the file cannot be mistaken for its declaration.
        Where several statements declare it, the first that no other procedure holds: a
        procedure contained in another may share a module procedure's name, and the module
        procedure is the one asked about. Only where every one is contained is it the first.
        """
        if name in self._procedure_lines:
            return self._procedure_lines[name]

        pattern = self._statement_pattern(name)
        found = [number for number, line in self._statement_lines() if pattern.match(line)]
        outermost = [number for number in found if not self._is_contained(number)]
        line = (outermost or found or [None])[0]
        self._procedure_lines[name] = line
        return line

    @staticmethod
    def _statement_pattern(name: str) -> re.Pattern:
        return re.compile(
            rf"^\s*{_PROCEDURE_PREFIX}\b(?:subroutine|function)\s+{_word(name)}\s*[(\s]",
            re.IGNORECASE,
        )

    def _statement_lines(self):
        """The lines that may be a subroutine or function statement, with their numbers."""
        for number, line in enumerate(self.lines, start=1):
            lowered = line.lower()
            if not any(word in lowered for word in _PROCEDURE_WORDS):
                continue
            if _COMMENT_RE.match(line) or _END_RE.match(line):
                continue
            yield number, line

    def _is_contained(self, number: int) -> bool:
        """Whether line `number` lies inside a procedure whose statement is above it."""
        if self._containers is None:
            self._containers = [
                (start, self._end_of(start, match.group("name")))
                for start, line in self._statement_lines()
                if (match := _ANY_PROCEDURE_RE.match(line)) is not None
            ]
        return any(start < number <= end for start, end in self._containers)

    def _end_of(self, start: int, name: str) -> int:
        """The line of the `end` that closes the procedure `name` declared on line `start`."""
        end_pattern = re.compile(
            rf"^\s*end\s+(?:subroutine|function)(?:\s+{_word(name)})?\s*$", re.IGNORECASE
        )
        for number in range(start + 1, len(self.lines) + 1):
            if end_pattern.match(self.lines[number - 1]):
                return number
        return len(self.lines)

    def procedure_prefix(self, name: str) -> str | None:
        """What the statement of `name` says before `subroutine`/`function`, as written.

        `M_OMP_PURE`, `pure`, `impure elemental`, or "" for nothing. None where the statement
        is not found, or where it is continued from the line above (`pure &`), whose prefix
        this line does not hold: a rule reading the prefix must not take a part for the whole.
        """
        number = self.procedure_line(name)
        if number is None or self.statement_start(number) != number:
            return None
        match = re.match(
            rf"^\s*(?P<prefix>{_PROCEDURE_PREFIX})\b(?:subroutine|function)\s+{_word(name)}",
            self.lines[number - 1],
            re.IGNORECASE,
        )
        return None if match is None else match.group("prefix").strip()

    def statement_start(self, number: int) -> int:
        """The first line of the statement on line `number`: the earliest of the lines above it,
        each ending in `&`, that it continues. Blank and comment lines in between are skipped."""
        start = number
        for previous in range(number - 1, 0, -1):
            line = self.lines[previous - 1]
            if not line.strip() or _COMMENT_RE.match(line):
                continue
            if not line.split("!")[0].rstrip().endswith("&"):
                break
            start = previous
        return start

    def module_span(self, name: str) -> tuple[int, int] | None:
        """The inclusive line range of a module, from its statement to its `end`."""
        start = self.module_line(name)
        if start is None:
            return None
        end_pattern = re.compile(
            rf"^\s*end\s+module(?:\s+{_word(name)})?\s*$", re.IGNORECASE
        )
        for number in range(start + 1, len(self.lines) + 1):
            if end_pattern.match(self.lines[number - 1]):
                return start, number
        return start, len(self.lines)

    def macro_lines(self, macro: str, span: tuple[int, int] | None = None) -> list[int]:
        """The lines within `span` whose code names `macro`: not a comment, not a directive."""
        pattern = re.compile(_word(macro))
        found = []
        for number in self._numbers(span):
            line = self.lines[number - 1]
            if _COMMENT_RE.match(line) or line.lstrip().startswith("#"):
                continue
            if pattern.search(line.split("!")[0]):
                found.append(number)
        return found

    def comment_macro_lines(self, macro: str, span: tuple[int, int] | None = None) -> list[int]:
        """The lines within `span` that open with `!<macro>(`: a directive a macro writes."""
        pattern = re.compile(rf"^\s*!{re.escape(macro)}\s*\(")
        return [number for number in self._numbers(span) if pattern.match(self.lines[number - 1])]

    def _numbers(self, span: tuple[int, int] | None) -> range:
        start, end = span if span is not None else (1, len(self.lines))
        return range(start, end + 1)

    def argument_line(self, procedure: str, argument: str) -> int | None:
        """The line declaring `argument` inside `procedure`.

        A declaration is a line with `::` whose right-hand side names the argument, so
        `intent(in) :: n_dims` matches while a use of `n_dims` in the body does not.
        """
        key = (procedure, argument)
        if key in self._argument_lines:
            return self._argument_lines[key]

        span = self.procedure_span(procedure)
        if span is None:
            self._argument_lines[key] = None
            return None

        start, end = span
        pattern = re.compile(_word(argument), re.IGNORECASE)
        for number in range(start + 1, end + 1):
            line = self.lines[number - 1]
            if _COMMENT_RE.match(line):
                continue
            _, separator, declared = line.partition("::")
            if not separator:
                continue
            declared = declared.split("!")[0]
            if pattern.search(declared):
                self._argument_lines[key] = number
                return number
        self._argument_lines[key] = None
        return None

    def argument_doc_line(self, procedure: str, argument: str) -> int | None:
        """Where an argument's `!!` documentation starts: the line after its declaration."""
        declaration = self.argument_line(procedure, argument)
        return None if declaration is None else declaration + 1

    def procedure_span(self, name: str) -> tuple[int, int] | None:
        """The inclusive line range of a procedure, from its statement to its `end`."""
        if name in self._procedure_spans:
            return self._procedure_spans[name]

        start = self.procedure_line(name)
        if start is None:
            self._procedure_spans[name] = None
            return None

        span = (start, self._end_of(start, name))
        self._procedure_spans[name] = span
        return span

    def variable_line(self, name: str) -> int | None:
        """A module-level declaration, e.g. `integer, parameter :: ERR_OK = 0`."""
        if name in self._variable_lines:
            return self._variable_lines[name]

        pattern = re.compile(_word(name), re.IGNORECASE)
        for number, line in enumerate(self.lines, start=1):
            if _COMMENT_RE.match(line):
                continue
            _, separator, declared = line.partition("::")
            if not separator:
                continue
            if pattern.search(declared.split("!")[0]):
                self._variable_lines[name] = number
                return number
        self._variable_lines[name] = None
        return None

    def variable_doc_line(self, name: str) -> int | None:
        declaration = self.variable_line(name)
        return None if declaration is None else declaration + 1

    def _first_match(self, pattern: re.Pattern) -> int | None:
        for number, line in enumerate(self.lines, start=1):
            if pattern.match(line):
                return number
        return None


class SourceIndex:
    """Caches `SourceFile`s, since every entity of a module asks about the same file."""

    def __init__(self, root: Path = Path()):
        self.root = Path(root)
        #: Keyed on the path as asked for. An `lru_cache` on the method would key on
        #: `self` too and keep every index -- and every file it read -- alive for the
        #: process, which is why these are instance dicts.
        self._files: dict[Path, SourceFile | None] = {}
        self._relatives: dict[Path, Path] = {}
        self._root_resolved: Path | None = None

    def file(self, path: Path) -> SourceFile | None:
        if path in self._files:
            return self._files[path]
        try:
            source = SourceFile(path, self._absolute(path).read_text(errors="replace"))
        except OSError:
            source = None
        self._files[path] = source
        return source

    def _absolute(self, path: Path) -> Path:
        """Where `path` is on disk. A relative one is relative to the root, not the cwd."""
        path = Path(path)
        return path if path.is_absolute() else self.root / path

    def module(self, path: Path, name: str) -> SourceLocation:
        return self._locate(path, lambda f: f.module_line(name))

    def procedure(self, path: Path, name: str) -> SourceLocation:
        return self._locate(path, lambda f: f.procedure_line(name))

    def procedure_prefix(self, path: Path, name: str) -> str | None:
        source = self.file(path)
        return None if source is None else source.procedure_prefix(name)

    def argument(self, path: Path, procedure: str, name: str) -> SourceLocation:
        return self._locate(path, lambda f: f.argument_line(procedure, name))

    def argument_doc(self, path: Path, procedure: str, name: str) -> int | None:
        source = self.file(path)
        return None if source is None else source.argument_doc_line(procedure, name)

    def variable(self, path: Path, name: str) -> SourceLocation:
        return self._locate(path, lambda f: f.variable_line(name))

    def variable_doc(self, path: Path, name: str) -> int | None:
        source = self.file(path)
        return None if source is None else source.variable_doc_line(name)

    def at(self, path: Path, line: int | None) -> SourceLocation:
        """A location in `path` at a line found some other way, the path as the author writes it."""
        return SourceLocation(self._relative(path), line)

    def _locate(self, path: Path, find) -> SourceLocation:
        source = self.file(path)
        if source is None:
            return SourceLocation(self._relative(path))
        return SourceLocation(self._relative(path), find(source))

    def _relative(self, path: Path) -> Path:
        """`path` as the author would write it, relative to the repository root.

        Ford is handed absolute source directories -- it has to be, since the generator
        chdirs -- so every path it reports back is absolute, and a diagnostic quoting one
        in full is both unreadable and different on every machine. A path outside the root
        is left alone rather than turned into a chain of `..`.

        Memoised, and the root resolved once: this is asked for every module, procedure,
        argument and parameter -- thousands of calls for a few dozen distinct answers --
        and `Path.resolve()` is a walk of `lstat` calls, not a string operation.
        """
        if path in self._relatives:
            return self._relatives[path]

        if self._root_resolved is None:
            try:
                self._root_resolved = self.root.resolve()
            except OSError:
                self._root_resolved = self.root

        try:
            relative = self._absolute(path).resolve().relative_to(self._root_resolved)
        except (ValueError, OSError):
            relative = Path(path)
        self._relatives[path] = relative
        return relative
