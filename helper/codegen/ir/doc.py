"""The documentation tree: Ford comments parsed into text, links and tables.

This layer knows Ford *syntax* and nothing else. It does not know that a two-column
table headed `Mode | Value` means something, or that certain prose encodes a default
value -- those are meanings, and they live in `directives.py` and `roles.py`. The
previous generator interpreted mode tables inside the table class itself, which put
knowledge of the argument model into the markdown parser.

Everything here is immutable. Transformations return new nodes, so an emitter that
rewrites a table (the C wrapper replaces mode parameters with their strings) cannot
disturb the doc every other emitter reads.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from typing import Iterator, Sequence


class DocParseError(Exception):
    """Malformed Ford documentation. Carries a relative line for the frontend to place."""

    def __init__(self, message: str, line_number: int | None = None, note: str | None = None):
        super().__init__(message)
        self.line_number = line_number
        self.note = note


_NAME = "[A-Za-z][A-Za-z_0-9]*"

#: The two type vocabularies are **not interchangeable**, and Ford says so: give a component an
#: item type (or the reverse) and it warns and generates no link. So they are spelled out
#: separately, from Ford's own documentation of `[[component(type):item(type)]]`, rather than
#: shared. Both lists were wrong until 2026-08-12 -- each carried two of the other's types and
#: was missing `interface`/`absinterface` -- and the way they were wrong was the bad way: a type
#: not in the list makes the *whole link* fail to match, so it was never a link at all. Not
#: resolved for Python or R, not checked by `validate.check_doc_links`, and rendered into every
#: binding as the raw `[[...]]` text the author typed.
#:
#: Longest alternative first where one is a prefix of another (`procedure` before `proc`): the
#: engine backtracks correctly either way, but the intent should not depend on that.
_COMPONENT_TYPES = (
    "procedure|proc|subroutine|function|absinterface|interface|block|type|file|"
    "submodule|module|program|namelist"
)
_ITEM_TYPES = (
    "absinterface|interface|bound|common|constructor|final|function|modproc|"
    "subroutine|type|variable"
)
#: A component type may carry an `ext` prefix, naming an entity in an *external* project.
#: Nothing in this project can resolve one, and `validate.check_doc_links` says so rather than
#: reporting it missing.
EXTERNAL_PREFIX = "ext"
_COMPONENT_TYPES = f"(?:{EXTERNAL_PREFIX})?(?:{_COMPONENT_TYPES})"


def _link_part(name: str, types: str) -> str:
    return rf"(?P<{name}>{_NAME})\s*(?:\(\s*(?P<{name}_type>{types})\s*\)\s*)?"


@dataclass(frozen=True)
class FordLink:
    """A Ford cross-reference, e.g. `[[tox_clustering(module):METHOD_WARD(variable)]]`."""

    RE = re.compile(
        rf"\[\[\s*{_link_part('component', _COMPONENT_TYPES)}"
        rf"(?::\s*{_link_part('item', _ITEM_TYPES)})?\]\]"
    )

    component: str
    component_type: str | None = None
    item: str | None = None
    item_type: str | None = None

    @classmethod
    def _from_match(cls, match: re.Match) -> FordLink:
        return cls(
            component=match.group("component"),
            component_type=match.group("component_type"),
            item=match.group("item"),
            item_type=match.group("item_type"),
        )

    @property
    def target(self) -> str:
        """The linked name: the item if there is one, otherwise the component."""
        return self.item if self.item is not None else self.component

    def __str__(self) -> str:
        text = self.component
        if self.component_type is not None:
            text += f"({self.component_type})"
        if self.item is not None:
            text += f":{self.item}"
            if self.item_type is not None:
                text += f"({self.item_type})"
        return f"[[{text}]]"


@dataclass(frozen=True)
class Text:
    """A run of plain prose."""

    text: str

    def __str__(self) -> str:
        return self.text


Span = Text | FordLink

#: A run of backticks, which opens an inline-code span and is closed by a run of the same length
_BACKTICK_RUN = re.compile(r"`+")


def _code_spans(text: str) -> list[tuple[int, int]]:
    """The `(start, end)` of every inline-code span, as Markdown delimits them.

    A run of *n* backticks opens a span that the next run of exactly *n* closes; an unclosed
    run is not a span at all. Ford has left the contents of these verbatim since v7, markup
    and cross-references alike, so nothing inside one is this parser's business either.
    """
    runs = list(_BACKTICK_RUN.finditer(text))
    spans: list[tuple[int, int]] = []
    index = 0
    while index < len(runs):
        opener = runs[index]
        for offset in range(index + 1, len(runs)):
            if runs[offset].group() == opener.group():
                spans.append((opener.start(), runs[offset].end()))
                index = offset
                break
        index += 1
    return spans


@dataclass(frozen=True)
class DocLine:
    """One line of documentation, split into prose runs and links."""

    spans: tuple[Span, ...] = ()
    line_number: int | None = None

    @classmethod
    def parse(cls, text: str, line_number: int | None = None) -> DocLine:
        text = text.strip()
        spans: list[Span] = []
        position = 0
        code = _code_spans(text)

        for match in FordLink.RE.finditer(text):
            start, end = match.span()
            # Ford 7 leaves a link inside inline code verbatim, like any other markup: an
            # author writing `` `call [[my_sub]]` `` means to *show* the link syntax. Parsing
            # it here would resolve something Ford does not, and -- since a dangling link is
            # an error now -- could refuse a build over prose that was never a link.
            if any(s <= start and end <= e for s, e in code):
                continue
            if start > position:
                spans.append(Text(text[position:start]))
            spans.append(FordLink._from_match(match))
            position = end

        if position < len(text):
            spans.append(Text(text[position:]))

        return cls(tuple(spans), line_number)

    @property
    def text(self) -> str:
        return "".join(str(span) for span in self.spans)

    @property
    def links(self) -> tuple[FordLink, ...]:
        return tuple(span for span in self.spans if isinstance(span, FordLink))

    @property
    def is_blank(self) -> bool:
        return not self.text.strip()

    def with_spans(self, spans: Sequence[Span]) -> DocLine:
        return DocLine(tuple(spans), self.line_number)

    def __str__(self) -> str:
        return self.text


class Alignment(Enum):
    """A markdown column alignment.

    DEFAULT (`---`) is kept distinct from LEFT (`:---`) so a table survives a
    parse/render round trip unchanged.
    """

    DEFAULT = "default"
    LEFT = "left"
    RIGHT = "right"
    CENTER = "center"

    def marker(self, width: int = 3) -> str:
        """Render back to markdown, padded to `width`."""
        left = ":" if self in (Alignment.LEFT, Alignment.CENTER) else ""
        right = ":" if self in (Alignment.RIGHT, Alignment.CENTER) else ""
        fill = max(width - len(left) - len(right), 1)
        return f"{left}{'-' * fill}{right}"

    @classmethod
    def parse(cls, text: str) -> Alignment:
        text = text.strip()
        if _ALIGNMENT_CELL_RE.match(text) is None:
            raise DocParseError(f"'{text}' is not a markdown table alignment")
        left = text.startswith(":")
        right = text.endswith(":")
        if left and right:
            return cls.CENTER
        if right:
            return cls.RIGHT
        if left:
            return cls.LEFT
        return cls.DEFAULT


_ALIGNMENT_CELL_RE = re.compile(r":?-+:?\Z")


@dataclass(frozen=True)
class DocTable:
    """A markdown table in a Ford comment."""

    header: tuple[DocLine, ...]
    alignment: tuple[Alignment, ...]
    rows: tuple[tuple[DocLine, ...], ...] = ()
    line_number: int | None = None

    def __post_init__(self):
        if not self.header:
            raise DocParseError("table has no columns", self.line_number)

        expected = len(self.header)
        if len(self.alignment) != expected:
            raise DocParseError(
                f"table has {expected} header columns but "
                f"{len(self.alignment)} alignment columns",
                self.line_number,
            )
        for index, row in enumerate(self.rows):
            if len(row) != expected:
                raise DocParseError(
                    f"table row {index + 1} has {len(row)} columns, expected {expected}",
                    self.line_number,
                    note="every row needs the same number of '|' separated cells",
                )

    @property
    def n_columns(self) -> int:
        return len(self.header)

    @property
    def header_text(self) -> tuple[str, ...]:
        return tuple(cell.text for cell in self.header)

    def column(self, index: int) -> tuple[DocLine, ...]:
        return tuple(row[index] for row in self.rows)

    def with_rows(self, rows) -> DocTable:
        return DocTable(self.header, self.alignment, tuple(rows), self.line_number)

    @property
    def is_blank(self) -> bool:
        return False


Block = DocLine | DocTable


@dataclass(frozen=True)
class Doc:
    """A parsed Ford comment: a sequence of lines and tables."""

    blocks: tuple[Block, ...] = ()

    @classmethod
    def parse(cls, lines: Sequence[str], first_line_number: int | None = None) -> Doc:
        """Parse raw documentation lines.

        `first_line_number`, when given, is the file line of `lines[0]`, so diagnostics
        can point at the offending line rather than at the procedure.
        """
        lines = list(lines)
        while lines and not lines[-1].strip():
            lines.pop()

        def number_of(index: int) -> int | None:
            return None if first_line_number is None else first_line_number + index

        blocks: list[Block] = []
        index = 0
        while index < len(lines):
            table, consumed = _parse_table(lines, index, number_of(index))
            if table is not None:
                blocks.append(table)
                index = consumed
            else:
                blocks.append(DocLine.parse(lines[index], number_of(index)))
                index += 1

        return cls(tuple(blocks))

    @property
    def lines(self) -> tuple[DocLine, ...]:
        return tuple(b for b in self.blocks if isinstance(b, DocLine))

    @property
    def tables(self) -> tuple[DocTable, ...]:
        return tuple(b for b in self.blocks if isinstance(b, DocTable))

    @property
    def links(self) -> tuple[FordLink, ...]:
        return tuple(link for line in self.all_lines for link in line.links)

    @property
    def all_lines(self) -> tuple[DocLine, ...]:
        """Every DocLine, including the cells of tables."""
        collected: list[DocLine] = []
        for block in self.blocks:
            if isinstance(block, DocTable):
                collected.extend(block.header)
                for row in block.rows:
                    collected.extend(row)
            else:
                collected.append(block)
        return tuple(collected)

    @property
    def text(self) -> str:
        """The prose, one block per line. Tables render as their source-ish form."""
        return "\n".join(
            block.text if isinstance(block, DocLine) else _table_text(block)
            for block in self.blocks
        )

    @property
    def summary(self) -> str:
        """The first non-blank line, for a one-line description."""
        for line in self.lines:
            if not line.is_blank:
                return line.text
        return ""

    def with_blocks(self, blocks: Sequence[Block]) -> Doc:
        return Doc(tuple(blocks))

    def replace(self, old: Block, new: Block) -> Doc:
        return Doc(tuple(new if block is old else block for block in self.blocks))

    def without_line(self, line_number: int | None) -> Doc:
        """A copy without the block at `line_number`.

        For prose one of the generator's own `DM_` macros wrote and a wrapper has since
        resolved: a mode-split wrapper *is* its mode, so the `DM_REQUIRED_IF_MODE` line it
        inherited states a condition that no longer exists, on an argument the split has
        already made mandatory. Only macro-written lines are ever removed this way -- author
        prose is never rewritten. `None` is a no-op, so a doc built without line numbers
        (every hand-made test fixture) is unaffected.
        """
        if line_number is None:
            return self
        return Doc(tuple(b for b in self.blocks if b.line_number != line_number))

    def __bool__(self) -> bool:
        return any(not block.is_blank for block in self.blocks)

    def __len__(self) -> int:
        return len(self.blocks)

    def __iter__(self) -> Iterator[Block]:
        return iter(self.blocks)


def _table_text(table: DocTable) -> str:
    """Render a table back to markdown, alignment row included.

    Faithful rather than pretty: an emitter that re-emits documentation must not turn a
    table into something that no longer parses as one. Column-width alignment is a
    presentation concern and belongs to the emitters.
    """
    def render_row(cells) -> str:
        return "| " + " | ".join(cell.text for cell in cells) + " |"

    alignment = "|" + "|".join(a.marker() for a in table.alignment) + "|"
    return "\n".join((render_row(table.header), alignment, *map(render_row, table.rows)))


def _split_row(line: str) -> list[str] | None:
    """Split `| a | b |` into cells, or return None if the line is not a table row."""
    line = line.strip()
    if len(line) < 2 or not line.startswith("|") or not line.endswith("|"):
        return None
    return [cell.strip() for cell in line[1:-1].split("|")]


def _is_alignment_row(cells: Sequence[str] | None) -> bool:
    return bool(cells) and all(_ALIGNMENT_CELL_RE.match(cell) for cell in cells)


def _parse_table(lines: Sequence[str], start: int, line_number: int | None):
    """Try to read a table at `start`. Returns (table, next_index) or (None, start).

    A row-shaped line is only a table when an alignment row follows it. Anything else is
    left as prose rather than reported as an error -- a lone `| note |` in a comment is
    not a mistake, and failing generation over it would be hostile.
    """
    header_cells = _split_row(lines[start])
    if header_cells is None:
        return None, start

    if start + 1 >= len(lines):
        return None, start

    alignment_cells = _split_row(lines[start + 1])
    if not _is_alignment_row(alignment_cells):
        return None, start

    header = tuple(DocLine.parse(cell, line_number) for cell in header_cells)
    alignment = tuple(Alignment.parse(cell) for cell in alignment_cells)

    rows: list[tuple[DocLine, ...]] = []
    index = start + 2
    while index < len(lines):
        cells = _split_row(lines[index])
        if cells is None:
            break
        row_number = None if line_number is None else line_number + (index - start)
        rows.append(tuple(DocLine.parse(cell, row_number) for cell in cells))
        index += 1

    return DocTable(header, alignment, tuple(rows), line_number), index
