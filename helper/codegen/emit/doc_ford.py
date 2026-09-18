"""Rendering a `Doc` back into a Ford comment.

The C wrapper inherits its documentation from the procedure it wraps, so the doc tree has
to go back out as Ford comments. Tables are rendered with their columns padded, which is
cosmetic but makes the generated source readable -- and readable generated source is what
lets a reviewer check the generator rather than trust it.
"""

from __future__ import annotations

from ..ir.doc import Doc, DocLine, DocTable
from ..render import Writer

#: How each kind of entity marks its documentation
MARKERS = {
    # `!!` documents the entity it follows
    "argument": ("!! ", "!! "),
    # `!>` opens a block before the entity, `!|` continues it
    "procedure": ("!> ", "!| "),
    "module": ("!> ", "!| "),
}


def render_doc(doc: Doc, kind: str, summary: str | None = None) -> str:
    """Render `doc` as a Ford comment.

    `summary` is prepended as a `summary:` meta tag, which is how a generated wrapper
    says what it is before repeating what it wraps.
    """
    opener, continuation = MARKERS[kind]
    writer = Writer()

    first = True
    if summary is not None:
        writer.line(f"{opener}summary: {summary}")
        first = False

    for block in doc:
        marker = opener if first else continuation
        first = False
        if isinstance(block, DocTable):
            for line in render_table(block).split("\n"):
                writer.line(f"{continuation}{line}")
        else:
            text = block.text
            writer.line(f"{marker}{text}" if text else continuation.rstrip())

    return writer.render()


def render_table(table: DocTable) -> str:
    """Render a table with padded columns, as a human would write it."""
    rows = [tuple(cell.text for cell in table.header)]
    rows.extend(tuple(cell.text for cell in row) for row in table.rows)

    widths = [
        max(len(row[column]) for row in rows) for column in range(table.n_columns)
    ]

    writer = Writer()
    writer.line("| " + " | ".join(t.ljust(w) for t, w in zip(rows[0], widths)) + " |")
    writer.line(
        "|" + "|".join(a.marker(w + 2) for a, w in zip(table.alignment, widths)) + "|"
    )
    for row in rows[1:]:
        writer.line("| " + " | ".join(t.ljust(w) for t, w in zip(row, widths)) + " |")
    return writer.render()
