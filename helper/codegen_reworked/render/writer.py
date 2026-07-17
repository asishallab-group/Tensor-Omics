"""Text building with indentation.

Replaces the old `Indentable(str)` + `>>` operator. A str subclass that overloads `>>`
to mean "indent" reads as a shift on every use site and loses its type through any plain
string operation, so indentation depended on remembering to re-wrap. Here indentation is
owned by the writer and applied as lines are added, which is also why emitters can nest
without threading a level through every call.

Blank runs collapse by default. Emitting an empty section -- a wrapper with no type
conversions, say -- otherwise leaves a run of blank lines behind, which is exactly why
the current generated wrappers have ragged vertical gaps.
"""

from __future__ import annotations

from contextlib import contextmanager


class Writer:
    """Accumulates lines of generated code at a current indentation level."""

    def __init__(self, indent: str = "    ", level: int = 0):
        self._indent = indent
        self._level = level
        self._lines: list[str] = []

    @property
    def level(self) -> int:
        return self._level

    @property
    def prefix(self) -> str:
        return self._indent * self._level

    def line(self, text: str = "") -> None:
        """Write one line at the current indent.

        Multi-line text is accepted and behaves like `block`, so callers never have to
        care whether an interpolated fragment happened to contain a newline.
        """
        if text == "":
            self.blank()
        else:
            self.block(text)

    def block(self, text: str) -> None:
        """Write multi-line text, preserving its relative indentation."""
        for raw in text.split("\n"):
            stripped = raw.rstrip()
            if stripped:
                self._lines.append(self.prefix + stripped)
            else:
                self.blank()

    def lines(self, texts) -> None:
        for text in texts:
            self.line(text)

    def blank(self, collapse: bool = True) -> None:
        """Write a blank line.

        Collapses by default: consecutive blanks, and blanks at the very top, are
        dropped. This lets an emitter unconditionally separate sections without knowing
        whether the neighbouring section produced any output.
        """
        if collapse and (not self._lines or self._lines[-1] == ""):
            return
        self._lines.append("")

    def extend(self, other: Writer) -> None:
        """Splice another writer's output in at this writer's current indent."""
        self.block(other.render())

    @contextmanager
    def indent(self, levels: int = 1):
        self._level += levels
        try:
            yield self
        finally:
            self._level -= levels

    def render(self, trailing_newline: bool = False) -> str:
        """Return the accumulated text, without trailing blank lines."""
        lines = list(self._lines)
        while lines and lines[-1] == "":
            lines.pop()
        text = "\n".join(lines)
        if trailing_newline and text:
            text += "\n"
        return text

    def __str__(self) -> str:
        return self.render()

    def __len__(self) -> int:
        return len(self._lines)

    def __bool__(self) -> bool:
        return any(self._lines)
