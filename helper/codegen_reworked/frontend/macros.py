"""Access to the C preprocessor macros defined in `src/macros.h`.

Ford preprocesses the sources, so by the time the generator sees a Ford comment the
`DM_*` doc macros have already expanded to prose. To recognise them again the generator
needs the *expanded* form -- and hardcoding that prose would mean any reword of a macro
silently stops the generator from seeing it.

So patterns are derived from the macro definitions themselves: expand a template such as
`DM_DEFAULT((?P<value>.*))` through a preprocessor whose macro bodies have been
regex-escaped, and the result is a regex matching the expanded prose with the capture
groups in place. The macro and its pattern cannot drift apart.

This trick is carried over from the previous generator; the wrapper around it is not.
That one used pcpp's *command-line* class, which preprocesses the whole header to stdout
as a side effect of collecting the macros -- at import time, twice.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

from pcpp import Preprocessor


class MacroTable:
    """The macros defined by a header, with expansion of arbitrary snippets."""

    def __init__(self, header: Path, include_paths: tuple[Path, ...] = ()):
        self.header = Path(header)
        self.include_paths = tuple(include_paths)

    @cached_property
    def _plain(self) -> Preprocessor:
        return self._load(escape_bodies=False)

    @cached_property
    def _escaped(self) -> Preprocessor:
        # A separate instance: escaping rewrites the macro bodies in place.
        return self._load(escape_bodies=True)

    def _load(self, escape_bodies: bool) -> Preprocessor:
        preprocessor = Preprocessor()
        for path in self.include_paths:
            preprocessor.add_path(str(path))

        text = self.header.read_text()
        preprocessor.parse(text, str(self.header))
        # Drain the token stream so every #define is actually processed. Nothing is
        # written anywhere -- the tokens themselves are of no interest here.
        for _ in iter(preprocessor.token, None):
            pass

        if escape_bodies:
            for macro in preprocessor.macros.values():
                for token in macro.value:
                    token.value = re.escape(token.value)

        return preprocessor

    @property
    def names(self) -> frozenset[str]:
        return frozenset(self._plain.macros)

    def __contains__(self, name: str) -> bool:
        return name in self._plain.macros

    def expand(self, text: str) -> str:
        """Expand macro invocations in `text` the way the Fortran sources see them."""
        return self._expand_with(self._plain, text)

    def pattern(self, template: str) -> str:
        """Turn a macro invocation template into a regex matching its expansion.

        `template` is a macro invocation whose arguments are regex fragments, e.g.
        `DM_DEFAULT((?P<value>.*))`. Everything contributed by the macro body is
        escaped, so only the fragments passed in stay live.
        """
        return self._expand_with(self._escaped, template)

    def compiled(self, template: str, flags: int = 0) -> re.Pattern:
        return re.compile(self.pattern(template), flags)

    @staticmethod
    def _expand_with(preprocessor: Preprocessor, text: str) -> str:
        tokens = preprocessor.tokenize(text)
        return "".join(token.value for token in preprocessor.expand_macros(tokens))


@dataclass(frozen=True)
class DocMacroNames:
    """Names of the documentation macros the generator understands.

    Names live here rather than being spelled inline at each use, so that adding a doc
    macro is a change in one place. The wording of each macro stays in `src/macros.h`.
    """

    default: str = "DM_DEFAULT"
    required_if_mode: str = "DM_REQUIRED_IF_MODE"
    optional_output: str = "DM_OPTIONAL_OUTPUT"
    result_size_is: str = "DM_RESULT_SIZE_IS"
    output_from: str = "DM_OUTPUT_FROM"

    def all(self) -> tuple[str, ...]:
        return (
            self.default,
            self.required_if_mode,
            self.optional_output,
            self.result_size_is,
            self.output_from,
        )


DOC_MACROS = DocMacroNames()
