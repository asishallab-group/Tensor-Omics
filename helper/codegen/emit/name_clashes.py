"""An argument named like something its generated Fortran wrapper uses.

A dummy of a generated wrapper hides whatever else that wrapper calls by its name. An argument
called `count` beside an `n_selected_<x>` count makes the wrapper's `count(mask, kind=int32)` a
call of an integer, an argument called `size` does the same to `size(x, 1)`, one called `int32`
to every `arg_pos=3_int32`, one called `handled` to the prologue flag the wrapper declares --
and each is a compile error in generated code, which fails the whole library.

Python and the R shim are out of this: every name they bind beside the arguments starts with
`_` (`names`), which no Fortran name can. A Fortran wrapper cannot do that, so it is the
argument that is refused -- but only where it really collides. `count` is a perfectly good
argument of `mask_chunk_count`, whose C wrapper counts nothing, and main publishes it there.

"Where it collides" is read off the wrapper as it is emitted, not off a second list of what
each emitter might write: the wrappers are rendered, their macros expanded, and every dummy is
looked for in the roles that make a clash -- called (`call d(`, or `d(` for a scalar), used as a
kind (`3_d`, `kind=d`, `integer(d)`), imported (`use ..., only: d`), or declared a second
time as a local. A role the emitters start writing is covered without anyone remembering to
list it.
"""

from __future__ import annotations

import re

from ..abi.model import CBinding
from ..config import CONVENTIONS, Conventions
from ..diagnostics import DiagnosticBag
from .fortran_c import FortranCEmitter
from .fortran_wrapper import FortranWrapperEmitter, WrapperInfo

#: Fortran's own words where they are followed by `(` or name a type, which are syntax and not
#: a use of a name: `intent(in)` does not call anything called `intent`
_SYNTAX = frozenset({
    "intent", "dimension", "character", "integer", "real", "logical", "complex", "type",
    "bind", "if", "allocate", "nullify", "select", "case", "subroutine", "function",
    "kind", "len",
})

_STRING_RE = re.compile(r'"[^"\n]*"|\'[^\'\n]*\'')
_COMMENT_RE = re.compile(r"!.*")
_USE_RE = re.compile(r"^\s*use\b[^\n]*?\bonly\s*:([^\n]*)", re.I | re.M)
_SUBROUTINE_RE = re.compile(r"\bsubroutine\s+(\w+)\s*\(.*?\bend\s+subroutine\s+\1\b", re.I | re.S)


def check_wrapper_names(binding: CBinding, synthesis, macros, diagnostics: DiagnosticBag,
                        conventions: Conventions = CONVENTIONS,
                        macros_header: str = "src/macros.h") -> None:
    """Refuse every argument named like something its C wrapper or Fortran wrapper uses."""
    reported: set[tuple] = set()

    c_emitter = FortranCEmitter(conventions, macros_header)
    for module in binding:
        wrappers = {wrapper.name.lower(): wrapper for wrapper in module}
        text = _scannable(c_emitter.module(module), macros)
        header, subroutines = _split(text)
        for name, body in subroutines:
            wrapper = wrappers.get(name.lower())
            if wrapper is None:
                continue
            procedure = wrapper.procedure
            dummies = [
                (a.name, a.is_scalar, a.source if a.source is not None else procedure)
                for a in wrapper
            ]
            _check(dummies, body, header, f"C wrapper '{wrapper.name}'", procedure,
                   diagnostics, reported)

    project = synthesis.project
    emitter = FortranWrapperEmitter(conventions, macros_header=macros_header, project=project)
    info: dict[str, dict[str, WrapperInfo]] = {}
    for spec in synthesis.specs:
        for wrapper in (spec.validating, spec.allocating):
            if wrapper is not None:
                info.setdefault(spec.module_name, {})[wrapper.name.lower()] = WrapperInfo(
                    spec.impl.name, spec.mode_fix
                )
    for module in project:
        if module.name not in synthesis.generated_names or not module.procedures:
            continue
        text = _scannable(emitter.module(module, info.get(module.name)), macros)
        header, subroutines = _split(text)
        for name, body in subroutines:
            procedure = module.procedure(name)
            if procedure is None:
                continue
            dummies = [(a.name, a.is_scalar, a) for a in procedure.arguments]
            _check(dummies, body, header, f"Fortran wrapper '{procedure.name}'", procedure,
                   diagnostics, reported)


def _scannable(text: str, macros) -> str:
    """The generated text with strings and comments gone and the macros expanded."""
    text = _STRING_RE.sub('""', text)
    text = _COMMENT_RE.sub("", text)
    lines = [line for line in text.splitlines() if not line.lstrip().startswith("#")]
    text = "\n".join(lines)
    return macros.expand(text) if macros is not None else text


def _split(text: str) -> tuple[str, list[tuple[str, str]]]:
    """The module's specification part, and each subroutine's text by name."""
    head, _, _ = text.partition("contains")
    return head, [(m.group(1), m.group(0)) for m in _SUBROUTINE_RE.finditer(text)]


def _imports(text: str) -> set[str]:
    names = set()
    for match in _USE_RE.finditer(text):
        names |= {n.strip().lower() for n in match.group(1).split(",") if n.strip()}
    return names


def _role(name: str, is_scalar: bool, body: str, imported: set[str]) -> str | None:
    """How the wrapper uses `name` other than as the dummy, or None."""
    d = re.escape(name)
    if re.search(rf"\bcall\s+{d}\s*\(", body, re.I):
        return f"calls '{name}'"
    if name.lower() in imported:
        return f"imports '{name}'"
    if (re.search(rf"(?:(?<!\w)\d+|\.)_{d}\b", body, re.I)
            or re.search(rf"\bkind\s*=\s*{d}\b", body, re.I)
            or re.search(rf"\b(?:integer|real|complex|logical)\s*\(\s*{d}\s*\)", body, re.I)):
        return f"uses '{name}' as a kind"
    if is_scalar and name.lower() not in _SYNTAX and re.search(rf"(?<![%\w]){d}\s*\(", body, re.I):
        return f"calls '{name}(...)'"
    if len(re.findall(rf"::\s*{d}\s*$", body, re.I | re.M)) > 1:
        return f"declares a local called '{name}'"
    return None


def _check(dummies, body: str, header: str, where: str, procedure, diagnostics: DiagnosticBag,
           reported: set) -> None:
    imported = _imports(header) | _imports(body)
    for name, is_scalar, entity in dummies:
        if name.lower() == procedure.name.lower():
            continue  # `validate` refuses an argument named like its own procedure
        role = _role(name, is_scalar, body, imported)
        if role is None:
            continue
        location = getattr(entity, "location", None)
        key = (getattr(location, "file", None), getattr(location, "line", None), name.lower())
        if key in reported:
            continue
        reported.add(key)
        diagnostics.error(
            f"argument '{name}' is a name the generated {where} uses: it {role}",
            entity=entity,
            note=(
                "a dummy hides what else the wrapper calls by that name, and the wrapper "
                "would not compile; rename the argument"
            ),
        )
