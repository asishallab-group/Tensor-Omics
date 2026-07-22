"""Emitting VS Code snippets for the interface.

A convenience layer over the same model the C/Python/R emitters read, so the snippets a
developer expands cannot drift from the functions that actually exist.

Output is split by language and by module root into `snippets/<language>_<root>_snippets.json`
(`fortran_f42_snippets.json`, `python_tox_snippets.json`, ...). The language is in the file
name rather than a per-snippet `scope`, and the root (`f42`/`tox`, taken from the module
name up to its first underscore) keeps the two namespaces in separate files.

What is emitted, and why:

- **Calling a procedure.** One snippet per exported procedure, in each language. Fortran
  gets the native `call` (or `result = f(...)` for a function), plus a variant that guards
  `ierr`; Python and R get the generated wrapper call. Arguments are keyword tabstops, so
  the template is self-documenting; a `mode`/`method` argument becomes a choice of its
  accepted values (the string API in Python/R, the Fortran parameters in Fortran).
- **Setting a module up.** Python gets a per-module import with a choice of that module's
  functions. R has no per-module unit -- everything is sourced through one loader -- so it
  gets a single global setup snippet rather than one per module.
- **A handful of generic aids**: a Fortran `use ..., only:` per module, an error-handling
  wrapper per interfacing language, and an error-code picker built from the catalogue.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

from ..abi.model import CArgument, CWrapper, CWrapperModule, CInterface
from ..ir.errors import ErrorCatalogue
from .python_ctypes import PythonEmitter, python_literal
from .r_wrapper import r_literal

#: The three languages, used internally as bucket keys
FORTRAN, PYTHON, R = "fortran", "python", "r"

#: The capitalised name each language takes in its file name, matching the existing
#: `snippets/` convention (`Fortran_f42_snippets.json`, ...)
_FILE_LANGUAGE = {FORTRAN: "Fortran", PYTHON: "Python", R: "R"}

#: The R loader a user sources before calling any wrapper, relative to the repo root
R_LOADER = "rcpp/load_tensor_omics.R"


@dataclass(frozen=True)
class _Style:
    """How one language lays a multi-line keyword call out."""

    sep: str  #: between an argument name and its value
    indent: str  #: one level of indentation
    cont: str  #: line continuation, appended to every line but the last
    trailing_comma: bool  #: whether the last argument keeps its comma
    close_on_own_line: bool  #: whether `)` sits on a line of its own


_STYLES = {
    FORTRAN: _Style(sep="=", indent="    ", cont=" &", trailing_comma=False,
                    close_on_own_line=False),
    PYTHON: _Style(sep="=", indent="    ", cont="", trailing_comma=True,
                   close_on_own_line=True),
    R: _Style(sep=" = ", indent="  ", cont="", trailing_comma=False,
              close_on_own_line=True),
}


class SnippetEmitter:
    def __init__(self):
        # The user-facing parameter list is the one Python asks for; the emitter-parity
        # test guarantees R asks for the very same set, so either serves the snippets.
        self._python = PythonEmitter()

    def snippets_files(self, interface: CInterface,
                       catalogue: ErrorCatalogue) -> dict[str, str]:
        """The generated files, keyed by bare file name (`fortran_tox_snippets.json`)."""
        buckets: dict[tuple[str, str], dict] = {}
        for module in interface:
            self._module_snippets(module, buckets)
        self._generic_snippets(catalogue, buckets)
        return {
            f"{_FILE_LANGUAGE[language]}_{root}_snippets.json":
                json.dumps(snippets, indent=2) + "\n"
            for (language, root), snippets in buckets.items()
            if snippets
        }

    @staticmethod
    def _add(buckets: dict, language: str, root: str, name: str,
             prefix: str, description: str, body: list[str]) -> None:
        bucket = buckets.setdefault((language, root), {})
        bucket[name] = {"prefix": prefix, "description": description, "body": body}

    # -- per module -------------------------------------------------------------

    def _module_snippets(self, module: CWrapperModule, buckets: dict) -> None:
        root, rest = _split_module(module.stripped_name)
        summary = module.doc.summary or f"the {module.stripped_name} module"

        names = [w.procedure.name for w in module]
        self._add(
            buckets, FORTRAN, root, f"{module.stripped_name} (use)",
            f"{root}:{rest}.use", f"use {module.stripped_name}, only: ...",
            [f"use {module.stripped_name}, only: {_choice(1, names)}"],
        )

        stripped = [w.stripped_name for w in module]
        self._add(
            buckets, PYTHON, root, f"{module.stripped_name} (import)",
            f"{root}:{rest}.setup", f"Import a function from {summary}",
            [
                f"from tensor_omics.{module.stripped_name} import (",
                f"    {_choice(1, stripped)},",
                ")",
            ],
        )

        for wrapper in module:
            self._call_snippets(wrapper, root, buckets)

    def _call_snippets(self, wrapper: CWrapper, root: str, buckets: dict) -> None:
        summary = wrapper.doc.summary or wrapper.stripped_name

        self._add(
            buckets, FORTRAN, root, f"{wrapper.procedure.name} (call)",
            f"{root}:{wrapper.procedure.name}", summary,
            self._fortran_call(wrapper),
        )
        checked = self._fortran_call(wrapper, checked=True)
        if checked is not None:
            self._add(
                buckets, FORTRAN, root, f"{wrapper.procedure.name} (call, checked)",
                f"{root}:{wrapper.procedure.name}.checked",
                f"{summary} (with error check)", checked,
            )

        self._add(
            buckets, PYTHON, root, f"{wrapper.stripped_name} (call)",
            f"{root}:{wrapper.stripped_name}", summary,
            self._wrapper_call(wrapper, PYTHON),
        )
        self._add(
            buckets, R, root, f"{wrapper.stripped_name} (call)",
            f"{root}:{wrapper.stripped_name}", summary,
            self._wrapper_call(wrapper, R),
        )

    # -- Fortran ----------------------------------------------------------------

    def _fortran_call(self, wrapper: CWrapper, checked: bool = False):
        """The native Fortran call. `checked` appends an `ierr` guard, or None if there
        is no error argument to guard on."""
        procedure = wrapper.procedure
        counter = _Counter()

        if procedure.result is not None:
            head = f"{_ph(counter.next(), 'result')} = {procedure.name}("
        else:
            head = f"call {procedure.name}("

        pairs = [(arg.name, self._fortran_value(arg, wrapper, counter))
                 for arg in procedure.arguments]
        lines = _render_call(head, pairs, _STYLES[FORTRAN])

        if not checked:
            return lines
        error = procedure.argument(procedure.conventions.error_arg)
        if error is None:
            return None
        lines.append(f"if (is_err({error.name})) return")
        return lines

    def _fortran_value(self, argument, wrapper: CWrapper, counter: _Counter) -> str:
        c_argument = wrapper.argument(argument.name)
        if c_argument is not None and c_argument.mode is not None:
            return _choice(counter.next(),
                           [v.parameter for v in c_argument.mode.values])
        return _ph(counter.next(), argument.name)

    # -- Python / R -------------------------------------------------------------

    def _wrapper_call(self, wrapper: CWrapper, language: str) -> list[str]:
        counter = _Counter()
        assign = self._python._outputs(wrapper)
        callee = wrapper.stripped_name

        if assign:
            operator = "<-" if language is R else "="
            head = f"{_ph(counter.next(), 'result')} {operator} {callee}("
        else:
            head = f"{callee}("

        pairs = [(arg.name, self._wrapper_value(arg, counter, language))
                 for arg in self._python._inputs(wrapper)]
        return _render_call(head, pairs, _STYLES[language])

    def _wrapper_value(self, argument: CArgument, counter: _Counter,
                       language: str) -> str:
        if argument.mode is not None:
            return _choice(counter.next(),
                           [v.string for v in argument.mode.values])
        if argument.optional:
            return _ph(counter.next(), "NULL" if language is R else "None")
        if argument.has_default:
            literal = r_literal if language is R else python_literal
            return _ph(counter.next(), literal(argument.default))
        return _ph(counter.next(), argument.name)

    # -- generic aids -----------------------------------------------------------

    def _generic_snippets(self, catalogue: ErrorCatalogue, buckets: dict) -> None:
        # all generic aids are tox-rooted, so they land in the tox files
        self._add(
            buckets, R, "tox", "tensor_omics (R setup)", "tox:setup",
            "Source the tensor-omics R wrappers",
            [f'source("{_ph(1, R_LOADER)}")'],
        )
        self._add(
            buckets, PYTHON, "tox", "tox_error (handler)", "tox:try",
            "Call a wrapper and handle a ToxError",
            [
                "try:",
                f"    {_ph(1, 'result')} = {_ph(2, 'call')}",
                f"except ToxError as {_ph(3, 'error')}:",
                f"    {_ph(0, 'raise')}",
            ],
        )
        self._add(
            buckets, R, "tox", "tox_error (handler)", "tox:trycatch",
            "Call a wrapper and handle a tox_error condition",
            [
                f"{_ph(1, 'result')} <- tryCatch(",
                f"  {_ph(2, 'call')},",
                f"  tox_error = function({_ph(3, 'error')}) {{",
                f"    {_ph(0, 'stop(error)')}",
                "  }",
                ")",
            ],
        )
        codes = [code.name for code in catalogue.codes]
        if codes:
            self._add(
                buckets, FORTRAN, "tox", "tox_errors code", "tox:err",
                "Pick a tox_errors code", [_choice(1, codes)],
            )


# -- helpers --------------------------------------------------------------------


class _Counter:
    """Hands out the ascending tabstop numbers a single snippet uses."""

    def __init__(self):
        self._n = 0

    def next(self) -> int:
        self._n += 1
        return self._n


def _ph(index: int, text: str) -> str:
    """A placeholder tabstop, `${1:text}`."""
    return f"${{{index}:{text}}}"


def _choice(index: int, values: list[str]) -> str:
    """A choice tabstop, `${1|a,b,c|}`. Our names never contain `,` or `|`."""
    return f"${{{index}|{','.join(values)}|}}"


def _render_call(head: str, pairs: list[tuple[str, str]], style: _Style) -> list[str]:
    """A keyword call spread one argument per line, in the given language's style."""
    if not pairs:
        return [f"{head})"]

    lines = [f"{head}{style.cont}"]
    last = len(pairs) - 1
    for index, (name, value) in enumerate(pairs):
        is_last = index == last
        comma = "" if is_last and not style.trailing_comma else ","
        if not is_last:
            tail = style.cont
        elif style.close_on_own_line:
            tail = ""
        else:
            tail = ")"
        lines.append(f"{style.indent}{name}{style.sep}{value}{comma}{tail}")
    if style.close_on_own_line:
        lines.append(")")
    return lines


def _split_module(name: str) -> tuple[str, str]:
    """A module name into its root (`f42`/`tox`) and the rest.

    `tox_clustering` -> `(tox, clustering)`; `tox_errors` -> `(tox, errors)`, which is why
    the error module needs no special-casing for its `tox:` prefix.
    """
    root, _, rest = name.partition("_")
    return root, rest or root
