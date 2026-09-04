"""Emitting the Python error module from the `tox_errors` catalogue.

Fortran reports failure by returning a code; Python reports it by raising. The generated
module is what turns one into the other, so that `ierr` never reaches a caller: a
successful call returns its results and nothing else, and a failed one raises.

The messages are `tox_errors`' own documentation, and the exception hierarchy follows the
numeric ranges it already groups its codes into. Nothing here is a second opinion about
what the errors are -- that is the point of generating it.
"""

from __future__ import annotations

from ..ir.errors import ErrorCatalogue, ErrorGroup
from ..render import Writer

#: The exception class for each group. `Tox` prefixed so `except ToxError` cannot catch
#: an unrelated RuntimeError, and one class per group rather than per code: 30-odd
#: exception types nobody would ever name individually is not a usable API.
GROUP_EXCEPTIONS = {
    ErrorGroup.IO: ("ToxIOError", "Reading or writing a file failed."),
    ErrorGroup.INPUT: ("ToxInputError", "An argument was not acceptable."),
    ErrorGroup.MEMORY: ("ToxMemoryError", "Allocation failed, or a pointer was null."),
    ErrorGroup.RUNTIME: ("ToxRuntimeError", "The Fortran runtime reported a problem."),
    ErrorGroup.INTERNAL: ("ToxInternalError", "A bug in tensor-omics."),
    ErrorGroup.OTHER: ("ToxError", "An error from tensor-omics."),
}

BASE_EXCEPTION = "ToxError"


class PythonErrorEmitter:
    def __init__(self, catalogue: ErrorCatalogue):
        self.catalogue = catalogue

    def module(self) -> str:
        writer = Writer()
        self._header(writer)
        self._exceptions(writer)
        self._codes(writer)
        self._tables(writer)
        self._check(writer)
        return writer.render(trailing_newline=True)

    def _header(self, writer: Writer) -> None:
        writer.block(
            '"""Errors raised by the tensor_omics binding.\n'
            "\n"
            "Generated from tox_errors. Do not edit.\n"
            '"""'
        )
        writer.blank()
        writer.line("#: How tox_errors packs the position of the offending argument into a code")
        writer.line(f"ARG_POS_FACTOR = {self.catalogue.arg_pos_factor}")
        writer.blank()

    def _exceptions(self, writer: Writer) -> None:
        writer.line(f"class {BASE_EXCEPTION}(RuntimeError):")
        with writer.indent():
            writer.block(
                '"""An error reported by tensor-omics.\n'
                "\n"
                "Attributes\n"
                "----------\n"
                "code : int\n"
                "    The tox_errors code, without the argument position packed into it.\n"
                "name : str or None\n"
                "    The name of the code, e.g. 'ERR_INVALID_INPUT'.\n"
                "argument : str or None\n"
                "    The argument the error was raised for, if it named one.\n"
                '"""'
            )
            writer.blank()
            writer.line("def __init__(self, message, code, name=None, argument=None):")
            with writer.indent():
                writer.line("super().__init__(message)")
                writer.line("self.code = code")
                writer.line("self.name = name")
                writer.line("self.argument = argument")
        writer.blank(collapse=False)
        writer.blank(collapse=False)

        for group in self.catalogue.groups():
            name, description = GROUP_EXCEPTIONS[group]
            if name == BASE_EXCEPTION:
                continue
            writer.line(f"class {name}({BASE_EXCEPTION}):")
            with writer.indent():
                writer.line(f'"""{description}"""')
            writer.blank(collapse=False)
            writer.blank(collapse=False)

    def _codes(self, writer: Writer) -> None:
        writer.line("# The codes themselves, so a caller can compare against a name")
        for code in self.catalogue:
            writer.line(f"{code.name} = {code.value}")
            if code.doc.summary:
                writer.line(f"#: {code.doc.summary}")
        writer.blank()

    def _tables(self, writer: Writer) -> None:
        errors = self.catalogue.errors
        writer.line("_MESSAGES = {")
        with writer.indent():
            for code in errors:
                writer.line(f"{code.value}: {_quote(code.message)},")
        writer.line("}")
        writer.blank()

        writer.line("_NAMES = {")
        with writer.indent():
            for code in errors:
                writer.line(f"{code.value}: {_quote(code.name)},")
        writer.line("}")
        writer.blank()

        writer.line("_EXCEPTIONS = {")
        with writer.indent():
            for code in errors:
                writer.line(f"{code.value}: {GROUP_EXCEPTIONS[code.group][0]},")
        writer.line("}")
        writer.blank()

        writer.line("#: Status codes are outcomes, not failures, and never raise")
        writer.line("_STATUSES = {")
        with writer.indent():
            for code in self.catalogue.statuses:
                writer.line(f"{code.value}: {_quote(code.name)},")
        writer.line("}")
        writer.blank()

    def _check(self, writer: Writer) -> None:
        writer.line("def check_err_code(ierr, arguments=(), sources=()):")
        with writer.indent():
            writer.block(
                '"""Raise if `ierr` reports an error.\n'
                "\n"
                "Parameters\n"
                "----------\n"
                "ierr : int\n"
                "    The code as tox_errors encoded it, argument position included.\n"
                "arguments : sequence of str, optional\n"
                "    The names of the wrapped procedure's arguments, in declaration\n"
                "    order, so the message can name the offending one rather than\n"
                "    give its number.\n"
                "sources : sequence of str or None, optional\n"
                "    Positionally alongside `arguments`: the argument the caller actually\n"
                "    passed, where the Fortran one was derived from it. An extent read off\n"
                "    an array names that array here, so the message can blame something the\n"
                "    caller wrote. None where the Fortran argument is the caller's own.\n"
                "\n"
                "Returns\n"
                "-------\n"
                "str or None\n"
                "    The name of the status code, if the call reported one. Status\n"
                "    codes are outcomes rather than failures and never raise.\n"
                "\n"
                "Raises\n"
                "------\n"
                "ToxError\n"
                "    If `ierr` reports an error. The subclass follows the kind of\n"
                "    error; every one of them is a ToxError.\n"
                '"""'
            )
            writer.line("code = ierr % ARG_POS_FACTOR")
            writer.line("if code == 0:")
            with writer.indent():
                writer.line("return None")
            writer.blank()
            writer.line("if code in _STATUSES:")
            with writer.indent():
                writer.line("return _STATUSES[code]")
            writer.blank()
            writer.line("arg_pos = ierr // ARG_POS_FACTOR")
            writer.line('message = _MESSAGES.get(code, f"unmapped error code {code}")')
            writer.line("argument = None")
            writer.line("if arg_pos > 0:")
            with writer.indent():
                writer.line("if arg_pos <= len(arguments):")
                with writer.indent():
                    writer.line("argument = arguments[arg_pos - 1]")
                    writer.line("source = sources[arg_pos - 1] if arg_pos <= len(sources) else None")
                    writer.line("if source:")
                    with writer.indent():
                        writer.line("# the caller's word first, then the Fortran argument it")
                        writer.line("# was derived from, which is what the signature calls it")
                        writer.line(
                            'message = f"{message} (argument \'{source}\', via \'{argument}\')"'
                        )
                        writer.line("argument = source")
                    writer.line("else:")
                    with writer.indent():
                        writer.line("message = f\"{message} (argument '{argument}')\"")
                writer.line("else:")
                with writer.indent():
                    writer.line('message = f"{message} (argument {arg_pos})"')
            writer.blank()
            writer.line(f"raise _EXCEPTIONS.get(code, {BASE_EXCEPTION})(")
            with writer.indent():
                writer.line("message, code=code, name=_NAMES.get(code), argument=argument")
            writer.line(")")
        return None


def _quote(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'
