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

#: What a generated binding raises for an argument it refuses itself -- a wrong shape, type or
#: value found before the library is called. It carries no tox code (`code` is None), because
#: none was reported: main's convention for a binding-side refusal (#221). An input error all
#: the same, so `except ToxInputError` catches it beside the Fortran's own; and a `TypeError`
#: and a `ValueError`, which is what these refusals were before, so an existing `except
#: ValueError` -- the numpy idiom -- keeps catching them.
ARGUMENT_EXCEPTION = "ToxArgumentError"

#: The code the C wrapper reports for an array that does not hold what its declaration computes
#: (#221). Its message is completed with the declared shape the binding passes along.
SIZE_MISMATCH = "ERR_SIZE_MISMATCH"

#: The constant the argument position is packed by, and the function that decodes a code, as
#: the module defines them and the package re-exports them
ARG_POS_FACTOR = "ARG_POS_FACTOR"
CHECK_ERR_CODE = "check_err_code"


def group_exceptions(catalogue: ErrorCatalogue) -> tuple[tuple[str, str], ...]:
    """The class of each group the catalogue has codes in, as (name, docstring), base aside."""
    return tuple(
        GROUP_EXCEPTIONS[group]
        for group in catalogue.groups()
        if GROUP_EXCEPTIONS[group][0] != BASE_EXCEPTION
    )


def public_names(catalogue: ErrorCatalogue) -> tuple[str, ...]:
    """Everything the error module defines for a caller, in the order it defines them.

    The one list of it: the module is emitted from its parts, and the package's `__init__`
    re-exports exactly this, so a code added to `tox_errors` -- an `ERR_*` or a `STAT_*` -- is
    importable from `tensor_omics` with no edit anywhere, and a caller never has to know which
    submodule holds it. The underscored tables are the module's own.
    """
    return (
        ARG_POS_FACTOR,
        BASE_EXCEPTION,
        *(name for name, _ in group_exceptions(catalogue)),
        ARGUMENT_EXCEPTION,
        *(code.name for code in catalogue),
        CHECK_ERR_CODE,
    )


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
        writer.line(f"{ARG_POS_FACTOR} = {self.catalogue.arg_pos_factor}")
        writer.blank()

    def _exceptions(self, writer: Writer) -> None:
        writer.line(f"class {BASE_EXCEPTION}(RuntimeError):")
        with writer.indent():
            writer.block(
                '"""An error reported by tensor-omics.\n'
                "\n"
                "Attributes\n"
                "----------\n"
                "code : int or None\n"
                "    The tox_errors code, without the argument position packed into it; None\n"
                "    for a ToxArgumentError, which the binding raised before any code existed.\n"
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

        for name, description in group_exceptions(self.catalogue):
            writer.line(f"class {name}({BASE_EXCEPTION}):")
            with writer.indent():
                writer.line(f'"""{description}"""')
            writer.blank(collapse=False)
            writer.blank(collapse=False)

        self._argument_exception(writer)

    def _argument_exception(self, writer: Writer) -> None:
        """`ToxArgumentError`, under the input group's class where the catalogue has one."""
        parent = (GROUP_EXCEPTIONS[ErrorGroup.INPUT][0]
                  if ErrorGroup.INPUT in self.catalogue.groups() else BASE_EXCEPTION)
        writer.line(f"class {ARGUMENT_EXCEPTION}({parent}, TypeError, ValueError):")
        with writer.indent():
            writer.block(
                '"""An argument the binding refused before calling the library.\n'
                "\n"
                "A wrong shape, type or value, found by the generated function itself, so it\n"
                "carries no tox code: `code` and `name` are None. It is a TypeError and a\n"
                "ValueError as well, so code catching either keeps working.\n"
                '"""'
            )
            writer.blank()
            writer.line("def __init__(self, message, argument=None):")
            with writer.indent():
                writer.line("super().__init__(message, code=None, argument=argument)")
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
        writer.line(f"def {CHECK_ERR_CODE}(ierr, arguments=(), sources=(), shapes=()):")
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
                "shapes : sequence of str or None, optional\n"
                "    Positionally alongside `arguments`: the declared shape of an array\n"
                "    whose size the C wrapper compares with its declaration, which the\n"
                f"    message of {SIZE_MISMATCH} then states. None elsewhere.\n"
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
                    if self.catalogue.code(SIZE_MISMATCH) is not None:
                        writer.line("shape = shapes[arg_pos - 1] if arg_pos <= len(shapes) else None")
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
                    size_mismatch = self.catalogue.code(SIZE_MISMATCH)
                    if size_mismatch is not None:
                        writer.line(f"if shape and code == {size_mismatch.value}:")
                        with writer.indent():
                            writer.line("# what the array must be, which no code can carry")
                            writer.line('message = f"{message}; \'{arguments[arg_pos - 1]}\' is '
                                        'declared with shape ({shape})"')
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
