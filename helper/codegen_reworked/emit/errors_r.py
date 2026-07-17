"""Emitting the R error module from the `tox_errors` catalogue.

The R mirror of `errors_python`: Fortran returns a code, R raises. `check_err_code` decodes
`ierr` and signals a **classed condition** so `tryCatch(tox_input_error = ...)` can select
on it. The message is the Fortran documentation, and the argument the position points at is
named.

Conditions rather than a bare `stop`, because a classed condition is what makes R errors
catchable by kind -- the R analogue of Python's exception subclasses.
"""

from __future__ import annotations

from ..ir.errors import ErrorCatalogue, ErrorGroup
from ..render import Writer

#: The condition class for each group. `tox_error` is the parent every one inherits, so
#: `tryCatch(tox_error = ...)` catches them all; the specific class allows finer nets.
GROUP_CONDITIONS = {
    ErrorGroup.IO: "tox_io_error",
    ErrorGroup.INPUT: "tox_input_error",
    ErrorGroup.MEMORY: "tox_memory_error",
    ErrorGroup.RUNTIME: "tox_runtime_error",
    ErrorGroup.INTERNAL: "tox_internal_error",
    ErrorGroup.OTHER: "tox_error",
}

BASE_CONDITION = "tox_error"


class RErrorEmitter:
    def __init__(self, catalogue: ErrorCatalogue):
        self.catalogue = catalogue

    def module(self) -> str:
        writer = Writer()
        writer.line("# Generated from tox_errors. Do not edit.")
        writer.blank()
        writer.line(f"ARG_POS_FACTOR <- {self.catalogue.arg_pos_factor}L")
        writer.blank()
        self._codes(writer)
        self._tables(writer)
        self._constructors(writer)
        self._check(writer)
        return writer.render(trailing_newline=True)

    def _codes(self, writer: Writer) -> None:
        writer.line("# the codes, exported so a caller can compare against a name")
        for code in self.catalogue:
            writer.line(f"{code.name} <- {code.value}L")
        writer.blank()

    def _tables(self, writer: Writer) -> None:
        errors = self.catalogue.errors

        writer.line(".tox_messages <- c(")
        with writer.indent():
            for index, code in enumerate(errors):
                comma = "" if index == len(errors) - 1 else ","
                writer.line(f'"{code.value}" = {_quote(code.message)}{comma}')
        writer.line(")")
        writer.blank()

        writer.line(".tox_names <- c(")
        with writer.indent():
            for index, code in enumerate(errors):
                comma = "" if index == len(errors) - 1 else ","
                writer.line(f'"{code.value}" = {_quote(code.name)}{comma}')
        writer.line(")")
        writer.blank()

        writer.line(".tox_classes <- c(")
        with writer.indent():
            for index, code in enumerate(errors):
                comma = "" if index == len(errors) - 1 else ","
                writer.line(f'"{code.value}" = "{GROUP_CONDITIONS[code.group]}"{comma}')
        writer.line(")")
        writer.blank()

        writer.line("# status codes are outcomes, not failures, and never raise")
        statuses = self.catalogue.statuses
        writer.line(".tox_statuses <- c(")
        with writer.indent():
            for index, code in enumerate(statuses):
                comma = "" if index == len(statuses) - 1 else ","
                writer.line(f'"{code.value}" = {_quote(code.name)}{comma}')
        writer.line(")")
        writer.blank()

    def _constructors(self, writer: Writer) -> None:
        writer.block(
            "#' Signal a tensor-omics error\n"
            "#'\n"
            "#' @param message the human-readable message\n"
            "#' @param class the specific condition class, a subclass of `tox_error`\n"
            "#' @param code the tox_errors code, position stripped\n"
            "#' @param name the code's name, e.g. `\"ERR_INVALID_INPUT\"`\n"
            "#' @param argument the offending argument, if the code named one\n"
            "#' @keywords internal"
        )
        writer.line(".tox_raise <- function(message, class, code = NA_integer_,")
        with writer.indent(5):
            writer.line("name = NA_character_, argument = NA_character_) {")
        with writer.indent():
            writer.line("cond <- structure(")
            with writer.indent():
                writer.line("class = c(class, \"tox_error\", \"error\", \"condition\"),")
                writer.line("list(message = message, call = sys.call(-1),")
            with writer.indent(5):
                writer.line("code = code, name = name, argument = argument)")
            writer.line(")")
            writer.line("stop(cond)")
        writer.line("}")
        writer.blank()

        # the type/shape/na conditions the R wrappers raise before the call
        writer.line("#' @keywords internal")
        writer.line(".tox_type_error <- function(argument, expected, value) {")
        with writer.indent():
            writer.line(".tox_raise(")
            with writer.indent():
                writer.line(
                    'sprintf("\'%s\' must be %s, not %s", argument, expected, '
                    'class(value)[1]),'
                )
                writer.line('"tox_type_error", argument = argument')
            writer.line(")")
        writer.line("}")
        writer.blank()

        writer.line("#' @keywords internal")
        writer.line(".tox_na_error <- function(argument) {")
        with writer.indent():
            writer.line(".tox_raise(")
            with writer.indent():
                writer.line(
                    'sprintf("\'%s\' contains NA, which tensor-omics cannot represent", '
                    "argument),"
                )
                writer.line('"tox_na_error", argument = argument')
            writer.line(")")
        writer.line("}")
        writer.blank()

        writer.line("#' @keywords internal")
        writer.line(".tox_shape_error <- function(argument, actual, other, expected) {")
        with writer.indent():
            writer.line(".tox_raise(")
            with writer.indent():
                writer.line(
                    'sprintf("\'%s\' has extent %d, but \'%s\' implies it should be %d",'
                )
                writer.line("argument, actual, other, expected),")
                writer.line('"tox_shape_error", argument = argument')
            writer.line(")")
        writer.line("}")
        writer.blank()

    def _check(self, writer: Writer) -> None:
        writer.block(
            "#' Raise if a tox_errors code reports an error\n"
            "#'\n"
            "#' @param ierr the encoded error code, argument position included\n"
            "#' @param arguments the wrapped procedure's argument names, in order, so the\n"
            "#'   message can name the offending one\n"
            "#' @return invisibly, the name of the status code if the call reported one;\n"
            "#'   status codes are outcomes, not failures, and never raise\n"
            "#' @keywords internal"
        )
        writer.line("check_err_code <- function(ierr, arguments = character()) {")
        with writer.indent():
            writer.line("code <- ierr %% ARG_POS_FACTOR")
            writer.line("if (code == 0L) return(invisible(NULL))")
            writer.blank()
            writer.line("key <- as.character(code)")
            writer.line("if (key %in% names(.tox_statuses))")
            with writer.indent():
                writer.line("return(invisible(.tox_statuses[[key]]))")
            writer.blank()
            writer.line("arg_pos <- ierr %/% ARG_POS_FACTOR")
            writer.line('message <- if (key %in% names(.tox_messages)) .tox_messages[[key]]')
            with writer.indent(11):
                writer.line('else sprintf("unmapped error code %d", code)')
            writer.line("argument <- NA_character_")
            writer.line("if (arg_pos > 0L) {")
            with writer.indent():
                writer.line("if (arg_pos <= length(arguments)) {")
                with writer.indent():
                    writer.line("argument <- arguments[[arg_pos]]")
                    writer.line('message <- sprintf("%s (argument \'%s\')", message, argument)')
                writer.line("} else {")
                with writer.indent():
                    writer.line('message <- sprintf("%s (argument %d)", message, arg_pos)')
                writer.line("}")
            writer.line("}")
            writer.blank()
            writer.line("class <- if (key %in% names(.tox_classes)) .tox_classes[[key]]")
            with writer.indent(9):
                writer.line('else "tox_error"')
            writer.line(".tox_raise(message, class, code = code, name = .tox_names[key],")
            with writer.indent(11):
                writer.line("argument = argument)")
        writer.line("}")


def _quote(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'
