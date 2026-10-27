"""Emitting the R half of the R binding.

The thinking layer. Per `design/language-layers.md`, this checks what every input is --
its kind and rank, and a single value in full -- coerces a double, cross-checks shapes,
applies defaults, calls the internal `.Call("<name>_call", ...)`, and turns `ierr` into a
raised condition. The one thing it leaves to the C shim is reading every element of an
integer, logical or character array, which the shim does in the pass that converts it
(#221): in R that check cost more than the cheapest call it guards. Both raise the same
classed conditions a caller can `tryCatch`, naming the argument in the caller's own words.
"""

from __future__ import annotations

from .error_arguments import shapes_of, sources_of
from ..abi.model import CArgument, Conversion, CWrapper, CWrapperModule, Origin
from ..ir.types import BaseType, Intent
from ..render import Writer
from .doc_literals import render as _render
from .doc_roxygen import render_roxygen
from .shape_checks import Position, declared_shape, extent_provider, reason, shape_checks


class RWrapperEmitter:
    def __init__(self, links=None, tiers=None):
        #: resolves a Ford link to what it is called here; None in a unit test, where a link
        #: then renders as plain code rather than a cross-reference
        self.links = links
        #: what each half of an expert/allocating pair says about the other, by published
        #: name (`emit.doc_tiers`); empty where a procedure reaches R only once
        self.tiers = tiers or {}

    def validators(self) -> str:
        """The R half of the input contract, in one file rather than repeated per module.

        Each checks an argument's kind (`is.numeric`, `is.logical`, `is.character`, never
        `typeof`) and rank, coerces a double, copying only when it must, and raises a
        classed `tox_type_error` naming the argument. A single value is checked in full
        here, NA included; an integer, logical or character array's elements are the C
        shim's to check, in the pass that converts them (`tox_marshal.h`); a double's are
        never scanned (an NA is a NaN payload Fortran itself catches). See
        `design/language-layers.md`.
        """
        return _VALIDATORS

    def module(self, module: CWrapperModule) -> str:
        # a module constant in an extent is its value here: this language cannot name it
        module = module.with_constant_sizes_evaluated()
        writer = Writer()
        writer.line("# Generated. Do not edit.")
        writer.blank()
        for wrapper in module:
            writer.block(self.function(wrapper))
            writer.blank(collapse=False)
        return writer.render(trailing_newline=True)

    # -- function ---------------------------------------------------------------

    def function(self, wrapper: CWrapper) -> str:
        writer = Writer()
        writer.block(render_roxygen(wrapper, self))

        params = ", ".join(self._param(a) for a in self._inputs(wrapper))
        writer.line(f"{wrapper.stripped_name} <- function({params}) {{")
        with writer.indent():
            self._validate(writer, wrapper)
            self._compute_auto(writer, wrapper)
            self._check_shapes(writer, wrapper)
            self._call(writer, wrapper)
            self._return(writer, wrapper)
        writer.line("}")
        return writer.render()

    def _compute_auto(self, writer: Writer, wrapper: CWrapper) -> None:
        """Fill each DM_OUTPUT_FROM(AUTO) argument by calling the producer's R wrapper.

        As in Python: the producer's own wrapper validates and error-checks, so this only
        calls it, with the consumer's already-validated inputs. One call serves every
        argument it supplies -- the producer behind `read_tox_data_into` unpacks a zip
        archive, so calling it once per extent would re-do that work each time.
        """
        from ..abi.c_abi import stripped_name

        wrote = False
        calls: dict[tuple, str] = {}
        for argument in wrapper:
            roles = argument.source.roles if argument.source else None
            plan = roles.computed_from if roles else None
            if plan is None or not plan.is_automatic:
                continue
            producer = stripped_name(plan.producer)
            call_args = ", ".join(
                f"{supply.name} = {self._producer_value(supply, wrapper)}"
                for supply in plan.inputs
            )
            if self._producer_returns_single(plan.producer):
                writer.line(f"{argument.name} <- {producer}({call_args})")
                wrote = True
                continue

            key = (producer, call_args)
            if key not in calls:
                suffix = "" if not calls else f"_{len(calls)}"
                calls[key] = f".{producer}_result{suffix}"
                writer.line(f"{calls[key]} <- {producer}({call_args})")
                for supply in plan.inout_feedback:
                    # the producer refined a value we handed it (e.g. capped it); adopt it
                    writer.line(f"{supply.argument} <- {calls[key]}${supply.name}")
            writer.line(f"{argument.name} <- {calls[key]}${plan.output.name}")
            wrote = True
        if wrote:
            writer.blank()

    def _producer_value(self, supply, wrapper: CWrapper) -> str:
        """What to pass for one producer input: an argument, or a constant outright."""
        if supply.is_constant:
            return r_literal(supply.constant)
        return self._r_value_of(supply.argument, wrapper)

    def _r_value_of(self, consumer_arg: str, wrapper: CWrapper) -> str:
        """The R expression for a consumer argument, to feed an AUTO producer call.

        A user input is just its name; a derived extent is not in scope as a variable, so
        it is spelled out (`length(values)`), which is what the caller would compute anyway.
        """
        argument = wrapper.argument(consumer_arg)
        roles = argument.source.roles if argument and argument.source else None
        if roles is not None and roles.is_extent:
            # the input it is read off, as everywhere (`shape_checks.extent_provider`); the
            # other inputs carrying it are what `_check_shapes` verifies against it
            position = extent_provider(argument, wrapper)
            if position is not None:
                return self._position(position)
        if roles is not None and roles.is_mask_count:
            return _mask_count(roles.mask_count_of.name)
        return consumer_arg

    @staticmethod
    def _producer_returns_single(producer) -> bool:
        outputs = [
            a for a in producer.arguments
            if a.intent and a.intent.is_output
            and a.name.lower() != "ierr"
            and not (a.roles and a.roles.is_temporary)
        ]
        if producer.result is not None:
            outputs.append(producer.result)
        return len(outputs) == 1
        return writer.render()

    def _inputs(self, wrapper: CWrapper) -> list[CArgument]:
        """What the user supplies -- the R function's parameters.

        Computed (AUTO) arguments are excluded: R works them out itself. They are not,
        however, excluded from the `.Call` shim, which C needs them for.
        """
        return [
            argument
            for argument in wrapper
            if argument.intent.is_input
            and not argument.is_temporary
            and not self._is_derived(argument)
            and (not argument.is_synthesised or self._must_be_supplied(argument, wrapper))
        ]

    def _call_inputs(self, wrapper: CWrapper) -> list[CArgument]:
        """What the `.Call` shim receives: the user inputs, plus computed args R fills."""
        return [
            argument
            for argument in wrapper
            if argument.intent.is_input
            and not argument.is_temporary
            and not self._is_derived_in_c(argument)
            and (not argument.is_synthesised or self._must_be_supplied(argument, wrapper))
        ]

    @staticmethod
    def _must_be_supplied(argument: CArgument, wrapper: CWrapper) -> bool:
        """Whether a synthesised extent or strlen has to come from the caller after all.

        Kept in step with the C emitter, which decides the same thing for the `.Call`
        signature: an extent or strlen that sizes only an *output* has nothing to be read
        off, so it stays a parameter on both sides.
        """
        from .c_call import CCallEmitter

        return CCallEmitter()._must_be_supplied(argument, wrapper)

    @staticmethod
    def _is_derived(argument: CArgument) -> bool:
        roles = argument.source.roles if argument.source else None
        return bool(roles and roles.is_derived)

    @staticmethod
    def _is_derived_in_c(argument: CArgument) -> bool:
        # C derives extents/shapes/counts, but not computed args (R passes those)
        roles = argument.source.roles if argument.source else None
        if roles is None:
            return False
        if roles.is_computed:
            return False
        return (roles.is_inferable_extent or roles.is_inferable_shape_arg
                or roles.is_mask_count)

    def _param(self, argument: CArgument) -> str:
        if argument.optional:
            return f"{argument.name} = NULL"
        if argument.has_default:
            return f"{argument.name} = {r_literal(argument.default)}"
        return argument.name

    # -- validation -------------------------------------------------------------

    def _validate(self, writer: Writer, wrapper: CWrapper) -> None:
        for argument in self._inputs(wrapper):
            coercion = self._coercion(argument)
            line = f'{argument.name} <- {coercion}'
            if argument.optional:
                # only coerce when given; NULL means absent
                writer.line(f"if (!is.null({argument.name}))")
                with writer.indent():
                    writer.line(line)
            else:
                writer.line(line)

    def _coercion(self, argument: CArgument) -> str:
        name = f'"{argument.name}"'
        value = argument.name

        if argument.shape_arg is not None:
            # its shape is read off it after this, so the coercion must not drop it
            base = {
                BaseType.INTEGER: "integer",
                BaseType.REAL: "double",
                BaseType.COMPLEX: "complex",
                BaseType.LOGICAL: "logical",
                BaseType.CHARACTER: "character",
            }[argument.type.base]
            return f".tox_as_{base}_shaped({value}, {name})"

        if argument.conversion is Conversion.MODE:
            choices = ", ".join(f'"{v.string}"' for v in argument.mode.values)
            return f".tox_as_mode({value}, {name}, c({choices}))"

        if argument.type.is_character:
            # a string's length is its leading C extent, so a single string has rank 1 here
            helper = ".tox_as_string" if argument.rank <= 1 else ".tox_as_character"
            return f"{helper}({value}, {name})"

        roles = argument.source.roles if argument.source else None
        if roles is not None and roles.is_shape_arg and argument.type.base is BaseType.INTEGER:
            # the shim uses this vector itself, as the dim of what it allocates, rather than
            # converting it for Fortran; a shape is a handful of elements, checked here
            return f".tox_as_integer_shape({value}, {name})"

        base = {
            BaseType.INTEGER: "integer",
            BaseType.REAL: "double",
            BaseType.COMPLEX: "complex",
            BaseType.LOGICAL: "logical",
        }[argument.type.base]

        if argument.is_scalar:
            return f".tox_as_{base}_scalar({value}, {name})"
        if argument.rank == 1:
            return f".tox_as_{base}_vector({value}, {name})"
        if argument.rank == 2:
            return f".tox_as_{base}_matrix({value}, {name})"
        return f".tox_as_{base}_array({value}, {name}, {argument.rank}L)"

    def _check_shapes(self, writer: Writer, wrapper: CWrapper) -> None:
        """Compare every input axis with what it must be, as Python does (`emit.shape_checks`).

        Raised as a `tox_shape_error` classed condition rather than a `ToxArgumentError`. R has
        no extent of its own in hand where the shim derives it, so the value read off the
        providing array is written out here as that array's extent.
        """
        wrote = False
        for check in shape_checks(wrapper):
            array = check.position.array
            actual = self._position(check.position)
            expected = self._expected(check.expected)
            condition = f"{actual} != {expected}"
            if array.optional:
                condition = f"!is.null({array.name}) && {condition}"
            dimension = "NA_integer_" if check.position.axis is None else f"{check.position.axis + 1}L"
            shape = _render(declared_shape(array), "r")
            writer.line(f"if ({condition})")
            with writer.indent():
                writer.line(
                    f'.tox_shape_error("{array.name}", {dimension}, {actual}, '
                    f'"{reason(check.expected)}", {expected}, "{shape}")'
                )
            wrote = True
        if wrote:
            writer.blank()

    def _expected(self, expected) -> str:
        """The value a compared axis must have, as R computes it before the call."""
        if expected.origin == "constant":
            return f"{expected.value}L"
        if expected.origin == "axis":
            return self._position(expected.position)
        if expected.origin == "mask":
            return _mask_count(expected.mask)
        # an argument the caller passed, or one the producer call above computed
        return expected.extent

    def _position(self, position: Position) -> str:
        """An axis of an input as R reads it, or its element count."""
        owner = position.array
        if position.axis is None:
            return f"length({owner.name})"
        return self._extent_expression(owner, position.axis)

    def _extent_expression(self, owner: CArgument, axis: int) -> str:
        if owner.rank <= 1 or (owner.type.is_character and owner.rank == 2):
            return f"length({owner.name})"
        return f"dim({owner.name})[{axis + 1}]"

    # -- call and return --------------------------------------------------------

    def _call(self, writer: Writer, wrapper: CWrapper) -> None:
        inputs = ", ".join(a.name for a in self._call_inputs(wrapper))
        prefix = f'"{wrapper.stripped_name}_call"'
        args = f"{prefix}, {inputs}" if inputs else prefix
        writer.line(f".result <- .Call({args})")
        names = ", ".join(f'"{a.name}"' for a in wrapper.procedure.arguments)
        writer.line(f".arguments <- c({names})")
        # and, positionally, the argument the caller passed that each was derived from
        sources = sources_of(wrapper.procedure)
        if any(sources):
            rendered = ", ".join(f'"{s}"' if s else "NA_character_" for s in sources)
            writer.line(f".sources <- c({rendered})")
        # and the declared shape of an array the C wrapper compares, which its error states
        shapes = shapes_of(wrapper, "r")
        if any(shapes):
            rendered = ", ".join(f'"{s}"' if s else "NA_character_" for s in shapes)
            writer.line(f".shapes <- c({rendered})")
        tables = [".arguments"]
        if any(sources):
            tables.append(".sources")
        if any(shapes):
            tables.append("shapes = .shapes")
        writer.line(f".status <- check_err_code(.result$ierr, {', '.join(tables)})")
        writer.blank()

    def _return(self, writer: Writer, wrapper: CWrapper) -> None:
        outputs = self._outputs(wrapper)
        if not outputs:
            writer.line("invisible(NULL)")
            return
        if len(outputs) == 1:
            writer.line(self._result_expression(outputs[0], wrapper))
            return
        writer.line("list(")
        with writer.indent():
            for index, argument in enumerate(outputs):
                comma = "" if index == len(outputs) - 1 else ","
                writer.line(f'{argument.name} = {self._result_expression(argument, wrapper)}{comma}')
        writer.line(")")

    def _outputs(self, wrapper: CWrapper) -> list[CArgument]:
        """What the R caller gets back: outputs and modified inputs, ierr and scaffolding aside."""
        consumed = {
            a.source.roles.result_size_arg.name
            for a in wrapper
            if a.source and a.source.roles and a.source.roles.result_size_arg
        }
        return [
            argument
            for argument in wrapper
            if argument.intent.is_output
            and not argument.is_error
            and not argument.is_temporary
            and argument.name not in consumed
        ]

    def _result_expression(self, argument: CArgument, wrapper: CWrapper) -> str:
        source = argument.source
        roles = source.roles if source else None
        value = f".result${argument.name}"

        if roles is not None and roles.result_size_arg is not None:
            count = f".result${roles.result_size_arg.name}"
            if argument.rank <= 1:
                return f"utils::head({value}, {count})"
            leading = ", " * (argument.rank - 1)
            return f"{value}[{leading}seq_len({count}), drop = FALSE]"

        if argument.type.is_character and argument.is_scalar:
            # a scalar string comes back as a length-1 character vector
            return f"{value}[[1]]"
        return value


def _mask_count(mask: str) -> str:
    """How many elements of a mask are set, as R counts them before the call.

    A mask reaches the shim unchecked -- its elements are the shim's to check, as it converts
    them (`tox_validate.R`) -- so it may still hold an NA or a number that is neither 0 nor 1
    here. Counted as the shim counts it, the exact TRUEs and 1s, it is always a whole number,
    which `.tox_shape_error` prints as one, and the element is then refused by name in C.
    """
    return f"sum({mask} == 1, na.rm = TRUE)"


def r_literal(value) -> str:
    """A Python value as R source, for a default."""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, int):
        return f"{value}L"
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, str):
        return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'
    return repr(value)


_VALIDATORS = r'''# Generated. Do not edit.
#
# The R half of the input contract (helper/codegen/design/language-layers.md). Each helper
# checks what an argument IS -- its kind, by is.numeric / is.logical / is.character and never
# by typeof (a factor's typeof is integer, and its codes are not its values), and its rank --
# and raises a classed tox_type_error naming it. A single value is checked here in full, NA
# included. The elements of an integer, logical or character array are not: the C shim checks
# every one while it converts it for Fortran, in the pass it makes anyway (tox_marshal.h), so
# an array reaches it unconverted. Doubles are never scanned: an NA is a NaN payload, which
# Fortran's own validation catches.

# The rank R sees: a plain vector, a 1-D array and a table are all rank 1
.tox_rank <- function(x) {
  dims <- dim(x)
  if (is.null(dims)) 1L else length(dims)
}

# A number is what is.numeric() calls one -- not a factor, a Date, a POSIXct or a difftime --
# except bit64's integer64: a double vector holding int64 bits, which is.numeric() calls
# numeric and every value of which would be read as a different number
.tox_check_number <- function(x, name, expected, logical = FALSE) {
  if (inherits(x, "integer64"))
    .tox_argument_error(
      sprintf("'%s' is a bit64 integer64, whose doubles hold int64 bits; convert it with as.integer() or as.double()", name),
      "tox_type_error", name)
  if (!is.numeric(x) && !(logical && is.logical(x))) .tox_type_error(name, expected, x)
}

# A data.frame whose every column is a number (or a logical, for a logical matrix) is a
# matrix, its rows staying rows, as a pandas DataFrame is in Python
.tox_frame_as_matrix <- function(x, name, expected, logical = FALSE) {
  if (!is.data.frame(x)) return(x)
  is_column <- function(column)
    !inherits(column, "integer64") && (is.numeric(column) || (logical && is.logical(column)))
  if (!all(vapply(x, is_column, NA))) .tox_type_error(name, expected, x)
  as.matrix(x)
}

# A single number that must be a whole number in R's integer range, which leaves out
# -2^31: that is NA_integer_. A NaN is refused for its value, an NA for being NA.
.tox_check_whole <- function(x, name) {
  if (is.na(x) && !is.nan(x)) .tox_na_error(name)
  if (is.double(x) && !isTRUE(x >= -.Machine$integer.max && x <= .Machine$integer.max && x == trunc(x)))
    .tox_argument_error(
      sprintf("'%s' must be a whole number in the integer range, not %s", name, format(x, digits = 15)),
      "tox_type_error", name)
}

.tox_as_double_vector <- function(x, name) {
  .tox_check_number(x, name, "a numeric vector")
  if (.tox_rank(x) != 1L) .tox_type_error(name, "a numeric vector", x)
  as.double(x)
}

.tox_as_double_scalar <- function(x, name) {
  .tox_check_number(x, name, "a numeric scalar")
  if (length(x) != 1L || .tox_rank(x) != 1L) .tox_type_error(name, "a numeric scalar", x)
  as.double(x)
}

.tox_as_double_matrix <- function(x, name) {
  x <- .tox_frame_as_matrix(x, name, "a numeric matrix")
  if (!is.matrix(x)) .tox_type_error(name, "a numeric matrix", x)
  .tox_check_number(x, name, "a numeric matrix")
  storage.mode(x) <- "double"
  x
}

.tox_as_double_array <- function(x, name, ndim) {
  expected <- sprintf("a numeric array of rank %d", ndim)
  if (!is.array(x) || length(dim(x)) != ndim) .tox_type_error(name, expected, x)
  .tox_check_number(x, name, expected)
  storage.mode(x) <- "double"
  x
}

# An array whose shape travels in a separate argument keeps that shape through
# coercion: as.double() and friends drop the dim attribute, and the shape is read off
# the array after this point. storage.mode() converts in place instead.
.tox_as_double_shaped <- function(x, name) {
  .tox_check_number(x, name, "a numeric array")
  storage.mode(x) <- "double"
  x
}

.tox_as_integer_shaped <- function(x, name) {
  .tox_check_number(x, name, "an integer array")
  x
}

.tox_as_complex_shaped <- function(x, name) {
  if (!is.numeric(x) && !is.complex(x)) .tox_type_error(name, "a complex array", x)
  storage.mode(x) <- "complex"
  x
}

.tox_as_logical_shaped <- function(x, name) {
  .tox_check_number(x, name, "a logical array", logical = TRUE)
  x
}

.tox_as_character_shaped <- function(x, name) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) .tox_type_error(name, "a character array", x)
  x
}

.tox_as_complex_vector <- function(x, name) {
  if (!is.numeric(x) && !is.complex(x)) .tox_type_error(name, "a complex vector", x)
  as.complex(x)
}

.tox_as_complex_scalar <- function(x, name) {
  if ((!is.numeric(x) && !is.complex(x)) || length(x) != 1L)
    .tox_type_error(name, "a complex scalar", x)
  as.complex(x)
}

.tox_as_complex_matrix <- function(x, name) {
  if (!is.matrix(x) || (!is.numeric(x) && !is.complex(x)))
    .tox_type_error(name, "a complex matrix", x)
  storage.mode(x) <- "complex"
  x
}

.tox_as_complex_array <- function(x, name, ndim) {
  if (!is.array(x) || (!is.numeric(x) && !is.complex(x)) || length(dim(x)) != ndim)
    .tox_type_error(name, sprintf("a complex array of %d dimensions", ndim), x)
  storage.mode(x) <- "complex"
  x
}

.tox_as_integer_vector <- function(x, name) {
  .tox_check_number(x, name, "an integer vector")
  if (.tox_rank(x) != 1L) .tox_type_error(name, "an integer vector", x)
  x
}

.tox_as_integer_scalar <- function(x, name) {
  .tox_check_number(x, name, "an integer scalar")
  if (length(x) != 1L || .tox_rank(x) != 1L) .tox_type_error(name, "an integer scalar", x)
  .tox_check_whole(x, name)
  as.integer(x)
}

.tox_as_integer_matrix <- function(x, name) {
  x <- .tox_frame_as_matrix(x, name, "an integer matrix")
  if (!is.matrix(x)) .tox_type_error(name, "an integer matrix", x)
  .tox_check_number(x, name, "an integer matrix")
  x
}

.tox_as_integer_array <- function(x, name, ndim) {
  expected <- sprintf("an integer array of rank %d", ndim)
  if (!is.array(x) || length(dim(x)) != ndim) .tox_type_error(name, expected, x)
  .tox_check_number(x, name, expected)
  x
}

# The shape of an array the shim allocates: it uses the vector itself, as the dim of the
# result, so it is checked and made an integer vector here. It is a handful of elements.
.tox_as_integer_shape <- function(x, name) {
  .tox_check_number(x, name, "an integer vector")
  if (.tox_rank(x) != 1L) .tox_type_error(name, "an integer vector", x)
  for (value in x) .tox_check_whole(value, name)
  as.integer(x)
}

.tox_as_logical_scalar <- function(x, name) {
  .tox_check_number(x, name, "a logical scalar", logical = TRUE)
  if (length(x) != 1L || .tox_rank(x) != 1L) .tox_type_error(name, "a logical scalar", x)
  if (is.na(x) && !is.nan(x)) .tox_na_error(name)
  if (!is.logical(x) && !isTRUE(x == 0 || x == 1))
    .tox_argument_error(
      sprintf("'%s' must be TRUE or FALSE, or the number 0 or 1, not %s", name, format(x, digits = 15)),
      "tox_type_error", name)
  as.logical(x)
}

.tox_as_logical_vector <- function(x, name) {
  .tox_check_number(x, name, "a logical vector", logical = TRUE)
  if (.tox_rank(x) != 1L) .tox_type_error(name, "a logical vector", x)
  x
}

.tox_as_logical_matrix <- function(x, name) {
  x <- .tox_frame_as_matrix(x, name, "a logical matrix", logical = TRUE)
  if (!is.matrix(x)) .tox_type_error(name, "a logical matrix", x)
  .tox_check_number(x, name, "a logical matrix", logical = TRUE)
  x
}

.tox_as_logical_array <- function(x, name, ndim) {
  expected <- sprintf("a logical array of rank %d", ndim)
  if (!is.array(x) || length(dim(x)) != ndim) .tox_type_error(name, expected, x)
  .tox_check_number(x, name, expected, logical = TRUE)
  x
}

# A factor is taken by its labels, which are its values; as.integer() would give its codes
.tox_as_character <- function(x, name) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x) || .tox_rank(x) != 1L) .tox_type_error(name, "a character vector", x)
  x
}

.tox_as_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || .tox_rank(x) != 1L)
    .tox_type_error(name, "a single string", x)
  x
}

# Trailing blanks are dropped first, as for every string: to Fortran "ward " is "ward"
.tox_as_mode <- function(x, name, choices) {
  if (!is.character(x) || length(x) != 1L) .tox_type_error(name, "a single string", x)
  x <- tolower(sub(" +$", "", x))
  if (!x %in% choices)
    .tox_argument_error(
      sprintf("'%s' must be one of %s, not \"%s\"",
              name, paste(sprintf('"%s"', choices), collapse = ", "), x),
      "tox_type_error", name)
  x
}
'''
