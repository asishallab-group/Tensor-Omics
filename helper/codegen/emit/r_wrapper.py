"""Emitting the R half of the R interface.

The thinking layer. Per `design/language-layers.md`, **R decides and raises, C++ marshals
and calls**: this validates and coerces every input, cross-checks shapes, applies
defaults, calls the internal `.<name>_rcpp`, and turns `ierr` into a raised condition. All
of it in one place, all raised as classed conditions a caller can `tryCatch`.

It is a thin R function over the C++ one, but the thinness is the point -- everything a
user can get wrong is caught here, in their language, with their argument names.
"""

from __future__ import annotations

from ..abi.model import CArgument, Conversion, CWrapper, CWrapperModule, Origin
from ..ir.types import BaseType, Intent
from ..render import Writer
from .doc_roxygen import render_roxygen


class RWrapperEmitter:
    def validators(self) -> str:
        """The coercion helpers, in one file rather than repeated per module.

        Each checks the R type and coerces, copying only when it must, and raises a
        classed `tox_type_error` naming the argument. NA is checked exactly where it is
        free: on integers (whose sentinel is an ordinary Fortran number), on logicals and
        characters (converted anyway), and never on doubles (a NaN payload Fortran itself
        catches). See `design/language-layers.md`.
        """
        return _VALIDATORS

    def module(self, module: CWrapperModule) -> str:
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
            for owner in roles.extent_of:
                c_owner = wrapper.argument(owner.name)
                if c_owner is not None and c_owner.is_temporary:
                    # a work array the wrapper allocates: not in scope, and sized *from*
                    # this extent rather than the other way round
                    continue
                if c_owner is not None and c_owner.intent.is_input and not c_owner.optional:
                    axis = self._axis_of(consumer_arg, c_owner)
                    if axis is not None:
                        return self._extent_expression(c_owner, axis)
        if roles is not None and roles.is_mask_count:
            return f"sum({roles.mask_count_of.name})"
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
        however, excluded from the `.rcpp` call, which C++ needs them for.
        """
        return [
            argument
            for argument in wrapper
            if argument.intent.is_input
            and not argument.is_temporary
            and not self._is_derived(argument)
            and (not argument.is_synthesised or self._must_be_supplied(argument, wrapper))
        ]

    def _rcpp_inputs(self, wrapper: CWrapper) -> list[CArgument]:
        """What the .rcpp function receives: the user inputs, plus computed args R fills."""
        return [
            argument
            for argument in wrapper
            if argument.intent.is_input
            and not argument.is_temporary
            and not self._is_derived_in_cpp(argument)
            and (not argument.is_synthesised or self._must_be_supplied(argument, wrapper))
        ]

    @staticmethod
    def _must_be_supplied(argument: CArgument, wrapper: CWrapper) -> bool:
        """Whether a synthesised extent or strlen has to come from the caller after all.

        Kept in step with the Rcpp emitter, which decides the same thing for the C++
        signature: an extent or strlen that sizes only an *output* has nothing to be read
        off, so it stays a parameter on both sides.
        """
        from .rcpp import RcppEmitter

        return RcppEmitter()._must_be_supplied(argument, wrapper)

    @staticmethod
    def _is_derived(argument: CArgument) -> bool:
        roles = argument.source.roles if argument.source else None
        return bool(roles and roles.is_derived)

    @staticmethod
    def _is_derived_in_cpp(argument: CArgument) -> bool:
        # C++ derives extents/shapes/counts, but not computed args (R passes those)
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

        if argument.conversion is Conversion.MODE:
            choices = ", ".join(f'"{v.string}"' for v in argument.mode.values)
            return f".tox_as_mode({value}, {name}, c({choices}))"

        if argument.type.is_character:
            helper = ".tox_as_string" if argument.rank <= 1 and argument.is_scalar \
                else ".tox_as_character"
            return f"{helper}({value}, {name})"

        base = {
            BaseType.INTEGER: "integer",
            BaseType.REAL: "double",
            BaseType.COMPLEX: "complex",
            BaseType.LOGICAL: "logical",
        }[argument.type.base]

        if argument.type.base is BaseType.LOGICAL:
            return f".tox_as_logical({value}, {name})"
        if argument.is_scalar:
            return f".tox_as_{base}_scalar({value}, {name})"
        if argument.rank == 1:
            return f".tox_as_{base}_vector({value}, {name})"
        if argument.rank == 2:
            return f".tox_as_{base}_matrix({value}, {name})"
        return f".tox_as_{base}_array({value}, {name}, {argument.rank}L)"

    def _check_shapes(self, writer: Writer, wrapper: CWrapper) -> None:
        """The cross-checks Fortran cannot make: shared extents must agree.

        The same idea as the Python emitter, raised as a `tox_shape_error` classed
        condition rather than a `ValueError`.
        """
        wrote = False
        for argument in wrapper:
            roles = argument.source.roles if argument.source else None
            if roles is None or not roles.is_extent:
                continue

            owners = []
            for owner in roles.extent_of:
                c_owner = wrapper.argument(owner.name)
                if c_owner is None or not c_owner.intent.is_input or c_owner.optional:
                    continue
                if c_owner.is_temporary:
                    # a work array C++ allocates itself always agrees, and is not a
                    # variable R has in scope
                    continue
                axis = self._axis_of(argument.name, c_owner)
                if axis is not None:
                    owners.append((c_owner, axis))
            if len(owners) < 2:
                continue

            first, first_axis = owners[0]
            reference = self._extent_expression(first, first_axis)
            for owner, axis in owners[1:]:
                actual = self._extent_expression(owner, axis)
                writer.line(f"if ({actual} != {reference})")
                with writer.indent():
                    writer.line(
                        f'.tox_shape_error("{owner.name}", {actual}, "{first.name}", '
                        f"{reference})"
                    )
                wrote = True
        if wrote:
            writer.blank()

    def _extent_expression(self, owner: CArgument, axis: int) -> str:
        if owner.rank <= 1 or (owner.type.is_character and owner.rank == 2):
            return f"length({owner.name})"
        return f"dim({owner.name})[{axis + 1}]"

    @staticmethod
    def _axis_of(extent: str, owner: CArgument) -> int | None:
        extents = list(owner.dimension.extents)
        if owner.type.is_character and extents:
            extents = extents[1:]
        try:
            return extents.index(extent)
        except ValueError:
            return None

    # -- call and return --------------------------------------------------------

    def _call(self, writer: Writer, wrapper: CWrapper) -> None:
        inputs = ", ".join(a.name for a in self._rcpp_inputs(wrapper))
        writer.line(f".result <- .{wrapper.stripped_name}_rcpp({inputs})")
        names = ", ".join(f'"{a.name}"' for a in wrapper.procedure.arguments)
        writer.line(f".arguments <- c({names})")
        writer.line(".status <- check_err_code(.result$ierr, .arguments)")
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
            return f"utils::head({value}, .result${roles.result_size_arg.name})"

        if argument.type.is_character and argument.is_scalar:
            # a scalar string comes back as a length-1 character vector
            return f"{value}[[1]]"
        return value


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
# Coercion helpers for the R interface. Each checks the type and coerces, copying only
# when it must, and raises a classed tox_type_error naming the argument. NA is checked
# only where the check is free: integers (an ordinary Fortran number), logicals and
# characters (converted anyway), never doubles (a NaN payload Fortran already catches).

.tox_as_double_vector <- function(x, name) {
  if (!is.numeric(x)) .tox_type_error(name, "a numeric vector", x)
  as.double(x)
}

.tox_as_double_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L) .tox_type_error(name, "a numeric scalar", x)
  as.double(x)
}

.tox_as_double_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x)) .tox_type_error(name, "a numeric matrix", x)
  storage.mode(x) <- "double"
  x
}

.tox_as_double_array <- function(x, name, ndim) {
  if (!is.array(x) || !is.numeric(x) || length(dim(x)) != ndim)
    .tox_type_error(name, sprintf("a numeric array of rank %d", ndim), x)
  storage.mode(x) <- "double"
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
  if (!is.numeric(x)) .tox_type_error(name, "an integer vector", x)
  x <- as.integer(x)
  if (anyNA(x)) .tox_na_error(name)
  x
}

.tox_as_integer_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L) .tox_type_error(name, "an integer scalar", x)
  x <- as.integer(x)
  if (anyNA(x)) .tox_na_error(name)
  x
}

.tox_as_integer_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x)) .tox_type_error(name, "an integer matrix", x)
  storage.mode(x) <- "integer"
  if (anyNA(x)) .tox_na_error(name)
  x
}

.tox_as_logical <- function(x, name) {
  if (!is.logical(x)) .tox_type_error(name, "a logical vector", x)
  if (anyNA(x)) .tox_na_error(name)
  x
}

.tox_as_character <- function(x, name) {
  if (!is.character(x)) .tox_type_error(name, "a character vector", x)
  if (anyNA(x)) .tox_na_error(name)
  x
}

.tox_as_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L) .tox_type_error(name, "a single string", x)
  if (anyNA(x)) .tox_na_error(name)
  x
}

.tox_as_mode <- function(x, name, choices) {
  if (!is.character(x) || length(x) != 1L) .tox_type_error(name, "a single string", x)
  x <- tolower(x)
  if (!x %in% choices)
    .tox_raise(
      sprintf("'%s' must be one of %s, not \"%s\"",
              name, paste(sprintf('"%s"', choices), collapse = ", "), x),
      "tox_type_error", argument = name)
  x
}
'''
