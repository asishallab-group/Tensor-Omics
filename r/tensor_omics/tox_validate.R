# Generated. Do not edit.
#
# Coercion helpers for the R binding. Each checks the type and coerces, copying only
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

# An array whose shape travels in a separate argument keeps that shape through
# coercion: as.integer() and friends drop the dim attribute, and the shape is read off
# the array after this point. storage.mode() converts in place instead.
.tox_as_double_shaped <- function(x, name) {
  if (!is.numeric(x)) .tox_type_error(name, "a numeric array", x)
  storage.mode(x) <- "double"
  x
}

.tox_as_integer_shaped <- function(x, name) {
  if (!is.numeric(x)) .tox_type_error(name, "an integer array", x)
  storage.mode(x) <- "integer"
  if (anyNA(x)) .tox_na_error(name)
  x
}

.tox_as_complex_shaped <- function(x, name) {
  if (!is.numeric(x) && !is.complex(x)) .tox_type_error(name, "a complex array", x)
  storage.mode(x) <- "complex"
  x
}

.tox_as_logical_shaped <- function(x, name) {
  if (!is.logical(x) && !is.numeric(x)) .tox_type_error(name, "a logical array", x)
  dims <- dim(x)
  x <- as.logical(x)
  if (anyNA(x)) .tox_na_error(name)
  dim(x) <- dims
  x
}

.tox_as_character_shaped <- function(x, name) {
  if (!is.character(x)) .tox_type_error(name, "a character array", x)
  if (anyNA(x)) .tox_na_error(name)
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
