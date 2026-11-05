# Generated. Do not edit.
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
