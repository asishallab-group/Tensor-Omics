# Generated from tox_errors. Do not edit.

ARG_POS_FACTOR <- 10000L

# the codes, exported so a caller can compare against a name
ERR_OK <- 0L
ERR_FILE_OPEN <- 101L
ERR_READ_MAGIC <- 102L
ERR_READ_TYPE <- 103L
ERR_READ_NDIMS <- 104L
ERR_READ_DIMS <- 105L
ERR_READ_CHARLEN <- 106L
ERR_READ_DATA <- 107L
ERR_WRITE_MAGIC <- 112L
ERR_WRITE_TYPE <- 113L
ERR_WRITE_NDIMS <- 114L
ERR_WRITE_DIMS <- 115L
ERR_WRITE_CHARLEN <- 116L
ERR_WRITE_DATA <- 117L
ERR_FILE_ADD <- 121L
ERR_FILE_EXTRACT <- 122L
ERR_FILE_CLOSE <- 123L
ERR_INVALID_FORMAT <- 200L
ERR_INVALID_INPUT <- 201L
ERR_EMPTY_INPUT <- 202L
ERR_DIM_MISMATCH <- 203L
ERR_NAN_INF <- 204L
ERR_UNSUPPORTED_TYPE <- 205L
ERR_SIZE_MISMATCH <- 206L
ERR_TYPE_MISMATCH <- 207L
ERR_STRING_TOO_LONG <- 208L
ERR_IDX_OUT_OF_BOUNDS <- 209L
ERR_DIVISION_BY_ZERO <- 210L
ERR_ALLOC_FAIL <- 301L
ERR_POINTER_NULL <- 302L
ERR_UNIT_NOT_CONNECTED <- 5002L
ERR_INTERNAL <- 9001L
ERR_UNKNOWN <- 9999L

.tox_messages <- c(
    "101" = "could not open file",
    "102" = "could not read magic number",
    "103" = "could not read array type error",
    "104" = "could not read number of dimensions",
    "105" = "could not read array dimensions",
    "106" = "could not read character length",
    "107" = "could not read array data",
    "112" = "could not write magic number",
    "113" = "could not write array type error",
    "114" = "could not write number of dimensions",
    "115" = "could not write array dimensions",
    "116" = "could not write character length",
    "117" = "could not write array data",
    "121" = "Could not add file to archive",
    "122" = "Could not extract file from archive",
    "123" = "Failed to close the file",
    "200" = "invalid format detected",
    "201" = "invalid input arguments",
    "202" = "empty input arrays",
    "203" = "dimensions do not match expected shape",
    "204" = "NaN or Inf found where not allowed",
    "205" = "unsupported data type encountered",
    "206" = "Array size mismatch",
    "207" = "Array type read does not match expected type",
    "208" = "String exceeds buffer size",
    "209" = "Array index out of bounds",
    "210" = "Division by zero encountered",
    "301" = "memory allocation failed",
    "302" = "null pointer dereference",
    "5002" = "Fortran runtime error: unit not connected",
    "9001" = "unexpected internal state or logic error",
    "9999" = "unknown error"
)

.tox_names <- c(
    "101" = "ERR_FILE_OPEN",
    "102" = "ERR_READ_MAGIC",
    "103" = "ERR_READ_TYPE",
    "104" = "ERR_READ_NDIMS",
    "105" = "ERR_READ_DIMS",
    "106" = "ERR_READ_CHARLEN",
    "107" = "ERR_READ_DATA",
    "112" = "ERR_WRITE_MAGIC",
    "113" = "ERR_WRITE_TYPE",
    "114" = "ERR_WRITE_NDIMS",
    "115" = "ERR_WRITE_DIMS",
    "116" = "ERR_WRITE_CHARLEN",
    "117" = "ERR_WRITE_DATA",
    "121" = "ERR_FILE_ADD",
    "122" = "ERR_FILE_EXTRACT",
    "123" = "ERR_FILE_CLOSE",
    "200" = "ERR_INVALID_FORMAT",
    "201" = "ERR_INVALID_INPUT",
    "202" = "ERR_EMPTY_INPUT",
    "203" = "ERR_DIM_MISMATCH",
    "204" = "ERR_NAN_INF",
    "205" = "ERR_UNSUPPORTED_TYPE",
    "206" = "ERR_SIZE_MISMATCH",
    "207" = "ERR_TYPE_MISMATCH",
    "208" = "ERR_STRING_TOO_LONG",
    "209" = "ERR_IDX_OUT_OF_BOUNDS",
    "210" = "ERR_DIVISION_BY_ZERO",
    "301" = "ERR_ALLOC_FAIL",
    "302" = "ERR_POINTER_NULL",
    "5002" = "ERR_UNIT_NOT_CONNECTED",
    "9001" = "ERR_INTERNAL",
    "9999" = "ERR_UNKNOWN"
)

.tox_classes <- c(
    "101" = "tox_io_error",
    "102" = "tox_io_error",
    "103" = "tox_io_error",
    "104" = "tox_io_error",
    "105" = "tox_io_error",
    "106" = "tox_io_error",
    "107" = "tox_io_error",
    "112" = "tox_io_error",
    "113" = "tox_io_error",
    "114" = "tox_io_error",
    "115" = "tox_io_error",
    "116" = "tox_io_error",
    "117" = "tox_io_error",
    "121" = "tox_io_error",
    "122" = "tox_io_error",
    "123" = "tox_io_error",
    "200" = "tox_input_error",
    "201" = "tox_input_error",
    "202" = "tox_input_error",
    "203" = "tox_input_error",
    "204" = "tox_input_error",
    "205" = "tox_input_error",
    "206" = "tox_input_error",
    "207" = "tox_input_error",
    "208" = "tox_input_error",
    "209" = "tox_input_error",
    "210" = "tox_input_error",
    "301" = "tox_memory_error",
    "302" = "tox_memory_error",
    "5002" = "tox_runtime_error",
    "9001" = "tox_internal_error",
    "9999" = "tox_internal_error"
)

# status codes are outcomes, not failures, and never raise
.tox_statuses <- c(
)

#' Signal a tensor-omics error
#'
#' @param message the human-readable message
#' @param class the specific condition class, a subclass of `tox_error`
#' @param code the tox_errors code, position stripped
#' @param name the code's name, e.g. `"ERR_INVALID_INPUT"`
#' @param argument the offending argument, if the code named one
#' @keywords internal
.tox_raise <- function(message, class, code = NA_integer_,
                    name = NA_character_, argument = NA_character_) {
    cond <- structure(
        class = c(class, "tox_error", "error", "condition"),
        list(message = message, call = sys.call(-1),
                        code = code, name = name, argument = argument)
    )
    stop(cond)
}

#' @keywords internal
.tox_type_error <- function(argument, expected, value) {
    .tox_raise(
        sprintf("'%s' must be %s, not %s", argument, expected, class(value)[1]),
        "tox_type_error", argument = argument
    )
}

#' @keywords internal
.tox_na_error <- function(argument) {
    .tox_raise(
        sprintf("'%s' contains NA, which tensor-omics cannot represent", argument),
        "tox_na_error", argument = argument
    )
}

#' @keywords internal
.tox_shape_error <- function(argument, actual, other, expected) {
    .tox_raise(
        sprintf("'%s' has extent %d, but '%s' implies it should be %d",
        argument, actual, other, expected),
        "tox_shape_error", argument = argument
    )
}

#' Raise if a tox_errors code reports an error
#'
#' @param ierr the encoded error code, argument position included
#' @param sources positionally alongside `arguments`: the argument the caller
#'   actually passed, where the Fortran one was derived from it -- an extent read
#'   off a matrix names that matrix here. NA where the Fortran argument is the
#'   caller's own.
#' @param arguments the wrapped procedure's argument names, in order, so the
#'   message can name the offending one
#' @return invisibly, the name of the status code if the call reported one;
#'   status codes are outcomes, not failures, and never raise
#' @keywords internal
check_err_code <- function(ierr, arguments = character(), sources = character()) {
    code <- ierr %% ARG_POS_FACTOR
    if (code == 0L) return(invisible(NULL))

    key <- as.character(code)
    if (key %in% names(.tox_statuses))
        return(invisible(.tox_statuses[[key]]))

    arg_pos <- ierr %/% ARG_POS_FACTOR
    message <- if (key %in% names(.tox_messages)) .tox_messages[[key]]
                                                else sprintf("unmapped error code %d", code)
    argument <- NA_character_
    if (arg_pos > 0L) {
        if (arg_pos <= length(arguments)) {
            argument <- arguments[[arg_pos]]
            source <- if (arg_pos <= length(sources)) sources[[arg_pos]]
                                                else NA_character_
            if (!is.na(source)) {
                # the caller's word first, then the Fortran argument it
                # was derived from, which is what the signature calls it
                message <- sprintf("%s (argument '%s', via '%s')", message,
                                                    source, argument)
                argument <- source
            } else {
                message <- sprintf("%s (argument '%s')", message, argument)
            }
        } else {
            message <- sprintf("%s (argument %d)", message, arg_pos)
        }
    }

    class <- if (key %in% names(.tox_classes)) .tox_classes[[key]]
                                        else "tox_error"
    .tox_raise(message, class, code = code, name = .tox_names[key],
                                                argument = argument)
}
