# Generated. Do not edit.

#' Subroutine to serialize a flat character array into a file
#'
#' @param arr a character vector. Array to be serialized
#' @param filename a string. Name of the file to write to
#' @return invisibly `NULL`; called for its effect.
#'
#' Generated from the Fortran module \code{f42_serde_arrays_serialize_char}.
#' @export
serialize_char_helper <- function(arr, filename) {
    arr <- .tox_as_character_shaped(arr, "arr")
    filename <- .tox_as_character(filename, "filename")
    .result <- .Call("serialize_char_helper_call", arr, filename)
    .arguments <- c("arr", "n_strings", "arr_shape", "filename", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}
