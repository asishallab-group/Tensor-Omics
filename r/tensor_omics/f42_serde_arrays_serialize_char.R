# Generated. Do not edit.

#' Subroutine to serialize a flat character array into a file
#'
#' Generated from the Fortran procedure \code{f42_serde_arrays_serialize_char::serialize_char_helper}, whose argument names
#' are the ones an error message reports.
#'
#' @param arr a character vector. Array to be serialized
#' @param filename a string. Name of the file to write to
#' @return invisibly `NULL`; called for its effect.
#' @export
serialize_char_helper <- function(arr, filename) {
    arr <- .tox_as_character_shaped(arr, "arr")
    filename <- .tox_as_character(filename, "filename")
    .result <- .Call("serialize_char_helper_call", arr, filename)
    .arguments <- c("arr", "n_strings", "arr_shape", "filename", "ierr")
    .sources <- c(NA_character_, "arr", "arr", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    invisible(NULL)
}
