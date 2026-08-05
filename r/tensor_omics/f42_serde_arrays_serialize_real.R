# Generated. Do not edit.

#' Subroutine to serialize a flat real array into a file
#'
#' Generated from the Fortran procedure \code{f42_serde_arrays_serialize_real::serialize_real_helper}, whose argument names
#' are the ones an error message reports.
#'
#' @param arr a numeric vector. Array to be serialized
#' @param filename a string. Name of the file to write to
#' @return invisibly `NULL`; called for its effect.
#' @export
serialize_real_helper <- function(arr, filename) {
    arr <- .tox_as_double_shaped(arr, "arr")
    filename <- .tox_as_character(filename, "filename")
    .result <- .Call("serialize_real_helper_call", arr, filename)
    .arguments <- c("arr", "n_elements", "arr_shape", "filename", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}
