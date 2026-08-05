# Generated. Do not edit.

#' Subroutine to serialize a flat logical array into a file
#'
#' Generated from the Fortran procedure \code{f42_serde_arrays_serialize_logical::serialize_logical_helper}, whose argument names
#' are the ones an error message reports.
#'
#' @param arr a logical vector. Array to be serialized
#' @param filename a string. Name of the file to write to
#' @return invisibly `NULL`; called for its effect.
#' @export
serialize_logical_helper <- function(arr, filename) {
    arr <- .tox_as_logical_shaped(arr, "arr")
    filename <- .tox_as_character(filename, "filename")
    .result <- .Call("serialize_logical_helper_call", arr, filename)
    .arguments <- c("arr", "n_elements", "arr_shape", "filename", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}
