# Generated. Do not edit.

#' Subroutine to serialize a flat complex array into a file
#'
#' @param arr a complex vector. Array to be serialized
#' @param filename a string. Name of the file to write to
#' @return invisibly `NULL`; called for its effect.
#'
#' Generated from the Fortran procedure \code{f42_serde_arrays_serialize_complex::serialize_complex_helper}.
#' @export
serialize_complex_helper <- function(arr, filename) {
    arr <- .tox_as_complex_shaped(arr, "arr")
    filename <- .tox_as_character(filename, "filename")
    .result <- .serialize_complex_helper_rcpp(arr, filename)
    .arguments <- c("arr", "n_elements", "arr_shape", "filename", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}
