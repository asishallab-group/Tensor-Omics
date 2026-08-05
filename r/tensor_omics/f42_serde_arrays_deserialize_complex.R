# Generated. Do not edit.

#' Deserialize a flat complex array from a file
#'
#' Generated from the Fortran procedure \code{f42_serde_arrays_deserialize_complex::deserialize_complex_helper}, whose argument names
#' are the ones an error message reports.
#'
#' @param filename a string. Name of the file
#' @return a complex vector. Pre-allocated array to read the data into
#' @export
deserialize_complex_helper <- function(filename) {
    filename <- .tox_as_character(filename, "filename")
    .get_array_metadata_result <- get_array_metadata(filename = filename, dims_out_capacity = 5L)
    arr_shape <- .get_array_metadata_result$dims_out

    .result <- .Call("deserialize_complex_helper_call", arr_shape, filename)
    .arguments <- c("arr", "n_elements", "arr_shape", "filename", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$arr
}
