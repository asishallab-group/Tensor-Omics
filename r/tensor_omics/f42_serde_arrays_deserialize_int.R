# Generated. Do not edit.

#' Deserialize a flat integer array from a file
#'
#' Generated from the Fortran procedure \code{f42_serde_arrays_deserialize_int::deserialize_int_helper}, whose argument names
#' are the ones an error message reports.
#'
#' @param filename a string. Name of the file
#' @return a integer vector. Pre-allocated array to read the data into
#' @export
deserialize_int_helper <- function(filename) {
    filename <- .tox_as_character(filename, "filename")
    .get_array_metadata_result <- get_array_metadata(filename = filename, dims_out_capacity = 5L)
    arr_shape <- .get_array_metadata_result$dims_out

    .result <- .Call("deserialize_int_helper_call", arr_shape, filename)
    .arguments <- c("arr", "n_elements", "arr_shape", "filename", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$arr
}
