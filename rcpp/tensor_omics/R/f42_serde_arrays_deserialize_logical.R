# Generated. Do not edit.

#' Deserialize a flat logical array from a file
#'
#' @param filename a string. Name of the file
#' @return Pre-allocated array to read the data into
#'
#' Generated from the Fortran procedure \code{f42_serde_arrays_deserialize_logical::deserialize_logical_helper}.
#' @export
deserialize_logical_helper <- function(filename) {
    filename <- .tox_as_character(filename, "filename")
    .get_array_metadata_result <- get_array_metadata(filename = filename, dims_out_capacity = 5L)
    arr_shape <- .get_array_metadata_result$dims_out

    .result <- .deserialize_logical_helper_rcpp(arr_shape, filename)
    .arguments <- c("arr", "n_elements", "arr_shape", "filename", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$arr
}
