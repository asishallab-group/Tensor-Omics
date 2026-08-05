# Generated. Do not edit.

#' Subroutine to deserialize a flat character array from a file
#'
#' Generated from the Fortran module \code{f42_serde_arrays_deserialize_char}.
#'
#' @param filename a string. Name of the file
#' @return a character vector. Pre-allocated array to read the data into
#' @export
deserialize_char_helper <- function(filename) {
    filename <- .tox_as_character(filename, "filename")
    .get_array_metadata_result <- get_array_metadata(filename = filename, dims_out_capacity = 5L)
    strlen <- .get_array_metadata_result$type_code
    arr_shape <- .get_array_metadata_result$dims_out

    .result <- .Call("deserialize_char_helper_call", strlen, arr_shape, filename)
    .arguments <- c("arr", "n_strings", "strlen", "arr_shape", "filename", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$arr
}
