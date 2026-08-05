# Generated. Do not edit.

#' Get the metadata of an array file
#'
#' Generated from the Fortran module \code{f42_serde_arrays_utils}.
#'
#' @param filename a string. Name of the file
#' @param dims_out_capacity a integer scalar. Capacity of the dims_out array
#' @return a named list with elements:
#'   \item{dims_out}{a integer vector. Array to store output dimensions
#'     The first `ndims` elements will hold the results.}
#'   \item{type_code}{a integer scalar. Type code of the serialized array}
#' @export
get_array_metadata <- function(filename, dims_out_capacity) {
    filename <- .tox_as_character(filename, "filename")
    dims_out_capacity <- .tox_as_integer_scalar(dims_out_capacity, "dims_out_capacity")
    .result <- .Call("get_array_metadata_call", filename, dims_out_capacity)
    .arguments <- c("filename", "dims_out", "dims_out_capacity", "ndims", "type_code", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        dims_out = utils::head(.result$dims_out, .result$ndims),
        type_code = .result$type_code
    )
}
