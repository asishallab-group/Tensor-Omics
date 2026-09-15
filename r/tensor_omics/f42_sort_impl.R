# Generated. Do not edit.

#' Sort a real vector, returning the permutation (ascending, NaN last)
#'
#' Generated from the Fortran procedure \code{f42_sort_impl::sort_real_get_perm}, whose argument names
#' are the ones an error message reports.
#'
#' @param array a numeric vector. values to sort
#' @return a integer vector. permutation that sorts `array` ascending, NaN last
#' @export
sort_real_get_perm <- function(array) {
    array <- .tox_as_double_vector(array, "array")
    .result <- .Call("sort_real_get_perm_call", array)
    .arguments <- c("array", "n", "perm", "ierr")
    .sources <- c(NA_character_, "array", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$perm
}
