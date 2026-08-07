# Generated. Do not edit.

#' Recommend workspace sizes based on Netlib exact formulas
#'
#' Computes the required sizes for integer and real workspace arrays.
#' These sizes depend on the dimensionality of the data and the maximum neighborhood size.
#'
#' Generated from the Fortran procedure \code{tox_loess_impl::tox_loess_required_workspace}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_dim a integer scalar. Dimensionality of the data
#' @param max_neighborhood_size a integer scalar. Maximum neighborhood size
#' @param save_factorization a logical scalar. Save matrix factorization flag
#' @return a named list with elements:
#'   \item{int_workspace_size}{a integer scalar. Required size of the integer workspace array}
#'   \item{real_workspace_size}{a integer scalar. Required size of the real workspace array}
#' @export
tox_loess_required_workspace <- function(n_dim, max_neighborhood_size, save_factorization) {
    n_dim <- .tox_as_integer_scalar(n_dim, "n_dim")
    max_neighborhood_size <- .tox_as_integer_scalar(max_neighborhood_size, "max_neighborhood_size")
    save_factorization <- .tox_as_logical(save_factorization, "save_factorization")
    .result <- .Call("tox_loess_required_workspace_call", n_dim, max_neighborhood_size, save_factorization)
    .arguments <- c("n_dim", "max_neighborhood_size", "int_workspace_size", "real_workspace_size", "save_factorization")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        int_workspace_size = .result$int_workspace_size,
        real_workspace_size = .result$real_workspace_size
    )
}
