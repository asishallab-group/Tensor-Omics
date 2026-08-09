# Generated. Do not edit.

#' Recommend LAPACK workspace sizes for estimate_stc_parameters' SVD calls
#'
#' Worst-case sizing for both `observable`'s dgesdd (an EA cloud can be as large as
#' n_vectors) and the pairwise principal-angle dgesvd (shared rank at most n_dimensions),
#' see `tox_stc_observable_svd_workspace` and `tox_stc_accept_ensemble_svd_workspace` for
#' the individual formulas this combines.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_parameter_estimation_kernel::tox_stc_estimate_parameters_svd_workspace}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_dimensions a integer scalar. Ambient dimension D
#' @param n_vectors a integer scalar. Number of input vectors N
#' @return a named list with elements:
#'   \item{lwork_observable}{a integer scalar. Recommended size of observable's real LAPACK workspace (worst case)}
#'   \item{iwork_size}{a integer scalar. Recommended size of observable's integer LAPACK workspace (worst case)}
#'   \item{lwork_angle}{a integer scalar. Recommended size of the pairwise principal-angle LAPACK workspace (worst case)}
#' @export
tox_stc_estimate_parameters_svd_workspace <- function(n_dimensions, n_vectors) {
    n_dimensions <- .tox_as_integer_scalar(n_dimensions, "n_dimensions")
    n_vectors <- .tox_as_integer_scalar(n_vectors, "n_vectors")
    .result <- .Call("tox_stc_estimate_parameters_svd_workspace_call", n_dimensions, n_vectors)
    .arguments <- c("n_dimensions", "n_vectors", "lwork_observable", "iwork_size", "lwork_angle")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        lwork_observable = .result$lwork_observable,
        iwork_size = .result$iwork_size,
        lwork_angle = .result$lwork_angle
    )
}
