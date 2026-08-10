# Generated. Do not edit.

#' Recommend LAPACK dgesvd workspace size for accept_ensemble's principal-angle SVDs
#'
#' Sized for the worst case across every reference comparison accept_ensemble performs (up
#' to n_dimensions-square, since a comparison's actual shared rank is always
#' <= n_dimensions): the documented minimum-workspace formula for a square M=N=n_dimensions
#' input with JOBU='N', JOBVT='N' (see `man dgesvd`): LWORK >= max(1, 5*n_dimensions). A
#' larger-than-required LWORK is always safe per LAPACK's own convention, so one size,
#' computed once, serves every one of accept_ensemble's (up to o+1) comparisons.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_accept_kernel::tox_stc_accept_ensemble_svd_workspace}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_dimensions a integer scalar. Ambient dimension D
#' @return a integer scalar. Recommended size of the real LAPACK workspace
#' @export
tox_stc_accept_ensemble_svd_workspace <- function(n_dimensions) {
    n_dimensions <- .tox_as_integer_scalar(n_dimensions, "n_dimensions")
    .result <- .Call("tox_stc_accept_ensemble_svd_workspace_call", n_dimensions)
    .arguments <- c("n_dimensions", "lwork")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$lwork
}
