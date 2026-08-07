# Generated. Do not edit.

#' Recommend LAPACK dgesvd workspace size for accept_ensemble's principal-angle SVD
#'
#' The documented minimum-workspace formula for a square M=N=min(d_t,d_tp1) input with
#' JOBU='N', JOBVT='N' (see `man dgesvd`): LWORK >= max(1, 5*min(M,N)).
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_accept_kernel::tox_stc_accept_ensemble_svd_workspace}, whose argument names
#' are the ones an error message reports.
#'
#' @param d_t a integer scalar. Ensemble's intrinsic dimension at t
#' @param d_tp1 a integer scalar. Ensemble's intrinsic dimension at t+1
#' @return a integer scalar. Recommended size of the real LAPACK workspace
#' @export
tox_stc_accept_ensemble_svd_workspace <- function(d_t, d_tp1) {
    d_t <- .tox_as_integer_scalar(d_t, "d_t")
    d_tp1 <- .tox_as_integer_scalar(d_tp1, "d_tp1")
    .result <- .Call("tox_stc_accept_ensemble_svd_workspace_call", d_t, d_tp1)
    .arguments <- c("d_t", "d_tp1", "lwork")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$lwork
}
