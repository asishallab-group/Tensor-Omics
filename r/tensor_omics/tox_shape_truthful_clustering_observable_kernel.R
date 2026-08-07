# Generated. Do not edit.

#' Recommend LAPACK dgesdd workspace sizes for observable's economy-mode SVD
#'
#' The documented minimum-workspace formula for JOBZ='S' (see `man dgesdd`):
#' LWORK >= 4*min(M,N)**2 + 7*min(M,N), IWORK size = 8*min(M,N), where M=n_dimensions and
#' N=n_selected_member.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_observable_kernel::tox_stc_observable_svd_workspace}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_dimensions a integer scalar. Ambient dimension D
#' @param n_selected_member a integer scalar. Number of selected ensemble members
#' @return a named list with elements:
#'   \item{lwork}{a integer scalar. Recommended size of the real LAPACK workspace}
#'   \item{iwork_size}{a integer scalar. Recommended size of the integer LAPACK workspace}
#' @export
tox_stc_observable_svd_workspace <- function(n_dimensions, n_selected_member) {
    n_dimensions <- .tox_as_integer_scalar(n_dimensions, "n_dimensions")
    n_selected_member <- .tox_as_integer_scalar(n_selected_member, "n_selected_member")
    .result <- .Call("tox_stc_observable_svd_workspace_call", n_dimensions, n_selected_member)
    .arguments <- c("n_dimensions", "n_selected_member", "lwork", "iwork_size")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        lwork = .result$lwork,
        iwork_size = .result$iwork_size
    )
}
