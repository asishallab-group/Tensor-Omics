# Generated. Do not edit.

#' Whether a grown ensemble at t+1 is still compatible with its own state at t
#'
#' Three criteria, all must hold: (1) principal angles between the d-dimensional tangent
#' bases, via `dgesvd` on M = U_t(:,1:d)^T U_tp1(:,1:d), whose singular values are
#' cos(alpha_i) directly -- but only when d_t == d_tp1: when the estimated intrinsic
#' dimension itself changed, the two tangent bases don't share a common dimension to
#' compare angles over at all, so this criterion is skipped (no SVD is computed) and
#' treated as vacuously satisfied; criterion (2) is what actually judges whether that
#' change in d is acceptable. (2) |d_tp1 - d_t| <= d_max. (3) |log(G_tp1/G_t)| <= G_max.
#' `ierr` is set only if the LAPACK SVD fails to converge -- not a condition any input
#' check could foresee.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_accept::accept_ensemble_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param U_t a numeric matrix. Ensemble's tangent+normal basis at t, see `observable`
#' @param d_t a integer scalar. Ensemble's intrinsic dimension at t
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param G_t a numeric scalar. Ensemble's spectral gap at t
#'   The minimum valid value is `above(0.0)`.
#' @param U_tp1 a numeric matrix. Ensemble's tangent+normal basis at t+1
#' @param d_tp1 a integer scalar. Ensemble's intrinsic dimension at t+1
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param G_tp1 a numeric scalar. Ensemble's spectral gap at t+1
#'   The minimum valid value is `above(0.0)`.
#' @param alpha_max a numeric scalar. Maximum tolerated principal angle (radians)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `2.0 * atan(1.0)`.
#' @param d_max a integer scalar. Maximum tolerated change in intrinsic dimension
#'   The minimum valid value is `0`.
#' @param G_max a numeric scalar. Maximum tolerated |log(G_tp1/G_t)|
#'   The minimum valid value is `0.0`.
#' @return a logical scalar. TRUE if all three acceptance criteria are satisfied
#' @export
accept_ensemble <- function(U_t, d_t, G_t, U_tp1, d_tp1, G_tp1, alpha_max, d_max, G_max) {
    U_t <- .tox_as_double_matrix(U_t, "U_t")
    d_t <- .tox_as_integer_scalar(d_t, "d_t")
    G_t <- .tox_as_double_scalar(G_t, "G_t")
    U_tp1 <- .tox_as_double_matrix(U_tp1, "U_tp1")
    d_tp1 <- .tox_as_integer_scalar(d_tp1, "d_tp1")
    G_tp1 <- .tox_as_double_scalar(G_tp1, "G_tp1")
    alpha_max <- .tox_as_double_scalar(alpha_max, "alpha_max")
    d_max <- .tox_as_integer_scalar(d_max, "d_max")
    G_max <- .tox_as_double_scalar(G_max, "G_max")
    if (dim(U_tp1)[1] != dim(U_t)[1])
        .tox_shape_error("U_tp1", dim(U_tp1)[1], "U_t", dim(U_t)[1])

    .result <- .Call("accept_ensemble_call", U_t, d_t, G_t, U_tp1, d_tp1, G_tp1, alpha_max, d_max, G_max)
    .arguments <- c("n_dimensions", "U_t", "d_t", "G_t", "U_tp1", "d_tp1", "G_tp1", "alpha_max", "d_max", "G_max", "is_accepted", "ierr")
    .sources <- c("U_t", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$is_accepted
}
