# Generated. Do not edit.

#' Mean squared residual of an ensemble's members off its tangent subspace
#'
#' No pass over the ensemble's member vectors is required; the sum is already implied by
#' the singular value decomposition \code{\link{observable}}
#' computes.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_observable::normal_error}, whose argument names
#' are the ones an error message reports.
#'
#' @param d a integer scalar. Intrinsic (tangent) dimension of the ensemble
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param eigenvalues a numeric vector. Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
#'   The minimum valid value is `0.0`.
#' @return a numeric scalar. Mean squared residual off the d-dimensional tangent subspace
#' @export
normal_error <- function(d, eigenvalues) {
    d <- .tox_as_integer_scalar(d, "d")
    eigenvalues <- .tox_as_double_vector(eigenvalues, "eigenvalues")
    .result <- .Call("normal_error_call", d, eigenvalues)
    .arguments <- c("d", "eigenvalues", "n_dimensions", "normal_error_value", "ierr")
    .sources <- c(NA_character_, NA_character_, "eigenvalues", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$normal_error_value
}

#' Extent along each tangent direction of an ensemble's tangent subspace
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_observable::tangent_scales}, whose argument names
#' are the ones an error message reports.
#'
#' @param d a integer scalar. Intrinsic (tangent) dimension of the ensemble
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param eigenvalues a numeric vector. Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
#'   The minimum valid value is `0.0`.
#' @return a numeric vector. Extent along each of the d tangent directions
#' @export
tangent_scales <- function(d, eigenvalues) {
    d <- .tox_as_integer_scalar(d, "d")
    eigenvalues <- .tox_as_double_vector(eigenvalues, "eigenvalues")
    .result <- .Call("tangent_scales_call", d, eigenvalues)
    .arguments <- c("d", "eigenvalues", "n_dimensions", "tangent_scales_value", "ierr")
    .sources <- c("tangent_scales_value", NA_character_, "eigenvalues", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$tangent_scales_value
}

#' The tuple (U, d, G, mu, normal_error, tangent_scales) for an ensemble
#'
#' `U` and `eigenvalues` are zero-padded to the full ambient dimension `n_dimensions`:
#' the economy SVD only yields `rank = min(n_dimensions, n_selected_member)` genuine
#' columns/values, less than `n_dimensions` whenever an ensemble is smaller than the
#' ambient space (typical early in growth). This keeps the output shape fixed regardless
#' of ensemble size, and slots directly into `normal_error`/`tangent_scales`'s existing
#' `n_dimensions`-length interface. `ierr` is set only if the LAPACK SVD fails to
#' converge -- not a condition any input check could foresee.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_observable::observable_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param member_selection_mask a logical vector. Ensemble membership over the full dataset
#' @return a named list with elements:
#'   \item{U}{a numeric matrix. Tangent+normal basis, zero-padded beyond rank}
#'   \item{eigenvalues}{a numeric vector. Covariance eigenvalues, descending, zero-padded beyond rank}
#'   \item{mu}{a numeric vector. Ensemble center}
#'   \item{d}{a integer scalar. Estimated intrinsic (tangent) dimension}
#'   \item{G}{a numeric scalar. Spectral gap at d}
#'   \item{normal_error_value}{a numeric scalar. Mean squared residual off the tangent subspace}
#'   \item{tangent_scales_value}{a numeric vector. Extent along each tangent direction, zero-padded beyond d}
#' @export
observable <- function(vectors, member_selection_mask) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    member_selection_mask <- .tox_as_logical(member_selection_mask, "member_selection_mask")
    if (length(member_selection_mask) != dim(vectors)[2])
        .tox_shape_error("member_selection_mask", length(member_selection_mask), "vectors", dim(vectors)[2])

    .result <- .Call("observable_call", vectors, member_selection_mask)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "member_selection_mask", "n_selected_member", "U", "eigenvalues", "mu", "d", "G", "normal_error_value", "tangent_scales_value", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, "member_selection_mask", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        U = .result$U,
        eigenvalues = .result$eigenvalues,
        mu = .result$mu,
        d = .result$d,
        G = .result$G,
        normal_error_value = .result$normal_error_value,
        tangent_scales_value = .result$tangent_scales_value
    )
}
