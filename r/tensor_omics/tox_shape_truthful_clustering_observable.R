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
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_observable::observable}, whose argument names
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

#' Each ensemble's final *accepted* growth-history state
#'
#' Not simply the last *populated* history column: `stc_push_ensemble_history`
#' (`tox_shape_truthful_clustering_impl`) also pushes a *rejected* final candidate right
#' before `ensemble_identification` halts growth via `STOP_REASON_REJECTED_IMMEDIATELY`/
#' `STOP_REASON_REJECTED_AFTER_STABLE`, so the last populated column is, in exactly those
#' two cases, the discarded candidate's geometry, not the ensemble's real final state.
#' Scans each ensemble's history backward for the last column that is both populated
#' (`ensemble_k_history /= 0`) and itself accepted (`ensemble_accepted_history`), and
#' slices `U`/`d`/`S`/`mu`/`G`/`k` out at that column. `ensemble_has_final` is `FALSE`
#' (all other outputs zero for that ensemble) only when no column qualifies at all --
#' possible for `STOP_REASON_MAX_SIZE` firing at the bootstrap step itself, before any
#' genuine SVD ever ran for that seed, and rarely when a small `o` lets a rejected push
#' evict every accepted entry the window ever held (`ensemble_U_first`/`ensemble_d_first`
#' still hold the bootstrap iteration in that case, but this kernel does not fall back to
#' them, since there is no `G_first`/`mu_first` counterpart to complete a fallback "final
#' state" from -- see `misc/mod_STC.md`, "Ensemble identification"). Consolidates what was
#' previously a private, single-consumer helper in `tox_stc_json.F90`
#' (`stc_last_accepted_history_index`) into a proper, independently testable kernel now
#' that `tox_shape_truthful_clustering_filter_impl` needs the exact same extraction --
#' placed here, not in the parent module, specifically so both of those (and any future
#' sibling) can depend on it without a circular module dependency (the parent module
#' already `use`s `tox_shape_truthful_clustering_reconciliation_impl`, which will in turn
#' `use` the filter kernel module, which needs this).
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_observable::ensemble_final_observable}, whose argument names
#' are the ones an error message reports.
#'
#' @param ensemble_U_history a numeric array of rank 4. Per-ensemble trailing tangent+normal bases, see `ensemble_identification`'s
#'   merged `ensemble_U_history`
#' @param ensemble_d_history a integer matrix. Per-ensemble trailing intrinsic dimensions
#' @param ensemble_S_history a numeric array of rank 3. Per-ensemble trailing singular values
#' @param ensemble_mu_history a numeric array of rank 3. Per-ensemble trailing centers
#' @param ensemble_G_history a numeric matrix. Per-ensemble trailing spectral gaps
#' @param ensemble_k_history a integer matrix. Per-ensemble trailing sizes; 0 marks an unpopulated column
#' @param ensemble_accepted_history a logical matrix. Whether the growth iteration retained in each history column was itself
#'   accepted -- see this kernel's own summary above
#' @return a named list with elements:
#'   \item{ensemble_U_final}{a numeric array of rank 3. Each ensemble's final accepted tangent+normal basis; zero when
#'     `ensemble_has_final` is `FALSE` for that ensemble}
#'   \item{ensemble_d_final}{a integer vector. Each ensemble's final accepted intrinsic dimension; zero when
#'     `ensemble_has_final` is `FALSE` for that ensemble}
#'   \item{ensemble_S_final}{a numeric matrix. Each ensemble's final accepted singular values; zero when
#'     `ensemble_has_final` is `FALSE` for that ensemble}
#'   \item{ensemble_mu_final}{a numeric matrix. Each ensemble's final accepted center; zero when `ensemble_has_final` is
#'     `FALSE` for that ensemble}
#'   \item{ensemble_G_final}{a numeric vector. Each ensemble's final accepted spectral gap; zero when `ensemble_has_final` is
#'     `FALSE` for that ensemble}
#'   \item{ensemble_k_final}{a integer vector. Each ensemble's final accepted size; zero when `ensemble_has_final` is
#'     `FALSE` for that ensemble}
#'   \item{ensemble_has_final}{a logical vector. Whether any history column at all qualifies as this ensemble's final accepted
#'     state -- see this kernel's own summary above for the (rare) `FALSE` cases}
#'   \item{ensemble_final_index}{a integer vector. The history column each `_final` output was sliced from (0 when
#'     `ensemble_has_final` is `FALSE`) -- also, since every column 1..this index is
#'     itself guaranteed accepted (only ever the single *last* populated column can be
#'     the rejected candidate this kernel's own summary describes), this doubles as the
#'     count of genuinely accepted, plottable history columns, for callers (e.g.
#'     `tox_stc_json`'s own `observable_history`) that need to iterate the whole
#'     trailing window, not just its final entry}
#' @export
ensemble_final_observable <- function(ensemble_U_history, ensemble_d_history, ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, ensemble_accepted_history) {
    ensemble_U_history <- .tox_as_double_array(ensemble_U_history, "ensemble_U_history", 4L)
    ensemble_d_history <- .tox_as_integer_matrix(ensemble_d_history, "ensemble_d_history")
    ensemble_S_history <- .tox_as_double_array(ensemble_S_history, "ensemble_S_history", 3L)
    ensemble_mu_history <- .tox_as_double_array(ensemble_mu_history, "ensemble_mu_history", 3L)
    ensemble_G_history <- .tox_as_double_matrix(ensemble_G_history, "ensemble_G_history")
    ensemble_k_history <- .tox_as_integer_matrix(ensemble_k_history, "ensemble_k_history")
    ensemble_accepted_history <- .tox_as_logical(ensemble_accepted_history, "ensemble_accepted_history")
    if (dim(ensemble_S_history)[1] != dim(ensemble_U_history)[1])
        .tox_shape_error("ensemble_S_history", dim(ensemble_S_history)[1], "ensemble_U_history", dim(ensemble_U_history)[1])
    if (dim(ensemble_mu_history)[1] != dim(ensemble_U_history)[1])
        .tox_shape_error("ensemble_mu_history", dim(ensemble_mu_history)[1], "ensemble_U_history", dim(ensemble_U_history)[1])
    if (dim(ensemble_d_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_d_history", dim(ensemble_d_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_S_history)[2] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_S_history", dim(ensemble_S_history)[2], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_mu_history)[2] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_mu_history", dim(ensemble_mu_history)[2], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_G_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_G_history", dim(ensemble_G_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_k_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_k_history", dim(ensemble_k_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_accepted_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_accepted_history", dim(ensemble_accepted_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_d_history)[2] != dim(ensemble_U_history)[4])
        .tox_shape_error("ensemble_d_history", dim(ensemble_d_history)[2], "ensemble_U_history", dim(ensemble_U_history)[4])
    if (dim(ensemble_S_history)[3] != dim(ensemble_U_history)[4])
        .tox_shape_error("ensemble_S_history", dim(ensemble_S_history)[3], "ensemble_U_history", dim(ensemble_U_history)[4])
    if (dim(ensemble_mu_history)[3] != dim(ensemble_U_history)[4])
        .tox_shape_error("ensemble_mu_history", dim(ensemble_mu_history)[3], "ensemble_U_history", dim(ensemble_U_history)[4])
    if (dim(ensemble_G_history)[2] != dim(ensemble_U_history)[4])
        .tox_shape_error("ensemble_G_history", dim(ensemble_G_history)[2], "ensemble_U_history", dim(ensemble_U_history)[4])
    if (dim(ensemble_k_history)[2] != dim(ensemble_U_history)[4])
        .tox_shape_error("ensemble_k_history", dim(ensemble_k_history)[2], "ensemble_U_history", dim(ensemble_U_history)[4])
    if (dim(ensemble_accepted_history)[2] != dim(ensemble_U_history)[4])
        .tox_shape_error("ensemble_accepted_history", dim(ensemble_accepted_history)[2], "ensemble_U_history", dim(ensemble_U_history)[4])

    .result <- .Call("ensemble_final_observable_call", ensemble_U_history, ensemble_d_history, ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, ensemble_accepted_history)
    .arguments <- c("n_dimensions", "o", "n_ensembles", "ensemble_U_history", "ensemble_d_history", "ensemble_S_history", "ensemble_mu_history", "ensemble_G_history", "ensemble_k_history", "ensemble_accepted_history", "ensemble_U_final", "ensemble_d_final", "ensemble_S_final", "ensemble_mu_final", "ensemble_G_final", "ensemble_k_final", "ensemble_has_final", "ensemble_final_index", "ierr")
    .sources <- c("ensemble_U_history", "ensemble_U_history", "ensemble_U_history", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        ensemble_U_final = .result$ensemble_U_final,
        ensemble_d_final = .result$ensemble_d_final,
        ensemble_S_final = .result$ensemble_S_final,
        ensemble_mu_final = .result$ensemble_mu_final,
        ensemble_G_final = .result$ensemble_G_final,
        ensemble_k_final = .result$ensemble_k_final,
        ensemble_has_final = .result$ensemble_has_final,
        ensemble_final_index = .result$ensemble_final_index
    )
}
