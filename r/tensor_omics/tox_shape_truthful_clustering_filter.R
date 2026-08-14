# Generated. Do not edit.

#' Ensemble eligibility by Stop Condition
#'
#' `eligible(e) = allowed_stop_reasons(ensemble_stop_reason(e))` when `allowed_stop_reasons`
#' is present; all `TRUE` (no filtering) when absent.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_filter::filter_ensembles_by_stop_condition}, whose argument names
#' are the ones an error message reports.
#'
#' @param ensemble_stop_reason a integer vector. Per-ensemble Stop Condition, see `ensemble_identification`'s merged
#'   `ensemble_stop_reason` -- an index 1..4 into `allowed_stop_reasons` below, in the
#'   order `tox_shape_truthful_clustering_kernel`'s own `STOP_REASON_MAX_SIZE` (1),
#'   `STOP_REASON_REJECTED_AFTER_STABLE` (2), `STOP_REASON_REJECTED_IMMEDIATELY` (3),
#'   `STOP_REASON_FIXED_POINT` (4) -- not imported by name here, to avoid a circular
#'   module dependency (the parent module already `use`s the reconciliation module,
#'   which `use`s this one)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `4`.
#' @param allowed_stop_reasons a logical vector. Per-Stop-Condition eligibility, indexed as documented on `ensemble_stop_reason`
#'   above. Absent means no filtering (every Stop Condition allowed) -- deliberately
#'   nullable, not annotated with a generated default: the generator only evaluates
#'   constant *scalar* expressions for that annotation (`codegen_guide.md` section 5.5)
#' @return a logical vector. Per-ensemble eligibility from this criterion alone
#' @export
filter_ensembles_by_stop_condition <- function(ensemble_stop_reason, allowed_stop_reasons = NULL) {
    ensemble_stop_reason <- .tox_as_integer_vector(ensemble_stop_reason, "ensemble_stop_reason")
    if (!is.null(allowed_stop_reasons))
        allowed_stop_reasons <- .tox_as_logical(allowed_stop_reasons, "allowed_stop_reasons")
    .result <- .Call("filter_ensembles_by_stop_condition_call", ensemble_stop_reason, allowed_stop_reasons)
    .arguments <- c("n_ensembles", "ensemble_stop_reason", "allowed_stop_reasons", "eligible", "ierr")
    .sources <- c("ensemble_stop_reason", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$eligible
}

#' Ensemble eligibility by final intrinsic dimension
#'
#' `d_min <= d <= d_max`, both inclusive, each independently optional (an absent bound
#' contributes no constraint on that side). An ensemble with no final accepted state at all
#' (`ensemble_has_final(e)` false) is never eligible under this criterion once at least one
#' of `d_min`/`d_max` is supplied -- there is no `d` to judge. Both bounds absent is a true
#' no-op (every ensemble eligible, `ensemble_has_final` not even consulted), matching
#' `filter_ensembles_by_stop_condition_kernel`'s own "omitted means unconstrained"
#' convention. `d_min` could in principle default to `0` (a genuine constant
#' expression), but is left nullable like `d_max` (whose own natural default, `n_dimensions`,
#' is a runtime value and so cannot be a generated default) rather than have the two bounds
#' of one range behave asymmetrically.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_filter::filter_ensembles_by_dimension}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_dimensions a integer scalar. Ambient dimension D
#'   The minimum valid value is `2`.
#' @param ensemble_d_final a integer vector. Each ensemble's final accepted intrinsic dimension, see
#'   `ensemble_final_observable`
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param ensemble_has_final a logical vector. Whether each ensemble has a final accepted state at all, see
#'   `ensemble_final_observable`
#' @param d_min a integer scalar. Minimum tolerated final intrinsic dimension, inclusive
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param d_max a integer scalar. Maximum tolerated final intrinsic dimension, inclusive
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @return a logical vector. Per-ensemble eligibility from this criterion alone
#' @export
filter_ensembles_by_dimension <- function(n_dimensions, ensemble_d_final, ensemble_has_final, d_min = NULL, d_max = NULL) {
    n_dimensions <- .tox_as_integer_scalar(n_dimensions, "n_dimensions")
    ensemble_d_final <- .tox_as_integer_vector(ensemble_d_final, "ensemble_d_final")
    ensemble_has_final <- .tox_as_logical(ensemble_has_final, "ensemble_has_final")
    if (!is.null(d_min))
        d_min <- .tox_as_integer_scalar(d_min, "d_min")
    if (!is.null(d_max))
        d_max <- .tox_as_integer_scalar(d_max, "d_max")
    if (length(ensemble_has_final) != length(ensemble_d_final))
        .tox_shape_error("ensemble_has_final", length(ensemble_has_final), "ensemble_d_final", length(ensemble_d_final))

    .result <- .Call("filter_ensembles_by_dimension_call", n_dimensions, ensemble_d_final, ensemble_has_final, d_min, d_max)
    .arguments <- c("n_dimensions", "n_ensembles", "ensemble_d_final", "ensemble_has_final", "d_min", "d_max", "eligible", "ierr")
    .sources <- c(NA_character_, "ensemble_d_final", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$eligible
}

#' Ensemble eligibility by final classical variance explained
#'
#' Variance explained $= \sum_{j=1}^{d}\lambda_j / (\sum_{j=1}^{d}\lambda_j +
#' \text{normal\_error})$, the familiar PCA energy-ratio (not the scale/sqrt-based ratio
#' considered and rejected during this module's own design -- the classical, squared-units
#' form was chosen specifically for being the measure data scientists already expect), with
#' eigenvalues recovered from the final singular values exactly as `observable` itself does,
#' $\lambda_j = s_j^2/(k-1)$, and `normal_error` reusing
#' \code{\link{normal_error}} directly
#' rather than re-deriving its sum. An ensemble with no final accepted state, or whose final
#' size is too small for the $k-1$ denominator to be meaningful ($k \leq 1$), is never
#' eligible under this criterion once `var_explained_min` is supplied. Absent
#' `var_explained_min` is a true no-op, matching this module's other two filters.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_filter::filter_ensembles_by_var_explained}, whose argument names
#' are the ones an error message reports.
#'
#' @param ensemble_S_final a numeric matrix. Each ensemble's final accepted singular values, see `ensemble_final_observable`
#' @param ensemble_d_final a integer vector. Each ensemble's final accepted intrinsic dimension
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param ensemble_k_final a integer vector. Each ensemble's final accepted size
#'   The minimum valid value is `0`.
#' @param ensemble_has_final a logical vector. Whether each ensemble has a final accepted state at all
#' @param var_explained_min a numeric scalar. Minimum tolerated fraction of variance explained by the tangent subspace,
#'   inclusive
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @return a logical vector. Per-ensemble eligibility from this criterion alone
#' @export
filter_ensembles_by_var_explained <- function(ensemble_S_final, ensemble_d_final, ensemble_k_final, ensemble_has_final, var_explained_min = NULL) {
    ensemble_S_final <- .tox_as_double_matrix(ensemble_S_final, "ensemble_S_final")
    ensemble_d_final <- .tox_as_integer_vector(ensemble_d_final, "ensemble_d_final")
    ensemble_k_final <- .tox_as_integer_vector(ensemble_k_final, "ensemble_k_final")
    ensemble_has_final <- .tox_as_logical(ensemble_has_final, "ensemble_has_final")
    if (!is.null(var_explained_min))
        var_explained_min <- .tox_as_double_scalar(var_explained_min, "var_explained_min")
    if (length(ensemble_d_final) != dim(ensemble_S_final)[2])
        .tox_shape_error("ensemble_d_final", length(ensemble_d_final), "ensemble_S_final", dim(ensemble_S_final)[2])
    if (length(ensemble_k_final) != dim(ensemble_S_final)[2])
        .tox_shape_error("ensemble_k_final", length(ensemble_k_final), "ensemble_S_final", dim(ensemble_S_final)[2])
    if (length(ensemble_has_final) != dim(ensemble_S_final)[2])
        .tox_shape_error("ensemble_has_final", length(ensemble_has_final), "ensemble_S_final", dim(ensemble_S_final)[2])

    .result <- .Call("filter_ensembles_by_var_explained_call", ensemble_S_final, ensemble_d_final, ensemble_k_final, ensemble_has_final, var_explained_min)
    .arguments <- c("n_dimensions", "n_ensembles", "ensemble_S_final", "ensemble_d_final", "ensemble_k_final", "ensemble_has_final", "var_explained_min", "eligible", "ierr")
    .sources <- c("ensemble_S_final", "ensemble_S_final", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$eligible
}

#' Combined ensemble eligibility for `merge_to_super_ensembles`
#'
#' Extracts each ensemble's final accepted state once (via `ensemble_final_observable`,
#' shared so this and `tox_stc_json`'s own reporting never derive it two different ways),
#' then calls each of this module's three per-criterion filters and combines their masks
#' with a plain logical AND. Also returns the three individual masks, not just the
#' combination -- so a caller (the report, in particular) can say *which* criterion excluded
#' a given ensemble, not merely that one did. Supplying none of `allowed_stop_reasons`/
#' `d_min`/`d_max`/`var_explained_min` makes every ensemble eligible (all four masks
#' all-`TRUE`), matching each individual filter's own no-op convention.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_filter::filter_ensembles}, whose argument names
#' are the ones an error message reports.
#'
#' @param ensemble_U_history a numeric array of rank 4. Per-ensemble trailing tangent+normal bases, see `ensemble_identification`'s
#'   merged output
#' @param ensemble_d_history a integer matrix. Per-ensemble trailing intrinsic dimensions
#' @param ensemble_S_history a numeric array of rank 3. Per-ensemble trailing singular values
#' @param ensemble_mu_history a numeric array of rank 3. Per-ensemble trailing centers
#' @param ensemble_G_history a numeric matrix. Per-ensemble trailing spectral gaps
#' @param ensemble_k_history a integer matrix. Per-ensemble trailing sizes
#' @param ensemble_accepted_history a logical matrix. Whether the growth iteration retained in each history column was itself accepted
#' @param ensemble_stop_reason a integer vector. Per-ensemble Stop Condition, see `filter_ensembles_by_stop_condition_kernel`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `4`.
#' @param allowed_stop_reasons a logical vector. See `filter_ensembles_by_stop_condition_kernel`
#' @param d_min a integer scalar. See `filter_ensembles_by_dimension_kernel`
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param d_max a integer scalar. See `filter_ensembles_by_dimension_kernel`
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param var_explained_min a numeric scalar. See `filter_ensembles_by_var_explained_kernel`
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @return a named list with elements:
#'   \item{eligible}{a logical vector. Combined eligibility: `TRUE` only where all three per-criterion masks are}
#'   \item{eligible_by_stop_condition}{a logical vector. See `filter_ensembles_by_stop_condition_kernel`}
#'   \item{eligible_by_dimension}{a logical vector. See `filter_ensembles_by_dimension_kernel`}
#'   \item{eligible_by_var_explained}{a logical vector. See `filter_ensembles_by_var_explained_kernel`}
#' @export
filter_ensembles <- function(ensemble_U_history, ensemble_d_history, ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, ensemble_accepted_history, ensemble_stop_reason, allowed_stop_reasons = NULL, d_min = NULL, d_max = NULL, var_explained_min = NULL) {
    ensemble_U_history <- .tox_as_double_array(ensemble_U_history, "ensemble_U_history", 4L)
    ensemble_d_history <- .tox_as_integer_matrix(ensemble_d_history, "ensemble_d_history")
    ensemble_S_history <- .tox_as_double_array(ensemble_S_history, "ensemble_S_history", 3L)
    ensemble_mu_history <- .tox_as_double_array(ensemble_mu_history, "ensemble_mu_history", 3L)
    ensemble_G_history <- .tox_as_double_matrix(ensemble_G_history, "ensemble_G_history")
    ensemble_k_history <- .tox_as_integer_matrix(ensemble_k_history, "ensemble_k_history")
    ensemble_accepted_history <- .tox_as_logical(ensemble_accepted_history, "ensemble_accepted_history")
    ensemble_stop_reason <- .tox_as_integer_vector(ensemble_stop_reason, "ensemble_stop_reason")
    if (!is.null(allowed_stop_reasons))
        allowed_stop_reasons <- .tox_as_logical(allowed_stop_reasons, "allowed_stop_reasons")
    if (!is.null(d_min))
        d_min <- .tox_as_integer_scalar(d_min, "d_min")
    if (!is.null(d_max))
        d_max <- .tox_as_integer_scalar(d_max, "d_max")
    if (!is.null(var_explained_min))
        var_explained_min <- .tox_as_double_scalar(var_explained_min, "var_explained_min")
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
    if (length(ensemble_stop_reason) != dim(ensemble_U_history)[4])
        .tox_shape_error("ensemble_stop_reason", length(ensemble_stop_reason), "ensemble_U_history", dim(ensemble_U_history)[4])

    .result <- .Call("filter_ensembles_call", ensemble_U_history, ensemble_d_history, ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, ensemble_accepted_history, ensemble_stop_reason, allowed_stop_reasons, d_min, d_max, var_explained_min)
    .arguments <- c("n_dimensions", "o", "n_ensembles", "ensemble_U_history", "ensemble_d_history", "ensemble_S_history", "ensemble_mu_history", "ensemble_G_history", "ensemble_k_history", "ensemble_accepted_history", "ensemble_stop_reason", "allowed_stop_reasons", "d_min", "d_max", "var_explained_min", "eligible", "eligible_by_stop_condition", "eligible_by_dimension", "eligible_by_var_explained", "ierr")
    .sources <- c("ensemble_U_history", "ensemble_U_history", "ensemble_U_history", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        eligible = .result$eligible,
        eligible_by_stop_condition = .result$eligible_by_stop_condition,
        eligible_by_dimension = .result$eligible_by_dimension,
        eligible_by_var_explained = .result$eligible_by_var_explained
    )
}
