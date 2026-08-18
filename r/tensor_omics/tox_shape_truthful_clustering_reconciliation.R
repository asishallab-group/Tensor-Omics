# Generated. Do not edit.

#' Filter eligible ensembles, then group/report their intersections
#'
#' See this module's own header comment for why this is a two-call orchestrator, and
#' `tox_shape_truthful_clustering_filter_impl`'s own `filter_ensembles_impl` for the
#' full eligibility-filtering algorithm (Stop Condition/final dimension/final variance
#' explained, each independently optional, combined by logical AND). `ierr` is set only if
#' `merge_to_super_ensembles_impl` discovers a component larger than `max_group_size` --
#' see that kernel's own doc comment.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_reconciliation::ensemble_reconciliation}, whose argument names
#' are the ones an error message reports.
#'
#' @param ensemble_masks a logical matrix. Per-ensemble membership, see Ensemble Identification's merged output
#' @param ensemble_stop_reason a integer vector. Per-ensemble Stop Condition, see `filter_ensembles_impl`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `4`.
#' @param ensemble_U_history a numeric array of rank 4. Per-ensemble trailing tangent+normal bases, see Ensemble Identification's merged
#'   output
#' @param ensemble_d_history a integer matrix. Per-ensemble trailing intrinsic dimensions
#' @param ensemble_S_history a numeric array of rank 3. Per-ensemble trailing singular values
#' @param ensemble_mu_history a numeric array of rank 3. Per-ensemble trailing centers
#' @param ensemble_G_history a numeric matrix. Per-ensemble trailing spectral gaps
#' @param ensemble_k_history a integer matrix. Per-ensemble trailing sizes
#' @param ensemble_accepted_history a logical matrix. Whether the growth iteration retained in each history column was itself accepted
#' @param mode a string, one of "report", "merge_overlap_coefficient", "merge_any". How intersections are processed
#'
#'   The default value is `"report"`.
#' @param min_overlap_coefficient a numeric scalar. Minimum Overlap Coefficient ($|\mathcal{E}_i \cap \mathcal{E}_j| /
#'   \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$) for an edge to qualify in mode
#'   \code{'merge_overlap_coefficient'};
#'   ignored in every other mode
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @param report_overlap_coefficient a logical scalar. Whether to compute and return `super_ensembles_overlap_coefficient` at all --
#'   see `merge_to_super_ensembles_impl`'s own note on this being guarded, not
#'   unconditional
#'   The default value is `FALSE`.
#' @param allowed_stop_reasons a logical vector. See `tox_shape_truthful_clustering_filter_impl`'s own
#'   `filter_ensembles_by_stop_condition_impl`
#' @param filter_dim_min a integer scalar. See `tox_shape_truthful_clustering_filter_impl`'s own
#'   `filter_ensembles_by_dimension_impl`
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param filter_dim_max a integer scalar. See `tox_shape_truthful_clustering_filter_impl`'s own
#'   `filter_ensembles_by_dimension_impl`
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param var_explained_min a numeric scalar. See `tox_shape_truthful_clustering_filter_impl`'s own
#'   `filter_ensembles_by_var_explained_impl`
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @param max_group_size a integer scalar. Maximum number of ensembles one super-ensemble (one column of `super_ensembles`)
#'   can hold; sizes its row dimension. `misc/mod_STC.md` suggests
#'   $\min(1024, N_{\mathcal{E}})$ as a sensible default -- always required, never
#'   optional with an auto-applied default here, for the same reason as
#'   `ensemble_identification`'s own `o`: a Fortran array bound cannot depend on a
#'   possibly-absent optional dummy, and a runtime-dependent value like
#'   $\min(1024, N_{\mathcal{E}})$ is not the constant expression an auto-applied
#'   default would need to be either.
#'   The minimum valid value is `2`.
#'   The maximum valid value is `n_ensembles`.
#' @return a named list with elements:
#'   \item{super_ensembles}{a integer matrix. One super-ensemble per column: the 1-indexed column indices of `ensemble_masks`
#'     belonging to that group, padded with 0 (invalid, ensembles are 1-indexed) below
#'     the group's actual size, and 0 in every row of an unused trailing column beyond
#'     `n_super_ensembles`. Sized at $N_{\mathcal{E}}(N_{\mathcal{E}}-1)$, twice mode
#'     \code{'report'}'s
#'     own true worst case ($N_{\mathcal{E}}(N_{\mathcal{E}}-1)/2$, every pair
#'     intersects) -- deliberately not divided by 2: the generator translates this
#'     specification expression close to verbatim into the Python/R bindings, where
#'     `/` on two integers is true division, not Fortran's own truncating integer
#'     division, so a literal `/2` here breaks the generated Python binding (a `float`
#'     where `np.empty`'s shape wants an `int`); see `misc/code_gen_footgun.md`. A
#'     safe, if looser, upper bound for modes 2 and 3 too, whose groups can never
#'     outnumber mode 1's own worst case.}
#'   \item{n_super_ensembles}{a integer scalar. Number of leading columns of `super_ensembles`/`super_ensembles_overlap_coefficient`
#'     actually filled}
#'   \item{super_ensembles_overlap_coefficient}{a numeric matrix. Column $l$, row $c_i$: the Overlap Coefficient between the ensembles in
#'     `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`. All zero unless
#'     `report_overlap_coefficient` was requested -- see the note above.}
#'   \item{eligible}{a logical vector. Combined per-ensemble eligibility actually used for merging above -- see
#'     `filter_ensembles_impl`. Ineligible ensembles are otherwise untouched: they
#'     are never removed from `ensemble_masks` or anything else this whole family
#'     reports, only excluded from contributing a pair here.}
#'   \item{eligible_by_stop_condition}{a logical vector. See `filter_ensembles_by_stop_condition_impl`}
#'   \item{eligible_by_dimension}{a logical vector. See `filter_ensembles_by_dimension_impl`}
#'   \item{eligible_by_var_explained}{a logical vector. See `filter_ensembles_by_var_explained_impl`}
#' @export
ensemble_reconciliation <- function(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, ensemble_accepted_history, mode = "report", min_overlap_coefficient = 0.9, report_overlap_coefficient = FALSE, allowed_stop_reasons = NULL, filter_dim_min = NULL, filter_dim_max = NULL, var_explained_min = NULL, max_group_size) {
    ensemble_masks <- .tox_as_logical(ensemble_masks, "ensemble_masks")
    ensemble_stop_reason <- .tox_as_integer_vector(ensemble_stop_reason, "ensemble_stop_reason")
    ensemble_U_history <- .tox_as_double_array(ensemble_U_history, "ensemble_U_history", 4L)
    ensemble_d_history <- .tox_as_integer_matrix(ensemble_d_history, "ensemble_d_history")
    ensemble_S_history <- .tox_as_double_array(ensemble_S_history, "ensemble_S_history", 3L)
    ensemble_mu_history <- .tox_as_double_array(ensemble_mu_history, "ensemble_mu_history", 3L)
    ensemble_G_history <- .tox_as_double_matrix(ensemble_G_history, "ensemble_G_history")
    ensemble_k_history <- .tox_as_integer_matrix(ensemble_k_history, "ensemble_k_history")
    ensemble_accepted_history <- .tox_as_logical(ensemble_accepted_history, "ensemble_accepted_history")
    mode <- .tox_as_mode(mode, "mode", c("report", "merge_overlap_coefficient", "merge_any"))
    min_overlap_coefficient <- .tox_as_double_scalar(min_overlap_coefficient, "min_overlap_coefficient")
    report_overlap_coefficient <- .tox_as_logical(report_overlap_coefficient, "report_overlap_coefficient")
    if (!is.null(allowed_stop_reasons))
        allowed_stop_reasons <- .tox_as_logical(allowed_stop_reasons, "allowed_stop_reasons")
    if (!is.null(filter_dim_min))
        filter_dim_min <- .tox_as_integer_scalar(filter_dim_min, "filter_dim_min")
    if (!is.null(filter_dim_max))
        filter_dim_max <- .tox_as_integer_scalar(filter_dim_max, "filter_dim_max")
    if (!is.null(var_explained_min))
        var_explained_min <- .tox_as_double_scalar(var_explained_min, "var_explained_min")
    max_group_size <- .tox_as_integer_scalar(max_group_size, "max_group_size")
    if (dim(ensemble_S_history)[1] != dim(ensemble_U_history)[1])
        .tox_shape_error("ensemble_S_history", dim(ensemble_S_history)[1], "ensemble_U_history", dim(ensemble_U_history)[1])
    if (dim(ensemble_mu_history)[1] != dim(ensemble_U_history)[1])
        .tox_shape_error("ensemble_mu_history", dim(ensemble_mu_history)[1], "ensemble_U_history", dim(ensemble_U_history)[1])
    if (length(ensemble_stop_reason) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_stop_reason", length(ensemble_stop_reason), "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_U_history)[4] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_U_history", dim(ensemble_U_history)[4], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_d_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_d_history", dim(ensemble_d_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_S_history)[3] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_S_history", dim(ensemble_S_history)[3], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_mu_history)[3] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_mu_history", dim(ensemble_mu_history)[3], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_G_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_G_history", dim(ensemble_G_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_k_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_k_history", dim(ensemble_k_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_accepted_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_accepted_history", dim(ensemble_accepted_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
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

    .result <- .Call("ensemble_reconciliation_call", ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, ensemble_accepted_history, mode, min_overlap_coefficient, report_overlap_coefficient, allowed_stop_reasons, filter_dim_min, filter_dim_max, var_explained_min, max_group_size)
    .arguments <- c("ensemble_masks", "ensemble_stop_reason", "n_dimensions", "n_vectors", "n_ensembles", "ensemble_U_history", "ensemble_d_history", "ensemble_S_history", "ensemble_mu_history", "ensemble_G_history", "ensemble_k_history", "ensemble_accepted_history", "o", "mode", "min_overlap_coefficient", "report_overlap_coefficient", "allowed_stop_reasons", "filter_dim_min", "filter_dim_max", "var_explained_min", "max_group_size", "super_ensembles", "n_super_ensembles", "super_ensembles_overlap_coefficient", "eligible", "eligible_by_stop_condition", "eligible_by_dimension", "eligible_by_var_explained", "ierr")
    .sources <- c(NA_character_, NA_character_, "ensemble_U_history", "ensemble_masks", "ensemble_masks", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "ensemble_U_history", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "super_ensembles", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        super_ensembles = .result$super_ensembles,
        n_super_ensembles = .result$n_super_ensembles,
        super_ensembles_overlap_coefficient = .result$super_ensembles_overlap_coefficient,
        eligible = .result$eligible,
        eligible_by_stop_condition = .result$eligible_by_stop_condition,
        eligible_by_dimension = .result$eligible_by_dimension,
        eligible_by_var_explained = .result$eligible_by_var_explained
    )
}

#' Group/report intersections among eligible ensembles into super-ensembles
#'
#' Detecting an intersection at all already requires $|\mathcal{E}_i \cap \mathcal{E}_j|$,
#' needed by every mode; the Overlap Coefficient itself is a single extra $O(1)$ step per
#' pair once each ensemble's own size is known --
#' $\text{OC}(\mathcal{E}_i, \mathcal{E}_j) = |\mathcal{E}_i \cap \mathcal{E}_j| /
#' \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$, cheaper even than the Jaccard Similarity Index
#' it replaces (no union to derive, just the smaller of the two already-precomputed sizes)
#' -- but modes 1 and 3 do not need it for their own decision, so its computation is
#' guarded behind `report_overlap_coefficient`, never unconditional (see `misc/mod_STC.md`'s
#' explicit note on this).
#'
#' An ineligible ensemble (`.not. eligible(i)`, see `filter_ensembles_impl`) never
#' contributes a pair here at all -- a plain `eligible(i) .and. eligible(j)` guard, no array
#' copying/compaction, the same "logical AND of masks" shape mode
#' \code{'merge_overlap_coefficient'}'s
#' own Overlap Coefficient threshold check already uses.
#'
#' Modes 2 and 3 group via a union-find over the qualifying-edge graph (`stc_uf_find`/
#' `stc_uf_union` below), unioning the smaller index under the larger's root so that a
#' component's root is always its own smallest member -- which is what makes the single
#' pass `do r = 1, n_ensembles` below both find every component exactly once and emit them
#' in ascending order of each group's smallest member, with no separate bookkeeping. A
#' discovered component larger than `max_group_size` is a genuine runtime condition no
#' static input check could foresee (it depends on the actual intersection pattern), so it
#' is reported via `ierr` rather than silently truncated -- see `codegen_guide.md` section
#' 5.14.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_reconciliation::merge_to_super_ensembles}, whose argument names
#' are the ones an error message reports.
#'
#' @param ensemble_masks a logical matrix. Per-ensemble membership, see Ensemble Identification's merged output
#' @param eligible a logical vector. Per-ensemble eligibility to contribute a pair here at all -- see
#'   `tox_shape_truthful_clustering_filter_impl`'s own `filter_ensembles_impl`,
#'   this kernel's own sibling in `ensemble_reconciliation`'s two-call orchestration
#' @param mode a string, one of "report", "merge_overlap_coefficient", "merge_any". How intersections are processed
#'
#'   The default value is `"report"`.
#' @param min_overlap_coefficient a numeric scalar. Minimum Overlap Coefficient ($|\mathcal{E}_i \cap \mathcal{E}_j| /
#'   \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$) for an edge to qualify in mode
#'   \code{'merge_overlap_coefficient'};
#'   ignored in every other mode
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @param report_overlap_coefficient a logical scalar. Whether to compute and return `super_ensembles_overlap_coefficient` at all --
#'   see the note above on this being guarded, not unconditional
#'   The default value is `FALSE`.
#' @param max_group_size a integer scalar. Maximum number of ensembles one super-ensemble (one column of `super_ensembles`)
#'   can hold; sizes its row dimension. `misc/mod_STC.md` suggests
#'   $\min(1024, N_{\mathcal{E}})$ as a sensible default -- always required, never
#'   optional with an auto-applied default here, for the same reason as
#'   `ensemble_identification`'s own `o`: a Fortran array bound cannot depend on a
#'   possibly-absent optional dummy, and a runtime-dependent value like
#'   $\min(1024, N_{\mathcal{E}})$ is not the constant expression an auto-applied
#'   default would need to be either.
#'   The minimum valid value is `2`.
#'   The maximum valid value is `n_ensembles`.
#' @return a named list with elements:
#'   \item{super_ensembles}{a integer matrix. One super-ensemble per column: the 1-indexed column indices of `ensemble_masks`
#'     belonging to that group, padded with 0 (invalid, ensembles are 1-indexed) below
#'     the group's actual size, and 0 in every row of an unused trailing column beyond
#'     `n_super_ensembles`. Sized at $N_{\mathcal{E}}(N_{\mathcal{E}}-1)$, twice mode
#'     \code{'report'}'s
#'     own true worst case ($N_{\mathcal{E}}(N_{\mathcal{E}}-1)/2$, every pair
#'     intersects) -- deliberately not divided by 2, see `ensemble_reconciliation_impl`'s
#'     own identical note.}
#'   \item{n_super_ensembles}{a integer scalar. Number of leading columns of `super_ensembles`/`super_ensembles_overlap_coefficient`
#'     actually filled}
#'   \item{super_ensembles_overlap_coefficient}{a numeric matrix. Column $l$, row $c_i$: the Overlap Coefficient between the ensembles in
#'     `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`. All zero unless
#'     `report_overlap_coefficient` was requested -- see the note above.}
#' @export
merge_to_super_ensembles <- function(ensemble_masks, eligible, mode = "report", min_overlap_coefficient = 0.9, report_overlap_coefficient = FALSE, max_group_size) {
    ensemble_masks <- .tox_as_logical(ensemble_masks, "ensemble_masks")
    eligible <- .tox_as_logical(eligible, "eligible")
    mode <- .tox_as_mode(mode, "mode", c("report", "merge_overlap_coefficient", "merge_any"))
    min_overlap_coefficient <- .tox_as_double_scalar(min_overlap_coefficient, "min_overlap_coefficient")
    report_overlap_coefficient <- .tox_as_logical(report_overlap_coefficient, "report_overlap_coefficient")
    max_group_size <- .tox_as_integer_scalar(max_group_size, "max_group_size")
    if (length(eligible) != dim(ensemble_masks)[2])
        .tox_shape_error("eligible", length(eligible), "ensemble_masks", dim(ensemble_masks)[2])

    .result <- .Call("merge_to_super_ensembles_call", ensemble_masks, eligible, mode, min_overlap_coefficient, report_overlap_coefficient, max_group_size)
    .arguments <- c("ensemble_masks", "eligible", "n_vectors", "n_ensembles", "mode", "min_overlap_coefficient", "report_overlap_coefficient", "max_group_size", "super_ensembles", "n_super_ensembles", "super_ensembles_overlap_coefficient", "ierr")
    .sources <- c(NA_character_, NA_character_, "ensemble_masks", "ensemble_masks", NA_character_, NA_character_, NA_character_, "super_ensembles", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        super_ensembles = .result$super_ensembles,
        n_super_ensembles = .result$n_super_ensembles,
        super_ensembles_overlap_coefficient = .result$super_ensembles_overlap_coefficient
    )
}
