# Generated. Do not edit.

#' Serializes an STC run's raw pipeline results as JSON
#'
#' Matches the schema consumed by misc/STC-experiments/interactive_template.html's D3 report
#'
#' Generated from the Fortran procedure \code{tox_stc_json::serialize_stc_results_as_json}, whose argument names
#' are the ones an error message reports.
#'
#' @param filename a string. Name of the JSON file to write
#' @param n_super_ensembles a integer scalar. Number of leading columns of `super_ensembles` actually filled
#' @param vectors a numeric matrix. Input data matrix
#' @param dim_names a character vector. Per-dimension display name
#' @param seed_selection_mask a logical vector. Seed selection, see `seeds`
#' @param ensemble_masks a logical matrix. Per-ensemble accepted membership, one column per seed
#' @param ensemble_stop_reason a integer vector. Per-ensemble Stop Condition
#' @param ensemble_growth_radii a numeric vector. Per-ensemble growth radius
#' @param ensemble_U_history a numeric array of rank 4. Per-ensemble trailing tangent+normal bases
#' @param ensemble_S_history a numeric array of rank 3. Per-ensemble trailing singular values
#' @param ensemble_d_history a integer matrix. Per-ensemble trailing intrinsic dimensions
#' @param ensemble_G_history a numeric matrix. Per-ensemble trailing spectral gaps
#' @param ensemble_mu_history a numeric array of rank 3. Per-ensemble trailing centers
#' @param ensemble_k_history a integer matrix. Per-ensemble trailing sizes
#' @param ensemble_accepted_history a logical matrix. Whether the growth iteration retained in each history column was itself accepted
#'   -- `stc_push_ensemble_history` also pushes a *rejected* final candidate before
#'   `ensemble_identification` halts growth via `STOP_REASON_REJECTED_IMMEDIATELY`/
#'   `STOP_REASON_REJECTED_AFTER_STABLE`, so the last populated column is not always
#'   the ensemble's actual last accepted state; this module uses this array to find
#'   the last column that genuinely is (see `stc_last_accepted_history_index`)
#' @param ensemble_member_added_at_step a integer matrix. Per-ensemble growth-iteration-joined bookkeeping, see `ensemble_identification`'s
#'   `member_added_at_step`; this module only ever reads its column max (= T, the
#'   final accepted growth iteration), not the per-vector values themselves
#' @param ensemble_low_confidence_masks a logical matrix. Per-ensemble iteration-1 fallback membership
#' @param ensemble_U_first a numeric array of rank 3. Per-ensemble tangent+normal basis at the bootstrap iteration (iteration 1)
#' @param ensemble_d_first a integer vector. Per-ensemble intrinsic dimension at the bootstrap iteration
#' @param super_ensembles a integer matrix. One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
#' @param k_min a integer scalar. This run's neighborhood size for each seed's growth radius
#' @param k_density a integer scalar. This run's density estimation neighborhood size
#' @param chordal_dist_max_as_prcnt_of_range a numeric scalar. This run's maximum tolerated chordal distance between tangent bases
#' @param d_max a integer scalar. This run's maximum tolerated change in intrinsic dimension
#' @param G_max a numeric scalar. This run's maximum tolerated |log(G_tp1/G_t)|
#' @param RMSE_change_max a numeric scalar. This run's maximum tolerated |log(RMSE_tp1/RMSE_t)|
#' @param f_max a numeric scalar. This run's ensemble size fraction of N above which growth is abandoned
#' @param a a integer scalar. This run's minimum accepted-iteration count for a stable rejection
#' @param exclusion_radius_percentile a numeric scalar. This run's seeding exclusion radius percentile
#' @param bandwidth_percentile a numeric scalar. This run's density-estimate kernel bandwidth percentile
#' @param reconciliation_mode a string, one of "report", "merge_overlap_coefficient", "merge_any". This run's `ensemble_reconciliation` mode
#' @param min_overlap_coefficient a numeric scalar. This run's minimum Overlap Coefficient for `MODE_MERGE_OVERLAP_COEFFICIENT`
#' @param allowed_stop_reasons a logical vector. This run's per-Stop-Condition eligibility actually used by
#'   `ensemble_reconciliation` -- reported here (as `params.excluded_stop_reasons`)
#'   for transparency only; this module no longer derives eligibility from it itself,
#'   see `ensemble_eligible` below
#' @param filter_d_min a integer scalar. This run's minimum tolerated final intrinsic dimension for reconciliation
#'   eligibility, inclusive -- see `tox_shape_truthful_clustering_filter_impl`'s own
#'   `d_min`; reported for transparency only, same as `allowed_stop_reasons` above
#' @param filter_d_max a integer scalar. This run's maximum tolerated final intrinsic dimension for reconciliation
#'   eligibility, inclusive -- see `tox_shape_truthful_clustering_filter_impl`'s own
#'   `d_max`; reported for transparency only, same as `allowed_stop_reasons` above
#' @param filter_var_explained_min a numeric scalar. This run's minimum tolerated final variance explained for reconciliation
#'   eligibility -- see `tox_shape_truthful_clustering_filter_impl`'s own
#'   `var_explained_min`; reported for transparency only, same as
#'   `allowed_stop_reasons` above
#' @param ensemble_eligible a logical vector. Per-ensemble combined reconciliation eligibility actually used by
#'   `ensemble_reconciliation`, see its own `eligible` output
#' @param ensemble_eligible_by_stop_condition a logical vector. See `ensemble_reconciliation`'s own `eligible_by_stop_condition`
#' @param ensemble_eligible_by_dimension a logical vector. See `ensemble_reconciliation`'s own `eligible_by_dimension`
#' @param ensemble_eligible_by_var_explained a logical vector. See `ensemble_reconciliation`'s own `eligible_by_var_explained`
#' @param estimated_k_min a integer scalar. `estimate_stc_parameters`'s proposed `k_min`, if estimation was used
#' @param estimated_k_density a integer scalar. `estimate_stc_parameters`'s proposed `k_density`, if estimation was used
#' @param estimated_density_quantile a numeric scalar. `estimate_stc_parameters`'s proposed density quantile, if estimation was used
#' @param estimated_chordal_dist_max_as_prcnt_of_range a numeric scalar. `estimate_stc_parameters`'s proposed `chordal_dist_max_as_prcnt_of_range`, if
#'   estimation was used
#' @param estimated_G_max a numeric scalar. `estimate_stc_parameters`'s proposed `G_max`, if estimation was used
#' @param estimated_d_max a integer scalar. `estimate_stc_parameters`'s proposed `d_max`, if estimation was used
#' @return invisibly `NULL`; called for its effect.
#' @export
serialize_stc_results_as_json <- function(filename, n_super_ensembles, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, k_min, k_density, chordal_dist_max_as_prcnt_of_range, d_max, G_max, RMSE_change_max, f_max, a, exclusion_radius_percentile, bandwidth_percentile, reconciliation_mode, min_overlap_coefficient, allowed_stop_reasons = NULL, filter_d_min = NULL, filter_d_max = NULL, filter_var_explained_min = NULL, ensemble_eligible, ensemble_eligible_by_stop_condition, ensemble_eligible_by_dimension, ensemble_eligible_by_var_explained, estimated_k_min = NULL, estimated_k_density = NULL, estimated_density_quantile = NULL, estimated_chordal_dist_max_as_prcnt_of_range = NULL, estimated_G_max = NULL, estimated_d_max = NULL) {
    filename <- .tox_as_character(filename, "filename")
    n_super_ensembles <- .tox_as_integer_scalar(n_super_ensembles, "n_super_ensembles")
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    dim_names <- .tox_as_character(dim_names, "dim_names")
    seed_selection_mask <- .tox_as_logical(seed_selection_mask, "seed_selection_mask")
    ensemble_masks <- .tox_as_logical(ensemble_masks, "ensemble_masks")
    ensemble_stop_reason <- .tox_as_integer_vector(ensemble_stop_reason, "ensemble_stop_reason")
    ensemble_growth_radii <- .tox_as_double_vector(ensemble_growth_radii, "ensemble_growth_radii")
    ensemble_U_history <- .tox_as_double_array(ensemble_U_history, "ensemble_U_history", 4L)
    ensemble_S_history <- .tox_as_double_array(ensemble_S_history, "ensemble_S_history", 3L)
    ensemble_d_history <- .tox_as_integer_matrix(ensemble_d_history, "ensemble_d_history")
    ensemble_G_history <- .tox_as_double_matrix(ensemble_G_history, "ensemble_G_history")
    ensemble_mu_history <- .tox_as_double_array(ensemble_mu_history, "ensemble_mu_history", 3L)
    ensemble_k_history <- .tox_as_integer_matrix(ensemble_k_history, "ensemble_k_history")
    ensemble_accepted_history <- .tox_as_logical(ensemble_accepted_history, "ensemble_accepted_history")
    ensemble_member_added_at_step <- .tox_as_integer_matrix(ensemble_member_added_at_step, "ensemble_member_added_at_step")
    ensemble_low_confidence_masks <- .tox_as_logical(ensemble_low_confidence_masks, "ensemble_low_confidence_masks")
    ensemble_U_first <- .tox_as_double_array(ensemble_U_first, "ensemble_U_first", 3L)
    ensemble_d_first <- .tox_as_integer_vector(ensemble_d_first, "ensemble_d_first")
    super_ensembles <- .tox_as_integer_matrix(super_ensembles, "super_ensembles")
    k_min <- .tox_as_integer_scalar(k_min, "k_min")
    k_density <- .tox_as_integer_scalar(k_density, "k_density")
    chordal_dist_max_as_prcnt_of_range <- .tox_as_double_scalar(chordal_dist_max_as_prcnt_of_range, "chordal_dist_max_as_prcnt_of_range")
    d_max <- .tox_as_integer_scalar(d_max, "d_max")
    G_max <- .tox_as_double_scalar(G_max, "G_max")
    RMSE_change_max <- .tox_as_double_scalar(RMSE_change_max, "RMSE_change_max")
    f_max <- .tox_as_double_scalar(f_max, "f_max")
    a <- .tox_as_integer_scalar(a, "a")
    exclusion_radius_percentile <- .tox_as_double_scalar(exclusion_radius_percentile, "exclusion_radius_percentile")
    bandwidth_percentile <- .tox_as_double_scalar(bandwidth_percentile, "bandwidth_percentile")
    reconciliation_mode <- .tox_as_mode(reconciliation_mode, "reconciliation_mode", c("report", "merge_overlap_coefficient", "merge_any"))
    min_overlap_coefficient <- .tox_as_double_scalar(min_overlap_coefficient, "min_overlap_coefficient")
    if (!is.null(allowed_stop_reasons))
        allowed_stop_reasons <- .tox_as_logical(allowed_stop_reasons, "allowed_stop_reasons")
    if (!is.null(filter_d_min))
        filter_d_min <- .tox_as_integer_scalar(filter_d_min, "filter_d_min")
    if (!is.null(filter_d_max))
        filter_d_max <- .tox_as_integer_scalar(filter_d_max, "filter_d_max")
    if (!is.null(filter_var_explained_min))
        filter_var_explained_min <- .tox_as_double_scalar(filter_var_explained_min, "filter_var_explained_min")
    ensemble_eligible <- .tox_as_logical(ensemble_eligible, "ensemble_eligible")
    ensemble_eligible_by_stop_condition <- .tox_as_logical(ensemble_eligible_by_stop_condition, "ensemble_eligible_by_stop_condition")
    ensemble_eligible_by_dimension <- .tox_as_logical(ensemble_eligible_by_dimension, "ensemble_eligible_by_dimension")
    ensemble_eligible_by_var_explained <- .tox_as_logical(ensemble_eligible_by_var_explained, "ensemble_eligible_by_var_explained")
    if (!is.null(estimated_k_min))
        estimated_k_min <- .tox_as_integer_scalar(estimated_k_min, "estimated_k_min")
    if (!is.null(estimated_k_density))
        estimated_k_density <- .tox_as_integer_scalar(estimated_k_density, "estimated_k_density")
    if (!is.null(estimated_density_quantile))
        estimated_density_quantile <- .tox_as_double_scalar(estimated_density_quantile, "estimated_density_quantile")
    if (!is.null(estimated_chordal_dist_max_as_prcnt_of_range))
        estimated_chordal_dist_max_as_prcnt_of_range <- .tox_as_double_scalar(estimated_chordal_dist_max_as_prcnt_of_range, "estimated_chordal_dist_max_as_prcnt_of_range")
    if (!is.null(estimated_G_max))
        estimated_G_max <- .tox_as_double_scalar(estimated_G_max, "estimated_G_max")
    if (!is.null(estimated_d_max))
        estimated_d_max <- .tox_as_integer_scalar(estimated_d_max, "estimated_d_max")
    if (length(dim_names) != dim(vectors)[1])
        .tox_shape_error("dim_names", length(dim_names), "vectors", dim(vectors)[1])
    if (dim(ensemble_U_history)[1] != dim(vectors)[1])
        .tox_shape_error("ensemble_U_history", dim(ensemble_U_history)[1], "vectors", dim(vectors)[1])
    if (dim(ensemble_S_history)[1] != dim(vectors)[1])
        .tox_shape_error("ensemble_S_history", dim(ensemble_S_history)[1], "vectors", dim(vectors)[1])
    if (dim(ensemble_mu_history)[1] != dim(vectors)[1])
        .tox_shape_error("ensemble_mu_history", dim(ensemble_mu_history)[1], "vectors", dim(vectors)[1])
    if (dim(ensemble_U_first)[1] != dim(vectors)[1])
        .tox_shape_error("ensemble_U_first", dim(ensemble_U_first)[1], "vectors", dim(vectors)[1])
    if (length(seed_selection_mask) != dim(vectors)[2])
        .tox_shape_error("seed_selection_mask", length(seed_selection_mask), "vectors", dim(vectors)[2])
    if (dim(ensemble_masks)[1] != dim(vectors)[2])
        .tox_shape_error("ensemble_masks", dim(ensemble_masks)[1], "vectors", dim(vectors)[2])
    if (dim(ensemble_member_added_at_step)[1] != dim(vectors)[2])
        .tox_shape_error("ensemble_member_added_at_step", dim(ensemble_member_added_at_step)[1], "vectors", dim(vectors)[2])
    if (dim(ensemble_low_confidence_masks)[1] != dim(vectors)[2])
        .tox_shape_error("ensemble_low_confidence_masks", dim(ensemble_low_confidence_masks)[1], "vectors", dim(vectors)[2])
    if (length(ensemble_stop_reason) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_stop_reason", length(ensemble_stop_reason), "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_growth_radii) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_growth_radii", length(ensemble_growth_radii), "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_U_history)[4] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_U_history", dim(ensemble_U_history)[4], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_S_history)[3] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_S_history", dim(ensemble_S_history)[3], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_d_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_d_history", dim(ensemble_d_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_G_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_G_history", dim(ensemble_G_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_mu_history)[3] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_mu_history", dim(ensemble_mu_history)[3], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_k_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_k_history", dim(ensemble_k_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_accepted_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_accepted_history", dim(ensemble_accepted_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_member_added_at_step)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_member_added_at_step", dim(ensemble_member_added_at_step)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_low_confidence_masks)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_low_confidence_masks", dim(ensemble_low_confidence_masks)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_U_first)[3] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_U_first", dim(ensemble_U_first)[3], "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_d_first) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_d_first", length(ensemble_d_first), "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_eligible) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_eligible", length(ensemble_eligible), "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_eligible_by_stop_condition) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_eligible_by_stop_condition", length(ensemble_eligible_by_stop_condition), "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_eligible_by_dimension) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_eligible_by_dimension", length(ensemble_eligible_by_dimension), "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_eligible_by_var_explained) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_eligible_by_var_explained", length(ensemble_eligible_by_var_explained), "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_S_history)[2] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_S_history", dim(ensemble_S_history)[2], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_d_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_d_history", dim(ensemble_d_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_G_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_G_history", dim(ensemble_G_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_mu_history)[2] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_mu_history", dim(ensemble_mu_history)[2], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_k_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_k_history", dim(ensemble_k_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_accepted_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_accepted_history", dim(ensemble_accepted_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])

    .result <- .Call("serialize_stc_results_as_json_call", filename, n_super_ensembles, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, k_min, k_density, chordal_dist_max_as_prcnt_of_range, d_max, G_max, RMSE_change_max, f_max, a, exclusion_radius_percentile, bandwidth_percentile, reconciliation_mode, min_overlap_coefficient, allowed_stop_reasons, filter_d_min, filter_d_max, filter_var_explained_min, ensemble_eligible, ensemble_eligible_by_stop_condition, ensemble_eligible_by_dimension, ensemble_eligible_by_var_explained, estimated_k_min, estimated_k_density, estimated_density_quantile, estimated_chordal_dist_max_as_prcnt_of_range, estimated_G_max, estimated_d_max)
    .arguments <- c("filename", "n_dimensions", "n_vectors", "n_selected_seed", "o", "max_group_size", "n_super_ensembles", "vectors", "dim_names", "seed_selection_mask", "ensemble_masks", "ensemble_stop_reason", "ensemble_growth_radii", "ensemble_U_history", "ensemble_S_history", "ensemble_d_history", "ensemble_G_history", "ensemble_mu_history", "ensemble_k_history", "ensemble_accepted_history", "ensemble_member_added_at_step", "ensemble_low_confidence_masks", "ensemble_U_first", "ensemble_d_first", "super_ensembles", "k_min", "k_density", "chordal_dist_max_as_prcnt_of_range", "d_max", "G_max", "RMSE_change_max", "f_max", "a", "exclusion_radius_percentile", "bandwidth_percentile", "reconciliation_mode", "min_overlap_coefficient", "allowed_stop_reasons", "filter_d_min", "filter_d_max", "filter_var_explained_min", "ensemble_eligible", "ensemble_eligible_by_stop_condition", "ensemble_eligible_by_dimension", "ensemble_eligible_by_var_explained", "estimated_k_min", "estimated_k_density", "estimated_density_quantile", "estimated_chordal_dist_max_as_prcnt_of_range", "estimated_G_max", "estimated_d_max", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", "ensemble_masks", "ensemble_U_history", "super_ensembles", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    invisible(NULL)
}

#' Writes a self-contained interactive HTML report for an STC run
#'
#' Concatenates the vendored D3 library, the D3 report template (both baked in at compile
#' time, see `tox_stc_html_assets`), and this run's own results as JSON into one file
#'
#' Generated from the Fortran procedure \code{tox_stc_json::write_stc_interactive_html_report}, whose argument names
#' are the ones an error message reports.
#'
#' @param filename a string. Name of the HTML file to write
#' @param n_super_ensembles a integer scalar. Number of leading columns of `super_ensembles` actually filled
#' @param vectors a numeric matrix. Input data matrix
#' @param dim_names a character vector. Per-dimension display name
#' @param seed_selection_mask a logical vector. Seed selection, see `seeds`
#' @param ensemble_masks a logical matrix. Per-ensemble accepted membership, one column per seed
#' @param ensemble_stop_reason a integer vector. Per-ensemble Stop Condition
#' @param ensemble_growth_radii a numeric vector. Per-ensemble growth radius
#' @param ensemble_U_history a numeric array of rank 4. Per-ensemble trailing tangent+normal bases
#' @param ensemble_S_history a numeric array of rank 3. Per-ensemble trailing singular values
#' @param ensemble_d_history a integer matrix. Per-ensemble trailing intrinsic dimensions
#' @param ensemble_G_history a numeric matrix. Per-ensemble trailing spectral gaps
#' @param ensemble_mu_history a numeric array of rank 3. Per-ensemble trailing centers
#' @param ensemble_k_history a integer matrix. Per-ensemble trailing sizes
#' @param ensemble_accepted_history a logical matrix. Whether the growth iteration retained in each history column was itself accepted
#'   -- `stc_push_ensemble_history` also pushes a *rejected* final candidate before
#'   `ensemble_identification` halts growth via `STOP_REASON_REJECTED_IMMEDIATELY`/
#'   `STOP_REASON_REJECTED_AFTER_STABLE`, so the last populated column is not always
#'   the ensemble's actual last accepted state; this module uses this array to find
#'   the last column that genuinely is (see `stc_last_accepted_history_index`)
#' @param ensemble_member_added_at_step a integer matrix. Per-ensemble growth-iteration-joined bookkeeping, see `ensemble_identification`'s
#'   `member_added_at_step`; this module only ever reads its column max (= T, the
#'   final accepted growth iteration), not the per-vector values themselves
#' @param ensemble_low_confidence_masks a logical matrix. Per-ensemble iteration-1 fallback membership
#' @param ensemble_U_first a numeric array of rank 3. Per-ensemble tangent+normal basis at the bootstrap iteration (iteration 1)
#' @param ensemble_d_first a integer vector. Per-ensemble intrinsic dimension at the bootstrap iteration
#' @param super_ensembles a integer matrix. One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
#' @param k_min a integer scalar. This run's neighborhood size for each seed's growth radius
#' @param k_density a integer scalar. This run's density estimation neighborhood size
#' @param chordal_dist_max_as_prcnt_of_range a numeric scalar. This run's maximum tolerated chordal distance between tangent bases
#' @param d_max a integer scalar. This run's maximum tolerated change in intrinsic dimension
#' @param G_max a numeric scalar. This run's maximum tolerated |log(G_tp1/G_t)|
#' @param RMSE_change_max a numeric scalar. This run's maximum tolerated |log(RMSE_tp1/RMSE_t)|
#' @param f_max a numeric scalar. This run's ensemble size fraction of N above which growth is abandoned
#' @param a a integer scalar. This run's minimum accepted-iteration count for a stable rejection
#' @param exclusion_radius_percentile a numeric scalar. This run's seeding exclusion radius percentile
#' @param bandwidth_percentile a numeric scalar. This run's density-estimate kernel bandwidth percentile
#' @param reconciliation_mode a string, one of "report", "merge_overlap_coefficient", "merge_any". This run's `ensemble_reconciliation` mode
#' @param min_overlap_coefficient a numeric scalar. This run's minimum Overlap Coefficient for `MODE_MERGE_OVERLAP_COEFFICIENT`
#' @param allowed_stop_reasons a logical vector. This run's per-Stop-Condition eligibility actually used by
#'   `ensemble_reconciliation` -- reported here (as `params.excluded_stop_reasons`)
#'   for transparency only; this module no longer derives eligibility from it itself,
#'   see `ensemble_eligible` below
#' @param filter_d_min a integer scalar. This run's minimum tolerated final intrinsic dimension for reconciliation
#'   eligibility, inclusive -- see `tox_shape_truthful_clustering_filter_impl`'s own
#'   `d_min`; reported for transparency only, same as `allowed_stop_reasons` above
#' @param filter_d_max a integer scalar. This run's maximum tolerated final intrinsic dimension for reconciliation
#'   eligibility, inclusive -- see `tox_shape_truthful_clustering_filter_impl`'s own
#'   `d_max`; reported for transparency only, same as `allowed_stop_reasons` above
#' @param filter_var_explained_min a numeric scalar. This run's minimum tolerated final variance explained for reconciliation
#'   eligibility -- see `tox_shape_truthful_clustering_filter_impl`'s own
#'   `var_explained_min`; reported for transparency only, same as
#'   `allowed_stop_reasons` above
#' @param ensemble_eligible a logical vector. Per-ensemble combined reconciliation eligibility actually used by
#'   `ensemble_reconciliation`, see its own `eligible` output
#' @param ensemble_eligible_by_stop_condition a logical vector. See `ensemble_reconciliation`'s own `eligible_by_stop_condition`
#' @param ensemble_eligible_by_dimension a logical vector. See `ensemble_reconciliation`'s own `eligible_by_dimension`
#' @param ensemble_eligible_by_var_explained a logical vector. See `ensemble_reconciliation`'s own `eligible_by_var_explained`
#' @param estimated_k_min a integer scalar. `estimate_stc_parameters`'s proposed `k_min`, if estimation was used
#' @param estimated_k_density a integer scalar. `estimate_stc_parameters`'s proposed `k_density`, if estimation was used
#' @param estimated_density_quantile a numeric scalar. `estimate_stc_parameters`'s proposed density quantile, if estimation was used
#' @param estimated_chordal_dist_max_as_prcnt_of_range a numeric scalar. `estimate_stc_parameters`'s proposed `chordal_dist_max_as_prcnt_of_range`, if
#'   estimation was used
#' @param estimated_G_max a numeric scalar. `estimate_stc_parameters`'s proposed `G_max`, if estimation was used
#' @param estimated_d_max a integer scalar. `estimate_stc_parameters`'s proposed `d_max`, if estimation was used
#' @return invisibly `NULL`; called for its effect.
#' @export
write_stc_interactive_html_report <- function(filename, n_super_ensembles, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, k_min, k_density, chordal_dist_max_as_prcnt_of_range, d_max, G_max, RMSE_change_max, f_max, a, exclusion_radius_percentile, bandwidth_percentile, reconciliation_mode, min_overlap_coefficient, allowed_stop_reasons = NULL, filter_d_min = NULL, filter_d_max = NULL, filter_var_explained_min = NULL, ensemble_eligible, ensemble_eligible_by_stop_condition, ensemble_eligible_by_dimension, ensemble_eligible_by_var_explained, estimated_k_min = NULL, estimated_k_density = NULL, estimated_density_quantile = NULL, estimated_chordal_dist_max_as_prcnt_of_range = NULL, estimated_G_max = NULL, estimated_d_max = NULL) {
    filename <- .tox_as_character(filename, "filename")
    n_super_ensembles <- .tox_as_integer_scalar(n_super_ensembles, "n_super_ensembles")
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    dim_names <- .tox_as_character(dim_names, "dim_names")
    seed_selection_mask <- .tox_as_logical(seed_selection_mask, "seed_selection_mask")
    ensemble_masks <- .tox_as_logical(ensemble_masks, "ensemble_masks")
    ensemble_stop_reason <- .tox_as_integer_vector(ensemble_stop_reason, "ensemble_stop_reason")
    ensemble_growth_radii <- .tox_as_double_vector(ensemble_growth_radii, "ensemble_growth_radii")
    ensemble_U_history <- .tox_as_double_array(ensemble_U_history, "ensemble_U_history", 4L)
    ensemble_S_history <- .tox_as_double_array(ensemble_S_history, "ensemble_S_history", 3L)
    ensemble_d_history <- .tox_as_integer_matrix(ensemble_d_history, "ensemble_d_history")
    ensemble_G_history <- .tox_as_double_matrix(ensemble_G_history, "ensemble_G_history")
    ensemble_mu_history <- .tox_as_double_array(ensemble_mu_history, "ensemble_mu_history", 3L)
    ensemble_k_history <- .tox_as_integer_matrix(ensemble_k_history, "ensemble_k_history")
    ensemble_accepted_history <- .tox_as_logical(ensemble_accepted_history, "ensemble_accepted_history")
    ensemble_member_added_at_step <- .tox_as_integer_matrix(ensemble_member_added_at_step, "ensemble_member_added_at_step")
    ensemble_low_confidence_masks <- .tox_as_logical(ensemble_low_confidence_masks, "ensemble_low_confidence_masks")
    ensemble_U_first <- .tox_as_double_array(ensemble_U_first, "ensemble_U_first", 3L)
    ensemble_d_first <- .tox_as_integer_vector(ensemble_d_first, "ensemble_d_first")
    super_ensembles <- .tox_as_integer_matrix(super_ensembles, "super_ensembles")
    k_min <- .tox_as_integer_scalar(k_min, "k_min")
    k_density <- .tox_as_integer_scalar(k_density, "k_density")
    chordal_dist_max_as_prcnt_of_range <- .tox_as_double_scalar(chordal_dist_max_as_prcnt_of_range, "chordal_dist_max_as_prcnt_of_range")
    d_max <- .tox_as_integer_scalar(d_max, "d_max")
    G_max <- .tox_as_double_scalar(G_max, "G_max")
    RMSE_change_max <- .tox_as_double_scalar(RMSE_change_max, "RMSE_change_max")
    f_max <- .tox_as_double_scalar(f_max, "f_max")
    a <- .tox_as_integer_scalar(a, "a")
    exclusion_radius_percentile <- .tox_as_double_scalar(exclusion_radius_percentile, "exclusion_radius_percentile")
    bandwidth_percentile <- .tox_as_double_scalar(bandwidth_percentile, "bandwidth_percentile")
    reconciliation_mode <- .tox_as_mode(reconciliation_mode, "reconciliation_mode", c("report", "merge_overlap_coefficient", "merge_any"))
    min_overlap_coefficient <- .tox_as_double_scalar(min_overlap_coefficient, "min_overlap_coefficient")
    if (!is.null(allowed_stop_reasons))
        allowed_stop_reasons <- .tox_as_logical(allowed_stop_reasons, "allowed_stop_reasons")
    if (!is.null(filter_d_min))
        filter_d_min <- .tox_as_integer_scalar(filter_d_min, "filter_d_min")
    if (!is.null(filter_d_max))
        filter_d_max <- .tox_as_integer_scalar(filter_d_max, "filter_d_max")
    if (!is.null(filter_var_explained_min))
        filter_var_explained_min <- .tox_as_double_scalar(filter_var_explained_min, "filter_var_explained_min")
    ensemble_eligible <- .tox_as_logical(ensemble_eligible, "ensemble_eligible")
    ensemble_eligible_by_stop_condition <- .tox_as_logical(ensemble_eligible_by_stop_condition, "ensemble_eligible_by_stop_condition")
    ensemble_eligible_by_dimension <- .tox_as_logical(ensemble_eligible_by_dimension, "ensemble_eligible_by_dimension")
    ensemble_eligible_by_var_explained <- .tox_as_logical(ensemble_eligible_by_var_explained, "ensemble_eligible_by_var_explained")
    if (!is.null(estimated_k_min))
        estimated_k_min <- .tox_as_integer_scalar(estimated_k_min, "estimated_k_min")
    if (!is.null(estimated_k_density))
        estimated_k_density <- .tox_as_integer_scalar(estimated_k_density, "estimated_k_density")
    if (!is.null(estimated_density_quantile))
        estimated_density_quantile <- .tox_as_double_scalar(estimated_density_quantile, "estimated_density_quantile")
    if (!is.null(estimated_chordal_dist_max_as_prcnt_of_range))
        estimated_chordal_dist_max_as_prcnt_of_range <- .tox_as_double_scalar(estimated_chordal_dist_max_as_prcnt_of_range, "estimated_chordal_dist_max_as_prcnt_of_range")
    if (!is.null(estimated_G_max))
        estimated_G_max <- .tox_as_double_scalar(estimated_G_max, "estimated_G_max")
    if (!is.null(estimated_d_max))
        estimated_d_max <- .tox_as_integer_scalar(estimated_d_max, "estimated_d_max")
    if (length(dim_names) != dim(vectors)[1])
        .tox_shape_error("dim_names", length(dim_names), "vectors", dim(vectors)[1])
    if (dim(ensemble_U_history)[1] != dim(vectors)[1])
        .tox_shape_error("ensemble_U_history", dim(ensemble_U_history)[1], "vectors", dim(vectors)[1])
    if (dim(ensemble_S_history)[1] != dim(vectors)[1])
        .tox_shape_error("ensemble_S_history", dim(ensemble_S_history)[1], "vectors", dim(vectors)[1])
    if (dim(ensemble_mu_history)[1] != dim(vectors)[1])
        .tox_shape_error("ensemble_mu_history", dim(ensemble_mu_history)[1], "vectors", dim(vectors)[1])
    if (dim(ensemble_U_first)[1] != dim(vectors)[1])
        .tox_shape_error("ensemble_U_first", dim(ensemble_U_first)[1], "vectors", dim(vectors)[1])
    if (length(seed_selection_mask) != dim(vectors)[2])
        .tox_shape_error("seed_selection_mask", length(seed_selection_mask), "vectors", dim(vectors)[2])
    if (dim(ensemble_masks)[1] != dim(vectors)[2])
        .tox_shape_error("ensemble_masks", dim(ensemble_masks)[1], "vectors", dim(vectors)[2])
    if (dim(ensemble_member_added_at_step)[1] != dim(vectors)[2])
        .tox_shape_error("ensemble_member_added_at_step", dim(ensemble_member_added_at_step)[1], "vectors", dim(vectors)[2])
    if (dim(ensemble_low_confidence_masks)[1] != dim(vectors)[2])
        .tox_shape_error("ensemble_low_confidence_masks", dim(ensemble_low_confidence_masks)[1], "vectors", dim(vectors)[2])
    if (length(ensemble_stop_reason) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_stop_reason", length(ensemble_stop_reason), "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_growth_radii) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_growth_radii", length(ensemble_growth_radii), "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_U_history)[4] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_U_history", dim(ensemble_U_history)[4], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_S_history)[3] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_S_history", dim(ensemble_S_history)[3], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_d_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_d_history", dim(ensemble_d_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_G_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_G_history", dim(ensemble_G_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_mu_history)[3] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_mu_history", dim(ensemble_mu_history)[3], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_k_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_k_history", dim(ensemble_k_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_accepted_history)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_accepted_history", dim(ensemble_accepted_history)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_member_added_at_step)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_member_added_at_step", dim(ensemble_member_added_at_step)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_low_confidence_masks)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_low_confidence_masks", dim(ensemble_low_confidence_masks)[2], "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_U_first)[3] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_U_first", dim(ensemble_U_first)[3], "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_d_first) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_d_first", length(ensemble_d_first), "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_eligible) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_eligible", length(ensemble_eligible), "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_eligible_by_stop_condition) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_eligible_by_stop_condition", length(ensemble_eligible_by_stop_condition), "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_eligible_by_dimension) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_eligible_by_dimension", length(ensemble_eligible_by_dimension), "ensemble_masks", dim(ensemble_masks)[2])
    if (length(ensemble_eligible_by_var_explained) != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_eligible_by_var_explained", length(ensemble_eligible_by_var_explained), "ensemble_masks", dim(ensemble_masks)[2])
    if (dim(ensemble_S_history)[2] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_S_history", dim(ensemble_S_history)[2], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_d_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_d_history", dim(ensemble_d_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_G_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_G_history", dim(ensemble_G_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_mu_history)[2] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_mu_history", dim(ensemble_mu_history)[2], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_k_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_k_history", dim(ensemble_k_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])
    if (dim(ensemble_accepted_history)[1] != dim(ensemble_U_history)[3])
        .tox_shape_error("ensemble_accepted_history", dim(ensemble_accepted_history)[1], "ensemble_U_history", dim(ensemble_U_history)[3])

    .result <- .Call("write_stc_interactive_html_report_call", filename, n_super_ensembles, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, k_min, k_density, chordal_dist_max_as_prcnt_of_range, d_max, G_max, RMSE_change_max, f_max, a, exclusion_radius_percentile, bandwidth_percentile, reconciliation_mode, min_overlap_coefficient, allowed_stop_reasons, filter_d_min, filter_d_max, filter_var_explained_min, ensemble_eligible, ensemble_eligible_by_stop_condition, ensemble_eligible_by_dimension, ensemble_eligible_by_var_explained, estimated_k_min, estimated_k_density, estimated_density_quantile, estimated_chordal_dist_max_as_prcnt_of_range, estimated_G_max, estimated_d_max)
    .arguments <- c("filename", "n_dimensions", "n_vectors", "n_selected_seed", "o", "max_group_size", "n_super_ensembles", "vectors", "dim_names", "seed_selection_mask", "ensemble_masks", "ensemble_stop_reason", "ensemble_growth_radii", "ensemble_U_history", "ensemble_S_history", "ensemble_d_history", "ensemble_G_history", "ensemble_mu_history", "ensemble_k_history", "ensemble_accepted_history", "ensemble_member_added_at_step", "ensemble_low_confidence_masks", "ensemble_U_first", "ensemble_d_first", "super_ensembles", "k_min", "k_density", "chordal_dist_max_as_prcnt_of_range", "d_max", "G_max", "RMSE_change_max", "f_max", "a", "exclusion_radius_percentile", "bandwidth_percentile", "reconciliation_mode", "min_overlap_coefficient", "allowed_stop_reasons", "filter_d_min", "filter_d_max", "filter_var_explained_min", "ensemble_eligible", "ensemble_eligible_by_stop_condition", "ensemble_eligible_by_dimension", "ensemble_eligible_by_var_explained", "estimated_k_min", "estimated_k_density", "estimated_density_quantile", "estimated_chordal_dist_max_as_prcnt_of_range", "estimated_G_max", "estimated_d_max", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", "ensemble_masks", "ensemble_U_history", "super_ensembles", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    invisible(NULL)
}
