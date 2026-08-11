# Generated. Do not edit.

#' Serializes each input vector's ensemble/super-ensemble membership as CSV
#'
#' `row` is the vector's own 1-based position in the input table, so this file joins back
#' to the caller's original data purely by row number. See this module's own doc comment
#' for the quoting convention.
#'
#' Generated from the Fortran procedure \code{tox_stc_csv::serialize_stc_points_as_csv}, whose argument names
#' are the ones an error message reports.
#'
#' @param filename a string. Name of the CSV file to write
#' @param n_super_ensembles a integer scalar. Number of leading columns of `super_ensembles` actually filled
#' @param seed_selection_mask a logical vector. Seed selection, see `seeds`
#' @param ensemble_masks a logical matrix. Per-ensemble accepted membership, one column per seed
#' @param ensemble_low_confidence_masks a logical matrix. Per-ensemble iteration-1 fallback membership
#' @param super_ensembles a integer matrix. One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
#' @return invisibly `NULL`; called for its effect.
#' @export
serialize_stc_points_as_csv <- function(filename, n_super_ensembles, seed_selection_mask, ensemble_masks, ensemble_low_confidence_masks, super_ensembles) {
    filename <- .tox_as_character(filename, "filename")
    n_super_ensembles <- .tox_as_integer_scalar(n_super_ensembles, "n_super_ensembles")
    seed_selection_mask <- .tox_as_logical(seed_selection_mask, "seed_selection_mask")
    ensemble_masks <- .tox_as_logical(ensemble_masks, "ensemble_masks")
    ensemble_low_confidence_masks <- .tox_as_logical(ensemble_low_confidence_masks, "ensemble_low_confidence_masks")
    super_ensembles <- .tox_as_integer_matrix(super_ensembles, "super_ensembles")
    if (dim(ensemble_masks)[1] != length(seed_selection_mask))
        .tox_shape_error("ensemble_masks", dim(ensemble_masks)[1], "seed_selection_mask", length(seed_selection_mask))
    if (dim(ensemble_low_confidence_masks)[1] != length(seed_selection_mask))
        .tox_shape_error("ensemble_low_confidence_masks", dim(ensemble_low_confidence_masks)[1], "seed_selection_mask", length(seed_selection_mask))
    if (dim(ensemble_low_confidence_masks)[2] != dim(ensemble_masks)[2])
        .tox_shape_error("ensemble_low_confidence_masks", dim(ensemble_low_confidence_masks)[2], "ensemble_masks", dim(ensemble_masks)[2])

    .result <- .Call("serialize_stc_points_as_csv_call", filename, n_super_ensembles, seed_selection_mask, ensemble_masks, ensemble_low_confidence_masks, super_ensembles)
    .arguments <- c("filename", "n_vectors", "n_selected_seed", "max_group_size", "n_super_ensembles", "seed_selection_mask", "ensemble_masks", "ensemble_low_confidence_masks", "super_ensembles", "ierr")
    .sources <- c(NA_character_, "seed_selection_mask", "ensemble_masks", "super_ensembles", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    invisible(NULL)
}

#' Serializes the full pairwise ensemble Overlap Coefficient matrix as CSV
#'
#' Only pairs with a nonempty intersection are written, matching `tox_stc_json`'s own
#' `overlap_coefficient_matrix` convention -- an absent pair means Overlap Coefficient 0.
#'
#' Generated from the Fortran procedure \code{tox_stc_csv::serialize_stc_ensemble_overlap_as_csv}, whose argument names
#' are the ones an error message reports.
#'
#' @param filename a string. Name of the CSV file to write
#' @param ensemble_masks a logical matrix. Per-ensemble accepted membership, one column per seed
#' @return invisibly `NULL`; called for its effect.
#' @export
serialize_stc_ensemble_overlap_as_csv <- function(filename, ensemble_masks) {
    filename <- .tox_as_character(filename, "filename")
    ensemble_masks <- .tox_as_logical(ensemble_masks, "ensemble_masks")
    .result <- .Call("serialize_stc_ensemble_overlap_as_csv_call", filename, ensemble_masks)
    .arguments <- c("filename", "n_vectors", "n_selected_seed", "ensemble_masks", "ierr")
    .sources <- c(NA_character_, "ensemble_masks", "ensemble_masks", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    invisible(NULL)
}

#' Serializes the super-ensembles as a gene-family-file-style TSV
#'
#' One line per super-ensemble: `<group_id>` TAB `<comma-separated member ensemble ids>`,
#' no header, no quoting -- see this module's own doc comment for why.
#'
#' Generated from the Fortran procedure \code{tox_stc_csv::serialize_stc_super_ensembles_as_tsv}, whose argument names
#' are the ones an error message reports.
#'
#' @param filename a string. Name of the TSV file to write
#' @param super_ensembles a integer matrix. One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
#' @return invisibly `NULL`; called for its effect.
#' @export
serialize_stc_super_ensembles_as_tsv <- function(filename, super_ensembles) {
    filename <- .tox_as_character(filename, "filename")
    super_ensembles <- .tox_as_integer_matrix(super_ensembles, "super_ensembles")
    .result <- .Call("serialize_stc_super_ensembles_as_tsv_call", filename, super_ensembles)
    .arguments <- c("filename", "max_group_size", "n_super_ensembles", "super_ensembles", "ierr")
    .sources <- c(NA_character_, "super_ensembles", "super_ensembles", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    invisible(NULL)
}
