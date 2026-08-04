# Generated. Do not edit.

#' Computes normalized tissue versatility for selected expression vectors.
#'
#' The metric is based on the angle between each gene expression vector and the space diagonal.
#' Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
#'
#' The selection-consistency checks (`n_selected_axes` as a dimension, and each selection count
#' matching its claimed total) live here: they compare a `count(mask)` against a claimed size, which
#' the generated wrapper's per-argument validators cannot express.
#'
#' @param expression_vectors a numeric matrix. 2D array (n_axes, n_vectors), each column is a gene expression vector
#' @param exp_vecs_selection_index a logical vector. Logical array (n_vectors), .TRUE. for vectors to process
#' @param n_selected_vectors a integer scalar. Number of selected expression vectors (count of .TRUE. in exp_vecs_selection_index)
#' @param axes_selection a logical vector. Logical array (n_axes), .TRUE. for axes to include in calculation
#' @param n_selected_axes a integer scalar. Number of selected axes (count of .TRUE. in axes_selection)
#' @return a named list with elements `tissue_versatilities`, `tissue_angles_deg`.
#'
#' Generated from the Fortran procedure \code{tox_tissue_versatility::compute_tissue_versatility}.
#' @export
compute_tissue_versatility <- function(expression_vectors, exp_vecs_selection_index, n_selected_vectors, axes_selection, n_selected_axes) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    exp_vecs_selection_index <- .tox_as_logical(exp_vecs_selection_index, "exp_vecs_selection_index")
    n_selected_vectors <- .tox_as_integer_scalar(n_selected_vectors, "n_selected_vectors")
    axes_selection <- .tox_as_logical(axes_selection, "axes_selection")
    n_selected_axes <- .tox_as_integer_scalar(n_selected_axes, "n_selected_axes")
    if (length(axes_selection) != dim(expression_vectors)[1])
        .tox_shape_error("axes_selection", length(axes_selection), "expression_vectors", dim(expression_vectors)[1])
    if (length(exp_vecs_selection_index) != dim(expression_vectors)[2])
        .tox_shape_error("exp_vecs_selection_index", length(exp_vecs_selection_index), "expression_vectors", dim(expression_vectors)[2])

    .result <- .Call("compute_tissue_versatility_call", expression_vectors, exp_vecs_selection_index, n_selected_vectors, axes_selection, n_selected_axes)
    .arguments <- c("n_axes", "n_vectors", "expression_vectors", "exp_vecs_selection_index", "n_selected_vectors", "axes_selection", "n_selected_axes", "tissue_versatilities", "tissue_angles_deg", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        tissue_versatilities = .result$tissue_versatilities,
        tissue_angles_deg = .result$tissue_angles_deg
    )
}
