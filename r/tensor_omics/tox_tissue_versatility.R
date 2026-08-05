# Generated. Do not edit.

#' Computes normalized tissue versatility for selected expression vectors.
#'
#' The metric is based on the angle between each gene expression vector and the space diagonal.
#' Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
#'
#' The masks follow the `n_selected_` convention, so the generated wrapper validates that each
#' selection count matches its mask; `n_selected_axes` (not an array extent) carries its own floor.
#'
#' Generated from the Fortran module \code{tox_tissue_versatility}.
#'
#' @param expression_vectors a numeric matrix. 2D array (n_axes, n_vectors), each column is a gene expression vector
#' @param vectors_selection_mask a logical vector. Logical array (n_vectors), TRUE for vectors to process
#' @param axes_selection_mask a logical vector. Logical array (n_axes), TRUE for axes to include in calculation
#' @return a named list with elements:
#'   \item{tissue_versatilities}{a numeric vector. Output, real array, length = n_selected_vectors, stores the calculated tissue versatilities}
#'   \item{tissue_angles_deg}{a numeric vector. Output, real array, length = n_selected_vectors, stores the calculated angles in degrees}
#' @export
compute_tissue_versatility <- function(expression_vectors, vectors_selection_mask, axes_selection_mask) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    vectors_selection_mask <- .tox_as_logical(vectors_selection_mask, "vectors_selection_mask")
    axes_selection_mask <- .tox_as_logical(axes_selection_mask, "axes_selection_mask")
    if (length(axes_selection_mask) != dim(expression_vectors)[1])
        .tox_shape_error("axes_selection_mask", length(axes_selection_mask), "expression_vectors", dim(expression_vectors)[1])
    if (length(vectors_selection_mask) != dim(expression_vectors)[2])
        .tox_shape_error("vectors_selection_mask", length(vectors_selection_mask), "expression_vectors", dim(expression_vectors)[2])

    .result <- .Call("compute_tissue_versatility_call", expression_vectors, vectors_selection_mask, axes_selection_mask)
    .arguments <- c("n_axes", "n_vectors", "expression_vectors", "vectors_selection_mask", "n_selected_vectors", "axes_selection_mask", "n_selected_axes", "tissue_versatilities", "tissue_angles_deg", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        tissue_versatilities = .result$tissue_versatilities,
        tissue_angles_deg = .result$tissue_angles_deg
    )
}
