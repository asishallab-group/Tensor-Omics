# Generated. Do not edit.

#' Compute per-gene mean expression, ignoring NaN values
#'
#' Generated from the Fortran module \code{tox_data_integration_preprocessing}.
#'
#' @param expr a numeric matrix. Expression matrix
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @return a numeric vector. Per-gene mean expression values
#' @export
compute_gene_means <- function(expr) {
    expr <- .tox_as_double_matrix(expr, "expr")
    .result <- .Call("compute_gene_means_call", expr)
    .arguments <- c("n_genes", "n_reps", "expr", "means", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$means
}

#' Compute signed residuals (centering by mean)
#'
#' Generated from the Fortran module \code{tox_data_integration_preprocessing}.
#'
#' @param expr a numeric matrix. Expression matrix
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param means a numeric vector. Per-gene mean expression values; NaN where every replicate of a gene was NaN
#'   NaN is permitted for this value.
#' @return a numeric matrix. Matrix of signed residuals
#' @export
compute_residuals <- function(expr, means) {
    expr <- .tox_as_double_matrix(expr, "expr")
    means <- .tox_as_double_vector(means, "means")
    if (length(means) != dim(expr)[2])
        .tox_shape_error("means", length(means), "expr", dim(expr)[2])

    .result <- .Call("compute_residuals_call", expr, means)
    .arguments <- c("n_genes", "n_reps", "expr", "means", "resid", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$resid
}

#' Turn a sorted pool of per-gene mean expression values into reference points
#'
#' This takes the pool already built; `pool_study_means` pools the means of two studies
#' first, if that is what is at hand.
#'
#' Generated from the Fortran module \code{tox_data_integration_preprocessing}.
#'
#' @param pooled_means a numeric vector. Pooled means
#'   NaN is permitted for this value.
#' @param pooled_means_perm a integer vector. Sorting permutation for `pooled_means`
#' @param n_points a integer scalar. Number of reference points to define
#' @return a named list with elements:
#'   \item{n_pool}{a integer scalar. Total number of included (non-NaN) pooled mean-expression values}
#'   \item{x_star}{a numeric vector. Mean-expression reference points}
#' @export
pool_means_expert <- function(pooled_means, pooled_means_perm, n_points) {
    pooled_means <- .tox_as_double_vector(pooled_means, "pooled_means")
    pooled_means_perm <- .tox_as_integer_vector(pooled_means_perm, "pooled_means_perm")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    if (length(pooled_means_perm) != length(pooled_means))
        .tox_shape_error("pooled_means_perm", length(pooled_means_perm), "pooled_means", length(pooled_means))

    .result <- .Call("pool_means_expert_call", pooled_means, pooled_means_perm, n_points)
    .arguments <- c("pooled_means", "pooled_means_perm", "pool_size", "n_points", "n_pool", "x_star", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_pool = .result$n_pool,
        x_star = .result$x_star
    )
}

#' Turn a sorted pool of per-gene mean expression values into reference points
#'
#' This takes the pool already built; `pool_study_means` pools the means of two studies
#' first, if that is what is at hand.
#'
#' Generated from the Fortran module \code{tox_data_integration_preprocessing}.
#'
#' @param pooled_means a numeric vector. Pooled means
#'   NaN is permitted for this value.
#' @param n_points a integer scalar. Number of reference points to define
#' @return a named list with elements:
#'   \item{n_pool}{a integer scalar. Total number of included (non-NaN) pooled mean-expression values}
#'   \item{x_star}{a numeric vector. Mean-expression reference points}
#' @export
pool_means <- function(pooled_means, n_points) {
    pooled_means <- .tox_as_double_vector(pooled_means, "pooled_means")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    .result <- .Call("pool_means_call", pooled_means, n_points)
    .arguments <- c("pooled_means", "pool_size", "n_points", "n_pool", "x_star", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_pool = .result$n_pool,
        x_star = .result$x_star
    )
}

#' Pool the per-gene mean expression values of two studies into reference points
#'
#' Concatenates the two studies' means, sorts the pool, and turns it into reference
#' points exactly as `pool_means` does.
#'
#' Generated from the Fortran module \code{tox_data_integration_preprocessing}.
#'
#' @param mean_S1 a numeric vector. Per-gene mean expression values of study S1
#'   NaN is permitted for this value.
#' @param mean_S2 a numeric vector. Per-gene mean expression values of study S2
#'   NaN is permitted for this value.
#' @param n_points a integer scalar. Number of reference points to define
#' @return a named list with elements:
#'   \item{n_pool}{a integer scalar. Total number of included (non-NaN) pooled mean-expression values}
#'   \item{x_star}{a numeric vector. Mean-expression reference points}
#' @export
pool_study_means_expert <- function(mean_S1, mean_S2, n_points) {
    mean_S1 <- .tox_as_double_vector(mean_S1, "mean_S1")
    mean_S2 <- .tox_as_double_vector(mean_S2, "mean_S2")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    .result <- .Call("pool_study_means_expert_call", mean_S1, mean_S2, n_points)
    .arguments <- c("n_genes_S1", "mean_S1", "n_genes_S2", "mean_S2", "n_points", "tmp_pooled_means", "tmp_pooled_means_perm", "n_pool", "x_star", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_pool = .result$n_pool,
        x_star = .result$x_star
    )
}

#' Pool the per-gene mean expression values of two studies into reference points
#'
#' Concatenates the two studies' means, sorts the pool, and turns it into reference
#' points exactly as `pool_means` does.
#'
#' Generated from the Fortran module \code{tox_data_integration_preprocessing}.
#'
#' @param mean_S1 a numeric vector. Per-gene mean expression values of study S1
#'   NaN is permitted for this value.
#' @param mean_S2 a numeric vector. Per-gene mean expression values of study S2
#'   NaN is permitted for this value.
#' @param n_points a integer scalar. Number of reference points to define
#' @return a named list with elements:
#'   \item{n_pool}{a integer scalar. Total number of included (non-NaN) pooled mean-expression values}
#'   \item{x_star}{a numeric vector. Mean-expression reference points}
#' @export
pool_study_means <- function(mean_S1, mean_S2, n_points) {
    mean_S1 <- .tox_as_double_vector(mean_S1, "mean_S1")
    mean_S2 <- .tox_as_double_vector(mean_S2, "mean_S2")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    .result <- .Call("pool_study_means_call", mean_S1, mean_S2, n_points)
    .arguments <- c("n_genes_S1", "mean_S1", "n_genes_S2", "mean_S2", "n_points", "n_pool", "x_star", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_pool = .result$n_pool,
        x_star = .result$x_star
    )
}

#' Construct neighborhood-based residual sets (kNN)
#'
#' Generated from the Fortran module \code{tox_data_integration_preprocessing}.
#'
#' @param x_star a numeric vector. Mean-expression reference points
#'   NaN is permitted for this value.
#' @param mean_S a numeric vector. Per-gene mean expression values
#'   NaN is permitted for this value.
#' @param resid_S a numeric matrix. Matrix of signed residuals
#'   NaN is permitted for this value.
#' @param n_neighbors a integer scalar. Number of neighbors; a gene whose mean is NaN can never be a neighbor, so this
#'   cannot exceed the number of genes with a defined mean
#'   The minimum valid value is `1`.
#'   The maximum valid value is `count(.not. ieee_is_nan(mean_S), kind=int32)`.
#'   It is recommended to compute this with
#'   [[tox_data_integration_preprocessing_kernel(module):calc_neighborhood_size(function)]].
#' @return a named list with elements:
#'   \item{neighborhood_residuals}{a numeric array of rank 3. Collection of residual vectors for each neighborhood}
#'   \item{neighborhood_indices}{a integer matrix. Indices of selected neighborhood genes}
#' @export
construct_neighborhoods_expert <- function(x_star, mean_S, resid_S, n_neighbors) {
    x_star <- .tox_as_double_vector(x_star, "x_star")
    mean_S <- .tox_as_double_vector(mean_S, "mean_S")
    resid_S <- .tox_as_double_matrix(resid_S, "resid_S")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    if (dim(resid_S)[2] != length(mean_S))
        .tox_shape_error("resid_S", dim(resid_S)[2], "mean_S", length(mean_S))

    .result <- .Call("construct_neighborhoods_expert_call", x_star, mean_S, resid_S, n_neighbors)
    .arguments <- c("n_points", "x_star", "n_genes_S", "mean_S", "n_reps_S", "resid_S", "tmp_distances", "tmp_distances_perm", "neighborhood_residuals", "neighborhood_indices", "n_neighbors", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        neighborhood_residuals = .result$neighborhood_residuals,
        neighborhood_indices = .result$neighborhood_indices
    )
}

#' Construct neighborhood-based residual sets (kNN)
#'
#' Generated from the Fortran module \code{tox_data_integration_preprocessing}.
#'
#' @param x_star a numeric vector. Mean-expression reference points
#'   NaN is permitted for this value.
#' @param mean_S a numeric vector. Per-gene mean expression values
#'   NaN is permitted for this value.
#' @param resid_S a numeric matrix. Matrix of signed residuals
#'   NaN is permitted for this value.
#' @param n_neighbors a integer scalar. Number of neighbors; a gene whose mean is NaN can never be a neighbor, so this
#'   cannot exceed the number of genes with a defined mean
#'   The minimum valid value is `1`.
#'   The maximum valid value is `count(.not. ieee_is_nan(mean_S), kind=int32)`.
#'   It is recommended to compute this with
#'   [[tox_data_integration_preprocessing_kernel(module):calc_neighborhood_size(function)]].
#' @return a named list with elements:
#'   \item{neighborhood_residuals}{a numeric array of rank 3. Collection of residual vectors for each neighborhood}
#'   \item{neighborhood_indices}{a integer matrix. Indices of selected neighborhood genes}
#' @export
construct_neighborhoods <- function(x_star, mean_S, resid_S, n_neighbors) {
    x_star <- .tox_as_double_vector(x_star, "x_star")
    mean_S <- .tox_as_double_vector(mean_S, "mean_S")
    resid_S <- .tox_as_double_matrix(resid_S, "resid_S")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    if (dim(resid_S)[2] != length(mean_S))
        .tox_shape_error("resid_S", dim(resid_S)[2], "mean_S", length(mean_S))

    .result <- .Call("construct_neighborhoods_call", x_star, mean_S, resid_S, n_neighbors)
    .arguments <- c("n_points", "x_star", "n_genes_S", "mean_S", "n_reps_S", "resid_S", "neighborhood_residuals", "neighborhood_indices", "n_neighbors", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        neighborhood_residuals = .result$neighborhood_residuals,
        neighborhood_indices = .result$neighborhood_indices
    )
}
