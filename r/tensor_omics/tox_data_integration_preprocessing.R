# Generated. Do not edit.

#' Compute per-gene mean expression, ignoring NaN values
#'
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing::compute_gene_means}, whose argument names
#' are the ones an error message reports.
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
    .sources <- c("expr", "expr", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$means
}

#' Compute signed residuals (centering by mean)
#'
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing::compute_residuals}, whose argument names
#' are the ones an error message reports.
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
    .sources <- c("expr", "expr", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$resid
}

#' Turn a sorted pool of per-gene mean expression values into reference points
#'
#' This takes the pool already built; `pool_study_means` pools the means of two studies
#' first, if that is what is at hand.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing::pool_means}, whose argument names
#' are the ones an error message reports.
#'
#' This entry point seeds \code{pooled_means_perm} and sorts it by \code{pooled_means}.
#' Call \code{pool_means_expert} to do that yourself.
#'
#' @param pooled_means a numeric vector. Pooled means
#'   NaN is permitted for this value.
#' @param n_points a integer scalar. Number of reference points to define
#' @return a named list with elements:
#'   \item{n_pool}{a integer scalar. Total number of included (non-NaN) pooled mean-expression values}
#'   \item{x_star}{a numeric vector. Mean-expression reference points. If every pooled mean is NaN, every reference
#'     point is NaN, and `construct_neighborhoods_ranged` rejects such points.}
#' @export
pool_means <- function(pooled_means, n_points) {
    pooled_means <- .tox_as_double_vector(pooled_means, "pooled_means")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    .result <- .Call("pool_means_call", pooled_means, n_points)
    .arguments <- c("pooled_means", "pool_size", "n_points", "n_pool", "x_star", "ierr")
    .sources <- c(NA_character_, "pooled_means", "x_star", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

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
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing::pool_means_expert}, whose argument names
#' are the ones an error message reports.
#'
#' The expert entry point: you supply \code{pooled_means_perm} yourself.
#' \code{pool_means} seeds \code{pooled_means_perm} and sorts it by \code{pooled_means}.
#'
#' @param pooled_means a numeric vector. Pooled means
#'   NaN is permitted for this value.
#' @param pooled_means_perm a integer vector. Sorting permutation for `pooled_means`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `pool_size`.
#' @param n_points a integer scalar. Number of reference points to define
#' @return a named list with elements:
#'   \item{n_pool}{a integer scalar. Total number of included (non-NaN) pooled mean-expression values}
#'   \item{x_star}{a numeric vector. Mean-expression reference points. If every pooled mean is NaN, every reference
#'     point is NaN, and `construct_neighborhoods_ranged` rejects such points.}
#' @export
pool_means_expert <- function(pooled_means, pooled_means_perm, n_points) {
    pooled_means <- .tox_as_double_vector(pooled_means, "pooled_means")
    pooled_means_perm <- .tox_as_integer_vector(pooled_means_perm, "pooled_means_perm")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    if (length(pooled_means_perm) != length(pooled_means))
        .tox_shape_error("pooled_means_perm", length(pooled_means_perm), "pooled_means", length(pooled_means))

    .result <- .Call("pool_means_expert_call", pooled_means, pooled_means_perm, n_points)
    .arguments <- c("pooled_means", "pooled_means_perm", "pool_size", "n_points", "n_pool", "x_star", "ierr")
    .sources <- c(NA_character_, NA_character_, "pooled_means", "x_star", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

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
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing::pool_study_means}, whose argument names
#' are the ones an error message reports.
#'
#' @param mean_S1 a numeric vector. Per-gene mean expression values of study S1
#'   NaN is permitted for this value.
#' @param mean_S2 a numeric vector. Per-gene mean expression values of study S2
#'   NaN is permitted for this value.
#' @param n_points a integer scalar. Number of reference points to define
#' @return a named list with elements:
#'   \item{n_pool}{a integer scalar. Total number of included (non-NaN) pooled mean-expression values}
#'   \item{x_star}{a numeric vector. Mean-expression reference points. If every mean of both studies is NaN, every
#'     reference point is NaN, and `construct_neighborhoods_ranged` rejects such points.}
#' @export
pool_study_means <- function(mean_S1, mean_S2, n_points) {
    mean_S1 <- .tox_as_double_vector(mean_S1, "mean_S1")
    mean_S2 <- .tox_as_double_vector(mean_S2, "mean_S2")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    .result <- .Call("pool_study_means_call", mean_S1, mean_S2, n_points)
    .arguments <- c("n_genes_S1", "mean_S1", "n_genes_S2", "mean_S2", "n_points", "n_pool", "x_star", "ierr")
    .sources <- c("mean_S1", NA_character_, "mean_S2", NA_character_, "x_star", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        n_pool = .result$n_pool,
        x_star = .result$x_star
    )
}

#' Construct neighborhood-based residual sets (kNN)
#'
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing::construct_neighborhoods}, whose argument names
#' are the ones an error message reports.
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
#'   \code{\link{calc_neighborhood_size}}.
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
    .sources <- c("x_star", NA_character_, "mean_S", NA_character_, "resid_S", NA_character_, NA_character_, NA_character_, "neighborhood_residuals", NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        neighborhood_residuals = .result$neighborhood_residuals,
        neighborhood_indices = .result$neighborhood_indices
    )
}

#' Construct neighborhood-based residual sets (kNN), with a `[min_idx, max_idx]` range per reference point
#'
#' Ported from 125-stabilize-jscomp's `construct_neighborhoods_helper`: a binary-search +
#' two-pointer kNN construction over a pre-sorted `mean_S`, distinct from
#' \code{\link{construct_neighborhoods}}'s
#' distance-sort algorithm above. In addition to the neighborhood gene indices, this also
#' reports each reference point's `[min_idx, max_idx]` neighborhood span (tie-extended, so
#' genes tied with the span's edge value are never split from it), which a candidate
#' admissibility gate elsewhere in the JSD-Comp-Test parameter search reasons about directly,
#' without needing the gathered residuals.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing::construct_neighborhoods_ranged}, whose argument names
#' are the ones an error message reports.
#'
#' This entry point seeds \code{mean_S_perm} and sorts it by \code{mean_S}.
#' Call \code{construct_neighborhoods_ranged_expert} to do that yourself.
#'
#' @param x_star a numeric vector. Mean-expression reference points
#' @param mean_S a numeric vector. Per-gene mean expression values
#'   NaN is permitted for this value.
#' @param n_neighbors a integer scalar. Number of neighbors to select per reference point
#'   The minimum valid value is `1`.
#' @return a named list with elements:
#'   \item{neighborhood_indices}{a integer matrix. Indices of selected neighborhood genes per reference point.
#'
#'     All indices are in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case
#'     `n_genes_S` is lower than `n_neighbors`, remaining indices are filled with the
#'     ones from `n_genes_S+1...n_neighbors` (a documented, deliberately-preserved
#'     limitation, ported as-is from 125-stabilize-jscomp).}
#'   \item{neighborhood_range}{a integer matrix. For each reference point, the `[min_idx, max_idx]` of the included genes. The
#'     index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))`
#'     would be the min value. In case of duplicate means, `min_idx` points to the first
#'     appearance of the value and `max_idx` to the last, so even though their related
#'     mean value is the min/max in the neighborhood, the actual gene might not be
#'     included. If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`}
#' @export
construct_neighborhoods_ranged <- function(x_star, mean_S, n_neighbors) {
    x_star <- .tox_as_double_vector(x_star, "x_star")
    mean_S <- .tox_as_double_vector(mean_S, "mean_S")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    .result <- .Call("construct_neighborhoods_ranged_call", x_star, mean_S, n_neighbors)
    .arguments <- c("n_points", "x_star", "n_genes_S", "mean_S", "n_neighbors", "neighborhood_indices", "neighborhood_range", "ierr")
    .sources <- c("x_star", NA_character_, "mean_S", NA_character_, "neighborhood_indices", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        neighborhood_indices = .result$neighborhood_indices,
        neighborhood_range = .result$neighborhood_range
    )
}

#' Construct neighborhood-based residual sets (kNN), with a `[min_idx, max_idx]` range per reference point
#'
#' Ported from 125-stabilize-jscomp's `construct_neighborhoods_helper`: a binary-search +
#' two-pointer kNN construction over a pre-sorted `mean_S`, distinct from
#' \code{\link{construct_neighborhoods}}'s
#' distance-sort algorithm above. In addition to the neighborhood gene indices, this also
#' reports each reference point's `[min_idx, max_idx]` neighborhood span (tie-extended, so
#' genes tied with the span's edge value are never split from it), which a candidate
#' admissibility gate elsewhere in the JSD-Comp-Test parameter search reasons about directly,
#' without needing the gathered residuals.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing::construct_neighborhoods_ranged_expert}, whose argument names
#' are the ones an error message reports.
#'
#' The expert entry point: you supply \code{mean_S_perm} yourself.
#' \code{construct_neighborhoods_ranged} seeds \code{mean_S_perm} and sorts it by \code{mean_S}.
#'
#' @param x_star a numeric vector. Mean-expression reference points
#' @param mean_S a numeric vector. Per-gene mean expression values
#'   NaN is permitted for this value.
#' @param mean_S_perm a integer vector. Sorting permutation for `mean_S`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_genes_S`.
#' @param n_neighbors a integer scalar. Number of neighbors to select per reference point
#'   The minimum valid value is `1`.
#' @return a named list with elements:
#'   \item{neighborhood_indices}{a integer matrix. Indices of selected neighborhood genes per reference point.
#'
#'     All indices are in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case
#'     `n_genes_S` is lower than `n_neighbors`, remaining indices are filled with the
#'     ones from `n_genes_S+1...n_neighbors` (a documented, deliberately-preserved
#'     limitation, ported as-is from 125-stabilize-jscomp).}
#'   \item{neighborhood_range}{a integer matrix. For each reference point, the `[min_idx, max_idx]` of the included genes. The
#'     index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))`
#'     would be the min value. In case of duplicate means, `min_idx` points to the first
#'     appearance of the value and `max_idx` to the last, so even though their related
#'     mean value is the min/max in the neighborhood, the actual gene might not be
#'     included. If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`}
#' @export
construct_neighborhoods_ranged_expert <- function(x_star, mean_S, mean_S_perm, n_neighbors) {
    x_star <- .tox_as_double_vector(x_star, "x_star")
    mean_S <- .tox_as_double_vector(mean_S, "mean_S")
    mean_S_perm <- .tox_as_integer_vector(mean_S_perm, "mean_S_perm")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    if (length(mean_S_perm) != length(mean_S))
        .tox_shape_error("mean_S_perm", length(mean_S_perm), "mean_S", length(mean_S))

    .result <- .Call("construct_neighborhoods_ranged_expert_call", x_star, mean_S, mean_S_perm, n_neighbors)
    .arguments <- c("n_points", "x_star", "n_genes_S", "mean_S", "mean_S_perm", "n_neighbors", "neighborhood_indices", "neighborhood_range", "ierr")
    .sources <- c("x_star", NA_character_, "mean_S", NA_character_, NA_character_, "neighborhood_indices", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        neighborhood_indices = .result$neighborhood_indices,
        neighborhood_range = .result$neighborhood_range
    )
}

#' Grow adaptive, dispersion-controlled neighborhoods over all studies' pooled gene means (Issue #217)
#'
#' Instead of one fixed neighbor count, each neighborhood grows while the mean absolute residual
#' of its pool stays stable, and the next reference point is placed beyond it; the number of
#' reference points emerges. `k_start`, `k_step` and `k_max` count pooled entries (one entry is
#' one gene of one study).
#'
#' **Pooling.** All studies' gene means are pooled into `N = max_n_genes_all_studies *
#' n_studies` entries (entry `(s - 1)*G + g` is gene `g` of study `s`) and sorted ascending; NaN
#' means take no part. A neighborhood is a contiguous run `a..b` of that order
#' (`pooled_neighborhood_range`, the input of
#' \code{\link{check_neighborhood_overlaps}}),
#' and its per-study counts `n_neighbors_per_point` add up to `b - a + 1`.
#'
#' **Growth** around the reference point `x_star` (the mean at the seed position): entries join
#' nearest first, ties to the lower mean; the dispersion `S` is the mean absolute non-NaN
#' residual of the pool.
#'
#' 1. The first `k_start` entries join unconditionally. With fewer than
#' \code{ADAPTIVE_MIN_VALID_RESIDUALS}
#' (`10`) non-NaN residuals the neighborhood stops with
#' \code{ADAPTIVE_STOP_TOO_FEW_RESIDUALS}
#' (`6`); with all residuals zero, with
#' \code{ADAPTIVE_STOP_ZERO_DISPERSION}
#' (`5`).
#' 2. Then rounds of the next `min(k_step, k_max - count)` entries: a round with
#' `(S_new - S_old) / S_old > tau` is discarded and growth stops
#' (\code{ADAPTIVE_STOP_TAU}, `1`);
#' otherwise it is committed. Growth also stops at `k_max` entries
#' (\code{ADAPTIVE_STOP_K_MAX}, `2`,
#' checked first) or when no entry is left
#' (\code{ADAPTIVE_STOP_EXHAUSTED}, `4`).
#' 3. With `max_pooled_residuals > 0`, a round stages whole entries only while the pool's
#' non-NaN residuals stay within that cap (the `k_start` phase is not capped); a capped round
#' is still evaluated, and growth then stops with
#' \code{ADAPTIVE_STOP_RESIDUAL_CAP} (`3`).
#'
#' **Reference points.** The first seed is the smallest pooled mean. After a neighborhood ends at
#' sorted position `b`, the next seed is the first position after `b` whose mean is nearest
#' `mean(b) + mad_distance_factor * MAD` (ties to the lower mean), `MAD` being the raw median
#' absolute deviation of the neighborhood's pooled means. The construction ends with the
#' neighborhood that reaches the largest pooled mean.
#'
#' **Construction status.**
#' \code{ADAPTIVE_STATUS_TOO_FEW_MEANS}
#' (`1`) when fewer than `k_start` pooled means are non-NaN: nothing
#' is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction always completes,
#' and reports, in this priority,
#' \code{ADAPTIVE_STATUS_TOO_FEW_RESIDUALS}
#' (`3`) if any neighborhood stopped with too few residuals,
#' \code{ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD}
#' (`2`) if any has fewer than `min_study_neighbors`
#' entries of some study, else
#' \code{ADAPTIVE_STATUS_OK} (`0`).
#'
#' Outputs sized by `N` are filled only for their first `n_points` reference points.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing::construct_adaptive_neighborhoods}, whose argument names
#' are the ones an error message reports.
#'
#' @param gene_means a numeric matrix. Mean expression of every gene in every study, NaN for a missing gene
#'   NaN is permitted for this value.
#' @param residuals a numeric array of rank 3. Signed residuals of every replicate of every gene in every study, NaN for a missing value
#'   NaN is permitted for this value.
#' @param k_start a integer scalar. Pooled entries every neighborhood takes unconditionally
#'   The minimum valid value is `1`.
#' @param k_step a integer scalar. Pooled entries staged per adaptive growth round
#'   The minimum valid value is `1`.
#' @param k_max a integer scalar. Largest number of pooled entries in one neighborhood
#'   The minimum valid value is `k_start`.
#' @param tau a numeric scalar. Largest relative dispersion increase a growth round may cause and still be committed
#'   The minimum valid value is `0.0`.
#'   The default value is `0.1`.
#' @param mad_distance_factor a numeric scalar. Multiple of a neighborhood's MAD that the next reference point's target lies beyond
#'   the neighborhood's largest mean
#'   The minimum valid value is `0.0`.
#'   The default value is `1.0`.
#' @param max_pooled_residuals a integer scalar. Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
#'   The minimum valid value is `0`.
#'   The default value is `0`.
#' @param min_study_neighbors a integer scalar. Fewest entries of each study every neighborhood must have for the construction
#'   status to stay ok
#'   The minimum valid value is `1`.
#'   The default value is `1`.
#' @return a named list with elements:
#'   \item{x_star}{a numeric vector. Reference point of each neighborhood: the pooled mean at its seed position
#'     The first `n_points` elements will hold the results.}
#'   \item{pooled_neighborhood_range}{a integer matrix. For each neighborhood, its first and last position in the ascending order of the
#'     pooled means
#'     The first `n_points` elements will hold the results.}
#'   \item{n_neighbors_per_point}{a integer matrix. For each neighborhood, how many of its pooled entries belong to each study
#'     The first `n_points` elements will hold the results.}
#'   \item{stop_reason}{a integer vector. Why each neighborhood stopped growing: `1` tau exceeded,
#'     `2` `k_max` reached, `3` residual cap
#'     reached, `4` pool exhausted, `5`
#'     zero dispersion, `6` too few residuals
#'     The first `n_points` elements will hold the results.}
#'   \item{neighborhood_dispersion}{a numeric vector. Mean absolute non-NaN residual of each final neighborhood; NaN without any
#'     The first `n_points` elements will hold the results.}
#'   \item{neighborhood_mad}{a numeric vector. Raw median absolute deviation of each final neighborhood's pooled means
#'     The first `n_points` elements will hold the results.}
#'   \item{max_n_neighbors}{a integer scalar. Largest number of entries one study has in one neighborhood; 0 without a neighborhood}
#'   \item{construction_status}{a integer scalar. Overall outcome: `0` ok, `1` too few
#'     means, `2` a study neighborhood below
#'     `min_study_neighbors`, `3` a neighborhood with too
#'     few residuals}
#' @export
construct_adaptive_neighborhoods <- function(gene_means, residuals, k_start, k_step, k_max, tau = 0.1, mad_distance_factor = 1.0, max_pooled_residuals = 0L, min_study_neighbors = 1L) {
    gene_means <- .tox_as_double_matrix(gene_means, "gene_means")
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    k_start <- .tox_as_integer_scalar(k_start, "k_start")
    k_step <- .tox_as_integer_scalar(k_step, "k_step")
    k_max <- .tox_as_integer_scalar(k_max, "k_max")
    tau <- .tox_as_double_scalar(tau, "tau")
    mad_distance_factor <- .tox_as_double_scalar(mad_distance_factor, "mad_distance_factor")
    max_pooled_residuals <- .tox_as_integer_scalar(max_pooled_residuals, "max_pooled_residuals")
    min_study_neighbors <- .tox_as_integer_scalar(min_study_neighbors, "min_study_neighbors")
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("construct_adaptive_neighborhoods_call", gene_means, residuals, k_start, k_step, k_max, tau, mad_distance_factor, max_pooled_residuals, min_study_neighbors)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "k_start", "k_step", "k_max", "n_points", "x_star", "pooled_neighborhood_range", "n_neighbors_per_point", "stop_reason", "neighborhood_dispersion", "neighborhood_mad", "max_n_neighbors", "construction_status", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        x_star = utils::head(.result$x_star, .result$n_points),
        pooled_neighborhood_range = .result$pooled_neighborhood_range[, seq_len(.result$n_points), drop = FALSE],
        n_neighbors_per_point = .result$n_neighbors_per_point[, seq_len(.result$n_points), drop = FALSE],
        stop_reason = utils::head(.result$stop_reason, .result$n_points),
        neighborhood_dispersion = utils::head(.result$neighborhood_dispersion, .result$n_points),
        neighborhood_mad = utils::head(.result$neighborhood_mad, .result$n_points),
        max_n_neighbors = .result$max_n_neighbors,
        construction_status = .result$construction_status
    )
}

#' Test whether every pair of consecutive neighborhoods overlaps by at least a minimum fraction
#'
#' The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
#' neighborhood spans of
#' \code{\link{construct_neighborhoods_ranged}}
#' (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
#' number of positions they share divided by the length of the earlier span; a reversed span
#' (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
#' `test_`-prefixed export as a test case.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing::check_neighborhood_overlaps}, whose argument names
#' are the ones an error message reports.
#'
#' @param neighborhood_range a integer matrix. For each reference point, the `[min_idx, max_idx]` neighborhood span
#'   The minimum valid value is `1`.
#' @param min_neighbor_overlap a numeric scalar. Minimum fractional overlap two consecutive neighborhoods must have
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @return a logical scalar. `TRUE` if every pair of consecutive neighborhoods overlaps by at least
#'   `min_neighbor_overlap`
#' @export
check_neighborhood_overlaps <- function(neighborhood_range, min_neighbor_overlap) {
    neighborhood_range <- .tox_as_integer_matrix(neighborhood_range, "neighborhood_range")
    min_neighbor_overlap <- .tox_as_double_scalar(min_neighbor_overlap, "min_neighbor_overlap")
    .result <- .Call("check_neighborhood_overlaps_call", neighborhood_range, min_neighbor_overlap)
    .arguments <- c("neighborhood_range", "n_points", "min_neighbor_overlap", "all_have_min_neighbor_overlap", "ierr")
    .sources <- c(NA_character_, "neighborhood_range", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$all_have_min_neighbor_overlap
}
