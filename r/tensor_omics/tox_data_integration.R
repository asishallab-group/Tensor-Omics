# Generated. Do not edit.

#' Compute per-gene mean expression, ignoring NaN values
#'
#' @param expr a numeric matrix. Expression matrix
#' @return Per-gene mean expression values
#'
#' Generated from the Fortran procedure \code{tox_data_integration::compute_gene_means}.
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
#' @param expr a numeric matrix. Expression matrix containing
#' @param means a numeric vector. Per-gene mean expression values
#' @return Matrix of signed residuals
#'
#' Generated from the Fortran procedure \code{tox_data_integration::compute_residuals}.
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

#' Pool per-gene mean expression values across studies
#'
#' @param mean_S1 a numeric vector. Per-gene mean expression values
#' @param mean_S2 a numeric vector. Per-gene mean expression values
#' @param n_points a integer scalar. Number of reference points to define
#' @return a named list with elements `n_pool`, `x_star`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration::pool_means_alloc}.
#' @export
pool_means <- function(mean_S1, mean_S2, n_points) {
    mean_S1 <- .tox_as_double_vector(mean_S1, "mean_S1")
    mean_S2 <- .tox_as_double_vector(mean_S2, "mean_S2")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    .result <- .Call("pool_means_call", mean_S1, mean_S2, n_points)
    .arguments <- c("n_genes_S1", "mean_S1", "n_genes_S2", "mean_S2", "n_points", "n_pool", "x_star", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_pool = .result$n_pool,
        x_star = .result$x_star
    )
}

#' Pool per-gene mean expression values across studies (expert entry point)
#'
#' @param pooled_means a numeric vector. Pooled means
#' @param pooled_means_perm a integer vector. Sorting permutation for `pooled_means`
#' @param n_points a integer scalar. Number of reference points to define
#' @return a named list with elements `n_pool`, `x_star`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration::pool_means}.
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

#' Calculate the number of neighbors to be used for constructing neighborhoods
#'
#' The `desired_size` works as upper limit, as the actual neighborhood size might be lower due to few genes with non-NaN mean.
#'
#' @param n_pool a integer scalar. Total number of pooled mean-expression values across both studies
#' @param n_points a integer scalar. Number of reference points
#' @param mean_S a numeric vector. Per-gene mean expression values
#' @param desired_size a integer scalar. Optional desired neighborhood size.
#' @return Calculated neighborhood size
#'
#' Generated from the Fortran procedure \code{tox_data_integration::calc_neighborhood_size}.
#' @export
calc_neighborhood_size <- function(n_pool, n_points, mean_S, desired_size = 1000L) {
    n_pool <- .tox_as_integer_scalar(n_pool, "n_pool")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    mean_S <- .tox_as_double_vector(mean_S, "mean_S")
    desired_size <- .tox_as_integer_scalar(desired_size, "desired_size")
    .result <- .Call("calc_neighborhood_size_call", n_pool, n_points, mean_S, desired_size)
    .arguments <- c("n_pool", "n_points", "n_genes_S", "mean_S", "desired_size")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$n_neighbors
}

#' Construct neighborhood-based residual sets (kNN)
#'
#' @param x_star a numeric vector. Mean-expression reference points
#' @param mean_S a numeric vector. Per-gene mean expression values
#' @param resid_S a numeric matrix. Matrix of signed residuals
#' @param n_neighbors a integer scalar. Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**
#' @return a named list with elements `neighborhood_residuals`, `neighborhood_indices`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration::construct_neighborhoods_alloc}.
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

#' Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
#'
#' @param neighborhood_residuals_S1 a numeric array of rank 3. Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
#' @param neighborhood_residuals_S2 a numeric array of rank 3. Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
#' @param global_jsd_observed a numeric scalar. Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
#' @param n_bins a integer scalar. Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
#' @param shared_residual_range a numeric scalar. Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
#' @param n_permutations a integer scalar. Number of permutations to perform
#' @param random_seed a integer scalar. Seed to use for shuffling
#' @param neighbor_mask_S1 a logical matrix. Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
#' @param neighbor_mask_S2 a logical matrix. Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
#' @return a named list with elements `jsd_null`, `p_value`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration::gjct_permutation_test_alloc}.
#' @export
gjct_permutation_test <- function(neighborhood_residuals_S1, neighborhood_residuals_S2, global_jsd_observed, n_bins, shared_residual_range, n_permutations, random_seed = NULL, neighbor_mask_S1 = NULL, neighbor_mask_S2 = NULL) {
    neighborhood_residuals_S1 <- .tox_as_double_array(neighborhood_residuals_S1, "neighborhood_residuals_S1", 3L)
    neighborhood_residuals_S2 <- .tox_as_double_array(neighborhood_residuals_S2, "neighborhood_residuals_S2", 3L)
    global_jsd_observed <- .tox_as_double_scalar(global_jsd_observed, "global_jsd_observed")
    n_bins <- .tox_as_integer_scalar(n_bins, "n_bins")
    shared_residual_range <- .tox_as_double_scalar(shared_residual_range, "shared_residual_range")
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    if (!is.null(random_seed))
        random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    if (!is.null(neighbor_mask_S1))
        neighbor_mask_S1 <- .tox_as_logical(neighbor_mask_S1, "neighbor_mask_S1")
    if (!is.null(neighbor_mask_S2))
        neighbor_mask_S2 <- .tox_as_logical(neighbor_mask_S2, "neighbor_mask_S2")
    if (dim(neighborhood_residuals_S2)[2] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_residuals_S2)[3] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[3], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])

    .result <- .Call("gjct_permutation_test_call", neighborhood_residuals_S1, neighborhood_residuals_S2, global_jsd_observed, n_bins, shared_residual_range, n_permutations, random_seed, neighbor_mask_S1, neighbor_mask_S2)
    .arguments <- c("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "global_jsd_observed", "n_bins", "shared_residual_range", "n_permutations", "jsd_null", "p_value", "ierr", "random_seed", "neighbor_mask_S1", "neighbor_mask_S2")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        jsd_null = .result$jsd_null,
        p_value = .result$p_value
    )
}

#' Estimates the permutation-test p-value (expert entry point with caller-provided work arrays)
#'
#' @param neighborhood_residuals_S1_copy a numeric array of rank 3. Copy (if wanted) of the computed neighborhood residuals for study 1, will be shuffled in-place
#' @param neighborhood_residuals_S2_copy a numeric array of rank 3. Copy (if wanted) of the computed neighborhood residuals for study 2, will be shuffled in-place
#' @param global_jsd_observed a numeric scalar. Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
#' @param n_bins a integer scalar. Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
#' @param shared_residual_range a numeric scalar. Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
#' @param n_permutations a integer scalar. Number of permutations to perform
#' @param random_seed a integer scalar. Seed to use for shuffling
#' @param neighbor_mask_S1 a logical matrix. Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
#' @param neighbor_mask_S2 a logical matrix. Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
#' @return a named list with elements `neighborhood_residuals_S1_copy`, `neighborhood_residuals_S2_copy`, `jsd_null`, `p_value`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration::gjct_permutation_test}.
#' @export
gjct_permutation_test_expert <- function(neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, global_jsd_observed, n_bins, shared_residual_range, n_permutations, random_seed = NULL, neighbor_mask_S1 = NULL, neighbor_mask_S2 = NULL) {
    neighborhood_residuals_S1_copy <- .tox_as_double_array(neighborhood_residuals_S1_copy, "neighborhood_residuals_S1_copy", 3L)
    neighborhood_residuals_S2_copy <- .tox_as_double_array(neighborhood_residuals_S2_copy, "neighborhood_residuals_S2_copy", 3L)
    global_jsd_observed <- .tox_as_double_scalar(global_jsd_observed, "global_jsd_observed")
    n_bins <- .tox_as_integer_scalar(n_bins, "n_bins")
    shared_residual_range <- .tox_as_double_scalar(shared_residual_range, "shared_residual_range")
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    if (!is.null(random_seed))
        random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    if (!is.null(neighbor_mask_S1))
        neighbor_mask_S1 <- .tox_as_logical(neighbor_mask_S1, "neighbor_mask_S1")
    if (!is.null(neighbor_mask_S2))
        neighbor_mask_S2 <- .tox_as_logical(neighbor_mask_S2, "neighbor_mask_S2")
    if (dim(neighborhood_residuals_S2_copy)[2] != dim(neighborhood_residuals_S1_copy)[2])
        .tox_shape_error("neighborhood_residuals_S2_copy", dim(neighborhood_residuals_S2_copy)[2], "neighborhood_residuals_S1_copy", dim(neighborhood_residuals_S1_copy)[2])
    if (dim(neighborhood_residuals_S2_copy)[3] != dim(neighborhood_residuals_S1_copy)[3])
        .tox_shape_error("neighborhood_residuals_S2_copy", dim(neighborhood_residuals_S2_copy)[3], "neighborhood_residuals_S1_copy", dim(neighborhood_residuals_S1_copy)[3])

    .result <- .Call("gjct_permutation_test_expert_call", neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, global_jsd_observed, n_bins, shared_residual_range, n_permutations, random_seed, neighbor_mask_S1, neighbor_mask_S2)
    .arguments <- c("neighborhood_residuals_S1_copy", "neighborhood_residuals_S2_copy", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "global_jsd_observed", "n_bins", "shared_residual_range", "n_permutations", "jsd_null", "p_value", "tmp_pool", "tmp_pmf_S1", "tmp_pmf_S2", "tmp_counts", "tmp_included_n_reps_S1", "tmp_included_n_reps_S2", "tmp_js_divergences", "tmp_weights", "ierr", "random_seed", "neighbor_mask_S1", "neighbor_mask_S2")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        neighborhood_residuals_S1_copy = .result$neighborhood_residuals_S1_copy,
        neighborhood_residuals_S2_copy = .result$neighborhood_residuals_S2_copy,
        jsd_null = .result$jsd_null,
        p_value = .result$p_value
    )
}

#' Computes the shared residual range [-R, R] from S1/S2 residuals (expert entry point)
#'
#' @param abs_residual_pool a numeric vector. The absolute residual values of the concatenated S1,S2 residuals
#' @param abs_residual_pool_perm a integer vector. The permutation vector that sorts `abs_residual_pool`
#' @param residual_range_quantile a numeric scalar. Quantile for determining the residual range.
#' @return Computed residual range (R)
#'
#' Generated from the Fortran procedure \code{tox_data_integration::determine_shared_residual_range}.
#' @export
determine_shared_residual_range_expert <- function(abs_residual_pool, abs_residual_pool_perm, residual_range_quantile = 95.0) {
    abs_residual_pool <- .tox_as_double_vector(abs_residual_pool, "abs_residual_pool")
    abs_residual_pool_perm <- .tox_as_integer_vector(abs_residual_pool_perm, "abs_residual_pool_perm")
    residual_range_quantile <- .tox_as_double_scalar(residual_range_quantile, "residual_range_quantile")
    if (length(abs_residual_pool_perm) != length(abs_residual_pool))
        .tox_shape_error("abs_residual_pool_perm", length(abs_residual_pool_perm), "abs_residual_pool", length(abs_residual_pool))

    .result <- .Call("determine_shared_residual_range_expert_call", abs_residual_pool, abs_residual_pool_perm, residual_range_quantile)
    .arguments <- c("abs_residual_pool", "abs_residual_pool_perm", "pool_size", "shared_residual_range", "ierr", "residual_range_quantile")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$shared_residual_range
}

#' Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
#'
#' @param neighborhood_residuals_S1 a numeric array of rank 3. Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
#' @param neighborhood_residuals_S2 a numeric array of rank 3. Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
#' @param residual_range_quantile a numeric scalar. Quantile for determining the residual range.
#' @return Computed residual range (R)
#'
#' Generated from the Fortran procedure \code{tox_data_integration::determine_shared_residual_range_alloc}.
#' @export
determine_shared_residual_range <- function(neighborhood_residuals_S1, neighborhood_residuals_S2, residual_range_quantile = 95.0) {
    neighborhood_residuals_S1 <- .tox_as_double_array(neighborhood_residuals_S1, "neighborhood_residuals_S1", 3L)
    neighborhood_residuals_S2 <- .tox_as_double_array(neighborhood_residuals_S2, "neighborhood_residuals_S2", 3L)
    residual_range_quantile <- .tox_as_double_scalar(residual_range_quantile, "residual_range_quantile")
    if (dim(neighborhood_residuals_S2)[2] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_residuals_S2)[3] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[3], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])

    .result <- .Call("determine_shared_residual_range_call", neighborhood_residuals_S1, neighborhood_residuals_S2, residual_range_quantile)
    .arguments <- c("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "shared_residual_range", "ierr", "residual_range_quantile")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$shared_residual_range
}

#' Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
#'
#' @param neighborhood_residuals a numeric array of rank 3. Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
#' @param shared_residual_range a numeric scalar. Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
#' @param n_bins a integer scalar. Number of equally sized histogram bins in range [-R,R]
#' @param neighbor_mask a logical matrix. Optional mask to exclude specific neighbors (e.g. for family-wise analysis)
#' @return a named list with elements `counts`, `pmf`, `included_n_reps`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration::build_residual_histograms}.
#' @export
build_residual_histograms <- function(neighborhood_residuals, shared_residual_range, n_bins, neighbor_mask = NULL) {
    neighborhood_residuals <- .tox_as_double_array(neighborhood_residuals, "neighborhood_residuals", 3L)
    shared_residual_range <- .tox_as_double_scalar(shared_residual_range, "shared_residual_range")
    n_bins <- .tox_as_integer_scalar(n_bins, "n_bins")
    if (!is.null(neighbor_mask))
        neighbor_mask <- .tox_as_logical(neighbor_mask, "neighbor_mask")
    .result <- .Call("build_residual_histograms_call", neighborhood_residuals, shared_residual_range, n_bins, neighbor_mask)
    .arguments <- c("neighborhood_residuals", "n_reps", "n_neighbors", "n_points", "shared_residual_range", "n_bins", "counts", "pmf", "included_n_reps", "ierr", "neighbor_mask")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        counts = .result$counts,
        pmf = .result$pmf,
        included_n_reps = .result$included_n_reps
    )
}

#' Computes the Jensen-Shannon divergence per reference point from the histogram pmfs
#'
#' @param pmf_S1 a numeric matrix. Computed normalized histogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
#' @param pmf_S2 a numeric matrix. Computed normalized histogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 2
#' @return Jensen-Shannon divergence per reference point
#'
#' Generated from the Fortran procedure \code{tox_data_integration::compute_divergence_per_reference_point}.
#' @export
compute_divergence_per_reference_point <- function(pmf_S1, pmf_S2) {
    pmf_S1 <- .tox_as_double_matrix(pmf_S1, "pmf_S1")
    pmf_S2 <- .tox_as_double_matrix(pmf_S2, "pmf_S2")
    if (dim(pmf_S2)[1] != dim(pmf_S1)[1])
        .tox_shape_error("pmf_S2", dim(pmf_S2)[1], "pmf_S1", dim(pmf_S1)[1])
    if (dim(pmf_S2)[2] != dim(pmf_S1)[2])
        .tox_shape_error("pmf_S2", dim(pmf_S2)[2], "pmf_S1", dim(pmf_S1)[2])

    .result <- .Call("compute_divergence_per_reference_point_call", pmf_S1, pmf_S2)
    .arguments <- c("pmf_S1", "pmf_S2", "n_points", "n_bins", "js_divergences", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$js_divergences
}

#' Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences
#'
#' @param js_divergences a numeric vector. Jensen-Shannon divergence per reference point, computed for studies S1 and S2
#' @param included_n_reps_S1 a integer vector. Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
#' @param included_n_reps_S2 a integer vector. Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
#' @return a named list with elements `global_js_divergence`, `weights`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration::compute_weighted_global_divergence}.
#' @export
compute_weighted_global_divergence <- function(js_divergences, included_n_reps_S1, included_n_reps_S2) {
    js_divergences <- .tox_as_double_vector(js_divergences, "js_divergences")
    included_n_reps_S1 <- .tox_as_integer_vector(included_n_reps_S1, "included_n_reps_S1")
    included_n_reps_S2 <- .tox_as_integer_vector(included_n_reps_S2, "included_n_reps_S2")
    if (length(included_n_reps_S1) != length(js_divergences))
        .tox_shape_error("included_n_reps_S1", length(included_n_reps_S1), "js_divergences", length(js_divergences))
    if (length(included_n_reps_S2) != length(js_divergences))
        .tox_shape_error("included_n_reps_S2", length(included_n_reps_S2), "js_divergences", length(js_divergences))

    .result <- .Call("compute_weighted_global_divergence_call", js_divergences, included_n_reps_S1, included_n_reps_S2)
    .arguments <- c("js_divergences", "n_points", "included_n_reps_S1", "included_n_reps_S2", "global_js_divergence", "weights", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        global_js_divergence = .result$global_js_divergence,
        weights = .result$weights
    )
}

#' Computes the family-level compatibility score for a single gene family (`family_idx`)
#'
#' Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricting residual samples to genes belonging to the specified family
#'
#' @param family_idx a integer scalar. Index of the family that should be analyzed
#' @param gene_to_family_S1 a integer vector. Mapping for study 1: Each index (gene) holds the index of its family
#' @param gene_to_family_S2 a integer vector. Mapping for study 2: Each index (gene) holds the index of its family
#' @param neighborhood_residuals_S1 a numeric array of rank 3. Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
#' @param neighborhood_residuals_S2 a numeric array of rank 3. Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
#' @param neighborhood_genes_S1 a integer matrix. Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
#' @param neighborhood_genes_S2 a integer matrix. Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
#' @param n_bins a integer scalar. Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
#' @param shared_residual_range a numeric scalar. Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
#' @return a named list with elements `js_divergences`, `included_n_reps_S1`, `included_n_reps_S2`, `total_included_n_reps`, `global_js_divergence`, `weights`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration::fjct_compute_jsd_alloc}.
#' @export
fjct_compute_jsd <- function(family_idx, gene_to_family_S1, gene_to_family_S2, neighborhood_residuals_S1, neighborhood_residuals_S2, neighborhood_genes_S1, neighborhood_genes_S2, n_bins, shared_residual_range) {
    family_idx <- .tox_as_integer_scalar(family_idx, "family_idx")
    gene_to_family_S1 <- .tox_as_integer_vector(gene_to_family_S1, "gene_to_family_S1")
    gene_to_family_S2 <- .tox_as_integer_vector(gene_to_family_S2, "gene_to_family_S2")
    neighborhood_residuals_S1 <- .tox_as_double_array(neighborhood_residuals_S1, "neighborhood_residuals_S1", 3L)
    neighborhood_residuals_S2 <- .tox_as_double_array(neighborhood_residuals_S2, "neighborhood_residuals_S2", 3L)
    neighborhood_genes_S1 <- .tox_as_integer_matrix(neighborhood_genes_S1, "neighborhood_genes_S1")
    neighborhood_genes_S2 <- .tox_as_integer_matrix(neighborhood_genes_S2, "neighborhood_genes_S2")
    n_bins <- .tox_as_integer_scalar(n_bins, "n_bins")
    shared_residual_range <- .tox_as_double_scalar(shared_residual_range, "shared_residual_range")
    if (dim(neighborhood_residuals_S2)[2] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_genes_S1)[1] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_genes_S1", dim(neighborhood_genes_S1)[1], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_genes_S2)[1] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_genes_S2", dim(neighborhood_genes_S2)[1], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_residuals_S2)[3] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[3], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])
    if (dim(neighborhood_genes_S1)[2] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_genes_S1", dim(neighborhood_genes_S1)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])
    if (dim(neighborhood_genes_S2)[2] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_genes_S2", dim(neighborhood_genes_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])

    .result <- .Call("fjct_compute_jsd_call", family_idx, gene_to_family_S1, gene_to_family_S2, neighborhood_residuals_S1, neighborhood_residuals_S2, neighborhood_genes_S1, neighborhood_genes_S2, n_bins, shared_residual_range)
    .arguments <- c("family_idx", "gene_to_family_S1", "gene_to_family_S2", "n_genes_S1", "n_genes_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_genes_S1", "neighborhood_genes_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "n_bins", "shared_residual_range", "js_divergences", "included_n_reps_S1", "included_n_reps_S2", "total_included_n_reps", "global_js_divergence", "weights", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        js_divergences = .result$js_divergences,
        included_n_reps_S1 = .result$included_n_reps_S1,
        included_n_reps_S2 = .result$included_n_reps_S2,
        total_included_n_reps = .result$total_included_n_reps,
        global_js_divergence = .result$global_js_divergence,
        weights = .result$weights
    )
}

#' Computes the compatibility score for a single sub-neighborhood/family (expert entry point with caller-provided masks and work arrays)
#'
#' Restricts residual samples to the neighbors selected by `neighbor_mask_S1`/`neighbor_mask_S2` (typically all neighbors belonging to one gene family; see [[tox_data_integration(module):fjct_compute_jsd_alloc(interface)]] for the family-index-based entry point that builds these masks)
#'
#' @param neighborhood_residuals_S1 a numeric array of rank 3. Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
#' @param neighborhood_residuals_S2 a numeric array of rank 3. Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
#' @param neighbor_mask_S1 a logical matrix. Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
#' @param neighbor_mask_S2 a logical matrix. Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
#' @param n_bins a integer scalar. Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
#' @param shared_residual_range a numeric scalar. Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
#' @return a named list with elements `js_divergences`, `included_n_reps_S1`, `included_n_reps_S2`, `total_included_n_reps`, `global_js_divergence`, `weights`, `pmf_S1`, `pmf_S2`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration::fjct_compute_jsd}.
#' @export
fjct_compute_jsd_expert <- function(neighborhood_residuals_S1, neighborhood_residuals_S2, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range) {
    neighborhood_residuals_S1 <- .tox_as_double_array(neighborhood_residuals_S1, "neighborhood_residuals_S1", 3L)
    neighborhood_residuals_S2 <- .tox_as_double_array(neighborhood_residuals_S2, "neighborhood_residuals_S2", 3L)
    neighbor_mask_S1 <- .tox_as_logical(neighbor_mask_S1, "neighbor_mask_S1")
    neighbor_mask_S2 <- .tox_as_logical(neighbor_mask_S2, "neighbor_mask_S2")
    n_bins <- .tox_as_integer_scalar(n_bins, "n_bins")
    shared_residual_range <- .tox_as_double_scalar(shared_residual_range, "shared_residual_range")
    if (dim(neighborhood_residuals_S2)[2] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighbor_mask_S1)[1] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighbor_mask_S1", dim(neighbor_mask_S1)[1], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighbor_mask_S2)[1] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighbor_mask_S2", dim(neighbor_mask_S2)[1], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_residuals_S2)[3] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[3], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])
    if (dim(neighbor_mask_S1)[2] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighbor_mask_S1", dim(neighbor_mask_S1)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])
    if (dim(neighbor_mask_S2)[2] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighbor_mask_S2", dim(neighbor_mask_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])

    .result <- .Call("fjct_compute_jsd_expert_call", neighborhood_residuals_S1, neighborhood_residuals_S2, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range)
    .arguments <- c("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "neighbor_mask_S1", "neighbor_mask_S2", "n_bins", "shared_residual_range", "js_divergences", "included_n_reps_S1", "included_n_reps_S2", "total_included_n_reps", "global_js_divergence", "weights", "pmf_S1", "pmf_S2", "tmp_counts", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        js_divergences = .result$js_divergences,
        included_n_reps_S1 = .result$included_n_reps_S1,
        included_n_reps_S2 = .result$included_n_reps_S2,
        total_included_n_reps = .result$total_included_n_reps,
        global_js_divergence = .result$global_js_divergence,
        weights = .result$weights,
        pmf_S1 = .result$pmf_S1,
        pmf_S2 = .result$pmf_S2
    )
}

#' Computes the per-family/per-sub-neighborhood contribution score
#'
#' Combines (1) how divergent the family is between the studies, and (2) how much residual support the family has overall, using the outputs from [[tox_data_integration(module):fjct_compute_jsd(interface)]], collected for the analyzed sub-neighborhoods.
#'
#' @param global_js_divergences a numeric vector. Per-sub-neighborhood weighted global JSD
#' @param total_included_n_reps_per_f a integer vector. Per-sub-neighborhood `total_included_n_reps`
#' @return a named list with elements `support_weights`, `contribution_scores`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration::fjct_compute_contribution_scores}.
#' @export
fjct_compute_contribution_scores <- function(global_js_divergences, total_included_n_reps_per_f) {
    global_js_divergences <- .tox_as_double_vector(global_js_divergences, "global_js_divergences")
    total_included_n_reps_per_f <- .tox_as_integer_vector(total_included_n_reps_per_f, "total_included_n_reps_per_f")
    if (length(total_included_n_reps_per_f) != length(global_js_divergences))
        .tox_shape_error("total_included_n_reps_per_f", length(total_included_n_reps_per_f), "global_js_divergences", length(global_js_divergences))

    .result <- .Call("fjct_compute_contribution_scores_call", global_js_divergences, total_included_n_reps_per_f)
    .arguments <- c("global_js_divergences", "total_included_n_reps_per_f", "k_families", "support_weights", "contribution_scores", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        support_weights = .result$support_weights,
        contribution_scores = .result$contribution_scores
    )
}
