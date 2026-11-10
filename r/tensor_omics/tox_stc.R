# Generated. Do not edit.

#' Compute the density-label radius from distances to the data centroid.
#'
#' Generated from the Fortran procedure \code{tox_stc::calculate_density_radius}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix (n_dimensions x n_vectors)
#' @param mean_to_other_vecs_dist_quant a numeric scalar. Optional quantile fraction (0.0 to 1.0), defaults to 0.15
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.15`.
#' @return a numeric scalar. The resulting density search radius
#' @export
calculate_density_radius <- function(vectors, mean_to_other_vecs_dist_quant = 0.15) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    mean_to_other_vecs_dist_quant <- .tox_as_double_scalar(mean_to_other_vecs_dist_quant, "mean_to_other_vecs_dist_quant")
    .result <- .Call("calculate_density_radius_call", vectors, mean_to_other_vecs_dist_quant)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "mean_to_other_vecs_dist_quant", "radius", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$radius
}

#' Count points in each Euclidean neighborhood to label local density.
#'
#' Generated from the Fortran procedure \code{tox_stc::calculate_labels_as_density}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix (n_dimensions x n_vectors)
#' @param r a numeric scalar. Label-sphere search radius
#'   The minimum valid value is `0.0`.
#' @param dimension_order a integer vector. Sequence array tracking tree split axes by variance
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param kd_indices a integer vector. KD-tree index sequence array computed via \code{\link{build_kd_index}}.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param n_tiles a integer scalar. Number of concurrent vector tiles, at least 1. Capped at the number of points.
#'   The minimum valid value is `1`.
#' @return a numeric vector. Output array storing generated density scalars
#' @export
calculate_labels_as_density <- function(vectors, r, dimension_order, kd_indices, n_tiles) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    r <- .tox_as_double_scalar(r, "r")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    n_tiles <- .tox_as_integer_scalar(n_tiles, "n_tiles")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", 1L, length(dimension_order), "'vectors' implies n_dimensions == ", dim(vectors)[1], "n_dimensions")
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", 1L, length(kd_indices), "'vectors' implies n_vectors == ", dim(vectors)[2], "n_vectors")

    .result <- .Call("calculate_labels_as_density_call", vectors, r, dimension_order, kd_indices, n_tiles)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "r", "dimension_order", "kd_indices", "n_tiles", "label_densities", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$label_densities
}

#' Select seeds in descending density order and compute their local radii.
#'
#' Generated from the Fortran procedure \code{tox_stc::identify_ensemble_seeds}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix (n_dimensions x n_vectors)
#' @param density_labels a numeric vector. Precalculated density labels for all ambient vectors
#'   The minimum valid value is `1.0`.
#' @param dimension_order a integer vector. Dimension split order array tracking the K-D tree structure
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param kd_indices a integer vector. K-D tree index sequence array
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param k_seeding a integer scalar. Number of nearest neighbors used to estimate each seed's local coverage radius
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors - 1`.
#' @param coverage_quant a numeric scalar. Optional percentile of k-nearest-neighbor distances used as the seed coverage and growth radius; defaults to 0.5.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.5`.
#' @return a named list with elements:
#'   \item{sorted_perm}{a integer vector. Vector indices sorted by density descending}
#'   \item{n_seeds}{a integer scalar. Number of selected seed vectors}
#'   \item{seed_mask}{a logical vector. Logical mask identifying selected seed vectors}
#'   \item{seed_radii}{a numeric vector. Per-seed growth radius, stored at each seed's own vector index; zero at non-seeds}
#' @export
identify_ensemble_seeds <- function(vectors, density_labels, dimension_order, kd_indices, k_seeding, coverage_quant = 0.5) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    density_labels <- .tox_as_double_vector(density_labels, "density_labels")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    k_seeding <- .tox_as_integer_scalar(k_seeding, "k_seeding")
    coverage_quant <- .tox_as_double_scalar(coverage_quant, "coverage_quant")
    if (length(density_labels) != dim(vectors)[2])
        .tox_shape_error("density_labels", 1L, length(density_labels), "'vectors' implies n_vectors == ", dim(vectors)[2], "n_vectors")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", 1L, length(dimension_order), "'vectors' implies n_dimensions == ", dim(vectors)[1], "n_dimensions")
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", 1L, length(kd_indices), "'vectors' implies n_vectors == ", dim(vectors)[2], "n_vectors")

    .result <- .Call("identify_ensemble_seeds_call", vectors, density_labels, dimension_order, kd_indices, k_seeding, coverage_quant)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "density_labels", "dimension_order", "kd_indices", "k_seeding", "sorted_perm", "n_seeds", "seed_mask", "seed_radii", "coverage_quant", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        sorted_perm = .result$sorted_perm,
        n_seeds = .result$n_seeds,
        seed_mask = .result$seed_mask,
        seed_radii = .result$seed_radii
    )
}

#' Compute the median and median absolute deviation of the density labels.
#'
#' Generated from the Fortran procedure \code{tox_stc::compute_ambient_density_stats}, whose argument names
#' are the ones an error message reports.
#'
#' @param density_labels a numeric vector. Density labels for all vectors the density-label routine.
#'   The minimum valid value is `1.0`.
#' @return a named list with elements:
#'   \item{median_ambient}{a numeric scalar. Ambient median density median_{i=1..N}(rho_i)}
#'   \item{mad_ambient}{a numeric scalar. Ambient median absolute deviation median_{i=1..N}(|rho_i - median_ambient|)}
#' @export
compute_ambient_density_stats <- function(density_labels) {
    density_labels <- .tox_as_double_vector(density_labels, "density_labels")
    .result <- .Call("compute_ambient_density_stats_call", density_labels)
    .arguments <- c("density_labels", "n_vectors", "median_ambient", "mad_ambient", "ierr")
    .sources <- c(NA_character_, "density_labels", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        median_ambient = .result$median_ambient,
        mad_ambient = .result$mad_ambient
    )
}

#' Merge intersecting ensembles by overlap coefficient and discard final singletons.
#'
#' Generated from the Fortran procedure \code{tox_stc::merge_ensembles}, whose argument names
#' are the ones an error message reports.
#'
#' @param raw_masks a logical matrix. Matrix of raw unmerged boolean ensemble masks
#' @param min_overlap_coefficient a numeric scalar. Minimum Overlap Coefficient required to merge two ensembles
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @return a logical matrix. Output matrix storing merged boolean ensemble masks
#'   The first `n_ensembles` elements will hold the results.
#' @export
merge_ensembles <- function(raw_masks, min_overlap_coefficient = 0.9) {
    raw_masks <- .tox_as_logical_matrix(raw_masks, "raw_masks")
    min_overlap_coefficient <- .tox_as_double_scalar(min_overlap_coefficient, "min_overlap_coefficient")
    .result <- .Call("merge_ensembles_call", raw_masks, min_overlap_coefficient)
    .arguments <- c("raw_masks", "n_vectors", "n_seeds", "min_overlap_coefficient", "merged_masks", "n_ensembles", "ierr")
    .sources <- c(NA_character_, "raw_masks", "raw_masks", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$merged_masks[, seq_len(.result$n_ensembles), drop = FALSE]
}

#' Grow one ensemble by one density-compatible surface step.
#'
#' Generated from the Fortran procedure \code{tox_stc::grow_ensemble}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix (n_dimensions x n_vectors)
#' @param ensemble_mask a logical vector. Membership before the growth step.
#' @param r a numeric scalar. Search radius threshold for surface growth
#'   The minimum valid value is `0.0`.
#' @param dimension_order a integer vector. Sequence array tracking split dimensions
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param kd_indices a integer vector. KD-tree index sequence array computed via \code{\link{build_kd_index}}.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param density_labels a numeric vector. Density labels for all vectors the density-label routine.
#'   The minimum valid value is `1.0`.
#' @param alpha_mad a numeric scalar. Multiplier factor for MAD density compatibility threshold
#'   The minimum valid value is `0.0`.
#'   The default value is `0.5`.
#' @return a logical vector. Membership after the growth step.
#' @export
grow_ensemble <- function(vectors, ensemble_mask, r, dimension_order, kd_indices, density_labels, alpha_mad = 0.5) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    ensemble_mask <- .tox_as_logical_vector(ensemble_mask, "ensemble_mask")
    r <- .tox_as_double_scalar(r, "r")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    density_labels <- .tox_as_double_vector(density_labels, "density_labels")
    alpha_mad <- .tox_as_double_scalar(alpha_mad, "alpha_mad")
    if (length(ensemble_mask) != dim(vectors)[2])
        .tox_shape_error("ensemble_mask", 1L, length(ensemble_mask), "'vectors' implies n_vectors == ", dim(vectors)[2], "n_vectors")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", 1L, length(dimension_order), "'vectors' implies n_dimensions == ", dim(vectors)[1], "n_dimensions")
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", 1L, length(kd_indices), "'vectors' implies n_vectors == ", dim(vectors)[2], "n_vectors")
    if (length(density_labels) != dim(vectors)[2])
        .tox_shape_error("density_labels", 1L, length(density_labels), "'vectors' implies n_vectors == ", dim(vectors)[2], "n_vectors")

    .result <- .Call("grow_ensemble_call", vectors, ensemble_mask, r, dimension_order, kd_indices, density_labels, alpha_mad)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "ensemble_mask", "r", "dimension_order", "kd_indices", "density_labels", "alpha_mad", "grown_mask", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$grown_mask
}

#' Update the five ensemble observables for a growth iteration.
#'
#' Generated from the Fortran procedure \code{tox_stc::compute_ensemble_observable}, whose argument names
#' are the ones an error message reports.
#'
#' @param ensemble_mask a logical vector. Selection mask for active ensemble members
#' @param density_labels a numeric vector. Density labels for all vectors the density-label routine.
#'   The minimum valid value is `1.0`.
#' @param candidate_count a integer scalar. Number of candidates added in current growth step
#'   The minimum valid value is `0`.
#' @param current_iter a integer scalar. Current growth iteration index (1-based)
#'   The minimum valid value is `1`.
#' @param observables a numeric matrix. 2D observable trajectory matrix where columns match iteration indices
#' @param t_observables a integer scalar. History tracking depth (defaults to 10; <= 0 represents infinity / all)
#'   The default value is `10`.
#' @return a numeric matrix. 2D observable trajectory matrix where columns match iteration indices
#' @export
compute_ensemble_observable <- function(ensemble_mask, density_labels, candidate_count, current_iter, observables, t_observables = 10L) {
    ensemble_mask <- .tox_as_logical_vector(ensemble_mask, "ensemble_mask")
    density_labels <- .tox_as_double_vector(density_labels, "density_labels")
    candidate_count <- .tox_as_integer_scalar(candidate_count, "candidate_count")
    current_iter <- .tox_as_integer_scalar(current_iter, "current_iter")
    observables <- .tox_as_double_matrix(observables, "observables")
    t_observables <- .tox_as_integer_scalar(t_observables, "t_observables")
    if (length(density_labels) != length(ensemble_mask))
        .tox_shape_error("density_labels", 1L, length(density_labels), "'ensemble_mask' implies n_vectors == ", length(ensemble_mask), "n_vectors")

    .result <- .Call("compute_ensemble_observable_call", ensemble_mask, density_labels, candidate_count, current_iter, observables, t_observables)
    .arguments <- c("ensemble_mask", "density_labels", "n_vectors", "candidate_count", "current_iter", "observables", "n_history", "n_observables", "t_observables", "ierr")
    .sources <- c(NA_character_, NA_character_, "ensemble_mask", NA_character_, NA_character_, NA_character_, "observables", "observables", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$observables
}

#' Test whether the change in ensemble observables permits a growth step.
#'
#' Generated from the Fortran procedure \code{tox_stc::accept_ensemble}, whose argument names
#' are the ones an error message reports.
#'
#' @param observables a numeric matrix. 2D observables trajectory matrix [5 x window_capacity]
#' @param current_iter a integer scalar. Current growth iteration index (1-based)
#'   The minimum valid value is `1`.
#' @param alpha_accept a numeric scalar. Optional log2 fold change acceptance threshold (defaults to 0.5)
#'   The minimum valid value is `0.0`.
#'   The default value is `0.5`.
#' @return a logical scalar. Output flag indicating step acceptance decision
#' @export
accept_ensemble <- function(observables, current_iter, alpha_accept = 0.5) {
    observables <- .tox_as_double_matrix(observables, "observables")
    current_iter <- .tox_as_integer_scalar(current_iter, "current_iter")
    alpha_accept <- .tox_as_double_scalar(alpha_accept, "alpha_accept")
    .result <- .Call("accept_ensemble_call", observables, current_iter, alpha_accept)
    .arguments <- c("observables", "n_history", "n_observables", "current_iter", "alpha_accept", "is_accepted", "ierr")
    .sources <- c(NA_character_, "observables", "observables", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$is_accepted
}

#' Grow one raw ensemble from each seed, using a global radius or per-seed radii.
#'
#' Generated from the Fortran procedure \code{tox_stc::obtain_ensembles_global_radius}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix (n_dimensions x n_vectors)
#' @param dimension_order a integer vector. Sequence array tracking split dimensions for KD-tree
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param kd_indices a integer vector. KD-tree index sequence array computed via \code{\link{build_kd_index}}.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param density_labels a numeric vector. Precalculated density labels for ambient vectors the density-label routine.
#'   The minimum valid value is `1.0`.
#' @param seed_indices a integer vector. Array of starting seed vector indices the seed-selection routine.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param r a numeric scalar. Uniform growth radius for every seed.
#'   The minimum valid value is `0.0`.
#' @param alpha_mad a numeric scalar. MAD multiplier for density compatibility threshold
#'   The minimum valid value is `0.0`.
#'   The default value is `0.5`.
#' @param alpha_accept a numeric scalar. Acceptance log2 fold-change threshold
#'   The minimum valid value is `0.0`.
#'   The default value is `0.5`.
#' @param t_observables a integer scalar. History window in iterations. Pass n_vectors to retain the full history.
#'   The minimum valid value is `1`.
#' @param n_tiles a integer scalar. Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
#'   The minimum valid value is `1`.
#' @return a named list with elements:
#'   \item{ensemble_matrix}{a logical matrix. Output matrix storing unmerged grown ensemble masks}
#'   \item{stop_reasons}{a integer vector. Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
#'     3: first step rejected, 4: no candidates.}
#'   \item{mad_ambient}{a numeric scalar. Ambient MAD of density labels; zero makes alpha_mad ineffective by allowing only exactly equal densities.}
#'   \item{n_ensembles}{a integer scalar. Total count of extracted raw ensembles (equals n_seeds)}
#' @export
obtain_ensembles_global_radius <- function(vectors, dimension_order, kd_indices, density_labels, seed_indices, r, alpha_mad = 0.5, alpha_accept = 0.5, t_observables, n_tiles) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    density_labels <- .tox_as_double_vector(density_labels, "density_labels")
    seed_indices <- .tox_as_integer_vector(seed_indices, "seed_indices")
    r <- .tox_as_double_scalar(r, "r")
    alpha_mad <- .tox_as_double_scalar(alpha_mad, "alpha_mad")
    alpha_accept <- .tox_as_double_scalar(alpha_accept, "alpha_accept")
    t_observables <- .tox_as_integer_scalar(t_observables, "t_observables")
    n_tiles <- .tox_as_integer_scalar(n_tiles, "n_tiles")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", 1L, length(dimension_order), "'vectors' implies n_dimensions == ", dim(vectors)[1], "n_dimensions")
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", 1L, length(kd_indices), "'vectors' implies n_vectors == ", dim(vectors)[2], "n_vectors")
    if (length(density_labels) != dim(vectors)[2])
        .tox_shape_error("density_labels", 1L, length(density_labels), "'vectors' implies n_vectors == ", dim(vectors)[2], "n_vectors")

    .result <- .Call("obtain_ensembles_global_radius_call", vectors, dimension_order, kd_indices, density_labels, seed_indices, r, alpha_mad, alpha_accept, t_observables, n_tiles)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "dimension_order", "kd_indices", "density_labels", "seed_indices", "n_seeds", "r", "alpha_mad", "alpha_accept", "t_observables", "n_tiles", "ensemble_matrix", "stop_reasons", "mad_ambient", "n_ensembles", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, "seed_indices", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        ensemble_matrix = .result$ensemble_matrix,
        stop_reasons = .result$stop_reasons,
        mad_ambient = .result$mad_ambient,
        n_ensembles = .result$n_ensembles
    )
}

#' Grow one raw ensemble from each seed, using a global radius or per-seed radii.
#'
#' Generated from the Fortran procedure \code{tox_stc::obtain_ensembles_seed_radius}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix (n_dimensions x n_vectors)
#' @param dimension_order a integer vector. Sequence array tracking split dimensions for KD-tree
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param kd_indices a integer vector. KD-tree index sequence array computed via \code{\link{build_kd_index}}.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param density_labels a numeric vector. Precalculated density labels for ambient vectors the density-label routine.
#'   The minimum valid value is `1.0`.
#' @param seed_indices a integer vector. Array of starting seed vector indices the seed-selection routine.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param seed_radii a numeric vector. Per-seed growth radius, allowing tighter growth in dense regions and wider growth in sparse regions.
#'   The minimum valid value is `0.0`.
#' @param alpha_mad a numeric scalar. MAD multiplier for density compatibility threshold
#'   The minimum valid value is `0.0`.
#'   The default value is `0.5`.
#' @param alpha_accept a numeric scalar. Acceptance log2 fold-change threshold
#'   The minimum valid value is `0.0`.
#'   The default value is `0.5`.
#' @param t_observables a integer scalar. History window in iterations. Pass n_vectors to retain the full history.
#'   The minimum valid value is `1`.
#' @param n_tiles a integer scalar. Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
#'   The minimum valid value is `1`.
#' @return a named list with elements:
#'   \item{ensemble_matrix}{a logical matrix. Output matrix storing unmerged grown ensemble masks}
#'   \item{stop_reasons}{a integer vector. Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
#'     3: first step rejected, 4: no candidates.}
#'   \item{mad_ambient}{a numeric scalar. Ambient MAD of density labels; zero makes alpha_mad ineffective by allowing only exactly equal densities.}
#'   \item{n_ensembles}{a integer scalar. Total count of extracted raw ensembles (equals n_seeds)}
#' @export
obtain_ensembles_seed_radius <- function(vectors, dimension_order, kd_indices, density_labels, seed_indices, seed_radii, alpha_mad = 0.5, alpha_accept = 0.5, t_observables, n_tiles) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    density_labels <- .tox_as_double_vector(density_labels, "density_labels")
    seed_indices <- .tox_as_integer_vector(seed_indices, "seed_indices")
    seed_radii <- .tox_as_double_vector(seed_radii, "seed_radii")
    alpha_mad <- .tox_as_double_scalar(alpha_mad, "alpha_mad")
    alpha_accept <- .tox_as_double_scalar(alpha_accept, "alpha_accept")
    t_observables <- .tox_as_integer_scalar(t_observables, "t_observables")
    n_tiles <- .tox_as_integer_scalar(n_tiles, "n_tiles")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", 1L, length(dimension_order), "'vectors' implies n_dimensions == ", dim(vectors)[1], "n_dimensions")
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", 1L, length(kd_indices), "'vectors' implies n_vectors == ", dim(vectors)[2], "n_vectors")
    if (length(density_labels) != dim(vectors)[2])
        .tox_shape_error("density_labels", 1L, length(density_labels), "'vectors' implies n_vectors == ", dim(vectors)[2], "n_vectors")
    if (length(seed_radii) != length(seed_indices))
        .tox_shape_error("seed_radii", 1L, length(seed_radii), "'seed_indices' implies n_seeds == ", length(seed_indices), "n_seeds")

    .result <- .Call("obtain_ensembles_seed_radius_call", vectors, dimension_order, kd_indices, density_labels, seed_indices, seed_radii, alpha_mad, alpha_accept, t_observables, n_tiles)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "dimension_order", "kd_indices", "density_labels", "seed_indices", "n_seeds", "seed_radii", "alpha_mad", "alpha_accept", "t_observables", "n_tiles", "ensemble_matrix", "stop_reasons", "mad_ambient", "n_ensembles", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, "seed_indices", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        ensemble_matrix = .result$ensemble_matrix,
        stop_reasons = .result$stop_reasons,
        mad_ambient = .result$mad_ambient,
        n_ensembles = .result$n_ensembles
    )
}
