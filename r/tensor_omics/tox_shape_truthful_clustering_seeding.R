# Generated. Do not edit.

#' Per-vector local density label, an adaptive-bandwidth kernel density estimate over each vector's own k_density nearest neighbors
#'
#' For each vector: find its `k_density` nearest neighbors (excluding itself), take the
#' `bandwidth_percentile` percentile of the distances to them as a per-vector local
#' bandwidth, then sum a Gaussian kernel over those same distances at that bandwidth,
#' normalized by `bandwidth**n_dimensions`. Unlike a single dataset-wide radius, this
#' bandwidth shrinks in dense regions and grows in sparse ones, so the resulting labels
#' reflect local, not global, density. See `misc/mod_STC.md`, SKG `density_labels`.
#'
#' `bandwidth_percentile` is a heuristic knob, not a calibrated standard deviation, and
#' deliberately documented as one: for a genuine 1-D Gaussian, 68.27% of its mass sits
#' within one SD, so the 68.27th percentile of *samples from that Gaussian* equals the SD
#' exactly -- which is where the default comes from -- but our distances are norms in
#' `n_dimensions` dimensions, not draws from a 1-D Gaussian, and the same "percentile that
#' equals the SD" shifts with dimension (it is ~39% at 2 dimensions, ~20% at 3, following
#' a chi distribution with `n_dimensions` degrees of freedom, not a plain half-normal).
#' Correcting for that would need the *local intrinsic* dimension, not the ambient one --
#' STC's whole premise is data concentrated near a lower-dimensional manifold, so the
#' ambient `n_dimensions` is typically the wrong dimension to correct with anyway, and the
#' local one is not yet known at this point in the pipeline (estimating it is `observable`'s
#' job, run later, on an actual ensemble -- not something to redo per point just to
#' calibrate a seeding bandwidth). So this is left an explicit, undisguised heuristic:
#' `bandwidth_percentile` is exposed for exactly this reason -- to be explored
#' empirically (see `misc/STC-experiments/README.md`) rather than settled by further
#' first-principles argument that does not actually resolve for manifold-concentrated data.
#'
#' The `bandwidth**n_dimensions` normalization is not optional, independent of which
#' percentile is chosen: a raw `sum(exp(-d**2/(2*bandwidth**2)))`, without it, is
#' scale-invariant -- scaling every distance (and therefore the bandwidth) by the same
#' constant leaves `d/bandwidth`, the only thing that enters the exponent, unchanged, so a
#' tight cluster and the same cluster stretched out 20x would score identically. The
#' standard adaptive-KDE normalization (divide by the bandwidth to the power of the
#' ambient dimension) ties the estimate back to an absolute scale, so a genuinely tighter
#' neighborhood outscores a genuinely looser one, not just a differently-shaped one.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_seeding::density_labels_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param kd_indices a integer vector. Pre-built k-d tree index over `vectors`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param k_density a integer scalar. Neighborhood size the local density estimate is taken over
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors - 1`.
#'   The default value is `30`.
#' @param bandwidth_percentile a numeric scalar. Percentile (0 to 100) of the k_density neighbor distances used as the local
#'   Gaussian bandwidth -- a heuristic choice, not a calibrated standard deviation,
#'   see above
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `68.27`.
#' @return a numeric vector. Per-vector local density label
#' @export
density_labels <- function(vectors, kd_indices, dimension_order, k_density = 30L, bandwidth_percentile = 68.27) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    k_density <- .tox_as_integer_scalar(k_density, "k_density")
    bandwidth_percentile <- .tox_as_double_scalar(bandwidth_percentile, "bandwidth_percentile")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "vectors", dim(vectors)[1])
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "vectors", dim(vectors)[2])

    .result <- .Call("density_labels_call", vectors, kd_indices, dimension_order, k_density, bandwidth_percentile)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "k_density", "bandwidth_percentile", "labels", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$labels
}

#' Select seed points via greedy, density-ranked, coverage-based selection
#'
#' Ranks vectors by density label, descending (see `density_labels`). Starting with the
#' highest-density unvisited vector, marks it a seed, marks every vector within its own
#' coverage radius as visited, and continues with the next-highest-density unvisited
#' vector until none remain -- so only genuinely uncovered regions can seed another
#' ensemble. The coverage radius is
#' \code{\link{calc_ensemble_growth_radius}}'s
#' own computation, called on the newly-selected seed with `k_density` in place of
#' `k_min` -- not a separate, dataset-wide radius: a fixed global radius can suppress
#' seed placement across a region much larger than what that seed's own ensemble will
#' ever actually grow into, leaving points "covered" by seed-exclusion but never reached
#' by any grown ensemble, see `misc/STC-experiments/README.md`.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_seeding::seeds_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param kd_indices a integer vector. Pre-built k-d tree index over `vectors`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param k_density a integer scalar. Neighborhood size for both the density estimate and the coverage radius, see
#'   `density_labels` and `calc_ensemble_growth_radius`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors - 1`.
#'   The default value is `30`.
#' @param bandwidth_percentile a numeric scalar. Percentile (0 to 100) of the k_density neighbor distances used as the local
#'   Gaussian bandwidth, see `density_labels`
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `68.27`.
#' @return a logical vector. TRUE for points selected as seeds
#' @export
seeds <- function(vectors, kd_indices, dimension_order, k_density = 30L, bandwidth_percentile = 68.27) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    k_density <- .tox_as_integer_scalar(k_density, "k_density")
    bandwidth_percentile <- .tox_as_double_scalar(bandwidth_percentile, "bandwidth_percentile")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "vectors", dim(vectors)[1])
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "vectors", dim(vectors)[2])

    .result <- .Call("seeds_call", vectors, kd_indices, dimension_order, k_density, bandwidth_percentile)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "k_density", "bandwidth_percentile", "is_seed_mask", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$is_seed_mask
}
