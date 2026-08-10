# Generated. Do not edit.

#' Pick n_anchors point indices at evenly-spaced percentiles of the density-sorted order
#'
#' Nearest-rank selection (not `calc_percentile_helper`'s own linear interpolation): these
#' must be genuine point indices, not interpolated values. Percentiles are
#' $100/n_{\text{anchors}}, 200/n_{\text{anchors}}, \ldots, 100$ -- e.g. $n_{\text{anchors}}=5$
#' gives 20/40/60/80/100%ile. Duplicate anchor indices are possible (not deduplicated) when
#' `n_anchors` is close to `n_vectors`; harmless -- a duplicated anchor's cloud just grows
#' redundantly, one estimator-anchor "slot" among several effectively wasted, not incorrect.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_parameter_estimation::sample_estimator_anchors_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param density_labels a numeric vector. Per-vector density label, see density_labels
#' @param n_anchors a integer scalar. Number of estimator anchors (EAs) to pick
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @return a integer vector. Point indices of the n_anchors estimator anchors, ascending-percentile order
#' @export
sample_estimator_anchors <- function(density_labels, n_anchors) {
    density_labels <- .tox_as_double_vector(density_labels, "density_labels")
    n_anchors <- .tox_as_integer_scalar(n_anchors, "n_anchors")
    .result <- .Call("sample_estimator_anchors_call", density_labels, n_anchors)
    .arguments <- c("density_labels", "n_vectors", "n_anchors", "anchor_indices", "ierr")
    .sources <- c(NA_character_, "density_labels", "anchor_indices", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$anchor_indices
}

#' Multi-source competitive region growth of the n_anchors estimator anchors, bounded by seed_max_set_size
#'
#' Each EA starts as its own single-point cloud (its anchor). Every round, the globally
#' closest (unclaimed point, proposing EA) pair -- closest meaning nearest to *any* member
#' of that EA's own current cloud, not just to the anchor -- is claimed, until either no
#' unclaimed point remains reachable or the total claimed across all EAs reaches
#' `seed_max_set_size` percent of `n_vectors`. Implemented as a brute-force rescan every
#' round directly on `vectors`, not via the k-d tree: `f42_kd_tree` has no "nearest
#' unclaimed point to a growing region" primitive, and with `n_anchors` small and total
#' growth capped small by design, the rescan is cheap in absolute terms regardless -- see
#' `misc/mod_STC.md` for the full complexity discussion and why this is a deliberate
#' simplicity trade, not an oversight.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_parameter_estimation::grow_estimator_anchor_clouds}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param anchor_indices a integer vector. Point indices of the estimator anchors, see sample_estimator_anchors
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param seed_max_set_size a numeric scalar. Percent (0 to 100) of n_vectors at which total growth across all EAs stops
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `5.0`.
#' @return a named list with elements:
#'   \item{cloud_masks}{a logical matrix. TRUE for members of each EA's (column) final cloud, including its own anchor}
#'   \item{cloud_sizes}{a integer vector. Final cloud size per EA}
#' @export
grow_estimator_anchor_clouds <- function(vectors, anchor_indices, seed_max_set_size = 5.0) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    anchor_indices <- .tox_as_integer_vector(anchor_indices, "anchor_indices")
    seed_max_set_size <- .tox_as_double_scalar(seed_max_set_size, "seed_max_set_size")
    .result <- .Call("grow_estimator_anchor_clouds_call", vectors, anchor_indices, seed_max_set_size)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "anchor_indices", "n_anchors", "seed_max_set_size", "cloud_masks", "cloud_sizes", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, "anchor_indices", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        cloud_masks = .result$cloud_masks,
        cloud_sizes = .result$cloud_sizes
    )
}

#' Estimate k_min, k_density, density_quantile, chordal_dist_max_as_prcnt_of_range, G_max, d_max from the data
#'
#' Orchestrates density_labels -> sample_estimator_anchors -> grow_estimator_anchor_clouds
#' -> observable (once per EA) -> pairwise EA comparisons -> aggregation. See
#' `misc/mod_STC.md`, SKG `estimate_stc_parameters`, for the full definition of every
#' output and the reasoning behind each. EAs whose final cloud has fewer than 2 members
#' (no meaningful SVD possible -- a documented, deliberately unguarded-against possibility
#' of `grow_estimator_anchor_clouds`'s own stop condition, see there) are excluded from
#' every statistic below; `ierr` is set if fewer than 2 EAs remain usable (no pairwise
#' comparison -- and therefore no G_max/d_max estimate -- is possible at all), or if every
#' usable pair has a zero shared rank (no chordal-distance estimate possible either, even
#' though G_max/d_max still are). Both are the one genuine, input-shape-dependent runtime
#' failure this SKG can hit that no simple per-argument DM_* annotation could foresee (it
#' depends on the data's own spatial distribution, not just
#' n_vectors/n_anchors/seed_max_set_size in isolation) -- see `codegen_guide.md` section
#' 5.14.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_parameter_estimation::estimate_stc_parameters_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param kd_indices a integer vector. Pre-built k-d tree index over `vectors`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param k_density a integer scalar. Passed through to density_labels
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors - 1`.
#' @param bandwidth_percentile a numeric scalar. Passed through to density_labels
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#' @param n_anchors a integer scalar. Number of estimator anchors (EAs), see sample_estimator_anchors
#'   The minimum valid value is `2`.
#'   The maximum valid value is `n_vectors`.
#'   The default value is `5`.
#' @param seed_max_set_size a numeric scalar. Passed through to grow_estimator_anchor_clouds
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `5.0`.
#' @param first_quartile_percentile a numeric scalar. Percentile (0 to 100) of the pairwise-EA-comparison distributions used for
#'   chordal_dist_max_as_prcnt_of_range/G_max/d_max, see estimate_stc_parameters
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `25.0`.
#' @return a named list with elements:
#'   \item{estimated_k_min}{a numeric scalar. Estimated k_min (real-valued; round for direct use as an integer argument)}
#'   \item{estimated_k_density}{a numeric scalar. Estimated k_density (equal to estimated_k_min, see estimate_stc_parameters)}
#'   \item{estimated_density_quantile}{a numeric scalar. Estimated density_quantile -- a literal radius (data units), not a percentile}
#'   \item{estimated_chordal_dist_max_as_prcnt_of_range}{a numeric scalar. Estimated chordal_dist_max_as_prcnt_of_range (0 to 1)}
#'   \item{estimated_G_max}{a numeric scalar. Estimated G_max}
#'   \item{estimated_d_max}{a numeric scalar. Estimated d_max (real-valued; round for direct use as an integer argument)}
#' @export
estimate_stc_parameters <- function(vectors, kd_indices, dimension_order, k_density = NULL, bandwidth_percentile = NULL, n_anchors = 5L, seed_max_set_size = 5.0, first_quartile_percentile = 25.0) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    if (!is.null(k_density))
        k_density <- .tox_as_integer_scalar(k_density, "k_density")
    if (!is.null(bandwidth_percentile))
        bandwidth_percentile <- .tox_as_double_scalar(bandwidth_percentile, "bandwidth_percentile")
    n_anchors <- .tox_as_integer_scalar(n_anchors, "n_anchors")
    seed_max_set_size <- .tox_as_double_scalar(seed_max_set_size, "seed_max_set_size")
    first_quartile_percentile <- .tox_as_double_scalar(first_quartile_percentile, "first_quartile_percentile")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "vectors", dim(vectors)[1])
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "vectors", dim(vectors)[2])

    .result <- .Call("estimate_stc_parameters_call", vectors, kd_indices, dimension_order, k_density, bandwidth_percentile, n_anchors, seed_max_set_size, first_quartile_percentile)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "k_density", "bandwidth_percentile", "n_anchors", "seed_max_set_size", "first_quartile_percentile", "estimated_k_min", "estimated_k_density", "estimated_density_quantile", "estimated_chordal_dist_max_as_prcnt_of_range", "estimated_G_max", "estimated_d_max", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        estimated_k_min = .result$estimated_k_min,
        estimated_k_density = .result$estimated_k_density,
        estimated_density_quantile = .result$estimated_density_quantile,
        estimated_chordal_dist_max_as_prcnt_of_range = .result$estimated_chordal_dist_max_as_prcnt_of_range,
        estimated_G_max = .result$estimated_G_max,
        estimated_d_max = .result$estimated_d_max
    )
}
