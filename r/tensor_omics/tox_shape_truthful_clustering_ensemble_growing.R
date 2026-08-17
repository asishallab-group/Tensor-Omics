# Generated. Do not edit.

#' Locally adapted ensemble growth radius, a percentile of the distances among a seed's own k_min nearest neighbors
#'
#' Matches LoManLe's `local_scale_i` exactly at the default `radius_percentile=50.0` (the
#' median) -- a per-seed, locally adaptive radius rather than a single dataset-wide one
#' (see `misc/STC_for_LoManLe.md` section 2.2). Work arrays are sized for the worst case
#' (`k_min = n_vectors - 1`) and sliced internally, since `k_min`'s resolved value is only
#' known once its default (if any) has been applied.
#'
#' `radius_percentile` generalizes what used to be a hardcoded median: `seeds_impl`
#' reuses this same SKG for its seed-exclusion radius (see `misc/mod_STC.md`, SKG
#' `seeds`), and a fixed dataset-wide-in-spirit median there was observed to suppress
#' entire uncovered regions (e.g. curvature extrema on a wavy manifold) whose own growth
#' never actually reaches that far -- exposing the percentile lets that specific call site
#' shrink its exclusion radius independently of this kernel's own default, without
#' touching the actual growth-radius computation any other caller relies on.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_ensemble_growing::calc_ensemble_growth_radius}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param kd_indices a integer vector. Pre-built k-d tree index over `vectors`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param seed_index a integer scalar. Index into `vectors`/`kd_indices` of the seed to compute the growth radius for
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param k_min a integer scalar. Neighborhood size the median distance is taken over
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors - 1`.
#'   The default value is `30`.
#' @param radius_percentile a numeric scalar. Percentile (0 to 100) of the k_min neighbor distances reported as the growth
#'   radius -- 50.0 (the default) is the median, matching this SKG's original,
#'   non-parameterized behavior
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `50.0`.
#' @return a numeric scalar. radius_percentile-th percentile of the distances among the seed's own k_min
#'   nearest neighbors
#' @export
calc_ensemble_growth_radius <- function(vectors, kd_indices, dimension_order, seed_index, k_min = 30L, radius_percentile = 50.0) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    seed_index <- .tox_as_integer_scalar(seed_index, "seed_index")
    k_min <- .tox_as_integer_scalar(k_min, "k_min")
    radius_percentile <- .tox_as_double_scalar(radius_percentile, "radius_percentile")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "vectors", dim(vectors)[1])
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "vectors", dim(vectors)[2])

    .result <- .Call("calc_ensemble_growth_radius_call", vectors, kd_indices, dimension_order, seed_index, k_min, radius_percentile)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "seed_index", "k_min", "radius_percentile", "growth_radius", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$growth_radius
}

#' Grow an ensemble by one step, the union of every current member's growth-radius neighborhood
#'
#' $\mathcal{C}_{\mathcal{E}_{t+1}} = \{x_k \mid \|x_k-x_i\|\le r_{\mathcal{E}} \;\exists\,
#' x_i\in\mathcal{E}\}$. Deliberately a plain sequential loop, not parallelized: this runs
#' once per growth iteration per ensemble, and outer-level parallelism across
#' ensembles/seeds (in `ensemble_identification`) is the right place for that -- a single
#' ensemble's member count is typically small, especially early in growth. An all-FALSE
#' `is_member_mask` (an empty ensemble) is a well-defined degenerate case: there is nothing
#' to grow from, so the result is all-FALSE too, not an error.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_ensemble_growing::grow_ensemble}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param kd_indices a integer vector. Pre-built k-d tree index over `vectors`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param is_member_mask a logical vector. Current ensemble membership
#' @param growth_radius a numeric scalar. This ensemble's growth radius, see `calc_ensemble_growth_radius`
#'   The minimum valid value is `0.0`.
#' @return a logical vector. Grown ensemble membership (superset of `is_member_mask`)
#' @export
grow_ensemble <- function(vectors, kd_indices, dimension_order, is_member_mask, growth_radius) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    is_member_mask <- .tox_as_logical(is_member_mask, "is_member_mask")
    growth_radius <- .tox_as_double_scalar(growth_radius, "growth_radius")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "vectors", dim(vectors)[1])
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "vectors", dim(vectors)[2])
    if (length(is_member_mask) != dim(vectors)[2])
        .tox_shape_error("is_member_mask", length(is_member_mask), "vectors", dim(vectors)[2])

    .result <- .Call("grow_ensemble_call", vectors, kd_indices, dimension_order, is_member_mask, growth_radius)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "is_member_mask", "growth_radius", "is_member_mask_next", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$is_member_mask_next
}
