# Generated. Do not edit.

#' Bootstrap an interval for each study's global JSD by resampling the studies' histograms
#'
#' In every bootstrap, each study's histogram at each reference point is drawn again with
#' replacement (multinomial), with as many residuals as the study observed
#' (`included_n_reps`); the consensus pmf, each study's per-point JSD and its weighted global
#' JSD are recomputed from the draw, with the divergence `divergence_mode` selects (by default
#' for exactly two studies their pairwise JSD, otherwise against the consensus). Each study is
#' drawn from its own counts (`counts`): the interval shows how much the observed JSD would move
#' if each study's own residuals were drawn again -- a robustness interval of the observed
#' value. Drawing again adds sampling noise, so the bootstrapped JSDs are larger than the
#' observed one on average: when the studies hardly differ, the interval can lie entirely above
#' the observed JSD.
#'
#' `confidence_interval` is input and output: its incoming `[lower, upper]` values seed every
#' slot of the bottom-k/top-k heaps, so they are only displaced by more extreme bootstrap
#' values. On return it holds `[largest of the n_bootstrapping_top_k_jsds smallest values,
#' smallest of the n_bootstrapping_top_k_jsds largest values]`.
#'
#' Errors (`ierr`), with `confidence_interval` left unchanged: `ERR_INVALID_INPUT`, before any
#' draw, when at some reference point the studies' `included_n_reps` do not add up to
#' `mean_pmf_included_n_reps` (as
#' \code{\link{create_mean_pmf}} produces
#' them), or a study's `counts` at some reference point do not add up to its
#' `included_n_reps`; `ERR_ALLOC_FAIL` when the random number generator cannot be
#' created. The count checks, together with the wrappers' bounds on the counts, establish every
#' draw's preconditions, so no draw can fail later on.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_stats::bootstrap_histogram}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bootstraps a integer scalar. Number of bootstrap resamples to perform
#'   The minimum valid value is `1`.
#' @param counts a integer array of rank 3. Every study's own residual counts per bin: what each study is drawn from; at each
#'   reference point they must add up to the study's `included_n_reps` there
#'   The minimum valid value is `0`.
#' @param mean_pmf_included_n_reps a integer vector. Count of non-NaN replicates per reference point of all studies together (the
#'   consensus's count, the denominator of the weights); must be the sum of the studies'
#'   `included_n_reps` there
#'   The minimum valid value is `0`.
#' @param included_n_reps a integer matrix. Count of non-NaN replicates per reference point, per study: the size of each study's
#'   draw. Every study needs at least one at every reference point
#'   The minimum valid value is `1`.
#' @param confidence_interval a numeric matrix. On entry the seed values of the heaps (see above); on return the bootstrapped
#'   `[lower, upper]` interval per study
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param divergence_mode a string, one of "divergence_symmetric", "divergence_pairwise_k2". How each study's per-point divergence of every bootstrap draw is computed, as in the
#'   routine whose interval this is
#'
#'   The default value is `"divergence_pairwise_k2"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @return a numeric matrix. On entry the seed values of the heaps (see above); on return the bootstrapped
#'   `[lower, upper]` interval per study
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @export
bootstrap_histogram <- function(n_bootstraps, counts, mean_pmf_included_n_reps, included_n_reps, confidence_interval, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, divergence_mode = "divergence_pairwise_k2") {
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    counts <- .tox_as_integer_array(counts, "counts", 3L)
    mean_pmf_included_n_reps <- .tox_as_integer_vector(mean_pmf_included_n_reps, "mean_pmf_included_n_reps")
    included_n_reps <- .tox_as_integer_matrix(included_n_reps, "included_n_reps")
    confidence_interval <- .tox_as_double_matrix(confidence_interval, "confidence_interval")
    two_sided_bootstrapping_significance_level <- .tox_as_double_scalar(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    divergence_mode <- .tox_as_mode(divergence_mode, "divergence_mode", c("divergence_symmetric", "divergence_pairwise_k2"))
    if (length(mean_pmf_included_n_reps) != dim(counts)[1])
        .tox_shape_error("mean_pmf_included_n_reps", 1L, length(mean_pmf_included_n_reps), "'counts' implies n_points == ", dim(counts)[1], "n_points")
    if (dim(included_n_reps)[1] != dim(counts)[1])
        .tox_shape_error("included_n_reps", 1L, dim(included_n_reps)[1], "'counts' implies n_points == ", dim(counts)[1], "n_points, n_studies")
    if (dim(included_n_reps)[2] != dim(counts)[3])
        .tox_shape_error("included_n_reps", 2L, dim(included_n_reps)[2], "'counts' implies n_studies == ", dim(counts)[3], "n_points, n_studies")
    if (dim(confidence_interval)[1] != 2L)
        .tox_shape_error("confidence_interval", 1L, dim(confidence_interval)[1], "its declared shape fixes it at ", 2L, "2, n_studies")
    if (dim(confidence_interval)[2] != dim(counts)[3])
        .tox_shape_error("confidence_interval", 2L, dim(confidence_interval)[2], "'counts' implies n_studies == ", dim(counts)[3], "2, n_studies")

    .result <- .Call("bootstrap_histogram_call", n_bootstraps, counts, mean_pmf_included_n_reps, included_n_reps, confidence_interval, two_sided_bootstrapping_significance_level, random_seed, divergence_mode)
    .arguments <- c("n_bootstraps", "n_bins", "n_points", "n_studies", "counts", "mean_pmf_included_n_reps", "included_n_reps", "confidence_interval", "two_sided_bootstrapping_significance_level", "random_seed", "divergence_mode", "ierr")
    .sources <- c(NA_character_, "counts", "counts", "counts", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$confidence_interval
}
