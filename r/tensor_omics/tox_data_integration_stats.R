# Generated. Do not edit.

#' Estimate how likely each study's observed weighted global JSD is to occur by chance
#'
#' In every permutation, each reference point's pooled consensus counts (`mean_pmf_counts`) are
#' drawn without replacement (multivariate hypergeometric) into the studies, each study taking
#' as many residuals as it observed (`included_n_reps`). The permuted pmfs, their equal-weight
#' consensus (built exactly as
#' \code{\link{create_mean_pmf}} builds the
#' observed one) and each study's weighted global JSD against that consensus are recomputed, and
#' compared with the observed `global_jsd_observed`. `p_values_observed_consensus` comes from the
#' same draws, taking each permuted study's JSD against the fixed observed `mean_pmf` instead.
#'
#' Both p-values are `(count + 1)/(n_permutations + 1)`, so never exactly 0. With
#' `n_permutations = 0` no permutation runs and both are NaN ("not computed").
#'
#' Errors (`ierr`): `ERR_INVALID_INPUT`, before any draw, when at some reference point the
#' studies' counts do not add up to the pool's total (`mean_pmf_included_n_reps`) or that total
#' is not the sum of the pool's bin counts; `ERR_ALLOC_FAIL` when the random number generator
#' cannot be created. The p-values are then NaN. A failed draw later on is recorded (first
#' failure only) without stopping the remaining permutations.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_stats::gjct_permutation_test}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_permutations a integer scalar. Number of permutations; `0` skips the test (both p-value outputs are then NaN)
#'   The minimum valid value is `0`.
#' @param mean_pmf_counts a integer matrix. Pooled residual counts per bin of the consensus pmf: the pool each permutation draws from
#'   The minimum valid value is `0`.
#' @param mean_pmf a numeric matrix. Observed consensus pmf; used only for `p_values_observed_consensus`
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @param mean_pmf_included_n_reps a integer vector. Count of non-NaN replicates (included ones) per reference point for the consensus pmf
#'   The minimum valid value is `0`.
#' @param included_n_reps a integer matrix. Count of non-NaN replicates per reference point, per study: the size of each study's
#'   draw. Every study needs at least one at every reference point
#'   The minimum valid value is `1`.
#' @param global_jsd_observed a numeric vector. Observed weighted global JSD of each study against the consensus pmf
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @return a named list with elements:
#'   \item{p_values}{a numeric vector. Empirical p-value per study, `(count + 1)/(n_permutations + 1)` over the permutations
#'     whose global JSD reached the observed one; NaN when `n_permutations` is 0}
#'   \item{p_values_observed_consensus}{a numeric vector. Like `p_values`, from the same permutations, but against the fixed observed consensus
#'     `mean_pmf`; kept for comparison only. NaN when `n_permutations` is 0}
#' @export
gjct_permutation_test <- function(n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps, included_n_reps, global_jsd_observed, random_seed = 42L) {
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    mean_pmf_counts <- .tox_as_integer_matrix(mean_pmf_counts, "mean_pmf_counts")
    mean_pmf <- .tox_as_double_matrix(mean_pmf, "mean_pmf")
    mean_pmf_included_n_reps <- .tox_as_integer_vector(mean_pmf_included_n_reps, "mean_pmf_included_n_reps")
    included_n_reps <- .tox_as_integer_matrix(included_n_reps, "included_n_reps")
    global_jsd_observed <- .tox_as_double_vector(global_jsd_observed, "global_jsd_observed")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    if (dim(mean_pmf)[2] != dim(mean_pmf_counts)[2])
        .tox_shape_error("mean_pmf", dim(mean_pmf)[2], "mean_pmf_counts", dim(mean_pmf_counts)[2])
    if (dim(mean_pmf)[1] != dim(mean_pmf_counts)[1])
        .tox_shape_error("mean_pmf", dim(mean_pmf)[1], "mean_pmf_counts", dim(mean_pmf_counts)[1])
    if (length(mean_pmf_included_n_reps) != dim(mean_pmf_counts)[1])
        .tox_shape_error("mean_pmf_included_n_reps", length(mean_pmf_included_n_reps), "mean_pmf_counts", dim(mean_pmf_counts)[1])
    if (dim(included_n_reps)[1] != dim(mean_pmf_counts)[1])
        .tox_shape_error("included_n_reps", dim(included_n_reps)[1], "mean_pmf_counts", dim(mean_pmf_counts)[1])
    if (length(global_jsd_observed) != dim(included_n_reps)[2])
        .tox_shape_error("global_jsd_observed", length(global_jsd_observed), "included_n_reps", dim(included_n_reps)[2])

    .result <- .Call("gjct_permutation_test_call", n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps, included_n_reps, global_jsd_observed, random_seed)
    .arguments <- c("n_permutations", "n_bins", "n_points", "n_studies", "mean_pmf_counts", "mean_pmf", "mean_pmf_included_n_reps", "included_n_reps", "global_jsd_observed", "p_values", "p_values_observed_consensus", "random_seed", "ierr")
    .sources <- c(NA_character_, "mean_pmf_counts", "mean_pmf_counts", "included_n_reps", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        p_values = .result$p_values,
        p_values_observed_consensus = .result$p_values_observed_consensus
    )
}

#' Bootstrap a confidence interval for each study's global JSD by resampling the pooled consensus histogram
#'
#' In every bootstrap, each study's histogram at each reference point is drawn with replacement
#' (multinomial) from the pooled consensus counts `mean_pmf_counts`, with as many residuals as
#' the study observed; the consensus pmf, each study's JSD against it and its weighted global
#' JSD are recomputed from the draw.
#'
#' `confidence_interval` is input and output: its incoming `[lower, upper]` values seed every
#' slot of the bottom-k/top-k heaps, so they are only displaced by more extreme bootstrap
#' values. On return it holds `[largest of the n_bootstrapping_top_k_jsds smallest values,
#' smallest of the n_bootstrapping_top_k_jsds largest values]`.
#'
#' Errors (`ierr`), with `confidence_interval` left unchanged: `ERR_INVALID_INPUT`, before any
#' draw, when at some reference point the studies' counts do not add up to
#' `mean_pmf_included_n_reps` or that total is not the sum of `mean_pmf_counts` (as
#' \code{\link{create_mean_pmf}} produces
#' them); `ERR_ALLOC_FAIL` when the random number generator cannot be created. A failed draw
#' later on is recorded (first failure only) without stopping the bootstrap.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_stats::bootstrap_histogram}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bootstraps a integer scalar. Number of bootstrap resamples to perform
#'   The minimum valid value is `1`.
#' @param mean_pmf_counts a integer matrix. Pooled residual counts per bin of the consensus pmf, the pool every bootstrap draws from
#'   The minimum valid value is `0`.
#' @param mean_pmf_included_n_reps a integer vector. Count of non-NaN replicates per reference point of the pool
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
#' @return a numeric matrix. On entry the seed values of the heaps (see above); on return the bootstrapped
#'   `[lower, upper]` interval per study
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @export
bootstrap_histogram <- function(n_bootstraps, mean_pmf_counts, mean_pmf_included_n_reps, included_n_reps, confidence_interval, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L) {
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    mean_pmf_counts <- .tox_as_integer_matrix(mean_pmf_counts, "mean_pmf_counts")
    mean_pmf_included_n_reps <- .tox_as_integer_vector(mean_pmf_included_n_reps, "mean_pmf_included_n_reps")
    included_n_reps <- .tox_as_integer_matrix(included_n_reps, "included_n_reps")
    confidence_interval <- .tox_as_double_matrix(confidence_interval, "confidence_interval")
    two_sided_bootstrapping_significance_level <- .tox_as_double_scalar(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    if (length(mean_pmf_included_n_reps) != dim(mean_pmf_counts)[1])
        .tox_shape_error("mean_pmf_included_n_reps", length(mean_pmf_included_n_reps), "mean_pmf_counts", dim(mean_pmf_counts)[1])
    if (dim(included_n_reps)[1] != dim(mean_pmf_counts)[1])
        .tox_shape_error("included_n_reps", dim(included_n_reps)[1], "mean_pmf_counts", dim(mean_pmf_counts)[1])
    if (dim(confidence_interval)[2] != dim(included_n_reps)[2])
        .tox_shape_error("confidence_interval", dim(confidence_interval)[2], "included_n_reps", dim(included_n_reps)[2])

    .result <- .Call("bootstrap_histogram_call", n_bootstraps, mean_pmf_counts, mean_pmf_included_n_reps, included_n_reps, confidence_interval, two_sided_bootstrapping_significance_level, random_seed)
    .arguments <- c("n_bootstraps", "n_bins", "n_points", "n_studies", "mean_pmf_counts", "mean_pmf_included_n_reps", "included_n_reps", "confidence_interval", "two_sided_bootstrapping_significance_level", "random_seed", "ierr")
    .sources <- c(NA_character_, "mean_pmf_counts", "mean_pmf_counts", "included_n_reps", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$confidence_interval
}
