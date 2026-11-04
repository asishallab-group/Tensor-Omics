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
#' This test permutes the pooled residual counts within fixed histograms (residual level); the
#' study-label Permutation-Test, which shuffles the samples' study labels and reruns the whole
#' JS-Comp-Test, is
#' \code{\link{run_js_comp_test_label_permutation_test}}.
#'
#' Errors (`ierr`): `ERR_INVALID_INPUT`, before any draw, when at some reference point the
#' studies' counts do not add up to the pool's total (`mean_pmf_included_n_reps`) or that total
#' is not the sum of the pool's bin counts; `ERR_ALLOC_FAIL` when the random number generator
#' cannot be created. The p-values are then NaN. The count check, together with the wrappers'
#' bounds on the counts, establishes every draw's preconditions, so no draw can fail later on.
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
    if (dim(mean_pmf)[1] != dim(mean_pmf_counts)[1])
        .tox_shape_error("mean_pmf", 1L, dim(mean_pmf)[1], "'mean_pmf_counts' implies n_points == ", dim(mean_pmf_counts)[1], "n_points, n_bins")
    if (dim(mean_pmf)[2] != dim(mean_pmf_counts)[2])
        .tox_shape_error("mean_pmf", 2L, dim(mean_pmf)[2], "'mean_pmf_counts' implies n_bins == ", dim(mean_pmf_counts)[2], "n_points, n_bins")
    if (length(mean_pmf_included_n_reps) != dim(mean_pmf_counts)[1])
        .tox_shape_error("mean_pmf_included_n_reps", 1L, length(mean_pmf_included_n_reps), "'mean_pmf_counts' implies n_points == ", dim(mean_pmf_counts)[1], "n_points")
    if (dim(included_n_reps)[1] != dim(mean_pmf_counts)[1])
        .tox_shape_error("included_n_reps", 1L, dim(included_n_reps)[1], "'mean_pmf_counts' implies n_points == ", dim(mean_pmf_counts)[1], "n_points, n_studies")
    if (length(global_jsd_observed) != dim(included_n_reps)[2])
        .tox_shape_error("global_jsd_observed", 1L, length(global_jsd_observed), "'included_n_reps' implies n_studies == ", dim(included_n_reps)[2], "n_studies")

    .result <- .Call("gjct_permutation_test_call", n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps, included_n_reps, global_jsd_observed, random_seed)
    .arguments <- c("n_permutations", "n_bins", "n_points", "n_studies", "mean_pmf_counts", "mean_pmf", "mean_pmf_included_n_reps", "included_n_reps", "global_jsd_observed", "p_values", "p_values_observed_consensus", "random_seed", "ierr")
    .sources <- c(NA_character_, "mean_pmf_counts", "mean_pmf_counts", "included_n_reps", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        p_values = .result$p_values,
        p_values_observed_consensus = .result$p_values_observed_consensus
    )
}

#' Bootstrap an interval for each study's global JSD by resampling the studies' histograms
#'
#' In every bootstrap, each study's histogram at each reference point is drawn again with
#' replacement (multinomial), with as many residuals as the study observed
#' (`included_n_reps`); the consensus pmf, each study's JSD against it and its weighted global
#' JSD are recomputed from the draw. What a study is drawn from is `bootstrap_mode`:
#'
#' - by default its own counts (`counts`): the interval shows how much the observed JSD would
#' move if each study's own residuals were drawn again -- a robustness interval of the
#' observed value. Drawing again adds sampling noise, so the bootstrapped JSDs are larger
#' than the observed one on average: when the studies hardly differ, the interval can lie
#' entirely above the observed JSD;
#' - or the pooled consensus counts of all studies (`mean_pmf_counts`): the interval shows the
#' spread of the JSD expected if the studies did not differ at all -- a noise band, not the
#' robustness of the observed value. This was the only bootstrap before the mode existed and
#' is kept for comparison.
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
#' them), or, in the per-study mode, when a study's `counts` at some reference point do not add
#' up to its `included_n_reps`; `ERR_ALLOC_FAIL` when the random number generator cannot be
#' created. The count checks, together with the wrappers' bounds on the counts, establish every
#' draw's preconditions, so no draw can fail later on.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_stats::bootstrap_histogram}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bootstraps a integer scalar. Number of bootstrap resamples to perform
#'   The minimum valid value is `1`.
#' @param mean_pmf_counts a integer matrix. Pooled residual counts per bin of the consensus pmf: what every study is drawn from in
#'   the pooled mode
#'   The minimum valid value is `0`.
#' @param counts a integer array of rank 3. Every study's own residual counts per bin: what each study is drawn from in the
#'   per-study mode (the default); at each reference point they must add up to the
#'   study's `included_n_reps` there. Not drawn from in the pooled mode, but still required
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
#' @param bootstrap_mode a string, one of "bootstrap_per_study", "bootstrap_pooled". What each study's histograms are drawn from (see above). By default each study's own
#'   counts: a robustness interval of the observed JSD. The pooled counts give a noise
#'   band of "no difference between the studies" instead, kept for comparison
#'
#'   The default value is `"bootstrap_per_study"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @return a numeric matrix. On entry the seed values of the heaps (see above); on return the bootstrapped
#'   `[lower, upper]` interval per study
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @export
bootstrap_histogram <- function(n_bootstraps, mean_pmf_counts, counts, mean_pmf_included_n_reps, included_n_reps, confidence_interval, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, bootstrap_mode = "bootstrap_per_study") {
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    mean_pmf_counts <- .tox_as_integer_matrix(mean_pmf_counts, "mean_pmf_counts")
    counts <- .tox_as_integer_array(counts, "counts", 3L)
    mean_pmf_included_n_reps <- .tox_as_integer_vector(mean_pmf_included_n_reps, "mean_pmf_included_n_reps")
    included_n_reps <- .tox_as_integer_matrix(included_n_reps, "included_n_reps")
    confidence_interval <- .tox_as_double_matrix(confidence_interval, "confidence_interval")
    two_sided_bootstrapping_significance_level <- .tox_as_double_scalar(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    bootstrap_mode <- .tox_as_mode(bootstrap_mode, "bootstrap_mode", c("bootstrap_per_study", "bootstrap_pooled"))
    if (dim(counts)[1] != dim(mean_pmf_counts)[1])
        .tox_shape_error("counts", 1L, dim(counts)[1], "'mean_pmf_counts' implies n_points == ", dim(mean_pmf_counts)[1], "n_points, n_bins, n_studies")
    if (dim(counts)[2] != dim(mean_pmf_counts)[2])
        .tox_shape_error("counts", 2L, dim(counts)[2], "'mean_pmf_counts' implies n_bins == ", dim(mean_pmf_counts)[2], "n_points, n_bins, n_studies")
    if (length(mean_pmf_included_n_reps) != dim(mean_pmf_counts)[1])
        .tox_shape_error("mean_pmf_included_n_reps", 1L, length(mean_pmf_included_n_reps), "'mean_pmf_counts' implies n_points == ", dim(mean_pmf_counts)[1], "n_points")
    if (dim(included_n_reps)[1] != dim(mean_pmf_counts)[1])
        .tox_shape_error("included_n_reps", 1L, dim(included_n_reps)[1], "'mean_pmf_counts' implies n_points == ", dim(mean_pmf_counts)[1], "n_points, n_studies")
    if (dim(included_n_reps)[2] != dim(counts)[3])
        .tox_shape_error("included_n_reps", 2L, dim(included_n_reps)[2], "'counts' implies n_studies == ", dim(counts)[3], "n_points, n_studies")
    if (dim(confidence_interval)[1] != 2L)
        .tox_shape_error("confidence_interval", 1L, dim(confidence_interval)[1], "its declared shape fixes it at ", 2L, "2, n_studies")
    if (dim(confidence_interval)[2] != dim(counts)[3])
        .tox_shape_error("confidence_interval", 2L, dim(confidence_interval)[2], "'counts' implies n_studies == ", dim(counts)[3], "2, n_studies")

    .result <- .Call("bootstrap_histogram_call", n_bootstraps, mean_pmf_counts, counts, mean_pmf_included_n_reps, included_n_reps, confidence_interval, two_sided_bootstrapping_significance_level, random_seed, bootstrap_mode)
    .arguments <- c("n_bootstraps", "n_bins", "n_points", "n_studies", "mean_pmf_counts", "counts", "mean_pmf_included_n_reps", "included_n_reps", "confidence_interval", "two_sided_bootstrapping_significance_level", "random_seed", "bootstrap_mode", "ierr")
    .sources <- c(NA_character_, "mean_pmf_counts", "mean_pmf_counts", "counts", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$confidence_interval
}
