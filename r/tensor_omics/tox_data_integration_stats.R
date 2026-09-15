# Generated. Do not edit.

#' Estimate how likely each study's observed weighted global JSD is to occur by chance
#'
#' Ported from 125-stabilize-jscomp's `gjct_permutation_test_helper`, generalized from its
#' hardcoded 2-study shuffle to a K-study GSL `random_multiv_hypergeom` resample from the
#' pooled `mean_pmf_counts` (without replacement, per reference point, per study). Impure:
#' draws random numbers, so it carries `ierr` for the one genuine runtime failure a validated
#' caller cannot foresee -- `create_rng` failing to allocate the GSL generator -- folded in as
#' `ERR_ALLOC_FAIL`; a later `random_multiv_hypergeom` failure (which validated,
#' internally-consistent inputs should never trigger) is likewise folded in, first failure
#' only, without stopping the resampling already in flight, matching
#' \code{\link{bootstrap_histogram}}'s own
#' precedent.
#'
#' Known limitation: ported verbatim from 125-stabilize-jscomp's p-value formula, WITHOUT the
#' `(1+count)/(n+1)` Laplace add-one correction -- `p_values(i) = anint(count(i))/n_permutations`,
#' so a study whose observed JSD is never reached by any permutation gets `p = 0.0` exactly,
#' not the `1/(n_permutations+1)` a Laplace-corrected formula would give. Deliberately ported
#' as-is; see the project's JSD-Comp-Test follow-up issue for the fix.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_stats::gjct_permutation_test}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_permutations a integer scalar. Number of permutations to perform
#'   The minimum valid value is `0`.
#' @param mean_pmf_counts a integer matrix. Absolute counts of a residual per bin for the consensus pmf -- the pool each
#'   permutation resamples from without replacement, per reference point
#'   The minimum valid value is `0`.
#' @param mean_pmf a numeric matrix. The consensus pmf built from all studies' pmfs
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @param mean_pmf_included_n_reps a integer vector. Count of non-NaN replicates (included ones) per reference point for the consensus pmf
#'   The minimum valid value is `0`.
#' @param included_n_reps a integer matrix. Count of non-NaN replicates (included ones) per reference point, per study -- how
#'   many elements are drawn (without replacement) from the pooled pool per study
#'   The minimum valid value is `0`.
#' @param global_jsd_observed a numeric vector. Observed weighted global JSD of each study against the consensus pmf
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @return a numeric vector. Empirical p-value per study: the fraction of permutations whose resampled global
#'   JSD reached or exceeded the observed value -- see the known-limitation note above
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
        .tox_shape_error("mean_pmf", dim(mean_pmf)[1], "mean_pmf_counts", dim(mean_pmf_counts)[1])
    if (dim(mean_pmf)[2] != dim(mean_pmf_counts)[2])
        .tox_shape_error("mean_pmf", dim(mean_pmf)[2], "mean_pmf_counts", dim(mean_pmf_counts)[2])
    if (length(mean_pmf_included_n_reps) != dim(mean_pmf_counts)[2])
        .tox_shape_error("mean_pmf_included_n_reps", length(mean_pmf_included_n_reps), "mean_pmf_counts", dim(mean_pmf_counts)[2])
    if (dim(included_n_reps)[1] != dim(mean_pmf_counts)[2])
        .tox_shape_error("included_n_reps", dim(included_n_reps)[1], "mean_pmf_counts", dim(mean_pmf_counts)[2])
    if (length(global_jsd_observed) != dim(included_n_reps)[2])
        .tox_shape_error("global_jsd_observed", length(global_jsd_observed), "included_n_reps", dim(included_n_reps)[2])

    .result <- .Call("gjct_permutation_test_call", n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps, included_n_reps, global_jsd_observed, random_seed)
    .arguments <- c("n_permutations", "n_bins", "n_points", "n_studies", "mean_pmf_counts", "mean_pmf", "mean_pmf_included_n_reps", "included_n_reps", "global_jsd_observed", "p_values", "random_seed", "ierr")
    .sources <- c(NA_character_, "mean_pmf_counts", "mean_pmf_counts", "included_n_reps", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$p_values
}
