# Generated. Do not edit.

#' Recommend the bootstrap top-/bottom-k heap size for a two-sided confidence interval
#'
#' Sizes the top-k/bottom-k heaps
#' \code{\link{bootstrap_histogram}} uses
#' for the interval's endpoints: `n_top_k = max(1, floor(level/100 * n_bootstraps))`, e.g. 25 for
#' the default 2.5 and `n_bootstraps=1000`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_stats_impl::calc_js_comp_test_n_top_k_jsds}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bootstraps a integer scalar. Number of bootstrap resamples that will be performed
#'   The minimum valid value is `1`.
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Two-sided significance level, as a percentage reserved at each end of the bootstrap
#'   distribution (e.g. 2.5 for a 95% two-sided interval)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @return a integer scalar. Number of elements to keep at each end of the bootstrap distribution, at least 1
#' @export
calc_js_comp_test_n_top_k_jsds <- function(n_bootstraps, two_sided_bootstrapping_significance_level = 2.5) {
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    two_sided_bootstrapping_significance_level <- .tox_as_double_scalar(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    .result <- .Call("calc_js_comp_test_n_top_k_jsds_call", n_bootstraps, two_sided_bootstrapping_significance_level)
    .arguments <- c("n_bootstraps", "two_sided_bootstrapping_significance_level", "n_top_k")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$n_top_k
}
