# Generated. Do not edit.

#' Estimate the histogram bin count for one (n_points, n_neighbors) candidate
#'
#' Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
#' (`R = shared_residual_range`), each clamped to
#' \code{MAX_N_BINS} and returned as
#' `sturges_bins`/`fd_bins`; `n_bins` is their maximum. Both take the sample size to be
#' `max_n_reps_all_studies*n_neighbors` (one study's padded slot count), so they can differ
#' from the textbook rules. The JS-Comp-Test uses them as diagnostics only, never to decide a
#' bin count.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_binning::estimate_bin_count}, whose argument names
#' are the ones an error message reports.
#'
#' This entry point seeds \code{residuals_perm} and sorts it by \code{residuals}.
#' Call \code{estimate_bin_count_expert} to do that yourself.
#'
#' @param residuals a numeric vector. Pooled signed residuals across all studies, reference points and neighbors
#'   NaN is permitted for this value.
#' @param max_n_reps_all_studies a integer scalar. Maximum number of replicates across all studies
#'   The minimum valid value is `1`.
#' @param n_neighbors a integer scalar. Neighborhood size of the candidate under test
#'   The minimum valid value is `1`.
#' @param shared_residual_range a numeric scalar. Computed residual range (R)
#'   The minimum valid value is `0.0`.
#' @return a named list with elements:
#'   \item{n_bins}{a integer scalar. Estimated number of histogram bins, at least 1 and at most MAX_N_BINS: max(sturges_bins, fd_bins)}
#'   \item{sturges_bins}{a integer scalar. Sturges' rule estimate alone, at least 1 and at most MAX_N_BINS}
#'   \item{fd_bins}{a integer scalar. Freedman-Diaconis rule estimate alone, at least 1 and at most MAX_N_BINS; falls back
#'     to sturges_bins when the interquartile range is (close to) zero or the ratio is NaN}
#' @export
estimate_bin_count <- function(residuals, max_n_reps_all_studies, n_neighbors, shared_residual_range) {
    residuals <- .tox_as_double_vector(residuals, "residuals")
    max_n_reps_all_studies <- .tox_as_integer_scalar(max_n_reps_all_studies, "max_n_reps_all_studies")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    shared_residual_range <- .tox_as_double_scalar(shared_residual_range, "shared_residual_range")
    .result <- .Call("estimate_bin_count_call", residuals, max_n_reps_all_studies, n_neighbors, shared_residual_range)
    .arguments <- c("residuals", "n_residuals", "max_n_reps_all_studies", "n_neighbors", "shared_residual_range", "n_bins", "sturges_bins", "fd_bins", "ierr")
    .sources <- c(NA_character_, "residuals", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        n_bins = .result$n_bins,
        sturges_bins = .result$sturges_bins,
        fd_bins = .result$fd_bins
    )
}

#' Estimate the histogram bin count for one (n_points, n_neighbors) candidate
#'
#' Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
#' (`R = shared_residual_range`), each clamped to
#' \code{MAX_N_BINS} and returned as
#' `sturges_bins`/`fd_bins`; `n_bins` is their maximum. Both take the sample size to be
#' `max_n_reps_all_studies*n_neighbors` (one study's padded slot count), so they can differ
#' from the textbook rules. The JS-Comp-Test uses them as diagnostics only, never to decide a
#' bin count.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_binning::estimate_bin_count_expert}, whose argument names
#' are the ones an error message reports.
#'
#' The expert entry point: you supply \code{residuals_perm} yourself.
#' \code{estimate_bin_count} seeds \code{residuals_perm} and sorts it by \code{residuals}.
#'
#' @param residuals a numeric vector. Pooled signed residuals across all studies, reference points and neighbors
#'   NaN is permitted for this value.
#' @param residuals_perm a integer vector. Sorting permutation for `residuals`, ascending, NaN last
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_residuals`.
#' @param max_n_reps_all_studies a integer scalar. Maximum number of replicates across all studies
#'   The minimum valid value is `1`.
#' @param n_neighbors a integer scalar. Neighborhood size of the candidate under test
#'   The minimum valid value is `1`.
#' @param shared_residual_range a numeric scalar. Computed residual range (R)
#'   The minimum valid value is `0.0`.
#' @return a named list with elements:
#'   \item{n_bins}{a integer scalar. Estimated number of histogram bins, at least 1 and at most MAX_N_BINS: max(sturges_bins, fd_bins)}
#'   \item{sturges_bins}{a integer scalar. Sturges' rule estimate alone, at least 1 and at most MAX_N_BINS}
#'   \item{fd_bins}{a integer scalar. Freedman-Diaconis rule estimate alone, at least 1 and at most MAX_N_BINS; falls back
#'     to sturges_bins when the interquartile range is (close to) zero or the ratio is NaN}
#' @export
estimate_bin_count_expert <- function(residuals, residuals_perm, max_n_reps_all_studies, n_neighbors, shared_residual_range) {
    residuals <- .tox_as_double_vector(residuals, "residuals")
    residuals_perm <- .tox_as_integer_vector(residuals_perm, "residuals_perm")
    max_n_reps_all_studies <- .tox_as_integer_scalar(max_n_reps_all_studies, "max_n_reps_all_studies")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    shared_residual_range <- .tox_as_double_scalar(shared_residual_range, "shared_residual_range")
    if (length(residuals_perm) != length(residuals))
        .tox_shape_error("residuals_perm", length(residuals_perm), "residuals", length(residuals))

    .result <- .Call("estimate_bin_count_expert_call", residuals, residuals_perm, max_n_reps_all_studies, n_neighbors, shared_residual_range)
    .arguments <- c("residuals", "residuals_perm", "n_residuals", "max_n_reps_all_studies", "n_neighbors", "shared_residual_range", "n_bins", "sturges_bins", "fd_bins", "ierr")
    .sources <- c(NA_character_, NA_character_, "residuals", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        n_bins = .result$n_bins,
        sturges_bins = .result$sturges_bins,
        fd_bins = .result$fd_bins
    )
}

#' Determine one neighborhood's occupancy-constrained histogram bin count (Issue #187)
#'
#' For one neighborhood's pooled residuals (all its neighbors, all studies), finds the largest
#' bin count `M` in `[m_min, m_max]` whose equal-width histogram over the neighborhood's own range
#' `[shared_residual_range_low, shared_residual_range_high]` (the
#' `lower_residual_range_quantile`/`upper_residual_range_quantile` quantiles of these residuals)
#' has every bin at or above `min_residuals_per_bin`. A geometric search (factor
#' `gamma_occupancy`) brackets the first failing `M`, then every integer below it is tested,
#' because occupancy is not monotonic in `M`.
#'
#' `occupancy_failed = TRUE` (even `m_min` bins fail, or every residual is NaN) means the caller
#' should reject the neighborhood; `selected_n_bins` is then `m_min` only so that it is defined.
#' `sturges_bins`/`fd_bins` come from
#' \code{\link{estimate_bin_count}} and are
#' diagnostics only.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_binning::determine_bin_count_occupancy}, whose argument names
#' are the ones an error message reports.
#'
#' This entry point seeds \code{pooled_residuals_perm} and sorts it by \code{pooled_residuals} and computes \code{tmp_bin_counts} for you.
#' Call \code{determine_bin_count_occupancy_expert} to do that yourself.
#'
#' @param pooled_residuals a numeric vector. Pooled signed residuals for one neighborhood, across all its neighbors and all studies
#'   NaN is permitted for this value.
#' @param max_n_reps_all_studies a integer scalar. Maximum number of replicates across all studies
#'   The minimum valid value is `1`.
#' @param n_neighbors a integer scalar. Neighborhood size of the candidate under test
#'   The minimum valid value is `1`.
#' @param m_min a integer scalar. Smallest candidate bin count the search will ever test (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest candidate bin count the search will ever test (M_max); raised to `m_min` if
#'   smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach for a candidate bin count to
#'   be admissible (n_min)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor for the coarse search stage; must exceed 1 or the search
#'   never advances
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] for this neighborhood's own lower residual-range bound
#'   (shared_residual_range_low)
#'   Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
#'   with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] for this neighborhood's own upper residual-range bound
#'   (shared_residual_range_high)
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @return a named list with elements:
#'   \item{selected_n_bins}{a integer scalar. The chosen M_j: the largest bin count in [m_min, m_max] whose pooled histogram has
#'     every bin at or above min_residuals_per_bin; m_min when occupancy_failed}
#'   \item{occupancy_failed}{a logical scalar. `TRUE` iff even m_min bins fail the occupancy criterion (or every pooled residual is
#'     NaN); the caller should then reject this neighborhood}
#'   \item{shared_residual_range_low}{a numeric scalar. This neighborhood's lower residual-range bound (R_low): the
#'     lower_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN}
#'   \item{shared_residual_range_high}{a numeric scalar. This neighborhood's upper residual-range bound (R_high): the
#'     upper_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN}
#'   \item{n_pooled_residuals}{a integer scalar. Count of non-NaN pooled residuals (N_j)}
#'   \item{min_bin_occupancy}{a integer scalar. Minimum bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{mean_bin_occupancy}{a numeric scalar. Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
#'     occupancy_failed}
#'   \item{max_bin_occupancy}{a integer scalar. Maximum bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{sturges_bins}{a integer scalar. Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)}
#'   \item{fd_bins}{a integer scalar. Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)}
#' @export
determine_bin_count_occupancy <- function(pooled_residuals, max_n_reps_all_studies, n_neighbors, m_min = 3L, m_max = 120L, min_residuals_per_bin = 10L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95) {
    pooled_residuals <- .tox_as_double_vector(pooled_residuals, "pooled_residuals")
    max_n_reps_all_studies <- .tox_as_integer_scalar(max_n_reps_all_studies, "max_n_reps_all_studies")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    gamma_occupancy <- .tox_as_double_scalar(gamma_occupancy, "gamma_occupancy")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    .result <- .Call("determine_bin_count_occupancy_call", pooled_residuals, max_n_reps_all_studies, n_neighbors, m_min, m_max, min_residuals_per_bin, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile)
    .arguments <- c("pooled_residuals", "n_residuals", "max_n_reps_all_studies", "n_neighbors", "selected_n_bins", "occupancy_failed", "shared_residual_range_low", "shared_residual_range_high", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "m_min", "m_max", "min_residuals_per_bin", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr")
    .sources <- c(NA_character_, "pooled_residuals", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        selected_n_bins = .result$selected_n_bins,
        occupancy_failed = .result$occupancy_failed,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        n_pooled_residuals = .result$n_pooled_residuals,
        min_bin_occupancy = .result$min_bin_occupancy,
        mean_bin_occupancy = .result$mean_bin_occupancy,
        max_bin_occupancy = .result$max_bin_occupancy,
        sturges_bins = .result$sturges_bins,
        fd_bins = .result$fd_bins
    )
}

#' Determine one neighborhood's occupancy-constrained histogram bin count (Issue #187)
#'
#' For one neighborhood's pooled residuals (all its neighbors, all studies), finds the largest
#' bin count `M` in `[m_min, m_max]` whose equal-width histogram over the neighborhood's own range
#' `[shared_residual_range_low, shared_residual_range_high]` (the
#' `lower_residual_range_quantile`/`upper_residual_range_quantile` quantiles of these residuals)
#' has every bin at or above `min_residuals_per_bin`. A geometric search (factor
#' `gamma_occupancy`) brackets the first failing `M`, then every integer below it is tested,
#' because occupancy is not monotonic in `M`.
#'
#' `occupancy_failed = TRUE` (even `m_min` bins fail, or every residual is NaN) means the caller
#' should reject the neighborhood; `selected_n_bins` is then `m_min` only so that it is defined.
#' `sturges_bins`/`fd_bins` come from
#' \code{\link{estimate_bin_count}} and are
#' diagnostics only.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_binning::determine_bin_count_occupancy_expert}, whose argument names
#' are the ones an error message reports.
#'
#' The expert entry point: you supply \code{pooled_residuals_perm} and \code{tmp_bin_counts} yourself.
#' \code{determine_bin_count_occupancy} seeds \code{pooled_residuals_perm} and sorts it by \code{pooled_residuals} and computes \code{tmp_bin_counts} for you.
#'
#' @param pooled_residuals a numeric vector. Pooled signed residuals for one neighborhood, across all its neighbors and all studies
#'   NaN is permitted for this value.
#' @param pooled_residuals_perm a integer vector. Sorting permutation for `pooled_residuals`, ascending, NaN last
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_residuals`.
#' @param max_n_reps_all_studies a integer scalar. Maximum number of replicates across all studies
#'   The minimum valid value is `1`.
#' @param n_neighbors a integer scalar. Neighborhood size of the candidate under test
#'   The minimum valid value is `1`.
#' @param m_min a integer scalar. Smallest candidate bin count the search will ever test (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest candidate bin count the search will ever test (M_max); raised to `m_min` if
#'   smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach for a candidate bin count to
#'   be admissible (n_min)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor for the coarse search stage; must exceed 1 or the search
#'   never advances
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] for this neighborhood's own lower residual-range bound
#'   (shared_residual_range_low)
#'   Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
#'   with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] for this neighborhood's own upper residual-range bound
#'   (shared_residual_range_high)
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @return a named list with elements:
#'   \item{selected_n_bins}{a integer scalar. The chosen M_j: the largest bin count in [m_min, m_max] whose pooled histogram has
#'     every bin at or above min_residuals_per_bin; m_min when occupancy_failed}
#'   \item{occupancy_failed}{a logical scalar. `TRUE` iff even m_min bins fail the occupancy criterion (or every pooled residual is
#'     NaN); the caller should then reject this neighborhood}
#'   \item{shared_residual_range_low}{a numeric scalar. This neighborhood's lower residual-range bound (R_low): the
#'     lower_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN}
#'   \item{shared_residual_range_high}{a numeric scalar. This neighborhood's upper residual-range bound (R_high): the
#'     upper_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN}
#'   \item{n_pooled_residuals}{a integer scalar. Count of non-NaN pooled residuals (N_j)}
#'   \item{min_bin_occupancy}{a integer scalar. Minimum bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{mean_bin_occupancy}{a numeric scalar. Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
#'     occupancy_failed}
#'   \item{max_bin_occupancy}{a integer scalar. Maximum bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{sturges_bins}{a integer scalar. Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)}
#'   \item{fd_bins}{a integer scalar. Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)}
#' @export
determine_bin_count_occupancy_expert <- function(pooled_residuals, pooled_residuals_perm, max_n_reps_all_studies, n_neighbors, m_min = 3L, m_max = 120L, min_residuals_per_bin = 10L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95) {
    pooled_residuals <- .tox_as_double_vector(pooled_residuals, "pooled_residuals")
    pooled_residuals_perm <- .tox_as_integer_vector(pooled_residuals_perm, "pooled_residuals_perm")
    max_n_reps_all_studies <- .tox_as_integer_scalar(max_n_reps_all_studies, "max_n_reps_all_studies")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    gamma_occupancy <- .tox_as_double_scalar(gamma_occupancy, "gamma_occupancy")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    if (length(pooled_residuals_perm) != length(pooled_residuals))
        .tox_shape_error("pooled_residuals_perm", length(pooled_residuals_perm), "pooled_residuals", length(pooled_residuals))

    .result <- .Call("determine_bin_count_occupancy_expert_call", pooled_residuals, pooled_residuals_perm, max_n_reps_all_studies, n_neighbors, m_min, m_max, min_residuals_per_bin, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile)
    .arguments <- c("pooled_residuals", "pooled_residuals_perm", "n_residuals", "max_n_reps_all_studies", "n_neighbors", "selected_n_bins", "occupancy_failed", "shared_residual_range_low", "shared_residual_range_high", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "tmp_bin_counts", "m_min", "m_max", "min_residuals_per_bin", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr")
    .sources <- c(NA_character_, NA_character_, "pooled_residuals", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        selected_n_bins = .result$selected_n_bins,
        occupancy_failed = .result$occupancy_failed,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        n_pooled_residuals = .result$n_pooled_residuals,
        min_bin_occupancy = .result$min_bin_occupancy,
        mean_bin_occupancy = .result$mean_bin_occupancy,
        max_bin_occupancy = .result$max_bin_occupancy,
        sturges_bins = .result$sturges_bins,
        fd_bins = .result$fd_bins
    )
}

#' Exhaustive brute-force reference implementation of Issue #187's occupancy search
#'
#' Tests every candidate bin count `M` in `[m_min, m_max]` and keeps the largest one whose
#' pooled histogram satisfies the occupancy criterion. A reference for validating
#' \code{\link{determine_bin_count_occupancy}},
#' whose faster search could in principle miss a larger admissible `M`. Takes the range and
#' `n_pooled_residuals` from a prior call of that routine, so only the `M` selection is compared.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_binning::determine_bin_count_occupancy_exhaustive}, whose argument names
#' are the ones an error message reports.
#'
#' This entry point seeds \code{pooled_residuals_perm} and sorts it by \code{pooled_residuals}.
#' Call \code{determine_bin_count_occupancy_exhaustive_expert} to do that yourself.
#'
#' @param pooled_residuals a numeric vector. Pooled signed residuals for one neighborhood, across all its neighbors and all studies
#'   NaN is permitted for this value.
#' @param n_pooled_residuals a integer scalar. Count of non-NaN pooled residuals (N_j), from a prior
#'   \code{\link{determine_bin_count_occupancy}}
#'   call
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_residuals`.
#' @param shared_residual_range_low a numeric scalar. Lower bound of the histogram range (R_low), from a prior
#'   \code{\link{determine_bin_count_occupancy}}
#'   call
#' @param shared_residual_range_high a numeric scalar. Upper bound of the histogram range (R_high), from the same prior call as
#'   `shared_residual_range_low`
#'   The minimum valid value is `shared_residual_range_low`.
#' @param m_min a integer scalar. Smallest candidate bin count tested (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest candidate bin count tested (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach for a candidate bin count to
#'   be admissible (n_min)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @return a named list with elements:
#'   \item{selected_n_bins}{a integer scalar. The largest bin count in [m_min, m_max] whose pooled histogram has every bin at or
#'     above min_residuals_per_bin, found by exhaustive search; m_min when occupancy_failed}
#'   \item{occupancy_failed}{a logical scalar. `TRUE` iff no candidate bin count in [m_min, m_max] satisfies the occupancy criterion}
#'   \item{min_bin_occupancy}{a integer scalar. Minimum bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{mean_bin_occupancy}{a numeric scalar. Mean bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{max_bin_occupancy}{a integer scalar. Maximum bin count at selected_n_bins; 0 when occupancy_failed}
#' @export
determine_bin_count_occupancy_exhaustive <- function(pooled_residuals, n_pooled_residuals, shared_residual_range_low, shared_residual_range_high, m_min = 3L, m_max = 120L, min_residuals_per_bin = 10L) {
    pooled_residuals <- .tox_as_double_vector(pooled_residuals, "pooled_residuals")
    n_pooled_residuals <- .tox_as_integer_scalar(n_pooled_residuals, "n_pooled_residuals")
    shared_residual_range_low <- .tox_as_double_scalar(shared_residual_range_low, "shared_residual_range_low")
    shared_residual_range_high <- .tox_as_double_scalar(shared_residual_range_high, "shared_residual_range_high")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    .result <- .Call("determine_bin_count_occupancy_exhaustive_call", pooled_residuals, n_pooled_residuals, shared_residual_range_low, shared_residual_range_high, m_min, m_max, min_residuals_per_bin)
    .arguments <- c("pooled_residuals", "n_residuals", "n_pooled_residuals", "shared_residual_range_low", "shared_residual_range_high", "selected_n_bins", "occupancy_failed", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "m_min", "m_max", "min_residuals_per_bin", "ierr")
    .sources <- c(NA_character_, "pooled_residuals", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        selected_n_bins = .result$selected_n_bins,
        occupancy_failed = .result$occupancy_failed,
        min_bin_occupancy = .result$min_bin_occupancy,
        mean_bin_occupancy = .result$mean_bin_occupancy,
        max_bin_occupancy = .result$max_bin_occupancy
    )
}

#' Exhaustive brute-force reference implementation of Issue #187's occupancy search
#'
#' Tests every candidate bin count `M` in `[m_min, m_max]` and keeps the largest one whose
#' pooled histogram satisfies the occupancy criterion. A reference for validating
#' \code{\link{determine_bin_count_occupancy}},
#' whose faster search could in principle miss a larger admissible `M`. Takes the range and
#' `n_pooled_residuals` from a prior call of that routine, so only the `M` selection is compared.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_binning::determine_bin_count_occupancy_exhaustive_expert}, whose argument names
#' are the ones an error message reports.
#'
#' The expert entry point: you supply \code{pooled_residuals_perm} yourself.
#' \code{determine_bin_count_occupancy_exhaustive} seeds \code{pooled_residuals_perm} and sorts it by \code{pooled_residuals}.
#'
#' @param pooled_residuals a numeric vector. Pooled signed residuals for one neighborhood, across all its neighbors and all studies
#'   NaN is permitted for this value.
#' @param pooled_residuals_perm a integer vector. Sorting permutation for `pooled_residuals`, ascending, NaN last
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_residuals`.
#' @param n_pooled_residuals a integer scalar. Count of non-NaN pooled residuals (N_j), from a prior
#'   \code{\link{determine_bin_count_occupancy}}
#'   call
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_residuals`.
#' @param shared_residual_range_low a numeric scalar. Lower bound of the histogram range (R_low), from a prior
#'   \code{\link{determine_bin_count_occupancy}}
#'   call
#' @param shared_residual_range_high a numeric scalar. Upper bound of the histogram range (R_high), from the same prior call as
#'   `shared_residual_range_low`
#'   The minimum valid value is `shared_residual_range_low`.
#' @param m_min a integer scalar. Smallest candidate bin count tested (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest candidate bin count tested (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach for a candidate bin count to
#'   be admissible (n_min)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @return a named list with elements:
#'   \item{selected_n_bins}{a integer scalar. The largest bin count in [m_min, m_max] whose pooled histogram has every bin at or
#'     above min_residuals_per_bin, found by exhaustive search; m_min when occupancy_failed}
#'   \item{occupancy_failed}{a logical scalar. `TRUE` iff no candidate bin count in [m_min, m_max] satisfies the occupancy criterion}
#'   \item{min_bin_occupancy}{a integer scalar. Minimum bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{mean_bin_occupancy}{a numeric scalar. Mean bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{max_bin_occupancy}{a integer scalar. Maximum bin count at selected_n_bins; 0 when occupancy_failed}
#' @export
determine_bin_count_occupancy_exhaustive_expert <- function(pooled_residuals, pooled_residuals_perm, n_pooled_residuals, shared_residual_range_low, shared_residual_range_high, m_min = 3L, m_max = 120L, min_residuals_per_bin = 10L) {
    pooled_residuals <- .tox_as_double_vector(pooled_residuals, "pooled_residuals")
    pooled_residuals_perm <- .tox_as_integer_vector(pooled_residuals_perm, "pooled_residuals_perm")
    n_pooled_residuals <- .tox_as_integer_scalar(n_pooled_residuals, "n_pooled_residuals")
    shared_residual_range_low <- .tox_as_double_scalar(shared_residual_range_low, "shared_residual_range_low")
    shared_residual_range_high <- .tox_as_double_scalar(shared_residual_range_high, "shared_residual_range_high")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    if (length(pooled_residuals_perm) != length(pooled_residuals))
        .tox_shape_error("pooled_residuals_perm", length(pooled_residuals_perm), "pooled_residuals", length(pooled_residuals))

    .result <- .Call("determine_bin_count_occupancy_exhaustive_expert_call", pooled_residuals, pooled_residuals_perm, n_pooled_residuals, shared_residual_range_low, shared_residual_range_high, m_min, m_max, min_residuals_per_bin)
    .arguments <- c("pooled_residuals", "pooled_residuals_perm", "n_residuals", "n_pooled_residuals", "shared_residual_range_low", "shared_residual_range_high", "selected_n_bins", "occupancy_failed", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "m_min", "m_max", "min_residuals_per_bin", "ierr")
    .sources <- c(NA_character_, NA_character_, "pooled_residuals", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        selected_n_bins = .result$selected_n_bins,
        occupancy_failed = .result$occupancy_failed,
        min_bin_occupancy = .result$min_bin_occupancy,
        mean_bin_occupancy = .result$mean_bin_occupancy,
        max_bin_occupancy = .result$max_bin_occupancy
    )
}

#' Choose one neighborhood's histogram bin count by Stone's cross-validation rule
#'
#' Chooses the number of equal-width histogram bins for one neighborhood's pooled residuals by
#' Stone's leave-one-out cross-validation (Duke STA 293.3 lecture notes, ch. 6, Theorem 2): for
#' `n` values over `m` bins of width `h`, the bin count minimising
#' `J(m) = 2/(h(n-1)) - (n+1)/(h(n-1)) * sum_j p_j^2` (`p_j` the share of values in bin `j`) is
#' chosen. An alternative to
#' \code{\link{determine_bin_count_occupancy}},
#' with the same range, floor and diagnostics.
#'
#' The range is the `lower_residual_range_quantile`/`upper_residual_range_quantile` quantiles of
#' all non-NaN pooled residuals, so `h = (R_high - R_low)/m`. The score counts only the residuals
#' **inside** the range (`n` is their number). Every `m` in `[m_min, m_max]` is tested, ties
#' going to the smaller `m`. A single score (pooled mode, or one study) is bit-exact for up to
#' about 200,000 residuals, so ties are exact; a sum over studies may differ in the last bit.
#'
#' `residual_to_study` says which study each residual comes from, and `bin_count_mode` how the
#' studies are combined into one shared bin count:
#' \code{'bin_count_stone_equal'} (the
#' default) minimises the sum of every study's own score `J_k(m)`,
#' \code{'bin_count_stone_size_weighted'}
#' the sum of `n_k * J_k(m)` with `n_k` the study's inside count, and
#' \code{'bin_count_stone_pooled'} the
#' score `J(m)` of all studies' residuals together. The two per-study modes leave out a study
#' with fewer than 2 residuals inside the range.
#'
#' `min_residuals_per_bin` stays a hard floor, with the occupancy search's definition: a bin
#' count is admissible only if the pooled histogram with the outside residuals clamped into the
#' edge bins has every bin at or above it (`0` turns the floor off). The minimum is taken over
#' the admissible bin counts only.
#'
#' `stone_failed = TRUE` means no bin count could be chosen: every pooled residual is NaN, the
#' range has width 0 or less, no study has 2 residuals inside the range (pooled mode: fewer than
#' 2 inside in total), or no bin count passes the floor. `selected_n_bins` is then `m_min`, only
#' so a caller ignoring the flag has a value to build with, and the occupancy diagnostics are 0.
#' Otherwise the diagnostics describe the pooled clamped histogram at `selected_n_bins`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_binning::determine_bin_count_stone}, whose argument names
#' are the ones an error message reports.
#'
#' This entry point seeds \code{pooled_residuals_perm} and sorts it by \code{pooled_residuals} and computes \code{tmp_bin_counts}, \code{tmp_residual_to_study_perm} and \code{tmp_study_block_sizes} for you.
#' Call \code{determine_bin_count_stone_expert} to do that yourself.
#'
#' @param pooled_residuals a numeric vector. Pooled signed residuals for one neighborhood, across all its neighbors and all studies
#'   NaN is permitted for this value.
#' @param residual_to_study a integer vector. Study each pooled residual comes from
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_studies`.
#' @param n_studies a integer scalar. Number of studies
#'   The minimum valid value is `1`.
#' @param max_n_reps_all_studies a integer scalar. Maximum number of replicates across all studies
#'   The minimum valid value is `1`.
#' @param n_neighbors a integer scalar. Neighborhood size of the candidate under test
#'   The minimum valid value is `1`.
#' @param bin_count_mode a string, one of "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How the studies' scores are combined into one bin count
#'
#'   The default value is `"bin_count_stone_equal"`.
#' @param m_min a integer scalar. Smallest candidate bin count tested (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest candidate bin count tested (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin of the pooled clamped histogram must
#'   reach for a candidate bin count to be admissible (0 turns the floor off)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] for this neighborhood's own lower residual-range bound
#'   (shared_residual_range_low)
#'   Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
#'   with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] for this neighborhood's own upper residual-range bound
#'   (shared_residual_range_high)
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @return a named list with elements:
#'   \item{selected_n_bins}{a integer scalar. The chosen bin count: the admissible bin count in [m_min, m_max] with the smallest
#'     Stone score; m_min when stone_failed}
#'   \item{stone_failed}{a logical scalar. `TRUE` iff no bin count could be chosen (see above); the caller should then reject
#'     this neighborhood rather than build a histogram from selected_n_bins}
#'   \item{shared_residual_range_low}{a numeric scalar. This neighborhood's lower residual-range bound (R_low); 0.0 if every residual is NaN}
#'   \item{shared_residual_range_high}{a numeric scalar. This neighborhood's upper residual-range bound (R_high); 0.0 if every residual is NaN}
#'   \item{n_pooled_residuals}{a integer scalar. Count of non-NaN pooled residuals (N_j)}
#'   \item{min_bin_occupancy}{a integer scalar. Minimum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed}
#'   \item{mean_bin_occupancy}{a numeric scalar. Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
#'     stone_failed}
#'   \item{max_bin_occupancy}{a integer scalar. Maximum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed}
#'   \item{sturges_bins}{a integer scalar. Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)}
#'   \item{fd_bins}{a integer scalar. Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)}
#' @export
determine_bin_count_stone <- function(pooled_residuals, residual_to_study, n_studies, max_n_reps_all_studies, n_neighbors, bin_count_mode = "bin_count_stone_equal", m_min = 3L, m_max = 120L, min_residuals_per_bin = 10L, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95) {
    pooled_residuals <- .tox_as_double_vector(pooled_residuals, "pooled_residuals")
    residual_to_study <- .tox_as_integer_vector(residual_to_study, "residual_to_study")
    n_studies <- .tox_as_integer_scalar(n_studies, "n_studies")
    max_n_reps_all_studies <- .tox_as_integer_scalar(max_n_reps_all_studies, "max_n_reps_all_studies")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    bin_count_mode <- .tox_as_mode(bin_count_mode, "bin_count_mode", c("bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled"))
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    if (length(residual_to_study) != length(pooled_residuals))
        .tox_shape_error("residual_to_study", length(residual_to_study), "pooled_residuals", length(pooled_residuals))

    .result <- .Call("determine_bin_count_stone_call", pooled_residuals, residual_to_study, n_studies, max_n_reps_all_studies, n_neighbors, bin_count_mode, m_min, m_max, min_residuals_per_bin, lower_residual_range_quantile, upper_residual_range_quantile)
    .arguments <- c("pooled_residuals", "n_residuals", "residual_to_study", "n_studies", "max_n_reps_all_studies", "n_neighbors", "selected_n_bins", "stone_failed", "shared_residual_range_low", "shared_residual_range_high", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "bin_count_mode", "m_min", "m_max", "min_residuals_per_bin", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr")
    .sources <- c(NA_character_, "pooled_residuals", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        selected_n_bins = .result$selected_n_bins,
        stone_failed = .result$stone_failed,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        n_pooled_residuals = .result$n_pooled_residuals,
        min_bin_occupancy = .result$min_bin_occupancy,
        mean_bin_occupancy = .result$mean_bin_occupancy,
        max_bin_occupancy = .result$max_bin_occupancy,
        sturges_bins = .result$sturges_bins,
        fd_bins = .result$fd_bins
    )
}

#' Choose one neighborhood's histogram bin count by Stone's cross-validation rule
#'
#' Chooses the number of equal-width histogram bins for one neighborhood's pooled residuals by
#' Stone's leave-one-out cross-validation (Duke STA 293.3 lecture notes, ch. 6, Theorem 2): for
#' `n` values over `m` bins of width `h`, the bin count minimising
#' `J(m) = 2/(h(n-1)) - (n+1)/(h(n-1)) * sum_j p_j^2` (`p_j` the share of values in bin `j`) is
#' chosen. An alternative to
#' \code{\link{determine_bin_count_occupancy}},
#' with the same range, floor and diagnostics.
#'
#' The range is the `lower_residual_range_quantile`/`upper_residual_range_quantile` quantiles of
#' all non-NaN pooled residuals, so `h = (R_high - R_low)/m`. The score counts only the residuals
#' **inside** the range (`n` is their number). Every `m` in `[m_min, m_max]` is tested, ties
#' going to the smaller `m`. A single score (pooled mode, or one study) is bit-exact for up to
#' about 200,000 residuals, so ties are exact; a sum over studies may differ in the last bit.
#'
#' `residual_to_study` says which study each residual comes from, and `bin_count_mode` how the
#' studies are combined into one shared bin count:
#' \code{'bin_count_stone_equal'} (the
#' default) minimises the sum of every study's own score `J_k(m)`,
#' \code{'bin_count_stone_size_weighted'}
#' the sum of `n_k * J_k(m)` with `n_k` the study's inside count, and
#' \code{'bin_count_stone_pooled'} the
#' score `J(m)` of all studies' residuals together. The two per-study modes leave out a study
#' with fewer than 2 residuals inside the range.
#'
#' `min_residuals_per_bin` stays a hard floor, with the occupancy search's definition: a bin
#' count is admissible only if the pooled histogram with the outside residuals clamped into the
#' edge bins has every bin at or above it (`0` turns the floor off). The minimum is taken over
#' the admissible bin counts only.
#'
#' `stone_failed = TRUE` means no bin count could be chosen: every pooled residual is NaN, the
#' range has width 0 or less, no study has 2 residuals inside the range (pooled mode: fewer than
#' 2 inside in total), or no bin count passes the floor. `selected_n_bins` is then `m_min`, only
#' so a caller ignoring the flag has a value to build with, and the occupancy diagnostics are 0.
#' Otherwise the diagnostics describe the pooled clamped histogram at `selected_n_bins`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_binning::determine_bin_count_stone_expert}, whose argument names
#' are the ones an error message reports.
#'
#' The expert entry point: you supply \code{pooled_residuals_perm}, \code{tmp_bin_counts}, \code{tmp_residual_to_study_perm} and \code{tmp_study_block_sizes} yourself.
#' \code{determine_bin_count_stone} seeds \code{pooled_residuals_perm} and sorts it by \code{pooled_residuals} and computes \code{tmp_bin_counts}, \code{tmp_residual_to_study_perm} and \code{tmp_study_block_sizes} for you.
#'
#' @param pooled_residuals a numeric vector. Pooled signed residuals for one neighborhood, across all its neighbors and all studies
#'   NaN is permitted for this value.
#' @param pooled_residuals_perm a integer vector. Sorting permutation for `pooled_residuals`, ascending, NaN last
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_residuals`.
#' @param residual_to_study a integer vector. Study each pooled residual comes from
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_studies`.
#' @param n_studies a integer scalar. Number of studies
#'   The minimum valid value is `1`.
#' @param max_n_reps_all_studies a integer scalar. Maximum number of replicates across all studies
#'   The minimum valid value is `1`.
#' @param n_neighbors a integer scalar. Neighborhood size of the candidate under test
#'   The minimum valid value is `1`.
#' @param bin_count_mode a string, one of "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How the studies' scores are combined into one bin count
#'
#'   The default value is `"bin_count_stone_equal"`.
#' @param m_min a integer scalar. Smallest candidate bin count tested (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest candidate bin count tested (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin of the pooled clamped histogram must
#'   reach for a candidate bin count to be admissible (0 turns the floor off)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] for this neighborhood's own lower residual-range bound
#'   (shared_residual_range_low)
#'   Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
#'   with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] for this neighborhood's own upper residual-range bound
#'   (shared_residual_range_high)
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @return a named list with elements:
#'   \item{selected_n_bins}{a integer scalar. The chosen bin count: the admissible bin count in [m_min, m_max] with the smallest
#'     Stone score; m_min when stone_failed}
#'   \item{stone_failed}{a logical scalar. `TRUE` iff no bin count could be chosen (see above); the caller should then reject
#'     this neighborhood rather than build a histogram from selected_n_bins}
#'   \item{shared_residual_range_low}{a numeric scalar. This neighborhood's lower residual-range bound (R_low); 0.0 if every residual is NaN}
#'   \item{shared_residual_range_high}{a numeric scalar. This neighborhood's upper residual-range bound (R_high); 0.0 if every residual is NaN}
#'   \item{n_pooled_residuals}{a integer scalar. Count of non-NaN pooled residuals (N_j)}
#'   \item{min_bin_occupancy}{a integer scalar. Minimum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed}
#'   \item{mean_bin_occupancy}{a numeric scalar. Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
#'     stone_failed}
#'   \item{max_bin_occupancy}{a integer scalar. Maximum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed}
#'   \item{sturges_bins}{a integer scalar. Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)}
#'   \item{fd_bins}{a integer scalar. Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)}
#' @export
determine_bin_count_stone_expert <- function(pooled_residuals, pooled_residuals_perm, residual_to_study, n_studies, max_n_reps_all_studies, n_neighbors, bin_count_mode = "bin_count_stone_equal", m_min = 3L, m_max = 120L, min_residuals_per_bin = 10L, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95) {
    pooled_residuals <- .tox_as_double_vector(pooled_residuals, "pooled_residuals")
    pooled_residuals_perm <- .tox_as_integer_vector(pooled_residuals_perm, "pooled_residuals_perm")
    residual_to_study <- .tox_as_integer_vector(residual_to_study, "residual_to_study")
    n_studies <- .tox_as_integer_scalar(n_studies, "n_studies")
    max_n_reps_all_studies <- .tox_as_integer_scalar(max_n_reps_all_studies, "max_n_reps_all_studies")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    bin_count_mode <- .tox_as_mode(bin_count_mode, "bin_count_mode", c("bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled"))
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    if (length(pooled_residuals_perm) != length(pooled_residuals))
        .tox_shape_error("pooled_residuals_perm", length(pooled_residuals_perm), "pooled_residuals", length(pooled_residuals))
    if (length(residual_to_study) != length(pooled_residuals))
        .tox_shape_error("residual_to_study", length(residual_to_study), "pooled_residuals", length(pooled_residuals))

    .result <- .Call("determine_bin_count_stone_expert_call", pooled_residuals, pooled_residuals_perm, residual_to_study, n_studies, max_n_reps_all_studies, n_neighbors, bin_count_mode, m_min, m_max, min_residuals_per_bin, lower_residual_range_quantile, upper_residual_range_quantile)
    .arguments <- c("pooled_residuals", "pooled_residuals_perm", "n_residuals", "residual_to_study", "n_studies", "max_n_reps_all_studies", "n_neighbors", "selected_n_bins", "stone_failed", "shared_residual_range_low", "shared_residual_range_high", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "tmp_bin_counts", "tmp_residual_to_study_perm", "tmp_study_block_sizes", "bin_count_mode", "m_min", "m_max", "min_residuals_per_bin", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr")
    .sources <- c(NA_character_, NA_character_, "pooled_residuals", NA_character_, "tmp_study_block_sizes", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        selected_n_bins = .result$selected_n_bins,
        stone_failed = .result$stone_failed,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        n_pooled_residuals = .result$n_pooled_residuals,
        min_bin_occupancy = .result$min_bin_occupancy,
        mean_bin_occupancy = .result$mean_bin_occupancy,
        max_bin_occupancy = .result$max_bin_occupancy,
        sturges_bins = .result$sturges_bins,
        fd_bins = .result$fd_bins
    )
}
