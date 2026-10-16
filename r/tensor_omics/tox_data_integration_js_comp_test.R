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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::estimate_bin_count}, whose argument names
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::estimate_bin_count_expert}, whose argument names
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::determine_bin_count_occupancy}, whose argument names
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::determine_bin_count_occupancy_expert}, whose argument names
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::determine_bin_count_occupancy_exhaustive}, whose argument names
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::determine_bin_count_occupancy_exhaustive_expert}, whose argument names
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::determine_bin_count_stone}, whose argument names
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::determine_bin_count_stone_expert}, whose argument names
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

#' Generate the GAMMA-decay (n_points, n_neighbors) candidate grid
#'
#' Starting from `n_points_high = clamp(ceiling(4*sqrt(G)), MIN_POINTS, MAX_POINTS)` (`G` =
#' `max_n_genes_all_studies`), multiplies it by GAMMA until it drops below
#' `n_points_low = max(MIN_POINTS, ceiling(0.2*n_points_high))`, and pairs each distinct
#' `n_points` with up to `size(KX_FACTORS)` distinct `n_neighbors = max(1, floor(G/(f*n_points_high)))`.
#' Duplicate values from clamping or rounding collapse, so a small enough `G` yields a single
#' candidate. Bin counts are not part of a candidate; they are chosen per reference point later.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::generate_js_comp_test_candidates}, whose argument names
#' are the ones an error message reports.
#'
#' @param max_n_genes_all_studies a integer scalar. Maximum number of genes across all studies
#'   The minimum valid value is `1`.
#' @return a integer matrix. Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
#'   The first `n_candidates` elements will hold the results.
#' @export
generate_js_comp_test_candidates <- function(max_n_genes_all_studies) {
    max_n_genes_all_studies <- .tox_as_integer_scalar(max_n_genes_all_studies, "max_n_genes_all_studies")
    .result <- .Call("generate_js_comp_test_candidates_call", max_n_genes_all_studies)
    .arguments <- c("max_n_genes_all_studies", "candidates_n_points_n_neighbors", "n_candidates", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$candidates_n_points_n_neighbors[, seq_len(.result$n_candidates), drop = FALSE]
}

#' Generate the ascending adaptive (k_start, k_step, k_max) candidate sequence
#'
#' The adaptive counterpart of
#' \code{\link{generate_js_comp_test_candidates}}.
#' An adaptive candidate fixes only its growth knobs (each neighborhood grows from `k_start`
#' pooled entries in rounds of `k_step` up to `k_max`); its number of reference points emerges.
#' `k_start` shrinks from one candidate to the next, so the candidates **ascend** in the number
#' of reference points they are expected to produce.
#'
#' With `N = max_n_genes_all_studies * n_studies` (the padded pool size) and
#' `floor_k = max(10, 2*n_studies)`:
#'
#' - `k_start_1 = max(floor_k, ceiling(0.02 * N))`
#' - `k_start_t = floor(k_start_1 * 0.8**(t-1))`
#' - stop once `k_start_t < floor_k`, or once 16 candidates are accepted
#' - `k_step_t = ceiling(0.25 * k_start_t)`
#' - `k_max_t = ceiling(4 * k_start_t)`
#'
#' The constants are \code{ADAPTIVE_START_FRACTION} (0.02),
#' \code{ADAPTIVE_GAMMA} (0.8),
#' \code{ADAPTIVE_K_STEP_FRACTION} (0.25),
#' \code{ADAPTIVE_K_MAX_FACTOR} (4),
#' \code{ADAPTIVE_K_START_MIN_ABS} (10) and
#' \code{MAX_CANDIDATE_PAIRS} (16).
#' The values strictly descend without deduplication. At least one candidate is always
#' returned, even for data smaller than `k_start_1`, whose construction then fails with a
#' too-few-means status in the search.
#'
#' Not checked here: the values fit 32-bit integers only while `k_start_1 <= huge(1)/4`;
#' \code{\link{calc_adaptive_js_comp_test_bounds}}
#' rejects inputs beyond that.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::generate_adaptive_js_comp_test_candidates}, whose argument names
#' are the ones an error message reports.
#'
#' @param max_n_genes_all_studies a integer scalar. Maximum number of genes across all studies
#'   The minimum valid value is `1`.
#' @param n_studies a integer scalar. Number of studies
#'   The minimum valid value is `1`.
#' @return a integer matrix. Candidate `[k_start, k_step, k_max]` triples, `k_start` strictly descending
#'   The first `n_candidates` elements will hold the results.
#' @export
generate_adaptive_js_comp_test_candidates <- function(max_n_genes_all_studies, n_studies) {
    max_n_genes_all_studies <- .tox_as_integer_scalar(max_n_genes_all_studies, "max_n_genes_all_studies")
    n_studies <- .tox_as_integer_scalar(n_studies, "n_studies")
    .result <- .Call("generate_adaptive_js_comp_test_candidates_call", max_n_genes_all_studies, n_studies)
    .arguments <- c("max_n_genes_all_studies", "n_studies", "candidates_k_start_k_step_k_max", "n_candidates", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$candidates_k_start_k_step_k_max[, seq_len(.result$n_candidates), drop = FALSE]
}

#' Grow adaptive, dispersion-controlled neighborhoods over all studies' pooled gene means (Issue #217)
#'
#' Instead of one fixed neighbor count, each neighborhood grows while the mean absolute residual
#' of its pool stays stable, and the next reference point is placed beyond it; the number of
#' reference points emerges. `k_start`, `k_step` and `k_max` count pooled entries (one entry is
#' one gene of one study).
#'
#' **Pooling.** All studies' gene means are pooled into `N = max_n_genes_all_studies *
#' n_studies` entries (entry `(s - 1)*G + g` is gene `g` of study `s`) and sorted ascending; NaN
#' means take no part. A neighborhood is a contiguous run `a..b` of that order
#' (`pooled_neighborhood_range`, the input of
#' \code{\link{check_neighborhood_overlaps}}),
#' and its per-study counts `n_neighbors_per_point` add up to `b - a + 1`.
#'
#' **Growth** around the reference point `x_star` (the mean at the seed position): entries join
#' nearest first, ties to the lower mean; the dispersion `S` is the mean absolute non-NaN
#' residual of the pool.
#'
#' 1. The first `k_start` entries join unconditionally. With fewer than
#' \code{ADAPTIVE_MIN_VALID_RESIDUALS}
#' (`10`) non-NaN residuals the neighborhood stops with
#' \code{ADAPTIVE_STOP_TOO_FEW_RESIDUALS}
#' (`6`); with all residuals zero, with
#' \code{ADAPTIVE_STOP_ZERO_DISPERSION}
#' (`5`).
#' 2. Then rounds of the next `min(k_step, k_max - count)` entries: a round with
#' `(S_new - S_old) / S_old > tau` is discarded and growth stops
#' (\code{ADAPTIVE_STOP_TAU}, `1`);
#' otherwise it is committed. Growth also stops at `k_max` entries
#' (\code{ADAPTIVE_STOP_K_MAX}, `2`,
#' checked first) or when no entry is left
#' (\code{ADAPTIVE_STOP_EXHAUSTED}, `4`).
#' 3. With `max_pooled_residuals > 0`, a round stages whole entries only while the pool's
#' non-NaN residuals stay within that cap (the `k_start` phase is not capped); a capped round
#' is still evaluated, and growth then stops with
#' \code{ADAPTIVE_STOP_RESIDUAL_CAP} (`3`).
#'
#' **Reference points.** The first seed is the smallest pooled mean. After a neighborhood ends at
#' sorted position `b`, the next seed is the first position after `b` whose mean is nearest
#' `mean(b) + mad_distance_factor * MAD` (ties to the lower mean), `MAD` being the raw median
#' absolute deviation of the neighborhood's pooled means. The construction ends with the
#' neighborhood that reaches the largest pooled mean.
#'
#' **Construction status.**
#' \code{ADAPTIVE_STATUS_TOO_FEW_MEANS}
#' (`1`) when fewer than `k_start` pooled means are non-NaN: nothing
#' is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction always completes,
#' and reports, in this priority,
#' \code{ADAPTIVE_STATUS_TOO_FEW_RESIDUALS}
#' (`3`) if any neighborhood stopped with too few residuals,
#' \code{ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD}
#' (`2`) if any has fewer than `min_study_neighbors`
#' entries of some study, else
#' \code{ADAPTIVE_STATUS_OK} (`0`).
#'
#' Outputs sized by `N` are filled only for their first `n_points` reference points.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::construct_adaptive_neighborhoods}, whose argument names
#' are the ones an error message reports.
#'
#' @param gene_means a numeric matrix. Mean expression of every gene in every study, NaN for a missing gene
#'   NaN is permitted for this value.
#' @param residuals a numeric array of rank 3. Signed residuals of every replicate of every gene in every study, NaN for a missing value
#'   NaN is permitted for this value.
#' @param k_start a integer scalar. Pooled entries every neighborhood takes unconditionally
#'   The minimum valid value is `1`.
#' @param k_step a integer scalar. Pooled entries staged per adaptive growth round
#'   The minimum valid value is `1`.
#' @param k_max a integer scalar. Largest number of pooled entries in one neighborhood
#'   The minimum valid value is `k_start`.
#' @param tau a numeric scalar. Largest relative dispersion increase a growth round may cause and still be committed
#'   The minimum valid value is `0.0`.
#'   The default value is `0.1`.
#' @param mad_distance_factor a numeric scalar. Multiple of a neighborhood's MAD that the next reference point's target lies beyond
#'   the neighborhood's largest mean
#'   The minimum valid value is `0.0`.
#'   The default value is `1.0`.
#' @param max_pooled_residuals a integer scalar. Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
#'   The minimum valid value is `0`.
#'   The default value is `0`.
#' @param min_study_neighbors a integer scalar. Fewest entries of each study every neighborhood must have for the construction
#'   status to stay ok
#'   The minimum valid value is `1`.
#'   The default value is `1`.
#' @return a named list with elements:
#'   \item{x_star}{a numeric vector. Reference point of each neighborhood: the pooled mean at its seed position
#'     The first `n_points` elements will hold the results.}
#'   \item{pooled_neighborhood_range}{a integer matrix. For each neighborhood, its first and last position in the ascending order of the
#'     pooled means
#'     The first `n_points` elements will hold the results.}
#'   \item{n_neighbors_per_point}{a integer matrix. For each neighborhood, how many of its pooled entries belong to each study
#'     The first `n_points` elements will hold the results.}
#'   \item{stop_reason}{a integer vector. Why each neighborhood stopped growing: `1` tau exceeded,
#'     `2` `k_max` reached, `3` residual cap
#'     reached, `4` pool exhausted, `5`
#'     zero dispersion, `6` too few residuals
#'     The first `n_points` elements will hold the results.}
#'   \item{neighborhood_dispersion}{a numeric vector. Mean absolute non-NaN residual of each final neighborhood; NaN without any
#'     The first `n_points` elements will hold the results.}
#'   \item{neighborhood_mad}{a numeric vector. Raw median absolute deviation of each final neighborhood's pooled means
#'     The first `n_points` elements will hold the results.}
#'   \item{max_n_neighbors}{a integer scalar. Largest number of entries one study has in one neighborhood; 0 without a neighborhood}
#'   \item{construction_status}{a integer scalar. Overall outcome: `0` ok, `1` too few
#'     means, `2` a study neighborhood below
#'     `min_study_neighbors`, `3` a neighborhood with too
#'     few residuals}
#' @export
construct_adaptive_neighborhoods <- function(gene_means, residuals, k_start, k_step, k_max, tau = 0.1, mad_distance_factor = 1.0, max_pooled_residuals = 0L, min_study_neighbors = 1L) {
    gene_means <- .tox_as_double_matrix(gene_means, "gene_means")
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    k_start <- .tox_as_integer_scalar(k_start, "k_start")
    k_step <- .tox_as_integer_scalar(k_step, "k_step")
    k_max <- .tox_as_integer_scalar(k_max, "k_max")
    tau <- .tox_as_double_scalar(tau, "tau")
    mad_distance_factor <- .tox_as_double_scalar(mad_distance_factor, "mad_distance_factor")
    max_pooled_residuals <- .tox_as_integer_scalar(max_pooled_residuals, "max_pooled_residuals")
    min_study_neighbors <- .tox_as_integer_scalar(min_study_neighbors, "min_study_neighbors")
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("construct_adaptive_neighborhoods_call", gene_means, residuals, k_start, k_step, k_max, tau, mad_distance_factor, max_pooled_residuals, min_study_neighbors)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "k_start", "k_step", "k_max", "n_points", "x_star", "pooled_neighborhood_range", "n_neighbors_per_point", "stop_reason", "neighborhood_dispersion", "neighborhood_mad", "max_n_neighbors", "construction_status", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        x_star = utils::head(.result$x_star, .result$n_points),
        pooled_neighborhood_range = .result$pooled_neighborhood_range[, seq_len(.result$n_points), drop = FALSE],
        n_neighbors_per_point = .result$n_neighbors_per_point[, seq_len(.result$n_points), drop = FALSE],
        stop_reason = utils::head(.result$stop_reason, .result$n_points),
        neighborhood_dispersion = utils::head(.result$neighborhood_dispersion, .result$n_points),
        neighborhood_mad = utils::head(.result$neighborhood_mad, .result$n_points),
        max_n_neighbors = .result$max_n_neighbors,
        construction_status = .result$construction_status
    )
}

#' Test whether every pair of consecutive neighborhoods overlaps by at least a minimum fraction
#'
#' The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
#' neighborhood spans of
#' \code{\link{construct_neighborhoods_ranged}}
#' (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
#' number of positions they share divided by the length of the earlier span; a reversed span
#' (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
#' `test_`-prefixed export as a test case.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::check_neighborhood_overlaps}, whose argument names
#' are the ones an error message reports.
#'
#' @param neighborhood_range a integer matrix. For each reference point, the `[min_idx, max_idx]` neighborhood span
#'   The minimum valid value is `1`.
#' @param min_neighbor_overlap a numeric scalar. Minimum fractional overlap two consecutive neighborhoods must have
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @return a logical scalar. `TRUE` if every pair of consecutive neighborhoods overlaps by at least
#'   `min_neighbor_overlap`
#' @export
check_neighborhood_overlaps <- function(neighborhood_range, min_neighbor_overlap) {
    neighborhood_range <- .tox_as_integer_matrix(neighborhood_range, "neighborhood_range")
    min_neighbor_overlap <- .tox_as_double_scalar(min_neighbor_overlap, "min_neighbor_overlap")
    .result <- .Call("check_neighborhood_overlaps_call", neighborhood_range, min_neighbor_overlap)
    .arguments <- c("neighborhood_range", "n_points", "min_neighbor_overlap", "all_have_min_neighbor_overlap", "ierr")
    .sources <- c(NA_character_, "neighborhood_range", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$all_have_min_neighbor_overlap
}

#' Test whether every bin of a mean pmf reaches a minimum absolute count
#'
#' The second admissibility gate of a parameter-search candidate, checked after
#' \code{\link{check_neighborhood_overlaps}}.
#' Each point is checked over its own `n_bins_per_point` bins only, so the zero padding beyond
#' them does not count as a failure.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::check_mean_pmf_min_counts}, whose argument names
#' are the ones an error message reports.
#'
#' @param mean_pmf_counts a integer matrix. Absolute counts of a residual per bin for the mean pmf
#'   The minimum valid value is `0`.
#' @param n_bins_per_point a integer vector. Each point's own bin count; only `mean_pmf_counts(i_point, 1:n_bins_per_point(i_point))`
#'   is inspected
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_bins`.
#' @param min_count a integer scalar. Minimum count each bin of the mean pmf must reach
#'   The minimum valid value is `0`.
#' @return a logical scalar. `TRUE` if every bin within each reference point's own `n_bins_per_point`, at every
#'   reference point, reaches at least `min_count`
#' @export
check_mean_pmf_min_counts <- function(mean_pmf_counts, n_bins_per_point, min_count) {
    mean_pmf_counts <- .tox_as_integer_matrix(mean_pmf_counts, "mean_pmf_counts")
    n_bins_per_point <- .tox_as_integer_vector(n_bins_per_point, "n_bins_per_point")
    min_count <- .tox_as_integer_scalar(min_count, "min_count")
    if (length(n_bins_per_point) != dim(mean_pmf_counts)[1])
        .tox_shape_error("n_bins_per_point", length(n_bins_per_point), "mean_pmf_counts", dim(mean_pmf_counts)[1])

    .result <- .Call("check_mean_pmf_min_counts_call", mean_pmf_counts, n_bins_per_point, min_count)
    .arguments <- c("mean_pmf_counts", "n_bins", "n_bins_per_point", "n_points", "min_count", "all_bins_have_min_count", "ierr")
    .sources <- c(NA_character_, "mean_pmf_counts", NA_character_, "mean_pmf_counts", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$all_bins_have_min_count
}

#' Test one candidate pair's bootstrapped confidence intervals against the current best, and detect a JSD plateau
#'
#' For each study, the overlap of the new candidate's interval with the best one so far is
#' measured as a fraction of the new interval's width, and counted when it reaches
#' `succeeding_ci_overlap`. A plateau is found when fewer studies pass than for the current
#' best (the best is then kept), or, after the new candidate became the best, when the count
#' meets the condition `join_method` names.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::check_plateau_condition}, whose argument names
#' are the ones an error message reports.
#'
#' @param confidence_interval a numeric matrix. Bootstrapped JSD confidence interval `[lower, upper]` of the candidate under test
#' @param best_candidate_pair_confidence_interval a numeric matrix. JSD confidence intervals for the current best candidate pair; overwritten with
#'   `confidence_interval` unless the new candidate is worse
#' @param best_candidate_index a integer scalar. Candidate-grid index of the current best candidate pair; overwritten with
#'   `candidate_index` unless the new candidate is worse
#' @param best_exceeded_ci_overlap_count a integer scalar. Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
#'   best candidate pair; overwritten unless the new candidate is worse
#'   The minimum valid value is `0`.
#' @param candidate_index a integer scalar. Candidate-grid index of the candidate pair that produced `confidence_interval`
#'   The minimum valid value is `1`.
#' @param join_method a string, one of "join_min", "join_max", "join_median". How the studies' overlaps are joined into the plateau condition; the median method
#'   requires `count > (n_studies - 1) / 2`
#' @param succeeding_ci_overlap a numeric scalar. Minimum fractional overlap with the best candidate's interval for a study to pass
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @return a named list with elements:
#'   \item{best_candidate_pair_confidence_interval}{a numeric matrix. JSD confidence intervals for the current best candidate pair; overwritten with
#'     `confidence_interval` unless the new candidate is worse}
#'   \item{best_candidate_index}{a integer scalar. Candidate-grid index of the current best candidate pair; overwritten with
#'     `candidate_index` unless the new candidate is worse}
#'   \item{best_exceeded_ci_overlap_count}{a integer scalar. Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
#'     best candidate pair; overwritten unless the new candidate is worse
#'     The minimum valid value is `0`.}
#'   \item{ci_overlap}{a numeric vector. Per-study overlap, in `[0, 1]`, of `confidence_interval` with the best interval as
#'     it was on entry, relative to the new interval's width; always filled. Meaningless
#'     for the very first candidate, whose best interval is only a seed}
#'   \item{plateau_found}{a logical scalar. `TRUE` once the new candidate is no better than the previous best, or once
#'     `join_method`'s overlap condition is met by the new candidate}
#' @export
check_plateau_condition <- function(confidence_interval, best_candidate_pair_confidence_interval, best_candidate_index, best_exceeded_ci_overlap_count, candidate_index, join_method, succeeding_ci_overlap) {
    confidence_interval <- .tox_as_double_matrix(confidence_interval, "confidence_interval")
    best_candidate_pair_confidence_interval <- .tox_as_double_matrix(best_candidate_pair_confidence_interval, "best_candidate_pair_confidence_interval")
    best_candidate_index <- .tox_as_integer_scalar(best_candidate_index, "best_candidate_index")
    best_exceeded_ci_overlap_count <- .tox_as_integer_scalar(best_exceeded_ci_overlap_count, "best_exceeded_ci_overlap_count")
    candidate_index <- .tox_as_integer_scalar(candidate_index, "candidate_index")
    join_method <- .tox_as_mode(join_method, "join_method", c("join_min", "join_max", "join_median"))
    succeeding_ci_overlap <- .tox_as_double_scalar(succeeding_ci_overlap, "succeeding_ci_overlap")
    if (dim(best_candidate_pair_confidence_interval)[2] != dim(confidence_interval)[2])
        .tox_shape_error("best_candidate_pair_confidence_interval", dim(best_candidate_pair_confidence_interval)[2], "confidence_interval", dim(confidence_interval)[2])

    .result <- .Call("check_plateau_condition_call", confidence_interval, best_candidate_pair_confidence_interval, best_candidate_index, best_exceeded_ci_overlap_count, candidate_index, join_method, succeeding_ci_overlap)
    .arguments <- c("confidence_interval", "best_candidate_pair_confidence_interval", "n_studies", "best_candidate_index", "best_exceeded_ci_overlap_count", "candidate_index", "join_method", "succeeding_ci_overlap", "ci_overlap", "plateau_found", "ierr")
    .sources <- c(NA_character_, NA_character_, "confidence_interval", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        best_candidate_pair_confidence_interval = .result$best_candidate_pair_confidence_interval,
        best_candidate_index = .result$best_candidate_index,
        best_exceeded_ci_overlap_count = .result$best_exceeded_ci_overlap_count,
        ci_overlap = .result$ci_overlap,
        plateau_found = .result$plateau_found
    )
}

#' Test one candidate's per-study JSD against the previous admissible candidate for a relative-effect-size plateau
#'
#' The relative-effect-size alternative to
#' \code{\link{check_plateau_condition}}:
#' for each study `i`, the relative change of the observed global JSD between successive
#' admissible candidates,
#' `delta(i) = |global_js_divergence(i) - prev_global_js_divergence(i)| / max(prev_global_js_divergence(i),
#' delta_epsilon)`, is summarized by its median (`delta_median`) and maximum (`delta_max`). A
#' plateau is found once both stay below their thresholds for
#' `delta_min_consecutive_transitions` consecutive transitions, counted across calls in
#' `n_consecutive_ok`.
#'
#' For the first admissible candidate (`has_previous = FALSE`) there is no transition:
#' `delta`, `delta_median` and `delta_max` are `-1.0` ("not computed").
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::check_effect_size_plateau_condition}, whose argument names
#' are the ones an error message reports.
#'
#' @param global_js_divergence a numeric vector. Current admissible candidate's observed global JSD per study
#'   The minimum valid value is `0.0`.
#' @param prev_global_js_divergence a numeric vector. Previous admissible candidate's observed global JSD per study; ignored when
#'   `has_previous` is `FALSE`
#'   The minimum valid value is `0.0`.
#' @param has_previous a logical scalar. `FALSE` for the first admissible candidate, which has no transition
#' @param delta_median_threshold a numeric scalar. Upper bound on the median relative change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#' @param delta_max_threshold a numeric scalar. Upper bound on the largest relative change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#' @param delta_epsilon a numeric scalar. Small constant preventing division by zero when a study's previous JSD was zero
#'   The minimum valid value is `above(0.0)`.
#' @param delta_min_consecutive_transitions a integer scalar. Number of consecutive qualifying transitions required to declare a plateau
#'   The minimum valid value is `1`.
#' @param n_consecutive_ok a integer scalar. Running count of consecutive qualifying transitions; reset to zero when this one does
#'   not qualify or `has_previous` is `FALSE`
#'   The minimum valid value is `0`.
#' @return a named list with elements:
#'   \item{n_consecutive_ok}{a integer scalar. Running count of consecutive qualifying transitions; reset to zero when this one does
#'     not qualify or `has_previous` is `FALSE`
#'     The minimum valid value is `0`.}
#'   \item{delta}{a numeric vector. Per-study relative JSD change from the previous admissible candidate; `-1.0`
#'     throughout iff `.not. has_previous`}
#'   \item{delta_median}{a numeric scalar. Median of `delta` across studies; `-1.0` iff `.not. has_previous`}
#'   \item{delta_max}{a numeric scalar. Maximum of `delta` across studies; `-1.0` iff `.not. has_previous`}
#'   \item{plateau_found}{a logical scalar. `TRUE` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`}
#' @export
check_effect_size_plateau_condition <- function(global_js_divergence, prev_global_js_divergence, has_previous, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, n_consecutive_ok) {
    global_js_divergence <- .tox_as_double_vector(global_js_divergence, "global_js_divergence")
    prev_global_js_divergence <- .tox_as_double_vector(prev_global_js_divergence, "prev_global_js_divergence")
    has_previous <- .tox_as_logical_scalar(has_previous, "has_previous")
    delta_median_threshold <- .tox_as_double_scalar(delta_median_threshold, "delta_median_threshold")
    delta_max_threshold <- .tox_as_double_scalar(delta_max_threshold, "delta_max_threshold")
    delta_epsilon <- .tox_as_double_scalar(delta_epsilon, "delta_epsilon")
    delta_min_consecutive_transitions <- .tox_as_integer_scalar(delta_min_consecutive_transitions, "delta_min_consecutive_transitions")
    n_consecutive_ok <- .tox_as_integer_scalar(n_consecutive_ok, "n_consecutive_ok")
    if (length(prev_global_js_divergence) != length(global_js_divergence))
        .tox_shape_error("prev_global_js_divergence", length(prev_global_js_divergence), "global_js_divergence", length(global_js_divergence))

    .result <- .Call("check_effect_size_plateau_condition_call", global_js_divergence, prev_global_js_divergence, has_previous, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, n_consecutive_ok)
    .arguments <- c("global_js_divergence", "prev_global_js_divergence", "n_studies", "has_previous", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "n_consecutive_ok", "delta", "delta_median", "delta_max", "plateau_found", "ierr")
    .sources <- c(NA_character_, NA_character_, "global_js_divergence", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        n_consecutive_ok = .result$n_consecutive_ok,
        delta = .result$delta,
        delta_median = .result$delta_median,
        delta_max = .result$delta_max,
        plateau_found = .result$plateau_found
    )
}

#' Build the consensus pmf and its histogram counts from all studies' pmfs
#'
#' `mean_pmf` is the equal-weight mean of all studies' pmfs, `mean_pmf_counts` and
#' `mean_pmf_included_n_reps` the sums of their counts. Known limitation: every study is
#' compared against a consensus that includes itself, not a leave-one-out background.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::create_mean_pmf}, whose argument names
#' are the ones an error message reports.
#'
#' @param pmfs a numeric array of rank 3. Per-study probabilities of each bin per reference point, from
#'   \code{\link{build_residual_histograms}}
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @param counts a integer array of rank 3. Absolute counts of a residual per bin for `pmfs`
#'   The minimum valid value is `0`.
#' @param included_n_reps a integer matrix. Count of non-NaN replicates per reference point, per study; every study needs at
#'   least one at every reference point
#'   The minimum valid value is `1`.
#' @return a named list with elements:
#'   \item{mean_pmf}{a numeric matrix. The consensus pmf, `sum(pmfs, dim=3) / n_studies`}
#'   \item{mean_pmf_included_n_reps}{a integer vector. Count of non-NaN replicates per reference point, summed over the studies}
#'   \item{mean_pmf_counts}{a integer matrix. Residual counts per bin of the consensus, `sum(counts, dim=3)`}
#' @export
create_mean_pmf <- function(pmfs, counts, included_n_reps) {
    pmfs <- .tox_as_double_array(pmfs, "pmfs", 3L)
    counts <- .tox_as_integer_array(counts, "counts", 3L)
    included_n_reps <- .tox_as_integer_matrix(included_n_reps, "included_n_reps")
    if (dim(counts)[2] != dim(pmfs)[2])
        .tox_shape_error("counts", dim(counts)[2], "pmfs", dim(pmfs)[2])
    if (dim(counts)[1] != dim(pmfs)[1])
        .tox_shape_error("counts", dim(counts)[1], "pmfs", dim(pmfs)[1])
    if (dim(included_n_reps)[1] != dim(pmfs)[1])
        .tox_shape_error("included_n_reps", dim(included_n_reps)[1], "pmfs", dim(pmfs)[1])
    if (dim(counts)[3] != dim(pmfs)[3])
        .tox_shape_error("counts", dim(counts)[3], "pmfs", dim(pmfs)[3])
    if (dim(included_n_reps)[2] != dim(pmfs)[3])
        .tox_shape_error("included_n_reps", dim(included_n_reps)[2], "pmfs", dim(pmfs)[3])

    .result <- .Call("create_mean_pmf_call", pmfs, counts, included_n_reps)
    .arguments <- c("pmfs", "counts", "n_bins", "n_points", "n_studies", "included_n_reps", "mean_pmf", "mean_pmf_included_n_reps", "mean_pmf_counts", "ierr")
    .sources <- c(NA_character_, NA_character_, "pmfs", "pmfs", "pmfs", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        mean_pmf = .result$mean_pmf,
        mean_pmf_included_n_reps = .result$mean_pmf_included_n_reps,
        mean_pmf_counts = .result$mean_pmf_counts
    )
}

#' Build only the consensus pmf from all studies' pmfs, without its histogram counts
#'
#' As \code{\link{create_mean_pmf}}, for
#' callers that need no counts (the bootstrap). Same known limitation: the consensus includes
#' the study compared against it.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::create_mean_pmf_only}, whose argument names
#' are the ones an error message reports.
#'
#' @param pmfs a numeric array of rank 3. Per-study probabilities of each bin per reference point, from
#'   \code{\link{build_residual_histograms}}
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @return a numeric matrix. The consensus pmf, `sum(pmfs, dim=3) / n_studies`
#' @export
create_mean_pmf_only <- function(pmfs) {
    pmfs <- .tox_as_double_array(pmfs, "pmfs", 3L)
    .result <- .Call("create_mean_pmf_only_call", pmfs)
    .arguments <- c("pmfs", "n_bins", "n_points", "n_studies", "mean_pmf", "ierr")
    .sources <- c(NA_character_, "pmfs", "pmfs", "pmfs", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$mean_pmf
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::bootstrap_histogram}, whose argument names
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

#' Run the JSD-Comp-Test pipeline for one fixed (n_points, n_neighbors) parameter setting
#'
#' The final JS-Comp-Test of K studies for one fixed `(n_points, n_neighbors)` setting:
#'
#' 1. Each study's `n_neighbors` genes nearest every reference point `x_star` (computed once by
#' the caller, e.g. with
#' \code{\link{pool_means}}) form its
#' neighborhoods
#' (\code{\link{construct_neighborhoods_ranged}}).
#' 2. Per reference point, all studies' residuals are pooled and the bin-count search
#' `bin_count_mode` names chooses the point's bin count and residual range (occupancy search
#' \code{\link{determine_bin_count_occupancy}}
#' by default, or Stone's rule
#' \code{\link{determine_bin_count_stone}}).
#' 3. Each study's histogram and pmf is built per point
#' (\code{\link{build_residual_histograms}}).
#' 4. The consensus pmf
#' (\code{\link{create_mean_pmf}}), each
#' study's per-point JSD against it and its weighted global JSD, and the permutation test
#' (\code{\link{gjct_permutation_test}}) give
#' `p_values` (primary) and `p_values_observed_consensus` (comparison only). With
#' `n_permutations = 0` both p-values are NaN.
#'
#' There is no admissibility gate: a point whose bin-count search fails (`occupancy_failed`)
#' still gets a histogram at `m_min` bins and still contributes to `global_js_divergence`. A
#' caller that needs such neighborhoods rejected should use
#' \code{\link{run_js_comp_test_parameter_search}}.
#'
#' The bin-sized arrays (`pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`) have the fixed
#' second extent \code{MAX_N_BINS}
#' (`256`); only their first `max_n_bins_per_point` bins are meaningful, the rest is
#' zero. A Python/R caller slices `[:, :max_n_bins_per_point, ...]` (R:
#' `[, seq_len(max_n_bins_per_point), ]`) itself.
#'
#' Errors (`ierr`): `ERR_INVALID_INPUT` when some study has no non-NaN residual at some
#' reference point (for example a study whose gene means are all missing), since the
#' equal-weight consensus needs every study; the consensus, JSD, weight and p-value outputs are
#' then not filled. `ERR_ALLOC_FAIL` when the permutation test cannot create its random number
#' generator; `pmfs`, `js_divergences`, `weights` and `global_js_divergence` then hold their
#' values from before the permutation test. An out-of-range neighbor index is also reported
#' as `ERR_INVALID_INPUT`, a safeguard valid input cannot reach.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::run_js_comp_test}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_neighbors a integer scalar. Number of neighbors per neighborhood
#'   The minimum valid value is `1`.
#'   The maximum valid value is `max_n_genes_all_studies`.
#' @param gene_means a numeric matrix. Per-gene mean expression values for all studies
#'   NaN is permitted for this value.
#' @param gene_means_perms a integer matrix. Per-study sorting permutation for `gene_means` (ascending, NaN last)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `max_n_genes_all_studies`.
#' @param residuals a numeric array of rank 3. Matrix of signed residuals per study
#'   NaN is permitted for this value.
#' @param x_star a numeric vector. Mean-expression reference points
#' @param n_permutations a integer scalar. Number of permutations of the permutation test; `0` skips it (p-values NaN)
#'   The minimum valid value is `0`.
#'   The default value is `1000`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach for a bin count to be
#'   admissible, in every `bin_count_mode`
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param m_min a integer scalar. Smallest bin count the bin-count search tests (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor of the occupancy search; unused in the Stone modes
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's lower residual-range bound
#'   Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
#'   with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's upper residual-range bound
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen; the Stone modes test every
#'   bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @return a named list with elements:
#'   \item{neighborhood_indices}{a integer array of rank 3. Gene indices of the selected neighborhood, per reference point, per study}
#'   \item{neighborhood_range}{a integer array of rank 3. For each reference point and study, the `[min_idx, max_idx]` neighborhood span}
#'   \item{n_bins_per_point}{a integer vector. Each reference point's histogram bin count `M_j`, from the bin-count search}
#'   \item{shared_residual_range_low}{a numeric vector. Each reference point's lower residual-range bound (R_low), from the bin-count search}
#'   \item{shared_residual_range_high}{a numeric vector. Each reference point's upper residual-range bound (R_high), from the bin-count search}
#'   \item{max_n_bins_per_point}{a integer scalar. `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
#'     `mean_pmf` and `mean_pmf_counts` (see above)}
#'   \item{occupancy_failed}{a logical vector. `TRUE` iff the bin-count search failed for this point; it then uses `m_min` bins and
#'     still contributes (see above)}
#'   \item{n_pooled_residuals}{a integer vector. Each reference point's non-NaN pooled residual count (N_j)}
#'   \item{min_bin_occupancy}{a integer vector. Each reference point's minimum bin occupancy at `n_bins_per_point`}
#'   \item{mean_bin_occupancy}{a numeric vector. Each reference point's mean bin occupancy at `n_bins_per_point`}
#'   \item{max_bin_occupancy}{a integer vector. Each reference point's maximum bin occupancy at `n_bins_per_point`}
#'   \item{sturges_bins}{a integer vector. Each reference point's Sturges' rule diagnostic (not used for the decision)}
#'   \item{fd_bins}{a integer vector. Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)}
#'   \item{pmfs}{a numeric array of rank 3. `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
#'     are meaningful}
#'   \item{counts}{a integer array of rank 3. Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful}
#'   \item{included_n_reps}{a integer matrix. Count of non-NaN replicates per reference point, per study}
#'   \item{mean_pmf}{a numeric matrix. The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful}
#'   \item{mean_pmf_counts}{a integer matrix. Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
#'     meaningful}
#'   \item{mean_pmf_included_n_reps}{a integer vector. Count of non-NaN replicates per reference point of the consensus pmf}
#'   \item{js_divergences}{a numeric matrix. Per-reference-point JSD of each study against the consensus pmf}
#'   \item{weights}{a numeric matrix. Per-reference-point weights for `global_js_divergence`}
#'   \item{global_js_divergence}{a numeric vector. Weighted global JSD of each study against the consensus pmf}
#'   \item{p_values}{a numeric vector. Empirical p-value per study, the primary result; NaN when `n_permutations` is 0}
#'   \item{p_values_observed_consensus}{a numeric vector. Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
#'     comparison only; NaN when `n_permutations` is 0}
#' @export
run_js_comp_test <- function(n_neighbors, gene_means, gene_means_perms, residuals, x_star, n_permutations = 1000L, random_seed = 42L, min_residuals_per_bin = 10L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, bin_count_mode = "bin_count_occupancy") {
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    gene_means <- .tox_as_double_matrix(gene_means, "gene_means")
    gene_means_perms <- .tox_as_integer_matrix(gene_means_perms, "gene_means_perms")
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    x_star <- .tox_as_double_vector(x_star, "x_star")
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    gamma_occupancy <- .tox_as_double_scalar(gamma_occupancy, "gamma_occupancy")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    bin_count_mode <- .tox_as_mode(bin_count_mode, "bin_count_mode", c("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled"))
    if (dim(gene_means_perms)[2] != dim(gene_means)[2])
        .tox_shape_error("gene_means_perms", dim(gene_means_perms)[2], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(gene_means_perms)[1] != dim(gene_means)[1])
        .tox_shape_error("gene_means_perms", dim(gene_means_perms)[1], "gene_means", dim(gene_means)[1])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_call", n_neighbors, gene_means, gene_means_perms, residuals, x_star, n_permutations, random_seed, min_residuals_per_bin, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, bin_count_mode)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "n_neighbors", "gene_means", "gene_means_perms", "residuals", "x_star", "neighborhood_indices", "neighborhood_range", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "p_values_observed_consensus", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "bin_count_mode", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", "x_star", "neighborhood_indices", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        neighborhood_indices = .result$neighborhood_indices,
        neighborhood_range = .result$neighborhood_range,
        n_bins_per_point = .result$n_bins_per_point,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        max_n_bins_per_point = .result$max_n_bins_per_point,
        occupancy_failed = .result$occupancy_failed,
        n_pooled_residuals = .result$n_pooled_residuals,
        min_bin_occupancy = .result$min_bin_occupancy,
        mean_bin_occupancy = .result$mean_bin_occupancy,
        max_bin_occupancy = .result$max_bin_occupancy,
        sturges_bins = .result$sturges_bins,
        fd_bins = .result$fd_bins,
        pmfs = .result$pmfs,
        counts = .result$counts,
        included_n_reps = .result$included_n_reps,
        mean_pmf = .result$mean_pmf,
        mean_pmf_counts = .result$mean_pmf_counts,
        mean_pmf_included_n_reps = .result$mean_pmf_included_n_reps,
        js_divergences = .result$js_divergences,
        weights = .result$weights,
        global_js_divergence = .result$global_js_divergence,
        p_values = .result$p_values,
        p_values_observed_consensus = .result$p_values_observed_consensus
    )
}

#' Run the JSD-Comp-Test pipeline on adaptive (Issue #217) neighborhoods
#'
#' The adaptive counterpart of
#' \code{\link{run_js_comp_test}}: the
#' same pipeline and outputs, on neighborhoods whose size varies per reference point and per
#' study. Call
#' \code{\link{construct_adaptive_neighborhoods}}
#' first and pass its `pooled_neighborhood_range` (its first `n_points` columns) with the same
#' `gene_means` and `residuals`. Each range is decoded back to per-study gene lists; their sizes
#' are returned as `n_neighbors_per_point`. Neighborhoods that decode to the same gene sets as a
#' fixed-k run give bit-identical results. There is no admissibility gate.
#'
#' Errors (`ierr`): `ERR_INVALID_INPUT` before any histogram work for the first range that is
#' reversed, reaches past the non-NaN pooled means, or leaves a study without a gene (every
#' output then undefined); otherwise as in `run_js_comp_test_impl` (a study without a non-NaN
#' residual at some point, or `ERR_ALLOC_FAIL` from the permutation test).
#'
#' The pooling work arrays hold `max_n_reps_all_studies * max_n_genes_all_studies * n_studies`
#' values each, the size of `residuals`, since one study may contribute all its genes to one
#' neighborhood. The bin-sized arrays have the fixed second extent
#' \code{MAX_N_BINS}; only their first
#' `max_n_bins_per_point` bins are meaningful (Python `[:, :max_n_bins_per_point, ...]`, R
#' `[, seq_len(max_n_bins_per_point), ]`).
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::run_js_comp_test_adaptive}, whose argument names
#' are the ones an error message reports.
#'
#' @param gene_means a numeric matrix. Mean expression of every gene in every study, NaN for a missing gene; the array the
#'   neighborhoods were constructed from
#'   NaN is permitted for this value.
#' @param residuals a numeric array of rank 3. Signed residuals of every replicate of every gene in every study, NaN for a missing value
#'   NaN is permitted for this value.
#' @param pooled_neighborhood_range a integer matrix. For each neighborhood, its first and last sorted pooled position, from
#'   construct_adaptive_neighborhoods
#'   The minimum valid value is `1`.
#'   The maximum valid value is `max_n_genes_all_studies*n_studies`.
#' @param n_permutations a integer scalar. Number of permutations of the permutation test; `0` skips it (p-values NaN)
#'   The minimum valid value is `0`.
#'   The default value is `1000`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach for a bin count to be
#'   admissible, in every `bin_count_mode`
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param m_min a integer scalar. Smallest bin count the bin-count search tests (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor of the occupancy search; unused in the Stone modes
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's lower residual-range bound
#'   Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
#'   with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's upper residual-range bound
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen; the Stone modes test every
#'   bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @return a named list with elements:
#'   \item{n_neighbors_per_point}{a integer matrix. For each neighborhood, how many of its genes belong to each study (at least 1 each)}
#'   \item{n_bins_per_point}{a integer vector. Each reference point's selected histogram bin count, from the bin-count search}
#'   \item{shared_residual_range_low}{a numeric vector. Each reference point's lower residual-range bound (R_low), from the bin-count search}
#'   \item{shared_residual_range_high}{a numeric vector. Each reference point's upper residual-range bound (R_high), from the bin-count search}
#'   \item{max_n_bins_per_point}{a integer scalar. `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
#'     `mean_pmf` and `mean_pmf_counts`}
#'   \item{occupancy_failed}{a logical vector. `TRUE` iff the bin-count search failed for this point; it then uses `m_min` bins and
#'     still contributes}
#'   \item{n_pooled_residuals}{a integer vector. Each reference point's non-NaN pooled residual count (N_j)}
#'   \item{min_bin_occupancy}{a integer vector. Each reference point's minimum bin occupancy at `n_bins_per_point`}
#'   \item{mean_bin_occupancy}{a numeric vector. Each reference point's mean bin occupancy at `n_bins_per_point`}
#'   \item{max_bin_occupancy}{a integer vector. Each reference point's maximum bin occupancy at `n_bins_per_point`}
#'   \item{sturges_bins}{a integer vector. Each reference point's Sturges' rule diagnostic (not used for the decision)}
#'   \item{fd_bins}{a integer vector. Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)}
#'   \item{pmfs}{a numeric array of rank 3. `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
#'     are meaningful}
#'   \item{counts}{a integer array of rank 3. Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful}
#'   \item{included_n_reps}{a integer matrix. Count of non-NaN residuals binned per reference point, per study}
#'   \item{mean_pmf}{a numeric matrix. The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful}
#'   \item{mean_pmf_counts}{a integer matrix. Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
#'     meaningful}
#'   \item{mean_pmf_included_n_reps}{a integer vector. Count of non-NaN residuals per reference point for the consensus pmf}
#'   \item{js_divergences}{a numeric matrix. Per-reference-point JSD of each study against the consensus pmf}
#'   \item{weights}{a numeric matrix. Per-reference-point weights for `global_js_divergence`}
#'   \item{global_js_divergence}{a numeric vector. Weighted global JSD of each study against the consensus pmf}
#'   \item{p_values}{a numeric vector. Empirical p-value per study, the primary result; NaN when `n_permutations` is 0}
#'   \item{p_values_observed_consensus}{a numeric vector. Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
#'     comparison only; NaN when `n_permutations` is 0}
#' @export
run_js_comp_test_adaptive <- function(gene_means, residuals, pooled_neighborhood_range, n_permutations = 1000L, random_seed = 42L, min_residuals_per_bin = 10L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, bin_count_mode = "bin_count_occupancy") {
    gene_means <- .tox_as_double_matrix(gene_means, "gene_means")
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    pooled_neighborhood_range <- .tox_as_integer_matrix(pooled_neighborhood_range, "pooled_neighborhood_range")
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    gamma_occupancy <- .tox_as_double_scalar(gamma_occupancy, "gamma_occupancy")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    bin_count_mode <- .tox_as_mode(bin_count_mode, "bin_count_mode", c("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled"))
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_adaptive_call", gene_means, residuals, pooled_neighborhood_range, n_permutations, random_seed, min_residuals_per_bin, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, bin_count_mode)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "gene_means", "residuals", "pooled_neighborhood_range", "n_neighbors_per_point", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "p_values_observed_consensus", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "bin_count_mode", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", "pooled_neighborhood_range", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        n_neighbors_per_point = .result$n_neighbors_per_point,
        n_bins_per_point = .result$n_bins_per_point,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        max_n_bins_per_point = .result$max_n_bins_per_point,
        occupancy_failed = .result$occupancy_failed,
        n_pooled_residuals = .result$n_pooled_residuals,
        min_bin_occupancy = .result$min_bin_occupancy,
        mean_bin_occupancy = .result$mean_bin_occupancy,
        max_bin_occupancy = .result$max_bin_occupancy,
        sturges_bins = .result$sturges_bins,
        fd_bins = .result$fd_bins,
        pmfs = .result$pmfs,
        counts = .result$counts,
        included_n_reps = .result$included_n_reps,
        mean_pmf = .result$mean_pmf,
        mean_pmf_counts = .result$mean_pmf_counts,
        mean_pmf_included_n_reps = .result$mean_pmf_included_n_reps,
        js_divergences = .result$js_divergences,
        weights = .result$weights,
        global_js_divergence = .result$global_js_divergence,
        p_values = .result$p_values,
        p_values_observed_consensus = .result$p_values_observed_consensus
    )
}

#' Search a GAMMA-decay (n_points, n_neighbors) candidate grid for a stable JSD parameter setting
#'
#' Walks the candidates of
#' \code{\link{generate_js_comp_test_candidates}}
#' from finest to coarsest resolution. For each candidate:
#'
#' 1. Reference points from the pooled gene means
#' (\code{\link{pool_means}}) and every
#' study's neighborhoods; gate 1: consecutive neighborhoods overlap enough
#' (\code{\link{check_neighborhood_overlaps}}).
#' 2. Per reference point, the bin-count search `bin_count_mode` names on all studies' pooled
#' residuals, then every study's histograms.
#' 3. Gate 2: every study has a non-NaN residual at every point, and every consensus bin reaches
#' `min_residuals_per_bin`
#' (\code{\link{check_mean_pmf_min_counts}}).
#' A candidate that passes both gates is *admissible*.
#' 4. The observed global JSD per study, its bootstrap confidence interval
#' (\code{\link{bootstrap_histogram}}),
#' and both plateau criteria: CI overlap with the best candidate so far
#' (\code{\link{check_plateau_condition}})
#' and relative effect-size stability
#' (\code{\link{check_effect_size_plateau_condition}}).
#' Both are always computed and traced; `plateau_mode` decides which one stops the search. A
#' stop by the effect-size criterion returns the candidate that triggered it.
#'
#' The returned candidate:
#'
#' 1. A plateau was found, or the grid has a single candidate and it is admissible: the best
#' candidate, `plateau_established = TRUE`.
#' 2. No plateau, but some candidate was admissible: the admissible candidate with the smallest
#' median confidence-interval width across studies, `plateau_established = FALSE`.
#' 3. No candidate was admissible: the first (finest) candidate, with
#' `best_candidate_pair_confidence_interval = -1.0`, every bin count `m_min`, both residual
#' ranges `0.0`, and `plateau_established = FALSE`.
#'
#' Per-point outputs and `trace_*` arrays are sized by the candidate bounds; see each argument
#' for its meaningful part. A failure to create the bootstrap's random number generator is
#' returned in `ierr` and stops the search; the outputs are then undefined.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::run_js_comp_test_parameter_search}, whose argument names
#' are the ones an error message reports.
#'
#' @param gene_means a numeric matrix. Per-gene mean expression values for all studies
#'   NaN is permitted for this value.
#' @param residuals a numeric array of rank 3. Matrix of signed residuals per study
#'   NaN is permitted for this value.
#' @param n_bootstraps a integer scalar. Number of bootstraps to perform for a candidate pair
#'   The minimum valid value is `1`.
#' @param join_method a string, one of "join_min", "join_max", "join_median". How the studies' confidence-interval overlaps are joined into the plateau condition
#' @param max_n_points_candidate a integer scalar. Upper bound on the grid's `n_points` candidates, which sizes the per-point outputs;
#'   at least what `calc_js_comp_test_candidate_bounds` returns, else invalid input
#'   It is recommended to compute this argument from the `max_n_points_candidate` output produced by \code{\link{calc_js_comp_test_candidate_bounds}}.
#'   The minimum valid value is `max(MIN_POINTS, min(ceiling(4.0*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS))`.
#' @param max_n_neighbors_candidate a integer scalar. Upper bound on the grid's `n_neighbors` candidates; at least what
#'   `calc_js_comp_test_candidate_bounds` returns, else invalid input
#'   It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by \code{\link{calc_js_comp_test_candidate_bounds}}.
#'   The minimum valid value is `max(1, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2*real(max(MIN_POINTS, min(ceiling(4.0*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32))`.
#' @param min_residuals_per_bin a integer scalar. Minimum count of every consensus bin (gate 2) and the floor of the bin-count search,
#'   one threshold for both, so a point the search marks as failed also fails the gate
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param min_neighbor_overlap a numeric scalar. Minimum fractional overlap of two consecutive neighborhoods (gate 1)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.1`.
#' @param succeeding_ci_overlap a numeric scalar. Minimum fractional overlap of a candidate's confidence interval with the best one so
#'   far, per `join_method`. `0` makes the first admissible candidate plateau under the
#'   CI-overlap criterion (for baseline sweeps)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @param plateau_mode a string, one of "plateau_ci_overlap", "plateau_effect_size", "plateau_both". Which plateau criterion decides when the search stops
#'
#'   The default value is `"plateau_ci_overlap"`.
#' @param delta_median_threshold a numeric scalar. Upper bound on the median relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.05`.
#' @param delta_max_threshold a numeric scalar. Upper bound on the largest relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.10`.
#' @param delta_epsilon a numeric scalar. Floor of the denominators of the relative JSD change and relative CI width
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `1.0e-10`.
#' @param delta_min_consecutive_transitions a integer scalar. Number of consecutive qualifying transitions for an effect-size plateau
#'   The minimum valid value is `1`.
#'   The default value is `2`.
#' @param m_min a integer scalar. Smallest bin count the bin-count search tests (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor of the occupancy search; unused in the Stone modes
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's lower residual-range bound
#'   Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
#'   with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's upper residual-range bound
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Two-sided significance level of the bootstrap interval, in percent per end
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen, for every candidate; the
#'   Stone modes test every bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @return a named list with elements:
#'   \item{n_points}{a integer scalar. The finally chosen candidate's `n_points`}
#'   \item{n_neighbors}{a integer scalar. The finally chosen candidate's `n_neighbors`}
#'   \item{n_bins_per_point}{a integer vector. The chosen candidate's bin count per reference point; the first `n_points` entries
#'     are meaningful}
#'   \item{shared_residual_range_low}{a numeric vector. The chosen candidate's lower residual-range bound (R_low) per reference point; the
#'     first `n_points` entries are meaningful, `0.0` when no candidate was admissible}
#'   \item{shared_residual_range_high}{a numeric vector. The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`}
#'   \item{best_candidate_pair_confidence_interval}{a numeric matrix. The chosen candidate's bootstrapped JSD confidence interval; `-1.0` when no candidate
#'     was admissible}
#'   \item{plateau_established}{a logical scalar. `TRUE` when a plateau was found (by the criterion `plateau_mode` selects) or the
#'     only candidate was admissible; `FALSE` otherwise, also when the smallest-uncertainty
#'     fallback returns a real candidate and interval}
#'   \item{trace_n_points}{a integer vector. Per-admissible-candidate `n_points`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_n_neighbors}{a integer vector. Per-admissible-candidate `n_neighbors`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_global_js_divergence}{a numeric matrix. Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_lower}{a numeric matrix. Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
#'     (`L_{i,t}`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_upper}{a numeric matrix. Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
#'     (`U_{i,t}`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_width}{a numeric matrix. Per-admissible-candidate, per-study confidence-interval width (`W_{i,t} = U_{i,t} -
#'     L_{i,t}`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_width_relative}{a numeric matrix. Per-admissible-candidate, per-study relative confidence-interval width
#'     (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_overlap}{a numeric matrix. Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
#'     interval with the best candidate's so far, relative to its own width (the values
#'     `succeeding_ci_overlap` is tested against); `-1.0` for the first admissible candidate
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta}{a numeric matrix. Per-admissible-candidate, per-study relative JSD change from the previous admissible
#'     candidate (`Delta_{i,t}`); `-1.0` for the first admissible candidate
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta_median}{a numeric vector. Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta_max}{a numeric vector. Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_selected_n_bins}{a integer matrix. Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
#'     `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
#'     `[:trace_n_points[t], t]`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_occupancy_failed}{a logical matrix. Per-admissible-candidate, per-reference-point flag: `TRUE` iff the bin-count search
#'     failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_n_pooled_residuals}{a integer matrix. Per-admissible-candidate, per-reference-point pooled residual count (`N_j`). Jagged
#'     like trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_min_bin_occupancy}{a integer matrix. Per-admissible-candidate, per-reference-point minimum bin occupancy. Jagged like
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_mean_bin_occupancy}{a numeric matrix. Per-admissible-candidate, per-reference-point mean bin occupancy. Jagged like
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_max_bin_occupancy}{a integer matrix. Per-admissible-candidate, per-reference-point maximum bin occupancy. Jagged like
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_sturges_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Sturges' rule diagnostic. Jagged like
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_fd_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic.
#'     Jagged like trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_low}{a numeric matrix. Per-admissible-candidate, per-reference-point lower residual-range bound (R_low).
#'     Jagged like trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_high}{a numeric matrix. Per-admissible-candidate, per-reference-point upper residual-range bound (R_high).
#'     Jagged like trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#' @export
run_js_comp_test_parameter_search <- function(gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, max_n_neighbors_candidate, min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_ci_overlap", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, bin_count_mode = "bin_count_occupancy") {
    gene_means <- .tox_as_double_matrix(gene_means, "gene_means")
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    join_method <- .tox_as_mode(join_method, "join_method", c("join_min", "join_max", "join_median"))
    max_n_points_candidate <- .tox_as_integer_scalar(max_n_points_candidate, "max_n_points_candidate")
    max_n_neighbors_candidate <- .tox_as_integer_scalar(max_n_neighbors_candidate, "max_n_neighbors_candidate")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    min_neighbor_overlap <- .tox_as_double_scalar(min_neighbor_overlap, "min_neighbor_overlap")
    succeeding_ci_overlap <- .tox_as_double_scalar(succeeding_ci_overlap, "succeeding_ci_overlap")
    plateau_mode <- .tox_as_mode(plateau_mode, "plateau_mode", c("plateau_ci_overlap", "plateau_effect_size", "plateau_both"))
    delta_median_threshold <- .tox_as_double_scalar(delta_median_threshold, "delta_median_threshold")
    delta_max_threshold <- .tox_as_double_scalar(delta_max_threshold, "delta_max_threshold")
    delta_epsilon <- .tox_as_double_scalar(delta_epsilon, "delta_epsilon")
    delta_min_consecutive_transitions <- .tox_as_integer_scalar(delta_min_consecutive_transitions, "delta_min_consecutive_transitions")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    gamma_occupancy <- .tox_as_double_scalar(gamma_occupancy, "gamma_occupancy")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    two_sided_bootstrapping_significance_level <- .tox_as_double_scalar(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    bin_count_mode <- .tox_as_mode(bin_count_mode, "bin_count_mode", c("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled"))
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_parameter_search_call", gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, max_n_neighbors_candidate, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, bin_count_mode)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "n_bootstraps", "join_method", "max_n_points_candidate", "max_n_neighbors_candidate", "n_points", "n_neighbors", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "best_candidate_pair_confidence_interval", "plateau_established", "n_admissible_evaluated", "trace_n_points", "trace_n_neighbors", "trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_ci_overlap", "trace_delta", "trace_delta_median", "trace_delta_max", "trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy", "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins", "trace_shared_residual_range_low", "trace_shared_residual_range_high", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", NA_character_, NA_character_, NA_character_, NA_character_, "n_bins_per_point", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        n_points = .result$n_points,
        n_neighbors = .result$n_neighbors,
        n_bins_per_point = .result$n_bins_per_point,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        best_candidate_pair_confidence_interval = .result$best_candidate_pair_confidence_interval,
        plateau_established = .result$plateau_established,
        trace_n_points = utils::head(.result$trace_n_points, .result$n_admissible_evaluated),
        trace_n_neighbors = utils::head(.result$trace_n_neighbors, .result$n_admissible_evaluated),
        trace_global_js_divergence = .result$trace_global_js_divergence[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_lower = .result$trace_ci_lower[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_upper = .result$trace_ci_upper[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_width = .result$trace_ci_width[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_width_relative = .result$trace_ci_width_relative[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_overlap = .result$trace_ci_overlap[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_delta = .result$trace_delta[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_delta_median = utils::head(.result$trace_delta_median, .result$n_admissible_evaluated),
        trace_delta_max = utils::head(.result$trace_delta_max, .result$n_admissible_evaluated),
        trace_selected_n_bins = .result$trace_selected_n_bins[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_occupancy_failed = .result$trace_occupancy_failed[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_n_pooled_residuals = .result$trace_n_pooled_residuals[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_min_bin_occupancy = .result$trace_min_bin_occupancy[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_mean_bin_occupancy = .result$trace_mean_bin_occupancy[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_max_bin_occupancy = .result$trace_max_bin_occupancy[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_sturges_bins = .result$trace_sturges_bins[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_fd_bins = .result$trace_fd_bins[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_shared_residual_range_low = .result$trace_shared_residual_range_low[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_shared_residual_range_high = .result$trace_shared_residual_range_high[, seq_len(.result$n_admissible_evaluated), drop = FALSE]
    )
}

#' Search the ascending adaptive (k_start, k_step, k_max) candidate sequence for a stable JSD parameter setting (Issue #217)
#'
#' The adaptive-neighborhood counterpart of
#' \code{\link{run_js_comp_test_parameter_search}}.
#' Instead of a fixed `(n_points, n_neighbors)` grid, it walks the growth-knob candidates of
#' \code{\link{generate_adaptive_js_comp_test_candidates}}
#' in order. Their `k_start` shrinks, so the number of reference points that emerges ascends.
#' The plateau criteria, `plateau_mode` and the fallbacks are exactly the fixed-k search's.
#'
#' **Per candidate**, in order:
#'
#' 1. Grow its neighborhoods exactly as
#' \code{\link{construct_adaptive_neighborhoods}}
#' does with the same `tau`, `mad_distance_factor`, `max_pooled_residuals` and
#' `min_study_neighbors`.
#' 2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
#' reference points emerged.
#' 3. Gate 1:
#' \code{\link{check_neighborhood_overlaps}}
#' on the pooled ranges.
#' 4. The pipeline of
#' \code{\link{run_js_comp_test_adaptive}}
#' up to the consensus pmf, then gate 2: every study has a non-NaN residual at every point
#' (otherwise logged as an empty-study candidate), and every consensus bin reaches the minimum
#' count
#' (\code{\link{check_mean_pmf_min_counts}}).
#' 5. As in the fixed-k search: traces, bootstrap, both plateau criteria, and the stop at the
#' first plateau `plateau_mode` selects.
#'
#' **Candidate log.** Every candidate tried gets an entry in `candidate_k_start`,
#' `candidate_n_points` and `candidate_status`, until the search stops. The status is the first
#' failing check, in the order above:
#' \code{ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED}
#' (too few means or residuals),
#' \code{ADAPTIVE_CANDIDATE_EMPTY_STUDY}
#' (a study neighborhood below `min_study_neighbors`, or a study without a non-NaN residual at
#' some point),
#' \code{ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED},
#' \code{ADAPTIVE_CANDIDATE_OVERLAP_FAILED},
#' \code{ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED},
#' or else \code{ADAPTIVE_CANDIDATE_EVALUATED}
#' (admissible, with a trace column).
#'
#' **Result**, one of three cases, as in the fixed-k search:
#'
#' 1. A plateau was found, or the sequence had a single candidate and it passed both gates:
#' that candidate's knobs, `plateau_established = TRUE`.
#' 2. No plateau but at least one admissible candidate: the admissible candidate with the
#' smallest median confidence-interval width, `plateau_established = FALSE`.
#' 3. Nothing admissible, including a single candidate that was rejected: candidate 1's knobs,
#' confidence interval `-1.0`, `plateau_established = FALSE`.
#'
#' `n_points` is the selected candidate's logged point count; in case 3 the bin counts are
#' `m_min` and the ranges `0.0`, as in the fixed-k search. To rebuild the selected neighborhoods, call
#' `construct_adaptive_neighborhoods` with the returned knobs and the same construction
#' options, then `run_js_comp_test_adaptive` with the same bin-count options
#' (`bin_count_mode` included): its `global_js_divergence` equals the selected candidate's
#' `trace_global_js_divergence` column.
#'
#' **Traces** are the fixed-k search's, with `trace_k_start`/`trace_k_step`/`trace_k_max` in
#' place of `trace_n_neighbors`; the per-point traces are jagged in the same way.
#'
#' **Memory.** The per-point arrays are sized by `max_n_points_candidate` (from
#' \code{\link{calc_adaptive_js_comp_test_bounds}});
#' the construction buffers by `N = max_n_genes_all_studies * n_studies` (the largest,
#' `tmp_n_neighbors_per_point`, `n_studies * N` integers); the pooling buffers hold
#' `max_n_reps_all_studies * N` values each, the size of `residuals`.
#'
#' A failure to create the bootstrap's random number generator is returned in `ierr` and
#' stops the search; the outputs are then undefined.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::run_js_comp_test_adaptive_parameter_search}, whose argument names
#' are the ones an error message reports.
#'
#' @param gene_means a numeric matrix. Mean expression of every gene in every study, NaN for a missing gene
#'   NaN is permitted for this value.
#' @param residuals a numeric array of rank 3. Signed residuals of every replicate of every gene in every study, NaN for a missing value
#'   NaN is permitted for this value.
#' @param n_bootstraps a integer scalar. Number of bootstraps to perform for a candidate
#'   The minimum valid value is `1`.
#' @param join_method a string, one of "join_min", "join_max", "join_median". How the studies' confidence-interval overlaps are joined into the plateau condition
#' @param max_n_points_candidate a integer scalar. Most reference points a candidate may emerge with; a candidate with more is rejected
#'   with the capacity status. The per-point outputs and traces below are sized by it, so
#'   the caller has to know it up front
#'   It is recommended to compute this argument from the `max_n_points_candidate` output produced by \code{\link{calc_adaptive_js_comp_test_bounds}}.
#'   The minimum valid value is `1`.
#' @param min_residuals_per_bin a integer scalar. Minimum count every bin must reach, both in the bin-count search (the floor in every
#'   `bin_count_mode`) and in the second admissibility gate (the consensus pmf)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param min_neighbor_overlap a numeric scalar. Minimum fractional overlap of two consecutive neighborhoods' pooled ranges (gate 1)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.1`.
#' @param succeeding_ci_overlap a numeric scalar. Minimum fractional overlap of a candidate's confidence interval with the best one so
#'   far, per `join_method`. `0` makes the first admissible candidate plateau under the
#'   CI-overlap criterion (for baseline sweeps)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @param plateau_mode a string, one of "plateau_ci_overlap", "plateau_effect_size", "plateau_both". Which plateau criterion decides when the search stops
#'
#'   The default value is `"plateau_ci_overlap"`.
#' @param delta_median_threshold a numeric scalar. Upper bound on the median relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.05`.
#' @param delta_max_threshold a numeric scalar. Upper bound on the largest relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.10`.
#' @param delta_epsilon a numeric scalar. Floor of the denominators of the relative JSD change and relative CI width
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `1.0e-10`.
#' @param delta_min_consecutive_transitions a integer scalar. Number of consecutive qualifying transitions for an effect-size plateau
#'   The minimum valid value is `1`.
#'   The default value is `2`.
#' @param m_min a integer scalar. Smallest bin count the bin-count search tests (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor of the occupancy search; unused in the Stone modes
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's lower residual-range bound
#'   Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
#'   with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's upper residual-range bound
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Two-sided significance level of the bootstrap interval, in percent per end
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param tau a numeric scalar. Largest relative dispersion increase a growth round may cause and still be committed
#'   The minimum valid value is `0.0`.
#'   The default value is `0.1`.
#' @param mad_distance_factor a numeric scalar. Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
#'   The minimum valid value is `0.0`.
#'   The default value is `1.0`.
#' @param max_pooled_residuals a integer scalar. Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
#'   The minimum valid value is `0`.
#'   The default value is `0`.
#' @param min_study_neighbors a integer scalar. Fewest entries of each study every neighborhood must have for a candidate to be
#'   admissible
#'   The minimum valid value is `1`.
#'   The default value is `1`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen, for every candidate; the
#'   Stone modes test every bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @return a named list with elements:
#'   \item{k_start}{a integer scalar. The finally chosen candidate's `k_start`}
#'   \item{k_step}{a integer scalar. The finally chosen candidate's `k_step`}
#'   \item{k_max}{a integer scalar. The finally chosen candidate's `k_max`}
#'   \item{n_points}{a integer scalar. Number of reference points the finally chosen candidate's construction produced}
#'   \item{n_bins_per_point}{a integer vector. The chosen candidate's bin count per reference point; the first `n_points` entries
#'     are meaningful. Without an admissible candidate, the first
#'     `min(n_points, max_n_points_candidate)` entries are `m_min`}
#'   \item{shared_residual_range_low}{a numeric vector. The chosen candidate's lower residual-range bound (R_low) per reference point; the
#'     first `n_points` entries are meaningful, `0.0` without an admissible candidate}
#'   \item{shared_residual_range_high}{a numeric vector. The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`}
#'   \item{best_candidate_confidence_interval}{a numeric matrix. The bootstrapped JSD confidence interval of the finally chosen candidate;
#'     `-1.0` throughout where no admissible candidate backs the result}
#'   \item{plateau_established}{a logical scalar. `TRUE` when a plateau was found (by the criterion `plateau_mode` selects) or the
#'     only candidate was admissible; `FALSE` otherwise}
#'   \item{trace_k_start}{a integer vector. Per-admissible-candidate `k_start`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_k_step}{a integer vector. Per-admissible-candidate `k_step`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_k_max}{a integer vector. Per-admissible-candidate `k_max`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_n_points}{a integer vector. Per-admissible-candidate number of reference points
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_global_js_divergence}{a numeric matrix. Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_lower}{a numeric matrix. Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_upper}{a numeric matrix. Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_width}{a numeric matrix. Per-admissible-candidate, per-study confidence-interval width
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_width_relative}{a numeric matrix. Per-admissible-candidate, per-study relative confidence-interval width
#'     (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_overlap}{a numeric matrix. Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
#'     interval with the best candidate's so far, relative to its own width; `-1.0` for the
#'     first admissible candidate
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta}{a numeric matrix. Per-admissible-candidate, per-study relative JSD change from the previous admissible
#'     candidate; `-1.0` for the first
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta_median}{a numeric vector. Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta_max}{a numeric vector. Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_selected_n_bins}{a integer matrix. Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
#'     `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
#'     `[:trace_n_points[t], t]`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_occupancy_failed}{a logical matrix. Per-admissible-candidate, per-reference-point flag: `TRUE` iff the bin-count search
#'     failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_n_pooled_residuals}{a integer matrix. Per-admissible-candidate, per-reference-point pooled non-NaN residual count (`N_j`);
#'     jagged as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_min_bin_occupancy}{a integer matrix. Per-admissible-candidate, per-reference-point minimum bin occupancy; jagged as
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_mean_bin_occupancy}{a numeric matrix. Per-admissible-candidate, per-reference-point mean bin occupancy; jagged as
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_max_bin_occupancy}{a integer matrix. Per-admissible-candidate, per-reference-point maximum bin occupancy; jagged as
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_sturges_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Sturges' rule diagnostic; jagged as
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_fd_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic; jagged
#'     as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_low}{a numeric matrix. Per-admissible-candidate, per-reference-point lower residual-range bound (R_low);
#'     jagged as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_high}{a numeric matrix. Per-admissible-candidate, per-reference-point upper residual-range bound (R_high);
#'     jagged as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{candidate_k_start}{a integer vector. `k_start` of every candidate tried, in sequence order
#'     The first `n_candidates_tried` elements will hold the results.}
#'   \item{candidate_n_points}{a integer vector. Number of reference points each tried candidate's construction produced (0 for too
#'     few means)
#'     The first `n_candidates_tried` elements will hold the results.}
#'   \item{candidate_status}{a integer vector. What happened to every tried candidate: `0` evaluated,
#'     `1` construction failed,
#'     `2` a study neighborhood below `min_study_neighbors`
#'     or a study without a non-NaN residual at some point,
#'     `3` capacity exceeded,
#'     `4` overlap gate failed,
#'     `5` consensus-pmf count gate failed
#'     The first `n_candidates_tried` elements will hold the results.}
#' @export
run_js_comp_test_adaptive_parameter_search <- function(gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_ci_overlap", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, tau = 0.1, mad_distance_factor = 1.0, max_pooled_residuals = 0L, min_study_neighbors = 1L, bin_count_mode = "bin_count_occupancy") {
    gene_means <- .tox_as_double_matrix(gene_means, "gene_means")
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    join_method <- .tox_as_mode(join_method, "join_method", c("join_min", "join_max", "join_median"))
    max_n_points_candidate <- .tox_as_integer_scalar(max_n_points_candidate, "max_n_points_candidate")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    min_neighbor_overlap <- .tox_as_double_scalar(min_neighbor_overlap, "min_neighbor_overlap")
    succeeding_ci_overlap <- .tox_as_double_scalar(succeeding_ci_overlap, "succeeding_ci_overlap")
    plateau_mode <- .tox_as_mode(plateau_mode, "plateau_mode", c("plateau_ci_overlap", "plateau_effect_size", "plateau_both"))
    delta_median_threshold <- .tox_as_double_scalar(delta_median_threshold, "delta_median_threshold")
    delta_max_threshold <- .tox_as_double_scalar(delta_max_threshold, "delta_max_threshold")
    delta_epsilon <- .tox_as_double_scalar(delta_epsilon, "delta_epsilon")
    delta_min_consecutive_transitions <- .tox_as_integer_scalar(delta_min_consecutive_transitions, "delta_min_consecutive_transitions")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    gamma_occupancy <- .tox_as_double_scalar(gamma_occupancy, "gamma_occupancy")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    two_sided_bootstrapping_significance_level <- .tox_as_double_scalar(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    tau <- .tox_as_double_scalar(tau, "tau")
    mad_distance_factor <- .tox_as_double_scalar(mad_distance_factor, "mad_distance_factor")
    max_pooled_residuals <- .tox_as_integer_scalar(max_pooled_residuals, "max_pooled_residuals")
    min_study_neighbors <- .tox_as_integer_scalar(min_study_neighbors, "min_study_neighbors")
    bin_count_mode <- .tox_as_mode(bin_count_mode, "bin_count_mode", c("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled"))
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_adaptive_parameter_search_call", gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, tau, mad_distance_factor, max_pooled_residuals, min_study_neighbors, bin_count_mode)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "n_bootstraps", "join_method", "max_n_points_candidate", "k_start", "k_step", "k_max", "n_points", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "best_candidate_confidence_interval", "plateau_established", "n_admissible_evaluated", "trace_k_start", "trace_k_step", "trace_k_max", "trace_n_points", "trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_ci_overlap", "trace_delta", "trace_delta_median", "trace_delta_max", "trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy", "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins", "trace_shared_residual_range_low", "trace_shared_residual_range_high", "n_candidates_tried", "candidate_k_start", "candidate_n_points", "candidate_status", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "bin_count_mode", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", NA_character_, NA_character_, NA_character_, NA_character_, "n_bins_per_point", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        k_start = .result$k_start,
        k_step = .result$k_step,
        k_max = .result$k_max,
        n_points = .result$n_points,
        n_bins_per_point = .result$n_bins_per_point,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        best_candidate_confidence_interval = .result$best_candidate_confidence_interval,
        plateau_established = .result$plateau_established,
        trace_k_start = utils::head(.result$trace_k_start, .result$n_admissible_evaluated),
        trace_k_step = utils::head(.result$trace_k_step, .result$n_admissible_evaluated),
        trace_k_max = utils::head(.result$trace_k_max, .result$n_admissible_evaluated),
        trace_n_points = utils::head(.result$trace_n_points, .result$n_admissible_evaluated),
        trace_global_js_divergence = .result$trace_global_js_divergence[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_lower = .result$trace_ci_lower[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_upper = .result$trace_ci_upper[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_width = .result$trace_ci_width[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_width_relative = .result$trace_ci_width_relative[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_overlap = .result$trace_ci_overlap[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_delta = .result$trace_delta[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_delta_median = utils::head(.result$trace_delta_median, .result$n_admissible_evaluated),
        trace_delta_max = utils::head(.result$trace_delta_max, .result$n_admissible_evaluated),
        trace_selected_n_bins = .result$trace_selected_n_bins[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_occupancy_failed = .result$trace_occupancy_failed[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_n_pooled_residuals = .result$trace_n_pooled_residuals[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_min_bin_occupancy = .result$trace_min_bin_occupancy[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_mean_bin_occupancy = .result$trace_mean_bin_occupancy[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_max_bin_occupancy = .result$trace_max_bin_occupancy[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_sturges_bins = .result$trace_sturges_bins[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_fd_bins = .result$trace_fd_bins[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_shared_residual_range_low = .result$trace_shared_residual_range_low[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_shared_residual_range_high = .result$trace_shared_residual_range_high[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        candidate_k_start = utils::head(.result$candidate_k_start, .result$n_candidates_tried),
        candidate_n_points = utils::head(.result$candidate_n_points, .result$n_candidates_tried),
        candidate_status = utils::head(.result$candidate_status, .result$n_candidates_tried)
    )
}
