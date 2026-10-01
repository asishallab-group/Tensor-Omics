# Generated. Do not edit.

#' Estimate the histogram bin count for one (n_points, n_neighbors) candidate
#'
#' Ported from 125-stabilize-jscomp's `estimate_bin_count_helper`. Computes Sturges' rule and
#' the Freedman-Diaconis rule (without doubling the bin width, since `shared_residual_range`
#' is already the one-sided half of the full `[-R, R]` histogram range, so dividing the full
#' range by the undoubled Freedman-Diaconis width already gives the doubled rule's bin count)
#' independently, each clamped on its own to at most
#' \code{MAX_N_BINS} bins and returned as
#' `sturges_bins`/`fd_bins` for diagnostics, with their maximum returned as `n_bins`. When the
#' interquartile range is too close to zero to safely divide by, the Freedman-Diaconis term is
#' skipped and `fd_bins` falls back to the (clamped) Sturges estimate instead of blowing up.
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
#'     to the (clamped) sturges_bins when the interquartile range is too close to zero to
#'     divide by (see the is_close guard below)}
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
#' Ported from 125-stabilize-jscomp's `estimate_bin_count_helper`. Computes Sturges' rule and
#' the Freedman-Diaconis rule (without doubling the bin width, since `shared_residual_range`
#' is already the one-sided half of the full `[-R, R]` histogram range, so dividing the full
#' range by the undoubled Freedman-Diaconis width already gives the doubled rule's bin count)
#' independently, each clamped on its own to at most
#' \code{MAX_N_BINS} bins and returned as
#' `sturges_bins`/`fd_bins` for diagnostics, with their maximum returned as `n_bins`. When the
#' interquartile range is too close to zero to safely divide by, the Freedman-Diaconis term is
#' skipped and `fd_bins` falls back to the (clamped) Sturges estimate instead of blowing up.
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
#'     to the (clamped) sturges_bins when the interquartile range is too close to zero to
#'     divide by (see the is_close guard below)}
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
#' Implements Issue #187's two-stage geometric-search-then-local-refinement algorithm for one
#' neighborhood's pooled residuals (`pooled_residuals`, across all its neighbors and all
#' studies): find the largest bin count `M` in `[m_min, m_max]` whose equal-width histogram
#' over `[shared_residual_range_low, shared_residual_range_high]` -- this neighborhood's own
#' asymmetric range, the `lower_residual_range_quantile`/`upper_residual_range_quantile`
#' percentiles of its own pooled signed residuals, rather than a single dataset-wide symmetric
#' range -- has every bin at or above `min_residuals_per_bin` (the occupancy criterion), rather
#' than the generic
#' Sturges/Freedman-Diaconis rule
#' \code{\link{estimate_bin_count}} alone
#' applies, which is why that routine is still called here too -- purely for the
#' `sturges_bins`/`fd_bins` diagnostic outputs, never for the decision itself.
#'
#' Stage 1 grows a candidate `trial_m` geometrically from `m_min` (`next_m =
#' ceiling(gamma_occupancy*trial_m)`, guaranteed to advance by at least 1 via
#' `max(trial_m + 1, next_m)`, clamped to `m_max`) until occupancy first fails, recording the
#' largest admissible `m_valid` and the first inadmissible `m_invalid`; reaching `m_max` while
#' still admissible returns it immediately, skipping stage 2 entirely. Stage 2 then tests every
#' integer strictly between `m_valid` and `m_invalid` -- not a binary search, since equal-width
#' bin boundaries are recomputed for every candidate `M` and occupancy is therefore not
#' guaranteed monotonic in `M` (Issue #187 is explicit about this) -- and keeps the largest one
#' that still passes. Both stages reuse the occupancy diagnostics (`min`/`max_bin_occupancy`)
#' computed for whichever `M` ends up selected, rather than recomputing them a second time
#' afterward; `mean_bin_occupancy` needs no such bookkeeping, since every non-NaN pooled
#' residual lands in exactly one bin at any `M`, so it is always exactly
#' `n_pooled_residuals / selected_n_bins`.
#'
#' Per the issue's own FAILURE policy, `occupancy_failed = TRUE` (even `m_min` bins could not
#' satisfy the occupancy criterion, or every pooled residual is NaN) means the neighborhood
#' should be rejected by the caller rather than built from `selected_n_bins` -- which is still
#' set to `m_min` in that case, purely so a caller ignoring `occupancy_failed` has *a* value to
#' build with, never as an indication the search actually found `m_min` admissible.
#'
#' The outer geometric search is a genuine sequential `do`/`exit` state machine, not `do
#' concurrent`: `m_valid`/`m_invalid` accumulate across iterations and each iteration's
#' continuation depends on the previous one's outcome, exactly the "loops with data-dependent
#' control flow across iterations" case Fortran_Coding_Guides.pdf Sec 10 carves out as the
#' deliberate exception to `do concurrent`.
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
#' @param m_max a integer scalar. Largest candidate bin count the search will ever test (M_max); if a caller passes
#'   `m_max < m_min`, the implementation clamps it up to `m_min` internally rather than
#'   relying on an unconfirmed generator capability to bound one optional argument by
#'   another
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
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] for this neighborhood's own upper residual-range bound
#'   (shared_residual_range_high)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @return a named list with elements:
#'   \item{selected_n_bins}{a integer scalar. The chosen M_j: the largest bin count in [m_min, m_max] whose pooled histogram has
#'     every bin at or above min_residuals_per_bin; m_min when occupancy_failed}
#'   \item{occupancy_failed}{a logical scalar. `TRUE` iff even m_min bins could not satisfy the occupancy criterion (including
#'     the case where every pooled residual is NaN) -- per Issue #187's FAILURE policy, the
#'     caller should reject this neighborhood rather than build a histogram from
#'     selected_n_bins}
#'   \item{shared_residual_range_low}{a numeric scalar. This neighborhood's own lower residual-range bound (R_low): the
#'     lower_residual_range_quantile percentile of its own pooled signed residuals --
#'     replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
#'     occupancy_failed because every pooled residual is NaN}
#'   \item{shared_residual_range_high}{a numeric scalar. This neighborhood's own upper residual-range bound (R_high): the
#'     upper_residual_range_quantile percentile of its own pooled signed residuals --
#'     replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
#'     occupancy_failed because every pooled residual is NaN}
#'   \item{n_pooled_residuals}{a integer scalar. Count of non-NaN pooled residuals (N_j)}
#'   \item{min_bin_occupancy}{a integer scalar. Minimum bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{mean_bin_occupancy}{a numeric scalar. Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
#'     occupancy_failed}
#'   \item{max_bin_occupancy}{a integer scalar. Maximum bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{sturges_bins}{a integer scalar. Sturges' rule estimate for this neighborhood's own pooled residuals, from
#'     estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision}
#'   \item{fd_bins}{a integer scalar. Freedman-Diaconis rule estimate for this neighborhood's own pooled residuals, from
#'     estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision}
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
#' Implements Issue #187's two-stage geometric-search-then-local-refinement algorithm for one
#' neighborhood's pooled residuals (`pooled_residuals`, across all its neighbors and all
#' studies): find the largest bin count `M` in `[m_min, m_max]` whose equal-width histogram
#' over `[shared_residual_range_low, shared_residual_range_high]` -- this neighborhood's own
#' asymmetric range, the `lower_residual_range_quantile`/`upper_residual_range_quantile`
#' percentiles of its own pooled signed residuals, rather than a single dataset-wide symmetric
#' range -- has every bin at or above `min_residuals_per_bin` (the occupancy criterion), rather
#' than the generic
#' Sturges/Freedman-Diaconis rule
#' \code{\link{estimate_bin_count}} alone
#' applies, which is why that routine is still called here too -- purely for the
#' `sturges_bins`/`fd_bins` diagnostic outputs, never for the decision itself.
#'
#' Stage 1 grows a candidate `trial_m` geometrically from `m_min` (`next_m =
#' ceiling(gamma_occupancy*trial_m)`, guaranteed to advance by at least 1 via
#' `max(trial_m + 1, next_m)`, clamped to `m_max`) until occupancy first fails, recording the
#' largest admissible `m_valid` and the first inadmissible `m_invalid`; reaching `m_max` while
#' still admissible returns it immediately, skipping stage 2 entirely. Stage 2 then tests every
#' integer strictly between `m_valid` and `m_invalid` -- not a binary search, since equal-width
#' bin boundaries are recomputed for every candidate `M` and occupancy is therefore not
#' guaranteed monotonic in `M` (Issue #187 is explicit about this) -- and keeps the largest one
#' that still passes. Both stages reuse the occupancy diagnostics (`min`/`max_bin_occupancy`)
#' computed for whichever `M` ends up selected, rather than recomputing them a second time
#' afterward; `mean_bin_occupancy` needs no such bookkeeping, since every non-NaN pooled
#' residual lands in exactly one bin at any `M`, so it is always exactly
#' `n_pooled_residuals / selected_n_bins`.
#'
#' Per the issue's own FAILURE policy, `occupancy_failed = TRUE` (even `m_min` bins could not
#' satisfy the occupancy criterion, or every pooled residual is NaN) means the neighborhood
#' should be rejected by the caller rather than built from `selected_n_bins` -- which is still
#' set to `m_min` in that case, purely so a caller ignoring `occupancy_failed` has *a* value to
#' build with, never as an indication the search actually found `m_min` admissible.
#'
#' The outer geometric search is a genuine sequential `do`/`exit` state machine, not `do
#' concurrent`: `m_valid`/`m_invalid` accumulate across iterations and each iteration's
#' continuation depends on the previous one's outcome, exactly the "loops with data-dependent
#' control flow across iterations" case Fortran_Coding_Guides.pdf Sec 10 carves out as the
#' deliberate exception to `do concurrent`.
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
#' @param m_max a integer scalar. Largest candidate bin count the search will ever test (M_max); if a caller passes
#'   `m_max < m_min`, the implementation clamps it up to `m_min` internally rather than
#'   relying on an unconfirmed generator capability to bound one optional argument by
#'   another
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
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] for this neighborhood's own upper residual-range bound
#'   (shared_residual_range_high)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @return a named list with elements:
#'   \item{selected_n_bins}{a integer scalar. The chosen M_j: the largest bin count in [m_min, m_max] whose pooled histogram has
#'     every bin at or above min_residuals_per_bin; m_min when occupancy_failed}
#'   \item{occupancy_failed}{a logical scalar. `TRUE` iff even m_min bins could not satisfy the occupancy criterion (including
#'     the case where every pooled residual is NaN) -- per Issue #187's FAILURE policy, the
#'     caller should reject this neighborhood rather than build a histogram from
#'     selected_n_bins}
#'   \item{shared_residual_range_low}{a numeric scalar. This neighborhood's own lower residual-range bound (R_low): the
#'     lower_residual_range_quantile percentile of its own pooled signed residuals --
#'     replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
#'     occupancy_failed because every pooled residual is NaN}
#'   \item{shared_residual_range_high}{a numeric scalar. This neighborhood's own upper residual-range bound (R_high): the
#'     upper_residual_range_quantile percentile of its own pooled signed residuals --
#'     replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
#'     occupancy_failed because every pooled residual is NaN}
#'   \item{n_pooled_residuals}{a integer scalar. Count of non-NaN pooled residuals (N_j)}
#'   \item{min_bin_occupancy}{a integer scalar. Minimum bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{mean_bin_occupancy}{a numeric scalar. Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
#'     occupancy_failed}
#'   \item{max_bin_occupancy}{a integer scalar. Maximum bin count at selected_n_bins; 0 when occupancy_failed}
#'   \item{sturges_bins}{a integer scalar. Sturges' rule estimate for this neighborhood's own pooled residuals, from
#'     estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision}
#'   \item{fd_bins}{a integer scalar. Freedman-Diaconis rule estimate for this neighborhood's own pooled residuals, from
#'     estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision}
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
#' Tests every candidate bin count `M` in `[m_min, m_max]` independently and keeps the largest
#' one whose pooled histogram satisfies the occupancy criterion, instead of
#' \code{\link{determine_bin_count_occupancy}}'s
#' own fast geometric-search-then-refinement. Exists purely to validate that routine's result:
#' occupancy is not guaranteed monotonic in `M` once bin boundaries are recomputed per candidate
#' (Issue #187 is explicit about this), so a search that stops at the first failure can in
#' principle miss a larger, independently-admissible `M` the geometric ladder never tries. Takes
#' `shared_residual_range_low`/`shared_residual_range_high`/`n_pooled_residuals` as direct
#' inputs, already produced by a prior
#' \code{\link{determine_bin_count_occupancy}}
#' call, rather than re-deriving them -- this isolates the comparison to just the `M`-selection
#' algorithm, uncontaminated by a second independent percentile computation.
#'
#' Every candidate is independent (no early exit, no state carried between iterations, unlike
#' the production routine's own Stage 1/Stage 2), so this is a genuine `do concurrent` with
#' `reduce(max:...)`, not a sequential search: `candidate_bin_counts` is declared local to the
#' loop (an ordinary MAX_N_BINS-sized local, not a `tmp_` dummy -- precedented by
#' `run_js_comp_test_impl`'s own `candidates_n_points_n_neighbors` local) so each concurrent
#' iteration gets its own private scratch instead of racing on a shared buffer sliced by
#' `trial_m`. `reduce(max:...)` has no "argmax" form, so the winning `M`'s own
#' min/mean/max_bin_occupancy are recovered with one extra, ordinary (non-concurrent) call to
#' `histogram_bin_counts` for `best_m` alone once the reduction is done -- trivial cost next to
#' the search itself.
#'
#' `n_pooled_residuals == 0` and a degenerate zero-width range need no special-case branch here:
#' `histogram_bin_counts` already guards the degenerate range internally, and an all-zero pool
#' naturally resolves to `occupancy_failed` (or a trivial pass at `min_residuals_per_bin=0`,
#' with `mean_bin_occupancy=0.0`, no divide-by-zero) -- intentional, not an oversight.
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
#' @param m_max a integer scalar. Largest candidate bin count tested (M_max); if a caller passes `m_max < m_min`, the
#'   implementation clamps it up to `m_min` internally rather than relying on an
#'   unconfirmed generator capability to bound one optional argument by another
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
#'   \item{occupancy_failed}{a logical scalar. `TRUE` iff no candidate bin count in [m_min, m_max] satisfies the occupancy
#'     criterion (including the case where n_pooled_residuals is 0 and min_residuals_per_bin
#'     is not itself 0)}
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
#' Tests every candidate bin count `M` in `[m_min, m_max]` independently and keeps the largest
#' one whose pooled histogram satisfies the occupancy criterion, instead of
#' \code{\link{determine_bin_count_occupancy}}'s
#' own fast geometric-search-then-refinement. Exists purely to validate that routine's result:
#' occupancy is not guaranteed monotonic in `M` once bin boundaries are recomputed per candidate
#' (Issue #187 is explicit about this), so a search that stops at the first failure can in
#' principle miss a larger, independently-admissible `M` the geometric ladder never tries. Takes
#' `shared_residual_range_low`/`shared_residual_range_high`/`n_pooled_residuals` as direct
#' inputs, already produced by a prior
#' \code{\link{determine_bin_count_occupancy}}
#' call, rather than re-deriving them -- this isolates the comparison to just the `M`-selection
#' algorithm, uncontaminated by a second independent percentile computation.
#'
#' Every candidate is independent (no early exit, no state carried between iterations, unlike
#' the production routine's own Stage 1/Stage 2), so this is a genuine `do concurrent` with
#' `reduce(max:...)`, not a sequential search: `candidate_bin_counts` is declared local to the
#' loop (an ordinary MAX_N_BINS-sized local, not a `tmp_` dummy -- precedented by
#' `run_js_comp_test_impl`'s own `candidates_n_points_n_neighbors` local) so each concurrent
#' iteration gets its own private scratch instead of racing on a shared buffer sliced by
#' `trial_m`. `reduce(max:...)` has no "argmax" form, so the winning `M`'s own
#' min/mean/max_bin_occupancy are recovered with one extra, ordinary (non-concurrent) call to
#' `histogram_bin_counts` for `best_m` alone once the reduction is done -- trivial cost next to
#' the search itself.
#'
#' `n_pooled_residuals == 0` and a degenerate zero-width range need no special-case branch here:
#' `histogram_bin_counts` already guards the degenerate range internally, and an all-zero pool
#' naturally resolves to `occupancy_failed` (or a trivial pass at `min_residuals_per_bin=0`,
#' with `mean_bin_occupancy=0.0`, no divide-by-zero) -- intentional, not an oversight.
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
#' @param m_max a integer scalar. Largest candidate bin count tested (M_max); if a caller passes `m_max < m_min`, the
#'   implementation clamps it up to `m_min` internally rather than relying on an
#'   unconfirmed generator capability to bound one optional argument by another
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
#'   \item{occupancy_failed}{a logical scalar. `TRUE` iff no candidate bin count in [m_min, m_max] satisfies the occupancy
#'     criterion (including the case where n_pooled_residuals is 0 and min_residuals_per_bin
#'     is not itself 0)}
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

#' Generate the GAMMA-decay (n_points, n_neighbors) candidate grid
#'
#' Ported from the grid-building half of 125-stabilize-jscomp's
#' `determine_js_comp_test_n_points_n_neighbors_helper`: starting from an initial
#' `n_points_high` (clamped between MIN_POINTS and MAX_POINTS), repeatedly multiplies by
#' GAMMA until it would drop below `n_points_low`, and for each distinct resulting
#' `n_points` candidate pairs it with up to `size(KX_FACTORS)` distinct `n_neighbors`
#' candidates. A duplicate `n_points` or `n_neighbors` value (from clamping or
#' integer rounding) collapses rather than repeating -- this is real, derived behavior the
#' grid depends on to avoid redundant candidates at small `max_n_genes_all_studies`, not a
#' bug: a small enough `max_n_genes_all_studies` collapses the whole grid down to exactly one
#' candidate.
#'
#' Issue #187: this routine no longer estimates a per-candidate histogram bin count as a side
#' effect -- both
#' \code{\link{run_js_comp_test}} and
#' \code{\link{run_js_comp_test_parameter_search}}
#' now compute real per-neighborhood bin counts via
#' \code{\link{determine_bin_count_occupancy}}
#' once a candidate has passed admissibility, superseding the old global-pool
#' `estimate_bin_count_impl` estimate this routine used to produce for every candidate
#' regardless of admissibility.
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
#' The adaptive-neighborhood counterpart (Issue #217) of
#' \code{\link{generate_js_comp_test_candidates}}.
#' An adaptive candidate does not fix its number of reference points: each neighborhood is
#' grown from `k_start` pooled entries in rounds of `k_step` up to `k_max` while its residual
#' dispersion stays stable, and the next reference point is placed beyond it, so the point
#' count emerges from the growth knobs. This routine only generates those knobs. Because
#' `k_start` shrinks from one candidate to the next, the candidates **ascend** in the number of
#' reference points they are expected to produce.
#'
#' With `N = max_n_genes_all_studies * n_studies` (the padded pool size, not the count of
#' non-NaN means, so that the capacity producer
#' \code{\link{calc_adaptive_js_comp_test_bounds}}
#' and the search agree on the same number) and `floor_k = max(10, 2*n_studies)`:
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
#' With them, one decay step lowers any `k_start >= 10` by at least 2, so the values are
#' strictly descending without any deduplication, every `k_step` is at least 1 and every
#' `k_max` at least `k_start`.
#'
#' Each `k_start_t` is computed in double precision directly from `k_start_1`, not by
#' rounding the previous candidate again, so rounding never accumulates along the sequence.
#'
#' Every knob scales with `N`, so the sequence does not collapse for small data the way the
#' fixed-k grid does below about 8,743 genes, where its absolute `n_points_low` floor leaves
#' room for only a single `n_points` value. The first candidate is always accepted, so even a
#' data set smaller than `k_start_1` (for example 3 genes in 1 study, where `k_start_1 = 10`)
#' yields exactly one candidate; growing its neighborhoods then fails with a too-few-means
#' status, which is the adaptive search's to report, not this routine's.
#'
#' `N` is formed in 64-bit integer arithmetic and never overflows. The generated values are
#' representable as long as `k_start_1` does not exceed `huge(1)/4 = 536870911` (so
#' that `k_max_1` fits a 32-bit integer), i.e. for `N` up to about 2.7e10 and `n_studies` up
#' to about 2.7e8 -- far beyond any real data set. This routine does not check that bound;
#' \code{\link{calc_adaptive_js_comp_test_bounds}},
#' which the adaptive search needs for its sizing anyway, rejects inputs beyond it.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::generate_adaptive_js_comp_test_candidates}, whose argument names
#' are the ones an error message reports.
#'
#' @param max_n_genes_all_studies a integer scalar. Maximum number of genes across all studies
#'   The minimum valid value is `1`.
#' @param n_studies a integer scalar. Number of studies
#'   The minimum valid value is `1`.
#' @return a integer matrix. Candidate `[k_start, k_step, k_max]` triples, `k_start` strictly descending (so the
#'   emerging number of reference points ascends)
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
#' The heteroscedastic neighborhood construction of Issue #217, after Aaron Schroeder's
#' noise-model growth: instead of one fixed neighbor count, each neighborhood grows while the
#' mean absolute residual of its pool stays stable, and the next reference point is placed
#' beyond it. The number of reference points emerges from the growth knobs `k_start`, `k_step`
#' and `k_max`, all counted in pooled entries (one entry is one gene of one study).
#'
#' **Pooling.** All studies' gene means are pooled into `N = max_n_genes_all_studies *
#' n_studies` flat entries (entry `(s - 1)*G + g` is gene `g` of study `s`) and sorted
#' ascending; NaN means sort last and take no part. A neighborhood is a contiguous run
#' `a..b` of that sorted order, returned as `pooled_neighborhood_range`, so its per-study
#' membership follows by decoding every position back to its study, and the studies' counts
#' `n_neighbors_per_point` add up to `b - a + 1` exactly. These ranges are what
#' \code{\link{check_neighborhood_overlaps}}
#' takes; unlike the fixed-k ranges they are not extended over tied means, because a range is
#' the exact membership here.
#'
#' **Growth of one neighborhood** around its reference point `x_star`, the mean at its seed
#' position. Entries join nearest first by `|mean - x_star|`, a tie going to the lower mean.
#' The dispersion `S` is the mean absolute value of every non-NaN residual of every entry in
#' the pool.
#'
#' 1. The first `k_start` entries join unconditionally.
#' 2. If they hold fewer than
#' \code{ADAPTIVE_MIN_VALID_RESIDUALS}
#' (`10`) non-NaN residuals, the neighborhood stops there with stop reason
#' \code{ADAPTIVE_STOP_TOO_FEW_RESIDUALS}
#' (`6`). If their residuals are all zero, it stops with
#' \code{ADAPTIVE_STOP_ZERO_DISPERSION}
#' (`5`), since no relative change is defined.
#' 3. Otherwise it grows in rounds: the next `min(k_step, k_max - count)` nearest entries are
#' staged and `S_new` of the pool with them is compared with the previous round's `S_old`.
#' If `(S_new - S_old) / S_old > tau` the round is discarded and growth stops
#' (\code{ADAPTIVE_STOP_TAU}, `1`);
#' otherwise the round is committed and `S_old = S_new`. A round without any non-NaN
#' residual leaves `S` unchanged and is committed.
#' 4. Growth also stops once the neighborhood holds `k_max` entries
#' (\code{ADAPTIVE_STOP_K_MAX}, `2`),
#' or once no entry is left on either side
#' (\code{ADAPTIVE_STOP_EXHAUSTED}, `4`);
#' `k_max` is checked first. A round cut short by the pool's end is still evaluated.
#' 5. With `max_pooled_residuals > 0`, an adaptive round stages whole entries only while the
#' pool's non-NaN residuals stay at or below that cap; NaN residuals do not count toward
#' it, and the `k_start` phase is not capped. A round cut short by the cap is evaluated as
#' usual and, if committed, ends growth; a round that could stage nothing ends it at once.
#' Both stop with
#' \code{ADAPTIVE_STOP_RESIDUAL_CAP} (`3`).
#'
#' **Reference points.** The first seed is the smallest pooled mean. After neighborhood `i`
#' ends at sorted position `b_i`, the next target is `mean(b_i) + mad_distance_factor *
#' MAD_i`, `MAD_i` being the raw median absolute deviation of neighborhood `i`'s pooled means
#' (see \code{calc_sorted_slice_mad}),
#' and the next seed is the position after `b_i` whose mean is nearest that target (a tie
#' going to the lower mean; the largest pooled mean if the target lies beyond it); where
#' several positions share that mean, the seed is the first of them. Seeds therefore
#' strictly advance, which bounds the number of reference points by `N`; the
#' construction ends with the neighborhood that reaches the largest pooled mean.
#'
#' **Construction status.**
#' \code{ADAPTIVE_STATUS_TOO_FEW_MEANS}
#' (`1`) when fewer than `k_start` pooled means are non-NaN (all-NaN means included): no
#' neighborhood is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction
#' always completes, so that it can be inspected, and the status reports
#' \code{ADAPTIVE_STATUS_TOO_FEW_RESIDUALS}
#' (`3`) if any neighborhood stopped with too few residuals, else
#' \code{ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD}
#' (`2`) if any neighborhood has fewer than `min_study_neighbors` entries of some study, else
#' \code{ADAPTIVE_STATUS_OK} (`0`). This
#' priority is fixed, whichever neighborhood failed first.
#'
#' All outputs sized by `N` are filled only for their first `n_points` reference points.
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
#' @param tau a numeric scalar. Largest relative increase of the dispersion an adaptive round may cause and still
#'   be committed
#'   The minimum valid value is `0.0`.
#'   The default value is `0.1`.
#' @param mad_distance_factor a numeric scalar. Multiple of a neighborhood's median absolute deviation that the next reference
#'   point's target lies beyond the neighborhood's largest mean
#'   The minimum valid value is `0.0`.
#'   The default value is `1.0`.
#' @param max_pooled_residuals a integer scalar. Cap on the non-NaN residuals adaptive rounds may grow a neighborhood's pool to; 0
#'   for no cap
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
#'   \item{neighborhood_dispersion}{a numeric vector. Mean absolute non-NaN residual of each final neighborhood (the last committed
#'     dispersion); NaN when the neighborhood holds no non-NaN residual at all
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
#' Ported from 125-stabilize-jscomp's `test_neighborhood_overlaps_helper`: the first
#' admissibility gate a candidate `(n_points, n_neighbors)` pair must pass, using the
#' `[min_idx, max_idx]` neighborhood spans
#' \code{\link{construct_neighborhoods_ranged}}
#' produces. Named `check_*` rather than 125's `test_*`, a deliberate deviation from the
#' verbatim port: the generated R binding is published under the Fortran name, and the R test
#' harness (`r/test_helpers.R`'s `run_all_tests`) discovers every `test_`-prefixed name in the
#' environment as a test case to run with no arguments -- a `test_`-prefixed export would be
#' swept up and fail every R suite that sources the package, not just this module's own.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::check_neighborhood_overlaps}, whose argument names
#' are the ones an error message reports.
#'
#' @param neighborhood_range a integer matrix. For each reference point, the `[min_idx, max_idx]` neighborhood span, as produced
#'   by construct_neighborhoods_ranged_impl
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
#' Ported from 125-stabilize-jscomp's `test_mean_pmf_min_counts_helper`: the second
#' admissibility gate a candidate `(n_points, n_neighbors)` pair must pass, checked once the
#' first gate
#' (\code{\link{check_neighborhood_overlaps}})
#' has already passed. Named `check_*` rather than 125's `test_*` for the same reason as its
#' sibling above: a `test_`-prefixed R export collides with the R test harness's own
#' test-discovery convention.
#' Issue #187: `n_bins_per_point` scopes the reduction to each point's own valid bin range, so
#' a point whose own bin count is narrower than `n_bins` (this array's second extent, i.e. the
#' widest bin count any point uses) has its legitimate zero-padded columns skipped rather than
#' mistaken for an occupancy failure.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::check_mean_pmf_min_counts}, whose argument names
#' are the ones an error message reports.
#'
#' @param mean_pmf_counts a integer matrix. Absolute counts of a residual per bin for the mean pmf
#'   The minimum valid value is `0`.
#' @param n_bins_per_point a integer vector. This point's own bin count -- only `mean_pmf_counts(1:n_bins_per_point(i_point), i_point)`
#'   is inspected; columns beyond it are legitimate zero-padding, not failures
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
    if (length(n_bins_per_point) != dim(mean_pmf_counts)[2])
        .tox_shape_error("n_bins_per_point", length(n_bins_per_point), "mean_pmf_counts", dim(mean_pmf_counts)[2])

    .result <- .Call("check_mean_pmf_min_counts_call", mean_pmf_counts, n_bins_per_point, min_count)
    .arguments <- c("mean_pmf_counts", "n_bins", "n_bins_per_point", "n_points", "min_count", "all_bins_have_min_count", "ierr")
    .sources <- c(NA_character_, "mean_pmf_counts", NA_character_, "mean_pmf_counts", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$all_bins_have_min_count
}

#' Test one candidate pair's bootstrapped confidence intervals against the current best, and detect a JSD plateau
#'
#' Ported from 125-stabilize-jscomp's `check_plateau_condition_helper`. A search over
#' candidates (finest resolution to coarsest) stops -- "plateaus" -- either when a new
#' candidate is no better than the previous best (the short-circuit below: keep the previous
#' best and stop searching), or once the new candidate's confidence-interval overlap with the
#' previous best meets the condition `join_method` names. `join_method` replaces
#' 125-stabilize-jscomp's hand-rolled join-method-range validation macro entirely: the mode
#' table below is itself the validation, checked against exactly the values it names.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::check_plateau_condition}, whose argument names
#' are the ones an error message reports.
#'
#' @param confidence_interval a numeric matrix. JSD confidence interval `[lower, upper]` from bootstrapping, for the candidate
#'   pair under test
#' @param best_candidate_pair_confidence_interval a numeric matrix. JSD confidence intervals for the current best candidate pair; overwritten with
#'   `confidence_interval` unless the new candidate is worse
#' @param best_candidate_index a integer scalar. Candidate-grid index of the current best candidate pair; overwritten with
#'   `candidate_index` unless the new candidate is worse
#' @param best_exceeded_ci_overlap_count a integer scalar. Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
#'   best candidate pair; overwritten unless the new candidate is worse
#'   The minimum valid value is `0`.
#' @param candidate_index a integer scalar. Candidate-grid index of the candidate pair that produced `confidence_interval`
#'   The minimum valid value is `1`.
#' @param join_method a string, one of "join_min", "join_max", "join_median". The way to evaluate all studies' confidence-interval overlaps for the plateau
#'   condition: METHOD_JOIN_MIN requires every study's overlap to exceed
#'   `succeeding_ci_overlap`, METHOD_JOIN_MAX requires only one study's overlap to
#'   exceed it, and METHOD_JOIN_MEDIAN requires a majority
#'   (`count > (n_studies - 1) / 2`) to exceed it
#' @param succeeding_ci_overlap a numeric scalar. Minimum fractional overlap an interval in `confidence_interval` must have with its
#'   respective interval in `best_candidate_pair_confidence_interval` to count as
#'   "exceeded"
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
    .arguments <- c("confidence_interval", "best_candidate_pair_confidence_interval", "n_studies", "best_candidate_index", "best_exceeded_ci_overlap_count", "candidate_index", "join_method", "succeeding_ci_overlap", "plateau_found", "ierr")
    .sources <- c(NA_character_, NA_character_, "confidence_interval", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        best_candidate_pair_confidence_interval = .result$best_candidate_pair_confidence_interval,
        best_candidate_index = .result$best_candidate_index,
        best_exceeded_ci_overlap_count = .result$best_exceeded_ci_overlap_count,
        plateau_found = .result$plateau_found
    )
}

#' Test one candidate's per-study JSD against the previous admissible candidate for a relative-effect-size plateau
#'
#' Implements Issue #178's relative-effect-size plateau criterion, complementary to
#' \code{\link{check_plateau_condition}}'s
#' CI-overlap one: for each study `i`, the relative change in observed JSD between successive
#' ADMISSIBLE parameter settings (both admissibility gates already passed),
#' `delta(i) = |global_js_divergence(i) - prev_global_js_divergence(i)| / max(prev_global_js_divergence(i),
#' delta_epsilon)`, summarized across studies by its median (`delta_median`, via the
#' already-shipped
#' \code{\link{calc_percentile}}) and maximum (`delta_max`). A
#' plateau is declared once both stay under their respective thresholds for
#' `delta_min_consecutive_transitions` consecutive transitions in a row -- tracked across calls
#' via `n_consecutive_ok`, reset the moment either threshold is missed.
#'
#' No transition exists for the very first admissible candidate a caller ever passes in
#' (`has_previous = FALSE`): `delta`/`delta_median`/`delta_max` are all set to
#' `-1.0` -- the same not-yet-computed sentinel
#' \code{\link{run_js_comp_test_parameter_search}}
#' already uses for its own confidence-interval fallback, usable here for the same reason:
#' every quantity this routine tracks is structurally non-negative.
#'
#' The 0.05/0.10 defaults `run_js_comp_test_parameter_search_impl` passes for
#' `delta_median_threshold`/`delta_max_threshold` are Issue #178's own suggested starting
#' point, explicitly not yet empirically validated -- see that routine's doc comment.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::check_effect_size_plateau_condition}, whose argument names
#' are the ones an error message reports.
#'
#' @param global_js_divergence a numeric vector. Current admissible candidate's observed global JSD per study
#'   The minimum valid value is `0.0`.
#' @param prev_global_js_divergence a numeric vector. Previous admissible candidate's observed global JSD per study; ignored when
#'   `has_previous` is `FALSE`
#'   The minimum valid value is `0.0`.
#' @param has_previous a logical scalar. `FALSE` for the very first admissible candidate a caller has ever passed in, where
#'   no transition exists to compute a relative change from
#' @param delta_median_threshold a numeric scalar. Upper bound the median relative change across studies must stay under for a
#'   transition to count toward a plateau
#'   The minimum valid value is `above(0.0)`.
#' @param delta_max_threshold a numeric scalar. Upper bound the largest relative change across studies must stay under for a
#'   transition to count toward a plateau
#'   The minimum valid value is `above(0.0)`.
#' @param delta_epsilon a numeric scalar. Small constant preventing division by zero when a study's previous JSD was zero
#'   The minimum valid value is `above(0.0)`.
#' @param delta_min_consecutive_transitions a integer scalar. Number of consecutive qualifying transitions required to declare a plateau
#'   The minimum valid value is `1`.
#' @param n_consecutive_ok a integer scalar. Running count of consecutive qualifying transitions; incremented when this
#'   transition qualifies, reset to zero otherwise (and whenever `has_previous` is
#'   `FALSE`)
#'   The minimum valid value is `0`.
#' @return a named list with elements:
#'   \item{n_consecutive_ok}{a integer scalar. Running count of consecutive qualifying transitions; incremented when this
#'     transition qualifies, reset to zero otherwise (and whenever `has_previous` is
#'     `FALSE`)
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
#' Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_helper`.
#'
#' Known limitation: averages over all n_studies including the study being compared against it,
#' rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
#' as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
#' the fix.
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
#' @param included_n_reps a integer matrix. Count of non-NaN replicates (included ones) per reference point, per study
#'   The minimum valid value is `0`.
#' @return a named list with elements:
#'   \item{mean_pmf}{a numeric matrix. The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
#'     known-limitation note above}
#'   \item{mean_pmf_included_n_reps}{a integer vector. Count of non-NaN replicates (included ones) per reference point for `mean_pmf`,
#'     summed across all n_studies}
#'   \item{mean_pmf_counts}{a integer matrix. Absolute counts of a residual per bin for the mean pmf -> `sum(counts, dim=3)`}
#' @export
create_mean_pmf <- function(pmfs, counts, included_n_reps) {
    pmfs <- .tox_as_double_array(pmfs, "pmfs", 3L)
    counts <- .tox_as_integer_array(counts, "counts", 3L)
    included_n_reps <- .tox_as_integer_matrix(included_n_reps, "included_n_reps")
    if (dim(counts)[1] != dim(pmfs)[1])
        .tox_shape_error("counts", dim(counts)[1], "pmfs", dim(pmfs)[1])
    if (dim(counts)[2] != dim(pmfs)[2])
        .tox_shape_error("counts", dim(counts)[2], "pmfs", dim(pmfs)[2])
    if (dim(included_n_reps)[1] != dim(pmfs)[2])
        .tox_shape_error("included_n_reps", dim(included_n_reps)[1], "pmfs", dim(pmfs)[2])
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
#' Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_only_helper`: useful where the
#' mean pmf's own counts don't matter, e.g. for the bootstrap confidence interval a later
#' stage of this port adds.
#'
#' Known limitation: averages over all n_studies including the study being compared against it,
#' rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
#' as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
#' the fix.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::create_mean_pmf_only}, whose argument names
#' are the ones an error message reports.
#'
#' @param pmfs a numeric array of rank 3. Per-study probabilities of each bin per reference point, from
#'   \code{\link{build_residual_histograms}}
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @return a numeric matrix. The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
#'   known-limitation note above
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
#' Ported from 125-stabilize-jscomp's `bootstrap_histogram_helper`. Resamples from the
#' POOLED/consensus histogram counts (`mean_pmf_counts`), not from each study's own histogram
#' -- this is 125's deliberate design, preserved as-is (see the plan for this port). Draws
#' random numbers via
#' \code{random_multinomial}, so this implementation is
#' deliberately impure, matching the project's existing precedent for the other
#' permutation-test module's purity
#' (\code{\link{perform_permutation_test}}).
#'
#' `confidence_interval` is both an input and an output: its incoming `[lower, upper]` values
#' seed every slot of the top-k/bottom-k heaps (`tmp_bootstrapping_top_k_jsds`) -- ported
#' verbatim from 125, which fills the whole heap with the *same* incoming reference value
#' rather than the usual plus/minus-infinity heap initialization, so the observed
#' (pre-bootstrap) value can only be displaced by a strictly more extreme bootstrap draw. On
#' return it holds `[largest of the n_bootstrapping_top_k_jsds smallest bootstrap draws,
#' smallest of the n_bootstrapping_top_k_jsds largest bootstrap draws]`.
#'
#' `mean_pmf_counts`/`tmp_pmfs`/`tmp_mean_pmf` are laid out bin-major (`(n_bins, n_points[,
#' n_studies])`), matching
#' \code{\link{create_mean_pmf}} and its
#' own `create_mean_pmf_only_impl` above -- both ported from 125, whose own convention this
#' is. The already-shipped
#' \code{\link{compute_divergence_per_reference_point}}
#' predates 125's own port and instead takes its pmf arguments point-major
#' (`(n_points, n_bins)`); the two calls below bridge the two conventions with an explicit
#' `transpose`, rather than picking one shape and silently reinterpreting the other's memory
#' under it (which would scramble every non-square `(n_bins, n_points)` histogram).
#'
#' A GSL allocation failure in `create_rng` is a genuine runtime error no input check could
#' have foreseen (codegen_guide.md Sec 5.14): every work array and `confidence_interval` are
#' then left untouched (arrays not yet written to keep their caller-visible defined state) and
#' `ierr` reports `ERR_ALLOC_FAIL`. A `random_multinomial` draw failing (which validated,
#' internally-consistent inputs should never trigger) is likewise folded into `ierr`, first
#' failure only, without stopping the resampling already in flight.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::bootstrap_histogram}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bootstraps a integer scalar. Number of bootstrap resamples to perform
#'   The minimum valid value is `1`.
#' @param mean_pmf_counts a integer matrix. Absolute counts of a residual per bin for the pooled/consensus pmf, from
#'   create_mean_pmf_impl -- resampled with replacement each bootstrap
#'   The minimum valid value is `0`.
#' @param mean_pmf_included_n_reps a integer vector. Count of non-NaN replicates (included ones) per reference point for the pooled pmf
#'   The minimum valid value is `0`.
#' @param included_n_reps a integer matrix. Count of non-NaN replicates (included ones) per reference point, per study --
#'   how many elements are drawn (with replacement) from the pooled pool per study
#'   The minimum valid value is `0`.
#' @param confidence_interval a numeric matrix. Confidence interval to be bootstrapped -- incoming values are the reference values
#'   that seed the top-k/bottom-k heaps (see above); overwritten with the bootstrapped
#'   `[lower, upper]` interval per study
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Forwarded to calc_js_comp_test_n_top_k_jsds to size n_bootstrapping_top_k_jsds; not
#'   otherwise used here
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @return a numeric matrix. Confidence interval to be bootstrapped -- incoming values are the reference values
#'   that seed the top-k/bottom-k heaps (see above); overwritten with the bootstrapped
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
    if (length(mean_pmf_included_n_reps) != dim(mean_pmf_counts)[2])
        .tox_shape_error("mean_pmf_included_n_reps", length(mean_pmf_included_n_reps), "mean_pmf_counts", dim(mean_pmf_counts)[2])
    if (dim(included_n_reps)[1] != dim(mean_pmf_counts)[2])
        .tox_shape_error("included_n_reps", dim(included_n_reps)[1], "mean_pmf_counts", dim(mean_pmf_counts)[2])
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
#' Ported from 125-stabilize-jscomp's `js_comp_test_helper`, restructured for Issue #187's
#' occupancy-constrained per-neighborhood histogram binning into three passes, mirroring
#' \code{\link{run_js_comp_test_parameter_search}}'s
#' own Pass A/B/C split (Issue #187's own Steps 2.5/2.6):
#'
#' - Pass A (per study): builds every study's neighborhoods
#' (\code{\link{construct_neighborhoods_ranged}}),
#' writing into `neighborhood_indices`/`neighborhood_range`, which already retain every
#' study's own values simultaneously (both are real `intent(out)` arguments sized
#' `(..., n_points, n_studies)` -- unlike `run_js_comp_test_parameter_search_impl`, no new
#' buffer was needed for this). Unlike that routine, there is no admissibility gate here, so
#' Pass A always runs to completion for every study.
#' - Pass B (per point, sequential -- see the implementation body's own comment for why): pools
#' every study's residuals for one reference point at a time
#' (\code{\link{gather_pooled_neighborhood_residuals}},
#' now published as its own entry point) and runs Issue #187's occupancy-constrained
#' bin-count search on the pooled result
#' (\code{\link{determine_bin_count_occupancy}}),
#' deciding `n_bins_per_point(i_point)` independently for every reference point, plus the
#' `occupancy_failed`/`n_pooled_residuals`/`min_bin_occupancy`/`mean_bin_occupancy`/
#' `max_bin_occupancy`/`sturges_bins`/`fd_bins` diagnostics. `max_n_bins_per_point`
#' (`maxval(n_bins_per_point(1:n_points))`) is derived once after Pass B and replaces the old
#' caller-supplied scalar `n_bins` everywhere downstream.
#' - Pass C (per study): re-gathers this study's residual values from the neighbor indices Pass
#' A already computed, then builds its residual histograms at the real per-point bin counts
#' (\code{\link{build_residual_histograms}}).
#'
#' `n_neighbors` may not exceed `max_n_genes_all_studies`, so every neighbor gene index Pass A
#' produces lies in `[1, max_n_genes_all_studies]`. Pass B still checks each index and, should
#' one ever fall outside, reports `ERR_INVALID_INPUT` and returns before Pass C could read
#' `residuals` out of bounds -- a safeguard that valid input cannot reach.
#'
#' After Pass C, the pipeline continues exactly as before: pools the per-study pmfs into the
#' consensus pmf
#' (\code{\link{create_mean_pmf}}), computes
#' each study's observed JSD against that consensus
#' (\code{\link{compute_divergence_per_reference_point}}/\code{\link{compute_weighted_global_divergence}},
#' called with the consensus pmf as the second argument), runs the permutation test
#' (\code{\link{gjct_permutation_test}}) -- which
#' returns the primary `p_values` (each permuted study against the consensus of the permuted
#' studies) and, from the same permutations, `p_values_observed_consensus` (each permuted study
#' against the observed consensus `mean_pmf`, kept for comparison) -- and
#' finally re-derives each study's pmf/JSD/weights/global JSD from its own UNTOUCHED `counts`
#' via
#' \code{\link{calc_pmf}} -- `mean_pmf`/`mean_pmf_counts`
#' are NOT re-derived, since the permutation test never modifies the observed consensus (it
#' resamples only its own scratch copies, never `mean_pmf_counts` itself, and builds each
#' permutation's consensus in its own scratch too), exactly as 125 relies on.
#'
#' **Behavioral asymmetry vs.
#' \code{\link{run_js_comp_test_parameter_search}}
#' -- read before using this entry point where inadequately-supported neighborhoods must be
#' rejected:** unlike that routine, THIS one has NO multi-candidate fallback and NO
#' admissibility gate at all (no
#' \code{\link{check_neighborhood_overlaps}},
#' no \code{\link{check_mean_pmf_min_counts}},
#' no early exit). A reference point whose Pass B occupancy search fails even at `m_min`
#' (`occupancy_failed(i_point) = TRUE_c_bool`) still gets a real histogram built, at
#' `n_bins_per_point(i_point) == m_min`, and that point still contributes to
#' `global_js_divergence` exactly like every other point -- its contribution is down-weighted
#' only by `included_n_reps` (an orthogonal quantity: how many non-NaN replicates it has), never
#' by bin sparsity. A caller that needs inadequately-supported neighborhoods rejected outright
#' should use `run_js_comp_test_parameter_search_impl` instead, which gates on exactly this via
#' `check_mean_pmf_min_counts_impl`.
#'
#' **A real, deliberate change to this routine's public array shapes (Issue #187):** the old
#' mandatory scalar input `n_bins` is gone -- there is no way for a caller to know the right bin
#' count in advance, since it is now genuinely computed inside this routine by Pass B's
#' occupancy search, independently per reference point. Every array whose bin-sized dimension
#' used to be sized by that input (`pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`,
#' `tmp_counts_point_major`, `tmp_pmf_point_major`, `tmp_permutation_mean_pmf_counts`,
#' `tmp_permutation_counts`, `tmp_permutation_pmfs`) is now sized to the fixed compile-time
#' ceiling \code{MAX_N_BINS} (`256`)
#' instead, exactly mirroring how `run_js_comp_test_parameter_search_impl`'s own
#' `tmp_counts_point_major`/`tmp_pmf_point_major`/`tmp_pmfs`/`tmp_counts` etc. have been sized
#' since Issue #187's earlier steps. The new `max_n_bins_per_point` output tells a caller how many of the
#' LEADING bins/rows of each of those arrays are actually meaningful
#' (`maxval(n_bins_per_point(1:n_points))`); the rest is unused padding. The generator's own
#' result-size trimming directive cannot express this trim, because it only ever trims an
#' array's LAST declared extent, and bins is the FIRST declared extent of every one of those
#' arrays -- so a Python/R caller must slice `[:max_n_bins_per_point, ...]` themselves, exactly as a
#' caller of `run_js_comp_test_parameter_search_impl`'s own jagged `trace_*` arrays already has
#' to.
#'
#' `x_star` is an ordinary input here, not computed by this routine -- 125's own
#' `js_comp_test_helper` takes it the same way, since a caller running several studies/several
#' parameter settings is expected to compute the reference points once
#' (\code{\link{pool_means}}) and reuse
#' them consistently.
#'
#' `construct_neighborhoods_ranged_impl` reports neighbor gene INDICES, not gathered residual
#' values (unlike its distance-sort sibling
#' \code{\link{construct_neighborhoods}}),
#' so Pass C gathers each neighbor's actual residual values from `residuals` itself
#' (one slice of `tmp_neighborhood_residuals_gathered` per reference point, reused for every
#' study) before building that point's histogram with `build_residual_histograms_impl` for
#' that single point, whose bins are copied straight into the matching column of
#' `counts`/`pmfs` -- no transpose is needed there. `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts`
#' here are BIN-major (`(256, n_points, n_studies)`) to match
#' \code{\link{create_mean_pmf}}'s own
#' convention, while `compute_divergence_per_reference_point_impl` and `calc_pmf_impl` are
#' POINT-major (`(n_points, max_n_bins_per_point)`): the calls to those two after Pass C bridge
#' with an explicit `transpose`, exactly as
#' \code{\link{bootstrap_histogram}} and
#' \code{\link{gjct_permutation_test}} already do.
#'
#' Impure: calls the impure `gjct_permutation_test_impl`. A GSL failure it reports is returned
#' in `ierr`, and the routine returns right there: `pmfs`, `js_divergences`, `weights` and
#' `global_js_divergence` then hold the pre-permutation values, not the final re-derived ones.
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
#'   NaN is permitted for this value.
#' @param n_permutations a integer scalar. Number of permutations, forwarded to gjct_permutation_test_impl
#'   The minimum valid value is `0`.
#'   The default value is `1000`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach for a candidate bin count to
#'   be admissible in Pass B's occupancy search, forwarded to
#'   determine_bin_count_occupancy_impl
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param m_min a integer scalar. Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
#'   forwarded to determine_bin_count_occupancy_impl
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest candidate bin count Pass B's occupancy search will ever test (M_max),
#'   forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
#'   determine_bin_count_occupancy_impl clamps it up to `m_min` internally
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor for Pass B's occupancy search's coarse search stage,
#'   forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
#'   advances
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] for each reference point's own lower residual-range bound,
#'   forwarded to determine_bin_count_occupancy_impl
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] for each reference point's own upper residual-range bound,
#'   forwarded to determine_bin_count_occupancy_impl
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @return a named list with elements:
#'   \item{neighborhood_indices}{a integer array of rank 3. Gene indices of the selected neighborhood, per reference point, per study (Pass A)}
#'   \item{neighborhood_range}{a integer array of rank 3. For each reference point and study, the `[min_idx, max_idx]` neighborhood span, as
#'     produced by construct_neighborhoods_ranged_impl (Pass A)}
#'   \item{n_bins_per_point}{a integer vector. This reference point's own selected histogram bin count (Issue #187's `M_j`), from
#'     Pass B's occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood
#'     may use a different bin count}
#'   \item{shared_residual_range_low}{a numeric vector. This reference point's own lower residual-range bound (R_low), from Pass B's
#'     occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood may use a
#'     different, asymmetric range (Step 3)}
#'   \item{shared_residual_range_high}{a numeric vector. This reference point's own upper residual-range bound (R_high), from Pass B's
#'     occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood may use a
#'     different, asymmetric range (Step 3)}
#'   \item{max_n_bins_per_point}{a integer scalar. The widest `n_bins_per_point` value across all `n_points` reference points
#'     (`maxval(n_bins_per_point(1:n_points))`), derived once after Pass B. The number of
#'     leading, meaningful bins/rows in `pmfs`, `counts`, `mean_pmf` and `mean_pmf_counts`
#'     below (and of the bin-sized work arrays) -- those are all declared
#'     with a fixed 256-bin ceiling (MAX_N_BINS) rather than a caller-supplied bin count,
#'     since `n_bins_per_point` can no longer be known by a caller in advance. A Python/R
#'     caller must slice `[:max_n_bins_per_point, ...]` themselves: the generator's own result-size
#'     trimming directive cannot express this trim, because it only ever trims an array's
#'     LAST declared extent, and bins is the FIRST declared extent of every one of those
#'     arrays}
#'   \item{occupancy_failed}{a logical vector. This reference point's `occupancy_failed` flag from Pass B
#'     (determine_bin_count_occupancy_impl) -- `TRUE` iff even `m_min` bins could not
#'     satisfy the occupancy criterion for it. See this routine's own doc block above for
#'     the behavioral asymmetry this implies vs. run_js_comp_test_parameter_search_impl: a
#'     `TRUE` point here still gets a real histogram and still contributes to
#'     `global_js_divergence`, it is never rejected}
#'   \item{n_pooled_residuals}{a integer vector. This reference point's pooled residual count (N_j) from Pass B
#'     (determine_bin_count_occupancy_impl)}
#'   \item{min_bin_occupancy}{a integer vector. This reference point's minimum bin occupancy at `n_bins_per_point`, from Pass B
#'     (determine_bin_count_occupancy_impl)}
#'   \item{mean_bin_occupancy}{a numeric vector. This reference point's mean bin occupancy at `n_bins_per_point`, from Pass B
#'     (determine_bin_count_occupancy_impl)}
#'   \item{max_bin_occupancy}{a integer vector. This reference point's maximum bin occupancy at `n_bins_per_point`, from Pass B
#'     (determine_bin_count_occupancy_impl)}
#'   \item{sturges_bins}{a integer vector. This reference point's Sturges' rule bin-count diagnostic from Pass B
#'     (determine_bin_count_occupancy_impl) -- never part of the occupancy search's own
#'     decision}
#'   \item{fd_bins}{a integer vector. This reference point's Freedman-Diaconis rule bin-count diagnostic from Pass B
#'     (determine_bin_count_occupancy_impl) -- never part of the occupancy search's own
#'     decision}
#'   \item{pmfs}{a numeric array of rank 3. `counts` normalized to `0 <= pmfs(:, :, i) <= 1` and `sum(pmfs(:, j, i)) == 1`. `256`
#'     = MAX_N_BINS, a fixed ceiling (see `max_n_bins_per_point` above) -- only rows `1:max_n_bins_per_point`
#'     are meaningful; a Python/R caller must slice `[:max_n_bins_per_point, ...]` themselves}
#'   \item{counts}{a integer array of rank 3. Absolute counts of a residual per bin for `pmfs`. `256` = MAX_N_BINS; only rows
#'     `1:max_n_bins_per_point` are meaningful -- see `pmfs` above}
#'   \item{included_n_reps}{a integer matrix. Count of non-NaN replicates (included ones) per reference point, per study}
#'   \item{mean_pmf}{a numeric matrix. The consensus pmf, from create_mean_pmf_impl. `256` = MAX_N_BINS; only rows
#'     `1:max_n_bins_per_point` are meaningful -- see `pmfs` above}
#'   \item{mean_pmf_counts}{a integer matrix. Absolute counts of a residual per bin for the consensus pmf. `256` = MAX_N_BINS;
#'     only rows `1:max_n_bins_per_point` are meaningful -- see `pmfs` above}
#'   \item{mean_pmf_included_n_reps}{a integer vector. Count of non-NaN replicates (included ones) per reference point for the consensus pmf}
#'   \item{js_divergences}{a numeric matrix. Per-reference-point JSD of each study against the consensus pmf}
#'   \item{weights}{a numeric matrix. Per-reference-point weights for `global_js_divergence`}
#'   \item{global_js_divergence}{a numeric vector. Weighted global JSD of each study against the consensus pmf}
#'   \item{p_values}{a numeric vector. Empirical p-value per study from gjct_permutation_test_impl, each permuted study
#'     compared against the consensus of the permuted studies -- the primary p-value}
#'   \item{p_values_observed_consensus}{a numeric vector. Empirical p-value per study from the same permutations, but with each permuted study
#'     compared against the fixed observed consensus `mean_pmf`, as the permutation test did
#'     before it recomputed the consensus per permutation. Kept for comparison only;
#'     `p_values` is the primary result}
#' @export
run_js_comp_test <- function(n_neighbors, gene_means, gene_means_perms, residuals, x_star, n_permutations = 1000L, random_seed = 42L, min_residuals_per_bin = 10L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95) {
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
    if (dim(gene_means_perms)[2] != dim(gene_means)[2])
        .tox_shape_error("gene_means_perms", dim(gene_means_perms)[2], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(gene_means_perms)[1] != dim(gene_means)[1])
        .tox_shape_error("gene_means_perms", dim(gene_means_perms)[1], "gene_means", dim(gene_means)[1])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_call", n_neighbors, gene_means, gene_means_perms, residuals, x_star, n_permutations, random_seed, min_residuals_per_bin, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "n_neighbors", "gene_means", "gene_means_perms", "residuals", "x_star", "neighborhood_indices", "neighborhood_range", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "p_values_observed_consensus", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", "x_star", "neighborhood_indices", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
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
#' same pipeline and the same outputs, on neighborhoods whose size varies per reference point
#' and per study. Call
#' \code{\link{construct_adaptive_neighborhoods}}
#' first, then pass its `pooled_neighborhood_range` (trimmed to its `n_points` reference
#' points) here together with the same `gene_means` and `residuals`.
#'
#' **Membership.** The pooled gene means are sorted again exactly as the construction sorts them
#' (the same routine on the same input), so a pooled position means the same entry here as it
#' did there. Every range `a..b` is decoded back to per-study gene lists by
#' \code{materialize_pooled_neighborhood};
#' the per-study counts are returned as `n_neighbors_per_point`, in the construction's own
#' `(n_studies, n_points)` orientation, and add up to `b - a + 1`.
#'
#' **Pipeline, per reference point.** Pass B pools the point's residuals across all its genes
#' of all studies and picks its bin count and residual range with Issue #187's occupancy
#' search
#' (\code{\link{determine_bin_count_occupancy}}),
#' its Sturges/Freedman-Diaconis diagnostics computed with the rounded mean per-study neighbor
#' count. Pass C builds each study's histogram from that study's own genes of the point. Then,
#' exactly as in `run_js_comp_test_impl`: the consensus pmf
#' (\code{\link{create_mean_pmf}}), each
#' study's JSD against it and the weighted global JSD, the permutation test
#' (\code{\link{gjct_permutation_test}}, with
#' both its `p_values` and its comparison `p_values_observed_consensus`) and the
#' final re-derivation from the untouched `counts`. Neighborhoods that decode to the same gene
#' sets as a fixed-k run's therefore give bit-identical results. Every point is weighted by its
#' non-NaN residual count, which under adaptive growth genuinely differs between points.
#' Like `run_js_comp_test_impl`, there is no admissibility gate: overlap and occupancy are the
#' parameter search's concern.
#'
#' **Runtime errors.** Before any histogram work, every range is checked in point order; the
#' first failing one sets `ERR_INVALID_INPUT` and the routine returns at once, every output
#' undefined. For each point, in this order: its first position lies after its last; its last
#' position lies beyond the non-NaN pooled means (the entries a construction can ever use); a
#' study has no gene in it, which would leave that study's pmf empty and the consensus pmf, an
#' average over all studies, meaningless. A range a construction produced passes the first two
#' checks by construction, and the third whenever its status was not an empty-study one.
#'
#' **Memory.** The work arrays are sized by the per-study gene bound `max_n_genes_all_studies`,
#' the largest number of genes one study can have in one neighborhood, not by a caller-supplied
#' neighbor count. The pooling buffers `tmp_pooled_residuals`/`tmp_pooled_residuals_perm`
#' therefore hold `max_n_reps_all_studies * max_n_genes_all_studies * n_studies` values each,
#' the size of `residuals` itself; Pass B and Pass C are sequential over points so one such
#' buffer serves every point.
#'
#' The bin-sized arrays have the same fixed 256-bin
#' (\code{MAX_N_BINS}) leading extent as in
#' `run_js_comp_test_impl`, of which only the first `max_n_bins_per_point` rows are
#' meaningful; a Python/R caller slices `[:max_n_bins_per_point, ...]` themselves.
#'
#' Impure: calls the impure `gjct_permutation_test_impl`. A GSL failure it reports is returned
#' in `ierr`, and the routine returns right there: `pmfs`, `js_divergences`, `weights` and
#' `global_js_divergence` then hold the pre-permutation values, not the final re-derived ones.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test::run_js_comp_test_adaptive}, whose argument names
#' are the ones an error message reports.
#'
#' @param gene_means a numeric matrix. Mean expression of every gene in every study, NaN for a missing gene; the same
#'   array the neighborhoods were constructed from
#'   NaN is permitted for this value.
#' @param residuals a numeric array of rank 3. Signed residuals of every replicate of every gene in every study, NaN for a missing value
#'   NaN is permitted for this value.
#' @param pooled_neighborhood_range a integer matrix. For each neighborhood, its first and last position in the ascending order of the
#'   pooled means, as returned by construct_adaptive_neighborhoods
#'   The minimum valid value is `1`.
#'   The maximum valid value is `max_n_genes_all_studies*n_studies`.
#' @param n_permutations a integer scalar. Number of permutations, forwarded to gjct_permutation_test_impl
#'   The minimum valid value is `0`.
#'   The default value is `1000`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach for a candidate bin count to
#'   be admissible in the occupancy search
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param m_min a integer scalar. Smallest candidate bin count the occupancy search tests (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest candidate bin count the occupancy search tests (M_max); raised to `m_min`
#'   when smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor of the occupancy search's coarse stage
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] for each reference point's lower residual-range bound
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] for each reference point's upper residual-range bound
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @return a named list with elements:
#'   \item{n_neighbors_per_point}{a integer matrix. For each neighborhood, how many of its genes belong to each study (at least 1 each)}
#'   \item{n_bins_per_point}{a integer vector. Each reference point's selected histogram bin count, from the occupancy search}
#'   \item{shared_residual_range_low}{a numeric vector. Each reference point's lower residual-range bound (R_low), from the occupancy search}
#'   \item{shared_residual_range_high}{a numeric vector. Each reference point's upper residual-range bound (R_high), from the occupancy search}
#'   \item{max_n_bins_per_point}{a integer scalar. The widest `n_bins_per_point` value, `maxval(n_bins_per_point)`: the number of leading,
#'     meaningful rows of `pmfs`, `counts`, `mean_pmf` and `mean_pmf_counts`}
#'   \item{occupancy_failed}{a logical vector. `TRUE` iff even `m_min` bins could not satisfy the occupancy criterion for this
#'     reference point; such a point still gets a histogram at `m_min` bins and still
#'     contributes to `global_js_divergence`}
#'   \item{n_pooled_residuals}{a integer vector. Each reference point's pooled non-NaN residual count (N_j), across all its genes of
#'     all studies}
#'   \item{min_bin_occupancy}{a integer vector. Each reference point's minimum bin occupancy at `n_bins_per_point`}
#'   \item{mean_bin_occupancy}{a numeric vector. Each reference point's mean bin occupancy at `n_bins_per_point`}
#'   \item{max_bin_occupancy}{a integer vector. Each reference point's maximum bin occupancy at `n_bins_per_point`}
#'   \item{sturges_bins}{a integer vector. Each reference point's Sturges' rule diagnostic, with the rounded mean per-study
#'     neighbor count; never part of the occupancy search's decision}
#'   \item{fd_bins}{a integer vector. Each reference point's Freedman-Diaconis rule diagnostic, with the rounded mean
#'     per-study neighbor count; never part of the occupancy search's decision}
#'   \item{pmfs}{a numeric array of rank 3. `counts` normalized per reference point and study. `256` = MAX_N_BINS; only rows
#'     `1:max_n_bins_per_point` are meaningful}
#'   \item{counts}{a integer array of rank 3. Absolute counts of a residual per bin for `pmfs`. `256` = MAX_N_BINS; only rows
#'     `1:max_n_bins_per_point` are meaningful}
#'   \item{included_n_reps}{a integer matrix. Count of non-NaN residuals binned per reference point, per study}
#'   \item{mean_pmf}{a numeric matrix. The consensus pmf. `256` = MAX_N_BINS; only rows `1:max_n_bins_per_point` are meaningful}
#'   \item{mean_pmf_counts}{a integer matrix. Absolute counts of a residual per bin for the consensus pmf. `256` = MAX_N_BINS; only
#'     rows `1:max_n_bins_per_point` are meaningful}
#'   \item{mean_pmf_included_n_reps}{a integer vector. Count of non-NaN residuals per reference point for the consensus pmf}
#'   \item{js_divergences}{a numeric matrix. Per-reference-point JSD of each study against the consensus pmf}
#'   \item{weights}{a numeric matrix. Per-reference-point weights for `global_js_divergence`}
#'   \item{global_js_divergence}{a numeric vector. Weighted global JSD of each study against the consensus pmf}
#'   \item{p_values}{a numeric vector. Empirical p-value per study from the permutation test, each permuted study compared
#'     against the consensus of the permuted studies -- the primary p-value}
#'   \item{p_values_observed_consensus}{a numeric vector. Empirical p-value per study from the same permutations, but with each permuted study
#'     compared against the fixed observed consensus `mean_pmf`, as the permutation test did
#'     before it recomputed the consensus per permutation. Kept for comparison only;
#'     `p_values` is the primary result}
#' @export
run_js_comp_test_adaptive <- function(gene_means, residuals, pooled_neighborhood_range, n_permutations = 1000L, random_seed = 42L, min_residuals_per_bin = 10L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95) {
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
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_adaptive_call", gene_means, residuals, pooled_neighborhood_range, n_permutations, random_seed, min_residuals_per_bin, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "gene_means", "residuals", "pooled_neighborhood_range", "n_neighbors_per_point", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "p_values_observed_consensus", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", "pooled_neighborhood_range", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
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
#' Ported from 125-stabilize-jscomp's `determine_js_comp_test_n_points_n_neighbors_helper` and
#' `_alloc`, merged into one implementation now that the new `_impl` rules leave no separate
#' hand-written allocation layer. Pools all studies' gene means, sorts them once
#' (\code{sort_real_heapsort_expl_size}), generates the candidate
#' grid from the gene count alone
#' (\code{\link{generate_js_comp_test_candidates}}
#' -- Issue #187, Step 2.7: this no longer needs the pooled residuals, since real per-neighborhood
#' bin counts are decided later, in Pass B below),
#' then walks it from finest to coarsest resolution: for each candidate, builds every study's
#' neighborhoods and checks the first admissibility gate
#' (\code{\link{check_neighborhood_overlaps}});
#' once every study passes, pools the consensus pmf and checks the second gate
#' (\code{\link{check_mean_pmf_min_counts}});
#' once that passes too, seeds a confidence interval with the observed JSD, bootstraps it
#' (\code{\link{bootstrap_histogram}}), and
#' tests it against the running best candidate for a plateau
#' (\code{\link{check_plateau_condition}}).
#' `plateau_mode` picks which of that CI-overlap criterion and Issue #178's complementary
#' relative-effect-size one
#' (\code{\link{check_effect_size_plateau_condition}})
#' governs the stop condition; both are always computed and traced (`trace_*` below) once a
#' candidate is admissible, regardless of `plateau_mode`, so a caller can compare what either
#' criterion would have decided. The search stops (`exit`) the moment the SELECTED criterion's
#' plateau is found -- see `plateau_mode`'s own mode table below for the accepted values.
#'
#' Issue #178 also names 2 blocking dependencies for validating the effect-size thresholds
#' empirically -- the KX_FACTORS default and the Freedman-Diaconis bin-count overestimate --
#' both deliberately left as-is here; see the project's JSD-Comp-Test follow-up issue. (A third
#' candidate blocker, a one-sided-vs-symmetric JSD formula question, was raised in the same
#' follow-up issue but confirmed by the issue's own author to be a mistake in the issue text,
#' not a real discrepancy -- the code's symmetric formula is correct as written.)
#' `delta_median_threshold`/`delta_max_threshold` default to the issue's own suggested (not yet
#' validated) 0.05/0.10.
#'
#' When `plateau_mode` selects the effect-size criterion (`MODE_PLATEAU_EFFECT_SIZE` or
#' `MODE_PLATEAU_BOTH`) and it plateaus independently of the CI-overlap criterion's own running
#' "best candidate" bookkeeping, `best_candidate_index`/`best_candidate_pair_confidence_interval`
#' are overridden to the candidate that actually triggered the effect-size plateau, so the
#' candidate this routine returns is always the one that stopped the search.
#'
#' The finally returned candidate is decided by one of three cases:
#'
#' 1. A plateau was found, or the grid collapsed to a single candidate (see
#' \code{\link{generate_js_comp_test_candidates}}'s
#' own small-N collapse note): the best candidate is returned and `plateau_established` is
#' `TRUE`. A lone candidate is used regardless of whether it plateaued or even passed
#' either gate -- the plateau machinery is bypassed entirely, exactly as 125 does. If that
#' lone candidate never passed both gates, no real per-point values exist for it, so every
#' point's bin count is set to `m_min`, both residual ranges to `0.0`, and
#' `best_candidate_pair_confidence_interval` stays `-1.0`.
#' 2. No plateau, but at least one candidate was admissible: the admissible candidate with the
#' smallest bootstrapped uncertainty (median confidence-interval width across studies) is
#' returned, with its real confidence interval, and `plateau_established` is `FALSE`. This
#' applies to every `plateau_mode`.
#' 3. No plateau, and no candidate was ever admissible: the search falls back to the FIRST
#' (finest-resolution) candidate, resets `best_candidate_pair_confidence_interval` to
#' `-1.0`, sets every point's bin count to `m_min` and both residual ranges to `0.0`, and
#' `plateau_established` is `FALSE`.
#'
#' Per the plan's work-array translation for this routine specifically: `max_n_bins_all_candidates`
#' (data-dependent, not cheaply closed-form in 125) is replaced by the fixed
#' \code{MAX_N_BINS} ceiling, so every
#' bin-dimensioned work array below is sized to MAX_N_BINS and sliced `(1:max_n_bins, ...)` per
#' candidate, rather than carrying a separate recommend-sized dimension argument for it.
#' `max_n_bins` is the widest per-point bin count Issue #187's occupancy search (Pass B below)
#' chose for the current candidate, `maxval(tmp_n_bins_per_point(1:n_points))` -- it replaces
#' the old single scalar `n_bins` that used to come from the global-pool Sturges/FD estimate.
#' `gene_means` is passed to
#' \code{sort_real_heapsort_expl_size}
#' as its own multi-dimensional self -- that callee declares its matching dummy with an
#' explicit shape, so standard Fortran sequence association reinterprets the contiguous actual
#' argument as the flat 1-D array it expects, exactly as 125's own `_alloc` layer did for the
#' same call.
#'
#' Impure: calls the impure
#' \code{\link{bootstrap_histogram}}. A GSL
#' failure it reports is folded into `ierr` (first failure only) without aborting the search,
#' matching that routine's own tolerant precedent.
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
#' @param join_method a string, one of "join_min", "join_max", "join_median". The way to evaluate all studies' confidence-interval overlaps for the plateau
#'   condition, forwarded to check_plateau_condition_impl
#' @param max_n_points_candidate a integer scalar. Exact upper bound on the grid's first (largest) `n_points` candidate. Issue #187's
#'   per-point outputs below (`n_bins_per_point`, `trace_selected_n_bins`, and the other
#'   jagged `trace_*` arrays) are sized by this argument, so unlike before Issue #187 it
#'   is no longer purely an internal sizing detail the plain wrapper can compute and
#'   hide -- the caller must know it up front to receive those arrays, hence JUST_INFO
#'   rather than AUTO here now
#'   It is recommended to compute this argument from the `max_n_points_candidate` output produced by \code{\link{calc_js_comp_test_candidate_bounds}}.
#'   The minimum valid value is `1`.
#' @param max_n_neighbors_candidate a integer scalar. Safe upper bound on the grid's largest `n_neighbors` candidate. Also JUST_INFO, not
#'   because anything returned is sized by it (nothing is), but because it comes from the
#'   same `calc_js_comp_test_candidate_bounds` call as `max_n_points_candidate` above --
#'   now that that call can no longer run automatically inside this wrapper, splitting
#'   this one back into an AUTO call would just be a second, redundant call to the same
#'   routine for no benefit
#'   It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by \code{\link{calc_js_comp_test_candidate_bounds}}.
#'   The minimum valid value is `1`.
#' @param min_residuals_per_bin a integer scalar. Minimum count each bin of the consensus pmf must reach to pass the second
#'   admissibility gate. Reuses Issue #187's occupancy-search default rather than an
#'   independently-tunable threshold of its own: once
#'   \code{\link{determine_bin_count_occupancy}}
#'   wires real per-neighborhood bin counts in, a separate laxer threshold here would
#'   silently let a candidate the occupancy search already marked `occupancy_failed`
#'   pass this gate anyway, defeating the FAILURE-detection mechanism
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param min_neighbor_overlap a numeric scalar. Minimum fractional overlap two consecutive neighborhoods must have to pass the first
#'   admissibility gate
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.1`.
#' @param succeeding_ci_overlap a numeric scalar. Minimum fractional overlap a candidate's confidence interval must have with the
#'   running best, per `join_method`, to plateau
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @param plateau_mode a string, one of "plateau_ci_overlap", "plateau_effect_size", "plateau_both". Which plateau criterion decides when the search stops
#'
#'   The default value is `"plateau_ci_overlap"`.
#' @param delta_median_threshold a numeric scalar. Upper bound the median relative JSD change across studies must stay under for a
#'   transition to count toward an effect-size plateau, forwarded to
#'   check_effect_size_plateau_condition_impl
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.05`.
#' @param delta_max_threshold a numeric scalar. Upper bound the largest relative JSD change across studies must stay under for a
#'   transition to count toward an effect-size plateau, forwarded to
#'   check_effect_size_plateau_condition_impl
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.10`.
#' @param delta_epsilon a numeric scalar. Small constant preventing division by zero when a study's previous admissible
#'   candidate's JSD was zero, forwarded to check_effect_size_plateau_condition_impl
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `1.0e-10`.
#' @param delta_min_consecutive_transitions a integer scalar. Number of consecutive qualifying transitions required to declare an effect-size
#'   plateau, forwarded to check_effect_size_plateau_condition_impl
#'   The minimum valid value is `1`.
#'   The default value is `2`.
#' @param m_min a integer scalar. Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
#'   forwarded to determine_bin_count_occupancy_impl
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest candidate bin count Pass B's occupancy search will ever test (M_max),
#'   forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
#'   determine_bin_count_occupancy_impl clamps it up to `m_min` internally
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor for Pass B's occupancy search's coarse search stage,
#'   forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
#'   advances
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] for each reference point's own lower residual-range bound,
#'   forwarded to determine_bin_count_occupancy_impl
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] for each reference point's own upper residual-range bound,
#'   forwarded to determine_bin_count_occupancy_impl
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
#'   to bootstrap_histogram_impl itself
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @return a named list with elements:
#'   \item{n_points}{a integer scalar. The finally chosen candidate's `n_points`}
#'   \item{n_neighbors}{a integer scalar. The finally chosen candidate's `n_neighbors`}
#'   \item{n_bins_per_point}{a integer vector. The finally chosen candidate's per-point histogram bin count, one per reference
#'     point (Issue #187: every neighborhood may use a different bin count). Only the
#'     leading `n_points` entries are meaningful, mirroring how `n_points`/`n_neighbors`
#'     above are the finally chosen candidate's own values}
#'   \item{shared_residual_range_low}{a numeric vector. The finally chosen candidate's per-point lower residual-range bound (R_low), one per
#'     reference point (Step 3: every neighborhood may use a different, asymmetric range).
#'     Only the leading `n_points` entries are meaningful, mirroring `n_bins_per_point`
#'     above; `0.0` throughout in the two genuinely-degenerate cases where Pass B
#'     never ran for the returned candidate (cases 1 and 3 of the routine description)}
#'   \item{shared_residual_range_high}{a numeric vector. The finally chosen candidate's per-point upper residual-range bound (R_high),
#'     mirroring `shared_residual_range_low` above in every respect}
#'   \item{best_candidate_pair_confidence_interval}{a numeric matrix. The bootstrapped JSD confidence interval for the finally chosen candidate pair;
#'     `-1.0` throughout when no candidate was ever admissible: either the lone
#'     candidate of a collapsed grid never passed both gates (`plateau_established` is
#'     `TRUE`), or no plateau was found and no smallest-bootstrap-uncertainty candidate
#'     could be substituted (`plateau_established` is `FALSE`)}
#'   \item{plateau_established}{a logical scalar. `TRUE` when a real plateau was found (by whichever criterion
#'     `plateau_mode` selected) or the candidate grid never had more than one candidate to
#'     begin with. `FALSE` when the search exhausted every admissible candidate
#'     without ever finding one -- Issue #178's own "report that parameter stability could
#'     not be established". When `FALSE` and at least one candidate was admissible, the
#'     routine still returns a real (non-`-1.0`) candidate and confidence interval: the
#'     admissible candidate with the smallest bootstrapped uncertainty, per the issue's own
#'     fallback recommendation -- `plateau_established` is what distinguishes that case from
#'     an actual plateau, not the confidence interval's sentinel value}
#'   \item{trace_n_points}{a integer vector. Per-admissible-candidate `n_points`, one entry per column of the other `trace_*`
#'     arrays. `16` = MAX_CANDIDATE_PAIRS, written as a literal for the same reason
#'     candidates_n_points_n_neighbors is in generate_js_comp_test_candidates_impl
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_n_neighbors}{a integer vector. Per-admissible-candidate `n_neighbors`, paired with trace_n_points above
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_global_js_divergence}{a numeric matrix. Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
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
#'     (`W_{i,t} / J_{i,t}`), denominator floored at `delta_epsilon` -- the issue's own
#'     formula omits this floor, but the same near-zero-JSD instability that motivates
#'     `delta_epsilon` in the `Delta_{i,t}` formula applies here too (a near-zero `J`
#'     destabilizes any ratio that divides by it, whichever candidate's `J` it is)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta}{a numeric matrix. Per-admissible-candidate, per-study relative JSD change from the previous admissible
#'     candidate (`Delta_{i,t}`), from check_effect_size_plateau_condition_impl;
#'     `-1.0` throughout at the first admissible candidate specifically (no
#'     predecessor to diff against) -- every other column within `1:n_admissible_evaluated`
#'     holds a real value
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta_median}{a numeric vector. Per-admissible-candidate median of trace_delta across studies (Delta-tilde_t);
#'     `-1.0` at the first admissible candidate, see trace_delta above
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta_max}{a numeric vector. Per-admissible-candidate maximum of trace_delta across studies (Delta^max_t);
#'     `-1.0` at the first admissible candidate, see trace_delta above
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_selected_n_bins}{a integer matrix. Per-admissible-candidate, per-reference-point selected histogram bin count
#'     (Issue #187's `M_j`), from determine_bin_count_occupancy_impl. Unlike every OTHER
#'     trace_* array above, whose first extent is a fixed thing like `n_studies`, this
#'     array's first extent is `max_n_points_candidate`, NOT `n_points`, because `n_points`
#'     itself varies per candidate (that is why `trace_n_points(16)` exists as its own
#'     array): this array is genuinely JAGGED per candidate column `t` -- only rows
#'     `1:trace_n_points(t)` are meaningful for that column, rows beyond that are undefined
#'     padding. The result-size directive below only trims the LAST extent (candidates, via
#'     `n_admissible_evaluated`), not this row dimension, so a Python/R caller must
#'     additionally slice `[:trace_n_points[t], t]` themselves
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_occupancy_failed}{a logical matrix. Per-admissible-candidate, per-reference-point `occupancy_failed` flag from
#'     determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
#'     trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
#'     column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_n_pooled_residuals}{a integer matrix. Per-admissible-candidate, per-reference-point pooled residual count (`N_j`) from
#'     determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
#'     trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
#'     column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_min_bin_occupancy}{a integer matrix. Per-admissible-candidate, per-reference-point minimum bin occupancy at
#'     trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
#'     column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
#'     meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
#'     themselves
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_mean_bin_occupancy}{a numeric matrix. Per-admissible-candidate, per-reference-point mean bin occupancy at
#'     trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
#'     column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
#'     meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
#'     themselves
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_max_bin_occupancy}{a integer matrix. Per-admissible-candidate, per-reference-point maximum bin occupancy at
#'     trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
#'     column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
#'     meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
#'     themselves
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_sturges_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Sturges' rule bin-count diagnostic
#'     from determine_bin_count_occupancy_impl (never part of the occupancy search's own
#'     decision). Jagged per candidate column exactly as trace_selected_n_bins above --
#'     only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R caller
#'     must slice `[:trace_n_points[t], t]` themselves
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_fd_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Freedman-Diaconis rule bin-count
#'     diagnostic from determine_bin_count_occupancy_impl (never part of the occupancy
#'     search's own decision). Jagged per candidate column exactly as trace_selected_n_bins
#'     above -- only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R
#'     caller must slice `[:trace_n_points[t], t]` themselves
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_low}{a numeric matrix. Per-admissible-candidate, per-reference-point lower residual-range bound (R_low)
#'     from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
#'     trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
#'     column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_high}{a numeric matrix. Per-admissible-candidate, per-reference-point upper residual-range bound (R_high)
#'     from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
#'     trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
#'     column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
#'     The first `n_admissible_evaluated` elements will hold the results.}
#' @export
run_js_comp_test_parameter_search <- function(gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, max_n_neighbors_candidate, min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_ci_overlap", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L) {
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
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_parameter_search_call", gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, max_n_neighbors_candidate, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "n_bootstraps", "join_method", "max_n_points_candidate", "max_n_neighbors_candidate", "n_points", "n_neighbors", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "best_candidate_pair_confidence_interval", "plateau_established", "n_admissible_evaluated", "trace_n_points", "trace_n_neighbors", "trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_delta", "trace_delta_median", "trace_delta_max", "trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy", "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins", "trace_shared_residual_range_low", "trace_shared_residual_range_high", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", NA_character_, NA_character_, NA_character_, NA_character_, "n_bins_per_point", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
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
#' 1. Grow its neighborhoods over the pooled, sorted gene means. The pooled means are sorted
#' once for the whole search, with the same routine on the same input as
#' \code{\link{construct_adaptive_neighborhoods}},
#' and the growth is that routine's own, with the same `tau`, `mad_distance_factor`,
#' `max_pooled_residuals` and `min_study_neighbors`. So a candidate's neighborhoods are
#' exactly what `construct_adaptive_neighborhoods` returns for its knobs.
#' 2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
#' reference points emerged.
#' 3. First admissibility gate:
#' \code{\link{check_neighborhood_overlaps}}
#' on the pooled ranges `[a_i, b_i]`.
#' 4. Per point, decode the range into per-study gene lists and run Issue #187's occupancy
#' search on the pooled residuals (Pass B). Then build each study's histogram per point
#' (Pass C), the consensus pmf, and apply the second gate,
#' \code{\link{check_mean_pmf_min_counts}}.
#' This is the pipeline of
#' \code{\link{run_js_comp_test_adaptive}}.
#' 5. From here on everything is the fixed-k search's post-gate bookkeeping, keyed on the
#' candidate index: the observed JSD and the traces, the effect-size check (before the
#' bootstrap, which reuses the observed-JSD buffer as scratch), the bootstrap, the
#' smallest-uncertainty tracking, the CI-overlap check, the `plateau_mode` selection, the
#' effect-size override, and the stop at the first plateau.
#'
#' **Candidate log.** Every candidate tried gets an entry in `candidate_k_start`,
#' `candidate_n_points` and `candidate_status`, in order, until the search stops. The status is
#'
#' - \code{ADAPTIVE_CANDIDATE_EVALUATED}
#' (`0`): passed both gates and has a trace column;
#' - \code{ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED}
#' (`1`): construction status too few means or too
#' few residuals, so the growth could not run as specified;
#' - \code{ADAPTIVE_CANDIDATE_EMPTY_STUDY}
#' (`2`): construction status empty study neighborhood;
#' - \code{ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED}
#' (`3`): more than `max_n_points_candidate` points;
#' - \code{ADAPTIVE_CANDIDATE_OVERLAP_FAILED}
#' (`4`): the first gate failed;
#' - \code{ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED}
#' (`5`): the second gate failed.
#'
#' The checks run in that order and the first failing one is logged. `candidate_n_points` is the
#' number of reference points the construction produced (0 for too few means), whatever the
#' status.
#'
#' **Result**, one of three cases, as in the fixed-k search:
#'
#' 1. A plateau was found, or the sequence had a single candidate: that candidate's knobs,
#' `plateau_established = TRUE`.
#' 2. No plateau but at least one admissible candidate: the admissible candidate with the
#' smallest median confidence-interval width, `plateau_established = FALSE`.
#' 3. Nothing admissible: candidate 1's knobs, confidence interval `-1.0`,
#' `plateau_established = FALSE`.
#'
#' `n_points` is the selected candidate's own logged point count. Where no admissible candidate
#' backs the result (case 3, or case 1 with a single, rejected candidate), every
#' `n_bins_per_point` entry up to `min(n_points, max_n_points_candidate)` is `m_min` and both
#' ranges are `0.0`, as in the fixed-k search. To rebuild the selected neighborhoods, call
#' `construct_adaptive_neighborhoods` with the returned knobs and the same construction
#' options, then `run_js_comp_test_adaptive` with the same occupancy options: its
#' `global_js_divergence` equals the selected candidate's `trace_global_js_divergence` column.
#'
#' **Traces** are the fixed-k search's, with `trace_k_start`/`trace_k_step`/`trace_k_max` in
#' place of `trace_n_neighbors`. The per-point traces are jagged per column exactly as there.
#'
#' **Memory.** Beyond the fixed-k search's per-point arrays, which are sized by
#' `max_n_points_candidate` (practical capacity from
#' \code{\link{calc_adaptive_js_comp_test_bounds}}),
#' the construction needs buffers over all `N = max_n_genes_all_studies * n_studies` pooled
#' entries, the largest being `tmp_n_neighbors_per_point` with `n_studies * N` integers. The
#' pooling buffers `tmp_pooled_residuals`/`tmp_pooled_residuals_perm` hold
#' `max_n_reps_all_studies * N` values each, the size of `residuals`, because one study may
#' have all its genes in one neighborhood. The histogram buffers hold `256 * max_n_points_candidate
#' * n_studies` values each.
#'
#' Impure: calls the impure
#' \code{\link{bootstrap_histogram}}. A GSL
#' failure it reports is folded into `ierr` (first failure only) without aborting the search,
#' as in the fixed-k search.
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
#' @param join_method a string, one of "join_min", "join_max", "join_median". The way to evaluate all studies' confidence-interval overlaps for the plateau
#'   condition, forwarded to check_plateau_condition_impl
#' @param max_n_points_candidate a integer scalar. Most reference points a candidate may emerge with; a candidate with more is rejected
#'   with the capacity status. The per-point outputs and traces below are sized by it, so
#'   the caller has to know it up front
#'   It is recommended to compute this argument from the `max_n_points_candidate` output produced by \code{\link{calc_adaptive_js_comp_test_bounds}}.
#'   The minimum valid value is `1`.
#' @param min_residuals_per_bin a integer scalar. Minimum count every bin must reach, both in the occupancy search and in the second
#'   admissibility gate (the consensus pmf)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param min_neighbor_overlap a numeric scalar. Minimum fractional overlap two consecutive neighborhoods' pooled ranges must have to
#'   pass the first admissibility gate
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.1`.
#' @param succeeding_ci_overlap a numeric scalar. Minimum fractional overlap a candidate's confidence interval must have with the
#'   running best, per `join_method`, to plateau
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @param plateau_mode a string, one of "plateau_ci_overlap", "plateau_effect_size", "plateau_both". Which plateau criterion decides when the search stops
#'
#'   The default value is `"plateau_ci_overlap"`.
#' @param delta_median_threshold a numeric scalar. Upper bound the median relative JSD change across studies must stay under for a
#'   transition to count toward an effect-size plateau
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.05`.
#' @param delta_max_threshold a numeric scalar. Upper bound the largest relative JSD change across studies must stay under for a
#'   transition to count toward an effect-size plateau
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.10`.
#' @param delta_epsilon a numeric scalar. Small constant preventing division by zero when a study's previous admissible JSD was
#'   zero
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `1.0e-10`.
#' @param delta_min_consecutive_transitions a integer scalar. Number of consecutive qualifying transitions required to declare an effect-size plateau
#'   The minimum valid value is `1`.
#'   The default value is `2`.
#' @param m_min a integer scalar. Smallest candidate bin count the occupancy search tests (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest candidate bin count the occupancy search tests (M_max); raised to `m_min`
#'   when smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor of the occupancy search's coarse stage
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] for each reference point's lower residual-range bound
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] for each reference point's upper residual-range bound
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
#'   to bootstrap_histogram_impl itself
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param tau a numeric scalar. Largest relative increase of the dispersion an adaptive growth round may cause and
#'   still be committed, as in construct_adaptive_neighborhoods
#'   The minimum valid value is `0.0`.
#'   The default value is `0.1`.
#' @param mad_distance_factor a numeric scalar. Multiple of a neighborhood's median absolute deviation that the next reference
#'   point's target lies beyond it, as in construct_adaptive_neighborhoods
#'   The minimum valid value is `0.0`.
#'   The default value is `1.0`.
#' @param max_pooled_residuals a integer scalar. Cap on the non-NaN residuals adaptive rounds may grow a neighborhood's pool to; 0
#'   for no cap, as in construct_adaptive_neighborhoods
#'   The minimum valid value is `0`.
#'   The default value is `0`.
#' @param min_study_neighbors a integer scalar. Fewest entries of each study every neighborhood must have for a candidate to be
#'   admissible, as in construct_adaptive_neighborhoods
#'   The minimum valid value is `1`.
#'   The default value is `1`.
#' @return a named list with elements:
#'   \item{k_start}{a integer scalar. The finally chosen candidate's `k_start`}
#'   \item{k_step}{a integer scalar. The finally chosen candidate's `k_step`}
#'   \item{k_max}{a integer scalar. The finally chosen candidate's `k_max`}
#'   \item{n_points}{a integer scalar. Number of reference points the finally chosen candidate's construction produced}
#'   \item{n_bins_per_point}{a integer vector. The finally chosen candidate's per-point histogram bin count; only the leading
#'     `n_points` entries are meaningful. Where no admissible candidate backs the result,
#'     the leading `min(n_points, max_n_points_candidate)` entries are `m_min`, since that
#'     candidate's point count may exceed the capacity}
#'   \item{shared_residual_range_low}{a numeric vector. The finally chosen candidate's per-point lower residual-range bound (R_low); only the
#'     leading `n_points` entries are meaningful, `0.0` where no admissible candidate
#'     backs the result}
#'   \item{shared_residual_range_high}{a numeric vector. The finally chosen candidate's per-point upper residual-range bound (R_high),
#'     mirroring `shared_residual_range_low`}
#'   \item{best_candidate_confidence_interval}{a numeric matrix. The bootstrapped JSD confidence interval of the finally chosen candidate;
#'     `-1.0` throughout where no admissible candidate backs the result}
#'   \item{plateau_established}{a logical scalar. `TRUE` when a plateau was found (by the criterion `plateau_mode` selects) or the
#'     sequence had a single candidate; `FALSE` when the search exhausted the sequence
#'     without one, in which case the smallest-uncertainty admissible candidate is
#'     returned if there is one}
#'   \item{trace_k_start}{a integer vector. Per-admissible-candidate `k_start`. `16` = MAX_CANDIDATE_PAIRS, the adaptive
#'     sequence's cap
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_k_step}{a integer vector. Per-admissible-candidate `k_step`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_k_max}{a integer vector. Per-admissible-candidate `k_max`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_n_points}{a integer vector. Per-admissible-candidate number of reference points
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_global_js_divergence}{a numeric matrix. Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_lower}{a numeric matrix. Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_upper}{a numeric matrix. Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_width}{a numeric matrix. Per-admissible-candidate, per-study confidence-interval width
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_width_relative}{a numeric matrix. Per-admissible-candidate, per-study confidence-interval width divided by the observed
#'     global JSD, the denominator floored at `delta_epsilon`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta}{a numeric matrix. Per-admissible-candidate, per-study relative JSD change from the previous admissible
#'     candidate; `-1.0` throughout at the first one
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta_median}{a numeric vector. Per-admissible-candidate median of trace_delta across studies; `-1.0` at the
#'     first one
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta_max}{a numeric vector. Per-admissible-candidate maximum of trace_delta across studies; `-1.0` at the
#'     first one
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_selected_n_bins}{a integer matrix. Per-admissible-candidate, per-reference-point selected histogram bin count. Jagged:
#'     only rows `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller
#'     must additionally slice `[:trace_n_points[t], t]` themselves
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_occupancy_failed}{a logical matrix. Per-admissible-candidate, per-reference-point `occupancy_failed` flag; jagged as
#'     trace_selected_n_bins
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
#'   \item{trace_sturges_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Sturges' rule diagnostic, with the
#'     rounded mean per-study neighbor count; jagged as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_fd_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic, with
#'     the rounded mean per-study neighbor count; jagged as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_low}{a numeric matrix. Per-admissible-candidate, per-reference-point lower residual-range bound (R_low);
#'     jagged as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_high}{a numeric matrix. Per-admissible-candidate, per-reference-point upper residual-range bound (R_high);
#'     jagged as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{candidate_k_start}{a integer vector. `k_start` of every candidate tried, in sequence order
#'     The first `n_candidates_tried` elements will hold the results.}
#'   \item{candidate_n_points}{a integer vector. Number of reference points every tried candidate's construction produced (0 for too
#'     few means), whether or not the candidate was admissible
#'     The first `n_candidates_tried` elements will hold the results.}
#'   \item{candidate_status}{a integer vector. What happened to every tried candidate: `0` evaluated,
#'     `1` construction failed,
#'     `2` a study neighborhood below `min_study_neighbors`,
#'     `3` capacity exceeded,
#'     `4` overlap gate failed,
#'     `5` consensus-pmf count gate failed
#'     The first `n_candidates_tried` elements will hold the results.}
#' @export
run_js_comp_test_adaptive_parameter_search <- function(gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_ci_overlap", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, tau = 0.1, mad_distance_factor = 1.0, max_pooled_residuals = 0L, min_study_neighbors = 1L) {
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
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_adaptive_parameter_search_call", gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, tau, mad_distance_factor, max_pooled_residuals, min_study_neighbors)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "n_bootstraps", "join_method", "max_n_points_candidate", "k_start", "k_step", "k_max", "n_points", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "best_candidate_confidence_interval", "plateau_established", "n_admissible_evaluated", "trace_k_start", "trace_k_step", "trace_k_max", "trace_n_points", "trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_delta", "trace_delta_median", "trace_delta_max", "trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy", "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins", "trace_shared_residual_range_low", "trace_shared_residual_range_high", "n_candidates_tried", "candidate_k_start", "candidate_n_points", "candidate_status", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", NA_character_, NA_character_, NA_character_, NA_character_, "n_bins_per_point", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
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
