# Generated. Do not edit.

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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_search::check_neighborhood_overlaps}, whose argument names
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_search::check_mean_pmf_min_counts}, whose argument names
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_search::check_plateau_condition}, whose argument names
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_search::check_effect_size_plateau_condition}, whose argument names
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
