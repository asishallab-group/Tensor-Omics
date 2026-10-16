#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter-Search Criteria
!|
!| What the two parameter searches of the JS-Comp-Test have in common: the fixed-k search
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
!| and the adaptive search
!| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]].
!| Both accept a candidate when its neighborhoods overlap enough and every consensus bin holds
!| enough residuals, and both stop at the first such candidate whose bootstrapped JSDs have
!| reached a plateau.
!|
!| The admissibility gates are
!| [[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl(interface)]]
!| (consecutive neighborhoods overlap by at least `min_neighbor_overlap`) and
!| [[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]]
!| (every consensus bin reaches the minimum count). The plateau checks are
!| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]
!| (the confidence-interval overlap with the best candidate so far, joined over the studies by
!| `join_method`) and
!| [[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]]
!| (the relative change of each study's observed JSD between successive candidates); a search's
!| `plateau_mode` selects which of them stops it. Both searches evaluate at most
!| [[tox_data_integration_js_comp_test_search_impl(module):MAX_CANDIDATE_PAIRS(variable)]] candidates.
!|
!| Generated from [[tox_data_integration_js_comp_test_search_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_js_comp_test_search
    use f42_safeguard
    use tox_data_integration_js_comp_test_search_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, check_effect_size_plateau_condition_impl
    use tox_data_integration_js_comp_test_search_impl, only: check_mean_pmf_min_counts_impl, check_neighborhood_overlaps_impl, check_plateau_condition_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_math_impl, only: above
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    use tox_errors, only: set_err, set_err_once, validate_all_in_range_int, validate_all_in_range_real
    use tox_errors, only: validate_dimension_size, validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: check_neighborhood_overlaps
    public :: check_mean_pmf_min_counts
    public :: check_plateau_condition
    public :: check_effect_size_plateau_condition
    public :: check_effect_size_plateau_condition_expert

contains

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl]].
    !| The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
    !| neighborhood spans of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
    !| number of positions they share divided by the length of the earlier span; a reversed span
    !| (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
    !| `test_`-prefixed export as a test case.
    pure subroutine check_neighborhood_overlaps(&
            neighborhood_range,&
            n_points,&
            min_neighbor_overlap,&
            all_have_min_neighbor_overlap,&
            ierr&
        )
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), dimension(2, n_points), intent(in) :: neighborhood_range
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span
            !! The minimum valid value is `1_int32`.
        real(real64), intent(in) :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        logical(c_bool), intent(out) :: all_have_min_neighbor_overlap
            !! `.true.` if every pair of consecutive neighborhoods overlaps by at least
            !! `min_neighbor_overlap`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_points, ierr, arg_pos=2_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=3_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_int(neighborhood_range, 2 * n_points, ierr, arg_pos=1_int32, min=1_int32)
        if (is_err(ierr)) return
#endif

        call check_neighborhood_overlaps_impl(&
            neighborhood_range = neighborhood_range,&
            n_points = n_points,&
            min_neighbor_overlap = min_neighbor_overlap,&
            all_have_min_neighbor_overlap = all_have_min_neighbor_overlap&
        )
    end subroutine check_neighborhood_overlaps

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl]].
    !| The second admissibility gate of a parameter-search candidate, checked after
    !| [[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl(interface)]].
    !| Each point is checked over its own `n_bins_per_point` bins only, so the zero padding beyond
    !| them does not count as a failure.
    pure subroutine check_mean_pmf_min_counts(&
            mean_pmf_counts,&
            n_bins,&
            n_bins_per_point,&
            n_points,&
            min_count,&
            all_bins_have_min_count,&
            ierr&
        )
        integer(int32), intent(in) :: n_bins
            !! Second extent of `mean_pmf_counts`: the widest bin count of any point
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), dimension(n_points, n_bins), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points), intent(in) :: n_bins_per_point
            !! Each point's own bin count; only `mean_pmf_counts(i_point, 1:n_bins_per_point(i_point))`
            !! is inspected
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_bins`.
        integer(int32), intent(in) :: min_count
            !! Minimum count each bin of the mean pmf must reach
            !! The minimum valid value is `0_int32`.
        logical(c_bool), intent(out) :: all_bins_have_min_count
            !! `.true.` if every bin within each reference point's own `n_bins_per_point`, at every
            !! reference point, reaches at least `min_count`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_bins, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_in_range_int(min_count, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_int(mean_pmf_counts, n_points * n_bins, ierr, arg_pos=1_int32, min=0_int32)
        call validate_all_in_range_int(n_bins_per_point, n_points, ierr, arg_pos=3_int32, min=1_int32, max=n_bins)
        if (is_err(ierr)) return
#endif

        call check_mean_pmf_min_counts_impl(&
            mean_pmf_counts = mean_pmf_counts,&
            n_bins = n_bins,&
            n_bins_per_point = n_bins_per_point,&
            n_points = n_points,&
            min_count = min_count,&
            all_bins_have_min_count = all_bins_have_min_count&
        )
    end subroutine check_mean_pmf_min_counts

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl]].
    !| For each study, the overlap of the new candidate's interval with the best one so far is
    !| measured as a fraction of the new interval's width, and counted when it reaches
    !| `succeeding_ci_overlap`. A plateau is found when fewer studies pass than for the current
    !| best (the best is then kept), or, after the new candidate became the best, when the count
    !| meets the condition `join_method` names.
    pure subroutine check_plateau_condition(&
            confidence_interval,&
            best_candidate_pair_confidence_interval,&
            n_studies,&
            best_candidate_index,&
            best_exceeded_ci_overlap_count,&
            candidate_index,&
            join_method,&
            succeeding_ci_overlap,&
            ci_overlap,&
            plateau_found,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(2, n_studies), intent(in) :: confidence_interval
            !! Bootstrapped JSD confidence interval `[lower, upper]` of the candidate under test
        real(real64), dimension(2, n_studies), intent(inout) :: best_candidate_pair_confidence_interval
            !! JSD confidence intervals for the current best candidate pair; overwritten with
            !! `confidence_interval` unless the new candidate is worse
        integer(int32), intent(inout) :: best_candidate_index
            !! Candidate-grid index of the current best candidate pair; overwritten with
            !! `candidate_index` unless the new candidate is worse
        integer(int32), intent(inout) :: best_exceeded_ci_overlap_count
            !! Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
            !! best candidate pair; overwritten unless the new candidate is worse
            !! The minimum valid value is `0_int32`.
        integer(int32), intent(in) :: candidate_index
            !! Candidate-grid index of the candidate pair that produced `confidence_interval`
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: join_method
            !! How the studies' overlaps are joined into the plateau condition; the median method
            !! requires `count > (n_studies - 1) / 2`
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        real(real64), intent(in) :: succeeding_ci_overlap
            !! Minimum fractional overlap with the best candidate's interval for a study to pass
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), dimension(n_studies), intent(out) :: ci_overlap
            !! Per-study overlap, in `[0, 1]`, of `confidence_interval` with the best interval as
            !! it was on entry, relative to the new interval's width; always filled. Meaningless
            !! for the very first candidate, whose best interval is only a seed
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` once the new candidate is no better than the previous best, or once
            !! `join_method`'s overlap condition is met by the new candidate
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_studies, ierr, arg_pos=3_int32)
        call validate_in_range_int(best_exceeded_ci_overlap_count, ierr, arg_pos=5_int32, min=0_int32)
        call validate_in_range_int(candidate_index, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=8_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(confidence_interval, 2 * n_studies, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(best_candidate_pair_confidence_interval, 2 * n_studies, ierr, arg_pos=2_int32)
        if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (is_err(ierr)) return
#endif

        call check_plateau_condition_impl(&
            confidence_interval = confidence_interval,&
            best_candidate_pair_confidence_interval = best_candidate_pair_confidence_interval,&
            n_studies = n_studies,&
            best_candidate_index = best_candidate_index,&
            best_exceeded_ci_overlap_count = best_exceeded_ci_overlap_count,&
            candidate_index = candidate_index,&
            join_method = join_method,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            ci_overlap = ci_overlap,&
            plateau_found = plateau_found&
        )
    end subroutine check_plateau_condition

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test_search(module):check_effect_size_plateau_condition_expert]] to prepare it yourself.
    !| The relative-effect-size alternative to
    !| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]:
    !| for each study `i`, the relative change of the observed global JSD between successive
    !| admissible candidates,
    !| `delta(i) = |global_js_divergence(i) - prev_global_js_divergence(i)| / max(prev_global_js_divergence(i),
    !| delta_epsilon)`, is summarized by its median (`delta_median`) and maximum (`delta_max`). A
    !| plateau is found once both stay below their thresholds for
    !| `delta_min_consecutive_transitions` consecutive transitions, counted across calls in
    !| `n_consecutive_ok`.
    !|
    !| For the first admissible candidate (`has_previous = .false.`) there is no transition:
    !| `delta`, `delta_median` and `delta_max` are `-1.0` ("not computed").
    pure subroutine check_effect_size_plateau_condition(&
            global_js_divergence,&
            prev_global_js_divergence,&
            n_studies,&
            has_previous,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            n_consecutive_ok,&
            delta,&
            delta_median,&
            delta_max,&
            plateau_found,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(n_studies), intent(in) :: global_js_divergence
            !! Current admissible candidate's observed global JSD per study
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_studies), intent(in) :: prev_global_js_divergence
            !! Previous admissible candidate's observed global JSD per study; ignored when
            !! `has_previous` is `.false.`
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), intent(in) :: has_previous
            !! `.false.` for the first admissible candidate, which has no transition
        real(real64), intent(in) :: delta_median_threshold
            !! Upper bound on the median relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_max_threshold
            !! Upper bound on the largest relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(int32), intent(in) :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(inout) :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; reset to zero when this one does
            !! not qualify or `has_previous` is `.false.`
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_studies), intent(out) :: delta
            !! Per-study relative JSD change from the previous admissible candidate; `-1.0_real64`
            !! throughout iff `.not. has_previous`
        real(real64), intent(out) :: delta_median
            !! Median of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        real(real64), intent(out) :: delta_max
            !! Maximum of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_delta_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_studies, ierr, arg_pos=3_int32)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=5_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=6_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=7_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_int(n_consecutive_ok, ierr, arg_pos=9_int32, min=0_int32)
        call validate_all_in_range_real(global_js_divergence, n_studies, ierr, arg_pos=1_int32, min=0.0_real64)
        call validate_all_in_range_real(prev_global_js_divergence, n_studies, ierr, arg_pos=2_int32, min=0.0_real64)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_delta_perm(n_studies))

        call check_effect_size_plateau_condition_impl(&
            global_js_divergence = global_js_divergence,&
            prev_global_js_divergence = prev_global_js_divergence,&
            n_studies = n_studies,&
            has_previous = has_previous,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            n_consecutive_ok = n_consecutive_ok,&
            delta = delta,&
            delta_median = delta_median,&
            delta_max = delta_max,&
            plateau_found = plateau_found,&
            tmp_delta_perm = tmp_delta_perm&
        )
    end subroutine check_effect_size_plateau_condition

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test_search(module):check_effect_size_plateau_condition]] does both.
    !| The relative-effect-size alternative to
    !| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]:
    !| for each study `i`, the relative change of the observed global JSD between successive
    !| admissible candidates,
    !| `delta(i) = |global_js_divergence(i) - prev_global_js_divergence(i)| / max(prev_global_js_divergence(i),
    !| delta_epsilon)`, is summarized by its median (`delta_median`) and maximum (`delta_max`). A
    !| plateau is found once both stay below their thresholds for
    !| `delta_min_consecutive_transitions` consecutive transitions, counted across calls in
    !| `n_consecutive_ok`.
    !|
    !| For the first admissible candidate (`has_previous = .false.`) there is no transition:
    !| `delta`, `delta_median` and `delta_max` are `-1.0` ("not computed").
    pure subroutine check_effect_size_plateau_condition_expert(&
            global_js_divergence,&
            prev_global_js_divergence,&
            n_studies,&
            has_previous,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            n_consecutive_ok,&
            delta,&
            delta_median,&
            delta_max,&
            plateau_found,&
            tmp_delta_perm,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(n_studies), intent(in) :: global_js_divergence
            !! Current admissible candidate's observed global JSD per study
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_studies), intent(in) :: prev_global_js_divergence
            !! Previous admissible candidate's observed global JSD per study; ignored when
            !! `has_previous` is `.false.`
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), intent(in) :: has_previous
            !! `.false.` for the first admissible candidate, which has no transition
        real(real64), intent(in) :: delta_median_threshold
            !! Upper bound on the median relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_max_threshold
            !! Upper bound on the largest relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(int32), intent(in) :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(inout) :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; reset to zero when this one does
            !! not qualify or `has_previous` is `.false.`
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_studies), intent(out) :: delta
            !! Per-study relative JSD change from the previous admissible candidate; `-1.0_real64`
            !! throughout iff `.not. has_previous`
        real(real64), intent(out) :: delta_median
            !! Median of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        real(real64), intent(out) :: delta_max
            !! Maximum of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`
        integer(int32), dimension(n_studies), intent(out) :: tmp_delta_perm
            !! Working array: sorting permutation for `delta`, used to compute `delta_median`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_studies, ierr, arg_pos=3_int32)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=5_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=6_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=7_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_int(n_consecutive_ok, ierr, arg_pos=9_int32, min=0_int32)
        call validate_all_in_range_real(global_js_divergence, n_studies, ierr, arg_pos=1_int32, min=0.0_real64)
        call validate_all_in_range_real(prev_global_js_divergence, n_studies, ierr, arg_pos=2_int32, min=0.0_real64)
        if (is_err(ierr)) return
#endif

        call check_effect_size_plateau_condition_impl(&
            global_js_divergence = global_js_divergence,&
            prev_global_js_divergence = prev_global_js_divergence,&
            n_studies = n_studies,&
            has_previous = has_previous,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            n_consecutive_ok = n_consecutive_ok,&
            delta = delta,&
            delta_median = delta_median,&
            delta_max = delta_max,&
            plateau_found = plateau_found,&
            tmp_delta_perm = tmp_delta_perm&
        )
    end subroutine check_effect_size_plateau_condition_expert

end module tox_data_integration_js_comp_test_search
