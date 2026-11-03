#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_search(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter-Search Criteria
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
module tox_data_integration_js_comp_test_search_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_as_view
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: check_neighborhood_overlaps_c
    public :: check_mean_pmf_min_counts_c
    public :: check_plateau_condition_c
    public :: check_effect_size_plateau_condition_c
    public :: check_effect_size_plateau_condition_expert_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_search(module):check_neighborhood_overlaps(subroutine)]]
    !| The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
    !| neighborhood spans of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
    !| number of positions they share divided by the length of the earlier span; a reversed span
    !| (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
    !| `test_`-prefixed export as a test case.
    subroutine check_neighborhood_overlaps_c(&
            neighborhood_range,&
            n_points,&
            min_neighbor_overlap,&
            all_have_min_neighbor_overlap,&
            ierr&
        ) bind(C, name="check_neighborhood_overlaps_c")
        use tox_data_integration_js_comp_test_search, only: check_neighborhood_overlaps

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), dimension(2, n_points), intent(in), target :: neighborhood_range
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        logical(c_bool), intent(out), target :: all_have_min_neighbor_overlap
            !! `.true.` if every pair of consecutive neighborhoods overlaps by at least
            !! `min_neighbor_overlap`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(min_neighbor_overlap)
        M_CHECK_NON_NULL(all_have_min_neighbor_overlap)
        M_CHECK_ARRAY_NON_NULL(neighborhood_range, 2 * n_points)

        call check_neighborhood_overlaps(&
            neighborhood_range = neighborhood_range,&
            n_points = n_points,&
            min_neighbor_overlap = min_neighbor_overlap,&
            all_have_min_neighbor_overlap = all_have_min_neighbor_overlap,&
            ierr = ierr&
        )
    end subroutine check_neighborhood_overlaps_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_search(module):check_mean_pmf_min_counts(subroutine)]]
    !| The second admissibility gate of a parameter-search candidate, checked after
    !| [[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl(interface)]].
    !| Each point is checked over its own `n_bins_per_point` bins only, so the zero padding beyond
    !| them does not count as a failure.
    subroutine check_mean_pmf_min_counts_c(&
            mean_pmf_counts,&
            n_bins,&
            n_bins_per_point,&
            n_points,&
            min_count,&
            all_bins_have_min_count,&
            ierr&
        ) bind(C, name="check_mean_pmf_min_counts_c")
        use tox_data_integration_js_comp_test_search, only: check_mean_pmf_min_counts

        integer(c_int), intent(in), target :: n_bins
            !! Second extent of `mean_pmf_counts`: the widest bin count of any point
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), dimension(n_points, n_bins), intent(in), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: n_bins_per_point
            !! Each point's own bin count; only `mean_pmf_counts(i_point, 1:n_bins_per_point(i_point))`
            !! is inspected
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_bins`.
        integer(c_int), intent(in), target :: min_count
            !! Minimum count each bin of the mean pmf must reach
            !! The minimum valid value is `0_int32`.
        logical(c_bool), intent(out), target :: all_bins_have_min_count
            !! `.true.` if every bin within each reference point's own `n_bins_per_point`, at every
            !! reference point, reaches at least `min_count`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(min_count)
        M_CHECK_NON_NULL(all_bins_have_min_count)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, n_points)

        call check_mean_pmf_min_counts(&
            mean_pmf_counts = mean_pmf_counts,&
            n_bins = n_bins,&
            n_bins_per_point = n_bins_per_point,&
            n_points = n_points,&
            min_count = min_count,&
            all_bins_have_min_count = all_bins_have_min_count,&
            ierr = ierr&
        )
    end subroutine check_mean_pmf_min_counts_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_search(module):check_plateau_condition(subroutine)]]
    !| For each study, the overlap of the new candidate's interval with the best one so far is
    !| measured as a fraction of the new interval's width, and counted when it reaches
    !| `succeeding_ci_overlap`. A plateau is found when fewer studies pass than for the current
    !| best (the best is then kept), or, after the new candidate became the best, when the count
    !| meets the condition `join_method` names.
    subroutine check_plateau_condition_c(&
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
        ) bind(C, name="check_plateau_condition_c")
        use tox_data_integration_js_comp_test_search, only: check_plateau_condition
        use tox_data_integration_js_comp_test_search_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(2, n_studies), intent(in), target :: confidence_interval
            !! Bootstrapped JSD confidence interval `[lower, upper]` of the candidate under test
        real(c_double), dimension(2, n_studies), intent(inout), target :: best_candidate_pair_confidence_interval
            !! JSD confidence intervals for the current best candidate pair; overwritten with
            !! `confidence_interval` unless the new candidate is worse
        integer(c_int), intent(inout), target :: best_candidate_index
            !! Candidate-grid index of the current best candidate pair; overwritten with
            !! `candidate_index` unless the new candidate is worse
        integer(c_int), intent(inout), target :: best_exceeded_ci_overlap_count
            !! Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
            !! best candidate pair; overwritten unless the new candidate is worse
            !! The minimum valid value is `0_int32`.
        integer(c_int), intent(in), target :: candidate_index
            !! Candidate-grid index of the candidate pair that produced `confidence_interval`
            !! The minimum valid value is `1_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
            !! How the studies' overlaps are joined into the plateau condition; the median method
            !! requires `count > (n_studies - 1) / 2`
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Minimum fractional overlap with the best candidate's interval for a study to pass
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), dimension(n_studies), intent(out), target :: ci_overlap
            !! Per-study overlap, in `[0, 1]`, of `confidence_interval` with the best interval as
            !! it was on entry, relative to the new interval's width; always filled. Meaningless
            !! for the very first candidate, whose best interval is only a seed
        logical(c_bool), intent(out), target :: plateau_found
            !! `.true.` once the new candidate is no better than the previous best, or once
            !! `join_method`'s overlap condition is met by the new candidate
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32) :: join_method_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(best_candidate_index)
        M_CHECK_NON_NULL(best_exceeded_ci_overlap_count)
        M_CHECK_NON_NULL(candidate_index)
        M_CHECK_NON_NULL(succeeding_ci_overlap)
        M_CHECK_NON_NULL(plateau_found)
        M_CHECK_ARRAY_NON_NULL(confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(best_candidate_pair_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)
        M_CHECK_ARRAY_NON_NULL(ci_overlap, n_studies)

        block
            character(len=:), pointer :: join_method_f
            join_method_f => c_char_as_view(join_method)

            select case (join_method_f)
                case ("join_min")
                    join_method_mode_f = METHOD_JOIN_MIN
                case ("join_max")
                    join_method_mode_f = METHOD_JOIN_MAX
                case ("join_median")
                    join_method_mode_f = METHOD_JOIN_MEDIAN
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call check_plateau_condition(&
            confidence_interval = confidence_interval,&
            best_candidate_pair_confidence_interval = best_candidate_pair_confidence_interval,&
            n_studies = n_studies,&
            best_candidate_index = best_candidate_index,&
            best_exceeded_ci_overlap_count = best_exceeded_ci_overlap_count,&
            candidate_index = candidate_index,&
            join_method = join_method_mode_f,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            ci_overlap = ci_overlap,&
            plateau_found = plateau_found,&
            ierr = ierr&
        )
    end subroutine check_plateau_condition_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_search(module):check_effect_size_plateau_condition(subroutine)]]
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
    subroutine check_effect_size_plateau_condition_c(&
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
        ) bind(C, name="check_effect_size_plateau_condition_c")
        use tox_data_integration_js_comp_test_search, only: check_effect_size_plateau_condition

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(n_studies), intent(in), target :: global_js_divergence
            !! Current admissible candidate's observed global JSD per study
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_studies), intent(in), target :: prev_global_js_divergence
            !! Previous admissible candidate's observed global JSD per study; ignored when
            !! `has_previous` is `.false.`
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), intent(in), target :: has_previous
            !! `.false.` for the first admissible candidate, which has no transition
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound on the median relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound on the largest relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(inout), target :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; reset to zero when this one does
            !! not qualify or `has_previous` is `.false.`
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_studies), intent(out), target :: delta
            !! Per-study relative JSD change from the previous admissible candidate; `-1.0_real64`
            !! throughout iff `.not. has_previous`
        real(c_double), intent(out), target :: delta_median
            !! Median of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        real(c_double), intent(out), target :: delta_max
            !! Maximum of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        logical(c_bool), intent(out), target :: plateau_found
            !! `.true.` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(has_previous)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(n_consecutive_ok)
        M_CHECK_NON_NULL(delta_median)
        M_CHECK_NON_NULL(delta_max)
        M_CHECK_NON_NULL(plateau_found)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(prev_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(delta, n_studies)

        call check_effect_size_plateau_condition(&
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
            ierr = ierr&
        )
    end subroutine check_effect_size_plateau_condition_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_search(module):check_effect_size_plateau_condition_expert(subroutine)]]
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
    subroutine check_effect_size_plateau_condition_expert_c(&
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
        ) bind(C, name="check_effect_size_plateau_condition_expert_c")
        use tox_data_integration_js_comp_test_search, only: check_effect_size_plateau_condition_expert

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(n_studies), intent(in), target :: global_js_divergence
            !! Current admissible candidate's observed global JSD per study
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_studies), intent(in), target :: prev_global_js_divergence
            !! Previous admissible candidate's observed global JSD per study; ignored when
            !! `has_previous` is `.false.`
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), intent(in), target :: has_previous
            !! `.false.` for the first admissible candidate, which has no transition
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound on the median relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound on the largest relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(inout), target :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; reset to zero when this one does
            !! not qualify or `has_previous` is `.false.`
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_studies), intent(out), target :: delta
            !! Per-study relative JSD change from the previous admissible candidate; `-1.0_real64`
            !! throughout iff `.not. has_previous`
        real(c_double), intent(out), target :: delta_median
            !! Median of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        real(c_double), intent(out), target :: delta_max
            !! Maximum of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        logical(c_bool), intent(out), target :: plateau_found
            !! `.true.` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_delta_perm
            !! Working array: sorting permutation for `delta`, used to compute `delta_median`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(has_previous)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(n_consecutive_ok)
        M_CHECK_NON_NULL(delta_median)
        M_CHECK_NON_NULL(delta_max)
        M_CHECK_NON_NULL(plateau_found)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(prev_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(delta, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_delta_perm, n_studies)

        call check_effect_size_plateau_condition_expert(&
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
            tmp_delta_perm = tmp_delta_perm,&
            ierr = ierr&
        )
    end subroutine check_effect_size_plateau_condition_expert_c

end module tox_data_integration_js_comp_test_search_c
#endif
