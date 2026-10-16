#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

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
!| [[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]] is the second of those
!| admissibility gates. The plateau checks are
!| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]
!| (the confidence-interval overlap with the best candidate so far, joined over the studies by
!| `join_method`) and
!| [[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]]
!| (the relative change of each study's observed JSD between successive candidates); a search's
!| `plateau_mode` selects which of them stops it. Both searches evaluate at most
!| [[tox_data_integration_js_comp_test_search_impl(module):MAX_CANDIDATE_PAIRS(variable)]] candidates.
module tox_data_integration_js_comp_test_search_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use f42_stats_impl, only: calc_quantile_impl
    use f42_sort_impl, only: init_perm, sort_array_heapsort
    use tox_data_integration_jsd_impl, only: compute_divergence_per_reference_point_impl, &
                                             compute_weighted_global_divergence_impl, create_mean_pmf_impl
    M_IMPLICIT_NONE
    private
    public :: check_mean_pmf_min_counts_impl, check_plateau_condition_impl, check_effect_size_plateau_condition_impl, &
             score_candidate_histograms, METHOD_JOIN_MIN, METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, &
             MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE, MODE_PLATEAU_BOTH, MAX_CANDIDATE_PAIRS

    ! `join_method`'s mode table; the `METHOD_` prefix follows from the argument name `join_method`.
    integer(int32), parameter :: METHOD_JOIN_MIN = 0_int32
        !! Join method: succeeds only once every study's confidence-interval overlap exceeds the
        !! threshold, in [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]
    integer(int32), parameter :: METHOD_JOIN_MAX = 1_int32
        !! Join method: succeeds once any one study's confidence-interval overlap exceeds the
        !! threshold, in [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]
    integer(int32), parameter :: METHOD_JOIN_MEDIAN = 2_int32
        !! Join method: succeeds once a majority of studies' confidence-interval overlaps exceed
        !! the threshold, in [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]

    ! M_MAX_CANDIDATE_PAIRS is defined in src/tox/data_integration/macros.h: both searches need the
    ! number itself as the extent of their candidate and trace arrays.
    integer(int32), parameter :: MAX_CANDIDATE_PAIRS = M_MAX_CANDIDATE_PAIRS
        !! = size(KX_FACTORS) * MAX_POINT_CANDIDATES = 2 * 8 of the fixed-k search
        !! ([[tox_data_integration_js_comp_test_fixed_k_search_impl(module):KX_FACTORS(variable)]],
        !! [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MAX_POINT_CANDIDATES(variable)]]);
        !! the adaptive search caps its candidate count at the same value

    ! `plateau_mode`'s mode table. M_MODE_PLATEAU_CI_OVERLAP, defined in
    ! src/tox/data_integration/macros.h, is also both searches' DM_DEFAULT, which cannot name a
    ! Fortran parameter.
    integer(int32), parameter :: MODE_PLATEAU_CI_OVERLAP = M_MODE_PLATEAU_CI_OVERLAP
        !! Plateau mode: CI overlap only
    integer(int32), parameter :: MODE_PLATEAU_EFFECT_SIZE = 1_int32
        !! Plateau mode: relative-effect-size stability only
    integer(int32), parameter :: MODE_PLATEAU_BOTH = 2_int32
        !! Plateau mode: either CI overlap or relative-effect-size stability

contains

    !> summary: Test whether every bin of a mean pmf reaches a minimum absolute count
    !| AUTHOR_LASZLO_LANG
    !| The second admissibility gate of a parameter-search candidate, checked after
    !| [[tox_data_integration_preprocessing_impl(module):check_neighborhood_overlaps_impl(interface)]].
    !| Each point is checked over its own `n_bins_per_point` bins only, so the zero padding beyond
    !| them does not count as a failure.
    pure subroutine check_mean_pmf_min_counts_impl(mean_pmf_counts, n_bins, n_bins_per_point, n_points, min_count, &
                                                   all_bins_have_min_count)
        integer(int32), intent(in) :: n_bins
            !! Second extent of `mean_pmf_counts`: the widest bin count of any point
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! Each point's own bin count; only `mean_pmf_counts(i_point, 1:n_bins_per_point(i_point))`
            !! is inspected
            !! DM_MIN(1_int32)
            !! DM_MAX(n_bins)
        integer(int32), intent(in) :: mean_pmf_counts(n_points, n_bins)
            !! Absolute counts of a residual per bin for the mean pmf
            !! DM_MIN(0_int32)
        integer(int32), intent(in) :: min_count
            !! Minimum count each bin of the mean pmf must reach
            !! DM_MIN(0_int32)
        logical(c_bool), intent(out) :: all_bins_have_min_count
            !! `.true.` if every bin within each reference point's own `n_bins_per_point`, at every
            !! reference point, reaches at least `min_count`

        integer(int32) :: i_point, i_bin
        logical(c_bool) :: all_pass

        all_pass = .true.
        do concurrent(i_bin=1:n_bins, i_point=1:n_points) shared(mean_pmf_counts, min_count, n_bins_per_point) &
                reduce(.and.:all_pass)
            if (i_bin <= n_bins_per_point(i_point)) then
                all_pass = all_pass .and. (mean_pmf_counts(i_point, i_bin) >= min_count)
            end if
        end do
        all_bins_have_min_count = logical(all_pass, kind=c_bool)
    end subroutine check_mean_pmf_min_counts_impl

    !> Calculates the fractional overlap between two intervals: the fraction of the overlap over
    !| the total range of the first interval, so for `a_min < b_min`,
    !| \[ \frac{\min(\texttt{a\_max}, \texttt{b\_max}) - \texttt{b\_min}}{\texttt{a\_max} - \texttt{a\_min}} \]
    !| For continuous intervals; inclusive integer index spans use `compute_index_span_overlap`.
    !| Private helper of
    !| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]].
    pure real(real64) function compute_fractional_overlap(a_min, a_max, b_min, b_max) result(overlap_percent)
        real(real64), intent(in) :: a_min
            !! Lower bound of the first interval
        real(real64), intent(in) :: a_max
            !! Upper bound of the first interval, assumed >= a_min
        real(real64), intent(in) :: b_min
            !! Lower bound of the second interval
        real(real64), intent(in) :: b_max
            !! Upper bound of the second interval, assumed >= b_min

        real(real64) :: left_max, right_min

        if (a_max == a_min) then
            ! A zero-width first interval: maximum overlap if it sits inside the second interval
            if (b_min <= a_max .and. b_max >= a_max) then
                overlap_percent = 1.0_real64
            else
                overlap_percent = 0.0_real64
            end if
        else
            left_max = min(a_max, b_max)
            right_min = max(a_min, b_min)
            ! Negative exactly when a_min < a_max < b_min < b_max, i.e. no overlap at all
            overlap_percent = max(0.0_real64, (left_max - right_min)/(a_max - a_min))
        end if
    end function compute_fractional_overlap

    !> summary: Test one candidate pair's bootstrapped confidence intervals against the current best, and detect a JSD plateau
    !| AUTHOR_LASZLO_LANG
    !| For each study, the overlap of the new candidate's interval with the best one so far is
    !| measured as a fraction of the new interval's width, and counted when it reaches
    !| `succeeding_ci_overlap`. A plateau is found when fewer studies pass than for the current
    !| best (the best is then kept), or, after the new candidate became the best, when the count
    !| meets the condition `join_method` names.
    pure subroutine check_plateau_condition_impl(confidence_interval, best_candidate_pair_confidence_interval, n_studies, &
                                                 best_candidate_index, best_exceeded_ci_overlap_count, candidate_index, &
                                                 join_method, succeeding_ci_overlap, ci_overlap, plateau_found)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), intent(in) :: confidence_interval(2, n_studies)
            !! Bootstrapped JSD confidence interval `[lower, upper]` of the candidate under test
        real(real64), intent(inout) :: best_candidate_pair_confidence_interval(2, n_studies)
            !! JSD confidence intervals for the current best candidate pair; overwritten with
            !! `confidence_interval` unless the new candidate is worse
        integer(int32), intent(inout) :: best_candidate_index
            !! Candidate-grid index of the current best candidate pair; overwritten with
            !! `candidate_index` unless the new candidate is worse
        integer(int32), intent(inout) :: best_exceeded_ci_overlap_count
            !! Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
            !! best candidate pair; overwritten unless the new candidate is worse
            !! DM_MIN(0_int32)
        integer(int32), intent(in) :: candidate_index
            !! Candidate-grid index of the candidate pair that produced `confidence_interval`
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: join_method
            !! How the studies' overlaps are joined into the plateau condition; the median method
            !! requires `count > (n_studies - 1) / 2`
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        real(real64), intent(in) :: succeeding_ci_overlap
            !! Minimum fractional overlap with the best candidate's interval for a study to pass
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), intent(out) :: ci_overlap(n_studies)
            !! Per-study overlap, in `[0, 1]`, of `confidence_interval` with the best interval as
            !! it was on entry, relative to the new interval's width; always filled. Meaningless
            !! for the very first candidate, whose best interval is only a seed
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` once the new candidate is no better than the previous best, or once
            !! `join_method`'s overlap condition is met by the new candidate

        integer(int32) :: exceeds_min_ci_overlap, i_study
        real(real64) :: overlap
        logical(c_bool) :: plateau

        exceeds_min_ci_overlap = 0_int32
        do concurrent(i_study=1:n_studies) local(overlap) shared(confidence_interval, &
                                                                  best_candidate_pair_confidence_interval, &
                                                                  succeeding_ci_overlap, ci_overlap) &
                                                           reduce(+:exceeds_min_ci_overlap)
            ! The new candidate goes first: its own width is the denominator, so an interval inside
            ! the best one scores 1.0 and a wider one less.
            overlap = compute_fractional_overlap(confidence_interval(1, i_study), confidence_interval(2, i_study), &
                                                 best_candidate_pair_confidence_interval(1, i_study), &
                                                 best_candidate_pair_confidence_interval(2, i_study))
            ci_overlap(i_study) = overlap
            if (overlap >= succeeding_ci_overlap) exceeds_min_ci_overlap = exceeds_min_ci_overlap + 1_int32
        end do

        if (exceeds_min_ci_overlap < best_exceeded_ci_overlap_count) then
            plateau = .true.
        else
            best_candidate_index = candidate_index
            best_exceeded_ci_overlap_count = exceeds_min_ci_overlap
            best_candidate_pair_confidence_interval = confidence_interval

            select case (join_method)
            case (METHOD_JOIN_MIN)
                plateau = exceeds_min_ci_overlap == n_studies
            case (METHOD_JOIN_MAX)
                plateau = exceeds_min_ci_overlap > 0_int32
            case default ! METHOD_JOIN_MEDIAN
                plateau = exceeds_min_ci_overlap >= (n_studies - 1_int32)/2_int32 + 1_int32
            end select
        end if
        plateau_found = logical(plateau, kind=c_bool)
    end subroutine check_plateau_condition_impl

    !> summary: Test one candidate's per-study JSD against the previous admissible candidate for a relative-effect-size plateau
    !| AUTHOR_LASZLO_LANG
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
    pure subroutine check_effect_size_plateau_condition_impl(global_js_divergence, prev_global_js_divergence, n_studies, &
                                                             has_previous, delta_median_threshold, delta_max_threshold, &
                                                             delta_epsilon, delta_min_consecutive_transitions, &
                                                             n_consecutive_ok, delta, delta_median, delta_max, &
                                                             plateau_found, tmp_delta_perm)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), intent(in) :: global_js_divergence(n_studies)
            !! Current admissible candidate's observed global JSD per study
            !! DM_MIN(0.0_real64)
        real(real64), intent(in) :: prev_global_js_divergence(n_studies)
            !! Previous admissible candidate's observed global JSD per study; ignored when
            !! `has_previous` is `.false.`
            !! DM_MIN(0.0_real64)
        logical(c_bool), intent(in) :: has_previous
            !! `.false.` for the first admissible candidate, which has no transition
        real(real64), intent(in) :: delta_median_threshold
            !! Upper bound on the median relative change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
        real(real64), intent(in) :: delta_max_threshold
            !! Upper bound on the largest relative change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
        real(real64), intent(in) :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! DM_MIN(above(0.0_real64))
        integer(int32), intent(in) :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! DM_MIN(1_int32)
        integer(int32), intent(inout) :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; reset to zero when this one does
            !! not qualify or `has_previous` is `.false.`
            !! DM_MIN(0_int32)
        real(real64), intent(out) :: delta(n_studies)
            !! Per-study relative JSD change from the previous admissible candidate; `-1.0_real64`
            !! throughout iff `.not. has_previous`
        real(real64), intent(out) :: delta_median
            !! Median of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        real(real64), intent(out) :: delta_max
            !! Maximum of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`
        integer(int32), intent(out) :: tmp_delta_perm(n_studies)
            !! Working array: sorting permutation for `delta`, used to compute `delta_median`

        integer(int32) :: i_study

        if (.not. has_previous) then
            delta = -1.0_real64
            delta_median = -1.0_real64
            delta_max = -1.0_real64
            n_consecutive_ok = 0_int32
            plateau_found = logical(.false., kind=c_bool)
            return
        end if

        do concurrent(i_study=1:n_studies) shared(delta, global_js_divergence, prev_global_js_divergence, delta_epsilon)
            delta(i_study) = abs(global_js_divergence(i_study) - prev_global_js_divergence(i_study)) &
                            /max(prev_global_js_divergence(i_study), delta_epsilon)
        end do

        call init_perm(tmp_delta_perm)
        call sort_array_heapsort(delta, tmp_delta_perm)
        call calc_quantile_impl(delta, n_studies, tmp_delta_perm, 0.5_real64, delta_median)
        delta_max = maxval(delta)

        if (delta_median < delta_median_threshold .and. delta_max < delta_max_threshold) then
            n_consecutive_ok = n_consecutive_ok + 1_int32
        else
            n_consecutive_ok = 0_int32
        end if

        plateau_found = logical(n_consecutive_ok >= delta_min_consecutive_transitions, kind=c_bool)
    end subroutine check_effect_size_plateau_condition_impl

    !> The second admissibility gate and the observed JSD of one parameter-search candidate,
    !| shared by both searches. First, every study must include at least one residual at every
    !| reference point (`has_empty_study` otherwise, and nothing else is computed): the
    !| equal-weight consensus divides by `n_studies`, so an empty study would leave it
    !| unnormalised. Then the consensus pmf
    !| ([[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]]) must
    !| reach `min_count` in every bin, scoped to each point's own bin count
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| Only if both pass are each study's per-point JSD, weights and observed global JSD computed,
    !| and the confidence interval seeded with the observed value. Internal helper, not published;
    !| public because both parameter searches call it. Every array dummy is a search buffer viewed
    !| packed at the candidate's sizes.
    pure subroutine score_candidate_histograms(n_points, max_n_bins, n_studies, n_bins_per_point, min_count, pmfs, counts, &
                                               included_n_reps, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps, &
                                               js_divergences, weights, global_js_divergence, confidence_interval, &
                                               has_empty_study, all_bins_have_min_count)
        integer(int32), intent(in) :: n_points
            !! The candidate's number of reference points
        integer(int32), intent(in) :: max_n_bins
            !! The candidate's widest per-point bin count
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! Bin count per reference point
        integer(int32), intent(in) :: min_count
            !! Minimum count every consensus bin must reach
        real(real64), intent(in) :: pmfs(n_points, max_n_bins, n_studies)
            !! Every study's pmfs
        integer(int32), intent(in) :: counts(n_points, max_n_bins, n_studies)
            !! Every study's histogram counts
        integer(int32), intent(in) :: included_n_reps(n_points, n_studies)
            !! Every study's included-replicate counts
        real(real64), intent(out) :: mean_pmf(n_points, max_n_bins)
            !! The consensus pmf
        integer(int32), intent(out) :: mean_pmf_counts(n_points, max_n_bins)
            !! The consensus histogram counts
        integer(int32), intent(out) :: mean_pmf_included_n_reps(n_points)
            !! The consensus included-replicate counts
        real(real64), intent(out) :: js_divergences(n_points, n_studies)
            !! Each study's per-point JSD against the consensus
        real(real64), intent(out) :: weights(n_points, n_studies)
            !! Each study's per-point weights
        real(real64), intent(out) :: global_js_divergence(n_studies)
            !! Each study's observed weighted global JSD
        real(real64), intent(out) :: confidence_interval(2, n_studies)
            !! `[J, J]` per study, the observed global JSD, to seed the bootstrap
        logical(c_bool), intent(out) :: has_empty_study
            !! `.true.` iff some study has no included residual at some reference point
        logical(c_bool), intent(out) :: all_bins_have_min_count
            !! `.true.` iff every consensus bin reaches `min_count`; `.false.` when `has_empty_study`

        integer(int32) :: i_study

        has_empty_study = logical(any(included_n_reps == 0_int32), kind=c_bool)
        all_bins_have_min_count = logical(.false., kind=c_bool)
        if (has_empty_study) return

        call create_mean_pmf_impl(pmfs, counts, max_n_bins, n_points, n_studies, included_n_reps, mean_pmf, &
                                  mean_pmf_included_n_reps, mean_pmf_counts)
        call check_mean_pmf_min_counts_impl(mean_pmf_counts, max_n_bins, n_bins_per_point, n_points, min_count, &
                                            all_bins_have_min_count)
        if (.not. all_bins_have_min_count) return

        do i_study = 1, n_studies
            call compute_divergence_per_reference_point_impl(pmfs(:, :, i_study), mean_pmf, n_points, max_n_bins, &
                                                             js_divergences(:, i_study))
            call compute_weighted_global_divergence_impl(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), &
                                                          mean_pmf_included_n_reps, global_js_divergence(i_study), &
                                                          weights(:, i_study))
            confidence_interval(1, i_study) = global_js_divergence(i_study)
            confidence_interval(2, i_study) = global_js_divergence(i_study)
        end do
    end subroutine score_candidate_histograms

end module tox_data_integration_js_comp_test_search_impl
