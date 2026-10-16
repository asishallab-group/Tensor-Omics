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
module tox_data_integration_js_comp_test_search_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use f42_stats_impl, only: calc_median
    use tox_errors, only: set_ok, set_err_once, is_err, get_err_code
    use tox_data_integration_jsd_impl, only: compute_consensus_divergences, create_mean_pmf_impl
    use tox_data_integration_stats_impl, only: bootstrap_histogram_impl
    M_IMPLICIT_NONE
    private
    public :: check_neighborhood_overlaps_impl, check_mean_pmf_min_counts_impl, check_plateau_condition_impl, check_effect_size_plateau_condition_impl, &
             score_candidate_histograms, evaluate_admissible_candidate, select_final_candidate, &
             METHOD_JOIN_MIN, METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE, &
             MODE_PLATEAU_BOTH, MAX_CANDIDATE_PAIRS

    ! `join_method`'s mode table; the `METHOD_` prefix follows from the argument name `join_method`.
    ! Plain parameters: no DM_DEFAULT and no other module needs the numbers themselves.
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
        !! Most candidates either search evaluates. The fixed-k search derives its
        !! [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MAX_POINT_CANDIDATES(variable)]]
        !! from it as `MAX_CANDIDATE_PAIRS / size(KX_FACTORS)` (16 / 2 = 8,
        !! [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):KX_FACTORS(variable)]]);
        !! the adaptive search caps its candidate count at the same value

    ! `plateau_mode`'s mode table, from the M_MODE_PLATEAU_* of src/tox/data_integration/macros.h:
    ! M_MODE_PLATEAU_CI_OVERLAP is also both searches' DM_DEFAULT, which cannot name a Fortran
    ! parameter, and the family is kept together.
    integer(int32), parameter :: MODE_PLATEAU_CI_OVERLAP = M_MODE_PLATEAU_CI_OVERLAP
        !! Plateau mode: CI overlap only
    integer(int32), parameter :: MODE_PLATEAU_EFFECT_SIZE = M_MODE_PLATEAU_EFFECT_SIZE
        !! Plateau mode: relative-effect-size stability only
    integer(int32), parameter :: MODE_PLATEAU_BOTH = M_MODE_PLATEAU_BOTH
        !! Plateau mode: either CI overlap or relative-effect-size stability

contains

    !> summary: Test whether every pair of consecutive neighborhoods overlaps by at least a minimum fraction
    !| AUTHOR_LASZLO_LANG
    !| The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
    !| neighborhood spans of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
    !| number of positions they share divided by the length of the earlier span; a reversed span
    !| (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
    !| `test_`-prefixed export as a test case.
    pure subroutine check_neighborhood_overlaps_impl(neighborhood_range, n_points, min_neighbor_overlap, &
                                                      all_have_min_neighbor_overlap)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: neighborhood_range(2, n_points)
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span
            !! DM_MIN(1_int32)
        real(real64), intent(in) :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        logical(c_bool), intent(out) :: all_have_min_neighbor_overlap
            !! `.true.` if every pair of consecutive neighborhoods overlaps by at least
            !! `min_neighbor_overlap`

        integer(int32) :: i_point
        real(real64) :: overlap
        logical(c_bool) :: all_pass

        all_pass = .true.
        do concurrent(i_point=1:n_points - 1) local(overlap) shared(neighborhood_range, min_neighbor_overlap) &
                reduce(.and.:all_pass)
            overlap = compute_index_span_overlap(neighborhood_range(1, i_point), neighborhood_range(2, i_point), &
                                                 neighborhood_range(1, i_point + 1), neighborhood_range(2, i_point + 1))
            all_pass = all_pass .and. (overlap >= min_neighbor_overlap)
        end do
        all_have_min_neighbor_overlap = logical(all_pass, kind=c_bool)
    end subroutine check_neighborhood_overlaps_impl

    !> Counts the positions two inclusive integer runs `[a_min, a_max]` and `[b_min, b_max]` share,
    !| relative to the number of positions in the first run:
    !| \[ \frac{\max(0, \min(\texttt{a\_max}, \texttt{b\_max}) - \max(\texttt{a\_min}, \texttt{b\_min}) + 1)}{\texttt{a\_max} - \texttt{a\_min} + 1} \]
    !| For example `[5,6]` and `[6,7]` give 0.5. A reversed run (`max < min`) is empty, so the
    !| result is 0 rather than a division by zero. Private helper of
    !| [[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl(interface)]].
    pure real(real64) function compute_index_span_overlap(a_min, a_max, b_min, b_max) result(overlap_fraction)
        integer(int32), intent(in) :: a_min
            !! First position of the first run
        integer(int32), intent(in) :: a_max
            !! Last position of the first run; the run is empty if `a_max < a_min`
        integer(int32), intent(in) :: b_min
            !! First position of the second run
        integer(int32), intent(in) :: b_max
            !! Last position of the second run; the run is empty if `b_max < b_min`

        integer(int32) :: n_shared

        ! An empty first run shares no position
        if (a_max < a_min) then
            overlap_fraction = 0.0_real64
            return
        end if

        n_shared = max(0_int32, min(a_max, b_max) - max(a_min, b_min) + 1_int32)
        overlap_fraction = real(n_shared, real64)/real(a_max - a_min + 1_int32, real64)
    end function compute_index_span_overlap

    !> summary: Test whether every bin of a mean pmf reaches a minimum absolute count
    !| AUTHOR_LASZLO_LANG
    !| The second admissibility gate of a parameter-search candidate, checked after
    !| [[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl(interface)]].
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

        call calc_median(delta, n_studies, tmp_delta_perm, delta_median)
        delta_max = maxval(delta)

        if (delta_median < delta_median_threshold .and. delta_max < delta_max_threshold) then
            n_consecutive_ok = n_consecutive_ok + 1_int32
        else
            n_consecutive_ok = 0_int32
        end if

        plateau_found = logical(n_consecutive_ok >= delta_min_consecutive_transitions, kind=c_bool)
    end subroutine check_effect_size_plateau_condition_impl

    !> summary: Apply the second admissibility gate to one parameter-search candidate and compute its observed JSD
    !| AUTHOR_LASZLO_LANG
    !| The second admissibility gate and the observed JSD of one parameter-search candidate,
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

        call compute_consensus_divergences(pmfs, mean_pmf, n_points, max_n_bins, max_n_bins, n_studies, included_n_reps, &
                                           mean_pmf_included_n_reps, js_divergences, global_js_divergence, weights)
        do i_study = 1, n_studies
            confidence_interval(1, i_study) = global_js_divergence(i_study)
            confidence_interval(2, i_study) = global_js_divergence(i_study)
        end do
    end subroutine score_candidate_histograms

    !> summary: Evaluate one admissible parameter-search candidate: traces, bootstrap and both plateau criteria
    !| AUTHOR_LASZLO_LANG
    !| The part of both parameter searches after a candidate passed both admissibility gates
    !| ([[tox_data_integration_js_comp_test_search_impl(module):score_candidate_histograms(subroutine)]])
    !| and the search counted it in `n_admissible_evaluated` and traced its identity. Fills trace
    !| column `n_admissible_evaluated`, tests the effect-size criterion
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]])
    !| before the bootstrap
    !| ([[tox_data_integration_stats_impl(module):bootstrap_histogram_impl(interface)]]) reuses
    !| the observed values as scratch, then the CI-overlap criterion
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]),
    !| and updates the running state: the best and the smallest-uncertainty candidate with their
    !| per-point snapshots, the previous candidate's JSD and the effect-size run. `plateau_found`
    !| is the criterion `plateau_mode` selects; an effect-size plateau it selects makes this
    !| candidate the best. Rows `n_points + 1` onwards of the per-point trace column are set to 0
    !| (`.false.` for the failure flags).
    !|
    !| Internal helper, not published; public because both parameter searches call it, and impure
    !| because the bootstrap is. Every candidate-sized dummy is a search buffer viewed packed at
    !| the candidate's sizes. A failed bootstrap returns its error code in `ierr` (argument
    !| position cleared) at once.
    subroutine evaluate_admissible_candidate(n_points, max_n_bins, n_studies, max_n_points_candidate, candidate_index, &
                                             n_admissible_evaluated, n_bootstraps, n_bootstrapping_top_k_jsds, &
                                             join_method, plateau_mode, succeeding_ci_overlap, delta_median_threshold, &
                                             delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, &
                                             n_bins_per_point, occupancy_failed, n_pooled_residuals, min_bin_occupancy, &
                                             mean_bin_occupancy, max_bin_occupancy, sturges_bins, fd_bins, &
                                             shared_residual_range_low, shared_residual_range_high, included_n_reps, &
                                             mean_pmf_counts, mean_pmf_included_n_reps, counts, pmfs, mean_pmf, &
                                             js_divergences, weights, global_js_divergence, confidence_interval, &
                                             trace_global_js_divergence, trace_ci_lower, trace_ci_upper, trace_ci_width, &
                                             trace_ci_width_relative, trace_ci_overlap, trace_delta, trace_delta_median, &
                                             trace_delta_max, trace_selected_n_bins, trace_occupancy_failed, &
                                             trace_n_pooled_residuals, trace_min_bin_occupancy, trace_mean_bin_occupancy, &
                                             trace_max_bin_occupancy, trace_sturges_bins, trace_fd_bins, &
                                             trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                             prev_global_js_divergence, n_consecutive_effect_size_ok, &
                                             best_candidate_confidence_interval, best_candidate_index, &
                                             best_exceeded_ci_overlap_count, best_n_bins_per_point, &
                                             best_shared_residual_range_low, best_shared_residual_range_high, &
                                             best_uncertainty_value, best_uncertainty_candidate_index, &
                                             best_uncertainty_confidence_interval, best_uncertainty_n_bins_per_point, &
                                             best_uncertainty_shared_residual_range_low, &
                                             best_uncertainty_shared_residual_range_high, plateau_found, &
                                             tmp_bootstrapping_top_k_jsds, tmp_point_pool, tmp_point_draw, tmp_delta_perm, &
                                             ierr, two_sided_bootstrapping_significance_level, random_seed)
        integer(int32), intent(in) :: n_points
            !! The candidate's number of reference points
        integer(int32), intent(in) :: max_n_bins
            !! The candidate's widest per-point bin count
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_points_candidate
            !! First extent of the per-point traces
        integer(int32), intent(in) :: candidate_index
            !! The candidate's index in the search's candidate sequence
        integer(int32), intent(in) :: n_admissible_evaluated
            !! Number of admissible candidates so far, this one included: its trace column
        integer(int32), intent(in) :: n_bootstraps
            !! Forwarded to the bootstrap
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Forwarded to the bootstrap (heap size)
        integer(int32), intent(in) :: join_method
            !! Forwarded to check_plateau_condition_impl
        integer(int32), intent(in) :: plateau_mode
            !! Which plateau criterion sets `plateau_found`, one of the `MODE_PLATEAU_*` values
        real(real64), intent(in) :: succeeding_ci_overlap
            !! Forwarded to check_plateau_condition_impl
        real(real64), intent(in) :: delta_median_threshold
            !! Forwarded to check_effect_size_plateau_condition_impl
        real(real64), intent(in) :: delta_max_threshold
            !! Forwarded to check_effect_size_plateau_condition_impl
        real(real64), intent(in) :: delta_epsilon
            !! Forwarded to check_effect_size_plateau_condition_impl; also the floor of the relative CI width
        integer(int32), intent(in) :: delta_min_consecutive_transitions
            !! Forwarded to check_effect_size_plateau_condition_impl
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! The candidate's bin count per reference point
        logical(c_bool), intent(in) :: occupancy_failed(n_points)
            !! The candidate's bin-count failure flag per reference point
        integer(int32), intent(in) :: n_pooled_residuals(n_points)
            !! The candidate's pooled residual count per reference point
        integer(int32), intent(in) :: min_bin_occupancy(n_points)
            !! The candidate's minimum bin occupancy per reference point
        real(real64), intent(in) :: mean_bin_occupancy(n_points)
            !! The candidate's mean bin occupancy per reference point
        integer(int32), intent(in) :: max_bin_occupancy(n_points)
            !! The candidate's maximum bin occupancy per reference point
        integer(int32), intent(in) :: sturges_bins(n_points)
            !! The candidate's Sturges' rule diagnostic per reference point
        integer(int32), intent(in) :: fd_bins(n_points)
            !! The candidate's Freedman-Diaconis rule diagnostic per reference point
        real(real64), intent(in) :: shared_residual_range_low(n_points)
            !! The candidate's R_low per reference point
        real(real64), intent(in) :: shared_residual_range_high(n_points)
            !! The candidate's R_high per reference point
        integer(int32), intent(in) :: included_n_reps(n_points, n_studies)
            !! Every study's included-replicate counts
        integer(int32), intent(in) :: mean_pmf_counts(n_points, max_n_bins)
            !! The consensus histogram counts, the bootstrap's pool
        integer(int32), intent(in) :: mean_pmf_included_n_reps(n_points)
            !! The consensus included-replicate counts
        integer(int32), intent(out) :: counts(n_points, max_n_bins, n_studies)
            !! Every study's histogram counts; bootstrap scratch here
        real(real64), intent(out) :: pmfs(n_points, max_n_bins, n_studies)
            !! Every study's pmfs; bootstrap scratch here
        real(real64), intent(out) :: mean_pmf(n_points, max_n_bins)
            !! The consensus pmf; bootstrap scratch here
        real(real64), intent(out) :: js_divergences(n_points, n_studies)
            !! Each study's per-point JSD; bootstrap scratch here
        real(real64), intent(out) :: weights(n_points, n_studies)
            !! Each study's per-point weights; bootstrap scratch here
        real(real64), intent(inout) :: global_js_divergence(n_studies)
            !! Each study's observed global JSD on entry; bootstrap scratch afterwards
        real(real64), intent(inout) :: confidence_interval(2, n_studies)
            !! The `[J, J]` seed on entry, the bootstrapped interval on return
        real(real64), intent(inout) :: trace_global_js_divergence(n_studies, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the observed global JSD
        real(real64), intent(inout) :: trace_ci_lower(n_studies, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the confidence-interval lower bound
        real(real64), intent(inout) :: trace_ci_upper(n_studies, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the confidence-interval upper bound
        real(real64), intent(inout) :: trace_ci_width(n_studies, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the confidence-interval width
        real(real64), intent(inout) :: trace_ci_width_relative(n_studies, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the relative confidence-interval width
        real(real64), intent(inout) :: trace_ci_overlap(n_studies, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the overlap with the best interval so far
        real(real64), intent(inout) :: trace_delta(n_studies, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the relative JSD change
        real(real64), intent(inout) :: trace_delta_median(MAX_CANDIDATE_PAIRS)
            !! The search's trace of the median relative JSD change
        real(real64), intent(inout) :: trace_delta_max(MAX_CANDIDATE_PAIRS)
            !! The search's trace of the largest relative JSD change
        integer(int32), intent(inout) :: trace_selected_n_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the per-point bin counts
        logical(c_bool), intent(inout) :: trace_occupancy_failed(max_n_points_candidate, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the per-point bin-count failure flags
        integer(int32), intent(inout) :: trace_n_pooled_residuals(max_n_points_candidate, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the per-point pooled residual counts
        integer(int32), intent(inout) :: trace_min_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the per-point minimum bin occupancy
        real(real64), intent(inout) :: trace_mean_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the per-point mean bin occupancy
        integer(int32), intent(inout) :: trace_max_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the per-point maximum bin occupancy
        integer(int32), intent(inout) :: trace_sturges_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the per-point Sturges' rule diagnostic
        integer(int32), intent(inout) :: trace_fd_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the per-point Freedman-Diaconis rule diagnostic
        real(real64), intent(inout) :: trace_shared_residual_range_low(max_n_points_candidate, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the per-point R_low
        real(real64), intent(inout) :: trace_shared_residual_range_high(max_n_points_candidate, MAX_CANDIDATE_PAIRS)
            !! The search's trace of the per-point R_high
        real(real64), intent(inout) :: prev_global_js_divergence(n_studies)
            !! The previous admissible candidate's observed global JSD; this candidate's on return
        integer(int32), intent(inout) :: n_consecutive_effect_size_ok
            !! Running count of consecutive qualifying effect-size transitions
        real(real64), intent(inout) :: best_candidate_confidence_interval(2, n_studies)
            !! The best candidate's confidence interval
        integer(int32), intent(inout) :: best_candidate_index
            !! The best candidate's index
        integer(int32), intent(inout) :: best_exceeded_ci_overlap_count
            !! Forwarded to check_plateau_condition_impl
        integer(int32), intent(inout) :: best_n_bins_per_point(n_points)
            !! The best candidate's bin counts, overwritten when this candidate becomes the best
        real(real64), intent(inout) :: best_shared_residual_range_low(n_points)
            !! The best candidate's R_low, like `best_n_bins_per_point`
        real(real64), intent(inout) :: best_shared_residual_range_high(n_points)
            !! The best candidate's R_high, like `best_n_bins_per_point`
        real(real64), intent(inout) :: best_uncertainty_value
            !! Smallest median confidence-interval width so far
        integer(int32), intent(inout) :: best_uncertainty_candidate_index
            !! Index of the candidate with that width
        real(real64), intent(inout) :: best_uncertainty_confidence_interval(2, n_studies)
            !! That candidate's confidence interval
        integer(int32), intent(inout) :: best_uncertainty_n_bins_per_point(n_points)
            !! That candidate's bin counts, overwritten when this candidate takes its place
        real(real64), intent(inout) :: best_uncertainty_shared_residual_range_low(n_points)
            !! That candidate's R_low, like `best_uncertainty_n_bins_per_point`
        real(real64), intent(inout) :: best_uncertainty_shared_residual_range_high(n_points)
            !! That candidate's R_high, like `best_uncertainty_n_bins_per_point`
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` iff the criterion `plateau_mode` selects found a plateau
        real(real64), intent(out) :: tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies)
            !! Working array: the bootstrap's bottom-k/top-k heaps
        integer(int32), intent(out) :: tmp_point_pool(M_MAX_N_BINS)
            !! Working array of the bootstrap: one point's pool
        integer(int32), intent(out) :: tmp_point_draw(M_MAX_N_BINS)
            !! Working array of the bootstrap: one study's draw at one point
        integer(int32), intent(out) :: tmp_delta_perm(n_studies)
            !! Working array: sorting permutation for the per-study medians
        integer(int32), intent(out) :: ierr
            !! Error code of the bootstrap if it fails (see above)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Forwarded to the bootstrap
        integer(int32), intent(in), optional :: random_seed
            !! Forwarded to the bootstrap

        integer(int32) :: bootstrap_ierr
        real(real64) :: median_ci_width
        logical(c_bool) :: ci_plateau_found, effect_size_plateau_found

        call set_ok(ierr)

        ! Trace the candidate and test the effect-size criterion, whatever plateau_mode is. Must
        ! come before the bootstrap, which reuses global_js_divergence as scratch.
        trace_global_js_divergence(1:n_studies, n_admissible_evaluated) = global_js_divergence
        trace_selected_n_bins(1:n_points, n_admissible_evaluated) = n_bins_per_point
        trace_occupancy_failed(1:n_points, n_admissible_evaluated) = occupancy_failed
        trace_n_pooled_residuals(1:n_points, n_admissible_evaluated) = n_pooled_residuals
        trace_min_bin_occupancy(1:n_points, n_admissible_evaluated) = min_bin_occupancy
        trace_mean_bin_occupancy(1:n_points, n_admissible_evaluated) = mean_bin_occupancy
        trace_max_bin_occupancy(1:n_points, n_admissible_evaluated) = max_bin_occupancy
        trace_sturges_bins(1:n_points, n_admissible_evaluated) = sturges_bins
        trace_fd_bins(1:n_points, n_admissible_evaluated) = fd_bins
        trace_shared_residual_range_low(1:n_points, n_admissible_evaluated) = shared_residual_range_low
        trace_shared_residual_range_high(1:n_points, n_admissible_evaluated) = shared_residual_range_high
        ! The rows beyond this candidate's points are padding: 0, as the final runs' unused bins.
        trace_selected_n_bins(n_points + 1:max_n_points_candidate, n_admissible_evaluated) = 0_int32
        trace_occupancy_failed(n_points + 1:max_n_points_candidate, n_admissible_evaluated) = .false.
        trace_n_pooled_residuals(n_points + 1:max_n_points_candidate, n_admissible_evaluated) = 0_int32
        trace_min_bin_occupancy(n_points + 1:max_n_points_candidate, n_admissible_evaluated) = 0_int32
        trace_mean_bin_occupancy(n_points + 1:max_n_points_candidate, n_admissible_evaluated) = 0.0_real64
        trace_max_bin_occupancy(n_points + 1:max_n_points_candidate, n_admissible_evaluated) = 0_int32
        trace_sturges_bins(n_points + 1:max_n_points_candidate, n_admissible_evaluated) = 0_int32
        trace_fd_bins(n_points + 1:max_n_points_candidate, n_admissible_evaluated) = 0_int32
        trace_shared_residual_range_low(n_points + 1:max_n_points_candidate, n_admissible_evaluated) = 0.0_real64
        trace_shared_residual_range_high(n_points + 1:max_n_points_candidate, n_admissible_evaluated) = 0.0_real64

        ! The first admissible candidate has no previous one.
        call check_effect_size_plateau_condition_impl(global_js_divergence, prev_global_js_divergence, &
                                                      n_studies, logical(n_admissible_evaluated > 1_int32, kind=c_bool), &
                                                      delta_median_threshold, delta_max_threshold, delta_epsilon, &
                                                      delta_min_consecutive_transitions, &
                                                      n_consecutive_effect_size_ok, &
                                                      trace_delta(1:n_studies, n_admissible_evaluated), &
                                                      trace_delta_median(n_admissible_evaluated), &
                                                      trace_delta_max(n_admissible_evaluated), &
                                                      effect_size_plateau_found, tmp_delta_perm)
        prev_global_js_divergence = global_js_divergence

        ! Same n_bootstraps/random_seed for every candidate, for comparability. The observed-value
        ! arrays above double as the bootstrap's work arrays, passed whole and viewed packed.
        call bootstrap_histogram_impl(n_bootstraps, max_n_bins, n_points, n_studies, mean_pmf_counts, &
                                      mean_pmf_included_n_reps, included_n_reps, &
                                      n_bootstrapping_top_k_jsds, confidence_interval, &
                                      tmp_bootstrapping_top_k_jsds, counts, pmfs, mean_pmf, &
                                      js_divergences, weights, global_js_divergence, tmp_point_pool, &
                                      tmp_point_draw, two_sided_bootstrapping_significance_level, random_seed, &
                                      bootstrap_ierr)
        ! Defensive: the counts add up by construction, so only a failed rng allocation can
        ! land here (untested). Without the return the interval would stay at its zero-width
        ! [J, J] seed and could win the smallest-uncertainty ranking.
        if (is_err(bootstrap_ierr)) then
            call set_err_once(ierr, get_err_code(bootstrap_ierr))
            return
        end if

        trace_ci_lower(1:n_studies, n_admissible_evaluated) = confidence_interval(1, 1:n_studies)
        trace_ci_upper(1:n_studies, n_admissible_evaluated) = confidence_interval(2, 1:n_studies)
        trace_ci_width(1:n_studies, n_admissible_evaluated) = &
            trace_ci_upper(1:n_studies, n_admissible_evaluated) - trace_ci_lower(1:n_studies, n_admissible_evaluated)
        trace_ci_width_relative(1:n_studies, n_admissible_evaluated) = &
            trace_ci_width(1:n_studies, n_admissible_evaluated) &
            / max(trace_global_js_divergence(1:n_studies, n_admissible_evaluated), delta_epsilon)

        ! Track the admissible candidate with the smallest median CI width (the no-plateau
        ! fallback). calc_median re-sorts tmp_delta_perm for trace_ci_width.
        call calc_median(trace_ci_width(1:n_studies, n_admissible_evaluated), n_studies, tmp_delta_perm, &
                         median_ci_width)
        if (median_ci_width < best_uncertainty_value) then
            best_uncertainty_value = median_ci_width
            best_uncertainty_candidate_index = candidate_index
            best_uncertainty_confidence_interval = confidence_interval
            best_uncertainty_n_bins_per_point = n_bins_per_point
            best_uncertainty_shared_residual_range_low = shared_residual_range_low
            best_uncertainty_shared_residual_range_high = shared_residual_range_high
        end if

        call check_plateau_condition_impl(confidence_interval, best_candidate_confidence_interval, n_studies, &
                                          best_candidate_index, best_exceeded_ci_overlap_count, candidate_index, join_method, &
                                          succeeding_ci_overlap, &
                                          trace_ci_overlap(1:n_studies, n_admissible_evaluated), ci_plateau_found)
        ! The first admissible candidate was compared against the -1.0 seed: record -1.0.
        if (n_admissible_evaluated == 1_int32) trace_ci_overlap(1:n_studies, 1) = -1.0_real64

        ! This candidate became the new best: snapshot its per-point values.
        if (best_candidate_index == candidate_index) then
            best_n_bins_per_point = n_bins_per_point
            best_shared_residual_range_low = shared_residual_range_low
            best_shared_residual_range_high = shared_residual_range_high
        end if

        select case (plateau_mode)
        case (MODE_PLATEAU_CI_OVERLAP)
            plateau_found = ci_plateau_found
        case (MODE_PLATEAU_EFFECT_SIZE)
            plateau_found = effect_size_plateau_found
        case default ! MODE_PLATEAU_BOTH
            plateau_found = logical(ci_plateau_found .or. effect_size_plateau_found, kind=c_bool)
        end select

        ! The best-candidate bookkeeping follows the CI-overlap criterion only; an effect-size
        ! plateau returns the candidate that triggered it.
        if (effect_size_plateau_found .and. plateau_mode /= MODE_PLATEAU_CI_OVERLAP) then
            best_candidate_index = candidate_index
            best_candidate_confidence_interval = confidence_interval
            best_n_bins_per_point = n_bins_per_point
            best_shared_residual_range_low = shared_residual_range_low
            best_shared_residual_range_high = shared_residual_range_high
        end if
    end subroutine evaluate_admissible_candidate

    !> summary: Select the candidate a parameter search returns, by its three cases
    !| AUTHOR_LASZLO_LANG
    !| The end of both parameter searches. Case 1, a plateau was found or the only candidate was
    !| admissible: the best candidate, `plateau_established = .true.`. Case 2, some candidate was
    !| admissible: the one with the smallest median confidence-interval width, with its interval.
    !| Case 3, nothing was admissible: candidate 1, every bin count `m_min`, both ranges `0.0`
    !| (for at most `max_n_points_candidate` points) and the interval `-1.0`. Cases 2 and 3 set
    !| `plateau_established = .false.`. In every case the per-point entries beyond those points
    !| are 0. The search reads the selected candidate's identity at
    !| `selected_candidate_index`. Internal helper, not published; public because both parameter
    !| searches call it.
    pure subroutine select_final_candidate(max_n_points_candidate, n_studies, n_candidates, n_admissible_evaluated, &
                                           plateau_found, m_min, best_candidate_index, best_n_points, &
                                           best_n_bins_per_point, best_shared_residual_range_low, &
                                           best_shared_residual_range_high, best_uncertainty_candidate_index, &
                                           best_uncertainty_n_points, best_uncertainty_confidence_interval, &
                                           best_uncertainty_n_bins_per_point, best_uncertainty_shared_residual_range_low, &
                                           best_uncertainty_shared_residual_range_high, first_n_points, &
                                           selected_candidate_index, n_points, n_bins_per_point, shared_residual_range_low, &
                                           shared_residual_range_high, best_candidate_confidence_interval, &
                                           plateau_established)
        integer(int32), intent(in) :: max_n_points_candidate
            !! Extent of the per-point arrays
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_candidates
            !! Number of candidates in the search's sequence, inadmissible ones included
        integer(int32), intent(in) :: n_admissible_evaluated
            !! Number of admissible candidates the search evaluated
        logical(c_bool), intent(in) :: plateau_found
            !! Whether the search stopped at a plateau
        integer(int32), intent(in) :: m_min
            !! The bin count of every point in case 3
        integer(int32), intent(in) :: best_candidate_index
            !! The best candidate's index
        integer(int32), intent(in) :: best_n_points
            !! The best candidate's number of reference points
        integer(int32), intent(in) :: best_n_bins_per_point(max_n_points_candidate)
            !! The best candidate's bin counts
        real(real64), intent(in) :: best_shared_residual_range_low(max_n_points_candidate)
            !! The best candidate's R_low
        real(real64), intent(in) :: best_shared_residual_range_high(max_n_points_candidate)
            !! The best candidate's R_high
        integer(int32), intent(in) :: best_uncertainty_candidate_index
            !! The smallest-uncertainty candidate's index
        integer(int32), intent(in) :: best_uncertainty_n_points
            !! The smallest-uncertainty candidate's number of reference points
        real(real64), intent(in) :: best_uncertainty_confidence_interval(2, n_studies)
            !! The smallest-uncertainty candidate's confidence interval
        integer(int32), intent(in) :: best_uncertainty_n_bins_per_point(max_n_points_candidate)
            !! The smallest-uncertainty candidate's bin counts
        real(real64), intent(in) :: best_uncertainty_shared_residual_range_low(max_n_points_candidate)
            !! The smallest-uncertainty candidate's R_low
        real(real64), intent(in) :: best_uncertainty_shared_residual_range_high(max_n_points_candidate)
            !! The smallest-uncertainty candidate's R_high
        integer(int32), intent(in) :: first_n_points
            !! Candidate 1's number of reference points
        integer(int32), intent(out) :: selected_candidate_index
            !! Index of the returned candidate
        integer(int32), intent(out) :: n_points
            !! The returned candidate's number of reference points
        integer(int32), intent(out) :: n_bins_per_point(max_n_points_candidate)
            !! The returned bin counts; the first `min(n_points, max_n_points_candidate)` are set,
            !! the rest are 0
        real(real64), intent(out) :: shared_residual_range_low(max_n_points_candidate)
            !! The returned R_low, like `n_bins_per_point`
        real(real64), intent(out) :: shared_residual_range_high(max_n_points_candidate)
            !! The returned R_high, like `n_bins_per_point`
        real(real64), intent(inout) :: best_candidate_confidence_interval(2, n_studies)
            !! The best candidate's interval on entry; the returned candidate's on return
        logical(c_bool), intent(out) :: plateau_established
            !! `.true.` in case 1 only

        integer(int32) :: n_filled

        ! In case 1 some candidate was admissible, so the best-candidate snapshot exists;
        ! n_candidates counts gate-failed candidates too.
        if (plateau_found .or. (n_candidates < 2_int32 .and. n_admissible_evaluated >= 1_int32)) then
            selected_candidate_index = best_candidate_index
            n_points = best_n_points
            n_bins_per_point(1:n_points) = best_n_bins_per_point(1:n_points)
            shared_residual_range_low(1:n_points) = best_shared_residual_range_low(1:n_points)
            shared_residual_range_high(1:n_points) = best_shared_residual_range_high(1:n_points)
            plateau_established = logical(.true., kind=c_bool)
        else if (n_admissible_evaluated >= 1_int32) then
            selected_candidate_index = best_uncertainty_candidate_index
            n_points = best_uncertainty_n_points
            n_bins_per_point(1:n_points) = best_uncertainty_n_bins_per_point(1:n_points)
            shared_residual_range_low(1:n_points) = best_uncertainty_shared_residual_range_low(1:n_points)
            shared_residual_range_high(1:n_points) = best_uncertainty_shared_residual_range_high(1:n_points)
            best_candidate_confidence_interval = best_uncertainty_confidence_interval
            plateau_established = logical(.false., kind=c_bool)
        else
            ! No bin count was kept for candidate 1, so m_min (as in a failed bin-count search) and
            ! a zero range (as for an all-NaN pool). An adaptive candidate 1 may have more points
            ! than the capacity, hence the clamp.
            selected_candidate_index = 1_int32
            n_points = first_n_points
            n_bins_per_point(1:min(n_points, max_n_points_candidate)) = m_min
            shared_residual_range_low(1:min(n_points, max_n_points_candidate)) = 0.0_real64
            shared_residual_range_high(1:min(n_points, max_n_points_candidate)) = 0.0_real64
            best_candidate_confidence_interval = -1.0_real64
            plateau_established = logical(.false., kind=c_bool)
        end if

        ! The entries beyond the returned points are padding: 0, as the final runs' unused bins.
        n_filled = min(n_points, max_n_points_candidate)
        n_bins_per_point(n_filled + 1:max_n_points_candidate) = 0_int32
        shared_residual_range_low(n_filled + 1:max_n_points_candidate) = 0.0_real64
        shared_residual_range_high(n_filled + 1:max_n_points_candidate) = 0.0_real64
    end subroutine select_final_candidate

end module tox_data_integration_js_comp_test_search_impl
