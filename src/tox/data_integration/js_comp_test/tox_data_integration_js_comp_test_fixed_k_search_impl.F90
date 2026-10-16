#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Fixed-k Parameter Search
!|
!| The data-driven choice of the JS-Comp-Test's fixed-k neighborhoods, `n_points` reference
!| points of `n_neighbors` genes each.
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
!| walks the `(n_points, n_neighbors)` candidate grid of
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]
!| from finest to coarsest and stops at the first admissible candidate whose bootstrapped JSDs
!| have reached a plateau, by the criteria of
!| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]].
!| `calc_js_comp_test_candidate_bounds` sizes the work arrays for a caller that allocates its own.
module tox_data_integration_js_comp_test_fixed_k_search_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use f42_math_impl, only: clamp
    use tox_errors, only: set_ok, is_err, validate_in_range_int
    use tox_data_integration_preprocessing_impl, only: construct_neighborhoods_ranged_impl, pool_means_impl, &
                                                       sort_gene_means_per_study, sort_pooled_gene_means
    use tox_data_integration_js_comp_test_binning_impl, only: determine_point_bin_count
    use tox_data_integration_js_comp_test_histograms_impl, only: build_fixed_k_histograms
    use tox_data_integration_js_comp_test_search_impl, only: check_neighborhood_overlaps_impl, score_candidate_histograms, &
                                                             evaluate_admissible_candidate, select_final_candidate, &
                                                             MAX_CANDIDATE_PAIRS
    M_IMPLICIT_NONE
    private
    public :: calc_js_comp_test_candidate_bounds, generate_js_comp_test_candidates_impl, &
             run_js_comp_test_parameter_search_impl, KX_FACTORS, MAX_POINTS, MIN_POINTS, GAMMA, MAX_POINT_CANDIDATES

    ! Tuning constants of the fixed-k candidate grid (generate_js_comp_test_candidates_impl).
    real(real64), parameter :: KX_FACTORS(2) = [0.25_real64, 0.5_real64]
        !! Neighbor-count factors, ascending: the first gives the larger `n_neighbors` candidate
    integer(int32), parameter :: MAX_POINTS = 1500_int32
    integer(int32), parameter :: MIN_POINTS = 300_int32
    real(real64), parameter :: GAMMA = 0.8_real64
    integer(int32), parameter :: MAX_POINT_CANDIDATES = MAX_CANDIDATE_PAIRS/size(KX_FACTORS, kind=int32)
        !! How often n_points_high may be multiplied by GAMMA before the grid gives up; derived so
        !! that the grid never exceeds MAX_CANDIDATE_PAIRS candidates

contains

    !> M_EXPORT_C
    !| summary: Recommend upper bounds for the js-comp-test candidate-grid work arrays
    !| AUTHOR_LASZLO_LANG
    !| For the candidate grid of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]:
    !| `max_n_points_candidate` is exactly its largest `n_points` candidate, and
    !| `max_n_neighbors_candidate` a safe (not necessarily tight) bound on its largest `n_neighbors`.
    !|
    !| Rejects a non-positive `max_n_genes_all_studies` with invalid input; both outputs are then
    !| left undefined.
    pure subroutine calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, &
                                                        max_n_neighbors_candidate, ierr)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: max_n_points_candidate
            !! Exact upper bound on the `n_points` candidate the grid ever produces
        integer(int32), intent(out) :: max_n_neighbors_candidate
            !! Safe upper bound on the `n_neighbors` candidate the grid ever produces
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure

        integer(int32) :: n_points_low

        call set_ok(ierr)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=1_int32, min=1_int32)
        if (is_err(ierr)) return

        max_n_points_candidate = clamp(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), &
                                       min_val=MIN_POINTS, max_val=MAX_POINTS)
        n_points_low = max(MIN_POINTS, ceiling(0.2_real64*real(max_n_points_candidate, real64), kind=int32))
        max_n_neighbors_candidate = max(1_int32, floor(real(max_n_genes_all_studies, real64)/ &
                                                        (KX_FACTORS(1)*real(n_points_low, real64)), kind=int32))
    end subroutine calc_js_comp_test_candidate_bounds

    !> summary: Generate the GAMMA-decay (n_points, n_neighbors) candidate grid
    !| AUTHOR_LASZLO_LANG
    !| Starting from `n_points_high = clamp(ceiling(4*sqrt(G)), MIN_POINTS, MAX_POINTS)` (`G` =
    !| `max_n_genes_all_studies`), multiplies it by GAMMA until it drops below
    !| `n_points_low = max(MIN_POINTS, ceiling(0.2*n_points_high))`, and pairs each distinct
    !| `n_points` with up to `size(KX_FACTORS)` distinct `n_neighbors = max(1, floor(G/(f*n_points_high)))`.
    !| Duplicate values from clamping or rounding collapse, so a small enough `G` yields a single
    !| candidate. Bin counts are not part of a candidate; they are chosen per reference point later.
    pure subroutine generate_js_comp_test_candidates_impl(max_n_genes_all_studies, candidates_n_points_n_neighbors, &
                                                          n_candidates)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: candidates_n_points_n_neighbors(2, M_MAX_CANDIDATE_PAIRS)
            !! Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
            !! DM_RESULT_SIZE_IS(n_candidates)
        integer(int32), intent(out) :: n_candidates
            !! Number of candidate pairs actually filled (at most M_MAX_CANDIDATE_PAIRS)

        real(real64) :: n_points_high, n_points_low
        integer(int32) :: i_point_candidate, point_candidate, prev_point_candidate
        integer(int32) :: i_neighbor_candidate, neighbor_candidate, prev_neighbor_candidate

        n_points_high = real(clamp(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), &
                                   min_val=MIN_POINTS, max_val=MAX_POINTS), kind=real64)
        n_points_low = real(max(MIN_POINTS, ceiling(0.2_real64*n_points_high, kind=int32)), kind=real64)

        prev_point_candidate = -1_int32
        n_candidates = 0_int32

        ! Sequential: each iteration dedups against the previous one.
        do i_point_candidate = 1, MAX_POINT_CANDIDATES
            if (n_points_high < n_points_low) exit

            point_candidate = nint(n_points_high, kind=int32)
            if (point_candidate /= prev_point_candidate) then
                prev_point_candidate = point_candidate
                prev_neighbor_candidate = -1_int32

                do i_neighbor_candidate = 1, size(KX_FACTORS, kind=int32)
                    neighbor_candidate = max(1_int32, floor(real(max_n_genes_all_studies, real64)/ &
                                                            (KX_FACTORS(i_neighbor_candidate)*n_points_high), kind=int32))

                    if (neighbor_candidate /= prev_neighbor_candidate) then
                        prev_neighbor_candidate = neighbor_candidate

                        n_candidates = n_candidates + 1_int32
                        candidates_n_points_n_neighbors(1, n_candidates) = point_candidate
                        candidates_n_points_n_neighbors(2, n_candidates) = neighbor_candidate
                    end if
                end do
            end if
            n_points_high = n_points_high*GAMMA
        end do
    end subroutine generate_js_comp_test_candidates_impl

    !> summary: Build one fixed-k candidate's neighborhoods for every study and check the first admissibility gate
    !| AUTHOR_LASZLO_LANG
    !| Pass A of the fixed-k parameter search for one candidate: builds every study's
    !| neighborhoods into `neighborhood_indices` and checks the first admissibility gate
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl(interface)]])
    !| study by study, stopping at the first study that fails it. Private helper; its
    !| explicit-shape dummies take the search's upper-bound-sized buffers whole and view their
    !| leading elements packed at the candidate's own sizes (see the parameter search).
    pure subroutine construct_fixed_k_candidate_neighborhoods(n_studies, max_n_genes_all_studies, n_points, n_neighbors, &
                                                              x_star, gene_means, gene_means_perms, min_neighbor_overlap, &
                                                              neighborhood_indices, neighborhood_range, &
                                                              all_have_min_neighbor_overlap)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_points
            !! The candidate's number of reference points
        integer(int32), intent(in) :: n_neighbors
            !! The candidate's number of neighbors per reference point
        real(real64), intent(in) :: x_star(n_points)
            !! The candidate's reference points
        real(real64), intent(in) :: gene_means(max_n_genes_all_studies, n_studies)
            !! Per-gene mean expression values for all studies
        integer(int32), intent(in) :: gene_means_perms(max_n_genes_all_studies, n_studies)
            !! Each study's sorting permutation for `gene_means`
        real(real64), intent(in) :: min_neighbor_overlap
            !! Forwarded to check_neighborhood_overlaps_impl
        integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points, n_studies)
            !! Every study's neighbor gene indices; only the studies up to the first failing one
            !! are filled
        integer(int32), intent(out) :: neighborhood_range(2, n_points)
            !! One study's neighborhood spans, reused per study
        logical(c_bool), intent(out) :: all_have_min_neighbor_overlap
            !! `.true.` iff every study passed the first admissibility gate

        integer(int32) :: i_study

        all_have_min_neighbor_overlap = logical(.true., kind=c_bool)
        do i_study = 1, n_studies
            call construct_neighborhoods_ranged_impl(n_points, x_star, max_n_genes_all_studies, gene_means(:, i_study), &
                                                     gene_means_perms(:, i_study), n_neighbors, &
                                                     neighborhood_indices(:, :, i_study), neighborhood_range)
            call check_neighborhood_overlaps_impl(neighborhood_range, n_points, min_neighbor_overlap, &
                                                   all_have_min_neighbor_overlap)
            if (.not. all_have_min_neighbor_overlap) exit
        end do
    end subroutine construct_fixed_k_candidate_neighborhoods

    !> summary: Search a GAMMA-decay (n_points, n_neighbors) candidate grid for a stable JSD parameter setting
    !| AUTHOR_LASZLO_LANG
    !| Walks the candidates of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]
    !| from finest to coarsest resolution. For each candidate:
    !|
    !| 1. Reference points from the pooled gene means
    !|    ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) and every
    !|    study's neighborhoods; gate 1: consecutive neighborhoods overlap enough
    !|    ([[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl(interface)]]).
    !| 2. Per reference point, the bin-count search `bin_count_mode` names on all studies' pooled
    !|    residuals, then every study's histograms.
    !| 3. Gate 2: every study has a non-NaN residual at every point, and every consensus bin reaches
    !|    `min_residuals_per_bin`
    !|    ([[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !|    A candidate that passes both gates is *admissible*.
    !| 4. The observed global JSD per study, its bootstrap confidence interval
    !|    ([[tox_data_integration_stats_impl(module):bootstrap_histogram_impl(interface)]]),
    !|    and both plateau criteria: CI overlap with the best candidate so far
    !|    ([[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]])
    !|    and relative effect-size stability
    !|    ([[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]]).
    !|    Both are always computed and traced; `plateau_mode` decides which one stops the search. A
    !|    stop by the effect-size criterion returns the candidate that triggered it.
    !|
    !| The returned candidate:
    !|
    !| 1. A plateau was found, or the grid has a single candidate and it is admissible: the best
    !|    candidate, `plateau_established = .true.`.
    !| 2. No plateau, but some candidate was admissible: the admissible candidate with the smallest
    !|    median confidence-interval width across studies, `plateau_established = .false.`.
    !| 3. No candidate was admissible: the first (finest) candidate, with
    !|    `best_candidate_pair_confidence_interval = -1.0`, every bin count `m_min`, both residual
    !|    ranges `0.0`, and `plateau_established = .false.`.
    !|
    !| Per-point outputs and `trace_*` arrays are sized by the candidate bounds; see each argument
    !| for its meaningful part. A failure to create the bootstrap's random number generator is
    !| returned in `ierr` and stops the search; the outputs are then undefined.
    subroutine run_js_comp_test_parameter_search_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                       gene_means, residuals, n_bootstraps, &
                                                       join_method, max_n_points_candidate, max_n_neighbors_candidate, &
                                                       n_bootstrapping_top_k_jsds, n_points, n_neighbors, n_bins_per_point, &
                                                       shared_residual_range_low, shared_residual_range_high, &
                                                       best_candidate_pair_confidence_interval, plateau_established, &
                                                       n_admissible_evaluated, &
                                                       trace_n_points, trace_n_neighbors, trace_global_js_divergence, &
                                                       trace_ci_lower, trace_ci_upper, trace_ci_width, &
                                                       trace_ci_width_relative, trace_ci_overlap, trace_delta, &
                                                       trace_delta_median, trace_delta_max, trace_selected_n_bins, &
                                                       trace_occupancy_failed, &
                                                       trace_n_pooled_residuals, trace_min_bin_occupancy, &
                                                       trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                                       trace_sturges_bins, trace_fd_bins, &
                                                       trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                                       tmp_gene_means_perms, &
                                                       tmp_gene_means_perm_all, &
                                                       tmp_x_star, tmp_neighborhood_indices_all_studies, &
                                                       tmp_neighborhood_range, tmp_neighborhood_residuals_gathered, &
                                                       tmp_point_pool, tmp_point_draw, tmp_n_bins_per_point, &
                                                       tmp_shared_residual_range_low, tmp_shared_residual_range_high, &
                                                       tmp_pmfs, tmp_counts, &
                                                       tmp_included_n_reps, tmp_mean_pmf, tmp_mean_pmf_counts, &
                                                       tmp_mean_pmf_included_n_reps, tmp_js_divergences, tmp_weights, &
                                                       tmp_global_js_divergence, tmp_confidence_interval, &
                                                       tmp_bootstrapping_top_k_jsds, tmp_prev_global_js_divergence, &
                                                       tmp_delta_perm, tmp_best_uncertainty_confidence_interval, &
                                                       tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                                       tmp_bin_counts_search, tmp_residual_to_study_perm, &
                                                       tmp_study_block_ends, tmp_study_block_sizes, tmp_occupancy_failed, &
                                                       tmp_n_pooled_residuals, tmp_min_bin_occupancy, &
                                                       tmp_mean_bin_occupancy, tmp_max_bin_occupancy, tmp_sturges_bins, &
                                                       tmp_fd_bins, tmp_best_n_bins_per_point, &
                                                       tmp_best_uncertainty_n_bins_per_point, &
                                                       tmp_best_shared_residual_range_low, tmp_best_shared_residual_range_high, &
                                                       tmp_best_uncertainty_shared_residual_range_low, &
                                                       tmp_best_uncertainty_shared_residual_range_high, &
                                                       min_residuals_per_bin, min_neighbor_overlap, &
                                                       succeeding_ci_overlap, plateau_mode, delta_median_threshold, &
                                                       delta_max_threshold, delta_epsilon, &
                                                       delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, &
                                                       lower_residual_range_quantile, upper_residual_range_quantile, &
                                                       two_sided_bootstrapping_significance_level, random_seed, &
                                                       bin_count_mode, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! DM_MIN(1_int32)
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
            !! DM_ALLOW_NAN
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps to perform for a candidate pair
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: join_method
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(int32), intent(in) :: max_n_points_candidate
            !! Upper bound on the grid's `n_points` candidates, which sizes the per-point outputs;
            !! at least what `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_js_comp_test_candidate_bounds, tox_data_integration_js_comp_test_fixed_k_search_impl, JUST_INFO)
            !! DM_MIN(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)))
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Upper bound on the grid's `n_neighbors` candidates; at least what
            !! `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! DM_OUTPUT_FROM(max_n_neighbors_candidate, calc_js_comp_test_candidate_bounds, tox_data_integration_js_comp_test_fixed_k_search_impl, JUST_INFO)
            !! DM_MIN(max(1_int32, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2_real64*real(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32)))
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! DM_OUTPUT_FROM(n_top_k, calc_js_comp_test_n_top_k_jsds, tox_data_integration_stats_impl, AUTO)
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: n_points
            !! The finally chosen candidate's `n_points`
        integer(int32), intent(out) :: n_neighbors
            !! The finally chosen candidate's `n_neighbors`
        integer(int32), dimension(max_n_points_candidate), intent(out) :: n_bins_per_point
            !! The chosen candidate's bin count per reference point; the first `n_points` entries
            !! are meaningful, the rest are 0
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_low
            !! The chosen candidate's lower residual-range bound (R_low) per reference point; the
            !! first `n_points` entries are meaningful, `0.0` when no candidate was admissible; the
            !! rest are `0.0`
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_high
            !! The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`
        real(real64), dimension(2, n_studies), intent(out) :: best_candidate_pair_confidence_interval
            !! The chosen candidate's bootstrapped JSD confidence interval; `-1.0` when no candidate
            !! was admissible
        logical(c_bool), intent(out) :: plateau_established
            !! `.true.` when a plateau was found (by the criterion `plateau_mode` selects) or the
            !! only candidate was admissible; `.false.` otherwise, also when the smallest-uncertainty
            !! fallback returns a real candidate and interval
        integer(int32), intent(out) :: n_admissible_evaluated
            !! Number of admissible candidates evaluated before the search stopped: the number of
            !! valid columns of every `trace_*` array (an inadmissible candidate has no entry)
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_n_points
            !! Per-admissible-candidate `n_points`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! (`L_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! (`U_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width (`W_{i,t} = U_{i,t} -
            !! L_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width (the values
            !! `succeeding_ci_overlap` is tested against); `-1.0` for the first admissible candidate
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`); `-1.0` for the first admissible candidate
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
            !! `1:trace_n_points(t)` of column `t` are meaningful (the rest are 0), so a Python/R
            !! caller slices `[:trace_n_points[t], t]`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        logical(c_bool), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins (`.false.`
            !! beyond)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled residual count (`N_j`). Jagged
            !! like trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic. Jagged like
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic.
            !! Jagged like trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low).
            !! Jagged like trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high).
            !! Jagged like trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means_perms
            !! Working array: each study's own sorting permutation for `gene_means`
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Working array: the current candidate's reference points
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: &
            tmp_neighborhood_indices_all_studies
            !! Working array: every study's neighbor gene indices for the current candidate
        integer(int32), dimension(2, max_n_points_candidate), intent(out) :: tmp_neighborhood_range
            !! Working array: one study's neighborhood spans, reused per study
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), &
            intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residuals, reused per study
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_pool
            !! Working array of the bootstrap: one point's pool
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_draw
            !! Working array of the bootstrap: one study's draw at one point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array: the current candidate's bin count per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array: the current candidate's R_low per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array: the current candidate's R_high per reference point
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_pmfs
            !! Working array: every study's pmf for the current candidate
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_counts
            !! Working array: every study's histogram counts for the current candidate
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array: every study's included-replicate counts for the current candidate
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf
            !! Working array: the current candidate's consensus pmf
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf_counts
            !! Working array: the current candidate's consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array: the current candidate's consensus included-replicate counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array: the current candidate's per-point JSD values, then bootstrap scratch
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array: the current candidate's per-point weights, then bootstrap scratch
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array: the current candidate's observed global JSD, then bootstrap scratch
        real(real64), dimension(2, n_studies), intent(out) :: tmp_confidence_interval
            !! Working array: the current candidate's confidence interval, bootstrapped in place
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array: the bootstrap's bottom-k/top-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array: the previous admissible candidate's observed global JSD per study
        integer(int32), dimension(n_studies), intent(out) :: tmp_delta_perm
            !! Working array: sorting permutation for the per-study medians; before that, every
            !! study's neighbor count for the bin-count search
        real(real64), dimension(2, n_studies), intent(out) :: tmp_best_uncertainty_confidence_interval
            !! Working array: confidence interval of the smallest-uncertainty admissible candidate
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array: the last index of each study's block in tmp_pooled_residuals
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in tmp_pooled_residuals
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
            !! Working array: the current candidate's occupancy_failed flag per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
            !! Working array: the current candidate's pooled residual count (N_j) per point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
            !! Working array: the current candidate's minimum bin occupancy per point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
            !! Working array: the current candidate's mean bin occupancy per point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
            !! Working array: the current candidate's maximum bin occupancy per point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
            !! Working array: the current candidate's Sturges' rule diagnostic per point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_fd_bins
            !! Working array: the current candidate's Freedman-Diaconis diagnostic per point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_n_bins_per_point
            !! Working array: bin counts per point of the current best candidate
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_n_bins_per_point
            !! Working array: bin counts per point of the smallest-uncertainty admissible candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_low
            !! Working array: R_low per point of the current best candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_high
            !! Working array: R_high per point of the current best candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array: R_low per point of the smallest-uncertainty admissible candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array: R_high per point of the smallest-uncertainty admissible candidate
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum count of every consensus bin (gate 2) and the floor of the bin-count search,
            !! one threshold for both, so a point the search marks as failed also fails the gate
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_MIN_NEIGHBOR_OVERLAP_DEFAULT)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap of a candidate's confidence interval with the best one so
            !! far, per `join_method`. `0` makes the first admissible candidate plateau under the
            !! CI-overlap criterion (for baseline sweeps)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_SUCCEEDING_CI_OVERLAP_DEFAULT)
        integer(int32), intent(in), optional :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(M_MODE_PLATEAU_CI_OVERLAP)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Upper bound on the median relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Upper bound on the largest relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Floor of the denominators of the relative JSD change and relative CI width
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions for an effect-size plateau
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(M_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! M_LOWER_RESIDUAL_RANGE_QUANTILE_ORDER_DOC
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! DM_MIN(lower_residual_range_quantile)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(M_RANDOM_SEED_DEFAULT)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen, for every candidate; the
            !! Stone modes test every bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(M_MODE_BIN_COUNT_OCCUPANCY)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if the bootstrap's random number generator could not be created

        integer(int32) :: candidates_n_points_n_neighbors(2, MAX_CANDIDATE_PAIRS)
        integer(int32) :: i_candidate, i_point, n_pool, max_n_bins
        integer(int32) :: prev_n_points, best_candidate_index, best_exceeded_ci_overlap_count, n_candidates
        integer(int32) :: pool_size, actual_min_residuals_per_bin, selected_candidate_index
        integer(int32) :: actual_plateau_mode, actual_delta_min_consecutive_transitions, n_consecutive_effect_size_ok
        integer(int32) :: best_uncertainty_candidate_index, actual_m_min
        real(real64) :: actual_min_neighbor_overlap, actual_succeeding_ci_overlap
        real(real64) :: actual_delta_median_threshold, actual_delta_max_threshold, actual_delta_epsilon
        real(real64) :: best_uncertainty_value
        logical(c_bool) :: all_have_min_neighbor_overlap, all_bins_have_min_count, plateau_found, has_empty_study

        call set_ok(ierr)
        M_DEFAULT_VAL(min_residuals_per_bin, actual_min_residuals_per_bin, M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        M_DEFAULT_VAL(min_neighbor_overlap, actual_min_neighbor_overlap, M_MIN_NEIGHBOR_OVERLAP_DEFAULT)
        M_DEFAULT_VAL(succeeding_ci_overlap, actual_succeeding_ci_overlap, M_SUCCEEDING_CI_OVERLAP_DEFAULT)
        M_DEFAULT_VAL(plateau_mode, actual_plateau_mode, M_MODE_PLATEAU_CI_OVERLAP)
        M_DEFAULT_VAL(delta_median_threshold, actual_delta_median_threshold, M_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        M_DEFAULT_VAL(delta_max_threshold, actual_delta_max_threshold, M_DELTA_MAX_THRESHOLD_DEFAULT)
        M_DEFAULT_VAL(delta_epsilon, actual_delta_epsilon, M_DELTA_EPSILON_DEFAULT)
        M_DEFAULT_VAL(delta_min_consecutive_transitions, actual_delta_min_consecutive_transitions, M_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        M_DEFAULT_VAL(m_min, actual_m_min, M_OCCUPANCY_M_MIN_DEFAULT)

        pool_size = max_n_genes_all_studies*n_studies

        call generate_js_comp_test_candidates_impl(max_n_genes_all_studies, candidates_n_points_n_neighbors, n_candidates)

        ! Sort each study's own gene means (for construct_neighborhoods_ranged_impl) and the
        ! flattened, all-studies-pooled gene means (for pool_means_impl's x_star).
        call sort_gene_means_per_study(gene_means, max_n_genes_all_studies, n_studies, tmp_gene_means_perms)
        call sort_pooled_gene_means(gene_means, max_n_genes_all_studies, n_studies, tmp_gene_means_perm_all)

        best_candidate_pair_confidence_interval = -1.0_real64
        best_candidate_index = 1_int32
        best_exceeded_ci_overlap_count = 0_int32
        best_uncertainty_value = huge(1.0_real64)
        best_uncertainty_candidate_index = 1_int32
        plateau_found = logical(.false., kind=c_bool)
        prev_n_points = -1_int32
        n_admissible_evaluated = 0_int32
        n_consecutive_effect_size_ok = 0_int32

        ! Test candidate pairs, finest resolution to coarsest, and stop at the first JSD plateau.
        do i_candidate = 1, n_candidates
            n_points = candidates_n_points_n_neighbors(1, i_candidate)
            n_neighbors = candidates_n_points_n_neighbors(2, i_candidate)

            if (prev_n_points /= n_points) then
                prev_n_points = n_points
                call pool_means_impl(gene_means, tmp_gene_means_perm_all, pool_size, n_points, n_pool, &
                                     tmp_x_star(1:n_points))
            end if

            ! ===== PASS A (per study): neighborhood construction + first admissibility gate.
            !
            ! The work arrays sized by an upper bound (max_n_neighbors_candidate,
            ! max_n_points_candidate, M_MAX_N_BINS) are used at the candidate's own sizes. They are
            ! always passed whole to explicit-shape dummies declared with the candidate's sizes, so
            ! every routine sees their leading elements packed (sequence association) and nothing is
            ! copied; every use of such a buffer goes through such a dummy, so all agree on the
            ! layout. This is why Passes A to C and the second gate live in separate, non-published helpers.
            call construct_fixed_k_candidate_neighborhoods(n_studies, max_n_genes_all_studies, n_points, n_neighbors, &
                                                           tmp_x_star, gene_means, tmp_gene_means_perms, &
                                                           actual_min_neighbor_overlap, &
                                                           tmp_neighborhood_indices_all_studies, tmp_neighborhood_range, &
                                                           all_have_min_neighbor_overlap)
            if (.not. all_have_min_neighbor_overlap) cycle

            ! ===== PASS B (per point): pool all studies' residuals and run the bin-count search, as
            ! run_js_comp_test_impl does. Sequential although the points are independent, so one
            ! tmp_pooled_residuals buffer serves every point (do concurrent would need one per point).
            ! tmp_delta_perm is free until the candidate is evaluated and holds every study's n_neighbors.
            tmp_delta_perm = n_neighbors
            do i_point = 1, n_points
                call determine_point_bin_count(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                               n_neighbors, n_points, tmp_neighborhood_indices_all_studies, i_point, &
                                               tmp_delta_perm, tmp_n_bins_per_point(i_point), tmp_occupancy_failed(i_point), &
                                               tmp_shared_residual_range_low(i_point), &
                                               tmp_shared_residual_range_high(i_point), tmp_n_pooled_residuals(i_point), &
                                               tmp_min_bin_occupancy(i_point), tmp_mean_bin_occupancy(i_point), &
                                               tmp_max_bin_occupancy(i_point), tmp_sturges_bins(i_point), &
                                               tmp_fd_bins(i_point), tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                               tmp_bin_counts_search, tmp_residual_to_study_perm, tmp_study_block_ends, &
                                               tmp_study_block_sizes, m_min=m_min, m_max=m_max, &
                                               min_residuals_per_bin=actual_min_residuals_per_bin, &
                                               gamma_occupancy=gamma_occupancy, &
                                               lower_residual_range_quantile=lower_residual_range_quantile, &
                                               upper_residual_range_quantile=upper_residual_range_quantile, &
                                               bin_count_mode=bin_count_mode)
            end do

            ! ===== PASS C (per study): every study's histograms at the chosen bin counts.
            max_n_bins = maxval(tmp_n_bins_per_point(1:n_points))
            call build_fixed_k_histograms(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                          n_neighbors, n_points, max_n_bins, max_n_bins, &
                                          tmp_neighborhood_indices_all_studies, tmp_shared_residual_range_low, &
                                          tmp_shared_residual_range_high, tmp_n_bins_per_point, tmp_counts, &
                                          tmp_pmfs, tmp_included_n_reps, tmp_neighborhood_residuals_gathered)

            ! Second admissibility gate; an empty study (an error in run_js_comp_test_impl) only
            ! makes the candidate inadmissible here.
            call score_candidate_histograms(n_points, max_n_bins, n_studies, tmp_n_bins_per_point, &
                                            actual_min_residuals_per_bin, tmp_pmfs, tmp_counts, tmp_included_n_reps, &
                                            tmp_mean_pmf, tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps, &
                                            tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                            tmp_confidence_interval, has_empty_study, all_bins_have_min_count)
            if (has_empty_study .or. .not. all_bins_have_min_count) cycle

            n_admissible_evaluated = n_admissible_evaluated + 1_int32
            trace_n_points(n_admissible_evaluated) = n_points
            trace_n_neighbors(n_admissible_evaluated) = n_neighbors
            call evaluate_admissible_candidate(n_points, max_n_bins, n_studies, max_n_points_candidate, i_candidate, &
                                               n_admissible_evaluated, n_bootstraps, n_bootstrapping_top_k_jsds, &
                                               join_method, actual_plateau_mode, actual_succeeding_ci_overlap, &
                                               actual_delta_median_threshold, actual_delta_max_threshold, &
                                               actual_delta_epsilon, actual_delta_min_consecutive_transitions, &
                                               tmp_n_bins_per_point, tmp_occupancy_failed, tmp_n_pooled_residuals, &
                                               tmp_min_bin_occupancy, tmp_mean_bin_occupancy, tmp_max_bin_occupancy, &
                                               tmp_sturges_bins, tmp_fd_bins, tmp_shared_residual_range_low, &
                                               tmp_shared_residual_range_high, tmp_included_n_reps, tmp_mean_pmf_counts, &
                                               tmp_mean_pmf_included_n_reps, tmp_counts, tmp_pmfs, tmp_mean_pmf, &
                                               tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                               tmp_confidence_interval, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_ci_overlap, &
                                               trace_delta, trace_delta_median, trace_delta_max, trace_selected_n_bins, &
                                               trace_occupancy_failed, trace_n_pooled_residuals, trace_min_bin_occupancy, &
                                               trace_mean_bin_occupancy, trace_max_bin_occupancy, trace_sturges_bins, &
                                               trace_fd_bins, trace_shared_residual_range_low, &
                                               trace_shared_residual_range_high, tmp_prev_global_js_divergence, &
                                               n_consecutive_effect_size_ok, best_candidate_pair_confidence_interval, &
                                               best_candidate_index, best_exceeded_ci_overlap_count, &
                                               tmp_best_n_bins_per_point, tmp_best_shared_residual_range_low, &
                                               tmp_best_shared_residual_range_high, best_uncertainty_value, &
                                               best_uncertainty_candidate_index, tmp_best_uncertainty_confidence_interval, &
                                               tmp_best_uncertainty_n_bins_per_point, &
                                               tmp_best_uncertainty_shared_residual_range_low, &
                                               tmp_best_uncertainty_shared_residual_range_high, plateau_found, &
                                               tmp_bootstrapping_top_k_jsds, tmp_point_pool, tmp_point_draw, &
                                               tmp_delta_perm, ierr, two_sided_bootstrapping_significance_level, &
                                               random_seed)
            if (is_err(ierr)) return
            if (plateau_found) exit
        end do

        call select_final_candidate(max_n_points_candidate, n_studies, n_candidates, n_admissible_evaluated, plateau_found, &
                                    actual_m_min, best_candidate_index, &
                                    candidates_n_points_n_neighbors(1, best_candidate_index), tmp_best_n_bins_per_point, &
                                    tmp_best_shared_residual_range_low, tmp_best_shared_residual_range_high, &
                                    best_uncertainty_candidate_index, &
                                    candidates_n_points_n_neighbors(1, best_uncertainty_candidate_index), &
                                    tmp_best_uncertainty_confidence_interval, tmp_best_uncertainty_n_bins_per_point, &
                                    tmp_best_uncertainty_shared_residual_range_low, &
                                    tmp_best_uncertainty_shared_residual_range_high, candidates_n_points_n_neighbors(1, 1), &
                                    selected_candidate_index, n_points, n_bins_per_point, shared_residual_range_low, &
                                    shared_residual_range_high, best_candidate_pair_confidence_interval, plateau_established)
        n_neighbors = candidates_n_points_n_neighbors(2, selected_candidate_index)
    end subroutine run_js_comp_test_parameter_search_impl

end module tox_data_integration_js_comp_test_fixed_k_search_impl
