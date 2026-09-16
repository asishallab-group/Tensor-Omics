#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search
!|
!| The data-driven `(n_points, n_neighbors)` parameter-stabilization search this pipeline runs
!| before the JSD-Comp-Test proper (Issue #126): a GAMMA-decay candidate grid
!| ([[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]],
!| each candidate's histogram bin count from
!| [[tox_data_integration_js_comp_test_impl(module):estimate_bin_count_impl(interface)]]), two
!| admissibility gates a candidate must pass before it is bootstrapped
!| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]],
!| [[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]),
!| and the plateau check that decides when the search has converged
!| ([[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]).
!| `calc_js_comp_test_candidate_bounds` sizes the candidate-grid work arrays for a caller that
!| allocates its own. Once a candidate has passed both gates, its bootstrap confidence interval
!| is resampled from the pooled consensus histogram by
!| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] (heap
!| size recommended by
!| [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds(interface)]]).
!| Later stages of the port add the K-study permutation test and the top-level orchestrator that
!| wire these building blocks together.
module tox_data_integration_js_comp_test_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_as_view
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: estimate_bin_count_c
    public :: estimate_bin_count_expert_c
    public :: generate_js_comp_test_candidates_c
    public :: generate_js_comp_test_candidates_expert_c
    public :: check_neighborhood_overlaps_c
    public :: check_mean_pmf_min_counts_c
    public :: check_plateau_condition_c
    public :: check_effect_size_plateau_condition_c
    public :: check_effect_size_plateau_condition_expert_c
    public :: create_mean_pmf_c
    public :: create_mean_pmf_only_c
    public :: bootstrap_histogram_c
    public :: bootstrap_histogram_expert_c
    public :: run_js_comp_test_c
    public :: run_js_comp_test_expert_c
    public :: run_js_comp_test_parameter_search_c
    public :: run_js_comp_test_parameter_search_expert_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):estimate_bin_count(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `estimate_bin_count_helper`. Takes the maximum of
    !| Sturges' rule and the Freedman-Diaconis rule (without doubling the bin width, since
    !| `shared_residual_range` is already the one-sided half of the full `[-R, R]` histogram
    !| range, so dividing the full range by the undoubled Freedman-Diaconis width already gives
    !| the doubled rule's bin count), clamped to at most
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] bins.
    subroutine estimate_bin_count_c(&
            residuals,&
            n_residuals,&
            max_n_reps_all_studies,&
            n_neighbors,&
            shared_residual_range,&
            n_bins,&
            ierr&
        ) bind(C, name="estimate_bin_count_c")
        use tox_data_integration_js_comp_test, only: estimate_bin_count

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        real(c_double), dimension(n_residuals), intent(in), target :: residuals
            !! Pooled signed residuals across all studies, reference points and neighbors
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), intent(out), target :: n_bins
            !! Estimated number of histogram bins, at least 1 and at most MAX_N_BINS
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_ARRAY_NON_NULL(residuals, n_residuals)

        call estimate_bin_count(&
            residuals = residuals,&
            n_residuals = n_residuals,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            shared_residual_range = shared_residual_range,&
            n_bins = n_bins,&
            ierr = ierr&
        )
    end subroutine estimate_bin_count_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):estimate_bin_count_expert(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `estimate_bin_count_helper`. Takes the maximum of
    !| Sturges' rule and the Freedman-Diaconis rule (without doubling the bin width, since
    !| `shared_residual_range` is already the one-sided half of the full `[-R, R]` histogram
    !| range, so dividing the full range by the undoubled Freedman-Diaconis width already gives
    !| the doubled rule's bin count), clamped to at most
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] bins.
    subroutine estimate_bin_count_expert_c(&
            residuals,&
            residuals_perm,&
            n_residuals,&
            max_n_reps_all_studies,&
            n_neighbors,&
            shared_residual_range,&
            n_bins,&
            ierr&
        ) bind(C, name="estimate_bin_count_expert_c")
        use tox_data_integration_js_comp_test, only: estimate_bin_count_expert

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        real(c_double), dimension(n_residuals), intent(in), target :: residuals
            !! Pooled signed residuals across all studies, reference points and neighbors
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_residuals), intent(in), target :: residuals_perm
            !! Sorting permutation for `residuals`, ascending, NaN last
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_residuals`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), intent(out), target :: n_bins
            !! Estimated number of histogram bins, at least 1 and at most MAX_N_BINS
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_ARRAY_NON_NULL(residuals, n_residuals)
        M_CHECK_ARRAY_NON_NULL(residuals_perm, n_residuals)

        call estimate_bin_count_expert(&
            residuals = residuals,&
            residuals_perm = residuals_perm,&
            n_residuals = n_residuals,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            shared_residual_range = shared_residual_range,&
            n_bins = n_bins,&
            ierr = ierr&
        )
    end subroutine estimate_bin_count_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):generate_js_comp_test_candidates(subroutine)]]
    !| Ported from the grid-building half of 125-stabilize-jscomp's
    !| `determine_js_comp_test_n_points_n_neighbors_helper`: starting from an initial
    !| `n_points_high` (clamped between MIN_POINTS and MAX_POINTS), repeatedly multiplies by
    !| GAMMA until it would drop below `n_points_low`, and for each distinct resulting
    !| `n_points` candidate pairs it with up to `size(KX_FACTORS)` distinct `n_neighbors`
    !| candidates, calling
    !| [[tox_data_integration_js_comp_test_impl(module):estimate_bin_count_impl(interface)]] for
    !| each pair's bin count. A duplicate `n_points` or `n_neighbors` value (from clamping or
    !| integer rounding) collapses rather than repeating -- this is real, derived behavior the
    !| grid depends on to avoid redundant candidates at small `max_n_genes_all_studies`, not a
    !| bug: a small enough `max_n_genes_all_studies` collapses the whole grid down to exactly one
    !| candidate.
    subroutine generate_js_comp_test_candidates_c(&
            max_n_genes_all_studies,&
            residuals,&
            n_residuals,&
            max_n_reps_all_studies,&
            shared_residual_range,&
            candidates_n_points_n_neighbors,&
            n_bins_candidates,&
            n_candidates,&
            ierr&
        ) bind(C, name="generate_js_comp_test_candidates_c")
        use tox_data_integration_js_comp_test, only: generate_js_comp_test_candidates

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_residuals), intent(in), target :: residuals
            !! Pooled signed residuals across all studies, reference points and neighbors
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(2, 16), intent(out), target :: candidates_n_points_n_neighbors
            !! Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
            !! The first `n_candidates` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: n_bins_candidates
            !! Per-candidate bin count from estimate_bin_count_impl, one per candidate pair
            !! The first `n_candidates` elements will hold the results.
        integer(c_int), intent(out), target :: n_candidates
            !! Number of candidate pairs actually filled (at most MAX_CANDIDATE_PAIRS = 16)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_candidates)
        M_CHECK_ARRAY_NON_NULL(residuals, n_residuals)
        M_CHECK_ARRAY_NON_NULL(candidates_n_points_n_neighbors, 2 * 16)
        M_CHECK_ARRAY_NON_NULL(n_bins_candidates, 16)

        call generate_js_comp_test_candidates(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            residuals = residuals,&
            n_residuals = n_residuals,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            shared_residual_range = shared_residual_range,&
            candidates_n_points_n_neighbors = candidates_n_points_n_neighbors,&
            n_bins_candidates = n_bins_candidates,&
            n_candidates = n_candidates,&
            ierr = ierr&
        )
    end subroutine generate_js_comp_test_candidates_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):generate_js_comp_test_candidates_expert(subroutine)]]
    !| Ported from the grid-building half of 125-stabilize-jscomp's
    !| `determine_js_comp_test_n_points_n_neighbors_helper`: starting from an initial
    !| `n_points_high` (clamped between MIN_POINTS and MAX_POINTS), repeatedly multiplies by
    !| GAMMA until it would drop below `n_points_low`, and for each distinct resulting
    !| `n_points` candidate pairs it with up to `size(KX_FACTORS)` distinct `n_neighbors`
    !| candidates, calling
    !| [[tox_data_integration_js_comp_test_impl(module):estimate_bin_count_impl(interface)]] for
    !| each pair's bin count. A duplicate `n_points` or `n_neighbors` value (from clamping or
    !| integer rounding) collapses rather than repeating -- this is real, derived behavior the
    !| grid depends on to avoid redundant candidates at small `max_n_genes_all_studies`, not a
    !| bug: a small enough `max_n_genes_all_studies` collapses the whole grid down to exactly one
    !| candidate.
    subroutine generate_js_comp_test_candidates_expert_c(&
            max_n_genes_all_studies,&
            residuals,&
            residuals_perm,&
            n_residuals,&
            max_n_reps_all_studies,&
            shared_residual_range,&
            candidates_n_points_n_neighbors,&
            n_bins_candidates,&
            n_candidates,&
            ierr&
        ) bind(C, name="generate_js_comp_test_candidates_expert_c")
        use tox_data_integration_js_comp_test, only: generate_js_comp_test_candidates_expert

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_residuals), intent(in), target :: residuals
            !! Pooled signed residuals across all studies, reference points and neighbors
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_residuals), intent(in), target :: residuals_perm
            !! Sorting permutation for `residuals`, ascending, NaN last
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_residuals`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(2, 16), intent(out), target :: candidates_n_points_n_neighbors
            !! Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
            !! The first `n_candidates` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: n_bins_candidates
            !! Per-candidate bin count from estimate_bin_count_impl, one per candidate pair
            !! The first `n_candidates` elements will hold the results.
        integer(c_int), intent(out), target :: n_candidates
            !! Number of candidate pairs actually filled (at most MAX_CANDIDATE_PAIRS = 16)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_candidates)
        M_CHECK_ARRAY_NON_NULL(residuals, n_residuals)
        M_CHECK_ARRAY_NON_NULL(residuals_perm, n_residuals)
        M_CHECK_ARRAY_NON_NULL(candidates_n_points_n_neighbors, 2 * 16)
        M_CHECK_ARRAY_NON_NULL(n_bins_candidates, 16)

        call generate_js_comp_test_candidates_expert(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            residuals = residuals,&
            residuals_perm = residuals_perm,&
            n_residuals = n_residuals,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            shared_residual_range = shared_residual_range,&
            candidates_n_points_n_neighbors = candidates_n_points_n_neighbors,&
            n_bins_candidates = n_bins_candidates,&
            n_candidates = n_candidates,&
            ierr = ierr&
        )
    end subroutine generate_js_comp_test_candidates_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_neighborhood_overlaps(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `test_neighborhood_overlaps_helper`: the first
    !| admissibility gate a candidate `(n_points, n_neighbors)` pair must pass, using the
    !| `[min_idx, max_idx]` neighborhood spans
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| produces. Named `check_*` rather than 125's `test_*`, a deliberate deviation from the
    !| verbatim port: the generated R binding is published under the Fortran name, and the R test
    !| harness (`r/test_helpers.R`'s `run_all_tests`) discovers every `test_`-prefixed name in the
    !| environment as a test case to run with no arguments -- a `test_`-prefixed export would be
    !| swept up and fail every R suite that sources the package, not just this module's own.
    subroutine check_neighborhood_overlaps_c(&
            neighborhood_range,&
            n_points,&
            min_neighbor_overlap,&
            all_have_min_neighbor_overlap,&
            ierr&
        ) bind(C, name="check_neighborhood_overlaps_c")
        use tox_data_integration_js_comp_test, only: check_neighborhood_overlaps

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), dimension(2, n_points), intent(in), target :: neighborhood_range
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span, as produced
            !! by construct_neighborhoods_ranged_impl
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_mean_pmf_min_counts(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `test_mean_pmf_min_counts_helper`: the second
    !| admissibility gate a candidate `(n_points, n_neighbors)` pair must pass, checked once the
    !| first gate
    !| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]])
    !| has already passed. Named `check_*` rather than 125's `test_*` for the same reason as its
    !| sibling above: a `test_`-prefixed R export collides with the R test harness's own
    !| test-discovery convention.
    subroutine check_mean_pmf_min_counts_c(&
            mean_pmf_counts,&
            n_bins,&
            n_points,&
            min_count,&
            all_bins_have_min_count,&
            ierr&
        ) bind(C, name="check_mean_pmf_min_counts_c")
        use tox_data_integration_js_comp_test, only: check_mean_pmf_min_counts

        integer(c_int), intent(in), target :: n_bins
            !! Number of histogram bins
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), dimension(n_bins, n_points), intent(in), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), intent(in), target :: min_count
            !! Minimum count each bin of the mean pmf must reach
            !! The minimum valid value is `0_int32`.
        logical(c_bool), intent(out), target :: all_bins_have_min_count
            !! `.true.` if every bin, at every reference point, reaches at least `min_count`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(min_count)
        M_CHECK_NON_NULL(all_bins_have_min_count)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)

        call check_mean_pmf_min_counts(&
            mean_pmf_counts = mean_pmf_counts,&
            n_bins = n_bins,&
            n_points = n_points,&
            min_count = min_count,&
            all_bins_have_min_count = all_bins_have_min_count,&
            ierr = ierr&
        )
    end subroutine check_mean_pmf_min_counts_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_plateau_condition(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `check_plateau_condition_helper`. A search over
    !| candidates (finest resolution to coarsest) stops -- "plateaus" -- either when a new
    !| candidate is no better than the previous best (the short-circuit below: keep the previous
    !| best and stop searching), or once the new candidate's confidence-interval overlap with the
    !| previous best meets the condition `join_method` names. `join_method` replaces
    !| 125-stabilize-jscomp's hand-rolled join-method-range validation macro entirely: the mode
    !| table below is itself the validation, checked against exactly the values it names.
    subroutine check_plateau_condition_c(&
            confidence_interval,&
            best_candidate_pair_confidence_interval,&
            n_studies,&
            best_candidate_index,&
            best_exceeded_ci_overlap_count,&
            candidate_index,&
            join_method,&
            succeeding_ci_overlap,&
            plateau_found,&
            ierr&
        ) bind(C, name="check_plateau_condition_c")
        use tox_data_integration_js_comp_test, only: check_plateau_condition
        use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(2, n_studies), intent(in), target :: confidence_interval
            !! JSD confidence interval `[lower, upper]` from bootstrapping, for the candidate
            !! pair under test
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
            !! The way to evaluate all studies' confidence-interval overlaps for the plateau
            !! condition: METHOD_JOIN_MIN requires every study's overlap to exceed
            !! `succeeding_ci_overlap`, METHOD_JOIN_MAX requires only one study's overlap to
            !! exceed it, and METHOD_JOIN_MEDIAN requires a majority
            !! (`count > (n_studies - 1) / 2`) to exceed it
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Minimum fractional overlap an interval in `confidence_interval` must have with its
            !! respective interval in `best_candidate_pair_confidence_interval` to count as
            !! "exceeded"
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
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
            plateau_found = plateau_found,&
            ierr = ierr&
        )
    end subroutine check_plateau_condition_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_effect_size_plateau_condition(subroutine)]]
    !| Implements Issue #178's relative-effect-size plateau criterion, complementary to
    !| [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]'s
    !| CI-overlap one: for each study `i`, the relative change in observed JSD between successive
    !| ADMISSIBLE parameter settings (both admissibility gates already passed),
    !| `delta(i) = |global_js_divergence(i) - prev_global_js_divergence(i)| / max(prev_global_js_divergence(i),
    !| delta_epsilon)`, summarized across studies by its median (`delta_median`, via the
    !| already-shipped
    !| [[f42_stats_impl(module):calc_percentile_impl(interface)]]) and maximum (`delta_max`). A
    !| plateau is declared once both stay under their respective thresholds for
    !| `delta_min_consecutive_transitions` consecutive transitions in a row -- tracked across calls
    !| via `n_consecutive_ok`, reset the moment either threshold is missed.
    !|
    !| No transition exists for the very first admissible candidate a caller ever passes in
    !| (`has_previous = .false.`): `delta`/`delta_median`/`delta_max` are all set to
    !| `-1.0_real64` -- the same not-yet-computed sentinel
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| already uses for its own confidence-interval fallback, usable here for the same reason:
    !| every quantity this routine tracks is structurally non-negative.
    !|
    !| The 0.05/0.10 defaults `run_js_comp_test_parameter_search_impl` passes for
    !| `delta_median_threshold`/`delta_max_threshold` are Issue #178's own suggested starting
    !| point, explicitly not yet empirically validated -- see that routine's doc comment.
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
        use tox_data_integration_js_comp_test, only: check_effect_size_plateau_condition

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
            !! `.false.` for the very first admissible candidate a caller has ever passed in, where
            !! no transition exists to compute a relative change from
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound the median relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound the largest relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(inout), target :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; incremented when this
            !! transition qualifies, reset to zero otherwise (and whenever `has_previous` is
            !! `.false.`)
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_effect_size_plateau_condition_expert(subroutine)]]
    !| Implements Issue #178's relative-effect-size plateau criterion, complementary to
    !| [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]'s
    !| CI-overlap one: for each study `i`, the relative change in observed JSD between successive
    !| ADMISSIBLE parameter settings (both admissibility gates already passed),
    !| `delta(i) = |global_js_divergence(i) - prev_global_js_divergence(i)| / max(prev_global_js_divergence(i),
    !| delta_epsilon)`, summarized across studies by its median (`delta_median`, via the
    !| already-shipped
    !| [[f42_stats_impl(module):calc_percentile_impl(interface)]]) and maximum (`delta_max`). A
    !| plateau is declared once both stay under their respective thresholds for
    !| `delta_min_consecutive_transitions` consecutive transitions in a row -- tracked across calls
    !| via `n_consecutive_ok`, reset the moment either threshold is missed.
    !|
    !| No transition exists for the very first admissible candidate a caller ever passes in
    !| (`has_previous = .false.`): `delta`/`delta_median`/`delta_max` are all set to
    !| `-1.0_real64` -- the same not-yet-computed sentinel
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| already uses for its own confidence-interval fallback, usable here for the same reason:
    !| every quantity this routine tracks is structurally non-negative.
    !|
    !| The 0.05/0.10 defaults `run_js_comp_test_parameter_search_impl` passes for
    !| `delta_median_threshold`/`delta_max_threshold` are Issue #178's own suggested starting
    !| point, explicitly not yet empirically validated -- see that routine's doc comment.
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
        use tox_data_integration_js_comp_test, only: check_effect_size_plateau_condition_expert

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
            !! `.false.` for the very first admissible candidate a caller has ever passed in, where
            !! no transition exists to compute a relative change from
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound the median relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound the largest relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(inout), target :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; incremented when this
            !! transition qualifies, reset to zero otherwise (and whenever `has_previous` is
            !! `.false.`)
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):create_mean_pmf(subroutine)]]
    !| Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_helper`.
    !|
    !| Known limitation: averages over all n_studies including the study being compared against it,
    !| rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
    !| as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
    !| the fix.
    subroutine create_mean_pmf_c(&
            pmfs,&
            counts,&
            n_bins,&
            n_points,&
            n_studies,&
            included_n_reps,&
            mean_pmf,&
            mean_pmf_included_n_reps,&
            mean_pmf_counts,&
            ierr&
        ) bind(C, name="create_mean_pmf_c")
        use tox_data_integration_js_comp_test, only: create_mean_pmf

        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(n_bins, n_points, n_studies), intent(in), target :: pmfs
            !! Per-study probabilities of each bin per reference point, from
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), dimension(n_bins, n_points, n_studies), intent(in), target :: counts
            !! Absolute counts of a residual per bin for `pmfs`
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_bins, n_points), intent(out), target :: mean_pmf
            !! The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
            !! known-limitation note above
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for `mean_pmf`,
            !! summed across all n_studies
        integer(c_int), dimension(n_bins, n_points), intent(out), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts, dim=3)`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)

        call create_mean_pmf(&
            pmfs = pmfs,&
            counts = counts,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            included_n_reps = included_n_reps,&
            mean_pmf = mean_pmf,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            mean_pmf_counts = mean_pmf_counts,&
            ierr = ierr&
        )
    end subroutine create_mean_pmf_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):create_mean_pmf_only(subroutine)]]
    !| Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_only_helper`: useful where the
    !| mean pmf's own counts don't matter, e.g. for the bootstrap confidence interval a later
    !| stage of this port adds.
    !|
    !| Known limitation: averages over all n_studies including the study being compared against it,
    !| rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
    !| as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
    !| the fix.
    subroutine create_mean_pmf_only_c(&
            pmfs,&
            n_bins,&
            n_points,&
            n_studies,&
            mean_pmf,&
            ierr&
        ) bind(C, name="create_mean_pmf_only_c")
        use tox_data_integration_js_comp_test, only: create_mean_pmf_only

        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(n_bins, n_points, n_studies), intent(in), target :: pmfs
            !! Per-study probabilities of each bin per reference point, from
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), dimension(n_bins, n_points), intent(out), target :: mean_pmf
            !! The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
            !! known-limitation note above
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_bins * n_points)

        call create_mean_pmf_only(&
            pmfs = pmfs,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            mean_pmf = mean_pmf,&
            ierr = ierr&
        )
    end subroutine create_mean_pmf_only_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):bootstrap_histogram(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `bootstrap_histogram_helper`. Resamples from the
    !| POOLED/consensus histogram counts (`mean_pmf_counts`), not from each study's own histogram
    !| -- this is 125's deliberate design, preserved as-is (see the plan for this port). Draws
    !| random numbers via
    !| [[f42_random_gsl(module):random_multinomial(subroutine)]], so this implementation is
    !| deliberately impure, matching the project's existing precedent for the other
    !| permutation-test module's purity
    !| ([[tox_trajectory_contribution_analysis_impl(module):perform_permutation_test_impl(interface)]]).
    !|
    !| `confidence_interval` is both an input and an output: its incoming `[lower, upper]` values
    !| seed every slot of the top-k/bottom-k heaps (`tmp_bootstrapping_top_k_jsds`) -- ported
    !| verbatim from 125, which fills the whole heap with the *same* incoming reference value
    !| rather than the usual plus/minus-infinity heap initialization, so the observed
    !| (pre-bootstrap) value can only be displaced by a strictly more extreme bootstrap draw. On
    !| return it holds `[largest of the n_bootstrapping_top_k_jsds smallest bootstrap draws,
    !| smallest of the n_bootstrapping_top_k_jsds largest bootstrap draws]`.
    !|
    !| `mean_pmf_counts`/`tmp_pmfs`/`tmp_mean_pmf` are laid out bin-major (`(n_bins, n_points[,
    !| n_studies])`), matching
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]] and its
    !| own `create_mean_pmf_only_impl` above -- both ported from 125, whose own convention this
    !| is. The already-shipped
    !| [[tox_data_integration_jsd_impl(module):compute_divergence_per_reference_point_impl(interface)]]
    !| predates 125's own port and instead takes its pmf arguments point-major
    !| (`(n_points, n_bins)`); the two calls below bridge the two conventions with an explicit
    !| `transpose`, rather than picking one shape and silently reinterpreting the other's memory
    !| under it (which would scramble every non-square `(n_bins, n_points)` histogram).
    !|
    !| A GSL allocation failure in `create_rng` is a genuine runtime error no input check could
    !| have foreseen (codegen_guide.md Sec 5.14): every work array and `confidence_interval` are
    !| then left untouched (arrays not yet written to keep their caller-visible defined state) and
    !| `ierr` reports `ERR_ALLOC_FAIL`. A `random_multinomial` draw failing (which validated,
    !| internally-consistent inputs should never trigger) is likewise folded into `ierr`, first
    !| failure only, without stopping the resampling already in flight.
    subroutine bootstrap_histogram_c(&
            n_bootstraps,&
            n_bins,&
            n_points,&
            n_studies,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            included_n_reps,&
            confidence_interval,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            ierr&
        ) bind(C, name="bootstrap_histogram_c")
        use tox_data_integration_js_comp_test, only: bootstrap_histogram

        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(n_bins, n_points), intent(in), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the pooled/consensus pmf, from
            !! create_mean_pmf_impl -- resampled with replacement each bootstrap
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the pooled pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study --
            !! how many elements are drawn (with replacement) from the pooled pool per study
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(2, n_studies), intent(inout), target :: confidence_interval
            !! Confidence interval to be bootstrapped -- incoming values are the reference values
            !! that seed the top-k/bottom-k heaps (see above); overwritten with the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds to size n_bootstrapping_top_k_jsds; not
            !! otherwise used here
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(confidence_interval, 2 * n_studies)

        call bootstrap_histogram(&
            n_bootstraps = n_bootstraps,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            included_n_reps = included_n_reps,&
            confidence_interval = confidence_interval,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine bootstrap_histogram_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):bootstrap_histogram_expert(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `bootstrap_histogram_helper`. Resamples from the
    !| POOLED/consensus histogram counts (`mean_pmf_counts`), not from each study's own histogram
    !| -- this is 125's deliberate design, preserved as-is (see the plan for this port). Draws
    !| random numbers via
    !| [[f42_random_gsl(module):random_multinomial(subroutine)]], so this implementation is
    !| deliberately impure, matching the project's existing precedent for the other
    !| permutation-test module's purity
    !| ([[tox_trajectory_contribution_analysis_impl(module):perform_permutation_test_impl(interface)]]).
    !|
    !| `confidence_interval` is both an input and an output: its incoming `[lower, upper]` values
    !| seed every slot of the top-k/bottom-k heaps (`tmp_bootstrapping_top_k_jsds`) -- ported
    !| verbatim from 125, which fills the whole heap with the *same* incoming reference value
    !| rather than the usual plus/minus-infinity heap initialization, so the observed
    !| (pre-bootstrap) value can only be displaced by a strictly more extreme bootstrap draw. On
    !| return it holds `[largest of the n_bootstrapping_top_k_jsds smallest bootstrap draws,
    !| smallest of the n_bootstrapping_top_k_jsds largest bootstrap draws]`.
    !|
    !| `mean_pmf_counts`/`tmp_pmfs`/`tmp_mean_pmf` are laid out bin-major (`(n_bins, n_points[,
    !| n_studies])`), matching
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]] and its
    !| own `create_mean_pmf_only_impl` above -- both ported from 125, whose own convention this
    !| is. The already-shipped
    !| [[tox_data_integration_jsd_impl(module):compute_divergence_per_reference_point_impl(interface)]]
    !| predates 125's own port and instead takes its pmf arguments point-major
    !| (`(n_points, n_bins)`); the two calls below bridge the two conventions with an explicit
    !| `transpose`, rather than picking one shape and silently reinterpreting the other's memory
    !| under it (which would scramble every non-square `(n_bins, n_points)` histogram).
    !|
    !| A GSL allocation failure in `create_rng` is a genuine runtime error no input check could
    !| have foreseen (codegen_guide.md Sec 5.14): every work array and `confidence_interval` are
    !| then left untouched (arrays not yet written to keep their caller-visible defined state) and
    !| `ierr` reports `ERR_ALLOC_FAIL`. A `random_multinomial` draw failing (which validated,
    !| internally-consistent inputs should never trigger) is likewise folded into `ierr`, first
    !| failure only, without stopping the resampling already in flight.
    subroutine bootstrap_histogram_expert_c(&
            n_bootstraps,&
            n_bins,&
            n_points,&
            n_studies,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            included_n_reps,&
            n_bootstrapping_top_k_jsds,&
            confidence_interval,&
            tmp_bootstrapping_top_k_jsds,&
            tmp_counts,&
            tmp_pmfs,&
            tmp_mean_pmf,&
            tmp_js_divergences,&
            tmp_weights,&
            tmp_global_js_divergence,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            ierr&
        ) bind(C, name="bootstrap_histogram_expert_c")
        use tox_data_integration_js_comp_test, only: bootstrap_histogram_expert

        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstrapping_top_k_jsds
            !! Number of elements kept at each end of the bootstrap distribution (top-k/bottom-k
            !! heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(n_bins, n_points), intent(in), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the pooled/consensus pmf, from
            !! create_mean_pmf_impl -- resampled with replacement each bootstrap
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the pooled pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study --
            !! how many elements are drawn (with replacement) from the pooled pool per study
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(2, n_studies), intent(inout), target :: confidence_interval
            !! Confidence interval to be bootstrapped -- incoming values are the reference values
            !! that seed the top-k/bottom-k heaps (see above); overwritten with the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out), target :: tmp_bootstrapping_top_k_jsds
            !! Working array used as top-k/bottom-k heaps for efficient percentile detection in
            !! the bootstrapped values -- `(:, 1, :)` the bottom-k (lower bound), `(:, 2, :)` the
            !! top-k (upper bound)
        integer(c_int), dimension(n_bins, n_points), intent(out), target :: tmp_counts
            !! Working array that holds one bootstrap's resampled histogram counts, one study at a time
        real(c_double), dimension(n_bins, n_points, n_studies), intent(out), target :: tmp_pmfs
            !! Working array that holds one bootstrap's resampled pmfs, all studies
        real(c_double), dimension(n_bins, n_points), intent(out), target :: tmp_mean_pmf
            !! Working array for one bootstrap's own (resampled) mean pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_js_divergences
            !! Working array for one bootstrap's per-point JSD values
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_weights
            !! Working array for one bootstrap's per-point global-JSD weights
        real(c_double), dimension(n_studies), intent(out), target :: tmp_global_js_divergence
            !! Working array for one bootstrap's global weighted JSD values
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds to size n_bootstrapping_top_k_jsds; not
            !! otherwise used here
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds * 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_global_js_divergence, n_studies)

        call bootstrap_histogram_expert(&
            n_bootstraps = n_bootstraps,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            included_n_reps = included_n_reps,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            confidence_interval = confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_counts = tmp_counts,&
            tmp_pmfs = tmp_pmfs,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine bootstrap_histogram_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `js_comp_test_helper`: for every study, builds its
    !| neighborhoods
    !| ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]])
    !| and residual histograms
    !| ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]), pools
    !| them into the consensus pmf
    !| ([[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]), computes
    !| each study's observed JSD against that consensus
    !| ([[tox_data_integration_jsd_impl(module):compute_divergence_per_reference_point_impl(interface)]]/[[tox_data_integration_jsd_impl(module):compute_weighted_global_divergence_impl(interface)]],
    !| called with the consensus pmf as the second argument), runs the permutation test
    !| ([[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]]), and
    !| finally re-derives each study's pmf/JSD/weights/global JSD from its own UNTOUCHED `counts`
    !| via
    !| [[tox_data_integration_jsd_impl(module):calc_pmf_impl(interface)]] -- `mean_pmf`/`mean_pmf_counts`
    !| are NOT re-derived, since they are invariant across permutations by construction (the
    !| permutation test above only resamples its own scratch copies, never `mean_pmf_counts`
    !| itself), exactly as 125 relies on.
    !|
    !| `x_star` is an ordinary input here, not computed by this routine -- 125's own
    !| `js_comp_test_helper` takes it the same way, since a caller running several studies/several
    !| parameter settings is expected to compute the reference points once
    !| ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) and reuse
    !| them consistently.
    !|
    !| `construct_neighborhoods_ranged_impl` reports neighbor gene INDICES, not gathered residual
    !| values (unlike its distance-sort sibling
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_impl(interface)]]),
    !| so this routine gathers each neighbor's actual residual values from `residuals` itself
    !| (`tmp_neighborhood_residuals_gathered`, a per-study scratch buffer) before calling
    !| `build_residual_histograms_impl`. `build_residual_histograms_impl`/`calc_pmf_impl` are
    !| POINT-major (`(n_points, n_bins)`), while `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts` here
    !| are BIN-major (`(n_bins, n_points, n_studies)`) to match
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]'s own
    !| convention -- every call across that boundary bridges with an explicit `transpose`, exactly
    !| as [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] and
    !| [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]] already do.
    !|
    !| Impure: calls the impure `gjct_permutation_test_impl`. A GSL failure it reports is folded
    !| into `ierr` (first failure only), matching that routine's own tolerant precedent.
    subroutine run_js_comp_test_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            n_points,&
            n_neighbors,&
            n_bins,&
            shared_residual_range,&
            gene_means,&
            gene_means_perms,&
            residuals,&
            x_star,&
            neighborhood_indices,&
            neighborhood_range,&
            pmfs,&
            counts,&
            included_n_reps,&
            mean_pmf,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            js_divergences,&
            weights,&
            global_js_divergence,&
            p_values,&
            n_permutations,&
            random_seed,&
            ierr&
        ) bind(C, name="run_js_comp_test_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points for neighborhoods
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors per neighborhood
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means_perms
            !! Per-study sorting permutation for `gene_means` (ascending, NaN last)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        real(c_double), dimension(n_points), intent(in), target :: x_star
            !! Mean-expression reference points
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_neighbors, n_points, n_studies), intent(out), target :: neighborhood_indices
            !! Gene indices of the selected neighborhood, per reference point, per study
        integer(c_int), dimension(2, n_points, n_studies), intent(out), target :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span, as
            !! produced by construct_neighborhoods_ranged_impl
        real(c_double), dimension(n_bins, n_points, n_studies), intent(out), target :: pmfs
            !! `counts` normalized to `0 <= pmfs(:, :, i) <= 1` and `sum(pmfs(:, j, i)) == 1`
        integer(c_int), dimension(n_bins, n_points, n_studies), intent(out), target :: counts
            !! Absolute counts of a residual per bin for `pmfs`
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study
        real(c_double), dimension(n_bins, n_points), intent(out), target :: mean_pmf
            !! The consensus pmf, from create_mean_pmf_impl
        integer(c_int), dimension(n_bins, n_points), intent(out), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study from gjct_permutation_test_impl
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations, forwarded to gjct_permutation_test_impl
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator for
            !! the permutation test

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(gene_means_perms, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_indices, n_neighbors * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(neighborhood_range, 2 * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)

        call run_js_comp_test(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            gene_means = gene_means,&
            gene_means_perms = gene_means_perms,&
            residuals = residuals,&
            x_star = x_star,&
            neighborhood_indices = neighborhood_indices,&
            neighborhood_range = neighborhood_range,&
            pmfs = pmfs,&
            counts = counts,&
            included_n_reps = included_n_reps,&
            mean_pmf = mean_pmf,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            js_divergences = js_divergences,&
            weights = weights,&
            global_js_divergence = global_js_divergence,&
            p_values = p_values,&
            n_permutations = n_permutations,&
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_expert(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `js_comp_test_helper`: for every study, builds its
    !| neighborhoods
    !| ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]])
    !| and residual histograms
    !| ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]), pools
    !| them into the consensus pmf
    !| ([[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]), computes
    !| each study's observed JSD against that consensus
    !| ([[tox_data_integration_jsd_impl(module):compute_divergence_per_reference_point_impl(interface)]]/[[tox_data_integration_jsd_impl(module):compute_weighted_global_divergence_impl(interface)]],
    !| called with the consensus pmf as the second argument), runs the permutation test
    !| ([[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]]), and
    !| finally re-derives each study's pmf/JSD/weights/global JSD from its own UNTOUCHED `counts`
    !| via
    !| [[tox_data_integration_jsd_impl(module):calc_pmf_impl(interface)]] -- `mean_pmf`/`mean_pmf_counts`
    !| are NOT re-derived, since they are invariant across permutations by construction (the
    !| permutation test above only resamples its own scratch copies, never `mean_pmf_counts`
    !| itself), exactly as 125 relies on.
    !|
    !| `x_star` is an ordinary input here, not computed by this routine -- 125's own
    !| `js_comp_test_helper` takes it the same way, since a caller running several studies/several
    !| parameter settings is expected to compute the reference points once
    !| ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) and reuse
    !| them consistently.
    !|
    !| `construct_neighborhoods_ranged_impl` reports neighbor gene INDICES, not gathered residual
    !| values (unlike its distance-sort sibling
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_impl(interface)]]),
    !| so this routine gathers each neighbor's actual residual values from `residuals` itself
    !| (`tmp_neighborhood_residuals_gathered`, a per-study scratch buffer) before calling
    !| `build_residual_histograms_impl`. `build_residual_histograms_impl`/`calc_pmf_impl` are
    !| POINT-major (`(n_points, n_bins)`), while `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts` here
    !| are BIN-major (`(n_bins, n_points, n_studies)`) to match
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]'s own
    !| convention -- every call across that boundary bridges with an explicit `transpose`, exactly
    !| as [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] and
    !| [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]] already do.
    !|
    !| Impure: calls the impure `gjct_permutation_test_impl`. A GSL failure it reports is folded
    !| into `ierr` (first failure only), matching that routine's own tolerant precedent.
    subroutine run_js_comp_test_expert_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            n_points,&
            n_neighbors,&
            n_bins,&
            shared_residual_range,&
            gene_means,&
            gene_means_perms,&
            residuals,&
            x_star,&
            neighborhood_indices,&
            neighborhood_range,&
            pmfs,&
            counts,&
            included_n_reps,&
            mean_pmf,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            js_divergences,&
            weights,&
            global_js_divergence,&
            p_values,&
            tmp_neighborhood_residuals_gathered,&
            tmp_counts_point_major,&
            tmp_pmf_point_major,&
            tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts,&
            tmp_permutation_pmfs,&
            tmp_permutation_js_divergences,&
            tmp_permutation_weights,&
            tmp_permutation_global_js_divergence,&
            n_permutations,&
            random_seed,&
            ierr&
        ) bind(C, name="run_js_comp_test_expert_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test_expert

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points for neighborhoods
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors per neighborhood
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means_perms
            !! Per-study sorting permutation for `gene_means` (ascending, NaN last)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        real(c_double), dimension(n_points), intent(in), target :: x_star
            !! Mean-expression reference points
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_neighbors, n_points, n_studies), intent(out), target :: neighborhood_indices
            !! Gene indices of the selected neighborhood, per reference point, per study
        integer(c_int), dimension(2, n_points, n_studies), intent(out), target :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span, as
            !! produced by construct_neighborhoods_ranged_impl
        real(c_double), dimension(n_bins, n_points, n_studies), intent(out), target :: pmfs
            !! `counts` normalized to `0 <= pmfs(:, :, i) <= 1` and `sum(pmfs(:, j, i)) == 1`
        integer(c_int), dimension(n_bins, n_points, n_studies), intent(out), target :: counts
            !! Absolute counts of a residual per bin for `pmfs`
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study
        real(c_double), dimension(n_bins, n_points), intent(out), target :: mean_pmf
            !! The consensus pmf, from create_mean_pmf_impl
        integer(c_int), dimension(n_bins, n_points), intent(out), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study from gjct_permutation_test_impl
        real(c_double), dimension(max_n_reps_all_studies, n_neighbors, n_points), intent(out), target :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residual values, reused per study
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: tmp_counts_point_major
            !! Working array: one study's point-major histogram counts from build_residual_histograms_impl
        real(c_double), dimension(n_points, n_bins), intent(out), target :: tmp_pmf_point_major
            !! Working array: one study's point-major pmf, reused both for
            !! build_residual_histograms_impl's output and for calc_pmf_impl's re-derived pmf
        integer(c_int), dimension(n_bins, n_points), intent(out), target :: tmp_permutation_mean_pmf_counts
            !! Working array forwarded to gjct_permutation_test_impl's own resampling pool
        integer(c_int), dimension(n_bins, n_points), intent(out), target :: tmp_permutation_counts
            !! Working array forwarded to gjct_permutation_test_impl's own per-study resampled counts
        real(c_double), dimension(n_bins, n_points, n_studies), intent(out), target :: tmp_permutation_pmfs
            !! Working array forwarded to gjct_permutation_test_impl's own resampled pmfs
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_permutation_js_divergences
            !! Working array forwarded to gjct_permutation_test_impl's own per-point JSD values
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_permutation_weights
            !! Working array forwarded to gjct_permutation_test_impl's own per-point weights
        real(c_double), dimension(n_studies), intent(out), target :: tmp_permutation_global_js_divergence
            !! Working array forwarded to gjct_permutation_test_impl's own resampled global JSD values
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations, forwarded to gjct_permutation_test_impl
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator for
            !! the permutation test

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(gene_means_perms, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_indices, n_neighbors * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(neighborhood_range, 2 * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_residuals_gathered, max_n_reps_all_studies * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_counts_point_major, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_pmf_point_major, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_mean_pmf_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_pmfs, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_global_js_divergence, n_studies)

        call run_js_comp_test_expert(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            gene_means = gene_means,&
            gene_means_perms = gene_means_perms,&
            residuals = residuals,&
            x_star = x_star,&
            neighborhood_indices = neighborhood_indices,&
            neighborhood_range = neighborhood_range,&
            pmfs = pmfs,&
            counts = counts,&
            included_n_reps = included_n_reps,&
            mean_pmf = mean_pmf,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            js_divergences = js_divergences,&
            weights = weights,&
            global_js_divergence = global_js_divergence,&
            p_values = p_values,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_counts_point_major = tmp_counts_point_major,&
            tmp_pmf_point_major = tmp_pmf_point_major,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_permutation_global_js_divergence = tmp_permutation_global_js_divergence,&
            n_permutations = n_permutations,&
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_parameter_search(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `determine_js_comp_test_n_points_n_neighbors_helper` and
    !| `_alloc`, merged into one implementation now that the new `_impl` rules leave no separate
    !| hand-written allocation layer. Pools all studies' residuals and gene means, sorts them once
    !| ([[f42_sort_impl(module):sort_real_heapsort_expl_size(interface)]]), generates the candidate
    !| grid
    !| ([[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]),
    !| then walks it from finest to coarsest resolution: for each candidate, builds every study's
    !| neighborhoods and checks the first admissibility gate
    !| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]);
    !| once every study passes, pools the consensus pmf and checks the second gate
    !| ([[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]);
    !| once that passes too, seeds a confidence interval with the observed JSD, bootstraps it
    !| ([[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]), and
    !| tests it against the running best candidate for a plateau
    !| ([[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]).
    !| `plateau_mode` picks which of that CI-overlap criterion and Issue #178's complementary
    !| relative-effect-size one
    !| ([[tox_data_integration_js_comp_test_impl(module):check_effect_size_plateau_condition_impl(interface)]])
    !| governs the stop condition; both are always computed and traced (`trace_*` below) once a
    !| candidate is admissible, regardless of `plateau_mode`, so a caller can compare what either
    !| criterion would have decided. The search stops (`exit`) the moment the SELECTED criterion's
    !| plateau is found -- see `plateau_mode`'s own mode table below for the accepted values.
    !|
    !| Issue #178 also names 3 blocking dependencies for validating the effect-size thresholds
    !| empirically -- the KX_FACTORS default, the Freedman-Diaconis bin-count overestimate, and the
    !| one-sided-vs-symmetric JSD formula question -- all deliberately left as-is here; see the
    !| project's JSD-Comp-Test follow-up issue. `delta_median_threshold`/`delta_max_threshold`
    !| default to the issue's own suggested (not yet validated) 0.05/0.10.
    !|
    !| When `plateau_mode` selects the effect-size criterion (`MODE_PLATEAU_EFFECT_SIZE` or
    !| `MODE_PLATEAU_BOTH`) and it plateaus independently of the CI-overlap criterion's own running
    !| "best candidate" bookkeeping, `best_candidate_index`/`best_candidate_pair_confidence_interval`
    !| are overridden to the candidate that actually triggered the effect-size plateau, so the
    !| candidate this routine returns is always the one that stopped the search.
    !|
    !| Ported verbatim, including the
    !| fallback 125 relies on: if no candidate ever plateaus, the search falls back to the FIRST
    !| (finest-resolution) candidate and resets `best_candidate_pair_confidence_interval` to
    !| `-1.0`; if the grid collapsed to a single candidate (see
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]'s
    !| own small-N collapse note), that one candidate is used regardless of whether it plateaued or
    !| even passed either gate -- the plateau machinery is bypassed entirely, exactly as 125 does.
    !|
    !| Per the plan's work-array translation for this routine specifically: `max_n_bins_all_candidates`
    !| (data-dependent, not cheaply closed-form in 125) is replaced by the fixed
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] ceiling, so every
    !| bin-dimensioned work array below is sized to MAX_N_BINS and sliced `(1:n_bins, ...)` per
    !| candidate, rather than carrying a separate recommend-sized dimension argument for it.
    !| `residuals`/`gene_means` are passed to
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]/[[f42_sort_impl(module):sort_real_heapsort_expl_size(interface)]]
    !| as their own multi-dimensional selves -- both callees declare their matching dummy with an
    !| explicit shape, so standard Fortran sequence association reinterprets the contiguous actual
    !| argument as the flat 1-D array they expect, exactly as 125's own `_alloc` layer did for the
    !| same calls.
    !|
    !| Impure: calls the impure
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]. A GSL
    !| failure it reports is folded into `ierr` (first failure only) without aborting the search,
    !| matching that routine's own tolerant precedent.
    subroutine run_js_comp_test_parameter_search_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            gene_means,&
            residuals,&
            shared_residual_range,&
            n_bootstraps,&
            join_method,&
            n_points,&
            n_neighbors,&
            n_bins,&
            best_candidate_pair_confidence_interval,&
            n_admissible_evaluated,&
            trace_n_points,&
            trace_n_neighbors,&
            trace_global_js_divergence,&
            trace_ci_lower,&
            trace_ci_upper,&
            trace_delta,&
            trace_delta_median,&
            trace_delta_max,&
            min_count_per_mean_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            ierr&
        ) bind(C, name="run_js_comp_test_parameter_search_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test_parameter_search
        use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps to perform for a candidate pair
            !! The minimum valid value is `1_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
            !! The way to evaluate all studies' confidence-interval overlaps for the plateau
            !! condition, forwarded to check_plateau_condition_impl
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(c_int), intent(out), target :: n_points
            !! The finally chosen candidate's `n_points`
        integer(c_int), intent(out), target :: n_neighbors
            !! The finally chosen candidate's `n_neighbors`
        integer(c_int), intent(out), target :: n_bins
            !! The finally chosen candidate's bin count
        real(c_double), dimension(2, n_studies), intent(out), target :: best_candidate_pair_confidence_interval
            !! The bootstrapped JSD confidence interval for the finally chosen candidate pair;
            !! `-1.0_real64` throughout if no candidate pair passed both admissibility gates and
            !! the search fell back to the finest-resolution candidate
        integer(c_int), intent(out), target :: n_admissible_evaluated
            !! Number of candidates that passed both admissibility gates and got a JSD/confidence
            !! interval computed before the search stopped (by plateau or grid exhaustion) -- the
            !! number of leading, valid columns/elements in every `trace_*` array below. This is
            !! Issue #178's own index `t` domain: "position in the ordered sequence of ADMISSIBLE
            !! parameter pairs" -- a candidate that failed either gate has no `trace_*` entry at
            !! all, rather than a zero-filled one
        integer(c_int), dimension(16), intent(out), target :: trace_n_points
            !! Per-admissible-candidate `n_points`, one entry per column of the other `trace_*`
            !! arrays. `16` = MAX_CANDIDATE_PAIRS, written as a literal for the same reason
            !! candidates_n_points_n_neighbors/n_bins_candidates are in
            !! generate_js_comp_test_candidates_impl
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`, paired with trace_n_points above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! (`L_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! (`U_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`), from check_effect_size_plateau_condition_impl;
            !! `-1.0_real64` throughout at the first admissible candidate specifically (no
            !! predecessor to diff against) -- every other column within `1:n_admissible_evaluated`
            !! holds a real value
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies (Delta-tilde_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies (Delta^max_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), intent(in), target :: min_count_per_mean_bin
            !! Minimum count each bin of the consensus pmf must reach to pass the second
            !! admissibility gate
            !! The minimum valid value is `0_int32`.
            !! The default value is `5_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have to pass the first
            !! admissibility gate
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Minimum fractional overlap a candidate's confidence interval must have with the
            !! running best, per `join_method`, to plateau
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        character(len=1, kind=c_char), dimension(19), intent(in), target :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode                                      | Value                                                                                 |
            !! |-------------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only (pre-Issue-#178 behavior) | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only       | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                          | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `'plateau_ci_overlap'`.
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound the median relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound the largest relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous admissible
            !! candidate's JSD was zero, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare an effect-size
            !! plateau, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
            !! to bootstrap_histogram_impl itself
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; folds any GSL allocation failure bootstrap_histogram_impl reports
        integer(int32) :: join_method_mode_f
        integer(int32) :: plateau_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_admissible_evaluated)
        M_CHECK_NON_NULL(min_count_per_mean_bin)
        M_CHECK_NON_NULL(min_neighbor_overlap)
        M_CHECK_NON_NULL(succeeding_ci_overlap)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)
        M_CHECK_ARRAY_NON_NULL(best_candidate_pair_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(trace_n_points, 16)
        M_CHECK_ARRAY_NON_NULL(trace_n_neighbors, 16)
        M_CHECK_ARRAY_NON_NULL(trace_global_js_divergence, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_lower, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_upper, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta_median, 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta_max, 16)
        M_CHECK_ARRAY_NON_NULL(plateau_mode, 19)

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
        block
            character(len=:), pointer :: plateau_mode_f
            plateau_mode_f => c_char_as_view(plateau_mode)

            select case (plateau_mode_f)
                case ("plateau_ci_overlap")
                    plateau_mode_mode_f = MODE_PLATEAU_CI_OVERLAP
                case ("plateau_effect_size")
                    plateau_mode_mode_f = MODE_PLATEAU_EFFECT_SIZE
                case ("plateau_both")
                    plateau_mode_mode_f = MODE_PLATEAU_BOTH
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call run_js_comp_test_parameter_search(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            shared_residual_range = shared_residual_range,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method_mode_f,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            n_bins = n_bins,&
            best_candidate_pair_confidence_interval = best_candidate_pair_confidence_interval,&
            n_admissible_evaluated = n_admissible_evaluated,&
            trace_n_points = trace_n_points,&
            trace_n_neighbors = trace_n_neighbors,&
            trace_global_js_divergence = trace_global_js_divergence,&
            trace_ci_lower = trace_ci_lower,&
            trace_ci_upper = trace_ci_upper,&
            trace_delta = trace_delta,&
            trace_delta_median = trace_delta_median,&
            trace_delta_max = trace_delta_max,&
            min_count_per_mean_bin = min_count_per_mean_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode_mode_f,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_parameter_search_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_parameter_search_expert(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `determine_js_comp_test_n_points_n_neighbors_helper` and
    !| `_alloc`, merged into one implementation now that the new `_impl` rules leave no separate
    !| hand-written allocation layer. Pools all studies' residuals and gene means, sorts them once
    !| ([[f42_sort_impl(module):sort_real_heapsort_expl_size(interface)]]), generates the candidate
    !| grid
    !| ([[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]),
    !| then walks it from finest to coarsest resolution: for each candidate, builds every study's
    !| neighborhoods and checks the first admissibility gate
    !| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]);
    !| once every study passes, pools the consensus pmf and checks the second gate
    !| ([[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]);
    !| once that passes too, seeds a confidence interval with the observed JSD, bootstraps it
    !| ([[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]), and
    !| tests it against the running best candidate for a plateau
    !| ([[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]).
    !| `plateau_mode` picks which of that CI-overlap criterion and Issue #178's complementary
    !| relative-effect-size one
    !| ([[tox_data_integration_js_comp_test_impl(module):check_effect_size_plateau_condition_impl(interface)]])
    !| governs the stop condition; both are always computed and traced (`trace_*` below) once a
    !| candidate is admissible, regardless of `plateau_mode`, so a caller can compare what either
    !| criterion would have decided. The search stops (`exit`) the moment the SELECTED criterion's
    !| plateau is found -- see `plateau_mode`'s own mode table below for the accepted values.
    !|
    !| Issue #178 also names 3 blocking dependencies for validating the effect-size thresholds
    !| empirically -- the KX_FACTORS default, the Freedman-Diaconis bin-count overestimate, and the
    !| one-sided-vs-symmetric JSD formula question -- all deliberately left as-is here; see the
    !| project's JSD-Comp-Test follow-up issue. `delta_median_threshold`/`delta_max_threshold`
    !| default to the issue's own suggested (not yet validated) 0.05/0.10.
    !|
    !| When `plateau_mode` selects the effect-size criterion (`MODE_PLATEAU_EFFECT_SIZE` or
    !| `MODE_PLATEAU_BOTH`) and it plateaus independently of the CI-overlap criterion's own running
    !| "best candidate" bookkeeping, `best_candidate_index`/`best_candidate_pair_confidence_interval`
    !| are overridden to the candidate that actually triggered the effect-size plateau, so the
    !| candidate this routine returns is always the one that stopped the search.
    !|
    !| Ported verbatim, including the
    !| fallback 125 relies on: if no candidate ever plateaus, the search falls back to the FIRST
    !| (finest-resolution) candidate and resets `best_candidate_pair_confidence_interval` to
    !| `-1.0`; if the grid collapsed to a single candidate (see
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]'s
    !| own small-N collapse note), that one candidate is used regardless of whether it plateaued or
    !| even passed either gate -- the plateau machinery is bypassed entirely, exactly as 125 does.
    !|
    !| Per the plan's work-array translation for this routine specifically: `max_n_bins_all_candidates`
    !| (data-dependent, not cheaply closed-form in 125) is replaced by the fixed
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] ceiling, so every
    !| bin-dimensioned work array below is sized to MAX_N_BINS and sliced `(1:n_bins, ...)` per
    !| candidate, rather than carrying a separate recommend-sized dimension argument for it.
    !| `residuals`/`gene_means` are passed to
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]/[[f42_sort_impl(module):sort_real_heapsort_expl_size(interface)]]
    !| as their own multi-dimensional selves -- both callees declare their matching dummy with an
    !| explicit shape, so standard Fortran sequence association reinterprets the contiguous actual
    !| argument as the flat 1-D array they expect, exactly as 125's own `_alloc` layer did for the
    !| same calls.
    !|
    !| Impure: calls the impure
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]. A GSL
    !| failure it reports is folded into `ierr` (first failure only) without aborting the search,
    !| matching that routine's own tolerant precedent.
    subroutine run_js_comp_test_parameter_search_expert_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            gene_means,&
            residuals,&
            shared_residual_range,&
            n_bootstraps,&
            join_method,&
            max_n_points_candidate,&
            max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds,&
            n_points,&
            n_neighbors,&
            n_bins,&
            best_candidate_pair_confidence_interval,&
            n_admissible_evaluated,&
            trace_n_points,&
            trace_n_neighbors,&
            trace_global_js_divergence,&
            trace_ci_lower,&
            trace_ci_upper,&
            trace_delta,&
            trace_delta_median,&
            trace_delta_max,&
            tmp_gene_means_perms,&
            tmp_gene_means_perm_all,&
            tmp_residuals_perm,&
            tmp_x_star,&
            tmp_neighborhood_indices,&
            tmp_neighborhood_range,&
            tmp_neighborhood_residuals_gathered,&
            tmp_counts_point_major,&
            tmp_pmf_point_major,&
            tmp_pmfs,&
            tmp_counts,&
            tmp_included_n_reps,&
            tmp_mean_pmf,&
            tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences,&
            tmp_weights,&
            tmp_global_js_divergence,&
            tmp_confidence_interval,&
            tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence,&
            tmp_delta_perm,&
            min_count_per_mean_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            ierr&
        ) bind(C, name="run_js_comp_test_parameter_search_expert_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test_parameter_search_expert
        use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_points_candidate
            !! Exact upper bound on the grid's first (largest) `n_points` candidate
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_neighbors_candidate
            !! Safe upper bound on the grid's largest `n_neighbors` candidate
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstrapping_top_k_jsds
            !! Number of elements kept at each end of the bootstrap distribution (top-k/bottom-k
            !! heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps to perform for a candidate pair
            !! The minimum valid value is `1_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
            !! The way to evaluate all studies' confidence-interval overlaps for the plateau
            !! condition, forwarded to check_plateau_condition_impl
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(c_int), intent(out), target :: n_points
            !! The finally chosen candidate's `n_points`
        integer(c_int), intent(out), target :: n_neighbors
            !! The finally chosen candidate's `n_neighbors`
        integer(c_int), intent(out), target :: n_bins
            !! The finally chosen candidate's bin count
        real(c_double), dimension(2, n_studies), intent(out), target :: best_candidate_pair_confidence_interval
            !! The bootstrapped JSD confidence interval for the finally chosen candidate pair;
            !! `-1.0_real64` throughout if no candidate pair passed both admissibility gates and
            !! the search fell back to the finest-resolution candidate
        integer(c_int), intent(out), target :: n_admissible_evaluated
            !! Number of candidates that passed both admissibility gates and got a JSD/confidence
            !! interval computed before the search stopped (by plateau or grid exhaustion) -- the
            !! number of leading, valid columns/elements in every `trace_*` array below. This is
            !! Issue #178's own index `t` domain: "position in the ordered sequence of ADMISSIBLE
            !! parameter pairs" -- a candidate that failed either gate has no `trace_*` entry at
            !! all, rather than a zero-filled one
        integer(c_int), dimension(16), intent(out), target :: trace_n_points
            !! Per-admissible-candidate `n_points`, one entry per column of the other `trace_*`
            !! arrays. `16` = MAX_CANDIDATE_PAIRS, written as a literal for the same reason
            !! candidates_n_points_n_neighbors/n_bins_candidates are in
            !! generate_js_comp_test_candidates_impl
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`, paired with trace_n_points above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! (`L_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! (`U_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`), from check_effect_size_plateau_condition_impl;
            !! `-1.0_real64` throughout at the first admissible candidate specifically (no
            !! predecessor to diff against) -- every other column within `1:n_admissible_evaluated`
            !! holds a real value
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies (Delta-tilde_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies (Delta^max_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(out), target :: tmp_gene_means_perms
            !! Working array: each study's own sorting permutation for `gene_means`
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_gene_means_perm_all
            !! Working array: sorting permutation for the flattened, all-studies-pooled `gene_means`
        integer(c_int), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out), target :: tmp_residuals_perm
            !! Working array: sorting permutation for the flattened, all-studies-pooled `residuals`
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_x_star
            !! Working array: reference points for the candidate whose `n_points` is current,
            !! recomputed only when `n_points` changes between candidates
        integer(c_int), dimension(max_n_neighbors_candidate, max_n_points_candidate), intent(out), target :: tmp_neighborhood_indices
            !! Working array: one study's neighbor gene indices for the current candidate, reused
            !! per study
        integer(c_int), dimension(2, max_n_points_candidate), intent(out), target :: tmp_neighborhood_range
            !! Working array: one study's `[min_idx, max_idx]` neighborhood spans for the current
            !! candidate, reused per study
        real(c_double), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), intent(out), target :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residual values for the current
            !! candidate, reused per study
        integer(c_int), dimension(max_n_points_candidate, 256), intent(out), target :: tmp_counts_point_major
            !! Working array: one study's point-major histogram counts for the current candidate,
            !! reused per study. The `256` is
            !! [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]], written as a
            !! literal because a dimension naming a module parameter the generated wrapper never
            !! `use`s would not compile there -- the same reason MAX_CANDIDATE_PAIRS is written as
            !! `16` on generate_js_comp_test_candidates_impl's own dummy arguments
        real(c_double), dimension(max_n_points_candidate, 256), intent(out), target :: tmp_pmf_point_major
            !! Working array: one study's point-major pmf for the current candidate, reused per
            !! study. `256` = MAX_N_BINS, see tmp_counts_point_major above
        real(c_double), dimension(256, max_n_points_candidate, n_studies), intent(out), target :: tmp_pmfs
            !! Working array: every study's bin-major pmf for the current candidate. `256` =
            !! MAX_N_BINS, see tmp_counts_point_major above
        integer(c_int), dimension(256, max_n_points_candidate, n_studies), intent(out), target :: tmp_counts
            !! Working array: every study's bin-major histogram counts for the current candidate.
            !! `256` = MAX_N_BINS, see tmp_counts_point_major above
        integer(c_int), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_included_n_reps
            !! Working array: every study's included-replicate counts for the current candidate
        real(c_double), dimension(256, max_n_points_candidate), intent(out), target :: tmp_mean_pmf
            !! Working array: the current candidate's consensus pmf. `256` = MAX_N_BINS, see
            !! tmp_counts_point_major above
        integer(c_int), dimension(256, max_n_points_candidate), intent(out), target :: tmp_mean_pmf_counts
            !! Working array: the current candidate's consensus histogram counts. `256` =
            !! MAX_N_BINS, see tmp_counts_point_major above
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_mean_pmf_included_n_reps
            !! Working array: the current candidate's consensus included-replicate counts
        real(c_double), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_js_divergences
            !! Working array: the current candidate's per-reference-point JSD values, reused as
            !! bootstrap_histogram_impl's own scratch once consumed
        real(c_double), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_weights
            !! Working array: the current candidate's per-reference-point weights, reused as
            !! bootstrap_histogram_impl's own scratch once consumed
        real(c_double), dimension(n_studies), intent(out), target :: tmp_global_js_divergence
            !! Working array: the current candidate's observed global JSD per study, reused as
            !! bootstrap_histogram_impl's own scratch once consumed
        real(c_double), dimension(2, n_studies), intent(out), target :: tmp_confidence_interval
            !! Working array: the current candidate's confidence interval, seeded with the observed
            !! global JSD and then bootstrapped in place
        real(c_double), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out), target :: tmp_bootstrapping_top_k_jsds
            !! Working array forwarded to bootstrap_histogram_impl's own top-k/bottom-k heaps
        real(c_double), dimension(n_studies), intent(out), target :: tmp_prev_global_js_divergence
            !! Working array: the previous admissible candidate's observed global JSD per study,
            !! forwarded to check_effect_size_plateau_condition_impl
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_delta_perm
            !! Working array forwarded to check_effect_size_plateau_condition_impl's own median
            !! computation
        integer(c_int), intent(in), target :: min_count_per_mean_bin
            !! Minimum count each bin of the consensus pmf must reach to pass the second
            !! admissibility gate
            !! The minimum valid value is `0_int32`.
            !! The default value is `5_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have to pass the first
            !! admissibility gate
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Minimum fractional overlap a candidate's confidence interval must have with the
            !! running best, per `join_method`, to plateau
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        character(len=1, kind=c_char), dimension(19), intent(in), target :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode                                      | Value                                                                                 |
            !! |-------------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only (pre-Issue-#178 behavior) | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only       | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                          | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `'plateau_ci_overlap'`.
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound the median relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound the largest relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous admissible
            !! candidate's JSD was zero, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare an effect-size
            !! plateau, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
            !! to bootstrap_histogram_impl itself
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; folds any GSL allocation failure bootstrap_histogram_impl reports
        integer(int32) :: join_method_mode_f
        integer(int32) :: plateau_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(max_n_neighbors_candidate)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_admissible_evaluated)
        M_CHECK_NON_NULL(min_count_per_mean_bin)
        M_CHECK_NON_NULL(min_neighbor_overlap)
        M_CHECK_NON_NULL(succeeding_ci_overlap)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)
        M_CHECK_ARRAY_NON_NULL(best_candidate_pair_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(trace_n_points, 16)
        M_CHECK_ARRAY_NON_NULL(trace_n_neighbors, 16)
        M_CHECK_ARRAY_NON_NULL(trace_global_js_divergence, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_lower, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_upper, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta_median, 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta_max, 16)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perms, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perm_all, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_residuals_perm, (max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_x_star, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_indices, max_n_neighbors_candidate * max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_range, 2 * max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_residuals_gathered, max_n_reps_all_studies * max_n_neighbors_candidate * max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_counts_point_major, max_n_points_candidate * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_pmf_point_major, max_n_points_candidate * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, 256 * max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, 256 * max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_included_n_reps, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf, 256 * max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf_counts, 256 * max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf_included_n_reps, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds * 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_prev_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_delta_perm, n_studies)
        M_CHECK_ARRAY_NON_NULL(plateau_mode, 19)

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
        block
            character(len=:), pointer :: plateau_mode_f
            plateau_mode_f => c_char_as_view(plateau_mode)

            select case (plateau_mode_f)
                case ("plateau_ci_overlap")
                    plateau_mode_mode_f = MODE_PLATEAU_CI_OVERLAP
                case ("plateau_effect_size")
                    plateau_mode_mode_f = MODE_PLATEAU_EFFECT_SIZE
                case ("plateau_both")
                    plateau_mode_mode_f = MODE_PLATEAU_BOTH
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call run_js_comp_test_parameter_search_expert(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            shared_residual_range = shared_residual_range,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method_mode_f,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            n_bins = n_bins,&
            best_candidate_pair_confidence_interval = best_candidate_pair_confidence_interval,&
            n_admissible_evaluated = n_admissible_evaluated,&
            trace_n_points = trace_n_points,&
            trace_n_neighbors = trace_n_neighbors,&
            trace_global_js_divergence = trace_global_js_divergence,&
            trace_ci_lower = trace_ci_lower,&
            trace_ci_upper = trace_ci_upper,&
            trace_delta = trace_delta,&
            trace_delta_median = trace_delta_median,&
            trace_delta_max = trace_delta_max,&
            tmp_gene_means_perms = tmp_gene_means_perms,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_residuals_perm = tmp_residuals_perm,&
            tmp_x_star = tmp_x_star,&
            tmp_neighborhood_indices = tmp_neighborhood_indices,&
            tmp_neighborhood_range = tmp_neighborhood_range,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_counts_point_major = tmp_counts_point_major,&
            tmp_pmf_point_major = tmp_pmf_point_major,&
            tmp_pmfs = tmp_pmfs,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps = tmp_included_n_reps,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps = tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_confidence_interval = tmp_confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence = tmp_prev_global_js_divergence,&
            tmp_delta_perm = tmp_delta_perm,&
            min_count_per_mean_bin = min_count_per_mean_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode_mode_f,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_parameter_search_expert_c

end module tox_data_integration_js_comp_test_c
#endif
