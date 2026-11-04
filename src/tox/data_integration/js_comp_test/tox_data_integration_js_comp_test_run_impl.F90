#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Final Run
!|
!| The JS-Comp-Test of K studies for neighborhood parameters already chosen: per-study
!| histograms at each reference point's bin count, each study's JSD (by default, for two studies,
!| their pairwise JSD; otherwise against the consensus pmf, see `divergence_mode`) and its weighting
!| into one global JSD per study. Its significance test is the study-label Permutation-Test,
!| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_permutation_impl(interface)]]
!| and its adaptive counterpart, which rerun these routines on reassigned samples.
!| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]] runs it
!| on fixed-k neighborhoods (`n_points` reference points of `n_neighbors` genes each),
!| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_adaptive_impl(interface)]]
!| on adaptive neighborhoods built by
!| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]].
!| The adaptive neighborhoods, chosen by
!| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]],
!| are the recommended way; the fixed-k ones stay available.
module tox_data_integration_js_comp_test_run_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT
    use f42_stats_impl, only: last_non_nan_position
    use tox_data_integration_jsd_impl, only: compute_consensus_divergences, create_mean_pmf_impl, &
                                             has_study_without_residuals
    use tox_data_integration_preprocessing_impl, only: construct_neighborhoods_ranged_impl, &
                                                       materialize_pooled_neighborhood, sort_gene_means_per_study, &
                                                       sort_pooled_gene_means
    use tox_data_integration_js_comp_test_binning_impl, only: MAX_N_BINS, determine_fixed_k_point_bin_counts, &
                                                              determine_adaptive_point_bin_counts
    use tox_data_integration_js_comp_test_histograms_impl, only: build_fixed_k_histograms, build_adaptive_histograms
    M_IMPLICIT_NONE
    private
    public :: run_js_comp_test_impl, run_js_comp_test_adaptive_impl, prepare_gene_means_perms

contains

    !> The shared tail of both final runs, after every study's histograms are built: the consensus
    !| pmf (create_mean_pmf_impl), each study's per-point JSD and its weighted global JSD
    !| (compute_consensus_divergences, with `divergence_mode`). Private.
    !|
    !| The consensus is built over the full `MAX_N_BINS` extent, so `pmfs`/`counts` must be zero
    !| beyond `max_n_bins_per_point`; `mean_pmf`/`mean_pmf_counts` are then zero there too.
    pure subroutine finish_js_comp_test(n_studies, n_points, max_n_bins_per_point, pmfs, counts, included_n_reps, &
                                        divergence_mode, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps, &
                                        js_divergences, weights, global_js_divergence)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: max_n_bins_per_point
            !! Widest bin count of any reference point; only bins `1:max_n_bins_per_point` of the
            !! bin-sized arrays below are meaningful
        real(real64), dimension(n_points, MAX_N_BINS, n_studies), intent(in) :: pmfs
            !! Per-study pmfs, zero beyond `max_n_bins_per_point`
        integer(int32), dimension(n_points, MAX_N_BINS, n_studies), intent(in) :: counts
            !! Per-study histogram counts, zero beyond `max_n_bins_per_point`
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        integer(int32), intent(in) :: divergence_mode
            !! One of the `MODE_DIVERGENCE_*` values, already resolved from its default by the caller
        real(real64), dimension(n_points, MAX_N_BINS), intent(out) :: mean_pmf
            !! The consensus pmf
        integer(int32), dimension(n_points, MAX_N_BINS), intent(out) :: mean_pmf_counts
            !! Absolute counts per bin for the consensus pmf
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point for the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study, the one `divergence_mode` selects
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study, the one `divergence_mode` selects

        ! Over the full extent: the studies' bins 1:max_n_bins_per_point are not one contiguous
        ! block, and the zero padding sums to zero.
        call create_mean_pmf_impl(pmfs, counts, MAX_N_BINS, n_points, n_studies, included_n_reps, mean_pmf, &
                                  mean_pmf_included_n_reps, mean_pmf_counts)

        call compute_consensus_divergences(pmfs, mean_pmf, n_points, max_n_bins_per_point, MAX_N_BINS, n_studies, &
                                           included_n_reps, mean_pmf_included_n_reps, divergence_mode, js_divergences, &
                                           global_js_divergence, weights)
    end subroutine finish_js_comp_test

    !> summary: Sort each study's gene means: one ascending, NaN-last permutation per study
    !| AUTHOR_LASZLO_LANG
    !| The prologue of the plain `run_js_comp_test`: column `i_study` of `gene_means_perms` is the
    !| heapsort permutation of `gene_means(:, i_study)`, built by
    !| [[tox_data_integration_preprocessing_impl(module):sort_gene_means_per_study(subroutine)]] exactly
    !| as the fixed-k parameter search builds its own. Never handles the call itself (`handled` is always `.false.`).
    pure subroutine prepare_gene_means_perms(gene_means, max_n_genes_all_studies, n_studies, gene_means_perms, &
                                             handled, ierr)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: gene_means_perms
            !! Per study (column), the permutation that sorts that study's `gene_means` ascending, NaN last
        ! Plain logical, not c_bool: the generated wrapper declares `handled` default-kind.
        logical, intent(out) :: handled
            !! Always `.false.`: the implementation always runs
        integer(int32), intent(out) :: ierr
            !! Error code; always success

        call sort_gene_means_per_study(gene_means, max_n_genes_all_studies, n_studies, gene_means_perms)
        handled = .false.
        call set_ok(ierr)
    end subroutine prepare_gene_means_perms

    !> summary: Run the JSD-Comp-Test pipeline for one fixed (n_points, n_neighbors) parameter setting
    !| AUTHOR_LASZLO_LANG
    !| The final JS-Comp-Test of K studies for one fixed `(n_points, n_neighbors)` setting:
    !|
    !| 1. Each study's `n_neighbors` genes nearest every reference point `x_star` (computed once by
    !|    the caller, e.g. with
    !|    [[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) form its
    !|    neighborhoods
    !|    ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]).
    !| 2. Per reference point, all studies' residuals are pooled and the bin-count search
    !|    `bin_count_mode` names chooses the point's bin count and residual range (occupancy search
    !|    [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !|    or Stone's rule
    !|    [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]],
    !|    by default with every study's score counted equally).
    !| 3. Each study's histogram and pmf is built per point
    !|    ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]).
    !| 4. The consensus pmf
    !|    ([[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]]), each study's
    !|    per-point JSD and its weighted global JSD: by default (`divergence_mode`), for exactly two
    !|    studies, the pairwise JSD of the two studies, the same for both; otherwise each study's JSD
    !|    against the consensus.
    !|
    !| The p-values of the study-label Permutation-Test come from
    !| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_permutation_impl(interface)]].
    !|
    !| There is no admissibility gate: a point whose bin-count search fails (`occupancy_failed`)
    !| still gets a histogram at `m_min` bins and still contributes to `global_js_divergence`. A
    !| caller that needs such neighborhoods rejected should use
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]].
    !|
    !| The bin-sized arrays (`pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`) have the fixed
    !| second extent [[tox_data_integration_js_comp_test_binning_impl(module):MAX_N_BINS(variable)]]
    !| (`M_MAX_N_BINS`); only their first `max_n_bins_per_point` bins are meaningful, the rest is
    !| zero. A Python/R caller slices `[:, :max_n_bins_per_point, ...]` (R:
    !| `[, seq_len(max_n_bins_per_point), ]`) itself.
    !|
    !| Errors (`ierr`): `ERR_INVALID_INPUT` when some study has no non-NaN residual at some
    !| reference point (for example a study whose gene means are all missing), since the
    !| equal-weight consensus needs every study; the consensus, JSD and weight outputs are then not
    !| filled.
    !|
    !| The plain `run_js_comp_test` sorts each study's gene means itself (heapsort, ascending, NaN
    !| last) and takes no permutation; `run_js_comp_test_expert` takes `gene_means_perms` from the
    !| caller instead. That order must still be ascending with NaN last (only the order among tied
    !| genes is the caller's choice); it is checked neither to be a permutation of
    !| `1..max_n_genes_all_studies` per column nor to be sorted.
    !| DM_PROLOGUE(prepare_gene_means_perms, tox_data_integration_js_comp_test_run_impl)
    pure subroutine run_js_comp_test_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, n_neighbors, &
                                     gene_means, gene_means_perms, residuals, x_star, &
                                     neighborhood_indices, neighborhood_range, n_bins_per_point, &
                                     shared_residual_range_low, shared_residual_range_high, max_n_bins_per_point, &
                                     occupancy_failed, n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                     max_bin_occupancy, sturges_bins, fd_bins, pmfs, counts, included_n_reps, mean_pmf, &
                                     mean_pmf_counts, mean_pmf_included_n_reps, js_divergences, weights, &
                                     global_js_divergence, tmp_neighborhood_residuals_gathered, tmp_pooled_residuals, &
                                     tmp_pooled_residuals_perm, tmp_bin_counts_search, tmp_residual_to_study_perm, &
                                     tmp_study_block_ends, tmp_study_block_sizes, tmp_point_n_neighbors, &
                                     min_residuals_per_bin, m_min, m_max, gamma_occupancy, &
                                     lower_residual_range_quantile, upper_residual_range_quantile, bin_count_mode, &
                                     divergence_mode, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_points
            !! Number of reference points for neighborhoods
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies)
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
            !! DM_ALLOW_NAN
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means_perms
            !! Per study (column), a permutation of that study's `gene_means` in ascending order, NaN last.
            !! The plain entry point builds it for you (heapsort); the expert one takes yours, which
            !! must be sorted that way too (only the order among tied genes is free). It is checked
            !! neither to be a permutation of `1..max_n_genes_all_studies` per column nor to be sorted.
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies)
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
            !! DM_ALLOW_NAN
        real(real64), dimension(n_points), intent(in) :: x_star
            !! Mean-expression reference points
        integer(int32), dimension(n_neighbors, n_points, n_studies), intent(out) :: neighborhood_indices
            !! Gene indices of the selected neighborhood, per reference point, per study
        integer(int32), dimension(2, n_points, n_studies), intent(out) :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span
        integer(int32), dimension(n_points), intent(out) :: n_bins_per_point
            !! Each reference point's histogram bin count `M_j`, from the bin-count search
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(int32), intent(out) :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts` (see above)
        logical(c_bool), dimension(n_points), intent(out) :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes (see above)
        integer(int32), dimension(n_points), intent(out) :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(int32), dimension(n_points), intent(out) :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(real64), dimension(n_points), intent(out) :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(int32), dimension(n_points), intent(out) :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(int32), dimension(n_points), intent(out) :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(int32), dimension(n_points), intent(out) :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(real64), dimension(n_points, M_MAX_N_BINS, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, M_MAX_N_BINS, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(real64), dimension(n_points, M_MAX_N_BINS), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, M_MAX_N_BINS), intent(out) :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study, the one `divergence_mode` selects
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study, the one `divergence_mode` selects
        real(real64), dimension(max_n_reps_all_studies, n_neighbors, n_points), intent(out) :: &
            tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residuals, reused per study
        real(real64), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(int32), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array: the last index of each study's block in tmp_pooled_residuals
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in tmp_pooled_residuals
        integer(int32), dimension(n_studies), intent(out) :: tmp_point_n_neighbors
            !! Working array: every study's neighbor count at a reference point (all `n_neighbors`)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
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
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(M_MODE_BIN_COUNT_STONE_EQUAL)
        integer(int32), intent(in), optional :: divergence_mode
            !! How each study's per-point divergence is computed
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Each study's JSD against the consensus of all studies | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]] |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! DM_DEFAULT(M_MODE_DIVERGENCE_PAIRWISE_K2)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT (see above)

        integer(int32) :: i_study, actual_divergence_mode

        call set_ok(ierr)
        M_DEFAULT_VAL(divergence_mode, actual_divergence_mode, M_MODE_DIVERGENCE_PAIRWISE_K2)

        ! ===== PASS A (per study): neighborhood construction.
        do i_study = 1, n_studies
            call construct_neighborhoods_ranged_impl(n_points, x_star, max_n_genes_all_studies, gene_means(:, i_study), &
                                                     gene_means_perms(:, i_study), n_neighbors, &
                                                     neighborhood_indices(:, :, i_study), neighborhood_range(:, :, i_study))
        end do

        ! ===== PASS B (per point): pool all studies' residuals and run the bin-count search.
        ! Sequential, so one tmp_pooled_residuals buffer serves every point (do concurrent would
        ! need one per point).
        call determine_fixed_k_point_bin_counts(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                                n_neighbors, n_points, neighborhood_indices, n_bins_per_point, &
                                                occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                                n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                                max_bin_occupancy, sturges_bins, fd_bins, tmp_point_n_neighbors, &
                                                tmp_pooled_residuals, tmp_pooled_residuals_perm, tmp_bin_counts_search, &
                                                tmp_residual_to_study_perm, tmp_study_block_ends, tmp_study_block_sizes, &
                                                m_min=m_min, m_max=m_max, min_residuals_per_bin=min_residuals_per_bin, &
                                                gamma_occupancy=gamma_occupancy, &
                                                lower_residual_range_quantile=lower_residual_range_quantile, &
                                                upper_residual_range_quantile=upper_residual_range_quantile, &
                                                bin_count_mode=bin_count_mode)

        ! ===== PASS C (per study): build every (point, study) histogram at the chosen bin counts.
        ! The bins beyond max_n_bins_per_point are zeroed, as finish_js_comp_test requires.
        max_n_bins_per_point = maxval(n_bins_per_point(1:n_points))
        call build_fixed_k_histograms(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, n_neighbors, &
                                      n_points, max_n_bins_per_point, MAX_N_BINS, neighborhood_indices, &
                                      shared_residual_range_low, shared_residual_range_high, n_bins_per_point, counts, &
                                      pmfs, included_n_reps, tmp_neighborhood_residuals_gathered)
        counts(:, max_n_bins_per_point + 1:, :) = 0_int32
        pmfs(:, max_n_bins_per_point + 1:, :) = 0.0_real64

        ! A study with no included residual at a point would leave that point's equal-weight
        ! consensus summing to less than 1. It depends on gene_means and residuals together, so
        ! no arg_pos.
        if (has_study_without_residuals(included_n_reps, n_points, n_studies)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        call finish_js_comp_test(n_studies, n_points, max_n_bins_per_point, pmfs, counts, included_n_reps, &
                                 actual_divergence_mode, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps, &
                                 js_divergences, weights, global_js_divergence)
    end subroutine run_js_comp_test_impl

    !> summary: Run the JSD-Comp-Test pipeline on adaptive (Issue #217) neighborhoods
    !| AUTHOR_LASZLO_LANG
    !| The adaptive counterpart of
    !| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]]: the
    !| same pipeline and outputs, on neighborhoods whose size varies per reference point and per
    !| study. Call
    !| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
    !| first and pass its `pooled_neighborhood_range` (its first `n_points` columns) with the same
    !| `gene_means` and `residuals`. Each range is decoded back to per-study gene lists; their sizes
    !| are returned as `n_neighbors_per_point`. Neighborhoods that decode to the same gene sets as a
    !| fixed-k run give bit-identical results. There is no admissibility gate.
    !|
    !| Errors (`ierr`): `ERR_INVALID_INPUT` before any histogram work for the first range that is
    !| reversed, reaches past the non-NaN pooled means, or leaves a study without a gene (every
    !| output then undefined); otherwise as in `run_js_comp_test_impl` (a study without a non-NaN
    !| residual at some point).
    !|
    !| The pooling work arrays hold `max_n_reps_all_studies * max_n_genes_all_studies * n_studies`
    !| values each, the size of `residuals`, since one study may contribute all its genes to one
    !| neighborhood. The bin-sized arrays have the fixed second extent
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MAX_N_BINS(variable)]]; only their first
    !| `max_n_bins_per_point` bins are meaningful (Python `[:, :max_n_bins_per_point, ...]`, R
    !| `[, seq_len(max_n_bins_per_point), ]`).
    pure subroutine run_js_comp_test_adaptive_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, &
                                              gene_means, residuals, pooled_neighborhood_range, n_neighbors_per_point, &
                                              n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                              max_n_bins_per_point, occupancy_failed, n_pooled_residuals, &
                                              min_bin_occupancy, mean_bin_occupancy, max_bin_occupancy, sturges_bins, &
                                              fd_bins, pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, &
                                              mean_pmf_included_n_reps, js_divergences, weights, global_js_divergence, &
                                              tmp_gene_means_perm_all, tmp_point_neighborhood_indices, &
                                              tmp_neighbor_residuals, tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                              tmp_bin_counts_search, tmp_residual_to_study_perm, tmp_study_block_ends, &
                                              tmp_study_block_sizes, min_residuals_per_bin, m_min, m_max, &
                                              gamma_occupancy, lower_residual_range_quantile, &
                                              upper_residual_range_quantile, bin_count_mode, divergence_mode, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (neighborhoods), the construction's own `n_points`
            !! DM_MIN(1_int32)
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene; the array the
            !! neighborhoods were constructed from
            !! DM_ALLOW_NAN
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! DM_ALLOW_NAN
        integer(int32), dimension(2, n_points), intent(in) :: pooled_neighborhood_range
            !! For each neighborhood, its first and last sorted pooled position, from
            !! construct_adaptive_neighborhoods
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies*n_studies)
        integer(int32), dimension(n_studies, n_points), intent(out) :: n_neighbors_per_point
            !! For each neighborhood, how many of its genes belong to each study (at least 1 each)
        integer(int32), dimension(n_points), intent(out) :: n_bins_per_point
            !! Each reference point's selected histogram bin count, from the bin-count search
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(int32), intent(out) :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts`
        logical(c_bool), dimension(n_points), intent(out) :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes
        integer(int32), dimension(n_points), intent(out) :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(int32), dimension(n_points), intent(out) :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(real64), dimension(n_points), intent(out) :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(int32), dimension(n_points), intent(out) :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(int32), dimension(n_points), intent(out) :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(int32), dimension(n_points), intent(out) :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(real64), dimension(n_points, M_MAX_N_BINS, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, M_MAX_N_BINS, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN residuals binned per reference point, per study
        real(real64), dimension(n_points, M_MAX_N_BINS), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, M_MAX_N_BINS), intent(out) :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN residuals per reference point for the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study, the one `divergence_mode` selects
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study, the one `divergence_mode` selects
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, as in the construction
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_point_neighborhood_indices
            !! Working array: one reference point's gene indices per study, reused per point
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out) :: tmp_neighbor_residuals
            !! Working array: one (reference point, study) pair's gathered residuals, reused per pair
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array: sorting permutation for `tmp_pooled_residuals`, reused per point
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array: the last index of each study's block in `tmp_pooled_residuals`
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in `tmp_pooled_residuals`
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
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
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(M_MODE_BIN_COUNT_STONE_EQUAL)
        integer(int32), intent(in), optional :: divergence_mode
            !! How each study's per-point divergence is computed
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Each study's JSD against the consensus of all studies | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]] |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! DM_DEFAULT(M_MODE_DIVERGENCE_PAIRWISE_K2)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT (see above)

        integer(int32) :: i_point, i_study, n_valid_means, actual_divergence_mode

        call set_ok(ierr)
        M_DEFAULT_VAL(divergence_mode, actual_divergence_mode, M_MODE_DIVERGENCE_PAIRWISE_K2)

        ! The construction's own sort, so every pooled position names the same entry as there.
        call sort_pooled_gene_means(gene_means, max_n_genes_all_studies, n_studies, tmp_gene_means_perm_all)
        ! NaN sorts last, so the non-NaN means occupy exactly the first n_valid_means positions.
        n_valid_means = last_non_nan_position(gene_means, max_n_genes_all_studies*n_studies, tmp_gene_means_perm_all)

        ! Runtime checks, all before any histogram work. Sequential with an early return: the
        ! first failing point ends the routine.
        do i_point = 1, n_points
            if (pooled_neighborhood_range(1, i_point) > pooled_neighborhood_range(2, i_point) .or. &
                pooled_neighborhood_range(2, i_point) > n_valid_means) then
                call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
                return
            end if
            call materialize_pooled_neighborhood(tmp_gene_means_perm_all, max_n_genes_all_studies, n_studies, &
                                                 pooled_neighborhood_range(1, i_point), &
                                                 pooled_neighborhood_range(2, i_point), max_n_genes_all_studies, &
                                                 tmp_point_neighborhood_indices, n_neighbors_per_point(:, i_point))
            do i_study = 1, n_studies
                if (n_neighbors_per_point(i_study, i_point) == 0_int32) then
                    call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
                    return
                end if
            end do
        end do

        ! ===== PASS B (per point, sequential: one pooling buffer serves every point).
        call determine_adaptive_point_bin_counts(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                                 n_points, tmp_gene_means_perm_all, pooled_neighborhood_range, &
                                                 n_neighbors_per_point, n_bins_per_point, occupancy_failed, &
                                                 shared_residual_range_low, shared_residual_range_high, n_pooled_residuals, &
                                                 min_bin_occupancy, mean_bin_occupancy, max_bin_occupancy, sturges_bins, &
                                                 fd_bins, tmp_point_neighborhood_indices, tmp_pooled_residuals, &
                                                 tmp_pooled_residuals_perm, tmp_bin_counts_search, &
                                                 tmp_residual_to_study_perm, tmp_study_block_ends, tmp_study_block_sizes, &
                                                 m_min=m_min, m_max=m_max, min_residuals_per_bin=min_residuals_per_bin, &
                                                 gamma_occupancy=gamma_occupancy, &
                                                 lower_residual_range_quantile=lower_residual_range_quantile, &
                                                 upper_residual_range_quantile=upper_residual_range_quantile, &
                                                 bin_count_mode=bin_count_mode)

        ! ===== PASS C (per point, per study, sequential: one gather buffer serves every pair).
        ! The bins beyond max_n_bins_per_point are zeroed, as finish_js_comp_test requires.
        max_n_bins_per_point = maxval(n_bins_per_point)
        counts(:, max_n_bins_per_point + 1:, :) = 0_int32
        pmfs(:, max_n_bins_per_point + 1:, :) = 0.0_real64
        call build_adaptive_histograms(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, n_points, &
                                       max_n_bins_per_point, MAX_N_BINS, tmp_gene_means_perm_all, &
                                       pooled_neighborhood_range, n_neighbors_per_point, shared_residual_range_low, &
                                       shared_residual_range_high, n_bins_per_point, counts, pmfs, included_n_reps, &
                                       tmp_point_neighborhood_indices, tmp_neighbor_residuals)

        ! The gene-count check cannot see a study whose genes at a point have only NaN residuals.
        ! Such a study would leave that point's equal-weight consensus summing to less than 1.
        ! It depends on gene_means and residuals together, so no arg_pos.
        if (has_study_without_residuals(included_n_reps, n_points, n_studies)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        call finish_js_comp_test(n_studies, n_points, max_n_bins_per_point, pmfs, counts, included_n_reps, &
                                 actual_divergence_mode, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps, &
                                 js_divergences, weights, global_js_divergence)
    end subroutine run_js_comp_test_adaptive_impl

end module tox_data_integration_js_comp_test_run_impl
