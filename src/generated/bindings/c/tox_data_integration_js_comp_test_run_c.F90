#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_run(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Final Run
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
module tox_data_integration_js_comp_test_run_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_as_view
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: run_js_comp_test_c
    public :: run_js_comp_test_expert_c
    public :: run_js_comp_test_adaptive_c
    public :: run_js_comp_test_adaptive_expert_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_run(module):run_js_comp_test(subroutine)]]
    !| The final JS-Comp-Test of K studies for one fixed `(n_points, n_neighbors)` setting:
    !|
    !| 1. Each study's `n_neighbors` genes nearest every reference point `x_star` (computed once by
    !| the caller, e.g. with
    !| [[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) form its
    !| neighborhoods
    !| ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]).
    !| 2. Per reference point, all studies' residuals are pooled and the bin-count search
    !| `bin_count_mode` names chooses the point's bin count and residual range (occupancy search
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !| or Stone's rule
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]],
    !| by default with every study's score counted equally).
    !| 3. Each study's histogram and pmf is built per point
    !| ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]).
    !| 4. The consensus pmf
    !| ([[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]]), each study's
    !| per-point JSD and its weighted global JSD: by default (`divergence_mode`), for exactly two
    !| studies, the pairwise JSD of the two studies, the same for both; otherwise each study's JSD
    !| against the consensus.
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
    !| (`256`); only their first `max_n_bins_per_point` bins are meaningful, the rest is
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
    !| The allocating wrapper runs [[tox_data_integration_js_comp_test_run_impl(module):prepare_gene_means_perms]] first, which may handle the call and skip this one.
    subroutine run_js_comp_test_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            n_points,&
            n_neighbors,&
            gene_means,&
            residuals,&
            x_star,&
            neighborhood_indices,&
            neighborhood_range,&
            n_bins_per_point,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            max_n_bins_per_point,&
            occupancy_failed,&
            n_pooled_residuals,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            sturges_bins,&
            fd_bins,&
            pmfs,&
            counts,&
            included_n_reps,&
            mean_pmf,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            js_divergences,&
            weights,&
            global_js_divergence,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            bin_count_mode,&
            divergence_mode,&
            ierr&
        ) bind(C, name="run_js_comp_test_c")
        use tox_data_integration_js_comp_test_run, only: run_js_comp_test
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC

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
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        real(c_double), dimension(n_points), intent(in), target :: x_star
            !! Mean-expression reference points
        integer(c_int), dimension(n_neighbors, n_points, n_studies), intent(out), target :: neighborhood_indices
            !! Gene indices of the selected neighborhood, per reference point, per study
        integer(c_int), dimension(2, n_points, n_studies), intent(out), target :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span
        integer(c_int), dimension(n_points), intent(out), target :: n_bins_per_point
            !! Each reference point's histogram bin count `M_j`, from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(c_int), intent(out), target :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts` (see above)
        logical(c_bool), dimension(n_points), intent(out), target :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes (see above)
        integer(c_int), dimension(n_points), intent(out), target :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(c_int), dimension(n_points), intent(out), target :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(c_double), dimension(n_points), intent(out), target :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(c_int), dimension(n_points), intent(out), target :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(c_double), dimension(n_points, 256, n_studies), intent(out), target :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(c_int), dimension(n_points, 256, n_studies), intent(out), target :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(c_double), dimension(n_points, 256), intent(out), target :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, 256), intent(out), target :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study, the one `divergence_mode` selects
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study, the one `divergence_mode` selects
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                | Value                                                                                                   |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin`    | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score (default)                                     | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                       | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_stone_equal'`.
        character(len=1, kind=c_char), dimension(22), intent(in), target :: divergence_mode
            !! How each study's per-point divergence is computed
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `'divergence_pairwise_k2'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT (see above)
        integer(int32) :: bin_count_mode_mode_f
        integer(int32) :: divergence_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(max_n_bins_per_point)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_indices, n_neighbors * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(neighborhood_range, 2 * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_low, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_high, n_points)
        M_CHECK_ARRAY_NON_NULL(occupancy_failed, n_points)
        M_CHECK_ARRAY_NON_NULL(n_pooled_residuals, n_points)
        M_CHECK_ARRAY_NON_NULL(min_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(max_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(sturges_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(fd_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)
        M_CHECK_ARRAY_NON_NULL(divergence_mode, 22)

        block
            character(len=:), pointer :: bin_count_mode_f
            bin_count_mode_f => c_char_as_view(bin_count_mode)

            select case (bin_count_mode_f)
                case ("bin_count_occupancy")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_OCCUPANCY
                case ("bin_count_stone_equal")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_EQUAL
                case ("bin_count_stone_size_weighted")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
                case ("bin_count_stone_pooled")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_POOLED
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        block
            character(len=:), pointer :: divergence_mode_f
            divergence_mode_f => c_char_as_view(divergence_mode)

            select case (divergence_mode_f)
                case ("divergence_symmetric")
                    divergence_mode_mode_f = MODE_DIVERGENCE_SYMMETRIC
                case ("divergence_pairwise_k2")
                    divergence_mode_mode_f = MODE_DIVERGENCE_PAIRWISE_K2
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call run_js_comp_test(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            gene_means = gene_means,&
            residuals = residuals,&
            x_star = x_star,&
            neighborhood_indices = neighborhood_indices,&
            neighborhood_range = neighborhood_range,&
            n_bins_per_point = n_bins_per_point,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            max_n_bins_per_point = max_n_bins_per_point,&
            occupancy_failed = occupancy_failed,&
            n_pooled_residuals = n_pooled_residuals,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            pmfs = pmfs,&
            counts = counts,&
            included_n_reps = included_n_reps,&
            mean_pmf = mean_pmf,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            js_divergences = js_divergences,&
            weights = weights,&
            global_js_divergence = global_js_divergence,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            bin_count_mode = bin_count_mode_mode_f,&
            divergence_mode = divergence_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_run(module):run_js_comp_test_expert(subroutine)]]
    !| The final JS-Comp-Test of K studies for one fixed `(n_points, n_neighbors)` setting:
    !|
    !| 1. Each study's `n_neighbors` genes nearest every reference point `x_star` (computed once by
    !| the caller, e.g. with
    !| [[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) form its
    !| neighborhoods
    !| ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]).
    !| 2. Per reference point, all studies' residuals are pooled and the bin-count search
    !| `bin_count_mode` names chooses the point's bin count and residual range (occupancy search
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !| or Stone's rule
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]],
    !| by default with every study's score counted equally).
    !| 3. Each study's histogram and pmf is built per point
    !| ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]).
    !| 4. The consensus pmf
    !| ([[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]]), each study's
    !| per-point JSD and its weighted global JSD: by default (`divergence_mode`), for exactly two
    !| studies, the pairwise JSD of the two studies, the same for both; otherwise each study's JSD
    !| against the consensus.
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
    !| (`256`); only their first `max_n_bins_per_point` bins are meaningful, the rest is
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
    !| The allocating wrapper runs [[tox_data_integration_js_comp_test_run_impl(module):prepare_gene_means_perms]] first, which may handle the call and skip this one.
    subroutine run_js_comp_test_expert_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            n_points,&
            n_neighbors,&
            gene_means,&
            gene_means_perms,&
            residuals,&
            x_star,&
            neighborhood_indices,&
            neighborhood_range,&
            n_bins_per_point,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            max_n_bins_per_point,&
            occupancy_failed,&
            n_pooled_residuals,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            sturges_bins,&
            fd_bins,&
            pmfs,&
            counts,&
            included_n_reps,&
            mean_pmf,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            js_divergences,&
            weights,&
            global_js_divergence,&
            tmp_neighborhood_residuals_gathered,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_bin_counts_search,&
            tmp_residual_to_study_perm,&
            tmp_study_block_ends,&
            tmp_study_block_sizes,&
            tmp_point_n_neighbors,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            bin_count_mode,&
            divergence_mode,&
            ierr&
        ) bind(C, name="run_js_comp_test_expert_c")
        use tox_data_integration_js_comp_test_run, only: run_js_comp_test_expert
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC

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
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means_perms
            !! Per study (column), a permutation of that study's `gene_means` in ascending order, NaN last.
            !! The plain entry point builds it for you (heapsort); the expert one takes yours, which
            !! must be sorted that way too (only the order among tied genes is free). It is checked
            !! neither to be a permutation of `1..max_n_genes_all_studies` per column nor to be sorted.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        real(c_double), dimension(n_points), intent(in), target :: x_star
            !! Mean-expression reference points
        integer(c_int), dimension(n_neighbors, n_points, n_studies), intent(out), target :: neighborhood_indices
            !! Gene indices of the selected neighborhood, per reference point, per study
        integer(c_int), dimension(2, n_points, n_studies), intent(out), target :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span
        integer(c_int), dimension(n_points), intent(out), target :: n_bins_per_point
            !! Each reference point's histogram bin count `M_j`, from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(c_int), intent(out), target :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts` (see above)
        logical(c_bool), dimension(n_points), intent(out), target :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes (see above)
        integer(c_int), dimension(n_points), intent(out), target :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(c_int), dimension(n_points), intent(out), target :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(c_double), dimension(n_points), intent(out), target :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(c_int), dimension(n_points), intent(out), target :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(c_double), dimension(n_points, 256, n_studies), intent(out), target :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(c_int), dimension(n_points, 256, n_studies), intent(out), target :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(c_double), dimension(n_points, 256), intent(out), target :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, 256), intent(out), target :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study, the one `divergence_mode` selects
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study, the one `divergence_mode` selects
        real(c_double), dimension(max_n_reps_all_studies, n_neighbors, n_points), intent(out), target :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residuals, reused per study
        real(c_double), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out), target :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out), target :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out), target :: tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_ends
            !! Working array: the last index of each study's block in tmp_pooled_residuals
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in tmp_pooled_residuals
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_point_n_neighbors
            !! Working array: every study's neighbor count at a reference point (all `n_neighbors`)
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                | Value                                                                                                   |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin`    | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score (default)                                     | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                       | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_stone_equal'`.
        character(len=1, kind=c_char), dimension(22), intent(in), target :: divergence_mode
            !! How each study's per-point divergence is computed
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `'divergence_pairwise_k2'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT (see above)
        integer(int32) :: bin_count_mode_mode_f
        integer(int32) :: divergence_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(max_n_bins_per_point)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(gene_means_perms, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_indices, n_neighbors * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(neighborhood_range, 2 * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_low, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_high, n_points)
        M_CHECK_ARRAY_NON_NULL(occupancy_failed, n_points)
        M_CHECK_ARRAY_NON_NULL(n_pooled_residuals, n_points)
        M_CHECK_ARRAY_NON_NULL(min_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(max_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(sturges_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(fd_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_residuals_gathered, max_n_reps_all_studies * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals, (max_n_reps_all_studies*n_neighbors*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals_perm, (max_n_reps_all_studies*n_neighbors*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_bin_counts_search, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_residual_to_study_perm, (max_n_reps_all_studies*n_neighbors*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_ends, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_sizes, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_point_n_neighbors, n_studies)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)
        M_CHECK_ARRAY_NON_NULL(divergence_mode, 22)

        block
            character(len=:), pointer :: bin_count_mode_f
            bin_count_mode_f => c_char_as_view(bin_count_mode)

            select case (bin_count_mode_f)
                case ("bin_count_occupancy")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_OCCUPANCY
                case ("bin_count_stone_equal")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_EQUAL
                case ("bin_count_stone_size_weighted")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
                case ("bin_count_stone_pooled")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_POOLED
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        block
            character(len=:), pointer :: divergence_mode_f
            divergence_mode_f => c_char_as_view(divergence_mode)

            select case (divergence_mode_f)
                case ("divergence_symmetric")
                    divergence_mode_mode_f = MODE_DIVERGENCE_SYMMETRIC
                case ("divergence_pairwise_k2")
                    divergence_mode_mode_f = MODE_DIVERGENCE_PAIRWISE_K2
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call run_js_comp_test_expert(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            gene_means = gene_means,&
            gene_means_perms = gene_means_perms,&
            residuals = residuals,&
            x_star = x_star,&
            neighborhood_indices = neighborhood_indices,&
            neighborhood_range = neighborhood_range,&
            n_bins_per_point = n_bins_per_point,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            max_n_bins_per_point = max_n_bins_per_point,&
            occupancy_failed = occupancy_failed,&
            n_pooled_residuals = n_pooled_residuals,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            pmfs = pmfs,&
            counts = counts,&
            included_n_reps = included_n_reps,&
            mean_pmf = mean_pmf,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            js_divergences = js_divergences,&
            weights = weights,&
            global_js_divergence = global_js_divergence,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_point_n_neighbors = tmp_point_n_neighbors,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            bin_count_mode = bin_count_mode_mode_f,&
            divergence_mode = divergence_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_run(module):run_js_comp_test_adaptive(subroutine)]]
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
    subroutine run_js_comp_test_adaptive_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            n_points,&
            gene_means,&
            residuals,&
            pooled_neighborhood_range,&
            n_neighbors_per_point,&
            n_bins_per_point,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            max_n_bins_per_point,&
            occupancy_failed,&
            n_pooled_residuals,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            sturges_bins,&
            fd_bins,&
            pmfs,&
            counts,&
            included_n_reps,&
            mean_pmf,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            js_divergences,&
            weights,&
            global_js_divergence,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            bin_count_mode,&
            divergence_mode,&
            ierr&
        ) bind(C, name="run_js_comp_test_adaptive_c")
        use tox_data_integration_js_comp_test_run, only: run_js_comp_test_adaptive
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC

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
            !! Number of reference points (neighborhoods), the construction's own `n_points`
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene; the array the
            !! neighborhoods were constructed from
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(c_int), dimension(2, n_points), intent(in), target :: pooled_neighborhood_range
            !! For each neighborhood, its first and last sorted pooled position, from
            !! construct_adaptive_neighborhoods
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies*n_studies`.
        integer(c_int), dimension(n_studies, n_points), intent(out), target :: n_neighbors_per_point
            !! For each neighborhood, how many of its genes belong to each study (at least 1 each)
        integer(c_int), dimension(n_points), intent(out), target :: n_bins_per_point
            !! Each reference point's selected histogram bin count, from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(c_int), intent(out), target :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts`
        logical(c_bool), dimension(n_points), intent(out), target :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes
        integer(c_int), dimension(n_points), intent(out), target :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(c_int), dimension(n_points), intent(out), target :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(c_double), dimension(n_points), intent(out), target :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(c_int), dimension(n_points), intent(out), target :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(c_double), dimension(n_points, 256, n_studies), intent(out), target :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(c_int), dimension(n_points, 256, n_studies), intent(out), target :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN residuals binned per reference point, per study
        real(c_double), dimension(n_points, 256), intent(out), target :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, 256), intent(out), target :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN residuals per reference point for the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study, the one `divergence_mode` selects
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study, the one `divergence_mode` selects
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                | Value                                                                                                   |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin`    | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score (default)                                     | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                       | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_stone_equal'`.
        character(len=1, kind=c_char), dimension(22), intent(in), target :: divergence_mode
            !! How each study's per-point divergence is computed
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `'divergence_pairwise_k2'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT (see above)
        integer(int32) :: bin_count_mode_mode_f
        integer(int32) :: divergence_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(max_n_bins_per_point)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(pooled_neighborhood_range, 2 * n_points)
        M_CHECK_ARRAY_NON_NULL(n_neighbors_per_point, n_studies * n_points)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_low, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_high, n_points)
        M_CHECK_ARRAY_NON_NULL(occupancy_failed, n_points)
        M_CHECK_ARRAY_NON_NULL(n_pooled_residuals, n_points)
        M_CHECK_ARRAY_NON_NULL(min_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(max_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(sturges_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(fd_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)
        M_CHECK_ARRAY_NON_NULL(divergence_mode, 22)

        block
            character(len=:), pointer :: bin_count_mode_f
            bin_count_mode_f => c_char_as_view(bin_count_mode)

            select case (bin_count_mode_f)
                case ("bin_count_occupancy")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_OCCUPANCY
                case ("bin_count_stone_equal")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_EQUAL
                case ("bin_count_stone_size_weighted")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
                case ("bin_count_stone_pooled")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_POOLED
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        block
            character(len=:), pointer :: divergence_mode_f
            divergence_mode_f => c_char_as_view(divergence_mode)

            select case (divergence_mode_f)
                case ("divergence_symmetric")
                    divergence_mode_mode_f = MODE_DIVERGENCE_SYMMETRIC
                case ("divergence_pairwise_k2")
                    divergence_mode_mode_f = MODE_DIVERGENCE_PAIRWISE_K2
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call run_js_comp_test_adaptive(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_points = n_points,&
            gene_means = gene_means,&
            residuals = residuals,&
            pooled_neighborhood_range = pooled_neighborhood_range,&
            n_neighbors_per_point = n_neighbors_per_point,&
            n_bins_per_point = n_bins_per_point,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            max_n_bins_per_point = max_n_bins_per_point,&
            occupancy_failed = occupancy_failed,&
            n_pooled_residuals = n_pooled_residuals,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            pmfs = pmfs,&
            counts = counts,&
            included_n_reps = included_n_reps,&
            mean_pmf = mean_pmf,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            js_divergences = js_divergences,&
            weights = weights,&
            global_js_divergence = global_js_divergence,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            bin_count_mode = bin_count_mode_mode_f,&
            divergence_mode = divergence_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_adaptive_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_run(module):run_js_comp_test_adaptive_expert(subroutine)]]
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
    subroutine run_js_comp_test_adaptive_expert_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            n_points,&
            gene_means,&
            residuals,&
            pooled_neighborhood_range,&
            n_neighbors_per_point,&
            n_bins_per_point,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            max_n_bins_per_point,&
            occupancy_failed,&
            n_pooled_residuals,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            sturges_bins,&
            fd_bins,&
            pmfs,&
            counts,&
            included_n_reps,&
            mean_pmf,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            js_divergences,&
            weights,&
            global_js_divergence,&
            tmp_gene_means_perm_all,&
            tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_bin_counts_search,&
            tmp_residual_to_study_perm,&
            tmp_study_block_ends,&
            tmp_study_block_sizes,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            bin_count_mode,&
            divergence_mode,&
            ierr&
        ) bind(C, name="run_js_comp_test_adaptive_expert_c")
        use tox_data_integration_js_comp_test_run, only: run_js_comp_test_adaptive_expert
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC

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
            !! Number of reference points (neighborhoods), the construction's own `n_points`
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene; the array the
            !! neighborhoods were constructed from
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(c_int), dimension(2, n_points), intent(in), target :: pooled_neighborhood_range
            !! For each neighborhood, its first and last sorted pooled position, from
            !! construct_adaptive_neighborhoods
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies*n_studies`.
        integer(c_int), dimension(n_studies, n_points), intent(out), target :: n_neighbors_per_point
            !! For each neighborhood, how many of its genes belong to each study (at least 1 each)
        integer(c_int), dimension(n_points), intent(out), target :: n_bins_per_point
            !! Each reference point's selected histogram bin count, from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(c_int), intent(out), target :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts`
        logical(c_bool), dimension(n_points), intent(out), target :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes
        integer(c_int), dimension(n_points), intent(out), target :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(c_int), dimension(n_points), intent(out), target :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(c_double), dimension(n_points), intent(out), target :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(c_int), dimension(n_points), intent(out), target :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(c_double), dimension(n_points, 256, n_studies), intent(out), target :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(c_int), dimension(n_points, 256, n_studies), intent(out), target :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN residuals binned per reference point, per study
        real(c_double), dimension(n_points, 256), intent(out), target :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, 256), intent(out), target :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN residuals per reference point for the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study, the one `divergence_mode` selects
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study, the one `divergence_mode` selects
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, as in the construction
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(out), target :: tmp_point_neighborhood_indices
            !! Working array: one reference point's gene indices per study, reused per point
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out), target :: tmp_neighbor_residuals
            !! Working array: one (reference point, study) pair's gathered residuals, reused per pair
        real(c_double), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out), target :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out), target :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for `tmp_pooled_residuals`, reused per point
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out), target :: tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_ends
            !! Working array: the last index of each study's block in `tmp_pooled_residuals`
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in `tmp_pooled_residuals`
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                | Value                                                                                                   |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin`    | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score (default)                                     | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                       | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_stone_equal'`.
        character(len=1, kind=c_char), dimension(22), intent(in), target :: divergence_mode
            !! How each study's per-point divergence is computed
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `'divergence_pairwise_k2'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT (see above)
        integer(int32) :: bin_count_mode_mode_f
        integer(int32) :: divergence_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(max_n_bins_per_point)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(pooled_neighborhood_range, 2 * n_points)
        M_CHECK_ARRAY_NON_NULL(n_neighbors_per_point, n_studies * n_points)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_low, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_high, n_points)
        M_CHECK_ARRAY_NON_NULL(occupancy_failed, n_points)
        M_CHECK_ARRAY_NON_NULL(n_pooled_residuals, n_points)
        M_CHECK_ARRAY_NON_NULL(min_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(max_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(sturges_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(fd_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perm_all, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_point_neighborhood_indices, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighbor_residuals, max_n_reps_all_studies * max_n_genes_all_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals, (max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals_perm, (max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_bin_counts_search, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_residual_to_study_perm, (max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_ends, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_sizes, n_studies)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)
        M_CHECK_ARRAY_NON_NULL(divergence_mode, 22)

        block
            character(len=:), pointer :: bin_count_mode_f
            bin_count_mode_f => c_char_as_view(bin_count_mode)

            select case (bin_count_mode_f)
                case ("bin_count_occupancy")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_OCCUPANCY
                case ("bin_count_stone_equal")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_EQUAL
                case ("bin_count_stone_size_weighted")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
                case ("bin_count_stone_pooled")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_POOLED
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        block
            character(len=:), pointer :: divergence_mode_f
            divergence_mode_f => c_char_as_view(divergence_mode)

            select case (divergence_mode_f)
                case ("divergence_symmetric")
                    divergence_mode_mode_f = MODE_DIVERGENCE_SYMMETRIC
                case ("divergence_pairwise_k2")
                    divergence_mode_mode_f = MODE_DIVERGENCE_PAIRWISE_K2
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call run_js_comp_test_adaptive_expert(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_points = n_points,&
            gene_means = gene_means,&
            residuals = residuals,&
            pooled_neighborhood_range = pooled_neighborhood_range,&
            n_neighbors_per_point = n_neighbors_per_point,&
            n_bins_per_point = n_bins_per_point,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            max_n_bins_per_point = max_n_bins_per_point,&
            occupancy_failed = occupancy_failed,&
            n_pooled_residuals = n_pooled_residuals,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            pmfs = pmfs,&
            counts = counts,&
            included_n_reps = included_n_reps,&
            mean_pmf = mean_pmf,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            js_divergences = js_divergences,&
            weights = weights,&
            global_js_divergence = global_js_divergence,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_point_neighborhood_indices = tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals = tmp_neighbor_residuals,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            bin_count_mode = bin_count_mode_mode_f,&
            divergence_mode = divergence_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_adaptive_expert_c

end module tox_data_integration_js_comp_test_run_c
#endif
