#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Final Run
!|
!| The JS-Comp-Test of K studies for neighborhood parameters already chosen: per-study
!| histograms at each reference point's bin count, each study's JSD against the consensus pmf,
!| their weighting into one global JSD per study, and the permutation test.
!| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]] runs it
!| on fixed-k neighborhoods (`n_points` reference points of `n_neighbors` genes each),
!| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_adaptive_impl(interface)]]
!| on adaptive neighborhoods built by
!| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]].
!|
!| Generated from [[tox_data_integration_js_comp_test_run_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_js_comp_test_run
    use f42_safeguard
    use tox_data_integration_js_comp_test_run_impl, only: prepare_gene_means_perms, run_js_comp_test_adaptive_impl, run_js_comp_test_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_math_impl, only: above
    use tox_data_integration_js_comp_test_binning_impl, only: MAX_N_BINS, MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED
    use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    use tox_errors, only: clear_err_arg_pos, set_err, set_err_once, validate_all_in_range_int
    use tox_errors, only: validate_all_in_range_real, validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: run_js_comp_test
    public :: run_js_comp_test_expert
    public :: run_js_comp_test_adaptive
    public :: run_js_comp_test_adaptive_expert

contains

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test_run(module):run_js_comp_test_expert]] to prepare it yourself.
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
    !| by default, or Stone's rule
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]).
    !| 3. Each study's histogram and pmf is built per point
    !| ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]).
    !| 4. The consensus pmf
    !| ([[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]]), each
    !| study's per-point JSD against it and its weighted global JSD, and the permutation test
    !| ([[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]]) give
    !| `p_values` (primary) and `p_values_observed_consensus` (comparison only). With
    !| `n_permutations = 0` both p-values are NaN.
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
    !| equal-weight consensus needs every study; the consensus, JSD, weight and p-value outputs are
    !| then not filled. `ERR_ALLOC_FAIL` when the permutation test cannot create its random number
    !| generator; `pmfs`, `js_divergences`, `weights` and `global_js_divergence` then hold their
    !| values from before the permutation test. An out-of-range neighbor index is also reported
    !| as `ERR_INVALID_INPUT`, a safeguard valid input cannot reach.
    !|
    !| The plain `run_js_comp_test` sorts each study's gene means itself (heapsort, ascending, NaN
    !| last) and takes no permutation; `run_js_comp_test_expert` takes `gene_means_perms` from the
    !| caller instead. That order must still be ascending with NaN last (only the order among tied
    !| genes is the caller's choice); it is not checked to be sorted.
    !| The allocating wrapper runs [[tox_data_integration_js_comp_test_run_impl(module):prepare_gene_means_perms]] first, which may handle the call and skip this one.
    subroutine run_js_comp_test(&
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
            p_values,&
            p_values_observed_consensus,&
            n_permutations,&
            random_seed,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            bin_count_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_points
            !! Number of reference points for neighborhoods
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
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
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, 256, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(real64), dimension(n_points, 256), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, 256), intent(out) :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            !! comparison only; NaN when `n_permutations` is 0
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                                   |
            !! |----------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)
        integer(int32), dimension(:, :), allocatable :: gene_means_perms
        real(real64), dimension(:, :, :), allocatable :: tmp_neighborhood_residuals_gathered
        real(real64), dimension(:, :), allocatable :: tmp_permutation_mean_pmf
        integer(int32), dimension(:), allocatable :: tmp_permutation_point_pool
        integer(int32), dimension(:), allocatable :: tmp_permutation_point_draw
        real(real64), dimension(:), allocatable :: tmp_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts_search
        integer(int32), dimension(:), allocatable :: tmp_residual_to_study_perm
        integer(int32), dimension(:), allocatable :: tmp_study_block_ends
        integer(int32), dimension(:), allocatable :: tmp_study_block_sizes
        integer(int32), dimension(:), allocatable :: tmp_point_n_neighbors
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_mean_pmf_counts
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_counts
        real(real64), dimension(:, :, :), allocatable :: tmp_permutation_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_permutation_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_permutation_weights
        real(real64), dimension(:), allocatable :: tmp_permutation_global_js_divergence
        logical :: handled

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=5_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_in_range_int(n_permutations, ierr, arg_pos=33_int32, min=0_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=35_int32, min=0_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=36_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=37_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=38_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=39_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=40_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=7_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(x_star, n_points, ierr, arg_pos=8_int32)
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=41_int32); end if
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(gene_means_perms(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, n_neighbors, n_points))
        M_ALLOCATE(tmp_permutation_mean_pmf(n_points, 256))
        M_ALLOCATE(tmp_permutation_point_pool(256))
        M_ALLOCATE(tmp_permutation_point_draw(256))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*n_neighbors*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*n_neighbors*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
        M_ALLOCATE(tmp_residual_to_study_perm(max_n_reps_all_studies*n_neighbors*n_studies))
        M_ALLOCATE(tmp_study_block_ends(n_studies))
        M_ALLOCATE(tmp_study_block_sizes(n_studies))
        M_ALLOCATE(tmp_point_n_neighbors(n_studies))
        M_ALLOCATE(tmp_permutation_mean_pmf_counts(n_points, 256))
        M_ALLOCATE(tmp_permutation_counts(n_points, 256))
        M_ALLOCATE(tmp_permutation_pmfs(n_points, 256, n_studies))
        M_ALLOCATE(tmp_permutation_js_divergences(n_points, n_studies))
        M_ALLOCATE(tmp_permutation_weights(n_points, n_studies))
        M_ALLOCATE(tmp_permutation_global_js_divergence(n_studies))

        call prepare_gene_means_perms(&
            gene_means = gene_means,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_studies = n_studies,&
            gene_means_perms = gene_means_perms,&
            handled = handled,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
        if (is_err(ierr)) return
        if (handled) return

        call run_js_comp_test_impl(&
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
            p_values = p_values,&
            p_values_observed_consensus = p_values_observed_consensus,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_point_pool = tmp_permutation_point_pool,&
            tmp_permutation_point_draw = tmp_permutation_point_draw,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_point_n_neighbors = tmp_point_n_neighbors,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_permutation_global_js_divergence = tmp_permutation_global_js_divergence,&
            n_permutations = n_permutations,&
            random_seed = random_seed,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test_run(module):run_js_comp_test]] does both.
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
    !| by default, or Stone's rule
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]).
    !| 3. Each study's histogram and pmf is built per point
    !| ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]).
    !| 4. The consensus pmf
    !| ([[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]]), each
    !| study's per-point JSD against it and its weighted global JSD, and the permutation test
    !| ([[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]]) give
    !| `p_values` (primary) and `p_values_observed_consensus` (comparison only). With
    !| `n_permutations = 0` both p-values are NaN.
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
    !| equal-weight consensus needs every study; the consensus, JSD, weight and p-value outputs are
    !| then not filled. `ERR_ALLOC_FAIL` when the permutation test cannot create its random number
    !| generator; `pmfs`, `js_divergences`, `weights` and `global_js_divergence` then hold their
    !| values from before the permutation test. An out-of-range neighbor index is also reported
    !| as `ERR_INVALID_INPUT`, a safeguard valid input cannot reach.
    !|
    !| The plain `run_js_comp_test` sorts each study's gene means itself (heapsort, ascending, NaN
    !| last) and takes no permutation; `run_js_comp_test_expert` takes `gene_means_perms` from the
    !| caller instead. That order must still be ascending with NaN last (only the order among tied
    !| genes is the caller's choice); it is not checked to be sorted.
    !| The allocating wrapper runs [[tox_data_integration_js_comp_test_run_impl(module):prepare_gene_means_perms]] first, which may handle the call and skip this one.
    subroutine run_js_comp_test_expert(&
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
            p_values,&
            p_values_observed_consensus,&
            tmp_neighborhood_residuals_gathered,&
            tmp_permutation_mean_pmf,&
            tmp_permutation_point_pool,&
            tmp_permutation_point_draw,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_bin_counts_search,&
            tmp_residual_to_study_perm,&
            tmp_study_block_ends,&
            tmp_study_block_sizes,&
            tmp_point_n_neighbors,&
            tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts,&
            tmp_permutation_pmfs,&
            tmp_permutation_js_divergences,&
            tmp_permutation_weights,&
            tmp_permutation_global_js_divergence,&
            n_permutations,&
            random_seed,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            bin_count_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_points
            !! Number of reference points for neighborhoods
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means_perms
            !! Per study (column), a permutation of that study's `gene_means` in ascending order, NaN last.
            !! The plain entry point builds it for you (heapsort); the expert one takes yours, which
            !! must be sorted that way too (only the order among tied genes is free) and is not
            !! checked to be sorted.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
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
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, 256, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(real64), dimension(n_points, 256), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, 256), intent(out) :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            !! comparison only; NaN when `n_permutations` is 0
        real(real64), dimension(max_n_reps_all_studies, n_neighbors, n_points), intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residuals, reused per study
        real(real64), dimension(n_points, 256), intent(out) :: tmp_permutation_mean_pmf
            !! Working array of the permutation test: consensus pmf of the permuted studies
        integer(int32), dimension(256), intent(out) :: tmp_permutation_point_pool
            !! Working array of the permutation test: one point's pool
        integer(int32), dimension(256), intent(out) :: tmp_permutation_point_draw
            !! Working array of the permutation test: one study's draw at one point
        real(real64), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
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
        integer(int32), dimension(n_points, 256), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array of the permutation test: the remaining pool
        integer(int32), dimension(n_points, 256), intent(out) :: tmp_permutation_counts
            !! Working array of the permutation test: one study's resampled counts
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array of the permutation test: resampled pmfs
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array of the permutation test: per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array of the permutation test: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_permutation_global_js_divergence
            !! Working array of the permutation test: global JSD values
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                                   |
            !! |----------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=5_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_in_range_int(n_permutations, ierr, arg_pos=51_int32, min=0_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=53_int32, min=0_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=54_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=55_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=56_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=57_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=58_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(gene_means_perms, max_n_genes_all_studies * n_studies, ierr, arg_pos=7_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=8_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(x_star, n_points, ierr, arg_pos=9_int32)
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=59_int32); end if
        if (is_err(ierr)) return
#endif

        call run_js_comp_test_impl(&
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
            p_values = p_values,&
            p_values_observed_consensus = p_values_observed_consensus,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_point_pool = tmp_permutation_point_pool,&
            tmp_permutation_point_draw = tmp_permutation_point_draw,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_point_n_neighbors = tmp_point_n_neighbors,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_permutation_global_js_divergence = tmp_permutation_global_js_divergence,&
            n_permutations = n_permutations,&
            random_seed = random_seed,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_adaptive_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test_run(module):run_js_comp_test_adaptive_expert]] to prepare it yourself.
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
    !| residual at some point, or `ERR_ALLOC_FAIL` from the permutation test).
    !|
    !| The pooling work arrays hold `max_n_reps_all_studies * max_n_genes_all_studies * n_studies`
    !| values each, the size of `residuals`, since one study may contribute all its genes to one
    !| neighborhood. The bin-sized arrays have the fixed second extent
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MAX_N_BINS(variable)]]; only their first
    !| `max_n_bins_per_point` bins are meaningful (Python `[:, :max_n_bins_per_point, ...]`, R
    !| `[, seq_len(max_n_bins_per_point), ]`).
    subroutine run_js_comp_test_adaptive(&
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
            p_values,&
            p_values_observed_consensus,&
            n_permutations,&
            random_seed,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            bin_count_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_points
            !! Number of reference points (neighborhoods), the construction's own `n_points`
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene; the array the
            !! neighborhoods were constructed from
            !! NaN is permitted for this value.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(int32), dimension(2, n_points), intent(in) :: pooled_neighborhood_range
            !! For each neighborhood, its first and last sorted pooled position, from
            !! construct_adaptive_neighborhoods
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies*n_studies`.
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
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, 256, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN residuals binned per reference point, per study
        real(real64), dimension(n_points, 256), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, 256), intent(out) :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN residuals per reference point for the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            !! comparison only; NaN when `n_permutations` is 0
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                                   |
            !! |----------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)
        integer(int32), dimension(:), allocatable :: tmp_gene_means_perm_all
        integer(int32), dimension(:, :), allocatable :: tmp_point_neighborhood_indices
        real(real64), dimension(:, :), allocatable :: tmp_neighbor_residuals
        real(real64), dimension(:, :), allocatable :: tmp_permutation_mean_pmf
        integer(int32), dimension(:), allocatable :: tmp_permutation_point_pool
        integer(int32), dimension(:), allocatable :: tmp_permutation_point_draw
        real(real64), dimension(:), allocatable :: tmp_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts_search
        integer(int32), dimension(:), allocatable :: tmp_residual_to_study_perm
        integer(int32), dimension(:), allocatable :: tmp_study_block_ends
        integer(int32), dimension(:), allocatable :: tmp_study_block_sizes
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_mean_pmf_counts
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_counts
        real(real64), dimension(:, :, :), allocatable :: tmp_permutation_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_permutation_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_permutation_weights
        real(real64), dimension(:), allocatable :: tmp_permutation_global_js_divergence

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_permutations, ierr, arg_pos=31_int32, min=0_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=33_int32, min=0_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=34_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=35_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=36_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=37_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=38_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(pooled_neighborhood_range, 2 * n_points, ierr, arg_pos=7_int32, min=1_int32, max=max_n_genes_all_studies*n_studies)
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=39_int32); end if
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_point_neighborhood_indices(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_neighbor_residuals(max_n_reps_all_studies, max_n_genes_all_studies))
        M_ALLOCATE(tmp_permutation_mean_pmf(n_points, 256))
        M_ALLOCATE(tmp_permutation_point_pool(256))
        M_ALLOCATE(tmp_permutation_point_draw(256))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
        M_ALLOCATE(tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_study_block_ends(n_studies))
        M_ALLOCATE(tmp_study_block_sizes(n_studies))
        M_ALLOCATE(tmp_permutation_mean_pmf_counts(n_points, 256))
        M_ALLOCATE(tmp_permutation_counts(n_points, 256))
        M_ALLOCATE(tmp_permutation_pmfs(n_points, 256, n_studies))
        M_ALLOCATE(tmp_permutation_js_divergences(n_points, n_studies))
        M_ALLOCATE(tmp_permutation_weights(n_points, n_studies))
        M_ALLOCATE(tmp_permutation_global_js_divergence(n_studies))

        call run_js_comp_test_adaptive_impl(&
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
            p_values = p_values,&
            p_values_observed_consensus = p_values_observed_consensus,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_point_neighborhood_indices = tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals = tmp_neighbor_residuals,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_point_pool = tmp_permutation_point_pool,&
            tmp_permutation_point_draw = tmp_permutation_point_draw,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_permutation_global_js_divergence = tmp_permutation_global_js_divergence,&
            n_permutations = n_permutations,&
            random_seed = random_seed,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_adaptive

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_adaptive_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test_run(module):run_js_comp_test_adaptive]] does both.
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
    !| residual at some point, or `ERR_ALLOC_FAIL` from the permutation test).
    !|
    !| The pooling work arrays hold `max_n_reps_all_studies * max_n_genes_all_studies * n_studies`
    !| values each, the size of `residuals`, since one study may contribute all its genes to one
    !| neighborhood. The bin-sized arrays have the fixed second extent
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MAX_N_BINS(variable)]]; only their first
    !| `max_n_bins_per_point` bins are meaningful (Python `[:, :max_n_bins_per_point, ...]`, R
    !| `[, seq_len(max_n_bins_per_point), ]`).
    subroutine run_js_comp_test_adaptive_expert(&
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
            p_values,&
            p_values_observed_consensus,&
            tmp_gene_means_perm_all,&
            tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals,&
            tmp_permutation_mean_pmf,&
            tmp_permutation_point_pool,&
            tmp_permutation_point_draw,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_bin_counts_search,&
            tmp_residual_to_study_perm,&
            tmp_study_block_ends,&
            tmp_study_block_sizes,&
            tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts,&
            tmp_permutation_pmfs,&
            tmp_permutation_js_divergences,&
            tmp_permutation_weights,&
            tmp_permutation_global_js_divergence,&
            n_permutations,&
            random_seed,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            bin_count_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_points
            !! Number of reference points (neighborhoods), the construction's own `n_points`
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene; the array the
            !! neighborhoods were constructed from
            !! NaN is permitted for this value.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(int32), dimension(2, n_points), intent(in) :: pooled_neighborhood_range
            !! For each neighborhood, its first and last sorted pooled position, from
            !! construct_adaptive_neighborhoods
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies*n_studies`.
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
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, 256, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN residuals binned per reference point, per study
        real(real64), dimension(n_points, 256), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, 256), intent(out) :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN residuals per reference point for the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            !! comparison only; NaN when `n_permutations` is 0
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, as in the construction
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_point_neighborhood_indices
            !! Working array: one reference point's gene indices per study, reused per point
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out) :: tmp_neighbor_residuals
            !! Working array: one (reference point, study) pair's gathered residuals, reused per pair
        real(real64), dimension(n_points, 256), intent(out) :: tmp_permutation_mean_pmf
            !! Working array forwarded to the permutation test
        integer(int32), dimension(256), intent(out) :: tmp_permutation_point_pool
            !! Working array forwarded to the permutation test
        integer(int32), dimension(256), intent(out) :: tmp_permutation_point_draw
            !! Working array forwarded to the permutation test
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for `tmp_pooled_residuals`, reused per point
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array: the last index of each study's block in `tmp_pooled_residuals`
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in `tmp_pooled_residuals`
        integer(int32), dimension(n_points, 256), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array forwarded to the permutation test
        integer(int32), dimension(n_points, 256), intent(out) :: tmp_permutation_counts
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_studies), intent(out) :: tmp_permutation_global_js_divergence
            !! Working array forwarded to the permutation test
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                                   |
            !! |----------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_permutations, ierr, arg_pos=49_int32, min=0_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=51_int32, min=0_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=52_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=53_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=54_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=55_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=56_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(pooled_neighborhood_range, 2 * n_points, ierr, arg_pos=7_int32, min=1_int32, max=max_n_genes_all_studies*n_studies)
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=57_int32); end if
        if (is_err(ierr)) return
#endif

        call run_js_comp_test_adaptive_impl(&
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
            p_values = p_values,&
            p_values_observed_consensus = p_values_observed_consensus,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_point_neighborhood_indices = tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals = tmp_neighbor_residuals,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_point_pool = tmp_permutation_point_pool,&
            tmp_permutation_point_draw = tmp_permutation_point_draw,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_permutation_global_js_divergence = tmp_permutation_global_js_divergence,&
            n_permutations = n_permutations,&
            random_seed = random_seed,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_adaptive_expert

end module tox_data_integration_js_comp_test_run
