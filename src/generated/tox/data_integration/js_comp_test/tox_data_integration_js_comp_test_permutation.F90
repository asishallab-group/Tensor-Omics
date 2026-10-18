#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Study-Label Permutation-Test
!|
!| Asks whether the assignment of samples to studies (batches) influences the JS-Comp-Test values
!| at all. Null hypothesis: it does not. The test shuffles the study labels of the samples
!| (without replacement), reruns the whole JS-Comp-Test on every shuffled data set, and compares
!| each study's observed global JSD with the values the shuffled data sets give.
!|
!| **Input layout.** All studies arrive in one array `expression(max_n_reps_all_studies,
!| max_n_genes_all_studies, n_studies)`: gene `g` is the same gene in every study, a gene missing in
!| a study is NaN there, and a replicate slot `(r, s)` whose values are all NaN is padding (a study
!| with fewer samples than `max_n_reps_all_studies`). Every other slot is a sample. Samples are
!| numbered in slot order, replicate fastest, then study.
!|
!| **Entry points.** One assignment of samples to slots, one JS-Comp-Test:
!| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]]
!| (k-NN neighborhoods, the fixed `(n_points, n_neighbors)` of `run_js_comp_test`) and
!| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]
!| (adaptive neighborhoods, `run_js_comp_test_adaptive`). The whole test, `n_permutations` random
!| assignments plus the observed one:
!| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_label_permutation_test_impl(interface)]]
!| and
!| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_label_permutation_test_impl(interface)]].
!| `calc_js_comp_test_permutation_bounds` and `calc_adaptive_js_comp_test_permutation_bounds` size
!| the work arrays for a caller that allocates its own.
!|
!| **What a p-value means.** Study `i`'s p-value is `(count + 1)/(n_valid + 1)`, where `count` is
!| how many of the `n_valid` valid shuffled data sets give study `i` a global JSD at least as large
!| as the observed one (the +1 counts the observed data itself, so it is never 0). A small p-value
!| (e.g. 0.01) says: hardly any random relabelling makes study `i` look as different from the
!| others as the real labels do, so its samples differ systematically (a batch effect). A large
!| one says the real labels are nothing special.
!|
!| This test is not the residual-level permutation test of `run_js_comp_test`
!| ([[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]]), which keeps
!| the data and the histograms fixed and only reshuffles the pooled histogram counts.
!|
!| Generated from [[tox_data_integration_js_comp_test_permutation_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_js_comp_test_permutation
    use f42_safeguard
    use tox_data_integration_js_comp_test_permutation_impl, only: calc_adaptive_js_comp_test_permutation_bounds, calc_js_comp_test_permutation_bounds, run_js_comp_test_adaptive_label_permutation_test_impl, run_js_comp_test_adaptive_on_sample_assignment_impl
    use tox_data_integration_js_comp_test_permutation_impl, only: run_js_comp_test_label_permutation_test_impl, run_js_comp_test_on_sample_assignment_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_math_impl, only: above
    use tox_data_integration_js_comp_test_binning_impl, only: MAX_N_BINS, MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED
    use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
    use tox_data_integration_js_comp_test_fixed_k_search_impl, only: MODE_NEIGHBORHOOD_SIZE_BOTH, MODE_NEIGHBORHOOD_SIZE_LARGE, MODE_NEIGHBORHOOD_SIZE_SMALL
    use tox_data_integration_js_comp_test_search_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH
    use tox_data_integration_js_comp_test_search_impl, only: MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE
    use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_ONE_SIDED, MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC, MODE_DIVERGENCE_WEIGHT_POOL
    use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL
    use tox_data_integration_stats_impl, only: MODE_BOOTSTRAP_PER_STUDY, MODE_BOOTSTRAP_POOLED
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    use tox_errors, only: clear_err_arg_pos, set_err, set_err_once, validate_all_in_range_int
    use tox_errors, only: validate_all_in_range_real, validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: run_js_comp_test_on_sample_assignment
    public :: run_js_comp_test_on_sample_assignment_expert
    public :: run_js_comp_test_adaptive_on_sample_assignment
    public :: run_js_comp_test_adaptive_on_sample_assignment_expert
    public :: run_js_comp_test_label_permutation_test
    public :: run_js_comp_test_label_permutation_test_expert
    public :: run_js_comp_test_adaptive_label_permutation_test
    public :: run_js_comp_test_adaptive_label_permutation_test_expert

contains

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test_permutation(module):run_js_comp_test_on_sample_assignment_expert]] to prepare it yourself.
    !| One building block of the study-label Permutation-Test: relabels the samples by
    !| `sample_assignment` and runs the whole JS-Comp-Test on the relabelled data set.
    !|
    !| **Input layout.** `expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)`
    !| holds every study: gene `g` is the same gene in every study, a gene missing in a study is NaN
    !| there, and a replicate slot whose values are all NaN is padding (a study with fewer samples).
    !| Every other slot is a sample, numbered in slot order: replicate fastest, then study (sample 1
    !| is the first non-padding slot of study 1).
    !|
    !| **Assignment.** `sample_assignment` is a permutation of `1..n_samples`: destination sample
    !| `k` (the `k`-th non-padding slot) receives all genes of source sample `sample_assignment(k)`.
    !| So every study keeps its slots, its sample count and its padding; only which samples sit in
    !| them changes. The identity `1..n_samples` is the observed data.
    !|
    !| **Steps.** Per study, the gene means and residuals of the relabelled data
    !| ([[tox_data_integration_preprocessing_impl(module):compute_gene_means_impl(interface)]],
    !| [[tox_data_integration_preprocessing_impl(module):compute_residuals_impl(interface)]]); the
    !| reference points from the pooled means
    !| ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]] on their
    !| ascending, NaN-last order); then
    !| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]] with
    !| `n_permutations = 0` (its own residual-level permutation test does not run) and the bin-count
    !| options as given. On the identity assignment the result is exactly that of a caller doing
    !| these steps on the observed data with `run_js_comp_test`.
    !|
    !| **Rerunning the parameter search.** With `rerun_parameter_search = .true.`, the relabelled
    !| data set first runs
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| with the search options below, and the run then uses the `(n_points, n_neighbors)` it returns
    !| (also its fallback when no candidate was admissible); `n_points` and `n_neighbors` are then
    !| not used for the run, only to size the work arrays. `used_n_points` and `used_n_neighbors`
    !| report what the run used. Run time: every call runs the whole JS-Comp-Test; with the rerun
    !| option it runs the whole parameter search first.
    !|
    !| **Known limitation.** A gene missing in a whole study (NaN in all its samples) travels with
    !| the samples, so after relabelling that study may have values for it.
    !|
    !| Errors (`ierr`): the call fails with invalid input (`ERR_INVALID_INPUT`) when some study has
    !| no sample at all, when `n_samples` is not the number of samples in `expression`, or when
    !| `sample_assignment` is not a permutation of `1..n_samples`. (Fortran expert tier only, where
    !| the caller supplies the work-array capacities: with the rerun option it also fails when a
    !| capacity is below the search's own bound, `calc_js_comp_test_candidate_bounds`.) Errors of
    !| the search and the run are returned with their code (for example `ERR_INVALID_INPUT` when some
    !| study has no residual at some reference point of the relabelled data). `global_js_divergence`
    !| is then NaN.
    subroutine run_js_comp_test_on_sample_assignment(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            expression,&
            n_samples,&
            sample_assignment,&
            n_points,&
            n_neighbors,&
            global_js_divergence,&
            used_n_points,&
            used_n_neighbors,&
            rerun_parameter_search,&
            n_bootstraps,&
            join_method,&
            min_residuals_per_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            bin_count_mode,&
            neighborhood_size_mode,&
            delta_abs_tolerance,&
            bootstrap_mode,&
            divergence_mode,&
            divergence_weight_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_samples
            !! Number of samples, the non-padding slots of `expression`
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
            !! A permutation of `1..n_samples`: destination sample `k` receives source sample
            !! `sample_assignment(k)`; the identity is the observed data
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_samples`.
        integer(int32), intent(in) :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus, for this assignment; NaN on error
        integer(int32), intent(out) :: used_n_points
            !! Number of reference points the run used (`n_points`, or the search's choice)
        integer(int32), intent(out) :: used_n_neighbors
            !! Number of neighbors per neighborhood the run used (`n_neighbors`, or the search's choice)
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! `.true.`: run the fixed-k parameter search on the relabelled data first and use its
            !! `(n_points, n_neighbors)`; `.false.`: use `n_points` and `n_neighbors`
            !! The default value is `.false.`.
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `1_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Rerun search: seed of the bootstrap's random number generator
            !! The default value is `42_int32`.
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
        integer(int32), intent(in), optional :: neighborhood_size_mode
            !! Rerun search: which neighborhood sizes the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are left out
            !! by default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode                                                       | Value                                                                                                    |
            !! |------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Large neighborhoods only, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods only, `KX_FACTORS(2)` = 0.5            | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! | Both, large then small for each `n_points`                 | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_BOTH(variable)]]  |
            !! The default value is `0_int32`.
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0_real64`.
        integer(int32), intent(in), optional :: bootstrap_mode
            !! Rerun search (used only with `rerun_parameter_search`): what its bootstrap draws each
            !! study's histograms from. By default each study's own counts: a robustness interval of
            !! the observed JSD. The pooled counts of all studies give a noise band of "no difference
            !! between the studies" instead, kept for comparison
            !!
            !! | Mode                                                                                             | Value                                                                          |
            !! |--------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
            !! | Each study's own counts, a robustness interval (default)                                         | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_PER_STUDY(variable)]] |
            !! | The pooled consensus counts of all studies, a noise band (the bootstrap before the mode existed) | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_POOLED(variable)]]    |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_mode
            !! Experiment option: how each study's per-point divergence is computed. By default the
            !! JSD of the study against the consensus pmf, the formula before the option existed
            !!
            !! | Mode                                                                                                                | Value                                                                           |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | JSD of each study against the consensus pmf (default)                                                               | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | For exactly two studies the pairwise JSD of the two, the same value for both (otherwise as the default)             | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! | One-sided Kullback-Leibler divergence of each study against the consensus pmf, at most log2(number of studies) bits | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_ONE_SIDED(variable)]]   |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_weight_mode
            !! Experiment option: how each reference point's weight in a study's global divergence is
            !! computed. By default the study's residual count there plus the count of all studies
            !! pooled, the weights before the option existed. For two studies in the pairwise
            !! `divergence_mode` the weights are always the two studies' counts
            !!
            !! | Mode                                                                   | Value                                                                                      |
            !! |------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
            !! | The study's count plus the pooled count of all studies (default)       | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL(variable)]] |
            !! | The pooled count of all studies only, the same weights for every study | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_POOL(variable)]]            |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)
        integer(int32) :: max_n_points_candidate
        integer(int32) :: max_n_neighbors_candidate
        integer(int32) :: n_bootstrapping_top_k_jsds
        integer(int32), dimension(:), allocatable :: tmp_sample_slots
        logical(c_bool), dimension(:), allocatable :: tmp_sample_mark
        real(real64), dimension(:, :, :), allocatable :: tmp_permuted_expression
        real(real64), dimension(:, :), allocatable :: tmp_gene_means
        real(real64), dimension(:, :, :), allocatable :: tmp_residuals
        integer(int32), dimension(:, :), allocatable :: tmp_gene_means_perms
        integer(int32), dimension(:), allocatable :: tmp_gene_means_perm_all
        real(real64), dimension(:), allocatable :: tmp_x_star
        integer(int32), dimension(:, :, :), allocatable :: tmp_neighborhood_indices
        integer(int32), dimension(:, :, :), allocatable :: tmp_neighborhood_range
        real(real64), dimension(:, :, :), allocatable :: tmp_neighborhood_residuals_gathered
        integer(int32), dimension(:), allocatable :: tmp_point_pool
        integer(int32), dimension(:), allocatable :: tmp_point_draw
        integer(int32), dimension(:), allocatable :: tmp_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_high
        logical(c_bool), dimension(:), allocatable :: tmp_occupancy_failed
        integer(int32), dimension(:), allocatable :: tmp_n_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_min_bin_occupancy
        real(real64), dimension(:), allocatable :: tmp_mean_bin_occupancy
        integer(int32), dimension(:), allocatable :: tmp_max_bin_occupancy
        integer(int32), dimension(:), allocatable :: tmp_sturges_bins
        integer(int32), dimension(:), allocatable :: tmp_fd_bins
        real(real64), dimension(:, :, :), allocatable :: tmp_pmfs
        integer(int32), dimension(:, :, :), allocatable :: tmp_counts
        integer(int32), dimension(:, :), allocatable :: tmp_included_n_reps
        real(real64), dimension(:, :), allocatable :: tmp_mean_pmf
        integer(int32), dimension(:, :), allocatable :: tmp_mean_pmf_counts
        integer(int32), dimension(:), allocatable :: tmp_mean_pmf_included_n_reps
        real(real64), dimension(:, :), allocatable :: tmp_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_weights
        real(real64), dimension(:), allocatable :: tmp_global_js_divergence
        real(real64), dimension(:), allocatable :: tmp_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_residual_to_study_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts_search
        integer(int32), dimension(:), allocatable :: tmp_study_block_ends
        integer(int32), dimension(:), allocatable :: tmp_study_block_sizes
        integer(int32), dimension(:), allocatable :: tmp_point_n_neighbors
        real(real64), dimension(:), allocatable :: tmp_p_values
        real(real64), dimension(:), allocatable :: tmp_p_values_observed_consensus
        real(real64), dimension(:, :), allocatable :: tmp_permutation_mean_pmf
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_mean_pmf_counts
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_counts
        real(real64), dimension(:, :, :), allocatable :: tmp_permutation_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_permutation_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_permutation_weights
        integer(int32), dimension(:), allocatable :: tmp_search_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_search_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_search_shared_residual_range_high
        real(real64), dimension(:, :), allocatable :: tmp_search_confidence_interval
        integer(int32), dimension(:), allocatable :: tmp_search_trace_n_points
        integer(int32), dimension(:), allocatable :: tmp_search_trace_n_neighbors
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_global_js_divergence
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_lower
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_upper
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_width
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_width_relative
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_overlap
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_delta
        real(real64), dimension(:), allocatable :: tmp_search_trace_delta_median
        real(real64), dimension(:), allocatable :: tmp_search_trace_delta_max
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_selected_n_bins
        logical(c_bool), dimension(:, :), allocatable :: tmp_search_trace_occupancy_failed
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_n_pooled_residuals
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_min_bin_occupancy
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_mean_bin_occupancy
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_max_bin_occupancy
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_sturges_bins
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_fd_bins
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_shared_residual_range_low
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_shared_residual_range_high
        real(real64), dimension(:, :), allocatable :: tmp_candidate_confidence_interval
        real(real64), dimension(:, :, :), allocatable :: tmp_bootstrapping_top_k_jsds
        real(real64), dimension(:), allocatable :: tmp_prev_global_js_divergence
        real(real64), dimension(:, :), allocatable :: tmp_best_uncertainty_confidence_interval
        integer(int32), dimension(:), allocatable :: tmp_best_n_bins_per_point
        integer(int32), dimension(:), allocatable :: tmp_best_uncertainty_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_best_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_best_shared_residual_range_high
        real(real64), dimension(:), allocatable :: tmp_best_uncertainty_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_best_uncertainty_shared_residual_range_high

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_samples, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=7_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=8_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=13_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=15_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=16_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=17_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=19_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=20_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=21_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=22_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=23_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=24_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=25_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=26_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=27_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=28_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(delta_abs_tolerance, ierr, arg_pos=32_int32, min=0.0_real64)
        call validate_all_in_range_real(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(sample_assignment, n_samples, ierr, arg_pos=6_int32, min=1_int32, max=n_samples)
        if (present(join_method)) then; if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=14_int32); end if
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=18_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=30_int32); end if
        if (present(neighborhood_size_mode)) then; if (neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_LARGE .and. neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_SMALL .and. neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=31_int32); end if
        if (present(bootstrap_mode)) then; if (bootstrap_mode /= MODE_BOOTSTRAP_PER_STUDY .and. bootstrap_mode /= MODE_BOOTSTRAP_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=33_int32); end if
        if (present(divergence_mode)) then; if (divergence_mode /= MODE_DIVERGENCE_SYMMETRIC .and. divergence_mode /= MODE_DIVERGENCE_PAIRWISE_K2 .and. divergence_mode /= MODE_DIVERGENCE_ONE_SIDED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=34_int32); end if
        if (present(divergence_weight_mode)) then; if (divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL .and. divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_POOL) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=35_int32); end if
        if (is_err(ierr)) return
#endif

        call calc_js_comp_test_permutation_bounds(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            ierr = ierr&
        )
        if (is_err(ierr)) return
        M_ALLOCATE(tmp_sample_slots(max_n_reps_all_studies*n_studies))
        M_ALLOCATE(tmp_sample_mark(n_samples))
        M_ALLOCATE(tmp_permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_gene_means(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_gene_means_perms(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_x_star(max_n_points_candidate))
        M_ALLOCATE(tmp_neighborhood_indices(max_n_neighbors_candidate, max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_neighborhood_range(2, max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate))
        M_ALLOCATE(tmp_point_pool(256))
        M_ALLOCATE(tmp_point_draw(256))
        M_ALLOCATE(tmp_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_occupancy_failed(max_n_points_candidate))
        M_ALLOCATE(tmp_n_pooled_residuals(max_n_points_candidate))
        M_ALLOCATE(tmp_min_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_mean_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_max_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_sturges_bins(max_n_points_candidate))
        M_ALLOCATE(tmp_fd_bins(max_n_points_candidate))
        M_ALLOCATE(tmp_pmfs(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_counts(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_included_n_reps(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_mean_pmf(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_mean_pmf_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        M_ALLOCATE(tmp_js_divergences(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_weights(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_ALLOCATE(tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
        M_ALLOCATE(tmp_study_block_ends(n_studies))
        M_ALLOCATE(tmp_study_block_sizes(n_studies))
        M_ALLOCATE(tmp_point_n_neighbors(n_studies))
        M_ALLOCATE(tmp_p_values(n_studies))
        M_ALLOCATE(tmp_p_values_observed_consensus(n_studies))
        M_ALLOCATE(tmp_permutation_mean_pmf(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_mean_pmf_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_pmfs(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_permutation_js_divergences(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_permutation_weights(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_search_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_search_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_search_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_search_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_search_trace_n_points(16))
        M_ALLOCATE(tmp_search_trace_n_neighbors(16))
        M_ALLOCATE(tmp_search_trace_global_js_divergence(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_lower(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_upper(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_width(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_width_relative(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_overlap(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_delta(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_delta_median(16))
        M_ALLOCATE(tmp_search_trace_delta_max(16))
        M_ALLOCATE(tmp_search_trace_selected_n_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_occupancy_failed(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_n_pooled_residuals(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_min_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_mean_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_max_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_sturges_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_fd_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_shared_residual_range_low(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_shared_residual_range_high(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_candidate_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        M_ALLOCATE(tmp_prev_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_best_uncertainty_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_best_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_best_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_best_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))

        call run_js_comp_test_on_sample_assignment_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            expression = expression,&
            n_samples = n_samples,&
            sample_assignment = sample_assignment,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            global_js_divergence = global_js_divergence,&
            used_n_points = used_n_points,&
            used_n_neighbors = used_n_neighbors,&
            tmp_sample_slots = tmp_sample_slots,&
            tmp_sample_mark = tmp_sample_mark,&
            tmp_permuted_expression = tmp_permuted_expression,&
            tmp_gene_means = tmp_gene_means,&
            tmp_residuals = tmp_residuals,&
            tmp_gene_means_perms = tmp_gene_means_perms,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_x_star = tmp_x_star,&
            tmp_neighborhood_indices = tmp_neighborhood_indices,&
            tmp_neighborhood_range = tmp_neighborhood_range,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            tmp_n_bins_per_point = tmp_n_bins_per_point,&
            tmp_shared_residual_range_low = tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high = tmp_shared_residual_range_high,&
            tmp_occupancy_failed = tmp_occupancy_failed,&
            tmp_n_pooled_residuals = tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy = tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy = tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy = tmp_max_bin_occupancy,&
            tmp_sturges_bins = tmp_sturges_bins,&
            tmp_fd_bins = tmp_fd_bins,&
            tmp_pmfs = tmp_pmfs,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps = tmp_included_n_reps,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps = tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_point_n_neighbors = tmp_point_n_neighbors,&
            tmp_p_values = tmp_p_values,&
            tmp_p_values_observed_consensus = tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_search_n_bins_per_point = tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low = tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high = tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval = tmp_search_confidence_interval,&
            tmp_search_trace_n_points = tmp_search_trace_n_points,&
            tmp_search_trace_n_neighbors = tmp_search_trace_n_neighbors,&
            tmp_search_trace_global_js_divergence = tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower = tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper = tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width = tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative = tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap = tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta = tmp_search_trace_delta,&
            tmp_search_trace_delta_median = tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max = tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins = tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed = tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals = tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy = tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy = tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy = tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins = tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins = tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low = tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high = tmp_search_trace_shared_residual_range_high,&
            tmp_candidate_confidence_interval = tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence = tmp_prev_global_js_divergence,&
            tmp_best_uncertainty_confidence_interval = tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point = tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point = tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low = tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high = tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low = tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high = tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search = rerun_parameter_search,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            bin_count_mode = bin_count_mode,&
            neighborhood_size_mode = neighborhood_size_mode,&
            delta_abs_tolerance = delta_abs_tolerance,&
            bootstrap_mode = bootstrap_mode,&
            divergence_mode = divergence_mode,&
            divergence_weight_mode = divergence_weight_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_on_sample_assignment

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test_permutation(module):run_js_comp_test_on_sample_assignment]] does both.
    !| One building block of the study-label Permutation-Test: relabels the samples by
    !| `sample_assignment` and runs the whole JS-Comp-Test on the relabelled data set.
    !|
    !| **Input layout.** `expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)`
    !| holds every study: gene `g` is the same gene in every study, a gene missing in a study is NaN
    !| there, and a replicate slot whose values are all NaN is padding (a study with fewer samples).
    !| Every other slot is a sample, numbered in slot order: replicate fastest, then study (sample 1
    !| is the first non-padding slot of study 1).
    !|
    !| **Assignment.** `sample_assignment` is a permutation of `1..n_samples`: destination sample
    !| `k` (the `k`-th non-padding slot) receives all genes of source sample `sample_assignment(k)`.
    !| So every study keeps its slots, its sample count and its padding; only which samples sit in
    !| them changes. The identity `1..n_samples` is the observed data.
    !|
    !| **Steps.** Per study, the gene means and residuals of the relabelled data
    !| ([[tox_data_integration_preprocessing_impl(module):compute_gene_means_impl(interface)]],
    !| [[tox_data_integration_preprocessing_impl(module):compute_residuals_impl(interface)]]); the
    !| reference points from the pooled means
    !| ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]] on their
    !| ascending, NaN-last order); then
    !| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]] with
    !| `n_permutations = 0` (its own residual-level permutation test does not run) and the bin-count
    !| options as given. On the identity assignment the result is exactly that of a caller doing
    !| these steps on the observed data with `run_js_comp_test`.
    !|
    !| **Rerunning the parameter search.** With `rerun_parameter_search = .true.`, the relabelled
    !| data set first runs
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| with the search options below, and the run then uses the `(n_points, n_neighbors)` it returns
    !| (also its fallback when no candidate was admissible); `n_points` and `n_neighbors` are then
    !| not used for the run, only to size the work arrays. `used_n_points` and `used_n_neighbors`
    !| report what the run used. Run time: every call runs the whole JS-Comp-Test; with the rerun
    !| option it runs the whole parameter search first.
    !|
    !| **Known limitation.** A gene missing in a whole study (NaN in all its samples) travels with
    !| the samples, so after relabelling that study may have values for it.
    !|
    !| Errors (`ierr`): the call fails with invalid input (`ERR_INVALID_INPUT`) when some study has
    !| no sample at all, when `n_samples` is not the number of samples in `expression`, or when
    !| `sample_assignment` is not a permutation of `1..n_samples`. (Fortran expert tier only, where
    !| the caller supplies the work-array capacities: with the rerun option it also fails when a
    !| capacity is below the search's own bound, `calc_js_comp_test_candidate_bounds`.) Errors of
    !| the search and the run are returned with their code (for example `ERR_INVALID_INPUT` when some
    !| study has no residual at some reference point of the relabelled data). `global_js_divergence`
    !| is then NaN.
    subroutine run_js_comp_test_on_sample_assignment_expert(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            expression,&
            n_samples,&
            sample_assignment,&
            n_points,&
            n_neighbors,&
            max_n_points_candidate,&
            max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds,&
            global_js_divergence,&
            used_n_points,&
            used_n_neighbors,&
            tmp_sample_slots,&
            tmp_sample_mark,&
            tmp_permuted_expression,&
            tmp_gene_means,&
            tmp_residuals,&
            tmp_gene_means_perms,&
            tmp_gene_means_perm_all,&
            tmp_x_star,&
            tmp_neighborhood_indices,&
            tmp_neighborhood_range,&
            tmp_neighborhood_residuals_gathered,&
            tmp_point_pool,&
            tmp_point_draw,&
            tmp_n_bins_per_point,&
            tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high,&
            tmp_occupancy_failed,&
            tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy,&
            tmp_sturges_bins,&
            tmp_fd_bins,&
            tmp_pmfs,&
            tmp_counts,&
            tmp_included_n_reps,&
            tmp_mean_pmf,&
            tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences,&
            tmp_weights,&
            tmp_global_js_divergence,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm,&
            tmp_bin_counts_search,&
            tmp_study_block_ends,&
            tmp_study_block_sizes,&
            tmp_point_n_neighbors,&
            tmp_p_values,&
            tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts,&
            tmp_permutation_pmfs,&
            tmp_permutation_js_divergences,&
            tmp_permutation_weights,&
            tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval,&
            tmp_search_trace_n_points,&
            tmp_search_trace_n_neighbors,&
            tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta,&
            tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high,&
            tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence,&
            tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search,&
            n_bootstraps,&
            join_method,&
            min_residuals_per_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            bin_count_mode,&
            neighborhood_size_mode,&
            delta_abs_tolerance,&
            bootstrap_mode,&
            divergence_mode,&
            divergence_weight_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_samples
            !! Number of samples, the non-padding slots of `expression`
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_permutation_impl(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `n_points`.
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Capacity of the per-neighbor work arrays
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_permutation_impl(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `n_neighbors`.
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! It is *VERY IMPORTANT* to compute this argument from the `n_bootstrapping_top_k_jsds` output produced by [[tox_data_integration_js_comp_test_permutation_impl(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
            !! A permutation of `1..n_samples`: destination sample `k` receives source sample
            !! `sample_assignment(k)`; the identity is the observed data
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_samples`.
        integer(int32), intent(in) :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus, for this assignment; NaN on error
        integer(int32), intent(out) :: used_n_points
            !! Number of reference points the run used (`n_points`, or the search's choice)
        integer(int32), intent(out) :: used_n_neighbors
            !! Number of neighbors per neighborhood the run used (`n_neighbors`, or the search's choice)
        integer(int32), dimension(max_n_reps_all_studies*n_studies), intent(out) :: tmp_sample_slots
            !! Working array: the flat slot index of every sample
        logical(c_bool), dimension(n_samples), intent(out) :: tmp_sample_mark
            !! Working array: the values of `sample_assignment` already seen
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_permuted_expression
            !! Working array: the relabelled data set
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means
            !! Working array: the relabelled data set's gene means per study
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_residuals
            !! Working array: the relabelled data set's residuals
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means_perms
            !! Working array of the search and the run: each study's sorting permutation of its gene means
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array of the search and the run: sorting permutation of the pooled gene means
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Working array of the search and the run: the reference points
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: tmp_neighborhood_indices
            !! Working array of the search and the run: every study's neighbor gene indices
        integer(int32), dimension(2, max_n_points_candidate, n_studies), intent(out) :: tmp_neighborhood_range
            !! Working array of the search and the run: neighborhood spans
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array of the search and the run: one study's gathered neighborhood residuals
        integer(int32), dimension(256), intent(out) :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(int32), dimension(256), intent(out) :: tmp_point_draw
            !! Working array of the search and the run: one study's draw at one point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array of the search and the run: bin count per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array of the search and the run: lower residual-range bound per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array of the search and the run: upper residual-range bound per reference point
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
            !! Working array of the search and the run: bin-count search failure flag per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
            !! Working array of the search and the run: pooled residual count per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
            !! Working array of the search and the run: minimum bin occupancy per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
            !! Working array of the search and the run: mean bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
            !! Working array of the search and the run: maximum bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
            !! Working array of the search and the run: Sturges' rule diagnostic per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_fd_bins
            !! Working array of the search and the run: Freedman-Diaconis diagnostic per reference point
        real(real64), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(int32), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(real64), dimension(max_n_points_candidate, 256), intent(out) :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        integer(int32), dimension(n_studies), intent(out) :: tmp_point_n_neighbors
            !! Working array of the search and the run: per-study neighbor counts, then sort scratch
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values
            !! Working array: the run's p-values (NaN, its permutation test does not run)
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values_observed_consensus
            !! Working array: the run's comparison p-values (NaN, its permutation test does not run)
        real(real64), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_mean_pmf
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_counts
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_search_n_bins_per_point
            !! Working array: the search's chosen bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_low
            !! Working array: the search's chosen lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_high
            !! Working array: the search's chosen upper residual-range bounds
        real(real64), dimension(2, n_studies), intent(out) :: tmp_search_confidence_interval
            !! Working array: the search's chosen confidence interval
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_n_neighbors
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(real64), dimension(16), intent(out) :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(real64), dimension(16), intent(out) :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_shared_residual_range_high
            !! Working array: the search's trace
        real(real64), dimension(2, n_studies), intent(out) :: tmp_candidate_confidence_interval
            !! Working array of the search: the current candidate's confidence interval
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array of the search: the bootstrap's bottom-k/top-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array of the search: the previous admissible candidate's global JSD
        real(real64), dimension(2, n_studies), intent(out) :: tmp_best_uncertainty_confidence_interval
            !! Working array of the search: the smallest-uncertainty candidate's confidence interval
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_n_bins_per_point
            !! Working array of the search: the best candidate's bin counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_n_bins_per_point
            !! Working array of the search: the smallest-uncertainty candidate's bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_low
            !! Working array of the search: the best candidate's lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_high
            !! Working array of the search: the best candidate's upper residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array of the search: the smallest-uncertainty candidate's lower bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array of the search: the smallest-uncertainty candidate's upper bounds
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! `.true.`: run the fixed-k parameter search on the relabelled data first and use its
            !! `(n_points, n_neighbors)`; `.false.`: use `n_points` and `n_neighbors`
            !! The default value is `.false.`.
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `1_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Rerun search: seed of the bootstrap's random number generator
            !! The default value is `42_int32`.
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
        integer(int32), intent(in), optional :: neighborhood_size_mode
            !! Rerun search: which neighborhood sizes the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are left out
            !! by default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode                                                       | Value                                                                                                    |
            !! |------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Large neighborhoods only, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods only, `KX_FACTORS(2)` = 0.5            | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! | Both, large then small for each `n_points`                 | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_BOTH(variable)]]  |
            !! The default value is `0_int32`.
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0_real64`.
        integer(int32), intent(in), optional :: bootstrap_mode
            !! Rerun search (used only with `rerun_parameter_search`): what its bootstrap draws each
            !! study's histograms from. By default each study's own counts: a robustness interval of
            !! the observed JSD. The pooled counts of all studies give a noise band of "no difference
            !! between the studies" instead, kept for comparison
            !!
            !! | Mode                                                                                             | Value                                                                          |
            !! |--------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
            !! | Each study's own counts, a robustness interval (default)                                         | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_PER_STUDY(variable)]] |
            !! | The pooled consensus counts of all studies, a noise band (the bootstrap before the mode existed) | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_POOLED(variable)]]    |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_mode
            !! Experiment option: how each study's per-point divergence is computed. By default the
            !! JSD of the study against the consensus pmf, the formula before the option existed
            !!
            !! | Mode                                                                                                                | Value                                                                           |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | JSD of each study against the consensus pmf (default)                                                               | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | For exactly two studies the pairwise JSD of the two, the same value for both (otherwise as the default)             | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! | One-sided Kullback-Leibler divergence of each study against the consensus pmf, at most log2(number of studies) bits | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_ONE_SIDED(variable)]]   |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_weight_mode
            !! Experiment option: how each reference point's weight in a study's global divergence is
            !! computed. By default the study's residual count there plus the count of all studies
            !! pooled, the weights before the option existed. For two studies in the pairwise
            !! `divergence_mode` the weights are always the two studies' counts
            !!
            !! | Mode                                                                   | Value                                                                                      |
            !! |------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
            !! | The study's count plus the pooled count of all studies (default)       | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL(variable)]] |
            !! | The pooled count of all studies only, the same weights for every study | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_POOL(variable)]]            |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_samples, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=7_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=8_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_in_range_int(max_n_points_candidate, ierr, arg_pos=9_int32, min=n_points)
        call validate_in_range_int(max_n_neighbors_candidate, ierr, arg_pos=10_int32, min=n_neighbors)
        call validate_in_range_int(n_bootstrapping_top_k_jsds, ierr, arg_pos=11_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=98_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=100_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=101_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=102_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=104_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=105_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=106_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=107_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=108_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=109_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=110_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=111_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=112_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=113_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(delta_abs_tolerance, ierr, arg_pos=117_int32, min=0.0_real64)
        call validate_all_in_range_real(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(sample_assignment, n_samples, ierr, arg_pos=6_int32, min=1_int32, max=n_samples)
        if (present(join_method)) then; if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=99_int32); end if
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=103_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=115_int32); end if
        if (present(neighborhood_size_mode)) then; if (neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_LARGE .and. neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_SMALL .and. neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=116_int32); end if
        if (present(bootstrap_mode)) then; if (bootstrap_mode /= MODE_BOOTSTRAP_PER_STUDY .and. bootstrap_mode /= MODE_BOOTSTRAP_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=118_int32); end if
        if (present(divergence_mode)) then; if (divergence_mode /= MODE_DIVERGENCE_SYMMETRIC .and. divergence_mode /= MODE_DIVERGENCE_PAIRWISE_K2 .and. divergence_mode /= MODE_DIVERGENCE_ONE_SIDED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=119_int32); end if
        if (present(divergence_weight_mode)) then; if (divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL .and. divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_POOL) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=120_int32); end if
        if (is_err(ierr)) return
#endif

        call run_js_comp_test_on_sample_assignment_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            expression = expression,&
            n_samples = n_samples,&
            sample_assignment = sample_assignment,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            global_js_divergence = global_js_divergence,&
            used_n_points = used_n_points,&
            used_n_neighbors = used_n_neighbors,&
            tmp_sample_slots = tmp_sample_slots,&
            tmp_sample_mark = tmp_sample_mark,&
            tmp_permuted_expression = tmp_permuted_expression,&
            tmp_gene_means = tmp_gene_means,&
            tmp_residuals = tmp_residuals,&
            tmp_gene_means_perms = tmp_gene_means_perms,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_x_star = tmp_x_star,&
            tmp_neighborhood_indices = tmp_neighborhood_indices,&
            tmp_neighborhood_range = tmp_neighborhood_range,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            tmp_n_bins_per_point = tmp_n_bins_per_point,&
            tmp_shared_residual_range_low = tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high = tmp_shared_residual_range_high,&
            tmp_occupancy_failed = tmp_occupancy_failed,&
            tmp_n_pooled_residuals = tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy = tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy = tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy = tmp_max_bin_occupancy,&
            tmp_sturges_bins = tmp_sturges_bins,&
            tmp_fd_bins = tmp_fd_bins,&
            tmp_pmfs = tmp_pmfs,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps = tmp_included_n_reps,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps = tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_point_n_neighbors = tmp_point_n_neighbors,&
            tmp_p_values = tmp_p_values,&
            tmp_p_values_observed_consensus = tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_search_n_bins_per_point = tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low = tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high = tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval = tmp_search_confidence_interval,&
            tmp_search_trace_n_points = tmp_search_trace_n_points,&
            tmp_search_trace_n_neighbors = tmp_search_trace_n_neighbors,&
            tmp_search_trace_global_js_divergence = tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower = tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper = tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width = tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative = tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap = tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta = tmp_search_trace_delta,&
            tmp_search_trace_delta_median = tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max = tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins = tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed = tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals = tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy = tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy = tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy = tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins = tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins = tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low = tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high = tmp_search_trace_shared_residual_range_high,&
            tmp_candidate_confidence_interval = tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence = tmp_prev_global_js_divergence,&
            tmp_best_uncertainty_confidence_interval = tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point = tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point = tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low = tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high = tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low = tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high = tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search = rerun_parameter_search,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            bin_count_mode = bin_count_mode,&
            neighborhood_size_mode = neighborhood_size_mode,&
            delta_abs_tolerance = delta_abs_tolerance,&
            bootstrap_mode = bootstrap_mode,&
            divergence_mode = divergence_mode,&
            divergence_weight_mode = divergence_weight_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_on_sample_assignment_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test_permutation(module):run_js_comp_test_adaptive_on_sample_assignment_expert]] to prepare it yourself.
    !| The adaptive counterpart of
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]],
    !| with the same input layout and assignment convention: `expression` holds every study (NaN for
    !| a missing gene, all-NaN replicate slots are padding, samples numbered replicate fastest, then
    !| study), and destination sample `k` receives all genes of source sample `sample_assignment(k)`;
    !| the identity is the observed data.
    !|
    !| **Steps.** Per study, the gene means and residuals of the relabelled data; the adaptive
    !| neighborhoods for the growth knobs `k_start`, `k_step`, `k_max` and the construction options
    !| ([[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]);
    !| then [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_adaptive_impl(interface)]]
    !| on its ranges with `n_permutations = 0`. `n_points` is the number of neighborhoods the
    !| construction built. On the identity assignment the result is exactly that of a caller doing
    !| these steps on the observed data.
    !|
    !| **Rerunning the parameter search.** With `rerun_parameter_search = .true.`, the relabelled
    !| data first runs
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]]
    !| (with its own capacity, `calc_adaptive_js_comp_test_bounds`, so it rejects the same candidates
    !| as a direct call), and the construction and run then use the knobs it returns (also its
    !| fallback when no candidate was admissible); `k_start`, `k_step` and `k_max` are then not used,
    !| except `k_start` to size the work arrays. `used_k_start`, `used_k_step` and `used_k_max` report
    !| what was used. Run time: every call runs the whole JS-Comp-Test; with the rerun option it runs
    !| the whole parameter search first.
    !|
    !| **Known limitation.** A gene missing in a whole study (NaN in all its samples) travels with
    !| the samples, so after relabelling that study may have values for it.
    !|
    !| Errors (`ierr`): the call fails with invalid input (`ERR_INVALID_INPUT`) when some study has
    !| no sample at all, when `n_samples` is not the number of samples in `expression`, when
    !| `sample_assignment` is not a permutation of `1..n_samples`, or when the construction builds no
    !| neighborhood or more neighborhoods than the work arrays hold. (Fortran expert tier only, where
    !| the caller supplies the work-array capacity: with the rerun option it also fails when the
    !| capacity is below the search's own.) Errors of the search and the run are returned with their
    !| code. `global_js_divergence` is then NaN.
    subroutine run_js_comp_test_adaptive_on_sample_assignment(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            expression,&
            n_samples,&
            sample_assignment,&
            k_start,&
            k_step,&
            k_max,&
            global_js_divergence,&
            used_k_start,&
            used_k_step,&
            used_k_max,&
            n_points,&
            rerun_parameter_search,&
            n_bootstraps,&
            join_method,&
            min_residuals_per_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            tau,&
            mad_distance_factor,&
            max_pooled_residuals,&
            min_study_neighbors,&
            bin_count_mode,&
            delta_abs_tolerance,&
            bootstrap_mode,&
            divergence_mode,&
            divergence_weight_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_samples
            !! Number of samples, the non-padding slots of `expression`
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
            !! A permutation of `1..n_samples`: destination sample `k` receives source sample
            !! `sample_assignment(k)`; the identity is the observed data
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_samples`.
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! The minimum valid value is `k_start`.
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus, for this assignment; NaN on error
        integer(int32), intent(out) :: used_k_start
            !! `k_start` the construction used (`k_start`, or the search's choice)
        integer(int32), intent(out) :: used_k_step
            !! `k_step` the construction used (`k_step`, or the search's choice)
        integer(int32), intent(out) :: used_k_max
            !! `k_max` the construction used (`k_max`, or the search's choice)
        integer(int32), intent(out) :: n_points
            !! Number of neighborhoods the construction built (0 when it built none)
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! `.true.`: run the adaptive parameter search on the relabelled data first and use its
            !! knobs; `.false.`: use `k_start`, `k_step` and `k_max`
            !! The default value is `.false.`.
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `1_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Rerun search: seed of the bootstrap's random number generator
            !! The default value is `42_int32`.
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have (construction status; in
            !! the rerun search, candidate admissibility)
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
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
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0_real64`.
        integer(int32), intent(in), optional :: bootstrap_mode
            !! Rerun search (used only with `rerun_parameter_search`): what its bootstrap draws each
            !! study's histograms from. By default each study's own counts: a robustness interval of
            !! the observed JSD. The pooled counts of all studies give a noise band of "no difference
            !! between the studies" instead, kept for comparison
            !!
            !! | Mode                                                                                             | Value                                                                          |
            !! |--------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
            !! | Each study's own counts, a robustness interval (default)                                         | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_PER_STUDY(variable)]] |
            !! | The pooled consensus counts of all studies, a noise band (the bootstrap before the mode existed) | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_POOLED(variable)]]    |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_mode
            !! Experiment option: how each study's per-point divergence is computed. By default the
            !! JSD of the study against the consensus pmf, the formula before the option existed
            !!
            !! | Mode                                                                                                                | Value                                                                           |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | JSD of each study against the consensus pmf (default)                                                               | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | For exactly two studies the pairwise JSD of the two, the same value for both (otherwise as the default)             | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! | One-sided Kullback-Leibler divergence of each study against the consensus pmf, at most log2(number of studies) bits | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_ONE_SIDED(variable)]]   |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_weight_mode
            !! Experiment option: how each reference point's weight in a study's global divergence is
            !! computed. By default the study's residual count there plus the count of all studies
            !! pooled, the weights before the option existed. For two studies in the pairwise
            !! `divergence_mode` the weights are always the two studies' counts
            !!
            !! | Mode                                                                   | Value                                                                                      |
            !! |------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
            !! | The study's count plus the pooled count of all studies (default)       | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL(variable)]] |
            !! | The pooled count of all studies only, the same weights for every study | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_POOL(variable)]]            |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)
        integer(int32) :: max_n_points_candidate
        integer(int32) :: n_bootstrapping_top_k_jsds
        integer(int32), dimension(:), allocatable :: tmp_sample_slots
        logical(c_bool), dimension(:), allocatable :: tmp_sample_mark
        real(real64), dimension(:, :, :), allocatable :: tmp_permuted_expression
        real(real64), dimension(:, :), allocatable :: tmp_gene_means
        real(real64), dimension(:, :, :), allocatable :: tmp_residuals
        integer(int32), dimension(:), allocatable :: tmp_gene_means_perm_all
        real(real64), dimension(:), allocatable :: tmp_x_star
        integer(int32), dimension(:, :), allocatable :: tmp_pooled_neighborhood_range
        integer(int32), dimension(:, :), allocatable :: tmp_n_neighbors_per_point
        integer(int32), dimension(:), allocatable :: tmp_stop_reason
        real(real64), dimension(:), allocatable :: tmp_neighborhood_dispersion
        real(real64), dimension(:), allocatable :: tmp_neighborhood_mad
        integer(int32), dimension(:, :), allocatable :: tmp_point_neighborhood_indices
        real(real64), dimension(:, :), allocatable :: tmp_neighbor_residuals
        integer(int32), dimension(:), allocatable :: tmp_point_pool
        integer(int32), dimension(:), allocatable :: tmp_point_draw
        integer(int32), dimension(:), allocatable :: tmp_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_high
        logical(c_bool), dimension(:), allocatable :: tmp_occupancy_failed
        integer(int32), dimension(:), allocatable :: tmp_n_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_min_bin_occupancy
        real(real64), dimension(:), allocatable :: tmp_mean_bin_occupancy
        integer(int32), dimension(:), allocatable :: tmp_max_bin_occupancy
        integer(int32), dimension(:), allocatable :: tmp_sturges_bins
        integer(int32), dimension(:), allocatable :: tmp_fd_bins
        real(real64), dimension(:, :, :), allocatable :: tmp_pmfs
        integer(int32), dimension(:, :, :), allocatable :: tmp_counts
        integer(int32), dimension(:, :), allocatable :: tmp_included_n_reps
        real(real64), dimension(:, :), allocatable :: tmp_mean_pmf
        integer(int32), dimension(:, :), allocatable :: tmp_mean_pmf_counts
        integer(int32), dimension(:), allocatable :: tmp_mean_pmf_included_n_reps
        real(real64), dimension(:, :), allocatable :: tmp_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_weights
        real(real64), dimension(:), allocatable :: tmp_global_js_divergence
        real(real64), dimension(:), allocatable :: tmp_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_residual_to_study_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts_search
        integer(int32), dimension(:), allocatable :: tmp_study_block_ends
        integer(int32), dimension(:), allocatable :: tmp_study_block_sizes
        real(real64), dimension(:), allocatable :: tmp_p_values
        real(real64), dimension(:), allocatable :: tmp_p_values_observed_consensus
        real(real64), dimension(:, :), allocatable :: tmp_permutation_mean_pmf
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_mean_pmf_counts
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_counts
        real(real64), dimension(:, :, :), allocatable :: tmp_permutation_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_permutation_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_permutation_weights
        integer(int32), dimension(:), allocatable :: tmp_search_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_search_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_search_shared_residual_range_high
        real(real64), dimension(:, :), allocatable :: tmp_search_confidence_interval
        integer(int32), dimension(:), allocatable :: tmp_search_trace_k_start
        integer(int32), dimension(:), allocatable :: tmp_search_trace_k_step
        integer(int32), dimension(:), allocatable :: tmp_search_trace_k_max
        integer(int32), dimension(:), allocatable :: tmp_search_trace_n_points
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_global_js_divergence
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_lower
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_upper
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_width
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_width_relative
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_overlap
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_delta
        real(real64), dimension(:), allocatable :: tmp_search_trace_delta_median
        real(real64), dimension(:), allocatable :: tmp_search_trace_delta_max
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_selected_n_bins
        logical(c_bool), dimension(:, :), allocatable :: tmp_search_trace_occupancy_failed
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_n_pooled_residuals
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_min_bin_occupancy
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_mean_bin_occupancy
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_max_bin_occupancy
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_sturges_bins
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_fd_bins
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_shared_residual_range_low
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_shared_residual_range_high
        integer(int32), dimension(:), allocatable :: tmp_search_candidate_k_start
        integer(int32), dimension(:), allocatable :: tmp_search_candidate_n_points
        integer(int32), dimension(:), allocatable :: tmp_search_candidate_status
        real(real64), dimension(:, :), allocatable :: tmp_candidate_confidence_interval
        real(real64), dimension(:, :, :), allocatable :: tmp_bootstrapping_top_k_jsds
        real(real64), dimension(:), allocatable :: tmp_prev_global_js_divergence
        integer(int32), dimension(:), allocatable :: tmp_delta_perm
        real(real64), dimension(:, :), allocatable :: tmp_best_uncertainty_confidence_interval
        integer(int32), dimension(:), allocatable :: tmp_best_n_bins_per_point
        integer(int32), dimension(:), allocatable :: tmp_best_uncertainty_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_best_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_best_shared_residual_range_high
        real(real64), dimension(:), allocatable :: tmp_best_uncertainty_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_best_uncertainty_shared_residual_range_high

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_samples, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(k_start, ierr, arg_pos=7_int32, min=1_int32)
        call validate_in_range_int(k_step, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_int(k_max, ierr, arg_pos=9_int32, min=k_start)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=16_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=18_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=19_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=20_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=22_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=23_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=24_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=25_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=26_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=27_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=28_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=29_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=30_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=31_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(tau, ierr, arg_pos=33_int32, min=0.0_real64)
        call validate_in_range_real(mad_distance_factor, ierr, arg_pos=34_int32, min=0.0_real64)
        call validate_in_range_int(max_pooled_residuals, ierr, arg_pos=35_int32, min=0_int32)
        call validate_in_range_int(min_study_neighbors, ierr, arg_pos=36_int32, min=1_int32)
        call validate_in_range_real(delta_abs_tolerance, ierr, arg_pos=38_int32, min=0.0_real64)
        call validate_all_in_range_real(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(sample_assignment, n_samples, ierr, arg_pos=6_int32, min=1_int32, max=n_samples)
        if (present(join_method)) then; if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=17_int32); end if
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=21_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=37_int32); end if
        if (present(bootstrap_mode)) then; if (bootstrap_mode /= MODE_BOOTSTRAP_PER_STUDY .and. bootstrap_mode /= MODE_BOOTSTRAP_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=39_int32); end if
        if (present(divergence_mode)) then; if (divergence_mode /= MODE_DIVERGENCE_SYMMETRIC .and. divergence_mode /= MODE_DIVERGENCE_PAIRWISE_K2 .and. divergence_mode /= MODE_DIVERGENCE_ONE_SIDED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=40_int32); end if
        if (present(divergence_weight_mode)) then; if (divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL .and. divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_POOL) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=41_int32); end if
        if (is_err(ierr)) return
#endif

        call calc_adaptive_js_comp_test_permutation_bounds(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_studies = n_studies,&
            k_start = k_start,&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            max_n_points_candidate = max_n_points_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            ierr = ierr&
        )
        if (is_err(ierr)) return
        M_ALLOCATE(tmp_sample_slots(max_n_reps_all_studies*n_studies))
        M_ALLOCATE(tmp_sample_mark(n_samples))
        M_ALLOCATE(tmp_permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_gene_means(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_x_star(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_pooled_neighborhood_range(2, max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_n_neighbors_per_point(n_studies, max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_stop_reason(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_neighborhood_dispersion(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_neighborhood_mad(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_point_neighborhood_indices(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_neighbor_residuals(max_n_reps_all_studies, max_n_genes_all_studies))
        M_ALLOCATE(tmp_point_pool(256))
        M_ALLOCATE(tmp_point_draw(256))
        M_ALLOCATE(tmp_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_occupancy_failed(max_n_points_candidate))
        M_ALLOCATE(tmp_n_pooled_residuals(max_n_points_candidate))
        M_ALLOCATE(tmp_min_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_mean_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_max_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_sturges_bins(max_n_points_candidate))
        M_ALLOCATE(tmp_fd_bins(max_n_points_candidate))
        M_ALLOCATE(tmp_pmfs(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_counts(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_included_n_reps(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_mean_pmf(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_mean_pmf_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        M_ALLOCATE(tmp_js_divergences(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_weights(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
        M_ALLOCATE(tmp_study_block_ends(n_studies))
        M_ALLOCATE(tmp_study_block_sizes(n_studies))
        M_ALLOCATE(tmp_p_values(n_studies))
        M_ALLOCATE(tmp_p_values_observed_consensus(n_studies))
        M_ALLOCATE(tmp_permutation_mean_pmf(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_mean_pmf_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_pmfs(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_permutation_js_divergences(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_permutation_weights(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_search_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_search_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_search_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_search_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_search_trace_k_start(16))
        M_ALLOCATE(tmp_search_trace_k_step(16))
        M_ALLOCATE(tmp_search_trace_k_max(16))
        M_ALLOCATE(tmp_search_trace_n_points(16))
        M_ALLOCATE(tmp_search_trace_global_js_divergence(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_lower(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_upper(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_width(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_width_relative(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_overlap(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_delta(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_delta_median(16))
        M_ALLOCATE(tmp_search_trace_delta_max(16))
        M_ALLOCATE(tmp_search_trace_selected_n_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_occupancy_failed(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_n_pooled_residuals(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_min_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_mean_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_max_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_sturges_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_fd_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_shared_residual_range_low(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_shared_residual_range_high(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_candidate_k_start(16))
        M_ALLOCATE(tmp_search_candidate_n_points(16))
        M_ALLOCATE(tmp_search_candidate_status(16))
        M_ALLOCATE(tmp_candidate_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        M_ALLOCATE(tmp_prev_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_delta_perm(n_studies))
        M_ALLOCATE(tmp_best_uncertainty_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_best_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_best_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_best_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))

        call run_js_comp_test_adaptive_on_sample_assignment_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            expression = expression,&
            n_samples = n_samples,&
            sample_assignment = sample_assignment,&
            k_start = k_start,&
            k_step = k_step,&
            k_max = k_max,&
            max_n_points_candidate = max_n_points_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            global_js_divergence = global_js_divergence,&
            used_k_start = used_k_start,&
            used_k_step = used_k_step,&
            used_k_max = used_k_max,&
            n_points = n_points,&
            tmp_sample_slots = tmp_sample_slots,&
            tmp_sample_mark = tmp_sample_mark,&
            tmp_permuted_expression = tmp_permuted_expression,&
            tmp_gene_means = tmp_gene_means,&
            tmp_residuals = tmp_residuals,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_x_star = tmp_x_star,&
            tmp_pooled_neighborhood_range = tmp_pooled_neighborhood_range,&
            tmp_n_neighbors_per_point = tmp_n_neighbors_per_point,&
            tmp_stop_reason = tmp_stop_reason,&
            tmp_neighborhood_dispersion = tmp_neighborhood_dispersion,&
            tmp_neighborhood_mad = tmp_neighborhood_mad,&
            tmp_point_neighborhood_indices = tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals = tmp_neighbor_residuals,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            tmp_n_bins_per_point = tmp_n_bins_per_point,&
            tmp_shared_residual_range_low = tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high = tmp_shared_residual_range_high,&
            tmp_occupancy_failed = tmp_occupancy_failed,&
            tmp_n_pooled_residuals = tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy = tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy = tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy = tmp_max_bin_occupancy,&
            tmp_sturges_bins = tmp_sturges_bins,&
            tmp_fd_bins = tmp_fd_bins,&
            tmp_pmfs = tmp_pmfs,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps = tmp_included_n_reps,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps = tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_p_values = tmp_p_values,&
            tmp_p_values_observed_consensus = tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_search_n_bins_per_point = tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low = tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high = tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval = tmp_search_confidence_interval,&
            tmp_search_trace_k_start = tmp_search_trace_k_start,&
            tmp_search_trace_k_step = tmp_search_trace_k_step,&
            tmp_search_trace_k_max = tmp_search_trace_k_max,&
            tmp_search_trace_n_points = tmp_search_trace_n_points,&
            tmp_search_trace_global_js_divergence = tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower = tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper = tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width = tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative = tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap = tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta = tmp_search_trace_delta,&
            tmp_search_trace_delta_median = tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max = tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins = tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed = tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals = tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy = tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy = tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy = tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins = tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins = tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low = tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high = tmp_search_trace_shared_residual_range_high,&
            tmp_search_candidate_k_start = tmp_search_candidate_k_start,&
            tmp_search_candidate_n_points = tmp_search_candidate_n_points,&
            tmp_search_candidate_status = tmp_search_candidate_status,&
            tmp_candidate_confidence_interval = tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence = tmp_prev_global_js_divergence,&
            tmp_delta_perm = tmp_delta_perm,&
            tmp_best_uncertainty_confidence_interval = tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point = tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point = tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low = tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high = tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low = tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high = tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search = rerun_parameter_search,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            bin_count_mode = bin_count_mode,&
            delta_abs_tolerance = delta_abs_tolerance,&
            bootstrap_mode = bootstrap_mode,&
            divergence_mode = divergence_mode,&
            divergence_weight_mode = divergence_weight_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_adaptive_on_sample_assignment

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test_permutation(module):run_js_comp_test_adaptive_on_sample_assignment]] does both.
    !| The adaptive counterpart of
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]],
    !| with the same input layout and assignment convention: `expression` holds every study (NaN for
    !| a missing gene, all-NaN replicate slots are padding, samples numbered replicate fastest, then
    !| study), and destination sample `k` receives all genes of source sample `sample_assignment(k)`;
    !| the identity is the observed data.
    !|
    !| **Steps.** Per study, the gene means and residuals of the relabelled data; the adaptive
    !| neighborhoods for the growth knobs `k_start`, `k_step`, `k_max` and the construction options
    !| ([[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]);
    !| then [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_adaptive_impl(interface)]]
    !| on its ranges with `n_permutations = 0`. `n_points` is the number of neighborhoods the
    !| construction built. On the identity assignment the result is exactly that of a caller doing
    !| these steps on the observed data.
    !|
    !| **Rerunning the parameter search.** With `rerun_parameter_search = .true.`, the relabelled
    !| data first runs
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]]
    !| (with its own capacity, `calc_adaptive_js_comp_test_bounds`, so it rejects the same candidates
    !| as a direct call), and the construction and run then use the knobs it returns (also its
    !| fallback when no candidate was admissible); `k_start`, `k_step` and `k_max` are then not used,
    !| except `k_start` to size the work arrays. `used_k_start`, `used_k_step` and `used_k_max` report
    !| what was used. Run time: every call runs the whole JS-Comp-Test; with the rerun option it runs
    !| the whole parameter search first.
    !|
    !| **Known limitation.** A gene missing in a whole study (NaN in all its samples) travels with
    !| the samples, so after relabelling that study may have values for it.
    !|
    !| Errors (`ierr`): the call fails with invalid input (`ERR_INVALID_INPUT`) when some study has
    !| no sample at all, when `n_samples` is not the number of samples in `expression`, when
    !| `sample_assignment` is not a permutation of `1..n_samples`, or when the construction builds no
    !| neighborhood or more neighborhoods than the work arrays hold. (Fortran expert tier only, where
    !| the caller supplies the work-array capacity: with the rerun option it also fails when the
    !| capacity is below the search's own.) Errors of the search and the run are returned with their
    !| code. `global_js_divergence` is then NaN.
    subroutine run_js_comp_test_adaptive_on_sample_assignment_expert(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            expression,&
            n_samples,&
            sample_assignment,&
            k_start,&
            k_step,&
            k_max,&
            max_n_points_candidate,&
            n_bootstrapping_top_k_jsds,&
            global_js_divergence,&
            used_k_start,&
            used_k_step,&
            used_k_max,&
            n_points,&
            tmp_sample_slots,&
            tmp_sample_mark,&
            tmp_permuted_expression,&
            tmp_gene_means,&
            tmp_residuals,&
            tmp_gene_means_perm_all,&
            tmp_x_star,&
            tmp_pooled_neighborhood_range,&
            tmp_n_neighbors_per_point,&
            tmp_stop_reason,&
            tmp_neighborhood_dispersion,&
            tmp_neighborhood_mad,&
            tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals,&
            tmp_point_pool,&
            tmp_point_draw,&
            tmp_n_bins_per_point,&
            tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high,&
            tmp_occupancy_failed,&
            tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy,&
            tmp_sturges_bins,&
            tmp_fd_bins,&
            tmp_pmfs,&
            tmp_counts,&
            tmp_included_n_reps,&
            tmp_mean_pmf,&
            tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences,&
            tmp_weights,&
            tmp_global_js_divergence,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm,&
            tmp_bin_counts_search,&
            tmp_study_block_ends,&
            tmp_study_block_sizes,&
            tmp_p_values,&
            tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts,&
            tmp_permutation_pmfs,&
            tmp_permutation_js_divergences,&
            tmp_permutation_weights,&
            tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval,&
            tmp_search_trace_k_start,&
            tmp_search_trace_k_step,&
            tmp_search_trace_k_max,&
            tmp_search_trace_n_points,&
            tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta,&
            tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high,&
            tmp_search_candidate_k_start,&
            tmp_search_candidate_n_points,&
            tmp_search_candidate_status,&
            tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence,&
            tmp_delta_perm,&
            tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search,&
            n_bootstraps,&
            join_method,&
            min_residuals_per_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            tau,&
            mad_distance_factor,&
            max_pooled_residuals,&
            min_study_neighbors,&
            bin_count_mode,&
            delta_abs_tolerance,&
            bootstrap_mode,&
            divergence_mode,&
            divergence_weight_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_samples
            !! Number of samples, the non-padding slots of `expression`
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays; a construction with more neighborhoods fails
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_permutation_impl(module):calc_adaptive_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! It is *VERY IMPORTANT* to compute this argument from the `n_bootstrapping_top_k_jsds` output produced by [[tox_data_integration_js_comp_test_permutation_impl(module):calc_adaptive_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
            !! A permutation of `1..n_samples`: destination sample `k` receives source sample
            !! `sample_assignment(k)`; the identity is the observed data
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_samples`.
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! The minimum valid value is `k_start`.
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus, for this assignment; NaN on error
        integer(int32), intent(out) :: used_k_start
            !! `k_start` the construction used (`k_start`, or the search's choice)
        integer(int32), intent(out) :: used_k_step
            !! `k_step` the construction used (`k_step`, or the search's choice)
        integer(int32), intent(out) :: used_k_max
            !! `k_max` the construction used (`k_max`, or the search's choice)
        integer(int32), intent(out) :: n_points
            !! Number of neighborhoods the construction built (0 when it built none)
        integer(int32), dimension(max_n_reps_all_studies*n_studies), intent(out) :: tmp_sample_slots
            !! Working array: the flat slot index of every sample
        logical(c_bool), dimension(n_samples), intent(out) :: tmp_sample_mark
            !! Working array: the values of `sample_assignment` already seen
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_permuted_expression
            !! Working array: the relabelled data set
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means
            !! Working array: the relabelled data set's gene means per study
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_residuals
            !! Working array: the relabelled data set's residuals
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array of the search, the construction and the run: sorting permutation of the
            !! pooled gene means
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_x_star
            !! Working array of the search and the construction: the reference points
        integer(int32), dimension(2, max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_neighborhood_range
            !! Working array of the search and the construction: each neighborhood's pooled range,
            !! read by the run
        integer(int32), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out) :: tmp_n_neighbors_per_point
            !! Working array of the search, the construction and the run: per-study entry count per
            !! neighborhood
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_stop_reason
            !! Working array of the search and the construction: stop reason per neighborhood
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_dispersion
            !! Working array of the search and the construction: dispersion per neighborhood
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_mad
            !! Working array of the search and the construction: MAD per neighborhood
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_point_neighborhood_indices
            !! Working array of the search and the run: one reference point's gene indices per study
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out) :: tmp_neighbor_residuals
            !! Working array of the search and the run: one (point, study) pair's gathered residuals
        integer(int32), dimension(256), intent(out) :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(int32), dimension(256), intent(out) :: tmp_point_draw
            !! Working array of the search and the run: one study's draw at one point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array of the search and the run: bin count per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array of the search and the run: lower residual-range bound per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array of the search and the run: upper residual-range bound per reference point
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
            !! Working array of the search and the run: bin-count search failure flag per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
            !! Working array of the search and the run: pooled residual count per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
            !! Working array of the search and the run: minimum bin occupancy per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
            !! Working array of the search and the run: mean bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
            !! Working array of the search and the run: maximum bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
            !! Working array of the search and the run: Sturges' rule diagnostic per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_fd_bins
            !! Working array of the search and the run: Freedman-Diaconis diagnostic per reference point
        real(real64), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(int32), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(real64), dimension(max_n_points_candidate, 256), intent(out) :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values
            !! Working array: the run's p-values (NaN, its permutation test does not run)
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values_observed_consensus
            !! Working array: the run's comparison p-values (NaN, its permutation test does not run)
        real(real64), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_mean_pmf
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_counts
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_search_n_bins_per_point
            !! Working array: the search's chosen bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_low
            !! Working array: the search's chosen lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_high
            !! Working array: the search's chosen upper residual-range bounds
        real(real64), dimension(2, n_studies), intent(out) :: tmp_search_confidence_interval
            !! Working array: the search's chosen confidence interval
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_k_start
            !! Working array: the search's trace
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_k_step
            !! Working array: the search's trace
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_k_max
            !! Working array: the search's trace
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(real64), dimension(16), intent(out) :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(real64), dimension(16), intent(out) :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_shared_residual_range_high
            !! Working array: the search's trace
        integer(int32), dimension(16), intent(out) :: tmp_search_candidate_k_start
            !! Working array: the search's candidate log
        integer(int32), dimension(16), intent(out) :: tmp_search_candidate_n_points
            !! Working array: the search's candidate log
        integer(int32), dimension(16), intent(out) :: tmp_search_candidate_status
            !! Working array: the search's candidate log
        real(real64), dimension(2, n_studies), intent(out) :: tmp_candidate_confidence_interval
            !! Working array of the search: the current candidate's confidence interval
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array of the search: the bootstrap's bottom-k/top-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array of the search: the previous admissible candidate's global JSD
        integer(int32), dimension(n_studies), intent(out) :: tmp_delta_perm
            !! Working array of the search: sorting permutation for the per-study medians
        real(real64), dimension(2, n_studies), intent(out) :: tmp_best_uncertainty_confidence_interval
            !! Working array of the search: the smallest-uncertainty candidate's confidence interval
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_n_bins_per_point
            !! Working array of the search: the best candidate's bin counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_n_bins_per_point
            !! Working array of the search: the smallest-uncertainty candidate's bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_low
            !! Working array of the search: the best candidate's lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_high
            !! Working array of the search: the best candidate's upper residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array of the search: the smallest-uncertainty candidate's lower bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array of the search: the smallest-uncertainty candidate's upper bounds
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! `.true.`: run the adaptive parameter search on the relabelled data first and use its
            !! knobs; `.false.`: use `k_start`, `k_step` and `k_max`
            !! The default value is `.false.`.
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `1_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Rerun search: seed of the bootstrap's random number generator
            !! The default value is `42_int32`.
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have (construction status; in
            !! the rerun search, candidate admissibility)
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
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
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0_real64`.
        integer(int32), intent(in), optional :: bootstrap_mode
            !! Rerun search (used only with `rerun_parameter_search`): what its bootstrap draws each
            !! study's histograms from. By default each study's own counts: a robustness interval of
            !! the observed JSD. The pooled counts of all studies give a noise band of "no difference
            !! between the studies" instead, kept for comparison
            !!
            !! | Mode                                                                                             | Value                                                                          |
            !! |--------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
            !! | Each study's own counts, a robustness interval (default)                                         | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_PER_STUDY(variable)]] |
            !! | The pooled consensus counts of all studies, a noise band (the bootstrap before the mode existed) | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_POOLED(variable)]]    |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_mode
            !! Experiment option: how each study's per-point divergence is computed. By default the
            !! JSD of the study against the consensus pmf, the formula before the option existed
            !!
            !! | Mode                                                                                                                | Value                                                                           |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | JSD of each study against the consensus pmf (default)                                                               | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | For exactly two studies the pairwise JSD of the two, the same value for both (otherwise as the default)             | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! | One-sided Kullback-Leibler divergence of each study against the consensus pmf, at most log2(number of studies) bits | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_ONE_SIDED(variable)]]   |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_weight_mode
            !! Experiment option: how each reference point's weight in a study's global divergence is
            !! computed. By default the study's residual count there plus the count of all studies
            !! pooled, the weights before the option existed. For two studies in the pairwise
            !! `divergence_mode` the weights are always the two studies' counts
            !!
            !! | Mode                                                                   | Value                                                                                      |
            !! |------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
            !! | The study's count plus the pooled count of all studies (default)       | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL(variable)]] |
            !! | The pooled count of all studies only, the same weights for every study | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_POOL(variable)]]            |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_samples, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(k_start, ierr, arg_pos=7_int32, min=1_int32)
        call validate_in_range_int(k_step, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_int(k_max, ierr, arg_pos=9_int32, min=k_start)
        call validate_in_range_int(max_n_points_candidate, ierr, arg_pos=10_int32, min=1_int32)
        call validate_in_range_int(n_bootstrapping_top_k_jsds, ierr, arg_pos=11_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=108_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=110_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=111_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=112_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=114_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=115_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=116_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=117_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=118_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=119_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=120_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=121_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=122_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=123_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(tau, ierr, arg_pos=125_int32, min=0.0_real64)
        call validate_in_range_real(mad_distance_factor, ierr, arg_pos=126_int32, min=0.0_real64)
        call validate_in_range_int(max_pooled_residuals, ierr, arg_pos=127_int32, min=0_int32)
        call validate_in_range_int(min_study_neighbors, ierr, arg_pos=128_int32, min=1_int32)
        call validate_in_range_real(delta_abs_tolerance, ierr, arg_pos=130_int32, min=0.0_real64)
        call validate_all_in_range_real(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(sample_assignment, n_samples, ierr, arg_pos=6_int32, min=1_int32, max=n_samples)
        if (present(join_method)) then; if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=109_int32); end if
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=113_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=129_int32); end if
        if (present(bootstrap_mode)) then; if (bootstrap_mode /= MODE_BOOTSTRAP_PER_STUDY .and. bootstrap_mode /= MODE_BOOTSTRAP_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=131_int32); end if
        if (present(divergence_mode)) then; if (divergence_mode /= MODE_DIVERGENCE_SYMMETRIC .and. divergence_mode /= MODE_DIVERGENCE_PAIRWISE_K2 .and. divergence_mode /= MODE_DIVERGENCE_ONE_SIDED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=132_int32); end if
        if (present(divergence_weight_mode)) then; if (divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL .and. divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_POOL) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=133_int32); end if
        if (is_err(ierr)) return
#endif

        call run_js_comp_test_adaptive_on_sample_assignment_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            expression = expression,&
            n_samples = n_samples,&
            sample_assignment = sample_assignment,&
            k_start = k_start,&
            k_step = k_step,&
            k_max = k_max,&
            max_n_points_candidate = max_n_points_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            global_js_divergence = global_js_divergence,&
            used_k_start = used_k_start,&
            used_k_step = used_k_step,&
            used_k_max = used_k_max,&
            n_points = n_points,&
            tmp_sample_slots = tmp_sample_slots,&
            tmp_sample_mark = tmp_sample_mark,&
            tmp_permuted_expression = tmp_permuted_expression,&
            tmp_gene_means = tmp_gene_means,&
            tmp_residuals = tmp_residuals,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_x_star = tmp_x_star,&
            tmp_pooled_neighborhood_range = tmp_pooled_neighborhood_range,&
            tmp_n_neighbors_per_point = tmp_n_neighbors_per_point,&
            tmp_stop_reason = tmp_stop_reason,&
            tmp_neighborhood_dispersion = tmp_neighborhood_dispersion,&
            tmp_neighborhood_mad = tmp_neighborhood_mad,&
            tmp_point_neighborhood_indices = tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals = tmp_neighbor_residuals,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            tmp_n_bins_per_point = tmp_n_bins_per_point,&
            tmp_shared_residual_range_low = tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high = tmp_shared_residual_range_high,&
            tmp_occupancy_failed = tmp_occupancy_failed,&
            tmp_n_pooled_residuals = tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy = tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy = tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy = tmp_max_bin_occupancy,&
            tmp_sturges_bins = tmp_sturges_bins,&
            tmp_fd_bins = tmp_fd_bins,&
            tmp_pmfs = tmp_pmfs,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps = tmp_included_n_reps,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps = tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_p_values = tmp_p_values,&
            tmp_p_values_observed_consensus = tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_search_n_bins_per_point = tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low = tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high = tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval = tmp_search_confidence_interval,&
            tmp_search_trace_k_start = tmp_search_trace_k_start,&
            tmp_search_trace_k_step = tmp_search_trace_k_step,&
            tmp_search_trace_k_max = tmp_search_trace_k_max,&
            tmp_search_trace_n_points = tmp_search_trace_n_points,&
            tmp_search_trace_global_js_divergence = tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower = tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper = tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width = tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative = tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap = tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta = tmp_search_trace_delta,&
            tmp_search_trace_delta_median = tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max = tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins = tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed = tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals = tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy = tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy = tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy = tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins = tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins = tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low = tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high = tmp_search_trace_shared_residual_range_high,&
            tmp_search_candidate_k_start = tmp_search_candidate_k_start,&
            tmp_search_candidate_n_points = tmp_search_candidate_n_points,&
            tmp_search_candidate_status = tmp_search_candidate_status,&
            tmp_candidate_confidence_interval = tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence = tmp_prev_global_js_divergence,&
            tmp_delta_perm = tmp_delta_perm,&
            tmp_best_uncertainty_confidence_interval = tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point = tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point = tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low = tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high = tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low = tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high = tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search = rerun_parameter_search,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            bin_count_mode = bin_count_mode,&
            delta_abs_tolerance = delta_abs_tolerance,&
            bootstrap_mode = bootstrap_mode,&
            divergence_mode = divergence_mode,&
            divergence_weight_mode = divergence_weight_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_adaptive_on_sample_assignment_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_label_permutation_test_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test_permutation(module):run_js_comp_test_label_permutation_test_expert]] to prepare it yourself.
    !| Tests whether the assignment of samples to studies influences each study's JS-Comp-Test value
    !| (null hypothesis: it does not). The same data, parameters and options as
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]],
    !| whose input layout (all studies in `expression`, NaN for a missing gene, all-NaN replicate
    !| slots are padding) and assignment convention apply here too.
    !|
    !| 1. The observed data (the identity assignment) give `global_js_divergence_observed`.
    !| 2. `n_permutations` times, the study labels of the samples are shuffled without replacement (a
    !| uniformly random permutation, unbiased Fisher-Yates, from one random number generator seeded
    !| with `random_seed`); every study keeps its sample count. The whole JS-Comp-Test runs on the
    !| shuffled data; row `i` of `jsd_null` holds its values.
    !| 3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`,
    !| the Laplace correction of the residual-level test, so it is never 0 and at least
    !| `1/(n_permutations + 1)`.
    !|
    !| **Rerunning the parameter search** (`rerun_parameter_search = .true.`): every data set, the
    !| observed one included, first runs the fixed-k parameter search and then the run with the
    !| `(n_points, n_neighbors)` it selects; `n_points` and `n_neighbors` then only size the work
    !| arrays. **Run time:** each permutation runs the whole JS-Comp-Test; with the rerun option it
    !| runs the whole parameter search first. The routine is sequential: a caller who wants
    !| permutations in parallel calls
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]]
    !| itself with its own assignments (e.g. from R), one call per assignment.
    !|
    !| **Failures.** A permutation whose JS-Comp-Test fails with invalid input (for example a study
    !| left without residuals at some reference point) gets a NaN row in `jsd_null`, is counted in
    !| `n_failed_permutations` and does not count in `n_valid`. No valid permutation gives NaN
    !| p-values (a direct call with `n_permutations = 0` too; the published entry points reject 0,
    !| since it is the extent of `jsd_null`). Any other error ends the test and is returned, with NaN
    !| p-values; the rows of the permutations run so far are kept. An error on the observed data is
    !| returned at once (`global_js_divergence_observed`, `p_values` and `jsd_null` NaN), as is a study
    !| without any sample (the call fails with invalid input). In the Fortran expert tier, where the
    !| caller supplies the work-array capacities, a capacity below the search's own bound fails the
    !| same way with the rerun option.
    !|
    !| **Known limitation.** A gene missing in a whole study travels with the samples, so after
    !| shuffling that study may have values for it.
    subroutine run_js_comp_test_label_permutation_test(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            expression,&
            n_points,&
            n_neighbors,&
            n_permutations,&
            global_js_divergence_observed,&
            p_values,&
            jsd_null,&
            n_failed_permutations,&
            rerun_parameter_search,&
            n_bootstraps,&
            join_method,&
            min_residuals_per_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            bin_count_mode,&
            neighborhood_size_mode,&
            delta_abs_tolerance,&
            bootstrap_mode,&
            divergence_mode,&
            divergence_weight_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_permutations
            !! Number of random assignments (rows of `jsd_null`)
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence_observed
            !! Weighted global JSD of each study on the observed data (the identity assignment)
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
            !! whose value reached the observed one; NaN when none is valid or on error
        real(real64), dimension(n_permutations, n_studies), intent(out) :: jsd_null
            !! Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation
        integer(int32), intent(out) :: n_failed_permutations
            !! Number of permutations whose JS-Comp-Test failed with invalid input
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! `.true.`: run the fixed-k parameter search on the relabelled data first and use its
            !! `(n_points, n_neighbors)`; `.false.`: use `n_points` and `n_neighbors`
            !! The default value is `.false.`.
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `1_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed of the random number generator of the label draws, and of the bootstrap of
            !! every rerun parameter search (a separate generator)
            !! The default value is `42_int32`.
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
        integer(int32), intent(in), optional :: neighborhood_size_mode
            !! Rerun search: which neighborhood sizes the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are left out
            !! by default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode                                                       | Value                                                                                                    |
            !! |------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Large neighborhoods only, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods only, `KX_FACTORS(2)` = 0.5            | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! | Both, large then small for each `n_points`                 | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_BOTH(variable)]]  |
            !! The default value is `0_int32`.
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0_real64`.
        integer(int32), intent(in), optional :: bootstrap_mode
            !! Rerun search (used only with `rerun_parameter_search`): what its bootstrap draws each
            !! study's histograms from. By default each study's own counts: a robustness interval of
            !! the observed JSD. The pooled counts of all studies give a noise band of "no difference
            !! between the studies" instead, kept for comparison
            !!
            !! | Mode                                                                                             | Value                                                                          |
            !! |--------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
            !! | Each study's own counts, a robustness interval (default)                                         | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_PER_STUDY(variable)]] |
            !! | The pooled consensus counts of all studies, a noise band (the bootstrap before the mode existed) | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_POOLED(variable)]]    |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_mode
            !! Experiment option: how each study's per-point divergence is computed. By default the
            !! JSD of the study against the consensus pmf, the formula before the option existed
            !!
            !! | Mode                                                                                                                | Value                                                                           |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | JSD of each study against the consensus pmf (default)                                                               | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | For exactly two studies the pairwise JSD of the two, the same value for both (otherwise as the default)             | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! | One-sided Kullback-Leibler divergence of each study against the consensus pmf, at most log2(number of studies) bits | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_ONE_SIDED(variable)]]   |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_weight_mode
            !! Experiment option: how each reference point's weight in a study's global divergence is
            !! computed. By default the study's residual count there plus the count of all studies
            !! pooled, the weights before the option existed. For two studies in the pairwise
            !! `divergence_mode` the weights are always the two studies' counts
            !!
            !! | Mode                                                                   | Value                                                                                      |
            !! |------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
            !! | The study's count plus the pooled count of all studies (default)       | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL(variable)]] |
            !! | The pooled count of all studies only, the same weights for every study | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_POOL(variable)]]            |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)
        integer(int32) :: max_n_points_candidate
        integer(int32) :: max_n_neighbors_candidate
        integer(int32) :: n_bootstrapping_top_k_jsds
        integer(int32), dimension(:), allocatable :: tmp_sample_slots
        integer(int32), dimension(:), allocatable :: tmp_sample_assignment
        real(real64), dimension(:), allocatable :: tmp_permuted_global_js_divergence
        logical(c_bool), dimension(:), allocatable :: tmp_sample_mark
        real(real64), dimension(:, :, :), allocatable :: tmp_permuted_expression
        real(real64), dimension(:, :), allocatable :: tmp_gene_means
        real(real64), dimension(:, :, :), allocatable :: tmp_residuals
        integer(int32), dimension(:, :), allocatable :: tmp_gene_means_perms
        integer(int32), dimension(:), allocatable :: tmp_gene_means_perm_all
        real(real64), dimension(:), allocatable :: tmp_x_star
        integer(int32), dimension(:, :, :), allocatable :: tmp_neighborhood_indices
        integer(int32), dimension(:, :, :), allocatable :: tmp_neighborhood_range
        real(real64), dimension(:, :, :), allocatable :: tmp_neighborhood_residuals_gathered
        integer(int32), dimension(:), allocatable :: tmp_point_pool
        integer(int32), dimension(:), allocatable :: tmp_point_draw
        integer(int32), dimension(:), allocatable :: tmp_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_high
        logical(c_bool), dimension(:), allocatable :: tmp_occupancy_failed
        integer(int32), dimension(:), allocatable :: tmp_n_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_min_bin_occupancy
        real(real64), dimension(:), allocatable :: tmp_mean_bin_occupancy
        integer(int32), dimension(:), allocatable :: tmp_max_bin_occupancy
        integer(int32), dimension(:), allocatable :: tmp_sturges_bins
        integer(int32), dimension(:), allocatable :: tmp_fd_bins
        real(real64), dimension(:, :, :), allocatable :: tmp_pmfs
        integer(int32), dimension(:, :, :), allocatable :: tmp_counts
        integer(int32), dimension(:, :), allocatable :: tmp_included_n_reps
        real(real64), dimension(:, :), allocatable :: tmp_mean_pmf
        integer(int32), dimension(:, :), allocatable :: tmp_mean_pmf_counts
        integer(int32), dimension(:), allocatable :: tmp_mean_pmf_included_n_reps
        real(real64), dimension(:, :), allocatable :: tmp_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_weights
        real(real64), dimension(:), allocatable :: tmp_global_js_divergence
        real(real64), dimension(:), allocatable :: tmp_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_residual_to_study_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts_search
        integer(int32), dimension(:), allocatable :: tmp_study_block_ends
        integer(int32), dimension(:), allocatable :: tmp_study_block_sizes
        integer(int32), dimension(:), allocatable :: tmp_point_n_neighbors
        real(real64), dimension(:), allocatable :: tmp_p_values
        real(real64), dimension(:), allocatable :: tmp_p_values_observed_consensus
        real(real64), dimension(:, :), allocatable :: tmp_permutation_mean_pmf
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_mean_pmf_counts
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_counts
        real(real64), dimension(:, :, :), allocatable :: tmp_permutation_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_permutation_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_permutation_weights
        integer(int32), dimension(:), allocatable :: tmp_search_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_search_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_search_shared_residual_range_high
        real(real64), dimension(:, :), allocatable :: tmp_search_confidence_interval
        integer(int32), dimension(:), allocatable :: tmp_search_trace_n_points
        integer(int32), dimension(:), allocatable :: tmp_search_trace_n_neighbors
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_global_js_divergence
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_lower
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_upper
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_width
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_width_relative
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_overlap
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_delta
        real(real64), dimension(:), allocatable :: tmp_search_trace_delta_median
        real(real64), dimension(:), allocatable :: tmp_search_trace_delta_max
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_selected_n_bins
        logical(c_bool), dimension(:, :), allocatable :: tmp_search_trace_occupancy_failed
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_n_pooled_residuals
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_min_bin_occupancy
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_mean_bin_occupancy
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_max_bin_occupancy
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_sturges_bins
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_fd_bins
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_shared_residual_range_low
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_shared_residual_range_high
        real(real64), dimension(:, :), allocatable :: tmp_candidate_confidence_interval
        real(real64), dimension(:, :, :), allocatable :: tmp_bootstrapping_top_k_jsds
        real(real64), dimension(:), allocatable :: tmp_prev_global_js_divergence
        real(real64), dimension(:, :), allocatable :: tmp_best_uncertainty_confidence_interval
        integer(int32), dimension(:), allocatable :: tmp_best_n_bins_per_point
        integer(int32), dimension(:), allocatable :: tmp_best_uncertainty_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_best_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_best_shared_residual_range_high
        real(real64), dimension(:), allocatable :: tmp_best_uncertainty_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_best_uncertainty_shared_residual_range_high

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=6_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_in_range_int(n_permutations, ierr, arg_pos=7_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=13_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=15_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=16_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=17_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=19_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=20_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=21_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=22_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=23_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=24_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=25_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=26_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=27_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=28_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(delta_abs_tolerance, ierr, arg_pos=32_int32, min=0.0_real64)
        call validate_all_in_range_real(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        if (present(join_method)) then; if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=14_int32); end if
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=18_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=30_int32); end if
        if (present(neighborhood_size_mode)) then; if (neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_LARGE .and. neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_SMALL .and. neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=31_int32); end if
        if (present(bootstrap_mode)) then; if (bootstrap_mode /= MODE_BOOTSTRAP_PER_STUDY .and. bootstrap_mode /= MODE_BOOTSTRAP_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=33_int32); end if
        if (present(divergence_mode)) then; if (divergence_mode /= MODE_DIVERGENCE_SYMMETRIC .and. divergence_mode /= MODE_DIVERGENCE_PAIRWISE_K2 .and. divergence_mode /= MODE_DIVERGENCE_ONE_SIDED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=34_int32); end if
        if (present(divergence_weight_mode)) then; if (divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL .and. divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_POOL) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=35_int32); end if
        if (is_err(ierr)) return
#endif

        call calc_js_comp_test_permutation_bounds(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            ierr = ierr&
        )
        if (is_err(ierr)) return
        M_ALLOCATE(tmp_sample_slots(max_n_reps_all_studies*n_studies))
        M_ALLOCATE(tmp_sample_assignment(max_n_reps_all_studies*n_studies))
        M_ALLOCATE(tmp_permuted_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_sample_mark(max_n_reps_all_studies*n_studies))
        M_ALLOCATE(tmp_permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_gene_means(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_gene_means_perms(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_x_star(max_n_points_candidate))
        M_ALLOCATE(tmp_neighborhood_indices(max_n_neighbors_candidate, max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_neighborhood_range(2, max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate))
        M_ALLOCATE(tmp_point_pool(256))
        M_ALLOCATE(tmp_point_draw(256))
        M_ALLOCATE(tmp_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_occupancy_failed(max_n_points_candidate))
        M_ALLOCATE(tmp_n_pooled_residuals(max_n_points_candidate))
        M_ALLOCATE(tmp_min_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_mean_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_max_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_sturges_bins(max_n_points_candidate))
        M_ALLOCATE(tmp_fd_bins(max_n_points_candidate))
        M_ALLOCATE(tmp_pmfs(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_counts(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_included_n_reps(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_mean_pmf(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_mean_pmf_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        M_ALLOCATE(tmp_js_divergences(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_weights(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_ALLOCATE(tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
        M_ALLOCATE(tmp_study_block_ends(n_studies))
        M_ALLOCATE(tmp_study_block_sizes(n_studies))
        M_ALLOCATE(tmp_point_n_neighbors(n_studies))
        M_ALLOCATE(tmp_p_values(n_studies))
        M_ALLOCATE(tmp_p_values_observed_consensus(n_studies))
        M_ALLOCATE(tmp_permutation_mean_pmf(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_mean_pmf_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_pmfs(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_permutation_js_divergences(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_permutation_weights(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_search_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_search_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_search_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_search_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_search_trace_n_points(16))
        M_ALLOCATE(tmp_search_trace_n_neighbors(16))
        M_ALLOCATE(tmp_search_trace_global_js_divergence(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_lower(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_upper(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_width(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_width_relative(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_overlap(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_delta(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_delta_median(16))
        M_ALLOCATE(tmp_search_trace_delta_max(16))
        M_ALLOCATE(tmp_search_trace_selected_n_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_occupancy_failed(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_n_pooled_residuals(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_min_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_mean_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_max_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_sturges_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_fd_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_shared_residual_range_low(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_shared_residual_range_high(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_candidate_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        M_ALLOCATE(tmp_prev_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_best_uncertainty_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_best_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_best_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_best_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))

        call run_js_comp_test_label_permutation_test_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            expression = expression,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            n_permutations = n_permutations,&
            global_js_divergence_observed = global_js_divergence_observed,&
            p_values = p_values,&
            jsd_null = jsd_null,&
            n_failed_permutations = n_failed_permutations,&
            tmp_sample_slots = tmp_sample_slots,&
            tmp_sample_assignment = tmp_sample_assignment,&
            tmp_permuted_global_js_divergence = tmp_permuted_global_js_divergence,&
            tmp_sample_mark = tmp_sample_mark,&
            tmp_permuted_expression = tmp_permuted_expression,&
            tmp_gene_means = tmp_gene_means,&
            tmp_residuals = tmp_residuals,&
            tmp_gene_means_perms = tmp_gene_means_perms,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_x_star = tmp_x_star,&
            tmp_neighborhood_indices = tmp_neighborhood_indices,&
            tmp_neighborhood_range = tmp_neighborhood_range,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            tmp_n_bins_per_point = tmp_n_bins_per_point,&
            tmp_shared_residual_range_low = tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high = tmp_shared_residual_range_high,&
            tmp_occupancy_failed = tmp_occupancy_failed,&
            tmp_n_pooled_residuals = tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy = tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy = tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy = tmp_max_bin_occupancy,&
            tmp_sturges_bins = tmp_sturges_bins,&
            tmp_fd_bins = tmp_fd_bins,&
            tmp_pmfs = tmp_pmfs,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps = tmp_included_n_reps,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps = tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_point_n_neighbors = tmp_point_n_neighbors,&
            tmp_p_values = tmp_p_values,&
            tmp_p_values_observed_consensus = tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_search_n_bins_per_point = tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low = tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high = tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval = tmp_search_confidence_interval,&
            tmp_search_trace_n_points = tmp_search_trace_n_points,&
            tmp_search_trace_n_neighbors = tmp_search_trace_n_neighbors,&
            tmp_search_trace_global_js_divergence = tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower = tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper = tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width = tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative = tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap = tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta = tmp_search_trace_delta,&
            tmp_search_trace_delta_median = tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max = tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins = tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed = tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals = tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy = tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy = tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy = tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins = tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins = tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low = tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high = tmp_search_trace_shared_residual_range_high,&
            tmp_candidate_confidence_interval = tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence = tmp_prev_global_js_divergence,&
            tmp_best_uncertainty_confidence_interval = tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point = tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point = tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low = tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high = tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low = tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high = tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search = rerun_parameter_search,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            bin_count_mode = bin_count_mode,&
            neighborhood_size_mode = neighborhood_size_mode,&
            delta_abs_tolerance = delta_abs_tolerance,&
            bootstrap_mode = bootstrap_mode,&
            divergence_mode = divergence_mode,&
            divergence_weight_mode = divergence_weight_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_label_permutation_test

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_label_permutation_test_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test_permutation(module):run_js_comp_test_label_permutation_test]] does both.
    !| Tests whether the assignment of samples to studies influences each study's JS-Comp-Test value
    !| (null hypothesis: it does not). The same data, parameters and options as
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]],
    !| whose input layout (all studies in `expression`, NaN for a missing gene, all-NaN replicate
    !| slots are padding) and assignment convention apply here too.
    !|
    !| 1. The observed data (the identity assignment) give `global_js_divergence_observed`.
    !| 2. `n_permutations` times, the study labels of the samples are shuffled without replacement (a
    !| uniformly random permutation, unbiased Fisher-Yates, from one random number generator seeded
    !| with `random_seed`); every study keeps its sample count. The whole JS-Comp-Test runs on the
    !| shuffled data; row `i` of `jsd_null` holds its values.
    !| 3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`,
    !| the Laplace correction of the residual-level test, so it is never 0 and at least
    !| `1/(n_permutations + 1)`.
    !|
    !| **Rerunning the parameter search** (`rerun_parameter_search = .true.`): every data set, the
    !| observed one included, first runs the fixed-k parameter search and then the run with the
    !| `(n_points, n_neighbors)` it selects; `n_points` and `n_neighbors` then only size the work
    !| arrays. **Run time:** each permutation runs the whole JS-Comp-Test; with the rerun option it
    !| runs the whole parameter search first. The routine is sequential: a caller who wants
    !| permutations in parallel calls
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]]
    !| itself with its own assignments (e.g. from R), one call per assignment.
    !|
    !| **Failures.** A permutation whose JS-Comp-Test fails with invalid input (for example a study
    !| left without residuals at some reference point) gets a NaN row in `jsd_null`, is counted in
    !| `n_failed_permutations` and does not count in `n_valid`. No valid permutation gives NaN
    !| p-values (a direct call with `n_permutations = 0` too; the published entry points reject 0,
    !| since it is the extent of `jsd_null`). Any other error ends the test and is returned, with NaN
    !| p-values; the rows of the permutations run so far are kept. An error on the observed data is
    !| returned at once (`global_js_divergence_observed`, `p_values` and `jsd_null` NaN), as is a study
    !| without any sample (the call fails with invalid input). In the Fortran expert tier, where the
    !| caller supplies the work-array capacities, a capacity below the search's own bound fails the
    !| same way with the rerun option.
    !|
    !| **Known limitation.** A gene missing in a whole study travels with the samples, so after
    !| shuffling that study may have values for it.
    subroutine run_js_comp_test_label_permutation_test_expert(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            expression,&
            n_points,&
            n_neighbors,&
            max_n_points_candidate,&
            max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds,&
            n_permutations,&
            global_js_divergence_observed,&
            p_values,&
            jsd_null,&
            n_failed_permutations,&
            tmp_sample_slots,&
            tmp_sample_assignment,&
            tmp_permuted_global_js_divergence,&
            tmp_sample_mark,&
            tmp_permuted_expression,&
            tmp_gene_means,&
            tmp_residuals,&
            tmp_gene_means_perms,&
            tmp_gene_means_perm_all,&
            tmp_x_star,&
            tmp_neighborhood_indices,&
            tmp_neighborhood_range,&
            tmp_neighborhood_residuals_gathered,&
            tmp_point_pool,&
            tmp_point_draw,&
            tmp_n_bins_per_point,&
            tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high,&
            tmp_occupancy_failed,&
            tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy,&
            tmp_sturges_bins,&
            tmp_fd_bins,&
            tmp_pmfs,&
            tmp_counts,&
            tmp_included_n_reps,&
            tmp_mean_pmf,&
            tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences,&
            tmp_weights,&
            tmp_global_js_divergence,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm,&
            tmp_bin_counts_search,&
            tmp_study_block_ends,&
            tmp_study_block_sizes,&
            tmp_point_n_neighbors,&
            tmp_p_values,&
            tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts,&
            tmp_permutation_pmfs,&
            tmp_permutation_js_divergences,&
            tmp_permutation_weights,&
            tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval,&
            tmp_search_trace_n_points,&
            tmp_search_trace_n_neighbors,&
            tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta,&
            tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high,&
            tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence,&
            tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search,&
            n_bootstraps,&
            join_method,&
            min_residuals_per_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            bin_count_mode,&
            neighborhood_size_mode,&
            delta_abs_tolerance,&
            bootstrap_mode,&
            divergence_mode,&
            divergence_weight_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_permutation_impl(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `n_points`.
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Capacity of the per-neighbor work arrays
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_permutation_impl(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `n_neighbors`.
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! It is *VERY IMPORTANT* to compute this argument from the `n_bootstrapping_top_k_jsds` output produced by [[tox_data_integration_js_comp_test_permutation_impl(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_permutations
            !! Number of random assignments (rows of `jsd_null`)
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence_observed
            !! Weighted global JSD of each study on the observed data (the identity assignment)
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
            !! whose value reached the observed one; NaN when none is valid or on error
        real(real64), dimension(n_permutations, n_studies), intent(out) :: jsd_null
            !! Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation
        integer(int32), intent(out) :: n_failed_permutations
            !! Number of permutations whose JS-Comp-Test failed with invalid input
        integer(int32), dimension(max_n_reps_all_studies*n_studies), intent(out) :: tmp_sample_slots
            !! Working array: the flat slot index of every sample
        integer(int32), dimension(max_n_reps_all_studies*n_studies), intent(out) :: tmp_sample_assignment
            !! Working array: the current assignment, in its first `n_samples` elements
        real(real64), dimension(n_studies), intent(out) :: tmp_permuted_global_js_divergence
            !! Working array: one permutation's global JSD per study
        logical(c_bool), dimension(max_n_reps_all_studies*n_studies), intent(out) :: tmp_sample_mark
            !! Working array: the values of the current assignment already seen
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_permuted_expression
            !! Working array: the relabelled data set
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means
            !! Working array: the relabelled data set's gene means per study
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_residuals
            !! Working array: the relabelled data set's residuals
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means_perms
            !! Working array of the search and the run: each study's sorting permutation of its gene means
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array of the search and the run: sorting permutation of the pooled gene means
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Working array of the search and the run: the reference points
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: tmp_neighborhood_indices
            !! Working array of the search and the run: every study's neighbor gene indices
        integer(int32), dimension(2, max_n_points_candidate, n_studies), intent(out) :: tmp_neighborhood_range
            !! Working array of the search and the run: neighborhood spans
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array of the search and the run: one study's gathered neighborhood residuals
        integer(int32), dimension(256), intent(out) :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(int32), dimension(256), intent(out) :: tmp_point_draw
            !! Working array of the search and the run: one study's draw at one point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array of the search and the run: bin count per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array of the search and the run: lower residual-range bound per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array of the search and the run: upper residual-range bound per reference point
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
            !! Working array of the search and the run: bin-count search failure flag per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
            !! Working array of the search and the run: pooled residual count per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
            !! Working array of the search and the run: minimum bin occupancy per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
            !! Working array of the search and the run: mean bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
            !! Working array of the search and the run: maximum bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
            !! Working array of the search and the run: Sturges' rule diagnostic per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_fd_bins
            !! Working array of the search and the run: Freedman-Diaconis diagnostic per reference point
        real(real64), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(int32), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(real64), dimension(max_n_points_candidate, 256), intent(out) :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        integer(int32), dimension(n_studies), intent(out) :: tmp_point_n_neighbors
            !! Working array of the search and the run: per-study neighbor counts, then sort scratch
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values
            !! Working array: the run's p-values (NaN, its permutation test does not run)
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values_observed_consensus
            !! Working array: the run's comparison p-values (NaN, its permutation test does not run)
        real(real64), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_mean_pmf
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_counts
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_search_n_bins_per_point
            !! Working array: the search's chosen bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_low
            !! Working array: the search's chosen lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_high
            !! Working array: the search's chosen upper residual-range bounds
        real(real64), dimension(2, n_studies), intent(out) :: tmp_search_confidence_interval
            !! Working array: the search's chosen confidence interval
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_n_neighbors
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(real64), dimension(16), intent(out) :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(real64), dimension(16), intent(out) :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_shared_residual_range_high
            !! Working array: the search's trace
        real(real64), dimension(2, n_studies), intent(out) :: tmp_candidate_confidence_interval
            !! Working array of the search: the current candidate's confidence interval
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array of the search: the bootstrap's bottom-k/top-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array of the search: the previous admissible candidate's global JSD
        real(real64), dimension(2, n_studies), intent(out) :: tmp_best_uncertainty_confidence_interval
            !! Working array of the search: the smallest-uncertainty candidate's confidence interval
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_n_bins_per_point
            !! Working array of the search: the best candidate's bin counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_n_bins_per_point
            !! Working array of the search: the smallest-uncertainty candidate's bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_low
            !! Working array of the search: the best candidate's lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_high
            !! Working array of the search: the best candidate's upper residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array of the search: the smallest-uncertainty candidate's lower bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array of the search: the smallest-uncertainty candidate's upper bounds
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! `.true.`: run the fixed-k parameter search on the relabelled data first and use its
            !! `(n_points, n_neighbors)`; `.false.`: use `n_points` and `n_neighbors`
            !! The default value is `.false.`.
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `1_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed of the random number generator of the label draws, and of the bootstrap of
            !! every rerun parameter search (a separate generator)
            !! The default value is `42_int32`.
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
        integer(int32), intent(in), optional :: neighborhood_size_mode
            !! Rerun search: which neighborhood sizes the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are left out
            !! by default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode                                                       | Value                                                                                                    |
            !! |------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Large neighborhoods only, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods only, `KX_FACTORS(2)` = 0.5            | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! | Both, large then small for each `n_points`                 | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_BOTH(variable)]]  |
            !! The default value is `0_int32`.
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0_real64`.
        integer(int32), intent(in), optional :: bootstrap_mode
            !! Rerun search (used only with `rerun_parameter_search`): what its bootstrap draws each
            !! study's histograms from. By default each study's own counts: a robustness interval of
            !! the observed JSD. The pooled counts of all studies give a noise band of "no difference
            !! between the studies" instead, kept for comparison
            !!
            !! | Mode                                                                                             | Value                                                                          |
            !! |--------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
            !! | Each study's own counts, a robustness interval (default)                                         | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_PER_STUDY(variable)]] |
            !! | The pooled consensus counts of all studies, a noise band (the bootstrap before the mode existed) | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_POOLED(variable)]]    |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_mode
            !! Experiment option: how each study's per-point divergence is computed. By default the
            !! JSD of the study against the consensus pmf, the formula before the option existed
            !!
            !! | Mode                                                                                                                | Value                                                                           |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | JSD of each study against the consensus pmf (default)                                                               | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | For exactly two studies the pairwise JSD of the two, the same value for both (otherwise as the default)             | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! | One-sided Kullback-Leibler divergence of each study against the consensus pmf, at most log2(number of studies) bits | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_ONE_SIDED(variable)]]   |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_weight_mode
            !! Experiment option: how each reference point's weight in a study's global divergence is
            !! computed. By default the study's residual count there plus the count of all studies
            !! pooled, the weights before the option existed. For two studies in the pairwise
            !! `divergence_mode` the weights are always the two studies' counts
            !!
            !! | Mode                                                                   | Value                                                                                      |
            !! |------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
            !! | The study's count plus the pooled count of all studies (default)       | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL(variable)]] |
            !! | The pooled count of all studies only, the same weights for every study | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_POOL(variable)]]            |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=6_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_in_range_int(max_n_points_candidate, ierr, arg_pos=7_int32, min=n_points)
        call validate_in_range_int(max_n_neighbors_candidate, ierr, arg_pos=8_int32, min=n_neighbors)
        call validate_in_range_int(n_bootstrapping_top_k_jsds, ierr, arg_pos=9_int32, min=1_int32)
        call validate_in_range_int(n_permutations, ierr, arg_pos=10_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=100_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=102_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=103_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=104_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=106_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=107_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=108_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=109_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=110_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=111_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=112_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=113_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=114_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=115_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(delta_abs_tolerance, ierr, arg_pos=119_int32, min=0.0_real64)
        call validate_all_in_range_real(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        if (present(join_method)) then; if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=101_int32); end if
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=105_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=117_int32); end if
        if (present(neighborhood_size_mode)) then; if (neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_LARGE .and. neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_SMALL .and. neighborhood_size_mode /= MODE_NEIGHBORHOOD_SIZE_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=118_int32); end if
        if (present(bootstrap_mode)) then; if (bootstrap_mode /= MODE_BOOTSTRAP_PER_STUDY .and. bootstrap_mode /= MODE_BOOTSTRAP_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=120_int32); end if
        if (present(divergence_mode)) then; if (divergence_mode /= MODE_DIVERGENCE_SYMMETRIC .and. divergence_mode /= MODE_DIVERGENCE_PAIRWISE_K2 .and. divergence_mode /= MODE_DIVERGENCE_ONE_SIDED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=121_int32); end if
        if (present(divergence_weight_mode)) then; if (divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL .and. divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_POOL) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=122_int32); end if
        if (is_err(ierr)) return
#endif

        call run_js_comp_test_label_permutation_test_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            expression = expression,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            n_permutations = n_permutations,&
            global_js_divergence_observed = global_js_divergence_observed,&
            p_values = p_values,&
            jsd_null = jsd_null,&
            n_failed_permutations = n_failed_permutations,&
            tmp_sample_slots = tmp_sample_slots,&
            tmp_sample_assignment = tmp_sample_assignment,&
            tmp_permuted_global_js_divergence = tmp_permuted_global_js_divergence,&
            tmp_sample_mark = tmp_sample_mark,&
            tmp_permuted_expression = tmp_permuted_expression,&
            tmp_gene_means = tmp_gene_means,&
            tmp_residuals = tmp_residuals,&
            tmp_gene_means_perms = tmp_gene_means_perms,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_x_star = tmp_x_star,&
            tmp_neighborhood_indices = tmp_neighborhood_indices,&
            tmp_neighborhood_range = tmp_neighborhood_range,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            tmp_n_bins_per_point = tmp_n_bins_per_point,&
            tmp_shared_residual_range_low = tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high = tmp_shared_residual_range_high,&
            tmp_occupancy_failed = tmp_occupancy_failed,&
            tmp_n_pooled_residuals = tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy = tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy = tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy = tmp_max_bin_occupancy,&
            tmp_sturges_bins = tmp_sturges_bins,&
            tmp_fd_bins = tmp_fd_bins,&
            tmp_pmfs = tmp_pmfs,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps = tmp_included_n_reps,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps = tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_point_n_neighbors = tmp_point_n_neighbors,&
            tmp_p_values = tmp_p_values,&
            tmp_p_values_observed_consensus = tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_search_n_bins_per_point = tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low = tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high = tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval = tmp_search_confidence_interval,&
            tmp_search_trace_n_points = tmp_search_trace_n_points,&
            tmp_search_trace_n_neighbors = tmp_search_trace_n_neighbors,&
            tmp_search_trace_global_js_divergence = tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower = tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper = tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width = tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative = tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap = tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta = tmp_search_trace_delta,&
            tmp_search_trace_delta_median = tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max = tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins = tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed = tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals = tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy = tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy = tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy = tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins = tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins = tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low = tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high = tmp_search_trace_shared_residual_range_high,&
            tmp_candidate_confidence_interval = tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence = tmp_prev_global_js_divergence,&
            tmp_best_uncertainty_confidence_interval = tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point = tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point = tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low = tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high = tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low = tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high = tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search = rerun_parameter_search,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            bin_count_mode = bin_count_mode,&
            neighborhood_size_mode = neighborhood_size_mode,&
            delta_abs_tolerance = delta_abs_tolerance,&
            bootstrap_mode = bootstrap_mode,&
            divergence_mode = divergence_mode,&
            divergence_weight_mode = divergence_weight_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_label_permutation_test_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_label_permutation_test_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test_permutation(module):run_js_comp_test_adaptive_label_permutation_test_expert]] to prepare it yourself.
    !| The adaptive counterpart of
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_label_permutation_test_impl(interface)]]:
    !| the same test (observed value, `n_permutations` shuffles of the study labels without
    !| replacement, p-value `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`), each
    !| data set run by
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]
    !| with the same data, knobs and options; its input layout and assignment convention apply here too.
    !|
    !| **Rerunning the parameter search** (`rerun_parameter_search = .true.`): every data set first
    !| runs the adaptive parameter search and then the construction and run with the knobs it
    !| selects. **Run time:** each permutation runs the whole JS-Comp-Test; with the rerun option it
    !| runs the whole parameter search first. The routine is sequential: a caller who wants
    !| permutations in parallel calls
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]
    !| itself with its own assignments (e.g. from R), one call per assignment.
    !|
    !| **Failures.** A permutation that fails with invalid input (a study left without residuals at
    !| some reference point, or a construction that builds no neighborhood or more neighborhoods
    !| than the work arrays hold) gets a NaN row in `jsd_null`, is counted in `n_failed_permutations`
    !| and does not count in `n_valid`. No valid permutation gives NaN p-values (a direct call with
    !| `n_permutations = 0` too; the published entry points reject 0, since it is the extent of
    !| `jsd_null`). Any other error ends the test and is returned, with NaN p-values; the rows of the
    !| permutations run so far are kept. An error on the observed data is returned at once
    !| (`global_js_divergence_observed`, `p_values` and `jsd_null` NaN), as is a study without any
    !| sample (the call fails with invalid input). In the Fortran expert tier, where the caller
    !| supplies the work-array capacity, a capacity below the search's own fails the same way with
    !| the rerun option.
    !|
    !| **Known limitation.** A gene missing in a whole study travels with the samples, so after
    !| shuffling that study may have values for it.
    subroutine run_js_comp_test_adaptive_label_permutation_test(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            expression,&
            k_start,&
            k_step,&
            k_max,&
            n_permutations,&
            global_js_divergence_observed,&
            p_values,&
            jsd_null,&
            n_failed_permutations,&
            rerun_parameter_search,&
            n_bootstraps,&
            join_method,&
            min_residuals_per_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            tau,&
            mad_distance_factor,&
            max_pooled_residuals,&
            min_study_neighbors,&
            bin_count_mode,&
            delta_abs_tolerance,&
            bootstrap_mode,&
            divergence_mode,&
            divergence_weight_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_permutations
            !! Number of random assignments (rows of `jsd_null`)
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! The minimum valid value is `k_start`.
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence_observed
            !! Weighted global JSD of each study on the observed data (the identity assignment)
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
            !! whose value reached the observed one; NaN when none is valid or on error
        real(real64), dimension(n_permutations, n_studies), intent(out) :: jsd_null
            !! Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation
        integer(int32), intent(out) :: n_failed_permutations
            !! Number of permutations whose JS-Comp-Test failed with invalid input
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! `.true.`: run the adaptive parameter search on the relabelled data first and use its
            !! knobs; `.false.`: use `k_start`, `k_step` and `k_max`
            !! The default value is `.false.`.
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `1_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed of the random number generator of the label draws, and of the bootstrap of
            !! every rerun parameter search (a separate generator)
            !! The default value is `42_int32`.
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have (construction status; in
            !! the rerun search, candidate admissibility)
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
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
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0_real64`.
        integer(int32), intent(in), optional :: bootstrap_mode
            !! Rerun search (used only with `rerun_parameter_search`): what its bootstrap draws each
            !! study's histograms from. By default each study's own counts: a robustness interval of
            !! the observed JSD. The pooled counts of all studies give a noise band of "no difference
            !! between the studies" instead, kept for comparison
            !!
            !! | Mode                                                                                             | Value                                                                          |
            !! |--------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
            !! | Each study's own counts, a robustness interval (default)                                         | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_PER_STUDY(variable)]] |
            !! | The pooled consensus counts of all studies, a noise band (the bootstrap before the mode existed) | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_POOLED(variable)]]    |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_mode
            !! Experiment option: how each study's per-point divergence is computed. By default the
            !! JSD of the study against the consensus pmf, the formula before the option existed
            !!
            !! | Mode                                                                                                                | Value                                                                           |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | JSD of each study against the consensus pmf (default)                                                               | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | For exactly two studies the pairwise JSD of the two, the same value for both (otherwise as the default)             | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! | One-sided Kullback-Leibler divergence of each study against the consensus pmf, at most log2(number of studies) bits | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_ONE_SIDED(variable)]]   |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_weight_mode
            !! Experiment option: how each reference point's weight in a study's global divergence is
            !! computed. By default the study's residual count there plus the count of all studies
            !! pooled, the weights before the option existed. For two studies in the pairwise
            !! `divergence_mode` the weights are always the two studies' counts
            !!
            !! | Mode                                                                   | Value                                                                                      |
            !! |------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
            !! | The study's count plus the pooled count of all studies (default)       | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL(variable)]] |
            !! | The pooled count of all studies only, the same weights for every study | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_POOL(variable)]]            |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)
        integer(int32) :: max_n_points_candidate
        integer(int32) :: n_bootstrapping_top_k_jsds
        integer(int32), dimension(:), allocatable :: tmp_sample_slots
        integer(int32), dimension(:), allocatable :: tmp_sample_assignment
        real(real64), dimension(:), allocatable :: tmp_permuted_global_js_divergence
        logical(c_bool), dimension(:), allocatable :: tmp_sample_mark
        real(real64), dimension(:, :, :), allocatable :: tmp_permuted_expression
        real(real64), dimension(:, :), allocatable :: tmp_gene_means
        real(real64), dimension(:, :, :), allocatable :: tmp_residuals
        integer(int32), dimension(:), allocatable :: tmp_gene_means_perm_all
        real(real64), dimension(:), allocatable :: tmp_x_star
        integer(int32), dimension(:, :), allocatable :: tmp_pooled_neighborhood_range
        integer(int32), dimension(:, :), allocatable :: tmp_n_neighbors_per_point
        integer(int32), dimension(:), allocatable :: tmp_stop_reason
        real(real64), dimension(:), allocatable :: tmp_neighborhood_dispersion
        real(real64), dimension(:), allocatable :: tmp_neighborhood_mad
        integer(int32), dimension(:, :), allocatable :: tmp_point_neighborhood_indices
        real(real64), dimension(:, :), allocatable :: tmp_neighbor_residuals
        integer(int32), dimension(:), allocatable :: tmp_point_pool
        integer(int32), dimension(:), allocatable :: tmp_point_draw
        integer(int32), dimension(:), allocatable :: tmp_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_high
        logical(c_bool), dimension(:), allocatable :: tmp_occupancy_failed
        integer(int32), dimension(:), allocatable :: tmp_n_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_min_bin_occupancy
        real(real64), dimension(:), allocatable :: tmp_mean_bin_occupancy
        integer(int32), dimension(:), allocatable :: tmp_max_bin_occupancy
        integer(int32), dimension(:), allocatable :: tmp_sturges_bins
        integer(int32), dimension(:), allocatable :: tmp_fd_bins
        real(real64), dimension(:, :, :), allocatable :: tmp_pmfs
        integer(int32), dimension(:, :, :), allocatable :: tmp_counts
        integer(int32), dimension(:, :), allocatable :: tmp_included_n_reps
        real(real64), dimension(:, :), allocatable :: tmp_mean_pmf
        integer(int32), dimension(:, :), allocatable :: tmp_mean_pmf_counts
        integer(int32), dimension(:), allocatable :: tmp_mean_pmf_included_n_reps
        real(real64), dimension(:, :), allocatable :: tmp_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_weights
        real(real64), dimension(:), allocatable :: tmp_global_js_divergence
        real(real64), dimension(:), allocatable :: tmp_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_residual_to_study_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts_search
        integer(int32), dimension(:), allocatable :: tmp_study_block_ends
        integer(int32), dimension(:), allocatable :: tmp_study_block_sizes
        real(real64), dimension(:), allocatable :: tmp_p_values
        real(real64), dimension(:), allocatable :: tmp_p_values_observed_consensus
        real(real64), dimension(:, :), allocatable :: tmp_permutation_mean_pmf
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_mean_pmf_counts
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_counts
        real(real64), dimension(:, :, :), allocatable :: tmp_permutation_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_permutation_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_permutation_weights
        integer(int32), dimension(:), allocatable :: tmp_search_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_search_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_search_shared_residual_range_high
        real(real64), dimension(:, :), allocatable :: tmp_search_confidence_interval
        integer(int32), dimension(:), allocatable :: tmp_search_trace_k_start
        integer(int32), dimension(:), allocatable :: tmp_search_trace_k_step
        integer(int32), dimension(:), allocatable :: tmp_search_trace_k_max
        integer(int32), dimension(:), allocatable :: tmp_search_trace_n_points
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_global_js_divergence
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_lower
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_upper
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_width
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_width_relative
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_ci_overlap
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_delta
        real(real64), dimension(:), allocatable :: tmp_search_trace_delta_median
        real(real64), dimension(:), allocatable :: tmp_search_trace_delta_max
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_selected_n_bins
        logical(c_bool), dimension(:, :), allocatable :: tmp_search_trace_occupancy_failed
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_n_pooled_residuals
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_min_bin_occupancy
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_mean_bin_occupancy
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_max_bin_occupancy
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_sturges_bins
        integer(int32), dimension(:, :), allocatable :: tmp_search_trace_fd_bins
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_shared_residual_range_low
        real(real64), dimension(:, :), allocatable :: tmp_search_trace_shared_residual_range_high
        integer(int32), dimension(:), allocatable :: tmp_search_candidate_k_start
        integer(int32), dimension(:), allocatable :: tmp_search_candidate_n_points
        integer(int32), dimension(:), allocatable :: tmp_search_candidate_status
        real(real64), dimension(:, :), allocatable :: tmp_candidate_confidence_interval
        real(real64), dimension(:, :, :), allocatable :: tmp_bootstrapping_top_k_jsds
        real(real64), dimension(:), allocatable :: tmp_prev_global_js_divergence
        integer(int32), dimension(:), allocatable :: tmp_delta_perm
        real(real64), dimension(:, :), allocatable :: tmp_best_uncertainty_confidence_interval
        integer(int32), dimension(:), allocatable :: tmp_best_n_bins_per_point
        integer(int32), dimension(:), allocatable :: tmp_best_uncertainty_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_best_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_best_shared_residual_range_high
        real(real64), dimension(:), allocatable :: tmp_best_uncertainty_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_best_uncertainty_shared_residual_range_high

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(k_start, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(k_step, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_int(k_max, ierr, arg_pos=7_int32, min=k_start)
        call validate_in_range_int(n_permutations, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=14_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=16_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=17_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=18_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=20_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=21_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=22_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=23_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=24_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=25_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=26_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=27_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=28_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=29_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(tau, ierr, arg_pos=31_int32, min=0.0_real64)
        call validate_in_range_real(mad_distance_factor, ierr, arg_pos=32_int32, min=0.0_real64)
        call validate_in_range_int(max_pooled_residuals, ierr, arg_pos=33_int32, min=0_int32)
        call validate_in_range_int(min_study_neighbors, ierr, arg_pos=34_int32, min=1_int32)
        call validate_in_range_real(delta_abs_tolerance, ierr, arg_pos=36_int32, min=0.0_real64)
        call validate_all_in_range_real(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        if (present(join_method)) then; if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=15_int32); end if
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=19_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=35_int32); end if
        if (present(bootstrap_mode)) then; if (bootstrap_mode /= MODE_BOOTSTRAP_PER_STUDY .and. bootstrap_mode /= MODE_BOOTSTRAP_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=37_int32); end if
        if (present(divergence_mode)) then; if (divergence_mode /= MODE_DIVERGENCE_SYMMETRIC .and. divergence_mode /= MODE_DIVERGENCE_PAIRWISE_K2 .and. divergence_mode /= MODE_DIVERGENCE_ONE_SIDED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=38_int32); end if
        if (present(divergence_weight_mode)) then; if (divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL .and. divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_POOL) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=39_int32); end if
        if (is_err(ierr)) return
#endif

        call calc_adaptive_js_comp_test_permutation_bounds(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_studies = n_studies,&
            k_start = k_start,&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            max_n_points_candidate = max_n_points_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            ierr = ierr&
        )
        if (is_err(ierr)) return
        M_ALLOCATE(tmp_sample_slots(max_n_reps_all_studies*n_studies))
        M_ALLOCATE(tmp_sample_assignment(max_n_reps_all_studies*n_studies))
        M_ALLOCATE(tmp_permuted_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_sample_mark(max_n_reps_all_studies*n_studies))
        M_ALLOCATE(tmp_permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_gene_means(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_x_star(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_pooled_neighborhood_range(2, max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_n_neighbors_per_point(n_studies, max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_stop_reason(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_neighborhood_dispersion(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_neighborhood_mad(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_point_neighborhood_indices(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_neighbor_residuals(max_n_reps_all_studies, max_n_genes_all_studies))
        M_ALLOCATE(tmp_point_pool(256))
        M_ALLOCATE(tmp_point_draw(256))
        M_ALLOCATE(tmp_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_occupancy_failed(max_n_points_candidate))
        M_ALLOCATE(tmp_n_pooled_residuals(max_n_points_candidate))
        M_ALLOCATE(tmp_min_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_mean_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_max_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_sturges_bins(max_n_points_candidate))
        M_ALLOCATE(tmp_fd_bins(max_n_points_candidate))
        M_ALLOCATE(tmp_pmfs(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_counts(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_included_n_reps(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_mean_pmf(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_mean_pmf_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        M_ALLOCATE(tmp_js_divergences(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_weights(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
        M_ALLOCATE(tmp_study_block_ends(n_studies))
        M_ALLOCATE(tmp_study_block_sizes(n_studies))
        M_ALLOCATE(tmp_p_values(n_studies))
        M_ALLOCATE(tmp_p_values_observed_consensus(n_studies))
        M_ALLOCATE(tmp_permutation_mean_pmf(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_mean_pmf_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_permutation_pmfs(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_permutation_js_divergences(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_permutation_weights(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_search_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_search_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_search_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_search_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_search_trace_k_start(16))
        M_ALLOCATE(tmp_search_trace_k_step(16))
        M_ALLOCATE(tmp_search_trace_k_max(16))
        M_ALLOCATE(tmp_search_trace_n_points(16))
        M_ALLOCATE(tmp_search_trace_global_js_divergence(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_lower(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_upper(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_width(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_width_relative(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_ci_overlap(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_delta(n_studies, 16))
        M_ALLOCATE(tmp_search_trace_delta_median(16))
        M_ALLOCATE(tmp_search_trace_delta_max(16))
        M_ALLOCATE(tmp_search_trace_selected_n_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_occupancy_failed(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_n_pooled_residuals(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_min_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_mean_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_max_bin_occupancy(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_sturges_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_fd_bins(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_shared_residual_range_low(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_trace_shared_residual_range_high(max_n_points_candidate, 16))
        M_ALLOCATE(tmp_search_candidate_k_start(16))
        M_ALLOCATE(tmp_search_candidate_n_points(16))
        M_ALLOCATE(tmp_search_candidate_status(16))
        M_ALLOCATE(tmp_candidate_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        M_ALLOCATE(tmp_prev_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_delta_perm(n_studies))
        M_ALLOCATE(tmp_best_uncertainty_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_best_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_best_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_best_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))

        call run_js_comp_test_adaptive_label_permutation_test_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            expression = expression,&
            k_start = k_start,&
            k_step = k_step,&
            k_max = k_max,&
            max_n_points_candidate = max_n_points_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            n_permutations = n_permutations,&
            global_js_divergence_observed = global_js_divergence_observed,&
            p_values = p_values,&
            jsd_null = jsd_null,&
            n_failed_permutations = n_failed_permutations,&
            tmp_sample_slots = tmp_sample_slots,&
            tmp_sample_assignment = tmp_sample_assignment,&
            tmp_permuted_global_js_divergence = tmp_permuted_global_js_divergence,&
            tmp_sample_mark = tmp_sample_mark,&
            tmp_permuted_expression = tmp_permuted_expression,&
            tmp_gene_means = tmp_gene_means,&
            tmp_residuals = tmp_residuals,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_x_star = tmp_x_star,&
            tmp_pooled_neighborhood_range = tmp_pooled_neighborhood_range,&
            tmp_n_neighbors_per_point = tmp_n_neighbors_per_point,&
            tmp_stop_reason = tmp_stop_reason,&
            tmp_neighborhood_dispersion = tmp_neighborhood_dispersion,&
            tmp_neighborhood_mad = tmp_neighborhood_mad,&
            tmp_point_neighborhood_indices = tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals = tmp_neighbor_residuals,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            tmp_n_bins_per_point = tmp_n_bins_per_point,&
            tmp_shared_residual_range_low = tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high = tmp_shared_residual_range_high,&
            tmp_occupancy_failed = tmp_occupancy_failed,&
            tmp_n_pooled_residuals = tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy = tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy = tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy = tmp_max_bin_occupancy,&
            tmp_sturges_bins = tmp_sturges_bins,&
            tmp_fd_bins = tmp_fd_bins,&
            tmp_pmfs = tmp_pmfs,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps = tmp_included_n_reps,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps = tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_p_values = tmp_p_values,&
            tmp_p_values_observed_consensus = tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_search_n_bins_per_point = tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low = tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high = tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval = tmp_search_confidence_interval,&
            tmp_search_trace_k_start = tmp_search_trace_k_start,&
            tmp_search_trace_k_step = tmp_search_trace_k_step,&
            tmp_search_trace_k_max = tmp_search_trace_k_max,&
            tmp_search_trace_n_points = tmp_search_trace_n_points,&
            tmp_search_trace_global_js_divergence = tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower = tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper = tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width = tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative = tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap = tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta = tmp_search_trace_delta,&
            tmp_search_trace_delta_median = tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max = tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins = tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed = tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals = tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy = tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy = tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy = tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins = tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins = tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low = tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high = tmp_search_trace_shared_residual_range_high,&
            tmp_search_candidate_k_start = tmp_search_candidate_k_start,&
            tmp_search_candidate_n_points = tmp_search_candidate_n_points,&
            tmp_search_candidate_status = tmp_search_candidate_status,&
            tmp_candidate_confidence_interval = tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence = tmp_prev_global_js_divergence,&
            tmp_delta_perm = tmp_delta_perm,&
            tmp_best_uncertainty_confidence_interval = tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point = tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point = tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low = tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high = tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low = tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high = tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search = rerun_parameter_search,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            bin_count_mode = bin_count_mode,&
            delta_abs_tolerance = delta_abs_tolerance,&
            bootstrap_mode = bootstrap_mode,&
            divergence_mode = divergence_mode,&
            divergence_weight_mode = divergence_weight_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_adaptive_label_permutation_test

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_label_permutation_test_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test_permutation(module):run_js_comp_test_adaptive_label_permutation_test]] does both.
    !| The adaptive counterpart of
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_label_permutation_test_impl(interface)]]:
    !| the same test (observed value, `n_permutations` shuffles of the study labels without
    !| replacement, p-value `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`), each
    !| data set run by
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]
    !| with the same data, knobs and options; its input layout and assignment convention apply here too.
    !|
    !| **Rerunning the parameter search** (`rerun_parameter_search = .true.`): every data set first
    !| runs the adaptive parameter search and then the construction and run with the knobs it
    !| selects. **Run time:** each permutation runs the whole JS-Comp-Test; with the rerun option it
    !| runs the whole parameter search first. The routine is sequential: a caller who wants
    !| permutations in parallel calls
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]
    !| itself with its own assignments (e.g. from R), one call per assignment.
    !|
    !| **Failures.** A permutation that fails with invalid input (a study left without residuals at
    !| some reference point, or a construction that builds no neighborhood or more neighborhoods
    !| than the work arrays hold) gets a NaN row in `jsd_null`, is counted in `n_failed_permutations`
    !| and does not count in `n_valid`. No valid permutation gives NaN p-values (a direct call with
    !| `n_permutations = 0` too; the published entry points reject 0, since it is the extent of
    !| `jsd_null`). Any other error ends the test and is returned, with NaN p-values; the rows of the
    !| permutations run so far are kept. An error on the observed data is returned at once
    !| (`global_js_divergence_observed`, `p_values` and `jsd_null` NaN), as is a study without any
    !| sample (the call fails with invalid input). In the Fortran expert tier, where the caller
    !| supplies the work-array capacity, a capacity below the search's own fails the same way with
    !| the rerun option.
    !|
    !| **Known limitation.** A gene missing in a whole study travels with the samples, so after
    !| shuffling that study may have values for it.
    subroutine run_js_comp_test_adaptive_label_permutation_test_expert(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            expression,&
            k_start,&
            k_step,&
            k_max,&
            max_n_points_candidate,&
            n_bootstrapping_top_k_jsds,&
            n_permutations,&
            global_js_divergence_observed,&
            p_values,&
            jsd_null,&
            n_failed_permutations,&
            tmp_sample_slots,&
            tmp_sample_assignment,&
            tmp_permuted_global_js_divergence,&
            tmp_sample_mark,&
            tmp_permuted_expression,&
            tmp_gene_means,&
            tmp_residuals,&
            tmp_gene_means_perm_all,&
            tmp_x_star,&
            tmp_pooled_neighborhood_range,&
            tmp_n_neighbors_per_point,&
            tmp_stop_reason,&
            tmp_neighborhood_dispersion,&
            tmp_neighborhood_mad,&
            tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals,&
            tmp_point_pool,&
            tmp_point_draw,&
            tmp_n_bins_per_point,&
            tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high,&
            tmp_occupancy_failed,&
            tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy,&
            tmp_sturges_bins,&
            tmp_fd_bins,&
            tmp_pmfs,&
            tmp_counts,&
            tmp_included_n_reps,&
            tmp_mean_pmf,&
            tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences,&
            tmp_weights,&
            tmp_global_js_divergence,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm,&
            tmp_bin_counts_search,&
            tmp_study_block_ends,&
            tmp_study_block_sizes,&
            tmp_p_values,&
            tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts,&
            tmp_permutation_pmfs,&
            tmp_permutation_js_divergences,&
            tmp_permutation_weights,&
            tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval,&
            tmp_search_trace_k_start,&
            tmp_search_trace_k_step,&
            tmp_search_trace_k_max,&
            tmp_search_trace_n_points,&
            tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta,&
            tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high,&
            tmp_search_candidate_k_start,&
            tmp_search_candidate_n_points,&
            tmp_search_candidate_status,&
            tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence,&
            tmp_delta_perm,&
            tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search,&
            n_bootstraps,&
            join_method,&
            min_residuals_per_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            tau,&
            mad_distance_factor,&
            max_pooled_residuals,&
            min_study_neighbors,&
            bin_count_mode,&
            delta_abs_tolerance,&
            bootstrap_mode,&
            divergence_mode,&
            divergence_weight_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays; a construction with more neighborhoods fails
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_permutation_impl(module):calc_adaptive_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! It is *VERY IMPORTANT* to compute this argument from the `n_bootstrapping_top_k_jsds` output produced by [[tox_data_integration_js_comp_test_permutation_impl(module):calc_adaptive_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_permutations
            !! Number of random assignments (rows of `jsd_null`)
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! The minimum valid value is `k_start`.
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence_observed
            !! Weighted global JSD of each study on the observed data (the identity assignment)
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
            !! whose value reached the observed one; NaN when none is valid or on error
        real(real64), dimension(n_permutations, n_studies), intent(out) :: jsd_null
            !! Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation
        integer(int32), intent(out) :: n_failed_permutations
            !! Number of permutations whose JS-Comp-Test failed with invalid input
        integer(int32), dimension(max_n_reps_all_studies*n_studies), intent(out) :: tmp_sample_slots
            !! Working array: the flat slot index of every sample
        integer(int32), dimension(max_n_reps_all_studies*n_studies), intent(out) :: tmp_sample_assignment
            !! Working array: the current assignment, in its first `n_samples` elements
        real(real64), dimension(n_studies), intent(out) :: tmp_permuted_global_js_divergence
            !! Working array: one permutation's global JSD per study
        logical(c_bool), dimension(max_n_reps_all_studies*n_studies), intent(out) :: tmp_sample_mark
            !! Working array: the values of the current assignment already seen
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_permuted_expression
            !! Working array: the relabelled data set
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means
            !! Working array: the relabelled data set's gene means per study
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_residuals
            !! Working array: the relabelled data set's residuals
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array of the search, the construction and the run: sorting permutation of the
            !! pooled gene means
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_x_star
            !! Working array of the search and the construction: the reference points
        integer(int32), dimension(2, max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_neighborhood_range
            !! Working array of the search and the construction: each neighborhood's pooled range,
            !! read by the run
        integer(int32), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out) :: tmp_n_neighbors_per_point
            !! Working array of the search, the construction and the run: per-study entry count per
            !! neighborhood
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_stop_reason
            !! Working array of the search and the construction: stop reason per neighborhood
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_dispersion
            !! Working array of the search and the construction: dispersion per neighborhood
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_mad
            !! Working array of the search and the construction: MAD per neighborhood
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_point_neighborhood_indices
            !! Working array of the search and the run: one reference point's gene indices per study
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out) :: tmp_neighbor_residuals
            !! Working array of the search and the run: one (point, study) pair's gathered residuals
        integer(int32), dimension(256), intent(out) :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(int32), dimension(256), intent(out) :: tmp_point_draw
            !! Working array of the search and the run: one study's draw at one point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array of the search and the run: bin count per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array of the search and the run: lower residual-range bound per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array of the search and the run: upper residual-range bound per reference point
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
            !! Working array of the search and the run: bin-count search failure flag per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
            !! Working array of the search and the run: pooled residual count per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
            !! Working array of the search and the run: minimum bin occupancy per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
            !! Working array of the search and the run: mean bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
            !! Working array of the search and the run: maximum bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
            !! Working array of the search and the run: Sturges' rule diagnostic per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_fd_bins
            !! Working array of the search and the run: Freedman-Diaconis diagnostic per reference point
        real(real64), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(int32), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(real64), dimension(max_n_points_candidate, 256), intent(out) :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values
            !! Working array: the run's p-values (NaN, its permutation test does not run)
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values_observed_consensus
            !! Working array: the run's comparison p-values (NaN, its permutation test does not run)
        real(real64), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_mean_pmf
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_permutation_counts
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_search_n_bins_per_point
            !! Working array: the search's chosen bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_low
            !! Working array: the search's chosen lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_high
            !! Working array: the search's chosen upper residual-range bounds
        real(real64), dimension(2, n_studies), intent(out) :: tmp_search_confidence_interval
            !! Working array: the search's chosen confidence interval
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_k_start
            !! Working array: the search's trace
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_k_step
            !! Working array: the search's trace
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_k_max
            !! Working array: the search's trace
        integer(int32), dimension(16), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(real64), dimension(n_studies, 16), intent(out) :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(real64), dimension(16), intent(out) :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(real64), dimension(16), intent(out) :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: tmp_search_trace_shared_residual_range_high
            !! Working array: the search's trace
        integer(int32), dimension(16), intent(out) :: tmp_search_candidate_k_start
            !! Working array: the search's candidate log
        integer(int32), dimension(16), intent(out) :: tmp_search_candidate_n_points
            !! Working array: the search's candidate log
        integer(int32), dimension(16), intent(out) :: tmp_search_candidate_status
            !! Working array: the search's candidate log
        real(real64), dimension(2, n_studies), intent(out) :: tmp_candidate_confidence_interval
            !! Working array of the search: the current candidate's confidence interval
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array of the search: the bootstrap's bottom-k/top-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array of the search: the previous admissible candidate's global JSD
        integer(int32), dimension(n_studies), intent(out) :: tmp_delta_perm
            !! Working array of the search: sorting permutation for the per-study medians
        real(real64), dimension(2, n_studies), intent(out) :: tmp_best_uncertainty_confidence_interval
            !! Working array of the search: the smallest-uncertainty candidate's confidence interval
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_n_bins_per_point
            !! Working array of the search: the best candidate's bin counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_n_bins_per_point
            !! Working array of the search: the smallest-uncertainty candidate's bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_low
            !! Working array of the search: the best candidate's lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_high
            !! Working array of the search: the best candidate's upper residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array of the search: the smallest-uncertainty candidate's lower bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array of the search: the smallest-uncertainty candidate's upper bounds
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! `.true.`: run the adaptive parameter search on the relabelled data first and use its
            !! knobs; `.false.`: use `k_start`, `k_step` and `k_max`
            !! The default value is `.false.`.
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `1_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed of the random number generator of the label draws, and of the bootstrap of
            !! every rerun parameter search (a separate generator)
            !! The default value is `42_int32`.
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have (construction status; in
            !! the rerun search, candidate admissibility)
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
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
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0_real64`.
        integer(int32), intent(in), optional :: bootstrap_mode
            !! Rerun search (used only with `rerun_parameter_search`): what its bootstrap draws each
            !! study's histograms from. By default each study's own counts: a robustness interval of
            !! the observed JSD. The pooled counts of all studies give a noise band of "no difference
            !! between the studies" instead, kept for comparison
            !!
            !! | Mode                                                                                             | Value                                                                          |
            !! |--------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
            !! | Each study's own counts, a robustness interval (default)                                         | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_PER_STUDY(variable)]] |
            !! | The pooled consensus counts of all studies, a noise band (the bootstrap before the mode existed) | [[tox_data_integration_stats_impl(module):MODE_BOOTSTRAP_POOLED(variable)]]    |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_mode
            !! Experiment option: how each study's per-point divergence is computed. By default the
            !! JSD of the study against the consensus pmf, the formula before the option existed
            !!
            !! | Mode                                                                                                                | Value                                                                           |
            !! |---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | JSD of each study against the consensus pmf (default)                                                               | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | For exactly two studies the pairwise JSD of the two, the same value for both (otherwise as the default)             | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! | One-sided Kullback-Leibler divergence of each study against the consensus pmf, at most log2(number of studies) bits | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_ONE_SIDED(variable)]]   |
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: divergence_weight_mode
            !! Experiment option: how each reference point's weight in a study's global divergence is
            !! computed. By default the study's residual count there plus the count of all studies
            !! pooled, the weights before the option existed. For two studies in the pairwise
            !! `divergence_mode` the weights are always the two studies' counts
            !!
            !! | Mode                                                                   | Value                                                                                      |
            !! |------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
            !! | The study's count plus the pooled count of all studies (default)       | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL(variable)]] |
            !! | The pooled count of all studies only, the same weights for every study | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_WEIGHT_POOL(variable)]]            |
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(k_start, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(k_step, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_int(k_max, ierr, arg_pos=7_int32, min=k_start)
        call validate_in_range_int(max_n_points_candidate, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_int(n_bootstrapping_top_k_jsds, ierr, arg_pos=9_int32, min=1_int32)
        call validate_in_range_int(n_permutations, ierr, arg_pos=10_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=108_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=110_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=111_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=112_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=114_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=115_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=116_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=117_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=118_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=119_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=120_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=121_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=122_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=123_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(tau, ierr, arg_pos=125_int32, min=0.0_real64)
        call validate_in_range_real(mad_distance_factor, ierr, arg_pos=126_int32, min=0.0_real64)
        call validate_in_range_int(max_pooled_residuals, ierr, arg_pos=127_int32, min=0_int32)
        call validate_in_range_int(min_study_neighbors, ierr, arg_pos=128_int32, min=1_int32)
        call validate_in_range_real(delta_abs_tolerance, ierr, arg_pos=130_int32, min=0.0_real64)
        call validate_all_in_range_real(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        if (present(join_method)) then; if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=109_int32); end if
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=113_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=129_int32); end if
        if (present(bootstrap_mode)) then; if (bootstrap_mode /= MODE_BOOTSTRAP_PER_STUDY .and. bootstrap_mode /= MODE_BOOTSTRAP_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=131_int32); end if
        if (present(divergence_mode)) then; if (divergence_mode /= MODE_DIVERGENCE_SYMMETRIC .and. divergence_mode /= MODE_DIVERGENCE_PAIRWISE_K2 .and. divergence_mode /= MODE_DIVERGENCE_ONE_SIDED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=132_int32); end if
        if (present(divergence_weight_mode)) then; if (divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_STUDY_PLUS_POOL .and. divergence_weight_mode /= MODE_DIVERGENCE_WEIGHT_POOL) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=133_int32); end if
        if (is_err(ierr)) return
#endif

        call run_js_comp_test_adaptive_label_permutation_test_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            expression = expression,&
            k_start = k_start,&
            k_step = k_step,&
            k_max = k_max,&
            max_n_points_candidate = max_n_points_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            n_permutations = n_permutations,&
            global_js_divergence_observed = global_js_divergence_observed,&
            p_values = p_values,&
            jsd_null = jsd_null,&
            n_failed_permutations = n_failed_permutations,&
            tmp_sample_slots = tmp_sample_slots,&
            tmp_sample_assignment = tmp_sample_assignment,&
            tmp_permuted_global_js_divergence = tmp_permuted_global_js_divergence,&
            tmp_sample_mark = tmp_sample_mark,&
            tmp_permuted_expression = tmp_permuted_expression,&
            tmp_gene_means = tmp_gene_means,&
            tmp_residuals = tmp_residuals,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_x_star = tmp_x_star,&
            tmp_pooled_neighborhood_range = tmp_pooled_neighborhood_range,&
            tmp_n_neighbors_per_point = tmp_n_neighbors_per_point,&
            tmp_stop_reason = tmp_stop_reason,&
            tmp_neighborhood_dispersion = tmp_neighborhood_dispersion,&
            tmp_neighborhood_mad = tmp_neighborhood_mad,&
            tmp_point_neighborhood_indices = tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals = tmp_neighbor_residuals,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            tmp_n_bins_per_point = tmp_n_bins_per_point,&
            tmp_shared_residual_range_low = tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high = tmp_shared_residual_range_high,&
            tmp_occupancy_failed = tmp_occupancy_failed,&
            tmp_n_pooled_residuals = tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy = tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy = tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy = tmp_max_bin_occupancy,&
            tmp_sturges_bins = tmp_sturges_bins,&
            tmp_fd_bins = tmp_fd_bins,&
            tmp_pmfs = tmp_pmfs,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps = tmp_included_n_reps,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps = tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_p_values = tmp_p_values,&
            tmp_p_values_observed_consensus = tmp_p_values_observed_consensus,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_search_n_bins_per_point = tmp_search_n_bins_per_point,&
            tmp_search_shared_residual_range_low = tmp_search_shared_residual_range_low,&
            tmp_search_shared_residual_range_high = tmp_search_shared_residual_range_high,&
            tmp_search_confidence_interval = tmp_search_confidence_interval,&
            tmp_search_trace_k_start = tmp_search_trace_k_start,&
            tmp_search_trace_k_step = tmp_search_trace_k_step,&
            tmp_search_trace_k_max = tmp_search_trace_k_max,&
            tmp_search_trace_n_points = tmp_search_trace_n_points,&
            tmp_search_trace_global_js_divergence = tmp_search_trace_global_js_divergence,&
            tmp_search_trace_ci_lower = tmp_search_trace_ci_lower,&
            tmp_search_trace_ci_upper = tmp_search_trace_ci_upper,&
            tmp_search_trace_ci_width = tmp_search_trace_ci_width,&
            tmp_search_trace_ci_width_relative = tmp_search_trace_ci_width_relative,&
            tmp_search_trace_ci_overlap = tmp_search_trace_ci_overlap,&
            tmp_search_trace_delta = tmp_search_trace_delta,&
            tmp_search_trace_delta_median = tmp_search_trace_delta_median,&
            tmp_search_trace_delta_max = tmp_search_trace_delta_max,&
            tmp_search_trace_selected_n_bins = tmp_search_trace_selected_n_bins,&
            tmp_search_trace_occupancy_failed = tmp_search_trace_occupancy_failed,&
            tmp_search_trace_n_pooled_residuals = tmp_search_trace_n_pooled_residuals,&
            tmp_search_trace_min_bin_occupancy = tmp_search_trace_min_bin_occupancy,&
            tmp_search_trace_mean_bin_occupancy = tmp_search_trace_mean_bin_occupancy,&
            tmp_search_trace_max_bin_occupancy = tmp_search_trace_max_bin_occupancy,&
            tmp_search_trace_sturges_bins = tmp_search_trace_sturges_bins,&
            tmp_search_trace_fd_bins = tmp_search_trace_fd_bins,&
            tmp_search_trace_shared_residual_range_low = tmp_search_trace_shared_residual_range_low,&
            tmp_search_trace_shared_residual_range_high = tmp_search_trace_shared_residual_range_high,&
            tmp_search_candidate_k_start = tmp_search_candidate_k_start,&
            tmp_search_candidate_n_points = tmp_search_candidate_n_points,&
            tmp_search_candidate_status = tmp_search_candidate_status,&
            tmp_candidate_confidence_interval = tmp_candidate_confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence = tmp_prev_global_js_divergence,&
            tmp_delta_perm = tmp_delta_perm,&
            tmp_best_uncertainty_confidence_interval = tmp_best_uncertainty_confidence_interval,&
            tmp_best_n_bins_per_point = tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point = tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low = tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high = tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low = tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high = tmp_best_uncertainty_shared_residual_range_high,&
            rerun_parameter_search = rerun_parameter_search,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            bin_count_mode = bin_count_mode,&
            delta_abs_tolerance = delta_abs_tolerance,&
            bootstrap_mode = bootstrap_mode,&
            divergence_mode = divergence_mode,&
            divergence_weight_mode = divergence_weight_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_adaptive_label_permutation_test_expert

end module tox_data_integration_js_comp_test_permutation
