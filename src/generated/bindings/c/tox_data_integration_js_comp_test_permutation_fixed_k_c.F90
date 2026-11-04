#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_permutation_fixed_k(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Study-Label Permutation-Test, k-NN
!|
!| Asks whether the assignment of samples to studies (batches) influences the JS-Comp-Test values
!| at all. Null hypothesis: it does not. The test shuffles the study labels of the samples
!| (without replacement), reruns the whole JS-Comp-Test on every shuffled data set, and compares
!| each study's observed global JSD with the values the shuffled data sets give. Study `i`'s p-value is
!| `(count + 1)/(n_valid + 1)`, where `count` is how many of the `n_valid` valid shuffled data sets give
!| study `i` a global JSD at least as large as the observed one. The +1 counts the observed data itself.
!|
!| **Input layout.** All studies arrive in one array `expression(max_n_reps_all_studies,
!| max_n_genes_all_studies, n_studies)`: gene `g` is the same gene in every study, a gene missing in
!| a study is NaN there, and a replicate slot `(r, s)` whose values are all NaN is padding (a study
!| with fewer samples than `max_n_reps_all_studies`). Every other slot is a sample. Samples are
!| numbered in slot order, replicate fastest, then study.
!|
!| **Entry points.** One assignment of samples to slots, one k-NN JS-Comp-Test (the fixed
!| `(n_points, n_neighbors)` of `run_js_comp_test`):
!| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]].
!| The whole test, `n_permutations` random assignments plus the observed one:
!| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_permutation_impl(interface)]].
!| The adaptive counterparts are in [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module)]], the shared helpers in [[tox_data_integration_js_comp_test_permutation_common_impl(module)]].
!| `calc_js_comp_test_permutation_bounds` and `calc_adaptive_js_comp_test_permutation_bounds`
!| ([[tox_data_integration_js_comp_test_permutation_bounds_impl(module)]]) size the work arrays for a caller that allocates its own.
!|
!| The adaptive routines are the recommended way, as for the JS-Comp-Test itself; the k-NN ones
!| stay available.
module tox_data_integration_js_comp_test_permutation_fixed_k_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_as_view
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: run_js_comp_test_on_sample_assignment_c
    public :: run_js_comp_test_on_sample_assignment_expert_c
    public :: run_js_comp_test_permutation_c
    public :: run_js_comp_test_permutation_expert_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_permutation_fixed_k(module):run_js_comp_test_on_sample_assignment(subroutine)]]
    !| One building block of the study-label Permutation-Test: relabels the samples by `sample_assignment` and runs
    !| the whole JS-Comp-Test on the relabelled data set. `expression` holds every study (input layout as in the module
    !| description). `sample_assignment` is a permutation of `1..n_samples`: destination sample `k` receives all genes of source
    !| sample `sample_assignment(k)`, so every study keeps its slots, sample count and padding. The identity is the observed
    !| data, and its result is exactly that of a caller doing the steps below on the observed data.
    !|
    !| Steps: per study, the gene means and residuals of the relabelled data
    !| ([[tox_data_integration_preprocessing_impl(module):compute_gene_means_impl(interface)]], [[tox_data_integration_preprocessing_impl(module):compute_residuals_impl(interface)]]); the reference points from
    !| the pooled means ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]); then [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]] with the
    !| bin-count options and `divergence_mode` as given.
    !|
    !| With `rerun_parameter_search = .true.` the relabelled data first run
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| with the search options below, and the run uses the `(n_points, n_neighbors)` it returns (also its fallback when no
    !| candidate was admissible). `n_points` and `n_neighbors` then only size the work arrays, and `used_n_points` and
    !| `used_n_neighbors` report what the run used. Every call runs the whole JS-Comp-Test, with the rerun option the
    !| whole parameter search first.
    !|
    !| Errors (`ierr`): invalid input (`ERR_INVALID_INPUT`) when some study has no sample at all, when `n_samples` is not the
    !| number of samples in `expression`, or when `sample_assignment` is not a permutation of `1..n_samples`. In the Fortran expert tier, with the rerun option, also when a capacity is below the
    !| search's own bound (`calc_js_comp_test_candidate_bounds_impl` for the same `grid_points_factor` and `grid_min_points`).
    !| Errors of the search and the run are returned with their code, and `global_js_divergence` is then NaN.
    !| Known limitation: a gene missing in a whole study travels with the samples, so after relabelling that study may have
    !| values for it.
    subroutine run_js_comp_test_on_sample_assignment_c(&
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
            divergence_mode,&
            grid_points_factor,&
            grid_min_points,&
            ierr&
        ) bind(C, name="run_js_comp_test_on_sample_assignment_c")
        use tox_data_integration_js_comp_test_permutation_fixed_k, only: run_js_comp_test_on_sample_assignment
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        use tox_data_integration_js_comp_test_fixed_k_search_impl, only: MODE_NEIGHBORHOOD_SIZE_LARGE, MODE_NEIGHBORHOOD_SIZE_SMALL
        use tox_data_integration_js_comp_test_search_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE
        use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_samples
            !! Number of samples, the non-padding slots of `expression`
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_samples), intent(in), target :: sample_assignment
            !! A permutation of `1..n_samples`: destination sample `k` receives source sample
            !! `sample_assignment(k)`; the identity is the observed data
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_samples`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study, the one `divergence_mode` selects (by default for exactly
            !! two studies their pairwise JSD, otherwise against the consensus), for this assignment; NaN on
            !! error
        integer(c_int), intent(out), target :: used_n_points
            !! Number of reference points the run used (`n_points`, or the search's choice)
        integer(c_int), intent(out), target :: used_n_neighbors
            !! Number of neighbors per neighborhood the run used (`n_neighbors`, or the search's choice)
        logical(c_bool), intent(in), target :: rerun_parameter_search
            !! `.true.`: run the fixed-k parameter search on the relabelled data first and use its
            !! `(n_points, n_neighbors)`; `.false.`: use `n_points` and `n_neighbors`
            !! The default value is `.false.`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                                  | Value                                                                                  |
            !! |---------------------------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass)                 | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)                  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (at least half of the studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `'join_min'`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        character(len=1, kind=c_char), dimension(19), intent(in), target :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `'plateau_effect_size'`.
        real(c_double), intent(in), target :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Rerun search: seed of the bootstrap's random number generator
            !! The default value is `42_int32`.
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
        character(len=1, kind=c_char), dimension(23), intent(in), target :: neighborhood_size_mode
            !! Rerun search: which neighborhood size the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are not the
            !! default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode                                                  | Value                                                                                                    |
            !! |-------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Large neighborhoods, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods, `KX_FACTORS(2)` = 0.5            | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! The default value is `'neighborhood_size_large'`.
        real(c_double), intent(in), target :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0001_real64`.
        character(len=1, kind=c_char), dimension(22), intent(in), target :: divergence_mode
            !! How each study's per-point divergence is computed, in the run (observed and every
            !! shuffled data set) and in the rerun search
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `'divergence_pairwise_k2'`.
        real(c_double), intent(in), target :: grid_points_factor
            !! Rerun search (used only with `rerun_parameter_search`): factor in the k-NN candidate
            !! grid's first reference-point count `ceiling(grid_points_factor * sqrt(max_n_genes_all_studies))`
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `4.0_real64`.
        integer(c_int), intent(in), target :: grid_min_points
            !! Rerun search (used only with `rerun_parameter_search`): the k-NN candidate grid's
            !! floor, its last and coarsest `n_points` candidate
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_POINTS`.
            !! The default value is `200_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)
        integer(int32) :: join_method_mode_f
        integer(int32) :: plateau_mode_mode_f
        integer(int32) :: bin_count_mode_mode_f
        integer(int32) :: neighborhood_size_mode_mode_f
        integer(int32) :: divergence_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_samples)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(used_n_points)
        M_CHECK_NON_NULL(used_n_neighbors)
        M_CHECK_NON_NULL(rerun_parameter_search)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(min_neighbor_overlap)
        M_CHECK_NON_NULL(succeeding_ci_overlap)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_NON_NULL(delta_abs_tolerance)
        M_CHECK_NON_NULL(grid_points_factor)
        M_CHECK_NON_NULL(grid_min_points)
        M_CHECK_ARRAY_NON_NULL(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(sample_assignment, n_samples)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)
        M_CHECK_ARRAY_NON_NULL(plateau_mode, 19)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)
        M_CHECK_ARRAY_NON_NULL(neighborhood_size_mode, 23)
        M_CHECK_ARRAY_NON_NULL(divergence_mode, 22)

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
            character(len=:), pointer :: neighborhood_size_mode_f
            neighborhood_size_mode_f => c_char_as_view(neighborhood_size_mode)

            select case (neighborhood_size_mode_f)
                case ("neighborhood_size_large")
                    neighborhood_size_mode_mode_f = MODE_NEIGHBORHOOD_SIZE_LARGE
                case ("neighborhood_size_small")
                    neighborhood_size_mode_mode_f = MODE_NEIGHBORHOOD_SIZE_SMALL
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

        call run_js_comp_test_on_sample_assignment(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            expression = expression,&
            n_samples = n_samples,&
            sample_assignment = sample_assignment,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            global_js_divergence = global_js_divergence,&
            used_n_points = used_n_points,&
            used_n_neighbors = used_n_neighbors,&
            rerun_parameter_search = rerun_parameter_search,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method_mode_f,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode_mode_f,&
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
            bin_count_mode = bin_count_mode_mode_f,&
            neighborhood_size_mode = neighborhood_size_mode_mode_f,&
            delta_abs_tolerance = delta_abs_tolerance,&
            divergence_mode = divergence_mode_mode_f,&
            grid_points_factor = grid_points_factor,&
            grid_min_points = grid_min_points,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_on_sample_assignment_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_permutation_fixed_k(module):run_js_comp_test_on_sample_assignment_expert(subroutine)]]
    !| One building block of the study-label Permutation-Test: relabels the samples by `sample_assignment` and runs
    !| the whole JS-Comp-Test on the relabelled data set. `expression` holds every study (input layout as in the module
    !| description). `sample_assignment` is a permutation of `1..n_samples`: destination sample `k` receives all genes of source
    !| sample `sample_assignment(k)`, so every study keeps its slots, sample count and padding. The identity is the observed
    !| data, and its result is exactly that of a caller doing the steps below on the observed data.
    !|
    !| Steps: per study, the gene means and residuals of the relabelled data
    !| ([[tox_data_integration_preprocessing_impl(module):compute_gene_means_impl(interface)]], [[tox_data_integration_preprocessing_impl(module):compute_residuals_impl(interface)]]); the reference points from
    !| the pooled means ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]); then [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]] with the
    !| bin-count options and `divergence_mode` as given.
    !|
    !| With `rerun_parameter_search = .true.` the relabelled data first run
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| with the search options below, and the run uses the `(n_points, n_neighbors)` it returns (also its fallback when no
    !| candidate was admissible). `n_points` and `n_neighbors` then only size the work arrays, and `used_n_points` and
    !| `used_n_neighbors` report what the run used. Every call runs the whole JS-Comp-Test, with the rerun option the
    !| whole parameter search first.
    !|
    !| Errors (`ierr`): invalid input (`ERR_INVALID_INPUT`) when some study has no sample at all, when `n_samples` is not the
    !| number of samples in `expression`, or when `sample_assignment` is not a permutation of `1..n_samples`. In the Fortran expert tier, with the rerun option, also when a capacity is below the
    !| search's own bound (`calc_js_comp_test_candidate_bounds_impl` for the same `grid_points_factor` and `grid_min_points`).
    !| Errors of the search and the run are returned with their code, and `global_js_divergence` is then NaN.
    !| Known limitation: a gene missing in a whole study travels with the samples, so after relabelling that study may have
    !| values for it.
    subroutine run_js_comp_test_on_sample_assignment_expert_c(&
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
            divergence_mode,&
            grid_points_factor,&
            grid_min_points,&
            ierr&
        ) bind(C, name="run_js_comp_test_on_sample_assignment_expert_c")
        use tox_data_integration_js_comp_test_permutation_fixed_k, only: run_js_comp_test_on_sample_assignment_expert
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        use tox_data_integration_js_comp_test_fixed_k_search_impl, only: MODE_NEIGHBORHOOD_SIZE_LARGE, MODE_NEIGHBORHOOD_SIZE_SMALL
        use tox_data_integration_js_comp_test_search_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE
        use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_samples
            !! Number of samples, the non-padding slots of `expression`
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_points_candidate
            !! Capacity of the per-point work arrays
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_permutation_bounds(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `n_points`.
        integer(c_int), intent(in), target :: max_n_neighbors_candidate
            !! Capacity of the per-neighbor work arrays
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_permutation_bounds(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `n_neighbors`.
        integer(c_int), intent(in), target :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! It is *VERY IMPORTANT* to compute this argument from the `n_bootstrapping_top_k_jsds` output produced by [[tox_data_integration_js_comp_test_permutation_bounds(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_samples), intent(in), target :: sample_assignment
            !! A permutation of `1..n_samples`: destination sample `k` receives source sample
            !! `sample_assignment(k)`; the identity is the observed data
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_samples`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study, the one `divergence_mode` selects (by default for exactly
            !! two studies their pairwise JSD, otherwise against the consensus), for this assignment; NaN on
            !! error
        integer(c_int), intent(out), target :: used_n_points
            !! Number of reference points the run used (`n_points`, or the search's choice)
        integer(c_int), intent(out), target :: used_n_neighbors
            !! Number of neighbors per neighborhood the run used (`n_neighbors`, or the search's choice)
        integer(c_int), dimension(max_n_reps_all_studies*n_studies), intent(out), target :: tmp_sample_slots
            !! Working array: the flat slot index of every sample
        logical(c_bool), dimension(n_samples), intent(out), target :: tmp_sample_mark
            !! Working array: the values of `sample_assignment` already seen
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out), target :: tmp_permuted_expression
            !! Working array: the relabelled data set
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(out), target :: tmp_gene_means
            !! Working array: the relabelled data set's gene means per study
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out), target :: tmp_residuals
            !! Working array: the relabelled data set's residuals
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(out), target :: tmp_gene_means_perms
            !! Working array of the search and the run: each study's sorting permutation of its gene means
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_gene_means_perm_all
            !! Working array of the search and the run: sorting permutation of the pooled gene means
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_x_star
            !! Working array of the search and the run: the reference points
        integer(c_int), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out), target :: tmp_neighborhood_indices
            !! Working array of the search and the run: every study's neighbor gene indices
        integer(c_int), dimension(2, max_n_points_candidate, n_studies), intent(out), target :: tmp_neighborhood_range
            !! Working array of the search and the run: neighborhood spans
        real(c_double), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), intent(out), target :: tmp_neighborhood_residuals_gathered
            !! Working array of the search and the run: one study's gathered neighborhood residuals
        integer(c_int), dimension(256), intent(out), target :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(c_int), dimension(256), intent(out), target :: tmp_point_draw
            !! Working array of the search and the run: one study's draw at one point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_n_bins_per_point
            !! Working array of the search and the run: bin count per reference point
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_shared_residual_range_low
            !! Working array of the search and the run: lower residual-range bound per reference point
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_shared_residual_range_high
            !! Working array of the search and the run: upper residual-range bound per reference point
        logical(c_bool), dimension(max_n_points_candidate), intent(out), target :: tmp_occupancy_failed
            !! Working array of the search and the run: bin-count search failure flag per reference point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_n_pooled_residuals
            !! Working array of the search and the run: pooled residual count per reference point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_min_bin_occupancy
            !! Working array of the search and the run: minimum bin occupancy per reference point
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_mean_bin_occupancy
            !! Working array of the search and the run: mean bin occupancy per reference point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_max_bin_occupancy
            !! Working array of the search and the run: maximum bin occupancy per reference point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_sturges_bins
            !! Working array of the search and the run: Sturges' rule diagnostic per reference point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_fd_bins
            !! Working array of the search and the run: Freedman-Diaconis diagnostic per reference point
        real(c_double), dimension(max_n_points_candidate, 256, n_studies), intent(out), target :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(c_int), dimension(max_n_points_candidate, 256, n_studies), intent(out), target :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(c_int), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(c_double), dimension(max_n_points_candidate, 256), intent(out), target :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(c_int), dimension(max_n_points_candidate, 256), intent(out), target :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(c_double), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(c_double), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(c_double), dimension(n_studies), intent(out), target :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(c_double), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out), target :: tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(c_int), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out), target :: tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(c_int), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out), target :: tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_point_n_neighbors
            !! Working array of the search and the run: per-study neighbor counts, then sort scratch
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_search_n_bins_per_point
            !! Working array: the search's chosen bin counts
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_search_shared_residual_range_low
            !! Working array: the search's chosen lower residual-range bounds
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_search_shared_residual_range_high
            !! Working array: the search's chosen upper residual-range bounds
        real(c_double), dimension(2, n_studies), intent(out), target :: tmp_search_confidence_interval
            !! Working array: the search's chosen confidence interval
        integer(c_int), dimension(16), intent(out), target :: tmp_search_trace_n_points
            !! Working array: the search's trace
        integer(c_int), dimension(16), intent(out), target :: tmp_search_trace_n_neighbors
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(c_double), dimension(16), intent(out), target :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(c_double), dimension(16), intent(out), target :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_shared_residual_range_high
            !! Working array: the search's trace
        real(c_double), dimension(2, n_studies), intent(out), target :: tmp_candidate_confidence_interval
            !! Working array of the search: the current candidate's confidence interval
        real(c_double), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out), target :: tmp_bootstrapping_top_k_jsds
            !! Working array of the search: the bootstrap's bottom-k/top-k heaps
        real(c_double), dimension(n_studies), intent(out), target :: tmp_prev_global_js_divergence
            !! Working array of the search: the previous admissible candidate's global JSD
        real(c_double), dimension(2, n_studies), intent(out), target :: tmp_best_uncertainty_confidence_interval
            !! Working array of the search: the smallest-uncertainty candidate's confidence interval
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_best_n_bins_per_point
            !! Working array of the search: the best candidate's bin counts
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_n_bins_per_point
            !! Working array of the search: the smallest-uncertainty candidate's bin counts
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_shared_residual_range_low
            !! Working array of the search: the best candidate's lower residual-range bounds
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_shared_residual_range_high
            !! Working array of the search: the best candidate's upper residual-range bounds
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array of the search: the smallest-uncertainty candidate's lower bounds
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array of the search: the smallest-uncertainty candidate's upper bounds
        logical(c_bool), intent(in), target :: rerun_parameter_search
            !! `.true.`: run the fixed-k parameter search on the relabelled data first and use its
            !! `(n_points, n_neighbors)`; `.false.`: use `n_points` and `n_neighbors`
            !! The default value is `.false.`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                                  | Value                                                                                  |
            !! |---------------------------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass)                 | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)                  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (at least half of the studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `'join_min'`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        character(len=1, kind=c_char), dimension(19), intent(in), target :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `'plateau_effect_size'`.
        real(c_double), intent(in), target :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Rerun search: seed of the bootstrap's random number generator
            !! The default value is `42_int32`.
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
        character(len=1, kind=c_char), dimension(23), intent(in), target :: neighborhood_size_mode
            !! Rerun search: which neighborhood size the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are not the
            !! default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode                                                  | Value                                                                                                    |
            !! |-------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Large neighborhoods, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods, `KX_FACTORS(2)` = 0.5            | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! The default value is `'neighborhood_size_large'`.
        real(c_double), intent(in), target :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0001_real64`.
        character(len=1, kind=c_char), dimension(22), intent(in), target :: divergence_mode
            !! How each study's per-point divergence is computed, in the run (observed and every
            !! shuffled data set) and in the rerun search
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `'divergence_pairwise_k2'`.
        real(c_double), intent(in), target :: grid_points_factor
            !! Rerun search (used only with `rerun_parameter_search`): factor in the k-NN candidate
            !! grid's first reference-point count `ceiling(grid_points_factor * sqrt(max_n_genes_all_studies))`
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `4.0_real64`.
        integer(c_int), intent(in), target :: grid_min_points
            !! Rerun search (used only with `rerun_parameter_search`): the k-NN candidate grid's
            !! floor, its last and coarsest `n_points` candidate
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_POINTS`.
            !! The default value is `200_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)
        integer(int32) :: join_method_mode_f
        integer(int32) :: plateau_mode_mode_f
        integer(int32) :: bin_count_mode_mode_f
        integer(int32) :: neighborhood_size_mode_mode_f
        integer(int32) :: divergence_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_samples)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(max_n_neighbors_candidate)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)
        M_CHECK_NON_NULL(used_n_points)
        M_CHECK_NON_NULL(used_n_neighbors)
        M_CHECK_NON_NULL(rerun_parameter_search)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(min_neighbor_overlap)
        M_CHECK_NON_NULL(succeeding_ci_overlap)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_NON_NULL(delta_abs_tolerance)
        M_CHECK_NON_NULL(grid_points_factor)
        M_CHECK_NON_NULL(grid_min_points)
        M_CHECK_ARRAY_NON_NULL(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(sample_assignment, n_samples)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_sample_slots, (max_n_reps_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_sample_mark, n_samples)
        M_CHECK_ARRAY_NON_NULL(tmp_permuted_expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perms, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perm_all, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_x_star, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_indices, max_n_neighbors_candidate * max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_range, 2 * max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_residuals_gathered, max_n_reps_all_studies * max_n_neighbors_candidate * max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_point_pool, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_point_draw, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_occupancy_failed, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_n_pooled_residuals, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_min_bin_occupancy, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_bin_occupancy, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_max_bin_occupancy, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_sturges_bins, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_fd_bins, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, max_n_points_candidate * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, max_n_points_candidate * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_included_n_reps, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf, max_n_points_candidate * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf_counts, max_n_points_candidate * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf_included_n_reps, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals, (max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals_perm, (max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_residual_to_study_perm, (max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_bin_counts_search, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_ends, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_sizes, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_point_n_neighbors, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_search_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_search_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_search_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_search_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_n_points, 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_n_neighbors, 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_global_js_divergence, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_ci_lower, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_ci_upper, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_ci_width, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_ci_width_relative, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_ci_overlap, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_delta, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_delta_median, 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_delta_max, 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_selected_n_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_occupancy_failed, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_n_pooled_residuals, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_min_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_mean_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_max_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_sturges_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_fd_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_shared_residual_range_low, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_shared_residual_range_high, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_candidate_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds * 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_prev_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_best_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)
        M_CHECK_ARRAY_NON_NULL(plateau_mode, 19)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)
        M_CHECK_ARRAY_NON_NULL(neighborhood_size_mode, 23)
        M_CHECK_ARRAY_NON_NULL(divergence_mode, 22)

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
            character(len=:), pointer :: neighborhood_size_mode_f
            neighborhood_size_mode_f => c_char_as_view(neighborhood_size_mode)

            select case (neighborhood_size_mode_f)
                case ("neighborhood_size_large")
                    neighborhood_size_mode_mode_f = MODE_NEIGHBORHOOD_SIZE_LARGE
                case ("neighborhood_size_small")
                    neighborhood_size_mode_mode_f = MODE_NEIGHBORHOOD_SIZE_SMALL
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

        call run_js_comp_test_on_sample_assignment_expert(&
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
            join_method = join_method_mode_f,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode_mode_f,&
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
            bin_count_mode = bin_count_mode_mode_f,&
            neighborhood_size_mode = neighborhood_size_mode_mode_f,&
            delta_abs_tolerance = delta_abs_tolerance,&
            divergence_mode = divergence_mode_mode_f,&
            grid_points_factor = grid_points_factor,&
            grid_min_points = grid_min_points,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_on_sample_assignment_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_permutation_fixed_k(module):run_js_comp_test_permutation(subroutine)]]
    !| Tests whether the assignment of samples to studies influences each study's JS-Comp-Test value (null
    !| hypothesis: it does not), with the data, parameters and options of
    !| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]].
    !|
    !| 1. The observed data (the identity assignment) give `global_js_divergence_observed`.
    !| 2. `n_permutations` times, the study labels are shuffled without replacement (a uniformly random permutation of
    !| `1..n_samples`, GSL's `gsl_ran_shuffle`, one generator seeded with `random_seed`), and the whole JS-Comp-Test runs on the
    !| shuffled data. Row `i` of `jsd_null` holds its values.
    !| 3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`, so it is never 0.
    !|
    !| With `rerun_parameter_search = .true.` every data set, the observed one included, first runs the fixed-k parameter search. `random_seed` then
    !| also seeds the bootstrap of every search: each search starts from the same seed, so every data set gets the same bootstrap
    !| draws. The routine is sequential: a caller who wants parallel permutations calls
    !| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]]
    !| itself, once per assignment.
    !|
    !| Failures: a permutation that fails with invalid input (a safeguard: the relabelling keeps every study's sample count, so the k-NN test does not normally trigger it) gets a NaN row in `jsd_null`, is
    !| counted in `n_failed_permutations` and is not in `n_valid`. No valid permutation gives NaN p-values (also for
    !| `n_permutations = 0` in a direct call; the published entry points reject 0). Any other error ends the test and is
    !| returned, with NaN p-values and the rows run so far kept. An error on the observed data, or a study without any sample,
    !| is returned at once with everything NaN. In the Fortran expert tier a capacity below the search's own bound
    !| (`calc_js_comp_test_candidate_bounds_impl` for the same `grid_points_factor` and `grid_min_points`) fails the same way with the rerun option.
    !| Known limitation: a gene missing in a whole study travels with the samples.
    subroutine run_js_comp_test_permutation_c(&
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
            divergence_mode,&
            grid_points_factor,&
            grid_min_points,&
            ierr&
        ) bind(C, name="run_js_comp_test_permutation_c")
        use tox_data_integration_js_comp_test_permutation_fixed_k, only: run_js_comp_test_permutation
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        use tox_data_integration_js_comp_test_fixed_k_search_impl, only: MODE_NEIGHBORHOOD_SIZE_LARGE, MODE_NEIGHBORHOOD_SIZE_SMALL
        use tox_data_integration_js_comp_test_search_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE
        use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_permutations
            !! Number of random assignments (rows of `jsd_null`)
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence_observed
            !! Weighted global JSD of each study on the observed data (the identity assignment)
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
            !! whose value reached the observed one; NaN when none is valid or on error
        real(c_double), dimension(n_permutations, n_studies), intent(out), target :: jsd_null
            !! Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation
        integer(c_int), intent(out), target :: n_failed_permutations
            !! Number of permutations whose JS-Comp-Test failed with invalid input
        logical(c_bool), intent(in), target :: rerun_parameter_search
            !! `.true.`: run the fixed-k parameter search on the relabelled data first and use its
            !! `(n_points, n_neighbors)`; `.false.`: use `n_points` and `n_neighbors`
            !! The default value is `.false.`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                                  | Value                                                                                  |
            !! |---------------------------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass)                 | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)                  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (at least half of the studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `'join_min'`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        character(len=1, kind=c_char), dimension(19), intent(in), target :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `'plateau_effect_size'`.
        real(c_double), intent(in), target :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed of the random number generator of the label draws, and of the bootstrap of
            !! every rerun parameter search (a separate generator)
            !! The default value is `42_int32`.
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
        character(len=1, kind=c_char), dimension(23), intent(in), target :: neighborhood_size_mode
            !! Rerun search: which neighborhood size the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are not the
            !! default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode                                                  | Value                                                                                                    |
            !! |-------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Large neighborhoods, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods, `KX_FACTORS(2)` = 0.5            | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! The default value is `'neighborhood_size_large'`.
        real(c_double), intent(in), target :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0001_real64`.
        character(len=1, kind=c_char), dimension(22), intent(in), target :: divergence_mode
            !! How each study's per-point divergence is computed, in the run (observed and every
            !! shuffled data set) and in the rerun search
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `'divergence_pairwise_k2'`.
        real(c_double), intent(in), target :: grid_points_factor
            !! Rerun search (used only with `rerun_parameter_search`): factor in the k-NN candidate
            !! grid's first reference-point count `ceiling(grid_points_factor * sqrt(max_n_genes_all_studies))`
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `4.0_real64`.
        integer(c_int), intent(in), target :: grid_min_points
            !! Rerun search (used only with `rerun_parameter_search`): the k-NN candidate grid's
            !! floor, its last and coarsest `n_points` candidate
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_POINTS`.
            !! The default value is `200_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)
        integer(int32) :: join_method_mode_f
        integer(int32) :: plateau_mode_mode_f
        integer(int32) :: bin_count_mode_mode_f
        integer(int32) :: neighborhood_size_mode_mode_f
        integer(int32) :: divergence_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(n_failed_permutations)
        M_CHECK_NON_NULL(rerun_parameter_search)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(min_neighbor_overlap)
        M_CHECK_NON_NULL(succeeding_ci_overlap)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_NON_NULL(delta_abs_tolerance)
        M_CHECK_NON_NULL(grid_points_factor)
        M_CHECK_NON_NULL(grid_min_points)
        M_CHECK_ARRAY_NON_NULL(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence_observed, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)
        M_CHECK_ARRAY_NON_NULL(jsd_null, n_permutations * n_studies)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)
        M_CHECK_ARRAY_NON_NULL(plateau_mode, 19)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)
        M_CHECK_ARRAY_NON_NULL(neighborhood_size_mode, 23)
        M_CHECK_ARRAY_NON_NULL(divergence_mode, 22)

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
            character(len=:), pointer :: neighborhood_size_mode_f
            neighborhood_size_mode_f => c_char_as_view(neighborhood_size_mode)

            select case (neighborhood_size_mode_f)
                case ("neighborhood_size_large")
                    neighborhood_size_mode_mode_f = MODE_NEIGHBORHOOD_SIZE_LARGE
                case ("neighborhood_size_small")
                    neighborhood_size_mode_mode_f = MODE_NEIGHBORHOOD_SIZE_SMALL
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

        call run_js_comp_test_permutation(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            expression = expression,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            n_permutations = n_permutations,&
            global_js_divergence_observed = global_js_divergence_observed,&
            p_values = p_values,&
            jsd_null = jsd_null,&
            n_failed_permutations = n_failed_permutations,&
            rerun_parameter_search = rerun_parameter_search,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method_mode_f,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode_mode_f,&
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
            bin_count_mode = bin_count_mode_mode_f,&
            neighborhood_size_mode = neighborhood_size_mode_mode_f,&
            delta_abs_tolerance = delta_abs_tolerance,&
            divergence_mode = divergence_mode_mode_f,&
            grid_points_factor = grid_points_factor,&
            grid_min_points = grid_min_points,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_permutation_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_permutation_fixed_k(module):run_js_comp_test_permutation_expert(subroutine)]]
    !| Tests whether the assignment of samples to studies influences each study's JS-Comp-Test value (null
    !| hypothesis: it does not), with the data, parameters and options of
    !| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]].
    !|
    !| 1. The observed data (the identity assignment) give `global_js_divergence_observed`.
    !| 2. `n_permutations` times, the study labels are shuffled without replacement (a uniformly random permutation of
    !| `1..n_samples`, GSL's `gsl_ran_shuffle`, one generator seeded with `random_seed`), and the whole JS-Comp-Test runs on the
    !| shuffled data. Row `i` of `jsd_null` holds its values.
    !| 3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`, so it is never 0.
    !|
    !| With `rerun_parameter_search = .true.` every data set, the observed one included, first runs the fixed-k parameter search. `random_seed` then
    !| also seeds the bootstrap of every search: each search starts from the same seed, so every data set gets the same bootstrap
    !| draws. The routine is sequential: a caller who wants parallel permutations calls
    !| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]]
    !| itself, once per assignment.
    !|
    !| Failures: a permutation that fails with invalid input (a safeguard: the relabelling keeps every study's sample count, so the k-NN test does not normally trigger it) gets a NaN row in `jsd_null`, is
    !| counted in `n_failed_permutations` and is not in `n_valid`. No valid permutation gives NaN p-values (also for
    !| `n_permutations = 0` in a direct call; the published entry points reject 0). Any other error ends the test and is
    !| returned, with NaN p-values and the rows run so far kept. An error on the observed data, or a study without any sample,
    !| is returned at once with everything NaN. In the Fortran expert tier a capacity below the search's own bound
    !| (`calc_js_comp_test_candidate_bounds_impl` for the same `grid_points_factor` and `grid_min_points`) fails the same way with the rerun option.
    !| Known limitation: a gene missing in a whole study travels with the samples.
    subroutine run_js_comp_test_permutation_expert_c(&
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
            divergence_mode,&
            grid_points_factor,&
            grid_min_points,&
            ierr&
        ) bind(C, name="run_js_comp_test_permutation_expert_c")
        use tox_data_integration_js_comp_test_permutation_fixed_k, only: run_js_comp_test_permutation_expert
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        use tox_data_integration_js_comp_test_fixed_k_search_impl, only: MODE_NEIGHBORHOOD_SIZE_LARGE, MODE_NEIGHBORHOOD_SIZE_SMALL
        use tox_data_integration_js_comp_test_search_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE
        use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_points_candidate
            !! Capacity of the per-point work arrays
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_permutation_bounds(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `n_points`.
        integer(c_int), intent(in), target :: max_n_neighbors_candidate
            !! Capacity of the per-neighbor work arrays
            !! It is *VERY IMPORTANT* to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_permutation_bounds(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `n_neighbors`.
        integer(c_int), intent(in), target :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! It is *VERY IMPORTANT* to compute this argument from the `n_bootstrapping_top_k_jsds` output produced by [[tox_data_integration_js_comp_test_permutation_bounds(module):calc_js_comp_test_permutation_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_permutations
            !! Number of random assignments (rows of `jsd_null`)
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence_observed
            !! Weighted global JSD of each study on the observed data (the identity assignment)
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
            !! whose value reached the observed one; NaN when none is valid or on error
        real(c_double), dimension(n_permutations, n_studies), intent(out), target :: jsd_null
            !! Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation
        integer(c_int), intent(out), target :: n_failed_permutations
            !! Number of permutations whose JS-Comp-Test failed with invalid input
        integer(c_int), dimension(max_n_reps_all_studies*n_studies), intent(out), target :: tmp_sample_slots
            !! Working array: the flat slot index of every sample
        integer(c_int), dimension(max_n_reps_all_studies*n_studies), intent(out), target :: tmp_sample_assignment
            !! Working array: the current assignment, in its first `n_samples` elements
        real(c_double), dimension(n_studies), intent(out), target :: tmp_permuted_global_js_divergence
            !! Working array: one permutation's global JSD per study
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out), target :: tmp_permuted_expression
            !! Working array: the relabelled data set
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(out), target :: tmp_gene_means
            !! Working array: the relabelled data set's gene means per study
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out), target :: tmp_residuals
            !! Working array: the relabelled data set's residuals
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(out), target :: tmp_gene_means_perms
            !! Working array of the search and the run: each study's sorting permutation of its gene means
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_gene_means_perm_all
            !! Working array of the search and the run: sorting permutation of the pooled gene means
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_x_star
            !! Working array of the search and the run: the reference points
        integer(c_int), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out), target :: tmp_neighborhood_indices
            !! Working array of the search and the run: every study's neighbor gene indices
        integer(c_int), dimension(2, max_n_points_candidate, n_studies), intent(out), target :: tmp_neighborhood_range
            !! Working array of the search and the run: neighborhood spans
        real(c_double), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), intent(out), target :: tmp_neighborhood_residuals_gathered
            !! Working array of the search and the run: one study's gathered neighborhood residuals
        integer(c_int), dimension(256), intent(out), target :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(c_int), dimension(256), intent(out), target :: tmp_point_draw
            !! Working array of the search and the run: one study's draw at one point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_n_bins_per_point
            !! Working array of the search and the run: bin count per reference point
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_shared_residual_range_low
            !! Working array of the search and the run: lower residual-range bound per reference point
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_shared_residual_range_high
            !! Working array of the search and the run: upper residual-range bound per reference point
        logical(c_bool), dimension(max_n_points_candidate), intent(out), target :: tmp_occupancy_failed
            !! Working array of the search and the run: bin-count search failure flag per reference point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_n_pooled_residuals
            !! Working array of the search and the run: pooled residual count per reference point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_min_bin_occupancy
            !! Working array of the search and the run: minimum bin occupancy per reference point
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_mean_bin_occupancy
            !! Working array of the search and the run: mean bin occupancy per reference point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_max_bin_occupancy
            !! Working array of the search and the run: maximum bin occupancy per reference point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_sturges_bins
            !! Working array of the search and the run: Sturges' rule diagnostic per reference point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_fd_bins
            !! Working array of the search and the run: Freedman-Diaconis diagnostic per reference point
        real(c_double), dimension(max_n_points_candidate, 256, n_studies), intent(out), target :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(c_int), dimension(max_n_points_candidate, 256, n_studies), intent(out), target :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(c_int), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(c_double), dimension(max_n_points_candidate, 256), intent(out), target :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(c_int), dimension(max_n_points_candidate, 256), intent(out), target :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(c_double), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(c_double), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(c_double), dimension(n_studies), intent(out), target :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(c_double), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out), target :: tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(c_int), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out), target :: tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(c_int), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out), target :: tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_point_n_neighbors
            !! Working array of the search and the run: per-study neighbor counts, then sort scratch
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_search_n_bins_per_point
            !! Working array: the search's chosen bin counts
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_search_shared_residual_range_low
            !! Working array: the search's chosen lower residual-range bounds
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_search_shared_residual_range_high
            !! Working array: the search's chosen upper residual-range bounds
        real(c_double), dimension(2, n_studies), intent(out), target :: tmp_search_confidence_interval
            !! Working array: the search's chosen confidence interval
        integer(c_int), dimension(16), intent(out), target :: tmp_search_trace_n_points
            !! Working array: the search's trace
        integer(c_int), dimension(16), intent(out), target :: tmp_search_trace_n_neighbors
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(c_double), dimension(n_studies, 16), intent(out), target :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(c_double), dimension(16), intent(out), target :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(c_double), dimension(16), intent(out), target :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: tmp_search_trace_shared_residual_range_high
            !! Working array: the search's trace
        real(c_double), dimension(2, n_studies), intent(out), target :: tmp_candidate_confidence_interval
            !! Working array of the search: the current candidate's confidence interval
        real(c_double), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out), target :: tmp_bootstrapping_top_k_jsds
            !! Working array of the search: the bootstrap's bottom-k/top-k heaps
        real(c_double), dimension(n_studies), intent(out), target :: tmp_prev_global_js_divergence
            !! Working array of the search: the previous admissible candidate's global JSD
        real(c_double), dimension(2, n_studies), intent(out), target :: tmp_best_uncertainty_confidence_interval
            !! Working array of the search: the smallest-uncertainty candidate's confidence interval
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_best_n_bins_per_point
            !! Working array of the search: the best candidate's bin counts
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_n_bins_per_point
            !! Working array of the search: the smallest-uncertainty candidate's bin counts
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_shared_residual_range_low
            !! Working array of the search: the best candidate's lower residual-range bounds
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_shared_residual_range_high
            !! Working array of the search: the best candidate's upper residual-range bounds
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array of the search: the smallest-uncertainty candidate's lower bounds
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array of the search: the smallest-uncertainty candidate's upper bounds
        logical(c_bool), intent(in), target :: rerun_parameter_search
            !! `.true.`: run the fixed-k parameter search on the relabelled data first and use its
            !! `(n_points, n_neighbors)`; `.false.`: use `n_points` and `n_neighbors`
            !! The default value is `.false.`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method                                                  | Value                                                                                  |
            !! |---------------------------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass)                 | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)                  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (at least half of the studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! The default value is `'join_min'`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        character(len=1, kind=c_char), dimension(19), intent(in), target :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode                                                                                 | Value                                                                                        |
            !! |--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                                                                      | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                                                                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `'plateau_effect_size'`.
        real(c_double), intent(in), target :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
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
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed of the random number generator of the label draws, and of the bootstrap of
            !! every rerun parameter search (a separate generator)
            !! The default value is `42_int32`.
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
        character(len=1, kind=c_char), dimension(23), intent(in), target :: neighborhood_size_mode
            !! Rerun search: which neighborhood size the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are not the
            !! default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode                                                  | Value                                                                                                    |
            !! |-------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Large neighborhoods, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods, `KX_FACTORS(2)` = 0.5            | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! The default value is `'neighborhood_size_large'`.
        real(c_double), intent(in), target :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.0001_real64`.
        character(len=1, kind=c_char), dimension(22), intent(in), target :: divergence_mode
            !! How each study's per-point divergence is computed, in the run (observed and every
            !! shuffled data set) and in the rerun search
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `'divergence_pairwise_k2'`.
        real(c_double), intent(in), target :: grid_points_factor
            !! Rerun search (used only with `rerun_parameter_search`): factor in the k-NN candidate
            !! grid's first reference-point count `ceiling(grid_points_factor * sqrt(max_n_genes_all_studies))`
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `4.0_real64`.
        integer(c_int), intent(in), target :: grid_min_points
            !! Rerun search (used only with `rerun_parameter_search`): the k-NN candidate grid's
            !! floor, its last and coarsest `n_points` candidate
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_POINTS`.
            !! The default value is `200_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)
        integer(int32) :: join_method_mode_f
        integer(int32) :: plateau_mode_mode_f
        integer(int32) :: bin_count_mode_mode_f
        integer(int32) :: neighborhood_size_mode_mode_f
        integer(int32) :: divergence_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(max_n_neighbors_candidate)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(n_failed_permutations)
        M_CHECK_NON_NULL(rerun_parameter_search)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(min_neighbor_overlap)
        M_CHECK_NON_NULL(succeeding_ci_overlap)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_NON_NULL(delta_abs_tolerance)
        M_CHECK_NON_NULL(grid_points_factor)
        M_CHECK_NON_NULL(grid_min_points)
        M_CHECK_ARRAY_NON_NULL(expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence_observed, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)
        M_CHECK_ARRAY_NON_NULL(jsd_null, n_permutations * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_sample_slots, (max_n_reps_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_sample_assignment, (max_n_reps_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_permuted_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permuted_expression, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perms, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perm_all, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_x_star, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_indices, max_n_neighbors_candidate * max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_range, 2 * max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_residuals_gathered, max_n_reps_all_studies * max_n_neighbors_candidate * max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_point_pool, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_point_draw, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_occupancy_failed, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_n_pooled_residuals, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_min_bin_occupancy, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_bin_occupancy, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_max_bin_occupancy, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_sturges_bins, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_fd_bins, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, max_n_points_candidate * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, max_n_points_candidate * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_included_n_reps, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf, max_n_points_candidate * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf_counts, max_n_points_candidate * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf_included_n_reps, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals, (max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals_perm, (max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_residual_to_study_perm, (max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_bin_counts_search, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_ends, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_sizes, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_point_n_neighbors, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_search_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_search_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_search_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_search_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_n_points, 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_n_neighbors, 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_global_js_divergence, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_ci_lower, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_ci_upper, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_ci_width, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_ci_width_relative, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_ci_overlap, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_delta, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_delta_median, 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_delta_max, 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_selected_n_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_occupancy_failed, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_n_pooled_residuals, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_min_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_mean_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_max_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_sturges_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_fd_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_shared_residual_range_low, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_search_trace_shared_residual_range_high, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_candidate_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds * 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_prev_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_best_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)
        M_CHECK_ARRAY_NON_NULL(plateau_mode, 19)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)
        M_CHECK_ARRAY_NON_NULL(neighborhood_size_mode, 23)
        M_CHECK_ARRAY_NON_NULL(divergence_mode, 22)

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
            character(len=:), pointer :: neighborhood_size_mode_f
            neighborhood_size_mode_f => c_char_as_view(neighborhood_size_mode)

            select case (neighborhood_size_mode_f)
                case ("neighborhood_size_large")
                    neighborhood_size_mode_mode_f = MODE_NEIGHBORHOOD_SIZE_LARGE
                case ("neighborhood_size_small")
                    neighborhood_size_mode_mode_f = MODE_NEIGHBORHOOD_SIZE_SMALL
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

        call run_js_comp_test_permutation_expert(&
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
            join_method = join_method_mode_f,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode_mode_f,&
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
            bin_count_mode = bin_count_mode_mode_f,&
            neighborhood_size_mode = neighborhood_size_mode_mode_f,&
            delta_abs_tolerance = delta_abs_tolerance,&
            divergence_mode = divergence_mode_mode_f,&
            grid_points_factor = grid_points_factor,&
            grid_min_points = grid_min_points,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_permutation_expert_c

end module tox_data_integration_js_comp_test_permutation_fixed_k_c
#endif
