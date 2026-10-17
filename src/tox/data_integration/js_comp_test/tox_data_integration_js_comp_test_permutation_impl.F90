#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

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
module tox_data_integration_js_comp_test_permutation_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use tox_errors, only: set_ok, set_err_once, is_err, get_err_code, validate_in_range_int, validate_in_range_real, &
                          ERR_INVALID_INPUT
    use f42_random_gsl, only: rng_t, create_rng, destroy_rng, random_uniform
    use f42_sort_impl, only: init_perm
    use tox_data_integration_preprocessing_impl, only: compute_gene_means_impl, compute_residuals_impl, pool_means_impl, &
                                                       sort_gene_means_per_study, sort_pooled_gene_means, &
                                                       construct_adaptive_neighborhoods_impl
    use tox_data_integration_stats_impl, only: calc_js_comp_test_n_top_k_jsds
    use tox_data_integration_js_comp_test_binning_impl, only: MAX_N_BINS
    use tox_data_integration_js_comp_test_fixed_k_search_impl, only: calc_js_comp_test_candidate_bounds, &
                                                                     run_js_comp_test_parameter_search_impl
    use tox_data_integration_js_comp_test_adaptive_search_impl, only: calc_adaptive_js_comp_test_bounds, &
                                                                      run_js_comp_test_adaptive_parameter_search_impl, &
                                                                      ADAPTIVE_POINT_CAPACITY_FACTOR
    use tox_data_integration_js_comp_test_run_impl, only: run_js_comp_test_impl, run_js_comp_test_adaptive_impl
    M_IMPLICIT_NONE
    private
    public :: calc_js_comp_test_permutation_bounds, calc_adaptive_js_comp_test_permutation_bounds, &
             run_js_comp_test_on_sample_assignment_impl, run_js_comp_test_adaptive_on_sample_assignment_impl, &
             run_js_comp_test_label_permutation_test_impl, run_js_comp_test_adaptive_label_permutation_test_impl

    ! `n_bootstraps`' default here (mandatory in both parameter searches), shared by DM_DEFAULT and
    ! M_DEFAULT_VAL of every routine below and of both sizing routines.
#define CM_N_BOOTSTRAPS_DEFAULT 1000_int32
    ! `rerun_parameter_search`'s default, shared by DM_DEFAULT and M_DEFAULT_VAL.
#define CM_RERUN_PARAMETER_SEARCH_DEFAULT .false.

contains

    !> M_EXPORT_C
    !| summary: Recommend the work-array capacities of the k-NN study-label Permutation-Test
    !| AUTHOR_LASZLO_LANG
    !| Sizes the per-point and per-neighbor work arrays of
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]]
    !| and
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_label_permutation_test_impl(interface)]]
    !| so that both the fixed `(n_points, n_neighbors)` run and a rerun of the parameter search fit:
    !| `max_n_points_candidate = max(n_points, P)` and `max_n_neighbors_candidate = max(n_neighbors,
    !| K)`, with `P` and `K` the bounds of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds(subroutine)]].
    !| `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
    !| [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds(subroutine)]] computes it.
    !|
    !| Rejects an out-of-range argument with invalid input (a non-finite significance level with
    !| `ERR_NAN_INF`); the outputs are then left undefined.
    pure subroutine calc_js_comp_test_permutation_bounds(max_n_genes_all_studies, n_points, n_neighbors, n_bootstraps, &
                                                         two_sided_bootstrapping_significance_level, max_n_points_candidate, &
                                                         max_n_neighbors_candidate, n_bootstrapping_top_k_jsds, ierr)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_points
            !! Number of reference points of the fixed-parameter run
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood of the fixed-parameter run
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies)
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps of the parameter search, when it is rerun
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_N_BOOTSTRAPS_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(out) :: max_n_points_candidate
            !! Capacity for reference points
        integer(int32), intent(out) :: max_n_neighbors_candidate
            !! Capacity for neighbors per neighborhood
        integer(int32), intent(out) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size (values kept at each end of the bootstrap distribution)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure

        integer(int32) :: actual_n_bootstraps, search_n_points, search_n_neighbors

        call set_ok(ierr)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=3_int32, min=1_int32, max=max(1_int32, max_n_genes_all_studies))
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=5_int32, &
                                    min=0.0_real64, max=100.0_real64)
        if (is_err(ierr)) return
        M_DEFAULT_VAL(n_bootstraps, actual_n_bootstraps, CM_N_BOOTSTRAPS_DEFAULT)

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, search_n_points, search_n_neighbors, ierr)
        if (is_err(ierr)) return
        max_n_points_candidate = max(n_points, search_n_points)
        max_n_neighbors_candidate = max(n_neighbors, search_n_neighbors)

        call calc_js_comp_test_n_top_k_jsds(actual_n_bootstraps, two_sided_bootstrapping_significance_level, &
                                            n_bootstrapping_top_k_jsds, ierr)
    end subroutine calc_js_comp_test_permutation_bounds

    !> M_EXPORT_C
    !| summary: Recommend the work-array capacity of the adaptive study-label Permutation-Test
    !| AUTHOR_LASZLO_LANG
    !| Sizes the per-point work arrays of
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]
    !| and
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_label_permutation_test_impl(interface)]]:
    !| `max_n_points_candidate = max(A, min(N, ceiling(F * N / k_start) + 1))`, with `N =
    !| max_n_genes_all_studies * n_studies` (the padded pool size), `F` the factor
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_POINT_CAPACITY_FACTOR(variable)]]
    !| and `A` the capacity of the adaptive parameter search,
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):calc_adaptive_js_comp_test_bounds(subroutine)]].
    !| As there, this is a **practical capacity, not a proven bound**: how many reference points a
    !| construction builds is only known once it has run, and a shuffled data set whose
    !| construction builds more fails with invalid input (it is then counted as a failed
    !| permutation). `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
    !| [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds(subroutine)]] computes it.
    !|
    !| Rejects an out-of-range argument with invalid input (a non-finite significance level with
    !| `ERR_NAN_INF`), and, naming neither argument, the joint bounds on `max_n_genes_all_studies`
    !| and `n_studies` that `calc_adaptive_js_comp_test_bounds` rejects; the outputs are then left
    !| undefined.
    pure subroutine calc_adaptive_js_comp_test_permutation_bounds(max_n_genes_all_studies, n_studies, k_start, n_bootstraps, &
                                                                  two_sided_bootstrapping_significance_level, &
                                                                  max_n_points_candidate, n_bootstrapping_top_k_jsds, ierr)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood of the fixed-knob run takes unconditionally
            !! DM_MIN(1_int32)
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps of the parameter search, when it is rerun
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_N_BOOTSTRAPS_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(out) :: max_n_points_candidate
            !! Practical capacity for reference points
        integer(int32), intent(out) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size (values kept at each end of the bootstrap distribution)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure

        integer(int32) :: actual_n_bootstraps, search_capacity
        integer(int64) :: n_pool, capacity

        call set_ok(ierr)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(k_start, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=5_int32, &
                                    min=0.0_real64, max=100.0_real64)
        if (is_err(ierr)) return
        M_DEFAULT_VAL(n_bootstraps, actual_n_bootstraps, CM_N_BOOTSTRAPS_DEFAULT)

        ! Its argument positions 1 and 2 are this routine's too; its joint bounds name none.
        call calc_adaptive_js_comp_test_bounds(max_n_genes_all_studies, n_studies, search_capacity, ierr)
        if (is_err(ierr)) return

        n_pool = int(max_n_genes_all_studies, int64)*int(n_studies, int64)
        capacity = ceiling(ADAPTIVE_POINT_CAPACITY_FACTOR*real(n_pool, real64)/real(k_start, real64), kind=int64) + 1_int64
        max_n_points_candidate = max(search_capacity, int(min(n_pool, capacity), kind=int32))

        call calc_js_comp_test_n_top_k_jsds(actual_n_bootstraps, two_sided_bootstrapping_significance_level, &
                                            n_bootstrapping_top_k_jsds, ierr)
    end subroutine calc_adaptive_js_comp_test_permutation_bounds

    !> summary: Find the sample slots of `expression`: every replicate slot that is not all NaN
    !| AUTHOR_LASZLO_LANG
    !| Private. Sample `k` is the `k`-th slot that has at least one non-NaN value, in slot order
    !| (replicate fastest, then study); `sample_slots(k)` is its flat slot index
    !| `(i_study - 1)*max_n_reps_all_studies + i_rep`. Sequential: the sample numbering is the order.
    pure subroutine locate_samples(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, sample_slots, &
                                   n_samples_found, has_empty_study)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), intent(in) :: expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Expression of every replicate slot of every gene in every study
        integer(int32), intent(out) :: sample_slots(max_n_reps_all_studies*n_studies)
            !! Flat slot index of every sample, in sample order; 0 beyond `n_samples_found`
        integer(int32), intent(out) :: n_samples_found
            !! Number of samples (slots that are not all NaN)
        logical(c_bool), intent(out) :: has_empty_study
            !! `.true.` iff some study has no sample at all

        integer(int32) :: i_study, i_rep, i_gene, n_in_study
        logical(c_bool) :: is_sample

        sample_slots = 0_int32
        n_samples_found = 0_int32
        has_empty_study = .false.
        do i_study = 1, n_studies
            n_in_study = 0_int32
            do i_rep = 1, max_n_reps_all_studies
                is_sample = .false.
                do i_gene = 1, max_n_genes_all_studies
                    if (.not. ieee_is_nan(expression(i_rep, i_gene, i_study))) then
                        is_sample = .true.
                        exit
                    end if
                end do
                if (is_sample) then
                    n_in_study = n_in_study + 1_int32
                    n_samples_found = n_samples_found + 1_int32
                    sample_slots(n_samples_found) = (i_study - 1_int32)*max_n_reps_all_studies + i_rep
                end if
            end do
            if (n_in_study == 0_int32) has_empty_study = .true.
        end do
    end subroutine locate_samples

    !> summary: Check that `sample_assignment` holds every value of `1..n_samples` exactly once
    !| AUTHOR_LASZLO_LANG
    !| Private. Sequential: `sample_mark` records the values seen so far.
    pure subroutine check_sample_assignment(n_samples, sample_assignment, sample_mark, is_permutation)
        integer(int32), intent(in) :: n_samples
            !! Number of samples
        integer(int32), intent(in) :: sample_assignment(n_samples)
            !! The assignment to check
        logical(c_bool), intent(out) :: sample_mark(n_samples)
            !! Working array: `.true.` for every value already seen
        logical(c_bool), intent(out) :: is_permutation
            !! `.true.` iff every value is in `1..n_samples` and none repeats

        integer(int32) :: i_sample, source

        sample_mark = .false.
        is_permutation = .true.
        do i_sample = 1, n_samples
            source = sample_assignment(i_sample)
            if (source < 1_int32 .or. source > n_samples) then
                is_permutation = .false.
                return
            end if
            if (sample_mark(source)) then
                is_permutation = .false.
                return
            end if
            sample_mark(source) = .true.
        end do
    end subroutine check_sample_assignment

    !> summary: Build the relabelled data set and every study's gene means and residuals
    !| AUTHOR_LASZLO_LANG
    !| Private. Destination sample `k` receives all genes of source sample `sample_assignment(k)`;
    !| padding slots stay as they are (all NaN). Then, per study, exactly what a caller does per study
    !| today: [[tox_data_integration_preprocessing_impl(module):compute_gene_means_impl(interface)]] and
    !| [[tox_data_integration_preprocessing_impl(module):compute_residuals_impl(interface)]] on that
    !| study's slice (NaN padding is skipped by both).
    pure subroutine assign_samples_and_center(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                              n_samples, sample_assignment, sample_slots, permuted_expression, gene_means, &
                                              residuals)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), intent(in) :: expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Expression of every replicate slot of every gene in every study
        integer(int32), intent(in) :: n_samples
            !! Number of samples
        integer(int32), intent(in) :: sample_assignment(n_samples)
            !! A permutation of `1..n_samples`, already checked
        integer(int32), intent(in) :: sample_slots(n_samples)
            !! Flat slot index of every sample, from locate_samples
        real(real64), intent(out) :: permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! The relabelled data set
        real(real64), intent(out) :: gene_means(max_n_genes_all_studies, n_studies)
            !! Per-gene mean of every study of the relabelled data set; NaN where a study has no value
        real(real64), intent(out) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Signed residuals of the relabelled data set; NaN where it has no value

        integer(int32) :: i_sample, i_gene, i_study, destination, source

        permuted_expression = expression
        ! Independent: the destination slots are distinct, and every source is read from `expression`.
        do concurrent(i_sample=1:n_samples, i_gene=1:max_n_genes_all_studies) local(destination, source) &
            shared(sample_slots, sample_assignment, permuted_expression, expression, max_n_reps_all_studies)
            destination = sample_slots(i_sample) - 1_int32
            source = sample_slots(sample_assignment(i_sample)) - 1_int32
            permuted_expression(mod(destination, max_n_reps_all_studies) + 1_int32, i_gene, &
                                destination/max_n_reps_all_studies + 1_int32) = &
                expression(mod(source, max_n_reps_all_studies) + 1_int32, i_gene, source/max_n_reps_all_studies + 1_int32)
        end do

        do concurrent(i_study=1:n_studies) shared(permuted_expression, gene_means, residuals, max_n_genes_all_studies, &
                                                  max_n_reps_all_studies)
            call compute_gene_means_impl(max_n_genes_all_studies, max_n_reps_all_studies, &
                                         permuted_expression(:, :, i_study), gene_means(:, i_study))
            call compute_residuals_impl(max_n_genes_all_studies, max_n_reps_all_studies, permuted_expression(:, :, i_study), &
                                        gene_means(:, i_study), residuals(:, :, i_study))
        end do
    end subroutine assign_samples_and_center

    !> summary: Check a rerun's capacities against the fixed-k parameter search's own bounds
    !| AUTHOR_LASZLO_LANG
    !| Private. `points_fit`/`neighbors_fit` are `.true.` iff the capacity reaches the bound of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds(subroutine)]].
    pure subroutine check_fixed_k_rerun_capacities(max_n_genes_all_studies, max_n_points_candidate, &
                                                   max_n_neighbors_candidate, points_fit, neighbors_fit)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies, at least 1
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity for reference points
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Capacity for neighbors per neighborhood
        logical(c_bool), intent(out) :: points_fit
            !! `.true.` iff `max_n_points_candidate` reaches the search's bound
        logical(c_bool), intent(out) :: neighbors_fit
            !! `.true.` iff `max_n_neighbors_candidate` reaches the search's bound

        integer(int32) :: search_n_points, search_n_neighbors, bounds_ierr

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, search_n_points, search_n_neighbors, bounds_ierr)
        points_fit = (.not. is_err(bounds_ierr)) .and. max_n_points_candidate >= search_n_points
        neighbors_fit = (.not. is_err(bounds_ierr)) .and. max_n_neighbors_candidate >= search_n_neighbors
    end subroutine check_fixed_k_rerun_capacities

    !> summary: Run the k-NN JS-Comp-Test on one assignment of samples to studies
    !| AUTHOR_LASZLO_LANG
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
    subroutine run_js_comp_test_on_sample_assignment_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                          expression, n_samples, sample_assignment, n_points, n_neighbors, &
                                                          max_n_points_candidate, max_n_neighbors_candidate, &
                                                          n_bootstrapping_top_k_jsds, global_js_divergence, used_n_points, &
                                                          used_n_neighbors, tmp_sample_slots, tmp_sample_mark, &
                                                          tmp_permuted_expression, tmp_gene_means, tmp_residuals, &
                                                          tmp_gene_means_perms, tmp_gene_means_perm_all, tmp_x_star, &
                                                          tmp_neighborhood_indices, tmp_neighborhood_range, &
                                                          tmp_neighborhood_residuals_gathered, tmp_point_pool, tmp_point_draw, &
                                                          tmp_n_bins_per_point, tmp_shared_residual_range_low, &
                                                          tmp_shared_residual_range_high, tmp_occupancy_failed, &
                                                          tmp_n_pooled_residuals, tmp_min_bin_occupancy, &
                                                          tmp_mean_bin_occupancy, tmp_max_bin_occupancy, tmp_sturges_bins, &
                                                          tmp_fd_bins, tmp_pmfs, tmp_counts, tmp_included_n_reps, tmp_mean_pmf, &
                                                          tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps, &
                                                          tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                                          tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                                          tmp_residual_to_study_perm, tmp_bin_counts_search, &
                                                          tmp_study_block_ends, tmp_study_block_sizes, tmp_point_n_neighbors, &
                                                          tmp_p_values, tmp_p_values_observed_consensus, &
                                                          tmp_permutation_mean_pmf, tmp_permutation_mean_pmf_counts, &
                                                          tmp_permutation_counts, tmp_permutation_pmfs, &
                                                          tmp_permutation_js_divergences, tmp_permutation_weights, &
                                                          tmp_search_n_bins_per_point, tmp_search_shared_residual_range_low, &
                                                          tmp_search_shared_residual_range_high, &
                                                          tmp_search_confidence_interval, tmp_search_trace_n_points, &
                                                          tmp_search_trace_n_neighbors, tmp_search_trace_global_js_divergence, &
                                                          tmp_search_trace_ci_lower, tmp_search_trace_ci_upper, &
                                                          tmp_search_trace_ci_width, tmp_search_trace_ci_width_relative, &
                                                          tmp_search_trace_ci_overlap, tmp_search_trace_delta, &
                                                          tmp_search_trace_delta_median, tmp_search_trace_delta_max, &
                                                          tmp_search_trace_selected_n_bins, tmp_search_trace_occupancy_failed, &
                                                          tmp_search_trace_n_pooled_residuals, &
                                                          tmp_search_trace_min_bin_occupancy, &
                                                          tmp_search_trace_mean_bin_occupancy, &
                                                          tmp_search_trace_max_bin_occupancy, tmp_search_trace_sturges_bins, &
                                                          tmp_search_trace_fd_bins, tmp_search_trace_shared_residual_range_low, &
                                                          tmp_search_trace_shared_residual_range_high, &
                                                          tmp_candidate_confidence_interval, tmp_bootstrapping_top_k_jsds, &
                                                          tmp_prev_global_js_divergence, &
                                                          tmp_best_uncertainty_confidence_interval, tmp_best_n_bins_per_point, &
                                                          tmp_best_uncertainty_n_bins_per_point, &
                                                          tmp_best_shared_residual_range_low, &
                                                          tmp_best_shared_residual_range_high, &
                                                          tmp_best_uncertainty_shared_residual_range_low, &
                                                          tmp_best_uncertainty_shared_residual_range_high, &
                                                          rerun_parameter_search, n_bootstraps, join_method, &
                                                          min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, &
                                                          plateau_mode, delta_median_threshold, delta_max_threshold, &
                                                          delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, &
                                                          gamma_occupancy, lower_residual_range_quantile, &
                                                          upper_residual_range_quantile, &
                                                          two_sided_bootstrapping_significance_level, random_seed, &
                                                          bin_count_mode, &
                                                          neighborhood_size_mode, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! DM_MIN(1_int32)
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: n_samples
            !! Number of samples, the non-padding slots of `expression`
            !! DM_MIN(1_int32)
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
            !! A permutation of `1..n_samples`: destination sample `k` receives source sample
            !! `sample_assignment(k)`; the identity is the observed data
            !! DM_MIN(1_int32)
            !! DM_MAX(n_samples)
        integer(int32), intent(in) :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies)
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_impl, AUTO)
            !! DM_MIN(n_points)
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Capacity of the per-neighbor work arrays
            !! DM_OUTPUT_FROM(max_n_neighbors_candidate, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_impl, AUTO)
            !! DM_MIN(n_neighbors)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! DM_OUTPUT_FROM(n_bootstrapping_top_k_jsds, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_impl, AUTO)
            !! DM_MIN(1_int32)
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
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: &
            tmp_permuted_expression
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
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: &
            tmp_neighborhood_indices
            !! Working array of the search and the run: every study's neighbor gene indices
        integer(int32), dimension(2, max_n_points_candidate, n_studies), intent(out) :: tmp_neighborhood_range
            !! Working array of the search and the run: neighborhood spans
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), &
            intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array of the search and the run: one study's gathered neighborhood residuals
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_draw
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
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
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
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_counts
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_permutation_pmfs
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
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_neighbors
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_high
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
            !! DM_DEFAULT(CM_RERUN_PARAMETER_SEARCH_DEFAULT)
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_N_BOOTSTRAPS_DEFAULT)
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! DM_DEFAULT(M_METHOD_JOIN_MIN)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_MIN_NEIGHBOR_OVERLAP_DEFAULT)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_SUCCEEDING_CI_OVERLAP_DEFAULT)
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(M_MODE_PLATEAU_EFFECT_SIZE)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
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
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(in), optional :: random_seed
            !! Rerun search: seed of the bootstrap's random number generator
            !! DM_DEFAULT(M_RANDOM_SEED_DEFAULT)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
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
        integer(int32), intent(in), optional :: neighborhood_size_mode
            !! Rerun search: which neighborhood sizes the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are left out
            !! by default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Large neighborhoods only, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods only, `KX_FACTORS(2)` = 0.5 | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! | Both, large then small for each `n_points` | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_BOTH(variable)]] |
            !! DM_DEFAULT(M_MODE_NEIGHBORHOOD_SIZE_LARGE)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        integer(int32) :: n_samples_found, n_pool, actual_n_bootstraps, actual_join_method, callee_ierr
        integer(int32) :: max_n_bins_per_point, search_n_admissible_evaluated
        logical(c_bool) :: actual_rerun, has_empty_study, is_permutation, points_fit, neighbors_fit
        logical(c_bool) :: search_plateau_established

        call set_ok(ierr)
        M_DEFAULT_VAL(rerun_parameter_search, actual_rerun, CM_RERUN_PARAMETER_SEARCH_DEFAULT)
        M_DEFAULT_VAL(n_bootstraps, actual_n_bootstraps, CM_N_BOOTSTRAPS_DEFAULT)
        M_DEFAULT_VAL(join_method, actual_join_method, M_METHOD_JOIN_MIN)
        global_js_divergence = M_NAN
        used_n_points = n_points
        used_n_neighbors = n_neighbors

        ! Runtime checks (codegen guide Sec 5.14), each naming the argument at fault.
        call locate_samples(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, tmp_sample_slots, &
                            n_samples_found, has_empty_study)
        if (has_empty_study) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
            return
        end if
        if (n_samples_found /= n_samples) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
            return
        end if
        call check_sample_assignment(n_samples, sample_assignment, tmp_sample_mark, is_permutation)
        if (.not. is_permutation) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=6_int32)
            return
        end if
        if (actual_rerun) then
            call check_fixed_k_rerun_capacities(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate, &
                                                points_fit, neighbors_fit)
            if (.not. points_fit) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=9_int32)
            if (.not. neighbors_fit) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=10_int32)
            if (is_err(ierr)) return
        end if

        call assign_samples_and_center(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, n_samples, &
                                       sample_assignment, tmp_sample_slots, tmp_permuted_expression, tmp_gene_means, &
                                       tmp_residuals)

        ! Every capacity-sized buffer is passed whole to explicit-shape dummies declared at the
        ! sizes actually used (sequence association, no copy), as the fixed-k search does.
        if (actual_rerun) then
            call run_js_comp_test_parameter_search_impl( &
                n_studies=n_studies, max_n_genes_all_studies=max_n_genes_all_studies, &
                max_n_reps_all_studies=max_n_reps_all_studies, gene_means=tmp_gene_means, residuals=tmp_residuals, &
                n_bootstraps=actual_n_bootstraps, join_method=actual_join_method, &
                max_n_points_candidate=max_n_points_candidate, max_n_neighbors_candidate=max_n_neighbors_candidate, &
                n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, n_points=used_n_points, &
                n_neighbors=used_n_neighbors, n_bins_per_point=tmp_search_n_bins_per_point, &
                shared_residual_range_low=tmp_search_shared_residual_range_low, &
                shared_residual_range_high=tmp_search_shared_residual_range_high, &
                best_candidate_pair_confidence_interval=tmp_search_confidence_interval, &
                plateau_established=search_plateau_established, n_admissible_evaluated=search_n_admissible_evaluated, &
                trace_n_points=tmp_search_trace_n_points, trace_n_neighbors=tmp_search_trace_n_neighbors, &
                trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
                trace_ci_lower=tmp_search_trace_ci_lower, trace_ci_upper=tmp_search_trace_ci_upper, &
                trace_ci_width=tmp_search_trace_ci_width, trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
                trace_ci_overlap=tmp_search_trace_ci_overlap, trace_delta=tmp_search_trace_delta, &
                trace_delta_median=tmp_search_trace_delta_median, trace_delta_max=tmp_search_trace_delta_max, &
                trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
                trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
                trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
                trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
                trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
                trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
                trace_sturges_bins=tmp_search_trace_sturges_bins, trace_fd_bins=tmp_search_trace_fd_bins, &
                trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
                trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
                tmp_gene_means_perms=tmp_gene_means_perms, tmp_gene_means_perm_all=tmp_gene_means_perm_all, &
                tmp_x_star=tmp_x_star, tmp_neighborhood_indices_all_studies=tmp_neighborhood_indices, &
                tmp_neighborhood_range=tmp_neighborhood_range, &
                tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, tmp_point_pool=tmp_point_pool, &
                tmp_point_draw=tmp_point_draw, tmp_n_bins_per_point=tmp_n_bins_per_point, &
                tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
                tmp_shared_residual_range_high=tmp_shared_residual_range_high, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
                tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
                tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
                tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
                tmp_global_js_divergence=tmp_global_js_divergence, &
                tmp_confidence_interval=tmp_candidate_confidence_interval, &
                tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
                tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, tmp_delta_perm=tmp_point_n_neighbors, &
                tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
                tmp_pooled_residuals=tmp_pooled_residuals, tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
                tmp_bin_counts_search=tmp_bin_counts_search, tmp_residual_to_study_perm=tmp_residual_to_study_perm, &
                tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
                tmp_occupancy_failed=tmp_occupancy_failed, tmp_n_pooled_residuals=tmp_n_pooled_residuals, &
                tmp_min_bin_occupancy=tmp_min_bin_occupancy, tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, &
                tmp_max_bin_occupancy=tmp_max_bin_occupancy, tmp_sturges_bins=tmp_sturges_bins, tmp_fd_bins=tmp_fd_bins, &
                tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
                tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
                tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
                tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
                tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
                tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
                min_residuals_per_bin=min_residuals_per_bin, min_neighbor_overlap=min_neighbor_overlap, &
                succeeding_ci_overlap=succeeding_ci_overlap, plateau_mode=plateau_mode, &
                delta_median_threshold=delta_median_threshold, delta_max_threshold=delta_max_threshold, &
                delta_epsilon=delta_epsilon, delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
                m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
                lower_residual_range_quantile=lower_residual_range_quantile, &
                upper_residual_range_quantile=upper_residual_range_quantile, &
                two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
                random_seed=random_seed, bin_count_mode=bin_count_mode, &
                neighborhood_size_mode=neighborhood_size_mode, ierr=callee_ierr)
            if (is_err(callee_ierr)) then
                call set_err_once(ierr, get_err_code(callee_ierr))
                return
            end if
        end if

        ! The reference points as the fixed-k search computes them, each study's order as the plain
        ! run_js_comp_test's prologue builds it.
        call sort_pooled_gene_means(tmp_gene_means, max_n_genes_all_studies, n_studies, tmp_gene_means_perm_all)
        call pool_means_impl(tmp_gene_means, tmp_gene_means_perm_all, max_n_genes_all_studies*n_studies, used_n_points, &
                             n_pool, tmp_x_star)
        call sort_gene_means_per_study(tmp_gene_means, max_n_genes_all_studies, n_studies, tmp_gene_means_perms)

        call run_js_comp_test_impl( &
            n_studies=n_studies, max_n_genes_all_studies=max_n_genes_all_studies, &
            max_n_reps_all_studies=max_n_reps_all_studies, n_points=used_n_points, n_neighbors=used_n_neighbors, &
            gene_means=tmp_gene_means, gene_means_perms=tmp_gene_means_perms, residuals=tmp_residuals, x_star=tmp_x_star, &
            neighborhood_indices=tmp_neighborhood_indices, neighborhood_range=tmp_neighborhood_range, &
            n_bins_per_point=tmp_n_bins_per_point, shared_residual_range_low=tmp_shared_residual_range_low, &
            shared_residual_range_high=tmp_shared_residual_range_high, max_n_bins_per_point=max_n_bins_per_point, &
            occupancy_failed=tmp_occupancy_failed, n_pooled_residuals=tmp_n_pooled_residuals, &
            min_bin_occupancy=tmp_min_bin_occupancy, mean_bin_occupancy=tmp_mean_bin_occupancy, &
            max_bin_occupancy=tmp_max_bin_occupancy, sturges_bins=tmp_sturges_bins, fd_bins=tmp_fd_bins, pmfs=tmp_pmfs, &
            counts=tmp_counts, included_n_reps=tmp_included_n_reps, mean_pmf=tmp_mean_pmf, &
            mean_pmf_counts=tmp_mean_pmf_counts, mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            js_divergences=tmp_js_divergences, weights=tmp_weights, global_js_divergence=global_js_divergence, &
            p_values=tmp_p_values, p_values_observed_consensus=tmp_p_values_observed_consensus, &
            tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, &
            tmp_permutation_mean_pmf=tmp_permutation_mean_pmf, tmp_permutation_point_pool=tmp_point_pool, &
            tmp_permutation_point_draw=tmp_point_draw, tmp_pooled_residuals=tmp_pooled_residuals, &
            tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
            tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_study_block_ends=tmp_study_block_ends, &
            tmp_study_block_sizes=tmp_study_block_sizes, tmp_point_n_neighbors=tmp_point_n_neighbors, &
            tmp_permutation_mean_pmf_counts=tmp_permutation_mean_pmf_counts, &
            tmp_permutation_counts=tmp_permutation_counts, tmp_permutation_pmfs=tmp_permutation_pmfs, &
            tmp_permutation_js_divergences=tmp_permutation_js_divergences, &
            tmp_permutation_weights=tmp_permutation_weights, &
            tmp_permutation_global_js_divergence=tmp_global_js_divergence, n_permutations=0_int32, &
            min_residuals_per_bin=min_residuals_per_bin, m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
            lower_residual_range_quantile=lower_residual_range_quantile, &
            upper_residual_range_quantile=upper_residual_range_quantile, bin_count_mode=bin_count_mode, ierr=callee_ierr)
        if (is_err(callee_ierr)) then
            call set_err_once(ierr, get_err_code(callee_ierr))
            global_js_divergence = M_NAN
            return
        end if
    end subroutine run_js_comp_test_on_sample_assignment_impl

    !> summary: Run the adaptive JS-Comp-Test on one assignment of samples to studies
    !| AUTHOR_LASZLO_LANG
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
    subroutine run_js_comp_test_adaptive_on_sample_assignment_impl(n_studies, max_n_genes_all_studies, &
                                                                   max_n_reps_all_studies, expression, n_samples, &
                                                                   sample_assignment, k_start, k_step, k_max, &
                                                                   max_n_points_candidate, n_bootstrapping_top_k_jsds, &
                                                                   global_js_divergence, used_k_start, used_k_step, &
                                                                   used_k_max, n_points, tmp_sample_slots, tmp_sample_mark, &
                                                                   tmp_permuted_expression, tmp_gene_means, tmp_residuals, &
                                                                   tmp_gene_means_perm_all, tmp_x_star, &
                                                                   tmp_pooled_neighborhood_range, tmp_n_neighbors_per_point, &
                                                                   tmp_stop_reason, tmp_neighborhood_dispersion, &
                                                                   tmp_neighborhood_mad, tmp_point_neighborhood_indices, &
                                                                   tmp_neighbor_residuals, tmp_point_pool, tmp_point_draw, &
                                                                   tmp_n_bins_per_point, tmp_shared_residual_range_low, &
                                                                   tmp_shared_residual_range_high, tmp_occupancy_failed, &
                                                                   tmp_n_pooled_residuals, tmp_min_bin_occupancy, &
                                                                   tmp_mean_bin_occupancy, tmp_max_bin_occupancy, &
                                                                   tmp_sturges_bins, tmp_fd_bins, tmp_pmfs, tmp_counts, &
                                                                   tmp_included_n_reps, tmp_mean_pmf, tmp_mean_pmf_counts, &
                                                                   tmp_mean_pmf_included_n_reps, tmp_js_divergences, &
                                                                   tmp_weights, tmp_global_js_divergence, &
                                                                   tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                                                   tmp_residual_to_study_perm, tmp_bin_counts_search, &
                                                                   tmp_study_block_ends, tmp_study_block_sizes, &
                                                                   tmp_p_values, tmp_p_values_observed_consensus, &
                                                                   tmp_permutation_mean_pmf, tmp_permutation_mean_pmf_counts, &
                                                                   tmp_permutation_counts, tmp_permutation_pmfs, &
                                                                   tmp_permutation_js_divergences, tmp_permutation_weights, &
                                                                   tmp_search_n_bins_per_point, &
                                                                   tmp_search_shared_residual_range_low, &
                                                                   tmp_search_shared_residual_range_high, &
                                                                   tmp_search_confidence_interval, tmp_search_trace_k_start, &
                                                                   tmp_search_trace_k_step, tmp_search_trace_k_max, &
                                                                   tmp_search_trace_n_points, &
                                                                   tmp_search_trace_global_js_divergence, &
                                                                   tmp_search_trace_ci_lower, tmp_search_trace_ci_upper, &
                                                                   tmp_search_trace_ci_width, &
                                                                   tmp_search_trace_ci_width_relative, &
                                                                   tmp_search_trace_ci_overlap, tmp_search_trace_delta, &
                                                                   tmp_search_trace_delta_median, tmp_search_trace_delta_max, &
                                                                   tmp_search_trace_selected_n_bins, &
                                                                   tmp_search_trace_occupancy_failed, &
                                                                   tmp_search_trace_n_pooled_residuals, &
                                                                   tmp_search_trace_min_bin_occupancy, &
                                                                   tmp_search_trace_mean_bin_occupancy, &
                                                                   tmp_search_trace_max_bin_occupancy, &
                                                                   tmp_search_trace_sturges_bins, tmp_search_trace_fd_bins, &
                                                                   tmp_search_trace_shared_residual_range_low, &
                                                                   tmp_search_trace_shared_residual_range_high, &
                                                                   tmp_search_candidate_k_start, &
                                                                   tmp_search_candidate_n_points, &
                                                                   tmp_search_candidate_status, &
                                                                   tmp_candidate_confidence_interval, &
                                                                   tmp_bootstrapping_top_k_jsds, &
                                                                   tmp_prev_global_js_divergence, tmp_delta_perm, &
                                                                   tmp_best_uncertainty_confidence_interval, &
                                                                   tmp_best_n_bins_per_point, &
                                                                   tmp_best_uncertainty_n_bins_per_point, &
                                                                   tmp_best_shared_residual_range_low, &
                                                                   tmp_best_shared_residual_range_high, &
                                                                   tmp_best_uncertainty_shared_residual_range_low, &
                                                                   tmp_best_uncertainty_shared_residual_range_high, &
                                                                   rerun_parameter_search, n_bootstraps, join_method, &
                                                                   min_residuals_per_bin, min_neighbor_overlap, &
                                                                   succeeding_ci_overlap, plateau_mode, &
                                                                   delta_median_threshold, delta_max_threshold, &
                                                                   delta_epsilon, delta_min_consecutive_transitions, m_min, &
                                                                   m_max, gamma_occupancy, lower_residual_range_quantile, &
                                                                   upper_residual_range_quantile, &
                                                                   two_sided_bootstrapping_significance_level, random_seed, &
                                                                   tau, mad_distance_factor, max_pooled_residuals, &
                                                                   min_study_neighbors, bin_count_mode, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! DM_MIN(1_int32)
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: n_samples
            !! Number of samples, the non-padding slots of `expression`
            !! DM_MIN(1_int32)
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
            !! A permutation of `1..n_samples`: destination sample `k` receives source sample
            !! `sample_assignment(k)`; the identity is the observed data
            !! DM_MIN(1_int32)
            !! DM_MAX(n_samples)
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! DM_MIN(k_start)
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays; a construction with more neighborhoods fails
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_adaptive_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_impl, AUTO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! DM_OUTPUT_FROM(n_bootstrapping_top_k_jsds, calc_adaptive_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_impl, AUTO)
            !! DM_MIN(1_int32)
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
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: &
            tmp_permuted_expression
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
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_draw
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
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values
            !! Working array: the run's p-values (NaN, its permutation test does not run)
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values_observed_consensus
            !! Working array: the run's comparison p-values (NaN, its permutation test does not run)
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_counts
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_permutation_pmfs
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
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_start
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_step
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_max
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_high
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_k_start
            !! Working array: the search's candidate log
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_n_points
            !! Working array: the search's candidate log
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_status
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
            !! DM_DEFAULT(CM_RERUN_PARAMETER_SEARCH_DEFAULT)
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_N_BOOTSTRAPS_DEFAULT)
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! DM_DEFAULT(M_METHOD_JOIN_MIN)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_MIN_NEIGHBOR_OVERLAP_DEFAULT)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_SUCCEEDING_CI_OVERLAP_DEFAULT)
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(M_MODE_PLATEAU_EFFECT_SIZE)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
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
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(in), optional :: random_seed
            !! Rerun search: seed of the bootstrap's random number generator
            !! DM_DEFAULT(M_RANDOM_SEED_DEFAULT)
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_TAU_DEFAULT)
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have (construction status; in
            !! the rerun search, candidate admissibility)
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
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
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        integer(int32) :: n_samples_found, actual_n_bootstraps, actual_join_method, callee_ierr, search_capacity
        integer(int32) :: max_n_bins_per_point, search_n_points, search_n_admissible_evaluated, search_n_candidates_tried
        integer(int32) :: max_n_neighbors, construction_status
        logical(c_bool) :: actual_rerun, has_empty_study, is_permutation, search_plateau_established

        call set_ok(ierr)
        M_DEFAULT_VAL(rerun_parameter_search, actual_rerun, CM_RERUN_PARAMETER_SEARCH_DEFAULT)
        M_DEFAULT_VAL(n_bootstraps, actual_n_bootstraps, CM_N_BOOTSTRAPS_DEFAULT)
        M_DEFAULT_VAL(join_method, actual_join_method, M_METHOD_JOIN_MIN)
        global_js_divergence = M_NAN
        used_k_start = k_start
        used_k_step = k_step
        used_k_max = k_max
        n_points = 0_int32

        ! Runtime checks (codegen guide Sec 5.14), each naming the argument at fault.
        call locate_samples(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, tmp_sample_slots, &
                            n_samples_found, has_empty_study)
        if (has_empty_study) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
            return
        end if
        if (n_samples_found /= n_samples) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
            return
        end if
        call check_sample_assignment(n_samples, sample_assignment, tmp_sample_mark, is_permutation)
        if (.not. is_permutation) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=6_int32)
            return
        end if
        if (actual_rerun) then
            ! The search gets the capacity a direct call would, so it rejects the same candidates.
            call calc_adaptive_js_comp_test_bounds(max_n_genes_all_studies, n_studies, search_capacity, callee_ierr)
            if (is_err(callee_ierr)) then
                call set_err_once(ierr, get_err_code(callee_ierr))
                return
            end if
            if (max_n_points_candidate < search_capacity) then
                call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=10_int32)
                return
            end if
        end if

        call assign_samples_and_center(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, n_samples, &
                                       sample_assignment, tmp_sample_slots, tmp_permuted_expression, tmp_gene_means, &
                                       tmp_residuals)

        ! Every capacity-sized buffer is passed whole to explicit-shape dummies declared at the
        ! sizes actually used (sequence association, no copy).
        if (actual_rerun) then
            call run_js_comp_test_adaptive_parameter_search_impl( &
                n_studies=n_studies, max_n_genes_all_studies=max_n_genes_all_studies, &
                max_n_reps_all_studies=max_n_reps_all_studies, gene_means=tmp_gene_means, residuals=tmp_residuals, &
                n_bootstraps=actual_n_bootstraps, join_method=actual_join_method, max_n_points_candidate=search_capacity, &
                n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, k_start=used_k_start, k_step=used_k_step, &
                k_max=used_k_max, n_points=search_n_points, n_bins_per_point=tmp_search_n_bins_per_point, &
                shared_residual_range_low=tmp_search_shared_residual_range_low, &
                shared_residual_range_high=tmp_search_shared_residual_range_high, &
                best_candidate_confidence_interval=tmp_search_confidence_interval, &
                plateau_established=search_plateau_established, n_admissible_evaluated=search_n_admissible_evaluated, &
                trace_k_start=tmp_search_trace_k_start, trace_k_step=tmp_search_trace_k_step, &
                trace_k_max=tmp_search_trace_k_max, trace_n_points=tmp_search_trace_n_points, &
                trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
                trace_ci_lower=tmp_search_trace_ci_lower, trace_ci_upper=tmp_search_trace_ci_upper, &
                trace_ci_width=tmp_search_trace_ci_width, trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
                trace_ci_overlap=tmp_search_trace_ci_overlap, trace_delta=tmp_search_trace_delta, &
                trace_delta_median=tmp_search_trace_delta_median, trace_delta_max=tmp_search_trace_delta_max, &
                trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
                trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
                trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
                trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
                trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
                trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
                trace_sturges_bins=tmp_search_trace_sturges_bins, trace_fd_bins=tmp_search_trace_fd_bins, &
                trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
                trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
                n_candidates_tried=search_n_candidates_tried, candidate_k_start=tmp_search_candidate_k_start, &
                candidate_n_points=tmp_search_candidate_n_points, candidate_status=tmp_search_candidate_status, &
                tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
                tmp_pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
                tmp_n_neighbors_per_point=tmp_n_neighbors_per_point, tmp_stop_reason=tmp_stop_reason, &
                tmp_neighborhood_dispersion=tmp_neighborhood_dispersion, tmp_neighborhood_mad=tmp_neighborhood_mad, &
                tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
                tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_n_bins_per_point=tmp_n_bins_per_point, &
                tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
                tmp_shared_residual_range_high=tmp_shared_residual_range_high, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
                tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
                tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
                tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
                tmp_global_js_divergence=tmp_global_js_divergence, tmp_point_pool=tmp_point_pool, &
                tmp_point_draw=tmp_point_draw, tmp_confidence_interval=tmp_candidate_confidence_interval, &
                tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
                tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, tmp_delta_perm=tmp_delta_perm, &
                tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
                tmp_pooled_residuals=tmp_pooled_residuals, tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
                tmp_bin_counts_search=tmp_bin_counts_search, tmp_residual_to_study_perm=tmp_residual_to_study_perm, &
                tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
                tmp_occupancy_failed=tmp_occupancy_failed, tmp_n_pooled_residuals=tmp_n_pooled_residuals, &
                tmp_min_bin_occupancy=tmp_min_bin_occupancy, tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, &
                tmp_max_bin_occupancy=tmp_max_bin_occupancy, tmp_sturges_bins=tmp_sturges_bins, tmp_fd_bins=tmp_fd_bins, &
                tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
                tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
                tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
                tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
                tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
                tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
                min_residuals_per_bin=min_residuals_per_bin, min_neighbor_overlap=min_neighbor_overlap, &
                succeeding_ci_overlap=succeeding_ci_overlap, plateau_mode=plateau_mode, &
                delta_median_threshold=delta_median_threshold, delta_max_threshold=delta_max_threshold, &
                delta_epsilon=delta_epsilon, delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
                m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
                lower_residual_range_quantile=lower_residual_range_quantile, &
                upper_residual_range_quantile=upper_residual_range_quantile, &
                two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
                random_seed=random_seed, tau=tau, mad_distance_factor=mad_distance_factor, &
                max_pooled_residuals=max_pooled_residuals, min_study_neighbors=min_study_neighbors, &
                bin_count_mode=bin_count_mode, ierr=callee_ierr)
            if (is_err(callee_ierr)) then
                call set_err_once(ierr, get_err_code(callee_ierr))
                return
            end if
        end if

        call construct_adaptive_neighborhoods_impl( &
            n_studies=n_studies, max_n_genes_all_studies=max_n_genes_all_studies, &
            max_n_reps_all_studies=max_n_reps_all_studies, gene_means=tmp_gene_means, residuals=tmp_residuals, &
            k_start=used_k_start, k_step=used_k_step, k_max=used_k_max, n_points=n_points, x_star=tmp_x_star, &
            pooled_neighborhood_range=tmp_pooled_neighborhood_range, n_neighbors_per_point=tmp_n_neighbors_per_point, &
            stop_reason=tmp_stop_reason, neighborhood_dispersion=tmp_neighborhood_dispersion, &
            neighborhood_mad=tmp_neighborhood_mad, max_n_neighbors=max_n_neighbors, &
            construction_status=construction_status, tmp_gene_means_perm_all=tmp_gene_means_perm_all, tau=tau, &
            mad_distance_factor=mad_distance_factor, max_pooled_residuals=max_pooled_residuals, &
            min_study_neighbors=min_study_neighbors)
        ! Not one argument's fault: the relabelled data decide how many neighborhoods emerge.
        if (n_points < 1_int32 .or. n_points > max_n_points_candidate) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        call run_js_comp_test_adaptive_impl( &
            n_studies=n_studies, max_n_genes_all_studies=max_n_genes_all_studies, &
            max_n_reps_all_studies=max_n_reps_all_studies, n_points=n_points, gene_means=tmp_gene_means, &
            residuals=tmp_residuals, pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
            n_neighbors_per_point=tmp_n_neighbors_per_point, n_bins_per_point=tmp_n_bins_per_point, &
            shared_residual_range_low=tmp_shared_residual_range_low, &
            shared_residual_range_high=tmp_shared_residual_range_high, max_n_bins_per_point=max_n_bins_per_point, &
            occupancy_failed=tmp_occupancy_failed, n_pooled_residuals=tmp_n_pooled_residuals, &
            min_bin_occupancy=tmp_min_bin_occupancy, mean_bin_occupancy=tmp_mean_bin_occupancy, &
            max_bin_occupancy=tmp_max_bin_occupancy, sturges_bins=tmp_sturges_bins, fd_bins=tmp_fd_bins, pmfs=tmp_pmfs, &
            counts=tmp_counts, included_n_reps=tmp_included_n_reps, mean_pmf=tmp_mean_pmf, &
            mean_pmf_counts=tmp_mean_pmf_counts, mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            js_divergences=tmp_js_divergences, weights=tmp_weights, global_js_divergence=global_js_divergence, &
            p_values=tmp_p_values, p_values_observed_consensus=tmp_p_values_observed_consensus, &
            tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
            tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_permutation_mean_pmf=tmp_permutation_mean_pmf, &
            tmp_permutation_point_pool=tmp_point_pool, tmp_permutation_point_draw=tmp_point_draw, &
            tmp_pooled_residuals=tmp_pooled_residuals, tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
            tmp_bin_counts_search=tmp_bin_counts_search, tmp_residual_to_study_perm=tmp_residual_to_study_perm, &
            tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
            tmp_permutation_mean_pmf_counts=tmp_permutation_mean_pmf_counts, &
            tmp_permutation_counts=tmp_permutation_counts, tmp_permutation_pmfs=tmp_permutation_pmfs, &
            tmp_permutation_js_divergences=tmp_permutation_js_divergences, &
            tmp_permutation_weights=tmp_permutation_weights, &
            tmp_permutation_global_js_divergence=tmp_global_js_divergence, n_permutations=0_int32, &
            min_residuals_per_bin=min_residuals_per_bin, m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
            lower_residual_range_quantile=lower_residual_range_quantile, &
            upper_residual_range_quantile=upper_residual_range_quantile, bin_count_mode=bin_count_mode, ierr=callee_ierr)
        if (is_err(callee_ierr)) then
            call set_err_once(ierr, get_err_code(callee_ierr))
            global_js_divergence = M_NAN
            return
        end if
    end subroutine run_js_comp_test_adaptive_on_sample_assignment_impl

    !> summary: Draw one uniformly random permutation of `1..n_samples` (unbiased Fisher-Yates)
    !| AUTHOR_LASZLO_LANG
    !| Private. `sample_assignment` is reset to the identity, then for `i = n_samples, ..., 2` element
    !| `i` is swapped with element `j = min(1 + int(u*i), i)`, `u` uniform in `[0, 1)`. A failure of the
    !| generator is returned in `ierr` and stops the draw.
    subroutine draw_sample_assignment(rng, n_samples, sample_assignment, ierr)
        type(rng_t), intent(in) :: rng
            !! The generator, created by the caller
        integer(int32), intent(in) :: n_samples
            !! Number of samples
        integer(int32), intent(out) :: sample_assignment(n_samples)
            !! The drawn permutation
        integer(int32), intent(out) :: ierr
            !! Error code of `random_uniform`

        integer(int32) :: i_sample, j_sample, swapped
        real(real64) :: uniform

        call set_ok(ierr)
        call init_perm(sample_assignment)
        ! Sequential: each swap depends on the previous ones.
        do i_sample = n_samples, 2, -1
            uniform = random_uniform(rng, ierr)
            if (is_err(ierr)) return
            j_sample = min(1_int32 + int(uniform*real(i_sample, real64), kind=int32), i_sample)
            swapped = sample_assignment(i_sample)
            sample_assignment(i_sample) = sample_assignment(j_sample)
            sample_assignment(j_sample) = swapped
        end do
    end subroutine draw_sample_assignment

    !> summary: Turn the permuted values into the study-label Permutation-Test's p-values
    !| AUTHOR_LASZLO_LANG
    !| Private. Per study, `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`, the
    !| existing test's Laplace correction and comparison; NaN rows are failed permutations and do not
    !| count. NaN when no permutation is valid.
    pure subroutine label_permutation_p_values(n_permutations, n_studies, global_js_divergence_observed, jsd_null, &
                                               p_values)
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations (rows of `jsd_null`)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), intent(in) :: global_js_divergence_observed(n_studies)
            !! Observed global JSD per study
        real(real64), intent(in) :: jsd_null(n_permutations, n_studies)
            !! Per-permutation global JSD per study; a NaN row is a failed permutation
        real(real64), intent(out) :: p_values(n_studies)
            !! Empirical p-value per study

        integer(int32) :: i_study, i_permutation, n_valid, n_at_least

        do concurrent(i_study=1:n_studies) local(i_permutation, n_valid, n_at_least) &
            shared(n_permutations, jsd_null, global_js_divergence_observed, p_values)
            n_valid = 0_int32
            n_at_least = 0_int32
            do i_permutation = 1, n_permutations
                if (ieee_is_nan(jsd_null(i_permutation, i_study))) cycle
                n_valid = n_valid + 1_int32
                if (jsd_null(i_permutation, i_study) >= global_js_divergence_observed(i_study)) &
                    n_at_least = n_at_least + 1_int32
            end do
            if (n_valid == 0_int32) then
                p_values(i_study) = M_NAN
            else
                p_values(i_study) = real(n_at_least + 1_int32, real64)/real(n_valid + 1_int32, real64)
            end if
        end do
    end subroutine label_permutation_p_values


    !> summary: Study-label Permutation-Test of the k-NN JS-Comp-Test
    !| AUTHOR_LASZLO_LANG
    !| Tests whether the assignment of samples to studies influences each study's JS-Comp-Test value
    !| (null hypothesis: it does not). The same data, parameters and options as
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]],
    !| whose input layout (all studies in `expression`, NaN for a missing gene, all-NaN replicate
    !| slots are padding) and assignment convention apply here too.
    !|
    !| 1. The observed data (the identity assignment) give `global_js_divergence_observed`.
    !| 2. `n_permutations` times, the study labels of the samples are shuffled without replacement (a
    !|    uniformly random permutation, unbiased Fisher-Yates, from one random number generator seeded
    !|    with `random_seed`); every study keeps its sample count. The whole JS-Comp-Test runs on the
    !|    shuffled data; row `i` of `jsd_null` holds its values.
    !| 3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`,
    !|    the Laplace correction of the residual-level test, so it is never 0 and at least
    !|    `1/(n_permutations + 1)`.
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
    subroutine run_js_comp_test_label_permutation_test_impl(n_studies, max_n_genes_all_studies, &
                                                            max_n_reps_all_studies, expression, n_points, &
                                                            n_neighbors, max_n_points_candidate, &
                                                            max_n_neighbors_candidate, n_bootstrapping_top_k_jsds, &
                                                            n_permutations, global_js_divergence_observed, p_values, &
                                                            jsd_null, n_failed_permutations, tmp_sample_slots, &
                                                            tmp_sample_assignment, tmp_permuted_global_js_divergence, &
                                                            tmp_sample_mark, tmp_permuted_expression, tmp_gene_means, &
                                                            tmp_residuals, tmp_gene_means_perms, &
                                                            tmp_gene_means_perm_all, tmp_x_star, &
                                                            tmp_neighborhood_indices, tmp_neighborhood_range, &
                                                            tmp_neighborhood_residuals_gathered, tmp_point_pool, &
                                                            tmp_point_draw, tmp_n_bins_per_point, &
                                                            tmp_shared_residual_range_low, &
                                                            tmp_shared_residual_range_high, tmp_occupancy_failed, &
                                                            tmp_n_pooled_residuals, tmp_min_bin_occupancy, &
                                                            tmp_mean_bin_occupancy, tmp_max_bin_occupancy, &
                                                            tmp_sturges_bins, tmp_fd_bins, tmp_pmfs, tmp_counts, &
                                                            tmp_included_n_reps, tmp_mean_pmf, tmp_mean_pmf_counts, &
                                                            tmp_mean_pmf_included_n_reps, tmp_js_divergences, &
                                                            tmp_weights, tmp_global_js_divergence, &
                                                            tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                                            tmp_residual_to_study_perm, tmp_bin_counts_search, &
                                                            tmp_study_block_ends, tmp_study_block_sizes, &
                                                            tmp_point_n_neighbors, tmp_p_values, &
                                                            tmp_p_values_observed_consensus, &
                                                            tmp_permutation_mean_pmf, &
                                                            tmp_permutation_mean_pmf_counts, tmp_permutation_counts, &
                                                            tmp_permutation_pmfs, tmp_permutation_js_divergences, &
                                                            tmp_permutation_weights, tmp_search_n_bins_per_point, &
                                                            tmp_search_shared_residual_range_low, &
                                                            tmp_search_shared_residual_range_high, &
                                                            tmp_search_confidence_interval, &
                                                            tmp_search_trace_n_points, tmp_search_trace_n_neighbors, &
                                                            tmp_search_trace_global_js_divergence, &
                                                            tmp_search_trace_ci_lower, tmp_search_trace_ci_upper, &
                                                            tmp_search_trace_ci_width, &
                                                            tmp_search_trace_ci_width_relative, &
                                                            tmp_search_trace_ci_overlap, tmp_search_trace_delta, &
                                                            tmp_search_trace_delta_median, &
                                                            tmp_search_trace_delta_max, &
                                                            tmp_search_trace_selected_n_bins, &
                                                            tmp_search_trace_occupancy_failed, &
                                                            tmp_search_trace_n_pooled_residuals, &
                                                            tmp_search_trace_min_bin_occupancy, &
                                                            tmp_search_trace_mean_bin_occupancy, &
                                                            tmp_search_trace_max_bin_occupancy, &
                                                            tmp_search_trace_sturges_bins, tmp_search_trace_fd_bins, &
                                                            tmp_search_trace_shared_residual_range_low, &
                                                            tmp_search_trace_shared_residual_range_high, &
                                                            tmp_candidate_confidence_interval, &
                                                            tmp_bootstrapping_top_k_jsds, &
                                                            tmp_prev_global_js_divergence, &
                                                            tmp_best_uncertainty_confidence_interval, &
                                                            tmp_best_n_bins_per_point, &
                                                            tmp_best_uncertainty_n_bins_per_point, &
                                                            tmp_best_shared_residual_range_low, &
                                                            tmp_best_shared_residual_range_high, &
                                                            tmp_best_uncertainty_shared_residual_range_low, &
                                                            tmp_best_uncertainty_shared_residual_range_high, &
                                                            rerun_parameter_search, n_bootstraps, join_method, &
                                                            min_residuals_per_bin, min_neighbor_overlap, &
                                                            succeeding_ci_overlap, plateau_mode, &
                                                            delta_median_threshold, delta_max_threshold, &
                                                            delta_epsilon, delta_min_consecutive_transitions, m_min, &
                                                            m_max, gamma_occupancy, lower_residual_range_quantile, &
                                                            upper_residual_range_quantile, &
                                                            two_sided_bootstrapping_significance_level, random_seed, &
                                                            bin_count_mode, &
                                                            neighborhood_size_mode, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! DM_MIN(1_int32)
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies)
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_impl, AUTO)
            !! DM_MIN(n_points)
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Capacity of the per-neighbor work arrays
            !! DM_OUTPUT_FROM(max_n_neighbors_candidate, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_impl, AUTO)
            !! DM_MIN(n_neighbors)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! DM_OUTPUT_FROM(n_bootstrapping_top_k_jsds, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_impl, AUTO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_permutations
            !! Number of random assignments (rows of `jsd_null`)
            !! DM_MIN(1_int32)
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
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: &
            tmp_permuted_expression
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
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: &
            tmp_neighborhood_indices
            !! Working array of the search and the run: every study's neighbor gene indices
        integer(int32), dimension(2, max_n_points_candidate, n_studies), intent(out) :: tmp_neighborhood_range
            !! Working array of the search and the run: neighborhood spans
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), &
            intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array of the search and the run: one study's gathered neighborhood residuals
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_draw
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
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
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
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_counts
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_permutation_pmfs
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
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_neighbors
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_high
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
            !! DM_DEFAULT(CM_RERUN_PARAMETER_SEARCH_DEFAULT)
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_N_BOOTSTRAPS_DEFAULT)
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! DM_DEFAULT(M_METHOD_JOIN_MIN)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_MIN_NEIGHBOR_OVERLAP_DEFAULT)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_SUCCEEDING_CI_OVERLAP_DEFAULT)
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(M_MODE_PLATEAU_EFFECT_SIZE)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
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
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(in), optional :: random_seed
            !! Seed of the random number generator of the label draws, and of the bootstrap of
            !! every rerun parameter search (a separate generator)
            !! DM_DEFAULT(M_RANDOM_SEED_DEFAULT)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
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
        integer(int32), intent(in), optional :: neighborhood_size_mode
            !! Rerun search: which neighborhood sizes the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are left out
            !! by default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Large neighborhoods only, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods only, `KX_FACTORS(2)` = 0.5 | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! | Both, large then small for each `n_points` | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_BOTH(variable)]] |
            !! DM_DEFAULT(M_MODE_NEIGHBORHOOD_SIZE_LARGE)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        integer(int32) :: n_samples, i_permutation, actual_random_seed, callee_ierr
        integer(int32) :: used_n_points, used_n_neighbors
        logical(c_bool) :: actual_rerun, has_empty_study
        logical(c_bool) :: points_fit, neighbors_fit
        type(rng_t) :: rng

        call set_ok(ierr)
        M_DEFAULT_VAL(rerun_parameter_search, actual_rerun, CM_RERUN_PARAMETER_SEARCH_DEFAULT)
        M_DEFAULT_VAL(random_seed, actual_random_seed, M_RANDOM_SEED_DEFAULT)
        global_js_divergence_observed = M_NAN
        p_values = M_NAN
        jsd_null = M_NAN
        n_failed_permutations = 0_int32

        ! The checks that do not depend on the assignment, here with this routine's argument
        ! positions; the level-(a) routine repeats them, and its errors are passed on by code only.
        call locate_samples(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, tmp_sample_slots, &
                            n_samples, has_empty_study)
        if (has_empty_study) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
            return
        end if
        if (actual_rerun) then
            call check_fixed_k_rerun_capacities(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate, &
                                                points_fit, neighbors_fit)
            if (.not. points_fit) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
            if (.not. neighbors_fit) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=8_int32)
            if (is_err(ierr)) return
        end if

        ! The observed data: the identity assignment. An error here ends the test at once.
        call init_perm(tmp_sample_assignment)
        call run_js_comp_test_on_sample_assignment_impl(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_samples=n_samples, sample_assignment=tmp_sample_assignment, n_points=n_points, &
            n_neighbors=n_neighbors, max_n_points_candidate=max_n_points_candidate, &
            max_n_neighbors_candidate=max_n_neighbors_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, &
            global_js_divergence=global_js_divergence_observed, used_n_points=used_n_points, &
            used_n_neighbors=used_n_neighbors, tmp_sample_slots=tmp_sample_slots, tmp_sample_mark=tmp_sample_mark, &
            tmp_permuted_expression=tmp_permuted_expression, tmp_gene_means=tmp_gene_means, &
            tmp_residuals=tmp_residuals, tmp_gene_means_perms=tmp_gene_means_perms, &
            tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_neighborhood_indices=tmp_neighborhood_indices, tmp_neighborhood_range=tmp_neighborhood_range, &
            tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, tmp_point_pool=tmp_point_pool, &
            tmp_point_draw=tmp_point_draw, tmp_n_bins_per_point=tmp_n_bins_per_point, &
            tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
            tmp_shared_residual_range_high=tmp_shared_residual_range_high, tmp_occupancy_failed=tmp_occupancy_failed, &
            tmp_n_pooled_residuals=tmp_n_pooled_residuals, tmp_min_bin_occupancy=tmp_min_bin_occupancy, &
            tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, tmp_max_bin_occupancy=tmp_max_bin_occupancy, &
            tmp_sturges_bins=tmp_sturges_bins, tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
            tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
            tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
            tmp_global_js_divergence=tmp_global_js_divergence, tmp_pooled_residuals=tmp_pooled_residuals, &
            tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
            tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
            tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
            tmp_point_n_neighbors=tmp_point_n_neighbors, tmp_p_values=tmp_p_values, &
            tmp_p_values_observed_consensus=tmp_p_values_observed_consensus, &
            tmp_permutation_mean_pmf=tmp_permutation_mean_pmf, &
            tmp_permutation_mean_pmf_counts=tmp_permutation_mean_pmf_counts, &
            tmp_permutation_counts=tmp_permutation_counts, tmp_permutation_pmfs=tmp_permutation_pmfs, &
            tmp_permutation_js_divergences=tmp_permutation_js_divergences, &
            tmp_permutation_weights=tmp_permutation_weights, tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_n_neighbors=tmp_search_trace_n_neighbors, &
            tmp_search_trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
            tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
            tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
            tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
            tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, tmp_search_trace_delta=tmp_search_trace_delta, &
            tmp_search_trace_delta_median=tmp_search_trace_delta_median, &
            tmp_search_trace_delta_max=tmp_search_trace_delta_max, &
            tmp_search_trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
            tmp_search_trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
            tmp_search_trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
            tmp_search_trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
            tmp_search_trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
            tmp_search_trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
            tmp_search_trace_sturges_bins=tmp_search_trace_sturges_bins, &
            tmp_search_trace_fd_bins=tmp_search_trace_fd_bins, &
            tmp_search_trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
            tmp_search_trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
            tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
            tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
            tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, &
            tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
            tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
            tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
            tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
            tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
            tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
            tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
            rerun_parameter_search=actual_rerun, n_bootstraps=n_bootstraps, join_method=join_method, &
            min_residuals_per_bin=min_residuals_per_bin, min_neighbor_overlap=min_neighbor_overlap, &
            succeeding_ci_overlap=succeeding_ci_overlap, plateau_mode=plateau_mode, &
            delta_median_threshold=delta_median_threshold, delta_max_threshold=delta_max_threshold, &
            delta_epsilon=delta_epsilon, delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
            m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
            lower_residual_range_quantile=lower_residual_range_quantile, &
            upper_residual_range_quantile=upper_residual_range_quantile, &
            two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
            random_seed=actual_random_seed, bin_count_mode=bin_count_mode, &
                neighborhood_size_mode=neighborhood_size_mode, ierr=callee_ierr)
        if (is_err(callee_ierr)) then
            call set_err_once(ierr, get_err_code(callee_ierr))
            global_js_divergence_observed = M_NAN
            return
        end if

        rng = create_rng(callee_ierr, actual_random_seed)
        if (is_err(callee_ierr)) then
            call set_err_once(ierr, get_err_code(callee_ierr))
            call destroy_rng(rng)
            return
        end if
        ! Sequential: one generator draws every assignment in turn.
        do i_permutation = 1, n_permutations
            call draw_sample_assignment(rng, n_samples, tmp_sample_assignment, callee_ierr)
            if (is_err(callee_ierr)) then
                call set_err_once(ierr, get_err_code(callee_ierr))
                exit
            end if
            call run_js_comp_test_on_sample_assignment_impl(n_studies=n_studies, &
                max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
                expression=expression, n_samples=n_samples, sample_assignment=tmp_sample_assignment, &
                n_points=n_points, n_neighbors=n_neighbors, max_n_points_candidate=max_n_points_candidate, &
                max_n_neighbors_candidate=max_n_neighbors_candidate, &
                n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, &
                global_js_divergence=tmp_permuted_global_js_divergence, used_n_points=used_n_points, &
                used_n_neighbors=used_n_neighbors, tmp_sample_slots=tmp_sample_slots, &
                tmp_sample_mark=tmp_sample_mark, tmp_permuted_expression=tmp_permuted_expression, &
                tmp_gene_means=tmp_gene_means, tmp_residuals=tmp_residuals, &
                tmp_gene_means_perms=tmp_gene_means_perms, tmp_gene_means_perm_all=tmp_gene_means_perm_all, &
                tmp_x_star=tmp_x_star, tmp_neighborhood_indices=tmp_neighborhood_indices, &
                tmp_neighborhood_range=tmp_neighborhood_range, &
                tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, &
                tmp_point_pool=tmp_point_pool, tmp_point_draw=tmp_point_draw, &
                tmp_n_bins_per_point=tmp_n_bins_per_point, &
                tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
                tmp_shared_residual_range_high=tmp_shared_residual_range_high, &
                tmp_occupancy_failed=tmp_occupancy_failed, tmp_n_pooled_residuals=tmp_n_pooled_residuals, &
                tmp_min_bin_occupancy=tmp_min_bin_occupancy, tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, &
                tmp_max_bin_occupancy=tmp_max_bin_occupancy, tmp_sturges_bins=tmp_sturges_bins, &
                tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
                tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
                tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
                tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
                tmp_global_js_divergence=tmp_global_js_divergence, tmp_pooled_residuals=tmp_pooled_residuals, &
                tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
                tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
                tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
                tmp_point_n_neighbors=tmp_point_n_neighbors, tmp_p_values=tmp_p_values, &
                tmp_p_values_observed_consensus=tmp_p_values_observed_consensus, &
                tmp_permutation_mean_pmf=tmp_permutation_mean_pmf, &
                tmp_permutation_mean_pmf_counts=tmp_permutation_mean_pmf_counts, &
                tmp_permutation_counts=tmp_permutation_counts, tmp_permutation_pmfs=tmp_permutation_pmfs, &
                tmp_permutation_js_divergences=tmp_permutation_js_divergences, &
                tmp_permutation_weights=tmp_permutation_weights, &
                tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
                tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
                tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
                tmp_search_confidence_interval=tmp_search_confidence_interval, &
                tmp_search_trace_n_points=tmp_search_trace_n_points, &
                tmp_search_trace_n_neighbors=tmp_search_trace_n_neighbors, &
                tmp_search_trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
                tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, &
                tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
                tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
                tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
                tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, &
                tmp_search_trace_delta=tmp_search_trace_delta, &
                tmp_search_trace_delta_median=tmp_search_trace_delta_median, &
                tmp_search_trace_delta_max=tmp_search_trace_delta_max, &
                tmp_search_trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
                tmp_search_trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
                tmp_search_trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
                tmp_search_trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
                tmp_search_trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
                tmp_search_trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
                tmp_search_trace_sturges_bins=tmp_search_trace_sturges_bins, &
                tmp_search_trace_fd_bins=tmp_search_trace_fd_bins, &
                tmp_search_trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
                tmp_search_trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
                tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
                tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
                tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, &
                tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
                tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
                tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
                tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
                tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
                tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
                tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
                rerun_parameter_search=actual_rerun, n_bootstraps=n_bootstraps, join_method=join_method, &
                min_residuals_per_bin=min_residuals_per_bin, min_neighbor_overlap=min_neighbor_overlap, &
                succeeding_ci_overlap=succeeding_ci_overlap, plateau_mode=plateau_mode, &
                delta_median_threshold=delta_median_threshold, delta_max_threshold=delta_max_threshold, &
                delta_epsilon=delta_epsilon, delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
                m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
                lower_residual_range_quantile=lower_residual_range_quantile, &
                upper_residual_range_quantile=upper_residual_range_quantile, &
                two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
                random_seed=actual_random_seed, bin_count_mode=bin_count_mode, &
                neighborhood_size_mode=neighborhood_size_mode, ierr=callee_ierr)
            if (is_err(callee_ierr)) then
                ! Invalid input of a relabelled data set (e.g. a study without residuals at some
                ! reference point) fails this permutation only; anything else ends the test.
                if (get_err_code(callee_ierr) == ERR_INVALID_INPUT) then
                    n_failed_permutations = n_failed_permutations + 1_int32
                    cycle
                end if
                call set_err_once(ierr, get_err_code(callee_ierr))
                exit
            end if
            jsd_null(i_permutation, :) = tmp_permuted_global_js_divergence
        end do
        call destroy_rng(rng)
        if (is_err(ierr)) return

        call label_permutation_p_values(n_permutations, n_studies, global_js_divergence_observed, jsd_null, p_values)
    end subroutine run_js_comp_test_label_permutation_test_impl

    !> summary: Study-label Permutation-Test of the adaptive JS-Comp-Test
    !| AUTHOR_LASZLO_LANG
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
    subroutine run_js_comp_test_adaptive_label_permutation_test_impl(n_studies, max_n_genes_all_studies, &
                                                                     max_n_reps_all_studies, expression, k_start, &
                                                                     k_step, k_max, max_n_points_candidate, &
                                                                     n_bootstrapping_top_k_jsds, n_permutations, &
                                                                     global_js_divergence_observed, p_values, &
                                                                     jsd_null, n_failed_permutations, &
                                                                     tmp_sample_slots, tmp_sample_assignment, &
                                                                     tmp_permuted_global_js_divergence, &
                                                                     tmp_sample_mark, tmp_permuted_expression, &
                                                                     tmp_gene_means, tmp_residuals, &
                                                                     tmp_gene_means_perm_all, tmp_x_star, &
                                                                     tmp_pooled_neighborhood_range, &
                                                                     tmp_n_neighbors_per_point, tmp_stop_reason, &
                                                                     tmp_neighborhood_dispersion, &
                                                                     tmp_neighborhood_mad, &
                                                                     tmp_point_neighborhood_indices, &
                                                                     tmp_neighbor_residuals, tmp_point_pool, &
                                                                     tmp_point_draw, tmp_n_bins_per_point, &
                                                                     tmp_shared_residual_range_low, &
                                                                     tmp_shared_residual_range_high, &
                                                                     tmp_occupancy_failed, tmp_n_pooled_residuals, &
                                                                     tmp_min_bin_occupancy, tmp_mean_bin_occupancy, &
                                                                     tmp_max_bin_occupancy, tmp_sturges_bins, &
                                                                     tmp_fd_bins, tmp_pmfs, tmp_counts, &
                                                                     tmp_included_n_reps, tmp_mean_pmf, &
                                                                     tmp_mean_pmf_counts, &
                                                                     tmp_mean_pmf_included_n_reps, &
                                                                     tmp_js_divergences, tmp_weights, &
                                                                     tmp_global_js_divergence, tmp_pooled_residuals, &
                                                                     tmp_pooled_residuals_perm, &
                                                                     tmp_residual_to_study_perm, &
                                                                     tmp_bin_counts_search, tmp_study_block_ends, &
                                                                     tmp_study_block_sizes, tmp_p_values, &
                                                                     tmp_p_values_observed_consensus, &
                                                                     tmp_permutation_mean_pmf, &
                                                                     tmp_permutation_mean_pmf_counts, &
                                                                     tmp_permutation_counts, tmp_permutation_pmfs, &
                                                                     tmp_permutation_js_divergences, &
                                                                     tmp_permutation_weights, &
                                                                     tmp_search_n_bins_per_point, &
                                                                     tmp_search_shared_residual_range_low, &
                                                                     tmp_search_shared_residual_range_high, &
                                                                     tmp_search_confidence_interval, &
                                                                     tmp_search_trace_k_start, &
                                                                     tmp_search_trace_k_step, tmp_search_trace_k_max, &
                                                                     tmp_search_trace_n_points, &
                                                                     tmp_search_trace_global_js_divergence, &
                                                                     tmp_search_trace_ci_lower, &
                                                                     tmp_search_trace_ci_upper, &
                                                                     tmp_search_trace_ci_width, &
                                                                     tmp_search_trace_ci_width_relative, &
                                                                     tmp_search_trace_ci_overlap, &
                                                                     tmp_search_trace_delta, &
                                                                     tmp_search_trace_delta_median, &
                                                                     tmp_search_trace_delta_max, &
                                                                     tmp_search_trace_selected_n_bins, &
                                                                     tmp_search_trace_occupancy_failed, &
                                                                     tmp_search_trace_n_pooled_residuals, &
                                                                     tmp_search_trace_min_bin_occupancy, &
                                                                     tmp_search_trace_mean_bin_occupancy, &
                                                                     tmp_search_trace_max_bin_occupancy, &
                                                                     tmp_search_trace_sturges_bins, &
                                                                     tmp_search_trace_fd_bins, &
                                                                     tmp_search_trace_shared_residual_range_low, &
                                                                     tmp_search_trace_shared_residual_range_high, &
                                                                     tmp_search_candidate_k_start, &
                                                                     tmp_search_candidate_n_points, &
                                                                     tmp_search_candidate_status, &
                                                                     tmp_candidate_confidence_interval, &
                                                                     tmp_bootstrapping_top_k_jsds, &
                                                                     tmp_prev_global_js_divergence, tmp_delta_perm, &
                                                                     tmp_best_uncertainty_confidence_interval, &
                                                                     tmp_best_n_bins_per_point, &
                                                                     tmp_best_uncertainty_n_bins_per_point, &
                                                                     tmp_best_shared_residual_range_low, &
                                                                     tmp_best_shared_residual_range_high, &
                                                                     tmp_best_uncertainty_shared_residual_range_low, &
                                                                     tmp_best_uncertainty_shared_residual_range_high, &
                                                                     rerun_parameter_search, n_bootstraps, &
                                                                     join_method, min_residuals_per_bin, &
                                                                     min_neighbor_overlap, succeeding_ci_overlap, &
                                                                     plateau_mode, delta_median_threshold, &
                                                                     delta_max_threshold, delta_epsilon, &
                                                                     delta_min_consecutive_transitions, m_min, m_max, &
                                                                     gamma_occupancy, lower_residual_range_quantile, &
                                                                     upper_residual_range_quantile, &
                                                                     two_sided_bootstrapping_significance_level, &
                                                                     random_seed, tau, mad_distance_factor, &
                                                                     max_pooled_residuals, min_study_neighbors, &
                                                                     bin_count_mode, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! DM_MIN(1_int32)
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! DM_MIN(k_start)
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays; a construction with more neighborhoods fails
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_adaptive_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_impl, AUTO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! DM_OUTPUT_FROM(n_bootstrapping_top_k_jsds, calc_adaptive_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_impl, AUTO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_permutations
            !! Number of random assignments (rows of `jsd_null`)
            !! DM_MIN(1_int32)
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
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: &
            tmp_permuted_expression
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
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_draw
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
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values
            !! Working array: the run's p-values (NaN, its permutation test does not run)
        real(real64), dimension(n_studies), intent(out) :: tmp_p_values_observed_consensus
            !! Working array: the run's comparison p-values (NaN, its permutation test does not run)
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array of the run's (skipped) permutation test
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_permutation_counts
            !! Working array of the run's (skipped) permutation test
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_permutation_pmfs
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
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_start
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_step
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_max
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_high
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_k_start
            !! Working array: the search's candidate log
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_n_points
            !! Working array: the search's candidate log
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_status
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
            !! DM_DEFAULT(CM_RERUN_PARAMETER_SEARCH_DEFAULT)
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_N_BOOTSTRAPS_DEFAULT)
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! DM_DEFAULT(M_METHOD_JOIN_MIN)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_MIN_NEIGHBOR_OVERLAP_DEFAULT)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_SUCCEEDING_CI_OVERLAP_DEFAULT)
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(M_MODE_PLATEAU_EFFECT_SIZE)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
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
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(in), optional :: random_seed
            !! Seed of the random number generator of the label draws, and of the bootstrap of
            !! every rerun parameter search (a separate generator)
            !! DM_DEFAULT(M_RANDOM_SEED_DEFAULT)
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_TAU_DEFAULT)
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have (construction status; in
            !! the rerun search, candidate admissibility)
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
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
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        integer(int32) :: n_samples, i_permutation, actual_random_seed, callee_ierr
        integer(int32) :: used_k_start, used_k_step, used_k_max, n_points, search_capacity
        logical(c_bool) :: actual_rerun, has_empty_study
        type(rng_t) :: rng

        call set_ok(ierr)
        M_DEFAULT_VAL(rerun_parameter_search, actual_rerun, CM_RERUN_PARAMETER_SEARCH_DEFAULT)
        M_DEFAULT_VAL(random_seed, actual_random_seed, M_RANDOM_SEED_DEFAULT)
        global_js_divergence_observed = M_NAN
        p_values = M_NAN
        jsd_null = M_NAN
        n_failed_permutations = 0_int32

        ! The checks that do not depend on the assignment, here with this routine's argument
        ! positions; the level-(a) routine repeats them, and its errors are passed on by code only.
        call locate_samples(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, tmp_sample_slots, &
                            n_samples, has_empty_study)
        if (has_empty_study) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
            return
        end if
        if (actual_rerun) then
            call calc_adaptive_js_comp_test_bounds(max_n_genes_all_studies, n_studies, search_capacity, callee_ierr)
            if (is_err(callee_ierr)) then
                call set_err_once(ierr, get_err_code(callee_ierr))
                return
            end if
            if (max_n_points_candidate < search_capacity) then
                call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=8_int32)
                return
            end if
        end if

        ! The observed data: the identity assignment. An error here ends the test at once.
        call init_perm(tmp_sample_assignment)
        call run_js_comp_test_adaptive_on_sample_assignment_impl(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_samples=n_samples, sample_assignment=tmp_sample_assignment, k_start=k_start, &
            k_step=k_step, k_max=k_max, max_n_points_candidate=max_n_points_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, &
            global_js_divergence=global_js_divergence_observed, used_k_start=used_k_start, used_k_step=used_k_step, &
            used_k_max=used_k_max, n_points=n_points, tmp_sample_slots=tmp_sample_slots, &
            tmp_sample_mark=tmp_sample_mark, tmp_permuted_expression=tmp_permuted_expression, &
            tmp_gene_means=tmp_gene_means, tmp_residuals=tmp_residuals, &
            tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
            tmp_n_neighbors_per_point=tmp_n_neighbors_per_point, tmp_stop_reason=tmp_stop_reason, &
            tmp_neighborhood_dispersion=tmp_neighborhood_dispersion, tmp_neighborhood_mad=tmp_neighborhood_mad, &
            tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
            tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_point_pool=tmp_point_pool, &
            tmp_point_draw=tmp_point_draw, tmp_n_bins_per_point=tmp_n_bins_per_point, &
            tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
            tmp_shared_residual_range_high=tmp_shared_residual_range_high, tmp_occupancy_failed=tmp_occupancy_failed, &
            tmp_n_pooled_residuals=tmp_n_pooled_residuals, tmp_min_bin_occupancy=tmp_min_bin_occupancy, &
            tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, tmp_max_bin_occupancy=tmp_max_bin_occupancy, &
            tmp_sturges_bins=tmp_sturges_bins, tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
            tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
            tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
            tmp_global_js_divergence=tmp_global_js_divergence, tmp_pooled_residuals=tmp_pooled_residuals, &
            tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
            tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
            tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
            tmp_p_values=tmp_p_values, tmp_p_values_observed_consensus=tmp_p_values_observed_consensus, &
            tmp_permutation_mean_pmf=tmp_permutation_mean_pmf, &
            tmp_permutation_mean_pmf_counts=tmp_permutation_mean_pmf_counts, &
            tmp_permutation_counts=tmp_permutation_counts, tmp_permutation_pmfs=tmp_permutation_pmfs, &
            tmp_permutation_js_divergences=tmp_permutation_js_divergences, &
            tmp_permutation_weights=tmp_permutation_weights, tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_k_start=tmp_search_trace_k_start, tmp_search_trace_k_step=tmp_search_trace_k_step, &
            tmp_search_trace_k_max=tmp_search_trace_k_max, tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
            tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
            tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
            tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
            tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, tmp_search_trace_delta=tmp_search_trace_delta, &
            tmp_search_trace_delta_median=tmp_search_trace_delta_median, &
            tmp_search_trace_delta_max=tmp_search_trace_delta_max, &
            tmp_search_trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
            tmp_search_trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
            tmp_search_trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
            tmp_search_trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
            tmp_search_trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
            tmp_search_trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
            tmp_search_trace_sturges_bins=tmp_search_trace_sturges_bins, &
            tmp_search_trace_fd_bins=tmp_search_trace_fd_bins, &
            tmp_search_trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
            tmp_search_trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
            tmp_search_candidate_k_start=tmp_search_candidate_k_start, &
            tmp_search_candidate_n_points=tmp_search_candidate_n_points, &
            tmp_search_candidate_status=tmp_search_candidate_status, &
            tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
            tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
            tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, tmp_delta_perm=tmp_delta_perm, &
            tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
            tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
            tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
            tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
            tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
            tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
            tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
            rerun_parameter_search=actual_rerun, n_bootstraps=n_bootstraps, join_method=join_method, &
            min_residuals_per_bin=min_residuals_per_bin, min_neighbor_overlap=min_neighbor_overlap, &
            succeeding_ci_overlap=succeeding_ci_overlap, plateau_mode=plateau_mode, &
            delta_median_threshold=delta_median_threshold, delta_max_threshold=delta_max_threshold, &
            delta_epsilon=delta_epsilon, delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
            m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
            lower_residual_range_quantile=lower_residual_range_quantile, &
            upper_residual_range_quantile=upper_residual_range_quantile, &
            two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
            random_seed=actual_random_seed, tau=tau, mad_distance_factor=mad_distance_factor, &
            max_pooled_residuals=max_pooled_residuals, min_study_neighbors=min_study_neighbors, &
            bin_count_mode=bin_count_mode, ierr=callee_ierr)
        if (is_err(callee_ierr)) then
            call set_err_once(ierr, get_err_code(callee_ierr))
            global_js_divergence_observed = M_NAN
            return
        end if

        rng = create_rng(callee_ierr, actual_random_seed)
        if (is_err(callee_ierr)) then
            call set_err_once(ierr, get_err_code(callee_ierr))
            call destroy_rng(rng)
            return
        end if
        ! Sequential: one generator draws every assignment in turn.
        do i_permutation = 1, n_permutations
            call draw_sample_assignment(rng, n_samples, tmp_sample_assignment, callee_ierr)
            if (is_err(callee_ierr)) then
                call set_err_once(ierr, get_err_code(callee_ierr))
                exit
            end if
            call run_js_comp_test_adaptive_on_sample_assignment_impl(n_studies=n_studies, &
                max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
                expression=expression, n_samples=n_samples, sample_assignment=tmp_sample_assignment, k_start=k_start, &
                k_step=k_step, k_max=k_max, max_n_points_candidate=max_n_points_candidate, &
                n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, &
                global_js_divergence=tmp_permuted_global_js_divergence, used_k_start=used_k_start, &
                used_k_step=used_k_step, used_k_max=used_k_max, n_points=n_points, tmp_sample_slots=tmp_sample_slots, &
                tmp_sample_mark=tmp_sample_mark, tmp_permuted_expression=tmp_permuted_expression, &
                tmp_gene_means=tmp_gene_means, tmp_residuals=tmp_residuals, &
                tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
                tmp_pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
                tmp_n_neighbors_per_point=tmp_n_neighbors_per_point, tmp_stop_reason=tmp_stop_reason, &
                tmp_neighborhood_dispersion=tmp_neighborhood_dispersion, tmp_neighborhood_mad=tmp_neighborhood_mad, &
                tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
                tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_point_pool=tmp_point_pool, &
                tmp_point_draw=tmp_point_draw, tmp_n_bins_per_point=tmp_n_bins_per_point, &
                tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
                tmp_shared_residual_range_high=tmp_shared_residual_range_high, &
                tmp_occupancy_failed=tmp_occupancy_failed, tmp_n_pooled_residuals=tmp_n_pooled_residuals, &
                tmp_min_bin_occupancy=tmp_min_bin_occupancy, tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, &
                tmp_max_bin_occupancy=tmp_max_bin_occupancy, tmp_sturges_bins=tmp_sturges_bins, &
                tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
                tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
                tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
                tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
                tmp_global_js_divergence=tmp_global_js_divergence, tmp_pooled_residuals=tmp_pooled_residuals, &
                tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
                tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
                tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
                tmp_p_values=tmp_p_values, tmp_p_values_observed_consensus=tmp_p_values_observed_consensus, &
                tmp_permutation_mean_pmf=tmp_permutation_mean_pmf, &
                tmp_permutation_mean_pmf_counts=tmp_permutation_mean_pmf_counts, &
                tmp_permutation_counts=tmp_permutation_counts, tmp_permutation_pmfs=tmp_permutation_pmfs, &
                tmp_permutation_js_divergences=tmp_permutation_js_divergences, &
                tmp_permutation_weights=tmp_permutation_weights, &
                tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
                tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
                tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
                tmp_search_confidence_interval=tmp_search_confidence_interval, &
                tmp_search_trace_k_start=tmp_search_trace_k_start, tmp_search_trace_k_step=tmp_search_trace_k_step, &
                tmp_search_trace_k_max=tmp_search_trace_k_max, tmp_search_trace_n_points=tmp_search_trace_n_points, &
                tmp_search_trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
                tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, &
                tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
                tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
                tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
                tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, &
                tmp_search_trace_delta=tmp_search_trace_delta, &
                tmp_search_trace_delta_median=tmp_search_trace_delta_median, &
                tmp_search_trace_delta_max=tmp_search_trace_delta_max, &
                tmp_search_trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
                tmp_search_trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
                tmp_search_trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
                tmp_search_trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
                tmp_search_trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
                tmp_search_trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
                tmp_search_trace_sturges_bins=tmp_search_trace_sturges_bins, &
                tmp_search_trace_fd_bins=tmp_search_trace_fd_bins, &
                tmp_search_trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
                tmp_search_trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
                tmp_search_candidate_k_start=tmp_search_candidate_k_start, &
                tmp_search_candidate_n_points=tmp_search_candidate_n_points, &
                tmp_search_candidate_status=tmp_search_candidate_status, &
                tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
                tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
                tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, tmp_delta_perm=tmp_delta_perm, &
                tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
                tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
                tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
                tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
                tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
                tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
                tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
                rerun_parameter_search=actual_rerun, n_bootstraps=n_bootstraps, join_method=join_method, &
                min_residuals_per_bin=min_residuals_per_bin, min_neighbor_overlap=min_neighbor_overlap, &
                succeeding_ci_overlap=succeeding_ci_overlap, plateau_mode=plateau_mode, &
                delta_median_threshold=delta_median_threshold, delta_max_threshold=delta_max_threshold, &
                delta_epsilon=delta_epsilon, delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
                m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
                lower_residual_range_quantile=lower_residual_range_quantile, &
                upper_residual_range_quantile=upper_residual_range_quantile, &
                two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
                random_seed=actual_random_seed, tau=tau, mad_distance_factor=mad_distance_factor, &
                max_pooled_residuals=max_pooled_residuals, min_study_neighbors=min_study_neighbors, &
                bin_count_mode=bin_count_mode, ierr=callee_ierr)
            if (is_err(callee_ierr)) then
                ! Invalid input of a relabelled data set (e.g. a study without residuals at some
                ! reference point) fails this permutation only; anything else ends the test.
                if (get_err_code(callee_ierr) == ERR_INVALID_INPUT) then
                    n_failed_permutations = n_failed_permutations + 1_int32
                    cycle
                end if
                call set_err_once(ierr, get_err_code(callee_ierr))
                exit
            end if
            jsd_null(i_permutation, :) = tmp_permuted_global_js_divergence
        end do
        call destroy_rng(rng)
        if (is_err(ierr)) return

        call label_permutation_p_values(n_permutations, n_studies, global_js_divergence_observed, jsd_null, p_values)
    end subroutine run_js_comp_test_adaptive_label_permutation_test_impl

end module tox_data_integration_js_comp_test_permutation_impl
