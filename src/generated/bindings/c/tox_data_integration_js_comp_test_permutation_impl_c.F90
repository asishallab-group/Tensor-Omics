#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_permutation_impl(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Study-Label Permutation-Test
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
module tox_data_integration_js_comp_test_permutation_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: calc_js_comp_test_permutation_bounds_c
    public :: calc_adaptive_js_comp_test_permutation_bounds_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_permutation_impl(module):calc_js_comp_test_permutation_bounds(subroutine)]]
    !| Sizes the per-point and per-neighbor work arrays of
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]]
    !| and
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_label_permutation_test_impl(interface)]]
    !| so that both the fixed `(n_points, n_neighbors)` run and a rerun of the parameter search fit:
    !| `max_n_points_candidate = max(n_points, P)` and `max_n_neighbors_candidate = max(n_neighbors,
    !| K)`, with `P` and `K` the bounds of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds_impl(subroutine)]].
    !| `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
    !| [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds(subroutine)]] computes it.
    !|
    !| Rejects an out-of-range argument with invalid input (a non-finite significance level with
    !| `ERR_NAN_INF`); the outputs are then left undefined.
    subroutine calc_js_comp_test_permutation_bounds_c(&
            max_n_genes_all_studies,&
            n_points,&
            n_neighbors,&
            n_bootstraps,&
            two_sided_bootstrapping_significance_level,&
            max_n_points_candidate,&
            max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds,&
            ierr&
        ) bind(C, name="calc_js_comp_test_permutation_bounds_c")
        use tox_data_integration_js_comp_test_permutation_impl, only: calc_js_comp_test_permutation_bounds

        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points of the fixed-parameter run
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors per neighborhood of the fixed-parameter run
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps of the parameter search, when it is rerun
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(out), target :: max_n_points_candidate
            !! Capacity for reference points
        integer(c_int), intent(out), target :: max_n_neighbors_candidate
            !! Capacity for neighbors per neighborhood
        integer(c_int), intent(out), target :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size (values kept at each end of the bootstrap distribution)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(max_n_neighbors_candidate)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)

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
    end subroutine calc_js_comp_test_permutation_bounds_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_permutation_impl(module):calc_adaptive_js_comp_test_permutation_bounds(subroutine)]]
    !| Sizes the per-point work arrays of
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]
    !| and
    !| [[tox_data_integration_js_comp_test_permutation_impl(module):run_js_comp_test_adaptive_label_permutation_test_impl(interface)]]:
    !| `max_n_points_candidate = max(A, min(N, ceiling(F * N / k_start) + 1))`, with `N =
    !| max_n_genes_all_studies * n_studies` (the padded pool size), `F` the factor
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_POINT_CAPACITY_FACTOR(variable)]]
    !| and `A` the capacity of the adaptive parameter search,
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):calc_adaptive_js_comp_test_bounds_impl(subroutine)]].
    !| As there, this is a **practical capacity, not a proven bound**: how many reference points a
    !| construction builds is only known once it has run, and a shuffled data set whose
    !| construction builds more fails with invalid input (it is then counted as a failed
    !| permutation). `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
    !| [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds(subroutine)]] computes it.
    !|
    !| Rejects an out-of-range argument with invalid input (a non-finite significance level with
    !| `ERR_NAN_INF`), and, naming neither argument, the joint bounds on `max_n_genes_all_studies`
    !| and `n_studies` that `calc_adaptive_js_comp_test_bounds_impl` rejects; the outputs are then left
    !| undefined.
    subroutine calc_adaptive_js_comp_test_permutation_bounds_c(&
            max_n_genes_all_studies,&
            n_studies,&
            k_start,&
            n_bootstraps,&
            two_sided_bootstrapping_significance_level,&
            max_n_points_candidate,&
            n_bootstrapping_top_k_jsds,&
            ierr&
        ) bind(C, name="calc_adaptive_js_comp_test_permutation_bounds_c")
        use tox_data_integration_js_comp_test_permutation_impl, only: calc_adaptive_js_comp_test_permutation_bounds

        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: k_start
            !! Pooled entries every neighborhood of the fixed-knob run takes unconditionally
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps of the parameter search, when it is rerun
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(out), target :: max_n_points_candidate
            !! Practical capacity for reference points
        integer(c_int), intent(out), target :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size (values kept at each end of the bootstrap distribution)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(k_start)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)

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
    end subroutine calc_adaptive_js_comp_test_permutation_bounds_c

end module tox_data_integration_js_comp_test_permutation_impl_c
#endif
