!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Pipeline
!|
!| The global JS-Comp-Test of K studies, with a data-driven choice of its neighborhood parameters.
!| Two kinds of neighborhoods are available; the adaptive ones are the recommended way, since they
!| follow the local density of the gene means.
!|
!| Adaptive neighborhoods, grown per reference point while the residual dispersion stays stable:
!| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]]
!| searches the `(k_start, k_step, k_max)` candidates of
!| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]],
!| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
!| builds the neighborhoods for the chosen knobs, and
!| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_adaptive_impl(interface)]] runs
!| the final test on them.
!|
!| Fixed-k neighborhoods (`n_points` reference points of `n_neighbors` genes each):
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
!| searches the `(n_points, n_neighbors)` candidates of
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]
!| for a bootstrap plateau, and
!| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]] runs the
!| final test (histograms, JSDs, weighting) for the chosen pair.
!|
!| Each study's divergence is, by default, for exactly two studies the pairwise JSD of the two
!| studies, and for any other number of studies the JSD of each study against the consensus of all
!| (`divergence_mode`). Each reference point's histogram bin count is chosen by Stone's
!| unbiased risk estimate
!| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]
!| (default: every study's score counted equally) or by the occupancy search
!| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
!| (`bin_count_mode`). The building blocks are published too: the admissibility gates
!| [[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]], the
!| plateau checks
!| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]],
!| the bootstrap
!| [[tox_data_integration_stats_impl(module):bootstrap_histogram_impl(interface)]], and the
!| brute-force reference
!| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]].
!| `calc_js_comp_test_candidate_bounds`, `calc_adaptive_js_comp_test_bounds` and
!| `calc_js_comp_test_n_top_k_jsds` size the work arrays for a caller that allocates its own.
!|
!| The significance test is the study-label Permutation-Test. It asks whether the assignment of
!| samples to studies influences the JS-Comp-Test values at all: it shuffles the study labels of the samples and reruns the whole
!| test on every shuffled data set.
!| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_permutation_impl(interface)]]
!| (k-NN) and
!| [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module):run_js_comp_test_adaptive_permutation_impl(interface)]]
!| (adaptive) run the whole test;
!| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]]
!| and
!| [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]
!| run one given assignment, for a caller who draws the assignments itself.
!| `calc_js_comp_test_permutation_bounds` and `calc_adaptive_js_comp_test_permutation_bounds` size
!| their work arrays.
!|
!| The pipeline is split over nine modules -- bin counts, histograms, the criteria both parameter
!| searches share, the fixed-k search, the adaptive search, the final run, the k-NN and the adaptive
!| study-label Permutation-Test and their sizing routines -- which this one gathers, so a caller reaches all of
!| them through it.
!|
!| Generated from [[tox_data_integration_js_comp_test_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_js_comp_test
    use tox_data_integration_js_comp_test_adaptive_search
    use tox_data_integration_js_comp_test_binning
    use tox_data_integration_js_comp_test_fixed_k_search
    use tox_data_integration_js_comp_test_permutation_adaptive
    use tox_data_integration_js_comp_test_permutation_bounds
    use tox_data_integration_js_comp_test_permutation_fixed_k
    use tox_data_integration_js_comp_test_run
    use tox_data_integration_js_comp_test_search
end module tox_data_integration_js_comp_test
