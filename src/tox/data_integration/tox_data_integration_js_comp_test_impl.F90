#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Pipeline
!|
!| The global JS-Comp-Test of K studies, with a data-driven choice of its neighborhood parameters.
!|
!| Fixed-k neighborhoods (`n_points` reference points of `n_neighbors` genes each):
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
!| searches the `(n_points, n_neighbors)` candidates of
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]
!| for a bootstrap plateau, and
!| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]] runs the
!| final test (histograms, JSDs, weighting, permutation test) for the chosen pair.
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
!| Each reference point's histogram bin count is chosen by the occupancy search
!| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
!| (default) or by Stone's cross-validation
!| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]
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
!| The pipeline is split over six modules -- bin counts, histograms, the criteria both parameter
!| searches share, the fixed-k search, the adaptive search and the final run -- which this one
!| gathers, so a caller reaches all of them through it.
module tox_data_integration_js_comp_test_impl
    use tox_data_integration_js_comp_test_binning_impl
    use tox_data_integration_js_comp_test_histograms_impl
    use tox_data_integration_js_comp_test_search_impl
    use tox_data_integration_js_comp_test_fixed_k_search_impl
    use tox_data_integration_js_comp_test_adaptive_search_impl
    use tox_data_integration_js_comp_test_run_impl
end module tox_data_integration_js_comp_test_impl
