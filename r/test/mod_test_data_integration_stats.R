source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# The data integration stats functions: the heap-size routine of the JS-Comp-Test's bootstrap. The
# bootstrap itself (bootstrap_histogram) is tested with the rest of the JS-Comp-Test in
# mod_test_data_integration_js_comp_test.R; the numbers are checked in the Fortran suite
# (test/mod_test_data_integration_js_comp_test.F90).

test_calc_js_comp_test_n_top_k_jsds_call <- function() {
  # max(1, floor(2.5/100 * 1000)) = 25 at the defaults.
  assert_true(calc_js_comp_test_n_top_k_jsds(1000L) == 25, "expected 25")
}

run_all_tests()
