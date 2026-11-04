"""
Python test suite for the data integration stats functions in tensoromics: the heap-size routine
of the JS-Comp-Test's bootstrap. The bootstrap itself (`bootstrap_histogram`) is tested with the
rest of the JS-Comp-Test in mod_test_data_integration_js_comp_test.py; the numbers are checked in
the Fortran suite (test/mod_test_data_integration_js_comp_test.F90).
"""
import numpy as np
import sys
import os

# Add parent directory to path to import tensor_omics
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from test_helpers import run_all_tests
from tensor_omics import calc_js_comp_test_n_top_k_jsds


def test_calc_js_comp_test_n_top_k_jsds_call():
    """The heap size is an int: max(1, floor(2.5/100 * 1000)) = 25 at the defaults."""
    n_top_k = calc_js_comp_test_n_top_k_jsds(1000)
    assert isinstance(n_top_k, (int, np.integer)) and n_top_k == 25, f"expected 25, got {n_top_k}"


if __name__ == '__main__':
    run_all_tests(globals().values())
