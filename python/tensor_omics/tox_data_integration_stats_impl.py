r"""tox_data_integration_stats_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test and Bootstrap

The permutation test of the JS-Comp-Test:
:func:`tensor_omics.gjct_permutation_test` estimates an
empirical p-value for each of K studies' weighted global JSD against the consensus pmf, by
resampling every reference point's pooled histogram counts into the studies without
replacement. `p_values` is the primary result; `p_values_observed_consensus` uses the older
null (a fixed observed consensus) and is kept for comparison.

The bootstrap
:func:`tensor_omics.bootstrap_histogram` gives each
study's weighted global JSD a confidence interval, by drawing each study's histograms with
replacement from the pooled consensus counts; `calc_js_comp_test_n_top_k_jsds` sizes its top-/bottom-k
heaps for a caller that allocates its own.

Python binding, generated from tox_data_integration_stats_impl. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.calc_js_comp_test_n_top_k_jsds_c.restype = None
_lib.calc_js_comp_test_n_top_k_jsds_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_JS_COMP_TEST_N_TOP_K_JSDS_ARGUMENTS = ("n_bootstraps", "two_sided_bootstrapping_significance_level", "n_top_k", "ierr",)

def calc_js_comp_test_n_top_k_jsds(
        n_bootstraps,
        two_sided_bootstrapping_significance_level=2.5,
):
    r"""Recommend the bootstrap top-/bottom-k heap size for a two-sided confidence interval

    Sizes the top-k/bottom-k heaps
    :func:`tensor_omics.bootstrap_histogram` uses
    for the interval's endpoints: `n_top_k = max(1, floor(level/100 * n_bootstraps))`, e.g. 25 for
    the default 2.5 and `n_bootstraps=1000`.

    Rejects an out-of-range argument with invalid input and a non-finite significance level with
    `ERR_NAN_INF`; `n_top_k` is then left undefined.

    Parameters
    ----------
    n_bootstraps : int
        Number of bootstrap resamples that will be performed
        The minimum valid value is `1`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Two-sided significance level, as a percentage reserved at each end of the bootstrap
        distribution (e.g. 2.5 for a 95% two-sided interval)
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.

    Returns
    -------
    n_top_k : int
        Number of elements to keep at each end of the bootstrap distribution, at least 1

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_stats_impl::calc_js_comp_test_n_top_k_jsds`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    n_top_k = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_js_comp_test_n_top_k_jsds_c(
        ctypes.byref(ctypes.c_int(n_bootstraps)),
        ctypes.byref(ctypes.c_double(two_sided_bootstrapping_significance_level)),
        ctypes.byref(n_top_k),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_JS_COMP_TEST_N_TOP_K_JSDS_ARGUMENTS)

    return n_top_k.value
