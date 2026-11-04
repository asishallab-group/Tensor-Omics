r"""tox_data_integration_stats_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Bootstrap

The bootstrap
:func:`tensor_omics.bootstrap_histogram` gives each
study's weighted global JSD an interval, by drawing each study's histograms again with
replacement from the study's own counts: a robustness interval, how much the observed value
would move if each study's residuals were drawn again. `calc_js_comp_test_n_top_k_jsds` sizes
its top-/bottom-k heaps for a caller that allocates its own. The significance test of the
JS-Comp-Test is the study-label Permutation-Test,
:func:`tensor_omics.run_js_comp_test_permutation`
and its adaptive counterpart.

Python binding, generated from tox_data_integration_stats_impl. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.calc_js_comp_test_n_top_k_jsds_c.restype = None
_lib.calc_js_comp_test_n_top_k_jsds_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_stats_impl::calc_js_comp_test_n_top_k_jsds`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_bootstraps = _contract.integer(n_bootstraps, "n_bootstraps")
    two_sided_bootstrapping_significance_level = _contract.real(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")

    # outputs and work arrays, which the caller never sees
    n_top_k = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.calc_js_comp_test_n_top_k_jsds_c(
        _ctypes.byref(_ctypes.c_int(n_bootstraps)),
        _ctypes.byref(_ctypes.c_double(two_sided_bootstrapping_significance_level)),
        _ctypes.byref(n_top_k),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CALC_JS_COMP_TEST_N_TOP_K_JSDS_ARGUMENTS)

    return n_top_k.value
