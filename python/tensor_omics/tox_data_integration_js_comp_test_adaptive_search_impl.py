r"""tox_data_integration_js_comp_test_adaptive_search_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Adaptive Parameter Search

The data-driven choice of the JS-Comp-Test's adaptive neighborhoods, each grown from `k_start`
pooled entries in rounds of `k_step` up to `k_max` while the residual dispersion stays stable.
:func:`tensor_omics.run_js_comp_test_adaptive_parameter_search`
walks the `(k_start, k_step, k_max)` candidates of
:func:`tensor_omics.generate_adaptive_js_comp_test_candidates`
and stops at the first admissible candidate whose bootstrapped JSDs have reached a plateau, by
the criteria of
:func:`tensor_omics.check_plateau_condition` and
:func:`tensor_omics.check_effect_size_plateau_condition`;
every candidate it tries is logged with an `ADAPTIVE_CANDIDATE_*` status.
:func:`tensor_omics.construct_adaptive_neighborhoods`
builds the neighborhoods for the chosen knobs. `calc_adaptive_js_comp_test_bounds` sizes the
work arrays for a caller that allocates its own.

Python binding, generated from tox_data_integration_js_comp_test_adaptive_search_impl. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.calc_adaptive_js_comp_test_bounds_c.restype = None
_lib.calc_adaptive_js_comp_test_bounds_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_ADAPTIVE_JS_COMP_TEST_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "n_studies", "max_n_points_candidate", "ierr",)

def calc_adaptive_js_comp_test_bounds(
        max_n_genes_all_studies,
        n_studies,
):
    r"""Recommend a reference-point capacity for the adaptive js-comp-test candidate sequence

    Sizes the per-point work arrays of the adaptive parameter search over the candidates of
    :func:`tensor_omics.generate_adaptive_js_comp_test_candidates`.
    A candidate's reference-point count is only known after its neighborhoods are grown, so this
    is a **practical capacity, not a proven bound**:
    `max_n_points_candidate = min(N, ceiling(2 * N / k_start_last) + 1)`, with the factor 2 being
    ``ADAPTIVE_POINT_CAPACITY_FACTOR``,
    `N = max_n_genes_all_studies * n_studies` the padded pool size and `k_start_last` the smallest
    candidate `k_start`. The search rejects a candidate that exceeds it with a capacity status.

    Rejects, with invalid input naming neither argument when the bound is joint, a
    `max_n_genes_all_studies` or `n_studies` below 1, a product `N` beyond the 32-bit range, or a
    first candidate `k_start_1 = max(10, 2*n_studies, ceiling(0.02 * N))` above
    `floor(huge(1)/ADAPTIVE_K_MAX_FACTOR)` (``ADAPTIVE_K_MAX_FACTOR``), so that its
    `k_max = ceiling(ADAPTIVE_K_MAX_FACTOR*k_start_1)` fits; the output is then left undefined.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.
    n_studies : int
        Number of studies
        The minimum valid value is `1`.

    Returns
    -------
    max_n_points_candidate : int
        Practical upper bound on the number of reference points any adaptive candidate is
        allowed to emerge with

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_adaptive_search_impl::calc_adaptive_js_comp_test_bounds`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    max_n_points_candidate = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_adaptive_js_comp_test_bounds_c(
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(max_n_points_candidate),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_ADAPTIVE_JS_COMP_TEST_BOUNDS_ARGUMENTS)

    return max_n_points_candidate.value
