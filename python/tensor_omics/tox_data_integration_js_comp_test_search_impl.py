r"""tox_data_integration_js_comp_test_search_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search

The data-driven choice of the JS-Comp-Test's neighborhood parameters: a candidate is admissible
when its neighborhoods overlap enough and every consensus bin holds enough residuals, and the
search stops at the first admissible candidate whose bootstrapped JSDs have reached a plateau.

Fixed-k neighborhoods:
:func:`tensor_omics.run_js_comp_test_parameter_search`
searches the `(n_points, n_neighbors)` candidates of
:func:`tensor_omics.generate_js_comp_test_candidates`.
Adaptive neighborhoods:
:func:`tensor_omics.run_js_comp_test_adaptive_parameter_search`
searches the `(k_start, k_step, k_max)` candidates of
:func:`tensor_omics.generate_adaptive_js_comp_test_candidates`.
The building blocks are published too: the gate
:func:`tensor_omics.check_mean_pmf_min_counts`
and the plateau checks
:func:`tensor_omics.check_plateau_condition` and
:func:`tensor_omics.check_effect_size_plateau_condition`.
`calc_js_comp_test_candidate_bounds` and `calc_adaptive_js_comp_test_bounds` size the work
arrays for a caller that allocates its own.

Python binding, generated from tox_data_integration_js_comp_test_search_impl. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.calc_js_comp_test_candidate_bounds_c.restype = None
_lib.calc_js_comp_test_candidate_bounds_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_JS_COMP_TEST_CANDIDATE_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "max_n_points_candidate", "max_n_neighbors_candidate", "ierr",)

_lib.calc_adaptive_js_comp_test_bounds_c.restype = None
_lib.calc_adaptive_js_comp_test_bounds_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_ADAPTIVE_JS_COMP_TEST_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "n_studies", "max_n_points_candidate", "ierr",)

def calc_js_comp_test_candidate_bounds(
        max_n_genes_all_studies,
):
    r"""Recommend upper bounds for the js-comp-test candidate-grid work arrays

    For the candidate grid of
    :func:`tensor_omics.generate_js_comp_test_candidates`:
    `max_n_points_candidate` is exactly its largest `n_points` candidate, and
    `max_n_neighbors_candidate` a safe (not necessarily tight) bound on its largest `n_neighbors`.

    Rejects a non-positive `max_n_genes_all_studies` with invalid input; both outputs are then
    left undefined.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.

    Returns
    -------
    dict
        with keys:

        max_n_points_candidate : int
            Exact upper bound on the `n_points` candidate the grid ever produces
        max_n_neighbors_candidate : int
            Safe upper bound on the `n_neighbors` candidate the grid ever produces

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_search_impl::calc_js_comp_test_candidate_bounds`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    max_n_points_candidate = ctypes.c_int(0)
    max_n_neighbors_candidate = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_js_comp_test_candidate_bounds_c(
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(max_n_points_candidate),
        ctypes.byref(max_n_neighbors_candidate),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_JS_COMP_TEST_CANDIDATE_BOUNDS_ARGUMENTS)

    return {
        "max_n_points_candidate": max_n_points_candidate.value,
        "max_n_neighbors_candidate": max_n_neighbors_candidate.value,
    }

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
    `huge(1)/4 = 536870911` (so that its `k_max = 4*k_start_1` fits); the output is then
    left undefined.

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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_search_impl::calc_adaptive_js_comp_test_bounds`, whose argument names are
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
