r"""tox_data_integration_js_comp_test_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search

The data-driven `(n_points, n_neighbors)` parameter-stabilization search this pipeline runs
before the JSD-Comp-Test proper (Issue #126): a GAMMA-decay candidate grid
(:func:`tensor_omics.generate_js_comp_test_candidates`,
each candidate's histogram bin count from
:func:`tensor_omics.estimate_bin_count`), two
admissibility gates a candidate must pass before it is bootstrapped
(:func:`tensor_omics.check_neighborhood_overlaps`,
:func:`tensor_omics.check_mean_pmf_min_counts`),
and the plateau check that decides when the search has converged
(:func:`tensor_omics.check_plateau_condition`).
`calc_js_comp_test_candidate_bounds` sizes the candidate-grid work arrays for a caller that
allocates its own. Once a candidate has passed both gates, its bootstrap confidence interval
is resampled from the pooled consensus histogram by
:func:`tensor_omics.bootstrap_histogram` (heap
size recommended by
:func:`tensor_omics.calc_js_comp_test_n_top_k_jsds`).
Later stages of the port add the K-study permutation test and the top-level orchestrator that
wire these building blocks together.

Python binding, generated from tox_data_integration_js_comp_test_impl. Do not edit.
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
_CALC_JS_COMP_TEST_CANDIDATE_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "max_n_points_candidate", "max_n_neighbors_candidate",)

_lib.calc_js_comp_test_n_top_k_jsds_c.restype = None
_lib.calc_js_comp_test_n_top_k_jsds_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_JS_COMP_TEST_N_TOP_K_JSDS_ARGUMENTS = ("n_bootstraps", "two_sided_bootstrapping_significance_level", "n_top_k",)

def calc_js_comp_test_candidate_bounds(
        max_n_genes_all_studies,
):
    r"""Recommend upper bounds for the js-comp-test candidate-grid work arrays

    Closed-form from the GAMMA-decay candidate-grid formula
    :func:`tensor_omics.generate_js_comp_test_candidates`
    uses: `max_n_points_candidate` is exactly the grid's first (largest) `n_points` candidate.
    `max_n_neighbors_candidate` is a safe, not necessarily tight, upper bound on the grid's
    largest `n_neighbors` candidate -- reached with the smallest `KX_FACTORS` entry at the
    smallest `n_points_high` the grid loop ever uses, which by construction never drops below
    `n_points_low`.

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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_impl::calc_js_comp_test_candidate_bounds`, whose argument names are
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

def calc_js_comp_test_n_top_k_jsds(
        n_bootstraps,
        two_sided_bootstrapping_significance_level=2.5,
):
    r"""Recommend the bootstrap top-/bottom-k heap size for a two-sided confidence interval

    Ported from 125-stabilize-jscomp's inline `n_bootstrapping_top_k_jsds` computation in
    `determine_js_comp_test_n_points_n_neighbors_alloc`: sizes the top-k/bottom-k heaps
    :func:`tensor_omics.bootstrap_histogram` uses
    to track a two-sided bootstrap confidence interval's endpoints. E.g. for a 95% two-sided
    interval (2.5% reserved at each end, the default) with `n_bootstraps=1000`,
    `n_top_k = max(1, floor(0.025*1000)) = 25`.

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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_impl::calc_js_comp_test_n_top_k_jsds`, whose argument names are
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
