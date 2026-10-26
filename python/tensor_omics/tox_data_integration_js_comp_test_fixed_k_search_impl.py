r"""tox_data_integration_js_comp_test_fixed_k_search_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Fixed-k Parameter Search

The data-driven choice of the JS-Comp-Test's fixed-k neighborhoods, `n_points` reference
points of `n_neighbors` genes each.
:func:`tensor_omics.run_js_comp_test_parameter_search`
walks the `(n_points, n_neighbors)` candidate grid of
:func:`tensor_omics.generate_js_comp_test_candidates`
from finest to coarsest and stops at the first admissible candidate whose bootstrapped JSDs
have reached a plateau, by the criteria of
:func:`tensor_omics.check_plateau_condition` and
:func:`tensor_omics.check_effect_size_plateau_condition`.
`calc_js_comp_test_candidate_bounds` sizes the work arrays for a caller that allocates its own.

Python binding, generated from tox_data_integration_js_comp_test_fixed_k_search_impl. Do not edit.
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_fixed_k_search_impl::calc_js_comp_test_candidate_bounds`, whose argument names are
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
