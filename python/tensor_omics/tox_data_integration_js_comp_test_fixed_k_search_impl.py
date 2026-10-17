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
By default the grid holds only the large neighborhoods
(``'neighborhood_size_large'``);
`neighborhood_size_mode` selects the small ones
(``'neighborhood_size_small'``)
or both
(``'neighborhood_size_both'``).
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
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_JS_COMP_TEST_CANDIDATE_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "n_studies", "max_n_points_candidate", "max_n_neighbors_candidate", "grid_sample_size_mode", "grid_points_factor", "grid_min_points", "ierr",)

def calc_js_comp_test_candidate_bounds(
        max_n_genes_all_studies,
        n_studies,
        grid_sample_size_mode='grid_sample_size_per_study',
        grid_points_factor=4.0,
        grid_min_points=200,
):
    r"""Recommend upper bounds for the js-comp-test candidate-grid work arrays

    For the candidate grid of
    :func:`tensor_omics.generate_js_comp_test_candidates`:
    `max_n_points_candidate` is exactly its largest `n_points` candidate, and
    `max_n_neighbors_candidate` a safe (not necessarily tight) bound on its largest `n_neighbors`,
    in every `neighborhood_size_mode` (the bound uses the large-neighborhood factor `KX_FACTORS(1)`)
    and every `grid_floor_mode` (the floor candidate of the clamp sits at the floor this bound
    already uses). Pass the same `grid_sample_size_mode`, `grid_points_factor` and
    `grid_min_points` as the parameter search that is sized with the result.

    Rejects a non-positive `max_n_genes_all_studies` or `n_studies`, an unknown
    `grid_sample_size_mode`, a non-positive `grid_points_factor` and a `grid_min_points` outside
    `[1, MAX_POINTS]` with invalid input (a non-finite factor with `ERR_NAN_INF`); both outputs
    are then left undefined.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.
    n_studies : int
        Number of studies; used only by the pooled grid sample size
        The minimum valid value is `1`.
    grid_sample_size_mode : str, one of 'grid_sample_size_per_study' | 'grid_sample_size_pooled', optional, default 'grid_sample_size_per_study'
        Sample size `N` in the grid's first reference-point count; must match the parameter
        search's own. An experiment switch

        The default value is `'grid_sample_size_per_study'`.
    grid_points_factor : float, optional, default 4.0
        Factor in the grid's first reference-point count `ceiling(grid_points_factor * sqrt(N))`;
        must match the parameter search's own. An experiment switch
        The minimum valid value is `above(0.0)`.
        The default value is `4.0`.
    grid_min_points : int, optional, default 200
        The grid's floor; must match the parameter search's own. An experiment switch: 300
        (`MIN_POINTS`) is the floor before it existed
        The minimum valid value is `1`.
        The maximum valid value is `MAX_POINTS`.
        The default value is `200`.

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
    # accept anything array-like, converting only when C needs it
    grid_sample_size_mode = np.array([str(grid_sample_size_mode).lower().encode().ljust(26)], dtype="S26")

    # outputs and work arrays, which the caller never sees
    max_n_points_candidate = ctypes.c_int(0)
    max_n_neighbors_candidate = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_js_comp_test_candidate_bounds_c(
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(max_n_points_candidate),
        ctypes.byref(max_n_neighbors_candidate),
        grid_sample_size_mode,
        ctypes.byref(ctypes.c_double(grid_points_factor)),
        ctypes.byref(ctypes.c_int(grid_min_points)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_JS_COMP_TEST_CANDIDATE_BOUNDS_ARGUMENTS)

    return {
        "max_n_points_candidate": max_n_points_candidate.value,
        "max_n_neighbors_candidate": max_n_neighbors_candidate.value,
    }
