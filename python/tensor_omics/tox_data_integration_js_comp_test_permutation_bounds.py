r"""tox_data_integration_js_comp_test_permutation_bounds

# Work-array capacities of the JS-Comp-Test study-label Permutation-Test

The sizing routines of the Permutation-Test that permutes the samples' study labels
(``tox_data_integration_js_comp_test_permutation_fixed_k_impl`` and ``tox_data_integration_js_comp_test_permutation_adaptive_impl``):
:func:`tensor_omics.calc_js_comp_test_permutation_bounds`
for the k-NN routines and
:func:`tensor_omics.calc_adaptive_js_comp_test_permutation_bounds`
for the adaptive ones. The plain entry points of the Permutation-Test call them for you; a
caller who allocates the work arrays itself (the expert entry points) calls them first.

Python binding, generated from tox_data_integration_js_comp_test_permutation_bounds. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.calc_js_comp_test_permutation_bounds_c.restype = None
_lib.calc_js_comp_test_permutation_bounds_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_JS_COMP_TEST_PERMUTATION_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "n_points", "n_neighbors", "n_bootstraps", "two_sided_bootstrapping_significance_level", "grid_points_factor", "grid_min_points", "max_n_points_candidate", "max_n_neighbors_candidate", "n_bootstrapping_top_k_jsds", "ierr",)

_lib.calc_adaptive_js_comp_test_permutation_bounds_c.restype = None
_lib.calc_adaptive_js_comp_test_permutation_bounds_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_ADAPTIVE_JS_COMP_TEST_PERMUTATION_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "n_studies", "k_start", "n_bootstraps", "two_sided_bootstrapping_significance_level", "max_n_points_candidate", "n_bootstrapping_top_k_jsds", "ierr",)

def calc_js_comp_test_permutation_bounds(
        max_n_genes_all_studies,
        n_points,
        n_neighbors,
        n_bootstraps=1000,
        two_sided_bootstrapping_significance_level=2.5,
        grid_points_factor=4.0,
        grid_min_points=200,
):
    r"""Recommend the work-array capacities of the k-NN study-label Permutation-Test

    Sizes the per-point and per-neighbor work arrays of
    :func:`tensor_omics.run_js_comp_test_on_sample_assignment`
    and
    :func:`tensor_omics.run_js_comp_test_permutation`
    so that both the fixed `(n_points, n_neighbors)` run and a rerun of the parameter search fit:
    `max_n_points_candidate = max(n_points, P)` and `max_n_neighbors_candidate = max(n_neighbors,
    K)`, with `P` and `K` the bounds of
    :func:`tensor_omics.calc_js_comp_test_candidate_bounds`
    for the same `grid_points_factor` and `grid_min_points`.
    `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
    :func:`tensor_omics.calc_js_comp_test_n_top_k_jsds` computes it.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.
    n_points : int
        Number of reference points of the fixed-parameter run
        The minimum valid value is `1`.
    n_neighbors : int
        Number of neighbors per neighborhood of the fixed-parameter run
        The minimum valid value is `1`.
        The maximum valid value is `max_n_genes_all_studies`.
    n_bootstraps : int, optional, default 1000
        Number of bootstraps of the parameter search, when it is rerun
        The minimum valid value is `1`.
        The default value is `1000`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Two-sided significance level of the bootstrap interval, in percent per end
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.
    grid_points_factor : float, optional, default 4.0
        The rerun parameter search's `grid_points_factor`
        The minimum valid value is `above(0.0)`.
        The default value is `4.0`.
    grid_min_points : int, optional, default 200
        The rerun parameter search's `grid_min_points`
        The minimum valid value is `1`.
        The maximum valid value is `MAX_POINTS`.
        The default value is `200`.

    Returns
    -------
    dict
        with keys:

        max_n_points_candidate : int
            Capacity for reference points
        max_n_neighbors_candidate : int
            Capacity for neighbors per neighborhood
        n_bootstrapping_top_k_jsds : int
            Bootstrap heap size (values kept at each end of the bootstrap distribution)

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation_bounds::calc_js_comp_test_permutation_bounds`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    max_n_genes_all_studies = _contract.integer(max_n_genes_all_studies, "max_n_genes_all_studies")
    n_points = _contract.integer(n_points, "n_points")
    n_neighbors = _contract.integer(n_neighbors, "n_neighbors")
    n_bootstraps = _contract.integer(n_bootstraps, "n_bootstraps")
    two_sided_bootstrapping_significance_level = _contract.real(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    grid_points_factor = _contract.real(grid_points_factor, "grid_points_factor")
    grid_min_points = _contract.integer(grid_min_points, "grid_min_points")

    # outputs and work arrays, which the caller never sees
    max_n_points_candidate = _ctypes.c_int(0)
    max_n_neighbors_candidate = _ctypes.c_int(0)
    n_bootstrapping_top_k_jsds = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.calc_js_comp_test_permutation_bounds_c(
        _ctypes.byref(_ctypes.c_int(max_n_genes_all_studies)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_neighbors)),
        _ctypes.byref(_ctypes.c_int(n_bootstraps)),
        _ctypes.byref(_ctypes.c_double(two_sided_bootstrapping_significance_level)),
        _ctypes.byref(_ctypes.c_double(grid_points_factor)),
        _ctypes.byref(_ctypes.c_int(grid_min_points)),
        _ctypes.byref(max_n_points_candidate),
        _ctypes.byref(max_n_neighbors_candidate),
        _ctypes.byref(n_bootstrapping_top_k_jsds),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CALC_JS_COMP_TEST_PERMUTATION_BOUNDS_ARGUMENTS)

    return {
        "max_n_points_candidate": max_n_points_candidate.value,
        "max_n_neighbors_candidate": max_n_neighbors_candidate.value,
        "n_bootstrapping_top_k_jsds": n_bootstrapping_top_k_jsds.value,
    }

def calc_adaptive_js_comp_test_permutation_bounds(
        max_n_genes_all_studies,
        n_studies,
        k_start,
        n_bootstraps=1000,
        two_sided_bootstrapping_significance_level=2.5,
):
    r"""Recommend the work-array capacity of the adaptive study-label Permutation-Test

    Sizes the per-point work arrays of
    :func:`tensor_omics.run_js_comp_test_adaptive_on_sample_assignment`
    and
    :func:`tensor_omics.run_js_comp_test_adaptive_permutation`:
    `max_n_points_candidate = max(A, min(N, ceiling(F * N / k_start) + 1))`, with `N =
    max_n_genes_all_studies * n_studies` (the padded pool size), `F` the factor
    ``ADAPTIVE_POINT_CAPACITY_FACTOR``
    and `A` the capacity of the adaptive parameter search,
    :func:`tensor_omics.calc_adaptive_js_comp_test_bounds`.
    As there, this is a **practical capacity, not a proven bound**: how many reference points a
    construction builds is only known once it has run, and a shuffled data set whose
    construction builds more fails with invalid input (it is then counted as a failed
    permutation). `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
    :func:`tensor_omics.calc_js_comp_test_n_top_k_jsds` computes it.

    Reports invalid input, naming neither argument, for the joint bounds on
    `max_n_genes_all_studies` and `n_studies` that the adaptive search's capacity rejects (a pool
    size or a first candidate's `k_max` beyond the 32-bit range); the outputs are then left
    undefined.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.
    n_studies : int
        Number of studies
        The minimum valid value is `1`.
    k_start : int
        Pooled entries every neighborhood of the fixed-knob run takes unconditionally
        The minimum valid value is `1`.
    n_bootstraps : int, optional, default 1000
        Number of bootstraps of the parameter search, when it is rerun
        The minimum valid value is `1`.
        The default value is `1000`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Two-sided significance level of the bootstrap interval, in percent per end
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.

    Returns
    -------
    dict
        with keys:

        max_n_points_candidate : int
            Practical capacity for reference points
        n_bootstrapping_top_k_jsds : int
            Bootstrap heap size (values kept at each end of the bootstrap distribution)

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation_bounds::calc_adaptive_js_comp_test_permutation_bounds`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    max_n_genes_all_studies = _contract.integer(max_n_genes_all_studies, "max_n_genes_all_studies")
    n_studies = _contract.integer(n_studies, "n_studies")
    k_start = _contract.integer(k_start, "k_start")
    n_bootstraps = _contract.integer(n_bootstraps, "n_bootstraps")
    two_sided_bootstrapping_significance_level = _contract.real(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")

    # outputs and work arrays, which the caller never sees
    max_n_points_candidate = _ctypes.c_int(0)
    n_bootstrapping_top_k_jsds = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.calc_adaptive_js_comp_test_permutation_bounds_c(
        _ctypes.byref(_ctypes.c_int(max_n_genes_all_studies)),
        _ctypes.byref(_ctypes.c_int(n_studies)),
        _ctypes.byref(_ctypes.c_int(k_start)),
        _ctypes.byref(_ctypes.c_int(n_bootstraps)),
        _ctypes.byref(_ctypes.c_double(two_sided_bootstrapping_significance_level)),
        _ctypes.byref(max_n_points_candidate),
        _ctypes.byref(n_bootstrapping_top_k_jsds),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CALC_ADAPTIVE_JS_COMP_TEST_PERMUTATION_BOUNDS_ARGUMENTS)

    return {
        "max_n_points_candidate": max_n_points_candidate.value,
        "n_bootstrapping_top_k_jsds": n_bootstrapping_top_k_jsds.value,
    }
