r"""tox_data_integration_preprocessing_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing

The step that turns expression vectors into the neighborhood residuals the rest of the test
consumes: gene-wise means, the signed deviation of each replicate from them, and the
neighborhoods of reference points those residuals are grouped into so the comparison is
conditioned on expression level rather than pooled across it.

Besides neighborhoods of a fixed size, the adaptive construction
:func:`tensor_omics.construct_adaptive_neighborhoods`
grows each reference point's neighborhood while its residual dispersion stays stable, and
:func:`tensor_omics.check_neighborhood_overlaps`
checks that consecutive neighborhoods overlap enough.

`calc_neighborhood_size` sizes a neighborhood for a caller that allocates its own, and
`gather_pooled_neighborhood_residuals` gathers a neighborhood's residuals of all studies.

Python binding, generated from tox_data_integration_preprocessing_impl. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.calc_neighborhood_size_c.restype = None
_lib.calc_neighborhood_size_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_NEIGHBORHOOD_SIZE_ARGUMENTS = ("n_pool", "n_points", "n_genes_S", "mean_S", "desired_size",)
#: For a derived argument, the one the caller passed it in
_CALC_NEIGHBORHOOD_SIZE_ARGUMENT_SOURCES = (None, None, "mean_S", None, None,)

_lib.gather_pooled_neighborhood_residuals_c.restype = None
_lib.gather_pooled_neighborhood_residuals_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GATHER_POOLED_NEIGHBORHOOD_RESIDUALS_ARGUMENTS = ("residuals", "max_n_reps_all_studies", "max_n_genes_all_studies", "n_neighbors", "n_studies", "neighborhood_indices_point", "pooled_residuals", "ierr",)
#: For a derived argument, the one the caller passed it in
_GATHER_POOLED_NEIGHBORHOOD_RESIDUALS_ARGUMENT_SOURCES = (None, "residuals", "residuals", "neighborhood_indices_point", "residuals", None, None, None,)

def calc_neighborhood_size(
        n_pool,
        n_points,
        mean_S,
        desired_size=1000,
):
    r"""Calculate the number of neighbors to be used for constructing neighborhoods

    The `desired_size` works as upper limit, as the actual neighborhood size might be lower
    due to few genes with non-NaN mean.

    Parameters
    ----------
    n_pool : int
        Total number of pooled mean-expression values across both studies
    n_points : int
        Number of reference points
    mean_S : np.ndarray[np.float64] of shape (n_genes_S,)
        Per-gene mean expression values
        NaN is permitted for this value.
    desired_size : int, optional, default 1000
        Optional desired neighborhood size
        The default value is `1000`.

    Returns
    -------
    n_neighbors : int
        Calculated neighborhood size

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing_impl::calc_neighborhood_size`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        mean_S = np.ascontiguousarray(mean_S, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_S' must be an array of np.float64: {error}") from None
    if mean_S.ndim != 1:
        raise ValueError(f"'mean_S' must have 1 dimension, but has {mean_S.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes_S = mean_S.shape[0]

    # outputs and work arrays, which the caller never sees
    n_neighbors = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_neighborhood_size_c(
        ctypes.byref(ctypes.c_int(n_pool)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_genes_S)),
        mean_S,
        ctypes.byref(ctypes.c_int(desired_size)),
        ctypes.byref(n_neighbors),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_NEIGHBORHOOD_SIZE_ARGUMENTS, _CALC_NEIGHBORHOOD_SIZE_ARGUMENT_SOURCES)

    return n_neighbors.value

def gather_pooled_neighborhood_residuals(
        residuals,
        neighborhood_indices_point,
):
    r"""Pool one reference point's residuals across every neighbor and every study

    Gathers one reference point's residuals from every neighbor gene of every study into one flat
    array, given that point's per-study neighbor gene indices (one column of the
    `neighborhood_indices_all_studies(n_neighbors, n_points, n_studies)` of
    :func:`tensor_omics.construct_neighborhoods_ranged`).
    Study `s` fills pool block `(s - 1)*max_n_reps_all_studies*n_neighbors + 1` to
    `s*max_n_reps_all_studies*n_neighbors`. This is the pooling
    :func:`tensor_omics.run_js_comp_test` and
    :func:`tensor_omics.run_js_comp_test_parameter_search`
    feed to the bin-count search, so a caller can pass the same input to
    :func:`tensor_omics.determine_bin_count_occupancy_exhaustive`
    on real data.

    Parameters
    ----------
    residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Matrix of signed residuals per study, NaN allowed for missing values; not validated
    neighborhood_indices_point : np.ndarray[np.int32] of shape (n_neighbors, n_studies,), column-major (order='F')
        Gene indices of one reference point's neighborhood, per study; each in
        [1, max_n_genes_all_studies]

    Returns
    -------
    pooled_residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies*n_neighbors*n_studies,), read-only
        This reference point's pooled residuals, laid out as a
        (max_n_reps_all_studies, n_neighbors, n_studies) array
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing_impl::gather_pooled_neighborhood_residuals`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        residuals = np.asfortranarray(residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'residuals' must be an array of np.float64: {error}") from None
    if residuals.ndim != 3:
        raise ValueError(f"'residuals' must have 3 dimensions, but has {residuals.ndim}")
    try:
        neighborhood_indices_point = np.asfortranarray(neighborhood_indices_point, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_indices_point' must be an array of np.int32: {error}") from None
    if neighborhood_indices_point.ndim != 2:
        raise ValueError(f"'neighborhood_indices_point' must have 2 dimensions, but has {neighborhood_indices_point.ndim}")

    # what the inputs already say, rather than asking for it again
    max_n_reps_all_studies = residuals.shape[0]
    max_n_genes_all_studies = residuals.shape[1]
    n_neighbors = neighborhood_indices_point.shape[0]
    n_studies = residuals.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if neighborhood_indices_point.shape[1] != n_studies:
        raise ValueError(f"'neighborhood_indices_point' has {neighborhood_indices_point.shape[1]} along axis 1, but "
            f"'residuals' implies n_studies == {n_studies}"
        )

    # outputs and work arrays, which the caller never sees
    pooled_residuals = np.empty((max_n_reps_all_studies*n_neighbors*n_studies,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.gather_pooled_neighborhood_residuals_c(
        residuals,
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_studies)),
        neighborhood_indices_point,
        pooled_residuals,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GATHER_POOLED_NEIGHBORHOOD_RESIDUALS_ARGUMENTS, _GATHER_POOLED_NEIGHBORHOOD_RESIDUALS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    pooled_residuals.flags.writeable = False

    return pooled_residuals
