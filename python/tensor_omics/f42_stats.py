"""f42_stats

Descriptive statistics: percentiles, empirical distribution functions, and 2-D LOESS smoothing.

Python binding, generated from f42_stats. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.loess_smooth_2d_c.restype = None
_lib.loess_smooth_2d_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LOESS_SMOOTH_2D_ARGUMENTS = ("n_total", "n_target", "x_ref", "y_ref", "indices_used", "n_used", "x_query", "kernel_sigma", "kernel_cutoff", "y_out", "ierr",)

_lib.compute_edf_expert_c.restype = None
_lib.compute_edf_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_EDF_EXPERT_ARGUMENTS = ("values", "n_values", "perm", "unique_values", "cdf_values", "n_unique", "ierr",)

_lib.compute_edf_c.restype = None
_lib.compute_edf_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_EDF_ARGUMENTS = ("values", "n_values", "unique_values", "cdf_values", "n_unique", "ierr",)

_lib.compute_scaled_distance_quantile_c.restype = None
_lib.compute_scaled_distance_quantile_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_SCALED_DISTANCE_QUANTILE_ARGUMENTS = ("n_genes", "rdi", "sorted_rdi", "perm", "quantile", "c_const",)

def loess_smooth_2d(
        x_ref,
        y_ref,
        indices_used,
        x_query,
        kernel_sigma,
        kernel_cutoff,
):
    r"""Performs LOESS smoothing on a set of data points

    Parameters
    ----------
    x_ref : np.ndarray[np.float64] of shape (n_total,)
        Reference x-coordinates.
    y_ref : np.ndarray[np.float64] of shape (n_total,)
        Reference y-coordinates (length n_total).
    indices_used : np.ndarray[np.int32] of shape (n_used,)
        Indices of reference points used for smoothing (only valid indices).
    x_query : np.ndarray[np.float64] of shape (n_target,)
        Target x-coordinates to smooth.
    kernel_sigma : float
        Bandwidth parameter for the kernel.
    kernel_cutoff : float
        Cutoff for the kernel, not used if zero

    Returns
    -------
    y_out : np.ndarray[np.float64] of shape (n_target,)
        Output smoothed values (length n_target).

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_stats::loess_smooth_2d`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        x_ref = np.ascontiguousarray(x_ref, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'x_ref' must be an array of np.float64: {error}") from None
    if x_ref.ndim != 1:
        raise ValueError(f"'x_ref' must have 1 dimension, but has {x_ref.ndim}")
    try:
        y_ref = np.ascontiguousarray(y_ref, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'y_ref' must be an array of np.float64: {error}") from None
    if y_ref.ndim != 1:
        raise ValueError(f"'y_ref' must have 1 dimension, but has {y_ref.ndim}")
    try:
        indices_used = np.ascontiguousarray(indices_used, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'indices_used' must be an array of np.int32: {error}") from None
    if indices_used.ndim != 1:
        raise ValueError(f"'indices_used' must have 1 dimension, but has {indices_used.ndim}")
    try:
        x_query = np.ascontiguousarray(x_query, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'x_query' must be an array of np.float64: {error}") from None
    if x_query.ndim != 1:
        raise ValueError(f"'x_query' must have 1 dimension, but has {x_query.ndim}")

    # what the inputs already say, rather than asking for it again
    n_total = x_ref.shape[0]
    n_target = x_query.shape[0]
    n_used = indices_used.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if y_ref.shape[0] != n_total:
        raise ValueError(f"'y_ref' has {y_ref.shape[0]} along axis 0, but "
            f"'x_ref' implies n_total == {n_total}"
        )

    # outputs and work arrays, which the caller never sees
    y_out = np.empty((n_target,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.loess_smooth_2d_c(
        ctypes.byref(ctypes.c_int(n_total)),
        ctypes.byref(ctypes.c_int(n_target)),
        x_ref,
        y_ref,
        indices_used,
        ctypes.byref(ctypes.c_int(n_used)),
        x_query,
        ctypes.byref(ctypes.c_double(kernel_sigma)),
        ctypes.byref(ctypes.c_double(kernel_cutoff)),
        y_out,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _LOESS_SMOOTH_2D_ARGUMENTS)

    return y_out

def compute_edf_expert(
        values,
        perm,
):
    r"""Compute the Empirical Distribution Function (EDF) from pre-sorted permutation

    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (n_values,)
        Array of observed data values (e.g., contributions or spikes).
    perm : np.ndarray[np.int32] of shape (n_values,)
        Pre-sorted permutation indices (must be sorted by values[perm]).

    Returns
    -------
    dict
        with keys:

        unique_values : np.ndarray[np.float64] of shape (n_values,)
            Sorted unique data values.
            The first `n_unique` elements will hold the results.
        cdf_values : np.ndarray[np.float64] of shape (n_values,)
            Corresponding cumulative frequencies between 0 and 1.
            The first `n_unique` elements will hold the results.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_stats::compute_edf`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        values = np.ascontiguousarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'values' must be an array of np.float64: {error}") from None
    if values.ndim != 1:
        raise ValueError(f"'values' must have 1 dimension, but has {values.ndim}")
    try:
        perm = np.ascontiguousarray(perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'perm' must be an array of np.int32: {error}") from None
    if perm.ndim != 1:
        raise ValueError(f"'perm' must have 1 dimension, but has {perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_values = values.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if perm.shape[0] != n_values:
        raise ValueError(f"'perm' has {perm.shape[0]} along axis 0, but "
            f"'values' implies n_values == {n_values}"
        )

    # outputs and work arrays, which the caller never sees
    unique_values = np.empty((n_values,), dtype=np.float64, order='C')
    cdf_values = np.empty((n_values,), dtype=np.float64, order='C')
    n_unique = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.compute_edf_expert_c(
        values,
        ctypes.byref(ctypes.c_int(n_values)),
        perm,
        unique_values,
        cdf_values,
        ctypes.byref(n_unique),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_EDF_EXPERT_ARGUMENTS)

    return {
        "unique_values": unique_values[..., :n_unique.value],
        "cdf_values": cdf_values[..., :n_unique.value],
    }

def compute_edf(
        values,
):
    r"""Sorts the values and computes the Empirical Distribution Function (EDF)

    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (n_values,)
        Array of observed data values (e.g., contributions or spikes).

    Returns
    -------
    dict
        with keys:

        unique_values : np.ndarray[np.float64] of shape (n_values,)
            Sorted unique data values.
            The first `n_unique` elements will hold the results.
        cdf_values : np.ndarray[np.float64] of shape (n_values,)
            Corresponding cumulative frequencies between 0 and 1.
            The first `n_unique` elements will hold the results.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_stats::compute_edf_alloc`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        values = np.ascontiguousarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'values' must be an array of np.float64: {error}") from None
    if values.ndim != 1:
        raise ValueError(f"'values' must have 1 dimension, but has {values.ndim}")

    # what the inputs already say, rather than asking for it again
    n_values = values.shape[0]

    # outputs and work arrays, which the caller never sees
    unique_values = np.empty((n_values,), dtype=np.float64, order='C')
    cdf_values = np.empty((n_values,), dtype=np.float64, order='C')
    n_unique = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.compute_edf_c(
        values,
        ctypes.byref(ctypes.c_int(n_values)),
        unique_values,
        cdf_values,
        ctypes.byref(n_unique),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_EDF_ARGUMENTS)

    return {
        "unique_values": unique_values[..., :n_unique.value],
        "cdf_values": cdf_values[..., :n_unique.value],
    }

def compute_scaled_distance_quantile(
        rdi,
        sorted_rdi,
        perm,
        c_const,
):
    r"""Calculate the empirical quantile (effect-size measure) of scaled expression distances (RDI)

    Parameters
    ----------
    rdi : np.ndarray[np.float64] of shape (n_genes,)
        empirical distribution D
    sorted_rdi : np.ndarray[np.float64] of shape (n_genes,)
        empirical distribution D with non negative values
    perm : np.ndarray[np.int32] of shape (n_genes,)
        Permutation array with sorted indices for sorted_rdi
    c_const : float
        Constant used in the computation, typically 1

    Returns
    -------
    quantile : np.ndarray[np.float64] of shape (n_genes,)
        Output array to store the computed quantile for each gene.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_stats::compute_scaled_distance_quantile`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        rdi = np.ascontiguousarray(rdi, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'rdi' must be an array of np.float64: {error}") from None
    if rdi.ndim != 1:
        raise ValueError(f"'rdi' must have 1 dimension, but has {rdi.ndim}")
    try:
        sorted_rdi = np.ascontiguousarray(sorted_rdi, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'sorted_rdi' must be an array of np.float64: {error}") from None
    if sorted_rdi.ndim != 1:
        raise ValueError(f"'sorted_rdi' must have 1 dimension, but has {sorted_rdi.ndim}")
    try:
        perm = np.ascontiguousarray(perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'perm' must be an array of np.int32: {error}") from None
    if perm.ndim != 1:
        raise ValueError(f"'perm' must have 1 dimension, but has {perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = rdi.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if sorted_rdi.shape[0] != n_genes:
        raise ValueError(f"'sorted_rdi' has {sorted_rdi.shape[0]} along axis 0, but "
            f"'rdi' implies n_genes == {n_genes}"
        )
    if perm.shape[0] != n_genes:
        raise ValueError(f"'perm' has {perm.shape[0]} along axis 0, but "
            f"'rdi' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    quantile = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_scaled_distance_quantile_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        rdi,
        sorted_rdi,
        perm,
        quantile,
        ctypes.byref(ctypes.c_double(c_const)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_SCALED_DISTANCE_QUANTILE_ARGUMENTS)

    return quantile
