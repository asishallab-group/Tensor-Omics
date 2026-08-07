"""tox_shape_truthful_clustering_observable

Generated from the kernel; do not edit -- regenerate instead.

Python binding, generated from tox_shape_truthful_clustering_observable. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.normal_error_c.restype = None
_lib.normal_error_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMAL_ERROR_ARGUMENTS = ("d", "eigenvalues", "n_dimensions", "normal_error_value", "ierr",)
#: For a derived argument, the one the caller passed it in
_NORMAL_ERROR_ARGUMENT_SOURCES = (None, None, "eigenvalues", None, None,)

_lib.tangent_scales_c.restype = None
_lib.tangent_scales_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_TANGENT_SCALES_ARGUMENTS = ("d", "eigenvalues", "n_dimensions", "tangent_scales_value", "ierr",)
#: For a derived argument, the one the caller passed it in
_TANGENT_SCALES_ARGUMENT_SOURCES = ("tangent_scales_value", None, "eigenvalues", None, None,)

_lib.observable_c.restype = None
_lib.observable_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_OBSERVABLE_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "member_selection_mask", "n_selected_member", "U", "eigenvalues", "mu", "d", "G", "normal_error_value", "tangent_scales_value", "ierr",)
#: For a derived argument, the one the caller passed it in
_OBSERVABLE_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, "member_selection_mask", None, None, None, None, None, None, None, None,)

def normal_error(
        d,
        eigenvalues,
):
    r"""Mean squared residual of an ensemble's members off its tangent subspace

    Parameters
    ----------
    d : int
        Intrinsic (tangent) dimension of the ensemble
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    eigenvalues : np.ndarray[np.float64] of shape (n_dimensions,)
        Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
        The minimum valid value is `0.0`.

    Returns
    -------
    normal_error_value : float
        Mean squared residual off the d-dimensional tangent subspace

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_observable::normal_error`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        eigenvalues = np.ascontiguousarray(eigenvalues, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'eigenvalues' must be an array of np.float64: {error}") from None
    if eigenvalues.ndim != 1:
        raise ValueError(f"'eigenvalues' must have 1 dimension, but has {eigenvalues.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = eigenvalues.shape[0]

    # outputs and work arrays, which the caller never sees
    normal_error_value = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.normal_error_c(
        ctypes.byref(ctypes.c_int(d)),
        eigenvalues,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(normal_error_value),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _NORMAL_ERROR_ARGUMENTS, _NORMAL_ERROR_ARGUMENT_SOURCES)

    return normal_error_value.value

def tangent_scales(
        d,
        eigenvalues,
):
    r"""Extent along each tangent direction of an ensemble's tangent subspace

    Parameters
    ----------
    d : int
        Intrinsic (tangent) dimension of the ensemble
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    eigenvalues : np.ndarray[np.float64] of shape (n_dimensions,)
        Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
        The minimum valid value is `0.0`.

    Returns
    -------
    tangent_scales_value : np.ndarray[np.float64] of shape (d,), read-only
        Extent along each of the d tangent directions
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_observable::tangent_scales`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        eigenvalues = np.ascontiguousarray(eigenvalues, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'eigenvalues' must be an array of np.float64: {error}") from None
    if eigenvalues.ndim != 1:
        raise ValueError(f"'eigenvalues' must have 1 dimension, but has {eigenvalues.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = eigenvalues.shape[0]

    # outputs and work arrays, which the caller never sees
    tangent_scales_value = np.empty((d,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.tangent_scales_c(
        ctypes.byref(ctypes.c_int(d)),
        eigenvalues,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        tangent_scales_value,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _TANGENT_SCALES_ARGUMENTS, _TANGENT_SCALES_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    tangent_scales_value.flags.writeable = False

    return tangent_scales_value

def observable(
        vectors,
        member_selection_mask,
):
    r"""The tuple (U, d, G, mu, normal_error, tangent_scales) for an ensemble

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix
    member_selection_mask : np.ndarray[np.bool_] of shape (n_vectors,)
        Ensemble membership over the full dataset

    Returns
    -------
    dict
        with keys:

        U : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions,), column-major (order='F'), read-only
            Tangent+normal basis, zero-padded beyond rank
            A result is a value; call `.copy()` to obtain a modifiable array.
        eigenvalues : np.ndarray[np.float64] of shape (n_dimensions,), read-only
            Covariance eigenvalues, descending, zero-padded beyond rank
            A result is a value; call `.copy()` to obtain a modifiable array.
        mu : np.ndarray[np.float64] of shape (n_dimensions,), read-only
            Ensemble center
            A result is a value; call `.copy()` to obtain a modifiable array.
        d : int
            Estimated intrinsic (tangent) dimension
        G : float
            Spectral gap at d
        normal_error_value : float
            Mean squared residual off the tangent subspace
        tangent_scales_value : np.ndarray[np.float64] of shape (n_dimensions,), read-only
            Extent along each tangent direction, zero-padded beyond d
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_observable::observable_alloc`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        vectors = np.asfortranarray(vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vectors' must be an array of np.float64: {error}") from None
    if vectors.ndim != 2:
        raise ValueError(f"'vectors' must have 2 dimensions, but has {vectors.ndim}")
    try:
        member_selection_mask = np.ascontiguousarray(member_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'member_selection_mask' must be an array of np.bool_: {error}") from None
    if member_selection_mask.ndim != 1:
        raise ValueError(f"'member_selection_mask' must have 1 dimension, but has {member_selection_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]
    n_selected_member = int(member_selection_mask.sum())

    # Fortran cannot check that shared extents agree; this can
    if member_selection_mask.shape[0] != n_vectors:
        raise ValueError(f"'member_selection_mask' has {member_selection_mask.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )

    # outputs and work arrays, which the caller never sees
    U = np.empty((n_dimensions, n_dimensions,), dtype=np.float64, order='F')
    eigenvalues = np.empty((n_dimensions,), dtype=np.float64, order='C')
    mu = np.empty((n_dimensions,), dtype=np.float64, order='C')
    d = ctypes.c_int(0)
    G = ctypes.c_double(0)
    normal_error_value = ctypes.c_double(0)
    tangent_scales_value = np.empty((n_dimensions,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.observable_c(
        vectors,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        member_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_member)),
        U,
        eigenvalues,
        mu,
        ctypes.byref(d),
        ctypes.byref(G),
        ctypes.byref(normal_error_value),
        tangent_scales_value,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _OBSERVABLE_ARGUMENTS, _OBSERVABLE_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    U.flags.writeable = False
    eigenvalues.flags.writeable = False
    mu.flags.writeable = False
    tangent_scales_value.flags.writeable = False

    return {
        "U": U,
        "eigenvalues": eigenvalues,
        "mu": mu,
        "d": d.value,
        "G": G.value,
        "normal_error_value": normal_error_value.value,
        "tangent_scales_value": tangent_scales_value,
    }
