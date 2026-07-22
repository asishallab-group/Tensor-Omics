"""Python interface to k-d tree spatial index over fixed-dimensional point sets.

Generated from f42_kd_tree. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.build_kd_index_c.restype = None
_lib.build_kd_index_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_KD_INDEX_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "ierr",)

_lib.build_spherical_kd_c.restype = None
_lib.build_spherical_kd_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_SPHERICAL_KD_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "ierr",)

def build_kd_index(
        points,
        dimension_order,
):
    r"""Build a k-d tree index using a stack-based, non-recursive approach

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Data points
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order (by variance)

    Returns
    -------
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        Output index array (k-d tree order)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::build_kd_index_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        points = np.asfortranarray(points, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'points' must be an array of np.float64: {error}") from None
    if points.ndim != 2:
        raise ValueError(f"'points' must have 2 dimensions, but has {points.ndim}")
    try:
        dimension_order = np.ascontiguousarray(dimension_order, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'dimension_order' must be an array of np.int32: {error}") from None
    if dimension_order.ndim != 1:
        raise ValueError(f"'dimension_order' must have 1 dimension, but has {dimension_order.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'points' implies n_dimensions == {n_dimensions}"
        )

    # outputs and work arrays, which the caller never sees
    kd_indices = np.empty((n_points,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.build_kd_index_c(
        points,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BUILD_KD_INDEX_ARGUMENTS)

    return kd_indices

def build_spherical_kd(
        points,
        dimension_order,
):
    r"""Build a k-d tree index over points on the unit sphere (unit vectors)

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Data points
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order (by variance)

    Returns
    -------
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        Output index array (k-d tree order)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::build_spherical_kd_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        points = np.asfortranarray(points, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'points' must be an array of np.float64: {error}") from None
    if points.ndim != 2:
        raise ValueError(f"'points' must have 2 dimensions, but has {points.ndim}")
    try:
        dimension_order = np.ascontiguousarray(dimension_order, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'dimension_order' must be an array of np.int32: {error}") from None
    if dimension_order.ndim != 1:
        raise ValueError(f"'dimension_order' must have 1 dimension, but has {dimension_order.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'points' implies n_dimensions == {n_dimensions}"
        )

    # outputs and work arrays, which the caller never sees
    kd_indices = np.empty((n_points,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.build_spherical_kd_c(
        points,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BUILD_SPHERICAL_KD_ARGUMENTS)

    return kd_indices
