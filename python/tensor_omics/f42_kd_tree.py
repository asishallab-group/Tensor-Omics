r"""f42_kd_tree

k-d tree spatial index over fixed-dimensional point sets.
Builds a k-d tree by recursively partitioning `kd_indices` around the median point along a
caller-supplied, cycling dimension order, using a stack-based (non-recursive) traversal so it
is safe to call from `pure` procedures. The tree is stored implicitly as an in-place-permuted
index array rather than as linked nodes.

Python binding, generated from f42_kd_tree. Do not edit.
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
#: For a derived argument, the one the caller passed it in
_BUILD_KD_INDEX_ARGUMENT_SOURCES = (None, "points", "points", None, None, None,)

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
#: For a derived argument, the one the caller passed it in
_BUILD_SPHERICAL_KD_ARGUMENT_SOURCES = (None, "points", "points", None, None, None,)

_lib.vicinity_vectors_c.restype = None
_lib.vicinity_vectors_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VICINITY_VECTORS_ARGUMENTS = ("query_point", "points", "n_dimensions", "n_points", "r", "dimension_order", "kd_indices", "vicinity_mask", "ierr",)
#: For a derived argument, the one the caller passed it in
_VICINITY_VECTORS_ARGUMENT_SOURCES = (None, None, "query_point", "points", None, None, None, None, None,)

_lib.vicinity_vectors_count_c.restype = None
_lib.vicinity_vectors_count_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VICINITY_VECTORS_COUNT_ARGUMENTS = ("query_point", "points", "n_dimensions", "n_points", "r", "dimension_order", "kd_indices", "n_neighbors", "ierr",)
#: For a derived argument, the one the caller passed it in
_VICINITY_VECTORS_COUNT_ARGUMENT_SOURCES = (None, None, "query_point", "points", None, None, None, None, None,)

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
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.

    Returns
    -------
    kd_indices : np.ndarray[np.int32] of shape (n_points,), read-only
        Output index array (k-d tree order)
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::build_kd_index`, whose argument names are
    the ones an error message reports.
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

    check_err_code(ierr.value, _BUILD_KD_INDEX_ARGUMENTS, _BUILD_KD_INDEX_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    kd_indices.flags.writeable = False

    return kd_indices

def build_spherical_kd(
        points,
        dimension_order,
):
    r"""Build a k-d tree index over points on the unit sphere (unit vectors)

    This is a thin, semantically-named wrapper: partitioning is identical to
    :func:`tensor_omics.build_kd_index` (plain per-axis median splits);
    callers are responsible for ensuring `points` are actually unit-normalized beforehand.

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Data points
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order (by variance)
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.

    Returns
    -------
    kd_indices : np.ndarray[np.int32] of shape (n_points,), read-only
        Output index array (k-d tree order)
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::build_spherical_kd`, whose argument names are
    the ones an error message reports.
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

    check_err_code(ierr.value, _BUILD_SPHERICAL_KD_ARGUMENTS, _BUILD_SPHERICAL_KD_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    kd_indices.flags.writeable = False

    return kd_indices

def vicinity_vectors(
        query_point,
        points,
        r,
        dimension_order,
        kd_indices,
):
    r"""Find reference points within a radius around a query point.

    Parameters
    ----------
    query_point : np.ndarray[np.float64] of shape (n_dimensions,)
        Coordinate vector used as the center of the search
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Ambient point matrix
    r : float
        Search radius
        The minimum valid value is `0.0`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Sequence of k-d tree split dimensions
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        K-d tree index sequence
        The minimum valid value is `1`.
        The maximum valid value is `n_points`.

    Returns
    -------
    vicinity_mask : np.ndarray[np.bool_] of shape (n_points,), read-only
        Mask indicating points within the search radius
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::vicinity_vectors`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        query_point = np.ascontiguousarray(query_point, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'query_point' must be an array of np.float64: {error}") from None
    if query_point.ndim != 1:
        raise ValueError(f"'query_point' must have 1 dimension, but has {query_point.ndim}")
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
    try:
        kd_indices = np.ascontiguousarray(kd_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'kd_indices' must be an array of np.int32: {error}") from None
    if kd_indices.ndim != 1:
        raise ValueError(f"'kd_indices' must have 1 dimension, but has {kd_indices.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = query_point.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if points.shape[0] != n_dimensions:
        raise ValueError(f"'points' has {points.shape[0]} along axis 0, but "
            f"'query_point' implies n_dimensions == {n_dimensions}"
        )
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'query_point' implies n_dimensions == {n_dimensions}"
        )
    if kd_indices.shape[0] != n_points:
        raise ValueError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but "
            f"'points' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    vicinity_mask = np.empty((n_points,), dtype=np.bool_, order='C')
    ierr = ctypes.c_int(0)

    _lib.vicinity_vectors_c(
        query_point,
        points,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_double(r)),
        dimension_order,
        kd_indices,
        vicinity_mask,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _VICINITY_VECTORS_ARGUMENTS, _VICINITY_VECTORS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    vicinity_mask.flags.writeable = False

    return vicinity_mask

def vicinity_vectors_count(
        query_point,
        points,
        r,
        dimension_order,
        kd_indices,
):
    r"""Count reference points within a radius around a query point.

    Parameters
    ----------
    query_point : np.ndarray[np.float64] of shape (n_dimensions,)
        Coordinate vector used as the center of the search
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Ambient point matrix
    r : float
        Search radius
        The minimum valid value is `0.0`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Sequence of k-d tree split dimensions
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        K-d tree index sequence
        The minimum valid value is `1`.
        The maximum valid value is `n_points`.

    Returns
    -------
    n_neighbors : int
        Number of points within the search radius

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::vicinity_vectors_count`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        query_point = np.ascontiguousarray(query_point, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'query_point' must be an array of np.float64: {error}") from None
    if query_point.ndim != 1:
        raise ValueError(f"'query_point' must have 1 dimension, but has {query_point.ndim}")
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
    try:
        kd_indices = np.ascontiguousarray(kd_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'kd_indices' must be an array of np.int32: {error}") from None
    if kd_indices.ndim != 1:
        raise ValueError(f"'kd_indices' must have 1 dimension, but has {kd_indices.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = query_point.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if points.shape[0] != n_dimensions:
        raise ValueError(f"'points' has {points.shape[0]} along axis 0, but "
            f"'query_point' implies n_dimensions == {n_dimensions}"
        )
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'query_point' implies n_dimensions == {n_dimensions}"
        )
    if kd_indices.shape[0] != n_points:
        raise ValueError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but "
            f"'points' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    n_neighbors = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.vicinity_vectors_count_c(
        query_point,
        points,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_double(r)),
        dimension_order,
        kd_indices,
        ctypes.byref(n_neighbors),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _VICINITY_VECTORS_COUNT_ARGUMENTS, _VICINITY_VECTORS_COUNT_ARGUMENT_SOURCES)

    return n_neighbors.value
