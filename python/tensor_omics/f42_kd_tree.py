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

_lib.kd_knn_query_c.restype = None
_lib.kd_knn_query_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_KD_KNN_QUERY_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "query_point", "k_neighbors", "neighbors", "distances", "ierr",)
#: For a derived argument, the one the caller passed it in
_KD_KNN_QUERY_ARGUMENT_SOURCES = (None, "points", "points", None, None, None, "neighbors", None, None, None,)

_lib.kd_range_query_mask_c.restype = None
_lib.kd_range_query_mask_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_KD_RANGE_QUERY_MASK_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "query_point", "radius", "in_radius_mask", "ierr",)
#: For a derived argument, the one the caller passed it in
_KD_RANGE_QUERY_MASK_ARGUMENT_SOURCES = (None, "points", "points", None, None, None, None, None, None,)

_lib.kd_range_query_list_c.restype = None
_lib.kd_range_query_list_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_KD_RANGE_QUERY_LIST_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "query_point", "radius", "neighbors", "n_found", "ierr",)
#: For a derived argument, the one the caller passed it in
_KD_RANGE_QUERY_LIST_ARGUMENT_SOURCES = (None, "points", "points", None, None, None, None, None, None, None,)

_lib.kd_range_query_count_c.restype = None
_lib.kd_range_query_count_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_KD_RANGE_QUERY_COUNT_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "query_point", "radius", "neighbor_count", "ierr",)
#: For a derived argument, the one the caller passed it in
_KD_RANGE_QUERY_COUNT_ARGUMENT_SOURCES = (None, "points", "points", None, None, None, None, None, None,)

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

def kd_knn_query(
        points,
        kd_indices,
        dimension_order,
        query_point,
        k_neighbors,
):
    r"""Find the k nearest neighbors of a query point in a pre-built k-d tree

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Original points dataset
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        Pre-built k-d tree index, see :func:`tensor_omics.build_kd_index`
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
    query_point : np.ndarray[np.float64] of shape (n_dimensions,)
        Query point coordinates
    k_neighbors : int
        Number of neighbors to find

    Returns
    -------
    dict
        with keys:

        neighbors : np.ndarray[np.int32] of shape (k_neighbors,), read-only
            guaranteed (max-heap order internally)
            A result is a value; call `.copy()` to obtain a modifiable array.
        distances : np.ndarray[np.float64] of shape (k_neighbors,), read-only
            Output: Euclidean distances to the k nearest neighbors
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::kd_knn_query_alloc`, whose argument names are
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
        kd_indices = np.ascontiguousarray(kd_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'kd_indices' must be an array of np.int32: {error}") from None
    if kd_indices.ndim != 1:
        raise ValueError(f"'kd_indices' must have 1 dimension, but has {kd_indices.ndim}")
    try:
        dimension_order = np.ascontiguousarray(dimension_order, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'dimension_order' must be an array of np.int32: {error}") from None
    if dimension_order.ndim != 1:
        raise ValueError(f"'dimension_order' must have 1 dimension, but has {dimension_order.ndim}")
    try:
        query_point = np.ascontiguousarray(query_point, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'query_point' must be an array of np.float64: {error}") from None
    if query_point.ndim != 1:
        raise ValueError(f"'query_point' must have 1 dimension, but has {query_point.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'points' implies n_dimensions == {n_dimensions}"
        )
    if query_point.shape[0] != n_dimensions:
        raise ValueError(f"'query_point' has {query_point.shape[0]} along axis 0, but "
            f"'points' implies n_dimensions == {n_dimensions}"
        )
    if kd_indices.shape[0] != n_points:
        raise ValueError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but "
            f"'points' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    neighbors = np.empty((k_neighbors,), dtype=np.int32, order='C')
    distances = np.empty((k_neighbors,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.kd_knn_query_c(
        points,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        query_point,
        ctypes.byref(ctypes.c_int(k_neighbors)),
        neighbors,
        distances,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _KD_KNN_QUERY_ARGUMENTS, _KD_KNN_QUERY_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    neighbors.flags.writeable = False
    distances.flags.writeable = False

    return {
        "neighbors": neighbors,
        "distances": distances,
    }

def kd_range_query_mask(
        points,
        kd_indices,
        dimension_order,
        query_point,
        radius,
):
    r"""Mark every point within `radius` of a query point in a pre-built k-d tree

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Original points dataset
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        Pre-built k-d tree index
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
    query_point : np.ndarray[np.float64] of shape (n_dimensions,)
        Query point coordinates
    radius : float
        Search radius

    Returns
    -------
    in_radius_mask : np.ndarray[np.bool_] of shape (n_points,), read-only
        Output: True for points within `radius` of `query_point`
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::kd_range_query_mask_alloc`, whose argument names are
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
        kd_indices = np.ascontiguousarray(kd_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'kd_indices' must be an array of np.int32: {error}") from None
    if kd_indices.ndim != 1:
        raise ValueError(f"'kd_indices' must have 1 dimension, but has {kd_indices.ndim}")
    try:
        dimension_order = np.ascontiguousarray(dimension_order, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'dimension_order' must be an array of np.int32: {error}") from None
    if dimension_order.ndim != 1:
        raise ValueError(f"'dimension_order' must have 1 dimension, but has {dimension_order.ndim}")
    try:
        query_point = np.ascontiguousarray(query_point, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'query_point' must be an array of np.float64: {error}") from None
    if query_point.ndim != 1:
        raise ValueError(f"'query_point' must have 1 dimension, but has {query_point.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'points' implies n_dimensions == {n_dimensions}"
        )
    if query_point.shape[0] != n_dimensions:
        raise ValueError(f"'query_point' has {query_point.shape[0]} along axis 0, but "
            f"'points' implies n_dimensions == {n_dimensions}"
        )
    if kd_indices.shape[0] != n_points:
        raise ValueError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but "
            f"'points' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    in_radius_mask = np.empty((n_points,), dtype=np.bool_, order='C')
    ierr = ctypes.c_int(0)

    _lib.kd_range_query_mask_c(
        points,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        query_point,
        ctypes.byref(ctypes.c_double(radius)),
        in_radius_mask,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _KD_RANGE_QUERY_MASK_ARGUMENTS, _KD_RANGE_QUERY_MASK_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    in_radius_mask.flags.writeable = False

    return in_radius_mask

def kd_range_query_list(
        points,
        kd_indices,
        dimension_order,
        query_point,
        radius,
):
    r"""List every point within `radius` of a query point in a pre-built k-d tree

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Original points dataset
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        Pre-built k-d tree index
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
    query_point : np.ndarray[np.float64] of shape (n_dimensions,)
        Query point coordinates
    radius : float
        Search radius

    Returns
    -------
    dict
        with keys:

        neighbors : np.ndarray[np.int32] of shape (n_points,), read-only
            Output: indices within `radius`, valid in `neighbors(1:n_found)`
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_found : int
            Output: number of points within `radius`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::kd_range_query_list_alloc`, whose argument names are
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
        kd_indices = np.ascontiguousarray(kd_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'kd_indices' must be an array of np.int32: {error}") from None
    if kd_indices.ndim != 1:
        raise ValueError(f"'kd_indices' must have 1 dimension, but has {kd_indices.ndim}")
    try:
        dimension_order = np.ascontiguousarray(dimension_order, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'dimension_order' must be an array of np.int32: {error}") from None
    if dimension_order.ndim != 1:
        raise ValueError(f"'dimension_order' must have 1 dimension, but has {dimension_order.ndim}")
    try:
        query_point = np.ascontiguousarray(query_point, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'query_point' must be an array of np.float64: {error}") from None
    if query_point.ndim != 1:
        raise ValueError(f"'query_point' must have 1 dimension, but has {query_point.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'points' implies n_dimensions == {n_dimensions}"
        )
    if query_point.shape[0] != n_dimensions:
        raise ValueError(f"'query_point' has {query_point.shape[0]} along axis 0, but "
            f"'points' implies n_dimensions == {n_dimensions}"
        )
    if kd_indices.shape[0] != n_points:
        raise ValueError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but "
            f"'points' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    neighbors = np.empty((n_points,), dtype=np.int32, order='C')
    n_found = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.kd_range_query_list_c(
        points,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        query_point,
        ctypes.byref(ctypes.c_double(radius)),
        neighbors,
        ctypes.byref(n_found),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _KD_RANGE_QUERY_LIST_ARGUMENTS, _KD_RANGE_QUERY_LIST_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    neighbors.flags.writeable = False

    return {
        "neighbors": neighbors,
        "n_found": n_found.value,
    }

def kd_range_query_count(
        points,
        kd_indices,
        dimension_order,
        query_point,
        radius,
):
    r"""Count the points within `radius` of a query point in a pre-built k-d tree

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Original points dataset
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        Pre-built k-d tree index
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
    query_point : np.ndarray[np.float64] of shape (n_dimensions,)
        Query point coordinates
    radius : float
        Search radius

    Returns
    -------
    neighbor_count : int
        Output: number of points within `radius`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::kd_range_query_count_alloc`, whose argument names are
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
        kd_indices = np.ascontiguousarray(kd_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'kd_indices' must be an array of np.int32: {error}") from None
    if kd_indices.ndim != 1:
        raise ValueError(f"'kd_indices' must have 1 dimension, but has {kd_indices.ndim}")
    try:
        dimension_order = np.ascontiguousarray(dimension_order, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'dimension_order' must be an array of np.int32: {error}") from None
    if dimension_order.ndim != 1:
        raise ValueError(f"'dimension_order' must have 1 dimension, but has {dimension_order.ndim}")
    try:
        query_point = np.ascontiguousarray(query_point, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'query_point' must be an array of np.float64: {error}") from None
    if query_point.ndim != 1:
        raise ValueError(f"'query_point' must have 1 dimension, but has {query_point.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'points' implies n_dimensions == {n_dimensions}"
        )
    if query_point.shape[0] != n_dimensions:
        raise ValueError(f"'query_point' has {query_point.shape[0]} along axis 0, but "
            f"'points' implies n_dimensions == {n_dimensions}"
        )
    if kd_indices.shape[0] != n_points:
        raise ValueError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but "
            f"'points' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    neighbor_count = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.kd_range_query_count_c(
        points,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        query_point,
        ctypes.byref(ctypes.c_double(radius)),
        ctypes.byref(neighbor_count),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _KD_RANGE_QUERY_COUNT_ARGUMENTS, _KD_RANGE_QUERY_COUNT_ARGUMENT_SOURCES)

    return neighbor_count.value
