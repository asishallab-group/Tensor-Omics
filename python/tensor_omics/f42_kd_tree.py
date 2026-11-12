r"""f42_kd_tree

k-d tree spatial index over fixed-dimensional point sets.
Builds a k-d tree by recursively partitioning `kd_indices` around the median point along a
caller-supplied, cycling dimension order, using a stack-based (non-recursive) traversal so it
is safe to call from `pure` procedures. The tree is stored implicitly as an in-place-permuted
index array rather than as linked nodes.

Python binding, generated from f42_kd_tree. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.build_kd_index_c.restype = None
_lib.build_kd_index_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_KD_INDEX_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "ierr",)
#: For a derived argument, the one the caller passed it in
_BUILD_KD_INDEX_ARGUMENT_SOURCES = (None, "points", "points", None, None, None,)

_lib.build_spherical_kd_c.restype = None
_lib.build_spherical_kd_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_SPHERICAL_KD_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "ierr",)
#: For a derived argument, the one the caller passed it in
_BUILD_SPHERICAL_KD_ARGUMENT_SOURCES = (None, "points", "points", None, None, None,)

_lib.kd_knn_query_c.restype = None
_lib.kd_knn_query_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_KD_KNN_QUERY_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "query_point", "k_neighbors", "neighbors", "distances", "ierr",)
#: For a derived argument, the one the caller passed it in
_KD_KNN_QUERY_ARGUMENT_SOURCES = (None, "points", "points", None, None, None, None, None, None, None,)

_lib.kd_range_query_mask_c.restype = None
_lib.kd_range_query_mask_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_KD_RANGE_QUERY_MASK_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "query_point", "radius", "in_radius_mask", "ierr",)
#: For a derived argument, the one the caller passed it in
_KD_RANGE_QUERY_MASK_ARGUMENT_SOURCES = (None, "points", "points", None, None, None, None, None, None,)

_lib.kd_range_query_list_c.restype = None
_lib.kd_range_query_list_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_KD_RANGE_QUERY_LIST_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "query_point", "radius", "neighbors", "n_found", "ierr",)
#: For a derived argument, the one the caller passed it in
_KD_RANGE_QUERY_LIST_ARGUMENT_SOURCES = (None, "points", "points", None, None, None, None, None, None, None,)

_lib.kd_range_query_count_c.restype = None
_lib.kd_range_query_count_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::build_kd_index`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    points = _contract.reals(points, "points", ndim=2, order="F")
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")

    # outputs and work arrays, which the caller never sees
    kd_indices = _np.empty((n_points,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.build_kd_index_c(
        points,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BUILD_KD_INDEX_ARGUMENTS, _BUILD_KD_INDEX_ARGUMENT_SOURCES)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::build_spherical_kd`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    points = _contract.reals(points, "points", ndim=2, order="F")
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")

    # outputs and work arrays, which the caller never sees
    kd_indices = _np.empty((n_points,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.build_spherical_kd_c(
        points,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BUILD_SPHERICAL_KD_ARGUMENTS, _BUILD_SPHERICAL_KD_ARGUMENT_SOURCES)

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

    Via a bounded max-heap kept directly in `neighbors`/`distances` and splitting-plane
    pruning. Requires `k_neighbors <= n_points`: every point is then guaranteed visited
    before the heap can still have room, so no fallback for "fewer than k found" is needed.
    Does not guarantee `neighbors`/`distances` are sorted nearest-to-farthest (max-heap
    order internally).

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Original points dataset
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        Pre-built k-d tree index, see :func:`tensor_omics.build_kd_index`
        The minimum valid value is `1`.
        The maximum valid value is `n_points`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    query_point : np.ndarray[np.float64] of shape (n_dimensions,)
        Query point coordinates
    k_neighbors : int
        Number of neighbors to find
        The minimum valid value is `1`.
        The maximum valid value is `n_points`.

    Returns
    -------
    dict
        with keys:

        neighbors : np.ndarray[np.int32] of shape (k_neighbors,), read-only
            Output: indices of the k nearest neighbors (nearest-to-farthest order not guaranteed, max-heap order internally)
            A result is a value; call `.copy()` to obtain a modifiable array.
        distances : np.ndarray[np.float64] of shape (k_neighbors,), read-only
            Output: Euclidean distances to the k nearest neighbors
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::kd_knn_query`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    points = _contract.reals(points, "points", ndim=2, order="F")
    kd_indices = _contract.integers(kd_indices, "kd_indices", ndim=1)
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)
    query_point = _contract.reals(query_point, "query_point", ndim=1)
    k_neighbors = _contract.integer(k_neighbors, "k_neighbors")

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if kd_indices.shape[0] != n_points:
        raise _ToxArgumentError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but 'points' implies n_points == {n_points}; 'kd_indices' is declared with shape (n_points)", argument="kd_indices")
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")
    if query_point.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'query_point' has {query_point.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'query_point' is declared with shape (n_dimensions)", argument="query_point")

    # outputs and work arrays, which the caller never sees
    neighbors = _np.empty((k_neighbors,), dtype=_np.int32, order='C')
    distances = _np.empty((k_neighbors,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.kd_knn_query_c(
        points,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        query_point,
        _ctypes.byref(_ctypes.c_int(k_neighbors)),
        neighbors,
        distances,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _KD_KNN_QUERY_ARGUMENTS, _KD_KNN_QUERY_ARGUMENT_SOURCES)

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

    Same iterative, stack-based traversal and splitting-plane pruning as
    :func:`tensor_omics.kd_knn_query` (the near side of each split is
    always explored, the far side only when it is still within `radius` of the splitting
    plane), with a fixed radius bound instead of a k-nearest-neighbor heap. Compares
    squared distances against a precomputed `radius**2` (no `sqrt` per node visited).

    Fits a caller that already does an O(n_points) pass over the result (e.g. merging it
    into an existing coverage mask via `.or.`). A caller issuing many independent range
    queries per outer step (e.g. one per candidate point in a greedy loop) should use
    :func:`tensor_omics.kd_range_query_list` or
    :func:`tensor_omics.kd_range_query_count` instead, since the
    `in_radius_mask = False` reset here costs O(n_points) on every call regardless of how
    few points are actually found.

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Original points dataset
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        Pre-built k-d tree index
        The minimum valid value is `1`.
        The maximum valid value is `n_points`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    query_point : np.ndarray[np.float64] of shape (n_dimensions,)
        Query point coordinates
    radius : float
        Search radius
        The minimum valid value is `0.0`.

    Returns
    -------
    in_radius_mask : np.ndarray[np.bool_] of shape (n_points,), read-only
        Output: True for points within `radius` of `query_point`
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::kd_range_query_mask`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    points = _contract.reals(points, "points", ndim=2, order="F")
    kd_indices = _contract.integers(kd_indices, "kd_indices", ndim=1)
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)
    query_point = _contract.reals(query_point, "query_point", ndim=1)
    radius = _contract.real(radius, "radius")

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if kd_indices.shape[0] != n_points:
        raise _ToxArgumentError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but 'points' implies n_points == {n_points}; 'kd_indices' is declared with shape (n_points)", argument="kd_indices")
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")
    if query_point.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'query_point' has {query_point.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'query_point' is declared with shape (n_dimensions)", argument="query_point")

    # outputs and work arrays, which the caller never sees
    in_radius_mask = _np.empty((n_points,), dtype=_np.bool_, order='C')
    ierr = _ctypes.c_int(0)

    _lib.kd_range_query_mask_c(
        points,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        query_point,
        _ctypes.byref(_ctypes.c_double(radius)),
        in_radius_mask,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _KD_RANGE_QUERY_MASK_ARGUMENTS, _KD_RANGE_QUERY_MASK_ARGUMENT_SOURCES)

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

    Same traversal and pruning as :func:`tensor_omics.kd_range_query_mask`,
    but writes matches into a caller-provided compact index buffer (`neighbors(1:n_found)`)
    instead of a full-size logical mask, so repeated calls -- e.g. once per candidate point
    in an outer greedy loop -- don't each pay an O(n_points) reset. A caller that only
    needs the count, not the identities, should use
    :func:`tensor_omics.kd_range_query_count` instead, to skip the
    index-buffer writes entirely.

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Original points dataset
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        Pre-built k-d tree index
        The minimum valid value is `1`.
        The maximum valid value is `n_points`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    query_point : np.ndarray[np.float64] of shape (n_dimensions,)
        Query point coordinates
    radius : float
        Search radius
        The minimum valid value is `0.0`.

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::kd_range_query_list`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    points = _contract.reals(points, "points", ndim=2, order="F")
    kd_indices = _contract.integers(kd_indices, "kd_indices", ndim=1)
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)
    query_point = _contract.reals(query_point, "query_point", ndim=1)
    radius = _contract.real(radius, "radius")

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if kd_indices.shape[0] != n_points:
        raise _ToxArgumentError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but 'points' implies n_points == {n_points}; 'kd_indices' is declared with shape (n_points)", argument="kd_indices")
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")
    if query_point.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'query_point' has {query_point.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'query_point' is declared with shape (n_dimensions)", argument="query_point")

    # outputs and work arrays, which the caller never sees
    neighbors = _np.empty((n_points,), dtype=_np.int32, order='C')
    n_found = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.kd_range_query_list_c(
        points,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        query_point,
        _ctypes.byref(_ctypes.c_double(radius)),
        neighbors,
        _ctypes.byref(n_found),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _KD_RANGE_QUERY_LIST_ARGUMENTS, _KD_RANGE_QUERY_LIST_ARGUMENT_SOURCES)

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

    Same traversal and pruning as
    :func:`tensor_omics.kd_range_query_mask`, but writes no index
    buffer at all -- only a scalar count. The right choice when the identities of the
    points found are never needed, e.g. a per-point local-density label.

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Original points dataset
    kd_indices : np.ndarray[np.int32] of shape (n_points,)
        Pre-built k-d tree index
        The minimum valid value is `1`.
        The maximum valid value is `n_points`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    query_point : np.ndarray[np.float64] of shape (n_dimensions,)
        Query point coordinates
    radius : float
        Search radius
        The minimum valid value is `0.0`.

    Returns
    -------
    neighbor_count : int
        Output: number of points within `radius`

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_kd_tree::kd_range_query_count`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    points = _contract.reals(points, "points", ndim=2, order="F")
    kd_indices = _contract.integers(kd_indices, "kd_indices", ndim=1)
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)
    query_point = _contract.reals(query_point, "query_point", ndim=1)
    radius = _contract.real(radius, "radius")

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if kd_indices.shape[0] != n_points:
        raise _ToxArgumentError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but 'points' implies n_points == {n_points}; 'kd_indices' is declared with shape (n_points)", argument="kd_indices")
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")
    if query_point.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'query_point' has {query_point.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'query_point' is declared with shape (n_dimensions)", argument="query_point")

    # outputs and work arrays, which the caller never sees
    neighbor_count = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.kd_range_query_count_c(
        points,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        query_point,
        _ctypes.byref(_ctypes.c_double(radius)),
        _ctypes.byref(neighbor_count),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _KD_RANGE_QUERY_COUNT_ARGUMENTS, _KD_RANGE_QUERY_COUNT_ARGUMENT_SOURCES)

    return neighbor_count.value
