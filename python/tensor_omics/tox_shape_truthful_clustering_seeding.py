"""tox_shape_truthful_clustering_seeding

Generated from the kernel; do not edit -- regenerate instead.

Python binding, generated from tox_shape_truthful_clustering_seeding. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.calculate_density_radius_c.restype = None
_lib.calculate_density_radius_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALCULATE_DENSITY_RADIUS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "mean_to_other_vecs_dist_quant", "radius", "ierr",)
#: For a derived argument, the one the caller passed it in
_CALCULATE_DENSITY_RADIUS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None,)

_lib.density_labels_c.restype = None
_lib.density_labels_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DENSITY_LABELS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "radius", "labels", "ierr",)
#: For a derived argument, the one the caller passed it in
_DENSITY_LABELS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, None,)

_lib.seeds_c.restype = None
_lib.seeds_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SEEDS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "mean_to_other_vecs_dist_quant", "is_seed_mask", "ierr",)
#: For a derived argument, the one the caller passed it in
_SEEDS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, None,)

def calculate_density_radius(
        vectors,
        mean_to_other_vecs_dist_quant=0.15,
):
    r"""Local-density search radius, a percentile of mean-to-vector distances

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix
    mean_to_other_vecs_dist_quant : float, optional, default 0.15
        Percentile (0.0 to 1.0) of mean-to-vector distances used as the density radius
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.15`.

    Returns
    -------
    radius : float
        Resulting density search radius

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_seeding::calculate_density_radius_alloc`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        vectors = np.asfortranarray(vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vectors' must be an array of np.float64: {error}") from None
    if vectors.ndim != 2:
        raise ValueError(f"'vectors' must have 2 dimensions, but has {vectors.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]

    # outputs and work arrays, which the caller never sees
    radius = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.calculate_density_radius_c(
        vectors,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        ctypes.byref(ctypes.c_double(mean_to_other_vecs_dist_quant)),
        ctypes.byref(radius),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALCULATE_DENSITY_RADIUS_ARGUMENTS, _CALCULATE_DENSITY_RADIUS_ARGUMENT_SOURCES)

    return radius.value

def density_labels(
        vectors,
        kd_indices,
        dimension_order,
        radius,
):
    r"""Per-vector density label, the count of vectors (including itself) within `radius`

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix
    kd_indices : np.ndarray[np.int32] of shape (n_vectors,)
        Pre-built k-d tree index over `vectors`
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    radius : float
        Density search radius, see `calculate_density_radius`
        The minimum valid value is `0.0`.

    Returns
    -------
    labels : np.ndarray[np.float64] of shape (n_vectors,), read-only
        Per-vector density label
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_seeding::density_labels_alloc`, whose argument names are
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

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if kd_indices.shape[0] != n_vectors:
        raise ValueError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )

    # outputs and work arrays, which the caller never sees
    labels = np.empty((n_vectors,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.density_labels_c(
        vectors,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        kd_indices,
        dimension_order,
        ctypes.byref(ctypes.c_double(radius)),
        labels,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DENSITY_LABELS_ARGUMENTS, _DENSITY_LABELS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    labels.flags.writeable = False

    return labels

def seeds(
        vectors,
        kd_indices,
        dimension_order,
        mean_to_other_vecs_dist_quant=0.15,
):
    r"""Select seed points via greedy, density-ranked, coverage-based selection

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix
    kd_indices : np.ndarray[np.int32] of shape (n_vectors,)
        Pre-built k-d tree index over `vectors`
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    mean_to_other_vecs_dist_quant : float, optional, default 0.15
        Density-radius percentile (0.0 to 1.0), see `calculate_density_radius`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.15`.

    Returns
    -------
    is_seed_mask : np.ndarray[np.bool_] of shape (n_vectors,), read-only
        True for points selected as seeds
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_seeding::seeds_alloc`, whose argument names are
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

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if kd_indices.shape[0] != n_vectors:
        raise ValueError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )

    # outputs and work arrays, which the caller never sees
    is_seed_mask = np.empty((n_vectors,), dtype=np.bool_, order='C')
    ierr = ctypes.c_int(0)

    _lib.seeds_c(
        vectors,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        kd_indices,
        dimension_order,
        ctypes.byref(ctypes.c_double(mean_to_other_vecs_dist_quant)),
        is_seed_mask,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SEEDS_ARGUMENTS, _SEEDS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    is_seed_mask.flags.writeable = False

    return is_seed_mask
