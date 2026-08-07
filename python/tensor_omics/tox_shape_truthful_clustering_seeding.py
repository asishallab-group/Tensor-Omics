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

_lib.density_labels_c.restype = None
_lib.density_labels_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DENSITY_LABELS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "k_density", "bandwidth_percentile", "labels", "ierr",)
#: For a derived argument, the one the caller passed it in
_DENSITY_LABELS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, None, None,)

_lib.seeds_c.restype = None
_lib.seeds_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SEEDS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "k_density", "bandwidth_percentile", "is_seed_mask", "ierr",)
#: For a derived argument, the one the caller passed it in
_SEEDS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, None, None,)

def density_labels(
        vectors,
        kd_indices,
        dimension_order,
        k_density=30,
        bandwidth_percentile=68.27,
):
    r"""Per-vector local density label, an adaptive-bandwidth kernel density estimate over each vector's own k_density nearest neighbors

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
    k_density : int, optional, default 30
        Neighborhood size the local density estimate is taken over
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors - 1`.
        The default value is `30`.
    bandwidth_percentile : float, optional, default 68.27
        Percentile (0 to 100) of the k_density neighbor distances used as the local
        Gaussian bandwidth -- a heuristic choice, not a calibrated standard deviation,
        see above
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `68.27`.

    Returns
    -------
    labels : np.ndarray[np.float64] of shape (n_vectors,), read-only
        Per-vector local density label
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
        ctypes.byref(ctypes.c_int(k_density)),
        ctypes.byref(ctypes.c_double(bandwidth_percentile)),
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
        k_density=30,
        bandwidth_percentile=68.27,
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
    k_density : int, optional, default 30
        Neighborhood size for both the density estimate and the coverage radius, see
        `density_labels` and `calc_ensemble_growth_radius`
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors - 1`.
        The default value is `30`.
    bandwidth_percentile : float, optional, default 68.27
        Percentile (0 to 100) of the k_density neighbor distances used as the local
        Gaussian bandwidth, see `density_labels`
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `68.27`.

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
        ctypes.byref(ctypes.c_int(k_density)),
        ctypes.byref(ctypes.c_double(bandwidth_percentile)),
        is_seed_mask,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SEEDS_ARGUMENTS, _SEEDS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    is_seed_mask.flags.writeable = False

    return is_seed_mask
