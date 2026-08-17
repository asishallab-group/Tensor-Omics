r"""tox_shape_truthful_clustering_ensemble_growing

# Shape Truthful Clustering (STC): Ensemble Growing

Kernels for growing an ensemble by one step: `calc_ensemble_growth_radius` (the
per-seed, locally adaptive growth radius, computed once) and `grow_ensemble` (the union,
over every current member, of the points within that radius). See `misc/mod_STC.md`,
sections "Local Radius Identification" and SKG `grow_ensemble`, for the full algorithm
definition. Both take an already-built k-d tree (`kd_indices`, `dimension_order`, see
:func:`tensor_omics.build_kd_index`) as input, the same one shared
with the seeding kernels, rather than building their own.

Python binding, generated from tox_shape_truthful_clustering_ensemble_growing. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.calc_ensemble_growth_radius_c.restype = None
_lib.calc_ensemble_growth_radius_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_ENSEMBLE_GROWTH_RADIUS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "seed_index", "k_min", "radius_percentile", "growth_radius", "ierr",)
#: For a derived argument, the one the caller passed it in
_CALC_ENSEMBLE_GROWTH_RADIUS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, None, None, None,)

_lib.grow_ensemble_c.restype = None
_lib.grow_ensemble_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GROW_ENSEMBLE_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "is_member_mask", "growth_radius", "is_member_mask_next", "ierr",)
#: For a derived argument, the one the caller passed it in
_GROW_ENSEMBLE_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, None, None,)

def calc_ensemble_growth_radius(
        vectors,
        kd_indices,
        dimension_order,
        seed_index,
        k_min=30,
        radius_percentile=50.0,
):
    r"""Locally adapted ensemble growth radius, a percentile of the distances among a seed's own k_min nearest neighbors

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
    seed_index : int
        Index into `vectors`/`kd_indices` of the seed to compute the growth radius for
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    k_min : int, optional, default 30
        Neighborhood size the median distance is taken over
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors - 1`.
        The default value is `30`.
    radius_percentile : float, optional, default 50.0
        Percentile (0 to 100) of the k_min neighbor distances reported as the growth
        radius -- 50.0 (the default) is the median, matching this SKG's original,
        non-parameterized behavior
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `50.0`.

    Returns
    -------
    growth_radius : float
        radius_percentile-th percentile of the distances among the seed's own k_min
        nearest neighbors

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_ensemble_growing::calc_ensemble_growth_radius`, whose argument names are
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
    growth_radius = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.calc_ensemble_growth_radius_c(
        vectors,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        kd_indices,
        dimension_order,
        ctypes.byref(ctypes.c_int(seed_index)),
        ctypes.byref(ctypes.c_int(k_min)),
        ctypes.byref(ctypes.c_double(radius_percentile)),
        ctypes.byref(growth_radius),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_ENSEMBLE_GROWTH_RADIUS_ARGUMENTS, _CALC_ENSEMBLE_GROWTH_RADIUS_ARGUMENT_SOURCES)

    return growth_radius.value

def grow_ensemble(
        vectors,
        kd_indices,
        dimension_order,
        is_member_mask,
        growth_radius,
):
    r"""Grow an ensemble by one step, the union of every current member's growth-radius neighborhood

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
    is_member_mask : np.ndarray[np.bool_] of shape (n_vectors,)
        Current ensemble membership
    growth_radius : float
        This ensemble's growth radius, see `calc_ensemble_growth_radius`
        The minimum valid value is `0.0`.

    Returns
    -------
    is_member_mask_next : np.ndarray[np.bool_] of shape (n_vectors,), read-only
        Grown ensemble membership (superset of `is_member_mask`)
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_ensemble_growing::grow_ensemble`, whose argument names are
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
    try:
        is_member_mask = np.ascontiguousarray(is_member_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'is_member_mask' must be an array of np.bool_: {error}") from None
    if is_member_mask.ndim != 1:
        raise ValueError(f"'is_member_mask' must have 1 dimension, but has {is_member_mask.ndim}")

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
    if is_member_mask.shape[0] != n_vectors:
        raise ValueError(f"'is_member_mask' has {is_member_mask.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )

    # outputs and work arrays, which the caller never sees
    is_member_mask_next = np.empty((n_vectors,), dtype=np.bool_, order='C')
    ierr = ctypes.c_int(0)

    _lib.grow_ensemble_c(
        vectors,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        kd_indices,
        dimension_order,
        is_member_mask,
        ctypes.byref(ctypes.c_double(growth_radius)),
        is_member_mask_next,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GROW_ENSEMBLE_ARGUMENTS, _GROW_ENSEMBLE_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    is_member_mask_next.flags.writeable = False

    return is_member_mask_next
