"""tox_shape_truthful_clustering_parameter_estimation

Generated from the kernel; do not edit -- regenerate instead.

Python binding, generated from tox_shape_truthful_clustering_parameter_estimation. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.sample_estimator_anchors_c.restype = None
_lib.sample_estimator_anchors_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SAMPLE_ESTIMATOR_ANCHORS_ARGUMENTS = ("density_labels", "n_vectors", "n_anchors", "anchor_indices", "ierr",)
#: For a derived argument, the one the caller passed it in
_SAMPLE_ESTIMATOR_ANCHORS_ARGUMENT_SOURCES = (None, "density_labels", "anchor_indices", None, None,)

_lib.grow_estimator_anchor_clouds_c.restype = None
_lib.grow_estimator_anchor_clouds_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GROW_ESTIMATOR_ANCHOR_CLOUDS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "anchor_indices", "n_anchors", "seed_max_set_size", "cloud_masks", "cloud_sizes", "ierr",)
#: For a derived argument, the one the caller passed it in
_GROW_ESTIMATOR_ANCHOR_CLOUDS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, "anchor_indices", None, None, None, None,)

_lib.estimate_stc_parameters_c.restype = None
_lib.estimate_stc_parameters_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(ctypes.POINTER(ctypes.c_double)),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ESTIMATE_STC_PARAMETERS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "k_density", "bandwidth_percentile", "n_anchors", "seed_max_set_size", "first_quartile_percentile", "estimated_k_min", "estimated_k_density", "estimated_density_quantile", "estimated_chordal_dist_max_as_prcnt_of_range", "estimated_G_max", "estimated_d_max", "ierr",)
#: For a derived argument, the one the caller passed it in
_ESTIMATE_STC_PARAMETERS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def sample_estimator_anchors(
        density_labels,
        n_anchors,
):
    r"""Pick n_anchors point indices at evenly-spaced percentiles of the density-sorted order

    Parameters
    ----------
    density_labels : np.ndarray[np.float64] of shape (n_vectors,)
        Per-vector density label, see density_labels
    n_anchors : int
        Number of estimator anchors (EAs) to pick
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.

    Returns
    -------
    anchor_indices : np.ndarray[np.int32] of shape (n_anchors,), read-only
        Point indices of the n_anchors estimator anchors, ascending-percentile order
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_parameter_estimation::sample_estimator_anchors_alloc`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        density_labels = np.ascontiguousarray(density_labels, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'density_labels' must be an array of np.float64: {error}") from None
    if density_labels.ndim != 1:
        raise ValueError(f"'density_labels' must have 1 dimension, but has {density_labels.ndim}")

    # what the inputs already say, rather than asking for it again
    n_vectors = density_labels.shape[0]

    # outputs and work arrays, which the caller never sees
    anchor_indices = np.empty((n_anchors,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.sample_estimator_anchors_c(
        density_labels,
        ctypes.byref(ctypes.c_int(n_vectors)),
        ctypes.byref(ctypes.c_int(n_anchors)),
        anchor_indices,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SAMPLE_ESTIMATOR_ANCHORS_ARGUMENTS, _SAMPLE_ESTIMATOR_ANCHORS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    anchor_indices.flags.writeable = False

    return anchor_indices

def grow_estimator_anchor_clouds(
        vectors,
        anchor_indices,
        seed_max_set_size=5.0,
):
    r"""Multi-source competitive region growth of the n_anchors estimator anchors, bounded by seed_max_set_size

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix
    anchor_indices : np.ndarray[np.int32] of shape (n_anchors,)
        Point indices of the estimator anchors, see sample_estimator_anchors
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    seed_max_set_size : float, optional, default 5.0
        Percent (0 to 100) of n_vectors at which total growth across all EAs stops
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `5.0`.

    Returns
    -------
    dict
        with keys:

        cloud_masks : np.ndarray[np.bool_] of shape (n_vectors, n_anchors,), column-major (order='F'), read-only
            True for members of each EA's (column) final cloud, including its own anchor
            A result is a value; call `.copy()` to obtain a modifiable array.
        cloud_sizes : np.ndarray[np.int32] of shape (n_anchors,), read-only
            Final cloud size per EA
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_parameter_estimation::grow_estimator_anchor_clouds`, whose argument names are
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
        anchor_indices = np.ascontiguousarray(anchor_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'anchor_indices' must be an array of np.int32: {error}") from None
    if anchor_indices.ndim != 1:
        raise ValueError(f"'anchor_indices' must have 1 dimension, but has {anchor_indices.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]
    n_anchors = anchor_indices.shape[0]

    # outputs and work arrays, which the caller never sees
    cloud_masks = np.empty((n_vectors, n_anchors,), dtype=np.bool_, order='F')
    cloud_sizes = np.empty((n_anchors,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.grow_estimator_anchor_clouds_c(
        vectors,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        anchor_indices,
        ctypes.byref(ctypes.c_int(n_anchors)),
        ctypes.byref(ctypes.c_double(seed_max_set_size)),
        cloud_masks,
        cloud_sizes,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GROW_ESTIMATOR_ANCHOR_CLOUDS_ARGUMENTS, _GROW_ESTIMATOR_ANCHOR_CLOUDS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    cloud_masks.flags.writeable = False
    cloud_sizes.flags.writeable = False

    return {
        "cloud_masks": cloud_masks,
        "cloud_sizes": cloud_sizes,
    }

def estimate_stc_parameters(
        vectors,
        kd_indices,
        dimension_order,
        k_density=None,
        bandwidth_percentile=None,
        n_anchors=5,
        seed_max_set_size=5.0,
        first_quartile_percentile=25.0,
):
    r"""Estimate k_min, k_density, density_quantile, chordal_dist_max_as_prcnt_of_range, G_max, d_max from the data

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
    k_density : int, optional
        Passed through to density_labels
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors - 1`.
    bandwidth_percentile : float, optional
        Passed through to density_labels
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
    n_anchors : int, optional, default 5
        Number of estimator anchors (EAs), see sample_estimator_anchors
        The minimum valid value is `2`.
        The maximum valid value is `n_vectors`.
        The default value is `5`.
    seed_max_set_size : float, optional, default 5.0
        Passed through to grow_estimator_anchor_clouds
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `5.0`.
    first_quartile_percentile : float, optional, default 25.0
        Percentile (0 to 100) of the pairwise-EA-comparison distributions used for
        chordal_dist_max_as_prcnt_of_range/G_max/d_max, see estimate_stc_parameters
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `25.0`.

    Returns
    -------
    dict
        with keys:

        estimated_k_min : float
            Estimated k_min (real-valued; round for direct use as an integer argument)
        estimated_k_density : float
            Estimated k_density (equal to estimated_k_min, see estimate_stc_parameters)
        estimated_density_quantile : float
            Estimated density_quantile -- a literal radius (data units), not a percentile
        estimated_chordal_dist_max_as_prcnt_of_range : float
            Estimated chordal_dist_max_as_prcnt_of_range (0 to 1)
        estimated_G_max : float
            Estimated G_max
        estimated_d_max : float
            Estimated d_max (real-valued; round for direct use as an integer argument)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_parameter_estimation::estimate_stc_parameters_alloc`, whose argument names are
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
    estimated_k_min = ctypes.c_double(0)
    estimated_k_density = ctypes.c_double(0)
    estimated_density_quantile = ctypes.c_double(0)
    estimated_chordal_dist_max_as_prcnt_of_range = ctypes.c_double(0)
    estimated_G_max = ctypes.c_double(0)
    estimated_d_max = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.estimate_stc_parameters_c(
        vectors,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        kd_indices,
        dimension_order,
        None if k_density is None else ctypes.byref(ctypes.c_int(k_density)),
        None if bandwidth_percentile is None else ctypes.byref(ctypes.c_double(bandwidth_percentile)),
        ctypes.byref(ctypes.c_int(n_anchors)),
        ctypes.byref(ctypes.c_double(seed_max_set_size)),
        ctypes.byref(ctypes.c_double(first_quartile_percentile)),
        ctypes.byref(estimated_k_min),
        ctypes.byref(estimated_k_density),
        ctypes.byref(estimated_density_quantile),
        ctypes.byref(estimated_chordal_dist_max_as_prcnt_of_range),
        ctypes.byref(estimated_G_max),
        ctypes.byref(estimated_d_max),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _ESTIMATE_STC_PARAMETERS_ARGUMENTS, _ESTIMATE_STC_PARAMETERS_ARGUMENT_SOURCES)

    return {
        "estimated_k_min": estimated_k_min.value,
        "estimated_k_density": estimated_k_density.value,
        "estimated_density_quantile": estimated_density_quantile.value,
        "estimated_chordal_dist_max_as_prcnt_of_range": estimated_chordal_dist_max_as_prcnt_of_range.value,
        "estimated_G_max": estimated_G_max.value,
        "estimated_d_max": estimated_d_max.value,
    }
