r"""tox_stc

Shape-truthful clustering with a global or per-seed growth radius.
Compute density labels, select seeds, grow their ensembles, and merge overlapping ensembles.
Both radius variants use the same density labels, seed ordering, growth rule and merge rule.

Python binding, generated from tox_stc. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.calculate_density_radius_c.restype = None
_lib.calculate_density_radius_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALCULATE_DENSITY_RADIUS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "mean_to_other_vecs_dist_quant", "radius", "ierr",)
#: For a derived argument, the one the caller passed it in
_CALCULATE_DENSITY_RADIUS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None,)

_lib.calculate_labels_as_density_c.restype = None
_lib.calculate_labels_as_density_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALCULATE_LABELS_AS_DENSITY_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "r", "dimension_order", "kd_indices", "n_tiles", "label_densities", "ierr",)
#: For a derived argument, the one the caller passed it in
_CALCULATE_LABELS_AS_DENSITY_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, None, None,)

_lib.identify_ensemble_seeds_c.restype = None
_lib.identify_ensemble_seeds_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_IDENTIFY_ENSEMBLE_SEEDS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "density_labels", "dimension_order", "kd_indices", "k_seeding", "sorted_perm", "n_seeds", "seed_mask", "seed_radii", "coverage_quant", "ierr",)
#: For a derived argument, the one the caller passed it in
_IDENTIFY_ENSEMBLE_SEEDS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, None, None, None, None, None, None,)

_lib.compute_ambient_density_stats_c.restype = None
_lib.compute_ambient_density_stats_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_AMBIENT_DENSITY_STATS_ARGUMENTS = ("density_labels", "n_vectors", "median_ambient", "mad_ambient", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_AMBIENT_DENSITY_STATS_ARGUMENT_SOURCES = (None, "density_labels", None, None, None,)

_lib.merge_ensembles_c.restype = None
_lib.merge_ensembles_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MERGE_ENSEMBLES_ARGUMENTS = ("raw_masks", "n_vectors", "n_seeds", "min_overlap_coefficient", "merged_masks", "n_ensembles", "ierr",)
#: For a derived argument, the one the caller passed it in
_MERGE_ENSEMBLES_ARGUMENT_SOURCES = (None, "raw_masks", "raw_masks", None, None, None, None,)

_lib.grow_ensemble_c.restype = None
_lib.grow_ensemble_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GROW_ENSEMBLE_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "ensemble_mask", "r", "dimension_order", "kd_indices", "density_labels", "alpha_mad", "grown_mask", "ierr",)
#: For a derived argument, the one the caller passed it in
_GROW_ENSEMBLE_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, None, None, None, None,)

_lib.compute_ensemble_observable_c.restype = None
_lib.compute_ensemble_observable_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ENSEMBLE_OBSERVABLE_ARGUMENTS = ("ensemble_mask", "density_labels", "n_vectors", "candidate_count", "current_iter", "observables", "n_history", "n_observables", "t_observables", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_ENSEMBLE_OBSERVABLE_ARGUMENT_SOURCES = (None, None, "ensemble_mask", None, None, None, "observables", "observables", None, None,)

_lib.accept_ensemble_c.restype = None
_lib.accept_ensemble_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ACCEPT_ENSEMBLE_ARGUMENTS = ("observables", "n_history", "n_observables", "current_iter", "alpha_accept", "is_accepted", "ierr",)
#: For a derived argument, the one the caller passed it in
_ACCEPT_ENSEMBLE_ARGUMENT_SOURCES = (None, "observables", "observables", None, None, None, None,)

_lib.obtain_ensembles_global_radius_c.restype = None
_lib.obtain_ensembles_global_radius_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_OBTAIN_ENSEMBLES_GLOBAL_RADIUS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "dimension_order", "kd_indices", "density_labels", "seed_indices", "n_seeds", "r", "alpha_mad", "alpha_accept", "t_observables", "n_tiles", "ensemble_matrix", "stop_reasons", "mad_ambient", "n_ensembles", "ierr",)
#: For a derived argument, the one the caller passed it in
_OBTAIN_ENSEMBLES_GLOBAL_RADIUS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, "seed_indices", None, None, None, None, None, None, None, None, None, None,)

_lib.obtain_ensembles_seed_radius_c.restype = None
_lib.obtain_ensembles_seed_radius_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_OBTAIN_ENSEMBLES_SEED_RADIUS_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "dimension_order", "kd_indices", "density_labels", "seed_indices", "n_seeds", "seed_radii", "alpha_mad", "alpha_accept", "t_observables", "n_tiles", "ensemble_matrix", "stop_reasons", "mad_ambient", "n_ensembles", "ierr",)
#: For a derived argument, the one the caller passed it in
_OBTAIN_ENSEMBLES_SEED_RADIUS_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, "seed_indices", None, None, None, None, None, None, None, None, None, None,)

def calculate_density_radius(
        vectors,
        mean_to_other_vecs_dist_quant=0.15,
):
    r"""Compute the density-label radius from distances to the data centroid.

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix (n_dimensions x n_vectors)
    mean_to_other_vecs_dist_quant : float, optional, default 0.15
        Optional quantile fraction (0.0 to 1.0), defaults to 0.15
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.15`.

    Returns
    -------
    radius : float
        The resulting density search radius

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_stc::calculate_density_radius`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    vectors = _contract.reals(vectors, "vectors", ndim=2, order="F")
    mean_to_other_vecs_dist_quant = _contract.real(mean_to_other_vecs_dist_quant, "mean_to_other_vecs_dist_quant")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]

    # outputs and work arrays, which the caller never sees
    radius = _ctypes.c_double(0)
    ierr = _ctypes.c_int(0)

    _lib.calculate_density_radius_c(
        vectors,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_vectors)),
        _ctypes.byref(_ctypes.c_double(mean_to_other_vecs_dist_quant)),
        _ctypes.byref(radius),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CALCULATE_DENSITY_RADIUS_ARGUMENTS, _CALCULATE_DENSITY_RADIUS_ARGUMENT_SOURCES)

    return radius.value

def calculate_labels_as_density(
        vectors,
        r,
        dimension_order,
        kd_indices,
        n_tiles,
):
    r"""Count points in each Euclidean neighborhood to label local density.

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix (n_dimensions x n_vectors)
    r : float
        Label-sphere search radius
        The minimum valid value is `0.0`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Sequence array tracking tree split axes by variance
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    kd_indices : np.ndarray[np.int32] of shape (n_vectors,)
        KD-tree index sequence array computed via :func:`tensor_omics.build_kd_index`.
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    n_tiles : int
        Number of concurrent vector tiles, at least 1. Capped at the number of points.
        The minimum valid value is `1`.

    Returns
    -------
    label_densities : np.ndarray[np.float64] of shape (n_vectors,), read-only
        Output array storing generated density scalars
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
    Generated from the Fortran procedure `tox_stc::calculate_labels_as_density`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    vectors = _contract.reals(vectors, "vectors", ndim=2, order="F")
    r = _contract.real(r, "r")
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)
    kd_indices = _contract.integers(kd_indices, "kd_indices", ndim=1)
    n_tiles = _contract.integer(n_tiles, "n_tiles")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'vectors' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")
    if kd_indices.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but 'vectors' implies n_vectors == {n_vectors}; 'kd_indices' is declared with shape (n_vectors)", argument="kd_indices")

    # outputs and work arrays, which the caller never sees
    label_densities = _np.empty((n_vectors,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.calculate_labels_as_density_c(
        vectors,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_vectors)),
        _ctypes.byref(_ctypes.c_double(r)),
        dimension_order,
        kd_indices,
        _ctypes.byref(_ctypes.c_int(n_tiles)),
        label_densities,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CALCULATE_LABELS_AS_DENSITY_ARGUMENTS, _CALCULATE_LABELS_AS_DENSITY_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    label_densities.flags.writeable = False

    return label_densities

def identify_ensemble_seeds(
        vectors,
        density_labels,
        dimension_order,
        kd_indices,
        k_seeding,
        coverage_quant=0.5,
):
    r"""Select seeds in descending density order and compute their local radii.

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix (n_dimensions x n_vectors)
    density_labels : np.ndarray[np.float64] of shape (n_vectors,)
        Precalculated density labels for all ambient vectors
        The minimum valid value is `1.0`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension split order array tracking the K-D tree structure
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    kd_indices : np.ndarray[np.int32] of shape (n_vectors,)
        K-D tree index sequence array
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    k_seeding : int
        Number of nearest neighbors used to estimate each seed's local coverage radius
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors - 1`.
    coverage_quant : float, optional, default 0.5
        Optional percentile of k-nearest-neighbor distances used as the seed coverage and growth radius; defaults to 0.5.
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.5`.

    Returns
    -------
    dict
        with keys:

        sorted_perm : np.ndarray[np.int32] of shape (n_vectors,), read-only
            Vector indices sorted by density descending
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_seeds : int
            Number of selected seed vectors
        seed_mask : np.ndarray[np.bool_] of shape (n_vectors,), read-only
            Logical mask identifying selected seed vectors
            A result is a value; call `.copy()` to obtain a modifiable array.
        seed_radii : np.ndarray[np.float64] of shape (n_vectors,), read-only
            Per-seed growth radius, stored at each seed's own vector index; zero at non-seeds
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
    Generated from the Fortran procedure `tox_stc::identify_ensemble_seeds`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    vectors = _contract.reals(vectors, "vectors", ndim=2, order="F")
    density_labels = _contract.reals(density_labels, "density_labels", ndim=1)
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)
    kd_indices = _contract.integers(kd_indices, "kd_indices", ndim=1)
    k_seeding = _contract.integer(k_seeding, "k_seeding")
    coverage_quant = _contract.real(coverage_quant, "coverage_quant")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if density_labels.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'density_labels' has {density_labels.shape[0]} along axis 0, but 'vectors' implies n_vectors == {n_vectors}; 'density_labels' is declared with shape (n_vectors)", argument="density_labels")
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'vectors' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")
    if kd_indices.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but 'vectors' implies n_vectors == {n_vectors}; 'kd_indices' is declared with shape (n_vectors)", argument="kd_indices")

    # outputs and work arrays, which the caller never sees
    sorted_perm = _np.empty((n_vectors,), dtype=_np.int32, order='C')
    n_seeds = _ctypes.c_int(0)
    seed_mask = _np.empty((n_vectors,), dtype=_np.bool_, order='C')
    seed_radii = _np.empty((n_vectors,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.identify_ensemble_seeds_c(
        vectors,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_vectors)),
        density_labels,
        dimension_order,
        kd_indices,
        _ctypes.byref(_ctypes.c_int(k_seeding)),
        sorted_perm,
        _ctypes.byref(n_seeds),
        seed_mask,
        seed_radii,
        _ctypes.byref(_ctypes.c_double(coverage_quant)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _IDENTIFY_ENSEMBLE_SEEDS_ARGUMENTS, _IDENTIFY_ENSEMBLE_SEEDS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    sorted_perm.flags.writeable = False
    seed_mask.flags.writeable = False
    seed_radii.flags.writeable = False

    return {
        "sorted_perm": sorted_perm,
        "n_seeds": n_seeds.value,
        "seed_mask": seed_mask,
        "seed_radii": seed_radii,
    }

def compute_ambient_density_stats(
        density_labels,
):
    r"""Compute the median and median absolute deviation of the density labels.

    Parameters
    ----------
    density_labels : np.ndarray[np.float64] of shape (n_vectors,)
        Density labels for all vectors the density-label routine.
        The minimum valid value is `1.0`.

    Returns
    -------
    dict
        with keys:

        median_ambient : float
            Ambient median density median_{i=1..N}(rho_i)
        mad_ambient : float
            Ambient median absolute deviation median_{i=1..N}(|rho_i - median_ambient|)

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_stc::compute_ambient_density_stats`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    density_labels = _contract.reals(density_labels, "density_labels", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_vectors = density_labels.shape[0]

    # outputs and work arrays, which the caller never sees
    median_ambient = _ctypes.c_double(0)
    mad_ambient = _ctypes.c_double(0)
    ierr = _ctypes.c_int(0)

    _lib.compute_ambient_density_stats_c(
        density_labels,
        _ctypes.byref(_ctypes.c_int(n_vectors)),
        _ctypes.byref(median_ambient),
        _ctypes.byref(mad_ambient),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_AMBIENT_DENSITY_STATS_ARGUMENTS, _COMPUTE_AMBIENT_DENSITY_STATS_ARGUMENT_SOURCES)

    return {
        "median_ambient": median_ambient.value,
        "mad_ambient": mad_ambient.value,
    }

def merge_ensembles(
        raw_masks,
        min_overlap_coefficient=0.9,
):
    r"""Merge intersecting ensembles by overlap coefficient and discard final singletons.

    Parameters
    ----------
    raw_masks : np.ndarray[np.bool_] of shape (n_vectors, n_seeds,), column-major (order='F')
        Matrix of raw unmerged boolean ensemble masks
    min_overlap_coefficient : float, optional, default 0.9
        Minimum Overlap Coefficient required to merge two ensembles
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.9`.

    Returns
    -------
    merged_masks : np.ndarray[np.bool_] of shape (n_vectors, n_seeds,), column-major (order='F'), read-only
        Output matrix storing merged boolean ensemble masks
        The first `n_ensembles` elements will hold the results.
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
    Generated from the Fortran procedure `tox_stc::merge_ensembles`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    raw_masks = _contract.logicals(raw_masks, "raw_masks", ndim=2, order="F")
    min_overlap_coefficient = _contract.real(min_overlap_coefficient, "min_overlap_coefficient")

    # what the inputs already say, rather than asking for it again
    n_vectors = raw_masks.shape[0]
    n_seeds = raw_masks.shape[1]

    # outputs and work arrays, which the caller never sees
    merged_masks = _np.empty((n_vectors, n_seeds,), dtype=_np.bool_, order='F')
    n_ensembles = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.merge_ensembles_c(
        raw_masks,
        _ctypes.byref(_ctypes.c_int(n_vectors)),
        _ctypes.byref(_ctypes.c_int(n_seeds)),
        _ctypes.byref(_ctypes.c_double(min_overlap_coefficient)),
        merged_masks,
        _ctypes.byref(n_ensembles),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _MERGE_ENSEMBLES_ARGUMENTS, _MERGE_ENSEMBLES_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    merged_masks.flags.writeable = False

    return merged_masks[..., :n_ensembles.value]

def grow_ensemble(
        vectors,
        ensemble_mask,
        r,
        dimension_order,
        kd_indices,
        density_labels,
        alpha_mad=0.5,
):
    r"""Grow one ensemble by one density-compatible surface step.

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix (n_dimensions x n_vectors)
    ensemble_mask : np.ndarray[np.bool_] of shape (n_vectors,)
        Membership before the growth step.
    r : float
        Search radius threshold for surface growth
        The minimum valid value is `0.0`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Sequence array tracking split dimensions
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    kd_indices : np.ndarray[np.int32] of shape (n_vectors,)
        KD-tree index sequence array computed via :func:`tensor_omics.build_kd_index`.
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    density_labels : np.ndarray[np.float64] of shape (n_vectors,)
        Density labels for all vectors the density-label routine.
        The minimum valid value is `1.0`.
    alpha_mad : float, optional, default 0.5
        Multiplier factor for MAD density compatibility threshold
        The minimum valid value is `0.0`.
        The default value is `0.5`.

    Returns
    -------
    grown_mask : np.ndarray[np.bool_] of shape (n_vectors,), read-only
        Membership after the growth step.
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
    Generated from the Fortran procedure `tox_stc::grow_ensemble`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    vectors = _contract.reals(vectors, "vectors", ndim=2, order="F")
    ensemble_mask = _contract.logicals(ensemble_mask, "ensemble_mask", ndim=1)
    r = _contract.real(r, "r")
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)
    kd_indices = _contract.integers(kd_indices, "kd_indices", ndim=1)
    density_labels = _contract.reals(density_labels, "density_labels", ndim=1)
    alpha_mad = _contract.real(alpha_mad, "alpha_mad")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if ensemble_mask.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'ensemble_mask' has {ensemble_mask.shape[0]} along axis 0, but 'vectors' implies n_vectors == {n_vectors}; 'ensemble_mask' is declared with shape (n_vectors)", argument="ensemble_mask")
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'vectors' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")
    if kd_indices.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but 'vectors' implies n_vectors == {n_vectors}; 'kd_indices' is declared with shape (n_vectors)", argument="kd_indices")
    if density_labels.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'density_labels' has {density_labels.shape[0]} along axis 0, but 'vectors' implies n_vectors == {n_vectors}; 'density_labels' is declared with shape (n_vectors)", argument="density_labels")

    # outputs and work arrays, which the caller never sees
    grown_mask = _np.empty((n_vectors,), dtype=_np.bool_, order='C')
    ierr = _ctypes.c_int(0)

    _lib.grow_ensemble_c(
        vectors,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_vectors)),
        ensemble_mask,
        _ctypes.byref(_ctypes.c_double(r)),
        dimension_order,
        kd_indices,
        density_labels,
        _ctypes.byref(_ctypes.c_double(alpha_mad)),
        grown_mask,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _GROW_ENSEMBLE_ARGUMENTS, _GROW_ENSEMBLE_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    grown_mask.flags.writeable = False

    return grown_mask

def compute_ensemble_observable(
        ensemble_mask,
        density_labels,
        candidate_count,
        current_iter,
        observables,
        t_observables=10,
):
    r"""Update the five ensemble observables for a growth iteration.

    Parameters
    ----------
    ensemble_mask : np.ndarray[np.bool_] of shape (n_vectors,)
        Selection mask for active ensemble members
    density_labels : np.ndarray[np.float64] of shape (n_vectors,)
        Density labels for all vectors the density-label routine.
        The minimum valid value is `1.0`.
    candidate_count : int
        Number of candidates added in current growth step
        The minimum valid value is `0`.
    current_iter : int
        Current growth iteration index (1-based)
        The minimum valid value is `1`.
    observables : np.ndarray[np.float64] of shape (n_observables, n_history,), column-major (order='F'), modified in place
        2D observable trajectory matrix where columns match iteration indices
    t_observables : int, optional, default 10
        History tracking depth (defaults to 10; <= 0 represents infinity / all)
        The default value is `10`.

    Returns
    -------
    None

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_stc::compute_ensemble_observable`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    ensemble_mask = _contract.logicals(ensemble_mask, "ensemble_mask", ndim=1)
    density_labels = _contract.reals(density_labels, "density_labels", ndim=1)
    candidate_count = _contract.integer(candidate_count, "candidate_count")
    current_iter = _contract.integer(current_iter, "current_iter")
    if not isinstance(observables, _np.ndarray) or observables.dtype != _np.float64:
        raise _ToxArgumentError("'observables' is modified in place, so it must already be a numpy array of {}".format(_np.float64), argument="observables")
    if observables.ndim != 2:
        raise _ToxArgumentError(f"'observables' must have 2 dimensions, but has {observables.ndim}", argument="observables")
    if not observables.flags.f_contiguous:
        raise _ToxArgumentError("'observables' is modified in place, so it must already be column-major (order='F')", argument="observables")
    t_observables = _contract.integer(t_observables, "t_observables")

    # what the inputs already say, rather than asking for it again
    n_vectors = ensemble_mask.shape[0]
    n_history = observables.shape[1]
    n_observables = observables.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if density_labels.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'density_labels' has {density_labels.shape[0]} along axis 0, but 'ensemble_mask' implies n_vectors == {n_vectors}; 'density_labels' is declared with shape (n_vectors)", argument="density_labels")

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.compute_ensemble_observable_c(
        ensemble_mask,
        density_labels,
        _ctypes.byref(_ctypes.c_int(n_vectors)),
        _ctypes.byref(_ctypes.c_int(candidate_count)),
        _ctypes.byref(_ctypes.c_int(current_iter)),
        observables,
        _ctypes.byref(_ctypes.c_int(n_history)),
        _ctypes.byref(_ctypes.c_int(n_observables)),
        _ctypes.byref(_ctypes.c_int(t_observables)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_ENSEMBLE_OBSERVABLE_ARGUMENTS, _COMPUTE_ENSEMBLE_OBSERVABLE_ARGUMENT_SOURCES)

    return None

def accept_ensemble(
        observables,
        current_iter,
        alpha_accept=0.5,
):
    r"""Test whether the change in ensemble observables permits a growth step.

    Parameters
    ----------
    observables : np.ndarray[np.float64] of shape (n_observables, n_history,), column-major (order='F')
        2D observables trajectory matrix [5 x window_capacity]
    current_iter : int
        Current growth iteration index (1-based)
        The minimum valid value is `1`.
    alpha_accept : float, optional, default 0.5
        Optional log2 fold change acceptance threshold (defaults to 0.5)
        The minimum valid value is `0.0`.
        The default value is `0.5`.

    Returns
    -------
    is_accepted : bool
        Output flag indicating step acceptance decision

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_stc::accept_ensemble`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    observables = _contract.reals(observables, "observables", ndim=2, order="F")
    current_iter = _contract.integer(current_iter, "current_iter")
    alpha_accept = _contract.real(alpha_accept, "alpha_accept")

    # what the inputs already say, rather than asking for it again
    n_history = observables.shape[1]
    n_observables = observables.shape[0]

    # outputs and work arrays, which the caller never sees
    is_accepted = _ctypes.c_bool(0)
    ierr = _ctypes.c_int(0)

    _lib.accept_ensemble_c(
        observables,
        _ctypes.byref(_ctypes.c_int(n_history)),
        _ctypes.byref(_ctypes.c_int(n_observables)),
        _ctypes.byref(_ctypes.c_int(current_iter)),
        _ctypes.byref(_ctypes.c_double(alpha_accept)),
        _ctypes.byref(is_accepted),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _ACCEPT_ENSEMBLE_ARGUMENTS, _ACCEPT_ENSEMBLE_ARGUMENT_SOURCES)

    return is_accepted.value

def obtain_ensembles_global_radius(
        vectors,
        dimension_order,
        kd_indices,
        density_labels,
        seed_indices,
        r,
        t_observables,
        n_tiles,
        alpha_mad=0.5,
        alpha_accept=0.5,
):
    r"""Grow one raw ensemble from each seed, using a global radius or per-seed radii.

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix (n_dimensions x n_vectors)
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Sequence array tracking split dimensions for KD-tree
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    kd_indices : np.ndarray[np.int32] of shape (n_vectors,)
        KD-tree index sequence array computed via :func:`tensor_omics.build_kd_index`.
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    density_labels : np.ndarray[np.float64] of shape (n_vectors,)
        Precalculated density labels for ambient vectors the density-label routine.
        The minimum valid value is `1.0`.
    seed_indices : np.ndarray[np.int32] of shape (n_seeds,)
        Array of starting seed vector indices the seed-selection routine.
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    r : float
        Uniform growth radius for every seed.
        The minimum valid value is `0.0`.
    alpha_mad : float, optional, default 0.5
        MAD multiplier for density compatibility threshold
        The minimum valid value is `0.0`.
        The default value is `0.5`.
    alpha_accept : float, optional, default 0.5
        Acceptance log2 fold-change threshold
        The minimum valid value is `0.0`.
        The default value is `0.5`.
    t_observables : int
        History window in iterations. Pass n_vectors to retain the full history.
        The minimum valid value is `1`.
    n_tiles : int
        Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
        The minimum valid value is `1`.

    Returns
    -------
    dict
        with keys:

        ensemble_matrix : np.ndarray[np.bool_] of shape (n_vectors, n_seeds,), column-major (order='F'), read-only
            Output matrix storing unmerged grown ensemble masks
            A result is a value; call `.copy()` to obtain a modifiable array.
        stop_reasons : np.ndarray[np.int32] of shape (n_seeds,), read-only
            Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
            3: first step rejected, 4: no candidates.
            A result is a value; call `.copy()` to obtain a modifiable array.
        mad_ambient : float
            Ambient MAD of density labels; zero makes alpha_mad ineffective by allowing only exactly equal densities.
        n_ensembles : int
            Total count of extracted raw ensembles (equals n_seeds)

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_stc::obtain_ensembles_global_radius`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    vectors = _contract.reals(vectors, "vectors", ndim=2, order="F")
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)
    kd_indices = _contract.integers(kd_indices, "kd_indices", ndim=1)
    density_labels = _contract.reals(density_labels, "density_labels", ndim=1)
    seed_indices = _contract.integers(seed_indices, "seed_indices", ndim=1)
    r = _contract.real(r, "r")
    alpha_mad = _contract.real(alpha_mad, "alpha_mad")
    alpha_accept = _contract.real(alpha_accept, "alpha_accept")
    t_observables = _contract.integer(t_observables, "t_observables")
    n_tiles = _contract.integer(n_tiles, "n_tiles")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]
    n_seeds = seed_indices.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'vectors' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")
    if kd_indices.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but 'vectors' implies n_vectors == {n_vectors}; 'kd_indices' is declared with shape (n_vectors)", argument="kd_indices")
    if density_labels.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'density_labels' has {density_labels.shape[0]} along axis 0, but 'vectors' implies n_vectors == {n_vectors}; 'density_labels' is declared with shape (n_vectors)", argument="density_labels")

    # outputs and work arrays, which the caller never sees
    ensemble_matrix = _np.empty((n_vectors, n_seeds,), dtype=_np.bool_, order='F')
    stop_reasons = _np.empty((n_seeds,), dtype=_np.int32, order='C')
    mad_ambient = _ctypes.c_double(0)
    n_ensembles = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.obtain_ensembles_global_radius_c(
        vectors,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_vectors)),
        dimension_order,
        kd_indices,
        density_labels,
        seed_indices,
        _ctypes.byref(_ctypes.c_int(n_seeds)),
        _ctypes.byref(_ctypes.c_double(r)),
        _ctypes.byref(_ctypes.c_double(alpha_mad)),
        _ctypes.byref(_ctypes.c_double(alpha_accept)),
        _ctypes.byref(_ctypes.c_int(t_observables)),
        _ctypes.byref(_ctypes.c_int(n_tiles)),
        ensemble_matrix,
        stop_reasons,
        _ctypes.byref(mad_ambient),
        _ctypes.byref(n_ensembles),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _OBTAIN_ENSEMBLES_GLOBAL_RADIUS_ARGUMENTS, _OBTAIN_ENSEMBLES_GLOBAL_RADIUS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    ensemble_matrix.flags.writeable = False
    stop_reasons.flags.writeable = False

    return {
        "ensemble_matrix": ensemble_matrix,
        "stop_reasons": stop_reasons,
        "mad_ambient": mad_ambient.value,
        "n_ensembles": n_ensembles.value,
    }

def obtain_ensembles_seed_radius(
        vectors,
        dimension_order,
        kd_indices,
        density_labels,
        seed_indices,
        seed_radii,
        t_observables,
        n_tiles,
        alpha_mad=0.5,
        alpha_accept=0.5,
):
    r"""Grow one raw ensemble from each seed, using a global radius or per-seed radii.

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix (n_dimensions x n_vectors)
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Sequence array tracking split dimensions for KD-tree
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    kd_indices : np.ndarray[np.int32] of shape (n_vectors,)
        KD-tree index sequence array computed via :func:`tensor_omics.build_kd_index`.
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    density_labels : np.ndarray[np.float64] of shape (n_vectors,)
        Precalculated density labels for ambient vectors the density-label routine.
        The minimum valid value is `1.0`.
    seed_indices : np.ndarray[np.int32] of shape (n_seeds,)
        Array of starting seed vector indices the seed-selection routine.
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    seed_radii : np.ndarray[np.float64] of shape (n_seeds,)
        Per-seed growth radius, allowing tighter growth in dense regions and wider growth in sparse regions.
        The minimum valid value is `0.0`.
    alpha_mad : float, optional, default 0.5
        MAD multiplier for density compatibility threshold
        The minimum valid value is `0.0`.
        The default value is `0.5`.
    alpha_accept : float, optional, default 0.5
        Acceptance log2 fold-change threshold
        The minimum valid value is `0.0`.
        The default value is `0.5`.
    t_observables : int
        History window in iterations. Pass n_vectors to retain the full history.
        The minimum valid value is `1`.
    n_tiles : int
        Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
        The minimum valid value is `1`.

    Returns
    -------
    dict
        with keys:

        ensemble_matrix : np.ndarray[np.bool_] of shape (n_vectors, n_seeds,), column-major (order='F'), read-only
            Output matrix storing unmerged grown ensemble masks
            A result is a value; call `.copy()` to obtain a modifiable array.
        stop_reasons : np.ndarray[np.int32] of shape (n_seeds,), read-only
            Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
            3: first step rejected, 4: no candidates.
            A result is a value; call `.copy()` to obtain a modifiable array.
        mad_ambient : float
            Ambient MAD of density labels; zero makes alpha_mad ineffective by allowing only exactly equal densities.
        n_ensembles : int
            Total count of extracted raw ensembles (equals n_seeds)

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_stc::obtain_ensembles_seed_radius`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    vectors = _contract.reals(vectors, "vectors", ndim=2, order="F")
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)
    kd_indices = _contract.integers(kd_indices, "kd_indices", ndim=1)
    density_labels = _contract.reals(density_labels, "density_labels", ndim=1)
    seed_indices = _contract.integers(seed_indices, "seed_indices", ndim=1)
    seed_radii = _contract.reals(seed_radii, "seed_radii", ndim=1)
    alpha_mad = _contract.real(alpha_mad, "alpha_mad")
    alpha_accept = _contract.real(alpha_accept, "alpha_accept")
    t_observables = _contract.integer(t_observables, "t_observables")
    n_tiles = _contract.integer(n_tiles, "n_tiles")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]
    n_seeds = seed_indices.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'vectors' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")
    if kd_indices.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but 'vectors' implies n_vectors == {n_vectors}; 'kd_indices' is declared with shape (n_vectors)", argument="kd_indices")
    if density_labels.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'density_labels' has {density_labels.shape[0]} along axis 0, but 'vectors' implies n_vectors == {n_vectors}; 'density_labels' is declared with shape (n_vectors)", argument="density_labels")
    if seed_radii.shape[0] != n_seeds:
        raise _ToxArgumentError(f"'seed_radii' has {seed_radii.shape[0]} along axis 0, but 'seed_indices' implies n_seeds == {n_seeds}; 'seed_radii' is declared with shape (n_seeds)", argument="seed_radii")

    # outputs and work arrays, which the caller never sees
    ensemble_matrix = _np.empty((n_vectors, n_seeds,), dtype=_np.bool_, order='F')
    stop_reasons = _np.empty((n_seeds,), dtype=_np.int32, order='C')
    mad_ambient = _ctypes.c_double(0)
    n_ensembles = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.obtain_ensembles_seed_radius_c(
        vectors,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_vectors)),
        dimension_order,
        kd_indices,
        density_labels,
        seed_indices,
        _ctypes.byref(_ctypes.c_int(n_seeds)),
        seed_radii,
        _ctypes.byref(_ctypes.c_double(alpha_mad)),
        _ctypes.byref(_ctypes.c_double(alpha_accept)),
        _ctypes.byref(_ctypes.c_int(t_observables)),
        _ctypes.byref(_ctypes.c_int(n_tiles)),
        ensemble_matrix,
        stop_reasons,
        _ctypes.byref(mad_ambient),
        _ctypes.byref(n_ensembles),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _OBTAIN_ENSEMBLES_SEED_RADIUS_ARGUMENTS, _OBTAIN_ENSEMBLES_SEED_RADIUS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    ensemble_matrix.flags.writeable = False
    stop_reasons.flags.writeable = False

    return {
        "ensemble_matrix": ensemble_matrix,
        "stop_reasons": stop_reasons,
        "mad_ambient": mad_ambient.value,
        "n_ensembles": n_ensembles.value,
    }
