r"""tox_euclidean_distance

Euclidean distances in gene expression space.

`euclidean_distance` between two expression vectors, and `distance_to_centroid` from every
gene of a family to that family's centroid -- the distance the outlier detection and the
relative-distance index are built on.

Python binding, generated from tox_euclidean_distance. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.euclidean_distance_c.restype = None
_lib.euclidean_distance_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_EUCLIDEAN_DISTANCE_ARGUMENTS = ("vec1", "vec2", "n_elements", "result", "ierr",)
#: For a derived argument, the one the caller passed it in
_EUCLIDEAN_DISTANCE_ARGUMENT_SOURCES = (None, None, "vec1", None, None,)

_lib.distance_to_centroid_c.restype = None
_lib.distance_to_centroid_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DISTANCE_TO_CENTROID_ARGUMENTS = ("n_genes", "n_families", "genes", "centroids", "gene_to_fam", "distances", "n_tissues", "ierr",)
#: For a derived argument, the one the caller passed it in
_DISTANCE_TO_CENTROID_ARGUMENT_SOURCES = ("genes", "centroids", None, None, None, None, "genes", None,)

def euclidean_distance(
        vec1,
        vec2,
):
    r"""Compute the Euclidean distance between two vectors.

    Calculates the L2 norm: `result = sqrt(sum((vec1_i - vec2_i)**2))`

    Parameters
    ----------
    vec1 : np.ndarray[np.float64] of shape (n_elements,)
        First expression vector
    vec2 : np.ndarray[np.float64] of shape (n_elements,)
        Second expression vector

    Returns
    -------
    result : float
        Output scalar distance

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_euclidean_distance::euclidean_distance`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        vec1 = _np.ascontiguousarray(vec1, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'vec1' must be an array of np.float64: {_error}", argument="vec1") from None
    if vec1.ndim != 1:
        raise _ToxArgumentError(f"'vec1' must have 1 dimension, but has {vec1.ndim}", argument="vec1")
    try:
        vec2 = _np.ascontiguousarray(vec2, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'vec2' must be an array of np.float64: {_error}", argument="vec2") from None
    if vec2.ndim != 1:
        raise _ToxArgumentError(f"'vec2' must have 1 dimension, but has {vec2.ndim}", argument="vec2")

    # what the inputs already say, rather than asking for it again
    n_elements = vec1.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if vec2.shape[0] != n_elements:
        raise _ToxArgumentError(f"'vec2' has {vec2.shape[0]} along axis 0, but 'vec1' implies n_elements == {n_elements}", argument="vec2")

    # outputs and work arrays, which the caller never sees
    result = _ctypes.c_double(0)
    ierr = _ctypes.c_int(0)

    _lib.euclidean_distance_c(
        vec1,
        vec2,
        _ctypes.byref(_ctypes.c_int(n_elements)),
        _ctypes.byref(result),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _EUCLIDEAN_DISTANCE_ARGUMENTS, _EUCLIDEAN_DISTANCE_ARGUMENT_SOURCES)

    return result.value

def distance_to_centroid(
        genes,
        centroids,
        gene_to_fam,
):
    r"""Compute distance from each gene to its corresponding family centroid.

    For each gene, extracts its expression vector and the centroid of its assigned family, then computes the Euclidean distance between them.

    Parameters
    ----------
    genes : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F')
        Gene expression matrix (n_tissues × n_genes), column-major
    centroids : np.ndarray[np.float64] of shape (n_tissues, n_families,), column-major (order='F')
        Family centroid matrix (n_tissues × n_families), column-major, `-1.0` for unassigned genes
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.

    Returns
    -------
    distances : np.ndarray[np.float64] of shape (n_genes,), read-only
        Output distances array
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
    Generated from the Fortran procedure `tox_euclidean_distance::distance_to_centroid`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        genes = _np.asfortranarray(genes, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'genes' must be an array of np.float64: {_error}", argument="genes") from None
    if genes.ndim != 2:
        raise _ToxArgumentError(f"'genes' must have 2 dimensions, but has {genes.ndim}", argument="genes")
    try:
        centroids = _np.asfortranarray(centroids, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'centroids' must be an array of np.float64: {_error}", argument="centroids") from None
    if centroids.ndim != 2:
        raise _ToxArgumentError(f"'centroids' must have 2 dimensions, but has {centroids.ndim}", argument="centroids")
    try:
        gene_to_fam = _np.ascontiguousarray(gene_to_fam, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'gene_to_fam' must be an array of np.int32: {_error}", argument="gene_to_fam") from None
    if gene_to_fam.ndim != 1:
        raise _ToxArgumentError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}", argument="gene_to_fam")

    # what the inputs already say, rather than asking for it again
    n_genes = genes.shape[1]
    n_families = centroids.shape[1]
    n_tissues = genes.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but 'genes' implies n_genes == {n_genes}", argument="gene_to_fam")
    if centroids.shape[0] != n_tissues:
        raise _ToxArgumentError(f"'centroids' has {centroids.shape[0]} along axis 0, but 'genes' implies n_tissues == {n_tissues}", argument="centroids")

    # outputs and work arrays, which the caller never sees
    distances = _np.empty((n_genes,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.distance_to_centroid_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_families)),
        genes,
        centroids,
        gene_to_fam,
        distances,
        _ctypes.byref(_ctypes.c_int(n_tissues)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _DISTANCE_TO_CENTROID_ARGUMENTS, _DISTANCE_TO_CENTROID_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    distances.flags.writeable = False

    return distances
