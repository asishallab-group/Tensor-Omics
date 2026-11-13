r"""tox_gene_centroids

Expression centroids of gene families.

`mean_vector` is the centroid of a set of expression vectors. `group_centroid_orthologs`
and `group_centroid_all` take the centroid of a family: over its orthologs only, or over
every gene in it -- two routines rather than one taking a flag, so which set a result is
over is visible at the call site.

Python binding, generated from tox_gene_centroids. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.mean_vector_c.restype = None
_lib.mean_vector_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MEAN_VECTOR_ARGUMENTS = ("expression_vectors", "n_axes", "n_genes", "genes_selection_mask", "n_selected_genes", "centroid", "ierr",)
#: For a derived argument, the one the caller passed it in
_MEAN_VECTOR_ARGUMENT_SOURCES = (None, "expression_vectors", "expression_vectors", None, "genes_selection_mask", None, None,)

_lib.group_centroid_orthologs_c.restype = None
_lib.group_centroid_orthologs_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GROUP_CENTROID_ORTHOLOGS_ARGUMENTS = ("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "n_omp_threads", "ortholog_set", "ierr",)
#: For a derived argument, the one the caller passed it in
_GROUP_CENTROID_ORTHOLOGS_ARGUMENT_SOURCES = (None, "expression_vectors", "expression_vectors", None, None, None, None, None, None,)

_lib.group_centroid_all_c.restype = None
_lib.group_centroid_all_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GROUP_CENTROID_ALL_ARGUMENTS = ("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "n_omp_threads", "ierr",)
#: For a derived argument, the one the caller passed it in
_GROUP_CENTROID_ALL_ARGUMENT_SOURCES = (None, "expression_vectors", "expression_vectors", None, None, None, None, None,)

def mean_vector(
        expression_vectors,
        genes_selection_mask,
):
    r"""Computes the element-wise mean of the selected gene vectors.

    A selection without genes gives the zero vector.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    genes_selection_mask : np.ndarray[np.bool_] of shape (n_genes,)
        `True` for the genes (columns of `expression_vectors`) to average

    Returns
    -------
    centroid : np.ndarray[np.float64] of shape (n_axes,), read-only
        The output vector representing the computed centroid.
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
    Generated from the Fortran procedure `tox_gene_centroids::mean_vector`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression_vectors = _contract.reals(expression_vectors, "expression_vectors", ndim=2, order="F")
    genes_selection_mask = _contract.logicals(genes_selection_mask, "genes_selection_mask", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]
    n_selected_genes = int(genes_selection_mask.sum())

    # Fortran cannot check that an array has the shape it is declared with; this can
    if genes_selection_mask.shape[0] != n_genes:
        raise _ToxArgumentError(f"'genes_selection_mask' has {genes_selection_mask.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'genes_selection_mask' is declared with shape (n_genes)", argument="genes_selection_mask")

    # outputs and work arrays, which the caller never sees
    centroid = _np.empty((n_axes,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.mean_vector_c(
        expression_vectors,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        _ctypes.byref(_ctypes.c_int(n_genes)),
        genes_selection_mask,
        _ctypes.byref(_ctypes.c_int(n_selected_genes)),
        centroid,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _MEAN_VECTOR_ARGUMENTS, _MEAN_VECTOR_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    centroid.flags.writeable = False

    return centroid

def group_centroid_orthologs(
        expression_vectors,
        gene_to_family,
        n_families,
        n_omp_threads,
        ortholog_set,
):
    r"""Computes one centroid per gene family, of all its genes or of its orthologs only.

    A family without selected genes gets the zero vector.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_to_family : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    n_families : int
        Total number of gene families to compute centroids for.
    n_omp_threads : int
        The most threads the parallel loop may use (OpenMP): it caps the team size.
        It also sizes the per-thread work arrays, where there are any, as their last extent: the plain entry point allocates that many columns, one in a build without OpenMP, and the expert one takes arrays with that many.
        It must be at least 1 in every build; without OpenMP the loop runs on one thread whatever is passed.
        Pass 1 when you parallelize around this call yourself, or call it in a forked worker process.
        The minimum valid value is `1`.
    ortholog_set : np.ndarray[np.bool_] of shape (n_genes,)
        A logical array indicating if a gene is part of a specific subset (e.g., orthologs).

    Returns
    -------
    centroid_matrix : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F'), read-only
        The output matrix (n_axes x n_families) to store the computed centroids.
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
    Generated from the Fortran procedure `tox_gene_centroids::group_centroid_orthologs`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression_vectors = _contract.reals(expression_vectors, "expression_vectors", ndim=2, order="F")
    gene_to_family = _contract.integers(gene_to_family, "gene_to_family", ndim=1)
    n_families = _contract.integer(n_families, "n_families")
    n_omp_threads = _contract.integer(n_omp_threads, "n_omp_threads")
    ortholog_set = _contract.logicals(ortholog_set, "ortholog_set", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if gene_to_family.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_family' has {gene_to_family.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'gene_to_family' is declared with shape (n_genes)", argument="gene_to_family")
    if ortholog_set.shape[0] != n_genes:
        raise _ToxArgumentError(f"'ortholog_set' has {ortholog_set.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'ortholog_set' is declared with shape (n_genes)", argument="ortholog_set")

    # outputs and work arrays, which the caller never sees
    centroid_matrix = _np.empty((n_axes, n_families,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.group_centroid_orthologs_c(
        expression_vectors,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        _ctypes.byref(_ctypes.c_int(n_genes)),
        gene_to_family,
        _ctypes.byref(_ctypes.c_int(n_families)),
        centroid_matrix,
        _ctypes.byref(_ctypes.c_int(n_omp_threads)),
        ortholog_set,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _GROUP_CENTROID_ORTHOLOGS_ARGUMENTS, _GROUP_CENTROID_ORTHOLOGS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    centroid_matrix.flags.writeable = False

    return centroid_matrix

def group_centroid_all(
        expression_vectors,
        gene_to_family,
        n_families,
        n_omp_threads,
):
    r"""Computes one centroid per gene family, of all its genes or of its orthologs only.

    A family without selected genes gets the zero vector.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_to_family : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    n_families : int
        Total number of gene families to compute centroids for.
    n_omp_threads : int
        The most threads the parallel loop may use (OpenMP): it caps the team size.
        It also sizes the per-thread work arrays, where there are any, as their last extent: the plain entry point allocates that many columns, one in a build without OpenMP, and the expert one takes arrays with that many.
        It must be at least 1 in every build; without OpenMP the loop runs on one thread whatever is passed.
        Pass 1 when you parallelize around this call yourself, or call it in a forked worker process.
        The minimum valid value is `1`.

    Returns
    -------
    centroid_matrix : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F'), read-only
        The output matrix (n_axes x n_families) to store the computed centroids.
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
    Generated from the Fortran procedure `tox_gene_centroids::group_centroid_all`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression_vectors = _contract.reals(expression_vectors, "expression_vectors", ndim=2, order="F")
    gene_to_family = _contract.integers(gene_to_family, "gene_to_family", ndim=1)
    n_families = _contract.integer(n_families, "n_families")
    n_omp_threads = _contract.integer(n_omp_threads, "n_omp_threads")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if gene_to_family.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_family' has {gene_to_family.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'gene_to_family' is declared with shape (n_genes)", argument="gene_to_family")

    # outputs and work arrays, which the caller never sees
    centroid_matrix = _np.empty((n_axes, n_families,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.group_centroid_all_c(
        expression_vectors,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        _ctypes.byref(_ctypes.c_int(n_genes)),
        gene_to_family,
        _ctypes.byref(_ctypes.c_int(n_families)),
        centroid_matrix,
        _ctypes.byref(_ctypes.c_int(n_omp_threads)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _GROUP_CENTROID_ALL_ARGUMENTS, _GROUP_CENTROID_ALL_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    centroid_matrix.flags.writeable = False

    return centroid_matrix
