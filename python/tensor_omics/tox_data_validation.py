r"""tox_data_validation

Semantic validation of TensorOmics data sets (dimensions, ID uniqueness, value ranges,
and cross-array consistency).

Complements ``tox_errors``'s generic per-argument validators
(`validate_dimension_size`, `validate_in_range_int/real`, ...) with checks specific to the
TensorOmics gene/family/expression/centroid/shift-vector data model, such as verifying that
`shift_vectors` was actually derived from `expression_vectors` and `family_centroids` as
expected. :func:`tensor_omics.validate_all_data` runs the full suite.

Python binding, generated from tox_data_validation. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.validate_data_structure_c.restype = None
_lib.validate_data_structure_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_DATA_STRUCTURE_ARGUMENTS = ("n_genes", "n_families", "n_samples", "gene_ids", "gene_family_ids", "gene_to_fam", "expression_vectors", "family_centroids", "shift_vectors", "ierr",)

_lib.validate_gene_to_family_mapping_c.restype = None
_lib.validate_gene_to_family_mapping_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_GENE_TO_FAMILY_MAPPING_ARGUMENTS = ("gene_to_fam", "n_families", "ierr",)

_lib.validate_expression_data_c.restype = None
_lib.validate_expression_data_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_EXPRESSION_DATA_ARGUMENTS = ("expression_vectors", "check_non_negative", "ierr",)

_lib.validate_family_centroids_c.restype = None
_lib.validate_family_centroids_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_FAMILY_CENTROIDS_ARGUMENTS = ("family_centroids", "ierr",)

_lib.validate_shift_vectors_c.restype = None
_lib.validate_shift_vectors_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_SHIFT_VECTORS_ARGUMENTS = ("shift_vectors", "expression_vectors", "family_centroids", "gene_to_fam", "n_samples", "ierr",)

_lib.validate_string_array_uniqueness_c.restype = None
_lib.validate_string_array_uniqueness_c.argtypes = (
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_STRING_ARRAY_UNIQUENESS_ARGUMENTS = ("str_arr", "ierr",)

_lib.validate_all_data_c.restype = None
_lib.validate_all_data_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_ALL_DATA_ARGUMENTS = ("n_genes", "n_families", "n_samples", "gene_ids", "gene_family_ids", "gene_to_fam", "expression_vectors", "family_centroids", "shift_vectors", "check_uniqueness", "check_shift_consistency", "ierr",)

def validate_data_structure(
        n_genes,
        n_families,
        n_samples,
        gene_ids,
        gene_family_ids,
        gene_to_fam,
        expression_vectors,
        family_centroids,
        shift_vectors,
):
    r"""Validate full data structure

    Parameters
    ----------
    n_genes : int
        Expected number of genes
    n_families : int
        Expected number of families
    n_samples : int
        Expected number of samples
    gene_ids : sequence of str, of length n_gene_ids_elements
        Gene ids
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    gene_family_ids : sequence of str, of length n_gene_family_ids_elements
        Gene family ids
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    gene_to_fam : np.ndarray[np.int32] of shape (n_gene_to_fam_elements,)
        gene to family mapping
    expression_vectors : np.ndarray[np.float64] of shape (n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2,), column-major (order='F')
        Expression vectors
    family_centroids : np.ndarray[np.float64] of shape (n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2,), column-major (order='F')
        Family centroids
    shift_vectors : np.ndarray[np.float64] of shape (n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2,), column-major (order='F')
        Shift vectors

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
    Generated from the Fortran procedure `tox_data_validation::validate_data_structure`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_genes = _contract.integer(n_genes, "n_genes")
    n_families = _contract.integer(n_families, "n_families")
    n_samples = _contract.integer(n_samples, "n_samples")
    _gene_ids_bytes = _contract.strings(gene_ids, "gene_ids")
    _gene_ids_width = max(map(len, _gene_ids_bytes), default=0) or 1
    gene_ids = _np.asarray([_b.ljust(_gene_ids_width) for _b in _gene_ids_bytes], dtype="S")
    _gene_family_ids_bytes = _contract.strings(gene_family_ids, "gene_family_ids")
    _gene_family_ids_width = max(map(len, _gene_family_ids_bytes), default=0) or 1
    gene_family_ids = _np.asarray([_b.ljust(_gene_family_ids_width) for _b in _gene_family_ids_bytes], dtype="S")
    gene_to_fam = _contract.integers(gene_to_fam, "gene_to_fam", ndim=1)
    expression_vectors = _contract.reals(expression_vectors, "expression_vectors", ndim=2, order="F")
    family_centroids = _contract.reals(family_centroids, "family_centroids", ndim=2, order="F")
    shift_vectors = _contract.reals(shift_vectors, "shift_vectors", ndim=2, order="F")

    # what the inputs already say, rather than asking for it again
    gene_ids_strlen = gene_ids.itemsize
    n_gene_ids_elements = gene_ids.shape[0]
    gene_family_ids_strlen = gene_family_ids.itemsize
    n_gene_family_ids_elements = gene_family_ids.shape[0]
    n_gene_to_fam_elements = gene_to_fam.shape[0]
    n_expression_vectors_elements_dim_1 = expression_vectors.shape[0]
    n_expression_vectors_elements_dim_2 = expression_vectors.shape[1]
    n_family_centroids_elements_dim_1 = family_centroids.shape[0]
    n_family_centroids_elements_dim_2 = family_centroids.shape[1]
    n_shift_vectors_elements_dim_1 = shift_vectors.shape[0]
    n_shift_vectors_elements_dim_2 = shift_vectors.shape[1]

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.validate_data_structure_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_families)),
        _ctypes.byref(_ctypes.c_int(n_samples)),
        gene_ids,
        _ctypes.byref(_ctypes.c_int(gene_ids_strlen)),
        _ctypes.byref(_ctypes.c_int(n_gene_ids_elements)),
        gene_family_ids,
        _ctypes.byref(_ctypes.c_int(gene_family_ids_strlen)),
        _ctypes.byref(_ctypes.c_int(n_gene_family_ids_elements)),
        gene_to_fam,
        _ctypes.byref(_ctypes.c_int(n_gene_to_fam_elements)),
        expression_vectors,
        _ctypes.byref(_ctypes.c_int(n_expression_vectors_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_expression_vectors_elements_dim_2)),
        family_centroids,
        _ctypes.byref(_ctypes.c_int(n_family_centroids_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_family_centroids_elements_dim_2)),
        shift_vectors,
        _ctypes.byref(_ctypes.c_int(n_shift_vectors_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_shift_vectors_elements_dim_2)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _VALIDATE_DATA_STRUCTURE_ARGUMENTS)

    return None

def validate_gene_to_family_mapping(
        gene_to_fam,
        n_families,
):
    r"""Validate gene to family mapping

    Parameters
    ----------
    gene_to_fam : np.ndarray[np.int32] of shape (n_gene_to_fam_elements,)
        gene to family mapping
    n_families : int
        number of families

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
    Generated from the Fortran procedure `tox_data_validation::validate_gene_to_family_mapping`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    gene_to_fam = _contract.integers(gene_to_fam, "gene_to_fam", ndim=1)
    n_families = _contract.integer(n_families, "n_families")

    # what the inputs already say, rather than asking for it again
    n_gene_to_fam_elements = gene_to_fam.shape[0]

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.validate_gene_to_family_mapping_c(
        gene_to_fam,
        _ctypes.byref(_ctypes.c_int(n_gene_to_fam_elements)),
        _ctypes.byref(_ctypes.c_int(n_families)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _VALIDATE_GENE_TO_FAMILY_MAPPING_ARGUMENTS)

    return None

def validate_expression_data(
        expression_vectors,
        check_non_negative,
):
    r"""Validate expresssion data

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2,), column-major (order='F')
        Expression vectors
    check_non_negative : bool
        Defines if non negative should be checked

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
    Generated from the Fortran procedure `tox_data_validation::validate_expression_data`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression_vectors = _contract.reals(expression_vectors, "expression_vectors", ndim=2, order="F")
    check_non_negative = _contract.logical(check_non_negative, "check_non_negative")

    # what the inputs already say, rather than asking for it again
    n_expression_vectors_elements_dim_1 = expression_vectors.shape[0]
    n_expression_vectors_elements_dim_2 = expression_vectors.shape[1]

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.validate_expression_data_c(
        expression_vectors,
        _ctypes.byref(_ctypes.c_int(n_expression_vectors_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_expression_vectors_elements_dim_2)),
        _ctypes.byref(_ctypes.c_bool(check_non_negative)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _VALIDATE_EXPRESSION_DATA_ARGUMENTS)

    return None

def validate_family_centroids(
        family_centroids,
):
    r"""Validate the family centroids

    Parameters
    ----------
    family_centroids : np.ndarray[np.float64] of shape (n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2,), column-major (order='F')
        Family centroids array

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
    Generated from the Fortran procedure `tox_data_validation::validate_family_centroids`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    family_centroids = _contract.reals(family_centroids, "family_centroids", ndim=2, order="F")

    # what the inputs already say, rather than asking for it again
    n_family_centroids_elements_dim_1 = family_centroids.shape[0]
    n_family_centroids_elements_dim_2 = family_centroids.shape[1]

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.validate_family_centroids_c(
        family_centroids,
        _ctypes.byref(_ctypes.c_int(n_family_centroids_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_family_centroids_elements_dim_2)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _VALIDATE_FAMILY_CENTROIDS_ARGUMENTS)

    return None

def validate_shift_vectors(
        shift_vectors,
        expression_vectors,
        family_centroids,
        gene_to_fam,
        n_samples,
):
    r"""Validates shift vectors

    Parameters
    ----------
    shift_vectors : np.ndarray[np.float64] of shape (n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2,), column-major (order='F')
        shift vectors
    expression_vectors : np.ndarray[np.float64] of shape (n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2,), column-major (order='F')
        expression vectors
    family_centroids : np.ndarray[np.float64] of shape (n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2,), column-major (order='F')
        family centroids
    gene_to_fam : np.ndarray[np.int32] of shape (n_gene_to_fam_elements,)
        gene to family mapping
    n_samples : int
        Number of samples

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
    Generated from the Fortran procedure `tox_data_validation::validate_shift_vectors`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    shift_vectors = _contract.reals(shift_vectors, "shift_vectors", ndim=2, order="F")
    expression_vectors = _contract.reals(expression_vectors, "expression_vectors", ndim=2, order="F")
    family_centroids = _contract.reals(family_centroids, "family_centroids", ndim=2, order="F")
    gene_to_fam = _contract.integers(gene_to_fam, "gene_to_fam", ndim=1)
    n_samples = _contract.integer(n_samples, "n_samples")

    # what the inputs already say, rather than asking for it again
    n_shift_vectors_elements_dim_1 = shift_vectors.shape[0]
    n_shift_vectors_elements_dim_2 = shift_vectors.shape[1]
    n_expression_vectors_elements_dim_1 = expression_vectors.shape[0]
    n_expression_vectors_elements_dim_2 = expression_vectors.shape[1]
    n_family_centroids_elements_dim_1 = family_centroids.shape[0]
    n_family_centroids_elements_dim_2 = family_centroids.shape[1]
    n_gene_to_fam_elements = gene_to_fam.shape[0]

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.validate_shift_vectors_c(
        shift_vectors,
        _ctypes.byref(_ctypes.c_int(n_shift_vectors_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_shift_vectors_elements_dim_2)),
        expression_vectors,
        _ctypes.byref(_ctypes.c_int(n_expression_vectors_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_expression_vectors_elements_dim_2)),
        family_centroids,
        _ctypes.byref(_ctypes.c_int(n_family_centroids_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_family_centroids_elements_dim_2)),
        gene_to_fam,
        _ctypes.byref(_ctypes.c_int(n_gene_to_fam_elements)),
        _ctypes.byref(_ctypes.c_int(n_samples)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _VALIDATE_SHIFT_VECTORS_ARGUMENTS)

    return None

def validate_string_array_uniqueness(
        str_arr,
):
    r"""Validate that no string appears more than once

    Parameters
    ----------
    str_arr : sequence of str, of length n_str_arr_elements
        string array
        Trailing blanks are not part of the value, as in Fortran: they are dropped.

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
    Generated from the Fortran procedure `tox_data_validation::validate_string_array_uniqueness`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    _str_arr_bytes = _contract.strings(str_arr, "str_arr")
    _str_arr_width = max(map(len, _str_arr_bytes), default=0) or 1
    str_arr = _np.asarray([_b.ljust(_str_arr_width) for _b in _str_arr_bytes], dtype="S")

    # what the inputs already say, rather than asking for it again
    str_arr_strlen = str_arr.itemsize
    n_str_arr_elements = str_arr.shape[0]

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.validate_string_array_uniqueness_c(
        str_arr,
        _ctypes.byref(_ctypes.c_int(str_arr_strlen)),
        _ctypes.byref(_ctypes.c_int(n_str_arr_elements)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _VALIDATE_STRING_ARRAY_UNIQUENESS_ARGUMENTS)

    return None

def validate_all_data(
        n_genes,
        n_families,
        n_samples,
        gene_ids,
        gene_family_ids,
        gene_to_fam,
        expression_vectors,
        family_centroids,
        shift_vectors,
        check_uniqueness=True,
        check_shift_consistency=True,
):
    r"""Comprehensive validation routine, combining all checks

    Parameters
    ----------
    n_genes : int
        Number of genes
    n_families : int
        Number of families
    n_samples : int
        Number of samples
    gene_ids : sequence of str, of length n_gene_ids_elements
        Gene ids array
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    gene_family_ids : sequence of str, of length n_gene_family_ids_elements
        gene family ids
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    gene_to_fam : np.ndarray[np.int32] of shape (n_gene_to_fam_elements,)
        gene to family mapping
    expression_vectors : np.ndarray[np.float64] of shape (n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2,), column-major (order='F')
        Expression vectors
    family_centroids : np.ndarray[np.float64] of shape (n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2,), column-major (order='F')
        family centroids
    shift_vectors : np.ndarray[np.float64] of shape (n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2,), column-major (order='F')
        shift vectors
    check_uniqueness : bool, optional, default True
        Check ID arrays for uniqueness.
        The default value is `True`.
    check_shift_consistency : bool, optional, default True
        Check consitency of shift array.
        The default value is `True`.

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
    Generated from the Fortran procedure `tox_data_validation::validate_all_data`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_genes = _contract.integer(n_genes, "n_genes")
    n_families = _contract.integer(n_families, "n_families")
    n_samples = _contract.integer(n_samples, "n_samples")
    _gene_ids_bytes = _contract.strings(gene_ids, "gene_ids")
    _gene_ids_width = max(map(len, _gene_ids_bytes), default=0) or 1
    gene_ids = _np.asarray([_b.ljust(_gene_ids_width) for _b in _gene_ids_bytes], dtype="S")
    _gene_family_ids_bytes = _contract.strings(gene_family_ids, "gene_family_ids")
    _gene_family_ids_width = max(map(len, _gene_family_ids_bytes), default=0) or 1
    gene_family_ids = _np.asarray([_b.ljust(_gene_family_ids_width) for _b in _gene_family_ids_bytes], dtype="S")
    gene_to_fam = _contract.integers(gene_to_fam, "gene_to_fam", ndim=1)
    expression_vectors = _contract.reals(expression_vectors, "expression_vectors", ndim=2, order="F")
    family_centroids = _contract.reals(family_centroids, "family_centroids", ndim=2, order="F")
    shift_vectors = _contract.reals(shift_vectors, "shift_vectors", ndim=2, order="F")
    check_uniqueness = _contract.logical(check_uniqueness, "check_uniqueness")
    check_shift_consistency = _contract.logical(check_shift_consistency, "check_shift_consistency")

    # what the inputs already say, rather than asking for it again
    gene_ids_strlen = gene_ids.itemsize
    n_gene_ids_elements = gene_ids.shape[0]
    gene_family_ids_strlen = gene_family_ids.itemsize
    n_gene_family_ids_elements = gene_family_ids.shape[0]
    n_gene_to_fam_elements = gene_to_fam.shape[0]
    n_expression_vectors_elements_dim_1 = expression_vectors.shape[0]
    n_expression_vectors_elements_dim_2 = expression_vectors.shape[1]
    n_family_centroids_elements_dim_1 = family_centroids.shape[0]
    n_family_centroids_elements_dim_2 = family_centroids.shape[1]
    n_shift_vectors_elements_dim_1 = shift_vectors.shape[0]
    n_shift_vectors_elements_dim_2 = shift_vectors.shape[1]

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.validate_all_data_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_families)),
        _ctypes.byref(_ctypes.c_int(n_samples)),
        gene_ids,
        _ctypes.byref(_ctypes.c_int(gene_ids_strlen)),
        _ctypes.byref(_ctypes.c_int(n_gene_ids_elements)),
        gene_family_ids,
        _ctypes.byref(_ctypes.c_int(gene_family_ids_strlen)),
        _ctypes.byref(_ctypes.c_int(n_gene_family_ids_elements)),
        gene_to_fam,
        _ctypes.byref(_ctypes.c_int(n_gene_to_fam_elements)),
        expression_vectors,
        _ctypes.byref(_ctypes.c_int(n_expression_vectors_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_expression_vectors_elements_dim_2)),
        family_centroids,
        _ctypes.byref(_ctypes.c_int(n_family_centroids_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_family_centroids_elements_dim_2)),
        shift_vectors,
        _ctypes.byref(_ctypes.c_int(n_shift_vectors_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_shift_vectors_elements_dim_2)),
        _ctypes.byref(_ctypes.c_bool(check_uniqueness)),
        _ctypes.byref(_ctypes.c_bool(check_shift_consistency)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _VALIDATE_ALL_DATA_ARGUMENTS)

    return None
