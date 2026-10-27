r"""tox_data_archive

Zip-archive backed persistence for TensorOmics data sets.

Wraps libzip (via C bindings declared in the interface block below) to create/extract zip
archives whose members are ``tox_data_read_write``-serialized arrays, indexed by a
plain-text `manifest.txt` mapping logical keys (e.g. `gene_ids`, `expression`) to member
filenames. :func:`tensor_omics.save_tox_data` and
``read_tox_data`` are the standard entry points for the
fixed TensorOmics data set schema; `create_zip_archive`/`extract_zip_archive` and the
`*_manifest*` routines below are the generic key/filename building blocks they are built on.

Python binding, generated from tox_data_archive. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.create_zip_archive_c.restype = None
_lib.create_zip_archive_c.argtypes = (
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CREATE_ZIP_ARCHIVE_ARGUMENTS = ("zip_filename", "keys", "filenames", "ierr",)

_lib.save_tox_data_c.restype = None
_lib.save_tox_data_c.argtypes = (
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(ndim=1)),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(ndim=1)),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS')),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(ndim=1)),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS')),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(ndim=1)),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(ndim=1)),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(ndim=1)),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS')),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(ndim=1)),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS')),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _nullable(_np.ctypeslib.ndpointer(ndim=1)),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SAVE_TOX_DATA_ARGUMENTS = ("zip_filename", "gene_ids", "gene_ids_file", "expression", "expression_file", "gene_to_family", "gene_to_family_file", "family_ids", "family_ids_file", "family_centroids", "family_centroids_file", "shift_vectors", "shift_vectors_file", "ierr",)

_lib.get_tox_data_dims_c.restype = None
_lib.get_tox_data_dims_c.argtypes = (
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GET_TOX_DATA_DIMS_ARGUMENTS = ("zip_filename", "n_gene_ids", "gene_id_len", "n_expression_rows", "n_expression_cols", "n_gene_to_family", "n_family_ids", "family_id_len", "n_family_centroids_rows", "n_family_centroids_cols", "n_shift_vectors_rows", "n_shift_vectors_cols", "ierr",)

_lib.read_tox_data_into_c.restype = None
_lib.read_tox_data_into_c.argtypes = (
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_READ_TOX_DATA_INTO_ARGUMENTS = ("zip_filename", "n_gene_ids", "gene_id_len", "gene_ids", "n_expression_rows", "n_expression_cols", "expression", "n_gene_to_family", "gene_to_family", "n_family_ids", "family_id_len", "family_ids", "n_family_centroids_rows", "n_family_centroids_cols", "family_centroids", "n_shift_vectors_rows", "n_shift_vectors_cols", "shift_vectors", "ierr",)
#: For a derived argument, the one the caller passed it in
_READ_TOX_DATA_INTO_ARGUMENT_SOURCES = (None, "gene_ids", None, None, "expression", "expression", None, "gene_to_family", None, "family_ids", None, None, "family_centroids", "family_centroids", None, "shift_vectors", "shift_vectors", None, None,)

def create_zip_archive(
        zip_filename,
        keys,
        filenames,
):
    r"""Creates a zip archive with generic file lists

    Parameters
    ----------
    zip_filename : str
        Name of the zip file to create
    keys : sequence of str, of length n_keys_elements
        Array of keys for manifest entries
    filenames : sequence of str, of length n_filenames_elements
        Array of filenames to add to zip

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
    Generated from the Fortran procedure `tox_data_archive::create_zip_archive`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    zip_filename = _np.array([str(zip_filename).encode().ljust(1)], dtype="S")
    try:
        _keys_bytes = [str(_s).encode() for _s in keys]
        _keys_width = max(map(len, _keys_bytes), default=0) or 1
        keys = _np.asarray([_b.ljust(_keys_width) for _b in _keys_bytes], dtype="S")
    except TypeError as _error:
        raise _ToxArgumentError(f"'keys' must be a sequence of strings: {_error}", argument="keys") from None
    if keys.ndim != 1:
        raise _ToxArgumentError(f"'keys' must have 1 dimension, but has {keys.ndim}", argument="keys")
    try:
        _filenames_bytes = [str(_s).encode() for _s in filenames]
        _filenames_width = max(map(len, _filenames_bytes), default=0) or 1
        filenames = _np.asarray([_b.ljust(_filenames_width) for _b in _filenames_bytes], dtype="S")
    except TypeError as _error:
        raise _ToxArgumentError(f"'filenames' must be a sequence of strings: {_error}", argument="filenames") from None
    if filenames.ndim != 1:
        raise _ToxArgumentError(f"'filenames' must have 1 dimension, but has {filenames.ndim}", argument="filenames")

    # what the inputs already say, rather than asking for it again
    zip_filename_strlen = zip_filename.itemsize
    keys_strlen = keys.itemsize
    n_keys_elements = keys.shape[0]
    filenames_strlen = filenames.itemsize
    n_filenames_elements = filenames.shape[0]

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.create_zip_archive_c(
        zip_filename,
        _ctypes.byref(_ctypes.c_int(zip_filename_strlen)),
        keys,
        _ctypes.byref(_ctypes.c_int(keys_strlen)),
        _ctypes.byref(_ctypes.c_int(n_keys_elements)),
        filenames,
        _ctypes.byref(_ctypes.c_int(filenames_strlen)),
        _ctypes.byref(_ctypes.c_int(n_filenames_elements)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CREATE_ZIP_ARCHIVE_ARGUMENTS)

    return None

def save_tox_data(
        zip_filename,
        gene_ids=None,
        gene_ids_file=None,
        expression=None,
        expression_file=None,
        gene_to_family=None,
        gene_to_family_file=None,
        family_ids=None,
        family_ids_file=None,
        family_centroids=None,
        family_centroids_file=None,
        shift_vectors=None,
        shift_vectors_file=None,
):
    r"""Save standard tox data

    Parameters
    ----------
    zip_filename : str
        Zip filename
    gene_ids : sequence of str, of length n_gene_ids_elements, optional
        Gene ids array, will be saved if provided
    gene_ids_file : str, optional
        Name of the gene ids file
    expression : np.ndarray[np.float64] of shape (n_expression_elements_dim_1, n_expression_elements_dim_2,), column-major (order='F'), optional
        Expression vectors array, will be saved if provided
    expression_file : str, optional
        Name of the expression file
    gene_to_family : np.ndarray[np.int32] of shape (n_gene_to_family_elements,), optional
        Gene to family mapping array, will be saved if provided
    gene_to_family_file : str, optional
        Name of the gene to family mapping file
    family_ids : sequence of str, of length n_family_ids_elements, optional
        Family ids array, will be saved if provided
    family_ids_file : str, optional
        Name of the family ids file
    family_centroids : np.ndarray[np.float64] of shape (n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2,), column-major (order='F'), optional
        Family centroids array, will be saved if provided
    family_centroids_file : str, optional
        Name of the family centroids file
    shift_vectors : np.ndarray[np.float64] of shape (n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2,), column-major (order='F'), optional
        Shift vectors array, will be saved if provided
    shift_vectors_file : str, optional
        Name of the shift vectors file

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
    Generated from the Fortran procedure `tox_data_archive::save_tox_data`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    zip_filename = _np.array([str(zip_filename).encode().ljust(1)], dtype="S")
    if gene_ids is not None:
        try:
            _gene_ids_bytes = [str(_s).encode() for _s in gene_ids]
            _gene_ids_width = max(map(len, _gene_ids_bytes), default=0) or 1
            gene_ids = _np.asarray([_b.ljust(_gene_ids_width) for _b in _gene_ids_bytes], dtype="S")
        except TypeError as _error:
            raise _ToxArgumentError(f"'gene_ids' must be a sequence of strings: {_error}", argument="gene_ids") from None
        if gene_ids.ndim != 1:
            raise _ToxArgumentError(f"'gene_ids' must have 1 dimension, but has {gene_ids.ndim}", argument="gene_ids")
    if gene_ids_file is not None:
        gene_ids_file = _np.array([str(gene_ids_file).encode().ljust(1)], dtype="S")
    if expression is not None:
        try:
            expression = _np.asfortranarray(expression, dtype=_np.float64)
        except (TypeError, ValueError) as _error:
            raise _ToxArgumentError(f"'expression' must be an array of np.float64: {_error}", argument="expression") from None
        if expression.ndim != 2:
            raise _ToxArgumentError(f"'expression' must have 2 dimensions, but has {expression.ndim}", argument="expression")
    if expression_file is not None:
        expression_file = _np.array([str(expression_file).encode().ljust(1)], dtype="S")
    if gene_to_family is not None:
        try:
            gene_to_family = _np.ascontiguousarray(gene_to_family, dtype=_np.int32)
        except (TypeError, ValueError) as _error:
            raise _ToxArgumentError(f"'gene_to_family' must be an array of np.int32: {_error}", argument="gene_to_family") from None
        if gene_to_family.ndim != 1:
            raise _ToxArgumentError(f"'gene_to_family' must have 1 dimension, but has {gene_to_family.ndim}", argument="gene_to_family")
    if gene_to_family_file is not None:
        gene_to_family_file = _np.array([str(gene_to_family_file).encode().ljust(1)], dtype="S")
    if family_ids is not None:
        try:
            _family_ids_bytes = [str(_s).encode() for _s in family_ids]
            _family_ids_width = max(map(len, _family_ids_bytes), default=0) or 1
            family_ids = _np.asarray([_b.ljust(_family_ids_width) for _b in _family_ids_bytes], dtype="S")
        except TypeError as _error:
            raise _ToxArgumentError(f"'family_ids' must be a sequence of strings: {_error}", argument="family_ids") from None
        if family_ids.ndim != 1:
            raise _ToxArgumentError(f"'family_ids' must have 1 dimension, but has {family_ids.ndim}", argument="family_ids")
    if family_ids_file is not None:
        family_ids_file = _np.array([str(family_ids_file).encode().ljust(1)], dtype="S")
    if family_centroids is not None:
        try:
            family_centroids = _np.asfortranarray(family_centroids, dtype=_np.float64)
        except (TypeError, ValueError) as _error:
            raise _ToxArgumentError(f"'family_centroids' must be an array of np.float64: {_error}", argument="family_centroids") from None
        if family_centroids.ndim != 2:
            raise _ToxArgumentError(f"'family_centroids' must have 2 dimensions, but has {family_centroids.ndim}", argument="family_centroids")
    if family_centroids_file is not None:
        family_centroids_file = _np.array([str(family_centroids_file).encode().ljust(1)], dtype="S")
    if shift_vectors is not None:
        try:
            shift_vectors = _np.asfortranarray(shift_vectors, dtype=_np.float64)
        except (TypeError, ValueError) as _error:
            raise _ToxArgumentError(f"'shift_vectors' must be an array of np.float64: {_error}", argument="shift_vectors") from None
        if shift_vectors.ndim != 2:
            raise _ToxArgumentError(f"'shift_vectors' must have 2 dimensions, but has {shift_vectors.ndim}", argument="shift_vectors")
    if shift_vectors_file is not None:
        shift_vectors_file = _np.array([str(shift_vectors_file).encode().ljust(1)], dtype="S")

    # what the inputs already say, rather than asking for it again
    zip_filename_strlen = zip_filename.itemsize
    gene_ids_strlen = 0 if gene_ids is None else gene_ids.itemsize
    n_gene_ids_elements = 0 if gene_ids is None else gene_ids.shape[0]
    gene_ids_file_strlen = 0 if gene_ids_file is None else gene_ids_file.itemsize
    n_expression_elements_dim_1 = 0 if expression is None else expression.shape[0]
    n_expression_elements_dim_2 = 0 if expression is None else expression.shape[1]
    expression_file_strlen = 0 if expression_file is None else expression_file.itemsize
    n_gene_to_family_elements = 0 if gene_to_family is None else gene_to_family.shape[0]
    gene_to_family_file_strlen = 0 if gene_to_family_file is None else gene_to_family_file.itemsize
    family_ids_strlen = 0 if family_ids is None else family_ids.itemsize
    n_family_ids_elements = 0 if family_ids is None else family_ids.shape[0]
    family_ids_file_strlen = 0 if family_ids_file is None else family_ids_file.itemsize
    n_family_centroids_elements_dim_1 = 0 if family_centroids is None else family_centroids.shape[0]
    n_family_centroids_elements_dim_2 = 0 if family_centroids is None else family_centroids.shape[1]
    family_centroids_file_strlen = 0 if family_centroids_file is None else family_centroids_file.itemsize
    n_shift_vectors_elements_dim_1 = 0 if shift_vectors is None else shift_vectors.shape[0]
    n_shift_vectors_elements_dim_2 = 0 if shift_vectors is None else shift_vectors.shape[1]
    shift_vectors_file_strlen = 0 if shift_vectors_file is None else shift_vectors_file.itemsize

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.save_tox_data_c(
        zip_filename,
        _ctypes.byref(_ctypes.c_int(zip_filename_strlen)),
        gene_ids,
        _ctypes.byref(_ctypes.c_int(gene_ids_strlen)),
        _ctypes.byref(_ctypes.c_int(n_gene_ids_elements)),
        gene_ids_file,
        _ctypes.byref(_ctypes.c_int(gene_ids_file_strlen)),
        expression,
        _ctypes.byref(_ctypes.c_int(n_expression_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_expression_elements_dim_2)),
        expression_file,
        _ctypes.byref(_ctypes.c_int(expression_file_strlen)),
        gene_to_family,
        _ctypes.byref(_ctypes.c_int(n_gene_to_family_elements)),
        gene_to_family_file,
        _ctypes.byref(_ctypes.c_int(gene_to_family_file_strlen)),
        family_ids,
        _ctypes.byref(_ctypes.c_int(family_ids_strlen)),
        _ctypes.byref(_ctypes.c_int(n_family_ids_elements)),
        family_ids_file,
        _ctypes.byref(_ctypes.c_int(family_ids_file_strlen)),
        family_centroids,
        _ctypes.byref(_ctypes.c_int(n_family_centroids_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_family_centroids_elements_dim_2)),
        family_centroids_file,
        _ctypes.byref(_ctypes.c_int(family_centroids_file_strlen)),
        shift_vectors,
        _ctypes.byref(_ctypes.c_int(n_shift_vectors_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_shift_vectors_elements_dim_2)),
        shift_vectors_file,
        _ctypes.byref(_ctypes.c_int(shift_vectors_file_strlen)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _SAVE_TOX_DATA_ARGUMENTS)

    return None

def get_tox_data_dims(
        zip_filename,
):
    r"""Report the shape of every member of a tox data archive

    Each count (and each string length) is 0 when the corresponding member is absent, so a
    caller can size all six output buffers up front. Character members report both an element
    count and a per-element string length. Pairs with
    :func:`tensor_omics.read_tox_data_into`.

    Parameters
    ----------
    zip_filename : str
        Name of the zip file

    Returns
    -------
    dict
        with keys:

        n_gene_ids : int
            Number of gene ids, 0 if absent
        gene_id_len : int
            String length of each gene id, 0 if absent
        n_expression_rows : int
            Rows (samples) of the expression matrix, 0 if absent
        n_expression_cols : int
            Columns (genes) of the expression matrix, 0 if absent
        n_gene_to_family : int
            Number of gene-to-family entries, 0 if absent
        n_family_ids : int
            Number of family ids, 0 if absent
        family_id_len : int
            String length of each family id, 0 if absent
        n_family_centroids_rows : int
            Rows (samples) of the family centroids matrix, 0 if absent
        n_family_centroids_cols : int
            Columns (families) of the family centroids matrix, 0 if absent
        n_shift_vectors_rows : int
            Rows of the shift vectors matrix, 0 if absent
        n_shift_vectors_cols : int
            Columns of the shift vectors matrix, 0 if absent

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_archive::get_tox_data_dims`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    zip_filename = _np.array([str(zip_filename).encode().ljust(1)], dtype="S")

    # what the inputs already say, rather than asking for it again
    zip_filename_strlen = zip_filename.itemsize

    # outputs and work arrays, which the caller never sees
    n_gene_ids = _ctypes.c_int(0)
    gene_id_len = _ctypes.c_int(0)
    n_expression_rows = _ctypes.c_int(0)
    n_expression_cols = _ctypes.c_int(0)
    n_gene_to_family = _ctypes.c_int(0)
    n_family_ids = _ctypes.c_int(0)
    family_id_len = _ctypes.c_int(0)
    n_family_centroids_rows = _ctypes.c_int(0)
    n_family_centroids_cols = _ctypes.c_int(0)
    n_shift_vectors_rows = _ctypes.c_int(0)
    n_shift_vectors_cols = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.get_tox_data_dims_c(
        zip_filename,
        _ctypes.byref(_ctypes.c_int(zip_filename_strlen)),
        _ctypes.byref(n_gene_ids),
        _ctypes.byref(gene_id_len),
        _ctypes.byref(n_expression_rows),
        _ctypes.byref(n_expression_cols),
        _ctypes.byref(n_gene_to_family),
        _ctypes.byref(n_family_ids),
        _ctypes.byref(family_id_len),
        _ctypes.byref(n_family_centroids_rows),
        _ctypes.byref(n_family_centroids_cols),
        _ctypes.byref(n_shift_vectors_rows),
        _ctypes.byref(n_shift_vectors_cols),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _GET_TOX_DATA_DIMS_ARGUMENTS)

    return {
        "n_gene_ids": n_gene_ids.value,
        "gene_id_len": gene_id_len.value,
        "n_expression_rows": n_expression_rows.value,
        "n_expression_cols": n_expression_cols.value,
        "n_gene_to_family": n_gene_to_family.value,
        "n_family_ids": n_family_ids.value,
        "family_id_len": family_id_len.value,
        "n_family_centroids_rows": n_family_centroids_rows.value,
        "n_family_centroids_cols": n_family_centroids_cols.value,
        "n_shift_vectors_rows": n_shift_vectors_rows.value,
        "n_shift_vectors_cols": n_shift_vectors_cols.value,
    }

def read_tox_data_into(
        zip_filename,
):
    r"""Read a tox data archive into caller-provided buffers

    Fills every buffer from the archive; size them from
    :func:`tensor_omics.get_tox_data_dims` first. A member that is absent
    has a zero extent and is left untouched.

    Parameters
    ----------
    zip_filename : str
        Name of the zip file

    Returns
    -------
    dict
        with keys:

        gene_ids : sequence of str, of length n_gene_ids
            Gene ids
        expression : np.ndarray[np.float64] of shape (n_expression_rows, n_expression_cols,), column-major (order='F'), read-only
            Expression vectors
            A result is a value; call `.copy()` to obtain a modifiable array.
        gene_to_family : np.ndarray[np.int32] of shape (n_gene_to_family,), read-only
            Gene to family mapping
            A result is a value; call `.copy()` to obtain a modifiable array.
        family_ids : sequence of str, of length n_family_ids
            Family ids
        family_centroids : np.ndarray[np.float64] of shape (n_family_centroids_rows, n_family_centroids_cols,), column-major (order='F'), read-only
            Family centroids
            A result is a value; call `.copy()` to obtain a modifiable array.
        shift_vectors : np.ndarray[np.float64] of shape (n_shift_vectors_rows, n_shift_vectors_cols,), column-major (order='F'), read-only
            Shift vectors
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
    Generated from the Fortran procedure `tox_data_archive::read_tox_data_into`, whose argument names are
    the ones an error message reports.
    """
    # kept before conversion, for the producers called below
    _zip_filename_raw = zip_filename

    # accept anything array-like, converting only when C needs it
    zip_filename = _np.array([str(zip_filename).encode().ljust(1)], dtype="S")

    # what the inputs already say, rather than asking for it again
    zip_filename_strlen = zip_filename.itemsize

    # work out what other procedures must supply, per DM_OUTPUT_FROM
    _get_tox_data_dims_result = get_tox_data_dims(zip_filename=_zip_filename_raw)
    n_gene_ids = _get_tox_data_dims_result["n_gene_ids"]
    gene_id_len = _get_tox_data_dims_result["gene_id_len"]
    n_expression_rows = _get_tox_data_dims_result["n_expression_rows"]
    n_expression_cols = _get_tox_data_dims_result["n_expression_cols"]
    n_gene_to_family = _get_tox_data_dims_result["n_gene_to_family"]
    n_family_ids = _get_tox_data_dims_result["n_family_ids"]
    family_id_len = _get_tox_data_dims_result["family_id_len"]
    n_family_centroids_rows = _get_tox_data_dims_result["n_family_centroids_rows"]
    n_family_centroids_cols = _get_tox_data_dims_result["n_family_centroids_cols"]
    n_shift_vectors_rows = _get_tox_data_dims_result["n_shift_vectors_rows"]
    n_shift_vectors_cols = _get_tox_data_dims_result["n_shift_vectors_cols"]

    # outputs and work arrays, which the caller never sees
    gene_ids = _np.zeros((n_gene_ids,), dtype=f"S{gene_id_len}")
    expression = _np.empty((n_expression_rows, n_expression_cols,), dtype=_np.float64, order='F')
    gene_to_family = _np.empty((n_gene_to_family,), dtype=_np.int32, order='C')
    family_ids = _np.zeros((n_family_ids,), dtype=f"S{family_id_len}")
    family_centroids = _np.empty((n_family_centroids_rows, n_family_centroids_cols,), dtype=_np.float64, order='F')
    shift_vectors = _np.empty((n_shift_vectors_rows, n_shift_vectors_cols,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.read_tox_data_into_c(
        zip_filename,
        _ctypes.byref(_ctypes.c_int(zip_filename_strlen)),
        _ctypes.byref(_ctypes.c_int(n_gene_ids)),
        _ctypes.byref(_ctypes.c_int(gene_id_len)),
        gene_ids,
        _ctypes.byref(_ctypes.c_int(n_expression_rows)),
        _ctypes.byref(_ctypes.c_int(n_expression_cols)),
        expression,
        _ctypes.byref(_ctypes.c_int(n_gene_to_family)),
        gene_to_family,
        _ctypes.byref(_ctypes.c_int(n_family_ids)),
        _ctypes.byref(_ctypes.c_int(family_id_len)),
        family_ids,
        _ctypes.byref(_ctypes.c_int(n_family_centroids_rows)),
        _ctypes.byref(_ctypes.c_int(n_family_centroids_cols)),
        family_centroids,
        _ctypes.byref(_ctypes.c_int(n_shift_vectors_rows)),
        _ctypes.byref(_ctypes.c_int(n_shift_vectors_cols)),
        shift_vectors,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _READ_TOX_DATA_INTO_ARGUMENTS, _READ_TOX_DATA_INTO_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    expression.flags.writeable = False
    gene_to_family.flags.writeable = False
    family_centroids.flags.writeable = False
    shift_vectors.flags.writeable = False

    return {
        "gene_ids": [_s.decode().rstrip(' ') for _s in gene_ids],
        "expression": expression,
        "gene_to_family": gene_to_family,
        "family_ids": [_s.decode().rstrip(' ') for _s in family_ids],
        "family_centroids": family_centroids,
        "shift_vectors": shift_vectors,
    }
