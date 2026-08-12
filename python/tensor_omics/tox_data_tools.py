r"""tox_data_tools

Parsers for the plain-text input formats TensorOmics data sets are built from (gene-expression
TSV/CSV files, OrthoFinder-style family files), plus small array-filtering helpers.

These populate the raw `gene_ids` / `expression` / `gene_to_family` / `family_ids` arrays that
:func:`tensor_omics.save_tox_data` later persists; unlike the archive/serde
layers, everything here works from delimited text rather than the library's binary array format.

Python binding, generated from tox_data_tools. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.read_expression_vectors_tsv_c.restype = None
_lib.read_expression_vectors_tsv_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
)

#: The wrapped procedure's arguments, so an error can name one
_READ_EXPRESSION_VECTORS_TSV_ARGUMENTS = ("file_list", "gene_ids", "expression_vectors", "n_header_rows", "gene_col", "value_cols", "start_row", "ierr", "delimiter",)

_lib.read_gene_ids_from_tsv_file_c.restype = None
_lib.read_gene_ids_from_tsv_file_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_READ_GENE_IDS_FROM_TSV_FILE_ARGUMENTS = ("filename", "gene_ids", "n_header_rows", "gene_col", "ierr",)

_lib.read_orthofinder_file_c.restype = None
_lib.read_orthofinder_file_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_READ_ORTHOFINDER_FILE_ARGUMENTS = ("filename", "gene_ids", "family_ids", "gene_to_fam", "ierr",)

_lib.get_unassigned_mask_c.restype = None
_lib.get_unassigned_mask_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GET_UNASSIGNED_MASK_ARGUMENTS = ("gene_to_fam", "mask", "n_genes_kept",)

def read_expression_vectors_tsv(
        file_list,
        gene_ids,
        expression_vectors,
        n_header_rows,
        gene_col,
        value_cols,
        start_row,
        delimiter='\t',
):
    r"""Read expression vectors from csv/tsv files

    Parameters
    ----------
    file_list : sequence of str, of length n_file_list_elements
        List of files to read from
    gene_ids : sequence of str, of length n_gene_ids_elements
        Array of gene IDS
    expression_vectors : np.ndarray[np.float64] of shape (n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2,), column-major (order='F'), modified in place
        Array of expression vectors
    n_header_rows : int
        Number of header rows to skip
    gene_col : int
        Index of column with gene_ids
    value_cols : np.ndarray[np.int32] of shape (n_value_cols_elements,)
        Indicies of columns containing values
    start_row : int
        Row in the expression vectors to start in
    delimiter : str, optional, default '\t'
        optional delimiter
        The default value is `char(9)`.

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_tools::read_expression_vectors_tsv`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        file_list = np.asarray([str(_s).encode() for _s in file_list], dtype="S")
    except TypeError as error:
        raise TypeError(f"'file_list' must be a sequence of strings: {error}") from None
    if file_list.ndim != 1:
        raise ValueError(f"'file_list' must have 1 dimension, but has {file_list.ndim}")
    try:
        gene_ids = np.asarray([str(_s).encode() for _s in gene_ids], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_ids' must be a sequence of strings: {error}") from None
    if gene_ids.ndim != 1:
        raise ValueError(f"'gene_ids' must have 1 dimension, but has {gene_ids.ndim}")
    if not isinstance(expression_vectors, np.ndarray) or expression_vectors.dtype != np.float64:
        raise TypeError("'expression_vectors' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")
    if not expression_vectors.flags.f_contiguous:
        raise ValueError("'expression_vectors' is modified in place, so it must already be column-major (order='F')")
    try:
        value_cols = np.ascontiguousarray(value_cols, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'value_cols' must be an array of np.int32: {error}") from None
    if value_cols.ndim != 1:
        raise ValueError(f"'value_cols' must have 1 dimension, but has {value_cols.ndim}")
    delimiter = np.array([str(delimiter).encode()], dtype="S1")

    # what the inputs already say, rather than asking for it again
    file_list_strlen = file_list.itemsize
    n_file_list_elements = file_list.shape[0]
    gene_ids_strlen = gene_ids.itemsize
    n_gene_ids_elements = gene_ids.shape[0]
    n_expression_vectors_elements_dim_1 = expression_vectors.shape[0]
    n_expression_vectors_elements_dim_2 = expression_vectors.shape[1]
    n_value_cols_elements = value_cols.shape[0]

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.read_expression_vectors_tsv_c(
        file_list,
        ctypes.byref(ctypes.c_int(file_list_strlen)),
        ctypes.byref(ctypes.c_int(n_file_list_elements)),
        gene_ids,
        ctypes.byref(ctypes.c_int(gene_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_gene_ids_elements)),
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_2)),
        ctypes.byref(ctypes.c_int(n_header_rows)),
        ctypes.byref(ctypes.c_int(gene_col)),
        value_cols,
        ctypes.byref(ctypes.c_int(n_value_cols_elements)),
        ctypes.byref(ctypes.c_int(start_row)),
        ctypes.byref(ierr),
        delimiter,
    )

    check_err_code(ierr.value, _READ_EXPRESSION_VECTORS_TSV_ARGUMENTS)

    return None

def read_gene_ids_from_tsv_file(
        filename,
        gene_ids_strlen,
        n_gene_ids_elements,
        n_header_rows,
        gene_col,
):
    r"""Only read the gene ids from a tsv file

    Parameters
    ----------
    filename : str
        Name of the file
    gene_ids_strlen : int
        length of the strings in `gene_ids`
    n_gene_ids_elements : int
        number of elements in `gene_ids`
    n_header_rows : int
        number of headers to skip
    gene_col : int
        Index of the column containing gene ids

    Returns
    -------
    gene_ids : sequence of str, of length n_gene_ids_elements
        gene ids array

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_tools::read_gene_ids_from_tsv_file`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    filename = np.array([str(filename).encode()], dtype="S")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize

    # outputs and work arrays, which the caller never sees
    gene_ids = np.zeros((n_gene_ids_elements,), dtype=f"S{gene_ids_strlen}")
    ierr = ctypes.c_int(0)

    _lib.read_gene_ids_from_tsv_file_c(
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        gene_ids,
        ctypes.byref(ctypes.c_int(gene_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_gene_ids_elements)),
        ctypes.byref(ctypes.c_int(n_header_rows)),
        ctypes.byref(ctypes.c_int(gene_col)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _READ_GENE_IDS_FROM_TSV_FILE_ARGUMENTS)

    return [_s.decode() for _s in gene_ids]

def read_orthofinder_file(
        filename,
        gene_ids,
        family_ids_strlen,
        n_family_ids_elements,
        n_gene_to_fam_elements,
):
    r"""Read a family file (Orthofinder)

    Parameters
    ----------
    filename : str
        Name of the file
    gene_ids : sequence of str, of length n_gene_ids_elements
        gene ids array
    family_ids_strlen : int
        length of the strings in `family_ids`
    n_family_ids_elements : int
        number of elements in `family_ids`
    n_gene_to_fam_elements : int
        number of elements in `gene_to_fam`

    Returns
    -------
    dict
        with keys:

        family_ids : sequence of str, of length n_family_ids_elements
            family ids array
        gene_to_fam : np.ndarray[np.int32] of shape (n_gene_to_fam_elements,), read-only
            gene to family mapping
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_tools::read_orthofinder_file`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    filename = np.array([str(filename).encode()], dtype="S")
    try:
        gene_ids = np.asarray([str(_s).encode() for _s in gene_ids], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_ids' must be a sequence of strings: {error}") from None
    if gene_ids.ndim != 1:
        raise ValueError(f"'gene_ids' must have 1 dimension, but has {gene_ids.ndim}")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize
    gene_ids_strlen = gene_ids.itemsize
    n_gene_ids_elements = gene_ids.shape[0]

    # outputs and work arrays, which the caller never sees
    family_ids = np.zeros((n_family_ids_elements,), dtype=f"S{family_ids_strlen}")
    gene_to_fam = np.empty((n_gene_to_fam_elements,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.read_orthofinder_file_c(
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        gene_ids,
        ctypes.byref(ctypes.c_int(gene_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_gene_ids_elements)),
        family_ids,
        ctypes.byref(ctypes.c_int(family_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_family_ids_elements)),
        gene_to_fam,
        ctypes.byref(ctypes.c_int(n_gene_to_fam_elements)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _READ_ORTHOFINDER_FILE_ARGUMENTS)

    # a result is a value: modify a copy, not this
    gene_to_fam.flags.writeable = False

    return {
        "family_ids": [_s.decode() for _s in family_ids],
        "gene_to_fam": gene_to_fam,
    }

def get_unassigned_mask(
        gene_to_fam,
):
    r"""Helper to create a mask of genes that are unassigned

    Parameters
    ----------
    gene_to_fam : np.ndarray[np.int32] of shape (n_gene_to_fam_elements,)
        gene to family mapping

    Returns
    -------
    dict
        with keys:

        mask : np.ndarray[np.bool_] of shape (size(gene_to_fam),), read-only
            mask for mapping
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_genes_kept : int
            number of genes kept

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_tools::get_unassigned_mask`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_gene_to_fam_elements = gene_to_fam.shape[0]

    # outputs and work arrays, which the caller never sees
    mask = np.empty((gene_to_fam.size,), dtype=np.bool_, order='C')
    n_genes_kept = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.get_unassigned_mask_c(
        gene_to_fam,
        ctypes.byref(ctypes.c_int(n_gene_to_fam_elements)),
        mask,
        ctypes.byref(n_genes_kept),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GET_UNASSIGNED_MASK_ARGUMENTS)

    # a result is a value: modify a copy, not this
    mask.flags.writeable = False

    return {
        "mask": mask,
        "n_genes_kept": n_genes_kept.value,
    }
