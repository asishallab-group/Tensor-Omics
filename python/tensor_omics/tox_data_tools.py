r"""tox_data_tools

Parsers for the plain-text input formats TensorOmics data sets are built from (gene-expression
TSV/CSV files, OrthoFinder-style family files), plus small array-filtering helpers.

These populate the raw `gene_ids` / `expression` / `gene_to_family` / `family_ids` arrays that
:func:`tensor_omics.save_tox_data` later persists; unlike the archive/serde
layers, everything here works from delimited text rather than the library's binary array format.

The readers are deliberately simple, not a general CSV reader: a line is read whole up to
4096 characters and refused with `ERR_STRING_TOO_LONG` past that; each delimiter separates two
fields, so two in a row enclose an empty one; blanks around a field are not part of it; there
is no quoting. Whatever a reader cannot read as it is meant, it refuses rather than guesses.
After an error, the contents of every output are undefined.

Python binding, generated from tox_data_tools. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.read_expression_vectors_tsv_c.restype = None
_lib.read_expression_vectors_tsv_c.argtypes = (
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_READ_EXPRESSION_VECTORS_TSV_ARGUMENTS = ("file_list", "gene_ids", "expression_vectors", "n_header_rows", "gene_col", "value_cols", "start_row", "delimiter", "ierr",)

_lib.read_gene_ids_from_tsv_file_c.restype = None
_lib.read_gene_ids_from_tsv_file_c.argtypes = (
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_READ_GENE_IDS_FROM_TSV_FILE_ARGUMENTS = ("filename", "gene_ids", "n_header_rows", "gene_col", "ierr",)

_lib.read_orthofinder_file_c.restype = None
_lib.read_orthofinder_file_c.argtypes = (
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_READ_ORTHOFINDER_FILE_ARGUMENTS = ("filename", "gene_ids", "family_ids", "gene_to_fam", "ierr",)

_lib.get_unassigned_mask_c.restype = None
_lib.get_unassigned_mask_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
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

    Each file's rows give the values of the gene their `gene_col` names: the `value_cols` of a
    row go to rows `start_row`, `start_row + 1`, ... of that gene's column, and each further file
    continues below the last. A row naming a gene not in `gene_ids` is skipped. Every row must
    have every column read, none of them empty; a value must be one complete, finite number.
    Blank lines at the end of a file are ignored. A line may be 4096 characters long.
    Values of genes no row names are left as they were; after an error, `expression_vectors`
    is undefined.

    Parameters
    ----------
    file_list : sequence of str, of length n_file_list_elements
        List of files to read from
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    gene_ids : sequence of str, of length n_gene_ids_elements
        Array of gene IDS
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
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
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
        At most 1 byte long in UTF-8: a longer value is refused, not cut short.

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
    Generated from the Fortran procedure `tox_data_tools::read_expression_vectors_tsv`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    _file_list_bytes = _contract.strings(file_list, "file_list")
    _file_list_width = max(map(len, _file_list_bytes), default=0) or 1
    file_list = _np.asarray([_b.ljust(_file_list_width) for _b in _file_list_bytes], dtype="S")
    _gene_ids_bytes = _contract.strings(gene_ids, "gene_ids")
    _gene_ids_width = max(map(len, _gene_ids_bytes), default=0) or 1
    gene_ids = _np.asarray([_b.ljust(_gene_ids_width) for _b in _gene_ids_bytes], dtype="S")
    if not isinstance(expression_vectors, _np.ndarray) or expression_vectors.dtype != _np.float64:
        raise _ToxArgumentError("'expression_vectors' is modified in place, so it must already be a numpy array of {}".format(_np.float64), argument="expression_vectors")
    if expression_vectors.ndim != 2:
        raise _ToxArgumentError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}", argument="expression_vectors")
    if not expression_vectors.flags.f_contiguous:
        raise _ToxArgumentError("'expression_vectors' is modified in place, so it must already be column-major (order='F')", argument="expression_vectors")
    n_header_rows = _contract.integer(n_header_rows, "n_header_rows")
    gene_col = _contract.integer(gene_col, "gene_col")
    value_cols = _contract.integers(value_cols, "value_cols", ndim=1)
    start_row = _contract.integer(start_row, "start_row")
    delimiter = _np.array([_contract.string(delimiter, "delimiter", 1).ljust(1)], dtype="S1")

    # what the inputs already say, rather than asking for it again
    file_list_strlen = file_list.itemsize
    n_file_list_elements = file_list.shape[0]
    gene_ids_strlen = gene_ids.itemsize
    n_gene_ids_elements = gene_ids.shape[0]
    n_expression_vectors_elements_dim_1 = expression_vectors.shape[0]
    n_expression_vectors_elements_dim_2 = expression_vectors.shape[1]
    n_value_cols_elements = value_cols.shape[0]

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.read_expression_vectors_tsv_c(
        file_list,
        _ctypes.byref(_ctypes.c_int(file_list_strlen)),
        _ctypes.byref(_ctypes.c_int(n_file_list_elements)),
        gene_ids,
        _ctypes.byref(_ctypes.c_int(gene_ids_strlen)),
        _ctypes.byref(_ctypes.c_int(n_gene_ids_elements)),
        expression_vectors,
        _ctypes.byref(_ctypes.c_int(n_expression_vectors_elements_dim_1)),
        _ctypes.byref(_ctypes.c_int(n_expression_vectors_elements_dim_2)),
        _ctypes.byref(_ctypes.c_int(n_header_rows)),
        _ctypes.byref(_ctypes.c_int(gene_col)),
        value_cols,
        _ctypes.byref(_ctypes.c_int(n_value_cols_elements)),
        _ctypes.byref(_ctypes.c_int(start_row)),
        delimiter,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _READ_EXPRESSION_VECTORS_TSV_ARGUMENTS)

    return None

def read_gene_ids_from_tsv_file(
        filename,
        gene_ids_strlen,
        n_gene_ids_elements,
        n_header_rows,
        gene_col,
):
    r"""Only read the gene ids from a tsv file

    Each tab-separated row after the header gives one id, from its `gene_col`, in file order.
    A row without that column, or with it empty, is refused, and so is an id longer than the
    strings of `gene_ids`. Blank lines at the end of the file are ignored. A line may be 4096
    characters long. With fewer rows than `gene_ids` has elements, the ones past the last row
    are left undefined, and the number of rows read is not returned: size `gene_ids` to the
    file. After an error, `gene_ids` is undefined.

    Parameters
    ----------
    filename : str
        Name of the file
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_tools::read_gene_ids_from_tsv_file`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    filename = _np.array([_contract.string(filename, "filename").ljust(1)], dtype="S")
    gene_ids_strlen = _contract.integer(gene_ids_strlen, "gene_ids_strlen")
    n_gene_ids_elements = _contract.integer(n_gene_ids_elements, "n_gene_ids_elements")
    n_header_rows = _contract.integer(n_header_rows, "n_header_rows")
    gene_col = _contract.integer(gene_col, "gene_col")

    # the outputs are allocated by a size the caller passes, which must not be negative
    _contract.extent(gene_ids_strlen, "gene_ids_strlen")
    _contract.extent(n_gene_ids_elements, "n_gene_ids_elements")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize

    # outputs and work arrays, which the caller never sees
    gene_ids = _np.zeros((n_gene_ids_elements,), dtype=f"S{gene_ids_strlen}")
    ierr = _ctypes.c_int(0)

    _lib.read_gene_ids_from_tsv_file_c(
        filename,
        _ctypes.byref(_ctypes.c_int(filename_strlen)),
        gene_ids,
        _ctypes.byref(_ctypes.c_int(gene_ids_strlen)),
        _ctypes.byref(_ctypes.c_int(n_gene_ids_elements)),
        _ctypes.byref(_ctypes.c_int(n_header_rows)),
        _ctypes.byref(_ctypes.c_int(gene_col)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _READ_GENE_IDS_FROM_TSV_FILE_ARGUMENTS)

    return [_s.decode().rstrip(' ') for _s in gene_ids]

def read_orthofinder_file(
        filename,
        gene_ids,
        family_ids_strlen,
        n_family_ids_elements,
        n_gene_to_fam_elements,
):
    r"""Read a family file (Orthofinder)

    After its header line, each line of an OrthoFinder `Orthogroups.tsv` is a family: its id,
    then a tab-separated column per species, each a comma-separated list of genes (empty
    columns allowed). Family k is the k-th such line, and each gene of `gene_ids` gets the
    index of the family that lists it, or 0 when none does; listed genes
    not in `gene_ids` are ignored. Blank lines are skipped. A family id longer than the strings
    of `family_ids`, a line without one, and more families than `family_ids` holds are refused.
    A line may be 4096 characters long. With fewer families than `family_ids` has elements, the
    ones past the last family are left undefined; after an error, all outputs are undefined.

    Parameters
    ----------
    filename : str
        Name of the file
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    gene_ids : sequence of str, of length n_gene_ids_elements
        gene ids array
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_tools::read_orthofinder_file`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    filename = _np.array([_contract.string(filename, "filename").ljust(1)], dtype="S")
    _gene_ids_bytes = _contract.strings(gene_ids, "gene_ids")
    _gene_ids_width = max(map(len, _gene_ids_bytes), default=0) or 1
    gene_ids = _np.asarray([_b.ljust(_gene_ids_width) for _b in _gene_ids_bytes], dtype="S")
    family_ids_strlen = _contract.integer(family_ids_strlen, "family_ids_strlen")
    n_family_ids_elements = _contract.integer(n_family_ids_elements, "n_family_ids_elements")
    n_gene_to_fam_elements = _contract.integer(n_gene_to_fam_elements, "n_gene_to_fam_elements")

    # the outputs are allocated by a size the caller passes, which must not be negative
    _contract.extent(family_ids_strlen, "family_ids_strlen")
    _contract.extent(n_family_ids_elements, "n_family_ids_elements")
    _contract.extent(n_gene_to_fam_elements, "n_gene_to_fam_elements")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize
    gene_ids_strlen = gene_ids.itemsize
    n_gene_ids_elements = gene_ids.shape[0]

    # outputs and work arrays, which the caller never sees
    family_ids = _np.zeros((n_family_ids_elements,), dtype=f"S{family_ids_strlen}")
    gene_to_fam = _np.empty((n_gene_to_fam_elements,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.read_orthofinder_file_c(
        filename,
        _ctypes.byref(_ctypes.c_int(filename_strlen)),
        gene_ids,
        _ctypes.byref(_ctypes.c_int(gene_ids_strlen)),
        _ctypes.byref(_ctypes.c_int(n_gene_ids_elements)),
        family_ids,
        _ctypes.byref(_ctypes.c_int(family_ids_strlen)),
        _ctypes.byref(_ctypes.c_int(n_family_ids_elements)),
        gene_to_fam,
        _ctypes.byref(_ctypes.c_int(n_gene_to_fam_elements)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _READ_ORTHOFINDER_FILE_ARGUMENTS)

    # a result is a value: modify a copy, not this
    gene_to_fam.flags.writeable = False

    return {
        "family_ids": [_s.decode().rstrip(' ') for _s in family_ids],
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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_tools::get_unassigned_mask`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    gene_to_fam = _contract.integers(gene_to_fam, "gene_to_fam", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_gene_to_fam_elements = gene_to_fam.shape[0]

    # outputs and work arrays, which the caller never sees
    mask = _np.empty((gene_to_fam.size,), dtype=_np.bool_, order='C')
    n_genes_kept = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.get_unassigned_mask_c(
        gene_to_fam,
        _ctypes.byref(_ctypes.c_int(n_gene_to_fam_elements)),
        mask,
        _ctypes.byref(n_genes_kept),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _GET_UNASSIGNED_MASK_ARGUMENTS)

    # a result is a value: modify a copy, not this
    mask.flags.writeable = False

    return {
        "mask": mask,
        "n_genes_kept": n_genes_kept.value,
    }
