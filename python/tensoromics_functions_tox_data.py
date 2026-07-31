"""
TensorOmics Functions Module
Python wrapper functions for Fortran routines via C interface
"""
from error_handling import check_err_code
import numpy as np
import ctypes
import os
from tensoromics_functions import (
    tox_deserialize_char_nd,
    tox_serialize_char_nd,
    tox_serialize_int_nd,
    tox_deserialize_int_nd,
    tox_serialize_real_nd,
    tox_deserialize_real_nd,
    _strings_to_c_char_matrix,
    _c_char_matrix_to_strings,
    _create_empty_c_char_matrix
)

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL(dll_path)


#> tox_data_read_write:read_gene_ids_from_tsv_file_c: Read gene ids from a tsv file
#' Function for read_gene_ids_from_tsv_file_c
def read_gene_ids_from_tsv_file(filename,n_genes, gene_ids_len, n_header_rows, gene_col):
    """
    Read gene ids from a tsv file

    Args:
        filename: Name of the TSV file
        n_genes: Number of genes to read
        gene_ids_len: Maximum length of each gene ID
        n_header_rows: Number of header rows to skip
        gene_col: Column index (1-based) for gene IDs
    """
    # Ensure filename is a string (single file)
    if isinstance(filename, list):
        if len(filename) > 0:
            filename = filename[0]
        else:
            raise ValueError("filename cannot be an empty list for read_gene_ids_from_tsv_file")

    gene_ids = _create_empty_c_char_matrix(n_genes, gene_ids_len)
    ierr = ctypes.c_int()

    # Example for read_gene_ids_from_tsv_file_c
    lib.read_gene_ids_from_tsv_file_c.argtypes = [
        ctypes.c_char_p,  # filename_raw
        ctypes.POINTER(ctypes.c_int),                  # fn_len
        np.ctypeslib.ndpointer(flags="F_CONTIGUOUS"),  # gene_ids_raw
        ctypes.POINTER(ctypes.c_int),                  # gene_ids_len
        ctypes.POINTER(ctypes.c_int),                  # n_genes
        ctypes.POINTER(ctypes.c_int),                  # n_header_rows
        ctypes.POINTER(ctypes.c_int),                  # gene_col
        ctypes.POINTER(ctypes.c_int)   # ierr
    ]
    lib.read_gene_ids_from_tsv_file_c.restype = None

    # Call C function
    lib.read_gene_ids_from_tsv_file_c(
        filename.encode("utf-8"),
        ctypes.byref(ctypes.c_int(len(filename))),
        gene_ids,
        ctypes.byref(ctypes.c_int(gene_ids_len)),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_header_rows)),
        ctypes.byref(ctypes.c_int(gene_col)),
        ctypes.byref(ierr)
    )

    check_err_code(ierr.value)

    return _c_char_matrix_to_strings(gene_ids, gene_ids_len)

#> tox_data_read_write:read_expression_vectors_tsv_c: Read expression vectors from given tabular (csv/tsv) files
# Function for read_expression_vectors_tsv_c
def read_expression_vectors_tsv(file_list, gene_ids, n_samples, n_header_rows,
                           gene_col, value_cols, delimiter='\t'):
    """
    Read expression vectors from given tabular (csv/tsv) files
    Args:
        file_list: List of filenames to read
        gene_ids: List of gene IDs to extract
        n_samples: Number of samples (files)
        n_header_rows: Number of header rows to skip in each file
        gene_col: Column index (1-based) for gene IDs
        value_cols: List of column indices (1-based) for expression values
        delimiter: Delimiter used in the files (default: tab)
    """
    # Ensure file_list is a list (multiple files)
    if not isinstance(file_list, list):
        file_list = [file_list]  # Convert single file to list

    # Ensure inputs are numpy arrays

    # Convert inputs to c_char matrices
    file_list_matrix, max_file_len = _strings_to_c_char_matrix(file_list)
    gene_ids_matrix, max_gene_len = _strings_to_c_char_matrix(gene_ids)

    # Prepare output arrays
    expression_vectors = np.empty((n_samples, len(gene_ids)), dtype=np.float64, order='F')
    ierr = ctypes.c_int()

    value_cols_ct = np.ascontiguousarray(value_cols, dtype=np.int32)

    # read_expression_vectors_tsv_c
    lib.read_expression_vectors_tsv_c.argtypes = [
        np.ctypeslib.ndpointer(flags="F_CONTIGUOUS"),  # file_list_raw
        ctypes.POINTER(ctypes.c_int),                  # file_list_len
        ctypes.POINTER(ctypes.c_int),                  # n_files
        np.ctypeslib.ndpointer(flags="F_CONTIGUOUS"),  # gene_ids_raw
        ctypes.POINTER(ctypes.c_int),                  # gene_ids_len
        ctypes.POINTER(ctypes.c_int),                  # n_genes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # expression_vectors_flat
        ctypes.POINTER(ctypes.c_int),                  # n_samples
        ctypes.POINTER(ctypes.c_int),                  # n_header_rows
        ctypes.POINTER(ctypes.c_int),                  # gene_col
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),  # value_cols
        ctypes.POINTER(ctypes.c_int),                  # n_value_cols
        ctypes.POINTER(ctypes.c_int),  # ierr
        ctypes.c_char_p  # delimiter_raw
    ]
    lib.read_expression_vectors_tsv_c.restype = None

    # Call C function
    lib.read_expression_vectors_tsv_c(
        file_list_matrix,
        ctypes.byref(ctypes.c_int(max_file_len)),
        ctypes.byref(ctypes.c_int(len(file_list))),
        gene_ids_matrix,
        ctypes.byref(ctypes.c_int(max_gene_len)),
        ctypes.byref(ctypes.c_int(len(gene_ids))),
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_header_rows)),
        ctypes.byref(ctypes.c_int(gene_col)),
        value_cols_ct,
        ctypes.byref(ctypes.c_int(len(value_cols))),
        ctypes.byref(ierr),
        delimiter.encode("utf-8")
    )

    check_err_code(ierr.value)

    return expression_vectors

#> tox_data_read_write:read_orthofinder_file_c: Read an orthofinder family file and map genes to families
# Function for read_orthofinder_file_c
def read_orthofinder_file(filename, gene_ids, family_ids_len, n_families):
    """
    Read an orthofinder family file and map genes to families
    Args:
        filename: Name of the orthofinder TSV file
        gene_ids: List of gene IDs to map
        family_ids_len: Maximum length of each family ID
        n_families: Number of families to read
    """
    # Ensure filename is a string (single file)
    if isinstance(filename, list):
        if len(filename) > 0:
            filename = filename[0]  # Take first element if it's a list
        else:
            raise ValueError("filename cannot be an empty list for read_orthofinder_file")

    # Ensure inputs are numpy arrays

    gene_ids_matrix, max_gene_len = _strings_to_c_char_matrix(gene_ids)

    # Prepare output arrays
    family_ids = _create_empty_c_char_matrix(n_families, family_ids_len)
    gene_to_fam = np.empty(len(gene_ids), dtype=np.int32, order='F')
    ierr = ctypes.c_int()

    # read_orthofinder_file_c
    lib.read_orthofinder_file_c.argtypes = [
        ctypes.c_char_p,  # filename_raw
        ctypes.POINTER(ctypes.c_int),                  # fn_len
        np.ctypeslib.ndpointer(flags="F_CONTIGUOUS"),  # gene_ids_raw
        ctypes.POINTER(ctypes.c_int),                  # gene_ids_len
        ctypes.POINTER(ctypes.c_int),                  # n_genes
        np.ctypeslib.ndpointer(flags="F_CONTIGUOUS"),  # family_ids_raw
        ctypes.POINTER(ctypes.c_int),                  # family_ids_len
        ctypes.POINTER(ctypes.c_int),                  # n_families
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),  # gene_to_fam
        ctypes.POINTER(ctypes.c_int)   # ierr
    ]
    lib.read_orthofinder_file_c.restype = None
    # Call C function
    lib.read_orthofinder_file_c(
        filename.encode("utf-8"),
        ctypes.byref(ctypes.c_int(len(filename))),
        gene_ids_matrix,
        ctypes.byref(ctypes.c_int(max_gene_len)),
        ctypes.byref(ctypes.c_int(len(gene_ids))),
        family_ids,
        ctypes.byref(ctypes.c_int(family_ids_len)),
        ctypes.byref(ctypes.c_int(n_families)),
        gene_to_fam,
        ctypes.byref(ierr)
    )

    check_err_code(ierr.value)

    return {
        'family_ids': _c_char_matrix_to_strings(family_ids, family_ids_len),
        'gene_to_fam': gene_to_fam
    }

#> tox_data_tools:filter_unassigned_genes_c: Filter out genes that are not assigned to any family (where gene_to_fam == 0).
def filter_unassigned_genes(gene_to_fam):
    """
    Filter out genes that are not assigned to any family (where gene_to_fam == 0).

    Args:
        gene_to_fam (list[int] or np.ndarray): Family assignment for each gene (0 means unassigned).

    Returns:
        dict: {
            'mask': list[int] of 1s (kept) and 0s (removed),
            'n_genes_kept': int, number of genes kept
        }
    """
    # Convert to numpy arrays for convenience
    gene_to_fam = np.array(gene_to_fam, dtype=int)

    # Logical mask: 1 if gene_to_fam != 0, else 0
    mask = (gene_to_fam != 0).astype(int)

    return {
        'mask': mask.tolist(),
        'n_genes_kept': int(np.sum(mask))
    }

# --- Python wrappers for validation ---

#> tox_data_validation:validate_gene_to_family_mapping_c: Validate gene to family mapping
def validate_gene_to_family_mapping(gene_to_fam, n_families):
    """
    Validate gene to family mapping
    Args:
        gene_to_fam: Array mapping each gene to a family index (0 if unassigned)
        n_families: Total number of families
    """
    n_genes = gene_to_fam.size
    ierr = ctypes.c_int()

    lib.validate_gene_to_family_mapping_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"), # gene_to_fam
        ctypes.POINTER(ctypes.c_int),                 # n_genes
        ctypes.POINTER(ctypes.c_int),                 # n_families
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    lib.validate_gene_to_family_mapping_c.restype = None

    lib.validate_gene_to_family_mapping_c(
        gene_to_fam,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

#> tox_data_validation:validate_expression_data_c: Validate expression data
def validate_expression_data(expression_vectors, check_non_negative=True):
    """
    Validate expression data
    Args:
        expression_vectors: 2D array of expression data (samples x genes)
        check_non_negative: Whether to check for non-negative values
    """
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    n_samples, n_genes = expression_vectors.shape
    ierr = ctypes.c_int()

    lib.validate_expression_data_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # expression_vectors
        ctypes.POINTER(ctypes.c_int),                    # n_genes
        ctypes.POINTER(ctypes.c_int),                    # n_samples
        ctypes.POINTER(ctypes.c_int),                    # check_non_negative
        ctypes.POINTER(ctypes.c_int)     # ierr
    ]
    lib.validate_expression_data_c.restype = None

    lib.validate_expression_data_c(
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(1 if check_non_negative else 0)),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

#> tox_data_validation:validate_family_centroids_c: Validate family centroids, checks for NaN/Inf
def validate_family_centroids(family_centroids):
    """
    Validate family centroids, checks for NaN/Inf
    Args:
        family_centroids: 2D array of family centroids (samples x families)
    """
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    n_samples, n_families = family_centroids.shape
    ierr = ctypes.c_int()

    lib.validate_family_centroids_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # family_centroids
        ctypes.POINTER(ctypes.c_int),                    # n_families
        ctypes.POINTER(ctypes.c_int),                    # n_samples
        ctypes.POINTER(ctypes.c_int)     # ierr
    ]
    lib.validate_family_centroids_c.restype = None

    lib.validate_family_centroids_c(
        family_centroids,
        ctypes.byref(ctypes.c_int(n_families)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

#> tox_data_validation:validate_shift_vectors_c: Validate shift vectors, checks if datatypes are correct and if the general structure matches
def validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_genes, n_samples, n_families):
    """
    Validate shift vectors, checks if datatypes are correct and if the general structure matches (first d rows = centroids, d+1 to 2d rows = shift)
    Args:
        shift_vectors: 2D array of shift vectors (2*samples x genes)
        expression_vectors: 2D array of expression data (samples x genes)
        family_centroids: 2D array of family centroids (samples x families)
        gene_to_fam: Array mapping each gene to a family index (0 if unassigned)
        n_genes: Number of genes
        n_samples: Number of samples
        n_families: Number of families
    """

    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    ierr = ctypes.c_int()

    lib.validate_shift_vectors_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # shift_vectors
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # expression_vectors
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # family_centroids
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),    # gene_to_fam
        ctypes.POINTER(ctypes.c_int),                    # n_genes
        ctypes.POINTER(ctypes.c_int),                    # n_samples
        ctypes.POINTER(ctypes.c_int),                    # n_families
        ctypes.POINTER(ctypes.c_int)     # ierr
    ]
    lib.validate_shift_vectors_c.restype = None

    lib.validate_shift_vectors_c(
        shift_vectors,
        expression_vectors,
        family_centroids,
        gene_to_fam,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_families)),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

#> tox_data_validation:validate_string_array_uniqueness_c: Validate uniqueness of strings
def validate_string_array_uniqueness(strings):
    """
    Validate uniqueness of strings - Note: Uses hashset internally which may increase memory usage temporarily for large datasets
    Args:
        strings: List of strings
    """
    n_strings = len(strings)
    gene_ids_raw, gene_ids_len = _strings_to_c_char_matrix(strings)
    ierr = ctypes.c_int()

    lib.validate_string_array_uniqueness_c.argtypes = [
        np.ctypeslib.ndpointer(flags="F_CONTIGUOUS"), # string_array_raw
        ctypes.POINTER(ctypes.c_int),                 # string_len
        ctypes.POINTER(ctypes.c_int),                 # n_strings
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    lib.validate_string_array_uniqueness_c.restype = None

    lib.validate_string_array_uniqueness_c(
        gene_ids_raw,
        ctypes.byref(ctypes.c_int(gene_ids_len)),
        ctypes.byref(ctypes.c_int(n_strings)),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

#> tox_data_validation:validate_data_structure_c: Validate overall data structure consistency. Confirms sizes and dependencies as far as possible.
def validate_data_structure(n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors):
    """
    Validate overall data structure consistency. Confirms sizes and dependencies as far as possible.
    Args:
        n_genes: Number of genes
        n_families: Number of families
        n_samples: Number of samples
        gene_ids: List of gene IDs
        gene_family_ids: List of family IDs
        gene_to_fam: Array mapping each gene to a family index (0 if unassigned)
        expression_vectors: 2D array of expression data (samples x genes)
        family_centroids: 2D array of family centroids (samples x families)
        shift_vectors: 2D array of shift vectors (2*samples x genes)
    """

    gene_ids_raw, gene_ids_len = _strings_to_c_char_matrix(gene_ids)
    gene_family_ids_raw, fam_len = _strings_to_c_char_matrix(gene_family_ids)

    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    ierr = ctypes.c_int()

    lib.validate_data_structure_c.argtypes = [
        ctypes.POINTER(ctypes.c_int),                 # n_genes
        ctypes.POINTER(ctypes.c_int),                 # n_families
        ctypes.POINTER(ctypes.c_int),                 # n_samples
        np.ctypeslib.ndpointer(flags="F_CONTIGUOUS"), # gene_ids_raw
        ctypes.POINTER(ctypes.c_int),                 # gene_ids_len
        np.ctypeslib.ndpointer(flags="F_CONTIGUOUS"), # gene_family_ids_raw
        ctypes.POINTER(ctypes.c_int),                 # fam_len
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"), # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # expression_vectors
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # family_centroids
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # shift_vectors
        ctypes.POINTER(ctypes.c_int)     # ierr
    ]
    lib.validate_data_structure_c.restype = None

    lib.validate_data_structure_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        ctypes.byref(ctypes.c_int(n_samples)),
        gene_ids_raw,
        ctypes.byref(ctypes.c_int(gene_ids_len)),
        gene_family_ids_raw,
        ctypes.byref(ctypes.c_int(fam_len)),
        gene_to_fam,
        expression_vectors,
        family_centroids,
        shift_vectors,
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

#> tox_data_validation:validate_all_data_c: Comprehensive validation of all data components. This function performs all individual validations in one go.
def validate_all_data(n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors):
    """
    Comprehensive validation of all data components. This function performs all individual validations in one go.
    Args:
        n_genes: Number of genes
        n_families: Number of families
        n_samples: Number of samples
        gene_ids: List of gene IDs
        gene_family_ids: List of family IDs
        gene_to_fam: Array mapping each gene to a family index (0 if unassigned)
        expression_vectors: 2D array of expression data (samples x genes)
        family_centroids: 2D array of family centroids (samples x families)
        shift_vectors: 2D array of shift vectors (genes x samples)
    """

    gene_ids_raw, gene_ids_len = _strings_to_c_char_matrix(gene_ids)
    gene_family_ids_raw, fam_len = _strings_to_c_char_matrix(gene_family_ids)

    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    ierr = ctypes.c_int()

    lib.validate_all_data_c.argtypes = [
        ctypes.POINTER(ctypes.c_int),                 # n_genes
        ctypes.POINTER(ctypes.c_int),                 # n_families
        ctypes.POINTER(ctypes.c_int),                 # n_samples
        np.ctypeslib.ndpointer(flags="F_CONTIGUOUS"), # gene_ids_raw
        ctypes.POINTER(ctypes.c_int),                 # gene_len
        np.ctypeslib.ndpointer(flags="F_CONTIGUOUS"), # gene_family_ids_raw
        ctypes.POINTER(ctypes.c_int),                 # fam_len
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"), # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # expression_vectors
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # family_centroids
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # shift_vectors
        ctypes.POINTER(ctypes.c_int)     # ierr
    ]
    lib.validate_all_data_c.restype = None

    lib.validate_all_data_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        ctypes.byref(ctypes.c_int(n_samples)),
        gene_ids_raw,
        ctypes.byref(ctypes.c_int(gene_ids_len)),
        gene_family_ids_raw,
        ctypes.byref(ctypes.c_int(fam_len)),
        gene_to_fam,
        expression_vectors,
        family_centroids,
        shift_vectors,
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)


#> tox_data_archive:create_zip_archive_c: Low-level function to create zip archive from keys and filenames.
def create_zip_archive(zip_filename: str, keys, filenames) -> None:
    """
    Low-level function to create zip archive from keys and filenames.
    Directly calls the Fortran function.

    Args:
        zip_filename: Name of the zip file to create
        keys: List of keys for the manifest
        filenames: List of filenames to include in archive
    """
    if len(keys) != len(filenames):
        raise ValueError("Keys and filenames must have the same length")

    # Convert keys and filenames to 2D numpy arrays with Fortran order
    n_keys = len(keys)

    # Create numpy arrays with Fortran order (column-major)
    # Shape: (string_length, string_count) - matching Fortran expectation
    keys_c, max_key_len = _strings_to_c_char_matrix(keys)
    filenames_c, max_filename_len = _strings_to_c_char_matrix(filenames)

    # Set up argument types to match Fortran subroutine
    lib.create_zip_archive_c.argtypes = [
        ctypes.c_char_p, ctypes.POINTER(ctypes.c_int),    # zip_filename, zip_len
        np.ctypeslib.ndpointer(flags='F_CONTIGUOUS'),  # keys
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),     # keys_len, keys_count
        np.ctypeslib.ndpointer(flags='F_CONTIGUOUS'),  # filenames
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),     # filenames_len, filenames_count
        ctypes.POINTER(ctypes.c_int)                                    # ierr
    ]

    ierr = ctypes.c_int()

    # Call Fortran function - arrays are already in Fortran order
    lib.create_zip_archive_c(
        zip_filename.encode("utf-8"),
        ctypes.byref(ctypes.c_int(len(zip_filename))),
        keys_c,
        ctypes.byref(ctypes.c_int(max_key_len)),
        ctypes.byref(ctypes.c_int(n_keys)),
        filenames_c,
        ctypes.byref(ctypes.c_int(max_filename_len)),
        ctypes.byref(ctypes.c_int(n_keys)),
        ctypes.byref(ierr)
    )

    check_err_code(ierr.value)
    print(f"Successfully created archive: {zip_filename}")


#> tox_data_archive:extract_zip_archive_c: Extract a zip archive created by create_zip_archive
def extract_zip_archive(zip_filename):
    """
    Extract a zip archive created by create_zip_archive

    Parameters:
    -----------
    zip_filename : str
        Path to the zip file to extract

    Returns:
    --------
    dict
        Dictionary mapping data keys to extracted filenames

    Raises:
    -------
    RuntimeError
        If extraction fails
    FileNotFoundError
        If zip file doesn't exist
    """

    # Check if zip file exists
    if not os.path.exists(zip_filename):
        raise FileNotFoundError(f"Zip file not found: {zip_filename}")

    lib.extract_zip_archive_c.argtypes = [
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    ]

    try:
        # Initialize error code
        ierr = ctypes.c_int()

        # Call the C extraction function
        lib.extract_zip_archive_c(
            zip_filename.encode("utf-8"),
            ctypes.byref(ctypes.c_int(len(zip_filename))),
            ctypes.byref(ierr)
        )

        check_err_code(ierr.value)

        # Read manifest file to get file mapping
        manifest_path = "manifest.txt"
        if not os.path.exists(manifest_path):
            raise RuntimeError("Manifest file not found after extraction")

        file_mapping = {}
        with open(manifest_path, 'r') as f:
            for line in f:
                parts = line.strip().split('=')
                if len(parts) == 2:
                    file_mapping[parts[0]] = parts[1]

        print(f"Successfully extracted {len(file_mapping)} files from {zip_filename}")

        return file_mapping
    except:
        print("Failed extracting archive")


#> f42_helper: Save tox data to zip archive
def save_tox_data(zip_filename: str,
                 gene_ids = None,
                 expression_vectors = None,
                 gene_to_fam = None,
                 family_ids = None,
                 family_centroids = None,
                 shift_vectors = None,
                 gene_ids_name = None,
                 expression_vectors_name = None,
                 gene_to_fam_name = None,
                 family_ids_name = None,
                 family_centroids_name = None,
                 shift_vectors_name = None) -> None:
    """
    High-level function to save TOX data to zip archive.
    Handles validation, serialization, and calls create_zip_archive.
    Use for standard conform tox gene data.

    Args:
        zip_filename: Name of the zip file to create
        gene_ids, expression_vectors, etc.: Data arrays
        gene_ids_name, expression_vectors_name, etc.: Temporary filenames
    """
    # First serialize the arrays to files
    temp_files = []
    keys = []
    filenames = []

    # Gene IDs
    if gene_ids is not None and gene_ids_name:
        tox_serialize_char_nd(gene_ids, gene_ids_name)
        temp_files.append(gene_ids_name)
        keys.append("gene_ids")
        filenames.append(gene_ids_name)

    # Expression vectors
    if expression_vectors is not None and expression_vectors_name:
        tox_serialize_real_nd(expression_vectors, expression_vectors_name)
        temp_files.append(expression_vectors_name)
        keys.append("expression")
        filenames.append(expression_vectors_name)

    # Gene to family mapping
    if gene_to_fam is not None and gene_to_fam_name:
        tox_serialize_int_nd(gene_to_fam, gene_to_fam_name)
        temp_files.append(gene_to_fam_name)
        keys.append("gene_to_family")
        filenames.append(gene_to_fam_name)

    # Family IDs
    if family_ids is not None and family_ids_name:
        tox_serialize_char_nd(family_ids, family_ids_name)
        temp_files.append(family_ids_name)
        keys.append("family_ids")
        filenames.append(family_ids_name)

    # Family centroids
    if family_centroids is not None and family_centroids_name:
        tox_serialize_real_nd(family_centroids, family_centroids_name)
        temp_files.append(family_centroids_name)
        keys.append("family_centroids")
        filenames.append(family_centroids_name)

    # Shift vectors
    if shift_vectors is not None and shift_vectors_name:
        tox_serialize_real_nd(shift_vectors, shift_vectors_name)
        temp_files.append(shift_vectors_name)
        keys.append("shift_vectors")
        filenames.append(shift_vectors_name)

    # Call the low-level function
    if keys:
        create_zip_archive(zip_filename, keys, filenames)
    else:
        print("No valid data provided to save - skipping archive creation")
        return

    # Clean up temporary files
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)
            print(f"Removed temporary file: {temp_file}")


#> f42_helper: Read tox data from zip archive
def read_tox_data(zip_filename: str,
                 load_gene_ids: bool = False,
                 load_expression_vectors: bool = False,
                 load_gene_to_fam: bool = False,
                 load_family_ids: bool = False,
                 load_family_centroids: bool = False,
                 load_shift_vectors: bool = False):
    """
    Read data from a zip archive created by save_tox_data. Use for standard conform tox gene data.
    Args:
        zip_filename: Name of the zip file to read
        load_gene_ids: Whether to load gene IDs
        load_expression_vectors: Whether to load expression vectors
        load_gene_to_fam: Whether to load gene to family mapping
        load_family_ids: Whether to load family IDs
        load_family_centroids: Whether to load family centroids
        load_shift_vectors: Whether to load shift vectors
    Returns:
        Dictionary with loaded data arrays
    """

    extract_zip_archive(zip_filename)

    result = {
        'gene_ids': None,
        'expression_vectors': None,
        'gene_to_fam': None,
        'family_ids': None,
        'family_centroids': None,
        'shift_vectors': None
    }

    # Read manifest file
    manifest_path = "manifest.txt"
    if not os.path.exists(manifest_path):
        raise FileNotFoundError("Manifest file not found in archive")

    # Parse manifest
    file_mapping = {}
    with open(manifest_path, 'r') as f:
        for line in f:
            parts = line.strip().split('=')
            if len(parts) == 2:
                file_mapping[parts[0]] = parts[1]

    print("Files found in archive:", file_mapping)

    # Load requested data using your existing deserialization functions
    try:
        if load_gene_ids and "gene_ids" in file_mapping:
            filename = file_mapping["gene_ids"]
            if os.path.exists(filename):
                result['gene_ids'] = tox_deserialize_char_nd(filename)
                print(f"Gene ids extracted from {filename}")

        if load_expression_vectors and "expression" in file_mapping:
            filename = file_mapping["expression"]
            if os.path.exists(filename):
                result['expression_vectors'] = tox_deserialize_real_nd(filename)
                print(f"Expression vectors extracted from {filename}")

        if load_gene_to_fam and "gene_to_family" in file_mapping:
            filename = file_mapping["gene_to_family"]
            if os.path.exists(filename):
                result['gene_to_fam'] = tox_deserialize_int_nd(filename)
                print(f"Gene to family mapping extracted from {filename}")

        if load_family_ids and "family_ids" in file_mapping:
            filename = file_mapping["family_ids"]
            if os.path.exists(filename):
                result['family_ids'] = tox_deserialize_char_nd(filename)
                print(f"Family IDs extracted from {filename}")

        if load_family_centroids and "family_centroids" in file_mapping:
            filename = file_mapping["family_centroids"]
            if os.path.exists(filename):
                result['family_centroids'] = tox_deserialize_real_nd(filename)
                print(f"Family centroids extracted from {filename}")

        if load_shift_vectors and "shift_vectors" in file_mapping:
            filename = file_mapping["shift_vectors"]
            if os.path.exists(filename):
                result['shift_vectors'] = tox_deserialize_real_nd(filename)
                print(f"Shift vectors extracted from {filename}")

    finally:
        # Cleanup extracted files
        files_to_remove = [manifest_path]
        for filename in file_mapping.values():
            if os.path.exists(filename):
                files_to_remove.append(filename)

        for file in files_to_remove:
            if os.path.exists(file):
                os.remove(file)
                print(f"Cleaned up: {file}")

    return result
