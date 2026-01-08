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
    tox_deserialize_real_nd
)

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL(dll_path)

# Example for read_gene_ids_from_tsv_file_C
lib.read_gene_ids_from_tsv_file_C.argtypes = [
    ctypes.POINTER(ctypes.c_char),  # filename_raw
    ctypes.c_int,                  # fn_len
    ctypes.POINTER(ctypes.c_char),  # gene_ids_raw
    ctypes.c_int,                  # gene_ids_len
    ctypes.c_int,                  # n_genes
    ctypes.c_int,                  # n_header_rows
    ctypes.c_int,                  # gene_col
    ctypes.POINTER(ctypes.c_int)   # ierr
]
lib.read_gene_ids_from_tsv_file_C.restype = None

# read_expression_vectors_tsv_C
lib.read_expression_vectors_tsv_C.argtypes = [
    ctypes.POINTER(ctypes.c_char),  # file_list_raw
    ctypes.c_int,                  # file_list_len
    ctypes.c_int,                  # n_files
    ctypes.POINTER(ctypes.c_char),  # gene_ids_raw
    ctypes.c_int,                  # gene_ids_len
    ctypes.c_int,                  # n_genes
    ctypes.POINTER(ctypes.c_double), # expression_vectors_flat
    ctypes.c_int,                  # n_samples
    ctypes.c_int,                  # n_header_rows
    ctypes.c_int,                  # gene_col
    ctypes.POINTER(ctypes.c_int),  # value_cols
    ctypes.c_int,                  # n_value_cols
    ctypes.POINTER(ctypes.c_int),  # ierr
    ctypes.POINTER(ctypes.c_char)  # delimiter_raw
]
lib.read_expression_vectors_tsv_C.restype = None

# read_orthofinder_file_C
lib.read_orthofinder_file_C.argtypes = [
    ctypes.POINTER(ctypes.c_char),  # filename_raw
    ctypes.c_int,                  # fn_len
    ctypes.POINTER(ctypes.c_char),  # gene_ids_raw
    ctypes.c_int,                  # gene_ids_len
    ctypes.c_int,                  # n_genes
    ctypes.POINTER(ctypes.c_char),  # family_ids_raw
    ctypes.c_int,                  # family_ids_len
    ctypes.c_int,                  # n_families
    ctypes.POINTER(ctypes.c_int),  # gene_to_fam
    ctypes.POINTER(ctypes.c_int)   # ierr
]
lib.read_orthofinder_file_C.restype = None

# --- Validation C function signatures ---

lib.validate_gene_to_family_mapping_C.argtypes = [
    ctypes.POINTER(ctypes.c_int), # gene_to_fam
    ctypes.c_int,                 # n_genes
    ctypes.c_int,                 # n_families
    ctypes.POINTER(ctypes.c_int)  # ierr
]
lib.validate_gene_to_family_mapping_C.restype = None

lib.validate_expression_data_C.argtypes = [
    ctypes.POINTER(ctypes.c_double), # expression_vectors
    ctypes.c_int,                    # n_genes
    ctypes.c_int,                    # n_samples
    ctypes.c_int,                    # check_non_negative
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_expression_data_C.restype = None

lib.validate_family_centroids_C.argtypes = [
    ctypes.POINTER(ctypes.c_double), # family_centroids
    ctypes.c_int,                    # n_families
    ctypes.c_int,                    # n_samples
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_family_centroids_C.restype = None

lib.validate_shift_vectors_C.argtypes = [
    ctypes.POINTER(ctypes.c_double), # shift_vectors
    ctypes.POINTER(ctypes.c_double), # expression_vectors
    ctypes.POINTER(ctypes.c_double), # family_centroids
    ctypes.POINTER(ctypes.c_int),    # gene_to_fam
    ctypes.c_int,                    # n_genes
    ctypes.c_int,                    # n_samples
    ctypes.c_int,                    # n_families
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_shift_vectors_C.restype = None

lib.validate_string_array_uniqueness_C.argtypes = [
    ctypes.POINTER(ctypes.c_char), # string_array_raw
    ctypes.c_int,                 # string_len
    ctypes.c_int,                 # n_strings
    ctypes.POINTER(ctypes.c_int)  # ierr
]
lib.validate_string_array_uniqueness_C.restype = None

lib.validate_data_structure_C.argtypes = [
    ctypes.c_int,                 # n_genes
    ctypes.c_int,                 # n_families
    ctypes.c_int,                 # n_samples
    ctypes.POINTER(ctypes.c_char), # gene_ids_raw
    ctypes.c_int,                 # gene_ids_len
    ctypes.POINTER(ctypes.c_char), # gene_family_ids_raw
    ctypes.c_int,                 # fam_len
    ctypes.POINTER(ctypes.c_int), # gene_to_fam
    ctypes.POINTER(ctypes.c_double), # expression_vectors
    ctypes.POINTER(ctypes.c_double), # family_centroids
    ctypes.POINTER(ctypes.c_double), # shift_vectors
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_data_structure_C.restype = None

lib.validate_all_data_C.argtypes = [
    ctypes.c_int,                 # n_genes
    ctypes.c_int,                 # n_families
    ctypes.c_int,                 # n_samples
    ctypes.POINTER(ctypes.c_char), # gene_ids_raw
    ctypes.c_int,                 # gene_len
    ctypes.POINTER(ctypes.c_char), # gene_family_ids_raw
    ctypes.c_int,                 # fam_len
    ctypes.POINTER(ctypes.c_int), # gene_to_fam
    ctypes.POINTER(ctypes.c_double), # expression_vectors
    ctypes.POINTER(ctypes.c_double), # family_centroids
    ctypes.POINTER(ctypes.c_double), # shift_vectors
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_all_data_C.restype = None

# converts a given string to a c_char array of given length
def string_to_c_char_array(s, length):
    """
    Convert a Python string to a ctypes c_char array of a given length, with null termination.

    Args:
        s (str): Input string.
        length (int): Length of the output c_char array.

    Returns:
        ctypes.Array: c_char array of specified length, null-terminated if space allows.
    """
    if s is None:
        s = ""
    # Create array of c_char with specified length
    arr = (ctypes.c_char * length)()
    # Encode string and copy to array
    encoded = s.encode('ascii')
    for i in range(min(length, len(encoded))):
        arr[i] = encoded[i]
    # Ensure null termination if there's space
    if len(encoded) < length:
        arr[len(encoded)] = b'\x00'
    return arr

# converts a c_char array back to a string
def c_char_array_to_string(c_array):
    """
    Convert a ctypes c_char array back to a Python string, stopping at the first null byte.

    Args:
        c_array (ctypes.Array): Input c_char array.

    Returns:
        str: Decoded string up to the first null byte.
    """
    # Find null terminator or use full length
    bytes_list = []
    for i in range(len(c_array)):
        byte = c_array[i]
        if byte == b'\x00':
            break
        bytes_list.append(byte)
    return b''.join(bytes_list).decode('ascii').strip()

def strings_to_c_char_matrix(strings, max_length):
    """
    Convert a list of strings to a flat c_char array (Fortran-compatible, NumPy-wrapped).

    Args:
        strings (list of str): List of strings to convert.
        max_length (int): Maximum length for each string (padded or truncated).

    Returns:
        np.ndarray: 2D NumPy array of shape (n_strings, max_length) in Fortran order.
    """
    import numpy as np
    n_strings = len(strings)
    total_size = n_strings * max_length
    
    # create flat c-types array
    matrix_type = ctypes.c_char * total_size
    matrix = matrix_type()

    # Initialize with all null bytes
    for i in range(total_size):
        matrix[i] = b'\x00'

    for i, s in enumerate(strings):
        encoded = s.encode('ascii')
        for j in range(min(max_length, len(encoded))):
            index = j + i * max_length
            matrix[index] = encoded[j:j+1]
        if len(encoded) < max_length:
            matrix[len(encoded) + i * max_length] = b'\x00'

    arr = np.ctypeslib.as_array(matrix)
    arr = arr.reshape((n_strings, max_length), order='F')

    return arr

def c_char_matrix_to_strings(matrix, max_length, n_strings):
    """
    Convert a 2D c_char matrix or NumPy array back to a list of strings (Fortran order).

    Args:
        matrix (np.ndarray or ctypes.Array): 2D c_char matrix or flat array.
        max_length (int): Maximum length of each string.
        n_strings (int): Number of strings.

    Returns:
        list of str: Decoded strings from the matrix.
    """
    import numpy as np

    if isinstance(matrix, np.ndarray):
        flat = matrix.ravel(order='F')
    else:
        flat = matrix

    strings = []
    for i in range(n_strings):
        chars = []
        for j in range(max_length):
            index = j + i * max_length
            if index >= len(flat):
                break

            char = flat[index]

            if isinstance(char, np.ndarray):
                char = char.item()
            
            if isinstance(char, (np.integer, int)):
                char = bytes([char])
            elif isinstance(char, (bytes, bytearray)):
                pass
            elif isinstance(char, str):
                char = char.encode('ascii')
            else:
                char = bytes([int(char)])

            if char == b'\x00':
                break

            chars.append(char)

        s = b''.join(chars).decode('ascii')
        strings.append(s)

    return strings


def _ensure_string_array(arr):
    """
    Ensure the input is a NumPy array of strings.

    Args:
        arr: Input array-like object.

    Returns:
        np.ndarray: Array of strings.
    """
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr, dtype=object)
    return arr

def _ensure_float_array(arr):
    """
    Ensure the input is a NumPy array of floats (np.float64).

    Args:
        arr: Input array-like object.

    Returns:
        np.ndarray: Array of floats.
    """
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr, dtype=np.float64)
    return arr

def _ensure_int_array(arr):
    """
    Ensure the input is a NumPy array of ints (np.int32).

    Args:
        arr: Input array-like object.

    Returns:
        np.ndarray: Array of ints.
    """
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr, dtype=np.int32)
    return arr

#' Function for read_gene_ids_from_tsv_file_C
def read_gene_ids_from_tsv_file(filename, n_genes, gene_ids_len, n_header_rows, gene_col):
    """
    Read gene IDs from a TSV file using the C interface.

    Args:
        filename (str): Name of the TSV file.
        n_genes (int): Number of genes to read.
        gene_ids_len (int): Maximum length of each gene ID.
        n_header_rows (int): Number of header rows to skip.
        gene_col (int): Column index (1-based) for gene IDs.

    Returns:
        np.ndarray: Array of gene IDs as strings.
    """
    # Ensure filename is a string (single file)
    if isinstance(filename, list):
        if len(filename) > 0:
            filename = filename[0]
        else:
            raise ValueError("filename cannot be an empty list for read_gene_ids_from_tsv_file")
    
    # Convert filename to c_char array
    fn_array = string_to_c_char_array(filename, len(filename))
    
    # Create output array for gene IDs - proper 2D array
    matrix_size = gene_ids_len * n_genes
    gene_ids_array = (ctypes.c_char * matrix_size)()
    
    ierr = ctypes.c_int()
    
    print(f"Debug: Creating gene_ids_array with size {matrix_size} = {gene_ids_len} * {n_genes}")
    
    # Call C function
    lib.read_gene_ids_from_tsv_file_C(
        ctypes.cast(ctypes.byref(fn_array), ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(len(filename)),
        ctypes.cast(ctypes.byref(gene_ids_array), ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(gene_ids_len),
        ctypes.c_int(n_genes),
        ctypes.c_int(n_header_rows),
        ctypes.c_int(gene_col),
        ctypes.byref(ierr)
    )
    
    check_err_code(ierr.value)
    
    # Convert result back to strings
    gene_ids = c_char_matrix_to_strings(gene_ids_array, gene_ids_len, n_genes)
    
    return np.array(gene_ids, dtype='U')

# Function for read_expression_vectors_tsv_C
def read_expression_vectors_tsv(file_list, gene_ids, n_samples, n_header_rows, 
                           gene_col, value_cols, delimiter='\t'):
    """
    Read expression vectors from given tabular (CSV/TSV) files using the C interface.

    Args:
        file_list (list of str): List of filenames to read.
        gene_ids (list of str): List of gene IDs to extract.
        n_samples (int): Number of samples (files).
        n_header_rows (int): Number of header rows to skip in each file.
        gene_col (int): Column index (1-based) for gene IDs.
        value_cols (list of int): List of column indices (1-based) for expression values.
        delimiter (str): Delimiter used in the files (default: tab).

    Returns:
        np.ndarray: 2D array of expression vectors (samples x genes).
    """
    # Ensure file_list is a list (multiple files)
    if not isinstance(file_list, list):
        file_list = [file_list]  # Convert single file to list
    
    # Ensure inputs are numpy arrays
    gene_ids = _ensure_string_array(gene_ids)
    
    # Convert inputs to c_char matrices
    max_file_len = max(len(f) for f in file_list)
    file_list_matrix = strings_to_c_char_matrix(file_list, max_file_len)
    
    max_gene_len = max(len(g) for g in gene_ids)
    gene_ids_matrix = strings_to_c_char_matrix(gene_ids, max_gene_len)
    
    # Prepare output arrays
    expression_vectors = np.zeros((n_samples, len(gene_ids)), dtype=np.float64, order='F')
    ierr = ctypes.c_int()
    delimiter_array = string_to_c_char_array(delimiter, len(delimiter))
    
    # Convert value_cols to ctypes array
    value_cols_ct = (ctypes.c_int * len(value_cols))(*value_cols)
    
    # Call C function
    lib.read_expression_vectors_tsv_C(
        file_list_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(max_file_len),
        ctypes.c_int(len(file_list)),
        gene_ids_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(max_gene_len),
        ctypes.c_int(len(gene_ids)),
        expression_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n_samples),
        ctypes.c_int(n_header_rows),
        ctypes.c_int(gene_col),
        value_cols_ct,
        ctypes.c_int(len(value_cols)),
        ctypes.byref(ierr),
        ctypes.cast(ctypes.byref(delimiter_array), ctypes.POINTER(ctypes.c_char))
    )
    
    check_err_code(ierr.value)
        
    return expression_vectors

# Function for read_orthofinder_file_C
def read_orthofinder_file(filename, gene_ids, family_ids_len, n_families):
    """
    Read an OrthoFinder family file and map genes to families using the C interface.

    Args:
        filename (str): Name of the OrthoFinder TSV file.
        gene_ids (list of str): List of gene IDs to map.
        family_ids_len (int): Maximum length of each family ID.
        n_families (int): Number of families to read.

    Returns:
        dict: {'family_ids': np.ndarray, 'gene_to_fam': np.ndarray}
    """
    # Ensure filename is a string (single file)
    if isinstance(filename, list):
        if len(filename) > 0:
            filename = filename[0]  # Take first element if it's a list
        else:
            raise ValueError("filename cannot be an empty list for read_orthofinder_file")
    
    # Ensure inputs are numpy arrays
    gene_ids = _ensure_string_array(gene_ids)
    
    # Convert inputs to c_char arrays
    fn_array = string_to_c_char_array(filename, len(filename))
    
    max_gene_len = max(len(g) for g in gene_ids)
    gene_ids_matrix = strings_to_c_char_matrix(gene_ids, max_gene_len)
    
    # Prepare output arrays
    family_ids_matrix = strings_to_c_char_matrix([""] * n_families, family_ids_len)
    gene_to_fam = np.zeros(len(gene_ids), dtype=np.int32, order='F')
    ierr = ctypes.c_int()
    
    # Call C function
    lib.read_orthofinder_file_C(
        ctypes.cast(ctypes.byref(fn_array), ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(len(filename)),
        gene_ids_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(max_gene_len),
        ctypes.c_int(len(gene_ids)),
        family_ids_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(family_ids_len),
        ctypes.c_int(n_families),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.byref(ierr)
    )

    
    check_err_code(ierr.value)
    
    # Convert result back to strings as numpy array
    family_ids = c_char_matrix_to_strings(family_ids_matrix, family_ids_len, n_families)
    
    return {
        'family_ids': np.array(family_ids, dtype='U'),
        'gene_to_fam': gene_to_fam
    }

def filter_unassigned_genes(gene_to_fam):
    """
    Filter out genes that are not assigned to any family (where gene_to_fam == 0).

    Args:
        gene_to_fam (list[int] or np.ndarray): Family assignment for each gene (0 means unassigned).

    Returns:
        dict: {'mask': list[int] of 1s (kept) and 0s (removed), 'n_genes_kept': int}
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

def validate_gene_to_family_mapping(gene_to_fam, n_families):
    """
    Validate gene to family mapping using the C interface.

    Args:
        gene_to_fam (array-like): Array mapping each gene to a family index (0 if unassigned).
        n_families (int): Total number of families.
    """
    gene_to_fam = _ensure_int_array(gene_to_fam)
    n_genes = gene_to_fam.size
    ierr = ctypes.c_int()
    lib.validate_gene_to_family_mapping_C(
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(n_genes),
        ctypes.c_int(n_families),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

def validate_expression_data(expression_vectors, check_non_negative=True):
    """
    Validate expression data using the C interface.

    Args:
        expression_vectors (np.ndarray): 2D array of expression data (samples x genes).
        check_non_negative (bool): Whether to check for non-negative values.
    """
    arr = _ensure_float_array(expression_vectors)
    arr = np.asfortranarray(arr, dtype=np.float64)
    n_samples, n_genes = arr.shape
    ierr = ctypes.c_int()
    lib.validate_expression_data_C(
        arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n_genes),
        ctypes.c_int(n_samples),
        ctypes.c_int(1 if check_non_negative else 0),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

def validate_family_centroids(family_centroids):
    """
    Validate family centroids for NaN/Inf using the C interface.

    Args:
        family_centroids (np.ndarray): 2D array of family centroids (samples x families).
    """
    arr = _ensure_float_array(family_centroids)
    arr = np.asfortranarray(arr, dtype=np.float64)
    n_samples, n_families = arr.shape
    ierr = ctypes.c_int()
    lib.validate_family_centroids_C(
        arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n_families),
        ctypes.c_int(n_samples),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

def validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_genes, n_samples, n_families):
    """
    Validate shift vectors, checking datatypes and structure using the C interface.

    Args:
        shift_vectors (np.ndarray): 2D array of shift vectors (2*samples x genes).
        expression_vectors (np.ndarray): 2D array of expression data (samples x genes).
        family_centroids (np.ndarray): 2D array of family centroids (samples x families).
        gene_to_fam (np.ndarray): Array mapping each gene to a family index (0 if unassigned).
        n_genes (int): Number of genes.
        n_samples (int): Number of samples.
        n_families (int): Number of families.
    """
    shift_vectors = _ensure_float_array(shift_vectors)
    expression_vectors = _ensure_float_array(expression_vectors)
    family_centroids = _ensure_float_array(family_centroids)
    gene_to_fam = _ensure_int_array(gene_to_fam)
    
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    ierr = ctypes.c_int()
    lib.validate_shift_vectors_C(
        shift_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        expression_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        family_centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(n_genes),
        ctypes.c_int(n_samples),
        ctypes.c_int(n_families),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

def validate_string_array_uniqueness(strings):
    """
    Validate uniqueness of strings using the C interface.

    Note: Uses a hashset internally, which may increase memory usage for large datasets.

    Args:
        strings (list of str): List of strings.
    """
    strings = _ensure_string_array(strings)
    string_len = max(len(g) for g in strings)
    n_strings = len(strings)
    gene_ids_raw = strings_to_c_char_matrix(strings, string_len)
    ierr = ctypes.c_int()
    lib.validate_string_array_uniqueness_C(
        gene_ids_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(string_len),
        ctypes.c_int(n_strings),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

def validate_data_structure(n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors):
    """
    Validate overall data structure consistency using the C interface.

    Args:
        n_genes (int): Number of genes.
        n_families (int): Number of families.
        n_samples (int): Number of samples.
        gene_ids (list of str): List of gene IDs.
        gene_family_ids (list of str): List of family IDs.
        gene_to_fam (np.ndarray): Array mapping each gene to a family index (0 if unassigned).
        expression_vectors (np.ndarray): 2D array of expression data (samples x genes).
        family_centroids (np.ndarray): 2D array of family centroids (samples x families).
        shift_vectors (np.ndarray): 2D array of shift vectors (2*samples x genes).
    """
    gene_ids = _ensure_string_array(gene_ids)
    gene_family_ids = _ensure_string_array(gene_family_ids)
    expression_vectors = _ensure_float_array(expression_vectors)
    family_centroids = _ensure_float_array(family_centroids)
    shift_vectors = _ensure_float_array(shift_vectors)
    gene_to_fam = _ensure_int_array(gene_to_fam)
    
    gene_ids_len = max(len(g) for g in gene_ids)
    fam_len = max(len(f) for f in gene_family_ids)
    gene_ids_raw = strings_to_c_char_matrix(gene_ids, gene_ids_len)
    gene_family_ids_raw = strings_to_c_char_matrix(gene_family_ids, fam_len)
    
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    ierr = ctypes.c_int()
    
    lib.validate_data_structure_C(
        ctypes.c_int(n_genes),
        ctypes.c_int(n_families),
        ctypes.c_int(n_samples),
        gene_ids_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(gene_ids_len),
        gene_family_ids_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(fam_len),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        expression_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        family_centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        shift_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

def validate_all_data(n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors):
    """
    Comprehensive validation of all data components using the C interface.
    This function performs all individual validations in one go.

    Args:
        n_genes (int): Number of genes.
        n_families (int): Number of families.
        n_samples (int): Number of samples.
        gene_ids (list of str): List of gene IDs.
        gene_family_ids (list of str): List of family IDs.
        gene_to_fam (np.ndarray): Array mapping each gene to a family index (0 if unassigned).
        expression_vectors (np.ndarray): 2D array of expression data (samples x genes).
        family_centroids (np.ndarray): 2D array of family centroids (samples x families).
        shift_vectors (np.ndarray): 2D array of shift vectors (genes x samples).
    """
    gene_ids = _ensure_string_array(gene_ids)
    gene_family_ids = _ensure_string_array(gene_family_ids)
    expression_vectors = _ensure_float_array(expression_vectors)
    family_centroids = _ensure_float_array(family_centroids)
    shift_vectors = _ensure_float_array(shift_vectors)
    gene_to_fam = _ensure_int_array(gene_to_fam)
    
    gene_ids_len = max(len(g) for g in gene_ids)
    fam_len = max(len(f) for f in gene_family_ids)
    gene_ids_raw = strings_to_c_char_matrix(gene_ids, gene_ids_len)
    gene_family_ids_raw = strings_to_c_char_matrix(gene_family_ids, fam_len)
    
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    ierr = ctypes.c_int()
    
    lib.validate_all_data_C(
        ctypes.c_int(n_genes),
        ctypes.c_int(n_families),
        ctypes.c_int(n_samples),
        gene_ids_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(gene_ids_len),
        gene_family_ids_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(fam_len),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        expression_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        family_centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        shift_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)
    
def _prepare_string(s) -> tuple:
    """
    Prepare a string for the C interface (returns bytes and length).

    Args:
        s (str): Input string.

    Returns:
        tuple: (bytes, int) for the C interface.
    """
    if s is None or s == "":
        # Return empty string with null terminator
        return b'\x00', 1
    # Include null terminator
    b = s.encode('utf-8') + b'\x00'
    return b, len(b)

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
    
    # Prepare arrays for C interface
    zip_b, zip_len = _prepare_string(zip_filename)
    
    # Convert keys and filenames to 2D numpy arrays with Fortran order
    max_key_len = max(len(key) for key in keys) + 1  # +1 for null terminator
    max_filename_len = max(len(fname) for fname in filenames) + 1  # +1 for null terminator
    count = len(keys)
    
    # Create numpy arrays with Fortran order (column-major)
    # Shape: (string_length, string_count) - matching Fortran expectation
    keys_array = np.zeros((max_key_len, count), dtype=np.byte, order='F')
    filenames_array = np.zeros((max_filename_len, count), dtype=np.byte, order='F')
    
    # Fill arrays in column-major order
    for i, (key, fname) in enumerate(zip(keys, filenames)):
        key_bytes = key.encode('utf-8')
        fname_bytes = fname.encode('utf-8')
        
        # Fill keys array - each column is one string
        for j in range(min(len(key_bytes), max_key_len - 1)):
            keys_array[j, i] = key_bytes[j]
        # Add null terminator
        if len(key_bytes) < max_key_len:
            keys_array[len(key_bytes), i] = 0  # null terminator
        else:
            keys_array[max_key_len - 1, i] = 0  # null terminator at last position
        
        # Fill filenames array - each column is one string  
        for j in range(min(len(fname_bytes), max_filename_len - 1)):
            filenames_array[j, i] = fname_bytes[j]
        # Add null terminator
        if len(fname_bytes) < max_filename_len:
            filenames_array[len(fname_bytes), i] = 0  # null terminator
        else:
            filenames_array[max_filename_len - 1, i] = 0  # null terminator at last position
    
    # Set up argument types to match Fortran subroutine
    lib.create_zip_archive_generic_c.argtypes = [
        ctypes.POINTER(ctypes.c_char), ctypes.c_int,                    # zip_filename, zip_len
        np.ctypeslib.ndpointer(dtype=np.byte, ndim=2, flags='F_CONTIGUOUS'),  # keys
        ctypes.c_int, ctypes.c_int,                                     # keys_len, keys_count
        np.ctypeslib.ndpointer(dtype=np.byte, ndim=2, flags='F_CONTIGUOUS'),  # filenames  
        ctypes.c_int, ctypes.c_int,                                     # filenames_len, filenames_count
        ctypes.POINTER(ctypes.c_int)                                    # ierr
    ]
    
    ierr = ctypes.c_int()
    
    # Call Fortran function - arrays are already in Fortran order
    lib.create_zip_archive_generic_c(
        zip_b, 
        ctypes.c_int(zip_len),
        keys_array,
        ctypes.c_int(max_key_len), 
        ctypes.c_int(count),
        filenames_array, 
        ctypes.c_int(max_filename_len), 
        ctypes.c_int(count),
        ctypes.byref(ierr)
    )
    
    check_err_code(ierr.value)
    print(f"Successfully created archive: {zip_filename}")

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

    lib.extract_zip_archive_generic_c.argtypes = [
        ctypes.POINTER(ctypes.c_char), 
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int)
    ]
    
    try:
        # Initialize error code
        ierr = ctypes.c_int()
        
        # Prepare filename for C function
        zip_b = (ctypes.c_char * len(zip_filename))(*zip_filename.encode('utf-8'))
        
        # Call the C extraction function
        lib.extract_zip_archive_generic_c(
            ctypes.cast(zip_b, ctypes.POINTER(ctypes.c_char)),
            ctypes.c_int(len(zip_filename)),
            ctypes.byref(ierr)
        )
        
        # Check for errors
        if ierr.value != 0:
            raise RuntimeError(f"Failed to extract zip archive. Error code: {ierr.value}")
        
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
    
    # Prepare the zip filename
    zip_b = string_to_c_char_array(zip_filename, len(zip_filename))
    zip_len = len(zip_filename)

    
    # Set up argument types for extract_zip_archive_c - using c_int
    lib.extract_zip_archive_generic_c.argtypes = [
        ctypes.POINTER(ctypes.c_char), ctypes.c_int,  # zip_filename, filename_len
        ctypes.POINTER(ctypes.c_int)                  # ierr
    ]
    lib.extract_zip_archive_generic_c.restype = None
    
    ierr = ctypes.c_int()
    
    # Call the C-bound Fortran subroutine to extract the zip archive
    lib.extract_zip_archive_generic_c(
        ctypes.cast(ctypes.byref(zip_b), ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(zip_len),
        ctypes.byref(ierr)
    )

    check_err_code(ierr.value)
    
    print(f"Successfully extracted archive: {zip_filename}")
    
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