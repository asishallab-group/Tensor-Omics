"""
TensorOmics Functions Module
Python wrapper functions for Fortran routines via C interface
"""

import numpy as np
import ctypes
import os
from error_handling import check_err_code

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL(dll_path)

# helper to convert a filename to ASCII chars to transfer it as integer to fortran
def _filename_to_ascii_array(filename):
    ascii_arr = np.array([ord(c) for c in filename], dtype=np.int32)
    return ascii_arr, np.int32(len(ascii_arr))

# Function to convert array of strings to integer ASCII array
def _string_array_to_ascii_matrix(strings: np.ndarray) -> tuple[np.ndarray, int]:
    """
    Transforms a char array to integer
    """
    if not isinstance(strings, np.ndarray) or strings.dtype.kind != 'U':
        raise ValueError("Input must be a numpy array of strings (dtype='U')")
    
    flat = strings.flatten(order='F')  # important: Fortran-order
    clen = max(len(s) for s in flat)
    total = flat.size

    mat = np.zeros((clen, total), dtype=np.int32, order='F')
    for i, s in enumerate(flat):
        codes = [ord(c) for c in s]
        mat[:len(codes), i] = codes

    return mat, clen

# Helper function to read dimensions of integer/real array
def tox_get_array_metadata(filename, max_dims=5, with_clen=False):
    """
    Reads dimensions (and optionally character length) of an array file.
    with_clen=True -> returns (dims, clen)
    with_clen=False -> returns dims only
    """
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    dims_out = np.zeros(max_dims, dtype=np.int32)
    ndims = ctypes.c_int()
    ierr = ctypes.c_int()
    clen = ctypes.c_int()  # always pass
    dims_out_capacity = ctypes.c_int(max_dims)

    # shared function
    lib.get_array_metadata_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"), # filename_ascii
        ctypes.c_int,                                                         # fn_len
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"), # dims_out
        ctypes.POINTER(ctypes.c_int),                                         # dims_out_capacity
        ctypes.POINTER(ctypes.c_int),                                         # ndims
        ctypes.POINTER(ctypes.c_int),                                         # ierr
        ctypes.POINTER(ctypes.c_int)                                          # clen
    ]
    lib.get_array_metadata_C.restype = None

    # call
    lib.get_array_metadata_C(
        filename_ascii,
        fn_len,
        dims_out,
        ctypes.byref(dims_out_capacity),
        ctypes.byref(ndims),
        ctypes.byref(ierr),
        ctypes.byref(clen)
    )

    check_err_code(ierr.value)

    if with_clen:
        return dims_out[:ndims.value], clen.value
    else:
        return dims_out[:ndims.value]

# serilization of an n-dimensional integer array
def tox_serialize_int_nd(arr: np.ndarray, filename: str):
    """
    Serializes an n-dimensional integer32 array to a binary file
    """
    if not isinstance(arr, np.ndarray) or arr.dtype != np.int32:
        raise ValueError("arr must be a numpy array of int32")

    # Make sure layout is fortran compatible
    arr_f = np.asfortranarray(arr)

    # dimensions
    dims = np.array(arr.shape, dtype=np.int32)
    ndim = arr.ndim
    ierr = ctypes.c_int()

    # flat array to pass to fortran
    flat = arr_f.ravel(order='F')

    # prepare filename
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    lib.serialize_int_nd_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # arr
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # dims
        ctypes.c_int,  # ndim
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,  # fn_len
        ctypes.POINTER(ctypes.c_int) 
    ]
    lib.serialize_int_nd_C.restype = None

    # call function
    lib.serialize_int_nd_C(
        flat,
        dims,
        ndim,
        filename_ascii,
        fn_len,
        ctypes.byref(ierr)
    )

    check_err_code(ierr.value)

# Deserialize an n dimensional integer array
def tox_deserialize_int_nd(filename):
    """
    Deserializes an n-dimensional int32-Array.
    """
    # read size of the array
    dims = tox_get_array_metadata(filename)
    print(f"Deserializing array with dimensions: {dims}")
    # create array with the proper size
    total_size = np.prod(dims)
    arr = np.zeros(total_size, dtype=np.int32, order='F')  # gets a 1D array
    ascii_arr, fn_len = _filename_to_ascii_array(filename)
    ierr = ctypes.c_int()

    lib.deserialize_int_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="F_CONTIGUOUS"),  # arr
        ctypes.c_int,                                                          # total size
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,                                                          # fn_len
        ctypes.POINTER(ctypes.c_int)                                           # ierr
    ]
    lib.deserialize_int_C.restype = None

    lib.deserialize_int_C(arr, total_size, ascii_arr, fn_len, ctypes.byref(ierr))
    check_err_code(ierr.value)
    return arr.reshape(dims, order='F')  # Reshape to original dimensions

def tox_serialize_real_nd(arr: np.ndarray, filename: str):
    """
    Serializes an n-dimensional real64 array to a binary file
    """
    if not isinstance(arr, np.ndarray) or arr.dtype != np.float64:
        raise ValueError("arr must be a numpy array of float64")

    # make sure layout is fortran compatible
    arr_f = np.asfortranarray(arr)

    # dimensions
    dims = np.array(arr.shape, dtype=np.int32)
    ndim = arr.ndim
    ierr = ctypes.c_int()

    # flat array with fortran order
    flat = arr_f.ravel(order='F')

    # ASCII-Filename preparation
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    # declare args
    lib.serialize_real_nd_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"), # arr
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # dims
        ctypes.c_int,  # ndim
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,  # fn_len
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    lib.serialize_real_nd_C.restype = None

    # call function
    lib.serialize_real_nd_C(
        flat,
        dims,
        ndim,
        filename_ascii,
        fn_len,
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

def tox_deserialize_real_nd(filename):
    """
    Deserializes an n-dimensional array of type real64
    """
    #read dimensions
    dims = tox_get_array_metadata(filename)
    print(f"Deserializing array with dimensions: {dims}")
    # create array with correct size
    total_size = np.prod(dims)
    arr = np.zeros(total_size, dtype=np.float64, order='F')  # accept flat array
    ascii_arr, fn_len = _filename_to_ascii_array(filename)
    ierr = ctypes.c_int()

    lib.deserialize_real_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # arr
        ctypes.c_int,                                                          # total size
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,                                                           # fn_len
        ctypes.POINTER(ctypes.c_int)                                           # ierr
    ]
    lib.deserialize_real_C.restype = None

    lib.deserialize_real_C(arr, total_size, ascii_arr, fn_len, ctypes.byref(ierr))
    check_err_code(ierr.value)
    return arr.reshape(dims, order='F')  # Reshape

def tox_serialize_char_nd(arr: np.ndarray, filename: str):
    """
    Serializes an n-dimensional character array to a binary file
    """
    if not isinstance(arr, np.ndarray) or arr.dtype.kind != 'U':
        raise ValueError("arr must be a numpy array of strings (dtype='U')")

    dims = np.array(arr.shape, dtype=np.int32)
    ndim = arr.ndim
    ierr = ctypes.c_int()

    ascii_mat, clen = _string_array_to_ascii_matrix(arr)
    ascii_ptr = np.ascontiguousarray(ascii_mat.ravel(order='F'))

    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    lib.serialize_char_flat_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),  # ascii_ptr
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),  # dims
        ctypes.c_int,                                                          # ndim
        ctypes.c_int,                                                          # clen
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),  # filename_ascii
        ctypes.c_int,                                                           # fn_len
        ctypes.POINTER(ctypes.c_int)                                           # ierr
    ]
    lib.serialize_char_flat_C.restype = None

    lib.serialize_char_flat_C(
        ascii_ptr,
        dims,
        ndim,
        clen,
        filename_ascii,
        fn_len,
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)


def tox_deserialize_char_nd(filename: str, ndim_max=5):
    """
    Deserializes an n-dimensional character array from a file
    """
    # convert filename to ASCII array
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    # Get metadata (dimensions + string length)
    dims, clen = tox_get_array_metadata(filename, max_dims=ndim_max, with_clen=True)
    ierr = ctypes.c_int()

    print(f"Deserializing char array with dimensions: {dims}, clen: {clen}")
    total = int(np.prod(dims))  # Anzahl Elemente

    # allocate flat ASCII array (1D buffer from Fortran)
    ascii_arr = np.zeros(total * clen, dtype=np.int32)

    # declare arguments for the flat Fortran routine
    lib.deserialize_char_flat_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # ascii_arr
        ctypes.c_int,           # clen
        ctypes.c_int,           # total
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,           # fn_len
        ctypes.POINTER(ctypes.c_int)                                           # ierr
    ]
    lib.deserialize_char_flat_C.restype = None

    lib.deserialize_char_flat_C(
        ascii_arr,
        clen,
        total,
        filename_ascii,
        fn_len,
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

    # empty array
    if total == 0:
        return np.empty(tuple(dims), dtype=f'U{clen}')

    # 1) per elemnt one block of "clen" codes
    codes_2d = ascii_arr.reshape((total, clen), order='C')

    # 2) ASCII -> Strings per Element
    strings_1d = np.array(
        [''.join(chr(c) for c in row if c > 0) for row in codes_2d],
        dtype=f'U{clen}'
    )

    # 3) reshape to target
    result = strings_1d.reshape(tuple(dims), order='F')
    return result


# Configure BST argument types
lib.build_bst_index_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # values
    ctypes.c_int32,                                                  # num_values
    np.ctypeslib.ndpointer(dtype=np.int32),                          # sorted_indices (out)
    np.ctypeslib.ndpointer(dtype=np.int32),                          # stack_left
    np.ctypeslib.ndpointer(dtype=np.int32),                          # stack_right
    ctypes.POINTER(ctypes.c_int)                                     # ierr   
]

lib.bst_range_query_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),  # values
    np.ctypeslib.ndpointer(dtype=np.int32),    # sorted_indices
    ctypes.c_int32,                            # num_values
    ctypes.c_double,                           # low
    ctypes.c_double,                           # high
    np.ctypeslib.ndpointer(dtype=np.int32),    # out_indices (out)
    ctypes.POINTER(ctypes.c_int32),            # number_matches (out)
    ctypes.POINTER(ctypes.c_int)
]

# Configure KD-Tree argument types
lib.build_kd_index_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),  # X_flat (col-major)
    ctypes.c_int32,                                                  # d
    ctypes.c_int32,                                                  # n
    np.ctypeslib.ndpointer(dtype=np.int32),                          # kd_ix (out)
    np.ctypeslib.ndpointer(dtype=np.int32),                          # dim_order
    np.ctypeslib.ndpointer(dtype=np.int32),                          # work
    np.ctypeslib.ndpointer(dtype=np.float64),                        # subarray
    np.ctypeslib.ndpointer(dtype=np.int32),                          # perm
    np.ctypeslib.ndpointer(dtype=np.int32),                          # stack_left
    np.ctypeslib.ndpointer(dtype=np.int32),                          # stack_right
    ctypes.POINTER(ctypes.c_int)
]

# --- BST Functions ---
def build_bst_index(values):
    """
    Build a BST index for the given values.
    
    Parameters:
    values (np.array): 1D array of values to index
    
    Returns:
    np.array: BST indices (0-based for Python)
    """
    n = len(values)
    indices = np.empty(n, dtype=np.int32)
    stack_left = np.empty(n, dtype=np.int32)
    stack_right = np.empty(n, dtype=np.int32)
    ierr = ctypes.c_int()
    
    # Build BST index
    lib.build_bst_index_C(values, n, indices, stack_left, stack_right, ctypes.byref(ierr))
    check_err_code(ierr.value)
    
    # Convert from Fortran 1-based to Python 0-based indexing
    return indices - 1

def bst_range_query(values, indices, lower_bound, upper_bound):
    """
    Perform a range query on BST-indexed values.
    
    Parameters:
    values (np.array): Original values array
    indices (np.array): BST indices from build_bst_index
    lower_bound (float): Lower bound of range (inclusive)
    upper_bound (float): Upper bound of range (inclusive)
    
    Returns:
    tuple: (matching_indices, count) where matching_indices are 0-based Python indices
    """
    n = len(values)
    output_indices = np.empty(n, dtype=np.int32)
    match_count = ctypes.c_int32(0)
    ierr = ctypes.c_int()
    
    # Convert indices back to 1-based for Fortran
    indices_1based = indices + 1
    
    # Perform range query
    lib.bst_range_query_C(values, indices_1based, n, lower_bound, upper_bound, 
                         output_indices, ctypes.byref(match_count), ctypes.byref(ierr))
    check_err_code(ierr.value)
    
    # Convert from Fortran 1-based to Python 0-based indexing
    matching_indices = output_indices[:match_count.value] - 1
    
    return matching_indices, match_count.value

# --- KD-Tree Functions ---
def build_kd_index(points, dimension_order=None):
    """
    Build a KD-Tree index for the given points.
    
    Parameters:
    points (np.array): 2D array of points (d x n, Fortran order)
    dimension_order (np.array): Order of dimensions for splitting (1-based)
    
    Returns:
    np.array: KD-Tree indices (0-based for Python)
    """
    d, n = points.shape
    
    # Use default dimension order if not provided
    if dimension_order is None:
        dimension_order = np.arange(1, d + 1, dtype=np.int32)  # 1-based dimensions
    
    # Ensure Fortran order (column-major)
    if not points.flags.f_contiguous:
        points = np.asfortranarray(points)

    
    # Initialize arrays
    kd_indices = np.empty(n, dtype=np.int32)
    workspace = np.empty(n, dtype=np.int32)
    value_buffer = np.empty(n, dtype=np.float64)
    permutation = np.empty(n, dtype=np.int32)
    stack_left = np.empty(n, dtype=np.int32)
    stack_right = np.empty(n, dtype=np.int32)
    ierr = ctypes.c_int()
    
    # Build KD-Tree index using the flat array
    lib.build_kd_index_C(points, d, n, kd_indices, dimension_order, workspace, 
                        value_buffer, permutation, stack_left, stack_right, ctypes.byref(ierr))
    check_err_code(ierr.value)
    
    # Convert from Fortran 1-based to Python 0-based indexing
    return kd_indices - 1

def build_spherical_kd(vectors, dimension_order=None):
    """
    Build a spherical KD-Tree index for the given unit vectors.
    
    Parameters:
    vectors (np.array): 2D array of unit vectors (d x n, Fortran order)
    dimension_order (np.array): Order of dimensions for splitting (1-based)
    
    Returns:
    np.array: Spherical KD-Tree indices (0-based for Python)
    """
    # For spherical KD-Tree, we use the same implementation as regular KD-Tree
    # but with a different name for clarity
    return build_kd_index(vectors, dimension_order)

def _readonly(*arrays: np.ndarray) -> None:
    """Mark all given NumPy arrays as read-only."""
    for a in arrays:
        if isinstance(a, np.ndarray):
            a.flags.writeable = False
            # NOTE: Returned NumPy arrays are read-only for safety.
            # If you need to modify them (e.g., for plotting), use `.copy()`.


def tox_vector_RAP_projection(vecs, vecs_selection_mask, axes_selection_mask):
    """
    Project selected vectors onto RAP constructed from selected axes.
    Args:
        vecs: Expression vectors (n_axes x n_vecs)
        vecs_selection_mask: Boolean/integer array (length n_vecs)
        axes_selection_mask: Boolean/integer array (length n_axes)
    Returns:
        np.ndarray: Projected vectors (n_selected_axes x n_selected_vecs)
    Raises:
        RuntimeError: If Fortran routine returns error
    """
    vecs = np.asfortranarray(vecs, dtype=np.float64)
    vecs_selection_mask = np.ascontiguousarray(vecs_selection_mask, dtype=np.int32)
    axes_selection_mask = np.ascontiguousarray(axes_selection_mask, dtype=np.int32)
    n_axes, n_vecs = vecs.shape
    if len(vecs_selection_mask) != n_vecs:
        raise ValueError("vecs_selection_mask length must match number of columns in vecs")
    if len(axes_selection_mask) != n_axes:
        raise ValueError("axes_selection_mask length must match number of rows in vecs")
    n_selected_vecs = int(np.sum(vecs_selection_mask))
    n_selected_axes = int(np.sum(axes_selection_mask))
    projections = np.empty((n_selected_axes, n_selected_vecs), order="F", dtype=np.float64)
    ierr = ctypes.c_int(0)
    omics_vector_RAP_projection = lib.omics_vector_RAP_projection_c
    omics_vector_RAP_projection.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # vecs
        ctypes.c_int,                                                    # n_axes
        ctypes.c_int,                                                    # n_vecs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),   # vecs_selection_mask
        ctypes.c_int,                                                    # n_selected_vecs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),   # axes_selection_mask
        ctypes.c_int,                                                    # n_selected_axes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # projections
        ctypes.POINTER(ctypes.c_int)                                     # ierr
    ]
    omics_vector_RAP_projection.restype = None
    omics_vector_RAP_projection(
        vecs, n_axes, n_vecs,
        vecs_selection_mask, n_selected_vecs,
        axes_selection_mask, n_selected_axes,
        projections, ctypes.byref(ierr)
    )
    check_err_code(ierr.value)
    _readonly(projections)
    return projections

def tox_field_RAP_projection(vecs, vecs_selection_mask, axes_selection_mask):
    """
    Project selected vector fields onto RAP constructed from selected axes.
    Args:
        vecs: Vector fields (2*n_axes x n_vecs)
        vecs_selection_mask: Boolean/integer array (length n_vecs)
        axes_selection_mask: Boolean/integer array (length n_axes)
    Returns:
        np.ndarray: Projected vectors (n_selected_axes x n_selected_vecs)
    Raises:
        RuntimeError: If Fortran routine returns error
    """
    vecs = np.asfortranarray(vecs, dtype=np.float64)
    vecs_selection_mask = np.ascontiguousarray(vecs_selection_mask, dtype=np.int32)
    axes_selection_mask = np.ascontiguousarray(axes_selection_mask, dtype=np.int32)
    n_axes = len(axes_selection_mask)
    n_vecs = vecs.shape[1]
    if len(vecs_selection_mask) != n_vecs:
        raise ValueError("vecs_selection_mask length must match number of columns in vecs")
    if vecs.shape[0] != 2 * n_axes:
        raise ValueError("vecs must have 2*n_axes rows")
    n_selected_vecs = int(np.sum(vecs_selection_mask))
    n_selected_axes = int(np.sum(axes_selection_mask))
    projections = np.empty((n_selected_axes, n_selected_vecs), order="F", dtype=np.float64)
    ierr = ctypes.c_int(0)
    omics_field_RAP_projection = lib.omics_field_RAP_projection_c
    omics_field_RAP_projection.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # vecs
        ctypes.c_int,                                                    # n_axes
        ctypes.c_int,                                                    # n_vecs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),   # vecs_selection_mask
        ctypes.c_int,                                                    # n_selected_vecs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),   # axes_selection_mask
        ctypes.c_int,                                                    # n_selected_axes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # projections
        ctypes.POINTER(ctypes.c_int)                                     # ierr
    ]
    omics_field_RAP_projection.restype = None
    omics_field_RAP_projection(
        vecs, n_axes, n_vecs,
        vecs_selection_mask, n_selected_vecs,
        axes_selection_mask, n_selected_axes,
        projections, ctypes.byref(ierr)
    )
    check_err_code(ierr.value)
    _readonly(projections)
    return projections


def tox_clock_hand_angle_between_vectors(v1, v2, selected_axes_for_signed):
    """
    Calculate clock hand angle between two vectors
    
    Args:
        v1: First vector (numpy array)
        v2: Second vector (numpy array)
        selected_axes_for_signed: Integer array of axes to use for signed angle (length n_dims)
    
    Returns:
        float: Signed angle between vectors in degrees
    Raises:
        RuntimeError: If Fortran routine returns error
    """
    # Input validation and conversion
    v1 = np.ascontiguousarray(v1, dtype=np.float64)  # First vector
    v2 = np.ascontiguousarray(v2, dtype=np.float64)  # Second vector
    selected_axes_for_signed = np.ascontiguousarray(selected_axes_for_signed, dtype=np.int32)  # Axes for signed angle
    n_dims = len(v1)
    if len(v2) != n_dims:
        raise ValueError("v1 and v2 must have same length")
    # Para 2D y 3D, Fortran ignora selected_axes_for_signed, pero requiere longitud 3
    if n_dims <= 3:
        selected_axes_for_signed = np.array([1, 2, 1], dtype=np.int32)
    else:
        selected_axes_for_signed = np.ascontiguousarray(selected_axes_for_signed, dtype=np.int32)
        if len(selected_axes_for_signed) != 3:
            raise ValueError("selected_axes_for_signed must have length 3 for n_dims > 3")
        if np.any(selected_axes_for_signed < 1) or np.any(selected_axes_for_signed > n_dims):
            raise ValueError("selected_axes_for_signed indices must be in [1, n_dims] for n_dims > 3")
    # Prepare output and error code
    signed_angle = ctypes.c_double(0.0)
    ierr = ctypes.c_int(0)
    # Setup C wrapper
    clock_hand_angle = lib.clock_hand_angle_between_vectors_c
    clock_hand_angle.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # v1
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # v2
        ctypes.c_int,  # n_dims
        ctypes.POINTER(ctypes.c_double),  # signed_angle
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # selected_axes_for_signed
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    clock_hand_angle.restype = None
    # Call Fortran routine
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(signed_angle), selected_axes_for_signed, ctypes.byref(ierr))
    # Check for errors
    check_err_code(ierr.value)

    _readonly(signed_angle)
    return signed_angle.value

def tox_clock_hand_angles_for_shift_vectors(origins, targets, vecs_selection_mask, selected_axes_for_signed):
    """
    Calculate clock hand angles for shift vectors
    
    Args:
        origins: Origin vectors (n_dims x n_vecs)
        targets: Target vectors (n_dims x n_vecs)
        vecs_selection_mask: Boolean array indicating which vectors to process
        selected_axes_for_signed: Integer array of axes to use for signed angle (length n_dims)
    
    Returns:
        numpy.ndarray: Signed angles for selected vectors in degrees
    Raises:
        RuntimeError: If Fortran routine returns error
    """
    # Input validation and conversion
    origins = np.asfortranarray(origins, dtype=np.float64)  # Origin vectors
    targets = np.asfortranarray(targets, dtype=np.float64)  # Target vectors
    vecs_selection_mask = np.ascontiguousarray(vecs_selection_mask, dtype=np.int32)  # Selection mask
    selected_axes_for_signed = np.ascontiguousarray(selected_axes_for_signed, dtype=np.int32)  # Axes for signed angle
    n_dims, n_vecs = origins.shape
    if targets.shape != (n_dims, n_vecs):
        raise ValueError("origins and targets must have same shape")
    if len(vecs_selection_mask) != n_vecs:
        raise ValueError("vecs_selection_mask must match number of vectors")
    if n_dims <= 3:
        selected_axes_for_signed = np.array([1, 2, 1], dtype=np.int32)
    else:
        selected_axes_for_signed = np.ascontiguousarray(selected_axes_for_signed, dtype=np.int32)
        if len(selected_axes_for_signed) != 3:
            raise ValueError("selected_axes_for_signed must have length 3 for n_dims > 3")
        if np.any(selected_axes_for_signed < 1) or np.any(selected_axes_for_signed > n_dims):
            raise ValueError("selected_axes_for_signed indices must be in [1, n_dims] for n_dims > 3")
    n_selected_vecs = int(np.sum(vecs_selection_mask))
    # Prepare output and error code
    signed_angles = np.zeros(n_selected_vecs, dtype=np.float64)
    ierr = ctypes.c_int(0)
    # Setup C wrapper
    clock_hand_angles = lib.clock_hand_angles_for_shift_vectors_c
    clock_hand_angles.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # origins
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # targets
        ctypes.c_int,  # n_dims
        ctypes.c_int,  # n_vecs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # vecs_selection_mask
        ctypes.c_int,  # n_selected_vecs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # selected_axes_for_signed
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # signed_angles
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    clock_hand_angles.restype = None
    # Call Fortran routine
    clock_hand_angles(origins, targets, n_dims, n_vecs, vecs_selection_mask, n_selected_vecs, selected_axes_for_signed, signed_angles, ctypes.byref(ierr))
    # Check for errors
    check_err_code(ierr.value)
    # Mark output as read-only
    _readonly(signed_angles)
    return signed_angles

def relative_axes_changes_from_shift_vector(shift_vector):
    """
    Compute relative axis contributions from a shift vector (RAP space).
    Args:
        shift_vector (array-like): Input vector (1D)
    Returns:
        np.ndarray: Relative axis contributions (sum to 1)
    Raises:
        RuntimeError: If Fortran routine returns error
    """
    # Input validation and conversion
    vec = np.ascontiguousarray(shift_vector, dtype=np.float64)  # Shift vector
    n_dims = len(vec)
    # Prepare output and error code
    contrib = np.zeros(n_dims, dtype=np.float64)
    ierr = ctypes.c_int(0)
    # Setup C wrapper
    relative_axes_changes = lib.relative_axes_changes_from_shift_vector_c
    relative_axes_changes.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # shift_vector
        ctypes.c_int,  # n_dims
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # contrib
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    relative_axes_changes.restype = None
    # Call Fortran routine
    relative_axes_changes(vec, n_dims, contrib, ctypes.byref(ierr))
    # Check for errors
    check_err_code(ierr.value)
    # Mark output as read-only
    _readonly(contrib)
    return contrib

def relative_axes_expression_from_expression_vector(expression_vector):
    """
    Compute relative axis contributions from an expression vector (RAP space).
    Args:
        expression_vector (array-like): Input vector (1D)
    Returns:
        np.ndarray: Relative axis contributions (sum to 1)
    Raises:
        RuntimeError: If Fortran routine returns error
    """
    # Input validation and conversion
    vec = np.ascontiguousarray(expression_vector, dtype=np.float64)  # Expression vector
    n_dims = len(vec)
    # Prepare output and error code
    contrib = np.zeros(n_dims, dtype=np.float64)
    ierr = ctypes.c_int(0)
    # Setup C wrapper
    relative_axes_changes = lib.relative_axes_expression_from_expression_vector_c
    relative_axes_changes.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # expression_vector
        ctypes.c_int,  # n_dims
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # contrib
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    relative_axes_changes.restype = None
    # Call Fortran routine
    relative_axes_changes(vec, n_dims, contrib, ctypes.byref(ierr))
    # Check for errors
    check_err_code(ierr.value)
    # Mark output as read-only
    _readonly(contrib)
    return contrib

def tox_normalize_by_std_dev(input_matrix):
    """
    Normalize gene expression values by standard deviation
    
    Args:
        input_matrix: A numeric matrix with genes as rows and tissues as columns
    
    Returns:
        numpy.ndarray: Normalized matrix with same dimensions as input
    """
    input_matrix = np.asarray(input_matrix, dtype=np.float64)
    n_genes, n_tissues = input_matrix.shape
    
    # Validate input data
    if np.any(np.isnan(input_matrix)):
        raise ValueError(f"Input matrix contains NaN values: {np.sum(np.isnan(input_matrix))}")
    if np.any(np.isinf(input_matrix)):
        raise ValueError(f"Input matrix contains infinite values: {np.sum(np.isinf(input_matrix))}")
    
    # Flatten input and prepare output
    input_flat = np.asfortranarray(input_matrix).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    ierr = ctypes.c_int(0)
    
    # Setup C wrapper
    normalize_c = lib.normalize_by_std_dev_c
    normalize_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_tissues
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # input
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # output
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    normalize_c.restype = None
    
    # Call Fortran routine
    normalize_c(n_genes, n_tissues, input_flat, output_flat, ctypes.byref(ierr))
    check_err_code(ierr.value)
    
    # Reshape and return
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    
    # Mark output as read-only
    _readonly(result)
    return result

def tox_quantile_normalization(input_matrix):
    """
    Quantile normalization of gene expression values
    
    Args:
        input_matrix: A numeric matrix with genes as rows and tissues as columns
    
    Returns:
        numpy.ndarray: Quantile-normalized matrix with same dimensions as input
    """
    input_matrix = np.asarray(input_matrix, dtype=np.float64)
    n_genes, n_tissues = input_matrix.shape
    
    # Flatten input and prepare output arrays
    input_flat = np.asfortranarray(input_matrix).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    temp_col = np.zeros(n_genes, dtype=np.float64)
    rank_means = np.zeros(n_genes, dtype=np.float64)
    perm = np.zeros(n_genes, dtype=np.int32)
    max_stack = max(2 * n_genes, 2)
    stack_left = np.zeros(max_stack, dtype=np.int32)
    stack_right = np.zeros(max_stack, dtype=np.int32)
    ierr = ctypes.c_int(0)
    
    # Setup C wrapper
    quantile_norm_c = lib.quantile_normalization_c
    quantile_norm_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_tissues
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # input
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # output
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # temp_col
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # rank_means
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right
        ctypes.c_int,  # max_stack
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    quantile_norm_c.restype = None
    
    # Call Fortran routine
    quantile_norm_c(n_genes, n_tissues, input_flat, output_flat,
                    temp_col, rank_means, perm, stack_left, stack_right, max_stack, ctypes.byref(ierr))
    check_err_code(ierr.value)
    
    # Reshape and return
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    
    # Mark output as read-only
    _readonly(result)
    return result

def tox_log2_transformation(input_matrix):
    """
    Apply log2(x + 1) transformation to gene expression values
    
    Args:
        input_matrix: A numeric matrix with genes as rows and tissues as columns
    
    Returns:
        numpy.ndarray: Log2-transformed matrix with same dimensions as input
    """
    input_matrix = np.asarray(input_matrix, dtype=np.float64)
    n_genes, n_tissues = input_matrix.shape
    
    # Flatten input and prepare output
    input_flat = np.asfortranarray(input_matrix).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    ierr = ctypes.c_int(0)
    
    # Setup C wrapper
    log2_transform_c = lib.log2_transformation_c
    log2_transform_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_tissues
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # input
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # output
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    log2_transform_c.restype = None
    
    # Call Fortran routine
    log2_transform_c(n_genes, n_tissues, input_flat, output_flat, ctypes.byref(ierr))
    check_err_code(ierr.value)
    
    # Reshape and return
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    
    # Mark output as read-only
    _readonly(result)
    return result

def tox_calculate_tissue_averages(input_matrix, group_starts, group_counts):
    """
    Calculate average expression across replicates for each tissue group
    
    Args:
        input_matrix: A numeric matrix with genes as rows and tissue replicates as columns
        group_starts: Array of starting column indices for each group (1-based for Fortran)
        group_counts: Array of counts for each group
    
    Returns:
        numpy.ndarray: Matrix with genes as rows and averaged tissues as columns
    """
    input_matrix = np.asarray(input_matrix, dtype=np.float64)
    group_starts = np.asarray(group_starts, dtype=np.int32)
    group_counts = np.asarray(group_counts, dtype=np.int32)
    
    n_genes, n_samples = input_matrix.shape
    n_groups = len(group_starts)
    
    if len(group_counts) != n_groups:
        raise ValueError("group_starts and group_counts must have same length")
    
    # Flatten input and prepare output
    input_flat = np.asfortranarray(input_matrix).ravel(order='F')
    output_flat = np.zeros(n_genes * n_groups, dtype=np.float64)
    ierr = ctypes.c_int(0)
    
    # Setup C wrapper
    tiss_avg_c = lib.calc_tiss_avg_c
    tiss_avg_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_groups
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # group_starts
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # group_counts
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # input
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # output
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    tiss_avg_c.restype = None
    
    # Call Fortran routine
    tiss_avg_c(n_genes, n_groups, group_starts, group_counts, input_flat, output_flat, ctypes.byref(ierr))
    check_err_code(ierr.value)
    
    # Reshape and return
    result = output_flat.reshape((n_genes, n_groups), order='F')
    
    # Mark output as read-only
    _readonly(result)
    return result

def tox_normalization_pipeline(input_matrix, group_starts, group_counts):
    """
    Complete normalization pipeline for gene expression data (up to log2(x+1))
    Mirrors Fortran normalization_pipeline (no fold change).

    Args:
        input_matrix: Numeric matrix (genes x tissues)
        group_starts: Integer array, start column index for each replicate group (1-based)
        group_counts: Integer array, number of columns per replicate group

    Returns:
        numpy.ndarray: log2(x+1) normalized expression (genes x groups)
    """
    input_matrix = np.asarray(input_matrix, dtype=np.float64)
    group_starts = np.asarray(group_starts, dtype=np.int32)
    group_counts = np.asarray(group_counts, dtype=np.int32)

    n_genes, n_tissues = input_matrix.shape
    n_grps = len(group_starts)

    # Flatten input and allocate workspace
    input_flat = np.asfortranarray(input_matrix).ravel(order='F')
    buf_stddev = np.zeros(n_genes * n_tissues, dtype=np.float64)
    buf_quant = np.zeros(n_genes * n_tissues, dtype=np.float64)
    buf_avg = np.zeros(n_genes * n_grps, dtype=np.float64)
    buf_log = np.zeros(n_genes * n_grps, dtype=np.float64)
    temp_col = np.zeros(n_genes, dtype=np.float64)
    rank_means = np.zeros(n_genes, dtype=np.float64)
    perm = np.zeros(n_genes, dtype=np.int32)
    max_stack = max(2 * n_genes, 2)
    stack_left = np.zeros(max_stack, dtype=np.int32)
    stack_right = np.zeros(max_stack, dtype=np.int32)
    ierr = ctypes.c_int(0)

    # Setup C wrapper
    normalization_pipeline_c = lib.normalization_pipeline_c
    normalization_pipeline_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_tissues
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # input_flat
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # buf_stddev
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # buf_quant
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # buf_avg
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # buf_log
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # temp_col
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # rank_means
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right
        ctypes.c_int,  # max_stack
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # group_starts
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # group_counts
        ctypes.c_int,  # n_grps
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    normalization_pipeline_c.restype = None

    # Call Fortran routine
    normalization_pipeline_c(
        n_genes, n_tissues, input_flat, buf_stddev, buf_quant, buf_avg, buf_log,
        temp_col, rank_means, perm, stack_left, stack_right, max_stack,
        group_starts, group_counts, n_grps, ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

    # Reshape and return log2(x+1) output
    result = buf_log.reshape((n_genes, n_grps), order='F')
    _readonly(result)
    return result

def tox_calculate_fold_changes(input_matrix, control_cols, condition_cols):
    """
    Calculate log2 fold changes between control and condition columns
    
    Args:
        input_matrix: A numeric matrix with genes as rows and tissues/conditions as columns
        control_cols: Array of control column indices (1-based for Fortran)
        condition_cols: Array of condition column indices (1-based for Fortran)
    
    Returns:
        numpy.ndarray: Matrix with genes as rows and fold change values as columns
    """
    input_matrix = np.asarray(input_matrix, dtype=np.float64)
    control_cols = np.asarray(control_cols, dtype=np.int32)
    condition_cols = np.asarray(condition_cols, dtype=np.int32)
    
    n_genes, n_samples = input_matrix.shape
    n_pairs = len(control_cols)
    
    if len(condition_cols) != n_pairs:
        raise ValueError("control_cols and condition_cols must have same length")
    
    # Flatten input and prepare output
    input_flat = np.asfortranarray(input_matrix).ravel(order='F')
    output_flat = np.zeros(n_genes * n_pairs, dtype=np.float64)
    ierr = ctypes.c_int(0)
    
    # Setup C wrapper
    fchange_c = lib.calc_fchange_c
    fchange_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_samples
        ctypes.c_int,  # n_pairs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # control_cols
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # condition_cols
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # input
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # output
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    fchange_c.restype = None
    
    # Call Fortran routine
    fchange_c(n_genes, n_samples, n_pairs, control_cols, condition_cols, input_flat, output_flat, ctypes.byref(ierr))
    check_err_code(ierr.value)
    
    # Reshape and return
    result = output_flat.reshape((n_genes, n_pairs), order='F')
    
    # Mark output as read-only
    _readonly(result)
    return result

def tox_calculate_tissue_versatility(expression_vectors, vector_selection, axis_selection):
    """
    Calculate Tissue Versatility
    
    Computes normalized tissue versatility for selected expression vectors.
    The metric is based on the angle between each gene expression vector and the space diagonal.
    Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
    This function automatically checks for errors and throws informative exceptions.
    
    Args:
        expression_vectors: Matrix where each column is a gene expression vector (n_axes x n_vectors)
        vector_selection: Boolean or integer array indicating which vectors to process (length n_vectors)
        axis_selection: Boolean or integer array indicating which axes to include in calculation (length n_axes)
    
    Returns:
        dict: Dictionary containing:
            - tissue_versatilities: Normalized tissue versatility values [0,1] for selected vectors
            - tissue_angles_deg: Angles in degrees [0,90] for selected vectors
            - n_selected_vectors: Number of vectors processed
            - n_selected_axes: Number of axes used in calculation

    """
    # Input validation
    if not isinstance(expression_vectors, np.ndarray):
        raise ValueError("expression_vectors must be a numpy array")
    
    if expression_vectors.ndim != 2:
        raise ValueError("expression_vectors must be a 2D array")
    
    # Convert inputs to numpy arrays
    expression_vectors = np.asarray(expression_vectors, dtype=np.float64)
    vector_selection = np.asarray(vector_selection, dtype=bool)
    axis_selection = np.asarray(axis_selection, dtype=bool)
    
    # Get dimensions
    n_axes, n_vectors = expression_vectors.shape
    
    # Validate dimensions
    if len(vector_selection) != n_vectors:
        raise ValueError("vector_selection length must match number of columns in expression_vectors")
    if len(axis_selection) != n_axes:
        raise ValueError("axis_selection length must match number of rows in expression_vectors")
    
    # Calculate counts
    n_selected_vectors = int(np.sum(vector_selection))
    n_selected_axes = int(np.sum(axis_selection))
    
    # Ensure arrays have correct dtype and memory layout
    expr_f = np.asfortranarray(expression_vectors, dtype=np.float64)
    select_vec = np.ascontiguousarray(vector_selection.astype(np.int32))
    select_axes = np.ascontiguousarray(axis_selection.astype(np.int32))
    
    # Prepare output arrays
    tissue_versatilities = np.empty(n_selected_vectors, dtype=np.float64)
    tissue_angles_deg = np.empty(n_selected_vectors, dtype=np.float64)
    ierr = ctypes.c_int(0)
    
    # Setup C/Fortran wrapper with proper type annotations
    tv = lib.compute_tissue_versatility_c
    tv.argtypes = [
        ctypes.c_int,  # n_axes
        ctypes.c_int,  # n_vectors
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # expression_vectors (Fortran order)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # exp_vecs_selection_index
        ctypes.c_int,  # n_selected_vectors
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # axes_selection
        ctypes.c_int,  # n_selected_axes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # tissue_versatilities
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # tissue_angles_deg
        ctypes.POINTER(ctypes.c_int),  # ierr
    ]
    tv.restype = None
    
    # Call the C/Fortran wrapper
    tv(
        n_axes,
        n_vectors,
        expr_f,
        select_vec,
        n_selected_vectors,
        select_axes,
        n_selected_axes,
        tissue_versatilities,
        tissue_angles_deg,
        ctypes.byref(ierr)
    )
    
    # Check for errors and throw informative messages
    check_err_code(ierr.value)
    
    # Mark outputs as read-only
    _readonly(tissue_versatilities, tissue_angles_deg)
    
    # Return structured result (no ierr since we checked for errors)
    return {
        'tissue_versatilities': tissue_versatilities,
        'tissue_angles_deg': tissue_angles_deg
    }


def tox_euclidean_distance(vec1, vec2):
    """
    Calculate Euclidean distance between two vectors
    
    Args:
        vec1: First vector (numpy array)
        vec2: Second vector (numpy array)
    
    Returns:
        float: Euclidean distance between the vectors
    """
    # Input validation
    if not isinstance(vec1, np.ndarray):
        vec1 = np.asarray(vec1, dtype=np.float64)
    if not isinstance(vec2, np.ndarray):
        vec2 = np.asarray(vec2, dtype=np.float64)
    
    vec1 = np.asarray(vec1, dtype=np.float64)
    vec2 = np.asarray(vec2, dtype=np.float64)
    
    if len(vec1) != len(vec2):
        raise ValueError("Vectors must have the same length")
    if len(vec1) == 0:
        raise ValueError("Vectors cannot be empty")
    if not np.issubdtype(vec1.dtype, np.number) or not np.issubdtype(vec2.dtype, np.number):
        raise ValueError("Vectors must be numeric")
    
    # Ensure contiguous arrays
    vec1 = np.ascontiguousarray(vec1, dtype=np.float64)
    vec2 = np.ascontiguousarray(vec2, dtype=np.float64)
    
    # Prepare output
    result = ctypes.c_double(0.0)
    
    # Setup C wrapper
    euclidean_distance_c = lib.euclidean_distance_c
    euclidean_distance_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double)
    ]
    euclidean_distance_c.restype = None
    
    # Call Fortran routine
    euclidean_distance_c(vec1, vec2, len(vec1), ctypes.byref(result))
    
    return result.value


def tox_distance_to_centroid(genes, centroids, gene_to_fam, d):
    """
    Calculate distance from each gene to its family centroid
    
    Args:
        genes: Gene expression data as flat array (n_genes * d elements)
        centroids: Family centroids as flat array (n_families * d elements)
        gene_to_fam: Gene-to-family mapping (0 = no family, >0 = family index)
        d: Number of dimensions
    
    Returns:
        numpy.ndarray: Distances from each gene to its centroid (-1 for invalid families)
    """
    # Input validation
    if not isinstance(genes, np.ndarray):
        genes = np.asarray(genes, dtype=np.float64)
    if not isinstance(centroids, np.ndarray):
        centroids = np.asarray(centroids, dtype=np.float64)
    if not isinstance(gene_to_fam, np.ndarray):
        gene_to_fam = np.asarray(gene_to_fam, dtype=np.int32)
    
    genes = np.asarray(genes, dtype=np.float64)
    centroids = np.asarray(centroids, dtype=np.float64)
    gene_to_fam = np.asarray(gene_to_fam, dtype=np.int32)
    d = int(d)
    
    # Calculate dimensions
    n_genes = len(genes) // d
    n_families = len(centroids) // d
    
    # Validate dimensions
    if len(genes) % d != 0:
        raise ValueError("Length of genes must be divisible by d")
    if len(centroids) % d != 0:
        raise ValueError("Length of centroids must be divisible by d")
    if len(gene_to_fam) != n_genes:
        raise ValueError("Length of gene_to_fam must equal number of genes")
    if np.any(gene_to_fam < 0):
        raise ValueError("gene_to_fam indices must be between 0 and n_families (0 = no family assignment)")
    
    # Ensure contiguous arrays
    genes = np.ascontiguousarray(genes, dtype=np.float64)
    centroids = np.ascontiguousarray(centroids, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    
    # Prepare output array
    distances = np.zeros(n_genes, dtype=np.float64)
    
    # Setup C wrapper
    distance_to_centroid_c = lib.distance_to_centroid_c
    distance_to_centroid_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # genes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # centroids
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        ctypes.c_int   # d
    ]
    distance_to_centroid_c.restype = None
    
    # Call Fortran routine
    distance_to_centroid_c(n_genes, n_families, genes, centroids, gene_to_fam, distances, d)
    
    # Mark output as read-only
    _readonly(distances)
    return distances


def tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff):
    """
    LOESS smoothing in 2D
    
    Args:
        x_ref: Reference x values
        y_ref: Reference y values  
        indices_used: Indices of points to use (1-based for Fortran)
        x_query: Query x values where to compute smoothed y
        kernel_sigma: Kernel bandwidth parameter
        kernel_cutoff: Kernel cutoff parameter
    
    Returns:
        numpy.ndarray: Smoothed y values at query points
    """
    # Input validation and conversion
    x_ref = np.ascontiguousarray(x_ref, dtype=np.float64)
    y_ref = np.ascontiguousarray(y_ref, dtype=np.float64)
    indices_used = np.ascontiguousarray(indices_used, dtype=np.int32)
    x_query = np.ascontiguousarray(x_query, dtype=np.float64)
    
    n_total = len(x_ref)
    n_target = len(x_query)
    n_used = len(indices_used)
    
    if len(y_ref) != n_total:
        raise ValueError("x_ref and y_ref must have same length")
    
    # Validate indices bounds
    if n_used > 0:
        if np.any(indices_used < 1) or np.any(indices_used > n_total):
            raise ValueError("indices_used must be 1-based and within bounds")
    
    # Validate parameters
    if kernel_sigma < 0:
        raise ValueError("kernel_sigma must be non-negative")
    if kernel_cutoff < 0:
        raise ValueError("kernel_cutoff must be non-negative")
    
    # Prepare output and error code
    y_out = np.zeros(n_target, dtype=np.float64)
    ierr = ctypes.c_int(0)
    
    # Setup C wrapper
    loess_smooth_2d_c = lib.loess_smooth_2d_c
    loess_smooth_2d_c.argtypes = [
        ctypes.c_int,  # n_total
        ctypes.c_int,  # n_target
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # x_ref
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # y_ref
        np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),    # indices_used
        ctypes.c_int,  # n_used
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # x_query
        ctypes.c_double,  # kernel_sigma
        ctypes.c_double,  # kernel_cutoff
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # y_out
        ctypes.POINTER(ctypes.c_int),  # ierr
    ]
    loess_smooth_2d_c.restype = None
    
    # Call Fortran routine
    loess_smooth_2d_c(
        n_total, n_target,
        x_ref, y_ref, indices_used, n_used, x_query,
        ctypes.c_double(kernel_sigma), ctypes.c_double(kernel_cutoff),
        y_out, ctypes.byref(ierr)
    )
    
    # Check for errors
    check_err_code(ierr.value)
    
    # Mark output as read-only
    _readonly(y_out)
    return y_out


def tox_compute_family_scaling(distances, gene_to_fam):
    """
    Compute family scaling factors for outlier detection
    
    Args:
        distances: Gene distances to family centroids
        gene_to_fam: Gene-to-family mapping
    
    Returns:
        dict: Dictionary containing scaling factors and intermediate results
    """
    distances = np.ascontiguousarray(distances, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    
    n_genes = len(distances)
    n_families = int(np.max(gene_to_fam)) if len(gene_to_fam) > 0 else 0
    
    if len(gene_to_fam) != n_genes:
        raise ValueError("distances and gene_to_fam must have same length")
    
    # Prepare output arrays
    dscale = np.zeros(n_families, dtype=np.float64)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    indices_used = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    # Setup C wrapper
    compute_family_scaling_c = lib.compute_family_scaling_c
    compute_family_scaling_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # dscale
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_x
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_y
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # indices_used
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # error_code
    ]
    compute_family_scaling_c.restype = None
    
    # Call Fortran routine
    compute_family_scaling_c(
        n_genes, n_families, distances, gene_to_fam,
        dscale, loess_x, loess_y, indices_used, error_code
    )
    
    # Check for errors
    check_err_code(error_code[0])
    
    # Mark outputs as read-only
    _readonly(dscale, loess_x, loess_y, indices_used)
    
    return {
        'dscale': dscale,
        'loess_x': loess_x,
        'loess_y': loess_y,
        'indices_used': indices_used
    }


def tox_compute_family_scaling_expert(distances, gene_to_fam, perm_tmp, stack_left_tmp, 
                                 stack_right_tmp, family_distances):
    """
    Expert version of compute_family_scaling with user-provided work arrays
    
    This version requires pre-allocated work arrays for maximum performance and control.
    Use this when you need fine-grained control over memory allocation or are calling
    this function many times in a tight loop.
    
    Args:
        distances: Gene distances to family centroids
        gene_to_fam: Gene-to-family mapping
        perm_tmp: Pre-allocated permutation array for sorting (n_genes)
        stack_left_tmp: Pre-allocated stack array for sorting (n_genes)
        stack_right_tmp: Pre-allocated stack array for sorting (n_genes)
        family_distances: Pre-allocated work array for family distances (n_genes)
    
    Returns:
        dict: Dictionary containing scaling factors and intermediate results
    """
    distances = np.ascontiguousarray(distances, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    perm_tmp = np.ascontiguousarray(perm_tmp, dtype=np.int32)
    stack_left_tmp = np.ascontiguousarray(stack_left_tmp, dtype=np.int32)
    stack_right_tmp = np.ascontiguousarray(stack_right_tmp, dtype=np.int32)
    family_distances = np.ascontiguousarray(family_distances, dtype=np.float64)
    
    n_genes = len(distances)
    n_families = int(np.max(gene_to_fam)) if len(gene_to_fam) > 0 else 0
    
    if len(gene_to_fam) != n_genes:
        raise ValueError("distances and gene_to_fam must have same length")
    if len(perm_tmp) != n_genes:
        raise ValueError("perm_tmp must have same length as distances")
    if len(stack_left_tmp) != n_genes:
        raise ValueError("stack_left_tmp must have same length as distances")
    if len(stack_right_tmp) != n_genes:
        raise ValueError("stack_right_tmp must have same length as distances")
    if len(family_distances) != n_genes:
        raise ValueError("family_distances must have same length as distances")
    
    # Prepare output arrays
    dscale = np.zeros(n_families, dtype=np.float64)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    indices_used = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    # Setup C wrapper
    compute_family_scaling_expert_c = lib.compute_family_scaling_expert_c
    compute_family_scaling_expert_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # dscale
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_x
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_y
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # indices_used
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm_tmp
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left_tmp
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right_tmp
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # family_distances
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # error_code
    ]
    compute_family_scaling_expert_c.restype = None
    
    # Call Fortran routine
    compute_family_scaling_expert_c(
        n_genes, n_families, distances, gene_to_fam,
        dscale, loess_x, loess_y, indices_used,
        perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code
    )
    
    # Check for errors
    check_err_code(error_code[0])
    
    # Mark outputs as read-only
    _readonly(dscale, loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances)
    
    return {
        'dscale': dscale,
        'loess_x': loess_x,
        'loess_y': loess_y,
        'indices_used': indices_used,
        'perm_tmp': perm_tmp,
        'stack_left_tmp': stack_left_tmp,
        'stack_right_tmp': stack_right_tmp,
        'family_distances': family_distances
    }


def tox_compute_rdi(distances, gene_to_fam, dscale):
    """
    Compute Relative Distance Index (RDI) for outlier detection
    
    Args:
        distances: Gene distances to centroids
        gene_to_fam: Gene-to-family mapping
        dscale: Family scaling factors
    
    Returns:
        numpy.ndarray: RDI values for each gene
    """
    distances = np.ascontiguousarray(distances, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    dscale = np.ascontiguousarray(dscale, dtype=np.float64)
    
    n_genes = len(distances)
    n_families = len(dscale)
    
    # Prepare output and work arrays
    rdi = np.zeros(n_genes, dtype=np.float64)
    sorted_rdi = np.zeros(n_genes, dtype=np.float64)
    perm = np.zeros(n_genes, dtype=np.int32)
    stack_left = np.zeros(n_genes, dtype=np.int32)
    stack_right = np.zeros(n_genes, dtype=np.int32)
    
    # Setup C wrapper
    compute_rdi_c = lib.compute_rdi_c
    compute_rdi_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # dscale
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # rdi
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # sorted_rdi
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right
    ]
    compute_rdi_c.restype = None
    
    # Call Fortran routine
    compute_rdi_c(n_genes, n_families, distances, gene_to_fam, dscale, 
                  rdi, sorted_rdi, perm, stack_left, stack_right)
    
    # Mark output as read-only
    _readonly(rdi)
    return rdi


def tox_identify_outliers(rdi, threshold=None, percentile=95.0):
    """
    Identify outliers based on RDI percentile or threshold
    
    Args:
        rdi: Relative Distance Index values
        threshold: Fixed RDI threshold (if None, uses percentile)
        percentile: Percentile threshold for outlier detection (default: 95 for top 5%)
    
    Returns:
        dict: Dictionary containing:
            - outliers: Boolean array indicating outliers
            - threshold: Threshold value used for detection
    """
    rdi = np.ascontiguousarray(rdi, dtype=np.float64)
    n_genes = len(rdi)
    
    # Prepare sorted RDI (copy and filter out negatives)
    sorted_rdi = rdi.copy()
    sorted_rdi[sorted_rdi < 0] = 0.0  # Filter out error values
    sorted_rdi.sort()  # Sort in ascending order
    
    # Prepare output arrays
    outliers_int = np.zeros(n_genes, dtype=np.int32)
    threshold_out = ctypes.c_double(0.0)
    
    # Setup C wrapper
    identify_outliers_c = lib.identify_outliers_c
    identify_outliers_c.argtypes = [
        ctypes.c_int,  # n_genes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # rdi
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # sorted_rdi
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # is_outlier_int
        ctypes.POINTER(ctypes.c_double),  # threshold (output)
        ctypes.c_double,  # percentile
    ]
    identify_outliers_c.restype = None
    
    # Call Fortran routine
    identify_outliers_c(n_genes, rdi, sorted_rdi, outliers_int, 
                        ctypes.byref(threshold_out), ctypes.c_double(percentile))
    
    # Mark output as read-only
    _readonly(outliers_int)
    
    return {
        'outliers': outliers_int,
        'threshold': threshold_out.value
    }


def tox_detect_outliers(distances, gene_to_fam, percentile=95.0):
    """
    Complete outlier detection pipeline
    
    Args:
        distances: Gene distances to centroids
        gene_to_fam: Gene-to-family mapping
        percentile: Percentile threshold for outlier detection (default: 95 for top 5%)
    
    Returns:
        dict: Dictionary containing outliers and intermediate results
    """
    distances = np.ascontiguousarray(distances, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    
    n_genes = len(distances)
    n_families = int(np.max(gene_to_fam)) if len(gene_to_fam) > 0 else 0
    
    # Prepare work arrays
    work_array = np.zeros(n_genes, dtype=np.float64)
    perm = np.zeros(n_genes, dtype=np.int32)
    stack_left = np.zeros(n_genes, dtype=np.int32)
    stack_right = np.zeros(n_genes, dtype=np.int32)
    
    # Prepare output arrays
    outliers_int = np.zeros(n_genes, dtype=np.int32)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    loess_n = np.zeros(n_families, dtype=np.int32)
    error_code = ctypes.c_int(0)
    
    # Setup C wrapper
    detect_outliers_c = lib.detect_outliers_c
    detect_outliers_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # work_array
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # outliers_int
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_x
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_y
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # loess_n
        ctypes.POINTER(ctypes.c_int),  # error_code
        ctypes.c_double,  # percentile
    ]
    detect_outliers_c.restype = None
    
    # Call Fortran routine
    detect_outliers_c(
        n_genes, n_families, distances, gene_to_fam,
        work_array, perm, stack_left, stack_right,
        outliers_int, loess_x, loess_y, loess_n,
        ctypes.byref(error_code), ctypes.c_double(percentile)
    )
    
    # Check for errors
    check_err_code(error_code.value)
    
    # Mark outputs as read-only
    _readonly(outliers_int, loess_x, loess_y, loess_n)
    
    return {
        'outliers': outliers_int,
        'loess_x': loess_x,
        'loess_y': loess_y,
        'loess_n': loess_n
    }


def tox_which(cond):
    """
    'which' utility for Python, like in R/MATLAB.
    Returns indices of TRUE elements in a logical array.
    
    Args:
        cond: Array of boolean/integer values (0/1)
    
    Returns:
        tuple: (idx_out, m_out) where idx_out contains 1-based indices and m_out is count
    """
    cond = np.ascontiguousarray(cond, dtype=np.int32)
    
    # Check for invalid input
    if np.any(np.isnan(cond.astype(float))):
        raise ValueError("Input contains NaN values")
    if np.any(np.isinf(cond.astype(float))):
        raise ValueError("Input contains infinite values")
    if not np.all((cond == 0) | (cond == 1)):
        raise ValueError("Input must contain only 0 and 1 values")
    
    n = cond.size
    idx_out = np.zeros(n, dtype=np.int32)
    m_max = n
    m_out = np.zeros(1, dtype=np.int32)
    error_code = ctypes.c_int(0)
    
    # Setup C wrapper
    which_c = lib.which_c
    which_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # mask
        ctypes.c_int,                 # n
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # idx_out
        ctypes.c_int,                 # m_max
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # m_out
        ctypes.POINTER(ctypes.c_int)  # error_code
    ]
    which_c.restype = None
    
    # Call Fortran subroutine
    which_c(
        cond, ctypes.c_int(n), idx_out, ctypes.c_int(m_max), m_out,
        ctypes.byref(error_code)
    )
    
    # Check for errors
    check_err_code(error_code.value)
    
    # Mark output as read-only
    _readonly(idx_out)

    return idx_out, int(m_out[0])

def tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid):
    """
    Calculate Shift Vector Field 

    Computes the shift vector field for each gene expression vector based on its family centroid.
    The shift vector is defined as the difference between the gene expression vector and its corresponding family centroid,
    starting at the expression vector and pointing to its family centroid.
    This function automatically checks for errors and throws informative exceptions.
    
    Args:
        expression_vectors: Matrix where each column is a gene expression vector (n_axes x n_vectors)
        family_centroids: Matrix where each column is a family centroid vector (n_axes x n_families)
        gene_to_centroid: Array mapping each gene to its corresponding family centroid index in family_centroids with length n_vectors (1 based for fortran)

    Returns:
        shift_vectors: The computed shift vectors for each gene expression vector

    """
    
    # Input validation
    if not isinstance(expression_vectors, np.ndarray):
        raise ValueError("expression_vectors must be a numpy array")
    
    if expression_vectors.ndim != 2:
        raise ValueError("expression_vectors must be a 2D array")

    if not isinstance(family_centroids, np.ndarray):
        raise ValueError("family_centroids must be a numpy array")

    if family_centroids.ndim != 2:
        raise ValueError("family_centroids must be a 2D array")
    
    # Convert inputs to numpy arrays
    expression_vectors = np.asarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asarray(family_centroids, dtype=np.float64)
    gene_to_centroid = np.asarray(gene_to_centroid, dtype=np.int32)
    
    # Get dimensions
    n_axes_genes, n_vectors = expression_vectors.shape
    n_axes_centroids, n_families = family_centroids.shape
    
    # Validate length of gene_to_centroid
    if len(gene_to_centroid) != n_vectors:
        raise ValueError("gene_to_centroid length must match number of columns in expression_vectors")
    
    # Validate dimensions
    if n_axes_genes != n_axes_centroids:
        raise ValueError("family_centroids must have the same number of axes as expression_vectors")

    # Validate correct mapping between genes and centroids
    if not np.all((gene_to_centroid > 0) & (gene_to_centroid <= n_families)):
        raise ValueError("gene_to_centroid contains invalid family IDs")

    # Ensure arrays have correct dtype and memory layout
    expr_v = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_c = np.asfortranarray(family_centroids, dtype=np.float64)
    gene_c = np.ascontiguousarray(gene_to_centroid.astype(np.int32))
    
    # Prepare output arrays
    shift_vectors = np.empty((2*n_axes_genes, n_vectors), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # Setup C/Fortran wrapper with proper type annotations
    sv = lib.compute_shift_vector_field_c
    sv.argtypes = [
        ctypes.c_int,  # n_axes_genes
        ctypes.c_int,  # n_vectors
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # expression_vectors (Fortran order)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # family_centroids (Fortran order)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_centroid
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # shift_vectors
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    sv.restype = None
    
    # Call the C/Fortran wrapper
    sv(
        n_axes_genes,
        n_vectors,
        n_families,
        expr_v,
        family_c,
        gene_c,
        shift_vectors,
        ctypes.byref(ierr)
    )

    # Check for errors and throw informative messages
    check_err_code(ierr.value)

    # Mark outputs as read-only
    _readonly(shift_vectors)
    
    # Return result (no ierr since we checked for errors)
    return shift_vectors


def tox_group_centroid(expression_vectors, gene_to_family, n_families, mode, ortholog_set = None):
    """
    Computes expression centroids for groups of genes.

    Computes the centroids for each gene family based on the expression vectors of its member genes.
    This function automatically checks for errors and throws informative exceptions.

    Args:
        vectors : np.ndarray
            A 2D NumPy array (n_axes x n_genes) of gene expression vectors.
        gene_to_family_map : np.ndarray
            A 1D NumPy array of length n_genes, mapping each gene to a family ID.
        n_families : int
            The total number of unique families.
        mode : str
            The calculation mode. 'all' or 'orthologs'.
        ortholog_set : np.ndarray
            (Optional) A 1D boolean NumPy array of length n_genes, indicating ortholog membership (only required in 'orthologs' mode).

    Returns:
        np.ndarray
            A read-only (n_axes x n_families) NumPy array containing the computed centroids.
    """

    # 1) Validate and prepare inputs
    if not isinstance(expression_vectors, np.ndarray) or expression_vectors.ndim != 2:
        raise ValueError("`vectors` must be a 2D NumPy array.")
    n_axes, n_genes = expression_vectors.shape

    if mode != 'all' and mode != 'orthologs':
        raise ValueError("'mode' must be either 'all' or 'orthologs'.")
    if mode == 'orthologs': 
        if ortholog_set is None:
            raise ValueError("`ortholog_set` must be provided when mode is 'orthologs'.")
    else:
        ortholog_set = np.ones(n_genes, dtype=np.int32, order="F")

    vecs_f = np.asarray(expression_vectors, dtype=np.float64, order="F")
    g2f_map_f = np.asarray(gene_to_family, dtype=np.int32, order="F")
    ortho_set_int_f = np.asarray(ortholog_set, dtype=np.int32, order="F")
    
    if g2f_map_f.size != n_genes:
        raise ValueError("`gene_to_family` must be a 1D NumPy array of size n_genes.")
    if ortho_set_int_f.size != n_genes:
        raise ValueError("`ortholog_set` must be a 1D NumPy array of size n_genes.")

    # 2) Prepare output buffers and mode flag
    centroids_out = np.zeros((n_axes, n_families), dtype=np.float64, order="F")
    selected_indices = np.zeros(n_genes, dtype=np.int32, order="F")
    ierr = ctypes.c_int(0)

    # 3) Setup C-interface signature
    group_centroid_c = lib.group_centroid_c
    group_centroid_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # expression_vectors
        ctypes.c_int,                                                   # n_axes
        ctypes.c_int,                                                   # n_genes
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),   # gene_to_family
        ctypes.c_int,                                                   # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # centroid_matrix (out)
        ctypes.c_char * 10,                                             # mode (character array)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),   # ortholog_set (as int array)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),   # selected_indices
        ctypes.c_int,                                                   # selected_indices_len
        ctypes.POINTER(ctypes.c_int)                                    # ierr
    ]
    group_centroid_c.restype = None

    # 4) Call the Fortran routine
    mode_buffer = ctypes.create_string_buffer(mode.encode('utf-8'), size=10)
    group_centroid_c(
        vecs_f,
        n_axes,
        n_genes,
        g2f_map_f,
        n_families,
        centroids_out,
        mode_buffer,
        ortho_set_int_f,
        selected_indices,
        n_genes,
        ctypes.byref(ierr)
    )

    # Check for errors and throw informative messages
    check_err_code(ierr.value)

    # 5) Mark output as read-only and return
    _readonly(centroids_out)
    return centroids_out

def tox_mean_vector(expression_vectors, gene_indices):
    """
    Compute the element-wise mean for a given set of gene expression vectors.

    This function wraps the Fortran subroutine `mean_vector_c`
    to compute the centroid (mean vector) for a selected set of genes.

    Args:
        expression_vectors: 2D numpy array (n_axes x n_genes) of gene expression vectors.
        gene_indices: 1D numpy array of column indices of selected genes (1-based).

    Returns:
        numpy.ndarray: 1D array of length n_axes representing the computed centroid.
    """
    # Validate inputs
    if not isinstance(expression_vectors, np.ndarray) or expression_vectors.ndim != 2:
        raise ValueError("expression_vectors must be a 2D numpy array.")
    n_axes, n_genes = expression_vectors.shape

    gene_indices = np.asarray(gene_indices, dtype=np.int32)
    n_selected_genes = len(gene_indices)
    if np.any(gene_indices < 1) or np.any(gene_indices > n_genes):
        raise ValueError("gene_indices must be integer indices between 1 and n_genes (1-based).")

    expr_f = np.asfortranarray(expression_vectors, dtype=np.float64)
    centroid_col = np.zeros(n_axes, dtype=np.float64)
    ierr = ctypes.c_int(0)

    # Setup C wrapper
    mean_vector_c = lib.mean_vector_c
    mean_vector_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # expression_vectors
        ctypes.c_int,                                                   # n_axes
        ctypes.c_int,                                                   # n_genes
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),   # gene_indices
        ctypes.c_int,                                                   # n_selected_genes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"), # centroid_col (out)
        ctypes.POINTER(ctypes.c_int)                                    # ierr
    ]
    mean_vector_c.restype = None

    # Call Fortran routine
    mean_vector_c(
        expr_f,
        n_axes,
        n_genes,
        gene_indices,
        n_selected_genes,
        centroid_col,
        ctypes.byref(ierr)
    )

    # Error handling
    check_err_code(ierr.value)

    # Mark output as read-only
    _readonly(centroid_col)
    return centroid_col
