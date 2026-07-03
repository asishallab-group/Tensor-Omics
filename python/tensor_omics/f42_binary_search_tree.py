from .error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


def get_sorted_value(
        values,
        sorted_indices,
        position
        ):
    """
    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (n_values_elements,) in column-major layout (order='F')
        Input real array
    sorted_indices : np.ndarray[np.int32] of shape (n_sorted_indices_elements,) in column-major layout (order='F')
        Permutation index array
    position : int
        Sorted position (1-based)

    Returns
    -------
    sorted_value : float


    Notes
    -----
    Get the value at the sorted position.
    """

    # ensure all array inputs are numpy arrays
    values = np.ascontiguousarray(values, dtype=np.float64)
    sorted_indices = np.ascontiguousarray(sorted_indices, dtype=np.int32)

    # extract dimension arguments
    n_values_elements = values.shape[0]
    n_sorted_indices_elements = sorted_indices.shape[0]


    # Create temporaries and/or outputs
    position = ctypes.c_int(position)
    ierr = ctypes.c_int(0)
    sorted_value = ctypes.c_double(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.get_sorted_value_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double)
    )
    tox.get_sorted_value_c.restype = None

    tox.get_sorted_value_c(
        values,
        ctypes.byref(ctypes.c_int(n_values_elements)),
        sorted_indices,
        ctypes.byref(ctypes.c_int(n_sorted_indices_elements)),
        ctypes.byref(position),
        ctypes.byref(ierr),
        ctypes.byref(sorted_value)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return sorted_value.value


def build_bst_index(
        values
        ):
    """
    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (num_values,) in column-major layout (order='F')
        Input real array to be indexed

    Returns
    -------
    results : dict
        sorted_indices : np.ndarray[np.int32] of shape (num_values,) in column-major layout (order='F')
            Output permutation index,
        tmp_left_stack : np.ndarray[np.int32] of shape (num_values,) in column-major layout (order='F')
            Manual stack for left indices,
        tmp_right_stack : np.ndarray[np.int32] of shape (num_values,) in column-major layout (order='F')
            Manual stack for right indices


    Notes
    -----
    Build the BST index by sorting indices using values in x.
    """

    # ensure all array inputs are numpy arrays
    values = np.ascontiguousarray(values, dtype=np.float64)

    # extract dimension arguments
    num_values = values.shape[0]


    # Create temporaries and/or outputs
    sorted_indices = np.empty((num_values,), dtype=np.int32, order='F')
    tmp_left_stack = np.empty((num_values,), dtype=np.int32, order='F')
    tmp_right_stack = np.empty((num_values,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.build_bst_index_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.build_bst_index_c.restype = None

    tox.build_bst_index_c(
        values,
        ctypes.byref(ctypes.c_int(num_values)),
        sorted_indices,
        tmp_left_stack,
        tmp_right_stack,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    sorted_indices.setflags(write=False)

    return sorted_indices


def bst_range_query(
        values,
        sorted_indices,
        lower_bound,
        upper_bound
        ):
    """
    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (num_values,) in column-major layout (order='F')
        Input real array
    sorted_indices : np.ndarray[np.int32] of shape (num_values,) in column-major layout (order='F')
        Permutation index array (sorted)
    lower_bound : float
        Lower bound of range (inclusive)
    upper_bound : float
        Upper bound of range (inclusive)

    Returns
    -------
    results : dict
        output_indices : np.ndarray[np.int32] of shape (num_matches,) in column-major layout (order='F')
            Output array of matching indicesThe first `num_matches` elements will hold the results.,
        num_matches : int
            Number of matches found


    Notes
    -----
    Perform a 1D range query over the sorted index.
    """

    # ensure all array inputs are numpy arrays
    values = np.ascontiguousarray(values, dtype=np.float64)
    sorted_indices = np.ascontiguousarray(sorted_indices, dtype=np.int32)

    # extract dimension arguments
    num_values = values.shape[0]


    # Create temporaries and/or outputs
    lower_bound = ctypes.c_double(lower_bound)
    upper_bound = ctypes.c_double(upper_bound)
    output_indices = np.empty((num_values,), dtype=np.int32, order='F')
    num_matches = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.bst_range_query_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.bst_range_query_c.restype = None

    tox.bst_range_query_c(
        values,
        sorted_indices,
        ctypes.byref(ctypes.c_int(num_values)),
        ctypes.byref(lower_bound),
        ctypes.byref(upper_bound),
        output_indices,
        ctypes.byref(num_matches),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    output_indices.setflags(write=False)

    return {
        "output_indices": output_indices[..., :num_matches.value].copy(),
        "num_matches": num_matches.value
    }
