from .error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


def build_kd_index(
        points,
        dimension_order
        ):
    """
    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (num_dimensions, num_points) in column-major layout (order='F')
        Data points
    dimension_order : np.ndarray[np.int32] of shape (num_dimensions,) in column-major layout (order='F')
        Dimension order (by variance)

    Returns
    -------
    results : dict
        kd_indices : np.ndarray[np.int32] of shape (num_points,) in column-major layout (order='F')
            Output index array (k-d tree order),
        tmp_workspace : np.ndarray[np.int32] of shape (num_points,) in column-major layout (order='F')
            Workspace array,
        tmp_value_buffer : np.ndarray[np.float64] of shape (num_points,) in column-major layout (order='F')
            Workspace for sorting,
        tmp_permutation : np.ndarray[np.int32] of shape (num_points,) in column-major layout (order='F')
            Workspace for sorting,
        tmp_left_stack : np.ndarray[np.int32] of shape (num_points,) in column-major layout (order='F')
            Workspace for sorting,
        tmp_right_stack : np.ndarray[np.int32] of shape (num_points,) in column-major layout (order='F')
            Workspace for sorting,
        tmp_recursion_stack : np.ndarray[np.int32] of shape (3, num_points) in column-major layout (order='F')
            Stack for l, r, depth


    Notes
    -----
    Build a k-d tree index using a stack-based, non-recursive approach.Initialize kd_indices to 1:num_points (original indices)Choose split dimension by cycling through dimension_orderFind median indexPartition kd_indices(left_idx:right_idx) by points(current_dim, kd_indices(:))Push right and left intervals onto stack
    """

    # ensure all array inputs are numpy arrays
    points = np.asfortranarray(points, dtype=np.float64)
    dimension_order = np.ascontiguousarray(dimension_order, dtype=np.int32)

    # extract dimension arguments
    num_dimensions = points.shape[0]
    num_points = points.shape[1]


    # Create temporaries and/or outputs
    kd_indices = np.empty((num_points,), dtype=np.int32, order='F')
    tmp_workspace = np.empty((num_points,), dtype=np.int32, order='F')
    tmp_value_buffer = np.empty((num_points,), dtype=np.float64, order='F')
    tmp_permutation = np.empty((num_points,), dtype=np.int32, order='F')
    tmp_left_stack = np.empty((num_points,), dtype=np.int32, order='F')
    tmp_right_stack = np.empty((num_points,), dtype=np.int32, order='F')
    tmp_recursion_stack = np.empty((3, num_points), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.build_kd_index_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.build_kd_index_c.restype = None

    tox.build_kd_index_c(
        points,
        ctypes.byref(ctypes.c_int(num_dimensions)),
        ctypes.byref(ctypes.c_int(num_points)),
        kd_indices,
        dimension_order,
        tmp_workspace,
        tmp_value_buffer,
        tmp_permutation,
        tmp_left_stack,
        tmp_right_stack,
        tmp_recursion_stack,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    kd_indices.setflags(write=False)

    return kd_indices


def build_spherical_kd(
        vectors,
        dimension_order
        ):
    """
    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (num_dimensions, num_vectors) in column-major layout (order='F')
        Input unit vectors
    dimension_order : np.ndarray[np.int32] of shape (num_dimensions,) in column-major layout (order='F')
        Dimension order

    Returns
    -------
    results : dict
        sphere_indices : np.ndarray[np.int32] of shape (num_vectors,) in column-major layout (order='F')
            Output index array,
        tmp_workspace : np.ndarray[np.int32] of shape (num_vectors,) in column-major layout (order='F')
            Workspace array,
        tmp_value_buffer : np.ndarray[np.float64] of shape (num_vectors,) in column-major layout (order='F')
            Value buffer,
        tmp_permutation : np.ndarray[np.int32] of shape (num_vectors,) in column-major layout (order='F')
            Permutation array,
        tmp_left_stack : np.ndarray[np.int32] of shape (num_vectors,) in column-major layout (order='F')
            Left stack,
        tmp_right_stack : np.ndarray[np.int32] of shape (num_vectors,) in column-major layout (order='F')
            Right stack,
        tmp_recursion_stack : np.ndarray[np.int32] of shape (3, num_vectors) in column-major layout (order='F')
            Stack for recursive calls


    Notes
    -----
    Build spherical k-d tree index
    """

    # ensure all array inputs are numpy arrays
    vectors = np.asfortranarray(vectors, dtype=np.float64)
    dimension_order = np.ascontiguousarray(dimension_order, dtype=np.int32)

    # extract dimension arguments
    num_dimensions = vectors.shape[0]
    num_vectors = vectors.shape[1]


    # Create temporaries and/or outputs
    sphere_indices = np.empty((num_vectors,), dtype=np.int32, order='F')
    tmp_workspace = np.empty((num_vectors,), dtype=np.int32, order='F')
    tmp_value_buffer = np.empty((num_vectors,), dtype=np.float64, order='F')
    tmp_permutation = np.empty((num_vectors,), dtype=np.int32, order='F')
    tmp_left_stack = np.empty((num_vectors,), dtype=np.int32, order='F')
    tmp_right_stack = np.empty((num_vectors,), dtype=np.int32, order='F')
    tmp_recursion_stack = np.empty((3, num_vectors), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.build_spherical_kd_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.build_spherical_kd_c.restype = None

    tox.build_spherical_kd_c(
        vectors,
        ctypes.byref(ctypes.c_int(num_dimensions)),
        ctypes.byref(ctypes.c_int(num_vectors)),
        sphere_indices,
        dimension_order,
        tmp_workspace,
        tmp_value_buffer,
        tmp_permutation,
        tmp_left_stack,
        tmp_right_stack,
        tmp_recursion_stack,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    sphere_indices.setflags(write=False)

    return sphere_indices


def get_kd_point(
        points,
        kd_indices,
        position,
        n_point_values_elements
        ):
    """
    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_points_elements_dim_1, n_points_elements_dim_2) in column-major layout (order='F')
        Input points
    kd_indices : np.ndarray[np.int32] of shape (n_kd_indices_elements,) in column-major layout (order='F')
        KD index array
    position : int
        Position in index

    Returns
    -------
    point_values : np.ndarray[np.float64] of shape (n_point_values_elements,) in column-major layout (order='F')
        Output point values

    Notes
    -----
    Get point from KD index
    """

    # ensure all array inputs are numpy arrays
    points = np.asfortranarray(points, dtype=np.float64)
    kd_indices = np.ascontiguousarray(kd_indices, dtype=np.int32)

    # extract dimension arguments
    n_points_elements_dim_1 = points.shape[0]
    n_points_elements_dim_2 = points.shape[1]
    n_kd_indices_elements = kd_indices.shape[0]


    # Create temporaries and/or outputs
    position = ctypes.c_int(position)
    point_values = np.empty((n_point_values_elements,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.get_kd_point_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.get_kd_point_c.restype = None

    tox.get_kd_point_c(
        points,
        ctypes.byref(ctypes.c_int(n_points_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_points_elements_dim_2)),
        kd_indices,
        ctypes.byref(ctypes.c_int(n_kd_indices_elements)),
        ctypes.byref(position),
        point_values,
        ctypes.byref(ctypes.c_int(n_point_values_elements)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    point_values.setflags(write=False)

    return point_values
