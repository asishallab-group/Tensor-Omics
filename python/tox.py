from error_handling import check_err_code

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
            Output array of matching indices
            The first `num_matches` elements will hold the results.,
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
        "output_indices": output_indices[..., :num_matches.value],
        "num_matches": num_matches.value
    }


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
    Build a k-d tree index using a stack-based, non-recursive approach.
    Initialize kd_indices to 1:num_points (original indices)
    Choose split dimension by cycling through dimension_order
    Find median index
    Partition kd_indices(left_idx:right_idx) by points(current_dim, kd_indices(:))
    Push right and left intervals onto stack
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


def which(
        mask,
        n,
        n_idx_out_elements,
        m_max
        ):
    """
    Parameters
    ----------
    mask : np.ndarray[np.int32] of shape (n_mask_elements,) in column-major layout (order='F')
        Logical array of size n.
    n : int
        Size of the mask.
    m_max : int
        Maximum size of idx_out.

    Returns
    -------
    results : dict
        idx_out : np.ndarray[np.int32] of shape (n_idx_out_elements,) in column-major layout (order='F')
            Integer array to store the indices of true values.,
        m_out : int
            Actual size of idx_out (number of true values found).

    Notes
    -----
    Finds the indices of the true values in a logical mask.
    """

    # ensure all array inputs are numpy arrays
    mask = np.ascontiguousarray(mask, dtype=np.int32)

    # extract dimension arguments
    n_mask_elements = mask.shape[0]


    # Create temporaries and/or outputs
    n = ctypes.c_int(n)
    idx_out = np.empty((n_idx_out_elements,), dtype=np.int32, order='F')
    m_max = ctypes.c_int(m_max)
    m_out = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.which_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.which_c.restype = None

    tox.which_c(
        mask,
        ctypes.byref(ctypes.c_int(n_mask_elements)),
        ctypes.byref(n),
        idx_out,
        ctypes.byref(ctypes.c_int(n_idx_out_elements)),
        ctypes.byref(m_max),
        ctypes.byref(m_out),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    idx_out.setflags(write=False)

    return {
        "idx_out": idx_out,
        "m_out": m_out.value
    }


def loess_smooth_2d(
        x_ref,
        y_ref,
        indices_used,
        x_query,
        kernel_sigma,
        kernel_cutoff
        ):
    """
    Parameters
    ----------
    x_ref : np.ndarray[np.float64] of shape (n_total,) in column-major layout (order='F')
        Reference x-coordinates.
    y_ref : np.ndarray[np.float64] of shape (n_total,) in column-major layout (order='F')
        Reference y-coordinates (length n_total).
    indices_used : np.ndarray[np.int32] of shape (n_used,) in column-major layout (order='F')
        Indices of reference points used for smoothing (only valid indices).
    x_query : np.ndarray[np.float64] of shape (n_target,) in column-major layout (order='F')
        Target x-coordinates to smooth.
    kernel_sigma : float
        Bandwidth parameter for the kernel.
    kernel_cutoff : float
        Cutoff for the kernel.

    Returns
    -------
    y_out : np.ndarray[np.float64] of shape (n_target,) in column-major layout (order='F')
        Output smoothed values (length n_target).

    Notes
    -----
    Performs LOESS smoothing on a set of data points.
    Smooths y_ref at x_query using reference points x_ref, y_ref, and kernel parameters.
    The user must pre-filter data and provide only valid indices in indices_used.
    """

    # ensure all array inputs are numpy arrays
    x_ref = np.ascontiguousarray(x_ref, dtype=np.float64)
    y_ref = np.ascontiguousarray(y_ref, dtype=np.float64)
    indices_used = np.ascontiguousarray(indices_used, dtype=np.int32)
    x_query = np.ascontiguousarray(x_query, dtype=np.float64)

    # extract dimension arguments
    n_total = x_ref.shape[0]
    n_target = x_query.shape[0]
    n_used = indices_used.shape[0]


    # Create temporaries and/or outputs
    kernel_sigma = ctypes.c_double(kernel_sigma)
    kernel_cutoff = ctypes.c_double(kernel_cutoff)
    y_out = np.empty((n_target,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.loess_smooth_2d_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.loess_smooth_2d_c.restype = None

    tox.loess_smooth_2d_c(
        ctypes.byref(ctypes.c_int(n_total)),
        ctypes.byref(ctypes.c_int(n_target)),
        x_ref,
        y_ref,
        indices_used,
        ctypes.byref(ctypes.c_int(n_used)),
        x_query,
        ctypes.byref(kernel_sigma),
        ctypes.byref(kernel_cutoff),
        y_out,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    y_out.setflags(write=False)

    return y_out


def compute_edf_expert(
        values,
        perm
        ):
    """
    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (n_values,) in column-major layout (order='F')
        Array of observed data values (e.g., contributions or spikes).
    perm : np.ndarray[np.int32] of shape (n_values,) in column-major layout (order='F')
        Pre-sorted permutation indices (must be sorted by values[perm]).

    Returns
    -------
    results : dict
        unique_values : np.ndarray[np.float64] of shape (n_values,) in column-major layout (order='F')
            Sorted unique data values.,
        cdf_values : np.ndarray[np.float64] of shape (n_values,) in column-major layout (order='F')
            Corresponding cumulative frequencies between 0 and 1.,
        n_unique : int
            Number of unique values found (actual size of output arrays)

    Notes
    -----
    Compute the Empirical Distribution Function (EDF) from pre-sorted permutation.
    Returns the sorted unique values and their cumulative frequencies in [0,1].
    Assumes perm is already sorted by values[perm]. Caller controls sorting algorithm.
    The number of unique values can be determined by finding the last non-zero cdf_value.
    """

    # ensure all array inputs are numpy arrays
    values = np.ascontiguousarray(values, dtype=np.float64)
    perm = np.ascontiguousarray(perm, dtype=np.int32)

    # extract dimension arguments
    n_values = values.shape[0]


    # Create temporaries and/or outputs
    unique_values = np.empty((n_values,), dtype=np.float64, order='F')
    cdf_values = np.empty((n_values,), dtype=np.float64, order='F')
    n_unique = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.compute_edf_expert_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.compute_edf_expert_c.restype = None

    tox.compute_edf_expert_c(
        values,
        ctypes.byref(ctypes.c_int(n_values)),
        perm,
        unique_values,
        cdf_values,
        ctypes.byref(n_unique),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    unique_values.setflags(write=False)
    cdf_values.setflags(write=False)

    return {
        "unique_values": unique_values,
        "cdf_values": cdf_values,
        "n_unique": n_unique.value
    }


def compute_edf(
        values
        ):
    """
    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (n_values,) in column-major layout (order='F')
        Array of observed data values (e.g., contributions or spikes).

    Returns
    -------
    results : dict
        unique_values : np.ndarray[np.float64] of shape (n_values,) in column-major layout (order='F')
            Sorted unique data values.,
        cdf_values : np.ndarray[np.float64] of shape (n_values,) in column-major layout (order='F')
            Corresponding cumulative frequencies between 0 and 1.,
        n_unique : int
            Number of unique values found (actual size of output arrays)

    Notes
    -----
    Helper routine that sorts and calls compute_edf.
    Allocates workspace internally and performs sorting before computing EDF.
    Use this for convenience; use compute_edf directly for custom sorting.
    """

    # ensure all array inputs are numpy arrays
    values = np.ascontiguousarray(values, dtype=np.float64)

    # extract dimension arguments
    n_values = values.shape[0]


    # Create temporaries and/or outputs
    unique_values = np.empty((n_values,), dtype=np.float64, order='F')
    cdf_values = np.empty((n_values,), dtype=np.float64, order='F')
    n_unique = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.compute_edf_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.compute_edf_c.restype = None

    tox.compute_edf_c(
        values,
        ctypes.byref(ctypes.c_int(n_values)),
        unique_values,
        cdf_values,
        ctypes.byref(n_unique),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    unique_values.setflags(write=False)
    cdf_values.setflags(write=False)

    return {
        "unique_values": unique_values,
        "cdf_values": cdf_values,
        "n_unique": n_unique.value
    }


def compute_empirical_p_values(
        rdi,
        sorted_rdi,
        perm,
        c_const
        ):
    """
    Parameters
    ----------
    rdi : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        Number of genes being processed.
    sorted_rdi : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        empirical distribution D
    perm : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        Constant used in the computation, typically 1
    c_const : float
        Output array to store the computed p-values for each gene.

    Returns
    -------
    p_values : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        empirical distribution D with non negative values

    Notes
    -----
    Calculate empirical p-values for scaled expression distances (RDI).
    Implements:
    P(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )
    Because distances are non-negative, a one-sided upper-tail empirical p-value is used.
    Assumptions / preconditions:
    - sorted_rdi(1:n_genes) contains the empirical distribution D.
    - If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
    """

    # ensure all array inputs are numpy arrays
    rdi = np.ascontiguousarray(rdi, dtype=np.float64)
    sorted_rdi = np.ascontiguousarray(sorted_rdi, dtype=np.float64)
    perm = np.ascontiguousarray(perm, dtype=np.int32)

    # extract dimension arguments
    n_genes = rdi.shape[0]


    # Create temporaries and/or outputs
    p_values = np.empty((n_genes,), dtype=np.float64, order='F')
    c_const = ctypes.c_double(c_const)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.compute_empirical_p_values_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.compute_empirical_p_values_c.restype = None

    tox.compute_empirical_p_values_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        rdi,
        sorted_rdi,
        perm,
        p_values,
        ctypes.byref(c_const),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    p_values.setflags(write=False)

    return p_values


def cluster_factor_trajectories_k_means(
        trajectories,
        centroids,
        max_iterations=300
        ):
    """
    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints) in column-major layout (order='F')
        matrix with data points to cluster
    centroids : np.ndarray[np.float64] of shape (n_factors, n_clusters) in column-major layout (order='F'), modified in-place
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
        The centroids should be unique. This is not checked in this routine.
        The final values will be the final centroids of the clusters
    max_iterations : int, optional
        number of maximum iterations of the clustering
        The default value is `300_int32`.

    Returns
    -------
    results : dict
        labels : np.ndarray[np.int32] of shape (n_samples * n_timepoints,) in column-major layout (order='F')
            array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            each label is the index of its related cluster -> `1<=label<=n_clusters=k`,
        label_counts : np.ndarray[np.int32] of shape (n_clusters,) in column-major layout (order='F')
            holds the number of points having the respective label assigned

    Notes
    -----
    Performs k-means clustering on factor trajectories, so factor evolution over time
    """

    # ensure all array inputs are numpy arrays
    trajectories = np.asfortranarray(trajectories, dtype=np.float64)
    centroids = np.asfortranarray(centroids, dtype=np.float64)

    # extract dimension arguments
    n_clusters = centroids.shape[1]
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]


    # Create temporaries and/or outputs
    labels = np.empty((n_samples * n_timepoints,), dtype=np.int32, order='F')
    label_counts = np.empty((n_clusters,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)
    max_iterations = ctypes.c_int(max_iterations)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.cluster_factor_trajectories_k_means_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=3, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.cluster_factor_trajectories_k_means_c.restype = None

    tox.cluster_factor_trajectories_k_means_c(
        ctypes.byref(ctypes.c_int(n_clusters)),
        trajectories,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        centroids,
        labels,
        label_counts,
        ctypes.byref(ierr),
        ctypes.byref(max_iterations)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    centroids.setflags(write=False)
    labels.setflags(write=False)
    label_counts.setflags(write=False)

    return {
        "labels": labels,
        "label_counts": label_counts
    }


def k_means_clustering(
        data_points,
        centroids,
        max_iterations=300
        ):
    """
    performs k-means clustering on `data_points`


    Parameters
    ----------
    data_points : np.ndarray[np.float64] of shape (n_dims, n_points) in column-major layout (order='F')
        matrix with data points to cluster
    centroids : np.ndarray[np.float64] of shape (n_dims, n_clusters) in column-major layout (order='F'), modified in-place
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
        The centroids should be unique. This is not checked in this routine.
        The final values will be the final centroids of the clusters
    max_iterations : int, optional
        number of maximum iterations of the clustering
        The default value is `300_int32`.

    Returns
    -------
    results : dict
        labels : np.ndarray[np.int32] of shape (n_points,) in column-major layout (order='F')
            array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            each label is the index of its related cluster -> `1<=label<=n_clusters=k`,
        label_counts : np.ndarray[np.int32] of shape (n_clusters,) in column-major layout (order='F')
            holds the number of points having the respective label assigned

    Notes
    -----
    k-means clustering algorithm:
    1. Assigns each data point to one of `k` clusters whose centroid is clostest
    2. Recalculates the centroids using the mean of its assigned points
    3. repeat 1-2 until assignment remains unchanged
    """

    # ensure all array inputs are numpy arrays
    data_points = np.asfortranarray(data_points, dtype=np.float64)
    centroids = np.asfortranarray(centroids, dtype=np.float64)

    # extract dimension arguments
    n_clusters = centroids.shape[1]
    n_points = data_points.shape[1]
    n_dims = data_points.shape[0]


    # Create temporaries and/or outputs
    labels = np.empty((n_points,), dtype=np.int32, order='F')
    label_counts = np.empty((n_clusters,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)
    max_iterations = ctypes.c_int(max_iterations)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.k_means_clustering_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.k_means_clustering_c.restype = None

    tox.k_means_clustering_c(
        ctypes.byref(ctypes.c_int(n_clusters)),
        data_points,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_dims)),
        centroids,
        labels,
        label_counts,
        ctypes.byref(ierr),
        ctypes.byref(max_iterations)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    centroids.setflags(write=False)
    labels.setflags(write=False)
    label_counts.setflags(write=False)

    return {
        "labels": labels,
        "label_counts": label_counts
    }


def linkage_clustering(
        distances,
        method
        ):
    """
    Perform linkage clustering on a distance matrix.


    Parameters
    ----------
    distances : np.ndarray[np.float64] of shape (n_points, n_points) in column-major layout (order='F'), modified in-place
        symmetric distance matrix, holding the positive distances between points. Distance of X->X is always zero.
        @note
        This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
        So there is no need to copy an existing distance matrix, just pass the original.
        @endnote
    method : str
        used algorithm
        |      Method      |                      Value                     |
        |------------------|------------------------------------------------|
        | Average / UPGMA  |    "average"      |
        | Weighted / WPGMA |    "weighted"     |
        |      Ward        |    "ward"         |

    Returns
    -------
    results : dict
        merge_i : np.ndarray[np.int32] of shape (n_points - 1,) in column-major layout (order='F')
            holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes,
        merge_j : np.ndarray[np.int32] of shape (n_points - 1,) in column-major layout (order='F')
            holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes,
        heights : np.ndarray[np.float64] of shape (n_points - 1,) in column-major layout (order='F')
            height of the shorter branch of the merge, e.g. if (A,B)+(C) merges to ((A,B),C), the branch to (A,B) is shorter,
        cluster_sizes : np.ndarray[np.int32] of shape (n_points - 1,) in column-major layout (order='F')
            size of cluster at iteration k

    Notes
    -----
    @note
    This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
    So there is no need to copy an existing distance matrix, just pass the original.
    @endnote
    """

    # ensure all array inputs are numpy arrays
    distances = np.asfortranarray(distances, dtype=np.float64)
    method = np.asarray(method)

    # extract dimension arguments
    n_points = distances.shape[0]


    # Create temporaries and/or outputs
    merge_i = np.empty((n_points - 1,), dtype=np.int32, order='F')
    merge_j = np.empty((n_points - 1,), dtype=np.int32, order='F')
    heights = np.empty((n_points - 1,), dtype=np.float64, order='F')
    cluster_sizes = np.empty((n_points - 1,), dtype=np.int32, order='F')
    method = method.astype(f"S{8}", order="F")
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.linkage_clustering_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{8}"),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.linkage_clustering_c.restype = None

    tox.linkage_clustering_c(
        distances,
        ctypes.byref(ctypes.c_int(n_points)),
        merge_i,
        merge_j,
        heights,
        cluster_sizes,
        method,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    distances.setflags(write=False)
    merge_i.setflags(write=False)
    merge_j.setflags(write=False)
    heights.setflags(write=False)
    cluster_sizes.setflags(write=False)

    return {
        "merge_i": merge_i,
        "merge_j": merge_j,
        "heights": heights,
        "cluster_sizes": cluster_sizes
    }


def mean_vector(
        expression_vectors,
        gene_indices
        ):
    """
    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes) in column-major layout (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_indices : np.ndarray[np.int32] of shape (n_selected_genes,) in column-major layout (order='F')
        An array containing the column indices of the selected genes in 'expression_vectors'.

    Returns
    -------
    centroid : np.ndarray[np.float64] of shape (n_axes,) in column-major layout (order='F')
        The output vector representing the computed centroid.

    Notes
    -----
    Computes the element-wise mean for a given set of vectors.
    """

    # ensure all array inputs are numpy arrays
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    gene_indices = np.ascontiguousarray(gene_indices, dtype=np.int32)

    # extract dimension arguments
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]
    n_selected_genes = gene_indices.shape[0]


    # Create temporaries and/or outputs
    centroid = np.empty((n_axes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.mean_vector_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.mean_vector_c.restype = None

    tox.mean_vector_c(
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        gene_indices,
        ctypes.byref(ctypes.c_int(n_selected_genes)),
        centroid,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    centroid.setflags(write=False)

    return centroid


def group_centroid(
        expression_vectors,
        gene_to_family,
        n_families,
        mode,
        ortholog_set=None
        ):
    """
    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes) in column-major layout (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_to_family : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        An array mapping each gene (by index) to a family ID.
    mode : str
        used mode for grouping
        |       Mode       |                             Value                               |
        |------------------|-----------------------------------------------------------------|
        | Group Orthologs  |   "group_orthologs"    |
        |    Group all     |      "group_all"       |
    ortholog_set : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F'), optional
        A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
        This optional argument needs to be passed if used mode (`mode`) is [[tox_gene_centroids(module):MODE_GROUP_ORTHOLOGS(variable)]].

    Returns
    -------
    results : dict
        centroid_matrix : np.ndarray[np.float64] of shape (n_axes, n_families) in column-major layout (order='F')
            The output matrix (n_axes x n_families) to store the computed centroids.,
        tmp_selected_indices : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
            An output array for storing indices.

    Notes
    -----
    Iterates over families, filters gene indices, and computes centroids.
    """

    # ensure all array inputs are numpy arrays
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    gene_to_family = np.ascontiguousarray(gene_to_family, dtype=np.int32)
    mode = np.asarray(mode)
    if ortholog_set is not None:
        ortholog_set = np.ascontiguousarray(ortholog_set, dtype=np.int32)

    # extract dimension arguments
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]


    # Create temporaries and/or outputs
    centroid_matrix = np.empty((n_axes, n_families), dtype=np.float64, order='F')
    mode = mode.astype(f"S{15}", order="F")
    tmp_selected_indices = np.empty((n_genes,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.group_centroid_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{15}"),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        nullable(np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32))
    )
    tox.group_centroid_c.restype = None

    tox.group_centroid_c(
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        gene_to_family,
        ctypes.byref(ctypes.c_int(n_families)),
        centroid_matrix,
        mode,
        tmp_selected_indices,
        ctypes.byref(ierr),
        ortholog_set
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    centroid_matrix.setflags(write=False)

    return centroid_matrix


def tox_loess_required_workspace(
        d,
        nvmax,
        setlf
        ):
    """
    Parameters
    ----------
    d : int
        Dimensionality of the data
    nvmax : int
        Maximum neighborhood size
    setlf : bool
        Save matrix factorization flag

    Returns
    -------
    results : dict
        int_workspace_size : int
            Required size of the integer workspace array,
        real_workspace_size : int
            Required size of the real workspace array

    Notes
    -----
    Recommend workspace sizes based on Netlib exact formulas.
    Computes the required sizes for integer and real workspace arrays.
    These sizes depend on the dimensionality of the data and the maximum neighborhood size.
    """

    # ensure all array inputs are numpy arrays


    # extract dimension arguments



    # Create temporaries and/or outputs
    d = ctypes.c_int(d)
    nvmax = ctypes.c_int(nvmax)
    int_workspace_size = ctypes.c_int(0)
    real_workspace_size = ctypes.c_int(0)
    setlf = ctypes.c_int(setlf)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.tox_loess_required_workspace_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.tox_loess_required_workspace_c.restype = None

    tox.tox_loess_required_workspace_c(
        ctypes.byref(d),
        ctypes.byref(nvmax),
        ctypes.byref(int_workspace_size),
        ctypes.byref(real_workspace_size),
        ctypes.byref(setlf),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return {
        "int_workspace_size": int_workspace_size.value,
        "real_workspace_size": real_workspace_size.value
    }


def loess_fit_plain(
        x,
        y,
        w,
        eval_points,
        span,
        degree,
        nvmax,
        infl,
        setlf,
        int_workspace,
        real_workspace,
        diagl
        ):
    """
    Parameters
    ----------
    x : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F')
        Predictor variable array
    y : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F')
        Response variable array
    w : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F')
        Weight array for data points
    eval_points : np.ndarray[np.float64] of shape (n, 1) in column-major layout (order='F')
        Evaluation points (x values at which the fitted curve is computed)
    span : float
        Smoothing parameter for LOESS
    degree : int
        Degree of the LOESS polynomial
    nvmax : int
        Maximum neighborhood size
    infl : bool
        Influence calculation flag
    setlf : bool
        Save matrix factorization flag
    int_workspace : np.ndarray[np.int32] of shape (int_workspace_size,) in column-major layout (order='F'), modified in-place
        Integer workspace array
    real_workspace : np.ndarray[np.float64] of shape (real_workspace_size,) in column-major layout (order='F'), modified in-place
        Real workspace array
    diagl : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F'), modified in-place
        Diagonal elements of the hat matrix

    Returns
    -------
    fitted_values : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F')
        Fitted (smoothed) values of y at the evaluation points

    Notes
    -----
    Perform plain LOESS fitting.
    Fits a LOESS model to the data using the specified smoothing parameter.
    Outputs the smoothed response variable array.
    """

    # ensure all array inputs are numpy arrays
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    w = np.ascontiguousarray(w, dtype=np.float64)
    eval_points = np.asfortranarray(eval_points, dtype=np.float64)
    int_workspace = np.ascontiguousarray(int_workspace, dtype=np.int32)
    real_workspace = np.ascontiguousarray(real_workspace, dtype=np.float64)
    diagl = np.ascontiguousarray(diagl, dtype=np.float64)

    # extract dimension arguments
    n = x.shape[0]
    int_workspace_size = int_workspace.shape[0]
    real_workspace_size = real_workspace.shape[0]


    # Create temporaries and/or outputs
    span = ctypes.c_double(span)
    degree = ctypes.c_int(degree)
    nvmax = ctypes.c_int(nvmax)
    infl = ctypes.c_int(infl)
    setlf = ctypes.c_int(setlf)
    fitted_values = np.empty((n,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.loess_fit_plain_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.loess_fit_plain_c.restype = None

    tox.loess_fit_plain_c(
        ctypes.byref(ctypes.c_int(n)),
        x,
        y,
        w,
        eval_points,
        ctypes.byref(span),
        ctypes.byref(degree),
        ctypes.byref(nvmax),
        ctypes.byref(infl),
        ctypes.byref(setlf),
        int_workspace,
        ctypes.byref(ctypes.c_int(int_workspace_size)),
        real_workspace,
        ctypes.byref(ctypes.c_int(real_workspace_size)),
        diagl,
        fitted_values,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    int_workspace.setflags(write=False)
    real_workspace.setflags(write=False)
    diagl.setflags(write=False)
    fitted_values.setflags(write=False)

    return fitted_values


def loess_fit_robust(
        x,
        y,
        w,
        eval_points,
        span,
        degree,
        nvmax,
        infl,
        setlf,
        n_iters,
        int_workspace,
        real_workspace,
        diagl,
        robust_weights,
        combined_weights,
        residuals,
        permutation_indices
        ):
    """
    Parameters
    ----------
    x : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F')
        Predictor variable array
    y : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F')
        Response variable array
    w : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F')
        Weight array for data points
    eval_points : np.ndarray[np.float64] of shape (n, 1) in column-major layout (order='F')
        Evaluation points (x values at which the fitted curve is computed)
    span : float
        Smoothing parameter for LOESS
    degree : int
        Degree of the LOESS polynomial
    nvmax : int
        Maximum neighborhood size
    infl : bool
        Influence calculation flag
    setlf : bool
        Save matrix factorization flag
    n_iters : int
        Number of robust iterations
    int_workspace : np.ndarray[np.int32] of shape (int_workspace_size,) in column-major layout (order='F'), modified in-place
        Integer workspace array
    real_workspace : np.ndarray[np.float64] of shape (real_workspace_size,) in column-major layout (order='F'), modified in-place
        Real workspace array
    diagl : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F'), modified in-place
        Diagonal elements of the hat matrix
    robust_weights : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F'), modified in-place
        Robust bisquare weights (updated each iteration, initialized to 1.0)
    combined_weights : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F'), modified in-place
        Combined weights: product of user weights and robust weights (w(i) * robust_weights(i))
    residuals : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F'), modified in-place
        Residuals (y - fitted_values), used to compute bisquare robust weights
    permutation_indices : np.ndarray[np.int32] of shape (n,) in column-major layout (order='F'), modified in-place
        Permutation indices array (from NetLib bisquare weight computation)

    Returns
    -------
    fitted_values : np.ndarray[np.float64] of shape (n,) in column-major layout (order='F')
        Fitted (smoothed) values of y at the evaluation points

    Notes
    -----
    Perform robust LOESS fitting with bisquare reweighting.
    Fits a LOESS model to the data using robust iterations to handle outliers.
    The robust fitting process iterates n_iters times, each iteration:
    - Combines original weights with robust weights (down-weights from previous iteration)
    - Runs LOESS fitting with combined weights
    - Computes residuals (y - fitted values)
    - Updates robust weights using bisquare function (suppresses large residuals)
    """

    # ensure all array inputs are numpy arrays
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    w = np.ascontiguousarray(w, dtype=np.float64)
    eval_points = np.asfortranarray(eval_points, dtype=np.float64)
    int_workspace = np.ascontiguousarray(int_workspace, dtype=np.int32)
    real_workspace = np.ascontiguousarray(real_workspace, dtype=np.float64)
    diagl = np.ascontiguousarray(diagl, dtype=np.float64)
    robust_weights = np.ascontiguousarray(robust_weights, dtype=np.float64)
    combined_weights = np.ascontiguousarray(combined_weights, dtype=np.float64)
    residuals = np.ascontiguousarray(residuals, dtype=np.float64)
    permutation_indices = np.ascontiguousarray(permutation_indices, dtype=np.int32)

    # extract dimension arguments
    n = x.shape[0]
    int_workspace_size = int_workspace.shape[0]
    real_workspace_size = real_workspace.shape[0]


    # Create temporaries and/or outputs
    span = ctypes.c_double(span)
    degree = ctypes.c_int(degree)
    nvmax = ctypes.c_int(nvmax)
    infl = ctypes.c_int(infl)
    setlf = ctypes.c_int(setlf)
    n_iters = ctypes.c_int(n_iters)
    fitted_values = np.empty((n,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.loess_fit_robust_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.loess_fit_robust_c.restype = None

    tox.loess_fit_robust_c(
        ctypes.byref(ctypes.c_int(n)),
        x,
        y,
        w,
        eval_points,
        ctypes.byref(span),
        ctypes.byref(degree),
        ctypes.byref(nvmax),
        ctypes.byref(infl),
        ctypes.byref(setlf),
        ctypes.byref(n_iters),
        int_workspace,
        ctypes.byref(ctypes.c_int(int_workspace_size)),
        real_workspace,
        ctypes.byref(ctypes.c_int(real_workspace_size)),
        diagl,
        robust_weights,
        combined_weights,
        residuals,
        permutation_indices,
        fitted_values,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    int_workspace.setflags(write=False)
    real_workspace.setflags(write=False)
    diagl.setflags(write=False)
    robust_weights.setflags(write=False)
    combined_weights.setflags(write=False)
    residuals.setflags(write=False)
    permutation_indices.setflags(write=False)
    fitted_values.setflags(write=False)

    return fitted_values


def loess(
        x,
        y,
        span,
        degree,
        mode,
        n_iters
        ):
    """
    Parameters
    ----------
    x : np.ndarray[np.float64] of shape (n_x_elements,) in column-major layout (order='F')
        Predictor variable array
    y : np.ndarray[np.float64] of shape (n_y_elements,) in column-major layout (order='F')
        Response variable array
    span : float
        Smoothing parameter for LOESS
    degree : int
        Degree of the LOESS polynomial
    mode : str
        Mode of operation
        |    Mode   |   Value  |
        |-----------|----------|
        | robust fitting |  "robust"   |
        | plain fitting  |  "plain"    |
    n_iters : int
        Number of robust iterations (only used when mode = 1)

    Returns
    -------
    fitted_values : np.ndarray[np.float64] of shape (size(y),) in column-major layout (order='F')
        Fitted (smoothed) values of y

    Notes
    -----
    Wrapper subroutine for LOESS fitting (plain or robust).
    This subroutine selects between plain and robust LOESS fitting based on the mode.
    It dynamically allocates the required arrays and computes workspace sizes.
    """

    # ensure all array inputs are numpy arrays
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    mode = np.asarray(mode)

    # extract dimension arguments
    n_x_elements = x.shape[0]
    n_y_elements = y.shape[0]


    # Create temporaries and/or outputs
    span = ctypes.c_double(span)
    degree = ctypes.c_int(degree)
    fitted_values = np.empty((size(y),), dtype=np.float64, order='F')
    mode = mode.astype(f"S{6}", order="F")
    n_iters = ctypes.c_int(n_iters)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.loess_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{6}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.loess_c.restype = None

    tox.loess_c(
        x,
        ctypes.byref(ctypes.c_int(n_x_elements)),
        y,
        ctypes.byref(ctypes.c_int(n_y_elements)),
        ctypes.byref(span),
        ctypes.byref(degree),
        fitted_values,
        mode,
        ctypes.byref(n_iters),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    fitted_values.setflags(write=False)

    return fitted_values


def normalize_unit_length(
        vector
        ):
    """
    Parameters
    ----------
    vector : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F'), modified in-place
        Vector that will be normalized to unit length

    Returns
    -------
    None

    Notes
    -----
    Normalizes an input vector to unit length in-place
    """

    # ensure all array inputs are numpy arrays
    vector = np.ascontiguousarray(vector, dtype=np.float64)

    # extract dimension arguments
    n_dims = vector.shape[0]


    # Create temporaries and/or outputs
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.normalize_unit_length_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.normalize_unit_length_c.restype = None

    tox.normalize_unit_length_c(
        vector,
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    vector.setflags(write=False)

    return None


def mask_get_first_successor_idx(
        bit_mask
        ):
    """
    Parameters
    ----------
    bit_mask : np.ndarray[np.int32] of shape (n_bit_mask_elements,) in column-major layout (order='F')
        chunked mask to mark active genes

    Returns
    -------
    idx : int
        index of last active gene

    Notes
    -----
    Helper function that returns the index after the last active gene in `bit_mask`, so the first succeeding gene.
    """

    # ensure all array inputs are numpy arrays
    bit_mask = np.ascontiguousarray(bit_mask, dtype=np.int32)

    # extract dimension arguments
    n_bit_mask_elements = bit_mask.shape[0]


    # Create temporaries and/or outputs
    idx = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.mask_get_first_successor_idx_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.mask_get_first_successor_idx_c.restype = None

    tox.mask_get_first_successor_idx_c(
        bit_mask,
        ctypes.byref(ctypes.c_int(n_bit_mask_elements)),
        ctypes.byref(idx),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return idx.value


def mask_check_state(
        bit_mask,
        i_gene
        ):
    """
    Parameters
    ----------
    bit_mask : np.ndarray[np.int32] of shape (n_bit_mask_elements,) in column-major layout (order='F')
        chunked mask to mark active paralogs
    i_gene : int
        index of paralog to be marked active

    Returns
    -------
    state : bool
        check result

    Notes
    -----
    Checks the state of a bit/paralog in `bit_mask` -> .true. if 1 else .false.
    """

    # ensure all array inputs are numpy arrays
    bit_mask = np.ascontiguousarray(bit_mask, dtype=np.int32)

    # extract dimension arguments
    n_bit_mask_elements = bit_mask.shape[0]


    # Create temporaries and/or outputs
    i_gene = ctypes.c_int(i_gene)
    state = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.mask_check_state_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.mask_check_state_c.restype = None

    tox.mask_check_state_c(
        bit_mask,
        ctypes.byref(ctypes.c_int(n_bit_mask_elements)),
        ctypes.byref(i_gene),
        ctypes.byref(state),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return state.value


def detect_neofunctionalization(
        ancestors,
        genes,
        gene_to_fam,
        thresholds
        ):
    """
    Parameters
    ----------
    ancestors : np.ndarray[np.float64] of shape (n_axes, n_families) in column-major layout (order='F')
        RAP projected unit length expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_axes, n_genes) in column-major layout (order='F')
        RAP projected unit length expression vectors of genes
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        mapping of gene index to family index
    thresholds : np.ndarray[np.float64] of shape (n_axes,) in column-major layout (order='F')
        threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis

    Returns
    -------
    neofunc : np.ndarray[np.int32] of shape (n_genes, n_axes) in column-major layout (order='F')
        `.true.` if neofunctionalization has been detected for the respective axes

    Notes
    -----
    Identifies neofunctionalization for genes by checking whether the difference of expression to its ancestor exceeds the threshold for the respective axis.
    """

    # ensure all array inputs are numpy arrays
    ancestors = np.asfortranarray(ancestors, dtype=np.float64)
    genes = np.asfortranarray(genes, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    thresholds = np.ascontiguousarray(thresholds, dtype=np.float64)

    # extract dimension arguments
    n_families = ancestors.shape[1]
    n_axes = ancestors.shape[0]
    n_genes = genes.shape[1]


    # Create temporaries and/or outputs
    neofunc = np.empty((n_genes, n_axes), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.detect_neofunctionalization_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.detect_neofunctionalization_c.restype = None

    tox.detect_neofunctionalization_c(
        ancestors,
        ctypes.byref(ctypes.c_int(n_families)),
        genes,
        ctypes.byref(ctypes.c_int(n_axes)),
        gene_to_fam,
        ctypes.byref(ctypes.c_int(n_genes)),
        thresholds,
        neofunc,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    neofunc.setflags(write=False)

    return neofunc


def detect_dosage_effect(
        ancestor,
        genes,
        filtered_paralogs_mask,
        max_subset_size,
        n_paralog_subsets,
        max_angle=3.141592653589793,
        gain_gamma=0.1
        ):
    """
    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes) in column-major layout (order='F')
        expression vectors of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
        bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
        M_DM_FROM_JUST_INFO to compute this argument using [[tox_paralog_analysis(module):filter_paralogs_by_pattern]]'s output `masks(:, family_idx)`.
    max_subset_size : int
        maximum subset size of checked gene subsets.
        M_DM_FROM_AUTO to compute this argument using [[tox_paralog_analysis(module):calc_work_arr_paralog_subsets_size]]'s output `max_subset_size`.
    max_angle : float, optional
        in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
        The default value is `pi`.
    gain_gamma : float, optional
        positive magnitude gain for dosage effect
        The default value is `0.1_real64`.

    Returns
    -------
    results : dict
        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`,
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_results) in column-major layout (order='F')
            working array to hold bitmask encoded subsets for detection.
            @note
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
            @endnote
            The first `n_results` elements will hold the results.,
        active_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
            working array to hold the extended subsets,
        temp_paralog_vector : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
            vector used for pruning subsets

    Notes
    -----
    Identifies subsets of paralogs with small angle to the `ancestor` (max_angle) and sum to a magnitude significantly exceeding `norm(ancestor)` (gain)
    """

    # ensure all array inputs are numpy arrays
    ancestor = np.ascontiguousarray(ancestor, dtype=np.float64)
    genes = np.asfortranarray(genes, dtype=np.float64)
    filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)

    # extract dimension arguments
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]


    # Create temporaries and/or outputs
    n_results = ctypes.c_int(0)
    max_subset_size = ctypes.c_int(max_subset_size)
    work_arr_paralog_subsets = np.empty((n_mask_chunks, n_paralog_subsets), dtype=np.int32, order='F')
    active_mask = np.empty((n_mask_chunks,), dtype=np.int32, order='F')
    temp_paralog_vector = np.empty((n_dims,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)
    max_angle = ctypes.c_double(max_angle)
    gain_gamma = ctypes.c_double(gain_gamma)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.detect_dosage_effect_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double)
    )
    tox.detect_dosage_effect_c.restype = None

    tox.detect_dosage_effect_c(
        ancestor,
        genes,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_dims)),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(n_results),
        ctypes.byref(max_subset_size),
        work_arr_paralog_subsets,
        ctypes.byref(ctypes.c_int(n_paralog_subsets)),
        active_mask,
        temp_paralog_vector,
        ctypes.byref(ierr),
        ctypes.byref(max_angle),
        ctypes.byref(gain_gamma)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    work_arr_paralog_subsets.setflags(write=False)
    active_mask.setflags(write=False)
    temp_paralog_vector.setflags(write=False)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets[..., :n_results.value],
        "active_mask": active_mask,
        "temp_paralog_vector": temp_paralog_vector
    }


def detect_subfunctionalization(
        ancestor,
        genes,
        rdi_threshold,
        filtered_paralogs_mask,
        max_subset_size,
        n_paralog_subsets,
        paralog_norms,
        sorted_paralog_norms_perm
        ):
    """
    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes) in column-major layout (order='F')
        expression vectors of genes
    rdi_threshold : float
        max allowed residual distance from `ancestor`
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
        bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    max_subset_size : int
        maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    paralog_norms : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` function from `f42_utils` function for this)
    sorted_paralog_norms_perm : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor

    Returns
    -------
    results : dict
        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`,
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_results) in column-major layout (order='F')
            working array to hold bitmask encoded subsets for detection.
            @note
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
            @endnote
            The first `n_results` elements will hold the results.,
        active_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
            working array to hold the extended subsets,
        temp_paralog_vector : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
            vector used for pruning subsets,
        temp_work_array : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
            needed for efficient check of minimum value after a certain index

    Notes
    -----
    Identifies subsets of paralogs exhibiting significant angles to the `ancestor`
    """

    # ensure all array inputs are numpy arrays
    ancestor = np.ascontiguousarray(ancestor, dtype=np.float64)
    genes = np.asfortranarray(genes, dtype=np.float64)
    filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)
    paralog_norms = np.ascontiguousarray(paralog_norms, dtype=np.float64)
    sorted_paralog_norms_perm = np.ascontiguousarray(sorted_paralog_norms_perm, dtype=np.int32)

    # extract dimension arguments
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]


    # Create temporaries and/or outputs
    rdi_threshold = ctypes.c_double(rdi_threshold)
    n_results = ctypes.c_int(0)
    max_subset_size = ctypes.c_int(max_subset_size)
    work_arr_paralog_subsets = np.empty((n_mask_chunks, n_paralog_subsets), dtype=np.int32, order='F')
    active_mask = np.empty((n_mask_chunks,), dtype=np.int32, order='F')
    temp_paralog_vector = np.empty((n_dims,), dtype=np.float64, order='F')
    temp_work_array = np.empty((n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.detect_subfunctionalization_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.detect_subfunctionalization_c.restype = None

    tox.detect_subfunctionalization_c(
        ancestor,
        genes,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(rdi_threshold),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(n_results),
        ctypes.byref(max_subset_size),
        work_arr_paralog_subsets,
        ctypes.byref(ctypes.c_int(n_paralog_subsets)),
        active_mask,
        temp_paralog_vector,
        paralog_norms,
        sorted_paralog_norms_perm,
        temp_work_array,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    work_arr_paralog_subsets.setflags(write=False)
    active_mask.setflags(write=False)
    temp_paralog_vector.setflags(write=False)
    temp_work_array.setflags(write=False)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets[..., :n_results.value],
        "active_mask": active_mask,
        "temp_paralog_vector": temp_paralog_vector,
        "temp_work_array": temp_work_array
    }


def detect_patterns(
        ancestor,
        genes,
        pattern_mode,
        filtered_paralogs_mask,
        max_subset_size,
        n_paralog_subsets,
        dosage_max_angle=3.141592653589793,
        dosage_gain_gamma=0.1,
        subfunc_rdi_threshold=None,
        subfunc_paralog_norms=None,
        subfunc_sorted_paralog_norms_perm=None
        ):
    """
    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes) in column-major layout (order='F')
        expression vectors of genes
    pattern_mode : str
        used pattern for detection
        |       Mode           |                               Value                                |
        |----------------------|--------------------------------------------------------------------|
        |    Dosage Effect     |    "dosage_pattern"     |
        | Subfunctionalization |    "subfunc_pattern"    |
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
        bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    max_subset_size : int
        maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    dosage_max_angle : float, optional
        in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
        The default value is `pi`.
    dosage_gain_gamma : float, optional
        in dosage mode required positive magnitude gain for dosage
        The default value is `0.1_real64`.
    subfunc_rdi_threshold : float, optional
        max allowed residual distance from `ancestor`
        This optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].
    subfunc_paralog_norms : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F'), optional
        in subfunctionalization mode needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` from `f42_utils` function for this)
        This optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].
    subfunc_sorted_paralog_norms_perm : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F'), optional
        in subfunctionalization mode needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor
        This optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].

    Returns
    -------
    results : dict
        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`,
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_results) in column-major layout (order='F')
            working array to hold bitmask encoded subsets for detection.
            @note
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
            @endnote
            The first `n_results` elements will hold the results.,
        active_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
            working array to hold the extended subsets,
        temp_paralog_vector : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
            vector used for pruning subsets,
        subfunc_temp_work_array : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
            in subfunctionalization mode needed for efficient check of minimum value after a certain index
            This optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].

    Notes
    -----
    Identifies subsets of paralogs where dosage effect or subfunctionalization applies, depending on `pattern`
    """

    # ensure all array inputs are numpy arrays
    ancestor = np.ascontiguousarray(ancestor, dtype=np.float64)
    genes = np.asfortranarray(genes, dtype=np.float64)
    pattern_mode = np.asarray(pattern_mode)
    filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)
    if subfunc_paralog_norms is not None:
        subfunc_paralog_norms = np.ascontiguousarray(subfunc_paralog_norms, dtype=np.float64)
    if subfunc_sorted_paralog_norms_perm is not None:
        subfunc_sorted_paralog_norms_perm = np.ascontiguousarray(subfunc_sorted_paralog_norms_perm, dtype=np.int32)

    # extract dimension arguments
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]


    # Create temporaries and/or outputs
    pattern_mode = pattern_mode.astype(f"S{15}", order="F")
    n_results = ctypes.c_int(0)
    max_subset_size = ctypes.c_int(max_subset_size)
    work_arr_paralog_subsets = np.empty((n_mask_chunks, n_paralog_subsets), dtype=np.int32, order='F')
    active_mask = np.empty((n_mask_chunks,), dtype=np.int32, order='F')
    temp_paralog_vector = np.empty((n_dims,), dtype=np.float64, order='F')
    dosage_max_angle = ctypes.c_double(dosage_max_angle)
    dosage_gain_gamma = ctypes.c_double(dosage_gain_gamma)
    subfunc_rdi_threshold = ctypes.c_double(subfunc_rdi_threshold)
    subfunc_temp_work_array = np.empty((n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.detect_patterns_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{15}"),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        nullable(ctypes.POINTER(ctypes.c_double)),
        nullable(np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64)),
        nullable(np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32)),
        nullable(np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64)),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.detect_patterns_c.restype = None

    tox.detect_patterns_c(
        ancestor,
        genes,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_dims)),
        pattern_mode,
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(n_results),
        ctypes.byref(max_subset_size),
        work_arr_paralog_subsets,
        ctypes.byref(ctypes.c_int(n_paralog_subsets)),
        active_mask,
        temp_paralog_vector,
        ctypes.byref(dosage_max_angle),
        ctypes.byref(dosage_gain_gamma),
        ctypes.byref(subfunc_rdi_threshold),
        subfunc_paralog_norms,
        subfunc_sorted_paralog_norms_perm,
        subfunc_temp_work_array,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    work_arr_paralog_subsets.setflags(write=False)
    active_mask.setflags(write=False)
    temp_paralog_vector.setflags(write=False)
    subfunc_temp_work_array.setflags(write=False)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets[..., :n_results.value],
        "active_mask": active_mask,
        "temp_paralog_vector": temp_paralog_vector,
        "subfunc_temp_work_array": subfunc_temp_work_array
    }


def mask_chunk_count(
        n_genes
        ):
    """
    Parameters
    ----------
    n_genes : int
        number of genes

    Returns
    -------
    count : int
        number of 32 bit chunks a mask needs to encode `n_genes` genes
        Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks

    Notes
    -----
    This subroutine easily determines the needed chunk count for subset bit masks, as an integer has only 32 bits.
    """

    # ensure all array inputs are numpy arrays


    # extract dimension arguments



    # Create temporaries and/or outputs
    n_genes = ctypes.c_int(n_genes)
    count = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.mask_chunk_count_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.mask_chunk_count_c.restype = None

    tox.mask_chunk_count_c(
        ctypes.byref(n_genes),
        ctypes.byref(count),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return count.value


def filter_paralogs_by_pattern_subfunctionalization(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks
        ):
    """
    Parameters
    ----------
    gene_angles : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
    threshold : float
        filter threshold
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families) in column-major layout (order='F')
        bit mask that will have indices of genes kept by pattern set to 1, else 0

    Notes
    -----
    This subroutine prefilters the genes for subfunctionalization,
    as genes that are already too close in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
    """

    # ensure all array inputs are numpy arrays
    gene_angles = np.ascontiguousarray(gene_angles, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)

    # extract dimension arguments
    n_genes = gene_angles.shape[0]


    # Create temporaries and/or outputs
    threshold = ctypes.c_double(threshold)
    masks = np.empty((n_mask_chunks, n_families), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.filter_paralogs_by_pattern_subfunctionalization_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.filter_paralogs_by_pattern_subfunctionalization_c.restype = None

    tox.filter_paralogs_by_pattern_subfunctionalization_c(
        gene_angles,
        ctypes.byref(threshold),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        gene_to_fam,
        masks,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    masks.setflags(write=False)

    return masks


def filter_paralogs_by_pattern_dosage_effect(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks
        ):
    """
    Parameters
    ----------
    gene_angles : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
    threshold : float
        filter threshold
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families) in column-major layout (order='F')
        bit mask that will have indices of genes kept by pattern set to 1, else 0

    Notes
    -----
    This subroutine prefilters the genes for dosage effect,
    as genes that are already too distant in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
    """

    # ensure all array inputs are numpy arrays
    gene_angles = np.ascontiguousarray(gene_angles, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)

    # extract dimension arguments
    n_genes = gene_angles.shape[0]


    # Create temporaries and/or outputs
    threshold = ctypes.c_double(threshold)
    masks = np.empty((n_mask_chunks, n_families), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.filter_paralogs_by_pattern_dosage_effect_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.filter_paralogs_by_pattern_dosage_effect_c.restype = None

    tox.filter_paralogs_by_pattern_dosage_effect_c(
        gene_angles,
        ctypes.byref(threshold),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        gene_to_fam,
        masks,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    masks.setflags(write=False)

    return masks


def filter_paralogs_by_pattern(
        pattern_mode,
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks
        ):
    """
    Parameters
    ----------
    pattern_mode : str
        used pattern for detection
        |       Mode           |                               Value                                |
        |----------------------|--------------------------------------------------------------------|
        |    Dosage Effect     |    "dosage_pattern"     |
        | Subfunctionalization |    "subfunc_pattern"    |
    gene_angles : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
    threshold : float
        filter threshold
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families) in column-major layout (order='F')
        bit mask that will have indices of genes kept by pattern set to 1, else 0

    Notes
    -----
    This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
    """

    # ensure all array inputs are numpy arrays
    pattern_mode = np.asarray(pattern_mode)
    gene_angles = np.ascontiguousarray(gene_angles, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)

    # extract dimension arguments
    n_genes = gene_angles.shape[0]


    # Create temporaries and/or outputs
    pattern_mode = pattern_mode.astype(f"S{15}", order="F")
    threshold = ctypes.c_double(threshold)
    masks = np.empty((n_mask_chunks, n_families), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.filter_paralogs_by_pattern_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{15}"),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.filter_paralogs_by_pattern_c.restype = None

    tox.filter_paralogs_by_pattern_c(
        pattern_mode,
        gene_angles,
        ctypes.byref(threshold),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        gene_to_fam,
        masks,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    masks.setflags(write=False)

    return masks


def calc_work_arr_paralog_subsets_size(
        max_subset_size,
        n_genes,
        filtered_paralogs_mask
        ):
    """
    Parameters
    ----------
    max_subset_size : int, modified in-place
        maximum size that a subset must not exceed.
        @warning
        If the desired size is too large and leads to an integer overflow, `max_subset_size` will be set to the maximum valid size.
        Also, size will be set to number of genes in `filtered_paralogs_mask` if larger.
        @endwarning
    n_genes : int
        number of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
        Output mask with all genes disabled that did not pass the filter
        M_DM_FROM_JUST_INFO to compute this argument using [[tox_paralog_analysis(module):filter_paralogs_by_pattern]]'s output `masks(:, family_idx)`.

    Returns
    -------
    work_array_size : int
        The calculated needed work array size in absolute worst case scenario. Look into source for details.

    Notes
    -----
    The `detect_*` subroutines need a work array for the to be tested subsets.
    In worst case, all need to be tried and subsets that cannot be extended will be kept as results.
    This is the reason why the work array holds the results as well, as all subsets that are stored in the array can be results as well.
    This subroutine calculates the needed size for the work array.
    """

    # ensure all array inputs are numpy arrays
    filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)

    # extract dimension arguments
    n_mask_chunks = filtered_paralogs_mask.shape[0]


    # Create temporaries and/or outputs
    max_subset_size = ctypes.c_int(max_subset_size)
    n_genes = ctypes.c_int(n_genes)
    work_array_size = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.calc_work_arr_paralog_subsets_size_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.calc_work_arr_paralog_subsets_size_c.restype = None

    tox.calc_work_arr_paralog_subsets_size_c(
        ctypes.byref(max_subset_size),
        ctypes.byref(n_genes),
        ctypes.byref(work_array_size),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return work_array_size.value


def mask_set_state(
        n_bit_mask_elements,
        i_gene,
        state
        ):
    """
    Parameters
    ----------
    i_gene : int
        index of paralog to be marked active
    state : bool
        state the bit should be set to

    Returns
    -------
    bit_mask : np.ndarray[np.int32] of shape (n_bit_mask_elements,) in column-major layout (order='F')
        chunked mask to mark active paralogs

    Notes
    -----
    Sets the state of a bit/gene in `bit_mask`
    """

    # ensure all array inputs are numpy arrays


    # extract dimension arguments



    # Create temporaries and/or outputs
    bit_mask = np.empty((n_bit_mask_elements,), dtype=np.int32, order='F')
    i_gene = ctypes.c_int(i_gene)
    state = ctypes.c_int(state)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.mask_set_state_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.mask_set_state_c.restype = None

    tox.mask_set_state_c(
        bit_mask,
        ctypes.byref(ctypes.c_int(n_bit_mask_elements)),
        ctypes.byref(i_gene),
        ctypes.byref(state),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    bit_mask.setflags(write=False)

    return bit_mask


def omics_vector_RAP_projection(
        vecs,
        vecs_selection_mask,
        axes_selection_mask
        ):
    """
    Parameters
    ----------
    vecs : np.ndarray[np.float64] of shape (n_axes, n_vecs) in column-major layout (order='F')
        matrix with expression vectors
    vecs_selection_mask : np.ndarray[np.int32] of shape (n_vecs,) in column-major layout (order='F')
        `.true.` for vectors where projection is to be computed
    axes_selection_mask : np.ndarray[np.int32] of shape (n_axes,) in column-major layout (order='F')
        `.true.` for axes to be included in RAP

    Returns
    -------
    projections : np.ndarray[np.float64] of shape (n_selected_axes, n_selected_vecs) in column-major layout (order='F')
        projected vectors

    Notes
    -----
    Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
    """

    # ensure all array inputs are numpy arrays
    vecs = np.asfortranarray(vecs, dtype=np.float64)
    vecs_selection_mask = np.ascontiguousarray(vecs_selection_mask, dtype=np.int32)
    axes_selection_mask = np.ascontiguousarray(axes_selection_mask, dtype=np.int32)

    # extract dimension arguments
    n_axes = vecs.shape[0]
    n_vecs = vecs.shape[1]
    n_selected_vecs = vecs_selection_mask.sum(axis=-1)
    n_selected_axes = axes_selection_mask.sum(axis=-1)

    # Create temporaries and/or outputs
    projections = np.empty((n_selected_axes, n_selected_vecs), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.omics_vector_RAP_projection_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.omics_vector_RAP_projection_c.restype = None

    tox.omics_vector_RAP_projection_c(
        vecs,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_vecs)),
        vecs_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_vecs)),
        axes_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_axes)),
        projections,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    projections.setflags(write=False)

    return projections


def omics_field_RAP_projection(
        vecs,
        vecs_selection_mask,
        axes_selection_mask
        ):
    """
    Parameters
    ----------
    vecs : np.ndarray[np.float64] of shape (2 * n_axes, n_vecs) in column-major layout (order='F')
        matrix with vector fields, first n rows mean vector origin, last n rows vector targets
    vecs_selection_mask : np.ndarray[np.int32] of shape (n_vecs,) in column-major layout (order='F')
        `.true.` for vectors where projection is to be computed
    axes_selection_mask : np.ndarray[np.int32] of shape (n_axes,) in column-major layout (order='F')
        `.true.` for axes to be included in RAP

    Returns
    -------
    projections : np.ndarray[np.float64] of shape (n_selected_axes, n_selected_vecs) in column-major layout (order='F')
        projected vectors

    Notes
    -----
    Project selected vector fields (e.g. shift vectors) onto the RAP constructed from a selected set of axes.
    """

    # ensure all array inputs are numpy arrays
    vecs = np.asfortranarray(vecs, dtype=np.float64)
    vecs_selection_mask = np.ascontiguousarray(vecs_selection_mask, dtype=np.int32)
    axes_selection_mask = np.ascontiguousarray(axes_selection_mask, dtype=np.int32)

    # extract dimension arguments
    n_axes = axes_selection_mask.shape[0]
    n_vecs = vecs.shape[1]
    n_selected_vecs = vecs_selection_mask.sum(axis=-1)
    n_selected_axes = axes_selection_mask.sum(axis=-1)

    # Create temporaries and/or outputs
    projections = np.empty((n_selected_axes, n_selected_vecs), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.omics_field_RAP_projection_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.omics_field_RAP_projection_c.restype = None

    tox.omics_field_RAP_projection_c(
        vecs,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_vecs)),
        vecs_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_vecs)),
        axes_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_axes)),
        projections,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    projections.setflags(write=False)

    return projections


def clock_hand_angle_between_vectors(
        v1,
        v2,
        selected_axes_for_signed
        ):
    """
    Parameters
    ----------
    v1 : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
        First normalized vector in RAP space
    v2 : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
        Second normalized vector in RAP space
    selected_axes_for_signed : np.ndarray[np.int32] of shape (3,) in column-major layout (order='F')
        Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)

    Returns
    -------
    signed_angle : float
        Signed angle between vectors in radians [-π, π]

    Notes
    -----
    Compute the signed clock hand angle between two RAP-projected and normalized vectors.
    Calculates the signed rotation angle between two normalized vectors in RAP space.
    For 2D/3D: automatic directionality calculation. For >3D: uses selected axes for directionality.
    """

    # ensure all array inputs are numpy arrays
    v1 = np.ascontiguousarray(v1, dtype=np.float64)
    v2 = np.ascontiguousarray(v2, dtype=np.float64)
    selected_axes_for_signed = np.ascontiguousarray(selected_axes_for_signed, dtype=np.int32)

    # extract dimension arguments
    n_dims = v1.shape[0]


    # Create temporaries and/or outputs
    signed_angle = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.clock_hand_angle_between_vectors_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.clock_hand_angle_between_vectors_c.restype = None

    tox.clock_hand_angle_between_vectors_c(
        v1,
        v2,
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(signed_angle),
        selected_axes_for_signed,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return signed_angle.value


def clock_hand_angles_for_shift_vectors(
        origins,
        targets,
        vecs_selection_mask,
        selected_axes_for_signed
        ):
    """
    Parameters
    ----------
    origins : np.ndarray[np.float64] of shape (n_dims, n_vecs) in column-major layout (order='F')
        First set of RAP-projected, normalized vectors (e.g. expression centroids)
    targets : np.ndarray[np.float64] of shape (n_dims, n_vecs) in column-major layout (order='F')
        Second set of RAP-projected, normalized vectors (e.g. paralogs)
    vecs_selection_mask : np.ndarray[np.int32] of shape (n_vecs,) in column-major layout (order='F')
        .true. for vector pairs where angle should be computed
    selected_axes_for_signed : np.ndarray[np.int32] of shape (3,) in column-major layout (order='F')
        Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)

    Returns
    -------
    signed_angles : np.ndarray[np.float64] of shape (n_selected_vecs,) in column-major layout (order='F')
        Signed rotation angles between vector pairs in radians [-π, π]

    Notes
    -----
    Compute signed rotation angles between RAP-projected and normalized vector pairs.
    Takes separate arrays of RAP-projected and normalized vectors (e.g. expression centroids and paralogs) and computes the signed rotation angle between corresponding pairs.
    This measures both magnitude and directionality of angular separation in RAP space.
    """

    # ensure all array inputs are numpy arrays
    origins = np.asfortranarray(origins, dtype=np.float64)
    targets = np.asfortranarray(targets, dtype=np.float64)
    vecs_selection_mask = np.ascontiguousarray(vecs_selection_mask, dtype=np.int32)
    selected_axes_for_signed = np.ascontiguousarray(selected_axes_for_signed, dtype=np.int32)

    # extract dimension arguments
    n_dims = origins.shape[0]
    n_vecs = origins.shape[1]
    n_selected_vecs = vecs_selection_mask.sum(axis=-1)

    # Create temporaries and/or outputs
    signed_angles = np.empty((n_selected_vecs,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.clock_hand_angles_for_shift_vectors_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.clock_hand_angles_for_shift_vectors_c.restype = None

    tox.clock_hand_angles_for_shift_vectors_c(
        origins,
        targets,
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(ctypes.c_int(n_vecs)),
        vecs_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_vecs)),
        selected_axes_for_signed,
        signed_angles,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    signed_angles.setflags(write=False)

    return signed_angles


def compute_relative_axis_contributions(
        vec
        ):
    """
    Parameters
    ----------
    vec : np.ndarray[np.float64] of shape (n_axes,) in column-major layout (order='F')
        RAP-projected and normalized vector (expression or shift)

    Returns
    -------
    contributions : np.ndarray[np.float64] of shape (n_axes,) in column-major layout (order='F')
        Fractional contribution of each axis (output), values in [0,1], sum to 1

    Notes
    -----
    Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
    Shared utility: computes fractional contribution of each axis to a RAP-projected and normalized vector.
    """

    # ensure all array inputs are numpy arrays
    vec = np.ascontiguousarray(vec, dtype=np.float64)

    # extract dimension arguments
    n_axes = vec.shape[0]


    # Create temporaries and/or outputs
    contributions = np.empty((n_axes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.compute_relative_axis_contributions_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.compute_relative_axis_contributions_c.restype = None

    tox.compute_relative_axis_contributions_c(
        vec,
        ctypes.byref(ctypes.c_int(n_axes)),
        contributions,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    contributions.setflags(write=False)

    return contributions


def relative_axes_changes_from_shift_vector(
        vec
        ):
    """
    Parameters
    ----------
    vec : np.ndarray[np.float64] of shape (n_axes,) in column-major layout (order='F')
        RAP-projected and normalized shift vector

    Returns
    -------
    contributions : np.ndarray[np.float64] of shape (n_axes,) in column-major layout (order='F')
        Fractional contribution of each axis (output), values in [0,1], sum to 1

    Notes
    -----
    Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
    Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
    """

    # ensure all array inputs are numpy arrays
    vec = np.ascontiguousarray(vec, dtype=np.float64)

    # extract dimension arguments
    n_axes = vec.shape[0]


    # Create temporaries and/or outputs
    contributions = np.empty((n_axes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.relative_axes_changes_from_shift_vector_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.relative_axes_changes_from_shift_vector_c.restype = None

    tox.relative_axes_changes_from_shift_vector_c(
        vec,
        ctypes.byref(ctypes.c_int(n_axes)),
        contributions,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    contributions.setflags(write=False)

    return contributions


def relative_axes_expression_from_expression_vector(
        vec
        ):
    """
    Parameters
    ----------
    vec : np.ndarray[np.float64] of shape (n_axes,) in column-major layout (order='F')
        RAP-projected and normalized expression vector

    Returns
    -------
    contributions : np.ndarray[np.float64] of shape (n_axes,) in column-major layout (order='F')
        Fractional contribution of each axis (output), values in [0,1], sum to 1

    Notes
    -----
    Compute fractional contribution of each axis to a RAP-projected and normalized expression vector.
    Wrapper for single RAP-projected expression vectors
    """

    # ensure all array inputs are numpy arrays
    vec = np.ascontiguousarray(vec, dtype=np.float64)

    # extract dimension arguments
    n_axes = vec.shape[0]


    # Create temporaries and/or outputs
    contributions = np.empty((n_axes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.relative_axes_expression_from_expression_vector_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.relative_axes_expression_from_expression_vector_c.restype = None

    tox.relative_axes_expression_from_expression_vector_c(
        vec,
        ctypes.byref(ctypes.c_int(n_axes)),
        contributions,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    contributions.setflags(write=False)

    return contributions
