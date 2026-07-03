from .error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


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
    Performs LOESS smoothing on a set of data points.Smooths y_ref at x_query using reference points x_ref, y_ref, and kernel parameters.The user must pre-filter data and provide only valid indices in indices_used.
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


def compute_edf(
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
    Compute the Empirical Distribution Function (EDF) from pre-sorted permutation.Returns the sorted unique values and their cumulative frequencies in [0,1].Assumes perm is already sorted by values[perm]. Caller controls sorting algorithm.The number of unique values can be determined by finding the last non-zero cdf_value.
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

    tox.compute_edf_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.compute_edf_c.restype = None

    tox.compute_edf_c(
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
    Helper routine that sorts and calls compute_edf.Allocates workspace internally and performs sorting before computing EDF.Use this for convenience; use compute_edf directly for custom sorting.
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
    Calculate empirical p-values for scaled expression distances (RDI).Implements:P(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )Because distances are non-negative, a one-sided upper-tail empirical p-value is used.Assumptions / preconditions:- sorted_rdi(1:n_genes) contains the empirical distribution D.- If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
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
