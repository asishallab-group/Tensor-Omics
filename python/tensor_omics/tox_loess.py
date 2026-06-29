from .error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


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
    Recommend workspace sizes based on Netlib exact formulas.Computes the required sizes for integer and real workspace arrays.These sizes depend on the dimensionality of the data and the maximum neighborhood size.
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
    Perform plain LOESS fitting.Fits a LOESS model to the data using the specified smoothing parameter.Outputs the smoothed response variable array.
    """

    # ensure all array inputs are numpy arrays
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    w = np.ascontiguousarray(w, dtype=np.float64)
    eval_points = np.asfortranarray(eval_points, dtype=np.float64)
    assert type(int_workspace) is np.ndarray and int_workspace.flags.f_contiguous and int_workspace.dtype == np.int32, "'int_workspace' must be column-major numpy array (order='F')"
    assert type(real_workspace) is np.ndarray and real_workspace.flags.f_contiguous and real_workspace.dtype == np.float64, "'real_workspace' must be column-major numpy array (order='F')"
    assert type(diagl) is np.ndarray and diagl.flags.f_contiguous and diagl.dtype == np.float64, "'diagl' must be column-major numpy array (order='F')"

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
    Perform robust LOESS fitting with bisquare reweighting.Fits a LOESS model to the data using robust iterations to handle outliers.The robust fitting process iterates n_iters times, each iteration:- Combines original weights with robust weights (down-weights from previous iteration)- Runs LOESS fitting with combined weights- Computes residuals (y - fitted values)- Updates robust weights using bisquare function (suppresses large residuals)
    """

    # ensure all array inputs are numpy arrays
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    w = np.ascontiguousarray(w, dtype=np.float64)
    eval_points = np.asfortranarray(eval_points, dtype=np.float64)
    assert type(int_workspace) is np.ndarray and int_workspace.flags.f_contiguous and int_workspace.dtype == np.int32, "'int_workspace' must be column-major numpy array (order='F')"
    assert type(real_workspace) is np.ndarray and real_workspace.flags.f_contiguous and real_workspace.dtype == np.float64, "'real_workspace' must be column-major numpy array (order='F')"
    assert type(diagl) is np.ndarray and diagl.flags.f_contiguous and diagl.dtype == np.float64, "'diagl' must be column-major numpy array (order='F')"
    assert type(robust_weights) is np.ndarray and robust_weights.flags.f_contiguous and robust_weights.dtype == np.float64, "'robust_weights' must be column-major numpy array (order='F')"
    assert type(combined_weights) is np.ndarray and combined_weights.flags.f_contiguous and combined_weights.dtype == np.float64, "'combined_weights' must be column-major numpy array (order='F')"
    assert type(residuals) is np.ndarray and residuals.flags.f_contiguous and residuals.dtype == np.float64, "'residuals' must be column-major numpy array (order='F')"
    assert type(permutation_indices) is np.ndarray and permutation_indices.flags.f_contiguous and permutation_indices.dtype == np.int32, "'permutation_indices' must be column-major numpy array (order='F')"

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
        Mode of operation|      Mode      |  Value   ||----------------|----------|| robust fitting | "robust" || plain fitting  | "plain"  |
    n_iters : int
        Number of robust iterations (only used when mode = 1)

    Returns
    -------
    fitted_values : np.ndarray[np.float64] of shape (size(y),) in column-major layout (order='F')
        Fitted (smoothed) values of y

    Notes
    -----
    Wrapper subroutine for LOESS fitting (plain or robust).This subroutine selects between plain and robust LOESS fitting based on the mode.It dynamically allocates the required arrays and computes workspace sizes.
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
