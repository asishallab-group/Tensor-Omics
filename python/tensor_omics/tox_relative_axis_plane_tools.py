from .error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


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
    Compute the signed clock hand angle between two RAP-projected and normalized vectors.Calculates the signed rotation angle between two normalized vectors in RAP space.For 2D/3D: automatic directionality calculation. For >3D: uses selected axes for directionality.
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
    Compute signed rotation angles between RAP-projected and normalized vector pairs.Takes separate arrays of RAP-projected and normalized vectors (e.g. expression centroids and paralogs) and computes the signed rotation angle between corresponding pairs.This measures both magnitude and directionality of angular separation in RAP space.
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
    Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.Shared utility: computes fractional contribution of each axis to a RAP-projected and normalized vector.
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
    Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
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
    Compute fractional contribution of each axis to a RAP-projected and normalized expression vector.Wrapper for single RAP-projected expression vectors
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
