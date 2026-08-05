"""Python binding to Generated from the kernel; do not edit -- regenerate instead.

Generated from tox_relative_axis_plane_tools. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.omics_vector_RAP_projection_c.restype = None
_lib.omics_vector_RAP_projection_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_OMICS_VECTOR_RAP_PROJECTION_ARGUMENTS = ("vecs", "n_axes", "n_vecs", "vecs_selection_mask", "n_selected_vecs", "axes_selection_mask", "n_selected_axes", "projections", "ierr",)

_lib.omics_field_RAP_projection_c.restype = None
_lib.omics_field_RAP_projection_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_OMICS_FIELD_RAP_PROJECTION_ARGUMENTS = ("fields", "n_axes", "n_fields", "fields_selection_mask", "n_selected_fields", "axes_selection_mask", "n_selected_axes", "projections", "ierr",)

_lib.clock_hand_angle_between_vectors_c.restype = None
_lib.clock_hand_angle_between_vectors_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CLOCK_HAND_ANGLE_BETWEEN_VECTORS_ARGUMENTS = ("v1", "v2", "n_dims", "signed_angle", "selected_axes_for_signed", "ierr",)

_lib.clock_hand_angles_for_shift_vectors_c.restype = None
_lib.clock_hand_angles_for_shift_vectors_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CLOCK_HAND_ANGLES_FOR_SHIFT_VECTORS_ARGUMENTS = ("fields", "n_dims", "n_fields", "fields_selection_mask", "n_selected_fields", "selected_axes_for_signed", "signed_angles", "ierr",)

_lib.compute_relative_axis_contributions_c.restype = None
_lib.compute_relative_axis_contributions_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_RELATIVE_AXIS_CONTRIBUTIONS_ARGUMENTS = ("vec", "n_axes", "contributions", "ierr",)

_lib.relative_axes_changes_from_shift_vector_c.restype = None
_lib.relative_axes_changes_from_shift_vector_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RELATIVE_AXES_CHANGES_FROM_SHIFT_VECTOR_ARGUMENTS = ("vec", "n_axes", "contributions", "ierr",)

_lib.relative_axes_expression_from_expression_vector_c.restype = None
_lib.relative_axes_expression_from_expression_vector_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RELATIVE_AXES_EXPRESSION_FROM_EXPRESSION_VECTOR_ARGUMENTS = ("vec", "n_axes", "contributions", "ierr",)

def omics_vector_RAP_projection(
        vecs,
        vecs_selection_mask,
        axes_selection_mask,
):
    r"""Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.

    Parameters
    ----------
    vecs : np.ndarray[np.float64] of shape (n_axes, n_vecs,), column-major (order='F')
        matrix with expression vectors
    vecs_selection_mask : np.ndarray[np.bool_] of shape (n_vecs,)
        `.true.` for vectors where projection is to be computed
    axes_selection_mask : np.ndarray[np.bool_] of shape (n_axes,)
        `.true.` for axes to be included in RAP

    Returns
    -------
    projections : np.ndarray[np.float64] of shape (n_selected_axes, n_selected_vecs,), column-major (order='F')
        projected vectors

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::omics_vector_RAP_projection`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        vecs = np.asfortranarray(vecs, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vecs' must be an array of np.float64: {error}") from None
    if vecs.ndim != 2:
        raise ValueError(f"'vecs' must have 2 dimensions, but has {vecs.ndim}")
    try:
        vecs_selection_mask = np.ascontiguousarray(vecs_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vecs_selection_mask' must be an array of np.bool_: {error}") from None
    if vecs_selection_mask.ndim != 1:
        raise ValueError(f"'vecs_selection_mask' must have 1 dimension, but has {vecs_selection_mask.ndim}")
    try:
        axes_selection_mask = np.ascontiguousarray(axes_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'axes_selection_mask' must be an array of np.bool_: {error}") from None
    if axes_selection_mask.ndim != 1:
        raise ValueError(f"'axes_selection_mask' must have 1 dimension, but has {axes_selection_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = vecs.shape[0]
    n_vecs = vecs.shape[1]
    n_selected_vecs = int(vecs_selection_mask.sum())
    n_selected_axes = int(axes_selection_mask.sum())

    # Fortran cannot check that shared extents agree; this can
    if axes_selection_mask.shape[0] != n_axes:
        raise ValueError(f"'axes_selection_mask' has {axes_selection_mask.shape[0]} along axis 0, but "
            f"'vecs' implies n_axes == {n_axes}"
        )
    if vecs_selection_mask.shape[0] != n_vecs:
        raise ValueError(f"'vecs_selection_mask' has {vecs_selection_mask.shape[0]} along axis 0, but "
            f"'vecs' implies n_vecs == {n_vecs}"
        )

    # outputs and work arrays, which the caller never sees
    projections = np.empty((n_selected_axes, n_selected_vecs,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.omics_vector_RAP_projection_c(
        vecs,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_vecs)),
        vecs_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_vecs)),
        axes_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_axes)),
        projections,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _OMICS_VECTOR_RAP_PROJECTION_ARGUMENTS)

    return projections

def omics_field_RAP_projection(
        fields,
        fields_selection_mask,
        axes_selection_mask,
):
    r"""Project selected vector fields (e.g. shift vectors) onto the RAP constructed from a selected set of axes.

    Parameters
    ----------
    fields : np.ndarray[np.float64] of shape (n_axes, 2, n_fields,), column-major (order='F')
        matrix with vector fields, `fields(:, 1, i_vec)` mean vector origin, `fields(:, 2, i_vec)` mean vector targets
    fields_selection_mask : np.ndarray[np.bool_] of shape (n_fields,)
        `.true.` for vectors where projection is to be computed
    axes_selection_mask : np.ndarray[np.bool_] of shape (n_axes,)
        `.true.` for axes to be included in RAP

    Returns
    -------
    projections : np.ndarray[np.float64] of shape (n_selected_axes, n_selected_fields,), column-major (order='F')
        projected vectors

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::omics_field_RAP_projection`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        fields = np.asfortranarray(fields, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'fields' must be an array of np.float64: {error}") from None
    if fields.ndim != 3:
        raise ValueError(f"'fields' must have 3 dimensions, but has {fields.ndim}")
    try:
        fields_selection_mask = np.ascontiguousarray(fields_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'fields_selection_mask' must be an array of np.bool_: {error}") from None
    if fields_selection_mask.ndim != 1:
        raise ValueError(f"'fields_selection_mask' must have 1 dimension, but has {fields_selection_mask.ndim}")
    try:
        axes_selection_mask = np.ascontiguousarray(axes_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'axes_selection_mask' must be an array of np.bool_: {error}") from None
    if axes_selection_mask.ndim != 1:
        raise ValueError(f"'axes_selection_mask' must have 1 dimension, but has {axes_selection_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = fields.shape[0]
    n_fields = fields.shape[2]
    n_selected_fields = int(fields_selection_mask.sum())
    n_selected_axes = int(axes_selection_mask.sum())

    # Fortran cannot check that shared extents agree; this can
    if axes_selection_mask.shape[0] != n_axes:
        raise ValueError(f"'axes_selection_mask' has {axes_selection_mask.shape[0]} along axis 0, but "
            f"'fields' implies n_axes == {n_axes}"
        )
    if fields_selection_mask.shape[0] != n_fields:
        raise ValueError(f"'fields_selection_mask' has {fields_selection_mask.shape[0]} along axis 0, but "
            f"'fields' implies n_fields == {n_fields}"
        )

    # outputs and work arrays, which the caller never sees
    projections = np.empty((n_selected_axes, n_selected_fields,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.omics_field_RAP_projection_c(
        fields,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_fields)),
        fields_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_fields)),
        axes_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_axes)),
        projections,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _OMICS_FIELD_RAP_PROJECTION_ARGUMENTS)

    return projections

def clock_hand_angle_between_vectors(
        v1,
        v2,
        selected_axes_for_signed,
):
    r"""Compute the signed clock hand angle between two RAP-projected and normalized vectors.

    Parameters
    ----------
    v1 : np.ndarray[np.float64] of shape (n_dims,)
        First normalized vector in RAP space
    v2 : np.ndarray[np.float64] of shape (n_dims,)
        Second normalized vector in RAP space
    selected_axes_for_signed : np.ndarray[np.int32] of shape (3,)
        Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)

    Returns
    -------
    signed_angle : float
        Signed angle between vectors in radians [-π, π]

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::clock_hand_angle_between_vectors`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        v1 = np.ascontiguousarray(v1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'v1' must be an array of np.float64: {error}") from None
    if v1.ndim != 1:
        raise ValueError(f"'v1' must have 1 dimension, but has {v1.ndim}")
    try:
        v2 = np.ascontiguousarray(v2, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'v2' must be an array of np.float64: {error}") from None
    if v2.ndim != 1:
        raise ValueError(f"'v2' must have 1 dimension, but has {v2.ndim}")
    try:
        selected_axes_for_signed = np.ascontiguousarray(selected_axes_for_signed, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'selected_axes_for_signed' must be an array of np.int32: {error}") from None
    if selected_axes_for_signed.ndim != 1:
        raise ValueError(f"'selected_axes_for_signed' must have 1 dimension, but has {selected_axes_for_signed.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dims = v1.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if v2.shape[0] != n_dims:
        raise ValueError(f"'v2' has {v2.shape[0]} along axis 0, but "
            f"'v1' implies n_dims == {n_dims}"
        )

    # outputs and work arrays, which the caller never sees
    signed_angle = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.clock_hand_angle_between_vectors_c(
        v1,
        v2,
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(signed_angle),
        selected_axes_for_signed,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CLOCK_HAND_ANGLE_BETWEEN_VECTORS_ARGUMENTS)

    return signed_angle.value

def clock_hand_angles_for_shift_vectors(
        fields,
        fields_selection_mask,
        selected_axes_for_signed,
):
    r"""Compute signed rotation angles between for shift vectors, so between their origin and target

    Parameters
    ----------
    fields : np.ndarray[np.float64] of shape (n_dims, 2, n_fields,), column-major (order='F')
        matrix with vector fields, `fields(:, 1, i_vec)` mean vector origin, `fields(:, 2, i_vec)` mean vector targets
    fields_selection_mask : np.ndarray[np.bool_] of shape (n_fields,)
        .true. for vector pairs where angle should be computed
    selected_axes_for_signed : np.ndarray[np.int32] of shape (3,)
        Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)

    Returns
    -------
    signed_angles : np.ndarray[np.float64] of shape (n_selected_fields,)
        Signed rotation angles between vector pairs in radians [-π, π]

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::clock_hand_angles_for_shift_vectors`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        fields = np.asfortranarray(fields, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'fields' must be an array of np.float64: {error}") from None
    if fields.ndim != 3:
        raise ValueError(f"'fields' must have 3 dimensions, but has {fields.ndim}")
    try:
        fields_selection_mask = np.ascontiguousarray(fields_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'fields_selection_mask' must be an array of np.bool_: {error}") from None
    if fields_selection_mask.ndim != 1:
        raise ValueError(f"'fields_selection_mask' must have 1 dimension, but has {fields_selection_mask.ndim}")
    try:
        selected_axes_for_signed = np.ascontiguousarray(selected_axes_for_signed, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'selected_axes_for_signed' must be an array of np.int32: {error}") from None
    if selected_axes_for_signed.ndim != 1:
        raise ValueError(f"'selected_axes_for_signed' must have 1 dimension, but has {selected_axes_for_signed.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dims = fields.shape[0]
    n_fields = fields.shape[2]
    n_selected_fields = int(fields_selection_mask.sum())

    # Fortran cannot check that shared extents agree; this can
    if fields_selection_mask.shape[0] != n_fields:
        raise ValueError(f"'fields_selection_mask' has {fields_selection_mask.shape[0]} along axis 0, but "
            f"'fields' implies n_fields == {n_fields}"
        )

    # outputs and work arrays, which the caller never sees
    signed_angles = np.empty((n_selected_fields,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.clock_hand_angles_for_shift_vectors_c(
        fields,
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(ctypes.c_int(n_fields)),
        fields_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_fields)),
        selected_axes_for_signed,
        signed_angles,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CLOCK_HAND_ANGLES_FOR_SHIFT_VECTORS_ARGUMENTS)

    return signed_angles

def compute_relative_axis_contributions(
        vec,
):
    r"""Compute the fractional contribution of each axis to a RAP-projected and normalized vector

    Parameters
    ----------
    vec : np.ndarray[np.float64] of shape (n_axes,)
        RAP-projected and normalized vector (expression or shift)

    Returns
    -------
    contributions : np.ndarray[np.float64] of shape (n_axes,)
        Fractional contribution of each axis (output), values in [0,1], sum to 1

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::compute_relative_axis_contributions`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        vec = np.ascontiguousarray(vec, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vec' must be an array of np.float64: {error}") from None
    if vec.ndim != 1:
        raise ValueError(f"'vec' must have 1 dimension, but has {vec.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = vec.shape[0]

    # outputs and work arrays, which the caller never sees
    contributions = np.empty((n_axes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_relative_axis_contributions_c(
        vec,
        ctypes.byref(ctypes.c_int(n_axes)),
        contributions,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_RELATIVE_AXIS_CONTRIBUTIONS_ARGUMENTS)

    return contributions

def relative_axes_changes_from_shift_vector(
        vec,
):
    r"""Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.

    Parameters
    ----------
    vec : np.ndarray[np.float64] of shape (n_axes,)
        RAP-projected and normalized shift vector

    Returns
    -------
    contributions : np.ndarray[np.float64] of shape (n_axes,)
        Fractional contribution of each axis (output), values in [0,1], sum to 1

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::relative_axes_changes_from_shift_vector`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        vec = np.ascontiguousarray(vec, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vec' must be an array of np.float64: {error}") from None
    if vec.ndim != 1:
        raise ValueError(f"'vec' must have 1 dimension, but has {vec.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = vec.shape[0]

    # outputs and work arrays, which the caller never sees
    contributions = np.empty((n_axes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.relative_axes_changes_from_shift_vector_c(
        vec,
        ctypes.byref(ctypes.c_int(n_axes)),
        contributions,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _RELATIVE_AXES_CHANGES_FROM_SHIFT_VECTOR_ARGUMENTS)

    return contributions

def relative_axes_expression_from_expression_vector(
        vec,
):
    r"""Compute fractional contribution of each axis to a RAP-projected and normalized expression vector.

    Parameters
    ----------
    vec : np.ndarray[np.float64] of shape (n_axes,)
        RAP-projected and normalized expression vector

    Returns
    -------
    contributions : np.ndarray[np.float64] of shape (n_axes,)
        Fractional contribution of each axis (output), values in [0,1], sum to 1

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::relative_axes_expression_from_expression_vector`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        vec = np.ascontiguousarray(vec, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vec' must be an array of np.float64: {error}") from None
    if vec.ndim != 1:
        raise ValueError(f"'vec' must have 1 dimension, but has {vec.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = vec.shape[0]

    # outputs and work arrays, which the caller never sees
    contributions = np.empty((n_axes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.relative_axes_expression_from_expression_vector_c(
        vec,
        ctypes.byref(ctypes.c_int(n_axes)),
        contributions,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _RELATIVE_AXES_EXPRESSION_FROM_EXPRESSION_VECTOR_ARGUMENTS)

    return contributions
