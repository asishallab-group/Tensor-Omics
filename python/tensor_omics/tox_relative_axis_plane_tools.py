r"""tox_relative_axis_plane_tools

Relative axis planes (RAPs): planes through higher-dimensional gene expression space, and
what can be read off a vector once it is projected onto one.

A RAP is picked by selecting axes (tissues) from the full expression space.
`vector_RAP_projection` projects a single vector onto it and `field_RAP_projection` a whole
field of them. Within the plane, `clock_hand_angle_between_vectors` measures the signed angle
between two vectors -- signed by an orientation reference, so the sign means the same thing in
every dimension -- and `clock_hand_angles_for_shift_vectors` does that for a whole shift
vector field at once.

`relative_axes_changes_from_shift_vector` and `relative_axes_expression_from_expression_vector`
give the per-axis breakdown instead of the angle: how much of a change, or of an expression
level, falls on each selected axis.

Python binding, generated from tox_relative_axis_plane_tools. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.omics_vector_RAP_projection_c.restype = None
_lib.omics_vector_RAP_projection_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_OMICS_VECTOR_RAP_PROJECTION_ARGUMENTS = ("vecs", "n_axes", "n_vecs", "vecs_selection_mask", "n_selected_vecs", "axes_selection_mask", "n_selected_axes", "projections", "ierr",)
#: For a derived argument, the one the caller passed it in
_OMICS_VECTOR_RAP_PROJECTION_ARGUMENT_SOURCES = (None, "vecs", "vecs", None, "projections", None, "projections", None, None,)

_lib.omics_field_RAP_projection_c.restype = None
_lib.omics_field_RAP_projection_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_OMICS_FIELD_RAP_PROJECTION_ARGUMENTS = ("fields", "n_axes", "n_fields", "fields_selection_mask", "n_selected_fields", "axes_selection_mask", "n_selected_axes", "projections", "ierr",)
#: For a derived argument, the one the caller passed it in
_OMICS_FIELD_RAP_PROJECTION_ARGUMENT_SOURCES = (None, "fields", "fields", None, "projections", None, "projections", None, None,)

_lib.clock_hand_angle_between_vectors_c.restype = None
_lib.clock_hand_angle_between_vectors_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CLOCK_HAND_ANGLE_BETWEEN_VECTORS_ARGUMENTS = ("v1", "v2", "n_dims", "orientation_reference", "signed_angle", "ierr",)
#: For a derived argument, the one the caller passed it in
_CLOCK_HAND_ANGLE_BETWEEN_VECTORS_ARGUMENT_SOURCES = (None, None, "v1", None, None, None,)

_lib.clock_hand_angles_for_shift_vectors_c.restype = None
_lib.clock_hand_angles_for_shift_vectors_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CLOCK_HAND_ANGLES_FOR_SHIFT_VECTORS_ARGUMENTS = ("fields", "n_dims", "n_fields", "fields_selection_mask", "n_selected_fields", "orientation_reference", "signed_angles", "ierr",)
#: For a derived argument, the one the caller passed it in
_CLOCK_HAND_ANGLES_FOR_SHIFT_VECTORS_ARGUMENT_SOURCES = (None, "fields", "fields", None, "signed_angles", None, None, None,)

_lib.compute_relative_axis_contributions_c.restype = None
_lib.compute_relative_axis_contributions_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_RELATIVE_AXIS_CONTRIBUTIONS_ARGUMENTS = ("vec", "n_axes", "contributions", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_RELATIVE_AXIS_CONTRIBUTIONS_ARGUMENT_SOURCES = (None, "vec", None, None,)

_lib.relative_axes_changes_from_shift_vector_c.restype = None
_lib.relative_axes_changes_from_shift_vector_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RELATIVE_AXES_CHANGES_FROM_SHIFT_VECTOR_ARGUMENTS = ("vec", "n_axes", "contributions", "ierr",)
#: For a derived argument, the one the caller passed it in
_RELATIVE_AXES_CHANGES_FROM_SHIFT_VECTOR_ARGUMENT_SOURCES = (None, "vec", None, None,)

_lib.relative_axes_expression_from_expression_vector_c.restype = None
_lib.relative_axes_expression_from_expression_vector_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RELATIVE_AXES_EXPRESSION_FROM_EXPRESSION_VECTOR_ARGUMENTS = ("vec", "n_axes", "contributions", "ierr",)
#: For a derived argument, the one the caller passed it in
_RELATIVE_AXES_EXPRESSION_FROM_EXPRESSION_VECTOR_ARGUMENT_SOURCES = (None, "vec", None, None,)

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
        `True` for vectors where projection is to be computed
    axes_selection_mask : np.ndarray[np.bool_] of shape (n_axes,)
        `True` for axes to be included in RAP

    Returns
    -------
    projections : np.ndarray[np.float64] of shape (n_selected_axes, n_selected_vecs,), column-major (order='F'), read-only
        projected vectors
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::omics_vector_RAP_projection`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    vecs = _contract.reals(vecs, "vecs", ndim=2, order="F")
    vecs_selection_mask = _contract.logicals(vecs_selection_mask, "vecs_selection_mask", ndim=1)
    axes_selection_mask = _contract.logicals(axes_selection_mask, "axes_selection_mask", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_axes = vecs.shape[0]
    n_vecs = vecs.shape[1]
    n_selected_vecs = int(vecs_selection_mask.sum())
    n_selected_axes = int(axes_selection_mask.sum())

    # Fortran cannot check that an array has the shape it is declared with; this can
    if vecs_selection_mask.shape[0] != n_vecs:
        raise _ToxArgumentError(f"'vecs_selection_mask' has {vecs_selection_mask.shape[0]} along axis 0, but 'vecs' implies n_vecs == {n_vecs}; 'vecs_selection_mask' is declared with shape (n_vecs)", argument="vecs_selection_mask")
    if axes_selection_mask.shape[0] != n_axes:
        raise _ToxArgumentError(f"'axes_selection_mask' has {axes_selection_mask.shape[0]} along axis 0, but 'vecs' implies n_axes == {n_axes}; 'axes_selection_mask' is declared with shape (n_axes)", argument="axes_selection_mask")

    # outputs and work arrays, which the caller never sees
    projections = _np.empty((n_selected_axes, n_selected_vecs,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.omics_vector_RAP_projection_c(
        vecs,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        _ctypes.byref(_ctypes.c_int(n_vecs)),
        vecs_selection_mask,
        _ctypes.byref(_ctypes.c_int(n_selected_vecs)),
        axes_selection_mask,
        _ctypes.byref(_ctypes.c_int(n_selected_axes)),
        projections,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _OMICS_VECTOR_RAP_PROJECTION_ARGUMENTS, _OMICS_VECTOR_RAP_PROJECTION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    projections.flags.writeable = False

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
        matrix with vector fields; each field holds two vectors, the origin first and the target second
    fields_selection_mask : np.ndarray[np.bool_] of shape (n_fields,)
        `True` for vectors where projection is to be computed
    axes_selection_mask : np.ndarray[np.bool_] of shape (n_axes,)
        `True` for axes to be included in RAP

    Returns
    -------
    projections : np.ndarray[np.float64] of shape (n_selected_axes, n_selected_fields,), column-major (order='F'), read-only
        projected vectors
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::omics_field_RAP_projection`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    fields = _contract.reals(fields, "fields", ndim=3, order="F")
    fields_selection_mask = _contract.logicals(fields_selection_mask, "fields_selection_mask", ndim=1)
    axes_selection_mask = _contract.logicals(axes_selection_mask, "axes_selection_mask", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_axes = fields.shape[0]
    n_fields = fields.shape[2]
    n_selected_fields = int(fields_selection_mask.sum())
    n_selected_axes = int(axes_selection_mask.sum())

    # Fortran cannot check that an array has the shape it is declared with; this can
    if fields.shape[1] != 2:
        raise _ToxArgumentError(f"'fields' has {fields.shape[1]} along axis 1, but its declared shape fixes it at {2}; 'fields' is declared with shape (n_axes, 2, n_fields)", argument="fields")
    if fields_selection_mask.shape[0] != n_fields:
        raise _ToxArgumentError(f"'fields_selection_mask' has {fields_selection_mask.shape[0]} along axis 0, but 'fields' implies n_fields == {n_fields}; 'fields_selection_mask' is declared with shape (n_fields)", argument="fields_selection_mask")
    if axes_selection_mask.shape[0] != n_axes:
        raise _ToxArgumentError(f"'axes_selection_mask' has {axes_selection_mask.shape[0]} along axis 0, but 'fields' implies n_axes == {n_axes}; 'axes_selection_mask' is declared with shape (n_axes)", argument="axes_selection_mask")

    # outputs and work arrays, which the caller never sees
    projections = _np.empty((n_selected_axes, n_selected_fields,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.omics_field_RAP_projection_c(
        fields,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        _ctypes.byref(_ctypes.c_int(n_fields)),
        fields_selection_mask,
        _ctypes.byref(_ctypes.c_int(n_selected_fields)),
        axes_selection_mask,
        _ctypes.byref(_ctypes.c_int(n_selected_axes)),
        projections,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _OMICS_FIELD_RAP_PROJECTION_ARGUMENTS, _OMICS_FIELD_RAP_PROJECTION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    projections.flags.writeable = False

    return projections

def clock_hand_angle_between_vectors(
        v1,
        v2,
        orientation_reference,
):
    r"""Compute the signed clock hand angle between two RAP-projected and normalized vectors.

    The unsigned angle is `acos(v1 . v2)`; `orientation_reference` supplies the sign by saying
    which way round the plane the two vectors span counts as positive. Reports
    `ERR_INVALID_INPUT` when the reference is orthogonal to the rotation and so orients nothing.

    Parameters
    ----------
    v1 : np.ndarray[np.float64] of shape (n_dims,)
        First normalized vector in RAP space
    v2 : np.ndarray[np.float64] of shape (n_dims,)
        Second normalized vector in RAP space
    orientation_reference : np.ndarray[np.float64] of shape (n_dims,)
        Orients the plane the rotation happens in, so the angle can carry a sign. A
        rotation from one vector to another has no inherent direction above two
        dimensions -- and in RAP space not even in two, since the axes are tissues or
        factors and carry no handedness -- so the caller states which way round counts
        as positive. The sign is that of this vector's component along the rotation.

    Returns
    -------
    signed_angle : float
        Signed angle between vectors in radians [-pi, pi]

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::clock_hand_angle_between_vectors`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    v1 = _contract.reals(v1, "v1", ndim=1)
    v2 = _contract.reals(v2, "v2", ndim=1)
    orientation_reference = _contract.reals(orientation_reference, "orientation_reference", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_dims = v1.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if v2.shape[0] != n_dims:
        raise _ToxArgumentError(f"'v2' has {v2.shape[0]} along axis 0, but 'v1' implies n_dims == {n_dims}; 'v2' is declared with shape (n_dims)", argument="v2")
    if orientation_reference.shape[0] != n_dims:
        raise _ToxArgumentError(f"'orientation_reference' has {orientation_reference.shape[0]} along axis 0, but 'v1' implies n_dims == {n_dims}; 'orientation_reference' is declared with shape (n_dims)", argument="orientation_reference")

    # outputs and work arrays, which the caller never sees
    signed_angle = _ctypes.c_double(0)
    ierr = _ctypes.c_int(0)

    _lib.clock_hand_angle_between_vectors_c(
        v1,
        v2,
        _ctypes.byref(_ctypes.c_int(n_dims)),
        orientation_reference,
        _ctypes.byref(signed_angle),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CLOCK_HAND_ANGLE_BETWEEN_VECTORS_ARGUMENTS, _CLOCK_HAND_ANGLE_BETWEEN_VECTORS_ARGUMENT_SOURCES)

    return signed_angle.value

def clock_hand_angles_for_shift_vectors(
        fields,
        fields_selection_mask,
        orientation_reference,
):
    r"""Compute signed rotation angles between for shift vectors, so between their origin and target

    Each selected field is angled by the rule of
    :func:`tensor_omics.clock_hand_angle_between_vectors`,
    with one `orientation_reference` shared by the whole batch. A single field whose rotation
    the reference fails to orient fails the call.

    Parameters
    ----------
    fields : np.ndarray[np.float64] of shape (n_dims, 2, n_fields,), column-major (order='F')
        matrix with vector fields; each field holds two vectors, the origin first and the target second
    fields_selection_mask : np.ndarray[np.bool_] of shape (n_fields,)
        True for vector pairs where angle should be computed
    orientation_reference : np.ndarray[np.float64] of shape (n_dims,)
        Orients the plane the rotation happens in, so the angle can carry a sign. A
        rotation from one vector to another has no inherent direction above two
        dimensions -- and in RAP space not even in two, since the axes are tissues or
        factors and carry no handedness -- so the caller states which way round counts
        as positive. The sign is that of this vector's component along the rotation.

    Returns
    -------
    signed_angles : np.ndarray[np.float64] of shape (n_selected_fields,), read-only
        Signed rotation angles between vector pairs in radians [-π, π]
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::clock_hand_angles_for_shift_vectors`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    fields = _contract.reals(fields, "fields", ndim=3, order="F")
    fields_selection_mask = _contract.logicals(fields_selection_mask, "fields_selection_mask", ndim=1)
    orientation_reference = _contract.reals(orientation_reference, "orientation_reference", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_dims = fields.shape[0]
    n_fields = fields.shape[2]
    n_selected_fields = int(fields_selection_mask.sum())

    # Fortran cannot check that an array has the shape it is declared with; this can
    if fields.shape[1] != 2:
        raise _ToxArgumentError(f"'fields' has {fields.shape[1]} along axis 1, but its declared shape fixes it at {2}; 'fields' is declared with shape (n_dims, 2, n_fields)", argument="fields")
    if fields_selection_mask.shape[0] != n_fields:
        raise _ToxArgumentError(f"'fields_selection_mask' has {fields_selection_mask.shape[0]} along axis 0, but 'fields' implies n_fields == {n_fields}; 'fields_selection_mask' is declared with shape (n_fields)", argument="fields_selection_mask")
    if orientation_reference.shape[0] != n_dims:
        raise _ToxArgumentError(f"'orientation_reference' has {orientation_reference.shape[0]} along axis 0, but 'fields' implies n_dims == {n_dims}; 'orientation_reference' is declared with shape (n_dims)", argument="orientation_reference")

    # outputs and work arrays, which the caller never sees
    signed_angles = _np.empty((n_selected_fields,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.clock_hand_angles_for_shift_vectors_c(
        fields,
        _ctypes.byref(_ctypes.c_int(n_dims)),
        _ctypes.byref(_ctypes.c_int(n_fields)),
        fields_selection_mask,
        _ctypes.byref(_ctypes.c_int(n_selected_fields)),
        orientation_reference,
        signed_angles,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CLOCK_HAND_ANGLES_FOR_SHIFT_VECTORS_ARGUMENTS, _CLOCK_HAND_ANGLES_FOR_SHIFT_VECTORS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    signed_angles.flags.writeable = False

    return signed_angles

def compute_relative_axis_contributions(
        vec,
):
    r"""Compute the fractional contribution of each axis to a RAP-projected and normalized vector

    Shared utility: the shift-vector and expression-vector entry points below both drive it.

    Parameters
    ----------
    vec : np.ndarray[np.float64] of shape (n_axes,)
        RAP-projected and normalized vector (expression or shift)

    Returns
    -------
    contributions : np.ndarray[np.float64] of shape (n_axes,), read-only
        Fractional contribution of each axis (output), values in [0,1], sum to 1
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::compute_relative_axis_contributions`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    vec = _contract.reals(vec, "vec", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_axes = vec.shape[0]

    # outputs and work arrays, which the caller never sees
    contributions = _np.empty((n_axes,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.compute_relative_axis_contributions_c(
        vec,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        contributions,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_RELATIVE_AXIS_CONTRIBUTIONS_ARGUMENTS, _COMPUTE_RELATIVE_AXIS_CONTRIBUTIONS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    contributions.flags.writeable = False

    return contributions

def relative_axes_changes_from_shift_vector(
        vec,
):
    r"""Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.

    Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)

    Parameters
    ----------
    vec : np.ndarray[np.float64] of shape (n_axes,)
        RAP-projected and normalized shift vector

    Returns
    -------
    contributions : np.ndarray[np.float64] of shape (n_axes,), read-only
        Fractional contribution of each axis (output), values in [0,1], sum to 1
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::relative_axes_changes_from_shift_vector`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    vec = _contract.reals(vec, "vec", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_axes = vec.shape[0]

    # outputs and work arrays, which the caller never sees
    contributions = _np.empty((n_axes,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.relative_axes_changes_from_shift_vector_c(
        vec,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        contributions,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _RELATIVE_AXES_CHANGES_FROM_SHIFT_VECTOR_ARGUMENTS, _RELATIVE_AXES_CHANGES_FROM_SHIFT_VECTOR_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    contributions.flags.writeable = False

    return contributions

def relative_axes_expression_from_expression_vector(
        vec,
):
    r"""Compute fractional contribution of each axis to a RAP-projected and normalized expression vector.

    Wrapper for single RAP-projected expression vectors

    Parameters
    ----------
    vec : np.ndarray[np.float64] of shape (n_axes,)
        RAP-projected and normalized expression vector

    Returns
    -------
    contributions : np.ndarray[np.float64] of shape (n_axes,), read-only
        Fractional contribution of each axis (output), values in [0,1], sum to 1
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_relative_axis_plane_tools::relative_axes_expression_from_expression_vector`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    vec = _contract.reals(vec, "vec", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_axes = vec.shape[0]

    # outputs and work arrays, which the caller never sees
    contributions = _np.empty((n_axes,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.relative_axes_expression_from_expression_vector_c(
        vec,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        contributions,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _RELATIVE_AXES_EXPRESSION_FROM_EXPRESSION_VECTOR_ARGUMENTS, _RELATIVE_AXES_EXPRESSION_FROM_EXPRESSION_VECTOR_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    contributions.flags.writeable = False

    return contributions
