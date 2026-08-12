r"""f42_stats

Descriptive statistics: percentiles, empirical distribution functions, and 2-D LOESS smoothing.

One of the modules the `f42_utils` family gathers.

Python binding, generated from f42_stats. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.loess_smooth_2d_c.restype = None
_lib.loess_smooth_2d_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LOESS_SMOOTH_2D_ARGUMENTS = ("n_total", "n_target", "x_ref", "y_ref", "indices_used", "n_used", "x_query", "kernel_sigma", "kernel_cutoff", "y_out", "ierr",)
#: For a derived argument, the one the caller passed it in
_LOESS_SMOOTH_2D_ARGUMENT_SOURCES = ("x_ref", "x_query", None, None, None, "indices_used", None, None, None, None, None,)

_lib.compute_edf_c.restype = None
_lib.compute_edf_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_EDF_ARGUMENTS = ("values", "n_values", "unique_values", "cdf_values", "n_unique", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_EDF_ARGUMENT_SOURCES = (None, "values", None, None, None, None,)

_lib.compute_edf_expert_c.restype = None
_lib.compute_edf_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_EDF_EXPERT_ARGUMENTS = ("values", "n_values", "values_perm", "unique_values", "cdf_values", "n_unique", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_EDF_EXPERT_ARGUMENT_SOURCES = (None, "values", None, None, None, None, None,)

_lib.calc_percentile_c.restype = None
_lib.calc_percentile_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_PERCENTILE_ARGUMENTS = ("array", "n_array", "percentile", "value", "n_considered", "ierr",)
#: For a derived argument, the one the caller passed it in
_CALC_PERCENTILE_ARGUMENT_SOURCES = (None, "array", None, None, None, None,)

_lib.calc_percentile_expert_c.restype = None
_lib.calc_percentile_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_PERCENTILE_EXPERT_ARGUMENTS = ("array", "n_array", "array_perm", "percentile", "value", "n_considered", "ierr",)
#: For a derived argument, the one the caller passed it in
_CALC_PERCENTILE_EXPERT_ARGUMENT_SOURCES = (None, "array", None, None, None, None, None,)

_lib.compute_scaled_distance_quantile_c.restype = None
_lib.compute_scaled_distance_quantile_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_SCALED_DISTANCE_QUANTILE_ARGUMENTS = ("n_genes", "rdi", "sorted_rdi", "quantile", "c_const", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_SCALED_DISTANCE_QUANTILE_ARGUMENT_SOURCES = ("rdi", None, None, None, None, None,)

_lib.compute_scaled_distance_quantile_expert_c.restype = None
_lib.compute_scaled_distance_quantile_expert_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_SCALED_DISTANCE_QUANTILE_EXPERT_ARGUMENTS = ("n_genes", "rdi", "sorted_rdi", "sorted_rdi_perm", "quantile", "c_const", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_SCALED_DISTANCE_QUANTILE_EXPERT_ARGUMENT_SOURCES = ("rdi", None, None, None, None, None, None,)

def loess_smooth_2d(
        x_ref,
        y_ref,
        indices_used,
        x_query,
        kernel_sigma,
        kernel_cutoff,
):
    r"""Performs LOESS smoothing on a set of data points

    Parameters
    ----------
    x_ref : np.ndarray[np.float64] of shape (n_total,)
        Reference x-coordinates.
    y_ref : np.ndarray[np.float64] of shape (n_total,)
        Reference y-coordinates (length n_total).
    indices_used : np.ndarray[np.int32] of shape (n_used,)
        Indices of reference points used for smoothing (only valid indices).
        The minimum valid value is `1`.
        The maximum valid value is `n_total`.
    x_query : np.ndarray[np.float64] of shape (n_target,)
        Target x-coordinates to smooth.
    kernel_sigma : float
        Bandwidth parameter for the kernel.
        The minimum valid value is `0.0`.
    kernel_cutoff : float
        Cutoff for the kernel, not used if zero
        The minimum valid value is `0.0`.

    Returns
    -------
    y_out : np.ndarray[np.float64] of shape (n_target,), read-only
        Output smoothed values (length n_target).
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_stats::loess_smooth_2d`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        x_ref = np.ascontiguousarray(x_ref, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'x_ref' must be an array of np.float64: {error}") from None
    if x_ref.ndim != 1:
        raise ValueError(f"'x_ref' must have 1 dimension, but has {x_ref.ndim}")
    try:
        y_ref = np.ascontiguousarray(y_ref, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'y_ref' must be an array of np.float64: {error}") from None
    if y_ref.ndim != 1:
        raise ValueError(f"'y_ref' must have 1 dimension, but has {y_ref.ndim}")
    try:
        indices_used = np.ascontiguousarray(indices_used, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'indices_used' must be an array of np.int32: {error}") from None
    if indices_used.ndim != 1:
        raise ValueError(f"'indices_used' must have 1 dimension, but has {indices_used.ndim}")
    try:
        x_query = np.ascontiguousarray(x_query, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'x_query' must be an array of np.float64: {error}") from None
    if x_query.ndim != 1:
        raise ValueError(f"'x_query' must have 1 dimension, but has {x_query.ndim}")

    # what the inputs already say, rather than asking for it again
    n_total = x_ref.shape[0]
    n_target = x_query.shape[0]
    n_used = indices_used.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if y_ref.shape[0] != n_total:
        raise ValueError(f"'y_ref' has {y_ref.shape[0]} along axis 0, but "
            f"'x_ref' implies n_total == {n_total}"
        )

    # outputs and work arrays, which the caller never sees
    y_out = np.empty((n_target,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.loess_smooth_2d_c(
        ctypes.byref(ctypes.c_int(n_total)),
        ctypes.byref(ctypes.c_int(n_target)),
        x_ref,
        y_ref,
        indices_used,
        ctypes.byref(ctypes.c_int(n_used)),
        x_query,
        ctypes.byref(ctypes.c_double(kernel_sigma)),
        ctypes.byref(ctypes.c_double(kernel_cutoff)),
        y_out,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _LOESS_SMOOTH_2D_ARGUMENTS, _LOESS_SMOOTH_2D_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    y_out.flags.writeable = False

    return y_out

def compute_edf(
        values,
):
    r"""Compute the Empirical Distribution Function (EDF) from a sorted permutation

    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (n_values,)
        Array of observed data values (e.g., contributions or spikes).

    Returns
    -------
    dict
        with keys:

        unique_values : np.ndarray[np.float64] of shape (n_values,), read-only
            Sorted unique data values.
            The first `n_unique` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        cdf_values : np.ndarray[np.float64] of shape (n_values,), read-only
            Corresponding cumulative frequencies between 0 and 1.
            The first `n_unique` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_stats::compute_edf`, whose argument names are
    the ones an error message reports.

    This entry point seeds `values_perm` and sorts it by `values`.
    Call `compute_edf_expert` to do that yourself.
    """
    # accept anything array-like, converting only when C needs it
    try:
        values = np.ascontiguousarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'values' must be an array of np.float64: {error}") from None
    if values.ndim != 1:
        raise ValueError(f"'values' must have 1 dimension, but has {values.ndim}")

    # what the inputs already say, rather than asking for it again
    n_values = values.shape[0]

    # outputs and work arrays, which the caller never sees
    unique_values = np.empty((n_values,), dtype=np.float64, order='C')
    cdf_values = np.empty((n_values,), dtype=np.float64, order='C')
    n_unique = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.compute_edf_c(
        values,
        ctypes.byref(ctypes.c_int(n_values)),
        unique_values,
        cdf_values,
        ctypes.byref(n_unique),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_EDF_ARGUMENTS, _COMPUTE_EDF_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    unique_values.flags.writeable = False
    cdf_values.flags.writeable = False

    return {
        "unique_values": unique_values[..., :n_unique.value],
        "cdf_values": cdf_values[..., :n_unique.value],
    }

def compute_edf_expert(
        values,
        values_perm,
):
    r"""Compute the Empirical Distribution Function (EDF) from a sorted permutation

    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (n_values,)
        Array of observed data values (e.g., contributions or spikes).
    values_perm : np.ndarray[np.int32] of shape (n_values,)
        Permutation of `values` in ascending order. The allocating entry point builds
        and heapsorts it for you; the expert one takes whatever order you supply.
        The minimum valid value is `1`.
        The maximum valid value is `n_values`.

    Returns
    -------
    dict
        with keys:

        unique_values : np.ndarray[np.float64] of shape (n_values,), read-only
            Sorted unique data values.
            The first `n_unique` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        cdf_values : np.ndarray[np.float64] of shape (n_values,), read-only
            Corresponding cumulative frequencies between 0 and 1.
            The first `n_unique` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_stats::compute_edf_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `values_perm` yourself.
    `compute_edf` seeds `values_perm` and sorts it by `values`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        values = np.ascontiguousarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'values' must be an array of np.float64: {error}") from None
    if values.ndim != 1:
        raise ValueError(f"'values' must have 1 dimension, but has {values.ndim}")
    try:
        values_perm = np.ascontiguousarray(values_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'values_perm' must be an array of np.int32: {error}") from None
    if values_perm.ndim != 1:
        raise ValueError(f"'values_perm' must have 1 dimension, but has {values_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_values = values.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if values_perm.shape[0] != n_values:
        raise ValueError(f"'values_perm' has {values_perm.shape[0]} along axis 0, but "
            f"'values' implies n_values == {n_values}"
        )

    # outputs and work arrays, which the caller never sees
    unique_values = np.empty((n_values,), dtype=np.float64, order='C')
    cdf_values = np.empty((n_values,), dtype=np.float64, order='C')
    n_unique = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.compute_edf_expert_c(
        values,
        ctypes.byref(ctypes.c_int(n_values)),
        values_perm,
        unique_values,
        cdf_values,
        ctypes.byref(n_unique),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_EDF_EXPERT_ARGUMENTS, _COMPUTE_EDF_EXPERT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    unique_values.flags.writeable = False
    cdf_values.flags.writeable = False

    return {
        "unique_values": unique_values[..., :n_unique.value],
        "cdf_values": cdf_values[..., :n_unique.value],
    }

def calc_percentile(
        array,
        percentile,
        n_considered=0,
):
    r"""Calculate the percentile of an array given a sorted permutation

    Parameters
    ----------
    array : np.ndarray[np.float64] of shape (n_array,)
        input array
    percentile : float
        desired percentile (0-100)
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
    n_considered : int, optional, default 0
        How many leading entries of `array_perm` the percentile is taken over, for a
        percentile of a subset -- the trailing entries are ignored rather than sliced
        off, so the permutation stays the shape the sort produced. Zero, the default,
        considers all `n_array` of them.
        The default value is `0`.
        The minimum valid value is `0`.
        The maximum valid value is `n_array`.

    Returns
    -------
    value : float
        output percentile value

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_stats::calc_percentile`, whose argument names are
    the ones an error message reports.

    This entry point seeds `array_perm` and sorts it by `array`.
    Call `calc_percentile_expert` to do that yourself.
    """
    # accept anything array-like, converting only when C needs it
    try:
        array = np.ascontiguousarray(array, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'array' must be an array of np.float64: {error}") from None
    if array.ndim != 1:
        raise ValueError(f"'array' must have 1 dimension, but has {array.ndim}")

    # what the inputs already say, rather than asking for it again
    n_array = array.shape[0]

    # outputs and work arrays, which the caller never sees
    value = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.calc_percentile_c(
        array,
        ctypes.byref(ctypes.c_int(n_array)),
        ctypes.byref(ctypes.c_double(percentile)),
        ctypes.byref(value),
        ctypes.byref(ctypes.c_int(n_considered)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_PERCENTILE_ARGUMENTS, _CALC_PERCENTILE_ARGUMENT_SOURCES)

    return value.value

def calc_percentile_expert(
        array,
        array_perm,
        percentile,
        n_considered=0,
):
    r"""Calculate the percentile of an array given a sorted permutation

    Parameters
    ----------
    array : np.ndarray[np.float64] of shape (n_array,)
        input array
    array_perm : np.ndarray[np.int32] of shape (n_array,)
        Permutation of `array` in ascending order. The allocating entry point builds and
        heapsorts it for you; the expert one takes whatever order you supply.
        The minimum valid value is `1`.
        The maximum valid value is `n_array`.
    percentile : float
        desired percentile (0-100)
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
    n_considered : int, optional, default 0
        How many leading entries of `array_perm` the percentile is taken over, for a
        percentile of a subset -- the trailing entries are ignored rather than sliced
        off, so the permutation stays the shape the sort produced. Zero, the default,
        considers all `n_array` of them.
        The default value is `0`.
        The minimum valid value is `0`.
        The maximum valid value is `n_array`.

    Returns
    -------
    value : float
        output percentile value

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_stats::calc_percentile_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `array_perm` yourself.
    `calc_percentile` seeds `array_perm` and sorts it by `array`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        array = np.ascontiguousarray(array, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'array' must be an array of np.float64: {error}") from None
    if array.ndim != 1:
        raise ValueError(f"'array' must have 1 dimension, but has {array.ndim}")
    try:
        array_perm = np.ascontiguousarray(array_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'array_perm' must be an array of np.int32: {error}") from None
    if array_perm.ndim != 1:
        raise ValueError(f"'array_perm' must have 1 dimension, but has {array_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_array = array.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if array_perm.shape[0] != n_array:
        raise ValueError(f"'array_perm' has {array_perm.shape[0]} along axis 0, but "
            f"'array' implies n_array == {n_array}"
        )

    # outputs and work arrays, which the caller never sees
    value = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.calc_percentile_expert_c(
        array,
        ctypes.byref(ctypes.c_int(n_array)),
        array_perm,
        ctypes.byref(ctypes.c_double(percentile)),
        ctypes.byref(value),
        ctypes.byref(ctypes.c_int(n_considered)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_PERCENTILE_EXPERT_ARGUMENTS, _CALC_PERCENTILE_EXPERT_ARGUMENT_SOURCES)

    return value.value

def compute_scaled_distance_quantile(
        rdi,
        sorted_rdi,
        c_const,
):
    r"""Calculate the empirical quantile (effect-size measure) of scaled expression distances (RDI)

    Parameters
    ----------
    rdi : np.ndarray[np.float64] of shape (n_genes,)
        empirical distribution D
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    sorted_rdi : np.ndarray[np.float64] of shape (n_genes,)
        empirical distribution D with non negative values
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    c_const : float
        Constant used in the computation, typically 1

    Returns
    -------
    quantile : np.ndarray[np.float64] of shape (n_genes,), read-only
        Output array to store the computed quantile for each gene.
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_stats::compute_scaled_distance_quantile`, whose argument names are
    the ones an error message reports.

    This entry point seeds `sorted_rdi_perm` and sorts it by `sorted_rdi`.
    Call `compute_scaled_distance_quantile_expert` to do that yourself.
    """
    # accept anything array-like, converting only when C needs it
    try:
        rdi = np.ascontiguousarray(rdi, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'rdi' must be an array of np.float64: {error}") from None
    if rdi.ndim != 1:
        raise ValueError(f"'rdi' must have 1 dimension, but has {rdi.ndim}")
    try:
        sorted_rdi = np.ascontiguousarray(sorted_rdi, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'sorted_rdi' must be an array of np.float64: {error}") from None
    if sorted_rdi.ndim != 1:
        raise ValueError(f"'sorted_rdi' must have 1 dimension, but has {sorted_rdi.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = rdi.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if sorted_rdi.shape[0] != n_genes:
        raise ValueError(f"'sorted_rdi' has {sorted_rdi.shape[0]} along axis 0, but "
            f"'rdi' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    quantile = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_scaled_distance_quantile_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        rdi,
        sorted_rdi,
        quantile,
        ctypes.byref(ctypes.c_double(c_const)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_SCALED_DISTANCE_QUANTILE_ARGUMENTS, _COMPUTE_SCALED_DISTANCE_QUANTILE_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    quantile.flags.writeable = False

    return quantile

def compute_scaled_distance_quantile_expert(
        rdi,
        sorted_rdi,
        sorted_rdi_perm,
        c_const,
):
    r"""Calculate the empirical quantile (effect-size measure) of scaled expression distances (RDI)

    Parameters
    ----------
    rdi : np.ndarray[np.float64] of shape (n_genes,)
        empirical distribution D
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    sorted_rdi : np.ndarray[np.float64] of shape (n_genes,)
        empirical distribution D with non negative values
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    sorted_rdi_perm : np.ndarray[np.int32] of shape (n_genes,)
        Permutation of `sorted_rdi` in ascending order. The allocating entry point builds
        and heapsorts it for you; the expert one takes whatever order you supply.
        The minimum valid value is `1`.
        The maximum valid value is `n_genes`.
    c_const : float
        Constant used in the computation, typically 1

    Returns
    -------
    quantile : np.ndarray[np.float64] of shape (n_genes,), read-only
        Output array to store the computed quantile for each gene.
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_stats::compute_scaled_distance_quantile_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `sorted_rdi_perm` yourself.
    `compute_scaled_distance_quantile` seeds `sorted_rdi_perm` and sorts it by `sorted_rdi`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        rdi = np.ascontiguousarray(rdi, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'rdi' must be an array of np.float64: {error}") from None
    if rdi.ndim != 1:
        raise ValueError(f"'rdi' must have 1 dimension, but has {rdi.ndim}")
    try:
        sorted_rdi = np.ascontiguousarray(sorted_rdi, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'sorted_rdi' must be an array of np.float64: {error}") from None
    if sorted_rdi.ndim != 1:
        raise ValueError(f"'sorted_rdi' must have 1 dimension, but has {sorted_rdi.ndim}")
    try:
        sorted_rdi_perm = np.ascontiguousarray(sorted_rdi_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'sorted_rdi_perm' must be an array of np.int32: {error}") from None
    if sorted_rdi_perm.ndim != 1:
        raise ValueError(f"'sorted_rdi_perm' must have 1 dimension, but has {sorted_rdi_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = rdi.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if sorted_rdi.shape[0] != n_genes:
        raise ValueError(f"'sorted_rdi' has {sorted_rdi.shape[0]} along axis 0, but "
            f"'rdi' implies n_genes == {n_genes}"
        )
    if sorted_rdi_perm.shape[0] != n_genes:
        raise ValueError(f"'sorted_rdi_perm' has {sorted_rdi_perm.shape[0]} along axis 0, but "
            f"'rdi' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    quantile = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_scaled_distance_quantile_expert_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        rdi,
        sorted_rdi,
        sorted_rdi_perm,
        quantile,
        ctypes.byref(ctypes.c_double(c_const)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_SCALED_DISTANCE_QUANTILE_EXPERT_ARGUMENTS, _COMPUTE_SCALED_DISTANCE_QUANTILE_EXPERT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    quantile.flags.writeable = False

    return quantile
