"""tox_get_outliers

Generated from the implementation; do not edit -- regenerate instead.

Python binding, generated from tox_get_outliers. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.compute_family_scaling_c.restype = None
_lib.compute_family_scaling_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_FAMILY_SCALING_ARGUMENTS = ("n_genes", "n_families", "distances", "gene_to_fam", "dscale", "loess_x", "loess_y", "indices_used", "span", "degree", "mode", "n_iters", "low_sd_cutoff", "excluded_low_sd", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_FAMILY_SCALING_ARGUMENT_SOURCES = ("distances", "dscale", None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.compute_rdi_c.restype = None
_lib.compute_rdi_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_RDI_ARGUMENTS = ("n_genes", "distances", "gene_to_fam", "dscale", "rdi", "sorted_rdi", "perm", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_RDI_ARGUMENT_SOURCES = ("distances", None, None, None, None, None, None, None,)

_lib.identify_outliers_c.restype = None
_lib.identify_outliers_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_IDENTIFY_OUTLIERS_ARGUMENTS = ("n_genes", "rdi", "sorted_rdi", "perm", "is_outlier", "threshold", "quantile", "percentile", "ierr",)
#: For a derived argument, the one the caller passed it in
_IDENTIFY_OUTLIERS_ARGUMENT_SOURCES = ("rdi", None, None, None, None, None, None, None, None,)

_lib.detect_outliers_c.restype = None
_lib.detect_outliers_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_OUTLIERS_ARGUMENTS = ("n_genes", "n_families", "distances", "gene_to_fam", "is_outlier", "loess_x", "loess_y", "loess_n", "quantile", "ierr", "percentile",)
#: For a derived argument, the one the caller passed it in
_DETECT_OUTLIERS_ARGUMENT_SOURCES = ("distances", "loess_x", None, None, None, None, None, None, None, None, None,)

def compute_family_scaling(
        n_families,
        distances,
        gene_to_fam,
        span=0.7,
        degree=2,
        mode='robust',
        n_iters=3,
):
    r"""Compute family scaling factors (dscale) to normalize distances

    Parameters
    ----------
    n_families : int
        Total number of gene families
    distances : np.ndarray[np.float64] of shape (n_genes,)
        Array of Euclidean distances for each gene
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Mapping of each gene to its family (1-based)
    span : float, optional, default 0.7
        Span parameter for LOESS smoothing, passed straight to
        :func:`tensor_omics.loess_fit_plain`, so it is held to that
        procedure's own range rather than to the NaN tolerance the distance data carries.
        The default value is `0.7`.
        The minimum valid value is `EPS_LOESS`.
        The maximum valid value is `1.0`.
    degree : int, optional, default 2
        Degree of the LOESS polynomial
        The default value is `2`.
    mode : str, one of 'plain' | 'robust', optional, default 'robust'
        Mode for LOESS fitting
        The default value is `1`.

    n_iters : int, optional, default 3
        Number of iterations for robust LOESS fitting
        The default value is `3`.

    Returns
    -------
    dict
        with keys:

        dscale : np.ndarray[np.float64] of shape (n_families,), read-only
            Array of scaling factors per family (output)
            A result is a value; call `.copy()` to obtain a modifiable array.
        loess_x : np.ndarray[np.float64] of shape (n_families,), read-only
            Reference x-coordinates for LOESS smoothing
            A result is a value; call `.copy()` to obtain a modifiable array.
        loess_y : np.ndarray[np.float64] of shape (n_families,), read-only
            Reference y-coordinates for LOESS smoothing
            A result is a value; call `.copy()` to obtain a modifiable array.
        indices_used : np.ndarray[np.int32] of shape (n_families,), read-only
            Indices of reference points used for smoothing
            A result is a value; call `.copy()` to obtain a modifiable array.
        low_sd_cutoff : float
            cutoff used to filter families with low std
        excluded_low_sd : np.ndarray[np.int32] of shape (n_families,), read-only
            Mask to save those families that have low sd
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers::compute_family_scaling`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        distances = np.ascontiguousarray(distances, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'distances' must be an array of np.float64: {error}") from None
    if distances.ndim != 1:
        raise ValueError(f"'distances' must have 1 dimension, but has {distances.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")
    mode = np.array([str(mode).lower().encode()], dtype="S6")

    # what the inputs already say, rather than asking for it again
    n_genes = distances.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'distances' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    dscale = np.empty((n_families,), dtype=np.float64, order='C')
    loess_x = np.empty((n_families,), dtype=np.float64, order='C')
    loess_y = np.empty((n_families,), dtype=np.float64, order='C')
    indices_used = np.empty((n_families,), dtype=np.int32, order='C')
    low_sd_cutoff = ctypes.c_double(0)
    excluded_low_sd = np.empty((n_families,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_family_scaling_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        distances,
        gene_to_fam,
        dscale,
        loess_x,
        loess_y,
        indices_used,
        ctypes.byref(ctypes.c_double(span)),
        ctypes.byref(ctypes.c_int(degree)),
        mode,
        ctypes.byref(ctypes.c_int(n_iters)),
        ctypes.byref(low_sd_cutoff),
        excluded_low_sd,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_FAMILY_SCALING_ARGUMENTS, _COMPUTE_FAMILY_SCALING_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    dscale.flags.writeable = False
    loess_x.flags.writeable = False
    loess_y.flags.writeable = False
    indices_used.flags.writeable = False
    excluded_low_sd.flags.writeable = False

    return {
        "dscale": dscale,
        "loess_x": loess_x,
        "loess_y": loess_y,
        "indices_used": indices_used,
        "low_sd_cutoff": low_sd_cutoff.value,
        "excluded_low_sd": excluded_low_sd,
    }

def compute_rdi(
        distances,
        gene_to_fam,
        dscale,
):
    r"""Compute the hybrid RDI (Relative Distance Index) for each gene

    Parameters
    ----------
    distances : np.ndarray[np.float64] of shape (n_genes,)
        Array of Euclidean distances for each gene to its centroid
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Gene-to-family mapping (1-based indexing)
    dscale : np.ndarray[np.float64] of shape (n_dscale_elements,)
        Array of scaling factors for each family
        NaN is permitted for this value.
        Infinite values are permitted for this value.

    Returns
    -------
    dict
        with keys:

        rdi : np.ndarray[np.float64] of shape (n_genes,), read-only
            Output array of RDI values for each gene
            A result is a value; call `.copy()` to obtain a modifiable array.
        sorted_rdi : np.ndarray[np.float64] of shape (n_genes,), read-only
            Work array for sorting (dimension n_genes)
            A result is a value; call `.copy()` to obtain a modifiable array.
        perm : np.ndarray[np.int32] of shape (n_genes,), read-only
            Permutation array for sorting (dimension n_genes, should be pre-initialized with 1:n_genes)
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers::compute_rdi`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        distances = np.ascontiguousarray(distances, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'distances' must be an array of np.float64: {error}") from None
    if distances.ndim != 1:
        raise ValueError(f"'distances' must have 1 dimension, but has {distances.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")
    try:
        dscale = np.ascontiguousarray(dscale, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'dscale' must be an array of np.float64: {error}") from None
    if dscale.ndim != 1:
        raise ValueError(f"'dscale' must have 1 dimension, but has {dscale.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = distances.shape[0]
    n_dscale_elements = dscale.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'distances' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    rdi = np.empty((n_genes,), dtype=np.float64, order='C')
    sorted_rdi = np.empty((n_genes,), dtype=np.float64, order='C')
    perm = np.empty((n_genes,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_rdi_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        distances,
        gene_to_fam,
        dscale,
        ctypes.byref(ctypes.c_int(n_dscale_elements)),
        rdi,
        sorted_rdi,
        perm,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_RDI_ARGUMENTS, _COMPUTE_RDI_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    rdi.flags.writeable = False
    sorted_rdi.flags.writeable = False
    perm.flags.writeable = False

    return {
        "rdi": rdi,
        "sorted_rdi": sorted_rdi,
        "perm": perm,
    }

def identify_outliers(
        rdi,
        sorted_rdi,
        perm,
        percentile=95.0,
):
    r"""Identify gene outliers based on the top percentile of RDI values

    Parameters
    ----------
    rdi : np.ndarray[np.float64] of shape (n_genes,)
        Array of RDI values for each gene
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    sorted_rdi : np.ndarray[np.float64] of shape (n_genes,)
        Sorted RDI array (must be filtered to remove negatives and sorted in ascending order before calling)
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    perm : np.ndarray[np.int32] of shape (n_genes,)
        Permutation array with sorted indices
    percentile : float, optional, default 95.0
        Percentile threshold (top 5% for the default).
        The default value is `95.0`.

    Returns
    -------
    dict
        with keys:

        is_outlier : np.ndarray[np.bool_] of shape (n_genes,), read-only
            Output boolean array indicating outliers
            A result is a value; call `.copy()` to obtain a modifiable array.
        threshold : float
            Output threshold value used for detection
        quantile : np.ndarray[np.float64] of shape (n_genes,), read-only
            Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            upper-tail quantile is used.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers::identify_outliers`, whose argument names are
    the ones an error message reports.
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
        perm = np.ascontiguousarray(perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'perm' must be an array of np.int32: {error}") from None
    if perm.ndim != 1:
        raise ValueError(f"'perm' must have 1 dimension, but has {perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = rdi.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if sorted_rdi.shape[0] != n_genes:
        raise ValueError(f"'sorted_rdi' has {sorted_rdi.shape[0]} along axis 0, but "
            f"'rdi' implies n_genes == {n_genes}"
        )
    if perm.shape[0] != n_genes:
        raise ValueError(f"'perm' has {perm.shape[0]} along axis 0, but "
            f"'rdi' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    is_outlier = np.empty((n_genes,), dtype=np.bool_, order='C')
    threshold = ctypes.c_double(0)
    quantile = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.identify_outliers_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        rdi,
        sorted_rdi,
        perm,
        is_outlier,
        ctypes.byref(threshold),
        quantile,
        ctypes.byref(ctypes.c_double(percentile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _IDENTIFY_OUTLIERS_ARGUMENTS, _IDENTIFY_OUTLIERS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    is_outlier.flags.writeable = False
    quantile.flags.writeable = False

    return {
        "is_outlier": is_outlier,
        "threshold": threshold.value,
        "quantile": quantile,
    }

def detect_outliers(
        n_families,
        distances,
        gene_to_fam,
        percentile=95.0,
):
    r"""Main routine to detect outliers using RDI and LOESS-based scaling

    Parameters
    ----------
    n_families : int
        Total number of gene families
    distances : np.ndarray[np.float64] of shape (n_genes,)
        Array of Euclidean distances for each gene to its centroid
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Gene-to-family mapping (1-based indexing)
    percentile : float, optional, default 95.0
        Percentile threshold for outlier detection.
        The default value is `95.0`.

    Returns
    -------
    dict
        with keys:

        is_outlier : np.ndarray[np.bool_] of shape (n_genes,), read-only
            Output boolean array indicating outliers
            A result is a value; call `.copy()` to obtain a modifiable array.
        loess_x : np.ndarray[np.float64] of shape (n_families,), read-only
            Reference x-coordinates.
            A result is a value; call `.copy()` to obtain a modifiable array.
        loess_y : np.ndarray[np.float64] of shape (n_families,), read-only
            Reference y-coordinates (length n_total).
            A result is a value; call `.copy()` to obtain a modifiable array.
        loess_n : np.ndarray[np.int32] of shape (n_families,), read-only
            Indices of reference points used for smoothing.
            A result is a value; call `.copy()` to obtain a modifiable array.
        quantile : np.ndarray[np.float64] of shape (n_genes,), read-only
            Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            upper-tail quantile is used.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers::detect_outliers`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        distances = np.ascontiguousarray(distances, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'distances' must be an array of np.float64: {error}") from None
    if distances.ndim != 1:
        raise ValueError(f"'distances' must have 1 dimension, but has {distances.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = distances.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'distances' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    is_outlier = np.empty((n_genes,), dtype=np.bool_, order='C')
    loess_x = np.empty((n_families,), dtype=np.float64, order='C')
    loess_y = np.empty((n_families,), dtype=np.float64, order='C')
    loess_n = np.empty((n_families,), dtype=np.int32, order='C')
    quantile = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.detect_outliers_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        distances,
        gene_to_fam,
        is_outlier,
        loess_x,
        loess_y,
        loess_n,
        quantile,
        ctypes.byref(ierr),
        ctypes.byref(ctypes.c_double(percentile)),
    )

    check_err_code(ierr.value, _DETECT_OUTLIERS_ARGUMENTS, _DETECT_OUTLIERS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    is_outlier.flags.writeable = False
    loess_x.flags.writeable = False
    loess_y.flags.writeable = False
    loess_n.flags.writeable = False
    quantile.flags.writeable = False

    return {
        "is_outlier": is_outlier,
        "loess_x": loess_x,
        "loess_y": loess_y,
        "loess_n": loess_n,
        "quantile": quantile,
    }
