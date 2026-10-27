r"""tox_get_outliers

Gene outliers, from how far each gene sits from its family's centroid.

The pipeline is three steps, each callable on its own. `compute_rdi` turns raw distances into
a relative distance index, scaled per family so families of different spread are comparable.
`compute_family_scaling` fits that scaling with LOESS against family size. `identify_outliers`
applies the threshold and reports which genes exceed it.

`detect_outliers` runs all three in one call, and is the entry point to reach for first.

Python binding, generated from tox_get_outliers. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.compute_family_scaling_c.restype = None
_lib.compute_family_scaling_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_FAMILY_SCALING_ARGUMENTS = ("n_genes", "n_families", "distances", "gene_to_fam", "dscale", "loess_x", "loess_y", "indices_used", "span", "degree", "mode", "n_iters", "low_sd_cutoff", "excluded_low_sd", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_FAMILY_SCALING_ARGUMENT_SOURCES = ("distances", "dscale", None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.compute_rdi_c.restype = None
_lib.compute_rdi_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_RDI_ARGUMENTS = ("n_genes", "distances", "gene_to_fam", "dscale", "rdi", "sorted_rdi", "perm", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_RDI_ARGUMENT_SOURCES = ("distances", None, None, None, None, None, None, None,)

_lib.identify_outliers_c.restype = None
_lib.identify_outliers_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_IDENTIFY_OUTLIERS_ARGUMENTS = ("n_genes", "rdi", "sorted_rdi", "perm", "is_outlier", "threshold", "tail_probability", "quantile_level", "ierr",)
#: For a derived argument, the one the caller passed it in
_IDENTIFY_OUTLIERS_ARGUMENT_SOURCES = ("rdi", None, None, None, None, None, None, None, None,)

_lib.detect_outliers_c.restype = None
_lib.detect_outliers_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_OUTLIERS_ARGUMENTS = ("n_genes", "n_families", "distances", "gene_to_fam", "is_outlier", "loess_x", "loess_y", "loess_n", "tail_probability", "quantile_level", "ierr",)
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

    Uses LOESS on the median/stddev of intra-family distances for scaling, regardless of orthologs.

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
        The minimum valid value is `0`.
        The maximum valid value is `2`.
    mode : str, one of 'plain' | 'robust', optional, default 'robust'
        Mode for LOESS fitting
        The default value is `'robust'`.

        Trailing blanks are not part of the value, as in Fortran: they are dropped.
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
            Mean distance of each family in the LOESS fit, packed at the front in the order of
            `indices_used`. The remaining slots hold `-1`.
            A result is a value; call `.copy()` to obtain a modifiable array.
        loess_y : np.ndarray[np.float64] of shape (n_families,), read-only
            Standard deviation of the distances of each family in the LOESS fit, packed like
            `loess_x`. The remaining slots hold `-1`.
            A result is a value; call `.copy()` to obtain a modifiable array.
        indices_used : np.ndarray[np.int32] of shape (n_families,), read-only
            Family index of each point in the LOESS fit, packed at the front. Families with a
            single member, and those with the lowest spread, are left out of the fit; the
            remaining slots hold `0`.
            A result is a value; call `.copy()` to obtain a modifiable array.
        low_sd_cutoff : float
            cutoff used to filter families with low std
        excluded_low_sd : np.ndarray[np.int32] of shape (n_families,), read-only
            Mask to save those families that have low sd
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
    Generated from the Fortran procedure `tox_get_outliers::compute_family_scaling`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_families = _contract.integer(n_families, "n_families")
    distances = _contract.reals(distances, "distances", ndim=1)
    gene_to_fam = _contract.integers(gene_to_fam, "gene_to_fam", ndim=1)
    span = _contract.real(span, "span")
    degree = _contract.integer(degree, "degree")
    mode = _np.array([_contract.mode(mode, "mode", ("plain", "robust",)).ljust(6)], dtype="S6")
    n_iters = _contract.integer(n_iters, "n_iters")

    # what the inputs already say, rather than asking for it again
    n_genes = distances.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if gene_to_fam.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but 'distances' implies n_genes == {n_genes}; 'gene_to_fam' is declared with shape (n_genes)", argument="gene_to_fam")

    # outputs and work arrays, which the caller never sees
    dscale = _np.empty((n_families,), dtype=_np.float64, order='C')
    loess_x = _np.empty((n_families,), dtype=_np.float64, order='C')
    loess_y = _np.empty((n_families,), dtype=_np.float64, order='C')
    indices_used = _np.empty((n_families,), dtype=_np.int32, order='C')
    low_sd_cutoff = _ctypes.c_double(0)
    excluded_low_sd = _np.empty((n_families,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.compute_family_scaling_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_families)),
        distances,
        gene_to_fam,
        dscale,
        loess_x,
        loess_y,
        indices_used,
        _ctypes.byref(_ctypes.c_double(span)),
        _ctypes.byref(_ctypes.c_int(degree)),
        mode,
        _ctypes.byref(_ctypes.c_int(n_iters)),
        _ctypes.byref(low_sd_cutoff),
        excluded_low_sd,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_FAMILY_SCALING_ARGUMENTS, _COMPUTE_FAMILY_SCALING_ARGUMENT_SOURCES)

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

    RDI = Euclidean distance / family scaling factor

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers::compute_rdi`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    distances = _contract.reals(distances, "distances", ndim=1)
    gene_to_fam = _contract.integers(gene_to_fam, "gene_to_fam", ndim=1)
    dscale = _contract.reals(dscale, "dscale", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_genes = distances.shape[0]
    n_dscale_elements = dscale.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if gene_to_fam.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but 'distances' implies n_genes == {n_genes}; 'gene_to_fam' is declared with shape (n_genes)", argument="gene_to_fam")

    # outputs and work arrays, which the caller never sees
    rdi = _np.empty((n_genes,), dtype=_np.float64, order='C')
    sorted_rdi = _np.empty((n_genes,), dtype=_np.float64, order='C')
    perm = _np.empty((n_genes,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.compute_rdi_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        distances,
        gene_to_fam,
        dscale,
        _ctypes.byref(_ctypes.c_int(n_dscale_elements)),
        rdi,
        sorted_rdi,
        perm,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_RDI_ARGUMENTS, _COMPUTE_RDI_ARGUMENT_SOURCES)

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
        quantile_level=0.95,
):
    r"""Identify gene outliers based on an upper quantile of the RDI values

    Expects sorted_rdi to be filtered (no negative values) and perm should be sorted in ascending order before calling.
    If sorted_rdi contains negatives or perm is not sorted, tmp_results may be invalid.

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
    quantile_level : float, optional, default 0.95
        Quantile level of the threshold, as a fraction in [0,1] (the top 5% for the default).
        The default value is `0.95`.
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    dict
        with keys:

        is_outlier : np.ndarray[np.bool_] of shape (n_genes,), read-only
            Output boolean array indicating outliers
            A result is a value; call `.copy()` to obtain a modifiable array.
        threshold : float
            Output threshold value used for detection
        tail_probability : np.ndarray[np.float64] of shape (n_genes,), read-only
            Empirical one-sided upper-tail probability (effect-size measure) for each gene, i.e. how extreme an
            observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            upper-tail probability is used.
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
    Generated from the Fortran procedure `tox_get_outliers::identify_outliers`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    rdi = _contract.reals(rdi, "rdi", ndim=1)
    sorted_rdi = _contract.reals(sorted_rdi, "sorted_rdi", ndim=1)
    perm = _contract.integers(perm, "perm", ndim=1)
    quantile_level = _contract.real(quantile_level, "quantile_level")

    # what the inputs already say, rather than asking for it again
    n_genes = rdi.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if sorted_rdi.shape[0] != n_genes:
        raise _ToxArgumentError(f"'sorted_rdi' has {sorted_rdi.shape[0]} along axis 0, but 'rdi' implies n_genes == {n_genes}; 'sorted_rdi' is declared with shape (n_genes)", argument="sorted_rdi")
    if perm.shape[0] != n_genes:
        raise _ToxArgumentError(f"'perm' has {perm.shape[0]} along axis 0, but 'rdi' implies n_genes == {n_genes}; 'perm' is declared with shape (n_genes)", argument="perm")

    # outputs and work arrays, which the caller never sees
    is_outlier = _np.empty((n_genes,), dtype=_np.bool_, order='C')
    threshold = _ctypes.c_double(0)
    tail_probability = _np.empty((n_genes,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.identify_outliers_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        rdi,
        sorted_rdi,
        perm,
        is_outlier,
        _ctypes.byref(threshold),
        tail_probability,
        _ctypes.byref(_ctypes.c_double(quantile_level)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _IDENTIFY_OUTLIERS_ARGUMENTS, _IDENTIFY_OUTLIERS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    is_outlier.flags.writeable = False
    tail_probability.flags.writeable = False

    return {
        "is_outlier": is_outlier,
        "threshold": threshold.value,
        "tail_probability": tail_probability,
    }

def detect_outliers(
        n_families,
        distances,
        gene_to_fam,
        quantile_level=0.95,
):
    r"""Main routine to detect outliers using RDI and LOESS-based scaling

    Orchestrates the full pipeline: per-family scaling via
    :func:`tensor_omics.compute_family_scaling`, the RDI per gene via
    :func:`tensor_omics.compute_rdi`, then flags outliers via
    :func:`tensor_omics.identify_outliers`.

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
    quantile_level : float, optional, default 0.95
        Quantile level of the threshold, as a fraction in [0,1], for outlier detection.
        The default value is `0.95`.
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    dict
        with keys:

        is_outlier : np.ndarray[np.bool_] of shape (n_genes,), read-only
            Output boolean array indicating outliers
            A result is a value; call `.copy()` to obtain a modifiable array.
        loess_x : np.ndarray[np.float64] of shape (n_families,), read-only
            Mean distance of each family in the LOESS fit of the family scaling, packed at the
            front in the order of `loess_n`. The remaining slots hold `-1`.
            A result is a value; call `.copy()` to obtain a modifiable array.
        loess_y : np.ndarray[np.float64] of shape (n_families,), read-only
            Standard deviation of the distances of each family in that fit, packed like
            `loess_x`. The remaining slots hold `-1`.
            A result is a value; call `.copy()` to obtain a modifiable array.
        loess_n : np.ndarray[np.int32] of shape (n_families,), read-only
            Family index of each point in that fit, packed at the front. Families with a single
            member, and those with the lowest spread, are left out of the fit; the remaining
            slots hold `0`.
            A result is a value; call `.copy()` to obtain a modifiable array.
        tail_probability : np.ndarray[np.float64] of shape (n_genes,), read-only
            Empirical one-sided upper-tail probability (effect-size measure) for each gene, i.e. how extreme an
            observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            upper-tail probability is used.
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
    Generated from the Fortran procedure `tox_get_outliers::detect_outliers`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_families = _contract.integer(n_families, "n_families")
    distances = _contract.reals(distances, "distances", ndim=1)
    gene_to_fam = _contract.integers(gene_to_fam, "gene_to_fam", ndim=1)
    quantile_level = _contract.real(quantile_level, "quantile_level")

    # what the inputs already say, rather than asking for it again
    n_genes = distances.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if gene_to_fam.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but 'distances' implies n_genes == {n_genes}; 'gene_to_fam' is declared with shape (n_genes)", argument="gene_to_fam")

    # outputs and work arrays, which the caller never sees
    is_outlier = _np.empty((n_genes,), dtype=_np.bool_, order='C')
    loess_x = _np.empty((n_families,), dtype=_np.float64, order='C')
    loess_y = _np.empty((n_families,), dtype=_np.float64, order='C')
    loess_n = _np.empty((n_families,), dtype=_np.int32, order='C')
    tail_probability = _np.empty((n_genes,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.detect_outliers_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_families)),
        distances,
        gene_to_fam,
        is_outlier,
        loess_x,
        loess_y,
        loess_n,
        tail_probability,
        _ctypes.byref(_ctypes.c_double(quantile_level)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _DETECT_OUTLIERS_ARGUMENTS, _DETECT_OUTLIERS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    is_outlier.flags.writeable = False
    loess_x.flags.writeable = False
    loess_y.flags.writeable = False
    loess_n.flags.writeable = False
    tail_probability.flags.writeable = False

    return {
        "is_outlier": is_outlier,
        "loess_x": loess_x,
        "loess_y": loess_y,
        "loess_n": loess_n,
        "tail_probability": tail_probability,
    }
