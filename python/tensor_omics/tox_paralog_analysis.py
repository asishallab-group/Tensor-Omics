"""Python interface to Module for detecting paralog-subset expression patterns (dosage effect and subfunctionalization) relative to an ancestral ortholog.

Generated from tox_paralog_analysis. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.mask_check_state_c.restype = None
_lib.mask_check_state_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MASK_CHECK_STATE_ARGUMENTS = ("bit_mask", "i_gene",)

_lib.detect_neofunctionalization_c.restype = None
_lib.detect_neofunctionalization_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_NEOFUNCTIONALIZATION_ARGUMENTS = ("ancestors", "n_families", "genes", "n_axes", "gene_to_fam", "n_genes", "thresholds", "neofunc", "ierr",)

_lib.detect_dosage_effect_c.restype = None
_lib.detect_dosage_effect_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    nullable(ctypes.POINTER(ctypes.c_double)),
    nullable(ctypes.POINTER(ctypes.c_double)),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_DOSAGE_EFFECT_ARGUMENTS = ("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "tmp_active_mask", "tmp_paralog_vector", "ierr", "max_angle", "gain_gamma",)

_lib.detect_subfunctionalization_c.restype = None
_lib.detect_subfunctionalization_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_SUBFUNCTIONALIZATION_ARGUMENTS = ("ancestor", "genes", "n_genes", "n_dims", "rdi_threshold", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "tmp_active_mask", "tmp_paralog_vector", "paralog_norms", "sorted_paralog_norms_perm", "tmp_work_array", "ierr",)

_lib.mask_chunk_count_c.restype = None
_lib.mask_chunk_count_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MASK_CHUNK_COUNT_ARGUMENTS = ("n_genes", "count",)

_lib.filter_paralogs_by_pattern_subfunctionalization_c.restype = None
_lib.filter_paralogs_by_pattern_subfunctionalization_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FILTER_PARALOGS_BY_PATTERN_SUBFUNCTIONALIZATION_ARGUMENTS = ("gene_angles", "threshold", "n_genes", "n_families", "gene_to_fam", "masks", "n_mask_chunks", "ierr",)

_lib.filter_paralogs_by_pattern_dosage_effect_c.restype = None
_lib.filter_paralogs_by_pattern_dosage_effect_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FILTER_PARALOGS_BY_PATTERN_DOSAGE_EFFECT_ARGUMENTS = ("gene_angles", "threshold", "n_genes", "n_families", "gene_to_fam", "masks", "n_mask_chunks", "ierr",)

_lib.calc_work_arr_paralog_subsets_size_c.restype = None
_lib.calc_work_arr_paralog_subsets_size_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_WORK_ARR_PARALOG_SUBSETS_SIZE_ARGUMENTS = ("max_subset_size", "n_genes", "work_array_size", "filtered_paralogs_mask", "n_mask_chunks", "ierr",)

def mask_check_state(
        bit_mask,
        i_gene,
):
    r"""Checks the state of a bit/paralog in `bit_mask` -> .true. if 1 else .false.

    Parameters
    ----------
    bit_mask : np.ndarray[np.int32] of shape (n_bit_mask_elements,)
        chunked mask to mark active paralogs
    i_gene : int
        index of paralog to be marked active

    Returns
    -------
    state : bool
        check result

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis::mask_check_state`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        bit_mask = np.ascontiguousarray(bit_mask, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'bit_mask' must be an array of np.int32: {error}") from None
    if bit_mask.ndim != 1:
        raise ValueError(f"'bit_mask' must have 1 dimension, but has {bit_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_bit_mask_elements = bit_mask.shape[0]

    # outputs and work arrays, which the caller never sees
    state = ctypes.c_bool(0)
    ierr = ctypes.c_int(0)

    _lib.mask_check_state_c(
        bit_mask,
        ctypes.byref(ctypes.c_int(n_bit_mask_elements)),
        ctypes.byref(ctypes.c_int(i_gene)),
        ctypes.byref(state),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _MASK_CHECK_STATE_ARGUMENTS)

    return state.value

def detect_neofunctionalization(
        ancestors,
        genes,
        gene_to_fam,
        thresholds,
):
    r"""Identifies neofunctionalization for genes by checking whether the difference of expression to its ancestor exceeds the threshold for the respective axis

    Parameters
    ----------
    ancestors : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F')
        RAP projected unit length expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        RAP projected unit length expression vectors of genes
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0_int32` for unassigned genes
    thresholds : np.ndarray[np.float64] of shape (n_axes,)
        threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis

    Returns
    -------
    neofunc : np.ndarray[np.bool_] of shape (n_genes, n_axes,), column-major (order='F')
        `.true.` if neofunctionalization has been detected for the respective axes, always `.false.` for unassigned genes

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis::detect_neofunctionalization`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ancestors = np.asfortranarray(ancestors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ancestors' must be an array of np.float64: {error}") from None
    if ancestors.ndim != 2:
        raise ValueError(f"'ancestors' must have 2 dimensions, but has {ancestors.ndim}")
    try:
        genes = np.asfortranarray(genes, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'genes' must be an array of np.float64: {error}") from None
    if genes.ndim != 2:
        raise ValueError(f"'genes' must have 2 dimensions, but has {genes.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")
    try:
        thresholds = np.ascontiguousarray(thresholds, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'thresholds' must be an array of np.float64: {error}") from None
    if thresholds.ndim != 1:
        raise ValueError(f"'thresholds' must have 1 dimension, but has {thresholds.ndim}")

    # what the inputs already say, rather than asking for it again
    n_families = ancestors.shape[1]
    n_axes = ancestors.shape[0]
    n_genes = genes.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if genes.shape[0] != n_axes:
        raise ValueError(f"'genes' has {genes.shape[0]} along axis 0, but "
            f"'ancestors' implies n_axes == {n_axes}"
        )
    if thresholds.shape[0] != n_axes:
        raise ValueError(f"'thresholds' has {thresholds.shape[0]} along axis 0, but "
            f"'ancestors' implies n_axes == {n_axes}"
        )
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'genes' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    neofunc = np.empty((n_genes, n_axes,), dtype=np.bool_, order='F')
    ierr = ctypes.c_int(0)

    _lib.detect_neofunctionalization_c(
        ancestors,
        ctypes.byref(ctypes.c_int(n_families)),
        genes,
        ctypes.byref(ctypes.c_int(n_axes)),
        gene_to_fam,
        ctypes.byref(ctypes.c_int(n_genes)),
        thresholds,
        neofunc,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETECT_NEOFUNCTIONALIZATION_ARGUMENTS)

    return neofunc

def detect_dosage_effect(
        ancestor,
        genes,
        filtered_paralogs_mask,
        max_subset_size,
        n_paralog_subsets,
        max_angle=None,
        gain_gamma=None,
):
    r"""Identifies subsets of paralogs with small angle to the `ancestor` (max_angle) and sum to a magnitude significantly exceeding `norm(ancestor)` (gain)

    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,)
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes,), column-major (order='F')
        expression vectors of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,)
        bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    max_subset_size : int
        maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    n_paralog_subsets : int
        number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    max_angle : float, optional
        in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned, default is Pi
    gain_gamma : float, optional
        positive magnitude gain for dosage effect, default 0.1

    Returns
    -------
    dict
        with keys:

        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_paralog_subsets,), column-major (order='F')
            working array to hold bitmask encoded subsets for detection.
            @note
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0_real64)` and represents the number of chunks
            @endnote

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis::detect_dosage_effect`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ancestor = np.ascontiguousarray(ancestor, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ancestor' must be an array of np.float64: {error}") from None
    if ancestor.ndim != 1:
        raise ValueError(f"'ancestor' must have 1 dimension, but has {ancestor.ndim}")
    try:
        genes = np.asfortranarray(genes, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'genes' must be an array of np.float64: {error}") from None
    if genes.ndim != 2:
        raise ValueError(f"'genes' must have 2 dimensions, but has {genes.ndim}")
    try:
        filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'filtered_paralogs_mask' must be an array of np.int32: {error}") from None
    if filtered_paralogs_mask.ndim != 1:
        raise ValueError(f"'filtered_paralogs_mask' must have 1 dimension, but has {filtered_paralogs_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if genes.shape[0] != n_dims:
        raise ValueError(f"'genes' has {genes.shape[0]} along axis 0, but "
            f"'ancestor' implies n_dims == {n_dims}"
        )

    # outputs and work arrays, which the caller never sees
    n_results = ctypes.c_int(0)
    work_arr_paralog_subsets = np.empty((n_mask_chunks, n_paralog_subsets,), dtype=np.int32, order='F')
    tmp_active_mask = np.empty((n_mask_chunks,), dtype=np.int32, order='C')
    tmp_paralog_vector = np.empty((n_dims,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.detect_dosage_effect_c(
        ancestor,
        genes,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_dims)),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(n_results),
        ctypes.byref(ctypes.c_int(max_subset_size)),
        work_arr_paralog_subsets,
        ctypes.byref(ctypes.c_int(n_paralog_subsets)),
        tmp_active_mask,
        tmp_paralog_vector,
        ctypes.byref(ierr),
        None if max_angle is None else ctypes.byref(ctypes.c_double(max_angle)),
        None if gain_gamma is None else ctypes.byref(ctypes.c_double(gain_gamma)),
    )

    check_err_code(ierr.value, _DETECT_DOSAGE_EFFECT_ARGUMENTS)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets,
    }

def detect_subfunctionalization(
        ancestor,
        genes,
        rdi_threshold,
        filtered_paralogs_mask,
        max_subset_size,
        n_paralog_subsets,
        paralog_norms,
        sorted_paralog_norms_perm,
):
    r"""Identifies subsets of paralogs exhibiting significant angles to the `ancestor`

    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,)
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes,), column-major (order='F')
        expression vectors of genes
    rdi_threshold : float
        max allowed residual distance from `ancestor`
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,)
        bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    max_subset_size : int
        maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    n_paralog_subsets : int
        number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    paralog_norms : np.ndarray[np.float64] of shape (n_genes,)
        needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` function from `f42_utils` function for this)
    sorted_paralog_norms_perm : np.ndarray[np.int32] of shape (n_genes,)
        needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor

    Returns
    -------
    dict
        with keys:

        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_paralog_subsets,), column-major (order='F')
            working array to hold bitmask encoded subsets for detection.
            @note
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0_real64)` and represents the number of chunks
            @endnote

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis::detect_subfunctionalization`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ancestor = np.ascontiguousarray(ancestor, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ancestor' must be an array of np.float64: {error}") from None
    if ancestor.ndim != 1:
        raise ValueError(f"'ancestor' must have 1 dimension, but has {ancestor.ndim}")
    try:
        genes = np.asfortranarray(genes, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'genes' must be an array of np.float64: {error}") from None
    if genes.ndim != 2:
        raise ValueError(f"'genes' must have 2 dimensions, but has {genes.ndim}")
    try:
        filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'filtered_paralogs_mask' must be an array of np.int32: {error}") from None
    if filtered_paralogs_mask.ndim != 1:
        raise ValueError(f"'filtered_paralogs_mask' must have 1 dimension, but has {filtered_paralogs_mask.ndim}")
    try:
        paralog_norms = np.ascontiguousarray(paralog_norms, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'paralog_norms' must be an array of np.float64: {error}") from None
    if paralog_norms.ndim != 1:
        raise ValueError(f"'paralog_norms' must have 1 dimension, but has {paralog_norms.ndim}")
    try:
        sorted_paralog_norms_perm = np.ascontiguousarray(sorted_paralog_norms_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'sorted_paralog_norms_perm' must be an array of np.int32: {error}") from None
    if sorted_paralog_norms_perm.ndim != 1:
        raise ValueError(f"'sorted_paralog_norms_perm' must have 1 dimension, but has {sorted_paralog_norms_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if paralog_norms.shape[0] != n_genes:
        raise ValueError(f"'paralog_norms' has {paralog_norms.shape[0]} along axis 0, but "
            f"'genes' implies n_genes == {n_genes}"
        )
    if sorted_paralog_norms_perm.shape[0] != n_genes:
        raise ValueError(f"'sorted_paralog_norms_perm' has {sorted_paralog_norms_perm.shape[0]} along axis 0, but "
            f"'genes' implies n_genes == {n_genes}"
        )
    if genes.shape[0] != n_dims:
        raise ValueError(f"'genes' has {genes.shape[0]} along axis 0, but "
            f"'ancestor' implies n_dims == {n_dims}"
        )

    # outputs and work arrays, which the caller never sees
    n_results = ctypes.c_int(0)
    work_arr_paralog_subsets = np.empty((n_mask_chunks, n_paralog_subsets,), dtype=np.int32, order='F')
    tmp_active_mask = np.empty((n_mask_chunks,), dtype=np.int32, order='C')
    tmp_paralog_vector = np.empty((n_dims,), dtype=np.float64, order='C')
    tmp_work_array = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.detect_subfunctionalization_c(
        ancestor,
        genes,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(ctypes.c_double(rdi_threshold)),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(n_results),
        ctypes.byref(ctypes.c_int(max_subset_size)),
        work_arr_paralog_subsets,
        ctypes.byref(ctypes.c_int(n_paralog_subsets)),
        tmp_active_mask,
        tmp_paralog_vector,
        paralog_norms,
        sorted_paralog_norms_perm,
        tmp_work_array,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETECT_SUBFUNCTIONALIZATION_ARGUMENTS)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets,
    }

def mask_chunk_count(
        n_genes,
):
    r"""Determines the needed chunk count for subset bit masks (an integer has only 32 bits)

    Parameters
    ----------
    n_genes : int
        number of genes

    Returns
    -------
    count : int
        number of 32 bit chunks a mask needs to encode `n_genes` genes

        Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0_real64)` and represents the number of chunks

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis::mask_chunk_count`.
    """
    # outputs and work arrays, which the caller never sees
    count = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.mask_chunk_count_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(count),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _MASK_CHUNK_COUNT_ARGUMENTS)

    return count.value

def filter_paralogs_by_pattern_subfunctionalization(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks,
):
    r"""Prefilters the genes for subfunctionalization

    Parameters
    ----------
    gene_angles : np.ndarray[np.float64] of shape (n_genes,)
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
    threshold : float
        filter threshold
    n_families : int
        number of families
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
    n_mask_chunks : int
        number of 32 bit chunks a mask needs to encode `n_genes` genes

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families,), column-major (order='F')
        bit mask that will have indices of genes kept by pattern set to 1, else 0

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis::filter_paralogs_by_pattern_subfunctionalization`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        gene_angles = np.ascontiguousarray(gene_angles, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_angles' must be an array of np.float64: {error}") from None
    if gene_angles.ndim != 1:
        raise ValueError(f"'gene_angles' must have 1 dimension, but has {gene_angles.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = gene_angles.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'gene_angles' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    masks = np.empty((n_mask_chunks, n_families,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    _lib.filter_paralogs_by_pattern_subfunctionalization_c(
        gene_angles,
        ctypes.byref(ctypes.c_double(threshold)),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        gene_to_fam,
        masks,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _FILTER_PARALOGS_BY_PATTERN_SUBFUNCTIONALIZATION_ARGUMENTS)

    return masks

def filter_paralogs_by_pattern_dosage_effect(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks,
):
    r"""Prefilters the genes for dosage effect

    Parameters
    ----------
    gene_angles : np.ndarray[np.float64] of shape (n_genes,)
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
    threshold : float
        filter threshold
    n_families : int
        number of families
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
    n_mask_chunks : int
        number of 32 bit chunks a mask needs to encode `n_genes` genes

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families,), column-major (order='F')
        bit mask that will have indices of genes kept by pattern set to 1, else 0

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis::filter_paralogs_by_pattern_dosage_effect`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        gene_angles = np.ascontiguousarray(gene_angles, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_angles' must be an array of np.float64: {error}") from None
    if gene_angles.ndim != 1:
        raise ValueError(f"'gene_angles' must have 1 dimension, but has {gene_angles.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = gene_angles.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'gene_angles' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    masks = np.empty((n_mask_chunks, n_families,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    _lib.filter_paralogs_by_pattern_dosage_effect_c(
        gene_angles,
        ctypes.byref(ctypes.c_double(threshold)),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        gene_to_fam,
        masks,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _FILTER_PARALOGS_BY_PATTERN_DOSAGE_EFFECT_ARGUMENTS)

    return masks

def calc_work_arr_paralog_subsets_size(
        max_subset_size,
        n_genes,
        filtered_paralogs_mask,
):
    r"""Calculates the needed size for the paralog-subsets work array

    Parameters
    ----------
    max_subset_size : int, modified in place
        maximum size that a subset must not exceed.
        @warning
        If the desired size is too large and leads to an integer overflow, `max_subset_size` will be set to the maximum valid size.

        Also, size will be set to number of genes in `filtered_paralogs_mask` if larger.
        @endwarning
    n_genes : int
        number of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,)
        Output mask with all genes disabled that did not pass the filter

    Returns
    -------
    dict
        with keys:

        max_subset_size : int
            maximum size that a subset must not exceed.
            @warning
            If the desired size is too large and leads to an integer overflow, `max_subset_size` will be set to the maximum valid size.

            Also, size will be set to number of genes in `filtered_paralogs_mask` if larger.
            @endwarning
        work_array_size : int
            The calculated needed work array size in absolute worst case scenario. Look into source for details.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis::calc_work_arr_paralog_subsets_size`.
    """
    # accept anything array-like, converting only when C needs it
    max_subset_size = ctypes.c_int(max_subset_size)
    try:
        filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'filtered_paralogs_mask' must be an array of np.int32: {error}") from None
    if filtered_paralogs_mask.ndim != 1:
        raise ValueError(f"'filtered_paralogs_mask' must have 1 dimension, but has {filtered_paralogs_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_mask_chunks = filtered_paralogs_mask.shape[0]

    # outputs and work arrays, which the caller never sees
    work_array_size = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_work_arr_paralog_subsets_size_c(
        ctypes.byref(max_subset_size),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(work_array_size),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_WORK_ARR_PARALOG_SUBSETS_SIZE_ARGUMENTS)

    return {
        "max_subset_size": max_subset_size.value,
        "work_array_size": work_array_size.value,
    }
