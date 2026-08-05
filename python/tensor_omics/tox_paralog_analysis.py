"""tox_paralog_analysis

Generated from the kernel; do not edit -- regenerate instead.

Python binding, generated from tox_paralog_analysis. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

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

_lib.detect_dosage_effect_expert_c.restype = None
_lib.detect_dosage_effect_expert_c.argtypes = (
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
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_DOSAGE_EFFECT_EXPERT_ARGUMENTS = ("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "tmp_active_mask", "tmp_paralog_vector", "max_angle", "gain_gamma", "ierr",)

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
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_DOSAGE_EFFECT_ARGUMENTS = ("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "max_angle", "gain_gamma", "ierr",)

_lib.detect_subfunctionalization_expert_c.restype = None
_lib.detect_subfunctionalization_expert_c.argtypes = (
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
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_SUBFUNCTIONALIZATION_EXPERT_ARGUMENTS = ("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "tmp_active_mask", "tmp_paralog_vector", "rdi_threshold", "paralog_norms", "sorted_paralog_norms_perm", "tmp_work_array", "ierr",)

_lib.detect_subfunctionalization_c.restype = None
_lib.detect_subfunctionalization_c.argtypes = (
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
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_SUBFUNCTIONALIZATION_ARGUMENTS = ("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "rdi_threshold", "paralog_norms", "sorted_paralog_norms_perm", "ierr",)

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
        The minimum valid value is `-1.0`.
        The maximum valid value is `1.0`.
    genes : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        RAP projected unit length expression vectors of genes
        The minimum valid value is `-1.0`.
        The maximum valid value is `1.0`.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    thresholds : np.ndarray[np.float64] of shape (n_axes,)
        threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis
        The minimum valid value is `-1.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    neofunc : np.ndarray[np.bool_] of shape (n_genes, n_axes,), column-major (order='F')
        `True` if neofunctionalization has been detected for the respective axes, always `False` for unassigned genes

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

def detect_dosage_effect_expert(
        ancestor,
        genes,
        filtered_paralogs_mask,
        max_subset_size,
        max_angle=3.141592653589793,
        gain_gamma=0.1,
):
    r"""Identifies subsets of paralogs matching this pattern

    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,)
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes,), column-major (order='F')
        expression vectors of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,)
        bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
    max_subset_size : int
        maximum subset size of checked gene subsets. Too large a value is capped to the
        maximum valid size. The bindings cap it automatically while sizing the work
        array; a Fortran caller caps it by calling
        [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size(subroutine)]] first.
        Zero is in range and means there is no subset to check -- the sizing routine reports
        it whenever the filtered families hold a single gene each. It reports a work array
        of zero slots along with it, which this routine does not accept, so a caller that
        gets zero back has nothing to detect and should not call here at all.
        The minimum valid value is `0`.
    max_angle : float, optional, default 3.141592653589793
        maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
        The default value is `4.0*atan(1.0)`.
        The minimum valid value is `0.0`.
        The maximum valid value is `PI`.
    gain_gamma : float, optional, default 0.1
        positive magnitude gain for dosage effect
        The default value is `0.1`.
        The minimum valid value is `above(0.0)`.

    Returns
    -------
    dict
        with keys:

        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_paralog_subsets,), column-major (order='F')
            working array to hold bitmask encoded subsets for detection.
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0)` and represents the number of chunks

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

    # work out what other procedures must supply, per DM_OUTPUT_FROM
    from .tox_paralog_analysis_kernel import calc_work_arr_paralog_subsets_size
    _calc_work_arr_paralog_subsets_size_result = calc_work_arr_paralog_subsets_size(max_subset_size=max_subset_size, n_genes=n_genes, filtered_paralogs_mask=filtered_paralogs_mask)
    max_subset_size = _calc_work_arr_paralog_subsets_size_result["max_subset_size"]
    n_paralog_subsets = _calc_work_arr_paralog_subsets_size_result["work_array_size"]

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

    _lib.detect_dosage_effect_expert_c(
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
        ctypes.byref(ctypes.c_double(max_angle)),
        ctypes.byref(ctypes.c_double(gain_gamma)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETECT_DOSAGE_EFFECT_EXPERT_ARGUMENTS)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets,
    }

def detect_dosage_effect(
        ancestor,
        genes,
        filtered_paralogs_mask,
        max_subset_size,
        max_angle=3.141592653589793,
        gain_gamma=0.1,
):
    r"""Identifies subsets of paralogs matching this pattern

    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,)
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes,), column-major (order='F')
        expression vectors of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,)
        bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
    max_subset_size : int
        maximum subset size of checked gene subsets. Too large a value is capped to the
        maximum valid size. The bindings cap it automatically while sizing the work
        array; a Fortran caller caps it by calling
        [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size(subroutine)]] first.
        Zero is in range and means there is no subset to check -- the sizing routine reports
        it whenever the filtered families hold a single gene each. It reports a work array
        of zero slots along with it, which this routine does not accept, so a caller that
        gets zero back has nothing to detect and should not call here at all.
        The minimum valid value is `0`.
    max_angle : float, optional, default 3.141592653589793
        maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
        The default value is `4.0*atan(1.0)`.
        The minimum valid value is `0.0`.
        The maximum valid value is `PI`.
    gain_gamma : float, optional, default 0.1
        positive magnitude gain for dosage effect
        The default value is `0.1`.
        The minimum valid value is `above(0.0)`.

    Returns
    -------
    dict
        with keys:

        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_paralog_subsets,), column-major (order='F')
            working array to hold bitmask encoded subsets for detection.
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0)` and represents the number of chunks

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis::detect_dosage_effect_alloc`.
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

    # work out what other procedures must supply, per DM_OUTPUT_FROM
    from .tox_paralog_analysis_kernel import calc_work_arr_paralog_subsets_size
    _calc_work_arr_paralog_subsets_size_result = calc_work_arr_paralog_subsets_size(max_subset_size=max_subset_size, n_genes=n_genes, filtered_paralogs_mask=filtered_paralogs_mask)
    max_subset_size = _calc_work_arr_paralog_subsets_size_result["max_subset_size"]
    n_paralog_subsets = _calc_work_arr_paralog_subsets_size_result["work_array_size"]

    # Fortran cannot check that shared extents agree; this can
    if genes.shape[0] != n_dims:
        raise ValueError(f"'genes' has {genes.shape[0]} along axis 0, but "
            f"'ancestor' implies n_dims == {n_dims}"
        )

    # outputs and work arrays, which the caller never sees
    n_results = ctypes.c_int(0)
    work_arr_paralog_subsets = np.empty((n_mask_chunks, n_paralog_subsets,), dtype=np.int32, order='F')
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
        ctypes.byref(ctypes.c_double(max_angle)),
        ctypes.byref(ctypes.c_double(gain_gamma)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETECT_DOSAGE_EFFECT_ARGUMENTS)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets,
    }

def detect_subfunctionalization_expert(
        ancestor,
        genes,
        filtered_paralogs_mask,
        max_subset_size,
        rdi_threshold,
        paralog_norms,
        sorted_paralog_norms_perm,
):
    r"""Identifies subsets of paralogs matching this pattern

    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,)
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes,), column-major (order='F')
        expression vectors of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,)
        bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
    max_subset_size : int
        maximum subset size of checked gene subsets. Too large a value is capped to the
        maximum valid size. The bindings cap it automatically while sizing the work
        array; a Fortran caller caps it by calling
        [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size(subroutine)]] first.
        Zero is in range and means there is no subset to check -- the sizing routine reports
        it whenever the filtered families hold a single gene each. It reports a work array
        of zero slots along with it, which this routine does not accept, so a caller that
        gets zero back has nothing to detect and should not call here at all.
        The minimum valid value is `0`.
    rdi_threshold : float
        max allowed residual distance from `ancestor`
        The minimum valid value is `0.0`.
    paralog_norms : np.ndarray[np.float64] of shape (n_genes,)
        euclidean norms of the genes, used for subset pruning (`norm` from `f42_utils` computes them)
        The minimum valid value is `0.0`.
    sorted_paralog_norms_perm : np.ndarray[np.int32] of shape (n_genes,)
        ascending permutation of the norms, for subset pruning: the smallest norm among the genes that could extend a subset must not fall below the subset's angle to the ancestor
        The minimum valid value is `1`.
        The maximum valid value is `n_genes`.

    Returns
    -------
    dict
        with keys:

        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_paralog_subsets,), column-major (order='F')
            working array to hold bitmask encoded subsets for detection.
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0)` and represents the number of chunks

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

    # work out what other procedures must supply, per DM_OUTPUT_FROM
    from .tox_paralog_analysis_kernel import calc_work_arr_paralog_subsets_size
    _calc_work_arr_paralog_subsets_size_result = calc_work_arr_paralog_subsets_size(max_subset_size=max_subset_size, n_genes=n_genes, filtered_paralogs_mask=filtered_paralogs_mask)
    max_subset_size = _calc_work_arr_paralog_subsets_size_result["max_subset_size"]
    n_paralog_subsets = _calc_work_arr_paralog_subsets_size_result["work_array_size"]

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

    _lib.detect_subfunctionalization_expert_c(
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
        ctypes.byref(ctypes.c_double(rdi_threshold)),
        paralog_norms,
        sorted_paralog_norms_perm,
        tmp_work_array,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETECT_SUBFUNCTIONALIZATION_EXPERT_ARGUMENTS)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets,
    }

def detect_subfunctionalization(
        ancestor,
        genes,
        filtered_paralogs_mask,
        max_subset_size,
        rdi_threshold,
        paralog_norms,
        sorted_paralog_norms_perm,
):
    r"""Identifies subsets of paralogs matching this pattern

    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,)
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes,), column-major (order='F')
        expression vectors of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,)
        bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
    max_subset_size : int
        maximum subset size of checked gene subsets. Too large a value is capped to the
        maximum valid size. The bindings cap it automatically while sizing the work
        array; a Fortran caller caps it by calling
        [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size(subroutine)]] first.
        Zero is in range and means there is no subset to check -- the sizing routine reports
        it whenever the filtered families hold a single gene each. It reports a work array
        of zero slots along with it, which this routine does not accept, so a caller that
        gets zero back has nothing to detect and should not call here at all.
        The minimum valid value is `0`.
    rdi_threshold : float
        max allowed residual distance from `ancestor`
        The minimum valid value is `0.0`.
    paralog_norms : np.ndarray[np.float64] of shape (n_genes,)
        euclidean norms of the genes, used for subset pruning (`norm` from `f42_utils` computes them)
        The minimum valid value is `0.0`.
    sorted_paralog_norms_perm : np.ndarray[np.int32] of shape (n_genes,)
        ascending permutation of the norms, for subset pruning: the smallest norm among the genes that could extend a subset must not fall below the subset's angle to the ancestor
        The minimum valid value is `1`.
        The maximum valid value is `n_genes`.

    Returns
    -------
    dict
        with keys:

        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_paralog_subsets,), column-major (order='F')
            working array to hold bitmask encoded subsets for detection.
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0)` and represents the number of chunks

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis::detect_subfunctionalization_alloc`.
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

    # work out what other procedures must supply, per DM_OUTPUT_FROM
    from .tox_paralog_analysis_kernel import calc_work_arr_paralog_subsets_size
    _calc_work_arr_paralog_subsets_size_result = calc_work_arr_paralog_subsets_size(max_subset_size=max_subset_size, n_genes=n_genes, filtered_paralogs_mask=filtered_paralogs_mask)
    max_subset_size = _calc_work_arr_paralog_subsets_size_result["max_subset_size"]
    n_paralog_subsets = _calc_work_arr_paralog_subsets_size_result["work_array_size"]

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
    ierr = ctypes.c_int(0)

    _lib.detect_subfunctionalization_c(
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
        ctypes.byref(ctypes.c_double(rdi_threshold)),
        paralog_norms,
        sorted_paralog_norms_perm,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETECT_SUBFUNCTIONALIZATION_ARGUMENTS)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets,
    }

def filter_paralogs_by_pattern_dosage_effect(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks,
):
    r"""Prefilters the genes for a pattern, so genes that cannot match it are not tried as subset extensions

    Parameters
    ----------
    gene_angles : np.ndarray[np.float64] of shape (n_genes,)
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
        The minimum valid value is `0.0`.
        The maximum valid value is `PI`.
    threshold : float
        filter threshold
    n_families : int
        number of families
        The minimum valid value is `1`.
        The maximum valid value is `n_genes`.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
    n_mask_chunks : int
        number of 32 bit chunks a mask needs to encode `n_genes` genes
        The minimum valid value is `(n_genes + 31) / 32`.

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families,), column-major (order='F')
        bit mask that will have the indices of genes kept by this pattern set to 1, else 0

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

def filter_paralogs_by_pattern_subfunctionalization(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks,
):
    r"""Prefilters the genes for a pattern, so genes that cannot match it are not tried as subset extensions

    Parameters
    ----------
    gene_angles : np.ndarray[np.float64] of shape (n_genes,)
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
        The minimum valid value is `0.0`.
        The maximum valid value is `PI`.
    threshold : float
        filter threshold
    n_families : int
        number of families
        The minimum valid value is `1`.
        The maximum valid value is `n_genes`.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
    n_mask_chunks : int
        number of 32 bit chunks a mask needs to encode `n_genes` genes
        The minimum valid value is `(n_genes + 31) / 32`.

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families,), column-major (order='F')
        bit mask that will have the indices of genes kept by this pattern set to 1, else 0

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
