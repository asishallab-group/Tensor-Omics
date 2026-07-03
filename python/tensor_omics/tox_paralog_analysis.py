from .error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


def mask_get_first_successor_idx(
        bit_mask
        ):
    """
    Parameters
    ----------
    bit_mask : np.ndarray[np.int32] of shape (n_bit_mask_elements,) in column-major layout (order='F')
        chunked mask to mark active genes

    Returns
    -------
    idx : int
        index of last active gene

    Notes
    -----
    Helper function that returns the index after the last active gene in `bit_mask`, so the first succeeding gene.
    """

    # ensure all array inputs are numpy arrays
    bit_mask = np.ascontiguousarray(bit_mask, dtype=np.int32)

    # extract dimension arguments
    n_bit_mask_elements = bit_mask.shape[0]


    # Create temporaries and/or outputs
    idx = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.mask_get_first_successor_idx_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.mask_get_first_successor_idx_c.restype = None

    tox.mask_get_first_successor_idx_c(
        bit_mask,
        ctypes.byref(ctypes.c_int(n_bit_mask_elements)),
        ctypes.byref(idx),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return idx.value


def mask_check_state(
        bit_mask,
        i_gene
        ):
    """
    Parameters
    ----------
    bit_mask : np.ndarray[np.int32] of shape (n_bit_mask_elements,) in column-major layout (order='F')
        chunked mask to mark active paralogs
    i_gene : int
        index of paralog to be marked active

    Returns
    -------
    state : bool
        check result

    Notes
    -----
    Checks the state of a bit/paralog in `bit_mask` -> .true. if 1 else .false.
    """

    # ensure all array inputs are numpy arrays
    bit_mask = np.ascontiguousarray(bit_mask, dtype=np.int32)

    # extract dimension arguments
    n_bit_mask_elements = bit_mask.shape[0]


    # Create temporaries and/or outputs
    i_gene = ctypes.c_int(i_gene)
    state = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.mask_check_state_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.mask_check_state_c.restype = None

    tox.mask_check_state_c(
        bit_mask,
        ctypes.byref(ctypes.c_int(n_bit_mask_elements)),
        ctypes.byref(i_gene),
        ctypes.byref(state),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return state.value


def detect_neofunctionalization(
        ancestors,
        genes,
        gene_to_fam,
        thresholds
        ):
    """
    Parameters
    ----------
    ancestors : np.ndarray[np.float64] of shape (n_axes, n_families) in column-major layout (order='F')
        RAP projected unit length expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_axes, n_genes) in column-major layout (order='F')
        RAP projected unit length expression vectors of genes
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        mapping of gene index to family index
    thresholds : np.ndarray[np.float64] of shape (n_axes,) in column-major layout (order='F')
        threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis

    Returns
    -------
    neofunc : np.ndarray[np.int32] of shape (n_genes, n_axes) in column-major layout (order='F')
        `.true.` if neofunctionalization has been detected for the respective axes

    Notes
    -----
    Identifies neofunctionalization for genes by checking whether the difference of expression to its ancestor exceeds the threshold for the respective axis.
    """

    # ensure all array inputs are numpy arrays
    ancestors = np.asfortranarray(ancestors, dtype=np.float64)
    genes = np.asfortranarray(genes, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    thresholds = np.ascontiguousarray(thresholds, dtype=np.float64)

    # extract dimension arguments
    n_families = ancestors.shape[1]
    n_axes = ancestors.shape[0]
    n_genes = genes.shape[1]


    # Create temporaries and/or outputs
    neofunc = np.empty((n_genes, n_axes), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.detect_neofunctionalization_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.detect_neofunctionalization_c.restype = None

    tox.detect_neofunctionalization_c(
        ancestors,
        ctypes.byref(ctypes.c_int(n_families)),
        genes,
        ctypes.byref(ctypes.c_int(n_axes)),
        gene_to_fam,
        ctypes.byref(ctypes.c_int(n_genes)),
        thresholds,
        neofunc,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    neofunc.setflags(write=False)

    return neofunc


def detect_dosage_effect(
        ancestor,
        genes,
        filtered_paralogs_mask,
        max_subset_size,
        n_paralog_subsets,
        max_angle=3.141592653589793,
        gain_gamma=0.1
        ):
    """
    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes) in column-major layout (order='F')
        expression vectors of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
        bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculationDM_FROM(masks(:, family_idx), filter_paralogs_by_pattern, tox_paralog_analysis, JUST_INFO)
    max_subset_size : int
        maximum subset size of checked gene subsets.DM_FROM(max_subset_size, calc_work_arr_paralog_subsets_size, tox_paralog_analysis, AUTO)
    max_angle : float, optional
        in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise prunedThe default value is `PI`.
    gain_gamma : float, optional
        positive magnitude gain for dosage effectThe default value is `0.1_real64`.

    Returns
    -------
    results : dict
        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`,
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_results) in column-major layout (order='F')
            working array to hold bitmask encoded subsets for detection.@noteEach bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks@endnoteThe first `n_results` elements will hold the results.,
        active_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
            working array to hold the extended subsets,
        temp_paralog_vector : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
            vector used for pruning subsets


    Notes
    -----
    Identifies subsets of paralogs with small angle to the `ancestor` (max_angle) and sum to a magnitude significantly exceeding `norm(ancestor)` (gain)
    """

    # ensure all array inputs are numpy arrays
    ancestor = np.ascontiguousarray(ancestor, dtype=np.float64)
    genes = np.asfortranarray(genes, dtype=np.float64)
    filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)

    # extract dimension arguments
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]


    # Create temporaries and/or outputs
    n_results = ctypes.c_int(0)
    max_subset_size = ctypes.c_int(max_subset_size)
    work_arr_paralog_subsets = np.empty((n_mask_chunks, n_paralog_subsets), dtype=np.int32, order='F')
    active_mask = np.empty((n_mask_chunks,), dtype=np.int32, order='F')
    temp_paralog_vector = np.empty((n_dims,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)
    max_angle = ctypes.c_double(max_angle)
    gain_gamma = ctypes.c_double(gain_gamma)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.detect_dosage_effect_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double)
    )
    tox.detect_dosage_effect_c.restype = None

    tox.detect_dosage_effect_c(
        ancestor,
        genes,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_dims)),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(n_results),
        ctypes.byref(max_subset_size),
        work_arr_paralog_subsets,
        ctypes.byref(ctypes.c_int(n_paralog_subsets)),
        active_mask,
        temp_paralog_vector,
        ctypes.byref(ierr),
        ctypes.byref(max_angle),
        ctypes.byref(gain_gamma)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    work_arr_paralog_subsets.setflags(write=False)
    active_mask.setflags(write=False)
    temp_paralog_vector.setflags(write=False)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets[..., :n_results.value].copy(),
        "active_mask": active_mask,
        "temp_paralog_vector": temp_paralog_vector
    }


def detect_subfunctionalization(
        ancestor,
        genes,
        rdi_threshold,
        filtered_paralogs_mask,
        max_subset_size,
        n_paralog_subsets,
        paralog_norms,
        sorted_paralog_norms_perm
        ):
    """
    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes) in column-major layout (order='F')
        expression vectors of genes
    rdi_threshold : float
        max allowed residual distance from `ancestor`
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
        bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    max_subset_size : int
        maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    paralog_norms : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` function from `f42_utils` function for this)
    sorted_paralog_norms_perm : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor

    Returns
    -------
    results : dict
        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`,
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_results) in column-major layout (order='F')
            working array to hold bitmask encoded subsets for detection.@noteEach bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks@endnoteThe first `n_results` elements will hold the results.,
        active_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
            working array to hold the extended subsets,
        temp_paralog_vector : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
            vector used for pruning subsets,
        temp_work_array : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
            needed for efficient check of minimum value after a certain index


    Notes
    -----
    Identifies subsets of paralogs exhibiting significant angles to the `ancestor`
    """

    # ensure all array inputs are numpy arrays
    ancestor = np.ascontiguousarray(ancestor, dtype=np.float64)
    genes = np.asfortranarray(genes, dtype=np.float64)
    filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)
    paralog_norms = np.ascontiguousarray(paralog_norms, dtype=np.float64)
    sorted_paralog_norms_perm = np.ascontiguousarray(sorted_paralog_norms_perm, dtype=np.int32)

    # extract dimension arguments
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]


    # Create temporaries and/or outputs
    rdi_threshold = ctypes.c_double(rdi_threshold)
    n_results = ctypes.c_int(0)
    max_subset_size = ctypes.c_int(max_subset_size)
    work_arr_paralog_subsets = np.empty((n_mask_chunks, n_paralog_subsets), dtype=np.int32, order='F')
    active_mask = np.empty((n_mask_chunks,), dtype=np.int32, order='F')
    temp_paralog_vector = np.empty((n_dims,), dtype=np.float64, order='F')
    temp_work_array = np.empty((n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.detect_subfunctionalization_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.detect_subfunctionalization_c.restype = None

    tox.detect_subfunctionalization_c(
        ancestor,
        genes,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(rdi_threshold),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(n_results),
        ctypes.byref(max_subset_size),
        work_arr_paralog_subsets,
        ctypes.byref(ctypes.c_int(n_paralog_subsets)),
        active_mask,
        temp_paralog_vector,
        paralog_norms,
        sorted_paralog_norms_perm,
        temp_work_array,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    work_arr_paralog_subsets.setflags(write=False)
    active_mask.setflags(write=False)
    temp_paralog_vector.setflags(write=False)
    temp_work_array.setflags(write=False)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets[..., :n_results.value].copy(),
        "active_mask": active_mask,
        "temp_paralog_vector": temp_paralog_vector,
        "temp_work_array": temp_work_array
    }


def detect_patterns(
        ancestor,
        genes,
        pattern_mode,
        filtered_paralogs_mask,
        max_subset_size,
        n_paralog_subsets,
        dosage_max_angle=3.141592653589793,
        dosage_gain_gamma=0.1,
        subfunc_rdi_threshold=None,
        subfunc_paralog_norms=None,
        subfunc_sorted_paralog_norms_perm=None
        ):
    """
    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes) in column-major layout (order='F')
        expression vectors of genes
    pattern_mode : int
        used pattern for detection|         Mode         |                              Value                              |
        |----------------------|-----------------------------------------------------------------|
        |    Dosage Effect     |  [[tox_paralog_analysis(module):MODE_DOSAGE_PATTERN(variable)]] |
        | Subfunctionalization | [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]] |
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
        bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    max_subset_size : int
        maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    dosage_max_angle : float, optional
        in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise prunedThe default value is `PI`.
    dosage_gain_gamma : float, optional
        in dosage mode required positive magnitude gain for dosageThe default value is `0.1_real64`.
    subfunc_rdi_threshold : float, optional
        max allowed residual distance from `ancestor`This optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].
    subfunc_paralog_norms : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F'), optional
        in subfunctionalization mode needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` from `f42_utils` function for this)This optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].
    subfunc_sorted_paralog_norms_perm : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F'), optional
        in subfunctionalization mode needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestorThis optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].

    Returns
    -------
    results : dict
        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`,
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_results) in column-major layout (order='F')
            working array to hold bitmask encoded subsets for detection.@noteEach bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks@endnoteThe first `n_results` elements will hold the results.,
        active_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
            working array to hold the extended subsets,
        temp_paralog_vector : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
            vector used for pruning subsets,
        subfunc_temp_work_array : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
            in subfunctionalization mode needed for efficient check of minimum value after a certain indexThis optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].


    Notes
    -----
    Identifies subsets of paralogs where dosage effect or subfunctionalization applies, depending on `pattern`
    """

    # ensure all array inputs are numpy arrays
    ancestor = np.ascontiguousarray(ancestor, dtype=np.float64)
    genes = np.asfortranarray(genes, dtype=np.float64)
    filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)
    if subfunc_paralog_norms is not None:
        subfunc_paralog_norms = np.ascontiguousarray(subfunc_paralog_norms, dtype=np.float64)

    if subfunc_sorted_paralog_norms_perm is not None:
        subfunc_sorted_paralog_norms_perm = np.ascontiguousarray(subfunc_sorted_paralog_norms_perm, dtype=np.int32)


    # extract dimension arguments
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]


    # Create temporaries and/or outputs
    pattern_mode = ctypes.c_int(pattern_mode)
    n_results = ctypes.c_int(0)
    max_subset_size = ctypes.c_int(max_subset_size)
    work_arr_paralog_subsets = np.empty((n_mask_chunks, n_paralog_subsets), dtype=np.int32, order='F')
    active_mask = np.empty((n_mask_chunks,), dtype=np.int32, order='F')
    temp_paralog_vector = np.empty((n_dims,), dtype=np.float64, order='F')
    dosage_max_angle = ctypes.c_double(dosage_max_angle)
    dosage_gain_gamma = ctypes.c_double(dosage_gain_gamma)
    subfunc_rdi_threshold = ctypes.c_double(subfunc_rdi_threshold)
    subfunc_temp_work_array = np.empty((n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.detect_patterns_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        nullable(ctypes.POINTER(ctypes.c_double)),
        nullable(np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64)),
        nullable(np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32)),
        nullable(np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64)),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.detect_patterns_c.restype = None

    tox.detect_patterns_c(
        ancestor,
        genes,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(pattern_mode),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(n_results),
        ctypes.byref(max_subset_size),
        work_arr_paralog_subsets,
        ctypes.byref(ctypes.c_int(n_paralog_subsets)),
        active_mask,
        temp_paralog_vector,
        ctypes.byref(dosage_max_angle),
        ctypes.byref(dosage_gain_gamma),
        ctypes.byref(subfunc_rdi_threshold),
        subfunc_paralog_norms,
        subfunc_sorted_paralog_norms_perm,
        subfunc_temp_work_array,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    work_arr_paralog_subsets.setflags(write=False)
    active_mask.setflags(write=False)
    temp_paralog_vector.setflags(write=False)
    subfunc_temp_work_array.setflags(write=False)

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets[..., :n_results.value].copy(),
        "active_mask": active_mask,
        "temp_paralog_vector": temp_paralog_vector,
        "subfunc_temp_work_array": subfunc_temp_work_array
    }


def mask_chunk_count(
        n_genes
        ):
    """
    Parameters
    ----------
    n_genes : int
        number of genes

    Returns
    -------
    count : int
        number of 32 bit chunks a mask needs to encode `n_genes` genesEach bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks

    Notes
    -----
    This subroutine easily determines the needed chunk count for subset bit masks, as an integer has only 32 bits.
    """

    # ensure all array inputs are numpy arrays


    # extract dimension arguments



    # Create temporaries and/or outputs
    n_genes = ctypes.c_int(n_genes)
    count = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.mask_chunk_count_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.mask_chunk_count_c.restype = None

    tox.mask_chunk_count_c(
        ctypes.byref(n_genes),
        ctypes.byref(count),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return count.value


def filter_paralogs_by_pattern_subfunctionalization(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks
        ):
    """
    Parameters
    ----------
    gene_angles : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
    threshold : float
        filter threshold
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families) in column-major layout (order='F')
        bit mask that will have indices of genes kept by pattern set to 1, else 0

    Notes
    -----
    This subroutine prefilters the genes for subfunctionalization,as genes that are already too close in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
    """

    # ensure all array inputs are numpy arrays
    gene_angles = np.ascontiguousarray(gene_angles, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)

    # extract dimension arguments
    n_genes = gene_angles.shape[0]


    # Create temporaries and/or outputs
    threshold = ctypes.c_double(threshold)
    masks = np.empty((n_mask_chunks, n_families), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.filter_paralogs_by_pattern_subfunctionalization_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.filter_paralogs_by_pattern_subfunctionalization_c.restype = None

    tox.filter_paralogs_by_pattern_subfunctionalization_c(
        gene_angles,
        ctypes.byref(threshold),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        gene_to_fam,
        masks,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    masks.setflags(write=False)

    return masks


def filter_paralogs_by_pattern_dosage_effect(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks
        ):
    """
    Parameters
    ----------
    gene_angles : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
    threshold : float
        filter threshold
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families) in column-major layout (order='F')
        bit mask that will have indices of genes kept by pattern set to 1, else 0

    Notes
    -----
    This subroutine prefilters the genes for dosage effect,as genes that are already too distant in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
    """

    # ensure all array inputs are numpy arrays
    gene_angles = np.ascontiguousarray(gene_angles, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)

    # extract dimension arguments
    n_genes = gene_angles.shape[0]


    # Create temporaries and/or outputs
    threshold = ctypes.c_double(threshold)
    masks = np.empty((n_mask_chunks, n_families), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.filter_paralogs_by_pattern_dosage_effect_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.filter_paralogs_by_pattern_dosage_effect_c.restype = None

    tox.filter_paralogs_by_pattern_dosage_effect_c(
        gene_angles,
        ctypes.byref(threshold),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        gene_to_fam,
        masks,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    masks.setflags(write=False)

    return masks


def filter_paralogs_by_pattern(
        pattern_mode,
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks
        ):
    """
    Parameters
    ----------
    pattern_mode : int
        used pattern for detection|         Mode         |                              Value                              |
        |----------------------|-----------------------------------------------------------------|
        |    Dosage Effect     |  [[tox_paralog_analysis(module):MODE_DOSAGE_PATTERN(variable)]] |
        | Subfunctionalization | [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]] |
    gene_angles : np.ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
    threshold : float
        filter threshold
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families) in column-major layout (order='F')
        bit mask that will have indices of genes kept by pattern set to 1, else 0

    Notes
    -----
    This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
    """

    # ensure all array inputs are numpy arrays
    gene_angles = np.ascontiguousarray(gene_angles, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)

    # extract dimension arguments
    n_genes = gene_angles.shape[0]


    # Create temporaries and/or outputs
    pattern_mode = ctypes.c_int(pattern_mode)
    threshold = ctypes.c_double(threshold)
    masks = np.empty((n_mask_chunks, n_families), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.filter_paralogs_by_pattern_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.filter_paralogs_by_pattern_c.restype = None

    tox.filter_paralogs_by_pattern_c(
        ctypes.byref(pattern_mode),
        gene_angles,
        ctypes.byref(threshold),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        gene_to_fam,
        masks,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    masks.setflags(write=False)

    return masks


def calc_work_arr_paralog_subsets_size(
        max_subset_size,
        n_genes,
        filtered_paralogs_mask
        ):
    """
    Parameters
    ----------
    max_subset_size : int, modified in-place
        maximum size that a subset must not exceed.@warningIf the desired size is too large and leads to an integer overflow, `max_subset_size` will be set to the maximum valid size.Also, size will be set to number of genes in `filtered_paralogs_mask` if larger.@endwarning
    n_genes : int
        number of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
        Output mask with all genes disabled that did not pass the filterDM_FROM(masks(:, family_idx), filter_paralogs_by_pattern, tox_paralog_analysis, JUST_INFO)

    Returns
    -------
    work_array_size : int
        The calculated needed work array size in absolute worst case scenario. Look into source for details.

    Notes
    -----
    The `detect_*` subroutines need a work array for the to be tested subsets.In worst case, all need to be tried and subsets that cannot be extended will be kept as results.This is the reason why the work array holds the results as well, as all subsets that are stored in the array can be results as well.This subroutine calculates the needed size for the work array.
    """

    # ensure all array inputs are numpy arrays
    filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)

    # extract dimension arguments
    n_mask_chunks = filtered_paralogs_mask.shape[0]


    # Create temporaries and/or outputs
    max_subset_size = ctypes.c_int(max_subset_size)
    n_genes = ctypes.c_int(n_genes)
    work_array_size = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.calc_work_arr_paralog_subsets_size_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.calc_work_arr_paralog_subsets_size_c.restype = None

    tox.calc_work_arr_paralog_subsets_size_c(
        ctypes.byref(max_subset_size),
        ctypes.byref(n_genes),
        ctypes.byref(work_array_size),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return work_array_size.value


def mask_set_state(
        n_bit_mask_elements,
        i_gene,
        state
        ):
    """
    Parameters
    ----------
    i_gene : int
        index of paralog to be marked active
    state : bool
        state the bit should be set to

    Returns
    -------
    bit_mask : np.ndarray[np.int32] of shape (n_bit_mask_elements,) in column-major layout (order='F')
        chunked mask to mark active paralogs

    Notes
    -----
    Sets the state of a bit/gene in `bit_mask`
    """

    # ensure all array inputs are numpy arrays


    # extract dimension arguments



    # Create temporaries and/or outputs
    bit_mask = np.empty((n_bit_mask_elements,), dtype=np.int32, order='F')
    i_gene = ctypes.c_int(i_gene)
    state = ctypes.c_int(state)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.mask_set_state_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.mask_set_state_c.restype = None

    tox.mask_set_state_c(
        bit_mask,
        ctypes.byref(ctypes.c_int(n_bit_mask_elements)),
        ctypes.byref(i_gene),
        ctypes.byref(state),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    bit_mask.setflags(write=False)

    return bit_mask
