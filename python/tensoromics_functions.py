"""
TensorOmics Functions Module
Python wrapper functions for Fortran routines via C interface
"""

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL(dll_path)


def tox_errors(ierr):
    """
    Check error code and raise informative exception if needed
    
    Args:
        ierr: Error code from Fortran routine
        context: Optional context string for debugging
    
    Raises:
        RuntimeError: If error code indicates failure
    """
    error_messages = {
        0: None,  # Success
        101: "File could not be opened",
        102: "Could not read magic number",
        103: "Could not read array type code", 
        104: "Could not read array dimension number",
        105: "Could not read array dimensions",
        106: "Could not read character length",
        107: "Could not read array data",
        200: "Invalid file format (magic number mismatch)",
        202: "No axes selected (empty input)",
        5002: "File not open or unit not connected",
        9999: "Unknown error"
    }
    
    msg = error_messages.get(ierr, f"Unknown Fortran error code: {ierr}")
    
    if msg is not None:
        raise RuntimeError(msg)

def _readonly(*arrays: np.ndarray) -> None:
    """Mark all given NumPy arrays as read-only."""
    for a in arrays:
        if isinstance(a, np.ndarray):
            a.flags.writeable = False
            # NOTE: Returned NumPy arrays are read-only for safety.
            # If you need to modify them (e.g., for plotting), use `.copy()`.

def tox_normalize_by_std_dev(input_matrix):
    """
    Normalize gene expression values by standard deviation
    
    Args:
        input_matrix: A numeric matrix with genes as rows and tissues as columns
    
    Returns:
        numpy.ndarray: Normalized matrix with same dimensions as input
    """
    input_matrix = np.asarray(input_matrix, dtype=np.float64)
    n_genes, n_tissues = input_matrix.shape
    
    # Validate input data
    if np.any(np.isnan(input_matrix)):
        raise ValueError(f"Input matrix contains NaN values: {np.sum(np.isnan(input_matrix))}")
    if np.any(np.isinf(input_matrix)):
        raise ValueError(f"Input matrix contains infinite values: {np.sum(np.isinf(input_matrix))}")
    
    # Flatten input and prepare output
    input_flat = np.asfortranarray(input_matrix).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    
    # Setup C wrapper
    normalize_c = lib.normalize_by_std_dev_c
    normalize_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_tissues
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # input
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # output
    ]
    normalize_c.restype = None
    
    # Call Fortran routine
    normalize_c(n_genes, n_tissues, input_flat, output_flat)
    
    # Reshape and return
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    
    # Mark output as read-only
    _readonly(result)
    return result

def tox_quantile_normalization(input_matrix):
    """
    Quantile normalization of gene expression values
    
    Args:
        input_matrix: A numeric matrix with genes as rows and tissues as columns
    
    Returns:
        numpy.ndarray: Quantile-normalized matrix with same dimensions as input
    """
    input_matrix = np.asarray(input_matrix, dtype=np.float64)
    n_genes, n_tissues = input_matrix.shape
    
    # Flatten input and prepare output arrays
    input_flat = np.asfortranarray(input_matrix).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    temp_col = np.zeros(n_genes, dtype=np.float64)
    rank_means = np.zeros(n_genes, dtype=np.float64)
    perm = np.zeros(n_genes, dtype=np.int32)
    
    # Prepare stack arrays for sorting
    max_stack = max(2 * n_genes, 2)
    stack_left = np.zeros(max_stack, dtype=np.int32)
    stack_right = np.zeros(max_stack, dtype=np.int32)
    
    # Setup C wrapper
    quantile_norm_c = lib.quantile_normalization_c
    quantile_norm_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_tissues
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # input
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # output
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # temp_col
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # rank_means
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right
        ctypes.c_int,  # max_stack
    ]
    quantile_norm_c.restype = None
    
    # Call Fortran routine
    quantile_norm_c(n_genes, n_tissues, input_flat, output_flat,
                    temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    
    # Reshape and return
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    
    # Mark output as read-only
    _readonly(result)
    return result

def tox_log2_transformation(input_matrix):
    """
    Apply log2(x + 1) transformation to gene expression values
    
    Args:
        input_matrix: A numeric matrix with genes as rows and tissues as columns
    
    Returns:
        numpy.ndarray: Log2-transformed matrix with same dimensions as input
    """
    input_matrix = np.asarray(input_matrix, dtype=np.float64)
    n_genes, n_tissues = input_matrix.shape
    
    # Flatten input and prepare output
    input_flat = np.asfortranarray(input_matrix).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    
    # Setup C wrapper
    log2_transform_c = lib.log2_transformation_c
    log2_transform_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_tissues
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # input
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # output
    ]
    log2_transform_c.restype = None
    
    # Call Fortran routine
    log2_transform_c(n_genes, n_tissues, input_flat, output_flat)
    
    # Reshape and return
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    
    # Mark output as read-only
    _readonly(result)
    return result

def tox_calculate_tissue_averages(input_matrix, group_starts, group_counts):
    """
    Calculate average expression across replicates for each tissue group
    
    Args:
        input_matrix: A numeric matrix with genes as rows and tissue replicates as columns
        group_starts: Array of starting column indices for each group (1-based for Fortran)
        group_counts: Array of counts for each group
    
    Returns:
        numpy.ndarray: Matrix with genes as rows and averaged tissues as columns
    """
    input_matrix = np.asarray(input_matrix, dtype=np.float64)
    group_starts = np.asarray(group_starts, dtype=np.int32)
    group_counts = np.asarray(group_counts, dtype=np.int32)
    
    n_genes, n_samples = input_matrix.shape
    n_groups = len(group_starts)
    
    if len(group_counts) != n_groups:
        raise ValueError("group_starts and group_counts must have same length")
    
    # Flatten input and prepare output
    input_flat = np.asfortranarray(input_matrix).ravel(order='F')
    output_flat = np.zeros(n_genes * n_groups, dtype=np.float64)
    
    # Setup C wrapper
    tiss_avg_c = lib.calc_tiss_avg_c
    tiss_avg_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_groups
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # group_starts
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # group_counts
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # input
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # output
    ]
    tiss_avg_c.restype = None
    
    # Call Fortran routine
    tiss_avg_c(n_genes, n_groups, group_starts, group_counts, input_flat, output_flat)
    
    # Reshape and return
    result = output_flat.reshape((n_genes, n_groups), order='F')
    
    # Mark output as read-only
    _readonly(result)
    return result

def tox_calculate_fold_changes(input_matrix, control_cols, condition_cols):
    """
    Calculate log2 fold changes between control and condition columns
    
    Args:
        input_matrix: A numeric matrix with genes as rows and tissues/conditions as columns
        control_cols: Array of control column indices (1-based for Fortran)
        condition_cols: Array of condition column indices (1-based for Fortran)
    
    Returns:
        numpy.ndarray: Matrix with genes as rows and fold change values as columns
    """
    input_matrix = np.asarray(input_matrix, dtype=np.float64)
    control_cols = np.asarray(control_cols, dtype=np.int32)
    condition_cols = np.asarray(condition_cols, dtype=np.int32)
    
    n_genes, n_samples = input_matrix.shape
    n_pairs = len(control_cols)
    
    if len(condition_cols) != n_pairs:
        raise ValueError("control_cols and condition_cols must have same length")
    
    # Flatten input and prepare output
    input_flat = np.asfortranarray(input_matrix).ravel(order='F')
    output_flat = np.zeros(n_genes * n_pairs, dtype=np.float64)
    
    # Setup C wrapper
    fchange_c = lib.calc_fchange_c
    fchange_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_samples
        ctypes.c_int,  # n_pairs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # control_cols
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # condition_cols
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # input
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # output
    ]
    fchange_c.restype = None
    
    # Call Fortran routine
    fchange_c(n_genes, n_samples, n_pairs, control_cols, condition_cols, input_flat, output_flat)
    
    # Reshape and return
    result = output_flat.reshape((n_genes, n_pairs), order='F')
    
    # Mark output as read-only
    _readonly(result)
    return result

def tox_calculate_tissue_versatility(expression_vectors, vector_selection, axis_selection):
    """
    Calculate Tissue Versatility
    
    Computes normalized tissue versatility for selected expression vectors.
    The metric is based on the angle between each gene expression vector and the space diagonal.
    Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
    This function automatically checks for errors and throws informative exceptions.
    
    Args:
        expression_vectors: Matrix where each column is a gene expression vector (n_axes x n_vectors)
        vector_selection: Boolean or integer array indicating which vectors to process (length n_vectors)
        axis_selection: Boolean or integer array indicating which axes to include in calculation (length n_axes)
    
    Returns:
        dict: Dictionary containing:
            - tissue_versatilities: Normalized tissue versatility values [0,1] for selected vectors
            - tissue_angles_deg: Angles in degrees [0,90] for selected vectors
            - n_selected_vectors: Number of vectors processed
            - n_selected_axes: Number of axes used in calculation

    """
    # Input validation
    if not isinstance(expression_vectors, np.ndarray):
        raise ValueError("expression_vectors must be a numpy array")
    
    if expression_vectors.ndim != 2:
        raise ValueError("expression_vectors must be a 2D array")
    
    # Convert inputs to numpy arrays
    expression_vectors = np.asarray(expression_vectors, dtype=np.float64)
    vector_selection = np.asarray(vector_selection, dtype=bool)
    axis_selection = np.asarray(axis_selection, dtype=bool)
    
    # Get dimensions
    n_axes, n_vectors = expression_vectors.shape
    
    # Validate dimensions
    if len(vector_selection) != n_vectors:
        raise ValueError("vector_selection length must match number of columns in expression_vectors")
    if len(axis_selection) != n_axes:
        raise ValueError("axis_selection length must match number of rows in expression_vectors")
    
    # Calculate counts
    n_selected_vectors = int(np.sum(vector_selection))
    n_selected_axes = int(np.sum(axis_selection))
    
    # Ensure arrays have correct dtype and memory layout
    expr_f = np.asfortranarray(expression_vectors, dtype=np.float64)
    select_vec = np.ascontiguousarray(vector_selection.astype(np.int32))
    select_axes = np.ascontiguousarray(axis_selection.astype(np.int32))
    
    # Prepare output arrays
    tissue_versatilities = np.empty(n_selected_vectors, dtype=np.float64)
    tissue_angles_deg = np.empty(n_selected_vectors, dtype=np.float64)
    ierr = ctypes.c_int(0)
    
    # Setup C/Fortran wrapper with proper type annotations
    tv = lib.compute_tissue_versatility_c
    tv.argtypes = [
        ctypes.c_int,  # n_axes
        ctypes.c_int,  # n_vectors
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # expression_vectors (Fortran order)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # exp_vecs_selection_index
        ctypes.c_int,  # n_selected_vectors
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # axes_selection
        ctypes.c_int,  # n_selected_axes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # tissue_versatilities
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # tissue_angles_deg
        ctypes.POINTER(ctypes.c_int),  # ierr
    ]
    tv.restype = None
    
    # Call the C/Fortran wrapper
    tv(
        n_axes,
        n_vectors,
        expr_f,
        select_vec,
        n_selected_vectors,
        select_axes,
        n_selected_axes,
        tissue_versatilities,
        tissue_angles_deg,
        ctypes.byref(ierr)
    )
    
    # Check for errors and throw informative messages
    tox_errors(ierr.value)
    
    # Mark outputs as read-only
    _readonly(tissue_versatilities, tissue_angles_deg)
    
    # Return structured result (no ierr since we checked for errors)
    return {
        'tissue_versatilities': tissue_versatilities,
        'tissue_angles_deg': tissue_angles_deg
    }


def tox_euclidean_distance(vec1, vec2):
    """
    Calculate Euclidean distance between two vectors
    
    Args:
        vec1: First vector (numpy array)
        vec2: Second vector (numpy array)
    
    Returns:
        float: Euclidean distance between the vectors
    """
    # Input validation
    if not isinstance(vec1, np.ndarray):
        vec1 = np.asarray(vec1, dtype=np.float64)
    if not isinstance(vec2, np.ndarray):
        vec2 = np.asarray(vec2, dtype=np.float64)
    
    vec1 = np.asarray(vec1, dtype=np.float64)
    vec2 = np.asarray(vec2, dtype=np.float64)
    
    if len(vec1) != len(vec2):
        raise ValueError("Vectors must have the same length")
    if len(vec1) == 0:
        raise ValueError("Vectors cannot be empty")
    if not np.issubdtype(vec1.dtype, np.number) or not np.issubdtype(vec2.dtype, np.number):
        raise ValueError("Vectors must be numeric")
    
    # Ensure contiguous arrays
    vec1 = np.ascontiguousarray(vec1, dtype=np.float64)
    vec2 = np.ascontiguousarray(vec2, dtype=np.float64)
    
    # Prepare output
    result = ctypes.c_double(0.0)
    
    # Setup C wrapper
    euclidean_distance_c = lib.euclidean_distance_c
    euclidean_distance_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double)
    ]
    euclidean_distance_c.restype = None
    
    # Call Fortran routine
    euclidean_distance_c(vec1, vec2, len(vec1), ctypes.byref(result))
    
    return result.value


def tox_distance_to_centroid(genes, centroids, gene_to_fam, d):
    """
    Calculate distance from each gene to its family centroid
    
    Args:
        genes: Gene expression data as flat array (n_genes * d elements)
        centroids: Family centroids as flat array (n_families * d elements)
        gene_to_fam: Gene-to-family mapping (0 = no family, >0 = family index)
        d: Number of dimensions
    
    Returns:
        numpy.ndarray: Distances from each gene to its centroid (-1 for invalid families)
    """
    # Input validation
    if not isinstance(genes, np.ndarray):
        genes = np.asarray(genes, dtype=np.float64)
    if not isinstance(centroids, np.ndarray):
        centroids = np.asarray(centroids, dtype=np.float64)
    if not isinstance(gene_to_fam, np.ndarray):
        gene_to_fam = np.asarray(gene_to_fam, dtype=np.int32)
    
    genes = np.asarray(genes, dtype=np.float64)
    centroids = np.asarray(centroids, dtype=np.float64)
    gene_to_fam = np.asarray(gene_to_fam, dtype=np.int32)
    d = int(d)
    
    # Calculate dimensions
    n_genes = len(genes) // d
    n_families = len(centroids) // d
    
    # Validate dimensions
    if len(genes) % d != 0:
        raise ValueError("Length of genes must be divisible by d")
    if len(centroids) % d != 0:
        raise ValueError("Length of centroids must be divisible by d")
    if len(gene_to_fam) != n_genes:
        raise ValueError("Length of gene_to_fam must equal number of genes")
    if np.any(gene_to_fam < 0):
        raise ValueError("gene_to_fam indices must be between 0 and n_families (0 = no family assignment)")
    
    # Ensure contiguous arrays
    genes = np.ascontiguousarray(genes, dtype=np.float64)
    centroids = np.ascontiguousarray(centroids, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    
    # Prepare output array
    distances = np.zeros(n_genes, dtype=np.float64)
    
    # Setup C wrapper
    distance_to_centroid_c = lib.distance_to_centroid_c
    distance_to_centroid_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # genes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # centroids
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        ctypes.c_int   # d
    ]
    distance_to_centroid_c.restype = None
    
    # Call Fortran routine
    distance_to_centroid_c(n_genes, n_families, genes, centroids, gene_to_fam, distances, d)
    
    # Mark output as read-only
    _readonly(distances)
    return distances


def tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff):
    """
    LOESS smoothing in 2D
    
    Args:
        x_ref: Reference x values
        y_ref: Reference y values  
        indices_used: Indices of points to use (1-based for Fortran)
        x_query: Query x values where to compute smoothed y
        kernel_sigma: Kernel bandwidth parameter
        kernel_cutoff: Kernel cutoff parameter
    
    Returns:
        numpy.ndarray: Smoothed y values at query points
    """
    # Input validation and conversion
    x_ref = np.ascontiguousarray(x_ref, dtype=np.float64)
    y_ref = np.ascontiguousarray(y_ref, dtype=np.float64)
    indices_used = np.ascontiguousarray(indices_used, dtype=np.int32)
    x_query = np.ascontiguousarray(x_query, dtype=np.float64)
    
    n_total = len(x_ref)
    n_target = len(x_query)
    n_used = len(indices_used)
    
    if len(y_ref) != n_total:
        raise ValueError("x_ref and y_ref must have same length")
    
    # Validate indices bounds
    if n_used > 0:
        if np.any(indices_used < 1) or np.any(indices_used > n_total):
            raise ValueError("indices_used must be 1-based and within bounds")
    
    # Validate parameters
    if kernel_sigma < 0:
        raise ValueError("kernel_sigma must be non-negative")
    if kernel_cutoff < 0:
        raise ValueError("kernel_cutoff must be non-negative")
    
    # Prepare output and error code
    y_out = np.zeros(n_target, dtype=np.float64)
    ierr = ctypes.c_int(0)
    
    # Setup C wrapper
    loess_smooth_2d_c = lib.loess_smooth_2d_c
    loess_smooth_2d_c.argtypes = [
        ctypes.c_int,  # n_total
        ctypes.c_int,  # n_target
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # x_ref
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # y_ref
        np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),    # indices_used
        ctypes.c_int,  # n_used
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # x_query
        ctypes.c_double,  # kernel_sigma
        ctypes.c_double,  # kernel_cutoff
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # y_out
        ctypes.POINTER(ctypes.c_int),  # ierr
    ]
    loess_smooth_2d_c.restype = None
    
    # Call Fortran routine
    loess_smooth_2d_c(
        n_total, n_target,
        x_ref, y_ref, indices_used, n_used, x_query,
        ctypes.c_double(kernel_sigma), ctypes.c_double(kernel_cutoff),
        y_out, ctypes.byref(ierr)
    )
    
    # Check for errors
    tox_errors(ierr.value)
    
    # Mark output as read-only
    _readonly(y_out)
    return y_out


def tox_compute_family_scaling(distances, gene_to_fam):
    """
    Compute family scaling factors for outlier detection
    
    Args:
        distances: Gene distances to family centroids
        gene_to_fam: Gene-to-family mapping
    
    Returns:
        dict: Dictionary containing scaling factors and intermediate results
    """
    distances = np.ascontiguousarray(distances, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    
    n_genes = len(distances)
    n_families = int(np.max(gene_to_fam)) if len(gene_to_fam) > 0 else 0
    
    if len(gene_to_fam) != n_genes:
        raise ValueError("distances and gene_to_fam must have same length")
    
    # Prepare output arrays
    dscale = np.zeros(n_families, dtype=np.float64)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    indices_used = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    # Setup C wrapper
    compute_family_scaling_c = lib.compute_family_scaling_c
    compute_family_scaling_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # dscale
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_x
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_y
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # indices_used
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # error_code
    ]
    compute_family_scaling_c.restype = None
    
    # Call Fortran routine
    compute_family_scaling_c(
        n_genes, n_families, distances, gene_to_fam,
        dscale, loess_x, loess_y, indices_used, error_code
    )
    
    # Check for errors
    tox_errors(error_code[0])
    
    # Mark outputs as read-only
    _readonly(dscale, loess_x, loess_y, indices_used)
    
    return {
        'dscale': dscale,
        'loess_x': loess_x,
        'loess_y': loess_y,
        'indices_used': indices_used
    }


def tox_compute_family_scaling_expert(distances, gene_to_fam, perm_tmp, stack_left_tmp, 
                                 stack_right_tmp, family_distances):
    """
    Expert version of compute_family_scaling with user-provided work arrays
    
    This version requires pre-allocated work arrays for maximum performance and control.
    Use this when you need fine-grained control over memory allocation or are calling
    this function many times in a tight loop.
    
    Args:
        distances: Gene distances to family centroids
        gene_to_fam: Gene-to-family mapping
        perm_tmp: Pre-allocated permutation array for sorting (n_genes)
        stack_left_tmp: Pre-allocated stack array for sorting (n_genes)
        stack_right_tmp: Pre-allocated stack array for sorting (n_genes)
        family_distances: Pre-allocated work array for family distances (n_genes)
    
    Returns:
        dict: Dictionary containing scaling factors and intermediate results
    """
    distances = np.ascontiguousarray(distances, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    perm_tmp = np.ascontiguousarray(perm_tmp, dtype=np.int32)
    stack_left_tmp = np.ascontiguousarray(stack_left_tmp, dtype=np.int32)
    stack_right_tmp = np.ascontiguousarray(stack_right_tmp, dtype=np.int32)
    family_distances = np.ascontiguousarray(family_distances, dtype=np.float64)
    
    n_genes = len(distances)
    n_families = int(np.max(gene_to_fam)) if len(gene_to_fam) > 0 else 0
    
    if len(gene_to_fam) != n_genes:
        raise ValueError("distances and gene_to_fam must have same length")
    if len(perm_tmp) != n_genes:
        raise ValueError("perm_tmp must have same length as distances")
    if len(stack_left_tmp) != n_genes:
        raise ValueError("stack_left_tmp must have same length as distances")
    if len(stack_right_tmp) != n_genes:
        raise ValueError("stack_right_tmp must have same length as distances")
    if len(family_distances) != n_genes:
        raise ValueError("family_distances must have same length as distances")
    
    # Prepare output arrays
    dscale = np.zeros(n_families, dtype=np.float64)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    indices_used = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    # Setup C wrapper
    compute_family_scaling_expert_c = lib.compute_family_scaling_expert_c
    compute_family_scaling_expert_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # dscale
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_x
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_y
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # indices_used
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm_tmp
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left_tmp
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right_tmp
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # family_distances
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # error_code
    ]
    compute_family_scaling_expert_c.restype = None
    
    # Call Fortran routine
    compute_family_scaling_expert_c(
        n_genes, n_families, distances, gene_to_fam,
        dscale, loess_x, loess_y, indices_used,
        perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code
    )
    
    # Check for errors
    tox_errors(error_code[0])
    
    # Mark outputs as read-only
    _readonly(dscale, loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances)
    
    return {
        'dscale': dscale,
        'loess_x': loess_x,
        'loess_y': loess_y,
        'indices_used': indices_used,
        'perm_tmp': perm_tmp,
        'stack_left_tmp': stack_left_tmp,
        'stack_right_tmp': stack_right_tmp,
        'family_distances': family_distances
    }


def tox_compute_rdi(distances, gene_to_fam, dscale):
    """
    Compute Relative Distance Index (RDI) for outlier detection
    
    Args:
        distances: Gene distances to centroids
        gene_to_fam: Gene-to-family mapping
        dscale: Family scaling factors
    
    Returns:
        numpy.ndarray: RDI values for each gene
    """
    distances = np.ascontiguousarray(distances, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    dscale = np.ascontiguousarray(dscale, dtype=np.float64)
    
    n_genes = len(distances)
    n_families = len(dscale)
    
    # Prepare output and work arrays
    rdi = np.zeros(n_genes, dtype=np.float64)
    sorted_rdi = np.zeros(n_genes, dtype=np.float64)
    perm = np.zeros(n_genes, dtype=np.int32)
    stack_left = np.zeros(n_genes, dtype=np.int32)
    stack_right = np.zeros(n_genes, dtype=np.int32)
    
    # Setup C wrapper
    compute_rdi_c = lib.compute_rdi_c
    compute_rdi_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # dscale
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # rdi
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # sorted_rdi
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right
    ]
    compute_rdi_c.restype = None
    
    # Call Fortran routine
    compute_rdi_c(n_genes, n_families, distances, gene_to_fam, dscale, 
                  rdi, sorted_rdi, perm, stack_left, stack_right)
    
    # Mark output as read-only
    _readonly(rdi)
    return rdi


def tox_identify_outliers(rdi, threshold=None, percentile=95.0):
    """
    Identify outliers based on RDI percentile or threshold
    
    Args:
        rdi: Relative Distance Index values
        threshold: Fixed RDI threshold (if None, uses percentile)
        percentile: Percentile threshold for outlier detection (default: 95 for top 5%)
    
    Returns:
        dict: Dictionary containing:
            - outliers: Boolean array indicating outliers
            - threshold: Threshold value used for detection
    """
    rdi = np.ascontiguousarray(rdi, dtype=np.float64)
    n_genes = len(rdi)
    
    # Prepare sorted RDI (copy and filter out negatives)
    sorted_rdi = rdi.copy()
    sorted_rdi[sorted_rdi < 0] = 0.0  # Filter out error values
    sorted_rdi.sort()  # Sort in ascending order
    
    # Prepare output arrays
    outliers_int = np.zeros(n_genes, dtype=np.int32)
    threshold_out = ctypes.c_double(0.0)
    
    # Setup C wrapper
    identify_outliers_c = lib.identify_outliers_c
    identify_outliers_c.argtypes = [
        ctypes.c_int,  # n_genes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # rdi
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # sorted_rdi
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # is_outlier_int
        ctypes.POINTER(ctypes.c_double),  # threshold (output)
        ctypes.c_double,  # percentile
    ]
    identify_outliers_c.restype = None
    
    # Call Fortran routine
    identify_outliers_c(n_genes, rdi, sorted_rdi, outliers_int, 
                        ctypes.byref(threshold_out), ctypes.c_double(percentile))
    
    # Mark output as read-only
    _readonly(outliers_int)
    
    return {
        'outliers': outliers_int,
        'threshold': threshold_out.value
    }


def tox_detect_outliers(distances, gene_to_fam, percentile=95.0):
    """
    Complete outlier detection pipeline
    
    Args:
        distances: Gene distances to centroids
        gene_to_fam: Gene-to-family mapping
        percentile: Percentile threshold for outlier detection (default: 95 for top 5%)
    
    Returns:
        dict: Dictionary containing outliers and intermediate results
    """
    distances = np.ascontiguousarray(distances, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    
    n_genes = len(distances)
    n_families = int(np.max(gene_to_fam)) if len(gene_to_fam) > 0 else 0
    
    # Prepare work arrays
    work_array = np.zeros(n_genes, dtype=np.float64)
    perm = np.zeros(n_genes, dtype=np.int32)
    stack_left = np.zeros(n_genes, dtype=np.int32)
    stack_right = np.zeros(n_genes, dtype=np.int32)
    
    # Prepare output arrays
    outliers_int = np.zeros(n_genes, dtype=np.int32)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    loess_n = np.zeros(n_families, dtype=np.int32)
    error_code = ctypes.c_int(0)
    
    # Setup C wrapper
    detect_outliers_c = lib.detect_outliers_c
    detect_outliers_c.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # work_array
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # outliers_int
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_x
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_y
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # loess_n
        ctypes.POINTER(ctypes.c_int),  # error_code
        ctypes.c_double,  # percentile
    ]
    detect_outliers_c.restype = None
    
    # Call Fortran routine
    detect_outliers_c(
        n_genes, n_families, distances, gene_to_fam,
        work_array, perm, stack_left, stack_right,
        outliers_int, loess_x, loess_y, loess_n,
        ctypes.byref(error_code), ctypes.c_double(percentile)
    )
    
    # Check for errors
    tox_errors(error_code.value)
    
    # Mark outputs as read-only
    _readonly(outliers_int, loess_x, loess_y, loess_n)
    
    return {
        'outliers': outliers_int,
        'loess_x': loess_x,
        'loess_y': loess_y,
        'loess_n': loess_n
    }


def tox_which(cond):
    """
    'which' utility for Python, like in R/MATLAB.
    Returns indices of TRUE elements in a logical array.
    
    Args:
        cond: Array of boolean/integer values (0/1)
    
    Returns:
        tuple: (idx_out, m_out) where idx_out contains 1-based indices and m_out is count
    """
    cond = np.ascontiguousarray(cond, dtype=np.int32)
    
    # Check for invalid input
    if np.any(np.isnan(cond.astype(float))):
        raise ValueError("Input contains NaN values")
    if np.any(np.isinf(cond.astype(float))):
        raise ValueError("Input contains infinite values")
    if not np.all((cond == 0) | (cond == 1)):
        raise ValueError("Input must contain only 0 and 1 values")
    
    n = cond.size
    idx_out = np.zeros(n, dtype=np.int32)
    m_max = n
    m_out = np.zeros(1, dtype=np.int32)
    error_code = ctypes.c_int(0)
    
    # Setup C wrapper
    which_c = lib.which_c
    which_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # mask
        ctypes.c_int,                 # n
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # idx_out
        ctypes.c_int,                 # m_max
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # m_out
        ctypes.POINTER(ctypes.c_int)  # error_code
    ]
    which_c.restype = None
    
    # Call Fortran subroutine
    which_c(
        cond, ctypes.c_int(n), idx_out, ctypes.c_int(m_max), m_out,
        ctypes.byref(error_code)
    )
    
    # Check for errors
    tox_errors(error_code.value)
    
    # Mark output as read-only
    _readonly(idx_out)

    return idx_out, int(m_out[0])

def tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid):
    """
    Calculate Shift Vector Field 

    Computes the shift vector field for each gene expression vector based on its family centroid.
    The shift vector is defined as the difference between the gene expression vector and its corresponding family centroid,
    starting at the expression vector and pointing to its family centroid.
    This function automatically checks for errors and throws informative exceptions.
    
    Args:
        expression_vectors: Matrix where each column is a gene expression vector (n_axes x n_vectors)
        family_centroids: Matrix where each column is a family centroid vector (n_axes x n_families)
        gene_to_centroid: Array mapping each gene to its corresponding family centroid index in family_centroids with length n_vectors (1 based for fortran)

    Returns:
        dict: Dictionary containing:
            - shift_vectors: The computed shift vectors for each gene expression vector

    """
    
    # Input validation
    if not isinstance(expression_vectors, np.ndarray):
        raise ValueError("expression_vectors must be a numpy array")
    
    if expression_vectors.ndim != 2:
        raise ValueError("expression_vectors must be a 2D array")

    if not isinstance(family_centroids, np.ndarray):
        raise ValueError("family_centroids must be a numpy array")

    if family_centroids.ndim != 2:
        raise ValueError("family_centroids must be a 2D array")
    
    # Convert inputs to numpy arrays
    expression_vectors = np.asarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asarray(family_centroids, dtype=np.float64)
    gene_to_centroid = np.asarray(gene_to_centroid, dtype=np.int32)
    
    # Get dimensions
    n_axes_genes, n_vectors = expression_vectors.shape
    n_axes_centroids, n_families = family_centroids.shape
    
    # Validate length of gene_to_centroid
    if len(gene_to_centroid) != n_vectors:
        raise ValueError("gene_to_centroid length must match number of columns in expression_vectors")
    
    # Validate dimensions
    if n_axes_genes != n_axes_centroids:
        raise ValueError("family_centroids must have the same number of axes as expression_vectors")

    # Validate correct mapping between genes and centroids
    if not np.all((gene_to_centroid > 0) & (gene_to_centroid <= n_families)):
        raise ValueError("gene_to_centroid contains invalid family IDs")

    # Ensure arrays have correct dtype and memory layout
    expr_v = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_c = np.asfortranarray(family_centroids, dtype=np.float64)
    gene_c = np.ascontiguousarray(gene_to_centroid.astype(np.int32))
    
    # Prepare output arrays
    shift_vectors = np.empty((2*n_axes_genes, n_vectors), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # Setup C/Fortran wrapper with proper type annotations
    sv = lib.compute_shift_vector_field_c
    sv.argtypes = [
        ctypes.c_int,  # n_axes_genes
        ctypes.c_int,  # n_vectors
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # expression_vectors (Fortran order)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # family_centroids (Fortran order)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_centroid
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # shift_vectors
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    sv.restype = None
    
    # Call the C/Fortran wrapper
    sv(
        n_axes_genes,
        n_vectors,
        n_families,
        expr_v,
        family_c,
        gene_c,
        shift_vectors,
        ctypes.byref(ierr)
    )

    # Check for errors and throw informative messages
    tox_errors(ierr.value)

    # Mark outputs as read-only
    _readonly(shift_vectors)
    
    # Return structured result (no ierr since we checked for errors)
    return {
        "shift_vectors": shift_vectors
    }

def _string_matrix_to_flat_ascii(str_matrix: np.ndarray) -> np.ndarray:
    """Converts a 2D numpy string array to a flat C-style ASCII integer array."""
    rows, cols = str_matrix.shape
    flat_ascii = np.zeros(rows * cols * MAX_FIELD_LEN, dtype=np.int32)
    for i in range(rows):
        for j in range(cols):
            s = str_matrix[i, j]
            for k, char in enumerate(s):
                if k >= MAX_FIELD_LEN: break
                idx = i * cols * MAX_FIELD_LEN + j * MAX_FIELD_LEN + k
                flat_ascii[idx] = ord(char)
    return flat_ascii

def _ascii_to_string_matrix(ascii_flat: np.ndarray, rows: int, cols: int) -> np.ndarray:
    """Converts a flat ASCII buffer from Fortran into a 2D numpy string array."""
    reshaped = ascii_flat.reshape((rows, cols, MAX_FIELD_LEN))
    str_matrix = np.empty((rows, cols), dtype=object)
    for i in range(rows):
        for j in range(cols):
            null_pos = np.where(reshaped[i, j] == 0)[0]
            end_pos = null_pos[0] if len(null_pos) > 0 else MAX_FIELD_LEN
            str_matrix[i, j] = "".join(chr(c) for c in reshaped[i, j, :end_pos])
    
    # Set the final array to read-only as required
    final_matrix = np.array(str_matrix, dtype=str)
    final_matrix.flags.writeable = False
    # NOTE: Returned NumPy arrays are read-only for safety.
    # If you need to modify them (e.g., for plotting), use .copy().
    return final_matrix

# --- User-Facing API Functions ---

def read_csv_dimensions(filepath: str, has_header: bool = True, delimiter: str = ',') -> tuple[int, int]:
    """
    Gets the dimensions (rows, columns) of a CSV file without reading its data.

    Parameters
    ----------
    filepath : str
        Path to the CSV file.
    has_header : bool, optional
        Whether the file has a header row, by default True.
    delimiter : str, optional
        The column delimiter, by default ','.

    Returns
    -------
    tuple[int, int]
        A tuple containing (number_of_data_rows, number_of_columns).
    """
    if not Path(filepath).exists():
        raise FileNotFoundError(f"File not found at '{filepath}'")
    if len(delimiter) != 1:
        raise ValueError("Delimiter must be a single character.")

    fname_ascii = np.array([ord(c) for c in filepath], dtype=np.int32)
    num_rows, num_cols, status = ctypes.c_int(0), ctypes.c_int(0), ctypes.c_int(0)

    f = fortran_lib.get_csv_dims_c
    f.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
        ctypes.c_int, ctypes.c_bool, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)
    ]
    f.restype = None
    
    f(fname_ascii, len(filepath), has_header, ord(delimiter),
      ctypes.byref(num_rows), ctypes.byref(num_cols), ctypes.byref(status))
    
    _handle_fortran_errors(status.value, "read_csv_dimensions")
    return num_rows.value, num_cols.value

def read_csv_as_strings(filepath: str, has_header: bool = True, delimiter: str = ',') -> tuple[np.ndarray | None, np.ndarray]:
    """
    Reads a CSV file into a NumPy array of strings.

    Parameters
    ----------
    filepath : str
        Path to the CSV file.
    has_header : bool, optional
        Whether the file has a header row, by default True.
    delimiter : str, optional
        The column delimiter, by default ','.

    Returns
    -------
    tuple[np.ndarray | None, np.ndarray]
        A tuple containing (header, data). Header is None if has_header is False.
        Data is a 2D NumPy array of strings. Returned arrays are read-only.
    """
    num_rows, num_cols = read_csv_dimensions(filepath, has_header, delimiter)
    fname_ascii = np.array([ord(c) for c in filepath], dtype=np.int32)
    status = ctypes.c_int(0)

    header_buffer = np.zeros(num_cols * MAX_FIELD_LEN, dtype=np.int32) if has_header else np.array([], dtype=np.int32)
    data_buffer = np.zeros(num_rows * num_cols * MAX_FIELD_LEN, dtype=np.int32)
    
    f = fortran_lib.read_csv_to_strings_c
    f.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        ctypes.c_bool, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), np.ctypeslib.ndpointer(dtype=np.int32),
        ctypes.POINTER(ctypes.c_int)
    ]
    f.restype = None
    
    f(fname_ascii, len(filepath), has_header, ord(delimiter),
      header_buffer, data_buffer, ctypes.byref(status))
    
    _handle_fortran_errors(status.value, "read_csv_as_strings")
    
    data_out = _ascii_to_string_matrix(data_buffer, num_rows, num_cols)
    header_out = _ascii_to_string_matrix(header_buffer, 1, num_cols).flatten() if has_header else None
    
    if header_out is not None:
        header_out.flags.writeable = False
        # NOTE: Returned NumPy arrays are read-only for safety.
        # If you need to modify them (e.g., for plotting), use .copy().
    
    return header_out, data_out

def _read_typed_columns(string_data: np.ndarray, columns: list[int], dtype: np.dtype, fortran_func_name: str) -> np.ndarray:
    """Generic internal function to convert string data to a numeric type."""
    if not isinstance(string_data, np.ndarray) or string_data.ndim != 2:
        raise TypeError("string_data must be a 2D NumPy array.")
    if not all(isinstance(c, int) and c > 0 for c in columns):
        raise ValueError("columns must be a list of 1-based positive integers.")

    num_rows, num_cols_in = string_data.shape
    num_cols_to_read = len(columns)
    
    flat_ascii_input = _string_matrix_to_flat_ascii(string_data)
    cols_to_read_f = np.array(columns, dtype=np.int32)
    output_buffer = np.zeros(num_rows * num_cols_to_read, dtype=dtype)
    status = ctypes.c_int(0)
    
    f = getattr(fortran_lib, fortran_func_name)
    f.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=dtype), ctypes.POINTER(ctypes.c_int)
    ]
    f.restype = None

    f(flat_ascii_input, num_rows, num_cols_in, cols_to_read_f, num_cols_to_read, output_buffer, ctypes.byref(status))
    _handle_fortran_errors(status.value, fortran_func_name)
    
    # Reshape to Fortran order, then transpose to get C order semantics
    reshaped_output = output_buffer.reshape((num_cols_to_read, num_rows)).T
    
    reshaped_output.flags.writeable = False
    # NOTE: Returned NumPy arrays are read-only for safety.
    # If you need to modify them (e.g., for plotting), use .copy().
    return reshaped_output

def read_integer_columns(string_data: np.ndarray, columns: list[int]) -> np.ndarray:
    """Extracts and converts specified columns to integers."""
    return _read_typed_columns(string_data, columns, np.int32, 'read_integer_columns_c')

def read_real_columns(string_data: np.ndarray, columns: list[int]) -> np.ndarray:
    """Extracts and converts specified columns to double precision reals."""
    return _read_typed_columns(string_data, columns, np.float64, 'read_real_columns_c')

def read_logical_columns(string_data: np.ndarray, columns: list[int]) -> np.ndarray:
    """Extracts and converts specified columns to booleans."""
    return _read_typed_columns(string_data, columns, np.bool_, 'read_logical_columns_c')
    
def read_complex_columns(string_data: np.ndarray, columns: list[int]) -> np.ndarray:
    """Extracts and converts specified columns to double precision complex numbers."""
    return _read_typed_columns(string_data, columns, np.complex128, 'read_complex_columns_c')

def read_character_columns(string_data: np.ndarray, columns: list[int]) -> np.ndarray:
    """Extracts specified string columns without type conversion."""
    if not isinstance(string_data, np.ndarray) or string_data.ndim != 2:
        raise TypeError("string_data must be a 2D NumPy array.")
    if not all(isinstance(c, int) and c > 0 for c in columns):
        raise ValueError("columns must be a list of 1-based positive integers.")

    # Fortran columns are 1-based, Python is 0-based.
    py_indices = [c - 1 for c in columns]
    output_array = string_data[:, py_indices]
    
    output_array.flags.writeable = False
    # NOTE: Returned NumPy arrays are read-only for safety.
    # If you need to modify them (e.g., for plotting), use .copy().
    return output_array