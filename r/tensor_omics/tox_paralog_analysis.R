# Generated. Do not edit.

#' Identifies neofunctionalization for genes by checking whether the difference of expression to its ancestor exceeds the threshold for the respective axis
#'
#' @param ancestors a numeric matrix. RAP projected unit length expression vector of ancestral ortholog
#' @param genes a numeric matrix. RAP projected unit length expression vectors of genes
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0` for unassigned genes
#' @param thresholds a numeric vector. threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis
#' @return `TRUE` if neofunctionalization has been detected for the respective axes, always `FALSE` for unassigned genes
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::detect_neofunctionalization}.
#' @export
detect_neofunctionalization <- function(ancestors, genes, gene_to_fam, thresholds) {
    ancestors <- .tox_as_double_matrix(ancestors, "ancestors")
    genes <- .tox_as_double_matrix(genes, "genes")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    thresholds <- .tox_as_double_vector(thresholds, "thresholds")
    if (dim(genes)[1] != dim(ancestors)[1])
        .tox_shape_error("genes", dim(genes)[1], "ancestors", dim(ancestors)[1])
    if (length(thresholds) != dim(ancestors)[1])
        .tox_shape_error("thresholds", length(thresholds), "ancestors", dim(ancestors)[1])
    if (length(gene_to_fam) != dim(genes)[2])
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "genes", dim(genes)[2])

    .result <- .Call("detect_neofunctionalization_call", ancestors, genes, gene_to_fam, thresholds)
    .arguments <- c("ancestors", "n_families", "genes", "n_axes", "gene_to_fam", "n_genes", "thresholds", "neofunc", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$neofunc
}

#' Identifies subsets of paralogs matching this pattern
#'
#' @param ancestor a numeric vector. expression vector of ancestral ortholog
#' @param genes a numeric matrix. expression vectors of genes
#' @param filtered_paralogs_mask a integer vector. bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
#' @param max_subset_size a integer scalar. maximum subset size of checked gene subsets. Too large a value is capped to the
#' @param max_angle a numeric scalar. maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
#' @param gain_gamma a numeric scalar. positive magnitude gain for dosage effect
#' @return a named list with elements `n_results`, `work_arr_paralog_subsets`.
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::detect_dosage_effect}.
#' @export
detect_dosage_effect_expert <- function(ancestor, genes, filtered_paralogs_mask, max_subset_size, max_angle = 3.141592653589793, gain_gamma = 0.1) {
    ancestor <- .tox_as_double_vector(ancestor, "ancestor")
    genes <- .tox_as_double_matrix(genes, "genes")
    filtered_paralogs_mask <- .tox_as_integer_vector(filtered_paralogs_mask, "filtered_paralogs_mask")
    max_subset_size <- .tox_as_integer_scalar(max_subset_size, "max_subset_size")
    max_angle <- .tox_as_double_scalar(max_angle, "max_angle")
    gain_gamma <- .tox_as_double_scalar(gain_gamma, "gain_gamma")
    .calc_work_arr_paralog_subsets_size_result <- calc_work_arr_paralog_subsets_size(max_subset_size = max_subset_size, n_genes = dim(genes)[2], filtered_paralogs_mask = filtered_paralogs_mask)
    max_subset_size <- .calc_work_arr_paralog_subsets_size_result$max_subset_size
    n_paralog_subsets <- .calc_work_arr_paralog_subsets_size_result$work_array_size

    if (dim(genes)[1] != length(ancestor))
        .tox_shape_error("genes", dim(genes)[1], "ancestor", length(ancestor))

    .result <- .Call("detect_dosage_effect_expert_call", ancestor, genes, filtered_paralogs_mask, max_subset_size, n_paralog_subsets, max_angle, gain_gamma)
    .arguments <- c("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "tmp_active_mask", "tmp_paralog_vector", "max_angle", "gain_gamma", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_results = .result$n_results,
        work_arr_paralog_subsets = .result$work_arr_paralog_subsets
    )
}

#' Identifies subsets of paralogs matching this pattern
#'
#' @param ancestor a numeric vector. expression vector of ancestral ortholog
#' @param genes a numeric matrix. expression vectors of genes
#' @param filtered_paralogs_mask a integer vector. bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
#' @param max_subset_size a integer scalar. maximum subset size of checked gene subsets. Too large a value is capped to the
#' @param max_angle a numeric scalar. maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
#' @param gain_gamma a numeric scalar. positive magnitude gain for dosage effect
#' @return a named list with elements `n_results`, `work_arr_paralog_subsets`.
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::detect_dosage_effect_alloc}.
#' @export
detect_dosage_effect <- function(ancestor, genes, filtered_paralogs_mask, max_subset_size, max_angle = 3.141592653589793, gain_gamma = 0.1) {
    ancestor <- .tox_as_double_vector(ancestor, "ancestor")
    genes <- .tox_as_double_matrix(genes, "genes")
    filtered_paralogs_mask <- .tox_as_integer_vector(filtered_paralogs_mask, "filtered_paralogs_mask")
    max_subset_size <- .tox_as_integer_scalar(max_subset_size, "max_subset_size")
    max_angle <- .tox_as_double_scalar(max_angle, "max_angle")
    gain_gamma <- .tox_as_double_scalar(gain_gamma, "gain_gamma")
    .calc_work_arr_paralog_subsets_size_result <- calc_work_arr_paralog_subsets_size(max_subset_size = max_subset_size, n_genes = dim(genes)[2], filtered_paralogs_mask = filtered_paralogs_mask)
    max_subset_size <- .calc_work_arr_paralog_subsets_size_result$max_subset_size
    n_paralog_subsets <- .calc_work_arr_paralog_subsets_size_result$work_array_size

    if (dim(genes)[1] != length(ancestor))
        .tox_shape_error("genes", dim(genes)[1], "ancestor", length(ancestor))

    .result <- .Call("detect_dosage_effect_call", ancestor, genes, filtered_paralogs_mask, max_subset_size, n_paralog_subsets, max_angle, gain_gamma)
    .arguments <- c("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "max_angle", "gain_gamma", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_results = .result$n_results,
        work_arr_paralog_subsets = .result$work_arr_paralog_subsets
    )
}

#' Identifies subsets of paralogs matching this pattern
#'
#' @param ancestor a numeric vector. expression vector of ancestral ortholog
#' @param genes a numeric matrix. expression vectors of genes
#' @param filtered_paralogs_mask a integer vector. bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
#' @param max_subset_size a integer scalar. maximum subset size of checked gene subsets. Too large a value is capped to the
#' @param rdi_threshold a numeric scalar. max allowed residual distance from `ancestor`
#' @param paralog_norms a numeric vector. euclidean norms of the genes, used for subset pruning (`norm` from `f42_utils` computes them)
#' @param sorted_paralog_norms_perm a integer vector. ascending permutation of the norms, for subset pruning: the smallest norm among the genes that could extend a subset must not fall below the subset's angle to the ancestor
#' @return a named list with elements `n_results`, `work_arr_paralog_subsets`.
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::detect_subfunctionalization}.
#' @export
detect_subfunctionalization_expert <- function(ancestor, genes, filtered_paralogs_mask, max_subset_size, rdi_threshold, paralog_norms, sorted_paralog_norms_perm) {
    ancestor <- .tox_as_double_vector(ancestor, "ancestor")
    genes <- .tox_as_double_matrix(genes, "genes")
    filtered_paralogs_mask <- .tox_as_integer_vector(filtered_paralogs_mask, "filtered_paralogs_mask")
    max_subset_size <- .tox_as_integer_scalar(max_subset_size, "max_subset_size")
    rdi_threshold <- .tox_as_double_scalar(rdi_threshold, "rdi_threshold")
    paralog_norms <- .tox_as_double_vector(paralog_norms, "paralog_norms")
    sorted_paralog_norms_perm <- .tox_as_integer_vector(sorted_paralog_norms_perm, "sorted_paralog_norms_perm")
    .calc_work_arr_paralog_subsets_size_result <- calc_work_arr_paralog_subsets_size(max_subset_size = max_subset_size, n_genes = dim(genes)[2], filtered_paralogs_mask = filtered_paralogs_mask)
    max_subset_size <- .calc_work_arr_paralog_subsets_size_result$max_subset_size
    n_paralog_subsets <- .calc_work_arr_paralog_subsets_size_result$work_array_size

    if (length(paralog_norms) != dim(genes)[2])
        .tox_shape_error("paralog_norms", length(paralog_norms), "genes", dim(genes)[2])
    if (length(sorted_paralog_norms_perm) != dim(genes)[2])
        .tox_shape_error("sorted_paralog_norms_perm", length(sorted_paralog_norms_perm), "genes", dim(genes)[2])
    if (dim(genes)[1] != length(ancestor))
        .tox_shape_error("genes", dim(genes)[1], "ancestor", length(ancestor))

    .result <- .Call("detect_subfunctionalization_expert_call", ancestor, genes, filtered_paralogs_mask, max_subset_size, n_paralog_subsets, rdi_threshold, paralog_norms, sorted_paralog_norms_perm)
    .arguments <- c("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "tmp_active_mask", "tmp_paralog_vector", "rdi_threshold", "paralog_norms", "sorted_paralog_norms_perm", "tmp_work_array", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_results = .result$n_results,
        work_arr_paralog_subsets = .result$work_arr_paralog_subsets
    )
}

#' Identifies subsets of paralogs matching this pattern
#'
#' @param ancestor a numeric vector. expression vector of ancestral ortholog
#' @param genes a numeric matrix. expression vectors of genes
#' @param filtered_paralogs_mask a integer vector. bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
#' @param max_subset_size a integer scalar. maximum subset size of checked gene subsets. Too large a value is capped to the
#' @param rdi_threshold a numeric scalar. max allowed residual distance from `ancestor`
#' @param paralog_norms a numeric vector. euclidean norms of the genes, used for subset pruning (`norm` from `f42_utils` computes them)
#' @param sorted_paralog_norms_perm a integer vector. ascending permutation of the norms, for subset pruning: the smallest norm among the genes that could extend a subset must not fall below the subset's angle to the ancestor
#' @return a named list with elements `n_results`, `work_arr_paralog_subsets`.
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::detect_subfunctionalization_alloc}.
#' @export
detect_subfunctionalization <- function(ancestor, genes, filtered_paralogs_mask, max_subset_size, rdi_threshold, paralog_norms, sorted_paralog_norms_perm) {
    ancestor <- .tox_as_double_vector(ancestor, "ancestor")
    genes <- .tox_as_double_matrix(genes, "genes")
    filtered_paralogs_mask <- .tox_as_integer_vector(filtered_paralogs_mask, "filtered_paralogs_mask")
    max_subset_size <- .tox_as_integer_scalar(max_subset_size, "max_subset_size")
    rdi_threshold <- .tox_as_double_scalar(rdi_threshold, "rdi_threshold")
    paralog_norms <- .tox_as_double_vector(paralog_norms, "paralog_norms")
    sorted_paralog_norms_perm <- .tox_as_integer_vector(sorted_paralog_norms_perm, "sorted_paralog_norms_perm")
    .calc_work_arr_paralog_subsets_size_result <- calc_work_arr_paralog_subsets_size(max_subset_size = max_subset_size, n_genes = dim(genes)[2], filtered_paralogs_mask = filtered_paralogs_mask)
    max_subset_size <- .calc_work_arr_paralog_subsets_size_result$max_subset_size
    n_paralog_subsets <- .calc_work_arr_paralog_subsets_size_result$work_array_size

    if (length(paralog_norms) != dim(genes)[2])
        .tox_shape_error("paralog_norms", length(paralog_norms), "genes", dim(genes)[2])
    if (length(sorted_paralog_norms_perm) != dim(genes)[2])
        .tox_shape_error("sorted_paralog_norms_perm", length(sorted_paralog_norms_perm), "genes", dim(genes)[2])
    if (dim(genes)[1] != length(ancestor))
        .tox_shape_error("genes", dim(genes)[1], "ancestor", length(ancestor))

    .result <- .Call("detect_subfunctionalization_call", ancestor, genes, filtered_paralogs_mask, max_subset_size, n_paralog_subsets, rdi_threshold, paralog_norms, sorted_paralog_norms_perm)
    .arguments <- c("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "rdi_threshold", "paralog_norms", "sorted_paralog_norms_perm", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_results = .result$n_results,
        work_arr_paralog_subsets = .result$work_arr_paralog_subsets
    )
}

#' Prefilters the genes for a pattern, so genes that cannot match it are not tried as subset extensions
#'
#' This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
#'
#' @param gene_angles a numeric vector. vector, holding the angles between ancestor and genes (0<=angle<=Pi)
#' @param threshold a numeric scalar. filter threshold
#' @param n_families a integer scalar. number of families
#' @param gene_to_fam a integer vector. a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
#' @param n_mask_chunks a integer scalar. number of 32 bit chunks a mask needs to encode `n_genes` genes
#' @return bit mask that will have the indices of genes kept by this pattern set to 1, else 0
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::filter_paralogs_by_pattern_dosage_effect}.
#' @export
filter_paralogs_by_pattern_dosage_effect <- function(gene_angles, threshold, n_families, gene_to_fam, n_mask_chunks) {
    gene_angles <- .tox_as_double_vector(gene_angles, "gene_angles")
    threshold <- .tox_as_double_scalar(threshold, "threshold")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    n_mask_chunks <- .tox_as_integer_scalar(n_mask_chunks, "n_mask_chunks")
    if (length(gene_to_fam) != length(gene_angles))
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "gene_angles", length(gene_angles))

    .result <- .Call("filter_paralogs_by_pattern_dosage_effect_call", gene_angles, threshold, n_families, gene_to_fam, n_mask_chunks)
    .arguments <- c("gene_angles", "threshold", "n_genes", "n_families", "gene_to_fam", "masks", "n_mask_chunks", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$masks
}

#' Prefilters the genes for a pattern, so genes that cannot match it are not tried as subset extensions
#'
#' This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
#'
#' @param gene_angles a numeric vector. vector, holding the angles between ancestor and genes (0<=angle<=Pi)
#' @param threshold a numeric scalar. filter threshold
#' @param n_families a integer scalar. number of families
#' @param gene_to_fam a integer vector. a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
#' @param n_mask_chunks a integer scalar. number of 32 bit chunks a mask needs to encode `n_genes` genes
#' @return bit mask that will have the indices of genes kept by this pattern set to 1, else 0
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::filter_paralogs_by_pattern_subfunctionalization}.
#' @export
filter_paralogs_by_pattern_subfunctionalization <- function(gene_angles, threshold, n_families, gene_to_fam, n_mask_chunks) {
    gene_angles <- .tox_as_double_vector(gene_angles, "gene_angles")
    threshold <- .tox_as_double_scalar(threshold, "threshold")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    n_mask_chunks <- .tox_as_integer_scalar(n_mask_chunks, "n_mask_chunks")
    if (length(gene_to_fam) != length(gene_angles))
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "gene_angles", length(gene_angles))

    .result <- .Call("filter_paralogs_by_pattern_subfunctionalization_call", gene_angles, threshold, n_families, gene_to_fam, n_mask_chunks)
    .arguments <- c("gene_angles", "threshold", "n_genes", "n_families", "gene_to_fam", "masks", "n_mask_chunks", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$masks
}
