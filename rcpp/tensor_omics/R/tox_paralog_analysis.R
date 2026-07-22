# Generated. Do not edit.

#' Checks the state of a bit/paralog in `bit_mask` -> .true. if 1 else .false.
#'
#' @param bit_mask a integer vector. chunked mask to mark active paralogs
#' @param i_gene a integer scalar. index of paralog to be marked active
#' @return check result
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::mask_check_state}.
#' @export
mask_check_state <- function(bit_mask, i_gene) {
    bit_mask <- .tox_as_integer_vector(bit_mask, "bit_mask")
    i_gene <- .tox_as_integer_scalar(i_gene, "i_gene")
    .result <- .mask_check_state_rcpp(bit_mask, i_gene)
    .arguments <- c("bit_mask", "i_gene")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$state
}

#' Identifies neofunctionalization for genes by checking whether the difference of expression to its ancestor exceeds the threshold for the respective axis
#'
#' @param ancestors a numeric matrix. RAP projected unit length expression vector of ancestral ortholog
#' @param genes a numeric matrix. RAP projected unit length expression vectors of genes
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0_int32` for unassigned genes
#' @param thresholds a numeric vector. threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis
#' @return `.true.` if neofunctionalization has been detected for the respective axes, always `.false.` for unassigned genes
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

    .result <- .detect_neofunctionalization_rcpp(ancestors, genes, gene_to_fam, thresholds)
    .arguments <- c("ancestors", "n_families", "genes", "n_axes", "gene_to_fam", "n_genes", "thresholds", "neofunc", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$neofunc
}

#' Identifies subsets of paralogs with small angle to the `ancestor` (max_angle) and sum to a magnitude significantly exceeding `norm(ancestor)` (gain)
#'
#' @param ancestor a numeric vector. expression vector of ancestral ortholog
#' @param genes a numeric matrix. expression vectors of genes
#' @param filtered_paralogs_mask a integer vector. bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
#' @param max_subset_size a integer scalar. maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
#' @param n_paralog_subsets a integer scalar. number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
#' @param max_angle a numeric scalar. in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned, default is Pi
#' @param gain_gamma a numeric scalar. positive magnitude gain for dosage effect, default 0.1
#' @return a named list with elements `n_results`, `work_arr_paralog_subsets`.
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::detect_dosage_effect}.
#' @export
detect_dosage_effect <- function(ancestor, genes, filtered_paralogs_mask, max_subset_size, n_paralog_subsets, max_angle = NULL, gain_gamma = NULL) {
    ancestor <- .tox_as_double_vector(ancestor, "ancestor")
    genes <- .tox_as_double_matrix(genes, "genes")
    filtered_paralogs_mask <- .tox_as_integer_vector(filtered_paralogs_mask, "filtered_paralogs_mask")
    max_subset_size <- .tox_as_integer_scalar(max_subset_size, "max_subset_size")
    n_paralog_subsets <- .tox_as_integer_scalar(n_paralog_subsets, "n_paralog_subsets")
    if (!is.null(max_angle))
        max_angle <- .tox_as_double_scalar(max_angle, "max_angle")
    if (!is.null(gain_gamma))
        gain_gamma <- .tox_as_double_scalar(gain_gamma, "gain_gamma")
    if (dim(genes)[1] != length(ancestor))
        .tox_shape_error("genes", dim(genes)[1], "ancestor", length(ancestor))

    .result <- .detect_dosage_effect_rcpp(ancestor, genes, filtered_paralogs_mask, max_subset_size, n_paralog_subsets, max_angle, gain_gamma)
    .arguments <- c("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "tmp_active_mask", "tmp_paralog_vector", "ierr", "max_angle", "gain_gamma")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_results = .result$n_results,
        work_arr_paralog_subsets = .result$work_arr_paralog_subsets
    )
}

#' Identifies subsets of paralogs exhibiting significant angles to the `ancestor`
#'
#' @param ancestor a numeric vector. expression vector of ancestral ortholog
#' @param genes a numeric matrix. expression vectors of genes
#' @param rdi_threshold a numeric scalar. max allowed residual distance from `ancestor`
#' @param filtered_paralogs_mask a integer vector. bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
#' @param max_subset_size a integer scalar. maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
#' @param n_paralog_subsets a integer scalar. number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
#' @param paralog_norms a numeric vector. needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` function from `f42_utils` function for this)
#' @param sorted_paralog_norms_perm a integer vector. needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor
#' @return a named list with elements `n_results`, `work_arr_paralog_subsets`.
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::detect_subfunctionalization}.
#' @export
detect_subfunctionalization <- function(ancestor, genes, rdi_threshold, filtered_paralogs_mask, max_subset_size, n_paralog_subsets, paralog_norms, sorted_paralog_norms_perm) {
    ancestor <- .tox_as_double_vector(ancestor, "ancestor")
    genes <- .tox_as_double_matrix(genes, "genes")
    rdi_threshold <- .tox_as_double_scalar(rdi_threshold, "rdi_threshold")
    filtered_paralogs_mask <- .tox_as_integer_vector(filtered_paralogs_mask, "filtered_paralogs_mask")
    max_subset_size <- .tox_as_integer_scalar(max_subset_size, "max_subset_size")
    n_paralog_subsets <- .tox_as_integer_scalar(n_paralog_subsets, "n_paralog_subsets")
    paralog_norms <- .tox_as_double_vector(paralog_norms, "paralog_norms")
    sorted_paralog_norms_perm <- .tox_as_integer_vector(sorted_paralog_norms_perm, "sorted_paralog_norms_perm")
    if (length(paralog_norms) != dim(genes)[2])
        .tox_shape_error("paralog_norms", length(paralog_norms), "genes", dim(genes)[2])
    if (length(sorted_paralog_norms_perm) != dim(genes)[2])
        .tox_shape_error("sorted_paralog_norms_perm", length(sorted_paralog_norms_perm), "genes", dim(genes)[2])
    if (dim(genes)[1] != length(ancestor))
        .tox_shape_error("genes", dim(genes)[1], "ancestor", length(ancestor))

    .result <- .detect_subfunctionalization_rcpp(ancestor, genes, rdi_threshold, filtered_paralogs_mask, max_subset_size, n_paralog_subsets, paralog_norms, sorted_paralog_norms_perm)
    .arguments <- c("ancestor", "genes", "n_genes", "n_dims", "rdi_threshold", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "tmp_active_mask", "tmp_paralog_vector", "paralog_norms", "sorted_paralog_norms_perm", "tmp_work_array", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_results = .result$n_results,
        work_arr_paralog_subsets = .result$work_arr_paralog_subsets
    )
}

#' Determines the needed chunk count for subset bit masks (an integer has only 32 bits)
#'
#' @param n_genes a integer scalar. number of genes
#' @return number of 32 bit chunks a mask needs to encode `n_genes` genes
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::mask_chunk_count}.
#' @export
mask_chunk_count <- function(n_genes) {
    n_genes <- .tox_as_integer_scalar(n_genes, "n_genes")
    .result <- .mask_chunk_count_rcpp(n_genes)
    .arguments <- c("n_genes", "count")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$count
}

#' Prefilters the genes for subfunctionalization
#'
#' as genes that are already too close in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
#'
#' @param gene_angles a numeric vector. vector, holding the angles between ancestor and genes (0<=angle<=Pi)
#' @param threshold a numeric scalar. filter threshold
#' @param n_families a integer scalar. number of families
#' @param gene_to_fam a integer vector. a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
#' @param n_mask_chunks a integer scalar. number of 32 bit chunks a mask needs to encode `n_genes` genes
#' @return bit mask that will have indices of genes kept by pattern set to 1, else 0
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

    .result <- .filter_paralogs_by_pattern_subfunctionalization_rcpp(gene_angles, threshold, n_families, gene_to_fam, n_mask_chunks)
    .arguments <- c("gene_angles", "threshold", "n_genes", "n_families", "gene_to_fam", "masks", "n_mask_chunks", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$masks
}

#' Prefilters the genes for dosage effect
#'
#' as genes that are already too distant in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
#'
#' @param gene_angles a numeric vector. vector, holding the angles between ancestor and genes (0<=angle<=Pi)
#' @param threshold a numeric scalar. filter threshold
#' @param n_families a integer scalar. number of families
#' @param gene_to_fam a integer vector. a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
#' @param n_mask_chunks a integer scalar. number of 32 bit chunks a mask needs to encode `n_genes` genes
#' @return bit mask that will have indices of genes kept by pattern set to 1, else 0
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

    .result <- .filter_paralogs_by_pattern_dosage_effect_rcpp(gene_angles, threshold, n_families, gene_to_fam, n_mask_chunks)
    .arguments <- c("gene_angles", "threshold", "n_genes", "n_families", "gene_to_fam", "masks", "n_mask_chunks", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$masks
}

#' Calculates the needed size for the paralog-subsets work array
#'
#' The `detect_*` subroutines need a work array for the to be tested subsets.
#' In worst case, all need to be tried and subsets that cannot be extended will be kept as results.
#' This is the reason why the work array holds the results as well, as all subsets that are stored in the array can be results as well.
#'
#' This subroutine calculates the needed size for the work array.
#'
#' @param max_subset_size a integer scalar. maximum size that a subset must not exceed.
#' @param n_genes a integer scalar. number of genes
#' @param filtered_paralogs_mask a integer vector. Output mask with all genes disabled that did not pass the filter
#' @return a named list with elements `max_subset_size`, `work_array_size`.
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::calc_work_arr_paralog_subsets_size}.
#' @export
calc_work_arr_paralog_subsets_size <- function(max_subset_size, n_genes, filtered_paralogs_mask) {
    max_subset_size <- .tox_as_integer_scalar(max_subset_size, "max_subset_size")
    n_genes <- .tox_as_integer_scalar(n_genes, "n_genes")
    filtered_paralogs_mask <- .tox_as_integer_vector(filtered_paralogs_mask, "filtered_paralogs_mask")
    .result <- .calc_work_arr_paralog_subsets_size_rcpp(max_subset_size, n_genes, filtered_paralogs_mask)
    .arguments <- c("max_subset_size", "n_genes", "work_array_size", "filtered_paralogs_mask", "n_mask_chunks", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        max_subset_size = .result$max_subset_size,
        work_array_size = .result$work_array_size
    )
}
