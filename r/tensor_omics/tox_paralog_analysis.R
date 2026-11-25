# Generated. Do not edit.

#' Identifies neofunctionalization for genes by checking whether the difference of expression to its ancestor exceeds the threshold for the respective axis
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::detect_neofunctionalization}, whose argument names
#' are the ones an error message reports.
#'
#' @param ancestors a numeric matrix. RAP projected unit length expression vector of ancestral ortholog
#'   The minimum valid value is `-1.0`.
#'   The maximum valid value is `1.0`.
#' @param genes a numeric matrix. RAP projected unit length expression vectors of genes
#'   The minimum valid value is `-1.0`.
#'   The maximum valid value is `1.0`.
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @param thresholds a numeric vector. threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis
#'   The minimum valid value is `-1.0`.
#'   The maximum valid value is `1.0`.
#' @return a logical matrix. `TRUE` if neofunctionalization has been detected for the respective axes, always `FALSE` for unassigned genes
#' @export
detect_neofunctionalization <- function(ancestors, genes, gene_to_fam, thresholds) {
    ancestors <- .tox_as_double_matrix(ancestors, "ancestors")
    genes <- .tox_as_double_matrix(genes, "genes")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    thresholds <- .tox_as_double_vector(thresholds, "thresholds")
    if (dim(genes)[1] != dim(ancestors)[1])
        .tox_shape_error("genes", 1L, dim(genes)[1], "'ancestors' implies n_axes == ", dim(ancestors)[1], "n_axes, n_genes")
    if (length(gene_to_fam) != dim(genes)[2])
        .tox_shape_error("gene_to_fam", 1L, length(gene_to_fam), "'genes' implies n_genes == ", dim(genes)[2], "n_genes")
    if (length(thresholds) != dim(ancestors)[1])
        .tox_shape_error("thresholds", 1L, length(thresholds), "'ancestors' implies n_axes == ", dim(ancestors)[1], "n_axes")

    .result <- .Call("detect_neofunctionalization_call", ancestors, genes, gene_to_fam, thresholds)
    .arguments <- c("ancestors", "n_families", "genes", "n_axes", "gene_to_fam", "n_genes", "thresholds", "neofunc", "ierr")
    .sources <- c(NA_character_, "ancestors", NA_character_, "ancestors", NA_character_, "genes", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$neofunc
}

#' Identifies subsets of paralogs matching this pattern
#'
#' A subset is grown from its genes in index order: it starts from one kept gene and is extended
#' only by kept genes after its last one. A single gene is never a result; a result has two genes
#' at least, and is not extended further. So which subsets are found depends on the genes' order,
#' and a result's superset is found too when the superset's own prefix did not meet the criterion
#' (issue #244).
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::detect_dosage_effect}, whose argument names
#' are the ones an error message reports.
#'
#' @param ancestor a numeric vector. expression vector of ancestral ortholog. Dosage effect only: it must not be all
#'   zeros, which have no direction to measure a subset's angle against, and such an
#'   ancestor is refused with `ERR_INVALID_INPUT`; subfunctionalization takes it.
#' @param genes a numeric matrix. expression vectors of genes. Dosage effect only: a subset whose summed vector has
#'   every component exactly 0.0, bit for bit, has no direction either, and is dropped
#'   from the search like one past `max_angle`. A sum that only rounds to near zero is
#'   not dropped: where its length times the ancestor's is at most 1e-12, the angle's
#'   computation refuses the whole detection with `ERR_DIVISION_BY_ZERO` (pending issue
#'   #223).
#' @param filtered_paralogs_bitmask a integer vector. bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine.
#'   A bit set past `n_genes` is refused with `ERR_INVALID_INPUT`; the sizing routine
#'   ignores such bits.
#' @param max_subset_size a integer scalar. maximum subset size of checked gene subsets. Too large a value is capped to the
#'   maximum valid size. The bindings cap it automatically while sizing the work
#'   array; a Fortran caller caps it by calling
#'   \code{\link{calc_paralog_subsets_bitmasks_size}} first.
#'   Zero is in range and means there is no subset to check; the sizing routine reports it
#'   when the filter kept no gene. With fewer than two kept genes there is no pair to form:
#'   the sizing then reports a work array of no subsets (and a size of 1 for one kept
#'   gene), which this routine takes, finding nothing (`n_results` 0).
#'   The minimum valid value is `0`.
#' @param max_angle a numeric scalar. maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
#'   The default value is `4.0*atan(1.0)`.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `PI`.
#' @param gain_gamma a numeric scalar. positive magnitude gain for dosage effect
#'   The default value is `0.1`.
#'   The minimum valid value is `above(0.0)`.
#' @return a named list with elements:
#'   \item{n_results}{a integer scalar. number of resulting subsets. They are stored as the first `n_results` elements of `paralog_subsets_bitmasks`}
#'   \item{paralog_subsets_bitmasks}{a integer matrix. working array to hold bitmask encoded subsets for detection.
#'     Each column is a bit mask of the `n_genes` genes, laid out as \code{f42_bit_masks_impl} describes.}
#' @export
detect_dosage_effect <- function(ancestor, genes, filtered_paralogs_bitmask, max_subset_size, max_angle = 3.141592653589793, gain_gamma = 0.1) {
    ancestor <- .tox_as_double_vector(ancestor, "ancestor")
    genes <- .tox_as_double_matrix(genes, "genes")
    filtered_paralogs_bitmask <- .tox_as_integer_vector(filtered_paralogs_bitmask, "filtered_paralogs_bitmask")
    max_subset_size <- .tox_as_integer_scalar(max_subset_size, "max_subset_size")
    max_angle <- .tox_as_double_scalar(max_angle, "max_angle")
    gain_gamma <- .tox_as_double_scalar(gain_gamma, "gain_gamma")
    .arguments <- c("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_bitmask", "n_mask_chunks", "n_results", "max_subset_size", "paralog_subsets_bitmasks", "n_paralog_subsets", "max_angle", "gain_gamma", "ierr")
    .sources <- c(NA_character_, NA_character_, "genes", "ancestor", NA_character_, "filtered_paralogs_bitmask", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .calc_paralog_subsets_bitmasks_size_result <- tryCatch(calc_paralog_subsets_bitmasks_size(max_subset_size = max_subset_size, n_genes = dim(genes)[2], filtered_paralogs_bitmask = filtered_paralogs_bitmask), tox_error = function(e) .tox_consumer_error(e, c(max_subset_size = 8L, n_genes = 3L, filtered_paralogs_bitmask = 5L), .arguments, .sources))
    max_subset_size <- .calc_paralog_subsets_bitmasks_size_result$max_subset_size
    n_paralog_subsets <- .calc_paralog_subsets_bitmasks_size_result$work_array_size

    if (dim(genes)[1] != length(ancestor))
        .tox_shape_error("genes", 1L, dim(genes)[1], "'ancestor' implies n_dims == ", length(ancestor), "n_dims, n_genes")

    .result <- .Call("detect_dosage_effect_call", ancestor, genes, filtered_paralogs_bitmask, max_subset_size, n_paralog_subsets, max_angle, gain_gamma)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        n_results = .result$n_results,
        paralog_subsets_bitmasks = .result$paralog_subsets_bitmasks
    )
}

#' Identifies subsets of paralogs matching this pattern
#'
#' A subset is grown from its genes in index order: it starts from one kept gene and is extended
#' only by kept genes after its last one. A single gene is never a result; a result has two genes
#' at least, and is not extended further. So which subsets are found depends on the genes' order,
#' and a result's superset is found too when the superset's own prefix did not meet the criterion
#' (issue #244).
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::detect_subfunctionalization}, whose argument names
#' are the ones an error message reports.
#'
#' @param ancestor a numeric vector. expression vector of ancestral ortholog. Dosage effect only: it must not be all
#'   zeros, which have no direction to measure a subset's angle against, and such an
#'   ancestor is refused with `ERR_INVALID_INPUT`; subfunctionalization takes it.
#' @param genes a numeric matrix. expression vectors of genes. Dosage effect only: a subset whose summed vector has
#'   every component exactly 0.0, bit for bit, has no direction either, and is dropped
#'   from the search like one past `max_angle`. A sum that only rounds to near zero is
#'   not dropped: where its length times the ancestor's is at most 1e-12, the angle's
#'   computation refuses the whole detection with `ERR_DIVISION_BY_ZERO` (pending issue
#'   #223).
#' @param filtered_paralogs_bitmask a integer vector. bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine.
#'   A bit set past `n_genes` is refused with `ERR_INVALID_INPUT`; the sizing routine
#'   ignores such bits.
#' @param max_subset_size a integer scalar. maximum subset size of checked gene subsets. Too large a value is capped to the
#'   maximum valid size. The bindings cap it automatically while sizing the work
#'   array; a Fortran caller caps it by calling
#'   \code{\link{calc_paralog_subsets_bitmasks_size}} first.
#'   Zero is in range and means there is no subset to check; the sizing routine reports it
#'   when the filter kept no gene. With fewer than two kept genes there is no pair to form:
#'   the sizing then reports a work array of no subsets (and a size of 1 for one kept
#'   gene), which this routine takes, finding nothing (`n_results` 0).
#'   The minimum valid value is `0`.
#' @param rdi_threshold a numeric scalar. max allowed residual distance from `ancestor`
#'   The minimum valid value is `0.0`.
#' @param paralog_norms a numeric vector. euclidean norms of the genes, used for subset pruning (`norm` from `f42_vector_impl` computes them).
#'   The routine trusts them to be the genes' own.
#'   The minimum valid value is `0.0`.
#' @param sorted_paralog_norms_perm a integer vector. ascending permutation of `paralog_norms`, for subset pruning: a subset is extended only
#'   while the smallest norm among all genes after its last one -- kept by the filter or
#'   not -- is at most the length of its residual. Only the entries' range is checked: a
#'   permutation that is not ascending gives wrong minima, and an order with a repeated
#'   entry (such as [1, 1, 1]) leaves part of the work array of minima unwritten, so the
#'   pruning reads uninitialised values.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_genes`.
#' @return a named list with elements:
#'   \item{n_results}{a integer scalar. number of resulting subsets. They are stored as the first `n_results` elements of `paralog_subsets_bitmasks`}
#'   \item{paralog_subsets_bitmasks}{a integer matrix. working array to hold bitmask encoded subsets for detection.
#'     Each column is a bit mask of the `n_genes` genes, laid out as \code{f42_bit_masks_impl} describes.}
#' @export
detect_subfunctionalization <- function(ancestor, genes, filtered_paralogs_bitmask, max_subset_size, rdi_threshold, paralog_norms, sorted_paralog_norms_perm) {
    ancestor <- .tox_as_double_vector(ancestor, "ancestor")
    genes <- .tox_as_double_matrix(genes, "genes")
    filtered_paralogs_bitmask <- .tox_as_integer_vector(filtered_paralogs_bitmask, "filtered_paralogs_bitmask")
    max_subset_size <- .tox_as_integer_scalar(max_subset_size, "max_subset_size")
    rdi_threshold <- .tox_as_double_scalar(rdi_threshold, "rdi_threshold")
    paralog_norms <- .tox_as_double_vector(paralog_norms, "paralog_norms")
    sorted_paralog_norms_perm <- .tox_as_integer_vector(sorted_paralog_norms_perm, "sorted_paralog_norms_perm")
    .arguments <- c("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_bitmask", "n_mask_chunks", "n_results", "max_subset_size", "paralog_subsets_bitmasks", "n_paralog_subsets", "rdi_threshold", "paralog_norms", "sorted_paralog_norms_perm", "ierr")
    .sources <- c(NA_character_, NA_character_, "genes", "ancestor", NA_character_, "filtered_paralogs_bitmask", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .calc_paralog_subsets_bitmasks_size_result <- tryCatch(calc_paralog_subsets_bitmasks_size(max_subset_size = max_subset_size, n_genes = dim(genes)[2], filtered_paralogs_bitmask = filtered_paralogs_bitmask), tox_error = function(e) .tox_consumer_error(e, c(max_subset_size = 8L, n_genes = 3L, filtered_paralogs_bitmask = 5L), .arguments, .sources))
    max_subset_size <- .calc_paralog_subsets_bitmasks_size_result$max_subset_size
    n_paralog_subsets <- .calc_paralog_subsets_bitmasks_size_result$work_array_size

    if (dim(genes)[1] != length(ancestor))
        .tox_shape_error("genes", 1L, dim(genes)[1], "'ancestor' implies n_dims == ", length(ancestor), "n_dims, n_genes")
    if (length(paralog_norms) != dim(genes)[2])
        .tox_shape_error("paralog_norms", 1L, length(paralog_norms), "'genes' implies n_genes == ", dim(genes)[2], "n_genes")
    if (length(sorted_paralog_norms_perm) != dim(genes)[2])
        .tox_shape_error("sorted_paralog_norms_perm", 1L, length(sorted_paralog_norms_perm), "'genes' implies n_genes == ", dim(genes)[2], "n_genes")

    .result <- .Call("detect_subfunctionalization_call", ancestor, genes, filtered_paralogs_bitmask, max_subset_size, n_paralog_subsets, rdi_threshold, paralog_norms, sorted_paralog_norms_perm)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        n_results = .result$n_results,
        paralog_subsets_bitmasks = .result$paralog_subsets_bitmasks
    )
}

#' Prefilters the genes for a pattern, so genes that cannot match it are not tried as subset extensions
#'
#' This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::filter_paralogs_by_pattern_dosage_effect}, whose argument names
#' are the ones an error message reports.
#'
#' @param gene_angles a numeric vector. vector, holding the angles between ancestor and genes (0<=angle<=Pi)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `PI`.
#' @param threshold a numeric scalar. filter threshold
#' @param n_families a integer scalar. number of families
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_genes`.
#' @param gene_to_fam a integer vector. a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#' @param n_mask_chunks a integer scalar. number of 32 bit chunks a mask needs to encode `n_genes` genes
#'   The minimum valid value is `((n_genes)/32 + min(1, mod((n_genes), 32)))`.
#' @return a integer matrix. bit mask that will have the indices of genes kept by this pattern set to 1, else 0
#' @export
filter_paralogs_by_pattern_dosage_effect <- function(gene_angles, threshold, n_families, gene_to_fam, n_mask_chunks) {
    gene_angles <- .tox_as_double_vector(gene_angles, "gene_angles")
    threshold <- .tox_as_double_scalar(threshold, "threshold")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    n_mask_chunks <- .tox_as_integer_scalar(n_mask_chunks, "n_mask_chunks")
    .tox_check_extent(n_families, "n_families")
    .tox_check_extent(n_mask_chunks, "n_mask_chunks")
    if (length(gene_to_fam) != length(gene_angles))
        .tox_shape_error("gene_to_fam", 1L, length(gene_to_fam), "'gene_angles' implies n_genes == ", length(gene_angles), "n_genes")

    .result <- .Call("filter_paralogs_by_pattern_dosage_effect_call", gene_angles, threshold, n_families, gene_to_fam, n_mask_chunks)
    .arguments <- c("gene_angles", "threshold", "n_genes", "n_families", "gene_to_fam", "family_bitmasks", "n_mask_chunks", "ierr")
    .sources <- c(NA_character_, NA_character_, "gene_angles", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$family_bitmasks
}

#' Prefilters the genes for a pattern, so genes that cannot match it are not tried as subset extensions
#'
#' This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::filter_paralogs_by_pattern_subfunctionalization}, whose argument names
#' are the ones an error message reports.
#'
#' @param gene_angles a numeric vector. vector, holding the angles between ancestor and genes (0<=angle<=Pi)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `PI`.
#' @param threshold a numeric scalar. filter threshold
#' @param n_families a integer scalar. number of families
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_genes`.
#' @param gene_to_fam a integer vector. a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#' @param n_mask_chunks a integer scalar. number of 32 bit chunks a mask needs to encode `n_genes` genes
#'   The minimum valid value is `((n_genes)/32 + min(1, mod((n_genes), 32)))`.
#' @return a integer matrix. bit mask that will have the indices of genes kept by this pattern set to 1, else 0
#' @export
filter_paralogs_by_pattern_subfunctionalization <- function(gene_angles, threshold, n_families, gene_to_fam, n_mask_chunks) {
    gene_angles <- .tox_as_double_vector(gene_angles, "gene_angles")
    threshold <- .tox_as_double_scalar(threshold, "threshold")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    n_mask_chunks <- .tox_as_integer_scalar(n_mask_chunks, "n_mask_chunks")
    .tox_check_extent(n_families, "n_families")
    .tox_check_extent(n_mask_chunks, "n_mask_chunks")
    if (length(gene_to_fam) != length(gene_angles))
        .tox_shape_error("gene_to_fam", 1L, length(gene_to_fam), "'gene_angles' implies n_genes == ", length(gene_angles), "n_genes")

    .result <- .Call("filter_paralogs_by_pattern_subfunctionalization_call", gene_angles, threshold, n_families, gene_to_fam, n_mask_chunks)
    .arguments <- c("gene_angles", "threshold", "n_genes", "n_families", "gene_to_fam", "family_bitmasks", "n_mask_chunks", "ierr")
    .sources <- c(NA_character_, NA_character_, "gene_angles", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$family_bitmasks
}

#' Calculates the needed size for the paralog-subsets work array
#'
#' The `detect_*` subroutines need a work array for the subsets they try. It holds, side by side,
#' the results found so far, the subsets still to be grown and the ones grown from them, so it is
#' sized for the worst case: every subset of the first F - 1 kept genes stays a candidate, and every
#' subset that adds the last kept gene becomes a result. With F genes kept and subsets of at most m
#' genes, that is one subset of the first F - 1 genes for each size 1 to m: the sum over
#' k = 1 .. min(m, F - 1) of C(F - 1, k), which is 2**(F - 1) - 1 for m = F, and 0 for fewer
#' than two kept genes.
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis::calc_paralog_subsets_bitmasks_size}, whose argument names
#' are the ones an error message reports.
#'
#' @param max_subset_size a integer scalar. maximum size that a subset must not exceed. Zero is in range and means there is
#'   nothing to size a work array for, which is reported back as a size of zero.
#'   The minimum valid value is `0`.
#'   It comes back lowered to the number of kept genes if it is larger, and lowered further
#'   to the largest size whose work array holds at most 2**31 - 1 (2,147,483,647) subsets, the
#'   maximum valid size, where the size asked for would need more.
#' @param n_genes a integer scalar. number of genes. It is checked before `n_mask_chunks`, so a mask of too few words, none
#'   included, is reported as `n_genes` out of range.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_mask_chunks * 32`.
#' @param filtered_paralogs_bitmask a integer vector. bit mask of the genes the filter kept. Only the first `n_genes` bits are counted;
#'   a bit past them is ignored here, while the detections refuse it.
#' @return a named list with elements:
#'   \item{max_subset_size}{a integer scalar. maximum size that a subset must not exceed. Zero is in range and means there is
#'     nothing to size a work array for, which is reported back as a size of zero.
#'     The minimum valid value is `0`.
#'     It comes back lowered to the number of kept genes if it is larger, and lowered further
#'     to the largest size whose work array holds at most 2**31 - 1 (2,147,483,647) subsets, the
#'     maximum valid size, where the size asked for would need more.}
#'   \item{work_array_size}{a integer scalar. The number of subsets the work array must hold in the worst case, for the
#'     `max_subset_size` handed back}
#' @export
calc_paralog_subsets_bitmasks_size <- function(max_subset_size, n_genes, filtered_paralogs_bitmask) {
    max_subset_size <- .tox_as_integer_scalar(max_subset_size, "max_subset_size")
    n_genes <- .tox_as_integer_scalar(n_genes, "n_genes")
    filtered_paralogs_bitmask <- .tox_as_integer_vector(filtered_paralogs_bitmask, "filtered_paralogs_bitmask")
    .result <- .Call("calc_paralog_subsets_bitmasks_size_call", max_subset_size, n_genes, filtered_paralogs_bitmask)
    .arguments <- c("max_subset_size", "n_genes", "work_array_size", "filtered_paralogs_bitmask", "n_mask_chunks", "ierr")
    .sources <- c(NA_character_, NA_character_, NA_character_, NA_character_, "filtered_paralogs_bitmask", NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        max_subset_size = .result$max_subset_size,
        work_array_size = .result$work_array_size
    )
}
