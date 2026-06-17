#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[tox_paralog_analysis(module)]]
module tox_paralog_analysis_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    implicit none
contains

    !> summary: C-wrapper for [[tox_paralog_analysis(module):mask_get_first_successor_idx(function)]]
    !| Helper function that returns the index after the last active gene in `bit_mask`, so the first succeeding gene.
    subroutine mask_get_first_successor_idx_c(bit_mask, n_bit_mask_elements, idx, ierr) bind(C, name="mask_get_first_successor_idx_c")
        use tox_paralog_analysis, only: mask_get_first_successor_idx
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_bit_mask_elements
            !!  Size of the 1. dimension/extent of `bit_mask`
        integer(c_int), intent(in), dimension(n_bit_mask_elements), target :: bit_mask
            !! chunked mask to mark active genes
        integer(c_int), intent(out), target :: idx
            !! index of last active gene
        integer(c_int), intent(out), target :: ierr
            !!  Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(bit_mask)
        M_CHECK_NON_NULL(n_bit_mask_elements)
        M_CHECK_NON_NULL(idx)
        idx = mask_get_first_successor_idx(bit_mask)
    end subroutine mask_get_first_successor_idx_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):detect_dosage_effect(subroutine)]]
    !| Identifies subsets of paralogs with small angle to the `ancestor` (max_angle) and sum to a magnitude significantly exceeding `norm(ancestor)` (gain)
    subroutine detect_dosage_effect_c(ancestor, genes, n_genes, n_dims, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, ierr, max_angle, gain_gamma) bind(C, name="detect_dosage_effect_c")
        use tox_paralog_analysis, only: detect_dosage_effect
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_genes
            !! number of vectors in `genes`
        integer(c_int), intent(in), target :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
            !! It is recommended to compute this argument using [[tox_paralog_analysis(module):mask_chunk_count]]'s output `count`.
        integer(c_int), intent(in), target :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`.
            !! It is *VERY IMPORTANT* to compute this argument using [[tox_paralog_analysis(module):calc_work_arr_paralog_subsets_size]]'s output `work_array_size`.
        real(c_double), intent(in), dimension(n_dims), target :: ancestor
            !! expression vector of ancestral ortholog
        real(c_double), intent(in), dimension(n_dims, n_genes), target :: genes
            !! expression vectors of genes
        integer(c_int), intent(in), dimension(n_mask_chunks), target :: filtered_paralogs_mask
            !! bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
            !! It is recommended to compute this argument using [[tox_paralog_analysis(module):filter_paralogs_by_pattern]]'s output `masks(:, family_idx)`.
        integer(c_int), intent(out), target :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(c_int), intent(in), target :: max_subset_size
            !! maximum subset size of checked gene subsets.
            !! It is *VERY IMPORTANT* to compute this argument using [[tox_paralog_analysis(module):calc_work_arr_paralog_subsets_size]]'s output `max_subset_size`.
        integer(c_int), intent(out), dimension(n_mask_chunks, n_paralog_subsets), target :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
            !! @endnote
        integer(c_int), intent(out), dimension(n_mask_chunks), target :: active_mask
            !! working array to hold the extended subsets
        real(c_double), intent(out), dimension(n_dims), target :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(c_int), intent(out), target :: ierr
            !! error code
        real(c_double), intent(in), target :: max_angle
            !! in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned, default is Pi
        real(c_double), intent(in), target :: gain_gamma
            !! positive magnitude gain for dosage effect, default 0.1
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(ancestor)
        M_CHECK_NON_NULL(genes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(filtered_paralogs_mask)
        M_CHECK_NON_NULL(n_mask_chunks)
        M_CHECK_NON_NULL(n_results)
        M_CHECK_NON_NULL(max_subset_size)
        M_CHECK_NON_NULL(work_arr_paralog_subsets)
        M_CHECK_NON_NULL(n_paralog_subsets)
        M_CHECK_NON_NULL(active_mask)
        M_CHECK_NON_NULL(temp_paralog_vector)
        M_CHECK_NON_NULL(max_angle)
        M_CHECK_NON_NULL(gain_gamma)
        call detect_dosage_effect(ancestor, genes, n_genes, n_dims, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, ierr, max_angle, gain_gamma)
    end subroutine detect_dosage_effect_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):filter_paralogs_by_pattern_dosage_effect(subroutine)]]
    !| This subroutine prefilters the genes for dosage effect,
    !| as genes that are already too distant in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
    subroutine filter_paralogs_by_pattern_dosage_effect_c(gene_angles, threshold, n_genes, n_families, gene_to_fam, masks, n_mask_chunks, ierr) bind(C, name="filter_paralogs_by_pattern_dosage_effect_c")
        use tox_paralog_analysis, only: filter_paralogs_by_pattern_dosage_effect
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_genes
            !! number of genes
        integer(c_int), intent(in), target :: n_families
            !! number of families
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        real(c_double), intent(in), dimension(n_genes), target :: gene_angles
            !! vector, holding the angles between ancestor and genes (0<=angle<=Pi)
        real(c_double), intent(in), target :: threshold
            !! filter threshold
        integer(c_int), intent(in), dimension(n_genes), target :: gene_to_fam
            !! a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
        integer(c_int), intent(out), dimension(n_mask_chunks, n_families), target :: masks
            !! bit mask that will have indices of genes kept by pattern set to 1, else 0
        integer(c_int), intent(out), target :: ierr
            !! error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(gene_angles)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(gene_to_fam)
        M_CHECK_NON_NULL(masks)
        M_CHECK_NON_NULL(n_mask_chunks)
        call filter_paralogs_by_pattern_dosage_effect(gene_angles, threshold, n_genes, n_families, gene_to_fam, masks, n_mask_chunks, ierr)
    end subroutine filter_paralogs_by_pattern_dosage_effect_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):calc_work_arr_paralog_subsets_size(subroutine)]]
    !| The `detect_*` subroutines need a work array for the to be tested subsets.
    !| In worst case, all need to be tried and subsets that cannot be extended will be kept as results.
    !| This is the reason why the work array holds the results as well, as all subsets that are stored in the array can be results as well.
    !| 
    !| This subroutine calculates the needed size for the work array.
    subroutine calc_work_arr_paralog_subsets_size_c(max_subset_size, n_genes, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr) bind(C, name="calc_work_arr_paralog_subsets_size_c")
        use tox_paralog_analysis, only: calc_work_arr_paralog_subsets_size
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        integer(c_int), intent(inout), target :: max_subset_size
            !! maximum size that a subset must not exceed.
            !! @warning
            !! If the desired size is too large and leads to an integer overflow, `max_subset_size` will be set to the maximum valid size.
            !! 
            !! Also, size will be set to number of genes in `filtered_paralogs_mask` if larger.
            !! @endwarning
        integer(c_int), intent(in), target :: n_genes
            !! number of genes
        integer(c_int), intent(out), target :: work_array_size
            !! The calculated needed work array size in absolute worst case scenario. Look into source for details.
        integer(c_int), intent(in), dimension(n_mask_chunks), target :: filtered_paralogs_mask
            !! Output mask with all genes disabled that did not pass the filter
            !! It is recommended to compute this argument using [[tox_paralog_analysis(module):filter_paralogs_by_pattern]]'s output `masks(:, family_idx)`.
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(max_subset_size)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(work_array_size)
        M_CHECK_NON_NULL(filtered_paralogs_mask)
        M_CHECK_NON_NULL(n_mask_chunks)
        call calc_work_arr_paralog_subsets_size(max_subset_size, n_genes, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
    end subroutine calc_work_arr_paralog_subsets_size_c

end module tox_paralog_analysis_c
#endif