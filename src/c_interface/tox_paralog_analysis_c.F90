#ifndef NO_C_INTERFACE
#include <src/macros.h>

!>  Module for C-wrappers for [[tox_paralog_analysis(module)]]
module tox_paralog_analysis_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL
contains

    !> C-wrapper for [[tox_paralog_analysis(module):mask_get_first_successor_idx(function)]]
    !| 
    !| Helper function that returns the index after the last active gene in `bit_mask`, so the first succeeding gene.
    subroutine mask_get_first_successor_idx_c(bit_mask, n_bit_mask_elements, idx, ierr) bind(C, name="mask_get_first_successor_idx_c")
        use tox_paralog_analysis, only: mask_get_first_successor_idx
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

    !> C-wrapper for [[tox_paralog_analysis(module):detect_patterns(subroutine)]]
    !| 
    !| Identifies subsets of paralogs where dosage effect or subfunctionalization applies, depending on `pattern`
    subroutine detect_patterns_c(ancestor, genes, n_genes, n_dims, pattern, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, dosage_max_angle, dosage_gain_gamma, subfunc_rdi_threshold, subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm, subfunc_temp_work_array, ierr) bind(C, name="detect_patterns_c")
        use tox_paralog_analysis, only: detect_patterns
        integer(c_int), intent(in), target :: n_genes
            !! number of vectors in `genes`
        integer(c_int), intent(in), target :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
        integer(c_int), intent(in), target :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        real(c_double), intent(in), dimension(n_dims), target :: ancestor
            !! expression vector of ancestral ortholog
        real(c_double), intent(in), dimension(n_dims, n_genes), target :: genes
            !! expression vectors of genes
        integer(c_int), intent(in), target :: pattern
            !! used pattern for detection
            !! 
            !! |       Pattern        | Value |
            !! |----------------------|-------|
            !! |    Dosage Effect     |   0   |
            !! | Subfunctionalization |   1   |
        integer(c_int), intent(in), dimension(n_mask_chunks), target :: filtered_paralogs_mask
            !! bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
        integer(c_int), intent(out), target :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(c_int), intent(in), target :: max_subset_size
            !! maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(c_int), intent(out), dimension(n_mask_chunks, n_paralog_subsets), target :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
            !! @endnote
        integer(c_int), intent(out), dimension(n_mask_chunks), target :: active_mask
            !! working array to hold the extended subsets
        real(c_double), intent(out), dimension(n_dims), target :: temp_paralog_vector
            !! vector used for pruning subsets
        real(c_double), intent(in), target :: dosage_max_angle
            !! in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned, default is Pi
        real(c_double), intent(in), target :: dosage_gain_gamma
            !! in dosage mode required positive magnitude gain for dosage, default 0.1
        real(c_double), intent(in), target :: subfunc_rdi_threshold
            !! max allowed residual distance from `ancestor`
        real(c_double), intent(in), dimension(n_genes), target :: subfunc_paralog_norms
            !! in subfunctionalization mode needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` from `f42_utils` function for this)
        integer(c_int), intent(in), dimension(n_genes), target :: subfunc_sorted_paralog_norms_perm
            !! in subfunctionalization mode needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor
            !! This optional argument needs to be passed if used mode is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].
        real(c_double), intent(out), dimension(n_genes), target :: subfunc_temp_work_array
            !! in subfunctionalization mode needed for efficient check of minimum value after a certain index
        integer(c_int), intent(out), target :: ierr
            !! error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(ancestor)
        M_CHECK_NON_NULL(genes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(pattern)
        M_CHECK_NON_NULL(filtered_paralogs_mask)
        M_CHECK_NON_NULL(n_mask_chunks)
        M_CHECK_NON_NULL(n_results)
        M_CHECK_NON_NULL(max_subset_size)
        M_CHECK_NON_NULL(work_arr_paralog_subsets)
        M_CHECK_NON_NULL(n_paralog_subsets)
        M_CHECK_NON_NULL(active_mask)
        M_CHECK_NON_NULL(temp_paralog_vector)
        M_CHECK_NON_NULL(dosage_max_angle)
        M_CHECK_NON_NULL(dosage_gain_gamma)
        M_CHECK_NON_NULL(subfunc_rdi_threshold)
        M_CHECK_NON_NULL(subfunc_paralog_norms)
        M_CHECK_NON_NULL(subfunc_sorted_paralog_norms_perm)
        M_CHECK_NON_NULL(subfunc_temp_work_array)
        call detect_patterns(ancestor, genes, n_genes, n_dims, pattern, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, dosage_max_angle, dosage_gain_gamma, subfunc_rdi_threshold, subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm, subfunc_temp_work_array, ierr)
    end subroutine detect_patterns_c

end module tox_paralog_analysis_c
#endif