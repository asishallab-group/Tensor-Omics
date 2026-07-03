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
    subroutine mask_get_first_successor_idx_c(&
            bit_mask,&
            n_bit_mask_elements,&
            idx,&
            ierr&
            ) bind(C, name="mask_get_first_successor_idx_c")
        use tox_paralog_analysis, only: mask_get_first_successor_idx
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_bit_mask_elements
            !! Size of the 1. dimension/extent of `bit_mask`
        integer(c_int), intent(in), dimension(n_bit_mask_elements), target :: bit_mask
            !! chunked mask to mark active genes
        integer(c_int), intent(out), target :: idx
            !! index of last active gene
        integer(c_int), intent(out), target :: ierr
            !! Error code


        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(bit_mask)
        M_CHECK_NON_NULL(n_bit_mask_elements)
        M_CHECK_NON_NULL(idx)


        idx = mask_get_first_successor_idx(&
            bit_mask = bit_mask&
        )

    end subroutine mask_get_first_successor_idx_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):mask_check_state(function)]]
    !| Checks the state of a bit/paralog in `bit_mask` -> .true. if 1 else .false.
    subroutine mask_check_state_c(&
            bit_mask,&
            n_bit_mask_elements,&
            i_gene,&
            state,&
            ierr&
            ) bind(C, name="mask_check_state_c")
        use tox_paralog_analysis, only: mask_check_state
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_bit_mask_elements
            !! Size of the 1. dimension/extent of `bit_mask`
        integer(c_int), intent(in), dimension(n_bit_mask_elements), target :: bit_mask
            !! chunked mask to mark active paralogs
        integer(c_int), intent(in), target :: i_gene
            !! index of paralog to be marked active
        integer(c_int), intent(out), target :: state
            !! check result
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: state_f

        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(bit_mask)
        M_CHECK_NON_NULL(n_bit_mask_elements)
        M_CHECK_NON_NULL(i_gene)
        M_CHECK_NON_NULL(state)

        call c_int_as_logical(state, state_f)
        state_f = mask_check_state(&
            bit_mask = bit_mask,&
            i_gene = i_gene&
        )
        call logical_as_c_int(state_f, state)
    end subroutine mask_check_state_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):detect_neofunctionalization(subroutine)]]
    !| Identifies neofunctionalization for genes by checking whether the difference of expression to its ancestor exceeds the threshold for the respective axis.
    subroutine detect_neofunctionalization_c(&
            ancestors,&
            n_families,&
            genes,&
            n_axes,&
            gene_to_fam,&
            n_genes,&
            thresholds,&
            neofunc,&
            ierr&
            ) bind(C, name="detect_neofunctionalization_c")
        use tox_paralog_analysis, only: detect_neofunctionalization
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_families
            !! number of vectors in `ancestors`
        integer(c_int), intent(in), target :: n_axes
            !! size of `ancestors` vector and vectors in `genes`
        integer(c_int), intent(in), target :: n_genes
            !! number of vectors in `genes`
        real(c_double), intent(in), dimension(n_axes, n_families), target :: ancestors
            !! RAP projected unit length expression vector of ancestral ortholog
        real(c_double), intent(in), dimension(n_axes, n_genes), target :: genes
            !! RAP projected unit length expression vectors of genes
        integer(c_int), intent(in), dimension(n_genes), target :: gene_to_fam
            !! mapping of gene index to family index
        real(c_double), intent(in), dimension(n_axes), target :: thresholds
            !! threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis
        integer(c_int), intent(out), dimension(n_genes, n_axes), target :: neofunc
            !! `.true.` if neofunctionalization has been detected for the respective axes
        integer(c_int), intent(out), target :: ierr
            !! error code
        logical, allocatable, dimension(:, :) :: neofunc_f

        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(ancestors)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(genes)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(gene_to_fam)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(thresholds)
        M_CHECK_NON_NULL(neofunc)

        M_ALLOCATE(neofunc_f(n_genes, n_axes))
        call detect_neofunctionalization(&
            ancestors = ancestors,&
            n_families = n_families,&
            genes = genes,&
            n_axes = n_axes,&
            gene_to_fam = gene_to_fam,&
            n_genes = n_genes,&
            thresholds = thresholds,&
            neofunc = neofunc_f,&
            ierr = ierr&
        )
        call logical_as_c_int(neofunc_f, neofunc)
    end subroutine detect_neofunctionalization_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):detect_dosage_effect(subroutine)]]
    !| Identifies subsets of paralogs with small angle to the `ancestor` (max_angle) and sum to a magnitude significantly exceeding `norm(ancestor)` (gain)
    subroutine detect_dosage_effect_c(&
            ancestor,&
            genes,&
            n_genes,&
            n_dims,&
            filtered_paralogs_mask,&
            n_mask_chunks,&
            n_results,&
            max_subset_size,&
            work_arr_paralog_subsets,&
            n_paralog_subsets,&
            active_mask,&
            temp_paralog_vector,&
            ierr,&
            max_angle,&
            gain_gamma&
            ) bind(C, name="detect_dosage_effect_c")
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
            !! DM_FROM(work_array_size, calc_work_arr_paralog_subsets_size, tox_paralog_analysis, AUTO)
        real(c_double), intent(in), dimension(n_dims), target :: ancestor
            !! expression vector of ancestral ortholog
        real(c_double), intent(in), dimension(n_dims, n_genes), target :: genes
            !! expression vectors of genes
        integer(c_int), intent(in), dimension(n_mask_chunks), target :: filtered_paralogs_mask
            !! bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
            !! DM_FROM(masks(:, family_idx), filter_paralogs_by_pattern, tox_paralog_analysis, JUST_INFO)
        integer(c_int), intent(out), target :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(c_int), intent(in), target :: max_subset_size
            !! maximum subset size of checked gene subsets.
            !! DM_FROM(max_subset_size, calc_work_arr_paralog_subsets_size, tox_paralog_analysis, AUTO)
        integer(c_int), intent(out), dimension(n_mask_chunks, n_paralog_subsets), target :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
            !! @endnote
            !! The first `n_results` elements will hold the results.
        integer(c_int), intent(out), dimension(n_mask_chunks), target :: active_mask
            !! working array to hold the extended subsets
        real(c_double), intent(out), dimension(n_dims), target :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(c_int), intent(out), target :: ierr
            !! error code
        real(c_double), intent(in), target :: max_angle
            !! in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
            !! The default value is `PI`.
        real(c_double), intent(in), target :: gain_gamma
            !! positive magnitude gain for dosage effect
            !! The default value is `0.1_real64`.


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


        call detect_dosage_effect(&
            ancestor = ancestor,&
            genes = genes,&
            n_genes = n_genes,&
            n_dims = n_dims,&
            filtered_paralogs_mask = filtered_paralogs_mask,&
            n_mask_chunks = n_mask_chunks,&
            n_results = n_results,&
            max_subset_size = max_subset_size,&
            work_arr_paralog_subsets = work_arr_paralog_subsets,&
            n_paralog_subsets = n_paralog_subsets,&
            active_mask = active_mask,&
            temp_paralog_vector = temp_paralog_vector,&
            ierr = ierr,&
            max_angle = max_angle,&
            gain_gamma = gain_gamma&
        )

    end subroutine detect_dosage_effect_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):detect_subfunctionalization(subroutine)]]
    !| Identifies subsets of paralogs exhibiting significant angles to the `ancestor`
    subroutine detect_subfunctionalization_c(&
            ancestor,&
            genes,&
            n_genes,&
            n_dims,&
            rdi_threshold,&
            filtered_paralogs_mask,&
            n_mask_chunks,&
            n_results,&
            max_subset_size,&
            work_arr_paralog_subsets,&
            n_paralog_subsets,&
            active_mask,&
            temp_paralog_vector,&
            paralog_norms,&
            sorted_paralog_norms_perm,&
            temp_work_array,&
            ierr&
            ) bind(C, name="detect_subfunctionalization_c")
        use tox_paralog_analysis, only: detect_subfunctionalization
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_genes
            !! number of vectors in `genes`
        integer(c_int), intent(in), target :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
            !! It is recommended to compute this argument using [[tox_paralog_analysis(module):mask_chunk_count]]'s output `count`.
        integer(c_int), intent(in), target :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        real(c_double), intent(in), dimension(n_dims), target :: ancestor
            !! expression vector of ancestral ortholog
        real(c_double), intent(in), dimension(n_dims, n_genes), target :: genes
            !! expression vectors of genes
        real(c_double), intent(in), target :: rdi_threshold
            !! max allowed residual distance from `ancestor`
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
            !! The first `n_results` elements will hold the results.
        integer(c_int), intent(out), dimension(n_mask_chunks), target :: active_mask
            !! working array to hold the extended subsets
        real(c_double), intent(out), dimension(n_dims), target :: temp_paralog_vector
            !! vector used for pruning subsets
        real(c_double), intent(in), dimension(n_genes), target :: paralog_norms
            !! needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` function from `f42_utils` function for this)
        integer(c_int), intent(in), dimension(n_genes), target :: sorted_paralog_norms_perm
            !! needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor
        real(c_double), intent(out), dimension(n_genes), target :: temp_work_array
            !! needed for efficient check of minimum value after a certain index
        integer(c_int), intent(out), target :: ierr
            !! error code


        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(ancestor)
        M_CHECK_NON_NULL(genes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(rdi_threshold)
        M_CHECK_NON_NULL(filtered_paralogs_mask)
        M_CHECK_NON_NULL(n_mask_chunks)
        M_CHECK_NON_NULL(n_results)
        M_CHECK_NON_NULL(max_subset_size)
        M_CHECK_NON_NULL(work_arr_paralog_subsets)
        M_CHECK_NON_NULL(n_paralog_subsets)
        M_CHECK_NON_NULL(active_mask)
        M_CHECK_NON_NULL(temp_paralog_vector)
        M_CHECK_NON_NULL(paralog_norms)
        M_CHECK_NON_NULL(sorted_paralog_norms_perm)
        M_CHECK_NON_NULL(temp_work_array)


        call detect_subfunctionalization(&
            ancestor = ancestor,&
            genes = genes,&
            n_genes = n_genes,&
            n_dims = n_dims,&
            rdi_threshold = rdi_threshold,&
            filtered_paralogs_mask = filtered_paralogs_mask,&
            n_mask_chunks = n_mask_chunks,&
            n_results = n_results,&
            max_subset_size = max_subset_size,&
            work_arr_paralog_subsets = work_arr_paralog_subsets,&
            n_paralog_subsets = n_paralog_subsets,&
            active_mask = active_mask,&
            temp_paralog_vector = temp_paralog_vector,&
            paralog_norms = paralog_norms,&
            sorted_paralog_norms_perm = sorted_paralog_norms_perm,&
            temp_work_array = temp_work_array,&
            ierr = ierr&
        )

    end subroutine detect_subfunctionalization_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):detect_patterns(subroutine)]]
    !| Identifies subsets of paralogs where dosage effect or subfunctionalization applies, depending on `pattern`
    subroutine detect_patterns_c(&
            ancestor,&
            genes,&
            n_genes,&
            n_dims,&
            pattern_mode,&
            filtered_paralogs_mask,&
            n_mask_chunks,&
            n_results,&
            max_subset_size,&
            work_arr_paralog_subsets,&
            n_paralog_subsets,&
            active_mask,&
            temp_paralog_vector,&
            dosage_max_angle,&
            dosage_gain_gamma,&
            subfunc_rdi_threshold,&
            subfunc_paralog_norms,&
            subfunc_sorted_paralog_norms_perm,&
            subfunc_temp_work_array,&
            ierr&
            ) bind(C, name="detect_patterns_c")
        use tox_paralog_analysis, only: detect_patterns
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_genes
            !! number of vectors in `genes`
        integer(c_int), intent(in), target :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
            !! It is recommended to compute this argument using [[tox_paralog_analysis(module):mask_chunk_count]]'s output `count`.
        integer(c_int), intent(in), target :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        real(c_double), intent(in), dimension(n_dims), target :: ancestor
            !! expression vector of ancestral ortholog
        real(c_double), intent(in), dimension(n_dims, n_genes), target :: genes
            !! expression vectors of genes
        integer(c_int), intent(in), target :: pattern_mode
            !! used pattern for detection
            !! 
            !! |         Mode         |                              Value                              |
            !! |----------------------|-----------------------------------------------------------------|
            !! |    Dosage Effect     |  [[tox_paralog_analysis(module):MODE_DOSAGE_PATTERN(variable)]] |
            !! | Subfunctionalization | [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]] |
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
            !! The first `n_results` elements will hold the results.
        integer(c_int), intent(out), dimension(n_mask_chunks), target :: active_mask
            !! working array to hold the extended subsets
        real(c_double), intent(out), dimension(n_dims), target :: temp_paralog_vector
            !! vector used for pruning subsets
        real(c_double), intent(in), target :: dosage_max_angle
            !! in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
            !! The default value is `PI`.
        real(c_double), intent(in), target :: dosage_gain_gamma
            !! in dosage mode required positive magnitude gain for dosage
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: subfunc_rdi_threshold
            !! max allowed residual distance from `ancestor`
            !! This optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].
        real(c_double), intent(in), dimension(n_genes), target :: subfunc_paralog_norms
            !! in subfunctionalization mode needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` from `f42_utils` function for this)
            !! This optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].
        integer(c_int), intent(in), dimension(n_genes), target :: subfunc_sorted_paralog_norms_perm
            !! in subfunctionalization mode needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor
            !! This optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].
        real(c_double), intent(out), dimension(n_genes), target :: subfunc_temp_work_array
            !! in subfunctionalization mode needed for efficient check of minimum value after a certain index
            !! This optional argument needs to be passed if used mode (`pattern_mode`) is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].
        integer(c_int), intent(out), target :: ierr
            !! error code


        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(ancestor)
        M_CHECK_NON_NULL(genes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(pattern_mode)
        M_CHECK_NON_NULL(filtered_paralogs_mask)
        M_CHECK_NON_NULL(n_mask_chunks)
        M_CHECK_NON_NULL(n_results)
        M_CHECK_NON_NULL(max_subset_size)
        M_CHECK_NON_NULL(work_arr_paralog_subsets)
        M_CHECK_NON_NULL(n_paralog_subsets)
        M_CHECK_NON_NULL(active_mask)
        M_CHECK_NON_NULL(temp_paralog_vector)


        if (pattern_mode_int_f == MODE_SUBFUNC_PATTERN) then
            call detect_patterns(&
                ancestor = ancestor,&
                genes = genes,&
                n_genes = n_genes,&
                n_dims = n_dims,&
                pattern_mode = pattern_mode,&
                filtered_paralogs_mask = filtered_paralogs_mask,&
                n_mask_chunks = n_mask_chunks,&
                n_results = n_results,&
                max_subset_size = max_subset_size,&
                work_arr_paralog_subsets = work_arr_paralog_subsets,&
                n_paralog_subsets = n_paralog_subsets,&
                active_mask = active_mask,&
                temp_paralog_vector = temp_paralog_vector,&
                dosage_max_angle = dosage_max_angle,&
                dosage_gain_gamma = dosage_gain_gamma,&
                ierr = ierr,&
                subfunc_rdi_threshold = subfunc_rdi_threshold,&
                subfunc_paralog_norms = subfunc_paralog_norms,&
                subfunc_sorted_paralog_norms_perm = subfunc_sorted_paralog_norms_perm,&
                subfunc_temp_work_array = subfunc_temp_work_array&
            )
        else
            call detect_patterns(&
                ancestor = ancestor,&
                genes = genes,&
                n_genes = n_genes,&
                n_dims = n_dims,&
                pattern_mode = pattern_mode,&
                filtered_paralogs_mask = filtered_paralogs_mask,&
                n_mask_chunks = n_mask_chunks,&
                n_results = n_results,&
                max_subset_size = max_subset_size,&
                work_arr_paralog_subsets = work_arr_paralog_subsets,&
                n_paralog_subsets = n_paralog_subsets,&
                active_mask = active_mask,&
                temp_paralog_vector = temp_paralog_vector,&
                dosage_max_angle = dosage_max_angle,&
                dosage_gain_gamma = dosage_gain_gamma,&
                ierr = ierr&
            )
        end if

    end subroutine detect_patterns_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):mask_chunk_count(subroutine)]]
    !| This subroutine easily determines the needed chunk count for subset bit masks, as an integer has only 32 bits.
    subroutine mask_chunk_count_c(&
            n_genes,&
            count,&
            ierr&
            ) bind(C, name="mask_chunk_count_c")
        use tox_paralog_analysis, only: mask_chunk_count
        use tox_paralog_analysis

        integer(c_int), intent(in), target :: n_genes
            !! number of genes
        integer(c_int), intent(out), target :: count
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
            !! 
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
        integer(c_int), intent(out), target :: ierr
            !! Error code


        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(count)


        call mask_chunk_count(&
            n_genes = n_genes,&
            count = count&
        )

    end subroutine mask_chunk_count_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):filter_paralogs_by_pattern_subfunctionalization(subroutine)]]
    !| This subroutine prefilters the genes for subfunctionalization,
    !| as genes that are already too close in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
    subroutine filter_paralogs_by_pattern_subfunctionalization_c(&
            gene_angles,&
            threshold,&
            n_genes,&
            n_families,&
            gene_to_fam,&
            masks,&
            n_mask_chunks,&
            ierr&
            ) bind(C, name="filter_paralogs_by_pattern_subfunctionalization_c")
        use tox_paralog_analysis, only: filter_paralogs_by_pattern_subfunctionalization
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_genes
            !! number of genes
        integer(c_int), intent(in), target :: n_families
            !! number of families
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
            !! It is *VERY IMPORTANT* to compute this argument using [[tox_paralog_analysis(module):mask_chunk_count]]'s output `count`.
            !! | Argument here | Argument there |
            !! |---------------|----------------|
            !! |    n_genes    |    n_genes     |
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


        call filter_paralogs_by_pattern_subfunctionalization(&
            gene_angles = gene_angles,&
            threshold = threshold,&
            n_genes = n_genes,&
            n_families = n_families,&
            gene_to_fam = gene_to_fam,&
            masks = masks,&
            n_mask_chunks = n_mask_chunks,&
            ierr = ierr&
        )

    end subroutine filter_paralogs_by_pattern_subfunctionalization_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):filter_paralogs_by_pattern_dosage_effect(subroutine)]]
    !| This subroutine prefilters the genes for dosage effect,
    !| as genes that are already too distant in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
    subroutine filter_paralogs_by_pattern_dosage_effect_c(&
            gene_angles,&
            threshold,&
            n_genes,&
            n_families,&
            gene_to_fam,&
            masks,&
            n_mask_chunks,&
            ierr&
            ) bind(C, name="filter_paralogs_by_pattern_dosage_effect_c")
        use tox_paralog_analysis, only: filter_paralogs_by_pattern_dosage_effect
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_genes
            !! number of genes
        integer(c_int), intent(in), target :: n_families
            !! number of families
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
            !! It is *VERY IMPORTANT* to compute this argument using [[tox_paralog_analysis(module):mask_chunk_count]]'s output `count`.
            !! | Argument here | Argument there |
            !! |---------------|----------------|
            !! |    n_genes    |    n_genes     |
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


        call filter_paralogs_by_pattern_dosage_effect(&
            gene_angles = gene_angles,&
            threshold = threshold,&
            n_genes = n_genes,&
            n_families = n_families,&
            gene_to_fam = gene_to_fam,&
            masks = masks,&
            n_mask_chunks = n_mask_chunks,&
            ierr = ierr&
        )

    end subroutine filter_paralogs_by_pattern_dosage_effect_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):filter_paralogs_by_pattern(subroutine)]]
    !| This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
    subroutine filter_paralogs_by_pattern_c(&
            pattern_mode,&
            gene_angles,&
            threshold,&
            n_genes,&
            n_families,&
            gene_to_fam,&
            masks,&
            n_mask_chunks,&
            ierr&
            ) bind(C, name="filter_paralogs_by_pattern_c")
        use tox_paralog_analysis, only: filter_paralogs_by_pattern
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_genes
            !! number of genes
        integer(c_int), intent(in), target :: n_families
            !! number of families
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
            !! It is *VERY IMPORTANT* to compute this argument using [[tox_paralog_analysis(module):mask_chunk_count]]'s output `count`.
            !! | Argument here | Argument there |
            !! |---------------|----------------|
            !! |    n_genes    |    n_genes     |
        integer(c_int), intent(in), target :: pattern_mode
            !! used pattern for detection
            !! 
            !! |         Mode         |                              Value                              |
            !! |----------------------|-----------------------------------------------------------------|
            !! |    Dosage Effect     |  [[tox_paralog_analysis(module):MODE_DOSAGE_PATTERN(variable)]] |
            !! | Subfunctionalization | [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]] |
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
        M_CHECK_NON_NULL(pattern_mode)
        M_CHECK_NON_NULL(gene_angles)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(gene_to_fam)
        M_CHECK_NON_NULL(masks)
        M_CHECK_NON_NULL(n_mask_chunks)


        call filter_paralogs_by_pattern(&
            pattern_mode = pattern_mode,&
            gene_angles = gene_angles,&
            threshold = threshold,&
            n_genes = n_genes,&
            n_families = n_families,&
            gene_to_fam = gene_to_fam,&
            masks = masks,&
            n_mask_chunks = n_mask_chunks,&
            ierr = ierr&
        )

    end subroutine filter_paralogs_by_pattern_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):calc_work_arr_paralog_subsets_size(subroutine)]]
    !| The `detect_*` subroutines need a work array for the to be tested subsets.
    !| In worst case, all need to be tried and subsets that cannot be extended will be kept as results.
    !| This is the reason why the work array holds the results as well, as all subsets that are stored in the array can be results as well.
    !| 
    !| This subroutine calculates the needed size for the work array.
    subroutine calc_work_arr_paralog_subsets_size_c(&
            max_subset_size,&
            n_genes,&
            work_array_size,&
            filtered_paralogs_mask,&
            n_mask_chunks,&
            ierr&
            ) bind(C, name="calc_work_arr_paralog_subsets_size_c")
        use tox_paralog_analysis, only: calc_work_arr_paralog_subsets_size
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
            !! It is recommended to compute this argument using [[tox_paralog_analysis(module):mask_chunk_count]]'s output `count`.
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
            !! DM_FROM(masks(:, family_idx), filter_paralogs_by_pattern, tox_paralog_analysis, JUST_INFO)
        integer(c_int), intent(out), target :: ierr
            !! Error code


        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(max_subset_size)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(work_array_size)
        M_CHECK_NON_NULL(filtered_paralogs_mask)
        M_CHECK_NON_NULL(n_mask_chunks)


        call calc_work_arr_paralog_subsets_size(&
            max_subset_size = max_subset_size,&
            n_genes = n_genes,&
            work_array_size = work_array_size,&
            filtered_paralogs_mask = filtered_paralogs_mask,&
            n_mask_chunks = n_mask_chunks,&
            ierr = ierr&
        )

    end subroutine calc_work_arr_paralog_subsets_size_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):mask_set_state(subroutine)]]
    !| Sets the state of a bit/gene in `bit_mask`
    subroutine mask_set_state_c(&
            bit_mask,&
            n_bit_mask_elements,&
            i_gene,&
            state,&
            ierr&
            ) bind(C, name="mask_set_state_c")
        use tox_paralog_analysis, only: mask_set_state
        use tox_paralog_analysis
        integer(c_int), intent(in), target :: n_bit_mask_elements
            !! Size of the 1. dimension/extent of `bit_mask`
        integer(c_int), intent(out), dimension(n_bit_mask_elements), target :: bit_mask
            !! chunked mask to mark active paralogs
        integer(c_int), intent(in), target :: i_gene
            !! index of paralog to be marked active
        integer(c_int), intent(in), target :: state
            !! state the bit should be set to
        integer(c_int), intent(out), target :: ierr
            !! error code
        logical :: state_f

        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(bit_mask)
        M_CHECK_NON_NULL(n_bit_mask_elements)
        M_CHECK_NON_NULL(i_gene)
        M_CHECK_NON_NULL(state)

        call c_int_as_logical(state, state_f)
        call mask_set_state(&
            bit_mask = bit_mask,&
            i_gene = i_gene,&
            state = state_f,&
            ierr = ierr&
        )

    end subroutine mask_set_state_c

end module tox_paralog_analysis_c
#endif