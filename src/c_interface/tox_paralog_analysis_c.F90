#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: C-wrappers for [[tox_paralog_analysis(module)]]
!| Module for detecting paralog-subset expression patterns (dosage effect and subfunctionalization) relative to an ancestral ortholog.
!|
!| Candidate paralog subsets are enumerated as bitmask-encoded gene sets, built up one gene at a time
!| starting from single genes. At every extension step the candidate is scored against the pattern-specific
!| criterion (small angle plus magnitude gain for dosage effect, or bounded residual distance for
!| subfunctionalization); subsets that can no longer satisfy the criterion are pruned instead of being
!| extended further, which keeps the combinatorial subset search tractable.
module tox_paralog_analysis_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: mask_check_state_c
    public :: detect_neofunctionalization_c
    public :: detect_dosage_effect_c
    public :: detect_subfunctionalization_c
    public :: mask_chunk_count_c
    public :: filter_paralogs_by_pattern_subfunctionalization_c
    public :: filter_paralogs_by_pattern_dosage_effect_c
    public :: calc_work_arr_paralog_subsets_size_c

contains

    !> summary: C-wrapper for [[tox_paralog_analysis(module):mask_check_state(function)]]
    subroutine mask_check_state_c(&
            bit_mask,&
            n_bit_mask_elements,&
            i_gene,&
            state,&
            ierr&
        ) bind(C, name="mask_check_state_c")
        use tox_paralog_analysis, only: mask_check_state

        integer(c_int), intent(in), target :: n_bit_mask_elements
            !! number of elements in `bit_mask`
        integer(c_int), dimension(n_bit_mask_elements), intent(in), target :: bit_mask
            !! chunked mask to mark active paralogs
        integer(c_int), intent(in), target :: i_gene
            !! index of paralog to be marked active
        logical(c_bool), intent(out), target :: state
            !! check result
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: state_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bit_mask_elements)
        M_CHECK_NON_NULL(i_gene)
        M_CHECK_NON_NULL(state)
        M_CHECK_ARRAY_NON_NULL(bit_mask, n_bit_mask_elements)

        state_f = mask_check_state(&
            bit_mask = bit_mask,&
            i_gene = i_gene&
        )

        state = state_f
    end subroutine mask_check_state_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):detect_neofunctionalization(subroutine)]]
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

        integer(c_int), intent(in), target :: n_families
            !! number of vectors in `ancestors`
        integer(c_int), intent(in), target :: n_axes
            !! size of `ancestors` vector and vectors in `genes`
        integer(c_int), intent(in), target :: n_genes
            !! number of vectors in `genes`
        real(c_double), dimension(n_axes, n_families), intent(in), target :: ancestors
            !! RAP projected unit length expression vector of ancestral ortholog
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: genes
            !! RAP projected unit length expression vectors of genes
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0_int32` for unassigned genes
        real(c_double), dimension(n_axes), intent(in), target :: thresholds
            !! threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis
        logical(c_bool), dimension(n_genes, n_axes), intent(out), target :: neofunc
            !! `.true.` if neofunctionalization has been detected for the respective axes, always `.false.` for unassigned genes
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical, dimension(n_genes, n_axes) :: neofunc_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_ARRAY_NON_NULL(ancestors, n_axes * n_families)
        M_CHECK_ARRAY_NON_NULL(genes, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(thresholds, n_axes)
        M_CHECK_ARRAY_NON_NULL(neofunc, n_genes * n_axes)

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

        neofunc = neofunc_f
    end subroutine detect_neofunctionalization_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):detect_dosage_effect(subroutine)]]
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
            tmp_active_mask,&
            tmp_paralog_vector,&
            ierr,&
            max_angle,&
            gain_gamma&
        ) bind(C, name="detect_dosage_effect_c")
        use tox_paralog_analysis, only: detect_dosage_effect

        integer(c_int), intent(in), target :: n_genes
            !! number of vectors in `genes`
        integer(c_int), intent(in), target :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
        integer(c_int), intent(in), target :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`.
            !! It is *VERY IMPORTANT* to compute this argument from the `work_array_size` output produced by [[tox_paralog_analysis(module):calc_work_arr_paralog_subsets_size]].
        real(c_double), dimension(n_dims), intent(in), target :: ancestor
            !! expression vector of ancestral ortholog
        real(c_double), dimension(n_dims, n_genes), intent(in), target :: genes
            !! expression vectors of genes
        integer(c_int), dimension(n_mask_chunks), intent(in), target :: filtered_paralogs_mask
            !! bit mask with genes' indices kept by pattern set to 1, else 0.
            !! It is recommended to compute this argument from the `masks` output produced by [[tox_paralog_analysis(module):filter_paralogs_by_pattern_dosage_effect]].
        integer(c_int), intent(out), target :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(c_int), intent(in), target :: max_subset_size
            !! maximum subset size of checked gene subsets. Too large a value is capped to the
            !! maximum valid size. The interfaces cap it automatically while sizing the work
            !! array; a Fortran caller caps it by calling
            !! [[tox_paralog_analysis(module):calc_work_arr_paralog_subsets_size(subroutine)]] first.
        integer(c_int), dimension(n_mask_chunks, n_paralog_subsets), intent(out), target :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0_real64)` and represents the number of chunks
            !! @endnote
        integer(c_int), dimension(n_mask_chunks), intent(out), target :: tmp_active_mask
            !! working array to hold the extended subsets
        real(c_double), dimension(n_dims), intent(out), target :: tmp_paralog_vector
            !! vector used for pruning subsets
        integer(c_int), intent(out), target :: ierr
            !! Error code
        real(c_double), intent(in), target :: max_angle
            !! in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
            !! The default value is `4.0_real64*atan(1.0_real64)`.
        real(c_double), intent(in), target :: gain_gamma
            !! positive magnitude gain for dosage effect
            !! The default value is `0.1_real64`.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(n_mask_chunks)
        M_CHECK_NON_NULL(n_results)
        M_CHECK_NON_NULL(max_subset_size)
        M_CHECK_NON_NULL(n_paralog_subsets)
        M_CHECK_NON_NULL(max_angle)
        M_CHECK_NON_NULL(gain_gamma)
        M_CHECK_ARRAY_NON_NULL(ancestor, n_dims)
        M_CHECK_ARRAY_NON_NULL(genes, n_dims * n_genes)
        M_CHECK_ARRAY_NON_NULL(filtered_paralogs_mask, n_mask_chunks)
        M_CHECK_ARRAY_NON_NULL(work_arr_paralog_subsets, n_mask_chunks * n_paralog_subsets)
        M_CHECK_ARRAY_NON_NULL(tmp_active_mask, n_mask_chunks)
        M_CHECK_ARRAY_NON_NULL(tmp_paralog_vector, n_dims)

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
            tmp_active_mask = tmp_active_mask,&
            tmp_paralog_vector = tmp_paralog_vector,&
            ierr = ierr,&
            max_angle = max_angle,&
            gain_gamma = gain_gamma&
        )
    end subroutine detect_dosage_effect_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):detect_subfunctionalization(subroutine)]]
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
            tmp_active_mask,&
            tmp_paralog_vector,&
            paralog_norms,&
            sorted_paralog_norms_perm,&
            tmp_work_array,&
            ierr&
        ) bind(C, name="detect_subfunctionalization_c")
        use tox_paralog_analysis, only: detect_subfunctionalization

        integer(c_int), intent(in), target :: n_genes
            !! number of vectors in `genes`
        integer(c_int), intent(in), target :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
        integer(c_int), intent(in), target :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`.
            !! It is *VERY IMPORTANT* to compute this argument from the `work_array_size` output produced by [[tox_paralog_analysis(module):calc_work_arr_paralog_subsets_size]].
        real(c_double), dimension(n_dims), intent(in), target :: ancestor
            !! expression vector of ancestral ortholog
        real(c_double), dimension(n_dims, n_genes), intent(in), target :: genes
            !! expression vectors of genes
        real(c_double), intent(in), target :: rdi_threshold
            !! max allowed residual distance from `ancestor`
        integer(c_int), dimension(n_mask_chunks), intent(in), target :: filtered_paralogs_mask
            !! bit mask with genes' indices kept by pattern set to 1, else 0.
            !! It is recommended to compute this argument from the `masks` output produced by [[tox_paralog_analysis(module):filter_paralogs_by_pattern_subfunctionalization]].
        integer(c_int), intent(out), target :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(c_int), intent(in), target :: max_subset_size
            !! maximum subset size of checked gene subsets. Too large a value is capped to the
            !! maximum valid size. The interfaces cap it automatically while sizing the work
            !! array; a Fortran caller caps it by calling
            !! [[tox_paralog_analysis(module):calc_work_arr_paralog_subsets_size(subroutine)]] first.
        integer(c_int), dimension(n_mask_chunks, n_paralog_subsets), intent(out), target :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0_real64)` and represents the number of chunks
            !! @endnote
        integer(c_int), dimension(n_mask_chunks), intent(out), target :: tmp_active_mask
            !! working array to hold the extended subsets
        real(c_double), dimension(n_dims), intent(out), target :: tmp_paralog_vector
            !! vector used for pruning subsets
        real(c_double), dimension(n_genes), intent(in), target :: paralog_norms
            !! needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` function from `f42_utils` function for this)
        integer(c_int), dimension(n_genes), intent(in), target :: sorted_paralog_norms_perm
            !! needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor
        real(c_double), dimension(n_genes), intent(out), target :: tmp_work_array
            !! needed for efficient check of minimum value after a certain index
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(rdi_threshold)
        M_CHECK_NON_NULL(n_mask_chunks)
        M_CHECK_NON_NULL(n_results)
        M_CHECK_NON_NULL(max_subset_size)
        M_CHECK_NON_NULL(n_paralog_subsets)
        M_CHECK_ARRAY_NON_NULL(ancestor, n_dims)
        M_CHECK_ARRAY_NON_NULL(genes, n_dims * n_genes)
        M_CHECK_ARRAY_NON_NULL(filtered_paralogs_mask, n_mask_chunks)
        M_CHECK_ARRAY_NON_NULL(work_arr_paralog_subsets, n_mask_chunks * n_paralog_subsets)
        M_CHECK_ARRAY_NON_NULL(tmp_active_mask, n_mask_chunks)
        M_CHECK_ARRAY_NON_NULL(tmp_paralog_vector, n_dims)
        M_CHECK_ARRAY_NON_NULL(paralog_norms, n_genes)
        M_CHECK_ARRAY_NON_NULL(sorted_paralog_norms_perm, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_work_array, n_genes)

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
            tmp_active_mask = tmp_active_mask,&
            tmp_paralog_vector = tmp_paralog_vector,&
            paralog_norms = paralog_norms,&
            sorted_paralog_norms_perm = sorted_paralog_norms_perm,&
            tmp_work_array = tmp_work_array,&
            ierr = ierr&
        )
    end subroutine detect_subfunctionalization_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):mask_chunk_count(subroutine)]]
    subroutine mask_chunk_count_c(&
            n_genes,&
            count,&
            ierr&
        ) bind(C, name="mask_chunk_count_c")
        use tox_paralog_analysis, only: mask_chunk_count

        integer(c_int), intent(in), target :: n_genes
            !! number of genes
        integer(c_int), intent(out), target :: count
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
            !!
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0_real64)` and represents the number of chunks
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(count)

        call mask_chunk_count(&
            n_genes = n_genes,&
            count = count&
        )
    end subroutine mask_chunk_count_c

    !> summary: C-wrapper for [[tox_paralog_analysis(module):filter_paralogs_by_pattern_subfunctionalization(subroutine)]]
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

        integer(c_int), intent(in), target :: n_genes
            !! number of genes
        integer(c_int), intent(in), target :: n_families
            !! number of families
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        real(c_double), dimension(n_genes), intent(in), target :: gene_angles
            !! vector, holding the angles between ancestor and genes (0<=angle<=Pi)
        real(c_double), intent(in), target :: threshold
            !! filter threshold
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
        integer(c_int), dimension(n_mask_chunks, n_families), intent(out), target :: masks
            !! bit mask that will have indices of genes kept by pattern set to 1, else 0
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(n_mask_chunks)
        M_CHECK_ARRAY_NON_NULL(gene_angles, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(masks, n_mask_chunks * n_families)

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

        integer(c_int), intent(in), target :: n_genes
            !! number of genes
        integer(c_int), intent(in), target :: n_families
            !! number of families
        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        real(c_double), dimension(n_genes), intent(in), target :: gene_angles
            !! vector, holding the angles between ancestor and genes (0<=angle<=Pi)
        real(c_double), intent(in), target :: threshold
            !! filter threshold
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
        integer(c_int), dimension(n_mask_chunks, n_families), intent(out), target :: masks
            !! bit mask that will have indices of genes kept by pattern set to 1, else 0
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(n_mask_chunks)
        M_CHECK_ARRAY_NON_NULL(gene_angles, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(masks, n_mask_chunks * n_families)

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
        integer(c_int), dimension(n_mask_chunks), intent(in), target :: filtered_paralogs_mask
            !! Output mask with all genes disabled that did not pass the filter
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_subset_size)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(work_array_size)
        M_CHECK_NON_NULL(n_mask_chunks)
        M_CHECK_ARRAY_NON_NULL(filtered_paralogs_mask, n_mask_chunks)

        call calc_work_arr_paralog_subsets_size(&
            max_subset_size = max_subset_size,&
            n_genes = n_genes,&
            work_array_size = work_array_size,&
            filtered_paralogs_mask = filtered_paralogs_mask,&
            n_mask_chunks = n_mask_chunks,&
            ierr = ierr&
        )
    end subroutine calc_work_arr_paralog_subsets_size_c

end module tox_paralog_analysis_c
#endif
