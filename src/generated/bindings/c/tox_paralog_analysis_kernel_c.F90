#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_paralog_analysis_kernel(module)]]
!| Kernels for detecting paralog-subset expression patterns (dosage effect and subfunctionalization) relative to an ancestral ortholog.
!|
!| Candidate paralog subsets are enumerated as bitmask-encoded gene sets, built up one gene at a time
!| starting from single genes. At every extension step the candidate is scored against the pattern-specific
!| criterion (small angle plus magnitude gain for dosage effect, or bounded residual distance for
!| subfunctionalization); subsets that can no longer satisfy the criterion are pruned instead of being
!| extended further, which keeps the combinatorial subset search tractable.
!|
!| Both pattern-taking kernels are **mode-split**: their `pattern` table names a procedure per value, so
!| the generator emits `detect_dosage_effect`/`detect_subfunctionalization` and
!| `filter_paralogs_by_pattern_dosage_effect`/`_subfunctionalization` into module `tox_paralog_analysis`,
!| replacing the hand-written per-mode wrappers -- and with them the `map_err_arg_pos` remapping, since
!| each generated wrapper validates its own arguments at their own positions. The kernels keep `ierr`
!| only for their runtime paths (bit-mask writes, the subset-search bookkeeping).
module tox_paralog_analysis_kernel_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: mask_check_state_c
    public :: mask_chunk_count_c
    public :: calc_work_arr_paralog_subsets_size_c

contains

    !> summary: C-wrapper for [[tox_paralog_analysis_kernel(module):mask_check_state(function)]]
    subroutine mask_check_state_c(&
            bit_mask,&
            n_bit_mask_elements,&
            i_gene,&
            state,&
            ierr&
        ) bind(C, name="mask_check_state_c")
        use tox_paralog_analysis_kernel, only: mask_check_state

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

    !> summary: C-wrapper for [[tox_paralog_analysis_kernel(module):mask_chunk_count(subroutine)]]
    subroutine mask_chunk_count_c(&
            n_genes,&
            count,&
            ierr&
        ) bind(C, name="mask_chunk_count_c")
        use tox_paralog_analysis_kernel, only: mask_chunk_count

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

    !> summary: C-wrapper for [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size(subroutine)]]
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
        use tox_paralog_analysis_kernel, only: calc_work_arr_paralog_subsets_size

        integer(c_int), intent(in), target :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        integer(c_int), intent(inout), target :: max_subset_size
            !! maximum size that a subset must not exceed. Zero is in range and means there is
            !! nothing to size a work array for, which is reported back as a size of zero.
            !! The minimum valid value is `0_int32`.
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

end module tox_paralog_analysis_kernel_c
#endif
