#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_paralog_analysis_impl(module)]]
!| Paralog-subset expression patterns -- dosage effect and subfunctionalization -- relative to an ancestral ortholog.
!|
!| Candidate paralog subsets are enumerated as bitmask-encoded gene sets, built up one gene at a time
!| starting from single genes. At every extension step the candidate is scored against the pattern-specific
!| criterion (small angle plus magnitude gain for dosage effect, or bounded residual distance for
!| subfunctionalization); subsets that can no longer satisfy the criterion are pruned instead of being
!| extended further, which keeps the combinatorial subset search tractable.
!|
!| There is a routine per pattern rather than one taking a pattern flag:
!| `detect_dosage_effect` and `detect_subfunctionalization` find the subsets, and
!| `filter_paralogs_by_pattern_dosage_effect` / `_subfunctionalization` reduce a set of genes to
!| those matching. Which criterion a result was scored against is therefore visible at the call
!| site, and each routine takes only the thresholds its own pattern needs.
module tox_paralog_analysis_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: mask_check_state_c
    public :: mask_chunk_count_c

contains

    !> summary: C-wrapper for [[tox_paralog_analysis_impl(module):mask_check_state(function)]]
    subroutine mask_check_state_c(&
            bit_mask,&
            n_bit_mask_elements,&
            i_gene,&
            state,&
            ierr&
        ) bind(C, name="mask_check_state_c")
        use tox_paralog_analysis_impl, only: mask_check_state

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

    !> summary: C-wrapper for [[tox_paralog_analysis_impl(module):mask_chunk_count(subroutine)]]
    subroutine mask_chunk_count_c(&
            n_genes,&
            count,&
            ierr&
        ) bind(C, name="mask_chunk_count_c")
        use tox_paralog_analysis_impl, only: mask_chunk_count

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

end module tox_paralog_analysis_impl_c
#endif
