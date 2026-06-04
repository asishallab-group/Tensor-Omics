#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> Module for C-wrappers for [[tox_paralog_analysis(module)]]
module tox_paralog_analysis_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err
contains

    !> C-wrapper for [[tox_paralog_analysis(module):mask_get_first_successor_idx(function)]]
    !| Helper function that returns the index after the last active gene in `bit_mask`, so the first succeeding gene.
    subroutine mask_get_first_successor_idx_c(n_bit_mask_elements, n_bit_mask_elements, bit_mask, idx) bind(C, name="mask_get_first_successor_idx_c")
        use tox_paralog_analysis, only: mask_get_first_successor_idx
        integer(c_int), intent(in), target :: n_bit_mask_elements
            !!  Size of the 1. dimension/extent of `bit_mask`
        integer(c_int), intent(in), target :: n_bit_mask_elements
            !!  Size of the 1. dimension/extent of `bit_mask`
        integer(c_int), dimension(n_bit_mask_elements), intent(in), target :: bit_mask
            !! chunked mask to mark active genes
        integer(c_int), intent(out), target :: idx
            !! index of last active gene
            !!
    
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(n_bit_mask_elements)
        M_CHECK_NON_NULL(n_bit_mask_elements)
        M_CHECK_NON_NULL(bit_mask)
        M_CHECK_NON_NULL(idx)
    
        idx = mask_get_first_successor_idx(bit_mask)
    end subroutine mask_get_first_successor_idx_c

end module tox_paralog_analysis_c
#endif