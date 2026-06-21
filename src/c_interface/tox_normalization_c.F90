#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[tox_normalization(module)]]
!| Module with normalization routines for tensor omics.
module tox_normalization_c
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

    !> summary: C-wrapper for [[tox_normalization(module):normalize_unit_length(subroutine)]]
    !| Normalizes an input vector to unit length in-place
    subroutine normalize_unit_length_c(&
            vector,&
            n_dims,&
            ierr&
            ) bind(C, name="normalize_unit_length_c")
        use tox_normalization, only: normalize_unit_length
        use tox_normalization
        integer(c_int), intent(in), target :: n_dims
            !! number of elements in `vector`
        real(c_double), intent(inout), dimension(n_dims), target :: vector
            !! Vector that will be normalized to unit length
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(vector)
        M_CHECK_NON_NULL(n_dims)
        call  normalize_unit_length(&
            vector = vector,&
            n_dims = n_dims,&
            ierr = ierr&
        )
    end subroutine normalize_unit_length_c

end module tox_normalization_c
#endif