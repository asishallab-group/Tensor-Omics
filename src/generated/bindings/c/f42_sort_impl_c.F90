#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_sort_impl(module)]]
!| Indirect sorting -- every routine reorders a permutation vector rather than the data -- plus the searches over a sorted array.
!|
!| One of the modules [[f42_utils_impl(module)]] gathers; `use f42_utils_impl` reaches all of them.
module f42_sort_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: sort_real_get_perm_c

contains

    !> summary: C-wrapper for [[f42_sort_impl(module):sort_real_get_perm(subroutine)]]
    subroutine sort_real_get_perm_c(&
            array,&
            n,&
            perm,&
            ierr&
        ) bind(C, name="sort_real_get_perm_c")
        use f42_sort_impl, only: sort_real_get_perm

        integer(c_int), intent(in), target :: n
            !! elements of `array`
        real(c_double), dimension(n), intent(in), target :: array
            !! values to sort
        integer(c_int), dimension(n), intent(out), target :: perm
            !! permutation that sorts `array` ascending, NaN last
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n)
        M_CHECK_ARRAY_NON_NULL(array, n)
        M_CHECK_ARRAY_NON_NULL(perm, n)

        call sort_real_get_perm(&
            array = array,&
            n = n,&
            perm = perm,&
            ierr = ierr&
        )
    end subroutine sort_real_get_perm_c

end module f42_sort_impl_c
#endif
