#include <src/macros.h>

!> summary: Implementation fixture for the extents the bindings compute and the ones they check
!| author: A Developer
!| An output sized by an integer division that Python must truncate as Fortran does. Compiled and
!| run by `test_end_to_end_fortran.py`, and through C, Python and R by `test_end_to_end_impl.py`.
module fx_extents_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private

    public :: truncated_count_impl

contains

    !> summary: The numbers one to `(n - 10)/3 + 5`
    !| author: A Developer
    !| An edge case: an output extent with an integer division of an operand that may be
    !| negative, which truncates toward zero in Fortran and C and floors in Python's `//`.
    pure subroutine truncated_count_impl(n, counted)
        integer(int32), intent(in) :: n
            !! decides the length
        integer(int32), dimension((n - 10)/3 + 5), intent(out) :: counted
            !! one to its own length
        integer(int32) :: i
        counted = [(i, i = 1, size(counted))]
    end subroutine truncated_count_impl

end module fx_extents_impl
