#include <src/macros.h>

!> summary: Implementation fixture for an optional mode with no default
!| author: A Developer
!| `fx_scale_by_method` copies its values, or doubles or halves them when asked to, and reports
!| whether it was asked. Its method is optional and has no default, so it is the one mode the C
!| wrapper may be handed as a null pointer -- and must then pass on as absent. Compiled and run by
!| `test_end_to_end_fortran.py`, and through Python and R by `test_end_to_end_impl.py`.
module fx_modes_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private

    public :: METHOD_DOUBLE, METHOD_HALVE
    public :: fx_scale_by_method_impl

    integer(int32), parameter :: METHOD_DOUBLE = 1_int32
        !! double every value
    integer(int32), parameter :: METHOD_HALVE = 2_int32
        !! halve every value

contains

    !> summary: The values as they are, or doubled or halved
    !| author: A Developer
    !| An edge case: an optional mode without a default. Omitted, the values are copied, and
    !| `method_given` says what `present()` saw, so a test can tell "absent" from "present with
    !| a value outside the table".
    pure subroutine fx_scale_by_method_impl(n, values, method, scaled, method_given)
        integer(int32), intent(in) :: n
            !! number of values
        real(real64), dimension(n), intent(in) :: values
            !! the values
        integer(int32), intent(in), optional :: method
            !! how to scale them; omitted, they are copied unchanged
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | double every value | [[fx_modes_impl(module):METHOD_DOUBLE(variable)]] |
            !! | halve every value | [[fx_modes_impl(module):METHOD_HALVE(variable)]] |
        real(real64), dimension(n), intent(out) :: scaled
            !! the scaled values
        integer(int32), intent(out) :: method_given
            !! 1 when a method was given, 0 when it was omitted
        scaled = values
        method_given = merge(1_int32, 0_int32, present(method))
        if (.not. present(method)) return
        select case (method)
            case (METHOD_DOUBLE)
                scaled = 2.0_real64*values
            case (METHOD_HALVE)
                scaled = 0.5_real64*values
        end select
    end subroutine fx_scale_by_method_impl

end module fx_modes_impl
