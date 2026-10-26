#include <src/macros.h>

!> summary: Implementation fixture for the extent forms and purities a wrapper must reproduce
!| author: A Developer
!| An assumed-shape input the wrapper validates, extents and character lengths named by a module
!| constant, and an `impure elemental` implementation. Compiled and run by
!| `test_end_to_end_fortran.py`, and through Python and R by `test_end_to_end_impl.py`.
module fx_shapes_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private

    public :: FX_N_AXES, FX_LABEL_LEN, n_noted
    public :: sum_values_impl, axis_sums_impl, label_length_impl, note_value_impl

    integer(int32), parameter :: FX_N_AXES = 2
        !! coordinates per point
    integer(int32), parameter :: FX_LABEL_LEN = 8
        !! characters in a label
    integer(int32), save :: n_noted = 0_int32
        !! how often `note_value_impl` ran

contains

    !> summary: The sum of an assumed-shape array
    !| author: A Developer
    !| An edge case: the array has no extent to name, and is finite by default, so the wrapper
    !| validates it with the size it asks of the array.
    pure subroutine sum_values_impl(values, total)
        real(real64), intent(in) :: values(:)
            !! the values
        real(real64), intent(out) :: total
            !! their sum
        total = sum(values)
    end subroutine sum_values_impl

    !> summary: The sum of the points along each axis
    !| author: A Developer
    !| An edge case: an extent named by a module constant, which the wrappers must import.
    pure subroutine axis_sums_impl(n_points, points, sums)
        integer(int32), intent(in) :: n_points
            !! number of points
        real(real64), dimension(FX_N_AXES, n_points), intent(in) :: points
            !! one point per column
        real(real64), dimension(FX_N_AXES), intent(out) :: sums
            !! the sum of every axis
        sums = sum(points, dim=2)
    end subroutine axis_sums_impl

    !> summary: The length of a label without its trailing blanks, and a fixed tag
    !| author: A Developer
    !| An edge case: a character length named by a module constant, on an input and an output.
    pure subroutine label_length_impl(label, n_chars, tag)
        character(len=FX_LABEL_LEN), intent(in) :: label
            !! the label
        integer(int32), intent(out) :: n_chars
            !! its length without trailing blanks
        character(len=FX_LABEL_LEN), intent(out) :: tag
            !! always `fx`
        n_chars = len_trim(label)
        tag = 'fx'
    end subroutine label_length_impl

    !> summary: Copies a value and counts that it ran
    !| author: A Developer
    !| An edge case: `impure elemental`, so the wrapper calling it must be impure too.
    impure elemental subroutine note_value_impl(value, copy)
        real(real64), intent(in) :: value
            !! the value
        real(real64), intent(out) :: copy
            !! the same value
        n_noted = n_noted + 1_int32
        copy = value
    end subroutine note_value_impl

end module fx_shapes_impl
