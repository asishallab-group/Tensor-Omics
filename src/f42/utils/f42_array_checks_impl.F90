#include <src/macros.h>

!> Questions about the elements of an array -- are they all equal to a value, is any of them below
!| or above a bound -- answered without the temporary array an expression such as
!| `all(array == value)` may build. Each stops at the first element that decides and takes any
!| section of an array, a strided one included, without a copy. Rank 1 and rank 2; a matrix is
!| checked column by column.
!|
!| The arrays and values hold no NaN, as everywhere in Tensor Omics: a caller that may hold one
!| handles it before the call.
!|
!| One of the modules [[f42_utils_impl(module)]] gathers; `use f42_utils_impl` reaches all of them.
module f42_array_checks_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    M_IMPLICIT_NONE

    !> `.true.` when every element of `array` equals `value`, as `all(array == value)` without its
    !| temporary; `.not. all_equal(array, value)` asks whether any element differs.
    interface all_equal
        module procedure all_equal_int, all_equal_real, all_equal_int_2d, all_equal_real_2d
    end interface all_equal

    !> `.true.` when some element of `array` is strictly less than `bound`, as `any(array < bound)`
    !| without its temporary.
    interface any_below
        module procedure any_below_int, any_below_real, any_below_int_2d, any_below_real_2d
    end interface any_below

    !> `.true.` when some element of `array` is strictly greater than `bound`, as
    !| `any(array > bound)` without its temporary.
    interface any_above
        module procedure any_above_int, any_above_real, any_above_int_2d, any_above_real_2d
    end interface any_above

contains

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether every element of `array` equals `value`; `.true.` for an empty array, as `all` gives.
    pure logical function all_equal_int(array, value) result(all_are_equal)
        integer(int32), dimension(:), intent(in) :: array
            !! The elements, any section of an array
        integer(int32), intent(in) :: value
            !! The value every element is compared with

        integer(int32) :: i_element

        all_are_equal = .true.
        do i_element = 1, size(array, kind=int32)
            if (array(i_element) /= value) then
                all_are_equal = .false.
                return
            end if
        end do
    end function all_equal_int

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether every element of `array` equals `value`, compared exactly with `==`, so `-0.0` equals
    !| `0.0`; `.true.` for an empty array, as `all` gives. No NaN in `array` or `value`.
    pure logical function all_equal_real(array, value) result(all_are_equal)
        real(real64), dimension(:), intent(in) :: array
            !! The elements, any section of an array
        real(real64), intent(in) :: value
            !! The value every element is compared with

        integer(int32) :: i_element

        all_are_equal = .true.
        do i_element = 1, size(array, kind=int32)
            if (array(i_element) /= value) then
                all_are_equal = .false.
                return
            end if
        end do
    end function all_equal_real

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether some element of `array` is strictly less than `bound`; `.false.` for an empty array.
    pure logical function any_below_int(array, bound) result(found)
        integer(int32), dimension(:), intent(in) :: array
            !! The elements, any section of an array
        integer(int32), intent(in) :: bound
            !! An element equal to it is not below it

        integer(int32) :: i_element

        found = .false.
        do i_element = 1, size(array, kind=int32)
            if (array(i_element) < bound) then
                found = .true.
                return
            end if
        end do
    end function any_below_int

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether some element of `array` is strictly less than `bound`; `.false.` for an empty array.
    !| `-0.0` is not below `0.0`, and `-Inf` is below every finite bound. No NaN in `array` or `bound`.
    pure logical function any_below_real(array, bound) result(found)
        real(real64), dimension(:), intent(in) :: array
            !! The elements, any section of an array
        real(real64), intent(in) :: bound
            !! An element equal to it is not below it

        integer(int32) :: i_element

        found = .false.
        do i_element = 1, size(array, kind=int32)
            if (array(i_element) < bound) then
                found = .true.
                return
            end if
        end do
    end function any_below_real

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether some element of `array` is strictly greater than `bound`; `.false.` for an empty array.
    pure logical function any_above_int(array, bound) result(found)
        integer(int32), dimension(:), intent(in) :: array
            !! The elements, any section of an array
        integer(int32), intent(in) :: bound
            !! An element equal to it is not above it

        integer(int32) :: i_element

        found = .false.
        do i_element = 1, size(array, kind=int32)
            if (array(i_element) > bound) then
                found = .true.
                return
            end if
        end do
    end function any_above_int

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether some element of `array` is strictly greater than `bound`; `.false.` for an empty array.
    !| `0.0` is not above `-0.0`, and `Inf` is above every finite bound. No NaN in `array` or `bound`.
    pure logical function any_above_real(array, bound) result(found)
        real(real64), dimension(:), intent(in) :: array
            !! The elements, any section of an array
        real(real64), intent(in) :: bound
            !! An element equal to it is not above it

        integer(int32) :: i_element

        found = .false.
        do i_element = 1, size(array, kind=int32)
            if (array(i_element) > bound) then
                found = .true.
                return
            end if
        end do
    end function any_above_real

    !> AUTHOR_FRANZ_ERIC_SILL
    !| [[f42_array_checks_impl(module):all_equal_int(function)]] over a matrix, column by column, so each
    !| extent, not the element count, has to fit `int32`.
    pure logical function all_equal_int_2d(array, value) result(all_are_equal)
        integer(int32), dimension(:, :), intent(in) :: array
            !! The elements, any section of a matrix
        integer(int32), intent(in) :: value
            !! The value every element is compared with

        integer(int32) :: i_column

        all_are_equal = .true.
        do i_column = 1, size(array, 2, kind=int32)
            if (.not. all_equal_int(array(:, i_column), value)) then
                all_are_equal = .false.
                return
            end if
        end do
    end function all_equal_int_2d

    !> AUTHOR_FRANZ_ERIC_SILL
    !| [[f42_array_checks_impl(module):all_equal_real(function)]] over a matrix, column by column, so each
    !| extent, not the element count, has to fit `int32`.
    pure logical function all_equal_real_2d(array, value) result(all_are_equal)
        real(real64), dimension(:, :), intent(in) :: array
            !! The elements, any section of a matrix
        real(real64), intent(in) :: value
            !! The value every element is compared with

        integer(int32) :: i_column

        all_are_equal = .true.
        do i_column = 1, size(array, 2, kind=int32)
            if (.not. all_equal_real(array(:, i_column), value)) then
                all_are_equal = .false.
                return
            end if
        end do
    end function all_equal_real_2d

    !> AUTHOR_FRANZ_ERIC_SILL
    !| [[f42_array_checks_impl(module):any_below_int(function)]] over a matrix, column by column, so each
    !| extent, not the element count, has to fit `int32`.
    pure logical function any_below_int_2d(array, bound) result(found)
        integer(int32), dimension(:, :), intent(in) :: array
            !! The elements, any section of a matrix
        integer(int32), intent(in) :: bound
            !! An element equal to it is not below it

        integer(int32) :: i_column

        found = .false.
        do i_column = 1, size(array, 2, kind=int32)
            if (any_below_int(array(:, i_column), bound)) then
                found = .true.
                return
            end if
        end do
    end function any_below_int_2d

    !> AUTHOR_FRANZ_ERIC_SILL
    !| [[f42_array_checks_impl(module):any_below_real(function)]] over a matrix, column by column, so each
    !| extent, not the element count, has to fit `int32`.
    pure logical function any_below_real_2d(array, bound) result(found)
        real(real64), dimension(:, :), intent(in) :: array
            !! The elements, any section of a matrix
        real(real64), intent(in) :: bound
            !! An element equal to it is not below it

        integer(int32) :: i_column

        found = .false.
        do i_column = 1, size(array, 2, kind=int32)
            if (any_below_real(array(:, i_column), bound)) then
                found = .true.
                return
            end if
        end do
    end function any_below_real_2d

    !> AUTHOR_FRANZ_ERIC_SILL
    !| [[f42_array_checks_impl(module):any_above_int(function)]] over a matrix, column by column, so each
    !| extent, not the element count, has to fit `int32`.
    pure logical function any_above_int_2d(array, bound) result(found)
        integer(int32), dimension(:, :), intent(in) :: array
            !! The elements, any section of a matrix
        integer(int32), intent(in) :: bound
            !! An element equal to it is not above it

        integer(int32) :: i_column

        found = .false.
        do i_column = 1, size(array, 2, kind=int32)
            if (any_above_int(array(:, i_column), bound)) then
                found = .true.
                return
            end if
        end do
    end function any_above_int_2d

    !> AUTHOR_FRANZ_ERIC_SILL
    !| [[f42_array_checks_impl(module):any_above_real(function)]] over a matrix, column by column, so each
    !| extent, not the element count, has to fit `int32`.
    pure logical function any_above_real_2d(array, bound) result(found)
        real(real64), dimension(:, :), intent(in) :: array
            !! The elements, any section of a matrix
        real(real64), intent(in) :: bound
            !! An element equal to it is not above it

        integer(int32) :: i_column

        found = .false.
        do i_column = 1, size(array, 2, kind=int32)
            if (any_above_real(array(:, i_column), bound)) then
                found = .true.
                return
            end if
        end do
    end function any_above_real_2d

end module f42_array_checks_impl
