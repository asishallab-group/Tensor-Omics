#include <src/macros.h>

!! Flat-index-based Binary Search Tree (BST) utilities for 1D range queries.
!!
!! This module provides routines to build a BST index (via sorting), access sorted values,
!! and perform efficient range queries over a real-valued array.
module f42_binary_search_tree
    use safeguard
    use f42_utils, only: sort_array_heapsort, init_perm
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: is_err, set_ok, validate_dimension_size, validate_in_range_int, validate_in_range_real
    implicit none
    public :: build_bst_index, get_sorted_value, bst_range_query
contains

    !> AUTHOR_AARON_SCHROEDER
    !| Build the BST index by sorting indices using values in x.
    subroutine build_bst_index(values, n_values, sorted_indices, ierr)
        integer(int32), intent(in) :: n_values
            !! Number of elements in values array
        real(real64), intent(in) :: values(n_values)
            !! Input real array to be indexed
        integer(int32), intent(out) :: sorted_indices(n_values)
            !! Output permutation index
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_values, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return

        call init_perm(sorted_indices)
        print "(*(G0, :, ', '))", sorted_indices
        print *
        call sort_array_heapsort(values, sorted_indices)
        print "(*(G0, :, ', '))", sorted_indices
        print *
        call sort_array_heapsort(values, sorted_indices, n_values)
        print "(*(G0, :, ', '))", sorted_indices
        print *
    end subroutine build_bst_index

    !> AUTHOR_AARON_SCHROEDER
    !| Get the value at the sorted position.
    function get_sorted_value(values, sorted_indices, position, ierr) result(sorted_value)
        real(real64), intent(in) :: values(:)
            !! Input real array
        integer(int32), intent(in) :: sorted_indices(:)
            !! Permutation index array
        integer(int32), intent(in) :: position
            !! Sorted position (1-based)
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64) :: sorted_value

        call set_ok(ierr)
        call validate_in_range_int(position, ierr, min=1_int32, max=size(sorted_indices, kind=int32), arg_pos=3_int32)
        if (is_err(ierr)) return
        call validate_in_range_int(sorted_indices(position), ierr, min=1_int32, max=size(values, kind=int32), arg_pos=2_int32)
        if (is_err(ierr)) return

        sorted_value = values(sorted_indices(position))
    end function get_sorted_value

    !> AUTHOR_AARON_SCHROEDER
    !| Perform a 1D range query over the sorted index.
    pure subroutine bst_range_query(values, sorted_indices, n_values, lower_bound, upper_bound, &
                                    output_indices, n_matches, ierr)

        integer(int32), intent(in) :: n_values
            !! Number of elements
        real(real64), intent(in) :: values(n_values)
            !! Input real array
        integer(int32), intent(in) :: sorted_indices(n_values)
            !! Permutation index array (sorted)
        real(real64), intent(in) :: lower_bound
            !! Lower bound of range (inclusive)
        real(real64), intent(in) :: upper_bound
            !! Upper bound of range (inclusive)
        integer(int32), intent(out) :: output_indices(n_values)
            !! Output array of matching indices
        integer(int32), intent(out) :: n_matches
            !! Number of matches found
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32) :: i_value

        call set_ok(ierr)

        call validate_dimension_size(n_values, ierr, arg_pos=3_int32)
        call validate_in_range_real(lower_bound, ierr, max=upper_bound, arg_pos=4_int32)
        if (is_err(ierr)) return

        n_matches = 0
        do i_value = 1, n_values
            if (values(sorted_indices(i_value)) >= lower_bound .and. &
                values(sorted_indices(i_value)) <= upper_bound) then
                n_matches = n_matches + 1
                output_indices(n_matches) = sorted_indices(i_value)
            else if (values(sorted_indices(i_value)) > upper_bound) then
                exit
            end if
        end do
    end subroutine bst_range_query

end module f42_binary_search_tree

!> Wrapper using C for getting range query usable by python
pure subroutine bst_range_query_c(values, sorted_indices, n_values, lower_bound, upper_bound, &
                                  output_indices, n_matches, ierr) bind(C, name='bst_range_query_c')
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use f42_binary_search_tree, only: bst_range_query
    M_USE_NULL_VALIDATION
    implicit none
    real(c_double), intent(in), target :: values(n_values)
        !! Input real array (C-style)
    integer(c_int), intent(in), target :: sorted_indices(n_values)
        !! Permutation index array (C-style)
    integer(c_int), intent(in), target :: n_values
        !! Number of elements
    real(c_double), intent(in), target :: lower_bound
        !! Lower bound of range
    real(c_double), intent(in), target :: upper_bound
        !! Upper bound of range
    integer(c_int), intent(out), target :: output_indices(n_values)
        !! Output array (C-style)
    integer(c_int), intent(out), target :: n_matches
        !! Number of matches found
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_values)
    M_CHECK_NON_NULL(lower_bound)
    M_CHECK_NON_NULL(upper_bound)
    M_CHECK_NON_NULL(values)
    M_CHECK_NON_NULL(sorted_indices)
    M_CHECK_NON_NULL(output_indices)
    M_CHECK_NON_NULL(n_matches)

    call bst_range_query(values, sorted_indices, n_values, lower_bound, upper_bound, &
                         output_indices, n_matches, ierr)
end subroutine bst_range_query_c

!> Wrapper using C for building BST index usable by python
subroutine build_bst_index_c(values, n_values, sorted_indices, ierr) &
    bind(C, name='build_bst_index_c')
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use f42_binary_search_tree
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: n_values
        !! Number of elements
    real(c_double), intent(in), target :: values(n_values)
        !! Input real array (C-style)
    integer(c_int), intent(out), target :: sorted_indices(n_values)
        !! Output permutation index (C-style)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_values)
    M_CHECK_NON_NULL(values)
    M_CHECK_NON_NULL(sorted_indices)

    call build_bst_index(values, n_values, sorted_indices, ierr)
end subroutine build_bst_index_c
