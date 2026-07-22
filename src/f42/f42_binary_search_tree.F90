#include <src/macros.h>

!> Flat-index-based Binary Search Tree (BST) utilities for 1D range queries.
!| This module provides routines to build a BST index (via sorting), access sorted values,
!| and perform range queries over a real-valued array using the sorted index.
module f42_binary_search_tree
    use safeguard
    use f42_utils, only: sort_array_heapsort, init_perm, binary_search_insertion
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: is_err, set_ok, validate_dimension_size, validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    public :: build_bst_index, get_sorted_value, bst_range_query
contains

    !> M_EXPORT_C
    !| summary: Build the BST index by sorting indices using values in x
    !| AUTHOR_AARON_SCHROEDER
    pure subroutine build_bst_index(values, n_values, sorted_indices, ierr)
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
        call sort_array_heapsort(values, sorted_indices)
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

    !> M_EXPORT_C
    !| summary: Perform a 1D range query over the sorted index
    !| AUTHOR_AARON_SCHROEDER
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
            !! Output array of matching indices.
            !! DM_RESULT_SIZE_IS(n_matches)
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

        ! find index of first value in range
        i_value = binary_search_insertion(values, sorted_indices, lower_bound)
        do while (i_value <= n_values)
            ! check if value is in range and add to output, else exit loop
            if (values(sorted_indices(i_value)) <= upper_bound) then
                n_matches = n_matches + 1
                output_indices(n_matches) = sorted_indices(i_value)
                i_value = i_value + 1
            else
                exit
            end if
        end do
    end subroutine bst_range_query

end module f42_binary_search_tree
