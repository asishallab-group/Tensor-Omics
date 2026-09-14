#include <src/macros.h>

!> Flat-index-based Binary Search Tree (BST) utilities for 1D range queries.
!| This module provides routines to build a BST index (via sorting), access sorted values,
!| and perform range queries over a real-valued array using the sorted index.
!|
!| Generated from [[f42_binary_search_tree_impl(module)]]; do not edit -- regenerate instead.
module f42_binary_search_tree
    use f42_binary_search_tree_impl, only: bst_range_query_impl, build_bst_index_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_sort_impl, only: init_perm, sort_array_heapsort
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: build_bst_index
    public :: bst_range_query
    public :: bst_range_query_expert

contains

    !> summary: Validates its inputs, then calls [[f42_binary_search_tree_impl(module):build_bst_index_impl]].
    pure subroutine build_bst_index(&
            values,&
            n_values,&
            sorted_indices,&
            ierr&
        )
        integer(int32), intent(in) :: n_values
            !! Number of elements in values array
        real(real64), dimension(n_values), intent(in) :: values
            !! Input real array to be indexed
        integer(int32), dimension(n_values), intent(out) :: sorted_indices
            !! Output permutation index
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_values, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(values, n_values, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call build_bst_index_impl(&
            values = values,&
            n_values = n_values,&
            sorted_indices = sorted_indices&
        )
    end subroutine build_bst_index

    !> summary: Validates its inputs, prepares what [[f42_binary_search_tree_impl(module):bst_range_query_impl]] needs, then calls it. The entry point to reach for first; see [[f42_binary_search_tree(module):bst_range_query_expert]] to prepare it yourself.
    pure subroutine bst_range_query(&
            values,&
            n_values,&
            lower_bound,&
            upper_bound,&
            output_indices,&
            n_matches,&
            ierr&
        )
        integer(int32), intent(in) :: n_values
            !! Number of elements
        real(real64), dimension(n_values), intent(in) :: values
            !! Input real array
        real(real64), intent(in) :: lower_bound
            !! Lower bound of range (inclusive)
            !! The maximum valid value is `upper_bound`.
        real(real64), intent(in) :: upper_bound
            !! Upper bound of range (inclusive)
        integer(int32), dimension(n_values), intent(out) :: output_indices
            !! Output array of matching indices.
            !! The first `n_matches` elements will hold the results.
        integer(int32), intent(out) :: n_matches
            !! Number of matches found
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: values_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_values, ierr, arg_pos=2_int32)
        call validate_in_range_real(lower_bound, ierr, arg_pos=3_int32, max=upper_bound)
        call validate_in_range_real(upper_bound, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(values, n_values, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(values_perm(n_values))
        call init_perm(values_perm)
        call sort_array_heapsort(values, values_perm)

        call bst_range_query_impl(&
            values = values,&
            values_perm = values_perm,&
            n_values = n_values,&
            lower_bound = lower_bound,&
            upper_bound = upper_bound,&
            output_indices = output_indices,&
            n_matches = n_matches&
        )
    end subroutine bst_range_query

    !> summary: Validates its inputs, then calls [[f42_binary_search_tree_impl(module):bst_range_query_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[f42_binary_search_tree(module):bst_range_query]] does both.
    pure subroutine bst_range_query_expert(&
            values,&
            values_perm,&
            n_values,&
            lower_bound,&
            upper_bound,&
            output_indices,&
            n_matches,&
            ierr&
        )
        integer(int32), intent(in) :: n_values
            !! Number of elements
        real(real64), dimension(n_values), intent(in) :: values
            !! Input real array
        integer(int32), dimension(n_values), intent(in) :: values_perm
            !! Permutation of `values` in ascending order -- the BST index. The allocating entry
            !! point builds and heapsorts it for you; the expert one takes whatever order you
            !! supply, so a caller that already holds one from
            !! [[f42_binary_search_tree(module):build_bst_index]] can reuse it across queries.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_values`.
        real(real64), intent(in) :: lower_bound
            !! Lower bound of range (inclusive)
            !! The maximum valid value is `upper_bound`.
        real(real64), intent(in) :: upper_bound
            !! Upper bound of range (inclusive)
        integer(int32), dimension(n_values), intent(out) :: output_indices
            !! Output array of matching indices.
            !! The first `n_matches` elements will hold the results.
        integer(int32), intent(out) :: n_matches
            !! Number of matches found
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_values, ierr, arg_pos=3_int32)
        call validate_in_range_real(lower_bound, ierr, arg_pos=4_int32, max=upper_bound)
        call validate_in_range_real(upper_bound, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(values, n_values, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(values_perm, n_values, ierr, arg_pos=2_int32, min=1_int32, max=n_values)
        if (is_err(ierr)) return
#endif

        call bst_range_query_impl(&
            values = values,&
            values_perm = values_perm,&
            n_values = n_values,&
            lower_bound = lower_bound,&
            upper_bound = upper_bound,&
            output_indices = output_indices,&
            n_matches = n_matches&
        )
    end subroutine bst_range_query_expert

end module f42_binary_search_tree
