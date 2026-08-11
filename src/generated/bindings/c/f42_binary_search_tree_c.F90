#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_binary_search_tree(module)]]
!| Generated from the implementation; do not edit -- regenerate instead.
module f42_binary_search_tree_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: build_bst_index_c
    public :: bst_range_query_c
    public :: bst_range_query_expert_c

contains

    !> summary: C-wrapper for [[f42_binary_search_tree(module):build_bst_index(subroutine)]]
    subroutine build_bst_index_c(&
            values,&
            n_values,&
            sorted_indices,&
            ierr&
        ) bind(C, name="build_bst_index_c")
        use f42_binary_search_tree, only: build_bst_index

        integer(c_int), intent(in), target :: n_values
            !! Number of elements in values array
        real(c_double), dimension(n_values), intent(in), target :: values
            !! Input real array to be indexed
        integer(c_int), dimension(n_values), intent(out), target :: sorted_indices
            !! Output permutation index
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_values)
        M_CHECK_ARRAY_NON_NULL(values, n_values)
        M_CHECK_ARRAY_NON_NULL(sorted_indices, n_values)

        call build_bst_index(&
            values = values,&
            n_values = n_values,&
            sorted_indices = sorted_indices,&
            ierr = ierr&
        )
    end subroutine build_bst_index_c

    !> summary: C-wrapper for [[f42_binary_search_tree(module):bst_range_query(subroutine)]]
    subroutine bst_range_query_c(&
            values,&
            n_values,&
            lower_bound,&
            upper_bound,&
            output_indices,&
            n_matches,&
            ierr&
        ) bind(C, name="bst_range_query_c")
        use f42_binary_search_tree, only: bst_range_query

        integer(c_int), intent(in), target :: n_values
            !! Number of elements
        real(c_double), dimension(n_values), intent(in), target :: values
            !! Input real array
        real(c_double), intent(in), target :: lower_bound
            !! Lower bound of range (inclusive)
            !! The maximum valid value is `upper_bound`.
        real(c_double), intent(in), target :: upper_bound
            !! Upper bound of range (inclusive)
        integer(c_int), dimension(n_values), intent(out), target :: output_indices
            !! Output array of matching indices.
            !! The first `n_matches` elements will hold the results.
        integer(c_int), intent(out), target :: n_matches
            !! Number of matches found
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_values)
        M_CHECK_NON_NULL(lower_bound)
        M_CHECK_NON_NULL(upper_bound)
        M_CHECK_NON_NULL(n_matches)
        M_CHECK_ARRAY_NON_NULL(values, n_values)
        M_CHECK_ARRAY_NON_NULL(output_indices, n_values)

        call bst_range_query(&
            values = values,&
            n_values = n_values,&
            lower_bound = lower_bound,&
            upper_bound = upper_bound,&
            output_indices = output_indices,&
            n_matches = n_matches,&
            ierr = ierr&
        )
    end subroutine bst_range_query_c

    !> summary: C-wrapper for [[f42_binary_search_tree(module):bst_range_query_expert(subroutine)]]
    subroutine bst_range_query_expert_c(&
            values,&
            values_perm,&
            n_values,&
            lower_bound,&
            upper_bound,&
            output_indices,&
            n_matches,&
            ierr&
        ) bind(C, name="bst_range_query_expert_c")
        use f42_binary_search_tree, only: bst_range_query_expert

        integer(c_int), intent(in), target :: n_values
            !! Number of elements
        real(c_double), dimension(n_values), intent(in), target :: values
            !! Input real array
        integer(c_int), dimension(n_values), intent(in), target :: values_perm
            !! Permutation of `values` in ascending order -- the BST index. The allocating entry
            !! point builds and heapsorts it for you; the expert one takes whatever order you
            !! supply, so a caller that already holds one from
            !! [[f42_binary_search_tree(module):build_bst_index]] can reuse it across queries.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_values`.
        real(c_double), intent(in), target :: lower_bound
            !! Lower bound of range (inclusive)
            !! The maximum valid value is `upper_bound`.
        real(c_double), intent(in), target :: upper_bound
            !! Upper bound of range (inclusive)
        integer(c_int), dimension(n_values), intent(out), target :: output_indices
            !! Output array of matching indices.
            !! The first `n_matches` elements will hold the results.
        integer(c_int), intent(out), target :: n_matches
            !! Number of matches found
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_values)
        M_CHECK_NON_NULL(lower_bound)
        M_CHECK_NON_NULL(upper_bound)
        M_CHECK_NON_NULL(n_matches)
        M_CHECK_ARRAY_NON_NULL(values, n_values)
        M_CHECK_ARRAY_NON_NULL(values_perm, n_values)
        M_CHECK_ARRAY_NON_NULL(output_indices, n_values)

        call bst_range_query_expert(&
            values = values,&
            values_perm = values_perm,&
            n_values = n_values,&
            lower_bound = lower_bound,&
            upper_bound = upper_bound,&
            output_indices = output_indices,&
            n_matches = n_matches,&
            ierr = ierr&
        )
    end subroutine bst_range_query_expert_c

end module f42_binary_search_tree_c
#endif
