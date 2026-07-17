#include <src/macros.h>

!> summary: Fixture module covering the awkward cases
!| Characters, shapes, masks, alloc/expert pairs and the directives that need resolving.
module fx_edges
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    implicit none

    integer(int32), parameter :: MODE_GROUP_ORTHOLOGS = 1
        !! group the values by their ortholog set
    integer(int32), parameter :: MODE_UNGROUPED = 2
        !! treat every value on its own

contains

    !> category: C-interface
    !| summary: Every character length form that is allowed
    !| author: A Developer
    subroutine fx_strings(assumed, fixed, sized, n_chars, names, n_names, ierr)
        character(len=*), intent(in) :: assumed
            !! assumed length, so C must be told how long it is
        character(len=8), intent(in) :: fixed
            !! a constant length, needing no extra argument
        integer(int32), intent(in) :: n_chars
            !! length of `sized`
        character(len=n_chars), intent(in) :: sized
            !! a length named by another argument
        integer(int32), intent(in) :: n_names
            !! elements of `names`
        character(len=16), dimension(n_names), intent(in) :: names
            !! a vector of strings, the deepest character rank supported
        integer(int32), intent(out) :: ierr
            !! Error code
    end subroutine fx_strings

    !> category: C-interface
    !| summary: A flat array whose shape travels separately
    !| author: A Developer
    subroutine fx_serialized(data, data_shape, ierr)
        real(real64), dimension(:), intent(in) :: data
            !! the values, flat, so any rank can be passed
        integer(int32), dimension(:), intent(in) :: data_shape
            !! the extents of `data`, one per dimension
        integer(int32), intent(out) :: ierr
            !! Error code
    end subroutine fx_serialized

    !> category: C-interface
    !| summary: A mask and the count derived from it
    !| author: A Developer
    subroutine fx_masked(genes_selection_mask, n_selected_genes, n_genes, results, n_results, ierr)
        integer(int32), intent(in) :: n_genes
            !! elements of the mask
        logical, dimension(n_genes), intent(in) :: genes_selection_mask
            !! which genes to use
        integer(int32), intent(in) :: n_selected_genes
            !! how many genes the mask selects, which the interfacing languages can count
        real(real64), dimension(n_genes), intent(out) :: results
            !! the results.
            !! DM_RESULT_SIZE_IS(n_results)
        integer(int32), intent(out) :: n_results
            !! how many leading elements of `results` were filled
        integer(int32), intent(out) :: ierr
            !! Error code
    end subroutine fx_masked

    !> category: C-interface
    !| summary: An already interoperable logical, which needs no conversion
    !| author: A Developer
    subroutine fx_c_bool_flag(flag, ierr)
        logical(c_bool), intent(in) :: flag
            !! a flag declared with the interoperable kind
        integer(int32), intent(out) :: ierr
            !! Error code
    end subroutine fx_c_bool_flag

    !> category: C-interface
    !| summary: The allocating half of a pair
    !| author: A Developer
    subroutine fx_cluster_alloc(values, n_values, n_clusters, ierr)
        integer(int32), intent(in) :: n_values
            !! elements of `values`
        real(real64), dimension(n_values), intent(in) :: values
            !! the values
        integer(int32), intent(out) :: n_clusters
            !! how many clusters were found
        integer(int32), intent(out) :: ierr
            !! Error code
    end subroutine fx_cluster_alloc

    !> category: C-interface
    !| summary: The expert half of a pair, which owns no allocation
    !| author: A Developer
    subroutine fx_cluster(values, n_values, tmp_work, n_work, n_clusters, ierr)
        integer(int32), intent(in) :: n_values
            !! elements of `values`
        real(real64), dimension(n_values), intent(in) :: values
            !! the values
        integer(int32), intent(in) :: n_work
            !! size of the work array.
            !! DM_OUTPUT_FROM(n_work, fx_work_size, fx_edges, AUTO)
        real(real64), dimension(n_work), intent(out) :: tmp_work
            !! scratch space
        integer(int32), intent(out) :: n_clusters
            !! how many clusters were found
        integer(int32), intent(out) :: ierr
            !! Error code
    end subroutine fx_cluster

    !> category: C-interface
    !| summary: Works out the work array size for fx_cluster
    !| author: A Developer
    subroutine fx_work_size(n_values, n_work, ierr)
        integer(int32), intent(in) :: n_values
            !! elements the caller intends to pass
        integer(int32), intent(out) :: n_work
            !! required size of the work array
        integer(int32), intent(out) :: ierr
            !! Error code
    end subroutine fx_work_size

    !> category: C-interface
    !| summary: An optional with no default, which C may pass a null pointer for
    !| author: A Developer
    subroutine fx_nullable(values, n_values, mode, ortholog_set, n_orthologs, ierr)
        integer(int32), intent(in) :: n_values
            !! elements of `values`
        real(real64), dimension(n_values), intent(in) :: values
            !! the values
        integer(int32), intent(in) :: mode
            !! how to group the values
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | group by ortholog set | [[fx_edges(module):MODE_GROUP_ORTHOLOGS(variable)]] |
            !! | do not group at all | [[fx_edges(module):MODE_UNGROUPED(variable)]] |
        integer(int32), intent(in) :: n_orthologs
            !! elements of `ortholog_set`
        integer(int32), dimension(n_orthologs), intent(in), optional :: ortholog_set
            !! which ortholog set each value belongs to.
            !! DM_REQUIRED_IF_MODE(mode, fx_edges, MODE_GROUP_ORTHOLOGS)
        integer(int32), intent(out) :: ierr
            !! Error code
    end subroutine fx_nullable

    !> category: C-interface
    !| summary: A procedure that reports no errors of its own
    !| author: A Developer
    pure subroutine fx_no_ierr(value, doubled)
        real(real64), intent(in) :: value
            !! the input
        real(real64), intent(out) :: doubled
            !! twice the input
    end subroutine fx_no_ierr

end module fx_edges
