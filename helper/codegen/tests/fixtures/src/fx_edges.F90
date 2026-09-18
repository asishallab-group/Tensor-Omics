#include <src/macros.h>

!> summary: Fixture module covering the awkward cases
!| Characters, shapes, masks, a two-tier pair and the directives that need resolving.
module fx_edges
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err
    implicit none

    integer(int32), parameter :: MODE_GROUP_ORTHOLOGS = 1
        !! group the values by their ortholog set
    integer(int32), parameter :: MODE_UNGROUPED = 2
        !! treat every value on its own

contains

    !> M_EXPORT_C
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

    !> M_EXPORT_C
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

    !> M_EXPORT_C
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

    !> M_EXPORT_C
    !| summary: A string coming back out, scalar and vector
    !| author: A Developer
    subroutine fx_labels(values, n_values, label, labels, ierr)
        integer(int32), intent(in) :: n_values
            !! elements of `values`
        real(real64), dimension(n_values), intent(in) :: values
            !! the values
        character(len=16), intent(out) :: label
            !! a single string coming back
        character(len=8), dimension(n_values), intent(out) :: labels
            !! one string per value, the deepest character rank supported
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i

        call set_ok(ierr)
        label = "summary"
        do i = 1, n_values
            if (values(i) > 0.0_real64) then
                labels(i) = "pos"
            else
                labels(i) = "nonpos"
            end if
        end do
    end subroutine fx_labels

    !> M_EXPORT_C
    !| summary: A string going in, scalar and vector
    !| author: A Developer
    subroutine fx_count_matching(names, n_names, wanted, n_matching, ierr)
        integer(int32), intent(in) :: n_names
            !! elements of `names`
        character(len=16), dimension(n_names), intent(in) :: names
            !! the names to search
        character(len=*), intent(in) :: wanted
            !! the name to look for
        integer(int32), intent(out) :: n_matching
            !! how many of `names` equal `wanted`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i

        call set_ok(ierr)
        n_matching = 0
        do i = 1, n_names
            if (trim(names(i)) == trim(wanted)) n_matching = n_matching + 1
        end do
    end subroutine fx_count_matching

    !> M_EXPORT_C
    !| summary: Strings the caller may omit, which C passes as a null pointer
    !| author: A Developer
    subroutine fx_optional_strings(tag, extras, n_given, ierr)
        character(len=*), intent(in), optional :: tag
            !! a scalar string the caller may leave out
        character(len=*), dimension(:), intent(in), optional :: extras
            !! a vector of strings the caller may leave out
        integer(int32), intent(out) :: n_given
            !! `tag`, if given and not blank, plus one per element of `extras`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        n_given = 0
        if (present(tag)) then
            if (len_trim(tag) > 0) n_given = n_given + 1
        end if
        if (present(extras)) n_given = n_given + int(size(extras), int32)
    end subroutine fx_optional_strings

    !> M_EXPORT_C
    !| summary: An already interoperable logical, which needs no conversion
    !| author: A Developer
    subroutine fx_c_bool_flag(flag, ierr)
        logical(c_bool), intent(in) :: flag
            !! a flag declared with the interoperable kind
        integer(int32), intent(out) :: ierr
            !! Error code
    end subroutine fx_c_bool_flag

    !> M_EXPORT_C
    !| summary: The entry point of a two-tier pair, which allocates
    !| author: A Developer
    subroutine fx_cluster(values, n_values, n_clusters, ierr)
        integer(int32), intent(in) :: n_values
            !! elements of `values`
        real(real64), dimension(n_values), intent(in) :: values
            !! the values
        integer(int32), intent(out) :: n_clusters
            !! how many clusters were found
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        n_clusters = count(values > 0.0_real64)
    end subroutine fx_cluster

    !> M_EXPORT_C
    !| summary: The expert tier of the pair, handed the work array instead of allocating it
    !| author: A Developer
    subroutine fx_cluster_expert(values, n_values, tmp_work, n_work, n_clusters, ierr)
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

        integer(int32) :: i

        call set_ok(ierr)
        ! use the work array so it is genuinely needed
        do i = 1, min(n_values, n_work)
            tmp_work(i) = values(i)
        end do
        n_clusters = count(values > 0.0_real64)
    end subroutine fx_cluster_expert

    !> M_EXPORT_C
    !| summary: A consumer whose producer names its input differently
    !| author: A Developer
    subroutine fx_renamed_input(samples, n_samples, tmp_work, n_work, ierr)
        integer(int32), intent(in) :: n_samples
            !! elements of `samples`
        real(real64), dimension(n_samples), intent(in) :: samples
            !! the values
        integer(int32), intent(in) :: n_work
            !! size of the work array.
            !! DM_OUTPUT_FROM(n_work, fx_work_size, fx_edges, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n_values       | n_samples   |
        real(real64), dimension(n_work), intent(out) :: tmp_work
            !! scratch space
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        tmp_work = real(n_work, real64)
    end subroutine fx_renamed_input

    !> M_EXPORT_C
    !| summary: A consumer whose producer takes an input it has no argument for
    !| author: A Developer
    subroutine fx_constant_input(samples, n_samples, tmp_work, n_work, ierr)
        integer(int32), intent(in) :: n_samples
            !! elements of `samples`
        real(real64), dimension(n_samples), intent(in) :: samples
            !! the values
        integer(int32), intent(in) :: n_work
            !! size of the work array.
            !! DM_OUTPUT_FROM(n_work, fx_scaled_size, fx_edges, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n_values       | n_samples   |
            !! | factor         | 3_int32     |
        real(real64), dimension(n_work), intent(out) :: tmp_work
            !! scratch space
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        tmp_work = real(n_work, real64)
    end subroutine fx_constant_input

    !> M_EXPORT_C
    !| summary: A sizer taking a factor the consumer has no argument for
    !| author: A Developer
    subroutine fx_scaled_size(n_values, factor, n_work, ierr)
        integer(int32), intent(in) :: n_values
            !! elements the caller intends to pass
        integer(int32), intent(in) :: factor
            !! how many work elements each value needs
        integer(int32), intent(out) :: n_work
            !! required size of the work array
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        n_work = factor*n_values
    end subroutine fx_scaled_size

    !> M_EXPORT_C
    !| summary: Works out the work array size for fx_cluster_expert
    !| author: A Developer
    subroutine fx_work_size(n_values, n_work, ierr)
        integer(int32), intent(in) :: n_values
            !! elements the caller intends to pass
        integer(int32), intent(out) :: n_work
            !! required size of the work array
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        n_work = 2*n_values
    end subroutine fx_work_size

    !> M_EXPORT_C
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

    !> M_EXPORT_C
    !| summary: A procedure that reports no errors of its own
    !| author: A Developer
    pure subroutine fx_no_ierr(value, doubled)
        real(real64), intent(in) :: value
            !! the input
        real(real64), intent(out) :: doubled
            !! twice the input
    end subroutine fx_no_ierr

end module fx_edges
