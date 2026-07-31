#include <src/macros.h>

!> summary: Fixture module covering the ordinary cases
!| Nothing here is compiled into the library. These procedures exist so the generator
!| has a stable, complete specimen of every construct it claims to support.
module fx_basics
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err, ERR_DIVISION_BY_ZERO, ERR_INVALID_INPUT
    implicit none

    integer(int32), parameter :: MODE_MEAN = 1
        !! average the values
    integer(int32), parameter :: MODE_MEDIAN = 2
        !! take the middle value

    integer(int32), parameter :: METHOD_WARD = 1
        !! minimises within-cluster variance
    integer(int32), parameter :: METHOD_SINGLE = 2
        !! nearest neighbour

contains

    !> M_EXPORT_C
    !| summary: Normalizes a vector to unit length in-place
    !| author: A Developer
    pure subroutine fx_normalize(vector, n_dims, ierr)
        integer(int32), intent(in) :: n_dims
            !! number of elements in `vector`
        real(real64), dimension(n_dims), intent(inout) :: vector
            !! Vector that will be normalized
        integer(int32), intent(out) :: ierr
            !! Error code

        real(real64) :: length

        call set_ok(ierr)
        length = sqrt(sum(vector**2))
        if (length <= 0.0_real64) then
            call set_err(ierr, ERR_DIVISION_BY_ZERO)
            return
        end if
        vector = vector/length
    end subroutine fx_normalize

    !> M_EXPORT_C
    !| summary: Sums a matrix, exercising a shared extent
    !| author: A Developer
    subroutine fx_sum_matrix(matrix, weights, n_rows, n_cols, total, ierr)
        integer(int32), intent(in) :: n_rows
            !! rows of `matrix`
        integer(int32), intent(in) :: n_cols
            !! columns of `matrix`, and elements of `weights`
        real(real64), dimension(n_rows, n_cols), intent(in) :: matrix
            !! the values
        real(real64), dimension(n_cols), intent(in) :: weights
            !! one weight per column, so its extent must match `matrix`
        real(real64), intent(out) :: total
            !! the weighted sum
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_col

        call set_ok(ierr)
        total = 0.0_real64
        do i_col = 1, n_cols
            total = total + weights(i_col)*sum(matrix(:, i_col))
        end do
    end subroutine fx_sum_matrix

    !> M_EXPORT_C
    !| summary: An output extent also named by a later input array, which must source it
    !| author: A Developer
    subroutine fx_grouped_output(n_cols, n_out, values, averaged, group_sizes, ierr)
        integer(int32), intent(in) :: n_cols
            !! columns of `values` and `averaged`
        integer(int32), intent(in) :: n_out
            !! output rows; also the number of `group_sizes`
        real(real64), dimension(n_cols), intent(in) :: values
            !! one value per column
        real(real64), dimension(n_out, n_cols), intent(out) :: averaged
            !! grouped output; its first extent `n_out` is shared with `group_sizes`
        integer(int32), dimension(n_out), intent(in) :: group_sizes
            !! one size per output group -- a *later* input carrying the same extent, so
            !! `n_out` must be read off it rather than the not-yet-allocated `averaged`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_out, i_col

        call set_ok(ierr)
        do i_col = 1, n_cols
            do i_out = 1, n_out
                averaged(i_out, i_col) = values(i_col)*real(group_sizes(i_out), real64)
            end do
        end do
    end subroutine fx_grouped_output

    !> M_EXPORT_C
    !| summary: An inout scalar the routine caps in place, handed back to the caller
    !| author: A Developer
    subroutine fx_cap_value(value, ierr)
        integer(int32), intent(inout) :: value
            !! capped at 10 and returned
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        if (value > 10_int32) value = 10_int32
    end subroutine fx_cap_value

    !> M_EXPORT_C
    !| summary: Every optional flavour at once
    !| author: A Developer
    subroutine fx_optionals(values, n_values, span, max_iter, use_quantile, tmp_work, ierr)
        integer(int32), intent(in) :: n_values
            !! elements of `values`
        real(real64), dimension(n_values), intent(in) :: values
            !! the values
        real(real64), intent(in), optional :: span
            !! smoothing span.
            !! DM_DEFAULT(0.1_real64)
        integer(int32), intent(in), optional :: max_iter
            !! iteration cap.
            !! DM_DEFAULT(300_int32)
        logical, intent(in), optional :: use_quantile
            !! whether to normalise by quantile.
            !! DM_DEFAULT(.false.)
        real(real64), dimension(n_values), intent(out) :: tmp_work
            !! scratch space, allocated by the caller
        integer(int32), intent(out) :: ierr
            !! Error code
    end subroutine fx_optionals

    !> M_EXPORT_C
    !| summary: A mode argument and a method argument
    !| author: A Developer
    subroutine fx_modes(values, n_values, mode, link_method, summary, ierr)
        integer(int32), intent(in) :: n_values
            !! elements of `values`
        real(real64), dimension(n_values), intent(in) :: values
            !! the values
        integer(int32), intent(in) :: mode
            !! how to summarise the values
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | average the values | [[fx_basics(module):MODE_MEAN(variable)]] |
            !! | take the middle value | [[fx_basics(module):MODE_MEDIAN(variable)]] |
        integer(int32), intent(in) :: link_method
            !! how to link clusters
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | minimises variance | [[fx_basics(module):METHOD_WARD(variable)]] |
            !! | nearest neighbour | [[fx_basics(module):METHOD_SINGLE(variable)]] |
        real(real64), intent(out) :: summary
            !! the summarised value
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        select case (mode)
            case (MODE_MEAN)
                summary = sum(values)/real(n_values, real64)
            case (MODE_MEDIAN)
                summary = values((n_values + 1)/2)
            case default
                call set_err(ierr, ERR_INVALID_INPUT)
                return
        end select
        ! link_method is accepted but unused: it is here to exercise a second mode
        summary = summary + 0.0_real64*real(link_method, real64)
    end subroutine fx_modes

    !> M_EXPORT_C
    !| summary: A function, which the C ABI turns into a subroutine
    !| author: A Developer
    pure function fx_count_positive(values, n_values) result(n_positive)
        integer(int32), intent(in) :: n_values
            !! elements of `values`
        real(real64), dimension(n_values), intent(in) :: values
            !! the values
        integer(int32) :: n_positive
            !! how many values are greater than zero

        n_positive = count(values > 0.0_real64)
    end function fx_count_positive

    !> M_EXPORT_C
    !| summary: A consumer whose work size is computed by a producer in another module
    !| author: A Developer
    !| The directive is documentation only, so naming fx_edges here creates no Fortran
    !| dependency: n_work is an ordinary input as far as the compiler is concerned.
    subroutine fx_cross_module(values, n_values, tmp_work, n_work, total, ierr)
        integer(int32), intent(in) :: n_values
            !! elements of `values`
        real(real64), dimension(n_values), intent(in) :: values
            !! the values
        integer(int32), intent(in) :: n_work
            !! size of the work array.
            !! DM_OUTPUT_FROM(n_work, fx_work_size, fx_edges, AUTO)
        real(real64), dimension(n_work), intent(out) :: tmp_work
            !! scratch space
        real(real64), intent(out) :: total
            !! the sum of the values
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i

        call set_ok(ierr)
        tmp_work = 0.0_real64
        do i = 1, min(n_values, n_work)
            tmp_work(i) = values(i)
        end do
        total = sum(values)
    end subroutine fx_cross_module

    !> summary: Not exported, and deliberately breaking the rules
    !| author: A Developer
    !| An internal routine is held to none of the export contract.
    subroutine fx_internal(text, mode)
        character(len=:), allocatable, intent(inout) :: text
            !! a deferred-length string, which an exported procedure could not have
        integer(int32), intent(in) :: mode
            !! a mode argument with no table, which an exported procedure could not have
    end subroutine fx_internal

end module fx_basics
