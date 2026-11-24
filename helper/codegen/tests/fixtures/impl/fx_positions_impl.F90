#include <src/macros.h>

!> summary: Implementation fixture for the argument position an implementation's `ierr` carries
!| author: A Developer
!| Every implementation and prologue here fails on request: it reports `ERR_INVALID_INPUT` at
!| `fx_report_pos` (or, for the prologue,
!| `fx_report_prologue_pos`), a position in its OWN dummy
!| list. The generated wrapper translates it into the caller's: a dummy it passes straight
!| through keeps its argument, everything it prepared itself is cleared to 0. Between them they
!| drop arguments ahead of kept ones in every way a wrapper can -- a permutation, a work array,
!| a recommend-sized value, a fixed mode, another mode's optional, a prologue-produced local, a
!| function result -- and keep a mask count and an optional. Compiled and run, every position
!| in turn, by `test_end_to_end_fortran.py`.
module fx_positions_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: fx_report_pos, fx_report_prologue_pos, MODE_FX_A, MODE_FX_B
    public :: fx_pos_work_size_impl, fx_pos_sum_impl, fx_pos_group_impl, fx_pos_clip_guard, fx_pos_clip_impl, fx_pos_share_impl

    integer(int32), save :: fx_report_pos = 0_int32
        !! the position each implementation here reports, 0 for none
    integer(int32), save :: fx_report_prologue_pos = 0_int32
        !! the position the prologue reports, 0 for none
    integer(int32), parameter :: MODE_FX_A = 1_int32
        !! mode a, which takes a weight
    integer(int32), parameter :: MODE_FX_B = 2_int32
        !! mode b, which takes none

contains

    !> summary: The size of `fx_pos_sum_impl`'s scratch: twice the number of values
    !| author: A Developer
    pure subroutine fx_pos_work_size_impl(n, work_size)
        integer(int32), intent(in) :: n
            !! number of values
        integer(int32), intent(out) :: work_size
            !! twice `n`
        work_size = 2_int32*n
    end subroutine fx_pos_work_size_impl

    !> summary: A sum over selected values, failing at `fx_report_pos` on request
    !| author: A Developer
    !| An edge case: a permutation, a work array and the recommend-sized value sizing it come
    !| before a mask, its count and an optional, so `fx_pos_sum` numbers those three otherwise.
    pure subroutine fx_pos_sum_impl(n, values, values_perm, work_size, tmp_work, values_mask, &
                                    n_selected_values, scale, total, ierr)
        integer(int32), intent(in) :: n
            !! number of values
        real(real64), dimension(n), intent(in) :: values
            !! the values
        integer(int32), dimension(n), intent(in) :: values_perm
            !! their ascending order
        integer(int32), intent(in) :: work_size
            !! size of the scratch
            !! DM_OUTPUT_FROM(work_size, fx_pos_work_size, fx_positions, AUTO)
        real(real64), dimension(work_size), intent(out) :: tmp_work
            !! scratch
        logical(c_bool), dimension(n), intent(in) :: values_mask
            !! which values count
        integer(int32), intent(in) :: n_selected_values
            !! how many count
        real(real64), intent(in), optional :: scale
            !! a factor; omitted, one
        real(real64), intent(out) :: total
            !! the sum of the selected values, times `scale`
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_INVALID_INPUT` at `fx_report_pos` when that is set
        call set_ok(ierr)
        tmp_work = 0.0_real64
        total = sum(values(values_perm), mask=values_mask(values_perm))
        if (present(scale)) total = total*scale
        if (fx_report_pos > 0) call set_err(ierr, ERR_INVALID_INPUT, fx_report_pos)
    end subroutine fx_pos_sum_impl

    !> summary: A sum in one of two modes, failing at `fx_report_pos` on request
    !| author: A Developer
    !| An edge case: a mode split, the mode and a work array ahead of an optional only mode a
    !| takes, each wrapper of each mode numbering what follows its own way.
    pure subroutine fx_pos_group_impl(n, values, mode, tmp_work, weight, total, ierr)
        integer(int32), intent(in) :: n
            !! number of values
        real(real64), dimension(n), intent(in) :: values
            !! the values
        integer(int32), intent(in) :: mode
            !! how to sum
            !!
            !! | Mode | Value | Procedure |
            !! |------|-------|-----------|
            !! | a | [[fx_positions_impl(module):MODE_FX_A(variable)]] | fx_pos_group_a |
            !! | b | [[fx_positions_impl(module):MODE_FX_B(variable)]] | fx_pos_group_b |
        real(real64), dimension(n), intent(out) :: tmp_work
            !! scratch
        real(real64), intent(in), optional :: weight
            !! a factor
            !! DM_REQUIRED_IF_MODE(mode, fx_positions_impl, MODE_FX_A)
        real(real64), intent(out) :: total
            !! the sum, times `weight` in mode a
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_INVALID_INPUT` at `fx_report_pos` when that is set
        call set_ok(ierr)
        tmp_work = values
        total = sum(tmp_work)
        if (mode == MODE_FX_A) total = total*weight
        if (fx_report_pos > 0) call set_err(ierr, ERR_INVALID_INPUT, fx_report_pos)
    end subroutine fx_pos_group_impl

    !> summary: The prologue of `fx_pos_clip_impl`, failing at `fx_report_prologue_pos` on request
    !| author: A Developer
    pure subroutine fx_pos_clip_guard(n, values, fraction, level, handled, ierr)
        integer(int32), intent(in) :: n
            !! number of values
        real(real64), dimension(n), intent(in) :: values
            !! the values
        real(real64), intent(in) :: fraction
            !! the level as a fraction of the largest value
        real(real64), intent(out) :: level
            !! the level to clip at
        logical, intent(out) :: handled
            !! never: the implementation always runs
        integer(int32), intent(inout) :: ierr
            !! Error code; `ERR_INVALID_INPUT` at `fx_report_prologue_pos` when that is set
        handled = .false.
        level = 0.0_real64
        if (n > 0) level = fraction*maxval(values)
        if (fx_report_prologue_pos > 0) call set_err(ierr, ERR_INVALID_INPUT, fx_report_prologue_pos)
    end subroutine fx_pos_clip_guard

    !> summary: Values clipped at a level, failing at `fx_report_pos` on request
    !| author: A Developer
    !| An edge case: a prologue, whose own argument joins the signature before `ierr`, produces
    !| the level ahead of a work array and the result.
    !| DM_PROLOGUE(fx_pos_clip_guard, fx_positions_impl)
    pure subroutine fx_pos_clip_impl(n, values, level, tmp_work, clipped, ierr)
        integer(int32), intent(in) :: n
            !! number of values
        real(real64), dimension(n), intent(in) :: values
            !! the values
        real(real64), intent(in) :: level
            !! the level to clip at
        real(real64), dimension(n), intent(out) :: tmp_work
            !! scratch
        real(real64), dimension(n), intent(out) :: clipped
            !! `values`, clipped at `level`
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_INVALID_INPUT` at `fx_report_pos` when that is set
        call set_ok(ierr)
        tmp_work = min(values, level)
        clipped = tmp_work
        if (fx_report_pos > 0) call set_err(ierr, ERR_INVALID_INPUT, fx_report_pos)
    end subroutine fx_pos_clip_impl

    !> summary: The sum of the values, failing at `fx_report_pos` on request
    !| author: A Developer
    !| An edge case: a function, whose result the wrappers publish where the implementation has
    !| no dummy, after a work array.
    function fx_pos_share_impl(n, values, tmp_copy, ierr) result(share)
        integer(int32), intent(in) :: n
            !! number of values
        real(real64), dimension(n), intent(in) :: values
            !! the values
        real(real64), dimension(n), intent(out) :: tmp_copy
            !! scratch
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_INVALID_INPUT` at `fx_report_pos` when that is set
        real(real64) :: share
            !! the sum
        call set_ok(ierr)
        tmp_copy = values
        share = sum(tmp_copy)
        if (fx_report_pos > 0) call set_err(ierr, ERR_INVALID_INPUT, fx_report_pos)
    end function fx_pos_share_impl
end module fx_positions_impl
