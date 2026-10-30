#include <src/macros.h>

!> summary: Implementation fixture for a prologue, beside an implementation's own `ierr`
!| author: A Developer
!| `fx_clip` clips values at a level its prologue derives from a fraction of the largest value.
!| The prologue asks for that fraction, which therefore joins `fx_clip`'s signature -- before
!| the implementation's own `ierr`, which still ends it. Read by `test_names.py` for the locals
!| an allocating wrapper with a prologue declares, and compiled by `test_end_to_end_fortran.py`.
module fx_prologue_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok
    M_IMPLICIT_NONE
    private

    public :: fx_clip_impl, fx_clip_guard

contains

    !> summary: The prologue of `fx_clip_impl`: derive the level from a fraction of the largest value
    !| author: A Developer
    pure subroutine fx_clip_guard(n, values, fraction, level, handled, ierr)
        integer(int32), intent(in) :: n
            !! number of values
        real(real64), dimension(n), intent(in) :: values
            !! the values
        real(real64), intent(in) :: fraction
            !! the level as a fraction of the largest value
        real(real64), intent(out) :: level
            !! the level to clip at
        logical, intent(out) :: handled
            !! never: the implementation always runs. A default logical, as the generated
            !! wrapper declares its own (the prologue contract in src/macros.h)
        integer(int32), intent(inout) :: ierr
            !! Error code, left as it is
        handled = .false.
        level = 0.0_real64
        if (n > 0) level = fraction * maxval(values)
        if (ierr /= 0) return
    end subroutine fx_clip_guard

    !> summary: Clip values at a level
    !| author: A Developer
    !| DM_PROLOGUE(fx_clip_guard, fx_prologue_impl)
    pure subroutine fx_clip_impl(n, values, level, tmp_work, clipped, ierr)
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
            !! Error code
        call set_ok(ierr)
        tmp_work = min(values, level)
        clipped = tmp_work
    end subroutine fx_clip_impl
end module fx_prologue_impl
