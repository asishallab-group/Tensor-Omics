#include <src/macros.h>

!> File I/O helpers.
!|
!| Not gathered by [[f42_utils_impl(module)]] -- `open_file` is internal Fortran-to-Fortran
!| plumbing (no `M_EXPORT_C`, not an `_impl`), used directly by the hand-written export-path
!| modules that need it (`tox_stc_json`, `tox_stc_csv`, `tox_flyer_json`), so it stays outside
!| the generated-utilities family rather than routing through it.
module f42_io
    use, intrinsic :: iso_fortran_env, only: int32
    use tox_errors, only: set_ok, set_err, is_err, ERR_FILE_OPEN
    M_IMPLICIT_NONE

    private
    public :: open_file

contains

    !> Opens `filename` on a new unit as a stream, replacing any existing file.
    !| Sets `ERR_FILE_OPEN` on failure. Use `formatted=.true.` for text output, `.false.` for binary.
    subroutine open_file(filename, unit, formatted, ierr)
        character(len=*), intent(in) :: filename
            !! Path of the file to open
        integer(int32), intent(out) :: unit
            !! Newly allocated unit connected to `filename`
        logical, intent(in) :: formatted
            !! `.true.` for formatted (text), `.false.` for unformatted (binary) access
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        if (formatted) then
            open(newunit=unit, file=filename, form='formatted', access='stream', status='replace', iostat=ierr)
        else
            open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace', iostat=ierr)
        end if
        if (is_err(ierr)) then
            call set_err(ierr, ERR_FILE_OPEN)
        end if
    end subroutine open_file

end module f42_io
