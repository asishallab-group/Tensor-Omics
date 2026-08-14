#include <src/macros.h>

!> Module with conversion helpers between Fortran-native types and their `iso_c_binding` counterparts.
!| Used by the C interface wrapper subroutines throughout tensor-omics to translate arguments (logicals,
!| characters/strings, and fixed-width integer kinds) crossing the Fortran/C boundary.
!|
!| It keeps its `tox_` name for continuity, but it lives in `src/f42/` because it is C-interop
!| glue with no exports of its own -- infrastructure the binding layer stands on.
module tox_conversions
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_null_char, c_char, c_size_t, c_int64_t
    use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

contains

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Converts a c_int value to logical, elemental -> any shape
    elemental subroutine c_int_as_logical(c_val, f_val)
        integer(c_int), intent(in) :: c_val
            !! the element containing the c variant of the number
        logical, intent(out) :: f_val
            !! the element that will hold the logical representation of the c_int: `0` means `.false.`, else `.true.`

        f_val = c_val /= 0
    end subroutine c_int_as_logical

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Converts a logical value to c_int, elemental -> any shape
    elemental subroutine logical_as_c_int(f_val, c_val)
        logical, intent(in) :: f_val
            !! the element containing the fortran variant of the number
        integer(c_int), intent(out) :: c_val
            !! the element that will hold the c_int representation of the logical: `0` if `.false.`, else `1`

        if (f_val) then
            c_val = 1
        else
            c_val = 0
        end if
    end subroutine logical_as_c_int

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Converts a c_char value to fortran character, elemental -> any shape
    elemental subroutine c_char_as_char(c_val, f_val)
        character(kind=c_char, len=1), intent(in) :: c_val
            !! the element containing the c variant of the character
        character(len=1), intent(out) :: f_val
            !! the element that will hold the fortran char representation of the c_char

        f_val = c_val
    end subroutine c_char_as_char

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Converts a character value to c_char, elemental -> any shape
    elemental subroutine char_as_c_char(f_val, c_val)
        character(len=1), intent(in) :: f_val
            !! the element containing the fortran variant of the character
        character(kind=c_char, len=1), intent(out) :: c_val
            !! the element that will hold the c_char (ASCII) representation of the character

        c_val = f_val
    end subroutine char_as_c_char

    !> AUTHOR_FRANZ_ERIC_SILL
    !| A zero-copy Fortran view of a `c_char` buffer, ending at the first NUL if there is one.
    !|
    !| Accepts both conventions a caller might use, which is the point of it: a C caller
    !| NUL-terminates (`ward` then four NULs), the generated bindings blank-pad (`ward    `),
    !| and both must select the same mode. A NUL is honoured where there is one, and where
    !| there is none the whole buffer is the string -- which compares equal to the shorter
    !| literal anyway, because Fortran blank-pads the shorter operand.
    !|
    !| The result is a **view**, not a copy: it is valid only while `buffer` is, and writing
    !| through `buffer` changes it. `buffer` must be contiguous, or `c_loc` takes the address
    !| of a compiler temporary and the view dangles on return. Every caller here passes an
    !| explicit-shape dummy or a `c_f_pointer`-mapped array, both of which are contiguous.
    !|
    !| Replaced an allocating version. Nothing here allocates now, which is what lets the
    !| module sit on `Conventions.impl_import_whitelist` and be provably allocation-free.
    function c_char_as_view(buffer) result(view)
        character(kind=c_char, len=1), dimension(:), intent(in), target, contiguous :: buffer
            !! the C buffer, NUL-terminated or blank-padded
        character(len=:), pointer :: view
            !! the buffer read as one Fortran string, up to the first NUL

        ! `size` of a dummy is a specification expression, so this length is legal even though
        ! it is not known until the call. What is NOT legal is a length computed inside the
        ! procedure -- hence the substring pointer assignment below rather than a second view.
        character(len=size(buffer)), pointer :: whole
        integer(int32) :: n

        call c_f_pointer(c_loc(buffer), whole)
        n = index(whole, c_null_char) - 1
        if (n < 0) n = len(whole)          ! no NUL: the whole buffer is the string
        ! pointer assignment, NOT c_f_pointer -- c_f_pointer onto a deferred length compiles
        ! and then yields an empty string on gfortran and a segfault on ifx
        view => whole(1:n)
    end function c_char_as_view

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Converts int32 to c_int64_t, elemental -> any shape
    elemental subroutine int32_as_c_int64(f_val, c_val)
        integer(int32), intent(in) :: f_val
            !! the element containing the fortran variant of the number
        integer(c_int64_t), intent(out) :: c_val
            !! the element that will hold the c_int64_t representation

        c_val = int(f_val, kind=c_int64_t)
    end subroutine int32_as_c_int64

    !> AUTHOR_AARON_SCHROEDER
    !| Converts c_int64_t to int32, elemental -> any shape
    elemental subroutine c_int64_as_int32(c_val, f_val)
        integer(c_int64_t), intent(in) :: c_val
            !! the element containing the c variant of the number
        integer(int32), intent(out) :: f_val
            !! the element that will hold the int32 representation

        f_val = int(c_val, kind=int32)
    end subroutine c_int64_as_int32

    !> AUTHOR_AARON_SCHROEDER
    !| Converts int32 to c_size_t, elemental -> any shape
    elemental subroutine int32_as_c_size(f_val, c_val)
        integer(int32), intent(in) :: f_val
            !! the element containing the fortran variant of the number
        integer(c_size_t), intent(out) :: c_val
            !! the element that will hold the c_size_t representation

        c_val = int(f_val, kind=c_size_t)
    end subroutine int32_as_c_size

    !> AUTHOR_AARON_SCHROEDER
    !| Converts c_size_t to int32, elemental -> any shape
    elemental subroutine c_size_as_int32(c_val, f_val)
        integer(c_size_t), intent(in) :: c_val
            !! the element containing the c variant of the number
        integer(int32), intent(out) :: f_val
            !! the element that will hold the int32 representation

        f_val = int(c_val, kind=int32)
    end subroutine c_size_as_int32

    !> AUTHOR_AARON_SCHROEDER
    !| Converts c_int64_t to c_size_t, elemental -> any shape
    elemental subroutine c_int64_as_c_size(c64_val, csize_val)
        integer(c_int64_t), intent(in) :: c64_val
            !! the element containing the c_int64_t value
        integer(c_size_t), intent(out) :: csize_val
            !! the element that will hold the c_size_t representation

        csize_val = int(c64_val, kind=c_size_t)
    end subroutine c_int64_as_c_size

    !> AUTHOR_AARON_SCHROEDER
    !| Converts c_size_t to c_int64_t, elemental -> any shape
    elemental subroutine c_size_as_c_int64(csize_val, c64_val)
        integer(c_size_t), intent(in) :: csize_val
            !! the element containing the c_size_t value
        integer(c_int64_t), intent(out) :: c64_val
            !! the element that will hold the c_int64_t representation

        c64_val = int(csize_val, kind=c_int64_t)
    end subroutine c_size_as_c_int64

end module tox_conversions
