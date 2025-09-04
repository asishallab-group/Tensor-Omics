module tox_conversions
    use iso_fortran_env, only: int32, real64
    use iso_c_binding, only: c_int, c_double, c_null_char, c_double_complex, c_char
    use tox_errors, only: ERR_ALLOC_FAIL, is_err, set_ok, set_err
    implicit none

contains

    !> Converts a c_double value to real64, elemental -> any shape
    elemental subroutine c_double_as_real64(c_val, f_val)
        real(c_double), intent(in) :: c_val
         !! the element containing the c variant of the number
        real(real64), intent(out) :: f_val
         !! the element that will hold the real64 representation of the c_double

        f_val = real(c_val, kind=real64)
    end subroutine c_double_as_real64

    !> Converts a real64 value to c_double, elemental -> any shape
    elemental subroutine real64_as_c_double(f_val, c_val)
        real(real64), intent(in) :: f_val
         !! the element containing the fortran variant of the number
        real(c_double), intent(out) :: c_val
         !! the element that will hold the c_double representation of the real64

        c_val = real(f_val, kind=c_double)
    end subroutine real64_as_c_double

    !> Converts a c_double_complex value to fortran complex(real64), elemental -> any shape
    elemental subroutine c_complex_as_complex(c_val, f_val)
        complex(c_double_complex), intent(in) :: c_val
         !! the element containing the c variant of the number
        complex(real64), intent(out) :: f_val
         !! the element that will hold the complex(real64) representation of the c_double_complex

        f_val = cmplx(real(c_val), aimag(c_val), kind=real64)
    end subroutine c_complex_as_complex

    !> Converts a fortran complex(real64) value to c_double_complex, elemental -> any shape
    elemental subroutine complex_as_c_complex(f_val, c_val)
        complex(real64), intent(in) :: f_val
         !! the element containing the fortran variant of the number
        complex(c_double_complex), intent(out) :: c_val
         !! the element that will hold the c_double_complex representation of the complex(real64)

        c_val = cmplx(real(f_val), aimag(f_val), kind=c_double_complex)
    end subroutine complex_as_c_complex

    !> Converts a c_int value to int32, elemental -> any shape
    elemental subroutine c_int_as_int32(c_val, f_val)
        integer(c_int), intent(in) :: c_val
         !! the element containing the c variant of the number
        integer(int32), intent(out) :: f_val
         !! the element that will hold the int32 representation of the c_int

        f_val = int(c_val, kind=int32)
    end subroutine c_int_as_int32

    !> Converts a int32 value to c_int, elemental -> any shape
    elemental subroutine int32_as_c_int(f_val, c_val)
        integer(int32), intent(in) :: f_val
         !! the element containing the fortran variant of the number
        integer(c_int), intent(out) :: c_val
         !! the element that will hold the c_int representation of the int32

        c_val = int(f_val, kind=c_int)
    end subroutine int32_as_c_int

    !> Converts a c_int value to logical, elemental -> any shape
    elemental subroutine c_int_as_logical(c_val, f_val)
        integer(c_int), intent(in) :: c_val
         !! the element containing the c variant of the number
        logical, intent(out) :: f_val
         !! the element that will hold the logical representation of the c_int: `0` means `.false.`, else `.true.`

        f_val = c_val /= 0
    end subroutine c_int_as_logical

    !> Converts a logical value to c_int, elemental -> any shape
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

    !> Converts a c_char value to fortran character, elemental -> any shape
    elemental subroutine c_char_as_char(c_val, f_val)
        character(kind=c_char, len=1), intent(in) :: c_val
         !! the element containing the c variant of the character
        character(len=1), intent(out) :: f_val
         !! the element that will hold the fortran char representation of the c_char

         f_val = c_val
    end subroutine c_char_as_char

    !> Converts a character value to c_char, elemental -> any shape
    elemental subroutine char_as_c_char(f_val, c_val)
        character(len=1), intent(in) :: f_val
         !! the element containing the fortran variant of the character
        character(kind=c_char, len=1), intent(out) :: c_val
         !! the element that will hold the c_char (ASCII) representation of the character

        c_val = f_val
    end subroutine char_as_c_char

    !> Converts a 1D c_char array to string
    pure subroutine c_char_1d_as_string(c_char_array, str_out, ierr)
        character(kind=c_char, len=1), dimension(:), intent(in) :: c_char_array
         !! c int array, representing characters
        character(len=:), allocatable, intent(out) :: str_out
         !! Fortran string, length determined by occuring null char in `c_char_array`
        integer(int32), intent(out) :: ierr
         !! Error code

        integer(int32) :: i, str_len, array_len

        call set_ok(ierr)

        array_len = size(c_char_array, 1)

        ! identify string length
        str_len = array_len
        do i = 1, array_len
            if (c_char_array(i) == c_null_char) then
                str_len = i - 1
                exit
            end if
        end do

        ! create string
        allocate(character(len=str_len) :: str_out, stat=ierr)
        if (is_err(ierr)) then
           call set_err(ierr, ERR_ALLOC_FAIL)
           return
        end if

        do i = 1, str_len
            call c_char_as_char(c_char_array(i), str_out(i:i))
        end do
    end subroutine c_char_1d_as_string

    !> Converts a string to 1D c_char array
    pure subroutine string_as_c_char_1d(str, c_char_array)
        character(len=*), intent(in) :: str
         !! Fortran string to be converted
        character(kind=c_char, len=1), dimension(:), intent(out) :: c_char_array
         !! c int array, representing chars of `str`, will always end with null char. If array too small, it will hold fitting trimmed `str`.

        integer(int32) :: i_str, str_len

        ! determine string length to be converted, +1 because of the null char in c
        str_len = min(len_trim(str) + 1, size(c_char_array, 1)) - 1

        do i_str = 1, str_len
           call c_char_as_char(str(i_str:i_str), c_char_array(i_str))
        end do

        if (size(c_char_array, 1) > 0) then
           c_char_array(str_len + 1) = c_null_char
        end if
    end subroutine string_as_c_char_1d

    !> Converts a 2D c_char array to 1D string array
    pure subroutine c_char_2d_as_string(c_char_array, str_out, ierr)
        character(kind=c_char, len=1), dimension(:, :), intent(in) :: c_char_array
         !! c int array, columns as ascii arrays
        character(len=:), dimension(:), allocatable, intent(out) :: str_out
         !! Fortran array of resulting strings
        integer(int32), intent(out) :: ierr
         !! Error code

        integer(int32) :: i_str, n_rows, n_strings
        character(len=:), allocatable :: string

        call set_ok(ierr)

        n_rows = size(c_char_array, 1)
        n_strings = size(c_char_array, 2)

        allocate(character(len=n_rows) :: str_out(n_strings), stat=ierr)
        if (is_err(ierr)) then
           call set_err(ierr, ERR_ALLOC_FAIL)
           return
        end if

        ! create strings
        do i_str = 1, n_strings
            call c_char_1d_as_string(c_char_array(:, i_str), string, ierr)
            if (is_err(ierr)) return
            str_out(i_str) = string
        end do
    end subroutine c_char_2d_as_string

    !> Converts a 1D string array to 2D c_char array
    pure subroutine string_as_c_char_2d(strings, c_char_array)
        character(len=*), dimension(:), intent(in) :: strings
         !! Fortran array of strings
        character(kind=c_char, len=1), dimension(:, :), intent(out) :: c_char_array
         !! c int array, columns as ascii arrays

        integer(int32) :: i_str, n_strings

        do i_str = 1, size(strings, 1)
           call string_as_c_char_1d(strings(i_str), c_char_array(:, i_str))
        end do
    end subroutine string_as_c_char_2d

end module tox_conversions
