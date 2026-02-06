module safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64

    ! safeguard to guarantee identity of c kinds and fortran kinds
    ! The preprocessor directives enforce a mismatch by overriding the C kinds
    ! Thus, in the final else-block all are used from iso_c_binding
    ! Using extra modules lowers the compilation priority of this module -> some other modules will be compiled first -> if they use c kinds but not safeguard, they fail first -> not wanted
#ifdef TEST_KIND_MISMATCH_C_INT
    use tox_conversions, only: c_char_as_char
    use config
    use, intrinsic :: iso_c_binding, only: c_double, c_double_complex
    implicit none
    integer(int32), parameter :: c_int = int32 * 2
#else
#ifdef TEST_KIND_MISMATCH_C_DOUBLE
    use tox_conversions, only: c_char_as_char
    use config
    use, intrinsic :: iso_c_binding, only: c_int, c_double_complex
    implicit none
    integer(int32), parameter ::  c_double = real64 * 2
#else
#ifdef TEST_KIND_MISMATCH_C_DOUBLE_COMPLEX
    use tox_conversions, only: c_char_as_char
    use config
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none
    integer(int32), parameter ::  c_double_complex = real64 * 2
#else
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
    implicit none
#endif
#endif
#endif


    ! type guards to guarantee kind identity between fortran and c for correct interop in the c wrapper routines
    logical, parameter :: THIS_FAILS_IF_C_INT_DOES_NOT_MATCH_INT32 = 1 == 1 / merge(1, 0, c_int == int32)
    logical, parameter :: THIS_FAILS_IF_C_DOUBLE_DOES_NOT_MATCH_REAL64 = 1 == 1 / merge(1, 0, c_double == real64)
    logical, parameter :: THIS_FAILS_IF_C_DOUBLE_COMPLEX_DOES_NOT_MATCH_REAL64 = 1 == 1 / merge(1, 0, c_double_complex == real64)
end module safeguard