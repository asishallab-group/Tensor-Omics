#include <src/macros.h>

!> AUTHOR_FRANZ_ERIC_SILL
!| This module ensures equivalence of used c types in this framework during compile time.
!| Thus, it needs to be used by every module to be compiled first, ***and*** it is highly recommended to compile the framework manually.
module f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64

    ! safeguard to guarantee identity of c kinds and fortran kinds
    ! The preprocessor directives enforce a mismatch by overriding the C kinds
    ! Thus, in the final else-block all are used from iso_c_binding
    ! Using extra modules lowers the compilation priority of this module -> some other modules will be compiled first -> if they use c kinds but not safeguard, they fail first -> not wanted
#ifdef TEST_KIND_MISMATCH_C_INT
    use tox_conversions, only: c_char_as_char
    use f42_config
    use, intrinsic :: iso_c_binding, only: c_double, c_double_complex
    M_IMPLICIT_NONE
    integer(int32), parameter :: c_int = int32*2
#else
#ifdef TEST_KIND_MISMATCH_C_DOUBLE
    use tox_conversions, only: c_char_as_char
    use f42_config
    use, intrinsic :: iso_c_binding, only: c_int, c_double_complex
    M_IMPLICIT_NONE
    integer(int32), parameter ::  c_double = real64*2
#else
#ifdef TEST_KIND_MISMATCH_C_DOUBLE_COMPLEX
    use tox_conversions, only: c_char_as_char
    use f42_config
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_IMPLICIT_NONE
    integer(int32), parameter ::  c_double_complex = real64*2
#else
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
    M_IMPLICIT_NONE
#endif
#endif
#endif

    ! type guards to guarantee kind identity between fortran and c for correct interop in the c wrapper routines
    ! NOTE: `1/merge(1, 0, cond)` is a compile-time-evaluated "static assert": when `cond` is `.false.`,
    ! `merge` yields 0 and the integer division by zero is a constant-expression error, which most
    ! compilers reject at compile time rather than at runtime -- turning a kind mismatch into a build
    ! failure instead of silent undefined behavior in the C interop wrappers.
    logical, parameter :: THIS_FAILS_IF_C_INT_DOES_NOT_MATCH_INT32 = 1 == 1/merge(1, 0, c_int == int32)
        !! Compile-time guard: fails to compile if `c_int` is not kind-identical to `int32`.
    logical, parameter :: THIS_FAILS_IF_C_DOUBLE_DOES_NOT_MATCH_REAL64 = 1 == 1/merge(1, 0, c_double == real64)
        !! Compile-time guard: fails to compile if `c_double` is not kind-identical to `real64`.
    logical, parameter :: THIS_FAILS_IF_C_DOUBLE_COMPLEX_DOES_NOT_MATCH_REAL64 = 1 == 1/merge(1, 0, c_double_complex == real64)
        !! Compile-time guard: fails to compile if the real component kind of `c_double_complex` is not kind-identical to `real64`.
end module f42_safeguard
