#include <authors.h>
#include <src/f42/serde/macros.h>
#include <src/macros.h>

!> AUTHOR_FRANZ_ERIC_SILL
!| Module for de-/serialization utilities.
module f42_serde_utils
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32
    M_IMPLICIT_NONE

    integer(int32), parameter :: INTEGER_TYPE_CODE = M_INTEGER_TYPE_CODE
        !! Type code for integer type
    integer(int32), parameter :: REAL_TYPE_CODE = M_REAL_TYPE_CODE
        !! Type code for real type
    integer(int32), parameter :: LOGICAL_TYPE_CODE = M_LOGICAL_TYPE_CODE
        !! Type code for logical type
    integer(int32), parameter :: COMPLEX_TYPE_CODE = M_COMPLEX_TYPE_CODE
        !! Type code for complex type
end module f42_serde_utils