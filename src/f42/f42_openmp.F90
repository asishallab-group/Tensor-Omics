#include <src/macros.h>

!> What this build of the library does with OpenMP.
!|
!| OpenMP is a build option (`./build.sh --openmp`), so whether a call runs on several threads
!| is decided when the library is built, not by the caller. A procedure that can run on several
!| threads takes their number as an argument; this module answers how many a build can use.
module f42_openmp
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32
    M_IMPLICIT_NONE
    private

    public :: max_omp_threads

contains

    !> M_EXPORT_C
    !| summary: The number of threads a parallel loop of this build can use: 1 without OpenMP.
    !| AUTHOR_FRANZ_ERIC_SILL
    !| A build with OpenMP answers the OpenMP runtime's `omp_get_max_threads()`, which follows
    !| the environment variable `OMP_NUM_THREADS` and is otherwise the runtime's own choice,
    !| usually the number of processors.
    !| A serial build answers 1, as one thread runs every loop there.
    function max_omp_threads() result(max_threads)
        !$ use omp_lib, only: omp_get_max_threads
        integer(int32) :: max_threads
            !! The number of threads a parallel loop can use; 1 without OpenMP

        max_threads = 1_int32
        !$ max_threads = int(omp_get_max_threads(), int32)
    end function max_omp_threads

end module f42_openmp
