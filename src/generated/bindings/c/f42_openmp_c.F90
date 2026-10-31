#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_openmp(module)]]
!| What this build of the library does with OpenMP.
!|
!| OpenMP is a build option (`./build.sh --openmp`), so whether a call runs on several threads
!| is decided when the library is built, not by the caller. A procedure that can run on several
!| threads takes their number as an argument; this module answers how many a build uses by default.
module f42_openmp_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: max_omp_threads_c

contains

    !> summary: C-wrapper for [[f42_openmp(module):max_omp_threads(function)]]
    !| A build with OpenMP answers the OpenMP runtime's `omp_get_max_threads()`, which follows
    !| the environment variable `OMP_NUM_THREADS` and is otherwise the runtime's own choice,
    !| usually the number of processors.
    !| It is a default, not a limit: a procedure's `n_omp_threads` sets the threads of its loop,
    !| and may ask for more.
    !| A serial build answers 1, as one thread runs every loop there.
    subroutine max_omp_threads_c(&
            max_threads,&
            ierr&
        ) bind(C, name="max_omp_threads_c")
        use f42_openmp, only: max_omp_threads

        integer(c_int), intent(out), target :: max_threads
            !! The default number of threads of a parallel loop; 1 without OpenMP
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_threads)

        max_threads = max_omp_threads(&
        )
    end subroutine max_omp_threads_c

end module f42_openmp_c
#endif
