#include <src/macros.h>

!> summary: Implementation fixture for OpenMP roots
!| author: A Developer
!| Two roots: implementations that take `n_omp_threads` and hold the parallel loop, declared pure
!| without OpenMP and impure with it. One has a per-thread work array, so its plain wrapper
!| allocates the thread axis with one column without OpenMP and hands the implementation the
!| same count; the other has none, so its only wrapper passes the count on. (No macro is named
!| here: a doc comment is preprocessed, so its name would be replaced.) Compiled and run by
!| `test_end_to_end_fortran.py`, and through C, Python and R by `test_end_to_end_impl.py`, which
!| builds it once more with OpenMP.
module fx_omp_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private

    public :: fx_column_sums_impl, fx_column_maxima_impl

contains

    !> summary: The sum of every column, one column per iteration of the parallel loop
    !| author: A Developer
    !| Each column is first copied into the calling thread's column of the work array, so a sum
    !| reads only what its own thread wrote. `n_work_columns` reports how many columns the work
    !| array has, which is how a test sees what the plain wrapper allocated.
    M_OMP_PURE subroutine fx_column_sums_impl(n_rows, n_cols, values, sums, n_omp_threads, &
        n_work_columns, tmp_column)
        !$ use omp_lib, only: omp_get_thread_num
        integer(int32), intent(in) :: n_rows
            !! rows of the matrix
        integer(int32), intent(in) :: n_cols
            !! columns of the matrix
        real(real64), intent(in) :: values(n_rows, n_cols)
            !! the matrix
        real(real64), intent(out) :: sums(n_cols)
            !! the sum of every column
        integer(int32), intent(in) :: n_omp_threads
        integer(int32), intent(out) :: n_work_columns
            !! the columns of the work array: the thread count the implementation was handed
        real(real64), intent(out) :: tmp_column(n_rows, n_omp_threads)
            !! one column of scratch per thread
        integer(int32) :: i_col, i_thread

        n_work_columns = n_omp_threads
        !M_OMP_PARALLEL(n_cols > 1) shared(n_rows, n_cols, values, sums, tmp_column) private(i_thread)
        do i_col = 1, n_cols
            i_thread = 1
            !$ i_thread = omp_get_thread_num() + 1
            tmp_column(:, i_thread) = values(:, i_col)
            sums(i_col) = sum(tmp_column(:, i_thread))
        end do
    end subroutine fx_column_sums_impl

    !> summary: The largest value of every column, one column per iteration of the parallel loop
    !| author: A Developer
    !| A root without work arrays: its one wrapper passes `n_omp_threads` on as it is.
    M_OMP_PURE subroutine fx_column_maxima_impl(n_rows, n_cols, values, maxima, n_omp_threads)
        integer(int32), intent(in) :: n_rows
            !! rows of the matrix
        integer(int32), intent(in) :: n_cols
            !! columns of the matrix
        real(real64), intent(in) :: values(n_rows, n_cols)
            !! the matrix
        real(real64), intent(out) :: maxima(n_cols)
            !! the largest value of every column
        integer(int32), intent(in) :: n_omp_threads
        integer(int32) :: i_col

        !M_OMP_PARALLEL(n_cols > 1) shared(n_cols, values, maxima)
        do i_col = 1, n_cols
            maxima(i_col) = maxval(values(:, i_col))
        end do
    end subroutine fx_column_maxima_impl

end module fx_omp_impl
