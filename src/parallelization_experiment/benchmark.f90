module c_time
   use iso_c_binding
   implicit none

   integer(c_int), parameter :: CLOCK_MONOTONIC = 1

   type, bind(C) :: timespec
      integer(c_long) :: tv_sec   ! seconds
      integer(c_long) :: tv_nsec  ! nanoseconds
   end type timespec

   interface
      subroutine clock_gettime(clk_id, tp) bind(C, name="clock_gettime")
         import :: c_int, timespec
         integer(c_int), intent(in), value :: clk_id
         type(timespec), intent(out) :: tp
      end subroutine
   end interface

contains

   function get_time() result(time)
      type(timespec) :: time

      call clock_gettime(CLOCK_MONOTONIC, time)   
   end function get_time
   
   function duration(start_time, stop_time) result(dur)
      use iso_fortran_env, only: int32, real64
      implicit none

      type(timespec), intent(in) :: start_time, stop_time
      real(real64) :: dur

      dur = real(stop_time%tv_sec - start_time%tv_sec, kind=real64)
      dur = dur + real(stop_time%tv_nsec - start_time%tv_nsec, kind=real64) / 1.0e9_real64
   end function duration

end module


program benchmark_parallelization
   use iso_fortran_env, only: int32, real64, stdout=>output_unit
   use c_time
   use parallelization_experiment
   implicit none

   integer(int32), parameter :: n_times_per_exp = 10, n(2) = [6000, 3000], m(2) = [4096, 16384]  ! matrix dimensions
   real(real64), allocatable :: x(:,:), y(:,:), out_mat(:,:)

   real(real64) :: time
   integer(int32) :: i, i_dim, ierr
   type(timespec) :: start, finish
   logical :: correct_result

   print "(A)", "Experiment,n_rows,n_cols,Time(sec),Correct_Result"
   flush(stdout)

   do i_dim = 1, 2
      ! Allocate and initialize matrices
      allocate(x(n(i_dim), m(i_dim)), y(n(i_dim), m(i_dim)), out_mat(n(i_dim), m(i_dim)), stat=ierr)
      if (ierr /= 0) error stop "could not allocate"
      x = 1.0_real64
      y = 0.5_real64


      !!!! run elemental experiment
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call elemental_calc(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "elemental_calc", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!

      !!!! run dc_collapsed
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call dc_collapsed(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "dc_collapsed", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run omp_collapsed
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call omp_collapsed(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "omp_collapsed", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run dc_dc
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call dc_dc(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "dc_dc", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run omp_dc
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call omp_dc(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "omp_dc", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run osimd_dc
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call osimd_dc(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "osimd_dc", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run ocon_dc
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call ocon_dc(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "ocon_dc", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run osimd_dc
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call osimd_dc(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "osimd_dc", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run dc_osimd
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call dc_osimd(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "dc_osimd", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run omp_osimd
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call omp_osimd(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "omp_osimd", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run ocon_osimd
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call ocon_osimd(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "ocon_osimd", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run osimd_osimd
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call osimd_osimd(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "osimd_osimd", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run dc_bare
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call dc_bare(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "dc_bare", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run omp_bare
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call omp_bare(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "omp_bare", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run ocon_bare
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call ocon_bare(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "ocon_bare", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run osimd_bare
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call osimd_bare(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "osimd_bare", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      !!!! run bare
         time = 0.0_real64
         correct_result = .true.
         do i = 1, n_times_per_exp
            start = get_time()
            call bare(x, y, out_mat)
            finish = get_time()

            correct_result = correct_result .and. all(out_mat == log(1.0_real64 + exp(0.5_real64)))

            time = time + duration(start, finish)

         end do
         time = time / real(n_times_per_exp, kind=real64)
         print '(A, ",", I0, ",", I0, ",", G0, ",", A)', "bare", size(out_mat, 1), size(out_mat, 2), time, merge("T", "F", correct_result)
         call flush(stdout)
      !!!!
      deallocate(x, y, out_mat, stat=ierr)
      if (ierr /= 0) error stop "could not deallocate"
   end do
end program benchmark_parallelization
