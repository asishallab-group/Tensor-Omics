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

   integer(int32), parameter :: n_times_per_exp = 10, n = 10000, m = 10000  ! matrix dimensions
   real(real64), allocatable :: x(:,:), y(:,:), out_mat(:,:)

   real(real64) :: time
   integer(int32) :: i
   type(timespec) :: start, finish
   logical :: correct_result

   ! Allocate and initialize matrices
   allocate(x(n, m), y(n, m), out_mat(n, m))
   x = 1.0_real64
   y = 0.5_real64

   print "(A)", "Experiment,Time(sec),Correct_Result"

   !!!! run elemental experiment
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call elemental_calc(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "elemental_calc", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!

   !!!! run dc_collapsed
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call dc_collapsed(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "dc_collapsed", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run omp_collapsed
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call omp_collapsed(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "omp_collapsed", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run dc_dc
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call dc_dc(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "dc_dc", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run omp_dc
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call omp_dc(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "omp_dc", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run osimd_dc
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call osimd_dc(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "osimd_dc", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run ocon_dc
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call ocon_dc(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "ocon_dc", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run osimd_dc
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call osimd_dc(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "osimd_dc", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run dc_osimd
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call dc_osimd(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "dc_osimd", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run omp_osimd
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call omp_osimd(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "omp_osimd", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run ocon_osimd
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call ocon_osimd(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "ocon_osimd", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run osimd_osimd
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call osimd_osimd(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "osimd_osimd", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run dc_bare
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call dc_bare(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "dc_bare", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run omp_bare
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call omp_bare(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "omp_bare", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run ocon_bare
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call ocon_bare(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "ocon_bare", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run osimd_bare
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call osimd_bare(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "osimd_bare", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
   !!!! run bare
      time = 0.0_real64
      correct_result = .true.
      do i = 1, n_times_per_exp
         start = get_time()
         call bare(x, y, out_mat)
         finish = get_time()

         correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

         time = time + duration(start, finish)

      end do
      time = time / real(n_times_per_exp, kind=real64)
      print '(A, ",", G0, ",", A)', "bare", time, merge("T", "F", correct_result)
      call flush(stdout)
   !!!!
end program benchmark_parallelization
