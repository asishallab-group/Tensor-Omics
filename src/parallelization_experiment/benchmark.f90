program benchmark_parallelization
   use iso_fortran_env, only: int32, real64
   use parallelization_experiment
   implicit none

   integer(int32), parameter :: n_times_per_exp = 10, n = 10000, m = 10000  ! matrix dimensions
   real(real64), allocatable :: x(:,:), y(:,:), out_mat(:,:)
   integer :: start, finish, rate
   integer :: i, j, k
   logical :: correct_result
   real(real64) :: time

   ! Procedure pointers for leaf and root routines
   procedure(leaf), pointer :: leaf_routine => null()

   interface
      subroutine root(x_mat, y_mat, out_mat, leaf_routine)
         use iso_fortran_env, only: real64
         use parallelization_experiment, only: leaf
         real(real64), dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat
         procedure(leaf) :: leaf_routine
      end subroutine root
   end interface

   type :: leaf_wrapper
      procedure(leaf), pointer, nopass :: p => null()
   end type
   type :: root_wrapper
      procedure(root), pointer, nopass :: p => null()
   end type

   type(leaf_wrapper) :: leaf_list(3)
   type(root_wrapper) :: root_list(4)
   character(len=10), parameter :: leaf_names(3) = ['dc_leaf   ', 'osimd_leaf', "bare_leaf "]
   character(len=10), parameter :: root_names(4) = ['dc_root   ', 'omp_root  ', 'ocon_root ', "osimd_root"]

   ! Allocate and initialize matrices
   allocate(x(n, m), y(n, m), out_mat(n, m))
   x = 1.0_real64
   y = 0.5_real64

   call system_clock(count_rate=rate)

   ! List of leaf routines
   leaf_list(1)%p => dc_leaf
   leaf_list(2)%p => osimd_leaf
   leaf_list(3)%p => bare_leaf


   ! List of root routines
   root_list(1)%p => dc_root
   root_list(2)%p => omp_root
   root_list(3)%p => ocon_root
   root_list(4)%p => osimd_root

   print "(A)", "Root,Leaf,Time(sec),Correct_Result"

   time = 0.0_real64
   correct_result = .true.
   do k = 1, n_times_per_exp
      call system_clock(start)
      call dc_collapsed(x, y, out_mat)
      call system_clock(finish)

      correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

      time = (time + real(finish - start) / real(rate)) / 2

   end do
   print '(A, ",", A, ",", G0, ",", A)', "dc_collapsed", "none", time, merge("T", "F", correct_result)

   time = 0.0_real64
   correct_result = .true.
   do k = 1, n_times_per_exp
      call system_clock(start)
      call omp_collapsed(x, y, out_mat)
      call system_clock(finish)

      correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

      time = (time + real(finish - start) / real(rate)) / 2

   end do
   print '(A, ",", A, ",", G0, ",", A)', "omp_collapsed", "none", time, merge("T", "F", correct_result)

   time = 0.0_real64
   correct_result = .true.
   do k = 1, n_times_per_exp
      call system_clock(start)
      call elemental_calc(x, y, out_mat)
      call system_clock(finish)

      correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

      time = (time + real(finish - start) / real(rate)) / 2

   end do
   print '(A, ",", A, ",", G0, ",", A)', "elemental_calc", "none", time, merge("T", "F", correct_result)

   ! Loop over all combinations
   do i = 1, size(root_list)
      do j = 1, size(leaf_list)
         time = 0.0_real64
         correct_result = .true.
         do k = 1, n_times_per_exp
            call system_clock(start)
            call root_list(i)%p(x, y, out_mat, leaf_list(j)%p)
            call system_clock(finish)

            correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

            time = (time + real(finish - start) / real(rate)) / 2

         end do
         print '(A, ",", A, ",", G0, ",", A)', trim(root_names(i)), trim(leaf_names(j)), time, merge("T", "F", correct_result)
      end do
   end do

   time = 0.0_real64
   correct_result = .true.
   do k = 1, n_times_per_exp
      call system_clock(start)
      call ocon_bare_inline(x, y, out_mat)
      call system_clock(finish)

      correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

      time = (time + real(finish - start) / real(rate)) / 2

   end do
   print '(A, ",", A, ",", G0, ",", A)', "ocon_bare_inline", "none", time, merge("T", "F", correct_result)

   time = 0.0_real64
   correct_result = .true.
   do k = 1, n_times_per_exp
      call system_clock(start)
      call bare(x, y, out_mat)
      call system_clock(finish)

      correct_result = correct_result .and. all(out_mat == sqrt(0.5_real64))

      time = (time + real(finish - start) / real(rate)) / 2

   end do
   print '(A, ",", A, ",", G0, ",", A)', "bare", "none", time, merge("T", "F", correct_result)
end program benchmark_parallelization
