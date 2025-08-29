module parallelization_experiment
   use iso_fortran_env, only: int32, real64
   implicit none

   interface
      pure subroutine leaf(x_arr, y_arr, out_arr)
         use iso_fortran_env, only: real64
         real(real64), dimension(:), intent(in) :: x_arr
         real(real64), dimension(size(x_arr, 1)), intent(in) :: y_arr
         real(real64), dimension(size(x_arr, 1)), intent(out) :: out_arr
      end subroutine
   end interface

contains

   elemental subroutine elemental_calc(x, y, out_val)
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      real(real64), intent(out) :: out_val

      out_val = sqrt(x - y)
   end subroutine elemental_calc

   pure subroutine dc_leaf(x_arr, y_arr, out_arr)
      real(real64), dimension(:), intent(in) :: x_arr
      real(real64), dimension(size(x_arr, 1)), intent(in) :: y_arr
      real(real64), dimension(size(x_arr, 1)), intent(out) :: out_arr

      integer(int32) :: i

      do concurrent (i = 1:size(x_arr, 1)) shared(x_arr, y_arr, out_arr)
         out_arr(i) = sqrt(x_arr(i) - y_arr(i))
      end do
   end subroutine dc_leaf

   subroutine dc_root(x_mat, y_mat, out_mat, leaf_routine)
      real(real64), dimension(:, :), intent(in) :: x_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat
      procedure(leaf) :: leaf_routine

      integer(int32) :: i

      do concurrent (i = 1:size(x_mat, 1)) shared(x_mat, y_mat)
         call leaf_routine(x_mat(:, i), y_mat(:, i), out_mat(:, i))
      end do
   end subroutine dc_root

   pure subroutine dc_collapsed(x_mat, y_mat, out_mat)
      real(real64), dimension(:, :), intent(in) :: x_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat
      procedure(leaf) :: leaf_routine

      integer(int32) :: i, j

      do concurrent (i = 1:size(x_mat, 1), j = 1:size(x_mat, 2)) shared(x_mat, y_mat, out_mat)
         out_mat(j, i) = sqrt(x_mat(j, i) - y_mat(j, i))
      end do
   end subroutine dc_collapsed

   pure subroutine bare(x_mat, y_mat, out_mat)
      real(real64), dimension(:, :), intent(in) :: x_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat
      procedure(leaf) :: leaf_routine

      integer(int32) :: i, j

      do j = 1, size(x_mat, 2)
         do i = 1, size(x_mat, 1)
            out_mat(j, i) = sqrt(x_mat(j, i) - y_mat(j, i))
         end do
      end do
   end subroutine bare

   subroutine omp_collapsed(x_mat, y_mat, out_mat)
      use iso_fortran_env, only: int32, real64
      use omp_lib
      implicit none

      real(real64), dimension(:, :), intent(in) :: x_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

      integer(int32) :: i, j

      !$omp parallel do collapse(2) default(none) shared(x_mat, y_mat, out_mat) private(i, j) schedule(static)
      do j = 1, size(x_mat, 2)
         do i = 1, size(x_mat, 1)
            out_mat(i, j) = sqrt(x_mat(i, j) - y_mat(i, j))
         end do
      end do
      !$omp end parallel do
   end subroutine omp_collapsed

   subroutine omp_root(x_mat, y_mat, out_mat, leaf_routine)
      use omp_lib
      real(real64), dimension(:, :), intent(in) :: x_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat
      procedure(leaf) :: leaf_routine

      integer(int32) :: i

      !$omp parallel do default(none) shared(x_mat, y_mat, out_mat) private(i) schedule(static)
      do i = 1, size(x_mat, 2)
         call leaf_routine(x_mat(:, i), y_mat(:, i), out_mat(:, i))
      end do
      !$omp end parallel do
   end subroutine omp_root

   pure subroutine osimd_leaf(x_arr, y_arr, out_arr)
      real(real64), dimension(:), intent(in) :: x_arr
      real(real64), dimension(size(x_arr, 1)), intent(in) :: y_arr
      real(real64), dimension(size(x_arr, 1)), intent(out) :: out_arr

      integer(int32) :: i

      !$omp simd
      do i = 1, size(x_arr, 1)
         out_arr(i) = sqrt(x_arr(i) - y_arr(i))
      end do
      !$omp end simd
   end subroutine osimd_leaf

   pure subroutine bare_leaf(x_arr, y_arr, out_arr)
      real(real64), dimension(:), intent(in) :: x_arr
      real(real64), dimension(size(x_arr, 1)), intent(in) :: y_arr
      real(real64), dimension(size(x_arr, 1)), intent(out) :: out_arr

      integer(int32) :: i

      do i = 1, size(x_arr, 1)
         out_arr(i) = sqrt(x_arr(i) - y_arr(i))
      end do
   end subroutine bare_leaf

   subroutine ocon_root(x_mat, y_mat, out_mat, leaf_routine)
      real(real64), dimension(:, :), intent(in) :: x_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat
      procedure(leaf) :: leaf_routine

      integer(int32) :: i

      !$omp parallel private(i) shared(x_mat, y_mat, out_mat)
      do concurrent (i = 1:size(x_mat, 2))
         call leaf_routine(x_mat(:, i), y_mat(:, i), out_mat(:, i))
      end do
      !$omp end parallel
   end subroutine ocon_root

   subroutine osimd_root(x_mat, y_mat, out_mat, leaf_routine)
      real(real64), dimension(:, :), intent(in) :: x_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat
      procedure(leaf) :: leaf_routine

      integer(int32) :: i

      !$omp simd
      do i = 1, size(x_mat, 2)
         call leaf_routine(x_mat(:, i), y_mat(:, i), out_mat(:, i))
      end do
      !$omp end simd
   end subroutine osimd_root

   subroutine ocon_bare_inline(x_mat, y_mat, out_mat)
      real(real64), dimension(:, :), intent(in) :: x_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

      integer(int32) :: i, j

      !$omp parallel private(i) shared(x_mat, y_mat, out_mat)
      do concurrent (i = 1:size(x_mat, 2))
         do j = 1, size(x_mat, 1)
            out_mat(j, i) = sqrt(x_mat(j, i) - y_mat(j, i))
         end do
      end do
      !$omp end parallel
   end subroutine ocon_bare_inline

end module parallelization_experiment