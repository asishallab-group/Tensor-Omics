module parallelization_experiment
   use iso_fortran_env, only: int32, real64
   use config, only: alignment
   implicit none

   interface
      subroutine experiment(x_mat, y_mat, out_mat)
         use iso_fortran_env, only: real64
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat
      end subroutine experiment
   end interface

contains

   elemental subroutine elemental_calc(x, y, out_val)
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      real(real64), intent(out) :: out_val

      out_val = log(1 + exp(x * y))
   end subroutine elemental_calc

   subroutine dc_collapsed(x_mat, y_mat, out_mat)
      real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

      integer(int32) :: i, j

      do concurrent (i = 1:size(x_mat, 2), j = 1:size(x_mat, 1)) shared(x_mat, y_mat, out_mat)
         out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
      end do
   end subroutine dc_collapsed

   subroutine bare(x_mat, y_mat, out_mat)
      real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

      integer(int32) :: i, j

      do j = 1, size(x_mat, 1)
         do i = 1, size(x_mat, 2)
            out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
         end do
      end do
   end subroutine bare

   subroutine omp_collapsed(x_mat, y_mat, out_mat)
      real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
      real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

      integer(int32) :: i, j

      !$omp parallel do collapse(2) default(none) shared(x_mat, y_mat, out_mat) private(i, j) schedule(static)
      do j = 1, size(x_mat, 1)
         do i = 1, size(x_mat, 2)
            out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
         end do
      end do
      !$omp end parallel do
   end subroutine omp_collapsed

!!!! Leaf: do concurrent
      subroutine dc_dc(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         do concurrent (i = 1:size(x_mat, 2)) shared(x_mat, y_mat)
            do concurrent (j = 1:size(x_mat, 1)) shared(x_mat, y_mat, out_mat)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
         end do
      end subroutine dc_dc

      subroutine omp_dc(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         !$omp parallel do shared(x_mat, y_mat, out_mat) private(i, j) schedule(static)
         do i = 1, size(x_mat, 2)
            do concurrent (j = 1:size(x_mat, 1)) shared(x_mat, y_mat, out_mat)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
         end do
         !$omp end parallel do
      end subroutine omp_dc

      subroutine ocon_dc(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         !$omp parallel private(i, j) shared(x_mat, y_mat, out_mat)
         do concurrent (i = 1:size(x_mat, 2))
            do concurrent (j = 1:size(x_mat, 1)) shared(x_mat, y_mat, out_mat)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
         end do
         !$omp end parallel
      end subroutine ocon_dc

      subroutine osimd_dc(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         !$omp simd
         do i = 1, size(x_mat, 2)
            do concurrent (j = 1:size(x_mat, 1)) shared(x_mat, y_mat, out_mat)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
         end do
         !$omp end simd
      end subroutine osimd_dc

!!!! Leaf: osimd
      subroutine dc_osimd(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         do concurrent (i = 1:size(x_mat, 2)) shared(x_mat, y_mat)
            !$omp simd
            do j = 1, size(x_mat, 1)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
            !$omp end simd
         end do
      end subroutine dc_osimd

      subroutine omp_osimd(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         !$omp parallel do default(none) shared(x_mat, y_mat, out_mat) private(i, j) schedule(static)
         do i = 1, size(x_mat, 2)
            !$omp simd
            do j = 1, size(x_mat, 1)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
            !$omp end simd
         end do
         !$omp end parallel do
      end subroutine omp_osimd

      subroutine ocon_osimd(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         !$omp parallel private(i, j) shared(x_mat, y_mat, out_mat)
         do concurrent (i = 1:size(x_mat, 2))
            !$omp simd
            do j = 1, size(x_mat, 1)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
            !$omp end simd
         end do
         !$omp end parallel
      end subroutine ocon_osimd

      subroutine osimd_osimd(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         !$omp simd
         do i = 1, size(x_mat, 2)
            !$omp simd
            do j = 1, size(x_mat, 1)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
            !$omp end simd
         end do
         !$omp end simd
      end subroutine osimd_osimd

!!!! Leaf: bare do-loop
      subroutine dc_bare(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         do concurrent (i = 1:size(x_mat, 2)) shared(x_mat, y_mat)
            do j = 1, size(x_mat, 1)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
         end do
      end subroutine dc_bare

      subroutine omp_bare(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         !$omp parallel do default(none) shared(x_mat, y_mat, out_mat) private(i, j) schedule(static)
         do i = 1, size(x_mat, 2)
            do j = 1, size(x_mat, 1)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
         end do
         !$omp end parallel do
      end subroutine omp_bare

      subroutine ocon_bare(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         !$omp parallel private(i, j) shared(x_mat, y_mat, out_mat)
         do concurrent (i = 1:size(x_mat, 2))
            do j = 1, size(x_mat, 1)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
         end do
         !$omp end parallel
      end subroutine ocon_bare

      subroutine osimd_bare(x_mat, y_mat, out_mat)
         real(real64), contiguous, dimension(:, :), intent(in) :: x_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(in) :: y_mat
         real(real64), dimension(size(x_mat, 1), size(x_mat, 2)), intent(out) :: out_mat

         integer(int32) :: i, j

         !$omp simd
         do i = 1, size(x_mat, 2)
            do j = 1, size(x_mat, 1)
               out_mat(j, i) = log(1 + exp(x_mat(j, i) * y_mat(j, i)))
            end do
         end do
         !$omp end simd
      end subroutine osimd_bare

end module parallelization_experiment