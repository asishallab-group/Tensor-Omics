module parallelization_experiment
   use iso_fortran_env, only: int32, real64
   use config, only: alignment
   use parallelized_normalization
   use tox_errors, only: set_ok, is_ok
   implicit none

contains
   subroutine bare(input_matrix, std_dev_out, log_out)
      real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
      real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


      integer(int32) :: i

      do i = 1, size(input_matrix, 2)
         call normalize_by_std_dev_bare(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
         call log2_transformation_bare(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
      end do
   end subroutine bare

!!!! Leaf: do concurrent
      subroutine dc_dc(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         do concurrent (i = 1:size(input_matrix, 2)) shared(input_matrix, std_dev_out)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
      end subroutine dc_dc

      subroutine omp_dc(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         !$omp parallel do shared(input_matrix, std_dev_out, log_out) private(i)
         do i = 1, size(input_matrix, 2)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
         !$omp end parallel do
      end subroutine omp_dc

      subroutine ocon_dc(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         !$omp parallel
         do concurrent (i = 1:size(input_matrix, 2)) shared(input_matrix, std_dev_out)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
         !$omp end parallel
      end subroutine ocon_dc

      subroutine osimd_dc(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         !$omp simd
         do i = 1, size(input_matrix, 2)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
         !$omp end simd
      end subroutine osimd_dc

!!!! Leaf: osimd
      subroutine dc_osimd(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         do concurrent (i = 1:size(input_matrix, 2)) shared(input_matrix, std_dev_out)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
      end subroutine dc_osimd

      subroutine omp_osimd(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         !$omp parallel do shared(input_matrix, std_dev_out, log_out) private(i)
         do i = 1, size(input_matrix, 2)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
         !$omp end parallel do
      end subroutine omp_osimd

      subroutine ocon_osimd(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         !$omp parallel
         do concurrent (i = 1:size(input_matrix, 2)) shared(input_matrix, std_dev_out)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
         !$omp end parallel
      end subroutine ocon_osimd

      subroutine osimd_osimd(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         !$omp simd
         do i = 1, size(input_matrix, 2)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
         !$omp end simd
      end subroutine osimd_osimd

!!!! Leaf: bare do-loop
      subroutine dc_bare(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         do concurrent (i = 1:size(input_matrix, 2)) shared(input_matrix, std_dev_out)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
      end subroutine dc_bare

      subroutine omp_bare(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         !$omp parallel do shared(input_matrix, std_dev_out, log_out) private(i)
         do i = 1, size(input_matrix, 2)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
         !$omp end parallel do
      end subroutine omp_bare

      subroutine ocon_bare(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         !$omp parallel
         do concurrent (i = 1:size(input_matrix, 2)) shared(input_matrix, std_dev_out)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
         !$omp end parallel
      end subroutine ocon_bare

      subroutine osimd_bare(input_matrix, std_dev_out, log_out)
         real(real64), contiguous, dimension(:, :), intent(in) :: input_matrix
         real(real64), dimension(size(input_matrix, 1), size(input_matrix, 2)), intent(out) :: std_dev_out, log_out


         integer(int32) :: i

         !$omp simd
         do i = 1, size(input_matrix, 2)
            call normalize_by_std_dev_dc(1, size(input_matrix, 1), input_matrix(:, i), std_dev_out(:, i))
            call log2_transformation_dc(1, size(input_matrix, 1), std_dev_out(:, i), log_out(:, i))
         end do
         !$omp end simd
      end subroutine osimd_bare

end module parallelization_experiment