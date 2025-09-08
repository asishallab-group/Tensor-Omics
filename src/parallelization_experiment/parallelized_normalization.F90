!> Module with normalization routines for tensor omics.
module parallelized_normalization
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
contains

  pure subroutine normalize_by_std_dev_bare(n_genes, n_tissues, input_matrix, output_matrix)
    integer(int32), intent(in) :: n_genes
    integer(int32), intent(in) :: n_tissues
    real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
    real(real64), intent(out) :: output_matrix(n_genes * n_tissues)

    integer(int32) :: i, j
    real(real64) :: std_dev, temp_sum

    do i = 1, n_genes
        temp_sum = 0.0d0
        do j = 1, n_tissues
            temp_sum = temp_sum + input_matrix((j-1)*n_genes + i)**2
        end do

        std_dev = sqrt(temp_sum / dble(n_tissues))
        if (std_dev == 0.0d0) std_dev = 1.0d0

        do j = 1, n_tissues
            output_matrix((j-1)*n_genes + i) = input_matrix((j-1)*n_genes + i) / std_dev
        end do
    end do      
  end subroutine normalize_by_std_dev_bare

  pure subroutine log2_transformation_bare(n_genes, n_tissues, input_matrix, output_matrix)
    implicit none

    !| Number of genes (rows)
    integer(int32), intent(in) :: n_genes
    !| Number of tissues (columns)
    integer(int32), intent(in) :: n_tissues
    !| Flattened input matrix (size: n_genes * n_tissues)
    real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
    !| Output matrix (same size as input)
    real(real64), intent(out) :: output_matrix(n_genes * n_tissues)

    ! Locals
    integer(int32) :: i
    real(real64), parameter :: LOG2 = log(2.0d0)

    ! Loop through all elements in the flattened input matrix
    do i = 1, n_genes * n_tissues
        ! Apply the log2(x + 1) transformation
        output_matrix(i) = log(input_matrix(i) + 1.0d0) / LOG2
    end do
  end subroutine log2_transformation_bare

  pure subroutine normalize_by_std_dev_dc(n_genes, n_tissues, input_matrix, output_matrix)
    integer(int32), intent(in) :: n_genes
    integer(int32), intent(in) :: n_tissues
    real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
    real(real64), intent(out) :: output_matrix(n_genes * n_tissues)

    integer(int32) :: i, j
    real(real64) :: std_dev, temp_sum

    do concurrent (i = 1:n_genes) shared(input_matrix, output_matrix) local(temp_sum, std_dev)
        temp_sum = 0.0d0
        do concurrent (j = 1:n_tissues) shared(input_matrix, i, n_genes) reduce(+:temp_sum)
            temp_sum = temp_sum + input_matrix((j-1)*n_genes + i)**2
        end do

        std_dev = sqrt(temp_sum / dble(n_tissues))
        if (std_dev == 0.0d0) std_dev = 1.0d0

        do concurrent (j = 1:n_tissues)
            output_matrix((j-1)*n_genes + i) = input_matrix((j-1)*n_genes + i) / std_dev
        end do
    end do      
  end subroutine normalize_by_std_dev_dc

  pure subroutine log2_transformation_dc(n_genes, n_tissues, input_matrix, output_matrix)
    implicit none

    !| Number of genes (rows)
    integer(int32), intent(in) :: n_genes
    !| Number of tissues (columns)
    integer(int32), intent(in) :: n_tissues
    !| Flattened input matrix (size: n_genes * n_tissues)
    real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
    !| Output matrix (same size as input)
    real(real64), intent(out) :: output_matrix(n_genes * n_tissues)

    ! Locals
    integer(int32) :: i
    real(real64), parameter :: LOG2 = log(2.0d0)

    ! Loop through all elements in the flattened input matrix
    do concurrent (i = 1:n_genes * n_tissues) shared(output_matrix, input_matrix)
        ! Apply the log2(x + 1) transformation
        output_matrix(i) = log(input_matrix(i) + 1.0d0) / LOG2
    end do
  end subroutine log2_transformation_dc

  subroutine normalize_by_std_dev_osimd(n_genes, n_tissues, input_matrix, output_matrix)
    integer(int32), intent(in) :: n_genes
    integer(int32), intent(in) :: n_tissues
    real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
    real(real64), intent(out) :: output_matrix(n_genes * n_tissues)

    integer(int32) :: i, j
    real(real64) :: std_dev, temp_sum

    !$omp simd
    do i = 1, n_genes
        temp_sum = 0.0d0
        !$omp simd
        do j = 1, n_tissues
            temp_sum = temp_sum + input_matrix((j-1)*n_genes + i)**2
        end do
        !$omp end simd

        std_dev = sqrt(temp_sum / dble(n_tissues))
        if (std_dev == 0.0d0) std_dev = 1.0d0

        !$omp simd
        do j = 1, n_tissues
            output_matrix((j-1)*n_genes + i) = input_matrix((j-1)*n_genes + i) / std_dev
        end do
        !$omp end simd
    end do
    !$omp end simd
  end subroutine normalize_by_std_dev_osimd

  subroutine log2_transformation_osimd(n_genes, n_tissues, input_matrix, output_matrix)
    implicit none

    !| Number of genes (rows)
    integer(int32), intent(in) :: n_genes
    !| Number of tissues (columns)
    integer(int32), intent(in) :: n_tissues
    !| Flattened input matrix (size: n_genes * n_tissues)
    real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
    !| Output matrix (same size as input)
    real(real64), intent(out) :: output_matrix(n_genes * n_tissues)

    ! Locals
    integer(int32) :: i
    real(real64), parameter :: LOG2 = log(2.0d0)

    ! Loop through all elements in the flattened input matrix
    !$omp simd
    do i = 1, n_genes * n_tissues
        ! Apply the log2(x + 1) transformation
        output_matrix(i) = log(input_matrix(i) + 1.0d0) / LOG2
    end do
    !$omp end simd
  end subroutine log2_transformation_osimd
end module parallelized_normalization
