!> Module with normalization routines for tensor omics.
module tox_normalization
  use, intrinsic :: iso_fortran_env, only: real64, int32
contains

  !> Normalizes each gene's expression vector using `sqrt(mean(x^2))`
  !| across tissues (not classical standard deviation).
  pure subroutine normalize_by_std_dev(n_genes, n_tissues, input_matrix, output_matrix)
      implicit none

      !| Number of genes (rows)
      integer(int32), intent(in) :: n_genes
      !| Number of tissues (columns)
      integer(int32), intent(in) :: n_tissues
      !| Flattened input matrix of gene expression values (column-major)
      real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
      !| Output normalized matrix (same shape as input)
      real(real64), intent(out) :: output_matrix(n_genes * n_tissues)

      ! Local variables
      integer(int32) :: i, j
      real(real64) :: std_dev, temp_sum

      ! Loop over each gene
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

      
  end subroutine normalize_by_std_dev

  !> Quantile normalization of a gene expression matrix (F42-compliant).
  !| Computes average expression per rank across tissues.

  pure subroutine quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix, &
                                      temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    use f42_utils, only: sort_array

    implicit none

    !| Number of genes (rows)
    integer(int32), intent(in) :: n_genes
    !| Number of tissues (columns)
    integer(int32), intent(in) :: n_tissues
    !| Flattened input matrix (column-major)
    real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
    !| Output normalized matrix (same shape as input)
    real(real64), intent(out) :: output_matrix(n_genes * n_tissues)
    !| Temporary vector for column sorting (size n_genes)
    real(real64), intent(inout) :: temp_col(n_genes)
    !| Preallocated vector to store rank means (size n_genes)
    real(real64), intent(inout) :: rank_means(n_genes)
    !| Permutation vector (size n_genes)
    integer(int32), intent(inout) :: perm(n_genes)
    !| Manual quicksort stack (≥ log2(n_genes) + 10)
    integer(int32), intent(inout) :: stack_left(max_stack)
    !| Manual quicksort stack (same size as stack_left)
    integer(int32), intent(inout) :: stack_right(max_stack)
    !| Stack size passed from R
    integer(int32), intent(in) :: max_stack

      ! Locals
      integer(int32) :: i, j

      ! Initialize rank means
      rank_means = 0.0d0

      ! === First pass: accumulate values by rank across tissues ===
      do j = 1, n_tissues
          ! Prepare current column and initialize permutation
          do i = 1, n_genes
              temp_col(i) = input_matrix((j - 1) * n_genes + i)
              perm(i) = i
          end do

          ! Sort current column with index tracking

          call sort_array(temp_col, perm, stack_left, stack_right)

          ! Accumulate values for each rank
          do i = 1, n_genes
              rank_means(i) = rank_means(i) + temp_col(perm(i))
          end do
      end do

      ! Average the rank values
      do i = 1, n_genes
          rank_means(i) = rank_means(i) / dble(n_tissues)
      end do

      ! === Second pass: assign averaged values by rank ===
      do j = 1, n_tissues
          ! Prepare column and reset permutation
          do i = 1, n_genes
              temp_col(i) = input_matrix((j - 1) * n_genes + i)
              perm(i) = i
          end do

          call sort_array(temp_col, perm, stack_left, stack_right)


          do i = 1, n_genes
              output_matrix((j - 1) * n_genes + perm(i)) = rank_means(i)
          end do
      end do
  end subroutine quantile_normalization

  !> Apply `log2(x + 1)` transformation to each element of the input matrix.
  !| This subroutine performs element-wise `log2(x + 1)` transformation on a
  !| matrix flattened in column-major order. The `log2` is computed via:
  !| `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
  !| non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).

  pure subroutine log2_transformation(n_genes, n_tissues, input_matrix, output_matrix)
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
  end subroutine log2_transformation


  !> Calculate tissue averages by averaging replicates within each group.
  !| For each group of tissue replicates, this subroutine computes the average
  !| expression per gene. The input matrix is column-major, flattened as a 1D array.

  pure subroutine calc_tiss_avg(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix)
      implicit none

      !| Number of genes (rows)
      integer(int32), intent(in) :: n_gene
      !| Number of tissue groups
      integer(int32), intent(in) :: n_grps
      !| Start column index for each group (length: n_grps)
      integer(int32), intent(in) :: group_s(n_grps)
      !| Number of columns per group (length: n_grps)
      integer(int32), intent(in) :: group_c(n_grps)
      !| Flattened input matrix (length: n_gene * n_col)
      real(real64), intent(in) :: input_matrix(n_gene * sum(group_c))
      !| Flattened output matrix (length: n_gene * n_grps)
      real(real64), intent(out) :: output_matrix(n_gene * n_grps)

      ! === Local variables ===
      integer(int32) :: i, j, g, col
      real(real64) :: sum_val
      integer(int32) :: start_idx, count_cols

      ! === Loop over each group ===
      do g = 1, n_grps
          start_idx = group_s(g)
          count_cols = group_c(g)

          do i = 1, n_gene
              sum_val = 0.0d0

              do j = 0, count_cols - 1
                  col = start_idx + j
                  sum_val = sum_val + input_matrix((col - 1) * n_gene + i)
              end do

              output_matrix((g - 1) * n_gene + i) = sum_val / dble(count_cols)
          end do
      end do
  end subroutine calc_tiss_avg


  !> Calculate `log2 fold changes` between condition and control columns.
  !| For each control-condition pair, this subroutine computes the `log2 fold change`
  !| by subtracting the expression value in the control column from the corresponding
  !| value in the condition column, for all genes.
  !|
  !| The input matrix must be column-major and flattened as a 1D array.

  pure subroutine calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix)

      implicit none

      ! === Arguments ===
      !| Number of genes (rows)
      integer(int32), intent(in) :: n_genes
      !| Number of columns in the input matrix
      integer(int32), intent(in) :: n_cols
      !| Number of control-condition pairs
      integer(int32), intent(in) :: n_pairs
      !| Control column indices (length n_pairs)
      integer(int32), intent(in) :: control_cols(n_pairs)
      !| Condition column indices (length n_pairs)
      integer(int32), intent(in) :: cond_cols(n_pairs)
      !| Input matrix, flattened (length: n_genes × n_cols)
      real(real64), intent(in) :: i_matrix(n_genes * n_cols)
      !| Output matrix for fold changes (length: n_genes × n_pairs)
      real(real64), intent(out) :: o_matrix(n_genes * n_pairs)

      ! === Locals ===
      integer(int32) :: i, p
      integer(int32) :: control_col, cond_col

      ! === Loop over each pair ===
      do p = 1, n_pairs
          control_col = control_cols(p)
          cond_col    = cond_cols(p)

          do i = 1, n_genes
              o_matrix((p - 1) * n_genes + i) = &
                  i_matrix((cond_col - 1) * n_genes + i) - &
                  i_matrix((control_col - 1) * n_genes + i)
          end do
      end do

  end subroutine calc_fchange


end module tox_normalization


!> R/Fortran wrapper for normalization by standard deviation.
!| Provides an interface for R (.Fortran) and Fortran code to call the normalization routine.
subroutine normalize_by_std_dev_r(n_genes, n_tissues, input_matrix, output_matrix)
  use tox_normalization
  !| Number of genes (rows)
  integer(int32), intent(in) :: n_genes
  !| Number of tissues (columns)
  integer(int32), intent(in) :: n_tissues
  !| Input matrix (n_genes x n_tissues)
  real(real64), intent(in)  :: input_matrix(n_genes, n_tissues)
  !| Output normalized matrix (same shape as input)
  real(real64), intent(out) :: output_matrix(n_genes, n_tissues)
  call normalize_by_std_dev(n_genes, n_tissues, input_matrix, output_matrix)
end subroutine normalize_by_std_dev_r

!> C/Python wrapper for normalization by standard deviation.
!| Provides a C/Python-compatible interface to the normalization routine.
subroutine normalize_by_std_dev_c(n_genes, n_tissues, input_matrix, output_matrix) bind(C, name="normalize_by_std_dev_c")
  use iso_c_binding, only : c_int, c_double, c_f_pointer, c_loc
  use tox_normalization
  !| Number of genes (rows)
  integer(c_int), value :: n_genes
  !| Number of tissues (columns)
  integer(c_int), value :: n_tissues
  !| Input matrix (flattened, n_genes * n_tissues)
  real(c_double), intent(in), target :: input_matrix(n_genes * n_tissues)
  !| Output normalized matrix (flattened, same shape as input)
  real(c_double), intent(out), target :: output_matrix(n_genes * n_tissues)
  real(c_double), pointer :: inmat(:,:), outmat(:,:)
  call c_f_pointer(c_loc(input_matrix(1)), inmat, [n_genes, n_tissues])
  call c_f_pointer(c_loc(output_matrix(1)), outmat, [n_genes, n_tissues])
  call normalize_by_std_dev(n_genes, n_tissues, inmat, outmat)
end subroutine normalize_by_std_dev_c


!> R/Fortran wrapper for quantile normalization.
!| Provides an interface for R (.Fortran) and Fortran code to call the quantile normalization routine.
subroutine quantile_normalization_r(n_genes, n_tissues, input_matrix, output_matrix, &
                                        temp_col, rank_means, perm, stack_left, stack_right, max_stack)
  use tox_normalization
  !| Number of genes (rows)
  integer(int32), intent(in) :: n_genes
  !| Number of tissues (columns)
  integer(int32), intent(in) :: n_tissues
  !| Stack size for sorting
  integer(int32), intent(in) :: max_stack
  !| Input matrix (n_genes x n_tissues)
  real(real64), intent(in)  :: input_matrix(n_genes, n_tissues)
  !| Output normalized matrix (same shape as input)
  real(real64), intent(out) :: output_matrix(n_genes, n_tissues)
  !| Temporary vector for column sorting (size n_genes)
  real(real64), intent(inout) :: temp_col(n_genes)
  !| Preallocated vector to store rank means (size n_genes)
  real(real64), intent(inout) :: rank_means(n_genes)
  !| Permutation vector (size n_genes)
  integer(int32), intent(inout) :: perm(n_genes)
  !| Manual quicksort stack (size max_stack)
  integer(int32), intent(inout) :: stack_left(max_stack)
  !| Manual quicksort stack (size max_stack)
  integer(int32), intent(inout) :: stack_right(max_stack)

  call quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix, &
                            temp_col, rank_means, perm, stack_left, stack_right, max_stack)
end subroutine quantile_normalization_r

!> C/Python wrapper for quantile normalization.
!| Provides a C/Python-compatible interface to the quantile normalization routine.
subroutine quantile_normalization_c(n_genes, n_tissues, input_matrix, output_matrix, &
                                    temp_col, rank_means, perm, stack_left, stack_right, max_stack) &
                                    bind(C, name="quantile_normalization_c")
  use iso_c_binding, only : c_int, c_double
  use tox_normalization
  !| Number of genes (rows)
  integer(c_int), intent(in), value :: n_genes
  !| Number of tissues (columns)
  integer(c_int), intent(in), value :: n_tissues
  !| Stack size for sorting
  integer(c_int), intent(in), value :: max_stack
  !| Input matrix (n_genes x n_tissues)
  real(c_double), intent(in), target :: input_matrix(n_genes, n_tissues)
  !| Output normalized matrix (same shape as input)
  real(c_double), intent(out), target :: output_matrix(n_genes, n_tissues)
  !| Temporary vector for column sorting (size n_genes)
  real(c_double), intent(inout), target :: temp_col(n_genes)
  !| Preallocated vector to store rank means (size n_genes)
  real(c_double), intent(inout), target :: rank_means(n_genes)
  !| Permutation vector (size n_genes)
  integer(c_int), intent(inout), target :: perm(n_genes)
  !| Manual quicksort stack (size max_stack)
  integer(c_int), intent(inout), target :: stack_left(max_stack)
  !| Manual quicksort stack (size max_stack)
  integer(c_int), intent(inout), target :: stack_right(max_stack)

  call quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix, &
                            temp_col, rank_means, perm, stack_left, stack_right, max_stack)
end subroutine quantile_normalization_c

!> R/Fortran wrapper for log2 transformation.
!| Provides an interface for R (.Fortran) and Fortran code to call the log2 transformation routine.
!| Applies log2(x+1) to each element of the input matrix. Arguments match R's .Fortran calling convention and expect flat arrays.
subroutine log2_transformation_r(n_genes, n_tissues, input_matrix, output_matrix)
  use tox_normalization
  !| Number of genes (rows)
  integer(int32), intent(in) :: n_genes
  !| Number of tissues (columns)
  integer(int32), intent(in) :: n_tissues
  !| Input matrix (flattened, n_genes * n_tissues)
  real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
  !| Output matrix (flattened, same shape as input)
  real(real64), intent(out) :: output_matrix(n_genes * n_tissues)
  call log2_transformation(n_genes, n_tissues, input_matrix, output_matrix)
end subroutine log2_transformation_r

!> C/Python wrapper for log2 transformation.
!| Provides a C/Python-compatible interface to the log2 transformation routine.
!| Expects flat arrays, matching C calling conventions. Suitable for use with ctypes.
!| Applies log2(x+1) to each element of the input matrix.
subroutine log2_transformation_c(n_genes, n_tissues, input_matrix, output_matrix) bind(C, name="log2_transformation_c")
  use iso_c_binding, only : c_int, c_double
  use tox_normalization
  !| Number of genes (rows)
  integer(c_int), intent(in), value :: n_genes
  !| Number of tissues (columns)
  integer(c_int), intent(in), value :: n_tissues
  !| Input matrix (flattened, n_genes * n_tissues)
  real(c_double), intent(in), target :: input_matrix(n_genes * n_tissues)
  !| Output matrix (flattened, same shape as input)
  real(c_double), intent(out), target :: output_matrix(n_genes * n_tissues)

  call log2_transformation(n_genes, n_tissues, input_matrix, output_matrix)
end subroutine log2_transformation_c

!> R/Fortran wrapper for tissue average calculation.
!| Provides an interface for R (.Fortran) and Fortran code to call the tissue average calculation routine.
!| Computes average expression per gene for each group of tissue replicates. Arguments match R's .Fortran calling convention.
subroutine calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix)
  use tox_normalization
  !| Number of genes (rows)
  integer(int32), intent(in) :: n_gene
  !| Number of tissue groups
  integer(int32), intent(in) :: n_grps
  !| Start column index for each group (length: n_grps)
  integer(int32), intent(in) :: group_s(n_grps)
  !| Number of columns per group (length: n_grps)
  integer(int32), intent(in) :: group_c(n_grps)
  !| Input matrix (flattened, n_gene * sum(group_c))
  real(real64), intent(in) :: input_matrix(n_gene * sum(group_c))
  !| Output matrix (flattened, n_gene * n_grps)
  real(real64), intent(out) :: output_matrix(n_gene * n_grps)
  call calc_tiss_avg(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix)
end subroutine calc_tiss_avg_r

!> C/Python wrapper for tissue average calculation.
!| Provides a C/Python-compatible interface to the tissue average calculation routine.
!| Suitable for use with ctypes. Computes average expression per gene for each group of tissue replicates.
subroutine calc_tiss_avg_c(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix) bind(C, name="calc_tiss_avg_c")
  use iso_c_binding, only : c_int, c_double
  use tox_normalization
  !| Number of genes (rows)
  integer(c_int), intent(in), value :: n_gene
  !| Number of tissue groups
  integer(c_int), intent(in), value :: n_grps
  !| Start column index for each group (length: n_grps)
  integer(c_int), intent(in), target :: group_s(n_grps)
  !| Number of columns per group (length: n_grps)
  integer(c_int), intent(in), target :: group_c(n_grps)
  !| Input matrix (flattened, n_gene * sum(group_c))
  real(c_double), intent(in), target :: input_matrix(n_gene * sum(group_c))
  !| Output matrix (flattened, n_gene * n_grps)
  real(c_double), intent(out) :: output_matrix(n_gene * n_grps)

  call calc_tiss_avg(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix)
end subroutine calc_tiss_avg_c

!> R/Fortran wrapper for fold change calculation.
!| Provides an interface for R (.Fortran) and Fortran code to call the fold change calculation routine.
!| Computes log2 fold changes between condition and control columns for all genes. Arguments match R's .Fortran calling convention.
subroutine calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix)
  use tox_normalization
  !| Number of genes (rows)
  integer(int32), intent(in) :: n_genes
  !| Number of columns in the input matrix
  integer(int32), intent(in) :: n_cols
  !| Number of control-condition pairs
  integer(int32), intent(in) :: n_pairs
  !| Control column indices (length n_pairs)
  integer(int32), intent(in) :: control_cols(n_pairs)
  !| Condition column indices (length n_pairs)
  integer(int32), intent(in) :: cond_cols(n_pairs)
  !| Input matrix (flattened, n_genes * n_cols)
  real(real64), intent(in) :: i_matrix(n_genes * n_cols)
  !| Output matrix for fold changes (flattened, n_genes * n_pairs)
  real(real64), intent(out) :: o_matrix(n_genes * n_pairs)

  call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix)
end subroutine calc_fchange_r

!> C/Python wrapper for fold change calculation.
!| Provides a C/Python-compatible interface to the fold change calculation routine.
!| Suitable for use with ctypes. Computes log2 fold changes between condition and control columns for all genes.
subroutine calc_fchange_c(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix) bind(C, name="calc_fchange_c")
  use iso_c_binding, only : c_int, c_double
  use tox_normalization
  !| Number of genes (rows)
  integer(c_int), intent(in), value :: n_genes
  !| Number of columns in the input matrix
  integer(c_int), intent(in), value :: n_cols
  !| Number of control-condition pairs
  integer(c_int), intent(in), value :: n_pairs
  !| Control column indices (length n_pairs)
  integer(c_int), intent(in), target :: control_cols(n_pairs)
  !| Condition column indices (length n_pairs)
  integer(c_int), intent(in), target :: cond_cols(n_pairs)
  !| Input matrix (flattened, n_genes * n_cols)
  real(c_double), intent(in), target :: i_matrix(n_genes * n_cols)
  !| Output matrix for fold changes (flattened, n_genes * n_pairs)
  real(c_double), intent(out) :: o_matrix(n_genes * n_pairs)

  call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix)
end subroutine calc_fchange_c