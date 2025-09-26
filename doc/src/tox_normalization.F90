!> Module with normalization routines for tensor omics.
module tox_normalization
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT, ERR_INVALID_INPUT, is_err
contains

  !> Complete normalization pipeline for gene expression data.
  !! Performs: std dev normalization, quantile normalization, replicate averaging, log2(x+1) transformation.
  !! Final result is in buf_log. If fold change is needed, call calc_fchange separately.
  pure subroutine normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, temp_col, rank_means, perm, stack_left, stack_right, max_stack, group_s, group_c, n_grps, ierr)

    !| Number of genes (rows)
    integer(int32), intent(in) :: n_genes
    !| Number of tissues (columns)
    integer(int32), intent(in) :: n_tissues
    !| Number of replicate groups
    integer(int32), intent(in) :: n_grps
    !| Stack size for quicksort
    integer(int32), intent(in) :: max_stack
    !| Flattened input matrix (n_genes * n_tissues), column-major
    real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
    !| Buffer for std dev normalization (n_genes * n_tissues)
    real(real64), intent(out) :: buf_stddev(n_genes * n_tissues)
    !| Buffer for quantile normalization (n_genes * n_tissues)
    real(real64), intent(out) :: buf_quant(n_genes * n_tissues)
    !| Buffer for replicate averaging (n_genes * n_grps)
    real(real64), intent(out) :: buf_avg(n_genes * n_grps)
    !| Buffer for log2(x+1) transformation (n_genes * n_grps)
    real(real64), intent(out) :: buf_log(n_genes * n_grps)
    !| Temporary column vector for sorting (n_genes)
    real(real64), intent(out) :: temp_col(n_genes)
    !| Buffer for rank means (n_genes)
    real(real64), intent(out) :: rank_means(n_genes)
    !| Permutation vector for sorting (n_genes)
    integer(int32), intent(out) :: perm(n_genes)
    !| Left stack for quicksort (max_stack)
    integer(int32), intent(out) :: stack_left(max_stack)
    !| Right stack for quicksort (max_stack)
    integer(int32), intent(out) :: stack_right(max_stack)
    !| Start column index for each replicate group (n_grps)
    integer(int32), intent(in) :: group_s(n_grps)
    !| Number of columns per replicate group (n_grps)
    integer(int32), intent(in) :: group_c(n_grps)
    !| Error code
    integer(int32), intent(out) :: ierr

    ! Error handling
    call set_ok(ierr)
    if (n_genes <= 0 .or. n_tissues <= 0 .or. n_grps <= 0 .or. max_stack <= 0) then
      call set_err(ierr, ERR_EMPTY_INPUT)
      return
    end if

    ! Step 1: Normalize per-gene by std dev
    call normalize_by_std_dev(n_genes, n_tissues, input_matrix, buf_stddev, ierr)
    if (is_err(ierr)) return

    ! Step 2: Quantile normalization
    call quantile_normalization(n_genes, n_tissues, buf_stddev, buf_quant, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    if (is_err(ierr)) return

    ! Step 3: Average replicates
    call calc_tiss_avg(n_genes, n_grps, group_s, group_c, buf_quant, buf_avg, ierr)
    if (is_err(ierr)) return

    ! Step 4: Log2(x+1) transformation
    call log2_transformation(n_genes, n_tissues, buf_avg, buf_log, ierr)
    if (is_err(ierr)) return

  end subroutine normalization_pipeline

  !> Normalizes each gene's expression vector using `sqrt(mean(x^2))`
  !| across tissues (not classical standard deviation).
  pure subroutine normalize_by_std_dev(n_genes, n_tissues, input_matrix, output_matrix, ierr)
      implicit none

      !| Number of genes (rows)
      integer(int32), intent(in) :: n_genes
      !| Number of tissues (columns)
      integer(int32), intent(in) :: n_tissues
      !| Flattened input matrix of gene expression values (column-major)
      real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
      !| Output normalized matrix (same shape as input)
      real(real64), intent(out) :: output_matrix(n_genes * n_tissues)
      !| Error code
      integer(int32), intent(out) :: ierr

      ! Local variables
      integer(int32) :: i_gene, i_tissue
      real(real64) :: std_dev, temp_sum

      ! Error handling
      call set_ok(ierr)
      if (n_genes <= 0 .or. n_tissues <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)
        return
      end if

      ! Loop over each gene
      do i_gene = 1, n_genes
          temp_sum = 0.0d0
          do i_tissue = 1, n_tissues
              temp_sum = temp_sum + input_matrix((i_tissue-1)*n_genes + i_gene)**2
          end do

          std_dev = sqrt(temp_sum / dble(n_tissues))
          if (std_dev == 0.0d0) std_dev = 1.0d0

          do i_tissue = 1, n_tissues
              output_matrix((i_tissue-1)*n_genes + i_gene) = input_matrix((i_tissue-1)*n_genes + i_gene) / std_dev
          end do
      end do

  end subroutine normalize_by_std_dev

  !> Quantile normalization of a gene expression matrix (F42-compliant).
  !| Computes average expression per rank across tissues.

  pure subroutine quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix, &
                                      temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
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
    real(real64), intent(out) :: temp_col(n_genes)
    !| Preallocated vector to store rank means (size n_genes)
    real(real64), intent(out) :: rank_means(n_genes)
    !| Permutation vector (size n_genes)
    integer(int32), intent(out) :: perm(n_genes)
    !| Stack size passed from R
    integer(int32), intent(in) :: max_stack
    !| Manual quicksort stack (≥ log2(n_genes) + 10)
    integer(int32), intent(out) :: stack_left(max_stack)
    !| Manual quicksort stack (same size as stack_left)
    integer(int32), intent(out) :: stack_right(max_stack)
    !| Error code
    integer(int32), intent(out) :: ierr

    ! Locals
    integer(int32) :: i_gene, i_tissue

    ! Error handling
    call set_ok(ierr)
    if (n_genes <= 0 .or. n_tissues <= 0 .or. max_stack <= 0) then
      call set_err(ierr, ERR_EMPTY_INPUT)
      return
    end if

    ! Initialize rank means
    rank_means = 0.0d0

    ! === First pass: accumulate values by rank across tissues ===
    do i_tissue = 1, n_tissues
        ! Prepare current column and initialize permutation
        do i_gene = 1, n_genes
            temp_col(i_gene) = input_matrix((i_tissue - 1) * n_genes + i_gene)
            perm(i_gene) = i_gene
        end do

        ! Sort current column with index tracking

        call sort_array(temp_col, perm, stack_left, stack_right)

        ! Accumulate values for each rank
        do i_gene = 1, n_genes
            rank_means(i_gene) = rank_means(i_gene) + temp_col(perm(i_gene))
        end do
    end do

    ! Average the rank values
    do i_gene = 1, n_genes
        rank_means(i_gene) = rank_means(i_gene) / dble(n_tissues)
    end do

    ! === Second pass: assign averaged values by rank ===
    do i_tissue = 1, n_tissues
        ! Prepare column and reset permutation
        do i_gene = 1, n_genes
            temp_col(i_gene) = input_matrix((i_tissue - 1) * n_genes + i_gene)
            perm(i_gene) = i_gene
        end do

        call sort_array(temp_col, perm, stack_left, stack_right)


        do i_gene = 1, n_genes
            output_matrix((i_tissue - 1) * n_genes + perm(i_gene)) = rank_means(i_gene)
        end do
    end do
  end subroutine quantile_normalization

  !> Apply `log2(x + 1)` transformation to each element of the input matrix.
  !| This subroutine performs element-wise `log2(x + 1)` transformation on a
  !| matrix flattened in column-major order. The `log2` is computed via:
  !| `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
  !| non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).

  pure subroutine log2_transformation(n_genes, n_tissues, input_matrix, output_matrix, ierr)
      implicit none

      !| Number of genes (rows)
      integer(int32), intent(in) :: n_genes
      !| Number of tissues (columns)
      integer(int32), intent(in) :: n_tissues
      !| Flattened input matrix (size: n_genes * n_tissues)
      real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
      !| Output matrix (same size as input)
      real(real64), intent(out) :: output_matrix(n_genes * n_tissues)
      !| Error code
      integer(int32), intent(out) :: ierr
      ! Locals
      integer(int32) :: i_elem
      real(real64), parameter :: LOG2 = log(2.0d0)

      ! Error handling
      call set_ok(ierr)
      if (n_genes <= 0 .or. n_tissues <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)
        return
      end if

      ! Loop through all elements in the flattened input matrix
      do i_elem = 1, n_genes * n_tissues
          ! Apply the log2(x + 1) transformation
          output_matrix(i_elem) = log(input_matrix(i_elem) + 1.0d0) / LOG2
      end do
  end subroutine log2_transformation


  !> Calculate tissue averages by averaging replicates within each group.
  !| For each group of tissue replicates, this subroutine computes the average
  !| expression per gene. The input matrix is column-major, flattened as a 1D array.

  pure subroutine calc_tiss_avg(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
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
      !| Error code
      integer(int32), intent(out) :: ierr

      ! === Local variables ===
      integer(int32) :: i_gene, i_group, i_col, j_col
      real(real64) :: sum_val
      integer(int32) :: start_idx, count_cols

      ! Error handling
      call set_ok(ierr)
      if (n_gene <= 0 .or. n_grps <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)
        return
      end if

      ! === Loop over each group ===
      do i_group = 1, n_grps
          start_idx = group_s(i_group)
          count_cols = group_c(i_group)

          do i_gene = 1, n_gene
              sum_val = 0.0d0

              do j_col = 0, count_cols - 1
                  i_col = start_idx + j_col
                  sum_val = sum_val + input_matrix((i_col - 1) * n_gene + i_gene)
              end do

              output_matrix((i_group - 1) * n_gene + i_gene) = sum_val / dble(count_cols)
          end do
      end do
  end subroutine calc_tiss_avg


  !> Calculate `log2 fold changes` between condition and control columns.
  !| For each control-condition pair, this subroutine computes the `log2 fold change`
  !| by subtracting the expression value in the control column from the corresponding
  !| value in the condition column, for all genes.
  !|
  !| The input matrix must be column-major and flattened as a 1D array.

  pure subroutine calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)

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
      !| Error code
      integer(int32), intent(out) :: ierr

      ! === Locals ===
      integer(int32) :: i_gene, i_pair
      integer(int32) :: control_col, cond_col

      ! Error handling
      call set_ok(ierr)
      if (n_genes <= 0 .or. n_cols <= 0 .or. n_pairs <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)
        return
      end if

      ! === Loop over each pair ===
      do i_pair = 1, n_pairs
          control_col = control_cols(i_pair)
          cond_col    = cond_cols(i_pair)

          do i_gene = 1, n_genes
              o_matrix((i_pair - 1) * n_genes + i_gene) = &
                  i_matrix((cond_col - 1) * n_genes + i_gene) - &
                  i_matrix((control_col - 1) * n_genes + i_gene)
          end do
      end do

  end subroutine calc_fchange


end module tox_normalization


!> R/Fortran wrapper for normalization by standard deviation.
!| Provides an interface for R (.Fortran) and Fortran code to call the normalization routine.
subroutine normalize_by_std_dev_r(n_genes, n_tissues, input_matrix, output_matrix, ierr)
  use tox_normalization
  !| Number of genes (rows)
  integer(int32), intent(in) :: n_genes
  !| Number of tissues (columns)
  integer(int32), intent(in) :: n_tissues
  !| Input matrix (n_genes x n_tissues)
  real(real64), intent(in)  :: input_matrix(n_genes, n_tissues)
  !| Output normalized matrix (same shape as input)
  real(real64), intent(out) :: output_matrix(n_genes, n_tissues)
  !| Error code
  integer(int32), intent(out) :: ierr
  
  call normalize_by_std_dev(n_genes, n_tissues, input_matrix, output_matrix, ierr)
  
end subroutine normalize_by_std_dev_r

!> C/Python wrapper for normalization by standard deviation.
!| Provides a C/Python-compatible interface to the normalization routine.
subroutine normalize_by_std_dev_c(n_genes, n_tissues, input_matrix, output_matrix, ierr) bind(C, name="normalize_by_std_dev_c")
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
  !| Error code
  integer(c_int), intent(out) :: ierr

  call normalize_by_std_dev(n_genes, n_tissues, input_matrix, output_matrix, ierr)

end subroutine normalize_by_std_dev_c


!> R/Fortran wrapper for quantile normalization.
!| Provides an interface for R (.Fortran) and Fortran code to call the quantile normalization routine.
subroutine quantile_normalization_r(n_genes, n_tissues, input_matrix, output_matrix, &
                                        temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
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
  real(real64), intent(out) :: temp_col(n_genes)
  !| Preallocated vector to store rank means (size n_genes)
  real(real64), intent(out) :: rank_means(n_genes)
  !| Permutation vector (size n_genes)
  integer(int32), intent(out) :: perm(n_genes)
  !| Manual quicksort stack (size max_stack)
  integer(int32), intent(out) :: stack_left(max_stack)
  !| Manual quicksort stack (size max_stack)
  integer(int32), intent(out) :: stack_right(max_stack)
  !| Error code
  integer(int32), intent(out) :: ierr

  call quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix, &
                            temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
end subroutine quantile_normalization_r

!> C/Python wrapper for quantile normalization.
!| Provides a C/Python-compatible interface to the quantile normalization routine.
subroutine quantile_normalization_c(n_genes, n_tissues, input_matrix, output_matrix, &
                                    temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr) &
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
  real(c_double), intent(out), target :: temp_col(n_genes)
  !| Preallocated vector to store rank means (size n_genes)
  real(c_double), intent(out), target :: rank_means(n_genes)
  !| Permutation vector (size n_genes)
  integer(c_int), intent(out), target :: perm(n_genes)
  !| Manual quicksort stack (size max_stack)
  integer(c_int), intent(out), target :: stack_left(max_stack)
  !| Manual quicksort stack (size max_stack)
  integer(c_int), intent(out), target :: stack_right(max_stack)
  !| Error code
  integer(c_int), intent(out) :: ierr

  call quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix, &
                            temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
end subroutine quantile_normalization_c

!> R/Fortran wrapper for log2 transformation.
!| Provides an interface for R (.Fortran) and Fortran code to call the log2 transformation routine.
!| Applies log2(x+1) to each element of the input matrix. Arguments match R's .Fortran calling convention and expect flat arrays.
subroutine log2_transformation_r(n_genes, n_tissues, input_matrix, output_matrix, ierr)
  use tox_normalization
  !| Number of genes (rows)
  integer(int32), intent(in) :: n_genes
  !| Number of tissues (columns)
  integer(int32), intent(in) :: n_tissues
  !| Input matrix (flattened, n_genes * n_tissues)
  real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
  !| Output matrix (flattened, same shape as input)
  real(real64), intent(out) :: output_matrix(n_genes * n_tissues)
  !| Error code
  integer(int32), intent(out) :: ierr

  call log2_transformation(n_genes, n_tissues, input_matrix, output_matrix, ierr)
  
end subroutine log2_transformation_r

!> C/Python wrapper for log2 transformation.
!| Provides a C/Python-compatible interface to the log2 transformation routine.
!| Expects flat arrays, matching C calling conventions. Suitable for use with ctypes.
!| Applies log2(x+1) to each element of the input matrix.
subroutine log2_transformation_c(n_genes, n_tissues, input_matrix, output_matrix, ierr) bind(C, name="log2_transformation_c")
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
  !| Error code
  integer(c_int), intent(out) :: ierr

  call log2_transformation(n_genes, n_tissues, input_matrix, output_matrix, ierr)
end subroutine log2_transformation_c

!> R/Fortran wrapper for tissue average calculation.
!| Provides an interface for R (.Fortran) and Fortran code to call the tissue average calculation routine.
!| Computes average expression per gene for each group of tissue replicates. Arguments match R's .Fortran calling convention.
subroutine calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
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
  !| Error code
  integer(int32), intent(out) :: ierr
  
  call calc_tiss_avg(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
  
end subroutine calc_tiss_avg_r

!> C/Python wrapper for tissue average calculation.
!| Provides a C/Python-compatible interface to the tissue average calculation routine.
!| Suitable for use with ctypes. Computes average expression per gene for each group of tissue replicates.
subroutine calc_tiss_avg_c(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr) bind(C, name="calc_tiss_avg_c")
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
  !| Error code
  integer(c_int), intent(out) :: ierr

  call calc_tiss_avg(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
end subroutine calc_tiss_avg_c

!> R/Fortran wrapper for fold change calculation.
!| Provides an interface for R (.Fortran) and Fortran code to call the fold change calculation routine.
!| Computes log2 fold changes between condition and control columns for all genes. Arguments match R's .Fortran calling convention.
subroutine calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
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
  !| Error code
  integer(int32), intent(out) :: ierr

  call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
end subroutine calc_fchange_r

!> C/Python wrapper for fold change calculation.
!| Provides a C/Python-compatible interface to the fold change calculation routine.
!| Suitable for use with ctypes. Computes log2 fold changes between condition and control columns for all genes.
subroutine calc_fchange_c(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr) bind(C, name="calc_fchange_c")
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
  !| Error code
  integer(c_int), intent(out) :: ierr

  call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
end subroutine calc_fchange_c

!> R/Fortran wrapper for normalization pipeline.
!| Provides an interface for R (.Fortran) and Fortran code to call the normalization pipeline routine.
!| Performs: std dev normalization, quantile normalization, replicate averaging, log2(x+1) transformation.
!| Final result is in buf_log. If fold change is needed, call calc_fchange separately.

subroutine normalization_pipeline_r(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, temp_col, rank_means, perm, stack_left, stack_right, max_stack, group_s, group_c, n_grps, ierr)
  use tox_normalization
  !| Number of genes (rows)
  integer(int32), intent(in) :: n_genes
  !| Number of tissues (columns)
  integer(int32), intent(in) :: n_tissues
  !| Number of replicate groups
  integer(int32), intent(in) :: n_grps
  !| Flattened input matrix (n_genes * n_tissues), column-major
  real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
  !| Buffer for std dev normalization (n_genes * n_tissues)
  real(real64), intent(out) :: buf_stddev(n_genes * n_tissues)
  !| Buffer for quantile normalization (n_genes * n_tissues)
  real(real64), intent(out) :: buf_quant(n_genes * n_tissues)
  !| Buffer for replicate averaging (n_genes * n_grps)
  real(real64), intent(out) :: buf_avg(n_genes * n_grps)
  !| Buffer for log2(x+1) transformation (n_genes * n_grps)
  real(real64), intent(out) :: buf_log(n_genes * n_grps)
  !| Temporary column vector for sorting (n_genes)
  real(real64), intent(out) :: temp_col(n_genes)
  !| Buffer for rank means (n_genes)
  real(real64), intent(out) :: rank_means(n_genes)
  !| Permutation vector for sorting (n_genes)
  integer(int32), intent(out) :: perm(n_genes)
  !| Stack size for quicksort
  integer(int32), intent(in) :: max_stack
  !| Left stack for quicksort (max_stack)
  integer(int32), intent(out) :: stack_left(max_stack)
  !| Right stack for quicksort (max_stack)
  integer(int32), intent(out) :: stack_right(max_stack)
  !| Start column index for each replicate group (n_grps)
  integer(int32), intent(in) :: group_s(n_grps)
  !| Number of columns per replicate group (n_grps)
  integer(int32), intent(in) :: group_c(n_grps)
  !| Error code
  integer(int32), intent(out) :: ierr

  call normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, temp_col, rank_means, perm, stack_left, stack_right, max_stack, group_s, group_c, n_grps, ierr)

end subroutine normalization_pipeline_r

!> C/Python wrapper for normalization pipeline.
!| Provides a C/Python-compatible interface to the normalization pipeline routine.
!| Suitable for use with ctypes. Performs: std dev normalization, quantile normalization, replicate averaging, log2(x+1) transformation.
!| Final result is in buf_log. If fold change is needed, call calc_fchange separately.

subroutine normalization_pipeline_c(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, temp_col, rank_means, perm, stack_left, stack_right, max_stack, group_s, group_c, n_grps, ierr) bind(C, name="normalization_pipeline_c")
  use iso_c_binding, only : c_int, c_double, c_bool
  use tox_normalization
  !| Number of genes (rows)
  integer(c_int), intent(in), value :: n_genes
  !| Number of tissues (columns)
  integer(c_int), intent(in), value :: n_tissues
  !| Number of replicate groups
  integer(c_int), intent(in), value :: n_grps
  !| Flattened input matrix (n_genes * n_tissues), column-major
  real(c_double), intent(in), target :: input_matrix(n_genes * n_tissues)
  !| Buffer for std dev normalization (n_genes * n_tissues)
  real(c_double), intent(out), target :: buf_stddev(n_genes * n_tissues)
  !| Buffer for quantile normalization (n_genes * n_tissues)
  real(c_double), intent(out), target :: buf_quant(n_genes * n_tissues)
  !| Buffer for replicate averaging (n_genes * n_grps)
  real(c_double), intent(out), target :: buf_avg(n_genes * n_grps)
  !| Buffer for log2(x+1) transformation (n_genes * n_grps)
  real(c_double), intent(out), target :: buf_log(n_genes * n_grps)
  !| Temporary column vector for sorting (n_genes)
  real(c_double), intent(out), target :: temp_col(n_genes)
  !| Buffer for rank means (n_genes)
  real(c_double), intent(out), target :: rank_means(n_genes)
  !| Permutation vector for sorting (n_genes)
  integer(c_int), intent(out), target :: perm(n_genes)
  !| Stack size for quicksort
  integer(c_int), intent(in), value :: max_stack
  !| Left stack for quicksort (max_stack)
  integer(c_int), intent(out), target :: stack_left(max_stack)
  !| Right stack for quicksort (max_stack)
  integer(c_int), intent(out), target :: stack_right(max_stack)
  !| Start column index for each replicate group (n_grps)
  integer(c_int), intent(in), target :: group_s(n_grps)
  !| Number of columns per replicate group (n_grps)
  integer(c_int), intent(in), target :: group_c(n_grps)
  !| Error code
  integer(c_int), intent(out) :: ierr

  call normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, temp_col, rank_means, perm, stack_left, stack_right, max_stack, group_s, group_c, n_grps, ierr)
end subroutine normalization_pipeline_c